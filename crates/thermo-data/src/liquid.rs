//! Pinned pure-liquid correlations and directed NRTL binary records (A7).
//! Source transcription: `data/source-extracts/liquid-nrtl.provisional.json`.
use serde::Deserialize;
use thermo_types::{CalculationError, ComponentId, DataMetadata, Temperature};

pub const LIQUID_DATA_VERSION: &str = "2026-11-08.1";

#[derive(Debug, Clone, Deserialize)]
struct Ancillary {
    form: String,
    #[serde(rename = "reducing_temperature_K")]
    t_r: f64,
    reducing_value: f64,
    using_tau_r: bool,
    n: Vec<f64>,
    t: Vec<f64>,
    #[serde(rename = "temperature_range_K")]
    range: [f64; 2],
}

#[derive(Debug, Clone, Deserialize)]
struct PureRecord {
    id: String,
    cas: String,
    vapor_pressure: Ancillary,
    saturated_liquid_density: Ancillary,
}

#[derive(Debug, Clone, Copy, Deserialize, PartialEq)]
pub struct TauTerms {
    pub t0: f64,
    #[serde(rename = "t1_K")]
    pub t1: f64,
    pub t2: f64,
    #[serde(rename = "t3_per_K")]
    pub t3: f64,
}

#[derive(Debug, Clone, Copy, Deserialize, PartialEq)]
pub struct AlphaTerms {
    pub a0: f64,
    #[serde(rename = "a1_per_K")]
    pub a1: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct BinaryRecord {
    i: String,
    j: String,
    source: String,
    source_row: String,
    tau_ij: TauTerms,
    tau_ji: TauTerms,
    alpha_ij: AlphaTerms,
    #[serde(rename = "regression_temperature_K")]
    regression_temperature: [f64; 2],
    #[serde(rename = "regression_pressure_Pa")]
    regression_pressure: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct Extract {
    schema: String,
    version: String,
    #[serde(rename = "gas_constant_J_per_mol_K")]
    gas_constant: f64,
    pure_liquid: Vec<PureRecord>,
    nrtl_binary: Vec<BinaryRecord>,
}

/// ln(psat/Pa) and its first two temperature derivatives (1/K, 1/K²).
#[derive(Debug, Clone, Copy)]
pub struct VaporPressure {
    pub pressure: f64,
    pub ln_pressure: f64,
    pub d_ln_pressure_dt: f64,
    pub d2_ln_pressure_dt2: f64,
}

/// Pure-liquid correlations for one component; no EOS coupling is implied.
#[derive(Debug, Clone)]
pub struct PureLiquid(PureRecord);

impl PureLiquid {
    pub fn id(&self) -> &str {
        &self.0.id
    }
    pub fn cas(&self) -> &str {
        &self.0.cas
    }
    /// Declared vapor-pressure correlation range in K (inclusive).
    pub fn vapor_pressure_range(&self) -> [f64; 2] {
        self.0.vapor_pressure.range
    }
    pub fn liquid_density_range(&self) -> [f64; 2] {
        self.0.saturated_liquid_density.range
    }

    /// CoolProp pV ancillary: psat = p_r exp((T_r/T) Σ n_k θ^t_k), θ = 1 − T/T_r.
    /// Analytical derivatives. Outside the declared range: `OutsideValidityDomain`.
    pub fn vapor_pressure(
        &self,
        temperature: Temperature,
    ) -> Result<VaporPressure, CalculationError> {
        let a = &self.0.vapor_pressure;
        let t = temperature.value();
        check_range(self.id(), "vapor pressure", a.range, t)?;
        let theta = 1.0 - t / a.t_r;
        // S = Σ n θ^e; dθ/dT = −1/T_r. θ > 0 inside the range (Tmax < T_r).
        let (mut s, mut ds, mut d2s) = (0.0, 0.0, 0.0);
        for (n, e) in a.n.iter().zip(&a.t) {
            s += n * theta.powf(*e);
            ds += -n * e * theta.powf(e - 1.0) / a.t_r;
            d2s += n * e * (e - 1.0) * theta.powf(e - 2.0) / (a.t_r * a.t_r);
        }
        // ln p = ln p_r + (T_r/T) S
        let q = a.t_r / t;
        let dq = -a.t_r / (t * t);
        let d2q = 2.0 * a.t_r / (t * t * t);
        let ln_p = a.reducing_value.ln() + q * s;
        let d1 = dq * s + q * ds;
        let d2 = d2q * s + 2.0 * dq * ds + q * d2s;
        let pressure = ln_p.exp();
        if ![pressure, ln_p, d1, d2].iter().all(|v| v.is_finite()) || pressure <= 0.0 {
            return Err(CalculationError::NumericalFailure(format!(
                "{} vapor pressure not representable",
                self.id()
            )));
        }
        Ok(VaporPressure {
            pressure,
            ln_pressure: ln_p,
            d_ln_pressure_dt: d1,
            d2_ln_pressure_dt2: d2,
        })
    }

    /// CoolProp rhoLnoexp ancillary: ρ' = ρ_r (1 + Σ n_k θ^t_k), mol/m³, saturated liquid.
    pub fn saturated_liquid_molar_density(
        &self,
        temperature: Temperature,
    ) -> Result<f64, CalculationError> {
        let a = &self.0.saturated_liquid_density;
        let t = temperature.value();
        check_range(self.id(), "saturated liquid density", a.range, t)?;
        let theta = 1.0 - t / a.t_r;
        let rho = a.reducing_value
            * (1.0
                + a.n
                    .iter()
                    .zip(&a.t)
                    .map(|(n, e)| n * theta.powf(*e))
                    .sum::<f64>());
        if !rho.is_finite() || rho <= 0.0 {
            return Err(CalculationError::NumericalFailure(format!(
                "{} liquid density not positive",
                self.id()
            )));
        }
        Ok(rho)
    }
}

fn check_range(id: &str, what: &str, range: [f64; 2], t: f64) -> Result<(), CalculationError> {
    if !(range[0]..=range[1]).contains(&t) {
        return Err(CalculationError::OutsideValidityDomain(format!(
            "{id} {what} correlation: {}..={} K",
            range[0], range[1]
        )));
    }
    Ok(())
}

/// Directed NRTL record for an ordered pair (i, j):
/// τ_ij = t0 + t1/T + t2 ln(T/K) + t3 T and symmetric α_ij = a0 + a1 T.
#[derive(Debug, Clone)]
pub struct NrtlBinary {
    pub tau_ij: TauTerms,
    pub tau_ji: TauTerms,
    pub alpha: AlphaTerms,
    /// Temperatures (K) spanned by the regression data; outside is extrapolation.
    pub regression_temperature_k: [f64; 2],
    pub regression_pressure_pa: f64,
    pub source: String,
}

/// Immutable checked liquid/NRTL dataset. Pairs absent from the file are missing:
/// no zero-interaction default is ever supplied.
#[derive(Debug, Clone)]
pub struct LiquidDatabase {
    pure: Vec<PureLiquid>,
    binaries: Vec<BinaryRecord>,
    version: String,
}

impl LiquidDatabase {
    pub fn embedded() -> Result<Self, CalculationError> {
        let extract: Extract = serde_json::from_str(include_str!(
            "../../../data/source-extracts/liquid-nrtl.provisional.json"
        ))
        .map_err(|e| CalculationError::MissingParameter(e.to_string()))?;
        let bad =
            |what: &str| CalculationError::MissingParameter(format!("invalid embedded {what}"));
        if extract.schema != "a7-liquid-nrtl-source-extract-v1"
            || extract.version != LIQUID_DATA_VERSION
            || extract.gas_constant != thermo_types::MOLAR_GAS_CONSTANT
        {
            return Err(bad("liquid extract header"));
        }
        for (k, p) in extract.pure_liquid.iter().enumerate() {
            ComponentId::new(&p.id)?;
            if extract.pure_liquid[..k].iter().any(|o| o.id == p.id) {
                return Err(bad("duplicate pure liquid"));
            }
            for (a, form, tau) in [
                (&p.vapor_pressure, "pV", true),
                (&p.saturated_liquid_density, "rhoLnoexp", false),
            ] {
                let finite = a.n.iter().chain(&a.t).all(|v| v.is_finite());
                if a.form != form
                    || a.using_tau_r != tau
                    || a.n.len() != a.t.len()
                    || a.n.is_empty()
                    || !finite
                    || !(a.t_r > 0.0 && a.reducing_value > 0.0)
                    || !(0.0 < a.range[0] && a.range[0] < a.range[1] && a.range[1] < a.t_r)
                {
                    return Err(bad("pure liquid ancillary"));
                }
            }
        }
        for (k, b) in extract.nrtl_binary.iter().enumerate() {
            ComponentId::new(&b.i)?;
            ComponentId::new(&b.j)?;
            let terms = [b.tau_ij, b.tau_ji];
            let finite = terms
                .iter()
                .flat_map(|t| [t.t0, t.t1, t.t2, t.t3])
                .chain([b.alpha_ij.a0, b.alpha_ij.a1])
                .all(f64::is_finite);
            let duplicate = extract.nrtl_binary[..k]
                .iter()
                .any(|o| (o.i == b.i && o.j == b.j) || (o.i == b.j && o.j == b.i));
            if b.i == b.j
                || duplicate
                || !finite
                || b.regression_temperature[0] >= b.regression_temperature[1]
                || b.source.trim().is_empty()
            {
                return Err(bad("NRTL binary"));
            }
        }
        Ok(Self {
            pure: extract.pure_liquid.into_iter().map(PureLiquid).collect(),
            binaries: extract.nrtl_binary,
            version: extract.version,
        })
    }

    pub fn pure(&self, id: &str) -> Result<&PureLiquid, CalculationError> {
        self.pure.iter().find(|p| p.id() == id).ok_or_else(|| {
            CalculationError::MissingParameter(format!("pure-liquid correlations for {id}"))
        })
    }

    pub fn pure_liquids(&self) -> &[PureLiquid] {
        &self.pure
    }

    /// Directed record oriented as (i, j): a stored (j, i) record is transposed.
    pub fn nrtl_binary(&self, i: &str, j: &str) -> Result<NrtlBinary, CalculationError> {
        for b in &self.binaries {
            let oriented = if b.i == i && b.j == j {
                Some((b.tau_ij, b.tau_ji))
            } else if b.i == j && b.j == i {
                Some((b.tau_ji, b.tau_ij))
            } else {
                None
            };
            if let Some((tau_ij, tau_ji)) = oriented {
                return Ok(NrtlBinary {
                    tau_ij,
                    tau_ji,
                    alpha: b.alpha_ij,
                    regression_temperature_k: b.regression_temperature,
                    regression_pressure_pa: b.regression_pressure,
                    source: format!("{}: {}", b.source, b.source_row),
                });
            }
        }
        Err(CalculationError::MissingParameter(format!(
            "NRTL interaction for pair {i}/{j}"
        )))
    }

    pub fn metadata(&self) -> DataMetadata {
        DataMetadata {
            identifier: "coolprop-7.2-ancillaries-chemsep-dechema-nrtl".into(),
            version: self.version.clone(),
            source: "data/source-extracts/liquid-nrtl.provisional.json; upstream URLs, hashes and notices retained".into(),
        }
    }
}
