//! Native GERG-2008 homogeneous Helmholtz evaluator.
//! Parameters derived from NIST AGA8 2.01; see data/source-extracts/gerg/LICENSE.
use serde::Deserialize;
use std::sync::OnceLock;
use thermo_types::{
    CalculationError as Error, ComponentId, Composition, CompositionBasis, InputError,
    NormalizationPolicy, Temperature, ValidityStatus,
};

/// GERG's formulation constant; deliberately distinct from the shared SI constant.
pub const GERG_GAS_CONSTANT: f64 = 8.314472;
pub const GERG_MODEL_VERSION: &str = "gerg-2008-nist-2.01-native-v1";
pub const GERG_DATA_REVISION: &str = "3bdb9ab8ff317c618b0b59d1b704c2c86ddc5fce";
/// Native GERG ideal reference, incompatible with cubic NASA7 anchored enthalpies.
pub const GERG_CALORIC_REFERENCE: &str =
    "GERG native ideal h=s=0 at 298.15 K, 101325 Pa (rounded NIST constants)";

#[derive(Deserialize)]
struct Dataset {
    components: Vec<String>,
    tables: Tables,
}

// Names retain upstream table notation to make the extraction auditable.
#[allow(non_snake_case)]
#[derive(Deserialize)]
struct Tables {
    Tc: Vec<f64>,
    MMiGERG: Vec<f64>,
    kpol: Vec<usize>,
    kexp: Vec<usize>,
    coik: Vec<Vec<i32>>,
    doik: Vec<Vec<i32>>,
    noik: Vec<Vec<f64>>,
    toik: Vec<Vec<f64>>,
    mNumb: Vec<Vec<i32>>,
    kpolij: Vec<usize>,
    kexpij: Vec<usize>,
    cijk: Vec<Vec<f64>>,
    eijk: Vec<Vec<f64>>,
    gijk: Vec<Vec<f64>>,
    nijk: Vec<Vec<f64>>,
    tijk: Vec<Vec<f64>>,
    dijk: Vec<Vec<i32>>,
    btij: Vec<Vec<f64>>,
    bvij: Vec<Vec<f64>>,
    gtij: Vec<Vec<f64>>,
    gvij: Vec<Vec<f64>>,
    fij: Vec<Vec<f64>>,
    th0i: Vec<Vec<f64>>,
    n0i: Vec<Vec<f64>>,
}

fn dataset() -> &'static Dataset {
    static DATA: OnceLock<Dataset> = OnceLock::new();
    DATA.get_or_init(|| {
        serde_json::from_str(include_str!("../../../data/gerg-2008.json"))
            .expect("checked embedded GERG parameters")
    })
}

/// Dimensionless Helmholtz derivatives scaled by powers of delta and tau.
/// Subscripts d/t denote differentiation with respect to delta/tau respectively.
#[derive(Debug, Clone, Copy, Default)]
pub struct GergDerivatives {
    pub ar: f64,
    pub ar_d: f64,
    pub ar_dd: f64,
    pub ar_t: f64,
    pub ar_tt: f64,
    pub ar_dt: f64,
    pub a0: f64,
    pub a0_t: f64,
    pub a0_tt: f64,
}

impl GergDerivatives {
    // A term n*delta^d*tau^t*exp(c*delta^2+e*delta+g-delta^q).
    // Logarithmic slopes give ordinary scaled second derivatives as L^2-L+L'.
    fn term(&mut self, value: f64, t: f64, slope: f64, curvature: f64) {
        self.ar += value;
        self.ar_d += value * slope;
        self.ar_dd += value * (slope * slope - slope + curvature);
        self.ar_t += value * t;
        self.ar_tt += value * t * (t - 1.0);
        self.ar_dt += value * t * slope;
    }
}

/// Homogeneous state in molar SI units. Mechanical/thermal checks are local;
/// neither phase stability nor equilibrium is established by this evaluator.
#[derive(Debug, Clone)]
pub struct GergState {
    /// Local admissibility never establishes physical phase validity.
    pub validity: ValidityStatus,
    pub temperature_k: f64,
    /// mol/m³
    pub molar_density: f64,
    pub mass_density: f64,
    pub pressure_pa: f64,
    pub compressibility: f64,
    /// Pa/(mol/m³), at fixed T and composition.
    pub dp_drho: f64,
    /// Pa/K at fixed molar density and composition.
    pub dp_dt: f64,
    pub internal_energy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub cv: f64,
    pub cp: f64,
    pub speed_of_sound: f64,
    pub helmholtz_energy: f64,
    pub gibbs_energy: f64,
    pub derivatives: GergDerivatives,
    pub model_version: &'static str,
    pub data_revision: &'static str,
    pub caloric_reference: &'static str,
}

/// Chemical-potential results in the original input component order, including
/// finite infinite-dilution coefficients for retained zero components.
#[derive(Debug, Clone)]
pub struct GergFugacityState {
    pub state: GergState,
    pub component_ids: Vec<String>,
    pub ln_fugacity_coefficients: Vec<f64>,
}

/// Immutable fixed-composition model; no mutable global calculation caches or FFI.
#[derive(Debug, Clone)]
pub struct GergModel {
    fractions: [f64; 22],
    component_indices: Vec<usize>,
    input_fractions: Vec<f64>,
    reducing_temperature: f64,
    reducing_density: f64,
    reducing_temperature_gradient: [f64; 22],
    reducing_volume_gradient: [f64; 22],
    molar_mass: f64,
}

impl GergModel {
    /// All 21 native GERG components. Methanol is unsupported, including at zero.
    /// Mass compositions must be converted explicitly with GERG molar masses.
    pub fn new(composition: &Composition) -> Result<Self, Error> {
        if composition.basis() != CompositionBasis::Mole {
            return Err(Error::UnsupportedCapability(
                "GERG requires mole fractions".into(),
            ));
        }
        let data = dataset();
        let p = &data.tables;
        let mut x = [0.0; 22];
        let mut component_indices = Vec::new();
        for (id, fraction) in composition.entries() {
            let i = data
                .components
                .iter()
                .position(|name| name == id.as_str())
                .ok_or_else(|| {
                    Error::MissingParameter(format!("GERG component {}", id.as_str()))
                })?;
            x[i + 1] = *fraction;
            component_indices.push(i + 1);
        }
        let mut tr = 0.0;
        let mut vr = 0.0;
        let mut mass = 0.0;
        let mut grad_t = [0.0; 22];
        let mut grad_v = [0.0; 22];
        for i in 1..=21 {
            mass += x[i] * p.MMiGERG[i] / 1000.0;
            for j in i..=21 {
                if x[i] + x[j] == 0.0 {
                    continue;
                }
                let factor = if i == j { 1.0 } else { 2.0 };
                let numerator = x[i] * x[j] * (x[i] + x[j]);
                for (b, g, value, gradient) in [
                    (p.bvij[i][j], p.gvij[i][j], &mut vr, &mut grad_v),
                    (p.btij[i][j], p.gtij[i][j], &mut tr, &mut grad_t),
                ] {
                    let denominator = b * x[i] + x[j];
                    let quotient = numerator / denominator;
                    *value += factor * g * quotient;
                    // i=j intentionally accumulates both partial contributions.
                    gradient[i] +=
                        factor * g * (x[j] * (2.0 * x[i] + x[j]) - b * quotient) / denominator;
                    gradient[j] +=
                        factor * g * (x[i] * (x[i] + 2.0 * x[j]) - quotient) / denominator;
                }
            }
        }
        Ok(Self {
            fractions: x,
            component_indices,
            input_fractions: composition.entries().iter().map(|(_, x)| *x).collect(),
            reducing_temperature: tr,
            reducing_density: 1000.0 / vr,
            reducing_temperature_gradient: grad_t,
            reducing_volume_gradient: grad_v,
            molar_mass: mass,
        })
    }

    pub fn supported_components() -> &'static [String] {
        &dataset().components
    }
    pub fn component_ids(&self) -> Vec<&str> {
        self.component_indices
            .iter()
            .map(|&i| dataset().components[i - 1].as_str())
            .collect()
    }
    pub fn fractions(&self) -> &[f64] {
        &self.input_fractions
    }
    /// Native critical temperature for one active species, including retained zeros.
    pub fn pure_critical_temperature(&self) -> Option<f64> {
        let active: Vec<_> = self
            .component_indices
            .iter()
            .filter(|&&i| self.fractions[i] > 0.0)
            .collect();
        if active.len() == 1 {
            Some(dataset().tables.Tc[*active[0]])
        } else {
            None
        }
    }
    /// Rebuild in the original component order, preserving inactive species.
    pub fn with_fractions(&self, fractions: &[f64]) -> Result<Self, Error> {
        if fractions.len() != self.component_indices.len() {
            return Err(InputError::ComponentMismatch.into());
        }
        let composition = Composition::new(
            CompositionBasis::Mole,
            self.component_ids()
                .iter()
                .zip(fractions)
                .map(|(id, x)| Ok((ComponentId::new(*id)?, *x)))
                .collect::<Result<_, InputError>>()?,
            NormalizationPolicy::RequireUnitSum,
        )?;
        Self::new(&composition)
    }
    pub fn reducing_temperature(&self) -> f64 {
        self.reducing_temperature
    }
    /// mol/m³
    pub fn reducing_density(&self) -> f64 {
        self.reducing_density
    }
    /// kg/mol, from GERG's own mass table.
    pub fn molar_mass(&self) -> f64 {
        self.molar_mass
    }

    /// Analytic residual composition derivatives and ln(phi), at specified T,rho.
    /// The original input order and zero components are retained. No equilibrium
    /// or phase-existence assessment is implied by finite chemical potentials.
    pub fn fugacity(&self, temperature: Temperature, rho: f64) -> Result<GergFugacityState, Error> {
        let state = self.evaluate(temperature, rho)?;
        let a = state.derivatives;
        let p = &dataset().tables;
        let delta = rho / self.reducing_density;
        let tau = self.reducing_temperature / temperature.value();
        let x = &self.fractions;
        let mut gradient = [0.0; 22];
        // Explicit mole-fraction derivatives at fixed reducing variables.
        // Terms for absent species are essential for infinite-dilution ln(phi).
        for i in 1..=21 {
            for k in 1..=p.kpol[i] + p.kexp[i] {
                let damping = if k > p.kpol[i] {
                    delta.powi(p.coik[i][k])
                } else {
                    0.0
                };
                gradient[i] += p.noik[i][k]
                    * delta.powi(p.doik[i][k])
                    * tau.powf(p.toik[i][k])
                    * (-damping).exp();
            }
            for j in i + 1..=21 {
                let m = p.mNumb[i][j];
                if m < 0 || (x[i] == 0.0 && x[j] == 0.0) {
                    continue;
                }
                let m = m as usize;
                for k in 1..=p.kpolij[m] + p.kexpij[m] {
                    let exponent = if k > p.kpolij[m] {
                        p.cijk[m][k] * delta * delta + p.eijk[m][k] * delta + p.gijk[m][k]
                    } else {
                        0.0
                    };
                    let term = p.fij[i][j]
                        * p.nijk[m][k]
                        * delta.powi(p.dijk[m][k])
                        * tau.powf(p.tijk[m][k])
                        * exponent.exp();
                    gradient[i] += x[j] * term;
                    gradient[j] += x[i] * term;
                }
            }
        }
        let vr = 1000.0 / self.reducing_density;
        for (i, derivative) in gradient.iter_mut().enumerate().skip(1) {
            *derivative += a.ar_d * self.reducing_volume_gradient[i] / vr
                + a.ar_t * self.reducing_temperature_gradient[i] / self.reducing_temperature;
        }
        let weighted: f64 = x.iter().zip(gradient).map(|(xi, gi)| xi * gi).sum();
        // d(n*ar)/dn_i at fixed T,V minus ln(Z) converts to the same-T,p
        // ideal reference used to define fugacity coefficients.
        let ln_phi: Vec<_> = self
            .component_indices
            .iter()
            .map(|&i| a.ar + a.ar_d + gradient[i] - weighted - state.compressibility.ln())
            .collect();
        if ln_phi.iter().any(|v| !v.is_finite()) {
            return Err(Error::NumericalFailure(
                "nonfinite GERG composition derivative".into(),
            ));
        }
        Ok(GergFugacityState {
            state,
            component_ids: self
                .component_indices
                .iter()
                .map(|&i| dataset().components[i - 1].clone())
                .collect(),
            ln_fugacity_coefficients: ln_phi,
        })
    }

    /// Evaluate at specified T and positive molar density, without root selection.
    /// Computational envelope: 200..=450 K, 1e-3..=30000 mol/m³.
    /// Derivatives may describe an unstable state; use `evaluate` for local checks.
    pub fn helmholtz(&self, temperature: Temperature, rho: f64) -> Result<GergDerivatives, Error> {
        let t = temperature.value();
        if !rho.is_finite() {
            return Err(InputError::NonFinite("GERG molar density").into());
        }
        if rho <= 0.0 {
            return Err(InputError::NotPositive("GERG molar density").into());
        }
        if !(200.0..=450.0).contains(&t) || !(1e-3..=30000.0).contains(&rho) {
            return Err(Error::OutsideValidityDomain(
                "GERG: 200..=450 K, 1e-3..=30000 mol/m³".into(),
            ));
        }
        let p = &dataset().tables;
        let delta = rho / self.reducing_density;
        let tau = self.reducing_temperature / t;
        let x = &self.fractions;
        let mut a = GergDerivatives::default();
        for i in 1..=21 {
            if x[i] == 0.0 {
                continue;
            }
            for k in 1..=p.kpol[i] + p.kexp[i] {
                let d = p.doik[i][k] as f64;
                let power = p.toik[i][k];
                let q = if k > p.kpol[i] { p.coik[i][k] } else { 0 };
                let damping = if q == 0 { 0.0 } else { delta.powi(q) };
                let v = x[i] * p.noik[i][k] * delta.powf(d) * tau.powf(power) * (-damping).exp();
                a.term(
                    v,
                    power,
                    d - f64::from(q) * damping,
                    -f64::from(q * q) * damping,
                );
            }
            let n = &p.n0i[i];
            // Separate logarithms preserve positive subnormal trace fractions.
            a.a0 += x[i] * ((rho / 1000.0).ln() + x[i].ln() + n[1] + n[2] / t - n[3] * t.ln());
            a.a0_t += x[i] * (n[3] + n[2] / t);
            a.a0_tt -= x[i] * n[3];
            for (j, &coefficient) in n.iter().enumerate().take(8).skip(4) {
                let theta = p.th0i[i][j];
                if theta == 0.0 {
                    continue;
                }
                let v = theta / t;
                if j == 4 || j == 6 {
                    a.a0 += x[i] * coefficient * v.sinh().ln();
                    a.a0_t += x[i] * coefficient * v / v.tanh();
                    a.a0_tt -= x[i] * coefficient * (v / v.sinh()).powi(2);
                } else {
                    a.a0 -= x[i] * coefficient * v.cosh().ln();
                    a.a0_t -= x[i] * coefficient * v * v.tanh();
                    a.a0_tt -= x[i] * coefficient * (v / v.cosh()).powi(2);
                }
            }
            for j in i + 1..=21 {
                if x[j] == 0.0 || p.mNumb[i][j] < 0 {
                    continue;
                }
                let m = p.mNumb[i][j] as usize;
                for k in 1..=p.kpolij[m] + p.kexpij[m] {
                    let d = p.dijk[m][k] as f64;
                    let power = p.tijk[m][k];
                    let (c, e, g) = if k > p.kpolij[m] {
                        (
                            p.cijk[m][k] * delta * delta,
                            p.eijk[m][k] * delta,
                            p.gijk[m][k],
                        )
                    } else {
                        (0.0, 0.0, 0.0)
                    };
                    let v = x[i]
                        * x[j]
                        * p.fij[i][j]
                        * p.nijk[m][k]
                        * delta.powf(d)
                        * tau.powf(power)
                        * (c + e + g).exp();
                    a.term(v, power, d + 2.0 * c + e, 4.0 * c + e);
                }
            }
        }
        if [
            a.ar, a.ar_d, a.ar_dd, a.ar_t, a.ar_tt, a.ar_dt, a.a0, a.a0_t, a.a0_tt,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(Error::NumericalFailure("GERG Helmholtz overflow".into()));
        }
        Ok(a)
    }

    /// Locally mechanically and thermally admissible homogeneous properties.
    /// Pressure must be positive and <=35 MPa. No phase/stability claim is made.
    pub fn evaluate(&self, temperature: Temperature, rho: f64) -> Result<GergState, Error> {
        let a = self.helmholtz(temperature, rho)?;
        let t = temperature.value();
        let r = GERG_GAS_CONSTANT;
        let z = 1.0 + a.ar_d;
        let pressure = rho * r * t * z;
        let pd = r * t * (1.0 + 2.0 * a.ar_d + a.ar_dd);
        let pt = rho * r * (1.0 + a.ar_d - a.ar_dt);
        let cv = -r * (a.a0_tt + a.ar_tt);
        if pd <= 0.0 || cv <= 0.0 {
            return Err(Error::OutsideValidityDomain(
                "GERG locally unstable homogeneous state".into(),
            ));
        }
        if !(0.0..=35e6).contains(&pressure) || pressure == 0.0 {
            return Err(Error::OutsideValidityDomain(
                "GERG pressure must be positive and <=35 MPa".into(),
            ));
        }
        let cp = cv + t * (pt / rho).powi(2) / pd;
        let u = r * t * (a.a0_t + a.ar_t);
        let helmholtz = r * t * (a.a0 + a.ar);
        let sound = (cp / cv * pd / self.molar_mass).sqrt();
        if [cp, u, helmholtz, sound, pd, pt]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(Error::NumericalFailure("GERG property overflow".into()));
        }
        Ok(GergState {
            validity: ValidityStatus::Unassessed,
            temperature_k: t,
            molar_density: rho,
            mass_density: rho * self.molar_mass,
            pressure_pa: pressure,
            compressibility: z,
            dp_drho: pd,
            dp_dt: pt,
            internal_energy: u,
            enthalpy: u + pressure / rho,
            entropy: r * (a.a0_t + a.ar_t - a.a0 - a.ar),
            cv,
            cp,
            speed_of_sound: sound,
            helmholtz_energy: helmholtz,
            gibbs_energy: helmholtz + pressure / rho,
            derivatives: a,
            model_version: GERG_MODEL_VERSION,
            data_revision: GERG_DATA_REVISION,
            caloric_reference: GERG_CALORIC_REFERENCE,
        })
    }
}
