//! NRTL excess Gibbs energy with directed, temperature-dependent parameters.
//!
//! g^E/(RT) = Σ_i x_i (Σ_j x_j τ_ji G_ji) / (Σ_k x_k G_ki), G_ji = exp(−α_ji τ_ji),
//! τ_ij = t0_ij + t1_ij/T + t2_ij ln(T/K) + t3_ij T, α_ij = α_ji = a0_ij + a1_ij T.
//! τ_ij ≠ τ_ji in general (directed); every diagonal term is zero. First and second
//! temperature derivatives are propagated analytically (second-order jets).
use std::ops::{Add, Div, Mul, Neg, Sub};
use thermo_data::{AlphaTerms, LiquidDatabase, TauTerms};
use thermo_types::{
    COMPOSITION_SUM_TOLERANCE, CalculationError as Error, ComponentId, InputError,
    MOLAR_GAS_CONSTANT as R, ModelMetadata, ValidityStatus,
};

/// Value with first and second temperature derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Jet {
    v: f64,
    d: f64,
    dd: f64,
}

impl Jet {
    const fn constant(v: f64) -> Self {
        Self { v, d: 0.0, dd: 0.0 }
    }
    fn exp(self) -> Self {
        let e = self.v.exp();
        Self {
            v: e,
            d: e * self.d,
            dd: e * (self.dd + self.d * self.d),
        }
    }
}

impl Add for Jet {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: self.d + o.d,
            dd: self.dd + o.dd,
        }
    }
}
impl Sub for Jet {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d: self.d - o.d,
            dd: self.dd - o.dd,
        }
    }
}
impl Neg for Jet {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: -self.d,
            dd: -self.dd,
        }
    }
}
impl Mul for Jet {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d: self.d * o.v + self.v * o.d,
            dd: self.dd * o.v + 2.0 * self.d * o.d + self.v * o.dd,
        }
    }
}
impl Mul<f64> for Jet {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        Self {
            v: self.v * o,
            d: self.d * o,
            dd: self.dd * o,
        }
    }
}
impl Div for Jet {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let q = self.v / o.v;
        let d = (self.d - q * o.d) / o.v;
        let dd = (self.dd - 2.0 * d * o.d - q * o.dd) / o.v;
        Self { v: q, d, dd }
    }
}

type JetMatrix = Vec<Vec<Jet>>;

/// Complete, directed NRTL parameter set in component order.
#[derive(Debug, Clone, PartialEq)]
pub struct NrtlParameters {
    ids: Vec<ComponentId>,
    tau: Vec<Vec<TauTerms>>,
    alpha: Vec<Vec<AlphaTerms>>,
    /// Per off-diagonal pair: regression temperature interval (K), when known.
    regression: Vec<Vec<Option<[f64; 2]>>>,
    source: String,
}

const ZERO_TAU: TauTerms = TauTerms {
    t0: 0.0,
    t1: 0.0,
    t2: 0.0,
    t3: 0.0,
};
const ZERO_ALPHA: AlphaTerms = AlphaTerms { a0: 0.0, a1: 0.0 };

impl NrtlParameters {
    /// Explicit matrices: `tau[i][j]` holds τ_ij terms. Diagonals must be exactly zero;
    /// α must be symmetric. A nonblank source or explicit approximation note is required.
    pub fn new(
        ids: Vec<ComponentId>,
        tau: Vec<Vec<TauTerms>>,
        alpha: Vec<Vec<AlphaTerms>>,
        source: impl Into<String>,
    ) -> Result<Self, Error> {
        let n = ids.len();
        let regression = vec![vec![None; n]; n];
        Self::checked(ids, tau, alpha, regression, source.into())
    }

    fn checked(
        ids: Vec<ComponentId>,
        tau: Vec<Vec<TauTerms>>,
        alpha: Vec<Vec<AlphaTerms>>,
        regression: Vec<Vec<Option<[f64; 2]>>>,
        source: String,
    ) -> Result<Self, Error> {
        let n = ids.len();
        if n == 0 {
            return Err(InputError::EmptyComposition.into());
        }
        if (0..n).any(|i| ids[..i].contains(&ids[i])) {
            return Err(InputError::DuplicateComponent.into());
        }
        if tau.len() != n
            || alpha.len() != n
            || tau.iter().any(|r| r.len() != n)
            || alpha.iter().any(|r| r.len() != n)
        {
            return Err(Error::MissingParameter(
                "complete n×n NRTL τ and α matrices required".into(),
            ));
        }
        if source.trim().is_empty() {
            return Err(Error::MissingParameter(
                "NRTL parameter source or explicit approximation description required".into(),
            ));
        }
        for i in 0..n {
            for j in 0..n {
                let t = tau[i][j];
                let a = alpha[i][j];
                if ![t.t0, t.t1, t.t2, t.t3, a.a0, a.a1]
                    .iter()
                    .all(|v| v.is_finite())
                {
                    return Err(InputError::NonFinite("NRTL parameter").into());
                }
                if i == j && (t != ZERO_TAU || a != ZERO_ALPHA) {
                    return Err(Error::UnsupportedCapability(
                        "NRTL diagonal τ_ii and α_ii must be zero".into(),
                    ));
                }
                if a != alpha[j][i] {
                    return Err(Error::UnsupportedCapability(
                        "NRTL α_ij must equal α_ji".into(),
                    ));
                }
            }
        }
        Ok(Self {
            ids,
            tau,
            alpha,
            regression,
            source,
        })
    }

    /// All ordered pairs from the dataset. Any absent pair is `MissingParameter`;
    /// a zero interaction is never substituted.
    pub fn from_database(database: &LiquidDatabase, ids: &[ComponentId]) -> Result<Self, Error> {
        let n = ids.len();
        let mut tau = vec![vec![ZERO_TAU; n]; n];
        let mut alpha = vec![vec![ZERO_ALPHA; n]; n];
        let mut regression = vec![vec![None; n]; n];
        let mut sources = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let b = database.nrtl_binary(ids[i].as_str(), ids[j].as_str())?;
                tau[i][j] = b.tau_ij;
                tau[j][i] = b.tau_ji;
                alpha[i][j] = b.alpha;
                alpha[j][i] = b.alpha;
                regression[i][j] = Some(b.regression_temperature_k);
                regression[j][i] = Some(b.regression_temperature_k);
                sources.push(b.source);
            }
        }
        let source = if sources.is_empty() {
            "single component: no NRTL interactions required".into()
        } else {
            format!(
                "{} data {}: {}",
                database.metadata().identifier,
                database.metadata().version,
                sources.join("; ")
            )
        };
        Self::checked(ids.to_vec(), tau, alpha, regression, source)
    }

    pub fn component_ids(&self) -> &[ComponentId] {
        &self.ids
    }
    pub fn source(&self) -> &str {
        &self.source
    }
    /// τ_ij terms (directed).
    pub fn tau_terms(&self, i: usize, j: usize) -> TauTerms {
        self.tau[i][j]
    }
    pub fn alpha_terms(&self, i: usize, j: usize) -> AlphaTerms {
        self.alpha[i][j]
    }

    /// Same parameters for the components in `order` (indices into this set).
    pub fn permuted(&self, order: &[usize]) -> Result<Self, Error> {
        let n = self.ids.len();
        let mut seen = vec![false; n];
        if order.len() != n
            || order
                .iter()
                .any(|&k| k >= n || std::mem::replace(&mut seen[k], true))
        {
            return Err(InputError::ComponentMismatch.into());
        }
        fn pick<T: Copy>(m: &[Vec<T>], order: &[usize]) -> Vec<Vec<T>> {
            order
                .iter()
                .map(|&i| order.iter().map(|&j| m[i][j]).collect())
                .collect()
        }
        Ok(Self {
            ids: order.iter().map(|&i| self.ids[i].clone()).collect(),
            tau: pick(&self.tau, order),
            alpha: pick(&self.alpha, order),
            regression: pick(&self.regression, order),
            source: self.source.clone(),
        })
    }
}

/// Excess and activity properties at T and x (model order). Units: J/mol, J/(mol K).
#[derive(Debug, Clone)]
pub struct ExcessState {
    pub temperature_k: f64,
    pub ln_activity_coefficients: Vec<f64>,
    /// ∂ln γ_i/∂T at fixed x, 1/K.
    pub d_ln_activity_coefficients_dt: Vec<f64>,
    /// g^E/(RT), dimensionless.
    pub reduced_excess_gibbs: f64,
    pub excess_gibbs: f64,
    /// h^E = −R T² ∂(g^E/RT)/∂T.
    pub excess_enthalpy: f64,
    pub excess_entropy: f64,
    pub excess_heat_capacity: f64,
    /// Within the regression temperature interval of every active pair, extrapolated
    /// otherwise, unassessed when explicit parameters carry no interval.
    pub validity: ValidityStatus,
}

/// NRTL liquid excess model. Pressure independent; no volumetric excess.
#[derive(Debug, Clone)]
pub struct NrtlModel {
    parameters: NrtlParameters,
}

impl NrtlModel {
    pub fn new(parameters: NrtlParameters) -> Self {
        Self { parameters }
    }

    pub fn parameters(&self) -> &NrtlParameters {
        &self.parameters
    }

    pub fn metadata(&self) -> ModelMetadata {
        let p = &self.parameters;
        ModelMetadata {
            identifier: "NRTL-Renon-Prausnitz-1968-directed-T".into(),
            implementation_version: "a7-v1".into(),
            formulation: format!(
                "tau_ij = t0 + t1/T + t2 ln(T/K) + t3 T; alpha_ij = a0 + a1 T (symmetric); G_ij = exp(-alpha_ij tau_ij); R={R}; ordered components={:?}; tau(t0,t1,t2,t3)={:?}; alpha(a0,a1)={:?}; source: {}",
                p.ids.iter().map(ComponentId::as_str).collect::<Vec<_>>(),
                p.tau
                    .iter()
                    .map(|r| r
                        .iter()
                        .map(|t| [t.t0, t.t1, t.t2, t.t3])
                        .collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                p.alpha
                    .iter()
                    .map(|r| r.iter().map(|a| [a.a0, a.a1]).collect::<Vec<_>>())
                    .collect::<Vec<_>>(),
                p.source
            ),
        }
    }

    fn matrices(&self, t: f64) -> Result<(JetMatrix, JetMatrix), Error> {
        let n = self.parameters.ids.len();
        let lnt = Jet {
            v: t.ln(),
            d: 1.0 / t,
            dd: -1.0 / (t * t),
        };
        let inv = Jet {
            v: 1.0 / t,
            d: -1.0 / (t * t),
            dd: 2.0 / (t * t * t),
        };
        let tj = Jet {
            v: t,
            d: 1.0,
            dd: 0.0,
        };
        let mut tau = vec![vec![Jet::constant(0.0); n]; n];
        let mut g = vec![vec![Jet::constant(1.0); n]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let p = self.parameters.tau[i][j];
                let a = self.parameters.alpha[i][j];
                let tij = Jet::constant(p.t0) + inv * p.t1 + lnt * p.t2 + tj * p.t3;
                let aij = Jet::constant(a.a0) + tj * a.a1;
                if aij.v.is_nan() || aij.v <= 0.0 {
                    return Err(Error::OutsideValidityDomain(format!(
                        "NRTL alpha_{i}{j} = {} is not positive at {t} K",
                        aij.v
                    )));
                }
                tau[i][j] = tij;
                g[i][j] = (-(aij * tij)).exp();
            }
        }
        Ok((tau, g))
    }

    /// Activity coefficients and excess properties. `x` must be finite, nonnegative and
    /// sum to one within 1e-10; zero entries give infinite-dilution ln γ.
    pub fn excess(&self, temperature: f64, x: &[f64]) -> Result<ExcessState, Error> {
        let n = self.parameters.ids.len();
        if x.len() != n {
            return Err(InputError::ComponentMismatch.into());
        }
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(InputError::NotPositive("temperature").into());
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(InputError::NonFinite("mole fraction").into());
        }
        if x.iter().any(|v| *v < 0.0) {
            return Err(InputError::Negative("mole fraction").into());
        }
        if (x.iter().sum::<f64>() - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
            return Err(InputError::NotNormalized.into());
        }
        let t = temperature;
        let (tau, g) = self.matrices(t)?;
        // Column sums over the first index: S_i = Σ_k x_k G_ki, C_i = Σ_j x_j τ_ji G_ji.
        let mut s = vec![Jet::constant(0.0); n];
        let mut c = vec![Jet::constant(0.0); n];
        for i in 0..n {
            for k in 0..n {
                s[i] = s[i] + g[k][i] * x[k];
                c[i] = c[i] + tau[k][i] * g[k][i] * x[k];
            }
        }
        let mut ge = Jet::constant(0.0);
        let mut ln_gamma = Vec::with_capacity(n);
        for i in 0..n {
            let ratio_i = c[i] / s[i];
            ge = ge + ratio_i * x[i];
            let mut l = ratio_i;
            for j in 0..n {
                l = l + (g[i][j] * x[j] / s[j]) * (tau[i][j] - c[j] / s[j]);
            }
            ln_gamma.push(l);
        }
        let rt = R * t;
        let excess_gibbs = rt * ge.v;
        let excess_enthalpy = -R * t * t * ge.d;
        let excess_heat_capacity = -R * (2.0 * t * ge.d + t * t * ge.dd);
        let excess_entropy = (excess_enthalpy - excess_gibbs) / t;
        let values: Vec<f64> = ln_gamma.iter().map(|j| j.v).collect();
        let slopes: Vec<f64> = ln_gamma.iter().map(|j| j.d).collect();
        if ![
            excess_gibbs,
            excess_enthalpy,
            excess_heat_capacity,
            excess_entropy,
        ]
        .iter()
        .chain(&values)
        .chain(&slopes)
        .all(|v| v.is_finite())
        {
            return Err(Error::NumericalFailure("nonfinite NRTL result".into()));
        }
        Ok(ExcessState {
            temperature_k: t,
            ln_activity_coefficients: values,
            d_ln_activity_coefficients_dt: slopes,
            reduced_excess_gibbs: ge.v,
            excess_gibbs,
            excess_enthalpy,
            excess_entropy,
            excess_heat_capacity,
            validity: self.validity(t, x),
        })
    }

    /// Regression-interval status over interacting (both present) pairs. With no
    /// interacting pair g^E = 0 and no interval applies, so the status is `Unassessed`.
    fn validity(&self, t: f64, x: &[f64]) -> ValidityStatus {
        let n = x.len();
        let mut status = None;
        for i in 0..n {
            for j in i + 1..n {
                if x[i] <= 0.0 || x[j] <= 0.0 {
                    continue;
                }
                match self.parameters.regression[i][j] {
                    None => return ValidityStatus::Unassessed,
                    Some([lo, hi]) if !(lo..=hi).contains(&t) => {
                        status = Some(ValidityStatus::Extrapolated)
                    }
                    Some(_) => {
                        status.get_or_insert(ValidityStatus::WithinDeclaredDomain);
                    }
                }
            }
        }
        status.unwrap_or(ValidityStatus::Unassessed)
    }
}
