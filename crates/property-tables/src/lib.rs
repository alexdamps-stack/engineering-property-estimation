//! A11 tabulation experiment for one fixed mixture of the exact PR cubic model.
//!
//! The table stores the molar Gibbs energy `g(T, ln p)` and its exact first and
//! cross derivatives from the direct model at every node, and interpolates with a
//! tensor-product cubic Hermite surface. All returned properties derive from that
//! one surface (`s = -∂g/∂T`, `v = ∂g/∂p`, `h = g + T s`, `cp = -T ∂²g/∂T²`, …),
//! so the Gibbs–Helmholtz, Maxwell (`∂v/∂T = -∂s/∂p`) and `h = g + Ts` relations
//! hold exactly within the table instead of being violated by independent
//! property splines.
//!
//! A table exists only for a domain that passed qualification on a refined
//! sample of the rectangle: exactly one real cubic root, stability in the finite
//! TPD trial set, and a reduced bulk modulus above a bound that excludes the
//! critical region. Outside the domain, or for another composition, evaluation
//! follows [`OutsidePolicy`]: a clear [`TableError`] or a direct, phase-qualified
//! TP flash whose result is labelled [`Source::DirectModel`]. Results never
//! silently mix sources.
//!
//! This crate is an experiment with measured, documented scope (`docs/a11-status.md`);
//! it does not replace the direct model in the flowsheet or project format.

use std::fmt;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant, HomogeneousState, RootSelection};
use thermo_solvers::{
    FlashOptions, PhFlashOptions, StabilityOptions, assess_stability, ph_flash, tp_flash,
};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, InputError,
    MOLAR_GAS_CONSTANT as R, NormalizationPolicy, Pressure, Temperature,
};

/// Rectangle in temperature and pressure with node counts. Temperature nodes are
/// uniform in T, pressure nodes uniform in ln p.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TableDomain {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub minimum_pressure_pa: f64,
    pub maximum_pressure_pa: f64,
    pub temperature_nodes: usize,
    pub pressure_nodes: usize,
}

impl TableDomain {
    /// Inclusive bounds; non-finite inputs are never inside.
    pub fn contains(&self, temperature_k: f64, pressure_pa: f64) -> bool {
        (self.minimum_temperature_k..=self.maximum_temperature_k).contains(&temperature_k)
            && (self.minimum_pressure_pa..=self.maximum_pressure_pa).contains(&pressure_pa)
    }

    fn validate(&self) -> Result<(), TableError> {
        let ok = [
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.minimum_pressure_pa,
            self.maximum_pressure_pa,
        ]
        .iter()
        .all(|v| v.is_finite())
            && self.minimum_temperature_k < self.maximum_temperature_k
            && 0.0 < self.minimum_pressure_pa
            && self.minimum_pressure_pa < self.maximum_pressure_pa
            && self.temperature_nodes >= 2
            && self.pressure_nodes >= 2
            && self.temperature_nodes.saturating_mul(self.pressure_nodes) <= 1_000_000;
        if ok {
            Ok(())
        } else {
            Err(TableError::InvalidDomain(
                "ordered finite bounds, positive pressures and 2..=1e6 total nodes (≥2 per axis) required".into(),
            ))
        }
    }

    fn temperature_step(&self) -> f64 {
        (self.maximum_temperature_k - self.minimum_temperature_k)
            / (self.temperature_nodes - 1) as f64
    }

    fn log_pressure_step(&self) -> f64 {
        (self.maximum_pressure_pa / self.minimum_pressure_pa).ln()
            / (self.pressure_nodes - 1) as f64
    }

    fn temperature(&self, i: usize) -> f64 {
        if i + 1 == self.temperature_nodes {
            self.maximum_temperature_k
        } else {
            self.minimum_temperature_k + i as f64 * self.temperature_step()
        }
    }

    fn pressure(&self, j: usize) -> f64 {
        if j + 1 == self.pressure_nodes {
            self.maximum_pressure_pa
        } else {
            self.minimum_pressure_pa * (j as f64 * self.log_pressure_step()).exp()
        }
    }

    fn describe(&self) -> String {
        format!(
            "{}..={} K, {}..={} Pa",
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.minimum_pressure_pa,
            self.maximum_pressure_pa
        )
    }
}

/// Conditions a domain must satisfy at every qualification point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualificationCriteria {
    /// Sample points per cell edge (≥1): 1 checks nodes only, 4 also checks three
    /// interior lines per cell in each direction.
    pub samples_per_cell: usize,
    /// Lower bound on the reduced bulk modulus `-v² (∂p/∂v)_T / (R T)` (1 for an
    /// ideal gas, 0 at a critical point). Excludes the critical region.
    pub minimum_reduced_bulk_modulus: f64,
    /// TPD stability trials used at every point.
    pub stability: StabilityOptions,
}

impl Default for QualificationCriteria {
    fn default() -> Self {
        Self {
            samples_per_cell: 4,
            minimum_reduced_bulk_modulus: 0.25,
            stability: StabilityOptions::default(),
        }
    }
}

/// Evidence retained from qualification.
#[derive(Debug, Clone, PartialEq)]
pub struct QualificationReport {
    pub points_checked: usize,
    /// Smallest reduced bulk modulus seen (closest approach to criticality).
    pub minimum_reduced_bulk_modulus: f64,
    /// Smallest TPD among nontrivial trials, if any trial was nontrivial.
    pub minimum_nontrivial_tpd: Option<f64>,
    pub criteria: QualificationCriteria,
}

/// What to do for a state outside the approved domain or composition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutsidePolicy {
    /// Return [`TableError::OutsideApprovedDomain`] / [`TableError::CompositionMismatch`].
    Error,
    /// Run the direct, stability-checked TP flash and return its single-phase
    /// properties labelled [`Source::DirectModel`]; a split is an error.
    DirectModel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Table,
    DirectModel,
}

/// Molar SI single-phase properties in one consistent set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// J/mol on the package caloric reference.
    pub gibbs_energy: f64,
    pub enthalpy: f64,
    /// J/(mol K)
    pub entropy: f64,
    pub internal_energy: f64,
    /// m³/mol
    pub molar_volume: f64,
    /// kg/m³
    pub mass_density: f64,
    pub compressibility: f64,
    /// J/(mol K)
    pub cp: f64,
    pub cv: f64,
    /// (∂v/∂T)_p, m³/(mol K)
    pub dv_dt_p: f64,
    /// (∂v/∂p)_T, m³/(mol Pa)
    pub dv_dp_t: f64,
    /// m/s
    pub speed_of_sound: f64,
}

impl Properties {
    /// Convert a direct homogeneous state into the same property set.
    pub fn from_direct(state: &HomogeneousState) -> Self {
        let t = state.temperature_k;
        Self {
            temperature_k: t,
            pressure_pa: state.pressure_pa,
            gibbs_energy: state.enthalpy - t * state.entropy,
            enthalpy: state.enthalpy,
            entropy: state.entropy,
            internal_energy: state.internal_energy,
            molar_volume: state.molar_volume,
            mass_density: state.mass_density,
            compressibility: state.compressibility,
            cp: state.cp,
            cv: state.cv,
            dv_dt_p: -state.dp_dt / state.dp_dv,
            dv_dp_t: 1.0 / state.dp_dv,
            speed_of_sound: state.speed_of_sound,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Evaluation {
    pub source: Source,
    pub properties: Properties,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TableError {
    InvalidDomain(String),
    /// Qualification failed at a sample point; the table is not built.
    QualificationFailed {
        temperature_k: f64,
        pressure_pa: f64,
        reason: String,
    },
    OutsideApprovedDomain {
        temperature_k: f64,
        pressure_pa: f64,
        domain: String,
    },
    /// PH pressure lies outside the approved pressure range (temperature is an
    /// output of PH and therefore not reported).
    PressureOutsideApprovedDomain {
        pressure_pa: f64,
        domain: String,
    },
    /// PH target enthalpy lies outside the table's enthalpy range at this pressure.
    EnthalpyOutsideApprovedDomain {
        pressure_pa: f64,
        enthalpy: f64,
        minimum: f64,
        maximum: f64,
    },
    CompositionMismatch,
    /// Direct fallback found more than one phase; single-phase properties are undefined.
    DirectResultNotSinglePhase {
        phases: usize,
    },
    Calculation(CalculationError),
}

impl fmt::Display for TableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDomain(m) => write!(f, "invalid table domain: {m}"),
            Self::QualificationFailed {
                temperature_k,
                pressure_pa,
                reason,
            } => write!(
                f,
                "table domain not qualified at T={temperature_k} K, p={pressure_pa} Pa: {reason}"
            ),
            Self::OutsideApprovedDomain {
                temperature_k,
                pressure_pa,
                domain,
            } => write!(
                f,
                "T={temperature_k} K, p={pressure_pa} Pa is outside the approved table domain ({domain}); use the direct model"
            ),
            Self::PressureOutsideApprovedDomain {
                pressure_pa,
                domain,
            } => write!(
                f,
                "p={pressure_pa} Pa is outside the approved table pressure range ({domain}); use the direct model"
            ),
            Self::EnthalpyOutsideApprovedDomain {
                pressure_pa,
                enthalpy,
                minimum,
                maximum,
            } => write!(
                f,
                "h={enthalpy} J/mol at p={pressure_pa} Pa is outside the approved table range {minimum}..={maximum} J/mol; use the direct model"
            ),
            Self::CompositionMismatch => {
                write!(f, "composition differs from the tabulated fixed mixture")
            }
            Self::DirectResultNotSinglePhase { phases } => write!(
                f,
                "direct TP flash found {phases} phases; single-phase table properties are undefined"
            ),
            Self::Calculation(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for TableError {}

impl From<CalculationError> for TableError {
    fn from(value: CalculationError) -> Self {
        Self::Calculation(value)
    }
}

impl From<InputError> for TableError {
    fn from(value: InputError) -> Self {
        Self::Calculation(value.into())
    }
}

/// Node data in (T, x = ln p): g, ∂g/∂T, ∂g/∂x, ∂²g/∂T∂x.
#[derive(Debug, Clone, Copy)]
struct Node([f64; 4]);

/// Tensor-product cubic Hermite surface of g(T, ln p) over a qualified domain.
#[derive(Debug, Clone)]
pub struct GibbsTable {
    domain: TableDomain,
    nodes: Vec<Node>,
    molar_mass: f64,
    ln_p_min: f64,
    report: QualificationReport,
}

/// Hermite basis on [0,1]: values, first and second derivatives of
/// [h00, h01, h10, h11] (value at 0, value at 1, slope at 0, slope at 1).
fn hermite(u: f64) -> [[f64; 4]; 3] {
    let u2 = u * u;
    let u3 = u2 * u;
    [
        [
            2.0 * u3 - 3.0 * u2 + 1.0,
            -2.0 * u3 + 3.0 * u2,
            u3 - 2.0 * u2 + u,
            u3 - u2,
        ],
        [
            6.0 * u2 - 6.0 * u,
            -6.0 * u2 + 6.0 * u,
            3.0 * u2 - 4.0 * u + 1.0,
            3.0 * u2 - 2.0 * u,
        ],
        [
            12.0 * u - 6.0,
            -12.0 * u + 6.0,
            6.0 * u - 4.0,
            6.0 * u - 2.0,
        ],
    ]
}

fn reduced_bulk_modulus(state: &HomogeneousState) -> f64 {
    -state.molar_volume * state.molar_volume * state.dp_dv / (R * state.temperature_k)
}

fn mixture_molar_mass(model: &CubicModel) -> f64 {
    model
        .components()
        .iter()
        .zip(model.fractions())
        .map(|(c, x)| x * c.molar_mass().value())
        .sum()
}

impl GibbsTable {
    /// Qualify the domain for the model's fixed composition, then tabulate.
    /// Only the exact PR variant is accepted: PR-BM's alpha second derivative is
    /// one-sided at each component Tc, which would put cp kinks inside cells.
    pub fn build(
        model: &CubicModel,
        domain: TableDomain,
        criteria: QualificationCriteria,
    ) -> Result<Self, TableError> {
        domain.validate()?;
        if model.variant() != CubicVariant::PengRobinson {
            return Err(TableError::Calculation(
                CalculationError::UnsupportedCapability(
                    "tabulation experiment qualified for exact PR only".into(),
                ),
            ));
        }
        if criteria.samples_per_cell == 0
            || criteria.samples_per_cell > 64
            || !criteria.minimum_reduced_bulk_modulus.is_finite()
            || criteria.minimum_reduced_bulk_modulus <= 0.0
        {
            return Err(TableError::InvalidDomain(
                "samples_per_cell in 1..=64 and a positive bulk-modulus bound required".into(),
            ));
        }
        let report = qualify(model, &domain, criteria)?;
        let ln_p_min = domain.minimum_pressure_pa.ln();
        let mut nodes = Vec::with_capacity(domain.temperature_nodes * domain.pressure_nodes);
        for i in 0..domain.temperature_nodes {
            for j in 0..domain.pressure_nodes {
                let (t, p) = (domain.temperature(i), domain.pressure(j));
                let s = model.evaluate(
                    Temperature::new(t)?,
                    Pressure::new(p)?,
                    RootSelection::MinimumGibbs,
                )?;
                let dv_dt = -s.dp_dt / s.dp_dv;
                nodes.push(Node([
                    s.enthalpy - t * s.entropy,
                    -s.entropy,
                    p * s.molar_volume,
                    p * dv_dt,
                ]));
            }
        }
        Ok(Self {
            domain,
            nodes,
            molar_mass: mixture_molar_mass(model),
            ln_p_min,
            report,
        })
    }

    pub fn domain(&self) -> &TableDomain {
        &self.domain
    }

    pub fn qualification(&self) -> &QualificationReport {
        &self.report
    }

    /// Stored node data in bytes (excluding the small fixed header).
    pub fn memory_bytes(&self) -> usize {
        self.nodes.len() * std::mem::size_of::<Node>()
    }

    /// Interpolate inside the approved domain; outside it, a clear error.
    pub fn evaluate(&self, temperature_k: f64, pressure_pa: f64) -> Result<Properties, TableError> {
        if !temperature_k.is_finite() || !pressure_pa.is_finite() {
            return Err(InputError::NonFinite("temperature or pressure").into());
        }
        if !self.domain.contains(temperature_k, pressure_pa) {
            return Err(TableError::OutsideApprovedDomain {
                temperature_k,
                pressure_pa,
                domain: self.domain.describe(),
            });
        }
        Ok(self.interpolate(temperature_k, pressure_pa))
    }

    /// Temperature and properties at which the table enthalpy equals `enthalpy`
    /// (J/mol) at `pressure_pa`. Safeguarded Newton on the table's own h(T) with
    /// cp as slope; converges to |Δh| ≤ 1e-9 max(|h|, 1 J/mol). Targets outside
    /// the table's h range at this pressure return a clear error.
    pub fn ph(&self, pressure_pa: f64, enthalpy: f64) -> Result<Properties, TableError> {
        if !pressure_pa.is_finite() || !enthalpy.is_finite() {
            return Err(InputError::NonFinite("pressure or enthalpy").into());
        }
        let d = &self.domain;
        if !(d.minimum_pressure_pa..=d.maximum_pressure_pa).contains(&pressure_pa) {
            return Err(TableError::PressureOutsideApprovedDomain {
                pressure_pa,
                domain: d.describe(),
            });
        }
        let (mut lo, mut hi) = (d.minimum_temperature_k, d.maximum_temperature_k);
        let h_lo = self.interpolate(lo, pressure_pa).enthalpy;
        let h_hi = self.interpolate(hi, pressure_pa).enthalpy;
        if !(h_lo..=h_hi).contains(&enthalpy) {
            return Err(TableError::EnthalpyOutsideApprovedDomain {
                pressure_pa,
                enthalpy,
                minimum: h_lo,
                maximum: h_hi,
            });
        }
        let tolerance = 1e-9 * enthalpy.abs().max(1.0);
        let mut t = lo + (enthalpy - h_lo) / (h_hi - h_lo) * (hi - lo);
        for _ in 0..100 {
            let state = self.interpolate(t, pressure_pa);
            let r = state.enthalpy - enthalpy;
            if r.abs() <= tolerance {
                return Ok(state);
            }
            if r > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
            let newton = t - r / state.cp;
            t = if newton > lo && newton < hi && state.cp > 0.0 {
                newton
            } else {
                0.5 * (lo + hi)
            };
        }
        Err(TableError::Calculation(CalculationError::NonConvergence(
            "table PH inversion did not converge in 100 iterations".into(),
        )))
    }

    fn interpolate(&self, t: f64, p: f64) -> Properties {
        let d = &self.domain;
        let ht = d.temperature_step();
        let hx = d.log_pressure_step();
        let x = p.ln();
        let ft = (t - d.minimum_temperature_k) / ht;
        let fx = (x - self.ln_p_min) / hx;
        let i = (ft.floor().max(0.0) as usize).min(d.temperature_nodes - 2);
        let j = (fx.floor().max(0.0) as usize).min(d.pressure_nodes - 2);
        let bu = hermite(ft - i as f64);
        let bw = hermite(fx - j as f64);
        let np = d.pressure_nodes;
        let corner = |a: usize, b: usize| self.nodes[(i + a) * np + (j + b)].0;
        // q[k][l]: ∂^(k+l) g / ∂T^k ∂x^l for k,l in {0,1,2} (k+l ≤ 3 used).
        let mut q = [[0.0; 3]; 3];
        for a in 0..2 {
            for b in 0..2 {
                let n = corner(a, b);
                for k in 0..3 {
                    let st = 1.0 / ht.powi(k as i32);
                    for l in 0..3 {
                        let sx = 1.0 / hx.powi(l as i32);
                        q[k][l] += st
                            * sx
                            * (bu[k][a] * bw[l][b] * n[0]
                                + bu[k][2 + a] * ht * bw[l][b] * n[1]
                                + bu[k][a] * bw[l][2 + b] * hx * n[2]
                                + bu[k][2 + a] * ht * bw[l][2 + b] * hx * n[3]);
                    }
                }
            }
        }
        let g = q[0][0];
        let entropy = -q[1][0];
        let v = q[0][1] / p;
        let dv_dt_p = q[1][1] / p;
        let dv_dp_t = (q[0][2] - q[0][1]) / (p * p);
        let cp = -t * q[2][0];
        let cv = cp + t * dv_dt_p * dv_dt_p / dv_dp_t;
        let enthalpy = g + t * entropy;
        Properties {
            temperature_k: t,
            pressure_pa: p,
            gibbs_energy: g,
            enthalpy,
            entropy,
            internal_energy: enthalpy - p * v,
            molar_volume: v,
            mass_density: self.molar_mass / v,
            compressibility: p * v / (R * t),
            cp,
            cv,
            dv_dt_p,
            dv_dp_t,
            speed_of_sound: (v * v / (-dv_dp_t) * cp / (self.molar_mass * cv)).sqrt(),
        }
    }
}

fn qualify(
    model: &CubicModel,
    domain: &TableDomain,
    criteria: QualificationCriteria,
) -> Result<QualificationReport, TableError> {
    let m = criteria.samples_per_cell;
    let nt = (domain.temperature_nodes - 1) * m + 1;
    let np = (domain.pressure_nodes - 1) * m + 1;
    let ht = (domain.maximum_temperature_k - domain.minimum_temperature_k) / (nt - 1) as f64;
    let hx = (domain.maximum_pressure_pa / domain.minimum_pressure_pa).ln() / (np - 1) as f64;
    let mut report = QualificationReport {
        points_checked: 0,
        minimum_reduced_bulk_modulus: f64::INFINITY,
        minimum_nontrivial_tpd: None,
        criteria,
    };
    for i in 0..nt {
        let t = if i + 1 == nt {
            domain.maximum_temperature_k
        } else {
            domain.minimum_temperature_k + i as f64 * ht
        };
        for j in 0..np {
            let p = if j + 1 == np {
                domain.maximum_pressure_pa
            } else {
                domain.minimum_pressure_pa * (j as f64 * hx).exp()
            };
            let fail = |reason: String| TableError::QualificationFailed {
                temperature_k: t,
                pressure_pa: p,
                reason,
            };
            let (tt, pp) = (Temperature::new(t)?, Pressure::new(p)?);
            let roots = model
                .compressibility_roots(tt, pp)
                .map_err(|e| fail(e.to_string()))?;
            if roots.len() != 1 {
                return Err(fail(format!(
                    "{} cubic roots: near a phase boundary",
                    roots.len()
                )));
            }
            let state = model
                .evaluate(tt, pp, RootSelection::MinimumGibbs)
                .map_err(|e| fail(e.to_string()))?;
            let kappa = reduced_bulk_modulus(&state);
            if kappa < criteria.minimum_reduced_bulk_modulus {
                return Err(fail(format!(
                    "reduced bulk modulus {kappa} below {}: critical region",
                    criteria.minimum_reduced_bulk_modulus
                )));
            }
            report.minimum_reduced_bulk_modulus = report.minimum_reduced_bulk_modulus.min(kappa);
            let stability = assess_stability(model, tt, pp, criteria.stability)
                .map_err(|e| fail(e.to_string()))?;
            if !stability.stable_in_trials {
                return Err(fail("TPD trial shows instability: two-phase region".into()));
            }
            for trial in &stability.trials {
                if trial.outcome == thermo_solvers::TrialOutcome::Converged {
                    let tm = trial.tangent_plane_distance;
                    report.minimum_nontrivial_tpd =
                        Some(report.minimum_nontrivial_tpd.map_or(tm, |m: f64| m.min(tm)));
                }
            }
            report.points_checked += 1;
        }
    }
    Ok(report)
}

/// A fixed-composition model with its qualified table and an explicit fallback rule.
#[derive(Debug, Clone)]
pub struct TabulatedMixture {
    model: CubicModel,
    table: GibbsTable,
    policy: OutsidePolicy,
    flash: FlashOptions,
}

impl TabulatedMixture {
    pub fn new(
        model: CubicModel,
        domain: TableDomain,
        criteria: QualificationCriteria,
        policy: OutsidePolicy,
    ) -> Result<Self, TableError> {
        let table = GibbsTable::build(&model, domain, criteria)?;
        Ok(Self {
            model,
            table,
            policy,
            flash: FlashOptions::default(),
        })
    }

    pub fn table(&self) -> &GibbsTable {
        &self.table
    }

    pub fn model(&self) -> &CubicModel {
        &self.model
    }

    pub fn policy(&self) -> OutsidePolicy {
        self.policy
    }

    /// The same qualified table with another outside-domain policy.
    pub fn with_policy(&self, policy: OutsidePolicy) -> Self {
        Self {
            policy,
            ..self.clone()
        }
    }

    /// Evaluate at the tabulated composition.
    pub fn evaluate(&self, temperature_k: f64, pressure_pa: f64) -> Result<Evaluation, TableError> {
        if !temperature_k.is_finite() || !pressure_pa.is_finite() {
            return Err(InputError::NonFinite("temperature or pressure").into());
        }
        if self.table.domain.contains(temperature_k, pressure_pa) {
            return Ok(Evaluation {
                source: Source::Table,
                properties: self.table.interpolate(temperature_k, pressure_pa),
            });
        }
        match self.policy {
            OutsidePolicy::Error => Err(TableError::OutsideApprovedDomain {
                temperature_k,
                pressure_pa,
                domain: self.table.domain.describe(),
            }),
            OutsidePolicy::DirectModel => {
                direct(&self.model, temperature_k, pressure_pa, self.flash)
            }
        }
    }

    /// PH at the tabulated composition: the table when p and h lie in the approved
    /// domain, otherwise per policy a clear error or the direct PH flash. A table
    /// inversion that does not converge is treated like a state outside the
    /// domain: an error under [`OutsidePolicy::Error`], the direct PH flash under
    /// [`OutsidePolicy::DirectModel`]. Input errors are never retried.
    pub fn ph(&self, pressure_pa: f64, enthalpy: f64) -> Result<Evaluation, TableError> {
        self.ph_from_table(self.table.ph(pressure_pa, enthalpy), pressure_pa, enthalpy)
    }

    /// Apply the outside policy to a table PH outcome.
    fn ph_from_table(
        &self,
        table: Result<Properties, TableError>,
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<Evaluation, TableError> {
        match table {
            Ok(properties) => Ok(Evaluation {
                source: Source::Table,
                properties,
            }),
            Err(
                e @ (TableError::OutsideApprovedDomain { .. }
                | TableError::PressureOutsideApprovedDomain { .. }
                | TableError::EnthalpyOutsideApprovedDomain { .. }
                | TableError::Calculation(CalculationError::NonConvergence(_))),
            ) => match self.policy {
                OutsidePolicy::Error => Err(e),
                OutsidePolicy::DirectModel => {
                    let flash = ph_flash(
                        &self.model,
                        Pressure::new(pressure_pa)?,
                        enthalpy,
                        PhFlashOptions {
                            tp: self.flash,
                            ..PhFlashOptions::default()
                        },
                    )?;
                    match flash.phases.as_slice() {
                        [single] => Ok(Evaluation {
                            source: Source::DirectModel,
                            properties: Properties::from_direct(&single.state),
                        }),
                        phases => Err(TableError::DirectResultNotSinglePhase {
                            phases: phases.len(),
                        }),
                    }
                }
            },
            Err(e) => Err(e),
        }
    }

    /// Evaluate at given mole fractions (model order). The table applies only to
    /// the tabulated composition (exact match within 1e-12 per entry).
    pub fn evaluate_at(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<Evaluation, TableError> {
        let same = fractions.len() == self.model.fractions().len()
            && fractions
                .iter()
                .zip(self.model.fractions())
                .all(|(a, b)| (a - b).abs() <= 1e-12);
        if same {
            return self.evaluate(temperature_k, pressure_pa);
        }
        match self.policy {
            OutsidePolicy::Error => Err(TableError::CompositionMismatch),
            OutsidePolicy::DirectModel => {
                let model = self.model.with_fractions(fractions)?;
                direct(&model, temperature_k, pressure_pa, self.flash)
            }
        }
    }
}

fn direct(
    model: &CubicModel,
    temperature_k: f64,
    pressure_pa: f64,
    options: FlashOptions,
) -> Result<Evaluation, TableError> {
    let flash = tp_flash(
        model,
        Temperature::new(temperature_k)?,
        Pressure::new(pressure_pa)?,
        options,
    )?;
    match flash.phases.as_slice() {
        [single] => Ok(Evaluation {
            source: Source::DirectModel,
            properties: Properties::from_direct(&single.state),
        }),
        phases => Err(TableError::DirectResultNotSinglePhase {
            phases: phases.len(),
        }),
    }
}

/// The A11 experiment definition shared by tests and examples.
pub mod experiment {
    use super::*;

    /// Dry natural-gas composition of `examples/homogeneous-natural-gas.json`.
    pub const COMPONENTS: [&str; 5] =
        ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"];
    pub const FRACTIONS: [f64; 5] = [0.88, 0.06, 0.03, 0.02, 0.01];

    /// Approved A11 domain: 280–420 K, 1–80 bar; 2.5 K steps in T, ln p steps of
    /// ln(80)/90 ≈ 0.0487 (grid chosen by the refinement study in docs/a11-status.md).
    pub const DOMAIN: TableDomain = TableDomain {
        minimum_temperature_k: 280.0,
        maximum_temperature_k: 420.0,
        minimum_pressure_pa: 1.0e5,
        maximum_pressure_pa: 8.0e6,
        temperature_nodes: 57,
        pressure_nodes: 91,
    };

    /// Exact PR with explicit synthetic zero kij, as in the natural-gas examples.
    pub fn model() -> CubicModel {
        let composition = Composition::new(
            CompositionBasis::Mole,
            COMPONENTS
                .iter()
                .zip(FRACTIONS)
                .map(|(id, x)| (ComponentId::new(*id).expect("valid id"), x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .expect("valid composition");
        CubicModel::new(
            &ComponentDatabase::embedded().expect("embedded data"),
            &composition,
            CubicVariant::PengRobinson,
            vec![vec![0.0; 5]; 5],
            "Explicit synthetic zero-kij approximation; not fitted VLE data",
        )
        .expect("valid model")
    }

    pub fn mixture(policy: OutsidePolicy) -> Result<TabulatedMixture, TableError> {
        TabulatedMixture::new(model(), DOMAIN, QualificationCriteria::default(), policy)
    }

    /// Fixed-seed SplitMix64 uniform points in (T, ln p) over the domain.
    pub fn sample_points(seed: u64, count: usize) -> Vec<(f64, f64)> {
        let mut state = seed;
        let mut next = || {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
        };
        let d = DOMAIN;
        let (lo, hi) = (d.minimum_pressure_pa.ln(), d.maximum_pressure_pa.ln());
        (0..count)
            .map(|_| {
                let t = d.minimum_temperature_k
                    + next() * (d.maximum_temperature_k - d.minimum_temperature_k);
                let p = (lo + next() * (hi - lo)).exp();
                (t, p)
            })
            .collect()
    }

    /// Maximum errors of table against direct evaluation over a point set.
    /// Values: absolute g/h/s (J/mol, J/(mol K)), relative v and density.
    /// Derivatives: relative cp, cv, (∂v/∂T)_p, (∂v/∂p)_T and speed of sound.
    #[derive(Debug, Clone, Copy, Default, PartialEq)]
    pub struct ErrorSummary {
        pub points: usize,
        pub gibbs_abs: f64,
        pub enthalpy_abs: f64,
        pub entropy_abs: f64,
        pub volume_rel: f64,
        pub density_rel: f64,
        pub cp_rel: f64,
        pub cv_rel: f64,
        pub dv_dt_rel: f64,
        pub dv_dp_rel: f64,
        pub speed_of_sound_rel: f64,
    }

    /// Cell centres in (T, ln p): the expected worst-case interpolation locations.
    pub fn cell_centres(domain: &TableDomain) -> Vec<(f64, f64)> {
        let mut points = Vec::new();
        for i in 0..domain.temperature_nodes - 1 {
            for j in 0..domain.pressure_nodes - 1 {
                let t = 0.5 * (domain.temperature(i) + domain.temperature(i + 1));
                let p = (domain.pressure(j) * domain.pressure(j + 1)).sqrt();
                points.push((t, p));
            }
        }
        points
    }

    pub fn compare(
        table: &GibbsTable,
        model: &CubicModel,
        points: &[(f64, f64)],
    ) -> Result<ErrorSummary, TableError> {
        let rel = |a: f64, b: f64| ((a - b) / b).abs();
        let mut e = ErrorSummary {
            points: points.len(),
            ..ErrorSummary::default()
        };
        for &(t, p) in points {
            let a = table.evaluate(t, p)?;
            let b = Properties::from_direct(&model.evaluate(
                Temperature::new(t)?,
                Pressure::new(p)?,
                RootSelection::MinimumGibbs,
            )?);
            e.gibbs_abs = e.gibbs_abs.max((a.gibbs_energy - b.gibbs_energy).abs());
            e.enthalpy_abs = e.enthalpy_abs.max((a.enthalpy - b.enthalpy).abs());
            e.entropy_abs = e.entropy_abs.max((a.entropy - b.entropy).abs());
            e.volume_rel = e.volume_rel.max(rel(a.molar_volume, b.molar_volume));
            e.density_rel = e.density_rel.max(rel(a.mass_density, b.mass_density));
            e.cp_rel = e.cp_rel.max(rel(a.cp, b.cp));
            e.cv_rel = e.cv_rel.max(rel(a.cv, b.cv));
            e.dv_dt_rel = e.dv_dt_rel.max(rel(a.dv_dt_p, b.dv_dt_p));
            e.dv_dp_rel = e.dv_dp_rel.max(rel(a.dv_dp_t, b.dv_dp_t));
            e.speed_of_sound_rel = e
                .speed_of_sound_rel
                .max(rel(a.speed_of_sound, b.speed_of_sound));
        }
        Ok(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small qualified table (nodes only) over part of the A11 domain.
    fn coarse(policy: OutsidePolicy) -> TabulatedMixture {
        let domain = TableDomain {
            minimum_temperature_k: 320.0,
            maximum_temperature_k: 380.0,
            minimum_pressure_pa: 1e6,
            maximum_pressure_pa: 5e6,
            temperature_nodes: 7,
            pressure_nodes: 5,
        };
        let criteria = QualificationCriteria {
            samples_per_cell: 1,
            ..Default::default()
        };
        TabulatedMixture::new(experiment::model(), domain, criteria, policy).unwrap()
    }

    fn non_convergence() -> Result<Properties, TableError> {
        Err(TableError::Calculation(CalculationError::NonConvergence(
            "table PH inversion did not converge in 100 iterations".into(),
        )))
    }

    #[test]
    fn table_ph_non_convergence_follows_the_outside_policy() {
        let fallback = coarse(OutsidePolicy::DirectModel);
        let p = 3e6;
        let h = fallback.table().evaluate(350.0, p).unwrap().enthalpy;

        let e = fallback.ph_from_table(non_convergence(), p, h).unwrap();
        assert_eq!(e.source, Source::DirectModel);
        let flash = ph_flash(
            fallback.model(),
            Pressure::new(p).unwrap(),
            h,
            PhFlashOptions::default(),
        )
        .unwrap();
        assert_eq!(e.properties.temperature_k, flash.temperature_k);

        let strict = fallback.with_policy(OutsidePolicy::Error);
        assert!(matches!(
            strict.ph_from_table(non_convergence(), p, h),
            Err(TableError::Calculation(CalculationError::NonConvergence(_)))
        ));
        // Input errors are not retried by the direct model.
        assert!(matches!(
            fallback.ph_from_table(Err(InputError::NonFinite("x").into()), p, h),
            Err(TableError::Calculation(CalculationError::InvalidInput(_)))
        ));
    }

    #[test]
    fn table_ph_pressure_outside_domain_reports_no_temperature() {
        let strict = coarse(OutsidePolicy::Error);
        let e = strict.table().ph(1.2e7, 0.0).unwrap_err();
        assert!(
            matches!(
                e,
                TableError::PressureOutsideApprovedDomain { pressure_pa, .. } if pressure_pa == 1.2e7
            ),
            "{e:?}"
        );
        assert!(!e.to_string().contains("NaN"), "{e}");
        assert!(matches!(
            strict.ph(1.2e7, 0.0),
            Err(TableError::PressureOutsideApprovedDomain { .. })
        ));
    }
}
