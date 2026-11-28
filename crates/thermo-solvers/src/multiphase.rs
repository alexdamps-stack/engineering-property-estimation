//! Multiphase (VLLE/LLE) stability and TP/PH flash (#224) for the PR/PR-BM cubic
//! models and the NRTL + ideal-gas (gamma–phi) activity package.
//!
//! **Stability.** Michelsen tangent-plane-distance trials from a finite,
//! documented set of starts: Wilson vapor- and liquid-like trials, an *aqueous*
//! trial (near-pure water, when water is present), a *hydrocarbon* trial (the
//! water-free feed restricted to its C/H components, or to all non-water
//! components when none is a hydrocarbon) and a near-pure trial per remaining
//! component. Unlike the two-phase assessment every trial runs, so all distinct
//! stationary points are available to initialise phases. Passing the trials is
//! evidence, not a proof, of stability.
//!
//! **Split.** Michelsen's multiphase formulation: for fixed fugacity
//! coefficients the phase amounts β ≥ 0 minimise the convex function
//! Q(β) = Σ_k β_k − Σ_i z_i ln Σ_k β_k/φ_ik (projected Newton with an active set).
//! A phase whose amount reaches zero is removed; the compositions
//! x_ik = z_i/(φ_ik Σ_l β_l/φ_il) update φ by successive substitution with
//! dominant-eigenvalue extrapolation and safeguarded Newton steps on ln K
//! relative to the first phase. After convergence the trials are repeated from
//! the equilibrium chemical potentials; a negative TPD adds that phase (at most
//! three phases). A fourth phase is reported as unsupported.
//!
//! **Plait points.** Two phases of the same kind that converge to the same
//! composition are merged. A split that keeps re-forming after a merge, or that
//! does not converge while two phases lie within 0.05 in mole fraction, fails
//! with an explicit plait/critical-point message. Converged phases closer than
//! 0.01 carry a warning diagnostic.
//!
//! **Labels.** Phases are reported vapor first, then liquids by increasing mass
//! density (light liquid, heavy liquid). For the cubic models the largest-volume
//! phase of a three-phase result is the vapor; of a two-phase result it is the
//! vapor unless the single-phase heuristic labels it liquid (LLE). Gamma–phi
//! phases carry the kind (liquid or ideal-gas vapor) they were solved with.
//! Liquid mass densities of the activity package come from its uncoupled
//! saturated-liquid density correlation.

mod adapters;
mod core;
mod ph;

pub use ph::{activity_vlle_ph_flash, vlle_ph_flash};

use thermo_models::{ActivityModel, ActivityPhaseState, CubicModel, HomogeneousState};
use thermo_types::{CalculationError as Error, Diagnostic, Phase, Pressure, Temperature};

use crate::stability::TrialOutcome;
use adapters::{ActivityAt, CubicAt};

/// Iteration limits, tolerances and the PH temperature bracket of the
/// multiphase solvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MultiphaseOptions {
    /// Successive-substitution limit of one split solve (including Newton steps).
    pub max_iterations: usize,
    /// Convergence: max_i (max_k − min_k) of ln(x_ik φ_ik) over active components.
    pub tolerance: f64,
    /// Successive-substitution limit of one TPD trial.
    pub trial_max_iterations: usize,
    /// TPD trial convergence: max |Δ ln W_i| between iterations.
    pub trial_tolerance: f64,
    /// Modified TPD below this value classifies a phase as unstable.
    pub tm_threshold: f64,
    /// Try safeguarded Newton steps on ln K once the substitution is close.
    pub newton: bool,
    /// PH temperature bracket, K.
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// PH: |h − h_target| / max(|h_target|, 1 J/mol).
    pub energy_tolerance: f64,
    /// PH outer temperature iteration limit.
    pub max_temperature_iterations: usize,
}

impl Default for MultiphaseOptions {
    fn default() -> Self {
        Self {
            max_iterations: 2000,
            tolerance: 1e-10,
            trial_max_iterations: 2000,
            trial_tolerance: 1e-10,
            tm_threshold: -1e-8,
            newton: true,
            minimum_temperature_k: thermo_models::CUBIC_TEMPERATURE_RANGE_K[0],
            maximum_temperature_k: thermo_models::CUBIC_TEMPERATURE_RANGE_K[1],
            energy_tolerance: 1e-10,
            max_temperature_iterations: 200,
        }
    }
}

/// Largest accepted split, trial or energy tolerance.
pub const MAX_MULTIPHASE_TOLERANCE: f64 = 1e-6;

impl MultiphaseOptions {
    /// Rejects tolerances outside (0, 1e-6], a TPD threshold outside [−1e-6, 0)
    /// and an unordered or nonfinite bracket (`UnsupportedCapability`), and zero
    /// iteration limits (`NonConvergence`).
    pub fn validate(&self) -> Result<(), Error> {
        let tolerance = |t: f64| t.is_finite() && t > 0.0 && t <= MAX_MULTIPHASE_TOLERANCE;
        if !tolerance(self.tolerance)
            || !tolerance(self.trial_tolerance)
            || !tolerance(self.energy_tolerance)
            || !self.tm_threshold.is_finite()
            || self.tm_threshold >= 0.0
            || self.tm_threshold < crate::MIN_TM_THRESHOLD
            || !self.minimum_temperature_k.is_finite()
            || !self.maximum_temperature_k.is_finite()
            || self.minimum_temperature_k <= 0.0
            || self.minimum_temperature_k >= self.maximum_temperature_k
        {
            return Err(Error::UnsupportedCapability(format!(
                "multiphase options require tolerances in (0, {MAX_MULTIPHASE_TOLERANCE:e}], a negative TPD threshold in [{:e}, 0) and an ordered finite temperature bracket",
                crate::MIN_TM_THRESHOLD
            )));
        }
        if self.max_iterations == 0
            || self.trial_max_iterations == 0
            || self.max_temperature_iterations == 0
        {
            return Err(Error::NonConvergence(
                "multiphase iteration limit is zero".into(),
            ));
        }
        Ok(())
    }
}

/// Phase kind a trial or phase is evaluated with. The cubic models use one
/// fluid kind (minimum-Gibbs root at each composition); the gamma–phi package
/// distinguishes its NRTL liquid and ideal-gas vapor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluidKind {
    /// Cubic EOS, minimum-Gibbs root at the phase composition.
    Fluid,
    /// Gamma–phi NRTL liquid.
    Liquid,
    /// Gamma–phi ideal-gas vapor.
    Vapor,
}

/// Initialisation of a multiphase TPD trial.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiphaseTrial {
    /// W_i = z_i K_i with Wilson K (gamma–phi: Raoult K = p_i^sat/p), vapor kind.
    WilsonVapor,
    /// W_i = z_i / K_i, liquid kind.
    WilsonLiquid,
    /// Aqueous trial phase: near-pure water (0.999), liquid kind.
    Aqueous,
    /// Hydrocarbon trial phase: 0.999 of the water-free feed over its C/H
    /// components (all non-water components when none is a hydrocarbon), liquid kind.
    Hydrocarbon,
    /// Near-pure trial of this component index (model order), liquid kind.
    PureComponent(usize),
}

#[derive(Debug, Clone)]
pub struct MultiphaseTrialResult {
    pub trial: MultiphaseTrial,
    pub kind: FluidKind,
    pub outcome: TrialOutcome,
    /// Modified TPD tm = 1 + Σ W_i (ln W_i + ln φ_i(w) − d_i − 1), d_i = ln(x_i φ_i) of the reference.
    pub tangent_plane_distance: f64,
    /// Normalised trial composition in model order (zeros for inactive components).
    pub composition: Vec<f64>,
    pub iterations: usize,
}

/// Multiphase stability from the finite trial set of the module documentation.
#[derive(Debug, Clone)]
pub struct MultiphaseStability {
    /// True when no trial located a stationary point with tm below the threshold.
    pub stable_in_trials: bool,
    /// Kind of the reference phase (the minimum-Gibbs kind of the tested composition).
    pub reference_kind: FluidKind,
    pub trials: Vec<MultiphaseTrialResult>,
}

impl MultiphaseStability {
    /// Trials that prove instability, most negative TPD first.
    pub fn unstable_trials(&self, threshold: f64) -> Vec<&MultiphaseTrialResult> {
        let mut t: Vec<_> = self
            .trials
            .iter()
            .filter(|t| {
                t.outcome.locates_stationary_point() && t.tangent_plane_distance < threshold
            })
            .collect();
        t.sort_by(|a, b| {
            a.tangent_plane_distance
                .total_cmp(&b.tangent_plane_distance)
        });
        t
    }
}

/// Role of a phase for unit outlets: vapor, the only liquid, or one of two
/// liquids ordered by mass density.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseRole {
    /// Vapor (also a single non-liquid phase, including one labelled supercritical).
    Vapor,
    /// The only liquid phase.
    Liquid,
    /// The liquid with the lower mass density of two.
    LightLiquid,
    /// The liquid with the higher mass density of two.
    HeavyLiquid,
}

#[derive(Debug, Clone)]
pub struct MultiphasePhase<S> {
    /// Label only (see the module documentation).
    pub phase: Phase,
    pub role: PhaseRole,
    pub kind: FluidKind,
    /// Molar phase fraction of the feed.
    pub phase_fraction: f64,
    /// Mole fractions in model order.
    pub composition: Vec<f64>,
    pub state: S,
}

#[derive(Debug, Clone)]
pub struct MultiphaseTpFlash<S> {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Vapor first, then liquids by increasing mass density.
    pub phases: Vec<MultiphasePhase<S>>,
    pub feed_stability: MultiphaseStability,
    /// Trials repeated from the equilibrium chemical potentials of the final
    /// split; `None` for a stable single-phase feed.
    pub final_stability: Option<MultiphaseStability>,
    /// Split iterations over all phase-addition rounds (0 for one phase).
    pub iterations: usize,
    /// max_i (max_k − min_k) ln(x_ik φ_ik) over active components (0 for one phase).
    pub max_ln_fugacity_residual: f64,
    /// max_i |Σ_k β_k x_ik − z_i|.
    pub max_component_balance_residual: f64,
    /// (G_split − G_feed)/(R T) per mole of feed; negative for a split, 0 for one phase.
    pub reduced_gibbs_change: f64,
    pub diagnostics: Vec<Diagnostic>,
}

impl<S> MultiphaseTpFlash<S> {
    /// Number of liquid phases.
    pub fn liquid_count(&self) -> usize {
        self.phases
            .iter()
            .filter(|p| p.phase == Phase::Liquid)
            .count()
    }
}

#[derive(Debug, Clone)]
pub struct MultiphasePhFlash<S> {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub target_enthalpy: f64,
    pub enthalpy: f64,
    /// (h − h_target) / max(|h_target|, 1 J/mol).
    pub scaled_energy_residual: f64,
    /// As [`MultiphaseTpFlash::phases`]. At an invariant temperature (a pure
    /// component, or the three-phase temperature of a binary) the phases come from
    /// both sides of the enthalpy step and their amounts from the lever rule.
    pub phases: Vec<MultiphasePhase<S>>,
    /// True when the target lies in an enthalpy step at an invariant temperature.
    pub invariant_temperature: bool,
    /// Outer temperature iterations.
    pub iterations: usize,
    /// max_i (max_k − min_k) ln(x_ik φ_ik) of the returned phases.
    pub max_ln_fugacity_residual: f64,
    /// max_i |Σ_k β_k x_ik − z_i|.
    pub max_component_balance_residual: f64,
    /// TP evidence at the solution temperature, or at the lower temperature
    /// bracket endpoint for an invariant solution.
    pub tp: MultiphaseTpFlash<S>,
    pub diagnostics: Vec<Diagnostic>,
}

/// Multiphase stability of the cubic model's composition at T and p.
pub fn vlle_stability(
    model: &CubicModel,
    temperature: Temperature,
    pressure: Pressure,
    options: MultiphaseOptions,
) -> Result<MultiphaseStability, Error> {
    options.validate()?;
    let a = CubicAt::new(model, temperature, pressure);
    core::stability(&a, model.fractions(), options)
}

/// Multiphase stability of the gamma–phi model's composition at T and p
/// (NRTL liquid and ideal-gas vapor trial phases).
pub fn activity_vlle_stability(
    model: &ActivityModel,
    temperature: Temperature,
    pressure: Pressure,
    options: MultiphaseOptions,
) -> Result<MultiphaseStability, Error> {
    options.validate()?;
    let a = ActivityAt::new(model, temperature, pressure)?;
    core::stability(&a, model.fractions(), options)
}

/// Isothermal–isobaric flash of the cubic model's composition with up to three
/// phases (vapor–liquid–liquid, liquid–liquid or vapor–liquid).
pub fn vlle_tp_flash(
    model: &CubicModel,
    temperature: Temperature,
    pressure: Pressure,
    options: MultiphaseOptions,
) -> Result<MultiphaseTpFlash<HomogeneousState>, Error> {
    options.validate()?;
    let a = CubicAt::new(model, temperature, pressure);
    core::tp_flash(&a, model.fractions(), options)
}

/// Isothermal–isobaric gamma–phi flash with up to three phases (one ideal-gas
/// vapor and up to two NRTL liquids, or three liquids without vapor).
pub fn activity_vlle_tp_flash(
    model: &ActivityModel,
    temperature: Temperature,
    pressure: Pressure,
    options: MultiphaseOptions,
) -> Result<MultiphaseTpFlash<ActivityPhaseState>, Error> {
    options.validate()?;
    let a = ActivityAt::new(model, temperature, pressure)?;
    core::tp_flash(&a, model.fractions(), options)
}
