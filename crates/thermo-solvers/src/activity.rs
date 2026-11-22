//! Gamma–phi (NRTL + ideal-gas vapor) equilibrium: bubble/dew pressure and
//! temperature, TP flash with an NRTL liquid tangent-plane check, and PH flash.
//!
//! With an ideal-gas vapor, liquid fugacities independent of pressure and a single
//! liquid phase, the feed is liquid for p ≥ p_bubble(T, z), vapor for
//! p ≤ p_dew(T, z), and two-phase in between. A liquid that fails the finite NRTL
//! tangent-plane trials (liquid–liquid split) is reported as unsupported; no VLLE.
//!
//! Laid out like the cubic solvers: `saturation` (bubble/dew points),
//! `flash` (TP flash), `ph` (PH flash) and `stability` (liquid TPD trials).
//! This module holds the shared options, result types and helpers.

mod flash;
mod ph;
mod saturation;
mod stability;

pub use flash::activity_tp_flash;
pub use ph::activity_ph_flash;
pub use saturation::{
    activity_bubble_pressure, activity_bubble_temperature, activity_dew_pressure,
    activity_dew_temperature,
};
pub use stability::liquid_stability;

use thermo_models::{
    ACTIVITY_TEMPERATURE_RANGE_K, ActivityModel, ActivityPhaseState, ExcessState,
    INFERRED_REGRESSION_INTERVAL_CODE, INFERRED_REGRESSION_INTERVAL_WARNING, PureLiquidState,
};
use thermo_types::{CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase, Temperature};

/// Iteration limits and tolerances of the gamma–phi solvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivityOptions {
    /// Successive-substitution limit (dew composition, TP split, TPD trials).
    pub max_iterations: usize,
    /// Convergence: max_i |ln(y_i p) − ln(x_i γ_i p_i^sat)| over active components.
    pub tolerance: f64,
    /// Outer temperature iteration limit (bubble/dew temperature, PH).
    pub max_temperature_iterations: usize,
    /// PH: |h − h_target| / energy scale (see `activity_ph_flash`), in (0, 1e-6].
    pub energy_tolerance: f64,
    /// Modified TPD below this value marks the liquid unstable (liquid–liquid split).
    pub tm_threshold: f64,
    /// Temperature bracket for bubble/dew temperatures and PH, within the package range.
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
}

impl Default for ActivityOptions {
    fn default() -> Self {
        Self {
            max_iterations: 500,
            tolerance: 1e-13,
            max_temperature_iterations: 200,
            energy_tolerance: 1e-10,
            tm_threshold: -1e-8,
            minimum_temperature_k: ACTIVITY_TEMPERATURE_RANGE_K[0],
            maximum_temperature_k: ACTIVITY_TEMPERATURE_RANGE_K[1],
        }
    }
}

impl ActivityOptions {
    fn validate(&self) -> Result<(), Error> {
        let [lo, hi] = ACTIVITY_TEMPERATURE_RANGE_K;
        let ok = self.tolerance.is_finite()
            && self.tolerance > 0.0
            && self.tolerance <= 1e-6
            && self.energy_tolerance.is_finite()
            && self.energy_tolerance > 0.0
            && self.energy_tolerance <= 1e-6
            && self.tm_threshold.is_finite()
            && self.tm_threshold < 0.0
            && self.minimum_temperature_k >= lo
            && self.maximum_temperature_k <= hi
            && self.minimum_temperature_k < self.maximum_temperature_k;
        if !ok {
            return Err(Error::UnsupportedCapability(format!(
                "activity options: tolerances in (0, 1e-6], negative TPD threshold, ordered bracket within {lo}..={hi} K"
            )));
        }
        if self.max_iterations == 0 || self.max_temperature_iterations == 0 {
            return Err(Error::NonConvergence(
                "activity iteration limit is zero".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaturationKind {
    Bubble,
    Dew,
}

/// Saturation point of the model composition; the incipient phase has zero amount.
#[derive(Debug, Clone)]
pub struct ActivitySaturation {
    pub kind: SaturationKind,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// Liquid mole fractions (feed for bubble, incipient for dew), model order.
    pub liquid: Vec<f64>,
    /// Vapor mole fractions (incipient for bubble, feed for dew), model order.
    pub vapor: Vec<f64>,
    pub max_ln_fugacity_residual: f64,
    /// Composition iterations at the final temperature plus outer temperature iterations.
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub struct LiquidStabilityTrial {
    /// Near-pure start for this component index (model order).
    pub component: usize,
    pub tangent_plane_distance: f64,
    pub composition: Vec<f64>,
    pub iterations: usize,
    /// Trial converged back to the tested composition.
    pub trivial: bool,
}

/// Finite near-pure NRTL TPD trials: evidence, not proof, of liquid stability.
#[derive(Debug, Clone)]
pub struct LiquidStability {
    pub stable_in_trials: bool,
    pub trials: Vec<LiquidStabilityTrial>,
}

#[derive(Debug, Clone)]
pub struct ActivityFlashPhase {
    pub phase: Phase,
    pub phase_fraction: f64,
    pub composition: Vec<f64>,
    pub state: ActivityPhaseState,
}

#[derive(Debug, Clone)]
pub struct ActivityTpFlash {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// One phase, or vapor first then liquid.
    pub phases: Vec<ActivityFlashPhase>,
    pub bubble_pressure_pa: f64,
    pub dew_pressure_pa: f64,
    pub iterations: usize,
    pub max_ln_fugacity_residual: f64,
    pub max_component_balance_residual: f64,
    /// NRTL liquid tangent-plane trials of the liquid phase, if one is present.
    pub liquid_stability: Option<LiquidStability>,
    /// True when any reported phase's NRTL validity rests on an inferred regression
    /// interval ([`ActivityPhaseState::inferred_regression_interval`]); a warning
    /// diagnostic is then included. Phase validity statuses are unchanged.
    pub inferred_regression_interval: bool,
    pub diagnostics: Vec<Diagnostic>,
}

impl ActivityTpFlash {
    /// Phase-fraction-weighted molar enthalpy, J/mol of feed.
    pub fn enthalpy(&self) -> f64 {
        self.phases
            .iter()
            .map(|p| p.phase_fraction * p.state.enthalpy)
            .sum()
    }
    pub fn vapor_fraction(&self) -> f64 {
        self.phases
            .iter()
            .filter(|p| p.phase == Phase::Vapor)
            .map(|p| p.phase_fraction)
            .sum()
    }
}
#[derive(Debug, Clone)]
pub struct ActivityPhFlash {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub target_enthalpy: f64,
    pub enthalpy: f64,
    /// (h − h_target) / energy_scale.
    pub scaled_energy_residual: f64,
    /// Energy scale of the branch that produced the result, J/mol (see `activity_ph_flash`).
    pub energy_scale: f64,
    /// Vapor first, then liquid. Pure saturation endpoints retain a zero-fraction phase.
    pub phases: Vec<ActivityFlashPhase>,
    /// Feed bubble temperature at p; `None` when it lies outside the option bracket.
    pub bubble_temperature_k: Option<f64>,
    /// Feed dew temperature at p; `None` when it lies outside the option bracket.
    pub dew_temperature_k: Option<f64>,
    pub iterations: usize,
    /// TP evidence for a two-phase mixture result.
    pub tp: Option<ActivityTpFlash>,
    /// As [`ActivityTpFlash::inferred_regression_interval`], over `phases`.
    pub inferred_regression_interval: bool,
    pub diagnostics: Vec<Diagnostic>,
}

/// Whether any phase rests on an inferred NRTL regression interval; if so, appends
/// the warning diagnostic. Only a flag is added: no result value changes.
fn flag_inferred_interval(
    phases: &[ActivityFlashPhase],
    diagnostics: &mut Vec<Diagnostic>,
) -> bool {
    let inferred = phases.iter().any(|p| p.state.inferred_regression_interval);
    if inferred {
        diagnostics.push(diagnostic(
            INFERRED_REGRESSION_INTERVAL_CODE,
            INFERRED_REGRESSION_INTERVAL_WARNING,
            DiagnosticSeverity::Warning,
        ));
    }
    inferred
}

fn diagnostic(code: &str, message: impl Into<String>, severity: DiagnosticSeverity) -> Diagnostic {
    Diagnostic {
        code: code.into(),
        severity,
        message: message.into(),
    }
}

fn active(z: &[f64]) -> Vec<usize> {
    (0..z.len()).filter(|i| z[*i] > 0.0).collect()
}

/// Pure-liquid states of every component at T.
fn pures(model: &ActivityModel, t: f64) -> Result<Vec<PureLiquidState>, Error> {
    model.pure_liquids(Temperature::new(t)?)
}

/// NRTL excess state at a solver temperature iterate in K.
fn excess(model: &ActivityModel, t: f64, x: &[f64]) -> Result<ExcessState, Error> {
    model.nrtl().excess_at(Temperature::new(t)?, x)
}

fn check_pressure(p: f64) -> Result<(), Error> {
    let [lo, hi] = thermo_models::ACTIVITY_PRESSURE_RANGE_PA;
    if !(lo..=hi).contains(&p) {
        return Err(Error::OutsideValidityDomain(format!(
            "gamma-phi saturation pressure {p} Pa outside {lo}..={hi} Pa"
        )));
    }
    Ok(())
}
