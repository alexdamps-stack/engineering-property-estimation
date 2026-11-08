//! Two-phase TP flash: TPD stability of the feed, safeguarded Newton / accelerated
//! successive substitution on ln K with a negative-flash Rachford–Rice inner solve, Gibbs
//! decrease check, and post-split stability of each phase.
use crate::equilibrium_model::EquilibriumModel;
use crate::rachford_rice::{RACHFORD_RICE_MAX_ITERATIONS, rachford_rice_negative};
use crate::stability::{
    StabilityAssessment, StabilityOptions, TrialOutcome, assess_stability_impl, normalized,
};
use thermo_models::{CubicModel, FugacityState, HomogeneousState, RootSelection};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase, Pressure, Temperature,
};

/// Critical Z/B of the exact-factor PR cubic (triple root): Zc = 0.30740130869...,
/// Bc = OMEGA_B. Used only for the single-phase label heuristic.
pub const PR_CRITICAL_VOLUME_OVER_COVOLUME: f64 = 0.307_401_308_698_7 / thermo_models::OMEGA_B;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlashOptions {
    /// Iteration limit for the phase split (including Newton steps).
    pub max_iterations: usize,
    /// Convergence: max_i |ln(y_i φ_i^V) − ln(x_i φ_i^L)| over active components.
    pub tolerance: f64,
    /// Try safeguarded Newton steps after three SS iterations; use dominant-
    /// eigenvalue extrapolation every fifth SS step when Newton is unavailable.
    pub acceleration: bool,
    /// Enable Newton within acceleration. Set false to retain eigenvalue-only SS.
    pub newton: bool,
    pub stability: StabilityOptions,
}

impl Default for FlashOptions {
    fn default() -> Self {
        Self {
            max_iterations: 2000,
            tolerance: 1e-10,
            acceleration: true,
            newton: true,
            stability: StabilityOptions::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct FlashPhase {
    /// Label only: two phases by molar volume; single phase by the documented heuristic.
    pub phase: Phase,
    /// Molar phase fraction of the feed.
    pub phase_fraction: f64,
    /// Mole fractions in model order.
    pub composition: Vec<f64>,
    /// Full homogeneous properties at the phase composition (minimum-Gibbs root).
    pub state: HomogeneousState,
}

#[derive(Debug, Clone)]
pub struct TpFlash {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    /// One phase, or vapor-labelled (larger molar volume) first, then liquid.
    pub phases: Vec<FlashPhase>,
    pub feed_stability: StabilityAssessment,
    /// Stability of each split phase, same order as `phases`; empty for one phase.
    pub phase_stability: Vec<StabilityAssessment>,
    /// Phase-split iterations, including Newton steps (0 for one phase).
    pub iterations: usize,
    /// max_i |ln f_i^V − ln f_i^L| over active components (0 for one phase).
    pub max_ln_fugacity_residual: f64,
    /// max_i |β y_i + (1−β) x_i − z_i| (mole-fraction units).
    pub max_component_balance_residual: f64,
    /// (G_split − G_feed)/(R T) per mole of feed; negative for a split, 0 for one phase.
    pub reduced_gibbs_change: f64,
    pub diagnostics: Vec<Diagnostic>,
}

impl TpFlash {
    /// Molar vapor fraction: the vapor-labelled phase fraction for a split, or 1/0 for
    /// a single phase labelled vapor/liquid. None for a supercritical single phase.
    pub fn vapor_fraction(&self) -> Option<f64> {
        match self.phases.as_slice() {
            [single] => match single.phase {
                Phase::Vapor => Some(1.0),
                Phase::Liquid => Some(0.0),
                Phase::Supercritical => None,
            },
            [vapor, _] => Some(vapor.phase_fraction),
            _ => None,
        }
    }
}

fn info(code: &str, message: impl Into<String>, severity: DiagnosticSeverity) -> Diagnostic {
    Diagnostic {
        code: code.into(),
        severity,
        message: message.into(),
    }
}

/// Isothermal–isobaric flash of the model's composition with at most two phases.
pub fn tp_flash(
    model: &CubicModel,
    temperature: Temperature,
    pressure: Pressure,
    options: FlashOptions,
) -> Result<TpFlash, Error> {
    tp_flash_impl(model, temperature, pressure, options)
}

pub(crate) fn tp_flash_impl<M: EquilibriumModel>(
    model: &M,
    temperature: Temperature,
    pressure: Pressure,
    options: FlashOptions,
) -> Result<TpFlash, Error> {
    let z = model.fractions().to_vec();
    let n = z.len();
    let active: Vec<usize> = (0..n).filter(|i| z[*i] > 0.0).collect();
    let feed_stability = assess_stability_impl(model, temperature, pressure, options.stability)?;
    let mut diagnostics = vec![info(
        "A4-STABILITY-FINITE-TRIALS",
        format!(
            "Feed stability from {} TPD trials (Wilson vapor/liquid{}); finite trials do not prove global stability.",
            feed_stability.trials.len(),
            if options.stability.pure_component_trials {
                ", near-pure components"
            } else {
                ""
            }
        ),
        DiagnosticSeverity::Information,
    )];
    if feed_stability.stable_in_trials {
        let state = model.evaluate(temperature, pressure, RootSelection::MinimumGibbs)?;
        let phase = model.phase_label(&state);
        diagnostics.push(info(
            "A4-SINGLE-PHASE-LABEL",
            "Single-phase label from the model adapter's documented pseudo-critical/density heuristic; a label, not a phase-boundary result.",
            DiagnosticSeverity::Information,
        ));
        return Ok(TpFlash {
            temperature_k: temperature.value(),
            pressure_pa: pressure.value(),
            phases: vec![FlashPhase {
                phase,
                phase_fraction: 1.0,
                composition: z,
                state,
            }],
            feed_stability,
            phase_stability: Vec::new(),
            iterations: 0,
            max_ln_fugacity_residual: 0.0,
            max_component_balance_residual: 0.0,
            reduced_gibbs_change: 0.0,
            diagnostics,
        });
    }

    // Initial K from the most negative trial: the trial is the lighter or heavier
    // incipient phase depending on its volume relative to the feed root.
    let trial = feed_stability
        .trials
        .iter()
        .filter(|t| t.outcome != TrialOutcome::TrivialSolution)
        .min_by(|a, b| {
            a.tangent_plane_distance
                .total_cmp(&b.tangent_plane_distance)
        })
        .ok_or_else(|| Error::NumericalFailure("unstable feed without a trial".into()))?;
    let trial_is_vapor = trial.compressibility > feed_stability.feed_compressibility;
    let mut ln_k: Vec<f64> = active
        .iter()
        .map(|&i| {
            let ratio = (trial.composition[i] / z[i]).ln();
            if trial_is_vapor { ratio } else { -ratio }
        })
        .collect();
    let z_active: Vec<f64> = active.iter().map(|&i| z[i]).collect();
    let g = |w: &[f64], f: &FugacityState| -> f64 {
        active
            .iter()
            .map(|&i| w[i] * (w[i].ln() + f.ln_fugacity_coefficients[i]))
            .sum()
    };
    let mut history: Vec<Vec<f64>> = Vec::new();
    let mut newton_steps = 0;
    let mut converged = None;
    for iteration in 1..=options.max_iterations {
        let k: Vec<f64> = ln_k.iter().map(|l| l.exp()).collect();
        let rr = rachford_rice_negative(&z_active, &k, RACHFORD_RICE_MAX_ITERATIONS)?;
        let (x, y) = (
            expand(&active, &rr.liquid, n),
            expand(&active, &rr.vapor, n),
        );
        let liquid = phase_fugacity(model, temperature, pressure, &x)?;
        let vapor = phase_fugacity(model, temperature, pressure, &y)?;
        let mut residual: f64 = 0.0;
        let mut residuals = Vec::with_capacity(active.len());
        let mut next = Vec::with_capacity(active.len());
        for &i in &active {
            let r = y[i].ln() + vapor.ln_fugacity_coefficients[i]
                - x[i].ln()
                - liquid.ln_fugacity_coefficients[i];
            residual = residual.max(r.abs());
            residuals.push(r);
            next.push(liquid.ln_fugacity_coefficients[i] - vapor.ln_fugacity_coefficients[i]);
        }
        if !residual.is_finite() {
            return Err(Error::NumericalFailure(
                "nonfinite fugacity residual".into(),
            ));
        }
        if residual <= options.tolerance {
            converged = Some((iteration, rr.vapor_fraction, x, y, liquid, vapor, residual));
            break;
        }
        if next.iter().all(|l| l.abs() < 1e-6) {
            return Err(Error::NumericalFailure(
                "successive substitution approached the trivial K=1 solution; near-critical or nearly identical phases".into(),
            ));
        }
        if options.acceleration
            && options.newton
            && iteration >= 3
            && rr.vapor_fraction > 0.0
            && rr.vapor_fraction < 1.0
        {
            let current_g =
                rr.vapor_fraction * g(&y, &vapor) + (1.0 - rr.vapor_fraction) * g(&x, &liquid);
            let evaluate = |candidate: &[f64]| -> Result<(Vec<f64>, f64), Error> {
                let k: Vec<f64> = candidate.iter().map(|l| l.exp()).collect();
                let rr = rachford_rice_negative(&z_active, &k, RACHFORD_RICE_MAX_ITERATIONS)?;
                if !(rr.vapor_fraction > 0.0 && rr.vapor_fraction < 1.0) {
                    return Err(Error::NumericalFailure(
                        "Newton trial has no physical phase fraction".into(),
                    ));
                }
                let x = expand(&active, &rr.liquid, n);
                let y = expand(&active, &rr.vapor, n);
                let l = phase_fugacity(model, temperature, pressure, &x)?;
                let v = phase_fugacity(model, temperature, pressure, &y)?;
                let candidate_g =
                    rr.vapor_fraction * g(&y, &v) + (1.0 - rr.vapor_fraction) * g(&x, &l);
                Ok((
                    active
                        .iter()
                        .map(|&i| {
                            y[i].ln() + v.ln_fugacity_coefficients[i]
                                - x[i].ln()
                                - l.ln_fugacity_coefficients[i]
                        })
                        .collect(),
                    candidate_g,
                ))
            };
            if let Some(candidate) = newton_step(&ln_k, &residuals, current_g, evaluate) {
                ln_k = candidate;
                newton_steps += 1;
                history.clear();
                continue;
            }
        }
        history.push(next.clone());
        ln_k = next;
        if options.acceleration && history.len() >= 3 && iteration % 5 == 0 {
            if let Some(extrapolated) = dominant_eigenvalue_step(&history) {
                ln_k = extrapolated;
            }
            history.clear();
        }
    }
    let Some((iterations, beta, x, y, liquid, vapor, residual)) = converged else {
        return Err(Error::NonConvergence(format!(
            "TP flash exhausted {} iterations (successive substitution / safeguarded Newton)",
            options.max_iterations
        )));
    };
    if !(beta > 0.0 && beta < 1.0) {
        return Err(Error::NumericalFailure(format!(
            "converged split has vapor fraction {beta} outside (0, 1) although the feed failed stability trials"
        )));
    }
    let feed = model.ln_fugacity(temperature, pressure, RootSelection::MinimumGibbs)?;
    let reduced_gibbs_change = beta * g(&y, &vapor) + (1.0 - beta) * g(&x, &liquid) - g(&z, &feed);
    if reduced_gibbs_change >= 0.0 {
        return Err(Error::NumericalFailure(format!(
            "converged split does not lower the Gibbs energy (Δg/RT = {reduced_gibbs_change:e})"
        )));
    }
    let max_component_balance_residual = (0..n)
        .map(|i| (beta * y[i] + (1.0 - beta) * x[i] - z[i]).abs())
        .fold(0.0, f64::max);
    let mut phases = Vec::new();
    for (composition, fraction) in [(y, beta), (x, 1.0 - beta)] {
        let state = model.with_fractions(&composition)?.evaluate(
            temperature,
            pressure,
            RootSelection::MinimumGibbs,
        )?;
        phases.push(FlashPhase {
            phase: Phase::Vapor,
            phase_fraction: fraction,
            composition,
            state,
        });
    }
    if phases[0].state.molar_volume < phases[1].state.molar_volume {
        phases.swap(0, 1);
    }
    phases[1].phase = Phase::Liquid;
    let (v_light, v_heavy) = (phases[0].state.molar_volume, phases[1].state.molar_volume);
    if (v_light - v_heavy).abs() <= 1e-3 * v_heavy {
        diagnostics.push(info(
            "A4-NEARLY-IDENTICAL-PHASES",
            "Phase molar volumes differ by less than 0.1%; vapor/liquid labels are not meaningful near a critical point.",
            DiagnosticSeverity::Warning,
        ));
    }
    let mut phase_stability = Vec::new();
    for phase in &phases {
        let assessment = assess_stability_impl(
            &model.with_fractions(&phase.composition)?,
            temperature,
            pressure,
            options.stability,
        )?;
        if !assessment.stable_in_trials {
            return Err(Error::UnsupportedCapability(
                "post-split stability indicates an additional phase; more than two phases are not supported".into(),
            ));
        }
        phase_stability.push(assessment);
    }
    diagnostics.push(info(
        "A4-NEWTON-STEPS",
        format!("Accepted {newton_steps} safeguarded Newton steps; rejected acceleration falls back to successive substitution."),
        DiagnosticSeverity::Information,
    ));
    diagnostics.push(info(
        "A4-PHASE-LABEL-BY-VOLUME",
        "Split phases are labelled vapor (larger molar volume) and liquid; liquid-liquid splits are not distinguished.",
        DiagnosticSeverity::Information,
    ));
    Ok(TpFlash {
        temperature_k: temperature.value(),
        pressure_pa: pressure.value(),
        phases,
        feed_stability,
        phase_stability,
        iterations,
        max_ln_fugacity_residual: residual,
        max_component_balance_residual,
        reduced_gibbs_change,
        diagnostics,
    })
}

fn expand(active: &[usize], values: &[f64], n: usize) -> Vec<f64> {
    let mut full = vec![0.0; n];
    for (p, &i) in active.iter().enumerate() {
        full[i] = values[p];
    }
    normalized(full)
}

fn phase_fugacity<M: EquilibriumModel>(
    model: &M,
    temperature: Temperature,
    pressure: Pressure,
    composition: &[f64],
) -> Result<FugacityState, Error> {
    model.with_fractions(composition)?.ln_fugacity(
        temperature,
        pressure,
        RootSelection::MinimumGibbs,
    )
}

/// Newton on active ln K, differentiating the complete residual (including the
/// implicit RR phase fraction) with central differences. Failed probes, singular
/// Jacobians and failed line searches leave the caller's SS update available.
fn newton_step(
    point: &[f64],
    residual: &[f64],
    current_g: f64,
    evaluate: impl Fn(&[f64]) -> Result<(Vec<f64>, f64), Error>,
) -> Option<Vec<f64>> {
    let n = point.len();
    let mut jacobian = vec![vec![0.0; n]; n];
    for j in 0..n {
        let h = f64::EPSILON.cbrt() * point[j].abs().max(1.0);
        let mut probe = point.to_vec();
        probe[j] += h;
        let (plus, _) = evaluate(&probe).ok()?;
        probe[j] = point[j] - h;
        let (minus, _) = evaluate(&probe).ok()?;
        for i in 0..n {
            jacobian[i][j] = (plus[i] - minus[i]) / (2.0 * h);
        }
    }
    let mut step: Vec<f64> = residual.iter().map(|r| -r).collect();
    // Small dense system with partial pivoting; no inverse is formed.
    for j in 0..n {
        let pivot =
            (j..n).max_by(|&a, &b| jacobian[a][j].abs().total_cmp(&jacobian[b][j].abs()))?;
        jacobian.swap(j, pivot);
        step.swap(j, pivot);
        if !jacobian[j][j].is_finite() || jacobian[j][j].abs() < 1e-12 {
            return None;
        }
        for i in j + 1..n {
            let factor = jacobian[i][j] / jacobian[j][j];
            let (upper, lower) = jacobian.split_at_mut(i);
            for (value, pivot_value) in lower[0].iter_mut().zip(&upper[j]).skip(j + 1) {
                *value -= factor * pivot_value;
            }
            step[i] -= factor * step[j];
        }
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            step[i] -= jacobian[i][j] * step[j];
        }
        step[i] /= jacobian[i][i];
    }
    if step.iter().any(|s| !s.is_finite()) {
        return None;
    }
    // Limit each ln K change to one, preserving direction. Backtrack for strict
    // sufficient decrease in the infinity norm of the fugacity residual and
    // no Gibbs increase beyond roundoff (dimensionless, per mole of feed).
    let norm = |r: &[f64]| r.iter().map(|r| r.abs()).fold(0.0, f64::max);
    let mut scale = 1.0 / norm(&step).max(1.0);
    for _ in 0..12 {
        let candidate: Vec<f64> = point
            .iter()
            .zip(&step)
            .map(|(p, s)| p + scale * s)
            .collect();
        if !candidate.iter().all(|k| k.abs() < 1e-6)
            && let Ok((r, g)) = evaluate(&candidate)
            && g.is_finite()
            && g <= current_g + 1e-14
            && r.iter().all(|r| r.is_finite())
            && norm(&r) < (1.0 - 1e-4 * scale) * norm(residual)
        {
            return Some(candidate);
        }
        scale *= 0.5;
    }
    None
}

/// Michelsen-style dominant-eigenvalue extrapolation from the last three SS iterates.
fn dominant_eigenvalue_step(history: &[Vec<f64>]) -> Option<Vec<f64>> {
    let [.., a, b, c] = history else {
        return None;
    };
    let d1: Vec<f64> = b.iter().zip(a).map(|(b, a)| b - a).collect();
    let d2: Vec<f64> = c.iter().zip(b).map(|(c, b)| c - b).collect();
    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(u, v)| u * v).sum::<f64>();
    let denominator = dot(&d1, &d1);
    if denominator <= 0.0 {
        return None;
    }
    let lambda = dot(&d2, &d1) / denominator;
    if !(lambda > 0.0 && lambda < 0.95) {
        return None;
    }
    Some(
        c.iter()
            .zip(&d2)
            .map(|(c, d)| c + d * lambda / (1.0 - lambda))
            .collect(),
    )
}

/// Single-phase label: supercritical above Kay's pseudo-critical T and p; otherwise
/// vapor when v/b exceeds the PR critical ratio Zc/Bc, liquid below it.
pub fn single_phase_label(model: &CubicModel, state: &HomogeneousState) -> Phase {
    let (tpc, ppc) =
        model
            .components()
            .iter()
            .zip(model.fractions())
            .fold((0.0, 0.0), |(t, p), (c, x)| {
                (
                    t + x * c.critical_temperature(),
                    p + x * c.critical_pressure(),
                )
            });
    if state.temperature_k > tpc && state.pressure_pa > ppc {
        Phase::Supercritical
    } else if state.molar_volume / model.covolume() >= PR_CRITICAL_VOLUME_OVER_COVOLUME {
        Phase::Vapor
    } else {
        Phase::Liquid
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newton_solves_a_coupled_system_with_row_pivoting() {
        let candidate = newton_step(&[2.0, 2.0], &[2.0, 1.0], 0.0, |x| {
            Ok((
                vec![2.0 * (x[1] - 1.0), 3.0 * (x[0] - 3.0) + 4.0 * (x[1] - 1.0)],
                0.0,
            ))
        })
        .unwrap();
        assert!((candidate[0] - 3.0).abs() < 1e-9);
        assert!((candidate[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn newton_backtracks_when_trial_leaves_model_domain() {
        let candidate = newton_step(&[2.0], &[1.0], 0.0, |x| {
            if x[0] < 1.8 {
                Err(Error::OutsideValidityDomain(
                    "trial outside test domain".into(),
                ))
            } else {
                Ok((vec![x[0] - 1.0], 0.0))
            }
        })
        .unwrap();
        assert!((candidate[0] - 1.875).abs() < 1e-10);
    }

    #[test]
    fn unusable_newton_steps_leave_successive_substitution_available() {
        // A smaller residual must not steer the split uphill toward K=1.
        assert!(newton_step(&[2.0], &[1.0], 0.0, |x| Ok((vec![x[0] - 1.0], 1.0))).is_none());
        assert!(newton_step(&[2.0], &[1.0], 0.0, |_| Ok((vec![1.0], 0.0))).is_none());
        assert!(newton_step(&[2.0], &[1.0], 0.0, |_| Ok((vec![f64::NAN], 0.0))).is_none());
        assert!(
            newton_step(&[2.0], &[1.0], 0.0, |_| {
                Err(Error::NumericalFailure("probe failed".into()))
            })
            .is_none()
        );
        // Probes are allowed, but every line-search candidate is out of domain.
        assert!(
            newton_step(&[2.0], &[1.0], 0.0, |x| {
                if (x[0] - 2.0).abs() > 1e-4 {
                    Err(Error::OutsideValidityDomain(
                        "trial outside test domain".into(),
                    ))
                } else {
                    Ok((vec![x[0] - 1.0], 0.0))
                }
            })
            .is_none()
        );
    }
}
