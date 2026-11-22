//! Gamma–phi TP flash with the NRTL liquid tangent-plane check.

use thermo_models::{ActivityModel, ActivityPhaseState};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, Phase, Pressure, Temperature,
};

use super::saturation::{bubble_at, dew_at};
use super::stability::liquid_stability;
use super::{
    ActivityFlashPhase, ActivityOptions, ActivityTpFlash, LiquidStability, active, diagnostic,
    excess, flag_inferred_interval, pures,
};
use crate::rachford_rice::{RACHFORD_RICE_MAX_ITERATIONS, rachford_rice_negative};
use crate::stability::normalized as normalize;

pub(super) fn phase_state(
    model: &ActivityModel,
    phase: Phase,
    t: f64,
    p: f64,
    composition: &[f64],
) -> Result<ActivityPhaseState, Error> {
    let m = model.with_fractions(composition)?;
    let (t, p) = (Temperature::new(t)?, Pressure::new(p)?);
    match phase {
        Phase::Liquid => m.liquid(t, p),
        _ => m.vapor(t, p),
    }
}

pub(super) fn require_stable_liquid(
    model: &ActivityModel,
    t: f64,
    x: &[f64],
    options: ActivityOptions,
) -> Result<LiquidStability, Error> {
    let s = liquid_stability(model, Temperature::new(t)?, x, options)?;
    if !s.stable_in_trials {
        return Err(Error::UnsupportedCapability(
            "NRTL liquid fails tangent-plane trials: liquid-liquid split (VLLE/LLE) is not supported".into(),
        ));
    }
    Ok(s)
}

/// Isothermal–isobaric gamma–phi flash of the model composition (at most one liquid).
pub fn activity_tp_flash(
    model: &ActivityModel,
    temperature: Temperature,
    pressure: Pressure,
    options: ActivityOptions,
) -> Result<ActivityTpFlash, Error> {
    options.validate()?;
    tp_flash(model, temperature, pressure, options, true)
}

/// Largest excursion of a converged vapor fraction outside (0, 1) accepted as
/// roundoff at a saturation boundary.
const SATURATION_BETA_ROUNDOFF: f64 = 1e-10;

/// TP flash; `check_stability = false` skips the liquid TPD trials (used inside the
/// PH search, whose accepted result is re-flashed with the trials).
pub(super) fn tp_flash(
    model: &ActivityModel,
    temperature: Temperature,
    pressure: Pressure,
    options: ActivityOptions,
    check_stability: bool,
) -> Result<ActivityTpFlash, Error> {
    model.check_domain(temperature, pressure)?;
    let (t, p) = (temperature.value(), pressure.value());
    let z = model.fractions().to_vec();
    let n = z.len();
    let act = active(&z);
    let pure = pures(model, t)?;
    let (p_bubble, _) = bubble_at(model, t, &z, &pure)?;
    let p_dew = if act.len() == 1 {
        p_bubble
    } else {
        dew_at(model, t, &z, &pure, options)?.pressure
    };
    let stability = |x: &[f64]| -> Result<Option<LiquidStability>, Error> {
        if check_stability {
            Ok(Some(require_stable_liquid(model, t, x, options)?))
        } else {
            Ok(None)
        }
    };
    let mut diagnostics = vec![diagnostic(
        "A7-GAMMA-PHI",
        "Modified Raoult with ideal-gas vapor; single-phase classification from exact bubble/dew pressures of the feed; at most one liquid phase.",
        DiagnosticSeverity::Information,
    )];
    let single = |phase: Phase,
                  stability: Option<LiquidStability>,
                  mut diagnostics: Vec<Diagnostic>|
     -> Result<ActivityTpFlash, Error> {
        let phases = vec![ActivityFlashPhase {
            phase,
            phase_fraction: 1.0,
            composition: z.clone(),
            state: phase_state(model, phase, t, p, &z)?,
        }];
        let inferred_regression_interval = flag_inferred_interval(&phases, &mut diagnostics);
        Ok(ActivityTpFlash {
            temperature_k: t,
            pressure_pa: p,
            phases,
            bubble_pressure_pa: p_bubble,
            dew_pressure_pa: p_dew,
            iterations: 0,
            max_ln_fugacity_residual: 0.0,
            max_component_balance_residual: 0.0,
            liquid_stability: stability,
            inferred_regression_interval,
            diagnostics,
        })
    };
    if p >= p_bubble {
        return single(Phase::Liquid, stability(&z)?, diagnostics);
    }
    if p <= p_dew {
        return single(Phase::Vapor, None, diagnostics);
    }
    let z_active: Vec<f64> = act.iter().map(|&i| z[i]).collect();
    let ln_ps: Vec<f64> = pure.iter().map(|s| s.vapor_pressure.ln()).collect();
    let mut gamma = excess(model, t, &z)?.ln_activity_coefficients;
    let mut converged = None;
    for iteration in 1..=options.max_iterations {
        let k: Vec<f64> = act
            .iter()
            .map(|&i| (gamma[i] + ln_ps[i] - p.ln()).exp())
            .collect();
        let rr = rachford_rice_negative(&z_active, &k, RACHFORD_RICE_MAX_ITERATIONS)?;
        let mut x = vec![0.0; n];
        let mut y = vec![0.0; n];
        for (q, &i) in act.iter().enumerate() {
            x[i] = rr.liquid[q];
            y[i] = rr.vapor[q];
        }
        let (x, y) = (normalize(x), normalize(y));
        let next = excess(model, t, &x)?.ln_activity_coefficients;
        let mut r: f64 = 0.0;
        for &i in &act {
            r = r.max((y[i].ln() + p.ln() - x[i].ln() - next[i] - ln_ps[i]).abs());
        }
        if !r.is_finite() {
            return Err(Error::NumericalFailure("nonfinite TP residual".into()));
        }
        gamma = next;
        if r <= options.tolerance {
            converged = Some((iteration, rr.vapor_fraction, x, y, r));
            break;
        }
    }
    let Some((iterations, beta, x, y, r)) = converged else {
        return Err(Error::NonConvergence(format!(
            "gamma-phi TP split exhausted {} iterations",
            options.max_iterations
        )));
    };
    // Within roundoff of p_dew or p_bubble the split can converge to β marginally
    // outside (0, 1); that is the saturated single phase. Larger excursions fail.
    if !(beta > 0.0 && beta < 1.0) && (beta.min(1.0 - beta)).abs() <= SATURATION_BETA_ROUNDOFF {
        diagnostics.push(diagnostic(
            "A7-SPLIT-AT-SATURATION-BOUNDARY",
            format!("Converged vapor fraction {beta:e} lies within {SATURATION_BETA_ROUNDOFF:e} of a saturation boundary; reported as the saturated single phase."),
            DiagnosticSeverity::Information,
        ));
        return if beta >= 1.0 {
            single(Phase::Vapor, None, diagnostics)
        } else {
            single(Phase::Liquid, stability(&z)?, diagnostics)
        };
    }
    if !(beta > 0.0 && beta < 1.0) {
        return Err(Error::NumericalFailure(format!(
            "converged split has vapor fraction {beta} outside (0, 1) between dew and bubble pressure"
        )));
    }
    let balance = (0..n)
        .map(|i| (beta * y[i] + (1.0 - beta) * x[i] - z[i]).abs())
        .fold(0.0, f64::max);
    let liquid_stability = stability(&x)?;
    if check_stability {
        diagnostics.push(diagnostic(
            "A7-LIQUID-TPD-FINITE-TRIALS",
            "Liquid stability from near-pure NRTL tangent-plane trials; finite trials do not prove global stability.",
            DiagnosticSeverity::Information,
        ));
    }
    let phases = vec![
        ActivityFlashPhase {
            phase: Phase::Vapor,
            phase_fraction: beta,
            state: phase_state(model, Phase::Vapor, t, p, &y)?,
            composition: y,
        },
        ActivityFlashPhase {
            phase: Phase::Liquid,
            phase_fraction: 1.0 - beta,
            state: phase_state(model, Phase::Liquid, t, p, &x)?,
            composition: x,
        },
    ];
    let inferred_regression_interval = flag_inferred_interval(&phases, &mut diagnostics);
    Ok(ActivityTpFlash {
        temperature_k: t,
        pressure_pa: p,
        phases,
        bubble_pressure_pa: p_bubble,
        dew_pressure_pa: p_dew,
        iterations,
        max_ln_fugacity_residual: r,
        max_component_balance_residual: balance,
        liquid_stability,
        inferred_regression_interval,
        diagnostics,
    })
}
