//! Gamma–phi PH flash on the option temperature bracket.

use thermo_models::{ActivityModel, ActivityPhaseState};
use thermo_types::{
    CalculationError as Error, DiagnosticSeverity, InputError, Phase, Pressure, Temperature,
};

use super::flash::{activity_tp_flash, phase_state, require_stable_liquid, tp_flash};
use super::saturation::{
    activity_bubble_temperature, activity_dew_temperature, bubble_at, dew_at, solve_temperature,
};
use super::{
    ActivityFlashPhase, ActivityOptions, ActivityPhFlash, ActivitySaturation, SaturationKind,
    active, check_pressure, diagnostic, pures,
};

/// Position of a feed saturation temperature at p relative to the option bracket.
enum Saturation {
    /// Below the bracket minimum: saturation pressure already exceeds p there.
    Below,
    At(ActivitySaturation),
    /// Above the bracket maximum: saturation pressure is still below p there.
    Above,
}

impl Saturation {
    fn temperature(&self) -> Option<f64> {
        match self {
            Self::At(s) => Some(s.temperature_k),
            _ => None,
        }
    }
}

fn locate_saturation(
    model: &ActivityModel,
    pressure: Pressure,
    options: ActivityOptions,
    kind: SaturationKind,
) -> Result<Saturation, Error> {
    let z = model.fractions();
    let ln_ratio = |t: f64| -> Result<f64, Error> {
        let pure = pures(model, t)?;
        let ps = match kind {
            SaturationKind::Bubble => bubble_at(model, t, z, &pure)?.0,
            SaturationKind::Dew => dew_at(model, t, z, &pure, options)?.pressure,
        };
        Ok(ps.ln() - pressure.value().ln())
    };
    if ln_ratio(options.minimum_temperature_k)? > options.tolerance {
        return Ok(Saturation::Below);
    }
    if ln_ratio(options.maximum_temperature_k)? < -options.tolerance {
        return Ok(Saturation::Above);
    }
    Ok(Saturation::At(match kind {
        SaturationKind::Bubble => activity_bubble_temperature(model, pressure, options)?,
        SaturationKind::Dew => activity_dew_temperature(model, pressure, options)?,
    }))
}

/// Pressure/enthalpy gamma–phi flash on the option temperature bracket.
///
/// The feed's bubble and dew temperatures at p split the bracket into a liquid branch
/// (below T_bubble), a two-phase range and a vapor branch (above T_dew); a saturation
/// temperature outside the bracket removes or truncates the corresponding branch.
/// Single-phase branches are solved by safeguarded Newton with exact c_p and the
/// energy scale max(|h_target|, 1 J/mol), as in A5. The two-phase mixture range uses
/// Illinois regula falsi on the split enthalpy with the scale
/// max(|h_target|, h(T_hi) − h(T_lo), 1 J/mol) over the two-phase range: the split's
/// composition tolerance limits the attainable enthalpy precision to a fraction of
/// that span, so a 1 J/mol floor would be unreachable near h = 0. A single active
/// component uses the lever rule at its saturation temperature.
pub fn activity_ph_flash(
    model: &ActivityModel,
    pressure: Pressure,
    target_enthalpy: f64,
    options: ActivityOptions,
) -> Result<ActivityPhFlash, Error> {
    if !target_enthalpy.is_finite() {
        return Err(InputError::NonFinite("target enthalpy").into());
    }
    options.validate()?;
    let p = pressure.value();
    check_pressure(p)?;
    let z = model.fractions().to_vec();
    let (t_min, t_max) = (options.minimum_temperature_k, options.maximum_temperature_k);
    let pure_feed = active(&z).len() == 1;
    let bubble = locate_saturation(model, pressure, options, SaturationKind::Bubble)?;
    let dew = if pure_feed {
        match &bubble {
            Saturation::Below => Saturation::Below,
            Saturation::At(s) => Saturation::At(s.clone()),
            Saturation::Above => Saturation::Above,
        }
    } else {
        locate_saturation(model, pressure, options, SaturationKind::Dew)?
    };
    let (tb, td) = (bubble.temperature(), dew.temperature());
    // Top of the single-liquid branch and bottom of the single-vapor branch.
    let liquid_top = match &bubble {
        Saturation::Below => None,
        Saturation::At(s) => Some(s.temperature_k),
        Saturation::Above => Some(t_max),
    };
    let vapor_bottom = match &dew {
        Saturation::Below => Some(t_min),
        Saturation::At(s) => Some(s.temperature_k),
        Saturation::Above => None,
    };
    let mut diagnostics = vec![diagnostic(
        "A7-PH-GAMMA-PHI",
        "PH energy uses phase-weighted gamma-phi enthalpy (NASA7 ideal gas; saturation-pressure slope for pure liquids; NRTL excess); single-phase scale max(|h_target|, 1 J/mol), two-phase mixture scale max(|h_target|, two-phase enthalpy span, 1 J/mol).",
        DiagnosticSeverity::Information,
    )];
    if tb.is_none() || td.is_none() {
        diagnostics.push(diagnostic(
            "A7-PH-SATURATION-OUTSIDE-BRACKET",
            format!(
                "Feed bubble and/or dew temperature at {p} Pa lies outside {t_min}..={t_max} K; only the phase branches inside the bracket are searched."
            ),
            DiagnosticSeverity::Information,
        ));
    }
    let finish =
        |t: f64, phases: Vec<ActivityFlashPhase>, iterations, tp, diagnostics, scale: f64| {
            let enthalpy: f64 = phases
                .iter()
                .map(|ph: &ActivityFlashPhase| ph.phase_fraction * ph.state.enthalpy)
                .sum();
            let scaled = (enthalpy - target_enthalpy) / scale;
            if !scaled.is_finite() || scaled.abs() > options.energy_tolerance {
                return Err(Error::NonConvergence(format!(
                    "PH energy residual {scaled:e} exceeds tolerance"
                )));
            }
            Ok(ActivityPhFlash {
                temperature_k: t,
                pressure_pa: p,
                target_enthalpy,
                enthalpy,
                scaled_energy_residual: scaled,
                energy_scale: scale,
                phases,
                bubble_temperature_k: tb,
                dew_temperature_k: td,
                iterations,
                tp,
                diagnostics,
            })
        };
    let one = |phase, state: ActivityPhaseState| {
        vec![ActivityFlashPhase {
            phase,
            phase_fraction: 1.0,
            composition: z.clone(),
            state,
        }]
    };
    let scale = target_enthalpy.abs().max(1.0);
    let budget = options.energy_tolerance * scale;
    // Single-phase branches.
    let branches = [
        liquid_top.map(|top| (Phase::Liquid, t_min, top)),
        vapor_bottom.map(|bottom| (Phase::Vapor, bottom, t_max)),
    ];
    for (phase, lo, hi) in branches.into_iter().flatten() {
        let state = |t: f64| phase_state(model, phase, t, p, &z);
        let edge_t = if phase == Phase::Liquid { hi } else { lo };
        let edge = state(edge_t)?;
        let inside = if phase == Phase::Liquid {
            target_enthalpy <= edge.enthalpy + budget
        } else {
            target_enthalpy >= edge.enthalpy - budget
        };
        if !inside {
            continue;
        }
        let (t, iterations, s) =
            if ((edge.enthalpy - target_enthalpy) / scale).abs() <= options.energy_tolerance {
                (edge_t, 0, edge)
            } else {
                let (t, iterations) = solve_temperature(
                    lo,
                    hi,
                    options,
                    options.energy_tolerance,
                    |t| {
                        let s = state(t)?;
                        Ok(((s.enthalpy - target_enthalpy) / scale, s.cp / scale))
                    },
                    "PH single-phase temperature",
                )?;
                (t, iterations, state(t)?)
            };
        if phase == Phase::Liquid {
            require_stable_liquid(model, t, &z, options)?;
        }
        return finish(t, one(phase, s), iterations, None, diagnostics, scale);
    }
    let outside = || {
        Error::OutsideValidityDomain(format!(
            "PH target {target_enthalpy} J/mol at {p} Pa is not reachable within {t_min}..={t_max} K"
        ))
    };
    if pure_feed {
        let Saturation::At(sat) = &bubble else {
            return Err(outside());
        };
        let t = sat.temperature_k;
        let lb = phase_state(model, Phase::Liquid, t, p, &z)?;
        let vb = phase_state(model, Phase::Vapor, t, p, &z)?;
        let latent = vb.enthalpy - lb.enthalpy;
        if latent <= 0.0 || !latent.is_finite() {
            return Err(Error::NumericalFailure("nonpositive latent heat".into()));
        }
        let beta = ((target_enthalpy - lb.enthalpy) / latent).clamp(0.0, 1.0);
        diagnostics.push(diagnostic(
            "A7-PURE-SATURATION-LEVER-RULE",
            "Single active component: saturation temperature from the vapor-pressure correlation; vapor fraction from the lever rule.",
            DiagnosticSeverity::Information,
        ));
        let phases = vec![
            ActivityFlashPhase {
                phase: Phase::Vapor,
                phase_fraction: beta,
                composition: z.clone(),
                state: vb,
            },
            ActivityFlashPhase {
                phase: Phase::Liquid,
                phase_fraction: 1.0 - beta,
                composition: z.clone(),
                state: lb,
            },
        ];
        return finish(t, phases, sat.iterations, None, diagnostics, scale);
    }
    // Two-phase mixture: the split enthalpy increases with T over [lo, hi].
    let (mut a, mut b) = (liquid_top.unwrap_or(t_min), vapor_bottom.unwrap_or(t_max));
    if a >= b {
        return Err(outside());
    }
    let split = |t: f64| -> Result<f64, Error> {
        Ok(
            tp_flash(model, Temperature::new(t)?, pressure, options, false)?.enthalpy()
                - target_enthalpy,
        )
    };
    let (mut fa, mut fb) = (split(a)?, split(b)?);
    let scale = target_enthalpy.abs().max(fb - fa).max(1.0);
    let tolerance = options.energy_tolerance * scale;
    if fa > tolerance || fb < -tolerance {
        return Err(outside());
    }
    let accept = |t: f64, iterations: usize, diagnostics| {
        let flash = activity_tp_flash(model, Temperature::new(t)?, pressure, options)?;
        finish(
            t,
            flash.phases.clone(),
            iterations,
            Some(flash),
            diagnostics,
            scale,
        )
    };
    if fa.abs() <= tolerance {
        return accept(a, 0, diagnostics);
    }
    if fb.abs() <= tolerance {
        return accept(b, 0, diagnostics);
    }
    let mut side = 0;
    for iteration in 1..=options.max_temperature_iterations {
        let secant = (a * fb - b * fa) / (fb - fa);
        let c = if secant > a && secant < b {
            secant
        } else {
            a + (b - a) / 2.0
        };
        if c <= a || c >= b {
            return Err(Error::NonConvergence(format!(
                "PH two-phase temperature reached floating-point resolution at {c} K without meeting the energy tolerance"
            )));
        }
        let fc = split(c)?;
        if fc.abs() <= tolerance {
            return accept(c, iteration, diagnostics);
        }
        if fc < 0.0 {
            a = c;
            fa = fc;
            if side == -1 {
                fb /= 2.0;
            }
            side = -1;
        } else {
            b = c;
            fb = fc;
            if side == 1 {
                fa /= 2.0;
            }
            side = 1;
        }
    }
    Err(Error::NonConvergence(format!(
        "PH two-phase temperature iteration did not reach the energy tolerance within {} iterations",
        options.max_temperature_iterations
    )))
}
