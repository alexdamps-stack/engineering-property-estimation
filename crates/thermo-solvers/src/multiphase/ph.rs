//! Multiphase PH flash: a bracketed search on h(T) of the multiphase TP flash,
//! with the lever rule at an invariant temperature.

use thermo_models::{ActivityModel, ActivityPhaseState, CubicModel, HomogeneousState};
use thermo_types::{CalculationError as Error, DiagnosticSeverity, Pressure, Temperature};

use super::adapters::{ActivityAt, CubicAt, VlleAdapter};
use super::core::{self, active, diagnostic, distance};
use super::{MultiphaseOptions, MultiphasePhFlash, MultiphasePhase, MultiphaseTpFlash};

/// Pressure/enthalpy flash of the cubic model's composition with up to three
/// phases, on the option temperature bracket (within 200..=600 K).
///
/// h(T) of the multiphase TP flash is continuous except at an invariant
/// temperature (a single component, or the three-phase temperature of a binary
/// at fixed pressure). The bracket is narrowed by Illinois regula falsi with a
/// bisection safeguard until |h − h_target| ≤ energy_tolerance · max(|h_target|,
/// 1 J/mol). If the bracket shrinks to a few ulp with an enthalpy step across
/// it, the phases from both sides are combined and their amounts follow from the
/// lever rule (component balances and the energy balance).
pub fn vlle_ph_flash(
    model: &CubicModel,
    pressure: Pressure,
    target_enthalpy: f64,
    options: MultiphaseOptions,
) -> Result<MultiphasePhFlash<HomogeneousState>, Error> {
    options.validate()?;
    check_bracket(options, thermo_models::CUBIC_TEMPERATURE_RANGE_K, "cubic")?;
    let z = model.fractions();
    ph_flash(
        |t| Ok(CubicAt::new(model, Temperature::new(t)?, pressure)),
        z,
        target_enthalpy,
        options,
    )
}

/// Pressure/enthalpy gamma–phi flash with up to three phases (see [`vlle_ph_flash`]),
/// on the option bracket within 273.16..=450 K.
pub fn activity_vlle_ph_flash(
    model: &ActivityModel,
    pressure: Pressure,
    target_enthalpy: f64,
    options: MultiphaseOptions,
) -> Result<MultiphasePhFlash<ActivityPhaseState>, Error> {
    options.validate()?;
    check_bracket(
        options,
        thermo_models::ACTIVITY_TEMPERATURE_RANGE_K,
        "gamma-phi",
    )?;
    let z = model.fractions();
    ph_flash(
        |t| ActivityAt::new(model, Temperature::new(t)?, pressure),
        z,
        target_enthalpy,
        options,
    )
}

fn check_bracket(options: MultiphaseOptions, range: [f64; 2], name: &str) -> Result<(), Error> {
    let [lo, hi] = range;
    if options.minimum_temperature_k < lo || options.maximum_temperature_k > hi {
        return Err(Error::UnsupportedCapability(format!(
            "{name} multiphase PH bracket must lie within {lo}..={hi} K"
        )));
    }
    Ok(())
}

fn flash_enthalpy<S>(flash: &MultiphaseTpFlash<S>, h: impl Fn(&S) -> f64) -> f64 {
    flash
        .phases
        .iter()
        .map(|p| p.phase_fraction * h(&p.state))
        .sum()
}

fn ph_flash<A: VlleAdapter>(
    make: impl Fn(f64) -> Result<A, Error>,
    z: &[f64],
    target: f64,
    options: MultiphaseOptions,
) -> Result<MultiphasePhFlash<A::State>, Error> {
    if !target.is_finite() {
        return Err(thermo_types::InputError::NonFinite("target enthalpy").into());
    }
    let scale = target.abs().max(1.0);
    let tolerance = options.energy_tolerance * scale;
    let evaluate = |t: f64| -> Result<(MultiphaseTpFlash<A::State>, f64), Error> {
        let flash = core::tp_flash(&make(t)?, z, options)?;
        let h = flash_enthalpy(&flash, A::enthalpy);
        Ok((flash, h - target))
    };
    let (mut a, mut b) = (options.minimum_temperature_k, options.maximum_temperature_k);
    let (mut flash_a, mut fa) = evaluate(a)?;
    let (mut flash_b, mut fb) = evaluate(b)?;
    if fa > tolerance || fb < -tolerance {
        return Err(Error::OutsideValidityDomain(format!(
            "PH target enthalpy {target} J/mol outside [{}, {}] J/mol of the {a}..={b} K bracket",
            fa + target,
            fb + target
        )));
    }
    let mut side = 0i8;
    let mut last_width = b - a;
    for iteration in 1..=options.max_temperature_iterations {
        for (t, f, flash) in [(a, fa, &flash_a), (b, fb, &flash_b)] {
            if f.abs() <= tolerance {
                return Ok(solution(t, target, f / scale, iteration, flash.clone()));
            }
        }
        if b - a <= 4.0 * f64::EPSILON * b {
            return invariant::<A>(z, target, scale, options, a, b, flash_a, flash_b, iteration);
        }
        // Illinois regula falsi; bisect when the bracket stopped halving.
        let (mut wa, mut wb) = (fa, fb);
        if side == 1 {
            wb *= 0.5;
        } else if side == -1 {
            wa *= 0.5;
        }
        let mut t = a - wa * (b - a) / (wb - wa);
        if (iteration % 3 == 0 && (b - a) > 0.5 * last_width) || !(t > a && t < b) {
            t = 0.5 * (a + b);
        }
        if iteration % 3 == 0 {
            last_width = b - a;
        }
        let (flash, f) = evaluate(t)?;
        if f < 0.0 {
            a = t;
            fa = f;
            flash_a = flash;
            side = -1;
        } else {
            b = t;
            fb = f;
            flash_b = flash;
            side = 1;
        }
    }
    Err(Error::NonConvergence(format!(
        "multiphase PH exhausted {} temperature iterations",
        options.max_temperature_iterations
    )))
}

fn solution<S>(
    t: f64,
    target: f64,
    scaled: f64,
    iterations: usize,
    tp: MultiphaseTpFlash<S>,
) -> MultiphasePhFlash<S>
where
    S: Clone,
{
    let mut diagnostics = tp.diagnostics.clone();
    diagnostics.push(diagnostic(
        "VLLE-PH-BRACKET",
        "Temperature from a bracketed search on h(T) of the multiphase TP flash.",
        DiagnosticSeverity::Information,
    ));
    MultiphasePhFlash {
        temperature_k: t,
        pressure_pa: tp.pressure_pa,
        target_enthalpy: target,
        enthalpy: target + scaled * target.abs().max(1.0),
        scaled_energy_residual: scaled,
        phases: tp.phases.clone(),
        invariant_temperature: false,
        iterations,
        max_ln_fugacity_residual: tp.max_ln_fugacity_residual,
        max_component_balance_residual: tp.max_component_balance_residual,
        tp,
        diagnostics,
    }
}

/// Lever rule over the phases on both sides of an enthalpy step at T*.
#[allow(clippy::too_many_arguments)]
fn invariant<A: VlleAdapter>(
    z: &[f64],
    target: f64,
    scale: f64,
    options: MultiphaseOptions,
    t_lo: f64,
    t_hi: f64,
    lower: MultiphaseTpFlash<A::State>,
    upper: MultiphaseTpFlash<A::State>,
    iterations: usize,
) -> Result<MultiphasePhFlash<A::State>, Error> {
    let t = 0.5 * (t_lo + t_hi);
    // Phases on both sides keep their own states: at T* a minimum-Gibbs
    // re-evaluation could not tell coexisting roots of one composition apart.
    let mut candidates: Vec<MultiphasePhase<A::State>> = Vec::new();
    for p in lower.phases.iter().chain(&upper.phases) {
        let same = |c: &MultiphasePhase<A::State>| {
            c.kind == p.kind
                && c.phase == p.phase
                && distance(&c.composition, &p.composition) < 1e-6
                && (A::molar_volume(&c.state) / A::molar_volume(&p.state) - 1.0).abs() < 1e-6
        };
        if !candidates.iter().any(same) {
            candidates.push(p.clone());
        }
    }
    if candidates.len() > core::MAX_PHASES {
        return Err(Error::UnsupportedCapability(
            "invariant PH state with more than three phases".into(),
        ));
    }
    let h: Vec<f64> = candidates.iter().map(|c| A::enthalpy(&c.state)).collect();
    let span = h.iter().fold(f64::NEG_INFINITY, |m, v| m.max(*v))
        - h.iter().fold(f64::INFINITY, |m, v| m.min(*v));
    let act = active(z);
    let np = candidates.len();
    // Rows: component balances, then the energy balance scaled by the step.
    let mut rows: Vec<(Vec<f64>, f64)> = act
        .iter()
        .map(|&i| (candidates.iter().map(|c| c.composition[i]).collect(), z[i]))
        .collect();
    let energy_scale = span.max(1.0);
    rows.push((
        h.iter().map(|v| v / energy_scale).collect(),
        target / energy_scale,
    ));
    let normal: Vec<Vec<f64>> = (0..np)
        .map(|k| {
            (0..np)
                .map(|l| rows.iter().map(|(r, _)| r[k] * r[l]).sum())
                .collect()
        })
        .collect();
    let rhs: Vec<f64> = (0..np)
        .map(|k| rows.iter().map(|(r, b)| r[k] * b).sum())
        .collect();
    let beta = core::solve_dense(normal, rhs).ok_or_else(|| {
        Error::NumericalFailure("singular lever-rule system at an invariant temperature".into())
    })?;
    if beta.iter().any(|b| *b < -1e-9) {
        return Err(Error::NumericalFailure(format!(
            "lever rule at the invariant temperature {t} K gives a negative phase amount {beta:?}"
        )));
    }
    let beta: Vec<f64> = beta.iter().map(|b| b.max(0.0)).collect();
    let enthalpy: f64 = beta.iter().zip(&h).map(|(b, h)| b * h).sum();
    let scaled = (enthalpy - target) / scale;
    let balance = (0..z.len())
        .map(|i| {
            (candidates
                .iter()
                .zip(&beta)
                .map(|(c, b)| b * c.composition[i])
                .sum::<f64>()
                - z[i])
                .abs()
        })
        .fold(0.0, f64::max);
    if scaled.abs() > options.energy_tolerance || balance > 1e-9 {
        return Err(Error::NonConvergence(format!(
            "PH target lies in an enthalpy step at {t} K, but the lever rule over {np} phases leaves energy residual {scaled:e} and component residual {balance:e}"
        )));
    }
    let kept: Vec<MultiphasePhase<A::State>> = candidates
        .into_iter()
        .zip(beta)
        .filter(|(_, b)| *b > 0.0)
        .map(|(mut c, b)| {
            c.phase_fraction = b;
            c
        })
        .collect();
    let residual = act
        .iter()
        .map(|&i| {
            let mu: Vec<f64> = kept
                .iter()
                .map(|c| c.composition[i].ln() + A::state_ln_phi(&c.state)[i])
                .collect();
            mu.iter().fold(f64::NEG_INFINITY, |m, v| m.max(*v))
                - mu.iter().fold(f64::INFINITY, |m, v| m.min(*v))
        })
        .fold(0.0, f64::max);
    if residual.is_nan() || residual > 1e-7 {
        return Err(Error::NumericalFailure(format!(
            "phases on both sides of the enthalpy step at {t} K are not in equilibrium (ln fugacity residual {residual:e})"
        )));
    }
    let phases = core::order::<A>(kept);
    let mut diagnostics = lower.diagnostics.clone();
    diagnostics.push(diagnostic(
        "VLLE-PH-INVARIANT",
        format!("PH target lies in the enthalpy step at the invariant temperature {t} K; phase amounts from the lever rule over the phases on both sides of the step."),
        DiagnosticSeverity::Information,
    ));
    Ok(MultiphasePhFlash {
        temperature_k: t,
        pressure_pa: lower.pressure_pa,
        target_enthalpy: target,
        enthalpy,
        scaled_energy_residual: scaled,
        phases,
        invariant_temperature: true,
        iterations,
        max_ln_fugacity_residual: residual,
        max_component_balance_residual: balance,
        tp: lower,
        diagnostics,
    })
}
