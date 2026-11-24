//! Gamma–phi PH and PS flashes on the option temperature bracket. Both use one
//! bracketed structure: liquid and vapor branches split at the feed bubble and dew
//! temperatures, a two-phase mixture range, and a lever rule at the saturation
//! temperature of a single active component.

use thermo_models::{ActivityModel, ActivityPhaseState};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, InputError, Phase, Pressure,
    Temperature,
};

use super::flash::{activity_tp_flash, phase_state, require_stable_liquid, tp_flash};
use super::saturation::{
    activity_bubble_temperature, activity_dew_temperature, bubble_at, dew_at, solve_temperature,
};
use super::{
    ActivityFlashPhase, ActivityOptions, ActivityPhFlash, ActivityPsFlash, ActivitySaturation,
    ActivityTpFlash, SaturationKind, active, check_pressure, diagnostic, flag_inferred_interval,
    pures,
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

/// The caloric quantity a flash specifies with pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Spec {
    Enthalpy,
    Entropy,
}

impl Spec {
    fn of(self, s: &ActivityPhaseState) -> f64 {
        match self {
            Self::Enthalpy => s.enthalpy,
            Self::Entropy => s.entropy,
        }
    }
    /// Temperature derivative at fixed p and composition: c_p or c_p/T.
    fn slope(self, s: &ActivityPhaseState) -> f64 {
        match self {
            Self::Enthalpy => s.cp,
            Self::Entropy => s.cp / s.temperature_k,
        }
    }
    fn mixture(self, phases: &[ActivityFlashPhase]) -> f64 {
        phases
            .iter()
            .map(|ph| ph.phase_fraction * self.of(&ph.state))
            .sum()
    }
    fn flash(self) -> &'static str {
        match self {
            Self::Enthalpy => "PH",
            Self::Entropy => "PS",
        }
    }
    fn residual(self) -> &'static str {
        match self {
            Self::Enthalpy => "PH energy",
            Self::Entropy => "PS entropy",
        }
    }
    fn unit(self) -> &'static str {
        match self {
            Self::Enthalpy => "J/mol",
            Self::Entropy => "J/(mol K)",
        }
    }
}

/// Result of the shared PH/PS structure before it is named for its specification.
struct CaloricFlash {
    temperature_k: f64,
    pressure_pa: f64,
    value: f64,
    scaled_residual: f64,
    scale: f64,
    phases: Vec<ActivityFlashPhase>,
    bubble_temperature_k: Option<f64>,
    dew_temperature_k: Option<f64>,
    iterations: usize,
    tp: Option<ActivityTpFlash>,
    inferred_regression_interval: bool,
    diagnostics: Vec<Diagnostic>,
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
    let r = caloric_flash(model, pressure, target_enthalpy, options, Spec::Enthalpy)?;
    Ok(ActivityPhFlash {
        temperature_k: r.temperature_k,
        pressure_pa: r.pressure_pa,
        target_enthalpy,
        enthalpy: r.value,
        scaled_energy_residual: r.scaled_residual,
        energy_scale: r.scale,
        phases: r.phases,
        bubble_temperature_k: r.bubble_temperature_k,
        dew_temperature_k: r.dew_temperature_k,
        iterations: r.iterations,
        tp: r.tp,
        inferred_regression_interval: r.inferred_regression_interval,
        diagnostics: r.diagnostics,
    })
}

/// Pressure/entropy gamma–phi flash with the structure of [`activity_ph_flash`]:
/// the same branches, saturation path and solvers, with molar entropy as the
/// specification, ∂s/∂T = c_p/T on single-phase branches, and the entropy scales
/// max(|s_target|, 1 J/(mol K)) and max(|s_target|, two-phase entropy span,
/// 1 J/(mol K)). `energy_tolerance` bounds the scaled entropy residual.
///
/// The gamma–phi liquid has no pressure dependence of enthalpy or entropy (no
/// Poynting term), so an isentropic liquid compression has zero work in this model.
pub fn activity_ps_flash(
    model: &ActivityModel,
    pressure: Pressure,
    target_entropy: f64,
    options: ActivityOptions,
) -> Result<ActivityPsFlash, Error> {
    let r = caloric_flash(model, pressure, target_entropy, options, Spec::Entropy)?;
    let enthalpy = Spec::Enthalpy.mixture(&r.phases);
    Ok(ActivityPsFlash {
        temperature_k: r.temperature_k,
        pressure_pa: r.pressure_pa,
        target_entropy,
        entropy: r.value,
        enthalpy,
        scaled_entropy_residual: r.scaled_residual,
        entropy_scale: r.scale,
        phases: r.phases,
        bubble_temperature_k: r.bubble_temperature_k,
        dew_temperature_k: r.dew_temperature_k,
        iterations: r.iterations,
        tp: r.tp,
        inferred_regression_interval: r.inferred_regression_interval,
        diagnostics: r.diagnostics,
    })
}

fn caloric_flash(
    model: &ActivityModel,
    pressure: Pressure,
    target: f64,
    options: ActivityOptions,
    spec: Spec,
) -> Result<CaloricFlash, Error> {
    if !target.is_finite() {
        return Err(InputError::NonFinite(match spec {
            Spec::Enthalpy => "target enthalpy",
            Spec::Entropy => "target entropy",
        })
        .into());
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
    let (flash, residual, unit) = (spec.flash(), spec.residual(), spec.unit());
    let mut diagnostics = vec![match spec {
        Spec::Enthalpy => diagnostic(
            "A7-PH-GAMMA-PHI",
            "PH energy uses phase-weighted gamma-phi enthalpy (NASA7 ideal gas; saturation-pressure slope for pure liquids; NRTL excess); single-phase scale max(|h_target|, 1 J/mol), two-phase mixture scale max(|h_target|, two-phase enthalpy span, 1 J/mol).",
            DiagnosticSeverity::Information,
        ),
        Spec::Entropy => diagnostic(
            "A14-PS-GAMMA-PHI",
            "PS entropy uses phase-weighted gamma-phi entropy (NASA7 ideal gas with ideal mixing; Clausius-Clapeyron pure liquids without pressure dependence; NRTL excess); single-phase scale max(|s_target|, 1 J/(mol K)), two-phase mixture scale max(|s_target|, two-phase entropy span, 1 J/(mol K)).",
            DiagnosticSeverity::Information,
        ),
    }];
    if tb.is_none() || td.is_none() {
        diagnostics.push(diagnostic(
            match spec {
                Spec::Enthalpy => "A7-PH-SATURATION-OUTSIDE-BRACKET",
                Spec::Entropy => "A14-PS-SATURATION-OUTSIDE-BRACKET",
            },
            format!(
                "Feed bubble and/or dew temperature at {p} Pa lies outside {t_min}..={t_max} K; only the phase branches inside the bracket are searched."
            ),
            DiagnosticSeverity::Information,
        ));
    }
    let finish = |t: f64,
                  phases: Vec<ActivityFlashPhase>,
                  iterations,
                  tp,
                  mut diagnostics: Vec<Diagnostic>,
                  scale: f64| {
        let value = spec.mixture(&phases);
        let scaled = (value - target) / scale;
        if !scaled.is_finite() || scaled.abs() > options.energy_tolerance {
            return Err(Error::NonConvergence(format!(
                "{residual} residual {scaled:e} exceeds tolerance"
            )));
        }
        let inferred_regression_interval = flag_inferred_interval(&phases, &mut diagnostics);
        Ok(CaloricFlash {
            temperature_k: t,
            pressure_pa: p,
            value,
            scaled_residual: scaled,
            scale,
            phases,
            bubble_temperature_k: tb,
            dew_temperature_k: td,
            iterations,
            tp,
            inferred_regression_interval,
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
    let scale = target.abs().max(1.0);
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
            target <= spec.of(&edge) + budget
        } else {
            target >= spec.of(&edge) - budget
        };
        if !inside {
            continue;
        }
        let (t, iterations, s) =
            if ((spec.of(&edge) - target) / scale).abs() <= options.energy_tolerance {
                (edge_t, 0, edge)
            } else {
                let (t, iterations) = solve_temperature(
                    lo,
                    hi,
                    options,
                    options.energy_tolerance,
                    |t| {
                        let s = state(t)?;
                        Ok(((spec.of(&s) - target) / scale, spec.slope(&s) / scale))
                    },
                    &format!("{flash} single-phase temperature"),
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
            "{flash} target {target} {unit} at {p} Pa is not reachable within {t_min}..={t_max} K"
        ))
    };
    if pure_feed {
        let Saturation::At(sat) = &bubble else {
            return Err(outside());
        };
        let t = sat.temperature_k;
        let lb = phase_state(model, Phase::Liquid, t, p, &z)?;
        let vb = phase_state(model, Phase::Vapor, t, p, &z)?;
        let latent = spec.of(&vb) - spec.of(&lb);
        if latent <= 0.0 || !latent.is_finite() {
            return Err(Error::NumericalFailure(match spec {
                Spec::Enthalpy => "nonpositive latent heat".into(),
                Spec::Entropy => "nonpositive entropy of vaporization".into(),
            }));
        }
        let beta = ((target - spec.of(&lb)) / latent).clamp(0.0, 1.0);
        diagnostics.push(match spec {
            Spec::Enthalpy => diagnostic(
                "A7-PURE-SATURATION-LEVER-RULE",
                "Single active component: saturation temperature from the vapor-pressure correlation; vapor fraction from the lever rule.",
                DiagnosticSeverity::Information,
            ),
            Spec::Entropy => diagnostic(
                "A14-PURE-SATURATION-LEVER-RULE",
                "Single active component: saturation temperature from the vapor-pressure correlation; vapor fraction from the entropy lever rule.",
                DiagnosticSeverity::Information,
            ),
        });
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
    // Two-phase mixture: the split enthalpy and entropy increase with T over [lo, hi].
    let (mut a, mut b) = (liquid_top.unwrap_or(t_min), vapor_bottom.unwrap_or(t_max));
    if a >= b {
        return Err(outside());
    }
    let split = |t: f64| -> Result<f64, Error> {
        Ok(spec.mixture(&tp_flash(model, Temperature::new(t)?, pressure, options, false)?.phases) - target)
    };
    let (mut fa, mut fb) = (split(a)?, split(b)?);
    let scale = target.abs().max(fb - fa).max(1.0);
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
                "{flash} two-phase temperature reached floating-point resolution at {c} K without meeting the {} tolerance",
                match spec {
                    Spec::Enthalpy => "energy",
                    Spec::Entropy => "entropy",
                }
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
        "{flash} two-phase temperature iteration did not reach the {} tolerance within {} iterations",
        match spec {
            Spec::Enthalpy => "energy",
            Spec::Entropy => "entropy",
        },
        options.max_temperature_iterations
    )))
}
