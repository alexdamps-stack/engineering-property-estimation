//! Gamma–phi (NRTL + ideal-gas vapor) equilibrium: bubble/dew pressure and
//! temperature, TP flash with an NRTL liquid tangent-plane check, and PH flash.
//!
//! With an ideal-gas vapor, liquid fugacities independent of pressure and a single
//! liquid phase, the feed is liquid for p ≥ p_bubble(T, z), vapor for
//! p ≤ p_dew(T, z), and two-phase in between. A liquid that fails the finite NRTL
//! tangent-plane trials (liquid–liquid split) is reported as unsupported; no VLLE.
use crate::rachford_rice::{RACHFORD_RICE_MAX_ITERATIONS, rachford_rice_negative};
use thermo_models::{
    ACTIVITY_TEMPERATURE_RANGE_K, ActivityModel, ActivityPhaseState, PureLiquidState,
};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, InputError, Phase, Pressure,
    Temperature,
};

/// Iteration limits and tolerances of the gamma–phi solvers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivityOptions {
    /// Successive-substitution limit (dew composition, TP split, TPD trials).
    pub max_iterations: usize,
    /// Convergence: max_i |ln(y_i p) − ln(x_i γ_i p_i^sat)| over active components.
    pub tolerance: f64,
    /// Outer temperature iteration limit (bubble/dew temperature, PH).
    pub max_temperature_iterations: usize,
    /// PH: |h − h_target| / max(|h_target|, 1 J/mol), in (0, 1e-6].
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
            tolerance: 1e-11,
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
    pub scaled_energy_residual: f64,
    /// Vapor first, then liquid. Pure saturation endpoints retain a zero-fraction phase.
    pub phases: Vec<ActivityFlashPhase>,
    pub bubble_temperature_k: f64,
    pub dew_temperature_k: f64,
    pub iterations: usize,
    /// TP evidence for a two-phase mixture result.
    pub tp: Option<ActivityTpFlash>,
    pub diagnostics: Vec<Diagnostic>,
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

fn normalize(mut x: Vec<f64>) -> Vec<f64> {
    let s: f64 = x.iter().sum();
    x.iter_mut().for_each(|v| *v /= s);
    x
}

/// Pure-liquid states of every component at T.
fn pures(model: &ActivityModel, t: f64) -> Result<Vec<PureLiquidState>, Error> {
    let temperature = Temperature::new(t)?;
    (0..model.fractions().len())
        .map(|i| model.pure_liquid(i, temperature))
        .collect()
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

/// Bubble pressure at T of the model composition: p = Σ x_i γ_i p_i^sat (explicit).
pub fn activity_bubble_pressure(
    model: &ActivityModel,
    temperature: Temperature,
) -> Result<ActivitySaturation, Error> {
    let t = temperature.value();
    model.check_temperature(temperature)?;
    let x = model.fractions().to_vec();
    let (p, y) = bubble_at(model, t, &x)?;
    check_pressure(p)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Bubble,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: residual(model, t, p, &x, &y)?,
        liquid: x,
        vapor: y,
        iterations: 0,
    })
}

fn bubble_at(model: &ActivityModel, t: f64, x: &[f64]) -> Result<(f64, Vec<f64>), Error> {
    let ex = model.nrtl().excess(t, x)?;
    let pure = pures(model, t)?;
    let partial: Vec<f64> = (0..x.len())
        .map(|i| {
            if x[i] > 0.0 {
                x[i] * ex.ln_activity_coefficients[i].exp() * pure[i].vapor_pressure
            } else {
                0.0
            }
        })
        .collect();
    let p: f64 = partial.iter().sum();
    if !p.is_finite() || p <= 0.0 {
        return Err(Error::NumericalFailure("nonfinite bubble pressure".into()));
    }
    Ok((p, partial.iter().map(|v| v / p).collect()))
}

/// max_i |ln(y_i p) − ln(x_i γ_i(x) p_i^sat)| over components present in both phases.
fn residual(model: &ActivityModel, t: f64, p: f64, x: &[f64], y: &[f64]) -> Result<f64, Error> {
    let ex = model.nrtl().excess(t, x)?;
    let pure = pures(model, t)?;
    let mut r: f64 = 0.0;
    for i in 0..x.len() {
        if x[i] > 0.0 && y[i] > 0.0 {
            r = r.max(
                (y[i].ln() + p.ln()
                    - x[i].ln()
                    - ex.ln_activity_coefficients[i]
                    - pure[i].vapor_pressure.ln())
                .abs(),
            );
        }
    }
    if !r.is_finite() {
        return Err(Error::NumericalFailure(
            "nonfinite fugacity residual".into(),
        ));
    }
    Ok(r)
}

/// Dew pressure at T: successive substitution on the incipient liquid composition.
pub fn activity_dew_pressure(
    model: &ActivityModel,
    temperature: Temperature,
    options: ActivityOptions,
) -> Result<ActivitySaturation, Error> {
    options.validate()?;
    model.check_temperature(temperature)?;
    let y = model.fractions().to_vec();
    let (p, x, iterations, r) = dew_at(model, temperature.value(), &y, options)?;
    check_pressure(p)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Dew,
        temperature_k: temperature.value(),
        pressure_pa: p,
        liquid: x,
        vapor: y,
        max_ln_fugacity_residual: r,
        iterations,
    })
}

fn dew_at(
    model: &ActivityModel,
    t: f64,
    y: &[f64],
    options: ActivityOptions,
) -> Result<(f64, Vec<f64>, usize, f64), Error> {
    let pure = pures(model, t)?;
    let n = y.len();
    let mut ln_gamma: Vec<f64> = vec![0.0; n];
    for iteration in 1..=options.max_iterations {
        let terms: Vec<f64> = (0..n)
            .map(|i| {
                if y[i] > 0.0 {
                    y[i] / (ln_gamma[i].exp() * pure[i].vapor_pressure)
                } else {
                    0.0
                }
            })
            .collect();
        let inv_p: f64 = terms.iter().sum();
        let p = 1.0 / inv_p;
        let x = normalize(terms);
        let ex = model.nrtl().excess(t, &x)?;
        let r = (0..n)
            .filter(|i| y[*i] > 0.0)
            .map(|i| (ex.ln_activity_coefficients[i] - ln_gamma[i]).abs())
            .fold(0.0, f64::max);
        if !r.is_finite() || !p.is_finite() {
            return Err(Error::NumericalFailure("nonfinite dew iteration".into()));
        }
        ln_gamma = ex.ln_activity_coefficients;
        if r <= options.tolerance {
            // Final consistent pressure from the converged γ.
            let terms: Vec<f64> = (0..n)
                .map(|i| {
                    if y[i] > 0.0 {
                        y[i] / (ln_gamma[i].exp() * pure[i].vapor_pressure)
                    } else {
                        0.0
                    }
                })
                .collect();
            let p = 1.0 / terms.iter().sum::<f64>();
            let x = normalize(terms);
            let residual = residual(model, t, p, &x, y)?;
            return Ok((p, x, iteration, residual));
        }
    }
    Err(Error::NonConvergence(format!(
        "dew-pressure composition iteration exhausted {} iterations",
        options.max_iterations
    )))
}

/// Safeguarded Newton on f(T) = 0 with increasing f over the bracket.
fn solve_temperature(
    mut lo: f64,
    mut hi: f64,
    options: ActivityOptions,
    tolerance: f64,
    mut f: impl FnMut(f64) -> Result<(f64, f64), Error>,
    what: &str,
) -> Result<(f64, usize), Error> {
    let (flo, _) = f(lo)?;
    let (fhi, _) = f(hi)?;
    if flo.abs() <= tolerance {
        return Ok((lo, 0));
    }
    if fhi.abs() <= tolerance {
        return Ok((hi, 0));
    }
    if flo > 0.0 || fhi < 0.0 {
        return Err(Error::OutsideValidityDomain(format!(
            "{what} not bracketed in {lo}..={hi} K"
        )));
    }
    let mut t = lo + (hi - lo) / 2.0;
    for iteration in 1..=options.max_temperature_iterations {
        let (value, slope) = f(t)?;
        if value.abs() <= tolerance {
            return Ok((t, iteration));
        }
        if value < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
        let newton = t - value / slope;
        let next = if slope > 0.0 && newton > lo && newton < hi {
            newton
        } else {
            lo + (hi - lo) / 2.0
        };
        if next == t || hi - lo <= 4.0 * f64::EPSILON * t {
            return Err(Error::NonConvergence(format!(
                "{what}: temperature bracket collapsed with residual {value:e}"
            )));
        }
        t = next;
    }
    Err(Error::NonConvergence(format!(
        "{what} exhausted {} temperature iterations",
        options.max_temperature_iterations
    )))
}

/// Bubble temperature at p: ln Σ x_i γ_i p_i^sat = ln p on the option bracket.
pub fn activity_bubble_temperature(
    model: &ActivityModel,
    pressure: Pressure,
    options: ActivityOptions,
) -> Result<ActivitySaturation, Error> {
    options.validate()?;
    let p = pressure.value();
    check_pressure(p)?;
    let x = model.fractions().to_vec();
    let f = |t: f64| -> Result<(f64, f64), Error> {
        let ex = model.nrtl().excess(t, &x)?;
        let pure = pures(model, t)?;
        let (mut sum, mut slope) = (0.0, 0.0);
        for i in 0..x.len() {
            if x[i] > 0.0 {
                let w = x[i] * ex.ln_activity_coefficients[i].exp() * pure[i].vapor_pressure;
                sum += w;
                slope += w * (ex.d_ln_activity_coefficients_dt[i] + pure[i].d_ln_vapor_pressure_dt);
            }
        }
        Ok((sum.ln() - p.ln(), slope / sum))
    };
    let (t, outer) = solve_temperature(
        options.minimum_temperature_k,
        options.maximum_temperature_k,
        options,
        options.tolerance,
        f,
        "bubble temperature",
    )?;
    let (pb, y) = bubble_at(model, t, &x)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Bubble,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: residual(model, t, p, &x, &y)?.max((pb / p).ln().abs()),
        liquid: x,
        vapor: y,
        iterations: outer,
    })
}

/// Dew temperature at p: ln p_dew(T) = ln p. The slope Σ x_i(∂lnγ_i/∂T + dln p_i^sat/dT)
/// is exact at fixed y by Gibbs–Duhem.
pub fn activity_dew_temperature(
    model: &ActivityModel,
    pressure: Pressure,
    options: ActivityOptions,
) -> Result<ActivitySaturation, Error> {
    options.validate()?;
    let p = pressure.value();
    check_pressure(p)?;
    let y = model.fractions().to_vec();
    let f = |t: f64| -> Result<(f64, f64), Error> {
        let (pd, x, _, _) = dew_at(model, t, &y, options)?;
        let ex = model.nrtl().excess(t, &x)?;
        let pure = pures(model, t)?;
        let slope = (0..x.len())
            .map(|i| x[i] * (ex.d_ln_activity_coefficients_dt[i] + pure[i].d_ln_vapor_pressure_dt))
            .sum::<f64>();
        Ok((pd.ln() - p.ln(), slope))
    };
    let (t, outer) = solve_temperature(
        options.minimum_temperature_k,
        options.maximum_temperature_k,
        options,
        options.tolerance,
        f,
        "dew temperature",
    )?;
    let (pd, x, it, _) = dew_at(model, t, &y, options)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Dew,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: residual(model, t, p, &x, &y)?.max((pd / p).ln().abs()),
        liquid: x,
        vapor: y,
        iterations: outer + it,
    })
}

/// Finite NRTL tangent-plane trials of a liquid at T: one near-pure start per active
/// component. A TPD below the threshold indicates a liquid–liquid split.
pub fn liquid_stability(
    model: &ActivityModel,
    temperature: Temperature,
    x: &[f64],
    options: ActivityOptions,
) -> Result<LiquidStability, Error> {
    options.validate()?;
    let t = temperature.value();
    let nrtl = model.nrtl();
    let ex = nrtl.excess(t, x)?;
    let act = active(x);
    let n = x.len();
    let mut result = LiquidStability {
        stable_in_trials: true,
        trials: Vec::new(),
    };
    if act.len() < 2 {
        return Ok(result);
    }
    let d: Vec<f64> = (0..n)
        .map(|i| {
            if x[i] > 0.0 {
                x[i].ln() + ex.ln_activity_coefficients[i]
            } else {
                0.0
            }
        })
        .collect();
    for (position, &j) in act.iter().enumerate() {
        let others = (act.len() - 1) as f64;
        let mut ln_w: Vec<f64> = (0..act.len())
            .map(|q| {
                if q == position {
                    (1.0 - 1e-3f64).ln()
                } else {
                    (1e-3 / others).ln()
                }
            })
            .collect();
        let mut done = None;
        for iteration in 1..=options.max_iterations {
            let big: Vec<f64> = ln_w.iter().map(|l| l.exp()).collect();
            let total: f64 = big.iter().sum();
            let mut w = vec![0.0; n];
            for (q, &i) in act.iter().enumerate() {
                w[i] = big[q] / total;
            }
            let w = normalize(w);
            let g = nrtl.excess(t, &w)?.ln_activity_coefficients;
            let mut tm = 1.0;
            let mut change: f64 = 0.0;
            let mut distance: f64 = 0.0;
            for (q, &i) in act.iter().enumerate() {
                tm += big[q] * (ln_w[q] + g[i] - d[i] - 1.0);
                let next = d[i] - g[i];
                change = change.max((next - ln_w[q]).abs());
                distance = distance.max((w[i] - x[i]).abs());
                ln_w[q] = next;
            }
            if !tm.is_finite() {
                return Err(Error::NumericalFailure("nonfinite liquid TPD".into()));
            }
            let trivial = distance <= 1e-6;
            if tm < options.tm_threshold && !trivial {
                done = Some((tm, w, iteration, false));
                break;
            }
            if change <= 1e-10 {
                done = Some((tm, w, iteration, trivial));
                break;
            }
        }
        let Some((tm, composition, iterations, trivial)) = done else {
            return Err(Error::NonConvergence(format!(
                "liquid TPD trial {j} exhausted {} iterations; liquid stability undetermined",
                options.max_iterations
            )));
        };
        if tm < options.tm_threshold && !trivial {
            result.stable_in_trials = false;
        }
        result.trials.push(LiquidStabilityTrial {
            component: j,
            tangent_plane_distance: tm,
            composition,
            iterations,
            trivial,
        });
    }
    Ok(result)
}

fn phase_state(
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

fn require_stable_liquid(
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
    model.check_domain(temperature, pressure)?;
    let (t, p) = (temperature.value(), pressure.value());
    let z = model.fractions().to_vec();
    let n = z.len();
    let act = active(&z);
    let (p_bubble, _) = bubble_at(model, t, &z)?;
    let (p_dew, _, _, _) = if act.len() == 1 {
        (p_bubble, z.clone(), 0, 0.0)
    } else {
        dew_at(model, t, &z, options)?
    };
    let mut diagnostics = vec![diagnostic(
        "A7-GAMMA-PHI",
        "Modified Raoult with ideal-gas vapor; single-phase classification from exact bubble/dew pressures of the feed; at most one liquid phase.",
        DiagnosticSeverity::Information,
    )];
    let single = |phase: Phase,
                  stability: Option<LiquidStability>,
                  diagnostics|
     -> Result<ActivityTpFlash, Error> {
        Ok(ActivityTpFlash {
            temperature_k: t,
            pressure_pa: p,
            phases: vec![ActivityFlashPhase {
                phase,
                phase_fraction: 1.0,
                composition: z.clone(),
                state: phase_state(model, phase, t, p, &z)?,
            }],
            bubble_pressure_pa: p_bubble,
            dew_pressure_pa: p_dew,
            iterations: 0,
            max_ln_fugacity_residual: 0.0,
            max_component_balance_residual: 0.0,
            liquid_stability: stability,
            diagnostics,
        })
    };
    if p >= p_bubble {
        let s = require_stable_liquid(model, t, &z, options)?;
        return single(Phase::Liquid, Some(s), diagnostics);
    }
    if p <= p_dew {
        return single(Phase::Vapor, None, diagnostics);
    }
    let z_active: Vec<f64> = act.iter().map(|&i| z[i]).collect();
    let pure = pures(model, t)?;
    let ln_ps: Vec<f64> = pure.iter().map(|s| s.vapor_pressure.ln()).collect();
    let mut gamma = model.nrtl().excess(t, &z)?.ln_activity_coefficients;
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
        let next = model.nrtl().excess(t, &x)?.ln_activity_coefficients;
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
    if !(beta > 0.0 && beta < 1.0) {
        return Err(Error::NumericalFailure(format!(
            "converged split has vapor fraction {beta} outside (0, 1) between dew and bubble pressure"
        )));
    }
    let balance = (0..n)
        .map(|i| (beta * y[i] + (1.0 - beta) * x[i] - z[i]).abs())
        .fold(0.0, f64::max);
    let stability = require_stable_liquid(model, t, &x, options)?;
    diagnostics.push(diagnostic(
        "A7-LIQUID-TPD-FINITE-TRIALS",
        "Liquid stability from near-pure NRTL tangent-plane trials; finite trials do not prove global stability.",
        DiagnosticSeverity::Information,
    ));
    Ok(ActivityTpFlash {
        temperature_k: t,
        pressure_pa: p,
        phases: vec![
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
        ],
        bubble_pressure_pa: p_bubble,
        dew_pressure_pa: p_dew,
        iterations,
        max_ln_fugacity_residual: r,
        max_component_balance_residual: balance,
        liquid_stability: Some(stability),
        diagnostics,
    })
}

/// Pressure/enthalpy gamma–phi flash. The feed's bubble and dew temperatures at p
/// must lie inside the option bracket. Liquid (h ≤ h_L(T_bubble)) and vapor
/// (h ≥ h_V(T_dew)) branches are solved by safeguarded Newton with exact c_p; the
/// two-phase branch by bisection on TP enthalpy, or by the lever rule at the
/// saturation temperature for a single active component.
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
    let scale = target_enthalpy.abs().max(1.0);
    let bubble = activity_bubble_temperature(model, pressure, options)?;
    let pure_feed = active(&z).len() == 1;
    let dew = if pure_feed {
        bubble.clone()
    } else {
        activity_dew_temperature(model, pressure, options)?
    };
    let (tb, td) = (bubble.temperature_k, dew.temperature_k);
    let h_liquid = |t: f64| phase_state(model, Phase::Liquid, t, p, &z);
    let h_vapor = |t: f64| phase_state(model, Phase::Vapor, t, p, &z);
    let lb = h_liquid(tb)?;
    let vd = h_vapor(td)?;
    let mut diagnostics = vec![diagnostic(
        "A7-PH-GAMMA-PHI",
        "PH energy uses phase-weighted gamma-phi enthalpy (NASA7 ideal gas; Clausius-Clapeyron pure liquids; NRTL excess); scale max(|h_target|, 1 J/mol).",
        DiagnosticSeverity::Information,
    )];
    let finish = |t: f64, phases: Vec<ActivityFlashPhase>, iterations, tp, diagnostics| {
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
    let budget = options.energy_tolerance * scale;
    // Single-phase branches.
    for (phase, lo, hi, edge, inside) in [
        (
            Phase::Liquid,
            options.minimum_temperature_k,
            tb,
            &lb,
            target_enthalpy <= lb.enthalpy + budget,
        ),
        (
            Phase::Vapor,
            td,
            options.maximum_temperature_k,
            &vd,
            target_enthalpy >= vd.enthalpy - budget,
        ),
    ] {
        if !inside {
            continue;
        }
        if ((edge.enthalpy - target_enthalpy) / scale).abs() <= options.energy_tolerance {
            let t = edge.temperature_k;
            return finish(t, one(phase, edge.clone()), 0, None, diagnostics);
        }
        let state = |t: f64| {
            if phase == Phase::Liquid {
                h_liquid(t)
            } else {
                h_vapor(t)
            }
        };
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
        if phase == Phase::Liquid {
            require_stable_liquid(model, t, &z, options)?;
        }
        return finish(t, one(phase, state(t)?), iterations, None, diagnostics);
    }
    if pure_feed {
        let t = tb;
        let vb = h_vapor(t)?;
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
        return finish(t, phases, bubble.iterations, None, diagnostics);
    }
    // Two-phase mixture: h(T) increases from h_L(T_b) to h_V(T_d).
    let (mut lo, mut hi) = (tb, td);
    for iteration in 1..=options.max_temperature_iterations {
        let mid = lo + (hi - lo) / 2.0;
        if mid <= lo || mid >= hi {
            break;
        }
        let flash = activity_tp_flash(model, Temperature::new(mid)?, pressure, options)?;
        let r = (flash.enthalpy() - target_enthalpy) / scale;
        if r.abs() <= options.energy_tolerance {
            let phases = flash.phases.clone();
            return finish(mid, phases, iteration, Some(flash), diagnostics);
        }
        if r < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(Error::NonConvergence(format!(
        "PH two-phase temperature iteration did not reach the energy tolerance within {} iterations",
        options.max_temperature_iterations
    )))
}
