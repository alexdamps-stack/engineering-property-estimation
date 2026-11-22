//! Gamma–phi bubble and dew pressure and temperature.

use thermo_models::{ActivityModel, ExcessState, PureLiquidState};
use thermo_types::{CalculationError as Error, Pressure, Temperature};

use super::{ActivityOptions, ActivitySaturation, SaturationKind, check_pressure, excess, pures};
use crate::stability::normalized as normalize;

/// Bubble pressure at T of the model composition: p = Σ x_i γ_i p_i^sat (explicit).
pub fn activity_bubble_pressure(
    model: &ActivityModel,
    temperature: Temperature,
) -> Result<ActivitySaturation, Error> {
    let t = temperature.value();
    model.check_temperature(temperature)?;
    let x = model.fractions().to_vec();
    let pure = pures(model, t)?;
    let (p, y) = bubble_at(model, t, &x, &pure)?;
    check_pressure(p)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Bubble,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: residual(model, t, p, &x, &y, &pure)?,
        liquid: x,
        vapor: y,
        iterations: 0,
    })
}

pub(super) fn bubble_at(
    model: &ActivityModel,
    t: f64,
    x: &[f64],
    pure: &[PureLiquidState],
) -> Result<(f64, Vec<f64>), Error> {
    let ex = excess(model, t, x)?;
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

/// max_i |ln(y_i p) − ln(x_i γ_i p_i^sat)| over components present in both phases.
fn residual_with(
    ln_gamma: &[f64],
    pure: &[PureLiquidState],
    p: f64,
    x: &[f64],
    y: &[f64],
) -> Result<f64, Error> {
    let mut r: f64 = 0.0;
    for i in 0..x.len() {
        if x[i] > 0.0 && y[i] > 0.0 {
            r = r.max(
                (y[i].ln() + p.ln() - x[i].ln() - ln_gamma[i] - pure[i].vapor_pressure.ln()).abs(),
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

fn residual(
    model: &ActivityModel,
    t: f64,
    p: f64,
    x: &[f64],
    y: &[f64],
    pure: &[PureLiquidState],
) -> Result<f64, Error> {
    let ex = excess(model, t, x)?;
    residual_with(&ex.ln_activity_coefficients, pure, p, x, y)
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
    let t = temperature.value();
    let d = dew_at(model, t, &y, &pures(model, t)?, options)?;
    check_pressure(d.pressure)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Dew,
        temperature_k: t,
        pressure_pa: d.pressure,
        liquid: d.liquid,
        vapor: y,
        max_ln_fugacity_residual: d.residual,
        iterations: d.iterations,
    })
}

/// Converged dew point at T; `excess` is evaluated at the returned liquid.
pub(super) struct DewPoint {
    pub(super) pressure: f64,
    liquid: Vec<f64>,
    iterations: usize,
    residual: f64,
    excess: ExcessState,
}

pub(super) fn dew_at(
    model: &ActivityModel,
    t: f64,
    y: &[f64],
    pure: &[PureLiquidState],
    options: ActivityOptions,
) -> Result<DewPoint, Error> {
    let n = y.len();
    let incipient = |ln_gamma: &[f64]| {
        let terms: Vec<f64> = (0..n)
            .map(|i| {
                if y[i] > 0.0 {
                    y[i] / (ln_gamma[i].exp() * pure[i].vapor_pressure)
                } else {
                    0.0
                }
            })
            .collect();
        (1.0 / terms.iter().sum::<f64>(), normalize(terms))
    };
    let mut ln_gamma: Vec<f64> = vec![0.0; n];
    for iteration in 1..=options.max_iterations {
        let (p, x) = incipient(&ln_gamma);
        let ex = excess(model, t, &x)?;
        let r = (0..n)
            .filter(|i| y[*i] > 0.0)
            .map(|i| (ex.ln_activity_coefficients[i] - ln_gamma[i]).abs())
            .fold(0.0, f64::max);
        if !r.is_finite() || !p.is_finite() {
            return Err(Error::NumericalFailure("nonfinite dew iteration".into()));
        }
        ln_gamma = ex.ln_activity_coefficients;
        if r <= options.tolerance {
            // Final consistent pressure and liquid from the converged γ.
            let (pressure, liquid) = incipient(&ln_gamma);
            let excess = excess(model, t, &liquid)?;
            let residual =
                residual_with(&excess.ln_activity_coefficients, pure, pressure, &liquid, y)?;
            return Ok(DewPoint {
                pressure,
                liquid,
                iterations: iteration,
                residual,
                excess,
            });
        }
    }
    Err(Error::NonConvergence(format!(
        "dew-pressure composition iteration exhausted {} iterations",
        options.max_iterations
    )))
}

/// Safeguarded Newton on f(T) = 0 with increasing f over the bracket.
pub(super) fn solve_temperature(
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
        let ex = excess(model, t, &x)?;
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
    let pure = pures(model, t)?;
    let (pb, y) = bubble_at(model, t, &x, &pure)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Bubble,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: residual(model, t, p, &x, &y, &pure)?.max((pb / p).ln().abs()),
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
    let mut last: Option<(f64, Vec<PureLiquidState>, DewPoint)> = None;
    let mut f = |t: f64| -> Result<(f64, f64), Error> {
        let pure = pures(model, t)?;
        let d = dew_at(model, t, &y, &pure, options)?;
        let slope = (0..y.len())
            .map(|i| {
                d.liquid[i]
                    * (d.excess.d_ln_activity_coefficients_dt[i] + pure[i].d_ln_vapor_pressure_dt)
            })
            .sum::<f64>();
        let value = d.pressure.ln() - p.ln();
        last = Some((t, pure, d));
        Ok((value, slope))
    };
    let (t, outer) = solve_temperature(
        options.minimum_temperature_k,
        options.maximum_temperature_k,
        options,
        options.tolerance,
        &mut f,
        "dew temperature",
    )?;
    let (pure, d) = match last {
        Some((at, pure, d)) if at == t => (pure, d),
        _ => {
            let pure = pures(model, t)?;
            let d = dew_at(model, t, &y, &pure, options)?;
            (pure, d)
        }
    };
    let r = residual_with(&d.excess.ln_activity_coefficients, &pure, p, &d.liquid, &y)?;
    Ok(ActivitySaturation {
        kind: SaturationKind::Dew,
        temperature_k: t,
        pressure_pa: p,
        max_ln_fugacity_residual: r.max((d.pressure / p).ln().abs()),
        liquid: d.liquid,
        vapor: y,
        iterations: outer + d.iterations,
    })
}
