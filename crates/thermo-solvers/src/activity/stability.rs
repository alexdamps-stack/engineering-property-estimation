//! Finite NRTL tangent-plane trials of a liquid (liquid–liquid split check).

use thermo_models::ActivityModel;
use thermo_types::{CalculationError as Error, Temperature};

use super::{ActivityOptions, LiquidStability, LiquidStabilityTrial, active};

/// Finite NRTL tangent-plane trials of a liquid at T: one near-pure start per active
/// component. A TPD below the threshold indicates a liquid–liquid split.
pub fn liquid_stability(
    model: &ActivityModel,
    temperature: Temperature,
    x: &[f64],
    options: ActivityOptions,
) -> Result<LiquidStability, Error> {
    options.validate()?;
    let nrtl = model.nrtl();
    let ex = nrtl.excess_at(temperature, x)?;
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
            let g = nrtl.excess_at(temperature, &w)?.ln_activity_coefficients;
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
