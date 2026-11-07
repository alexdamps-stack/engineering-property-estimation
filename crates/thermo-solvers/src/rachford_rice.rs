//! Rachford–Rice phase split for prescribed equilibrium ratios. This is not a
//! thermodynamic flash: K values are inputs, not results.
use thermo_types::{CalculationError as Error, InputError};

/// Iteration limit used by the TP flash for each Rachford–Rice solve.
pub const RACHFORD_RICE_MAX_ITERATIONS: usize = 200;

/// Root of the Rachford–Rice equation. Vectors follow the input order.
#[derive(Debug, Clone, PartialEq)]
pub struct RachfordRiceSolution {
    /// Molar vapor fraction. The negative-flash solver may return values outside [0, 1].
    pub vapor_fraction: f64,
    pub liquid: Vec<f64>,
    pub vapor: Vec<f64>,
    /// |g(beta)| at the returned root, dimensionless.
    pub residual: f64,
    pub iterations: usize,
}

/// Constrained (0 <= beta <= 1) outcome for prescribed K values.
#[derive(Debug, Clone, PartialEq)]
pub enum RachfordRiceOutcome {
    /// g(0) <= 0: the prescribed K values place no vapor in [0, 1].
    LiquidOnly,
    /// g(1) >= 0: the prescribed K values place no liquid in [0, 1].
    VaporOnly,
    TwoPhase(RachfordRiceSolution),
}

/// g(beta) = sum z_i (K_i - 1) / (1 + beta (K_i - 1)); strictly decreasing between poles.
pub fn rachford_rice_residual(z: &[f64], k: &[f64], vapor_fraction: f64) -> f64 {
    z.iter()
        .zip(k)
        .filter(|(z, _)| **z > 0.0)
        .map(|(z, k)| z * (k - 1.0) / (1.0 + vapor_fraction * (k - 1.0)))
        .sum()
}

fn validate(z: &[f64], k: &[f64]) -> Result<(), Error> {
    if z.is_empty() {
        return Err(InputError::EmptyComposition.into());
    }
    if z.len() != k.len() {
        return Err(InputError::ComponentMismatch.into());
    }
    if z.iter().chain(k).any(|v| !v.is_finite()) {
        return Err(InputError::NonFinite("feed fraction or K value").into());
    }
    if z.iter().any(|v| *v < 0.0) {
        return Err(InputError::Negative("feed fraction").into());
    }
    if k.iter().any(|v| *v <= 0.0) {
        return Err(InputError::NotPositive("K value").into());
    }
    if (z.iter().sum::<f64>() - 1.0).abs() > 1e-10 {
        return Err(InputError::NotNormalized.into());
    }
    Ok(())
}

/// Constrained Rachford–Rice: classify single-phase outcomes from g(0) and g(1),
/// otherwise solve for 0 < beta < 1.
pub fn rachford_rice(
    z: &[f64],
    k: &[f64],
    max_iterations: usize,
) -> Result<RachfordRiceOutcome, Error> {
    validate(z, k)?;
    if rachford_rice_residual(z, k, 0.0) <= 0.0 {
        return Ok(RachfordRiceOutcome::LiquidOnly);
    }
    if rachford_rice_residual(z, k, 1.0) >= 0.0 {
        return Ok(RachfordRiceOutcome::VaporOnly);
    }
    rachford_rice_negative(z, k, max_iterations).map(RachfordRiceOutcome::TwoPhase)
}

/// Negative-flash Rachford–Rice: solve g(beta)=0 between the poles
/// 1/(1-Kmax) < beta < 1/(1-Kmin), where phase compositions stay nonnegative.
/// Requires active (z>0) K values on both sides of one; otherwise no root exists.
/// Safeguarded Newton with bisection; `max_iterations` is a hard limit.
pub fn rachford_rice_negative(
    z: &[f64],
    k: &[f64],
    max_iterations: usize,
) -> Result<RachfordRiceSolution, Error> {
    validate(z, k)?;
    let active = || z.iter().zip(k).filter(|(z, _)| **z > 0.0).map(|(_, k)| *k);
    let k_max = active().fold(f64::NEG_INFINITY, f64::max);
    let k_min = active().fold(f64::INFINITY, f64::min);
    if k_max <= 1.0 || k_min >= 1.0 {
        return Err(Error::UnsupportedCapability(
            "no Rachford-Rice root: all active K values lie on one side of 1".into(),
        ));
    }
    let (mut lo, mut hi) = (1.0 / (1.0 - k_max), 1.0 / (1.0 - k_min));
    let mut beta = 0.5;
    let derivative = |beta: f64| -> f64 {
        -z.iter()
            .zip(k)
            .filter(|(z, _)| **z > 0.0)
            .map(|(z, k)| z * (k - 1.0).powi(2) / (1.0 + beta * (k - 1.0)).powi(2))
            .sum::<f64>()
    };
    for iteration in 1..=max_iterations {
        let g = rachford_rice_residual(z, k, beta);
        if !g.is_finite() {
            return Err(Error::NumericalFailure(
                "nonfinite Rachford-Rice residual".into(),
            ));
        }
        if g > 0.0 {
            lo = beta;
        } else {
            hi = beta;
        }
        let newton = beta - g / derivative(beta);
        let next = if newton > lo && newton < hi && newton.is_finite() {
            newton
        } else {
            lo + (hi - lo) / 2.0
        };
        let step = (next - beta).abs();
        if g == 0.0 || step <= 4.0 * f64::EPSILON * beta.abs().max(1.0) || next == lo || next == hi
        {
            let beta = if g == 0.0 { beta } else { next };
            return Ok(solution(z, k, beta, iteration));
        }
        beta = next;
    }
    Err(Error::NonConvergence(format!(
        "Rachford-Rice exhausted {max_iterations} iterations"
    )))
}

fn solution(z: &[f64], k: &[f64], beta: f64, iterations: usize) -> RachfordRiceSolution {
    let mut liquid: Vec<f64> = z
        .iter()
        .zip(k)
        .map(|(z, k)| {
            if *z > 0.0 {
                z / (1.0 + beta * (k - 1.0))
            } else {
                0.0
            }
        })
        .collect();
    let mut vapor: Vec<f64> = liquid.iter().zip(k).map(|(x, k)| x * k).collect();
    // At the root both sums are one; normalize only the last-bit roundoff.
    for phase in [&mut liquid, &mut vapor] {
        let sum: f64 = phase.iter().sum();
        phase.iter_mut().for_each(|v| *v /= sum);
    }
    RachfordRiceSolution {
        vapor_fraction: beta,
        liquid,
        vapor,
        residual: rachford_rice_residual(z, k, beta).abs(),
        iterations,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
    }

    #[test]
    fn single_phase_classification_for_one_sided_k() {
        let z = [0.3, 0.7];
        assert_eq!(
            rachford_rice(&z, &[2.0, 1.5], 50).unwrap(),
            RachfordRiceOutcome::VaporOnly
        );
        assert_eq!(
            rachford_rice(&z, &[0.2, 0.9], 50).unwrap(),
            RachfordRiceOutcome::LiquidOnly
        );
        // K on both sides of 1 but g(0) <= 0: subcooled liquid.
        assert_eq!(
            rachford_rice(&z, &[1.1, 0.5], 50).unwrap(),
            RachfordRiceOutcome::LiquidOnly
        );
        assert!(matches!(
            rachford_rice_negative(&z, &[2.0, 1.5], 50),
            Err(Error::UnsupportedCapability(_))
        ));
    }

    #[test]
    fn binary_matches_closed_form() {
        let (z1, k1, k2) = (0.4, 3.0, 0.2);
        let x1 = (1.0 - k2) / (k1 - k2);
        let beta = (z1 - x1) / (k1 * x1 - x1);
        let RachfordRiceOutcome::TwoPhase(s) =
            rachford_rice(&[z1, 1.0 - z1], &[k1, k2], 100).unwrap()
        else {
            panic!("expected two phases")
        };
        close(s.vapor_fraction, beta, 1e-15);
        close(s.liquid[0], x1, 1e-15);
        close(s.vapor[0], k1 * x1, 1e-15);
        assert!(s.residual <= 1e-15);
    }

    #[test]
    fn multicomponent_balances_and_monotone_residual() {
        let z = [0.5, 0.2, 0.15, 0.1, 0.05];
        let k = [6.0, 1.8, 0.7, 0.2, 0.01];
        let RachfordRiceOutcome::TwoPhase(s) = rachford_rice(&z, &k, 100).unwrap() else {
            panic!("expected two phases")
        };
        assert!(s.vapor_fraction > 0.0 && s.vapor_fraction < 1.0);
        for i in 0..z.len() {
            let balance = s.vapor_fraction * s.vapor[i] + (1.0 - s.vapor_fraction) * s.liquid[i];
            close(balance, z[i], 1e-15);
            close(s.vapor[i] / s.liquid[i], k[i], 1e-13 * k[i]);
        }
        let (lo, hi) = (1.0 / (1.0 - 6.0), 1.0 / (1.0 - 0.01));
        let mut previous = f64::INFINITY;
        for step in 1..1000 {
            let beta = lo + (hi - lo) * step as f64 / 1000.0;
            let g = rachford_rice_residual(&z, &k, beta);
            assert!(g < previous, "residual not decreasing at beta={beta}");
            previous = g;
        }
    }

    #[test]
    fn negative_flash_outside_unit_interval() {
        // Vapor-side negative flash: g(1) > 0 with K on both sides of one.
        let z = [0.9, 0.1];
        let k = [1.5, 0.5];
        assert_eq!(
            rachford_rice(&z, &k, 100).unwrap(),
            RachfordRiceOutcome::VaporOnly
        );
        let s = rachford_rice_negative(&z, &k, 100).unwrap();
        assert!(s.vapor_fraction > 1.0 && s.vapor_fraction < 1.0 / (1.0 - 0.5));
        assert!(s.residual <= 1e-14);
        assert!(s.liquid.iter().chain(&s.vapor).all(|v| *v >= 0.0));
        let x1 = (1.0 - k[1]) / (k[0] - k[1]);
        close(s.liquid[0], x1, 1e-14);
    }

    #[test]
    fn k_close_to_one_and_trace_components() {
        // K within 1e-7 of one: beta is ill-conditioned in the inputs (about 0.2 for
        // these nominal values), so compare with the closed form evaluated from the
        // same stored floats. The K=1 trace keeps x3 = y3 = z3:
        // x1 = (1-z3)(1-K2)/(K1-K2), beta = (z1-x1)/(x1 (K1-1)).
        let z = [0.5 + 1e-8, 0.5 - 1e-8 - 1e-12, 1e-12, 0.0];
        let k = [1.0 + 1e-7, 1.0 - 1e-7, 1.0, 50.0];
        let RachfordRiceOutcome::TwoPhase(s) = rachford_rice(&z, &k, 200).unwrap() else {
            panic!("expected two phases")
        };
        let x1 = (1.0 - z[2]) * (1.0 - k[1]) / (k[0] - k[1]);
        let beta = (z[0] - x1) / (x1 * (k[0] - 1.0));
        assert!((0.1..0.3).contains(&beta));
        close(s.vapor_fraction, beta, 1e-7 * beta);
        assert!(s.residual <= 1e-20);
        close(s.liquid[2], 1e-12, 1e-24);
        assert_eq!(s.liquid[3], 0.0);
        assert_eq!(s.vapor[3], 0.0);
    }

    #[test]
    fn invalid_inputs_and_iteration_exhaustion() {
        assert!(rachford_rice(&[], &[], 10).is_err());
        assert!(rachford_rice(&[0.5, 0.5], &[2.0], 10).is_err());
        assert!(rachford_rice(&[0.5, 0.5], &[f64::NAN, 0.5], 10).is_err());
        assert!(rachford_rice(&[0.5, 0.5], &[0.0, 0.5], 10).is_err());
        assert!(rachford_rice(&[1.5, -0.5], &[2.0, 0.5], 10).is_err());
        assert!(rachford_rice(&[0.5, 0.6], &[2.0, 0.5], 10).is_err());
        assert!(matches!(
            rachford_rice(&[0.5, 0.5], &[3.0, 0.2], 0),
            Err(Error::NonConvergence(_))
        ));
        assert!(matches!(
            rachford_rice(&[0.5, 0.5], &[3.0, 0.2], 1),
            Err(Error::NonConvergence(_))
        ));
    }
}
