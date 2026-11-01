//! Homogeneous PR76 and PR76 with Boston–Mathias supercritical continuation.
//! No phase stability, equilibrium flash, or ideal caloric contribution yet.

use serde::Serialize;
use thermo_types::{
    CubicModel, GAS_CONSTANT as R, Mixture, Pressure, RootSelection, Temperature, ThermoError,
};

pub const MODEL_VERSION: &str = "pr76-exact-constants-bm-v1";
const OMEGA_A: f64 = 0.457_235_528_921_382_2;
const OMEGA_B: f64 = 0.077_796_073_903_888_46;

#[derive(Debug, Serialize)]
pub struct HomogeneousProperties {
    pub roots_z: Vec<f64>,
    pub z: f64,
    pub molar_volume_m3_mol: f64,
    pub molar_density_mol_m3: f64,
    pub mass_density_kg_m3: f64,
    pub ln_fugacity_coefficients: Vec<f64>,
    /// Departure from an ideal gas of the same T, P, and composition.
    pub residual_enthalpy_j_mol: f64,
    /// Departure from an ideal gas of the same T, P, and composition.
    pub residual_entropy_j_mol_k: f64,
    pub pressure_relative_residual: f64,
}

/// Alpha and d(alpha)/dT. Kappa is the PR76 quadratic for both model variants.
pub fn alpha(
    model: CubicModel,
    temperature: Temperature,
    tc: Temperature,
    omega: f64,
) -> Result<(f64, f64), ThermoError> {
    if !omega.is_finite() || !(-0.5..=1.0).contains(&omega) {
        return Err(ThermoError::Unsupported(
            "initial alpha evaluator requires -0.5 <= omega <= 1.0".into(),
        ));
    }
    let t = temperature.value();
    let tr = t / tc.value();
    let k = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let (a, derivative) = if model == CubicModel::PengRobinsonBostonMathias && tr > 1.0 {
        let d = 1.0 + k / 2.0;
        let c = 1.0 - 1.0 / d;
        let a = (2.0 * c * (1.0 - tr.powf(d))).exp();
        (a, -2.0 * c * d * tr.powf(d - 1.0) / tc.value() * a)
    } else {
        let q = 1.0 + k * (1.0 - tr.sqrt());
        if q <= 0.0 {
            return Err(ThermoError::Unsupported(
                "PR alpha beyond its positive square-root branch".into(),
            ));
        }
        (q * q, -k * q / (tr.sqrt() * tc.value()))
    };
    if a <= 0.0 || !a.is_finite() || !derivative.is_finite() {
        return Err(ThermoError::NumericalFailure(
            "nonfinite or vanishing alpha".into(),
        ));
    }
    Ok((a, derivative))
}

/// Evaluate an explicitly selected homogeneous root, without assuming equilibrium.
pub fn evaluate(
    model: CubicModel,
    mixture: &Mixture,
    temperature: Temperature,
    pressure: Pressure,
    root: RootSelection,
) -> Result<HomogeneousProperties, ThermoError> {
    let t = temperature.value();
    let p = pressure.value();
    // Numerical guardrails, not an experimentally validated thermodynamic domain.
    if !(80.0..=1000.0).contains(&t) || !(1e3..=1e8).contains(&p) {
        return Err(ThermoError::Unsupported("initial evaluator guardrails: 80–1000 K and 1e3–1e8 Pa absolute; no full-domain accuracy claim".into()));
    }
    let cs = mixture.components();
    let x = mixture.mole_fractions();
    let n = cs.len();
    let mut ai = Vec::with_capacity(n);
    let mut dai = Vec::with_capacity(n);
    let mut bi = Vec::with_capacity(n);
    for c in cs {
        let tc = c.critical_temperature().value();
        let pc = c.critical_pressure().value();
        let (alpha, derivative) = alpha(
            model,
            temperature,
            c.critical_temperature(),
            c.acentric_factor(),
        )?;
        let a0 = OMEGA_A * (R * tc).powi(2) / pc;
        ai.push(a0 * alpha);
        dai.push(a0 * derivative);
        bi.push(OMEGA_B * R * tc / pc);
    }
    let b: f64 = bi.iter().zip(x).map(|(b, x)| b * x).sum();
    let mut a = 0.0;
    let mut da = 0.0;
    let mut row_sums = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            let aij = (ai[i] * ai[j]).sqrt() * (1.0 - mixture.binary_interactions()[i][j]);
            row_sums[i] += x[j] * aij;
            a += x[i] * x[j] * aij;
            da += x[i] * x[j] * aij * 0.5 * (dai[i] / ai[i] + dai[j] / ai[j]);
        }
    }
    let reduced_a = a * p / (R * t).powi(2);
    let reduced_b = b * p / (R * t);
    if ![a, b, da, reduced_a, reduced_b]
        .iter()
        .all(|v| v.is_finite())
        || a <= 0.0
        || b <= 0.0
    {
        return Err(ThermoError::NumericalFailure(
            "invalid mixed EOS parameters".into(),
        ));
    }
    let roots_z = physical_roots(reduced_a, reduced_b)?;
    let z = match root {
        RootSelection::Smallest => roots_z[0],
        RootSelection::Largest => roots_z[roots_z.len() - 1],
    };
    let v = z * R * t / p;
    let denominator = v * v + 2.0 * b * v - b * b;
    let dpdv = -R * t / (v - b).powi(2) + a * (2.0 * v + 2.0 * b) / denominator.powi(2);
    if !dpdv.is_finite() || dpdv >= 0.0 {
        return Err(ThermoError::Unsupported(
            "selected root is mechanically unstable or singular".into(),
        ));
    }
    let sqrt2 = 2.0_f64.sqrt();
    let log_ratio = (2.0 * sqrt2 * reduced_b / (z + (1.0 - sqrt2) * reduced_b)).ln_1p();
    let attraction = a / (2.0 * sqrt2 * b * R * t) * log_ratio;
    let ln_fugacity_coefficients: Vec<_> = (0..n)
        .map(|i| {
            let beta = bi[i] / b;
            beta * (z - 1.0) - (z - reduced_b).ln() - attraction * (2.0 * row_sums[i] / a - beta)
        })
        .collect();
    let h = R * t * (z - 1.0) + (t * da - a) / (2.0 * sqrt2 * b) * log_ratio;
    let s = R * (z - reduced_b).ln() + da / (2.0 * sqrt2 * b) * log_ratio;
    let rho = 1.0 / v;
    let calculated_pressure = R * t / (v - b) - a / denominator;
    let pressure_relative_residual = (calculated_pressure - p).abs() / p;
    let mass_density = rho * mixture.molar_mass_kg_mol();
    if ![z, v, rho, mass_density, h, s, pressure_relative_residual]
        .iter()
        .all(|v| v.is_finite())
        || !ln_fugacity_coefficients.iter().all(|v| v.is_finite())
        || pressure_relative_residual > 1e-8
    {
        return Err(ThermoError::NumericalFailure(
            "nonfinite result or EOS pressure residual above 1e-8".into(),
        ));
    }
    Ok(HomogeneousProperties {
        roots_z,
        z,
        molar_volume_m3_mol: v,
        molar_density_mol_m3: rho,
        mass_density_kg_m3: mass_density,
        ln_fugacity_coefficients,
        residual_enthalpy_j_mol: h,
        residual_entropy_j_mol_k: s,
        pressure_relative_residual,
    })
}

// Bracket the monic cubic between its turning points, then bisect each interval.
// This avoids cancellation in closed-form small-root expressions.
fn physical_roots(a: f64, b: f64) -> Result<Vec<f64>, ThermoError> {
    let c2 = b - 1.0;
    let c1 = a - 3.0 * b * b - 2.0 * b;
    let c0 = -a * b + b * b + b * b * b;
    let f = |z: f64| ((z + c2) * z + c1) * z + c0;
    let upper = 1.0 + c2.abs().max(c1.abs()).max(c0.abs());
    if !upper.is_finite() || b >= upper {
        return Err(ThermoError::NumericalFailure(
            "invalid cubic root bound".into(),
        ));
    }
    let mut bounds = vec![b];
    let discriminant = c2 * c2 - 3.0 * c1;
    // A triple root can have a zero or slightly negative discriminant after
    // rounding, so testing only distinct turning points misses this case.
    if discriminant.abs() < 1e-12 && f(-c2 / 3.0).abs() < 1e-14 {
        return Err(ThermoError::Unsupported(
            "near-critical triple root requires critical-region handling".into(),
        ));
    }
    if discriminant > 0.0 {
        for point in [
            (-c2 - discriminant.sqrt()) / 3.0,
            (-c2 + discriminant.sqrt()) / 3.0,
        ] {
            if point > b && point < upper {
                if f(point).abs() < 1e-14 {
                    return Err(ThermoError::Unsupported(
                        "near-multiple cubic root requires critical-region handling".into(),
                    ));
                }
                bounds.push(point);
            }
        }
    }
    bounds.push(upper);
    let mut roots = Vec::new();
    for interval in bounds.windows(2) {
        let (mut lo, mut hi) = (interval[0], interval[1]);
        let mut flo = f(lo);
        if flo.signum() == f(hi).signum() {
            continue;
        }
        for _ in 0..160 {
            let mid = lo + (hi - lo) / 2.0;
            if mid == lo || mid == hi {
                break;
            }
            let fm = f(mid);
            if fm == 0.0 {
                lo = mid;
                hi = mid;
                break;
            }
            if fm.signum() == flo.signum() {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }
        let root = (lo + hi) / 2.0;
        if root.is_finite() && root > b {
            roots.push(root);
        }
    }
    if roots.is_empty() {
        return Err(ThermoError::NumericalFailure(
            "no physical Z > B root found".into(),
        ));
    }
    Ok(roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_pr_critical_triple_root_is_explicitly_unsupported() {
        assert!(matches!(
            physical_roots(OMEGA_A, OMEGA_B),
            Err(ThermoError::Unsupported(_))
        ));
    }
}
