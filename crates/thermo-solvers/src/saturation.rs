//! Pure-component saturation pressure from fugacity equality of the outer cubic roots.
use thermo_models::{CubicModel, RootSelection};
use thermo_types::{CalculationError as Error, Pressure, Temperature};

/// Iteration limit for the bracketed Newton solve in ln p.
pub const SATURATION_MAX_ITERATIONS: usize = 200;
/// Declared upper reduced temperature; closer to Tc returns a near-critical error.
pub const SATURATION_MAX_REDUCED_TEMPERATURE: f64 = 0.95;

#[derive(Debug, Clone)]
pub struct PureSaturation {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub liquid_compressibility: f64,
    pub vapor_compressibility: f64,
    pub liquid_molar_volume: f64,
    pub vapor_molar_volume: f64,
    pub ln_fugacity_coefficient: f64,
    /// |ln φ^L − ln φ^V| at the returned pressure.
    pub fugacity_residual: f64,
    pub iterations: usize,
}

/// Canonicalize the single active fraction, retaining all inactive slots.
pub(crate) fn canonical_pure_model(model: &CubicModel) -> Result<Option<CubicModel>, Error> {
    let mut active = model
        .fractions()
        .iter()
        .enumerate()
        .filter(|(_, x)| **x > 0.0);
    if let (Some((index, _)), None) = (active.next(), active.next()) {
        let mut fractions = vec![0.0; model.fractions().len()];
        fractions[index] = 1.0;
        return model.with_fractions(&fractions).map(Some);
    }
    Ok(None)
}

/// Saturation pressure of a pure model (exactly one active component). The
/// solve is bracketed by the spinodal pressures, where ln φ^L − ln φ^V is strictly
/// decreasing in p (its ln p derivative is Z^L − Z^V < 0).
pub fn pure_saturation_pressure(
    model: &CubicModel,
    temperature: Temperature,
    max_iterations: usize,
) -> Result<PureSaturation, Error> {
    let Some(pure_model) = canonical_pure_model(model)? else {
        return Err(Error::UnsupportedCapability(
            "pure saturation requires exactly one active component".into(),
        ));
    };
    let model = &pure_model;
    let index = model.fractions().iter().position(|v| *v == 1.0).unwrap();
    let t = temperature.value();
    let tc = model.components()[index].critical_temperature();
    if t / tc > SATURATION_MAX_REDUCED_TEMPERATURE {
        return Err(Error::OutsideValidityDomain(format!(
            "pure saturation supported for T/Tc <= {SATURATION_MAX_REDUCED_TEMPERATURE}; near-critical states are not supported"
        )));
    }
    if !(200.0..=600.0).contains(&t) {
        return Err(Error::OutsideValidityDomain(
            "pure saturation follows the 200..=600 K homogeneous EOS envelope".into(),
        ));
    }
    let (v_liquid_spinodal, v_vapor_spinodal) = spinodal_volumes(model, temperature)?;
    let p_low = model
        .pressure_derivatives(temperature, v_liquid_spinodal)?
        .0;
    let p_high = model.pressure_derivatives(temperature, v_vapor_spinodal)?.0;
    let mut lo = (p_low.max(0.01) * (1.0 + 1e-9)).ln();
    let mut hi = (p_high * (1.0 - 1e-9)).ln();
    if hi <= lo {
        return Err(Error::OutsideValidityDomain(
            "vapor spinodal pressure below the 0.01 Pa envelope".into(),
        ));
    }
    let difference = |ln_p: f64| -> Result<(f64, f64, f64, f64), Error> {
        let p = Pressure::new(ln_p.exp())?;
        let l = model.ln_fugacity(temperature, p, RootSelection::SmallestVolume)?;
        let v = model.ln_fugacity(temperature, p, RootSelection::LargestVolume)?;
        Ok((
            l.ln_fugacity_coefficients[index] - v.ln_fugacity_coefficients[index],
            l.compressibility,
            v.compressibility,
            v.ln_fugacity_coefficients[index],
        ))
    };
    if difference(lo)?.0 <= 0.0 {
        return Err(Error::OutsideValidityDomain(
            "saturation pressure below the 0.01 Pa envelope".into(),
        ));
    }
    let mut ln_p = lo + (hi - lo) / 2.0;
    for iteration in 1..=max_iterations {
        let (f, zl, zv, ln_phi) = difference(ln_p)?;
        if f > 0.0 {
            lo = ln_p;
        } else {
            hi = ln_p;
        }
        let newton = ln_p - f / (zl - zv);
        let next = if newton > lo && newton < hi {
            newton
        } else {
            lo + (hi - lo) / 2.0
        };
        if f == 0.0 || (next - ln_p).abs() <= 1e-15 * ln_p.abs().max(1.0) {
            let p = ln_p.exp();
            return Ok(PureSaturation {
                temperature_k: t,
                pressure_pa: p,
                liquid_compressibility: zl,
                vapor_compressibility: zv,
                liquid_molar_volume: zl * thermo_types::MOLAR_GAS_CONSTANT * t / p,
                vapor_molar_volume: zv * thermo_types::MOLAR_GAS_CONSTANT * t / p,
                ln_fugacity_coefficient: ln_phi,
                fugacity_residual: f.abs(),
                iterations: iteration,
            });
        }
        ln_p = next;
    }
    Err(Error::NonConvergence(format!(
        "pure saturation exhausted {max_iterations} iterations"
    )))
}

/// Volumes where (dp/dv)_T = 0 below Tc: a log-spaced scan of v/b followed by bisection.
fn spinodal_volumes(model: &CubicModel, temperature: Temperature) -> Result<(f64, f64), Error> {
    let b = model.covolume();
    let slope = |v: f64| model.pressure_derivatives(temperature, v).map(|d| d.1);
    let grid: Vec<f64> = (0..=2000)
        .map(|i| b * (1.0 + 10f64.powf(-4.0 + 9.0 * i as f64 / 2000.0)))
        .collect();
    let mut changes = Vec::new();
    let mut previous = slope(grid[0])?;
    for pair in grid.windows(2) {
        let current = slope(pair[1])?;
        if previous.signum() != current.signum() {
            changes.push((pair[0], pair[1]));
        }
        previous = current;
    }
    let [first, second] = changes.as_slice() else {
        return Err(Error::NumericalFailure(format!(
            "expected two spinodal volumes, found {}",
            changes.len()
        )));
    };
    let bisect = |(mut lo, mut hi): (f64, f64)| -> Result<f64, Error> {
        let s_lo = slope(lo)?;
        for _ in 0..200 {
            let mid = lo + (hi - lo) / 2.0;
            if mid == lo || mid == hi {
                break;
            }
            if slope(mid)?.signum() == s_lo.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Ok(lo + (hi - lo) / 2.0)
    };
    Ok((bisect(*first)?, bisect(*second)?))
}
