//! Assay distillation data to pseudo-component inputs (#229).
//!
//! Distillation temperatures are in K and volume percentages in liquid
//! volume percent distilled. Coefficients are transcribed from Fahim,
//! Al-Sahhaf and Elkilani, *Fundamentals of Petroleum Refining* (2010),
//! chapter 3: the API (Riazi–Daubert 1980) ASTM D86 to TBP conversion,
//! Eq. (3.3) and Table 3.2, and the D86 volume and mean average boiling
//! points, Eqs. (3.6)–(3.7). Cuts follow the equal-temperature-interval
//! breakup of section 3.3.1 with constant Watson K (section 3.3.2,
//! Example E3.5). Nothing is extrapolated beyond the supplied curve.
use thermo_types::{CalculationError, InputError};

use crate::{PseudoComponentInput, specific_gravity_from_watson_k};

/// `(volume percent, a, b)` of the API D86→TBP method TBP = a·D86^b
/// (temperatures in °R); Fahim et al. (2010) Table 3.2.
pub const D86_TO_TBP_API_CONSTANTS: [(f64, f64, f64); 7] = [
    (0.0, 0.9167, 1.0019),
    (10.0, 0.5277, 1.0900),
    (30.0, 0.7429, 1.0425),
    (50.0, 0.8920, 1.0176),
    (70.0, 0.8705, 1.0226),
    (90.0, 0.9490, 1.0110),
    (95.0, 0.8008, 1.0355),
];

fn check_curve(points: &[(f64, f64)], what: &str) -> Result<(), CalculationError> {
    if points.len() < 2 {
        return Err(CalculationError::MissingParameter(format!(
            "{what}: at least two distillation points are required"
        )));
    }
    if points.iter().any(|(v, t)| !v.is_finite() || !t.is_finite()) {
        return Err(InputError::NonFinite("distillation point").into());
    }
    if points
        .iter()
        .any(|(v, t)| !(0.0..=100.0).contains(v) || *t <= 0.0)
    {
        return Err(CalculationError::OutsideValidityDomain(format!(
            "{what}: volume percent must lie in 0..=100 and temperatures must be positive"
        )));
    }
    if points
        .windows(2)
        .any(|w| w[1].0 <= w[0].0 || w[1].1 <= w[0].1)
    {
        return Err(CalculationError::UnsupportedCapability(format!(
            "{what}: volume percent and temperature must both increase strictly"
        )));
    }
    Ok(())
}

/// Converts ASTM D86 points `(volume %, T/K)` to TBP with the API
/// (Riazi–Daubert 1980) method. Every volume percent must be one of the
/// tabulated 0, 10, 30, 50, 70, 90 or 95; other points are refused.
pub fn d86_to_tbp(d86: &[(f64, f64)]) -> Result<Vec<(f64, f64)>, CalculationError> {
    check_curve(d86, "ASTM D86 curve")?;
    d86.iter()
        .map(|&(volume, t)| {
            let (_, a, b) = D86_TO_TBP_API_CONSTANTS
                .iter()
                .find(|(v, _, _)| *v == volume)
                .ok_or_else(|| {
                    CalculationError::UnsupportedCapability(format!(
                        "API D86 to TBP conversion is tabulated at 0, 10, 30, 50, 70, 90 and 95 vol%, not {volume}"
                    ))
                })?;
            Ok((volume, a * (1.8 * t).powf(*b) / 1.8))
        })
        .collect()
}

fn d86_at(d86: &[(f64, f64)], volume: f64) -> Result<f64, CalculationError> {
    d86.iter()
        .find(|(v, _)| *v == volume)
        .map(|(_, t)| *t)
        .ok_or_else(|| {
            CalculationError::MissingParameter(format!("ASTM D86 temperature at {volume} vol%"))
        })
}

/// Volume average boiling point, K, from the D86 10, 30, 50, 70 and 90 vol%
/// temperatures; Fahim et al. (2010) Eq. (3.6).
pub fn d86_volume_average_boiling_point(d86: &[(f64, f64)]) -> Result<f64, CalculationError> {
    check_curve(d86, "ASTM D86 curve")?;
    let mut sum = 0.0;
    for v in [10.0, 30.0, 50.0, 70.0, 90.0] {
        sum += d86_at(d86, v)?;
    }
    Ok(sum / 5.0)
}

/// Mean average boiling point, K, from D86 data: MeABP = VABP − D with
/// ln D = −0.94402 − 0.00865 (VABP − 32)^0.6667 + 2.99791 SL^0.333, VABP in °F
/// and SL = (T90 − T10)/80 in °F per vol%; Fahim et al. (2010) Eq. (3.7).
pub fn d86_mean_average_boiling_point(d86: &[(f64, f64)]) -> Result<f64, CalculationError> {
    let fahrenheit = |k: f64| k * 1.8 - 459.67;
    let vabp_f = fahrenheit(d86_volume_average_boiling_point(d86)?);
    let slope = (fahrenheit(d86_at(d86, 90.0)?) - fahrenheit(d86_at(d86, 10.0)?)) / 80.0;
    if vabp_f <= 32.0 || slope <= 0.0 {
        return Err(CalculationError::OutsideValidityDomain(
            "MeABP correction needs VABP above 32 °F and a rising D86 curve".into(),
        ));
    }
    let d =
        (-0.944_02 - 0.008_65 * (vabp_f - 32.0).powf(0.6667) + 2.997_91 * slope.powf(0.333)).exp();
    Ok((vabp_f - d + 459.67) / 1.8)
}

/// One narrow cut of a TBP curve: its pseudo-component input and its liquid
/// volume fraction of the whole curve range.
#[derive(Debug, Clone, PartialEq)]
pub struct AssayCut {
    pub input: PseudoComponentInput,
    pub liquid_volume_fraction: f64,
    /// Lower and upper cut temperatures, K.
    pub boiling_range_k: [f64; 2],
}

fn volume_at(tbp: &[(f64, f64)], t: f64) -> f64 {
    let i = tbp
        .windows(2)
        .position(|w| t <= w[1].1)
        .unwrap_or(tbp.len() - 2);
    let ((v0, t0), (v1, t1)) = (tbp[i], tbp[i + 1]);
    v0 + (v1 - v0) * (t - t0) / (t1 - t0)
}

/// Splits a TBP curve `(volume %, T/K)` into `count` cuts of equal
/// temperature width between its first and last points. Each cut's Tb is
/// its mid boiling point, its SG follows from the constant Watson factor
/// of the whole fraction, and its volume is read from the curve by linear
/// interpolation. Volume fractions are relative to the supplied curve range
/// (first to last point), which must therefore cover the material to be
/// represented; the curve is never extrapolated. Identifiers are
/// `{prefix}{n}` with n from 1, and `prefix` must start with `pseudo:`.
pub fn cut_tbp_curve(
    tbp: &[(f64, f64)],
    count: usize,
    watson_k: f64,
    prefix: &str,
) -> Result<Vec<AssayCut>, CalculationError> {
    check_curve(tbp, "TBP curve")?;
    if count == 0 {
        return Err(InputError::NotPositive("number of cuts").into());
    }
    if !watson_k.is_finite() || watson_k <= 0.0 {
        return Err(InputError::NotPositive("Watson K").into());
    }
    let (v_first, t_first) = tbp[0];
    let (v_last, t_last) = tbp[tbp.len() - 1];
    let width = (t_last - t_first) / count as f64;
    let span = v_last - v_first;
    (0..count)
        .map(|i| {
            let lo = t_first + width * i as f64;
            let hi = if i + 1 == count {
                t_last
            } else {
                t_first + width * (i + 1) as f64
            };
            let tb = 0.5 * (lo + hi);
            let input = PseudoComponentInput::new(
                format!("{prefix}{}", i + 1),
                tb,
                specific_gravity_from_watson_k(tb, watson_k),
            );
            crate::check_pseudo_component_id(&input.id)?;
            Ok(AssayCut {
                input,
                liquid_volume_fraction: (volume_at(tbp, hi) - volume_at(tbp, lo)) / span,
                boiling_range_k: [lo, hi],
            })
        })
        .collect()
}

/// Mole fractions of cuts from their liquid volume fractions:
/// nᵢ ∝ Vᵢ·SGᵢ/Mᵢ, with each Mᵢ from the cut's characterisation.
pub fn cut_mole_fractions(cuts: &[AssayCut]) -> Result<Vec<f64>, CalculationError> {
    let moles = cuts
        .iter()
        .map(|c| {
            let m = c.input.characterise()?.molar_mass_kg_per_mol;
            Ok(c.liquid_volume_fraction * c.input.specific_gravity / m)
        })
        .collect::<Result<Vec<f64>, CalculationError>>()?;
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(InputError::ZeroTotal.into());
    }
    Ok(moles.into_iter().map(|n| n / total).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tbp() -> Vec<(f64, f64)> {
        vec![(0.0, 400.0), (50.0, 480.0), (100.0, 560.0)]
    }

    #[test]
    fn cuts_partition_the_curve() {
        let cuts = cut_tbp_curve(&tbp(), 4, 11.8, "pseudo:kero-").unwrap();
        assert_eq!(cuts.len(), 4);
        assert_eq!(cuts[0].input.id, "pseudo:kero-1");
        assert!((cuts[0].input.normal_boiling_point_k - 420.0).abs() < 1e-12);
        let total: f64 = cuts.iter().map(|c| c.liquid_volume_fraction).sum();
        assert!((total - 1.0).abs() < 1e-14);
        for c in &cuts {
            assert!((c.liquid_volume_fraction - 0.25).abs() < 1e-14);
            assert!(
                (crate::watson_k(c.input.normal_boiling_point_k, c.input.specific_gravity) - 11.8)
                    .abs()
                    < 1e-12
            );
        }
        let x = cut_mole_fractions(&cuts).unwrap();
        assert!((x.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        // Equal volumes: the lighter cuts have more moles.
        assert!(x.windows(2).all(|w| w[0] > w[1]));
    }

    #[test]
    fn malformed_curves_are_refused() {
        assert!(cut_tbp_curve(&tbp(), 0, 11.8, "pseudo:a").is_err());
        assert!(cut_tbp_curve(&tbp(), 2, 11.8, "cut-").is_err());
        assert!(cut_tbp_curve(&[(0.0, 400.0)], 2, 11.8, "pseudo:a").is_err());
        let falling = [(0.0, 500.0), (50.0, 480.0)];
        assert!(matches!(
            cut_tbp_curve(&falling, 2, 11.8, "pseudo:a"),
            Err(CalculationError::UnsupportedCapability(_))
        ));
        assert!(matches!(
            d86_to_tbp(&[(0.0, 300.0), (20.0, 320.0)]),
            Err(CalculationError::UnsupportedCapability(_))
        ));
        assert!(matches!(
            d86_mean_average_boiling_point(&[(10.0, 300.0), (30.0, 320.0)]),
            Err(CalculationError::MissingParameter(_))
        ));
    }
}
