//! Product distillation curves from a composition (#190).
//!
//! A product made of components with known normal boiling points and
//! standard liquid volumes (petroleum pseudo-components, #229) has a TBP
//! curve, cumulative liquid volume percent distilled against temperature,
//! built here by the inverse of the equal-temperature cutting of
//! [`crate::cut_tbp_curve`]:
//!
//! - the components are ordered by normal boiling point Tb₁ ≤ … ≤ Tbₙ;
//! - component i boils over the interval between the midpoints to its
//!   neighbours, [(Tbᵢ₋₁ + Tbᵢ)/2, (Tbᵢ + Tbᵢ₊₁)/2]; the end intervals reach
//!   half the neighbouring spacing beyond the end components,
//!   Tb₁ − (Tb₂ − Tb₁)/2 and Tbₙ + (Tbₙ − Tbₙ₋₁)/2 (a single component boils
//!   at its Tb);
//! - within its interval the temperature rises linearly with the cumulative
//!   liquid volume, and each component's share of the curve is its share of
//!   the product's standard liquid volume.
//!
//! The intervals are fixed by the whole component list, including components
//! with zero amount, so a curve changes continuously with the composition.
//! The initial (0 %) and final (100 %) points are the lower end of the
//! lightest and the upper end of the heaviest component present. For a
//! product composed of equal-temperature cuts of a TBP curve in their assay
//! proportions the construction returns that curve, linearly interpolated,
//! at every cut boundary.
//!
//! ASTM D86 temperatures follow from the TBP curve by the exact inverse of
//! the API (Riazi–Daubert) correlation TBP = a·D86^b already used for assay
//! input ([`crate::D86_TO_TBP_API_CONSTANTS`]; M. R. Riazi and T. E. Daubert,
//! Hydrocarbon Processing 59(3) 115–116, 1980; Riazi and Daubert, Ind. Eng.
//! Chem. Process Des. Dev. 25, 1009–1015, 1986; constants as tabulated in
//! Fahim, Al-Sahhaf and Elkilani, *Fundamentals of Petroleum Refining*, 2010,
//! Eq. (3.3) and Table 3.2), D86 = (TBP/a)^(1/b) with temperatures in °R,
//! applied point by point at the tabulated 0, 10, 30, 50, 70, 90 and 95 vol%.
//! D86 temperatures between tabulated points are interpolated linearly in
//! volume percent; above 95 vol% the method has no constants and the D86
//! curve is not defined. The correlation was developed for whole petroleum
//! fractions; for a narrow cut its D86 points need not increase with volume,
//! and they are reported as calculated.
use thermo_types::{CalculationError, InputError};

use crate::D86_TO_TBP_API_CONSTANTS;

/// Volume percents reported for product curves.
pub const REPORTED_VOLUME_PERCENTS: [f64; 7] = [5.0, 10.0, 30.0, 50.0, 70.0, 90.0, 95.0];

/// Highest volume percent of the ASTM D86 curve.
pub const D86_MAX_VOLUME_PERCENT: f64 = 95.0;

/// One component of a product for its TBP curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveComponent {
    /// Normal boiling point, K.
    pub normal_boiling_point_k: f64,
    /// Standard liquid volume of the component in the product, in any
    /// consistent unit (for example m³ per mole of product); zero for a
    /// component that is absent but part of the component list.
    pub liquid_volume: f64,
}

/// One component's linear part of a TBP curve.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Segment {
    volume_percent: [f64; 2],
    temperature_k: [f64; 2],
}

/// TBP curve of a product, cumulative liquid volume percent against K.
#[derive(Debug, Clone, PartialEq)]
pub struct TbpCurve {
    /// Components present, in boiling order, with positive volume.
    segments: Vec<Segment>,
}

impl TbpCurve {
    /// Builds the curve of a product from its components (module
    /// documentation). Components may be given in any order.
    pub fn from_components(components: &[CurveComponent]) -> Result<Self, CalculationError> {
        if components.is_empty() {
            return Err(CalculationError::MissingParameter(
                "TBP curve: at least one component is required".into(),
            ));
        }
        for c in components {
            if !c.normal_boiling_point_k.is_finite() || !c.liquid_volume.is_finite() {
                return Err(InputError::NonFinite("TBP curve component").into());
            }
            if c.normal_boiling_point_k <= 0.0 {
                return Err(InputError::NotPositive("normal boiling point").into());
            }
            if c.liquid_volume < 0.0 {
                return Err(InputError::Negative("component liquid volume").into());
            }
        }
        let mut sorted = components.to_vec();
        sorted.sort_by(|a, b| {
            a.normal_boiling_point_k
                .total_cmp(&b.normal_boiling_point_k)
        });
        let total: f64 = sorted.iter().map(|c| c.liquid_volume).sum();
        if total <= 0.0 || !total.is_finite() {
            return Err(InputError::ZeroTotal.into());
        }
        let tb: Vec<f64> = sorted.iter().map(|c| c.normal_boiling_point_k).collect();
        let n = tb.len();
        // Interval boundaries: n + 1 temperatures.
        let mut bounds = Vec::with_capacity(n + 1);
        if n == 1 {
            bounds.extend([tb[0], tb[0]]);
        } else {
            bounds.push(tb[0] - 0.5 * (tb[1] - tb[0]));
            bounds.extend(tb.windows(2).map(|w| 0.5 * (w[0] + w[1])));
            bounds.push(tb[n - 1] + 0.5 * (tb[n - 1] - tb[n - 2]));
        }
        if bounds[0] <= 0.0 {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "TBP curve: the lightest component's interval reaches {} K; the boiling points are too widely spaced",
                bounds[0]
            )));
        }
        let mut cumulative = 0.0;
        let mut segments = Vec::new();
        for (i, c) in sorted.iter().enumerate() {
            if c.liquid_volume <= 0.0 {
                continue;
            }
            let lo = 100.0 * cumulative / total;
            cumulative += c.liquid_volume;
            segments.push(Segment {
                volume_percent: [lo, 100.0 * cumulative / total],
                temperature_k: [bounds[i], bounds[i + 1]],
            });
        }
        let last = segments.last_mut().expect("positive total volume");
        last.volume_percent[1] = 100.0;
        Ok(Self { segments })
    }

    /// TBP temperature, K, at `volume_percent` (0 to 100) distilled.
    pub fn temperature_at(&self, volume_percent: f64) -> Result<f64, CalculationError> {
        if !volume_percent.is_finite() {
            return Err(InputError::NonFinite("volume percent").into());
        }
        if !(0.0..=100.0).contains(&volume_percent) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "TBP curve is defined from 0 to 100 vol%, not {volume_percent}"
            )));
        }
        let segment = self
            .segments
            .iter()
            .find(|s| volume_percent <= s.volume_percent[1])
            .unwrap_or_else(|| self.segments.last().expect("non-empty curve"));
        let [v0, v1] = segment.volume_percent;
        let [t0, t1] = segment.temperature_k;
        let share = ((volume_percent - v0) / (v1 - v0)).clamp(0.0, 1.0);
        Ok(t0 + share * (t1 - t0))
    }

    /// Curve knots `(volume %, K)`: the initial point and the upper end of
    /// each component present.
    pub fn points(&self) -> Vec<(f64, f64)> {
        let first = &self.segments[0];
        std::iter::once((first.volume_percent[0], first.temperature_k[0]))
            .chain(
                self.segments
                    .iter()
                    .map(|s| (s.volume_percent[1], s.temperature_k[1])),
            )
            .collect()
    }

    /// ASTM D86 points at the tabulated volume percents of the API method.
    pub fn d86_points(&self) -> Result<Vec<(f64, f64)>, CalculationError> {
        let tbp = D86_TO_TBP_API_CONSTANTS
            .iter()
            .map(|(v, _, _)| Ok((*v, self.temperature_at(*v)?)))
            .collect::<Result<Vec<_>, CalculationError>>()?;
        tbp_to_d86(&tbp)
    }

    /// ASTM D86 temperature, K, at `volume_percent` (0 to 95), interpolated
    /// linearly between the tabulated points.
    pub fn d86_temperature_at(&self, volume_percent: f64) -> Result<f64, CalculationError> {
        if !volume_percent.is_finite() {
            return Err(InputError::NonFinite("volume percent").into());
        }
        if !(0.0..=D86_MAX_VOLUME_PERCENT).contains(&volume_percent) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "ASTM D86 curve from TBP is defined from 0 to 95 vol%, not {volume_percent}"
            )));
        }
        let d86 = self.d86_points()?;
        let i = d86
            .windows(2)
            .position(|w| volume_percent <= w[1].0)
            .expect("checked volume range");
        let ((v0, t0), (v1, t1)) = (d86[i], d86[i + 1]);
        Ok(t0 + (t1 - t0) * (volume_percent - v0) / (v1 - v0))
    }
}

/// Converts TBP points `(volume %, T/K)` to ASTM D86 by inverting the API
/// (Riazi–Daubert) method, D86 = (TBP/a)^(1/b) in °R. Every volume percent
/// must be one of the tabulated 0, 10, 30, 50, 70, 90 or 95. The points are
/// converted independently and need not be monotonic.
pub fn tbp_to_d86(tbp: &[(f64, f64)]) -> Result<Vec<(f64, f64)>, CalculationError> {
    tbp.iter()
        .map(|&(volume, t)| {
            if !volume.is_finite() || !t.is_finite() {
                return Err(InputError::NonFinite("distillation point").into());
            }
            if t <= 0.0 {
                return Err(InputError::NotPositive("TBP temperature").into());
            }
            let (_, a, b) = D86_TO_TBP_API_CONSTANTS
                .iter()
                .find(|(v, _, _)| *v == volume)
                .ok_or_else(|| {
                    CalculationError::UnsupportedCapability(format!(
                        "API TBP to D86 conversion is tabulated at 0, 10, 30, 50, 70, 90 and 95 vol%, not {volume}"
                    ))
                })?;
            Ok((volume, (1.8 * t / a).powf(1.0 / b) / 1.8))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{cut_tbp_curve, d86_to_tbp};

    fn component(tb: f64, volume: f64) -> CurveComponent {
        CurveComponent {
            normal_boiling_point_k: tb,
            liquid_volume: volume,
        }
    }

    #[test]
    fn equal_cuts_in_assay_proportions_return_the_assay_curve() {
        let assay = [(0.0, 340.0), (20.0, 400.0), (50.0, 460.0), (100.0, 550.0)];
        let cuts = cut_tbp_curve(&assay, 6, 11.8, "pseudo:c").unwrap();
        // Any common volume unit: the fractions themselves, in reverse order.
        let components: Vec<_> = cuts
            .iter()
            .rev()
            .map(|c| component(c.input.normal_boiling_point_k, c.liquid_volume_fraction))
            .collect();
        let curve = TbpCurve::from_components(&components).unwrap();
        let interpolate = |t: f64| {
            let i = assay.windows(2).position(|w| t <= w[1].1).unwrap();
            let ((v0, t0), (v1, t1)) = (assay[i], assay[i + 1]);
            v0 + (v1 - v0) * (t - t0) / (t1 - t0)
        };
        for (v, t) in curve.points() {
            assert!((interpolate(t) - v).abs() < 1e-12, "{v} {t}");
            assert!((curve.temperature_at(v).unwrap() - t).abs() < 1e-9);
        }
        assert_eq!(curve.points().len(), 7);
        assert!((curve.temperature_at(0.0).unwrap() - 340.0).abs() < 1e-12);
        assert!((curve.temperature_at(100.0).unwrap() - 550.0).abs() < 1e-12);
    }

    #[test]
    fn absent_components_keep_their_intervals() {
        // The middle component is absent: the curve jumps over its interval.
        let curve = TbpCurve::from_components(&[
            component(400.0, 1.0),
            component(420.0, 0.0),
            component(440.0, 1.0),
        ])
        .unwrap();
        assert_eq!(
            curve.points(),
            [(0.0, 390.0), (50.0, 410.0), (100.0, 450.0)]
        );
        assert!((curve.temperature_at(25.0).unwrap() - 400.0).abs() < 1e-12);
        assert!((curve.temperature_at(75.0).unwrap() - 440.0).abs() < 1e-12);
        // A trace of the middle component changes the curve continuously.
        let traced = TbpCurve::from_components(&[
            component(400.0, 1.0),
            component(420.0, 1e-12),
            component(440.0, 1.0),
        ])
        .unwrap();
        for v in [5.0, 49.0, 51.0, 95.0] {
            let gap = traced.temperature_at(v).unwrap() - curve.temperature_at(v).unwrap();
            assert!(gap.abs() < 1e-6, "{v}: {gap}");
        }
        // An absent lightest component moves the initial point up.
        let heavy =
            TbpCurve::from_components(&[component(400.0, 0.0), component(440.0, 2.0)]).unwrap();
        assert_eq!(heavy.points(), [(0.0, 420.0), (100.0, 460.0)]);
        let single = TbpCurve::from_components(&[component(400.0, 3.0)]).unwrap();
        assert_eq!(single.temperature_at(37.0).unwrap(), 400.0);
    }

    #[test]
    fn d86_inverts_the_api_tbp_conversion() {
        let d86 = [
            (0.0, 309.65),
            (10.0, 327.15),
            (30.0, 350.15),
            (50.0, 374.65),
            (70.0, 404.15),
            (90.0, 444.15),
            (95.0, 459.65),
        ];
        let back = tbp_to_d86(&d86_to_tbp(&d86).unwrap()).unwrap();
        for ((v, t), (vb, tb)) in d86.iter().zip(&back) {
            assert_eq!(v, vb);
            assert!((t - tb).abs() < 1e-10, "{v}: {t} {tb}");
        }
        let curve = TbpCurve::from_components(&[
            component(330.0, 1.0),
            component(370.0, 2.0),
            component(410.0, 2.0),
            component(450.0, 1.0),
        ])
        .unwrap();
        let points = curve.d86_points().unwrap();
        assert_eq!(points.len(), 7);
        let at_50 = curve.d86_temperature_at(50.0).unwrap();
        assert!((at_50 - points[3].1).abs() < 1e-12);
        let at_5 = curve.d86_temperature_at(5.0).unwrap();
        assert!((at_5 - 0.5 * (points[0].1 + points[1].1)).abs() < 1e-12);
    }

    #[test]
    fn malformed_inputs_are_refused() {
        assert!(TbpCurve::from_components(&[]).is_err());
        assert!(TbpCurve::from_components(&[component(400.0, 0.0)]).is_err());
        assert!(TbpCurve::from_components(&[component(400.0, -1.0)]).is_err());
        assert!(TbpCurve::from_components(&[component(f64::NAN, 1.0)]).is_err());
        assert!(matches!(
            TbpCurve::from_components(&[component(100.0, 1.0), component(400.0, 1.0)]),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        let curve = TbpCurve::from_components(&[component(400.0, 1.0)]).unwrap();
        assert!(curve.temperature_at(101.0).is_err());
        assert!(matches!(
            curve.d86_temperature_at(96.0),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(matches!(
            tbp_to_d86(&[(5.0, 400.0)]),
            Err(CalculationError::UnsupportedCapability(_))
        ));
    }
}
