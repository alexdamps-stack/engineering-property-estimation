//! Petroleum pseudo-component characterisation (#229).
//!
//! A pseudo-component is a narrow boiling cut described by its normal boiling
//! point Tb and specific gravity SG (60 °F/60 °F), and optionally a measured
//! molar mass. Its cubic-EOS constants and ideal-gas heat capacity come from
//! documented correlations:
//!
//! | Quantity | Correlation |
//! |---|---|
//! | Tc, Pc | Kesler and Lee (1976), SI form |
//! | ω | Lee and Kesler (1975) for Tb/Tc < 0.8; Kesler and Lee (1976) above |
//! | M (unless measured) | Riazi and Daubert (1987) |
//! | ideal-gas cp | Kesler and Lee (1976), quadratic in T |
//!
//! Coefficients are transcribed from Fahim, Al-Sahhaf and Elkilani,
//! *Fundamentals of Petroleum Refining*, Elsevier (2010), chapter 3,
//! Eqs. (3.5), (3.11), (3.25)–(3.28) and (3.38)–(3.41). The published worked
//! examples used to check them are transcribed in
//! `validation/pseudo-components/fahim-2010-worked-examples.json`. See
//! `docs/pseudo-components-status.md` for the declared domain and evidence.
//!
//! Pseudo-components have no molecular formula and no formation properties:
//! reactors reject them, and [`crate::Component::standard_gibbs_energy`]
//! returns `MissingParameter`.
use std::sync::Arc;

use serde::{Deserialize, Serialize};
use thermo_types::{CalculationError, ComponentId, MOLAR_GAS_CONSTANT as R};

use crate::{Component, Constants, Nasa7, Record};

/// Prefix of every pseudo-component identifier. Embedded database
/// identifiers never start with it, so the two cannot collide.
pub const PSEUDO_COMPONENT_PREFIX: &str = "pseudo:";
/// Identifier of the characterisation method set reported in data provenance.
pub const PSEUDO_CHARACTERISATION_ID: &str =
    "petroleum-pseudo-kesler-lee-1976-lee-kesler-1975-riazi-daubert-1987";
/// Version of the characterisation formulation. A change to any coefficient,
/// domain or construction bumps it, so saved results become stale.
pub const PSEUDO_CHARACTERISATION_VERSION: &str = "2026-11-24.1";
/// Declared normal-boiling-point domain of the characterisation, K.
pub const PSEUDO_BOILING_POINT_RANGE_K: [f64; 2] = [300.0, 850.0];
/// Declared specific-gravity (60 °F/60 °F) domain of the characterisation.
pub const PSEUDO_SPECIFIC_GRAVITY_RANGE: [f64; 2] = [0.60, 1.00];
/// Declared molar-mass domain, kg/mol (70–700 g/mol).
pub const PSEUDO_MOLAR_MASS_RANGE_KG_PER_MOL: [f64; 2] = [0.070, 0.700];
/// Declared temperature domain of the Kesler–Lee ideal-gas heat capacity, K.
/// Ideal-gas properties of a pseudo-component outside it are
/// `OutsideValidityDomain`; packages must keep their search interval inside.
pub const PSEUDO_IDEAL_GAS_RANGE_K: [f64; 2] = [250.0, 900.0];
/// Pressure at which Tb is defined (1 atm), Pa.
pub const NORMAL_BOILING_PRESSURE_PA: f64 = 101_325.0;

/// Stored characterisation input of one pseudo-component. Projects persist
/// exactly this record; every derived constant is recomputed from it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PseudoComponentInput {
    /// Stable identifier starting with [`PSEUDO_COMPONENT_PREFIX`].
    pub id: String,
    /// Normal (1 atm) boiling point of the cut, K. For a narrow TBP cut this
    /// is its mid boiling point; for a whole fraction, its mean average
    /// boiling point.
    pub normal_boiling_point_k: f64,
    /// Specific gravity at 60 °F/60 °F.
    pub specific_gravity: f64,
    /// Measured molar mass, kg/mol. `None` selects Riazi–Daubert (1987).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub molar_mass_kg_per_mol: Option<f64>,
    /// Critical constants from a stated source (for example a defined
    /// surrogate compound) instead of the Kesler–Lee correlations.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specified_constants: Option<SpecifiedConstants>,
}

/// Critical constants supplied with their source instead of correlated.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpecifiedConstants {
    pub critical_temperature_k: f64,
    pub critical_pressure_pa: f64,
    pub acentric_factor: f64,
    /// Where the values come from. Must not be empty.
    pub source: String,
}

/// Method that produced a characterised quantity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConstantSource {
    /// A named published correlation.
    Correlation(&'static str),
    /// A value supplied with the input, with its stated source.
    Specified(String),
}

/// Constants derived from a [`PseudoComponentInput`], with their methods.
#[derive(Debug, Clone, PartialEq)]
pub struct PseudoCharacterisation {
    pub input: PseudoComponentInput,
    /// Watson characterisation factor K = (1.8 Tb)^(1/3) / SG (Tb in K).
    pub watson_k: f64,
    pub critical_temperature_k: f64,
    pub critical_pressure_pa: f64,
    pub acentric_factor: f64,
    pub molar_mass_kg_per_mol: f64,
    pub critical_constants_source: ConstantSource,
    pub molar_mass_source: ConstantSource,
    /// Molar ideal-gas heat capacity cp/R = c0 + c1 T + c2 T², T in K
    /// (Kesler–Lee with this component's K, ω and M).
    pub ideal_gas_cp_over_r: [f64; 3],
}

/// Watson (UOP) characterisation factor from Tb in K and SG; Fahim et al.
/// (2010) Eq. (3.5), with Tb converted to °R.
pub fn watson_k(normal_boiling_point_k: f64, specific_gravity: f64) -> f64 {
    (1.8 * normal_boiling_point_k).cbrt() / specific_gravity
}

/// Specific gravity at constant Watson factor; inverse of [`watson_k`].
pub fn specific_gravity_from_watson_k(normal_boiling_point_k: f64, watson_k: f64) -> f64 {
    (1.8 * normal_boiling_point_k).cbrt() / watson_k
}

/// Specific gravity 60 °F/60 °F from API gravity: 141.5 / (API + 131.5).
pub fn specific_gravity_from_api(api_gravity: f64) -> f64 {
    141.5 / (api_gravity + 131.5)
}

/// Kesler–Lee (1976) critical temperature, K, from Tb in K and SG; Fahim et
/// al. (2010) Eq. (3.25). The Tb coefficient is 0.4244, the exact conversion
/// of the original °R form; the 0.422 printed in Eq. (3.25) is a misprint
/// (see `docs/pseudo-components-status.md`).
pub fn kesler_lee_critical_temperature(tb: f64, sg: f64) -> f64 {
    189.8 + 450.6 * sg + (0.4244 + 0.1174 * sg) * tb + (14_410.0 - 100_688.0 * sg) / tb
}

/// Kesler–Lee (1976) critical pressure, Pa, from Tb in K and SG; Fahim et al.
/// (2010) Eq. (3.26), which gives ln(Pc / bar).
pub fn kesler_lee_critical_pressure(tb: f64, sg: f64) -> f64 {
    let ln_pc_bar =
        5.689 - 0.0566 / sg - 1e-3 * tb * (0.436_392 + 4.121_64 / sg + 0.213_426 / (sg * sg))
            + 1e-7 * tb * tb * (4.757_94 + 11.819 / sg + 1.530_15 / (sg * sg))
            - 1e-10 * tb.powi(3) * (2.450_55 + 9.901 / (sg * sg));
    ln_pc_bar.exp() * 1e5
}

/// Acentric factor from Tb, Tc (K), Pc (Pa) and the Watson factor; Fahim et
/// al. (2010) Eqs. (3.27)–(3.28): Lee and Kesler (1975) for Tb/Tc < 0.8 and
/// Kesler and Lee (1976) otherwise. Pb is 1 atm.
pub fn lee_kesler_acentric_factor(tb: f64, tc: f64, pc: f64, watson_k: f64) -> f64 {
    let tbr = tb / tc;
    if tbr < 0.8 {
        // ln(Pb/Pc): Pb = 1 atm is the pressure at which Tb is defined.
        let ln_pbr = (NORMAL_BOILING_PRESSURE_PA / pc).ln();
        (ln_pbr - 5.927_14 + 6.096_48 / tbr + 1.288_62 * tbr.ln() - 0.169_347 * tbr.powi(6))
            / (15.2518 - 15.6875 / tbr - 13.4721 * tbr.ln() + 0.435_77 * tbr.powi(6))
    } else {
        let k = watson_k;
        -7.904 + 0.1352 * k - 0.007_465 * k * k + 8.359 * tbr + (1.408 - 0.010_63 * k) / tbr
    }
}

/// Riazi–Daubert (1987) molar mass, kg/mol, from Tb in K and SG; Fahim et al.
/// (2010) Eq. (3.11) (there attributed to Pedersen et al. 1989).
pub fn riazi_daubert_molar_mass(tb: f64, sg: f64) -> f64 {
    42.965
        * (2.097e-4 * tb - 7.787_12 * sg + 2.084_76e-3 * tb * sg).exp()
        * tb.powf(1.260_07)
        * sg.powf(4.983_08)
        / 1000.0
}

/// Kesler–Lee (1976) ideal-gas specific heat capacity coefficients
/// `[a0, a1, a2]` with cp = a0 + a1 T + a2 T² in kJ/(kg K) and T in K:
///
/// ```text
/// a0 = −1.41779 + 0.11828 K − CF (1.09223 − 2.48245 ω)
/// a1 = −(6.99724 − 8.69326 K + 0.27715 K²)·1e-4 + CF (3.434 − 7.14 ω)·1e-3
/// a2 = −2.2582e-6 + CF (7.2661 − 9.2561 ω)·1e-7
/// CF = ((12.8 − K)(10 − K) / (10 ω))²
/// ```
///
/// Fahim et al. (2010) Eqs. (3.38)–(3.41) print the same correlation per
/// kmol as A = 4.1843 M (−0.33886 + 0.02827 K − 0.26105 CF + 0.59332 ω CF),
/// B = −7.5318 M [(0.9291 − 1.1543 K + 0.0368 K²) − CF (4.56 − 9.48 ω)]·1e-4
/// and C = 13.5573 M [−1.6658 + CF (0.536 − 0.6828 ω)]·1e-7, which is this
/// form with the conversion factors 4.1843 and 1.8 folded into the
/// coefficients (every coefficient agrees to 2e-4 relative or better).
pub fn kesler_lee_ideal_gas_cp(watson_k: f64, omega: f64) -> [f64; 3] {
    let k = watson_k;
    let cf = ((12.8 - k) * (10.0 - k) / (10.0 * omega)).powi(2);
    [
        -1.417_79 + 0.118_28 * k - cf * (1.092_23 - 2.482_45 * omega),
        -(6.997_24 - 8.693_26 * k + 0.277_15 * k * k) * 1e-4 + cf * (3.434 - 7.14 * omega) * 1e-3,
        -2.2582e-6 + cf * (7.2661 - 9.2561 * omega) * 1e-7,
    ]
}

fn outside(what: &str, value: f64, range: [f64; 2], id: &str) -> CalculationError {
    CalculationError::OutsideValidityDomain(format!(
        "pseudo-component {id}: {what} {value} outside the declared characterisation domain {}..={}",
        range[0], range[1]
    ))
}

fn check_range(what: &str, value: f64, range: [f64; 2], id: &str) -> Result<(), CalculationError> {
    if !value.is_finite() {
        return Err(thermo_types::InputError::NonFinite("pseudo-component input").into());
    }
    if !(range[0]..=range[1]).contains(&value) {
        return Err(outside(what, value, range, id));
    }
    Ok(())
}

/// Validates the identifier of a pseudo-component.
pub fn check_pseudo_component_id(id: &str) -> Result<(), CalculationError> {
    ComponentId::new(id)?;
    if !id.starts_with(PSEUDO_COMPONENT_PREFIX) || id.len() == PSEUDO_COMPONENT_PREFIX.len() {
        return Err(CalculationError::InvalidInput(
            thermo_types::InputError::InvalidIdentifier,
        ));
    }
    Ok(())
}

impl PseudoComponentInput {
    /// A characterised cut: Tb in K, SG 60 °F/60 °F, molar mass from
    /// Riazi–Daubert and critical constants from Kesler–Lee.
    pub fn new(id: impl Into<String>, normal_boiling_point_k: f64, specific_gravity: f64) -> Self {
        Self {
            id: id.into(),
            normal_boiling_point_k,
            specific_gravity,
            molar_mass_kg_per_mol: None,
            specified_constants: None,
        }
    }

    /// Derives the constants and ideal-gas heat capacity. Inputs outside the
    /// declared domain are `OutsideValidityDomain`; derived constants that are
    /// not physical (Tc ≤ Tb, Pc ≤ 1 atm, nonpositive cp in the ideal-gas
    /// domain) are `NumericalFailure`. Nothing is clipped or extrapolated.
    pub fn characterise(&self) -> Result<PseudoCharacterisation, CalculationError> {
        let id = self.id.as_str();
        check_pseudo_component_id(id)?;
        let tb = self.normal_boiling_point_k;
        let sg = self.specific_gravity;
        check_range("Tb/K", tb, PSEUDO_BOILING_POINT_RANGE_K, id)?;
        check_range("SG", sg, PSEUDO_SPECIFIC_GRAVITY_RANGE, id)?;
        let k = watson_k(tb, sg);
        let (molar_mass, molar_mass_source) = match self.molar_mass_kg_per_mol {
            Some(m) => (m, ConstantSource::Specified("measured molar mass".into())),
            None => (
                riazi_daubert_molar_mass(tb, sg),
                ConstantSource::Correlation("Riazi-Daubert 1987"),
            ),
        };
        check_range(
            "M/(kg/mol)",
            molar_mass,
            PSEUDO_MOLAR_MASS_RANGE_KG_PER_MOL,
            id,
        )?;
        let (tc, pc, omega, critical_source) = match &self.specified_constants {
            Some(s) => {
                if s.source.trim().is_empty() {
                    return Err(CalculationError::MissingParameter(format!(
                        "pseudo-component {id}: specified critical constants need a source"
                    )));
                }
                (
                    s.critical_temperature_k,
                    s.critical_pressure_pa,
                    s.acentric_factor,
                    ConstantSource::Specified(s.source.clone()),
                )
            }
            None => {
                let tc = kesler_lee_critical_temperature(tb, sg);
                let pc = kesler_lee_critical_pressure(tb, sg);
                (
                    tc,
                    pc,
                    lee_kesler_acentric_factor(tb, tc, pc, k),
                    ConstantSource::Correlation("Kesler-Lee 1976 Tc/Pc; Lee-Kesler 1975 omega"),
                )
            }
        };
        if ![tc, pc, omega, molar_mass].iter().all(|v| v.is_finite())
            || tc <= tb
            || pc <= NORMAL_BOILING_PRESSURE_PA
            || omega <= 0.0
        {
            return Err(CalculationError::NumericalFailure(format!(
                "pseudo-component {id}: unphysical constants Tc={tc} K, Pc={pc} Pa, omega={omega}"
            )));
        }
        let per_mass = kesler_lee_ideal_gas_cp(k, omega);
        // kJ/(kg K) × kg/mol × 1000 J/kJ = J/(mol K); divide by R.
        let scale = molar_mass * 1000.0 / R;
        let cp = per_mass.map(|a| a * scale);
        let [lo, hi] = PSEUDO_IDEAL_GAS_RANGE_K;
        let at = |t: f64| cp[0] + t * (cp[1] + t * cp[2]);
        let vertex = -cp[1] / (2.0 * cp[2]);
        let mut checks = vec![lo, hi];
        if cp[2] > 0.0 && (lo..=hi).contains(&vertex) {
            checks.push(vertex);
        }
        if cp.iter().any(|c| !c.is_finite()) || checks.iter().any(|t| at(*t) <= 0.0) {
            return Err(CalculationError::NumericalFailure(format!(
                "pseudo-component {id}: Kesler-Lee ideal-gas cp is not positive on {lo}..={hi} K"
            )));
        }
        Ok(PseudoCharacterisation {
            input: self.clone(),
            watson_k: k,
            critical_temperature_k: tc,
            critical_pressure_pa: pc,
            acentric_factor: omega,
            molar_mass_kg_per_mol: molar_mass,
            critical_constants_source: critical_source,
            molar_mass_source,
            ideal_gas_cp_over_r: cp,
        })
    }

    /// The characterised component, ready to extend a component database.
    pub fn component(&self) -> Result<Component, CalculationError> {
        let c = self.characterise()?;
        let [lo, hi] = PSEUDO_IDEAL_GAS_RANGE_K;
        let [c0, c1, c2] = c.ideal_gas_cp_over_r;
        // One quadratic cp/R on the whole domain: both NASA7-shaped intervals
        // carry the same coefficients and the join is the upper bound, so the
        // integration constants a6/a7 are unused (no formation basis).
        let coefficients = [c0, c1, c2, 0.0, 0.0, 0.0, 0.0];
        Ok(Component(Record {
            id: c.input.id.clone(),
            cas: String::new(),
            formula: String::new(),
            constants: Constants {
                tc: c.critical_temperature_k,
                pc: c.critical_pressure_pa,
                acentric_factor: c.acentric_factor,
                molar_mass_kg_per_mol: c.molar_mass_kg_per_mol,
            },
            ideal_gas: Nasa7 {
                breaks: [lo, hi, hi],
                coefficients: [coefficients; 2],
            },
            extension: false,
            pseudo: Some(Arc::new(c)),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use thermo_types::Temperature;

    fn kerosene_cut() -> PseudoComponentInput {
        PseudoComponentInput::new("pseudo:cut-1", 480.0, 0.80)
    }

    #[test]
    fn identifiers_must_carry_the_prefix() {
        for bad in ["methane", "pseudo:", "", "pseudo:a b"] {
            let mut input = kerosene_cut();
            input.id = bad.into();
            assert!(
                matches!(input.characterise(), Err(CalculationError::InvalidInput(_))),
                "{bad}"
            );
        }
    }

    #[test]
    fn inputs_outside_the_domain_are_rejected_not_clipped() {
        for (tb, sg) in [(290.0, 0.8), (900.0, 0.8), (480.0, 0.55), (480.0, 1.05)] {
            let input = PseudoComponentInput::new("pseudo:x", tb, sg);
            assert!(
                matches!(
                    input.characterise(),
                    Err(CalculationError::OutsideValidityDomain(_))
                ),
                "{tb} {sg}"
            );
        }
        let mut nan = kerosene_cut();
        nan.specific_gravity = f64::NAN;
        assert!(matches!(
            nan.characterise(),
            Err(CalculationError::InvalidInput(_))
        ));
        let mut heavy = kerosene_cut();
        heavy.molar_mass_kg_per_mol = Some(0.9);
        assert!(matches!(
            heavy.characterise(),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }

    #[test]
    fn specified_constants_need_a_source() {
        let mut input = kerosene_cut();
        input.specified_constants = Some(SpecifiedConstants {
            critical_temperature_k: 650.0,
            critical_pressure_pa: 2e6,
            acentric_factor: 0.5,
            source: " ".into(),
        });
        assert!(matches!(
            input.characterise(),
            Err(CalculationError::MissingParameter(_))
        ));
        input.specified_constants.as_mut().unwrap().source = "test".into();
        let c = input.characterise().unwrap();
        assert_eq!(c.critical_temperature_k, 650.0);
        assert_eq!(
            c.critical_constants_source,
            ConstantSource::Specified("test".into())
        );
    }

    /// The NASA7-shaped record integrates the quadratic cp exactly: h and s
    /// agree with closed-form integrals from 298.15 K, and cp is continuous
    /// up to the upper bound.
    #[test]
    fn ideal_gas_record_integrates_the_quadratic_exactly() {
        let input = kerosene_cut();
        let c = input.characterise().unwrap();
        let component = input.component().unwrap();
        let [c0, c1, c2] = c.ideal_gas_cp_over_r;
        let t0 = crate::REFERENCE_TEMPERATURE;
        for t in [250.0, 298.15, 450.0, 700.0, 900.0] {
            let ig = component.ideal_gas(Temperature::new(t).unwrap()).unwrap();
            let h = R
                * (c0 * (t - t0)
                    + c1 / 2.0 * (t * t - t0 * t0)
                    + c2 / 3.0 * (t.powi(3) - t0.powi(3)));
            let s = R * (c0 * (t / t0).ln() + c1 * (t - t0) + c2 / 2.0 * (t * t - t0 * t0));
            let cp = R * (c0 + c1 * t + c2 * t * t);
            assert!((ig.enthalpy - h).abs() <= 1e-9 * h.abs().max(1.0), "{t}");
            assert!((ig.entropy - s).abs() <= 1e-9 * s.abs().max(1.0), "{t}");
            assert!((ig.heat_capacity - cp).abs() <= 1e-12 * cp, "{t}");
        }
        for t in [249.0, 901.0] {
            assert!(matches!(
                component.ideal_gas(Temperature::new(t).unwrap()),
                Err(CalculationError::OutsideValidityDomain(_))
            ));
        }
    }

    #[test]
    fn pseudo_components_have_no_formation_basis() {
        let component = kerosene_cut().component().unwrap();
        assert!(component.is_pseudo());
        assert!(component.elements().is_empty());
        assert!(component.standard_formation_enthalpy().is_nan());
        assert!(matches!(
            component.standard_gibbs_energy(Temperature::new(400.0).unwrap()),
            Err(CalculationError::MissingParameter(_))
        ));
    }

    #[test]
    fn stored_input_round_trips_and_refuses_unknown_fields() {
        let mut input = kerosene_cut();
        input.molar_mass_kg_per_mol = Some(0.17);
        let text = serde_json::to_string(&input).unwrap();
        assert_eq!(
            serde_json::from_str::<PseudoComponentInput>(&text).unwrap(),
            input
        );
        let minimal = r#"{"id":"pseudo:a","normal_boiling_point_k":450.0,"specific_gravity":0.78}"#;
        assert!(serde_json::from_str::<PseudoComponentInput>(minimal).is_ok());
        let unknown =
            r#"{"id":"pseudo:a","normal_boiling_point_k":450.0,"specific_gravity":0.78,"tc":1}"#;
        assert!(serde_json::from_str::<PseudoComponentInput>(unknown).is_err());
    }
}
