//! Engineering liquid density (#228): COSTALD saturated-liquid volumes
//! (Hankinson and Thomson 1979) with the COSTALD-Tait compressed-liquid
//! correction and generalized Riedel pseudo saturation pressure (Thomson,
//! Brobst and Hankinson 1982; API Technical Data Book procedure 6A3.4).
//!
//! This is a separately labelled *engineering* property. It is not an
//! equation of state: it has no Helmholtz potential and never replaces the EOS
//! volume in fugacities, caloric properties or derivative identities
//! (`docs/calculation-contract.md`, *Liquid density and caloric coupling*).
//!
//! Domain and qualification (`docs/liquid-density-status.md`):
//! - mixture reduced temperature 0.25 ≤ T/T_cm ≤ 0.95 and p ≤ 69 MPa,
//!   otherwise `OutsideValidityDomain`; nothing is extrapolated;
//! - the components in [`COSTALD_QUALIFIED_COMPONENTS`] were compared with
//!   reference equations of state; water, methanol and hydrogen exceed the
//!   declared budgets and are refused above a total trace mole fraction of
//!   [`COSTALD_TRACE_MOLE_FRACTION`];
//! - mixtures are `WithinDeclaredDomain` when every component above trace is
//!   in [`COSTALD_MIXTURE_QUALIFIED_COMPONENTS`], otherwise `Unassessed`
//!   (CO₂ + hydrocarbon liquids deviated by 3–9 %).
use thermo_data::{ComponentDatabase, CostaldDatabase};
use thermo_types::{
    CalculationError, DataMetadata, InputError, MOLAR_GAS_CONSTANT as R, ModelMetadata,
    ValidityStatus,
};

/// Model identifier recorded with every engineering density.
pub const LIQUID_DENSITY_MODEL_ID: &str = "COSTALD-HT1979+Tait-TBH1982";
/// Implementation version; bump when the formulation or the policy changes.
pub const LIQUID_DENSITY_IMPLEMENTATION_VERSION: &str = "1";
/// Declared mixture reduced-temperature interval (inclusive).
pub const COSTALD_REDUCED_TEMPERATURE_RANGE: [f64; 2] = [0.25, 0.95];
/// Declared upper pressure, Pa.
pub const COSTALD_MAXIMUM_PRESSURE_PA: f64 = 69e6;
/// Components whose pure-fluid densities agree with reference equations of
/// state within the declared budgets (`validation/fixtures/liquid-density-v1.json`).
pub const COSTALD_QUALIFIED_COMPONENTS: [&str; 7] = [
    "methane",
    "ethane",
    "propane",
    "nitrogen",
    "carbon-dioxide",
    "oxygen",
    "carbon-monoxide",
];
/// Components whose mixtures (bubble-point liquids) agree within the budgets.
pub const COSTALD_MIXTURE_QUALIFIED_COMPONENTS: [&str; 6] = [
    "methane",
    "ethane",
    "propane",
    "nitrogen",
    "oxygen",
    "carbon-monoxide",
];
/// Total mole fraction of non-qualified components that is tolerated, and the
/// threshold above which a component counts for the mixture qualification.
pub const COSTALD_TRACE_MOLE_FRACTION: f64 = 1e-3;
/// Allowed deviation of the mole-fraction sum from one.
const FRACTION_SUM_TOLERANCE: f64 = 1e-9;

#[derive(Debug, Clone, PartialEq)]
struct Pure {
    id: String,
    critical_temperature: f64,
    characteristic_volume: f64,
    omega_srk: f64,
    molar_mass: f64,
}

/// COSTALD engineering liquid density for a fixed, ordered component list.
#[derive(Debug, Clone, PartialEq)]
pub struct CostaldLiquidDensity {
    components: Vec<Pure>,
    data: DataMetadata,
}

/// One engineering liquid-density evaluation. Reported beside, never
/// instead of, the EOS density.
#[derive(Debug, Clone, PartialEq)]
pub struct EngineeringLiquidDensity {
    /// Engineering liquid molar volume, m³/mol.
    pub molar_volume: f64,
    /// Engineering liquid mass density, kg/m³.
    pub mass_density: f64,
    /// COSTALD saturated-liquid molar volume at T, m³/mol.
    pub saturated_molar_volume: f64,
    /// Generalized Riedel pseudo saturation pressure of the mixture, Pa. It
    /// is a correlation input, not the package's bubble pressure.
    pub pseudo_saturation_pressure_pa: f64,
    /// Whether the Tait correction was applied (p above the pseudo
    /// saturation pressure); otherwise the saturated value is reported.
    pub pressure_corrected: bool,
    /// Mixture reduced temperature T/T_cm.
    pub reduced_temperature: f64,
    pub model: ModelMetadata,
    pub data: DataMetadata,
    /// `WithinDeclaredDomain` or `Unassessed`; states outside the domain fail.
    pub validity: ValidityStatus,
}

impl CostaldLiquidDensity {
    /// Parameters for `ids` in that order. Every component needs database
    /// constants (T_c, molar mass) and COSTALD parameters; there is no
    /// generalized fallback.
    pub fn from_databases(
        components: &ComponentDatabase,
        costald: &CostaldDatabase,
        ids: &[&str],
    ) -> Result<Self, CalculationError> {
        if ids.is_empty() {
            return Err(InputError::EmptyComposition.into());
        }
        let components = ids
            .iter()
            .map(|id| {
                let c = components.get(id)?;
                let p = costald.get(id)?;
                Ok(Pure {
                    id: (*id).to_string(),
                    critical_temperature: c.critical_temperature(),
                    characteristic_volume: p.characteristic_volume,
                    omega_srk: p.omega_srk,
                    molar_mass: c.molar_mass().value(),
                })
            })
            .collect::<Result<Vec<_>, CalculationError>>()?;
        Ok(Self {
            components,
            data: costald.metadata(),
        })
    }

    /// The embedded databases for `ids`.
    pub fn embedded(ids: &[&str]) -> Result<Self, CalculationError> {
        Self::from_databases(
            &ComponentDatabase::embedded()?,
            &CostaldDatabase::embedded()?,
            ids,
        )
    }

    pub fn component_ids(&self) -> Vec<&str> {
        self.components.iter().map(|c| c.id.as_str()).collect()
    }

    pub fn metadata() -> ModelMetadata {
        ModelMetadata {
            identifier: LIQUID_DENSITY_MODEL_ID.into(),
            implementation_version: LIQUID_DENSITY_IMPLEMENTATION_VERSION.into(),
            formulation: "COSTALD saturated volume with Hankinson-Thomson mixing rules; COSTALD-Tait compressed-liquid correction above the generalized Riedel pseudo saturation pressure (API TDB 6A3.4); engineering property only, no EOS coupling".into(),
        }
    }

    pub fn data(&self) -> &DataMetadata {
        &self.data
    }

    /// Engineering liquid density at T (K), p (Pa) and liquid mole fractions
    /// in model order.
    pub fn evaluate(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EngineeringLiquidDensity, CalculationError> {
        if fractions.len() != self.components.len() {
            return Err(InputError::ComponentMismatch.into());
        }
        if fractions.iter().any(|x| !x.is_finite()) {
            return Err(InputError::NonFinite("mole fraction").into());
        }
        if fractions.iter().any(|x| *x < 0.0) {
            return Err(InputError::Negative("mole fraction").into());
        }
        let sum: f64 = fractions.iter().sum();
        if (sum - 1.0).abs() > FRACTION_SUM_TOLERANCE {
            return Err(InputError::NotNormalized.into());
        }
        for (value, what) in [(temperature_k, "temperature"), (pressure_pa, "pressure")] {
            if !value.is_finite() {
                return Err(InputError::NonFinite(what).into());
            }
            if value <= 0.0 {
                return Err(InputError::NotPositive(what).into());
            }
        }
        let x: Vec<f64> = fractions.iter().map(|v| v / sum).collect();
        let pairs = || self.components.iter().zip(&x);

        let unqualified: f64 = pairs()
            .filter(|(c, _)| !COSTALD_QUALIFIED_COMPONENTS.contains(&c.id.as_str()))
            .map(|(_, xi)| xi)
            .sum();
        if unqualified > COSTALD_TRACE_MOLE_FRACTION {
            let names: Vec<&str> = pairs()
                .filter(|(c, xi)| {
                    **xi > 0.0 && !COSTALD_QUALIFIED_COMPONENTS.contains(&c.id.as_str())
                })
                .map(|(c, _)| c.id.as_str())
                .collect();
            return Err(CalculationError::OutsideValidityDomain(format!(
                "COSTALD engineering liquid density is not qualified for {} \
                 (mole fraction {unqualified:.3e} above the trace limit {COSTALD_TRACE_MOLE_FRACTION})",
                names.join(", ")
            )));
        }

        // Hankinson-Thomson mixing rules.
        let (mut s1, mut s13, mut s23, mut st, mut omega, mut molar_mass) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for (c, xi) in pairs() {
            let v = c.characteristic_volume;
            let cbrt = v.cbrt();
            s1 += xi * v;
            s13 += xi * cbrt;
            s23 += xi * cbrt * cbrt;
            st += xi * (v * c.critical_temperature).sqrt();
            omega += xi * c.omega_srk;
            molar_mass += xi * c.molar_mass;
        }
        let vm = 0.25 * (s1 + 3.0 * s23 * s13);
        let tcm = st * st / vm;
        let tr = temperature_k / tcm;
        let [tr_lo, tr_hi] = COSTALD_REDUCED_TEMPERATURE_RANGE;
        if !(tr_lo..=tr_hi).contains(&tr) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "COSTALD reduced temperature {tr:.4} outside {tr_lo}..={tr_hi} \
                 (T = {temperature_k} K, pseudo-critical {tcm:.3} K)"
            )));
        }
        if pressure_pa > COSTALD_MAXIMUM_PRESSURE_PA {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "COSTALD-Tait pressure {pressure_pa} Pa above {COSTALD_MAXIMUM_PRESSURE_PA} Pa"
            )));
        }

        let tau = 1.0 - tr;
        let t13 = tau.cbrt();
        let v0 = 1.0 + t13 * (-1.52816 + t13 * (1.43907 + t13 * (-0.81446 + 0.190454 * t13)));
        let vdelta =
            (-0.296123 + tr * (0.386914 + tr * (-0.0427258 - 0.0480645 * tr))) / (tr - 1.00001);
        let saturated = vm * v0 * (1.0 - omega * vdelta);

        // Generalized Riedel pseudo saturation pressure (log10 form).
        let pcm = (0.291 - 0.080 * omega) * R * tcm / vm;
        let log_tr = tr.log10();
        let alpha = 35.0 - 36.0 / tr - 96.736 * log_tr + tr.powi(6);
        let beta = log_tr + 0.03721754 * alpha;
        let log_pr = 5.8031817 * log_tr + 0.07608141 * alpha + omega * 4.86601 * beta;
        let psat = 10f64.powf(log_pr) * pcm;

        let corrected = pressure_pa > psat;
        let molar_volume = if corrected {
            let e = (4.79594 + omega * (0.250047 + 1.14188 * omega)).exp();
            let c = 0.0861488 + 0.0344483 * omega;
            let b =
                pcm * (-1.0 + t13 * (-9.070217 + t13 * (62.45326 + t13 * (-135.1102 + e * t13))));
            let denominator = b + psat;
            if denominator.is_nan() || denominator <= 0.0 {
                return Err(CalculationError::NumericalFailure(
                    "COSTALD-Tait parameter B + p_sat is not positive".into(),
                ));
            }
            saturated * (1.0 - c * ((pressure_pa - psat) / denominator).ln_1p())
        } else {
            saturated
        };
        if ![saturated, molar_volume, psat]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(CalculationError::NumericalFailure(
                "COSTALD liquid volume is not finite and positive".into(),
            ));
        }

        let major: Vec<&str> = pairs()
            .filter(|(_, xi)| **xi > COSTALD_TRACE_MOLE_FRACTION)
            .map(|(c, _)| c.id.as_str())
            .collect();
        let validity = if major.len() <= 1
            || major
                .iter()
                .all(|id| COSTALD_MIXTURE_QUALIFIED_COMPONENTS.contains(id))
        {
            ValidityStatus::WithinDeclaredDomain
        } else {
            ValidityStatus::Unassessed
        };
        Ok(EngineeringLiquidDensity {
            molar_volume,
            mass_density: molar_mass / molar_volume,
            saturated_molar_volume: saturated,
            pseudo_saturation_pressure_pa: psat,
            pressure_corrected: corrected,
            reduced_temperature: tr,
            model: Self::metadata(),
            data: self.data.clone(),
            validity,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_inputs_are_rejected() {
        let m = CostaldLiquidDensity::embedded(&["propane", "ethane"]).unwrap();
        assert!(matches!(
            m.evaluate(&[1.0], 250.0, 1e6),
            Err(CalculationError::InvalidInput(
                InputError::ComponentMismatch
            ))
        ));
        assert!(matches!(
            m.evaluate(&[0.5, 0.4], 250.0, 1e6),
            Err(CalculationError::InvalidInput(InputError::NotNormalized))
        ));
        assert!(matches!(
            m.evaluate(&[1.5, -0.5], 250.0, 1e6),
            Err(CalculationError::InvalidInput(InputError::Negative(_)))
        ));
        assert!(matches!(
            m.evaluate(&[0.5, 0.5], f64::NAN, 1e6),
            Err(CalculationError::InvalidInput(InputError::NonFinite(_)))
        ));
        assert!(matches!(
            m.evaluate(&[0.5, 0.5], 250.0, 0.0),
            Err(CalculationError::InvalidInput(InputError::NotPositive(_)))
        ));
        assert!(matches!(
            CostaldLiquidDensity::embedded(&["propane", "unobtainium"]),
            Err(CalculationError::MissingParameter(_))
        ));
        assert!(matches!(
            CostaldLiquidDensity::embedded(&[]),
            Err(CalculationError::InvalidInput(InputError::EmptyComposition))
        ));
    }

    #[test]
    fn trace_policy_and_mixture_validity() {
        let m = CostaldLiquidDensity::embedded(&["propane", "water", "carbon-dioxide"]).unwrap();
        // Trace water is tolerated; the propane-rich liquid stays qualified.
        let r = m.evaluate(&[0.9995, 0.0005, 0.0], 300.0, 2e6).unwrap();
        assert_eq!(r.validity, ValidityStatus::WithinDeclaredDomain);
        // Water above the trace limit is refused, not extrapolated.
        assert!(matches!(
            m.evaluate(&[0.99, 0.01, 0.0], 300.0, 2e6),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        // CO2 + propane is computed but unassessed.
        let r = m.evaluate(&[0.8, 0.0, 0.2], 260.0, 2e6).unwrap();
        assert_eq!(r.validity, ValidityStatus::Unassessed);
        // Pure CO2 is qualified.
        let r = m.evaluate(&[0.0, 0.0, 1.0], 260.0, 5e6).unwrap();
        assert_eq!(r.validity, ValidityStatus::WithinDeclaredDomain);
        assert_eq!(r.model.identifier, LIQUID_DENSITY_MODEL_ID);
    }

    #[test]
    fn pressure_correction_is_monotonic_and_vanishes_at_saturation() {
        let m = CostaldLiquidDensity::embedded(&["propane"]).unwrap();
        let low = m.evaluate(&[1.0], 300.0, 1e5).unwrap();
        assert!(!low.pressure_corrected);
        assert_eq!(low.molar_volume, low.saturated_molar_volume);
        let mut last = low.molar_volume;
        for p in [2e6, 1e7, 3e7, 6e7] {
            let r = m.evaluate(&[1.0], 300.0, p).unwrap();
            assert!(r.pressure_corrected);
            assert!(r.molar_volume < last, "{p}");
            last = r.molar_volume;
        }
        let at = m
            .evaluate(
                &[1.0],
                300.0,
                low.pseudo_saturation_pressure_pa * (1.0 + 1e-12),
            )
            .unwrap();
        assert!((at.molar_volume / at.saturated_molar_volume - 1.0).abs() < 1e-12);
    }

    #[test]
    fn component_permutation_is_invariant() {
        let a = CostaldLiquidDensity::embedded(&["methane", "ethane", "propane"]).unwrap();
        let b = CostaldLiquidDensity::embedded(&["propane", "methane", "ethane"]).unwrap();
        let ra = a.evaluate(&[0.2, 0.3, 0.5], 200.0, 5e6).unwrap();
        let rb = b.evaluate(&[0.5, 0.2, 0.3], 200.0, 5e6).unwrap();
        assert!((ra.molar_volume / rb.molar_volume - 1.0).abs() < 1e-14);
        assert!((ra.mass_density / rb.mass_density - 1.0).abs() < 1e-14);
    }
}
