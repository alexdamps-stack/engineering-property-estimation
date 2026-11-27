//! ChemSep v8.3 saturated-liquid density correlations (#228) for the
//! components that COSTALD does not qualify (water, methanol). Source
//! transcription: `data/source-extracts/liquid-density-chemsep.provisional.json`
//! (Artistic License 2.0, via the pinned DWSIM commit of #225). These feed
//! only the engineering liquid density; they never enter an EOS.
use crate::{ComponentDatabase, liquid::check_cas};
use serde::Deserialize;
use thermo_types::{CalculationError, ComponentId, DataMetadata};

/// Version of the embedded ChemSep saturated-liquid density extract.
pub const CHEMSEP_LIQUID_DENSITY_VERSION: &str = "2026-11-27.1";
/// Identifier of the embedded ChemSep saturated-liquid density dataset.
pub const CHEMSEP_LIQUID_DENSITY_ID: &str = "chemsep-8.3-liquid-density";
const SCHEMA: &str = "liquid-density-chemsep-source-extract-v1";
const EMBEDDED_EXTRACT: &str =
    include_str!("../../../data/source-extracts/liquid-density-chemsep.provisional.json");

#[derive(Debug, Clone, Deserialize)]
struct Record {
    id: String,
    cas: String,
    equation: u32,
    coefficients: [f64; 5],
    #[serde(rename = "temperature_range_K")]
    range: [f64; 2],
    #[serde(rename = "chemsep_critical_temperature_K")]
    critical_temperature: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct Extract {
    schema: String,
    version: String,
    components: Vec<Record>,
}

/// ChemSep saturated-liquid density correlation of one component (DIPPR
/// equation 105 or 106).
#[derive(Debug, Clone, PartialEq)]
pub struct SaturatedLiquidDensityCorrelation {
    pub equation: u32,
    pub coefficients: [f64; 5],
    /// Declared source range, K (inclusive).
    pub temperature_range_k: [f64; 2],
    /// ChemSep critical temperature used by equation 106, K.
    pub critical_temperature_k: f64,
}

impl SaturatedLiquidDensityCorrelation {
    /// Saturated-liquid molar density, mol/m³. Outside the source range:
    /// `OutsideValidityDomain`.
    pub fn molar_density(&self, id: &str, temperature_k: f64) -> Result<f64, CalculationError> {
        let [lo, hi] = self.temperature_range_k;
        if !(lo..=hi).contains(&temperature_k) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "ChemSep liquid density of {id} is defined for {lo}..={hi} K, not {temperature_k} K"
            )));
        }
        let [a, b, c, d, e] = self.coefficients;
        let kmol_m3 = match self.equation {
            105 => a / b.powf(1.0 + (1.0 - temperature_k / c).powf(d)),
            106 => {
                let tr = temperature_k / self.critical_temperature_k;
                a * (1.0 - tr).powf(b + tr * (c + tr * (d + e * tr)))
            }
            other => {
                return Err(CalculationError::UnsupportedCapability(format!(
                    "ChemSep equation {other}"
                )));
            }
        };
        let rho = 1000.0 * kmol_m3;
        if rho.is_finite() && rho > 0.0 {
            Ok(rho)
        } else {
            Err(CalculationError::NumericalFailure(format!(
                "{id} ChemSep liquid density is not positive at {temperature_k} K"
            )))
        }
    }
}

/// Checked ChemSep saturated-liquid density records. Components absent from
/// the extract have no component-specific correlation.
#[derive(Debug, Clone)]
pub struct ChemsepLiquidDensityDatabase {
    records: Vec<Record>,
}

impl ChemsepLiquidDensityDatabase {
    /// The embedded dataset, parsed and checked once per process.
    pub fn embedded() -> Result<Self, CalculationError> {
        static EMBEDDED: std::sync::OnceLock<
            Result<ChemsepLiquidDensityDatabase, CalculationError>,
        > = std::sync::OnceLock::new();
        EMBEDDED
            .get_or_init(|| Self::from_extract(EMBEDDED_EXTRACT, &ComponentDatabase::embedded()?))
            .clone()
    }

    fn from_extract(text: &str, components: &ComponentDatabase) -> Result<Self, CalculationError> {
        let extract: Extract = serde_json::from_str(text)
            .map_err(|e| CalculationError::MissingParameter(e.to_string()))?;
        let bad = |what: String| {
            CalculationError::MissingParameter(format!(
                "invalid embedded ChemSep liquid density {what}"
            ))
        };
        if extract.schema != SCHEMA || extract.version != CHEMSEP_LIQUID_DENSITY_VERSION {
            return Err(bad("extract header".into()));
        }
        for (k, r) in extract.components.iter().enumerate() {
            ComponentId::new(&r.id)?;
            if extract.components[..k].iter().any(|o| o.id == r.id) {
                return Err(bad(format!("duplicate record {}", r.id)));
            }
            check_cas(&r.id, &r.cas)?;
            let component = components
                .get(&r.id)
                .map_err(|_| bad(format!("record {}: not in the component database", r.id)))?;
            if component.cas() != r.cas {
                return Err(bad(format!(
                    "record {}: CAS differs from component database",
                    r.id
                )));
            }
            if !r.coefficients.iter().all(|v| v.is_finite())
                || !matches!(r.equation, 105 | 106)
                || !(0.0 < r.range[0] && r.range[0] < r.range[1])
                || !(r.critical_temperature.is_finite() && r.critical_temperature > 0.0)
            {
                return Err(bad(format!("record {}: correlation", r.id)));
            }
        }
        Ok(Self {
            records: extract.components,
        })
    }

    /// The correlation of `id`, if the extract has one.
    pub fn get(&self, id: &str) -> Option<SaturatedLiquidDensityCorrelation> {
        self.records
            .iter()
            .find(|r| r.id == id)
            .map(|r| SaturatedLiquidDensityCorrelation {
                equation: r.equation,
                coefficients: r.coefficients,
                temperature_range_k: r.range,
                critical_temperature_k: r.critical_temperature,
            })
    }

    pub fn metadata(&self) -> DataMetadata {
        DataMetadata {
            identifier: CHEMSEP_LIQUID_DENSITY_ID.into(),
            version: CHEMSEP_LIQUID_DENSITY_VERSION.into(),
            source: "data/source-extracts/liquid-density-chemsep.provisional.json; ChemSep v8.3 (Artistic License 2.0) via pinned DWSIM, URL, hash and notice retained".into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_and_methanol_correlations_load_and_evaluate() {
        let db = ChemsepLiquidDensityDatabase::embedded().unwrap();
        let water = db.get("water").unwrap();
        assert_eq!(water.equation, 106);
        // About 55.3 kmol/m³ (997 kg/m³) at 298.15 K.
        let rho = water.molar_density("water", 298.15).unwrap();
        assert!((rho / 55_350.0 - 1.0).abs() < 0.005, "{rho}");
        let methanol = db.get("methanol").unwrap();
        assert_eq!(methanol.equation, 105);
        assert!(matches!(
            methanol.molar_density("methanol", 600.0),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(db.get("propane").is_none());
        let edited = EMBEDDED_EXTRACT.replacen("\"equation\": 106", "\"equation\": 107", 1);
        assert_ne!(edited, EMBEDDED_EXTRACT);
        assert!(
            ChemsepLiquidDensityDatabase::from_extract(
                &edited,
                &ComponentDatabase::embedded().unwrap()
            )
            .is_err()
        );
    }
}
