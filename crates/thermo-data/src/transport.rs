//! Transport-property source data (#225): Chung gas-method constants and
//! ChemSep saturated-liquid correlations. Source transcription:
//! `data/source-extracts/transport.provisional.json`, produced only by
//! `validation/prepare-transport-source-data.py`.
use serde::Deserialize;
use thermo_types::{CalculationError, ComponentId, DataMetadata, Temperature};

use crate::ComponentDatabase;

/// Version of the embedded transport extract; reported in [`TransportDatabase::metadata`].
pub const TRANSPORT_DATA_VERSION: &str = "2026-11-24.1";
/// Identifier of the embedded transport extract.
pub const TRANSPORT_DATASET_ID: &str = "chemsep-8.3-transport";
const TRANSPORT_SCHEMA: &str = "transport-source-extract-v1";
const EMBEDDED_EXTRACT: &str =
    include_str!("../../../data/source-extracts/transport.provisional.json");

/// A ChemSep temperature correlation, T in K. Only the two equation forms the
/// transcription accepts exist; the generator rejects any other number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CorrelationForm {
    /// ChemSep 101: `exp(A + B/T + C ln T + D T^E)`.
    Eq101,
    /// ChemSep 16: `A + exp(B/T + C + D T + E T²)`.
    Eq16,
}

#[derive(Debug, Clone, Deserialize)]
struct CorrelationRecord {
    equation: u32,
    coefficients: [f64; 5],
    units: String,
    #[serde(rename = "temperature_range_K")]
    range: [f64; 2],
}

/// One pure-component temperature correlation with its source range.
#[derive(Debug, Clone, PartialEq)]
pub struct Correlation {
    pub form: CorrelationForm,
    /// A..E in source order.
    pub coefficients: [f64; 5],
    /// Source temperature range in K (inclusive).
    pub range: [f64; 2],
}

impl Correlation {
    /// Value at `t` K. Outside the source range: `OutsideValidityDomain`
    /// naming the component and the property. No extrapolation.
    pub fn evaluate(&self, id: &str, property: &str, t: f64) -> Result<f64, CalculationError> {
        let [lo, hi] = self.range;
        if !(lo..=hi).contains(&t) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "{id} {property} correlation: {lo}..={hi} K"
            )));
        }
        let [a, b, c, d, e] = self.coefficients;
        let y = match self.form {
            CorrelationForm::Eq101 => (a + b / t + c * t.ln() + d * t.powf(e)).exp(),
            CorrelationForm::Eq16 => a + (b / t + c + d * t + e * t * t).exp(),
        };
        if !y.is_finite() || y <= 0.0 {
            return Err(CalculationError::NumericalFailure(format!(
                "{id} {property} correlation is not positive at {t} K"
            )));
        }
        Ok(y)
    }
}

#[derive(Debug, Clone, Deserialize)]
struct Record {
    id: String,
    cas: String,
    #[serde(rename = "critical_volume_m3_per_mol")]
    critical_volume: f64,
    #[serde(rename = "dipole_moment_debye")]
    dipole_debye: f64,
    association_factor: f64,
    liquid_viscosity: CorrelationRecord,
    liquid_thermal_conductivity: CorrelationRecord,
    surface_tension: CorrelationRecord,
}

/// Transport data of one component.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportComponent {
    id: String,
    cas: String,
    critical_volume: f64,
    dipole_debye: f64,
    association_factor: f64,
    liquid_viscosity: Correlation,
    liquid_thermal_conductivity: Correlation,
    surface_tension: Correlation,
}

impl TransportComponent {
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn cas(&self) -> &str {
        &self.cas
    }
    /// ChemSep critical volume, m³/mol (Chung reducing volume).
    pub fn critical_volume(&self) -> f64 {
        self.critical_volume
    }
    /// ChemSep dipole moment in debye.
    pub fn dipole_moment_debye(&self) -> f64 {
        self.dipole_debye
    }
    /// Chung et al. (1988) association factor κ (zero for non-associating species).
    pub fn association_factor(&self) -> f64 {
        self.association_factor
    }
    /// Saturated-liquid viscosity correlation, Pa s.
    pub fn liquid_viscosity(&self) -> &Correlation {
        &self.liquid_viscosity
    }
    /// Saturated-liquid thermal conductivity correlation, W/(m K).
    pub fn liquid_thermal_conductivity(&self) -> &Correlation {
        &self.liquid_thermal_conductivity
    }
    /// Surface-tension correlation (liquid against its vapor), N/m.
    pub fn surface_tension(&self) -> &Correlation {
        &self.surface_tension
    }
    /// Saturated-liquid viscosity at `temperature`, Pa s.
    pub fn liquid_viscosity_at(&self, temperature: Temperature) -> Result<f64, CalculationError> {
        self.liquid_viscosity
            .evaluate(&self.id, "liquid viscosity", temperature.value())
    }
    /// Saturated-liquid thermal conductivity at `temperature`, W/(m K).
    pub fn liquid_thermal_conductivity_at(
        &self,
        temperature: Temperature,
    ) -> Result<f64, CalculationError> {
        self.liquid_thermal_conductivity.evaluate(
            &self.id,
            "liquid thermal conductivity",
            temperature.value(),
        )
    }
    /// Surface tension at `temperature`, N/m.
    pub fn surface_tension_at(&self, temperature: Temperature) -> Result<f64, CalculationError> {
        self.surface_tension
            .evaluate(&self.id, "surface tension", temperature.value())
    }
}

/// Checked embedded transport dataset. Components are matched to the component
/// database by ID and CAS number; a component without a record is missing data.
#[derive(Debug, Clone)]
pub struct TransportDatabase {
    components: Vec<TransportComponent>,
}

impl TransportDatabase {
    /// The embedded dataset, parsed and checked once per process.
    pub fn embedded() -> Result<Self, CalculationError> {
        static EMBEDDED: std::sync::OnceLock<Result<TransportDatabase, CalculationError>> =
            std::sync::OnceLock::new();
        EMBEDDED
            .get_or_init(|| Self::from_extract(EMBEDDED_EXTRACT))
            .clone()
    }

    fn from_extract(text: &str) -> Result<Self, CalculationError> {
        #[derive(Deserialize)]
        struct Extract {
            schema: String,
            version: String,
            components: Vec<Record>,
        }
        let bad = |m: String| CalculationError::MissingParameter(m);
        let extract: Extract = serde_json::from_str(text).map_err(|e| bad(e.to_string()))?;
        if extract.schema != TRANSPORT_SCHEMA || extract.version != TRANSPORT_DATA_VERSION {
            return Err(bad(format!(
                "invalid embedded transport extract header: schema {} version {} (expected {TRANSPORT_SCHEMA} {TRANSPORT_DATA_VERSION})",
                extract.schema, extract.version
            )));
        }
        let database = ComponentDatabase::embedded()?;
        let mut components: Vec<TransportComponent> = Vec::new();
        for r in extract.components {
            ComponentId::new(&r.id)?;
            let invalid = |what: &str| bad(format!("invalid transport record {}: {what}", r.id));
            if components.iter().any(|c| c.id == r.id) {
                return Err(invalid("duplicate"));
            }
            if database.get(&r.id)?.cas() != r.cas {
                return Err(invalid("CAS number differs from the component database"));
            }
            if !(r.critical_volume.is_finite() && r.critical_volume > 0.0)
                || !(r.dipole_debye.is_finite() && r.dipole_debye >= 0.0)
                || !(r.association_factor.is_finite() && r.association_factor >= 0.0)
            {
                return Err(invalid("constants"));
            }
            let correlation = |c: &CorrelationRecord, units: &str| {
                let form = match c.equation {
                    101 => CorrelationForm::Eq101,
                    16 => CorrelationForm::Eq16,
                    _ => return Err(invalid("unsupported correlation equation")),
                };
                let [lo, hi] = c.range;
                if c.units != units
                    || c.coefficients.iter().any(|v| !v.is_finite())
                    || !(lo.is_finite() && hi.is_finite() && 0.0 < lo && lo < hi)
                {
                    return Err(invalid("correlation"));
                }
                Ok(Correlation {
                    form,
                    coefficients: c.coefficients,
                    range: c.range,
                })
            };
            components.push(TransportComponent {
                liquid_viscosity: correlation(&r.liquid_viscosity, "Pa s")?,
                liquid_thermal_conductivity: correlation(
                    &r.liquid_thermal_conductivity,
                    "W/(m K)",
                )?,
                surface_tension: correlation(&r.surface_tension, "N/m")?,
                id: r.id,
                cas: r.cas,
                critical_volume: r.critical_volume,
                dipole_debye: r.dipole_debye,
                association_factor: r.association_factor,
            });
        }
        Ok(Self { components })
    }

    /// Transport record of `id`, or `MissingParameter`.
    pub fn get(&self, id: &str) -> Result<&TransportComponent, CalculationError> {
        self.components.iter().find(|c| c.id == id).ok_or_else(|| {
            CalculationError::MissingParameter(format!("no transport data for component {id}"))
        })
    }

    pub fn components(&self) -> &[TransportComponent] {
        &self.components
    }

    pub fn metadata(&self) -> DataMetadata {
        DataMetadata {
            identifier: TRANSPORT_DATASET_ID.into(),
            version: TRANSPORT_DATA_VERSION.into(),
            source: "data/source-extracts/transport.provisional.json (ChemSep 8.3 via pinned DWSIM, Artistic License 2.0; Chung 1988 association factors); upstream URL and hash retained".into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(k: f64) -> Temperature {
        Temperature::new(k).unwrap()
    }

    #[test]
    fn embedded_extract_loads_the_eight_qualified_components() {
        let db = TransportDatabase::embedded().unwrap();
        let ids: Vec<_> = db.components().iter().map(|c| c.id()).collect();
        assert_eq!(
            ids,
            [
                "methane",
                "ethane",
                "propane",
                "nitrogen",
                "carbon-dioxide",
                "oxygen",
                "water",
                "methanol"
            ]
        );
        assert!(matches!(
            db.get("hydrogen"),
            Err(CalculationError::MissingParameter(_))
        ));
        assert_eq!(db.metadata().version, TRANSPORT_DATA_VERSION);
        let water = db.get("water").unwrap();
        assert_eq!(water.association_factor(), 0.075908);
        assert!((water.dipole_moment_debye() - 1.8497).abs() < 1e-3);
        assert!((water.critical_volume() - 5.595e-5).abs() < 1e-12);
    }

    /// Hand-evaluated source equations at one temperature each.
    #[test]
    fn correlations_evaluate_the_source_forms_inside_their_ranges() {
        let db = TransportDatabase::embedded().unwrap();
        let water = db.get("water").unwrap();
        // ChemSep 101 with A..E = -133.7, 6785.7, 18.47, -1.4736e-5, 2 at 300 K.
        let ln_mu = -133.7 + 6785.7 / 300.0 + 18.47 * 300f64.ln() - 1.4736e-5 * 300.0 * 300.0;
        let mu = water.liquid_viscosity_at(t(300.0)).unwrap();
        assert!((mu - ln_mu.exp()).abs() <= 1e-15 + 1e-13 * mu);
        // ChemSep 16 for the surface tension of water at 300 K.
        let sigma_ref = -0.031819
            + (167.09f64 / 300.0 - 3.6781 + 0.0053717 * 300.0 - 8.4188e-6 * 300.0 * 300.0).exp();
        let sigma = water.surface_tension_at(t(300.0)).unwrap();
        assert!((sigma - sigma_ref).abs() <= 1e-15 + 1e-13 * sigma);
        // Plausibility only (accuracy is qualified in the model tests).
        assert!((8.0e-4..9.0e-4).contains(&mu), "{mu}");
        assert!((0.070..0.073).contains(&sigma), "{sigma}");
        assert!(matches!(
            water.liquid_viscosity_at(t(250.0)),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(matches!(
            water.liquid_thermal_conductivity_at(t(640.0)),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
    }

    #[test]
    fn stale_version_schema_and_unknown_equations_are_rejected() {
        let stale = EMBEDDED_EXTRACT.replacen(TRANSPORT_DATA_VERSION, "1999-01-01.1", 1);
        assert!(TransportDatabase::from_extract(&stale).is_err());
        let schema = EMBEDDED_EXTRACT.replacen(TRANSPORT_SCHEMA, "transport-source-extract-v0", 1);
        assert!(TransportDatabase::from_extract(&schema).is_err());
        let equation = EMBEDDED_EXTRACT.replacen("\"equation\": 101", "\"equation\": 102", 1);
        assert_ne!(equation, EMBEDDED_EXTRACT);
        assert!(matches!(
            TransportDatabase::from_extract(&equation),
            Err(CalculationError::MissingParameter(_))
        ));
        let cas = EMBEDDED_EXTRACT.replacen("\"cas\": \"74-82-8\"", "\"cas\": \"74-84-0\"", 1);
        assert!(TransportDatabase::from_extract(&cas).is_err());
    }
}
