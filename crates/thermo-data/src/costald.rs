//! COSTALD characteristic volumes and SRK acentric factors (#228).
//! Source transcription: `data/source-extracts/liquid-density-costald.provisional.json`
//! (Hankinson and Thomson 1979, Table 1, via the MIT-licensed `chemicals` 1.5.2).
//! These parameters serve only the engineering liquid density; they never enter
//! an equation of state.
use crate::{ComponentDatabase, liquid::check_cas};
use serde::Deserialize;
use thermo_types::{CalculationError, ComponentId, DataMetadata};

/// Version of the embedded COSTALD extract.
pub const COSTALD_DATA_VERSION: &str = "2026-11-24.1";
/// Identifier of the embedded COSTALD dataset.
pub const COSTALD_DATASET_ID: &str = "costald-hankinson-thomson-1979";
const SCHEMA: &str = "issue-228-costald-source-extract-v1";
const EMBEDDED_EXTRACT: &str =
    include_str!("../../../data/source-extracts/liquid-density-costald.provisional.json");

#[derive(Debug, Clone, Deserialize)]
struct Record {
    id: String,
    cas: String,
    #[serde(rename = "characteristic_volume_m3_per_mol")]
    characteristic_volume: f64,
    omega_srk: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct Extract {
    schema: String,
    version: String,
    #[serde(rename = "gas_constant_J_per_mol_K")]
    gas_constant: f64,
    components: Vec<Record>,
}

/// COSTALD parameters of one component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CostaldParameters {
    /// Characteristic volume V*, m³/mol.
    pub characteristic_volume: f64,
    /// SRK acentric factor ω_SRK (dimensionless).
    pub omega_srk: f64,
}

/// Immutable, checked COSTALD dataset. Components absent from the file are
/// missing; no generalized estimate is substituted.
#[derive(Debug, Clone)]
pub struct CostaldDatabase {
    records: Vec<Record>,
}

impl CostaldDatabase {
    /// The embedded dataset, parsed and checked once per process.
    pub fn embedded() -> Result<Self, CalculationError> {
        static EMBEDDED: std::sync::OnceLock<Result<CostaldDatabase, CalculationError>> =
            std::sync::OnceLock::new();
        EMBEDDED
            .get_or_init(|| Self::from_extract(EMBEDDED_EXTRACT, &ComponentDatabase::embedded()?))
            .clone()
    }

    /// Identifiers and CAS numbers must match `components`; V* must be finite
    /// and positive and ω_SRK finite within (-1, 2).
    fn from_extract(text: &str, components: &ComponentDatabase) -> Result<Self, CalculationError> {
        let extract: Extract = serde_json::from_str(text)
            .map_err(|e| CalculationError::MissingParameter(e.to_string()))?;
        let bad = |what: String| {
            CalculationError::MissingParameter(format!("invalid embedded COSTALD {what}"))
        };
        if extract.schema != SCHEMA
            || extract.version != COSTALD_DATA_VERSION
            || extract.gas_constant != thermo_types::MOLAR_GAS_CONSTANT
        {
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
                    "record {}: CAS {} differs from component database CAS {}",
                    r.id,
                    r.cas,
                    component.cas()
                )));
            }
            if !(r.characteristic_volume.is_finite() && r.characteristic_volume > 0.0) {
                return Err(bad(format!("record {}: characteristic volume", r.id)));
            }
            if !(r.omega_srk.is_finite() && r.omega_srk > -1.0 && r.omega_srk < 2.0) {
                return Err(bad(format!("record {}: SRK acentric factor", r.id)));
            }
        }
        Ok(Self {
            records: extract.components,
        })
    }

    /// Parameters of `id`, or `MissingParameter`.
    pub fn get(&self, id: &str) -> Result<CostaldParameters, CalculationError> {
        self.records
            .iter()
            .find(|r| r.id == id)
            .map(|r| CostaldParameters {
                characteristic_volume: r.characteristic_volume,
                omega_srk: r.omega_srk,
            })
            .ok_or_else(|| {
                CalculationError::MissingParameter(format!("COSTALD parameters for {id}"))
            })
    }

    /// Identifiers with parameters, in extract order.
    pub fn component_ids(&self) -> Vec<&str> {
        self.records.iter().map(|r| r.id.as_str()).collect()
    }

    pub fn metadata(&self) -> DataMetadata {
        DataMetadata {
            identifier: COSTALD_DATASET_ID.into(),
            version: COSTALD_DATA_VERSION.into(),
            source: "data/source-extracts/liquid-density-costald.provisional.json; Hankinson and Thomson (1979) Table 1 via chemicals 1.5.2 (MIT), URL, hash and notice retained".into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rejected(from: &str, to: &str) -> String {
        let edited = EMBEDDED_EXTRACT.replacen(from, to, 1);
        assert_ne!(edited, EMBEDDED_EXTRACT, "edit {from:?} did not apply");
        match CostaldDatabase::from_extract(&edited, &ComponentDatabase::embedded().unwrap()) {
            Err(CalculationError::MissingParameter(m)) => m,
            other => panic!("{from:?} -> {to:?} should be rejected, got {other:?}"),
        }
    }

    #[test]
    fn committed_extract_covers_every_database_component() {
        let db = CostaldDatabase::embedded().unwrap();
        let components = ComponentDatabase::embedded().unwrap();
        for c in components.components() {
            let p = db.get(c.id()).unwrap();
            assert!(p.characteristic_volume > 0.0);
        }
        // Transcribed values (Hankinson and Thomson 1979, Table 1).
        let water = db.get("water").unwrap();
        assert_eq!(water.characteristic_volume, 4.357e-5);
        assert_eq!(water.omega_srk, 0.3852);
        assert!(matches!(
            db.get("ethanol"),
            Err(CalculationError::MissingParameter(_))
        ));
        assert_eq!(db.metadata().version, COSTALD_DATA_VERSION);
    }

    #[test]
    fn malformed_records_are_rejected() {
        let m = rejected(r#""cas": "74-98-6""#, r#""cas": "74-98-7""#);
        assert!(m.contains("CAS"), "{m}");
        let m = rejected(r#""cas": "74-98-6""#, r#""cas": "74-84-0""#);
        assert!(m.contains("differs from component database"), "{m}");
        let m = rejected(r#""id": "propane""#, r#""id": "butane""#);
        assert!(m.contains("not in the component database"), "{m}");
        let m = rejected(
            r#""characteristic_volume_m3_per_mol": 0.0002001"#,
            r#""characteristic_volume_m3_per_mol": -0.0002001"#,
        );
        assert!(m.contains("characteristic volume"), "{m}");
        let m = rejected(r#""omega_srk": 0.1532"#, r#""omega_srk": 2.5"#);
        assert!(m.contains("acentric"), "{m}");
        let m = rejected(r#""version": "2026-11-24.1""#, r#""version": "0""#);
        assert!(m.contains("header"), "{m}");
    }
}
