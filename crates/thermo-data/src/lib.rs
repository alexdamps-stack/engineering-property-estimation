//! Embedded, versioned constants. No binary interaction parameters are defaulted.

use serde::Deserialize;
use thermo_types::{Component, ThermoError};

#[derive(Debug, Deserialize)]
pub struct Catalog {
    pub version: String,
    pub source: String,
    pub source_url: String,
    pub limitations: String,
    components: Vec<ComponentRecord>,
}

#[derive(Debug, Deserialize)]
struct ComponentRecord {
    id: String,
    critical_temperature_k: f64,
    critical_pressure_pa: f64,
    molar_mass_kg_mol: f64,
    acentric_factor: f64,
}

impl Catalog {
    pub fn bundled() -> Result<Self, ThermoError> {
        serde_json::from_str(include_str!("../../../data/components.json"))
            .map_err(|e| ThermoError::MissingData(format!("invalid bundled catalog: {e}")))
    }

    pub fn component(&self, id: &str) -> Result<Component, ThermoError> {
        let c = self
            .components
            .iter()
            .find(|c| c.id == id)
            .ok_or_else(|| ThermoError::MissingData(format!("unknown component ID: {id}")))?;
        Component::new(
            &c.id,
            c.critical_temperature_k,
            c.critical_pressure_pa,
            c.molar_mass_kg_mol,
            c.acentric_factor,
        )
    }

    pub fn component_ids(&self) -> impl Iterator<Item = &str> {
        self.components.iter().map(|c| c.id.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bundled_records_are_valid_and_unknown_ids_are_explicit() {
        let catalog = Catalog::bundled().unwrap();
        assert_eq!(catalog.component_ids().count(), 8);
        for id in catalog.component_ids() {
            assert!(catalog.component(id).is_ok());
        }
        assert!(matches!(
            catalog.component("unobtainium"),
            Err(ThermoError::MissingData(_))
        ));
    }
}
