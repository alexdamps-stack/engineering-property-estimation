//! Pinned component constants and ideal-gas NASA7 properties.
//! Source coefficients and notices are retained in `data/source-extracts`.
mod liquid;
pub use liquid::*;

use serde::Deserialize;
use thermo_types::{
    CalculationError, ComponentId, DataMetadata, MOLAR_GAS_CONSTANT as R, MolarMass, Temperature,
};

pub const DATA_VERSION: &str = "2026-11-07.1";
pub const REFERENCE_TEMPERATURE: f64 = 298.15;
pub const REFERENCE_PRESSURE: f64 = 100_000.0;

#[derive(Debug, Clone, Deserialize)]
struct Constants {
    #[serde(rename = "critical_temperature_K")]
    tc: f64,
    #[serde(rename = "critical_pressure_Pa")]
    pc: f64,
    acentric_factor: f64,
    molar_mass_kg_per_mol: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct Nasa7 {
    #[serde(rename = "temperature_breaks_K")]
    breaks: [f64; 3],
    coefficients: [[f64; 7]; 2],
}

#[derive(Debug, Clone, Deserialize)]
struct Record {
    id: String,
    cas: String,
    constants: Constants,
    ideal_gas: Nasa7,
}

/// Immutable component from the checked embedded dataset.
#[derive(Debug, Clone)]
pub struct Component(Record);

/// Molar ideal-gas quantities, h and s anchored to zero at 298.15 K.
/// Pure entropy excludes pressure and mixture terms. Units: J/mol, J/(mol K).
#[derive(Debug, Clone, Copy)]
pub struct IdealGasProperties {
    pub enthalpy: f64,
    pub entropy: f64,
    pub heat_capacity: f64,
}

impl Component {
    pub fn id(&self) -> &str {
        &self.0.id
    }
    pub fn cas(&self) -> &str {
        &self.0.cas
    }
    pub fn critical_temperature(&self) -> f64 {
        self.0.constants.tc
    }
    pub fn critical_pressure(&self) -> f64 {
        self.0.constants.pc
    }
    pub fn acentric_factor(&self) -> f64 {
        self.0.constants.acentric_factor
    }
    pub fn molar_mass(&self) -> MolarMass {
        MolarMass::new(self.0.constants.molar_mass_kg_per_mol).expect("validated embedded mass")
    }

    pub fn ideal_gas(
        &self,
        temperature: Temperature,
    ) -> Result<IdealGasProperties, CalculationError> {
        let t = temperature.value();
        let [lo, join, hi] = self.0.ideal_gas.breaks;
        if !(lo..=hi).contains(&t) {
            return Err(CalculationError::OutsideValidityDomain(format!(
                "{} NASA7: {lo}..={hi} K",
                self.id()
            )));
        }
        // Upper interval owns the join. Integrating cp across it avoids a jump
        // from the independently fitted source integration constants a6/a7.
        let low = &self.0.ideal_gas.coefficients[0];
        let high = &self.0.ideal_gas.coefficients[1];
        let integral = |entropy| {
            if t < join {
                primitive(low, t, entropy) - primitive(low, REFERENCE_TEMPERATURE, entropy)
            } else {
                primitive(low, join, entropy) - primitive(low, REFERENCE_TEMPERATURE, entropy)
                    + primitive(high, t, entropy)
                    - primitive(high, join, entropy)
            }
        };
        let a = if t < join { low } else { high };
        Ok(IdealGasProperties {
            enthalpy: R * integral(false),
            entropy: R * integral(true),
            heat_capacity: R * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])))),
        })
    }
}

fn primitive(a: &[f64; 7], t: f64, entropy: bool) -> f64 {
    if entropy {
        a[0] * t.ln()
            + a[1] * t
            + a[2] * t * t / 2.0
            + a[3] * t.powi(3) / 3.0
            + a[4] * t.powi(4) / 4.0
    } else {
        a[0] * t
            + a[1] * t * t / 2.0
            + a[2] * t.powi(3) / 3.0
            + a[3] * t.powi(4) / 4.0
            + a[4] * t.powi(5) / 5.0
    }
}

/// Immutable, checked embedded dataset. This does not qualify experimental accuracy.
#[derive(Debug, Clone)]
pub struct ComponentDatabase {
    components: Vec<Component>,
}

impl ComponentDatabase {
    pub fn embedded() -> Result<Self, CalculationError> {
        #[derive(Deserialize)]
        struct Extract {
            components: Vec<Record>,
        }
        let extract: Extract = serde_json::from_str(include_str!(
            "../../../data/source-extracts/components.provisional.json"
        ))
        .map_err(|e| CalculationError::MissingParameter(e.to_string()))?;
        for (i, c) in extract.components.iter().enumerate() {
            ComponentId::new(&c.id)?;
            if extract.components[..i].iter().any(|other| other.id == c.id)
                || !c.constants.tc.is_finite()
                || c.constants.tc <= 0.0
                || !c.constants.pc.is_finite()
                || c.constants.pc <= 0.0
                || !c.constants.acentric_factor.is_finite()
                || c.ideal_gas.breaks != [200.0, 1000.0, 6000.0]
                || c.ideal_gas
                    .coefficients
                    .iter()
                    .flatten()
                    .any(|v| !v.is_finite())
            {
                return Err(CalculationError::MissingParameter(format!(
                    "invalid embedded record {}",
                    c.id
                )));
            }
            MolarMass::new(c.constants.molar_mass_kg_per_mol)?;
        }
        Ok(Self {
            components: extract.components.into_iter().map(Component).collect(),
        })
    }

    pub fn get(&self, id: &str) -> Result<&Component, CalculationError> {
        self.components
            .iter()
            .find(|c| c.id() == id)
            .ok_or_else(|| CalculationError::MissingParameter(format!("component {id}")))
    }
    pub fn components(&self) -> &[Component] {
        &self.components
    }
    pub fn metadata(&self) -> DataMetadata {
        DataMetadata { identifier: "coolprop-7.2-cubic-cantera-3.2-nasa7".into(), version: DATA_VERSION.into(),
            source: "data/source-extracts/components.provisional.json; upstream URLs, hashes and notices retained".into() }
    }
}
