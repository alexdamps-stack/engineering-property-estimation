//! JSON-compatible homogeneous requests routed through validated models and data.

use serde::{Deserialize, Serialize};
use thermo_data::Catalog;
use thermo_models::{HomogeneousProperties, MODEL_VERSION};
use thermo_types::{
    CompositionBasis, CubicModel, Mixture, Pressure, RootSelection, Temperature, ThermoError,
};

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HomogeneousRequest {
    pub model: CubicModel,
    pub components: Vec<String>,
    pub fractions: Vec<f64>,
    pub basis: CompositionBasis,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub binary_interactions: Vec<Vec<f64>>,
    pub root: RootSelection,
}

#[derive(Debug, Serialize)]
pub struct HomogeneousResult {
    pub schema_version: u32,
    pub model: CubicModel,
    pub model_version: &'static str,
    pub data_version: String,
    pub data_source: String,
    pub components: Vec<String>,
    pub mole_fractions: Vec<f64>,
    pub binary_interactions: Vec<Vec<f64>>,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub selected_root: RootSelection,
    pub phase_stability: &'static str,
    pub validity: &'static str,
    pub warnings: Vec<&'static str>,
    pub properties: HomogeneousProperties,
}

pub fn calculate(request: HomogeneousRequest) -> Result<HomogeneousResult, ThermoError> {
    let temperature = Temperature::new(request.temperature_k)?;
    let pressure = Pressure::new(request.pressure_pa)?;
    let catalog = Catalog::bundled()?;
    let components = request
        .components
        .iter()
        .map(|id| catalog.component(id))
        .collect::<Result<Vec<_>, _>>()?;
    let mixture = Mixture::new(
        components,
        request.fractions,
        request.basis,
        request.binary_interactions,
    )?;
    let properties =
        thermo_models::evaluate(request.model, &mixture, temperature, pressure, request.root)?;
    Ok(HomogeneousResult {
        schema_version: 1,
        model: request.model,
        model_version: MODEL_VERSION,
        data_version: catalog.version,
        data_source: catalog.source,
        components: request.components,
        mole_fractions: mixture.mole_fractions().to_vec(),
        binary_interactions: mixture.binary_interactions().to_vec(),
        temperature_k: temperature.value(),
        pressure_pa: pressure.value(),
        selected_root: request.root,
        phase_stability: "not_assessed",
        validity: "implementation_checked_at_reference_cases_only",
        warnings: vec![
            "Homogeneous algebraic root; not an equilibrium flash or phase identification.",
            "Residual enthalpy/entropy only, relative to an ideal gas at the same T, P, and composition; no total caloric properties.",
            "Input kij values are explicit user choices, not qualified fitted interaction data.",
            "No experimental accuracy or full-domain validation claim.",
        ],
        properties,
    })
}
