//! Phase transport properties for flowsheet packages (#225).
//!
//! The methods are package-independent side properties: the gas method uses
//! the phase's own molar volume from the package model, the liquid method
//! uses temperature, pressure and composition. Nothing here changes an
//! equilibrium or caloric result.
pub use thermo_models::{TransportModel, TransportProperties, TransportRegime};
use thermo_types::{CalculationError, Phase, Temperature};

use crate::{EquilibriumPhase, EquilibriumState, FlowsheetPackage};

/// Transport properties of `phase`, whose composition is in the order of
/// `component_ids`. Liquid phases use the liquid correlations; vapor and
/// supercritical phases use the Chung method. A component without transport
/// data is `MissingParameter`.
pub fn phase_transport(
    component_ids: &[&str],
    phase: &EquilibriumPhase,
) -> Result<TransportProperties, CalculationError> {
    let model = TransportModel::new(component_ids)?;
    let temperature = Temperature::new(phase.state.temperature_k)?;
    match phase.phase {
        Phase::Liquid => model.liquid(&phase.composition, temperature, phase.state.pressure_pa),
        Phase::Vapor | Phase::Supercritical => {
            model.gas(&phase.composition, temperature, phase.state.molar_volume)
        }
    }
}

/// Evaluates [`FlowsheetPackage::transport`] for every phase of `state` and
/// stores the result in `PhaseProperties::transport`, with the available
/// viscosity and thermal conductivity in their own fields. Unavailable
/// properties stay `None`; the state is otherwise unchanged.
pub fn annotate_transport<P: FlowsheetPackage + ?Sized>(package: &P, state: &mut EquilibriumState) {
    for phase in &mut state.phases {
        let result = package.transport(phase);
        let (viscosity, conductivity) = match &result {
            Ok(t) => (
                t.viscosity.clone().ok(),
                t.thermal_conductivity.clone().ok(),
            ),
            Err(_) => (None, None),
        };
        phase.state.dynamic_viscosity = viscosity;
        phase.state.thermal_conductivity = conductivity;
        phase.state.transport = Some(result);
    }
}
