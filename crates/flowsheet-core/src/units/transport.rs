//! Transport-property access for unit calculations (#225).
//!
//! Units that need a viscosity, thermal conductivity or surface tension ask
//! the context for it and fail explicitly when the selected package cannot
//! supply it at that phase state; nothing is estimated in its place.

use property_packages::{
    EquilibriumPhase, EquilibriumState, TransportProperties, annotate_transport,
};

use crate::units::CalcContext;

impl CalcContext<'_> {
    /// Transport properties of one phase through the selected package.
    pub fn transport(&self, phase: &EquilibriumPhase) -> Result<TransportProperties, String> {
        self.package()
            .transport(phase)
            .map_err(|e| format!("transport properties: {e}"))
    }

    /// Dynamic viscosity of `phase`, Pa s, or the reason it is unavailable.
    pub fn viscosity(&self, phase: &EquilibriumPhase) -> Result<f64, String> {
        self.transport(phase)?
            .viscosity
            .map_err(|e| format!("viscosity: {e}"))
    }

    /// Thermal conductivity of `phase`, W/(m K), or the reason it is unavailable.
    pub fn thermal_conductivity(&self, phase: &EquilibriumPhase) -> Result<f64, String> {
        self.transport(phase)?
            .thermal_conductivity
            .map_err(|e| format!("thermal conductivity: {e}"))
    }

    /// Surface tension of a liquid `phase`, N/m, or the reason it is unavailable.
    pub fn surface_tension(&self, phase: &EquilibriumPhase) -> Result<f64, String> {
        self.transport(phase)?
            .surface_tension
            .map_err(|e| format!("surface tension: {e}"))
    }

    /// `state` with every phase's transport properties evaluated (see
    /// [`annotate_transport`]). The engine applies this to every stored outlet.
    pub fn with_transport(&self, mut state: EquilibriumState) -> EquilibriumState {
        annotate_transport(self.package(), &mut state);
        state
    }
}
