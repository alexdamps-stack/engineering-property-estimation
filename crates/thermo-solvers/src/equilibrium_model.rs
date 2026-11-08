//! Internal contract shared by cubic and native GERG equilibrium adapters.
use thermo_models::{CubicModel, FugacityState, HomogeneousState, RootSelection};
use thermo_types::{CalculationError, Phase, Pressure, Temperature};

pub(crate) trait EquilibriumModel: Sized {
    fn fractions(&self) -> &[f64];
    fn with_fractions(&self, fractions: &[f64]) -> Result<Self, CalculationError>;
    fn ln_fugacity(
        &self,
        t: Temperature,
        p: Pressure,
        selection: RootSelection,
    ) -> Result<FugacityState, CalculationError>;
    fn evaluate(
        &self,
        t: Temperature,
        p: Pressure,
        selection: RootSelection,
    ) -> Result<HomogeneousState, CalculationError>;
    fn wilson_values(&self, t: Temperature, p: Pressure) -> Vec<f64>;
    fn phase_label(&self, state: &HomogeneousState) -> Phase;
}

impl EquilibriumModel for CubicModel {
    fn fractions(&self) -> &[f64] {
        self.fractions()
    }
    fn with_fractions(&self, x: &[f64]) -> Result<Self, CalculationError> {
        self.with_fractions(x)
    }
    fn ln_fugacity(
        &self,
        t: Temperature,
        p: Pressure,
        s: RootSelection,
    ) -> Result<FugacityState, CalculationError> {
        self.ln_fugacity(t, p, s)
    }
    fn evaluate(
        &self,
        t: Temperature,
        p: Pressure,
        s: RootSelection,
    ) -> Result<HomogeneousState, CalculationError> {
        self.evaluate(t, p, s)
    }
    fn wilson_values(&self, t: Temperature, p: Pressure) -> Vec<f64> {
        crate::wilson_k(self, t, p)
    }
    fn phase_label(&self, s: &HomogeneousState) -> Phase {
        crate::single_phase_label(self, s)
    }
}
