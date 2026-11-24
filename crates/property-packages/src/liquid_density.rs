//! Engineering liquid density of the cubic packages (#228).
//!
//! The PR, PR-BM and SRK packages keep their untranslated EOS volume for
//! every equilibrium, fugacity and caloric result. A package may in addition
//! carry a separately labelled COSTALD engineering density
//! ([`thermo_models::CostaldLiquidDensity`]); it is evaluated on request for
//! a liquid phase and never feeds back into the EOS results
//! (`docs/calculation-contract.md`, *Liquid density and caloric coupling*).
use std::sync::Arc;

use thermo_models::{CostaldLiquidDensity, EngineeringLiquidDensity};
use thermo_types::{CalculationError, InputError};

use crate::CubicPackage;

impl CubicPackage {
    /// The same package with a COSTALD engineering liquid density for its
    /// components. The correlation must list the package components in
    /// package order. Equilibrium and caloric results are unaffected.
    pub fn with_engineering_liquid_density(
        mut self,
        density: CostaldLiquidDensity,
    ) -> Result<Self, CalculationError> {
        if density.component_ids() != self.model().component_ids() {
            return Err(InputError::ComponentMismatch.into());
        }
        self.liquid_density = Some(Arc::new(density));
        Ok(self)
    }

    /// Engineering liquid density of a liquid phase with mole fractions
    /// `fractions` (package order) at T and p. `UnsupportedCapability` when
    /// the package carries no engineering density; states outside the
    /// correlation's qualified domain return `OutsideValidityDomain`.
    pub fn engineering_liquid_density(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EngineeringLiquidDensity, CalculationError> {
        self.liquid_density
            .as_ref()
            .ok_or_else(|| {
                CalculationError::UnsupportedCapability(
                    "this cubic package was built without an engineering liquid density".into(),
                )
            })?
            .evaluate(fractions, temperature_k, pressure_pa)
    }
}
