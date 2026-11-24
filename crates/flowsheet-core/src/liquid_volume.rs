//! Engineering liquid volumes of streams and phases (#228).
//!
//! Volumetric liquid quantities reported by units (pump hydraulics,
//! volumetric flow) use the package's separately labelled engineering liquid
//! density, never the EOS volume, and are *Unavailable* when the package has
//! none or the state is outside its domain. The EOS density stays on each
//! phase (`PhaseProperties::mass_density`) and is reported beside it.
use property_packages::{EngineeringLiquidDensity, EquilibriumState, FlowsheetPackage};
use thermo_types::{CalculationError, Phase};

use crate::engine::MaterialStream;

/// Engineering density of every phase of `state`, in phase order: `None` for
/// vapor and supercritical phases, the package's result for liquid phases.
pub fn phase_engineering_densities(
    package: &dyn FlowsheetPackage,
    state: &EquilibriumState,
) -> Vec<Option<Result<EngineeringLiquidDensity, CalculationError>>> {
    state
        .phases
        .iter()
        .map(|p| {
            (p.phase == Phase::Liquid).then(|| {
                package.engineering_liquid_density(
                    &p.composition,
                    p.state.temperature_k,
                    p.state.pressure_pa,
                )
            })
        })
        .collect()
}

/// Engineering liquid molar volume per mole of the whole state, m³/mol:
/// Σ over liquid phases of phase fraction × engineering molar volume.
/// `Ok(None)` when the state has no liquid phase with a positive fraction.
pub fn engineering_liquid_volume(
    package: &dyn FlowsheetPackage,
    state: &EquilibriumState,
) -> Result<Option<f64>, CalculationError> {
    let mut volume = None;
    for (p, density) in state
        .phases
        .iter()
        .zip(phase_engineering_densities(package, state))
    {
        if let Some(density) = density
            && p.phase_fraction > 0.0
        {
            *volume.get_or_insert(0.0) += p.phase_fraction * density?.molar_volume;
        }
    }
    Ok(volume)
}

impl MaterialStream {
    /// Volumetric flow of the liquid phases at stream conditions on the
    /// engineering liquid density, m³/s. `Ok(None)` without a liquid phase
    /// (or without flow state); an error when the package cannot supply the
    /// engineering density of a liquid phase.
    pub fn engineering_liquid_volume_flow_m3_s(
        &self,
        package: &dyn FlowsheetPackage,
    ) -> Result<Option<f64>, CalculationError> {
        let Some(state) = &self.state else {
            return Ok(None);
        };
        Ok(engineering_liquid_volume(package, state)?.map(|v| self.molar_flow_mol_s * v))
    }
}
