//! One immutable model/data/interaction selection for a connected flowsheet.
//! The initial adapter supports exact PR and PR-BM on their existing TP/PH domain.

use thermo_models::CubicModel;
use thermo_solvers::{FlashPhase, PhFlashOptions, ph_flash, tp_flash};
use thermo_types::{CalculationError, Diagnostic, Pressure, Temperature};

/// Equilibrium state in package component order. Phase states retain model,
/// dataset, caloric-reference and validity metadata, densities, Z and entropy.
#[derive(Debug, Clone)]
pub struct EquilibriumState {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub mole_fractions: Vec<f64>,
    pub phases: Vec<FlashPhase>,
    pub iterations: usize,
    pub convergence: Convergence,
    pub diagnostics: Vec<Diagnostic>,
}

/// Residuals retained from the actual equilibrium solve. PH pure coexistence
/// has a saturation fugacity residual and no TP split component residual.
#[derive(Debug, Clone)]
pub struct Convergence {
    pub max_ln_fugacity_residual: f64,
    pub max_component_balance_residual: Option<f64>,
    pub scaled_energy_residual: Option<f64>,
}

impl EquilibriumState {
    /// Phase-weighted molar enthalpy, J/mol, on the package caloric reference.
    pub fn enthalpy(&self) -> f64 {
        self.phases
            .iter()
            .map(|p| p.phase_fraction * p.state.enthalpy)
            .sum()
    }

    /// Phase-weighted molar entropy, J/(mol K).
    pub fn entropy(&self) -> f64 {
        self.phases
            .iter()
            .map(|p| p.phase_fraction * p.state.entropy)
            .sum()
    }
}

/// Fixed cubic package. Component order, interactions and source are supplied
/// explicitly through the validated model; there are no implicit zero pairs.
#[derive(Debug, Clone)]
pub struct CubicPackage {
    model: CubicModel,
    options: PhFlashOptions,
    identity: std::sync::Arc<()>,
}

impl CubicPackage {
    pub fn new(model: CubicModel, options: PhFlashOptions) -> Self {
        Self {
            model,
            options,
            identity: std::sync::Arc::new(()),
        }
    }

    pub fn model(&self) -> &CubicModel {
        &self.model
    }

    /// Clones share an immutable configuration; any newly selected package
    /// conservatively invalidates old calculation results.
    pub fn same_configuration(&self, other: &Self) -> bool {
        std::sync::Arc::ptr_eq(&self.identity, &other.identity)
    }

    pub fn validate_composition(&self, fractions: &[f64]) -> Result<(), CalculationError> {
        self.validate_options()?;
        self.model.with_fractions(fractions).map(|_| ())
    }

    fn validate_options(&self) -> Result<(), CalculationError> {
        let lo = self.options.minimum_temperature_k;
        let hi = self.options.maximum_temperature_k;
        if !lo.is_finite() || !hi.is_finite() || lo < 200.0 || hi > 600.0 || lo >= hi {
            return Err(CalculationError::UnsupportedCapability(
                "package requires ordered finite temperature bounds within 200..=600 K".into(),
            ));
        }
        Ok(())
    }

    /// TP equilibrium, with all solver errors and diagnostics propagated.
    pub fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.validate_options()?;
        if temperature_k < self.options.minimum_temperature_k
            || temperature_k > self.options.maximum_temperature_k
        {
            return Err(CalculationError::UnsupportedCapability(
                "temperature outside package search interval".into(),
            ));
        }
        let r = tp_flash(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            self.options.tp,
        )?;
        Ok(EquilibriumState {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r.phases,
            iterations: r.iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: r.max_ln_fugacity_residual,
                max_component_balance_residual: Some(r.max_component_balance_residual),
                scaled_energy_residual: None,
            },
            diagnostics: r.diagnostics,
        })
    }

    /// PH equilibrium on the same reference and bounded search interval.
    pub fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.validate_options()?;
        let r = ph_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            enthalpy,
            self.options,
        )?;
        Ok(EquilibriumState {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r.phases,
            iterations: r.iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: r.tp.as_ref().map_or_else(
                    || {
                        r.pure_saturation
                            .as_ref()
                            .map_or(0.0, |s| s.fugacity_residual)
                    },
                    |tp| tp.max_ln_fugacity_residual,
                ),
                max_component_balance_residual: r
                    .tp
                    .as_ref()
                    .map(|tp| tp.max_component_balance_residual),
                scaled_energy_residual: Some(r.scaled_energy_residual),
            },
            diagnostics: r.diagnostics,
        })
    }
}

impl CubicPackage {
    /// Bubble/dew pressure in the package's temperature interval and explicitly
    /// bounded A13 dry-fluid domain. Specified fractions belong to the material phase.
    pub fn saturation_pressure(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        kind: thermo_solvers::CubicSaturationKind,
        bounds: thermo_solvers::SaturationOptions,
    ) -> Result<thermo_solvers::SaturationPoint, CalculationError> {
        self.validate_options()?;
        let bounds = self.saturation_bounds(bounds)?;
        thermo_solvers::saturation_pressure(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            kind,
            bounds,
        )
    }
    pub fn saturation_temperature(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        kind: thermo_solvers::CubicSaturationKind,
        bounds: thermo_solvers::SaturationOptions,
    ) -> Result<thermo_solvers::SaturationPoint, CalculationError> {
        self.validate_options()?;
        let bounds = self.saturation_bounds(bounds)?;
        thermo_solvers::saturation_temperature(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            kind,
            bounds,
        )
    }
    pub fn phase_envelope(
        &self,
        fractions: &[f64],
        kind: thermo_solvers::CubicSaturationKind,
        step_k: f64,
        maximum_points: usize,
        bounds: thermo_solvers::SaturationOptions,
    ) -> Result<thermo_solvers::PhaseEnvelope, CalculationError> {
        self.validate_options()?;
        let bounds = self.saturation_bounds(bounds)?;
        thermo_solvers::trace_phase_envelope(
            &self.model.with_fractions(fractions)?,
            kind,
            step_k,
            maximum_points,
            bounds,
        )
    }
    fn saturation_bounds(
        &self,
        bounds: thermo_solvers::SaturationOptions,
    ) -> Result<thermo_solvers::SaturationOptions, CalculationError> {
        bounds.validate()?;
        if bounds.minimum_temperature_k < self.options.minimum_temperature_k
            || bounds.maximum_temperature_k > self.options.maximum_temperature_k
        {
            return Err(CalculationError::OutsideValidityDomain(
                "A13 search bounds outside package temperature interval".into(),
            ));
        }
        Ok(bounds)
    }
}
