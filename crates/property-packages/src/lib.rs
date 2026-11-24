//! One immutable model/data/interaction selection for a connected flowsheet.
//! The flowsheet package interface shares TP/PH states across EOS backends.

mod eos_cg;
pub use eos_cg::EosCgPackage;
mod flowsheet;
pub use flowsheet::{FlowsheetPackage, GergPackage, NrtlPackage};
mod pseudo;
pub use pseudo::{
    PSEUDO_CUBIC_TEMPERATURE_RANGE_K, cubic_model_with_pseudo_components,
    cubic_package_with_pseudo_components, pseudo_component_database,
};

use thermo_models::{CubicModel, HomogeneousState};
// The cubic homogeneous envelope, re-exported for interfaces that build packages
// without depending on `thermo-models` directly.
pub use thermo_models::{CUBIC_PRESSURE_RANGE_PA, CUBIC_TEMPERATURE_RANGE_K};
use thermo_solvers::{FlashPhase, PhFlashOptions, ph_flash, ps_flash, tp_flash};
use thermo_types::{
    CalculationError, CaloricReference, DataMetadata, Diagnostic, ModelMetadata, Phase, Pressure,
    Temperature, ValidityStatus,
};

/// Properties shared by homogeneous EOS and activity-package phases.
#[derive(Debug, Clone)]
pub struct PhaseProperties {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub compressibility: f64,
    pub molar_volume: f64,
    pub mass_density: f64,
    pub ln_fugacity_coefficients: Vec<f64>,
    pub residual_enthalpy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub cp: f64,
    /// Molar isochoric heat capacity in J/(mol K), when supplied by the model.
    pub cv: Option<f64>,
    /// Phase speed of sound in m/s, when supplied by the model.
    pub speed_of_sound: Option<f64>,
    /// Dynamic viscosity in Pa s, when supplied by the model.
    pub dynamic_viscosity: Option<f64>,
    /// Thermal conductivity in W/(m K), when supplied by the model.
    pub thermal_conductivity: Option<f64>,
    pub model: ModelMetadata,
    pub data: DataMetadata,
    pub caloric_reference: CaloricReference,
    pub validity: ValidityStatus,
    pub inferred_regression_interval: bool,
}

impl From<HomogeneousState> for PhaseProperties {
    fn from(s: HomogeneousState) -> Self {
        Self {
            temperature_k: s.temperature_k,
            pressure_pa: s.pressure_pa,
            compressibility: s.compressibility,
            molar_volume: s.molar_volume,
            mass_density: s.mass_density,
            ln_fugacity_coefficients: s.ln_fugacity_coefficients,
            residual_enthalpy: s.residual_enthalpy,
            enthalpy: s.enthalpy,
            entropy: s.entropy,
            cp: s.cp,
            cv: Some(s.cv),
            speed_of_sound: Some(s.speed_of_sound),
            dynamic_viscosity: None,
            thermal_conductivity: None,
            model: s.model,
            data: s.data,
            caloric_reference: s.caloric_reference,
            validity: s.validity,
            inferred_regression_interval: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EquilibriumPhase {
    pub phase: Phase,
    pub phase_fraction: f64,
    pub composition: Vec<f64>,
    pub state: PhaseProperties,
}

impl From<FlashPhase> for EquilibriumPhase {
    fn from(p: FlashPhase) -> Self {
        Self {
            phase: p.phase,
            phase_fraction: p.phase_fraction,
            composition: p.composition,
            state: p.state.into(),
        }
    }
}

/// Equilibrium state in package component order. Phase states retain model,
/// dataset, caloric-reference and validity metadata, densities, Z and entropy.
#[derive(Debug, Clone)]
pub struct EquilibriumState {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub mole_fractions: Vec<f64>,
    pub phases: Vec<EquilibriumPhase>,
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
    /// PS only: (s − s_target) / max(|s_target|, 1 J/(mol K)).
    pub scaled_entropy_residual: Option<f64>,
}

impl EquilibriumState {
    fn from_tp(r: thermo_solvers::TpFlash, fractions: &[f64]) -> Self {
        Self {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r.phases.into_iter().map(Into::into).collect(),
            iterations: r.iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: r.max_ln_fugacity_residual,
                max_component_balance_residual: Some(r.max_component_balance_residual),
                scaled_energy_residual: None,
                scaled_entropy_residual: None,
            },
            diagnostics: r.diagnostics,
        }
    }

    fn from_ph(r: thermo_solvers::PhFlash, fractions: &[f64]) -> Self {
        Self {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r.phases.into_iter().map(Into::into).collect(),
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
                scaled_entropy_residual: None,
            },
            diagnostics: r.diagnostics,
        }
    }

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

    /// Molar vapor fraction, as [`thermo_solvers::TpFlash::vapor_fraction`]: the
    /// first phase's fraction for a split, 1/0 for a single vapor/liquid phase and
    /// `None` for a single supercritical phase.
    pub fn vapor_fraction(&self) -> Option<f64> {
        (!self.phases.is_empty() && self.phases.iter().all(|p| p.phase != Phase::Supercritical))
            .then(|| {
                self.phases
                    .iter()
                    .filter(|p| p.phase == Phase::Vapor)
                    .map(|p| p.phase_fraction)
                    .sum()
            })
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
        let [t_lo, t_hi] = CUBIC_TEMPERATURE_RANGE_K;
        if !lo.is_finite() || !hi.is_finite() || lo < t_lo || hi > t_hi || lo >= hi {
            return Err(CalculationError::UnsupportedCapability(format!(
                "package requires ordered finite temperature bounds within {t_lo}..={t_hi} K"
            )));
        }
        Ok(())
    }

    /// Temperature search interval of the package, K (PH searches and the
    /// allowed TP range).
    pub fn temperature_interval(&self) -> [f64; 2] {
        [
            self.options.minimum_temperature_k,
            self.options.maximum_temperature_k,
        ]
    }

    /// Logarithmic fugacity coefficients of a single homogeneous phase at
    /// T, p and composition, on the minimum-Gibbs root. No phase split is
    /// assessed; reactors use this for the reacting fluid phase.
    pub fn ln_fugacity_coefficients(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<Vec<f64>, CalculationError> {
        self.validate_options()?;
        Ok(self
            .model
            .with_fractions(fractions)?
            .ln_fugacity(
                Temperature::new(temperature_k)?,
                Pressure::new(pressure_pa)?,
                thermo_models::RootSelection::MinimumGibbs,
            )?
            .ln_fugacity_coefficients)
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
            return Err(CalculationError::OutsideValidityDomain(
                "temperature outside package search interval".into(),
            ));
        }
        let r = tp_flash(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            self.options.tp,
        )?;
        Ok(EquilibriumState::from_tp(r, fractions))
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
        Ok(EquilibriumState::from_ph(r, fractions))
    }

    /// PS equilibrium (A14) on the same reference and bounded search interval.
    pub fn ps(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        entropy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.validate_options()?;
        let r = ps_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            entropy,
            self.options,
        )?;
        Ok(EquilibriumState {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r.phases.into_iter().map(Into::into).collect(),
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
                scaled_energy_residual: None,
                scaled_entropy_residual: Some(r.scaled_entropy_residual),
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
