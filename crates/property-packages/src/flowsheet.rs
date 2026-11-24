use std::{fmt::Debug, sync::Arc};

use thermo_models::{ACTIVITY_TEMPERATURE_RANGE_K, ActivityModel, GergModel};
use thermo_solvers::{
    ActivityFlashPhase, ActivityOptions, GergPhFlashOptions, QUALIFIED_DRY_COMPONENTS,
    activity_ph_flash, activity_tp_flash, gerg_ph_flash, gerg_tp_flash,
};
use thermo_types::{CalculationError, MOLAR_GAS_CONSTANT, Phase, Pressure, Temperature};

use crate::{
    Convergence, CubicPackage, EngineeringLiquidDensity, EquilibriumPhase, EquilibriumState,
    PhaseProperties,
};

/// Equilibrium operations required by a sequential flowsheet.
pub trait FlowsheetPackage: Debug + Send + Sync {
    fn identity(&self) -> &Arc<()>;
    fn component_ids(&self) -> Vec<&str>;
    fn temperature_interval(&self) -> [f64; 2];
    fn validate_composition(&self, fractions: &[f64]) -> Result<(), CalculationError>;
    fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, CalculationError>;
    fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, CalculationError>;

    /// Access for legacy unit models that depend on cubic-specific operations.
    fn cubic(&self) -> Result<&CubicPackage, CalculationError>;

    /// Engineering liquid density of one liquid phase at T, p and phase mole
    /// fractions (#228). It is reported beside the phase's EOS density and
    /// never enters equilibrium or caloric results. Packages without one
    /// return `UnsupportedCapability`; nothing is substituted.
    fn engineering_liquid_density(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _pressure_pa: f64,
    ) -> Result<EngineeringLiquidDensity, CalculationError> {
        Err(CalculationError::UnsupportedCapability(
            "this property package provides no engineering liquid density".into(),
        ))
    }
}

/// The bounded NRTL liquid and ideal-gas vapor package.
#[derive(Debug, Clone)]
pub struct NrtlPackage {
    model: ActivityModel,
    options: ActivityOptions,
    identity: Arc<()>,
}

impl NrtlPackage {
    pub const TEMPERATURE_RANGE_K: [f64; 2] = ACTIVITY_TEMPERATURE_RANGE_K;
    pub const COMPONENTS: [&str; 2] = ["methanol", "water"];

    pub fn new(model: ActivityModel, interval: [f64; 2]) -> Result<Self, CalculationError> {
        let [lo, hi] = interval;
        let range = Self::TEMPERATURE_RANGE_K[0]..=Self::TEMPERATURE_RANGE_K[1];
        if !range.contains(&lo) || !range.contains(&hi) || lo >= hi {
            return Err(CalculationError::OutsideValidityDomain(
                "NRTL-IG requires ordered temperature bounds within 273.16..=450 K".into(),
            ));
        }
        let options = ActivityOptions {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            ..Default::default()
        };
        Ok(Self {
            model,
            options,
            identity: Arc::new(()),
        })
    }

    fn phase(&self, p: ActivityFlashPhase) -> Result<EquilibriumPhase, CalculationError> {
        let s = p.state;
        let (cv, speed_of_sound) = if p.phase == Phase::Vapor {
            // Ideal-gas acoustic derivative at fixed composition.
            let cv = s.cp - MOLAR_GAS_CONSTANT;
            let sound = (s.cp / cv * s.pressure_pa / s.mass_density).sqrt();
            if !cv.is_finite() || cv <= 0.0 || !sound.is_finite() || sound <= 0.0 {
                return Err(CalculationError::NumericalFailure(
                    "invalid ideal-gas heat capacity or speed of sound".into(),
                ));
            }
            (Some(cv), Some(sound))
        } else {
            (None, None)
        };
        let residual_enthalpy = if p.phase == Phase::Vapor {
            0.0
        } else {
            let ideal = self.model.with_fractions(&p.composition)?.vapor(
                Temperature::new(s.temperature_k)?,
                Pressure::new(s.pressure_pa)?,
            )?;
            s.enthalpy - ideal.enthalpy
        };
        Ok(EquilibriumPhase {
            phase: p.phase,
            phase_fraction: p.phase_fraction,
            composition: p.composition,
            state: PhaseProperties {
                temperature_k: s.temperature_k,
                pressure_pa: s.pressure_pa,
                compressibility: s.pressure_pa * s.molar_volume
                    / (MOLAR_GAS_CONSTANT * s.temperature_k),
                molar_volume: s.molar_volume,
                mass_density: s.mass_density,
                ln_fugacity_coefficients: s.ln_fugacity_coefficients,
                residual_enthalpy,
                enthalpy: s.enthalpy,
                entropy: s.entropy,
                cp: s.cp,
                cv,
                speed_of_sound,
                dynamic_viscosity: None,
                thermal_conductivity: None,
                model: s.model,
                data: s.data,
                caloric_reference: s.caloric_reference,
                validity: s.validity,
                inferred_regression_interval: s.inferred_regression_interval,
            },
        })
    }
}

impl FlowsheetPackage for NrtlPackage {
    fn identity(&self) -> &Arc<()> {
        &self.identity
    }

    fn component_ids(&self) -> Vec<&str> {
        self.model.component_ids()
    }

    fn temperature_interval(&self) -> [f64; 2] {
        [
            self.options.minimum_temperature_k,
            self.options.maximum_temperature_k,
        ]
    }

    fn validate_composition(&self, fractions: &[f64]) -> Result<(), CalculationError> {
        self.model.with_fractions(fractions).map(|_| ())
    }

    fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let [lo, hi] = self.temperature_interval();
        if !(lo..=hi).contains(&temperature_k) {
            return Err(CalculationError::OutsideValidityDomain(
                "temperature outside NRTL-IG package search interval".into(),
            ));
        }
        let r = activity_tp_flash(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            self.options,
        )?;
        Ok(EquilibriumState {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r
                .phases
                .into_iter()
                .map(|p| self.phase(p))
                .collect::<Result<_, _>>()?,
            iterations: r.iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: r.max_ln_fugacity_residual,
                max_component_balance_residual: Some(r.max_component_balance_residual),
                scaled_energy_residual: None,
                scaled_entropy_residual: None,
            },
            diagnostics: r.diagnostics,
        })
    }

    fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let r = activity_ph_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            enthalpy,
            self.options,
        )?;
        Ok(EquilibriumState {
            temperature_k: r.temperature_k,
            pressure_pa: r.pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: r
                .phases
                .into_iter()
                .map(|p| self.phase(p))
                .collect::<Result<_, _>>()?,
            iterations: r.iterations,
            convergence: Convergence {
                max_ln_fugacity_residual: r
                    .tp
                    .as_ref()
                    .map_or(0.0, |tp| tp.max_ln_fugacity_residual),
                max_component_balance_residual: r
                    .tp
                    .as_ref()
                    .map(|tp| tp.max_component_balance_residual),
                scaled_energy_residual: Some(r.scaled_energy_residual),
                scaled_entropy_residual: None,
            },
            diagnostics: r.diagnostics,
        })
    }

    fn cubic(&self) -> Result<&CubicPackage, CalculationError> {
        Err(CalculationError::UnsupportedCapability(
            "this unit's thermodynamic calculation requires a cubic package; selected NRTL-IG"
                .into(),
        ))
    }
}

impl FlowsheetPackage for CubicPackage {
    fn identity(&self) -> &Arc<()> {
        &self.identity
    }

    fn component_ids(&self) -> Vec<&str> {
        self.model().component_ids()
    }

    fn temperature_interval(&self) -> [f64; 2] {
        self.temperature_interval()
    }

    fn validate_composition(&self, fractions: &[f64]) -> Result<(), CalculationError> {
        self.validate_composition(fractions)
    }

    fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.tp(fractions, temperature_k, pressure_pa)
    }

    fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.ph(fractions, pressure_pa, enthalpy)
    }

    fn cubic(&self) -> Result<&CubicPackage, CalculationError> {
        Ok(self)
    }

    fn engineering_liquid_density(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EngineeringLiquidDensity, CalculationError> {
        CubicPackage::engineering_liquid_density(self, fractions, temperature_k, pressure_pa)
    }
}

/// GERG-2008 on its qualified dry-gas equilibrium domain and native caloric reference.
#[derive(Debug, Clone)]
pub struct GergPackage {
    model: GergModel,
    options: GergPhFlashOptions,
    identity: Arc<()>,
}

impl GergPackage {
    pub const TEMPERATURE_RANGE_K: [f64; 2] = [230.0, 350.0];
    pub const COMPONENTS: [&str; 5] = QUALIFIED_DRY_COMPONENTS;

    pub fn new(model: GergModel, temperature_interval: [f64; 2]) -> Result<Self, CalculationError> {
        if model
            .component_ids()
            .iter()
            .any(|id| !Self::COMPONENTS.contains(id))
        {
            return Err(CalculationError::UnsupportedCapability(format!(
                "GERG-2008 flowsheet components must be selected from {}",
                Self::COMPONENTS.join(", ")
            )));
        }
        let [lo, hi] = temperature_interval;
        let range = Self::TEMPERATURE_RANGE_K[0]..=Self::TEMPERATURE_RANGE_K[1];
        if !range.contains(&lo) || !range.contains(&hi) || lo >= hi {
            return Err(CalculationError::OutsideValidityDomain(
                "GERG-2008 requires ordered temperature bounds within 230..=350 K".into(),
            ));
        }
        let options = GergPhFlashOptions {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            ..Default::default()
        };
        Ok(Self {
            model,
            options,
            identity: Arc::new(()),
        })
    }
}

impl FlowsheetPackage for GergPackage {
    fn identity(&self) -> &Arc<()> {
        &self.identity
    }

    fn component_ids(&self) -> Vec<&str> {
        self.model.component_ids()
    }

    fn temperature_interval(&self) -> [f64; 2] {
        [
            self.options.minimum_temperature_k,
            self.options.maximum_temperature_k,
        ]
    }

    fn validate_composition(&self, fractions: &[f64]) -> Result<(), CalculationError> {
        self.model.with_fractions(fractions).map(|_| ())
    }

    fn tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let [lo, hi] = self.temperature_interval();
        if !(lo..=hi).contains(&temperature_k) {
            return Err(CalculationError::OutsideValidityDomain(
                "temperature outside GERG-2008 package search interval".into(),
            ));
        }
        let result = gerg_tp_flash(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            self.options.tp,
        )?;
        Ok(EquilibriumState::from_tp(result, fractions))
    }

    fn ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let result = gerg_ph_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            enthalpy,
            self.options,
        )?;
        Ok(EquilibriumState::from_ph(result, fractions))
    }

    fn cubic(&self) -> Result<&CubicPackage, CalculationError> {
        Err(CalculationError::UnsupportedCapability(
            "this unit's thermodynamic calculation requires a cubic package; selected GERG-2008"
                .into(),
        ))
    }
}
