use std::{fmt::Debug, sync::Arc};

use thermo_data::Component;
use thermo_models::{ACTIVITY_TEMPERATURE_RANGE_K, ActivityModel, GergModel, RootSelection};
use thermo_solvers::{
    ActivityFlashPhase, ActivityOptions, ActivitySaturation, GergPhFlashOptions,
    QUALIFIED_DRY_COMPONENTS, SaturationOptions, activity_bubble_pressure,
    activity_bubble_temperature, activity_dew_pressure, activity_dew_temperature,
    activity_ph_flash, activity_ps_flash, activity_tp_flash, gerg_homogeneous, gerg_ln_fugacity,
    gerg_ph_flash, gerg_ps_flash, gerg_saturation_pressure, gerg_saturation_temperature,
    gerg_tp_flash, liquid_stability,
};
use thermo_types::{
    CalculationError, Diagnostic, DiagnosticSeverity, MOLAR_GAS_CONSTANT, Phase, Pressure,
    Temperature,
};

use crate::{
    Convergence, CubicPackage, EngineeringLiquidDensity, EquilibriumPhase, EquilibriumState,
    PhaseProperties, SaturationKind, SaturationState,
};

mod multiphase;
pub use multiphase::{MultiphaseEquilibrium, PhaseRole};

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
    /// No unit uses it; new code must call package operations instead.
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
        Err(unsupported("an engineering liquid density"))
    }

    /// Viscosity, thermal conductivity and surface tension of one phase of a
    /// state returned by this package (#225). The provided method applies the
    /// package-independent transport methods to the phase's own state
    /// (molar volume from this package's model); a package without a usable
    /// phase volume must override it and return `UnsupportedCapability`.
    fn transport(
        &self,
        phase: &EquilibriumPhase,
    ) -> Result<crate::TransportProperties, CalculationError> {
        crate::phase_transport(&self.component_ids(), phase)
    }

    // --- Operations beyond TP/PH (#226). Each package overrides the ones it
    // supports; the defaults fail explicitly and never substitute a model.

    /// PS equilibrium at p and molar entropy on the package caloric reference,
    /// with the structure of PH (including the pure saturation lever rule).
    fn ps(
        &self,
        _fractions: &[f64],
        _pressure_pa: f64,
        _entropy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        Err(unsupported("PS flash"))
    }

    /// Properties of one homogeneous phase of composition `fractions` at T, p
    /// on the selected root. No phase-stability claim.
    fn homogeneous(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _pressure_pa: f64,
        _root: RootSelection,
    ) -> Result<PhaseProperties, CalculationError> {
        Err(unsupported("homogeneous phase evaluation"))
    }

    /// ln φ per component (package order) of one homogeneous phase on the
    /// selected root, from the package's own phase model.
    fn ln_fugacity_coefficients(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _pressure_pa: f64,
        _root: RootSelection,
    ) -> Result<Vec<f64>, CalculationError> {
        Err(unsupported("fugacity coefficients"))
    }

    /// Equilibrium ratios K_i = φ_i(x, liquid root) / φ_i(y, vapor root) at T
    /// and p, from [`Self::ln_fugacity_coefficients`]. A successive-substitution
    /// building block; it does not check phase existence or stability.
    fn k_values(
        &self,
        x: &[f64],
        y: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<Vec<f64>, CalculationError> {
        let liquid = self.ln_fugacity_coefficients(
            x,
            temperature_k,
            pressure_pa,
            RootSelection::SmallestVolume,
        )?;
        let vapor = self.ln_fugacity_coefficients(
            y,
            temperature_k,
            pressure_pa,
            RootSelection::LargestVolume,
        )?;
        let k: Vec<f64> = liquid
            .iter()
            .zip(&vapor)
            .map(|(l, v)| (l - v).exp())
            .collect();
        if k.iter().all(|v| v.is_finite() && *v > 0.0) {
            Ok(k)
        } else {
            Err(CalculationError::NumericalFailure(
                "nonfinite equilibrium ratio".into(),
            ))
        }
    }

    /// Temperature interval, K, searched by the package's bubble/dew solver.
    fn saturation_temperature_interval(&self) -> Result<[f64; 2], CalculationError> {
        Err(unsupported("bubble/dew calculation"))
    }

    /// Bubble or dew pressure at T. `fractions` are the material phase: the
    /// liquid for a bubble point, the vapor for a dew point.
    fn saturation_pressure(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        Err(unsupported("bubble/dew pressure"))
    }

    /// Bubble or dew temperature at p, inside
    /// [`Self::saturation_temperature_interval`].
    fn saturation_temperature(
        &self,
        _fractions: &[f64],
        _pressure_pa: f64,
        _kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        Err(unsupported("bubble/dew temperature"))
    }

    /// Components with element and formation data, in package order, when the
    /// package enthalpy is on the formation-compatible ideal-gas reference
    /// (pure ideal gas h = s = 0 at 298.15 K and 1 bar). Reactors add standard
    /// formation enthalpies to stream enthalpies on this basis.
    fn formation_components(&self) -> Result<&[Component], CalculationError> {
        Err(unsupported("reaction formation basis"))
    }

    /// Three-phase (VLLE/LLE) TP flash with multiphase stability (#224). Packages
    /// without the capability return `UnsupportedCapability`.
    fn vlle_tp(
        &self,
        _fractions: &[f64],
        _temperature_k: f64,
        _pressure_pa: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        Err(multiphase::unsupported("TP"))
    }

    /// Three-phase (VLLE/LLE) PH flash (#224). Packages without the capability
    /// return `UnsupportedCapability`.
    fn vlle_ph(
        &self,
        _fractions: &[f64],
        _pressure_pa: f64,
        _enthalpy: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        Err(multiphase::unsupported("PH"))
    }
}

fn unsupported(operation: &str) -> CalculationError {
    CalculationError::UnsupportedCapability(format!(
        "the selected property package does not support {operation}"
    ))
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

    pub(crate) fn phase(&self, p: ActivityFlashPhase) -> Result<EquilibriumPhase, CalculationError> {
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
                transport: None,
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

    fn vlle_tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        multiphase::nrtl_tp(self, fractions, temperature_k, pressure_pa)
    }

    fn vlle_ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        multiphase::nrtl_ph(self, fractions, pressure_pa, enthalpy)
    }

    fn cubic(&self) -> Result<&CubicPackage, CalculationError> {
        Err(CalculationError::UnsupportedCapability(
            "this unit's thermodynamic calculation requires a cubic package; selected NRTL-IG"
                .into(),
        ))
    }

    fn ps(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        entropy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let r = activity_ps_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            entropy,
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
                scaled_energy_residual: None,
                scaled_entropy_residual: Some(r.scaled_entropy_residual),
            },
            diagnostics: r.diagnostics,
        })
    }

    fn homogeneous(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<PhaseProperties, CalculationError> {
        let phase = self.homogeneous_phase(fractions, temperature_k, pressure_pa, root)?;
        Ok(self.phase(phase)?.state)
    }

    fn ln_fugacity_coefficients(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<Vec<f64>, CalculationError> {
        Ok(self
            .homogeneous_phase(fractions, temperature_k, pressure_pa, root)?
            .state
            .ln_fugacity_coefficients)
    }

    fn saturation_temperature_interval(&self) -> Result<[f64; 2], CalculationError> {
        Ok(self.temperature_interval())
    }

    fn saturation_pressure(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let [lo, hi] = self.temperature_interval();
        if !(lo..=hi).contains(&temperature_k) {
            return Err(CalculationError::OutsideValidityDomain(
                "temperature outside NRTL-IG package search interval".into(),
            ));
        }
        let model = self.model.with_fractions(fractions)?;
        let t = Temperature::new(temperature_k)?;
        let point = match kind {
            SaturationKind::Bubble => activity_bubble_pressure(&model, t)?,
            SaturationKind::Dew => activity_dew_pressure(&model, t, self.options)?,
        };
        self.saturation_state(point, kind)
    }

    fn saturation_temperature(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let model = self.model.with_fractions(fractions)?;
        let p = Pressure::new(pressure_pa)?;
        let point = match kind {
            SaturationKind::Bubble => activity_bubble_temperature(&model, p, self.options)?,
            SaturationKind::Dew => activity_dew_temperature(&model, p, self.options)?,
        };
        self.saturation_state(point, kind)
    }

    fn formation_components(&self) -> Result<&[Component], CalculationError> {
        // NRTL-IG shares the NASA7 ideal-gas caloric zero of the cubic packages.
        Ok(self.model.components())
    }
}

impl NrtlPackage {
    /// The gamma–phi phase on the selected root: the NRTL liquid for the smallest
    /// volume, the ideal-gas vapor for the largest, and the one of lower molar
    /// Gibbs energy for minimum Gibbs. No liquid–liquid stability claim.
    fn homogeneous_phase(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<ActivityFlashPhase, CalculationError> {
        let model = self.model.with_fractions(fractions)?;
        let (t, p) = (
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
        );
        let state = match root {
            RootSelection::SmallestVolume => model.liquid(t, p)?,
            RootSelection::LargestVolume => model.vapor(t, p)?,
            RootSelection::MinimumGibbs => {
                let (liquid, vapor) = (model.liquid(t, p)?, model.vapor(t, p)?);
                if liquid.gibbs_energy <= vapor.gibbs_energy {
                    liquid
                } else {
                    vapor
                }
            }
        };
        Ok(ActivityFlashPhase {
            phase: state.phase,
            phase_fraction: 1.0,
            composition: fractions.to_vec(),
            state,
        })
    }

    /// Phase states of a gamma–phi bubble/dew point. The liquid must pass the
    /// NRTL tangent-plane trials: a liquid–liquid split is unsupported (#224).
    fn saturation_state(
        &self,
        point: ActivitySaturation,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let (t, p) = (point.temperature_k, point.pressure_pa);
        let stability = liquid_stability(
            &self.model,
            Temperature::new(t)?,
            &point.liquid,
            self.options,
        )?;
        if !stability.stable_in_trials {
            return Err(CalculationError::UnsupportedCapability(
                "NRTL saturation liquid fails tangent-plane trials: liquid-liquid split (VLLE/LLE) is not supported".into(),
            ));
        }
        let liquid = self.homogeneous_phase(&point.liquid, t, p, RootSelection::SmallestVolume)?;
        let vapor = self.homogeneous_phase(&point.vapor, t, p, RootSelection::LargestVolume)?;
        Ok(SaturationState {
            kind,
            temperature_k: t,
            pressure_pa: p,
            liquid: self.phase(liquid)?.state,
            vapor: self.phase(vapor)?.state,
            liquid_composition: point.liquid,
            vapor_composition: point.vapor,
            max_ln_fugacity_residual: point.max_ln_fugacity_residual,
            iterations: point.iterations,
            diagnostics: vec![Diagnostic {
                code: "A7-GAMMA-PHI-SATURATION".into(),
                severity: DiagnosticSeverity::Information,
                message: "Modified Raoult bubble/dew point: NRTL liquid, ideal-gas vapor, no Poynting factor; liquid stable in finite NRTL tangent-plane trials.".into(),
            }],
        })
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

    fn vlle_tp(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        multiphase::cubic_tp(self, fractions, temperature_k, pressure_pa)
    }

    fn vlle_ph(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        enthalpy: f64,
    ) -> Result<MultiphaseEquilibrium, CalculationError> {
        multiphase::cubic_ph(self, fractions, pressure_pa, enthalpy)
    }

    fn cubic(&self) -> Result<&CubicPackage, CalculationError> {
        Ok(self)
    }

    fn ps(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        entropy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        self.ps(fractions, pressure_pa, entropy)
    }

    fn homogeneous(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<PhaseProperties, CalculationError> {
        CubicPackage::homogeneous(self, fractions, temperature_k, pressure_pa, root).map(Into::into)
    }

    fn ln_fugacity_coefficients(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<Vec<f64>, CalculationError> {
        self.validate_composition(fractions)?;
        Ok(self
            .model()
            .with_fractions(fractions)?
            .ln_fugacity(
                Temperature::new(temperature_k)?,
                Pressure::new(pressure_pa)?,
                root,
            )?
            .ln_fugacity_coefficients)
    }

    fn saturation_temperature_interval(&self) -> Result<[f64; 2], CalculationError> {
        let bounds = self.default_saturation_bounds()?;
        Ok([bounds.minimum_temperature_k, bounds.maximum_temperature_k])
    }

    fn saturation_pressure(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let point = CubicPackage::saturation_pressure(
            self,
            fractions,
            temperature_k,
            kind.into(),
            self.default_saturation_bounds()?,
        )?;
        Ok(SaturationState::from_point(point, "Cubic"))
    }

    fn saturation_temperature(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let point = CubicPackage::saturation_temperature(
            self,
            fractions,
            pressure_pa,
            kind.into(),
            self.default_saturation_bounds()?,
        )?;
        Ok(SaturationState::from_point(point, "Cubic"))
    }

    fn formation_components(&self) -> Result<&[Component], CalculationError> {
        Ok(self.model().components())
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

    fn ps(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        entropy: f64,
    ) -> Result<EquilibriumState, CalculationError> {
        let result = gerg_ps_flash(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            entropy,
            self.options,
        )?;
        Ok(EquilibriumState::from_ps(result, fractions))
    }

    fn homogeneous(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<PhaseProperties, CalculationError> {
        gerg_homogeneous(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            root,
        )
        .map(Into::into)
    }

    fn ln_fugacity_coefficients(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        pressure_pa: f64,
        root: RootSelection,
    ) -> Result<Vec<f64>, CalculationError> {
        Ok(gerg_ln_fugacity(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            Pressure::new(pressure_pa)?,
            root,
        )?
        .ln_fugacity_coefficients)
    }

    fn saturation_temperature_interval(&self) -> Result<[f64; 2], CalculationError> {
        Ok(self.temperature_interval())
    }

    fn saturation_pressure(
        &self,
        fractions: &[f64],
        temperature_k: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let point = gerg_saturation_pressure(
            &self.model.with_fractions(fractions)?,
            Temperature::new(temperature_k)?,
            kind.into(),
            self.saturation_bounds(),
        )?;
        Ok(SaturationState::from_point(point, "GERG-2008"))
    }

    fn saturation_temperature(
        &self,
        fractions: &[f64],
        pressure_pa: f64,
        kind: SaturationKind,
    ) -> Result<SaturationState, CalculationError> {
        let point = gerg_saturation_temperature(
            &self.model.with_fractions(fractions)?,
            Pressure::new(pressure_pa)?,
            kind.into(),
            self.saturation_bounds(),
        )?;
        Ok(SaturationState::from_point(point, "GERG-2008"))
    }

    fn formation_components(&self) -> Result<&[Component], CalculationError> {
        Err(CalculationError::UnsupportedCapability(
            "reactors need the ideal-gas formation basis (h = s = 0 at 298.15 K and 1 bar); GERG-2008 uses its native caloric reference and supplies no conversion offset".into(),
        ))
    }
}

impl GergPackage {
    /// A13 bubble/dew bounds: the package temperature interval and the GERG
    /// equilibrium pressure envelope.
    fn saturation_bounds(&self) -> SaturationOptions {
        SaturationOptions {
            minimum_temperature_k: self.options.minimum_temperature_k,
            maximum_temperature_k: self.options.maximum_temperature_k,
            minimum_pressure_pa: 1e5,
            maximum_pressure_pa: 1e7,
            ..SaturationOptions::default()
        }
    }
}
