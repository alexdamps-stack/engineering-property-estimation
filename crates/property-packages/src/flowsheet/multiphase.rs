//! Three-phase (VLLE/LLE) TP and PH flash of the flowsheet packages (#224).
//!
//! PR/PR-BM (cubic) and NRTL-IG implement it; GERG-2008 and EOS-CG-2021 keep the
//! `FlowsheetPackage` default and return `UnsupportedCapability`.

use thermo_solvers::{
    ActivityFlashPhase, MultiphaseOptions, MultiphasePhase, activity_vlle_ph_flash,
    activity_vlle_tp_flash, vlle_ph_flash, vlle_tp_flash,
};
use thermo_types::{CalculationError, Diagnostic, Pressure, Temperature};

pub use thermo_solvers::PhaseRole;

use super::NrtlPackage;
use crate::{Convergence, CubicPackage, EquilibriumPhase, EquilibriumState};

/// Multiphase equilibrium in package component order.
///
/// `state.phases` lists the vapor first, then the liquids by increasing mass
/// density; `roles` gives the role of each phase in the same order. Unit outlets
/// map phases by role: vapor to the vapor outlet, the light liquid to the light
/// liquid outlet and the heavy liquid to the heavy liquid outlet. A single liquid
/// has the role [`PhaseRole::Liquid`]; its outlet is the unit's decision.
#[derive(Debug, Clone)]
pub struct MultiphaseEquilibrium {
    pub state: EquilibriumState,
    pub roles: Vec<PhaseRole>,
    /// PH only: the target fell in the enthalpy step at an invariant temperature
    /// (pure component or binary three-phase temperature) and the phase amounts
    /// come from the lever rule.
    pub invariant_temperature: bool,
}

impl MultiphaseEquilibrium {
    /// The phase with `role`, if present.
    pub fn phase_with_role(&self, role: PhaseRole) -> Option<&EquilibriumPhase> {
        self.roles
            .iter()
            .position(|r| *r == role)
            .map(|k| &self.state.phases[k])
    }
}

pub(super) fn unsupported(operation: &str) -> CalculationError {
    CalculationError::UnsupportedCapability(format!(
        "three-phase (VLLE/LLE) {operation} flash is not supported by the selected property package"
    ))
}

fn options([lo, hi]: [f64; 2]) -> MultiphaseOptions {
    MultiphaseOptions {
        minimum_temperature_k: lo,
        maximum_temperature_k: hi,
        ..Default::default()
    }
}

fn check_temperature(interval: [f64; 2], t: f64, name: &str) -> Result<(), CalculationError> {
    if !(interval[0]..=interval[1]).contains(&t) {
        return Err(CalculationError::OutsideValidityDomain(format!(
            "temperature outside {name} package search interval"
        )));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn assemble<S>(
    fractions: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
    phases: Vec<MultiphasePhase<S>>,
    convert: impl Fn(MultiphasePhase<S>) -> Result<EquilibriumPhase, CalculationError>,
    iterations: usize,
    convergence: Convergence,
    diagnostics: Vec<Diagnostic>,
    invariant_temperature: bool,
) -> Result<MultiphaseEquilibrium, CalculationError> {
    let roles = phases.iter().map(|p| p.role).collect();
    Ok(MultiphaseEquilibrium {
        state: EquilibriumState {
            temperature_k,
            pressure_pa,
            mole_fractions: fractions.to_vec(),
            phases: phases.into_iter().map(convert).collect::<Result<_, _>>()?,
            iterations,
            convergence,
            diagnostics,
        },
        roles,
        invariant_temperature,
    })
}

fn cubic_phase(
    p: MultiphasePhase<thermo_models::HomogeneousState>,
) -> Result<EquilibriumPhase, CalculationError> {
    Ok(EquilibriumPhase {
        phase: p.phase,
        phase_fraction: p.phase_fraction,
        composition: p.composition,
        state: p.state.into(),
    })
}

pub(super) fn cubic_tp(
    package: &CubicPackage,
    fractions: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<MultiphaseEquilibrium, CalculationError> {
    package.validate_options()?;
    check_temperature(package.temperature_interval(), temperature_k, "cubic")?;
    let r = vlle_tp_flash(
        &package.model.with_fractions(fractions)?,
        Temperature::new(temperature_k)?,
        Pressure::new(pressure_pa)?,
        options(package.temperature_interval()),
    )?;
    let convergence = Convergence {
        max_ln_fugacity_residual: r.max_ln_fugacity_residual,
        max_component_balance_residual: Some(r.max_component_balance_residual),
        scaled_energy_residual: None,
        scaled_entropy_residual: None,
    };
    assemble(
        fractions,
        r.temperature_k,
        r.pressure_pa,
        r.phases,
        cubic_phase,
        r.iterations,
        convergence,
        r.diagnostics,
        false,
    )
}

pub(super) fn cubic_ph(
    package: &CubicPackage,
    fractions: &[f64],
    pressure_pa: f64,
    enthalpy: f64,
) -> Result<MultiphaseEquilibrium, CalculationError> {
    package.validate_options()?;
    let r = vlle_ph_flash(
        &package.model.with_fractions(fractions)?,
        Pressure::new(pressure_pa)?,
        enthalpy,
        options(package.temperature_interval()),
    )?;
    let convergence = Convergence {
        max_ln_fugacity_residual: r.max_ln_fugacity_residual,
        max_component_balance_residual: Some(r.max_component_balance_residual),
        scaled_energy_residual: Some(r.scaled_energy_residual),
        scaled_entropy_residual: None,
    };
    assemble(
        fractions,
        r.temperature_k,
        r.pressure_pa,
        r.phases,
        cubic_phase,
        r.iterations + r.tp.iterations,
        convergence,
        r.diagnostics,
        r.invariant_temperature,
    )
}

fn nrtl_phase(
    package: &NrtlPackage,
    p: MultiphasePhase<thermo_models::ActivityPhaseState>,
) -> Result<EquilibriumPhase, CalculationError> {
    package.phase(ActivityFlashPhase {
        phase: p.phase,
        phase_fraction: p.phase_fraction,
        composition: p.composition,
        state: p.state,
    })
}

pub(super) fn nrtl_tp(
    package: &NrtlPackage,
    fractions: &[f64],
    temperature_k: f64,
    pressure_pa: f64,
) -> Result<MultiphaseEquilibrium, CalculationError> {
    let interval = [
        package.options.minimum_temperature_k,
        package.options.maximum_temperature_k,
    ];
    check_temperature(interval, temperature_k, "NRTL-IG")?;
    let r = activity_vlle_tp_flash(
        &package.model.with_fractions(fractions)?,
        Temperature::new(temperature_k)?,
        Pressure::new(pressure_pa)?,
        options(interval),
    )?;
    let convergence = Convergence {
        max_ln_fugacity_residual: r.max_ln_fugacity_residual,
        max_component_balance_residual: Some(r.max_component_balance_residual),
        scaled_energy_residual: None,
        scaled_entropy_residual: None,
    };
    assemble(
        fractions,
        r.temperature_k,
        r.pressure_pa,
        r.phases,
        |p| nrtl_phase(package, p),
        r.iterations,
        convergence,
        r.diagnostics,
        false,
    )
}

pub(super) fn nrtl_ph(
    package: &NrtlPackage,
    fractions: &[f64],
    pressure_pa: f64,
    enthalpy: f64,
) -> Result<MultiphaseEquilibrium, CalculationError> {
    let interval = [
        package.options.minimum_temperature_k,
        package.options.maximum_temperature_k,
    ];
    let r = activity_vlle_ph_flash(
        &package.model.with_fractions(fractions)?,
        Pressure::new(pressure_pa)?,
        enthalpy,
        options(interval),
    )?;
    let convergence = Convergence {
        max_ln_fugacity_residual: r.max_ln_fugacity_residual,
        max_component_balance_residual: Some(r.max_component_balance_residual),
        scaled_energy_residual: Some(r.scaled_energy_residual),
        scaled_entropy_residual: None,
    };
    assemble(
        fractions,
        r.temperature_k,
        r.pressure_pa,
        r.phases,
        |p| nrtl_phase(package, p),
        r.iterations + r.tp.iterations,
        convergence,
        r.diagnostics,
        r.invariant_temperature,
    )
}
