//! Model adapters bound to one temperature and pressure.

use thermo_data::Component;
use thermo_models::{
    ActivityModel, ActivityPhaseState, CubicModel, HomogeneousState, RootSelection,
};
use thermo_types::{CalculationError as Error, Phase, Pressure, Temperature};

use super::FluidKind;

/// Equilibrium model at fixed T and p as seen by the multiphase solver.
pub(crate) trait VlleAdapter {
    type State: Clone;
    fn temperature(&self) -> Temperature;
    fn pressure(&self) -> Pressure;
    /// Kinds compared for the minimum-Gibbs reference phase of a composition.
    fn kinds(&self) -> &'static [FluidKind];
    fn vapor_kind(&self) -> FluidKind;
    fn liquid_kind(&self) -> FluidKind;
    /// ln φ_i of composition `x` (normalised, model order) in phase `kind`.
    fn ln_phi(&self, x: &[f64], kind: FluidKind) -> Result<Vec<f64>, Error>;
    /// Wilson-type K values used for vapor/liquid trial starts.
    fn wilson_k(&self) -> Vec<f64>;
    fn components(&self) -> &[Component];
    fn state(&self, x: &[f64], kind: FluidKind) -> Result<Self::State, Error>;
    fn molar_volume(state: &Self::State) -> f64;
    fn mass_density(state: &Self::State) -> f64;
    fn enthalpy(state: &Self::State) -> f64;
    /// ln φ_i stored with a phase state (same definition as `ln_phi`).
    fn state_ln_phi(state: &Self::State) -> &[f64];
    /// Single-phase label of a phase.
    fn label(&self, x: &[f64], kind: FluidKind, state: &Self::State) -> Result<Phase, Error>;
    /// ln φ on the root matching a phase label (cubic: largest volume for vapor,
    /// smallest for liquid), for coexisting roots at an invariant temperature.
    fn ln_phi_labelled(&self, x: &[f64], kind: FluidKind, phase: Phase) -> Result<Vec<f64>, Error>;
    /// Full state on the root matching a phase label.
    fn state_labelled(
        &self,
        x: &[f64],
        kind: FluidKind,
        phase: Phase,
    ) -> Result<Self::State, Error>;
}

/// Index of water (by component identifier) among the components.
pub(crate) fn water_index(components: &[Component]) -> Option<usize> {
    components.iter().position(|c| c.id() == "water")
}

/// True for a component whose formula contains only carbon and hydrogen.
pub(crate) fn is_hydrocarbon(component: &Component) -> bool {
    let elements = component.elements();
    elements.iter().any(|(e, _)| e == "C") && elements.iter().all(|(e, _)| e == "C" || e == "H")
}

/// Cubic root of a labelled phase.
fn root_for(phase: Phase) -> RootSelection {
    if phase == Phase::Liquid {
        RootSelection::SmallestVolume
    } else {
        RootSelection::LargestVolume
    }
}

pub(crate) struct CubicAt<'a> {
    model: &'a CubicModel,
    t: Temperature,
    p: Pressure,
}

impl<'a> CubicAt<'a> {
    pub(crate) fn new(model: &'a CubicModel, t: Temperature, p: Pressure) -> Self {
        Self { model, t, p }
    }
}

impl VlleAdapter for CubicAt<'_> {
    type State = HomogeneousState;
    fn temperature(&self) -> Temperature {
        self.t
    }
    fn pressure(&self) -> Pressure {
        self.p
    }
    fn kinds(&self) -> &'static [FluidKind] {
        &[FluidKind::Fluid]
    }
    fn vapor_kind(&self) -> FluidKind {
        FluidKind::Fluid
    }
    fn liquid_kind(&self) -> FluidKind {
        FluidKind::Fluid
    }
    fn ln_phi(&self, x: &[f64], _kind: FluidKind) -> Result<Vec<f64>, Error> {
        Ok(self
            .model
            .with_fractions(x)?
            .ln_fugacity(self.t, self.p, RootSelection::MinimumGibbs)?
            .ln_fugacity_coefficients)
    }
    fn wilson_k(&self) -> Vec<f64> {
        crate::wilson_k(self.model, self.t, self.p)
    }
    fn components(&self) -> &[Component] {
        self.model.components()
    }
    fn state(&self, x: &[f64], _kind: FluidKind) -> Result<HomogeneousState, Error> {
        self.model
            .with_fractions(x)?
            .evaluate(self.t, self.p, RootSelection::MinimumGibbs)
    }
    fn molar_volume(s: &HomogeneousState) -> f64 {
        s.molar_volume
    }
    fn mass_density(s: &HomogeneousState) -> f64 {
        s.mass_density
    }
    fn enthalpy(s: &HomogeneousState) -> f64 {
        s.enthalpy
    }
    fn state_ln_phi(s: &HomogeneousState) -> &[f64] {
        &s.ln_fugacity_coefficients
    }
    fn label(&self, x: &[f64], _kind: FluidKind, s: &HomogeneousState) -> Result<Phase, Error> {
        Ok(crate::single_phase_label(&self.model.with_fractions(x)?, s))
    }
    fn ln_phi_labelled(
        &self,
        x: &[f64],
        _kind: FluidKind,
        phase: Phase,
    ) -> Result<Vec<f64>, Error> {
        Ok(self
            .model
            .with_fractions(x)?
            .ln_fugacity(self.t, self.p, root_for(phase))?
            .ln_fugacity_coefficients)
    }
    fn state_labelled(
        &self,
        x: &[f64],
        _kind: FluidKind,
        phase: Phase,
    ) -> Result<HomogeneousState, Error> {
        self.model
            .with_fractions(x)?
            .evaluate(self.t, self.p, root_for(phase))
    }
}

pub(crate) struct ActivityAt<'a> {
    model: &'a ActivityModel,
    t: Temperature,
    p: Pressure,
    /// ln(p_i^sat / p) per component.
    ln_psat_over_p: Vec<f64>,
}

impl<'a> ActivityAt<'a> {
    pub(crate) fn new(
        model: &'a ActivityModel,
        t: Temperature,
        p: Pressure,
    ) -> Result<Self, Error> {
        model.check_domain(t, p)?;
        let ln_psat_over_p = model
            .pure_liquids(t)?
            .iter()
            .map(|s| s.vapor_pressure.ln() - p.value().ln())
            .collect();
        Ok(Self {
            model,
            t,
            p,
            ln_psat_over_p,
        })
    }
}

impl VlleAdapter for ActivityAt<'_> {
    type State = ActivityPhaseState;
    fn temperature(&self) -> Temperature {
        self.t
    }
    fn pressure(&self) -> Pressure {
        self.p
    }
    fn kinds(&self) -> &'static [FluidKind] {
        &[FluidKind::Liquid, FluidKind::Vapor]
    }
    fn vapor_kind(&self) -> FluidKind {
        FluidKind::Vapor
    }
    fn liquid_kind(&self) -> FluidKind {
        FluidKind::Liquid
    }
    fn ln_phi(&self, x: &[f64], kind: FluidKind) -> Result<Vec<f64>, Error> {
        match kind {
            FluidKind::Vapor => Ok(vec![0.0; x.len()]),
            _ => {
                let ex = self.model.nrtl().excess_at(self.t, x)?;
                Ok(ex
                    .ln_activity_coefficients
                    .iter()
                    .zip(&self.ln_psat_over_p)
                    .map(|(g, s)| g + s)
                    .collect())
            }
        }
    }
    fn wilson_k(&self) -> Vec<f64> {
        self.ln_psat_over_p.iter().map(|l| l.exp()).collect()
    }
    fn components(&self) -> &[Component] {
        self.model.components()
    }
    fn state(&self, x: &[f64], kind: FluidKind) -> Result<ActivityPhaseState, Error> {
        let m = self.model.with_fractions(x)?;
        match kind {
            FluidKind::Vapor => m.vapor(self.t, self.p),
            _ => m.liquid(self.t, self.p),
        }
    }
    fn molar_volume(s: &ActivityPhaseState) -> f64 {
        s.molar_volume
    }
    fn mass_density(s: &ActivityPhaseState) -> f64 {
        s.mass_density
    }
    fn enthalpy(s: &ActivityPhaseState) -> f64 {
        s.enthalpy
    }
    fn state_ln_phi(s: &ActivityPhaseState) -> &[f64] {
        &s.ln_fugacity_coefficients
    }
    fn ln_phi_labelled(
        &self,
        x: &[f64],
        kind: FluidKind,
        _phase: Phase,
    ) -> Result<Vec<f64>, Error> {
        self.ln_phi(x, kind)
    }
    fn state_labelled(
        &self,
        x: &[f64],
        kind: FluidKind,
        _phase: Phase,
    ) -> Result<ActivityPhaseState, Error> {
        self.state(x, kind)
    }
    fn label(&self, _x: &[f64], kind: FluidKind, _s: &ActivityPhaseState) -> Result<Phase, Error> {
        Ok(match kind {
            FluidKind::Vapor => Phase::Vapor,
            _ => Phase::Liquid,
        })
    }
}
