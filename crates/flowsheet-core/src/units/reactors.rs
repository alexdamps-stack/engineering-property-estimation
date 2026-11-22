//! Reactors group (#191–#197): conversion, yield, equilibrium, Gibbs, CSTR,
//! plug-flow and batch reactors.
//!
//! Owned by the reactors work package (`docs/a17-architecture.md`). The
//! reaction-set data model and formation properties belong to this group;
//! packed (catalytic) beds and the thermal modes are described in
//! `docs/a17-implementation-brief.md`.

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::unit::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, UnitKind};

/// Metadata of a reactor kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    let reactor =
        |label, issue| kind_info(label, "R", PaletteCategory::Reactors, &SINGLE, &DUTY, issue);
    match kind {
        UnitKind::ConversionReactor => reactor("Conversion reactor", 191),
        UnitKind::YieldReactor => reactor("Yield reactor", 192),
        UnitKind::EquilibriumReactor => reactor("Equilibrium reactor", 193),
        UnitKind::GibbsReactor => reactor("Gibbs reactor", 194),
        UnitKind::Cstr => reactor("CSTR", 195),
        UnitKind::PlugFlowReactor => reactor("Plug-flow reactor", 196),
        UnitKind::BatchReactor => reactor("Batch reactor", 197),
        _ => unreachable!("{kind:?} is not a reactor kind"),
    }
}

/// Calculation status of a reactor kind.
pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    DRAWN_ONLY
}

/// Specifications of the reactor kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactorSpec {
    /// Conversion reactor (#191).
    Conversion(ConversionReactorSpec),
    /// Yield reactor (#192).
    Yield(YieldReactorSpec),
    /// Equilibrium reactor (#193).
    Equilibrium(EquilibriumReactorSpec),
    /// Gibbs reactor (#194).
    Gibbs(GibbsReactorSpec),
    /// CSTR (#195).
    Cstr(CstrSpec),
    /// Plug-flow reactor (#196).
    PlugFlow(PlugFlowReactorSpec),
    /// Batch reactor (#197).
    Batch(BatchReactorSpec),
}

/// Conversion reactor specification (#191).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConversionReactorSpec {}

/// Yield reactor specification (#192).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct YieldReactorSpec {}

/// Equilibrium reactor specification (#193).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EquilibriumReactorSpec {}

/// Gibbs reactor specification (#194).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GibbsReactorSpec {}

/// CSTR specification (#195).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CstrSpec {}

/// Plug-flow reactor specification (#196).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PlugFlowReactorSpec {}

/// Batch reactor specification (#197).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BatchReactorSpec {}

/// Default specification of a reactor kind.
pub(crate) fn default_spec(kind: UnitKind) -> ReactorSpec {
    match kind {
        UnitKind::ConversionReactor => ReactorSpec::Conversion(Default::default()),
        UnitKind::YieldReactor => ReactorSpec::Yield(Default::default()),
        UnitKind::EquilibriumReactor => ReactorSpec::Equilibrium(Default::default()),
        UnitKind::GibbsReactor => ReactorSpec::Gibbs(Default::default()),
        UnitKind::Cstr => ReactorSpec::Cstr(Default::default()),
        UnitKind::PlugFlowReactor => ReactorSpec::PlugFlow(Default::default()),
        UnitKind::BatchReactor => ReactorSpec::Batch(Default::default()),
        _ => unreachable!("{kind:?} has no ReactorSpec"),
    }
}

impl GroupSpec for ReactorSpec {
    fn kind(&self) -> UnitKind {
        match self {
            ReactorSpec::Conversion(_) => UnitKind::ConversionReactor,
            ReactorSpec::Yield(_) => UnitKind::YieldReactor,
            ReactorSpec::Equilibrium(_) => UnitKind::EquilibriumReactor,
            ReactorSpec::Gibbs(_) => UnitKind::GibbsReactor,
            ReactorSpec::Cstr(_) => UnitKind::Cstr,
            ReactorSpec::PlugFlow(_) => UnitKind::PlugFlowReactor,
            ReactorSpec::Batch(_) => UnitKind::BatchReactor,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        Vec::new()
    }
}

/// Calculates a reactor unit.
pub(crate) fn calculate(
    kind: UnitKind,
    _spec: &ReactorSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(super::not_implemented(kind))
}

#[cfg(test)]
mod tests {
    use crate::UnitGroup;

    #[test]
    fn reactor_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Reactors);
    }
}
