//! Registry metadata for the reactor calculations (#191–#197).
//!
//! The reactor specifications retain the flat variants used by existing
//! projects. Their calculations share the engine's reaction balance.

use super::kind_info;
use crate::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, UnitKind};

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

pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    CalculationStatus::Calculated
}
