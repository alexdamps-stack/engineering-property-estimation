//! Multiphase group: VLLE flash with three-phase separator (#176), decanter
//! (#177), three-phase distillation (#188) and liquid–liquid extractor
//! (#189).
//!
//! Owned by the multiphase work package (`docs/a17-architecture.md`). The
//! VLLE / LLE flash is a prerequisite this group adds (the package exposes
//! TP/PH two-phase flashes and fugacity coefficients through
//! [`CalcContext`]).

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::unit::{CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, UnitKind};

/// Metadata of a multiphase kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::{Columns as C, Separators as S};
    const THREE_PHASE: [PortSpec; 4] = [
        PortSpec::inlet("in"),
        PortSpec::outlet("vapor"),
        PortSpec::outlet("light liquid"),
        PortSpec::outlet("heavy liquid"),
    ];
    const DECANTER: [PortSpec; 3] = [
        PortSpec::inlet("in"),
        PortSpec::outlet("light liquid"),
        PortSpec::outlet("heavy liquid"),
    ];
    const COLUMN: [PortSpec; 4] = [
        PortSpec::inlet("feed"),
        PortSpec::outlet("distillate"),
        PortSpec::outlet("decant"),
        PortSpec::outlet("bottoms"),
    ];
    const EXTRACTOR: [PortSpec; 4] = [
        PortSpec::inlet("feed"),
        PortSpec::inlet("solvent"),
        PortSpec::outlet("extract"),
        PortSpec::outlet("raffinate"),
    ];
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    const CONDENSER_REBOILER: [PortSpec; 2] = [
        PortSpec::energy_outlet("condenser duty"),
        PortSpec::energy_inlet("reboiler duty"),
    ];
    match kind {
        UnitKind::ThreePhaseSeparator => {
            kind_info("Three-phase separator", "D", S, &THREE_PHASE, &DUTY, 176)
        }
        UnitKind::Decanter => kind_info("Decanter", "D", S, &DECANTER, &DUTY, 177),
        UnitKind::ThreePhaseColumn => kind_info(
            "Three-phase column",
            "T",
            C,
            &COLUMN,
            &CONDENSER_REBOILER,
            188,
        ),
        UnitKind::LiquidExtractor => {
            kind_info("Liquid-liquid extractor", "T", C, &EXTRACTOR, &[], 189)
        }
        _ => unreachable!("{kind:?} is not a multiphase kind"),
    }
}

/// Calculation status of a multiphase kind.
pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    DRAWN_ONLY
}

/// Specifications of the multiphase kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum MultiphaseSpec {
    /// Three-phase separator (#176).
    ThreePhaseSeparator(ThreePhaseSeparatorSpec),
    /// Decanter (#177).
    Decanter(DecanterSpec),
    /// Three-phase distillation (#188).
    ThreePhaseColumn(ThreePhaseColumnSpec),
    /// Liquid–liquid extractor (#189).
    LiquidExtractor(LiquidExtractorSpec),
}

/// Three-phase separator specification (#176).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThreePhaseSeparatorSpec {}

/// Decanter specification (#177).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DecanterSpec {}

/// Three-phase distillation specification (#188).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThreePhaseColumnSpec {}

/// Liquid–liquid extractor specification (#189).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct LiquidExtractorSpec {}

/// Default specification of a multiphase kind.
pub(crate) fn default_spec(kind: UnitKind) -> MultiphaseSpec {
    match kind {
        UnitKind::ThreePhaseSeparator => MultiphaseSpec::ThreePhaseSeparator(Default::default()),
        UnitKind::Decanter => MultiphaseSpec::Decanter(Default::default()),
        UnitKind::ThreePhaseColumn => MultiphaseSpec::ThreePhaseColumn(Default::default()),
        UnitKind::LiquidExtractor => MultiphaseSpec::LiquidExtractor(Default::default()),
        _ => unreachable!("{kind:?} has no MultiphaseSpec"),
    }
}

impl GroupSpec for MultiphaseSpec {
    fn kind(&self) -> UnitKind {
        match self {
            MultiphaseSpec::ThreePhaseSeparator(_) => UnitKind::ThreePhaseSeparator,
            MultiphaseSpec::Decanter(_) => UnitKind::Decanter,
            MultiphaseSpec::ThreePhaseColumn(_) => UnitKind::ThreePhaseColumn,
            MultiphaseSpec::LiquidExtractor(_) => UnitKind::LiquidExtractor,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        Vec::new()
    }
}

/// Calculates a multiphase unit.
pub(crate) fn calculate(
    kind: UnitKind,
    _spec: &MultiphaseSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(super::not_implemented(kind))
}

#[cfg(test)]
mod tests {
    use crate::UnitGroup;

    #[test]
    fn multiphase_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Multiphase);
    }
}
