//! Columns group: shortcut distillation (#185), rigorous distillation
//! (#186), absorber / stripper (#187) and petroleum fractionator (#190),
//! with equilibrium-stage and packed / rate-based sections
//! (`docs/a17-implementation-brief.md`).
//!
//! Owned by the columns work package (`docs/a17-architecture.md`). Column
//! feeds and side draws are variable, numbered ports.

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info, numbered};
use crate::unit::{
    CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, SpecProblem, UnitKind,
};

/// Condenser duty leaves the column; reboiler duty enters it.
const CONDENSER_REBOILER: [PortSpec; 2] = [
    PortSpec::energy_outlet("condenser duty"),
    PortSpec::energy_inlet("reboiler duty"),
];

/// Metadata of a column kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::Columns as C;
    const SHORTCUT: [PortSpec; 3] = [
        PortSpec::inlet("feed"),
        PortSpec::outlet("distillate"),
        PortSpec::outlet("bottoms"),
    ];
    // Feeds and side draws are variable ports (see `DistillationColumnSpec`).
    const PRODUCTS: [PortSpec; 2] = [PortSpec::outlet("distillate"), PortSpec::outlet("bottoms")];
    const ABSORBER: [PortSpec; 4] = [
        PortSpec::inlet("gas in"),
        PortSpec::inlet("liquid in"),
        PortSpec::outlet("gas out"),
        PortSpec::outlet("liquid out"),
    ];
    const FRACTIONATOR: [PortSpec; 3] = [
        PortSpec::inlet("feed"),
        PortSpec::outlet("overhead"),
        PortSpec::outlet("bottoms"),
    ];
    const CONDENSER: [PortSpec; 1] = [PortSpec::energy_outlet("condenser duty")];
    match kind {
        UnitKind::ShortcutColumn => kind_info(
            "Shortcut column",
            "T",
            C,
            &SHORTCUT,
            &CONDENSER_REBOILER,
            185,
        ),
        UnitKind::DistillationColumn => kind_info(
            "Distillation column",
            "T",
            C,
            &PRODUCTS,
            &CONDENSER_REBOILER,
            186,
        ),
        UnitKind::Absorber => kind_info("Absorber", "T", C, &ABSORBER, &[], 187),
        UnitKind::PetroleumFractionator => kind_info(
            "Petroleum fractionator",
            "T",
            C,
            &FRACTIONATOR,
            &CONDENSER,
            190,
        ),
        _ => unreachable!("{kind:?} is not a column kind"),
    }
}

/// Calculation status of a column kind.
pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    DRAWN_ONLY
}

/// Specifications of the column kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum ColumnSpec {
    /// Shortcut column (#185).
    Shortcut(ShortcutColumnSpec),
    /// Rigorous distillation column (#186).
    Distillation(DistillationColumnSpec),
    /// Absorber / stripper (#187).
    Absorber(AbsorberSpec),
    /// Petroleum fractionator (#190).
    PetroleumFractionator(PetroleumFractionatorSpec),
}

/// Shortcut column specification (#185).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ShortcutColumnSpec {}

/// Largest number of feeds or side draws of a column.
pub const MAX_COLUMN_PORTS: u16 = 20;

/// Rigorous distillation column specification (#186). Feed `k` is the
/// variable port `feed k`, side draw `k` the port `side draw k`.
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnSpec {
    /// Number of feeds, 1 to [`MAX_COLUMN_PORTS`].
    pub feeds: u16,
    /// Number of side draws, 0 to [`MAX_COLUMN_PORTS`].
    pub side_draws: u16,
}

impl Default for DistillationColumnSpec {
    fn default() -> Self {
        Self {
            feeds: 1,
            side_draws: 0,
        }
    }
}

/// Absorber / stripper specification (#187).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AbsorberSpec {}

/// Petroleum fractionator specification (#190). Side draw `k` is the
/// variable port `side draw k`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PetroleumFractionatorSpec {
    /// Number of side draws, 0 to [`MAX_COLUMN_PORTS`].
    pub side_draws: u16,
}

/// Default specification of a column kind.
pub(crate) fn default_spec(kind: UnitKind) -> ColumnSpec {
    match kind {
        UnitKind::ShortcutColumn => ColumnSpec::Shortcut(Default::default()),
        UnitKind::DistillationColumn => ColumnSpec::Distillation(Default::default()),
        UnitKind::Absorber => ColumnSpec::Absorber(Default::default()),
        UnitKind::PetroleumFractionator => ColumnSpec::PetroleumFractionator(Default::default()),
        _ => unreachable!("{kind:?} has no ColumnSpec"),
    }
}

fn count_issue(out: &mut Vec<SpecIssue>, field: &'static str, n: u16, min: u16) {
    if !(min..=MAX_COLUMN_PORTS).contains(&n) {
        out.push(SpecIssue {
            field,
            problem: SpecProblem::CountOutOfRange,
        });
    }
}

impl GroupSpec for ColumnSpec {
    fn kind(&self) -> UnitKind {
        match self {
            ColumnSpec::Shortcut(_) => UnitKind::ShortcutColumn,
            ColumnSpec::Distillation(_) => UnitKind::DistillationColumn,
            ColumnSpec::Absorber(_) => UnitKind::Absorber,
            ColumnSpec::PetroleumFractionator(_) => UnitKind::PetroleumFractionator,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        let mut out = Vec::new();
        match self {
            ColumnSpec::Distillation(s) => {
                count_issue(&mut out, "feeds", s.feeds, 1);
                count_issue(&mut out, "side draws", s.side_draws, 0);
            }
            ColumnSpec::PetroleumFractionator(s) => {
                count_issue(&mut out, "side draws", s.side_draws, 0)
            }
            _ => {}
        }
        out
    }

    fn variable_ports(&self) -> Vec<PortSpec> {
        let side_draws = |n| numbered(PortSpec::outlet("side draw"), n);
        match self {
            ColumnSpec::Distillation(s) => numbered(PortSpec::inlet("feed"), s.feeds)
                .chain(side_draws(s.side_draws))
                .collect(),
            ColumnSpec::PetroleumFractionator(s) => side_draws(s.side_draws).collect(),
            _ => Vec::new(),
        }
    }
}

/// Calculates a column unit.
pub(crate) fn calculate(
    kind: UnitKind,
    _spec: &ColumnSpec,
    _inputs: &UnitInputs,
    _ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    Err(super::not_implemented(kind))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UnitGroup;

    #[test]
    fn column_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::Columns);
    }

    #[test]
    fn feeds_and_side_draws_are_numbered_ports() {
        let spec = ColumnSpec::Distillation(DistillationColumnSpec {
            feeds: 2,
            side_draws: 1,
        });
        let labels: Vec<_> = spec.variable_ports().iter().map(|p| p.label()).collect();
        assert_eq!(labels, ["feed 1", "feed 2", "side draw 1"]);
        let none = ColumnSpec::Distillation(DistillationColumnSpec {
            feeds: 0,
            side_draws: 0,
        });
        assert_eq!(none.issues()[0].field, "feeds");
    }
}
