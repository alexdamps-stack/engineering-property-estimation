//! Simple and logical group: component splitter (#172), multiplier /
//! duplicator (#173), selector (#174), stream cutter (#175), tank (#178),
//! adjust (#212), set (#213), spreadsheet / calculator (#214), balance
//! (#215), recycle block (#216), control blocks (#217) and user / CAPE-OPEN
//! unit (#218).
//!
//! Owned by the simple/logical work package (`docs/a17-architecture.md`).
//! Logical blocks (adjust, set, spreadsheet, control) have no material
//! ports; the variable links they need are part of this group's design.

mod component_splitter;

pub use component_splitter::ComponentSplitterSpec;

use super::{CalcContext, DRAWN_ONLY, GroupSpec, UnitInputs, UnitOutput, kind_info, numbered};
use crate::unit::{
    CalculationStatus, KindInfo, PaletteCategory, PortSpec, SpecIssue, SpecProblem, UnitKind,
};

const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];

/// Metadata of a simple/logical kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::{
        Extensibility as X, Logical as L, MixingSplitting as M, Separators as S,
    };
    const SPLIT: [PortSpec; 3] = [
        PortSpec::inlet("in"),
        PortSpec::outlet("overhead"),
        PortSpec::outlet("bottoms"),
    ];
    const DUTY: [PortSpec; 1] = [PortSpec::energy_inlet("duty")];
    // Selector inlets are variable ports (see `SelectorSpec`).
    const SELECTED: [PortSpec; 1] = [PortSpec::outlet("out")];
    match kind {
        UnitKind::ComponentSplitter => kind_info("Component splitter", "X", M, &SPLIT, &DUTY, 172),
        UnitKind::Multiplier => kind_info("Stream multiplier", "MULT", M, &SINGLE, &[], 173),
        UnitKind::Selector => kind_info("Stream selector", "SEL", M, &SELECTED, &[], 174),
        UnitKind::StreamCutter => kind_info("Stream cutter", "CUT", M, &SINGLE, &[], 175),
        UnitKind::Tank => kind_info("Tank", "TK", S, &SINGLE, &DUTY, 178),
        UnitKind::Adjust => kind_info("Adjust", "ADJ", L, &[], &[], 212),
        UnitKind::Set => kind_info("Set", "SET", L, &[], &[], 213),
        UnitKind::Spreadsheet => kind_info("Spreadsheet", "SPRDSHT", L, &[], &[], 214),
        UnitKind::Balance => kind_info("Balance", "BAL", L, &SINGLE, &[], 215),
        UnitKind::RecycleBlock => kind_info("Recycle", "RCY", L, &SINGLE, &[], 216),
        UnitKind::ControlBlock => kind_info("Control block", "CTL", L, &[], &[], 217),
        UnitKind::UserUnit => kind_info("User unit", "UOP", X, &SINGLE, &DUTY, 218),
        _ => unreachable!("{kind:?} is not a simple/logical kind"),
    }
}

/// Calculation status of a simple/logical kind.
pub(crate) fn status(kind: UnitKind) -> CalculationStatus {
    match kind {
        UnitKind::ComponentSplitter => CalculationStatus::Calculated,
        _ => DRAWN_ONLY,
    }
}

/// Specifications of the simple/logical kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum SimpleLogicalSpec {
    /// Component splitter (#172).
    ComponentSplitter(ComponentSplitterSpec),
    /// Stream multiplier / duplicator (#173).
    Multiplier(MultiplierSpec),
    /// Stream selector (#174).
    Selector(SelectorSpec),
    /// Stream cutter (#175).
    StreamCutter(StreamCutterSpec),
    /// Tank (#178).
    Tank(TankSpec),
    /// Adjust (#212).
    Adjust(AdjustSpec),
    /// Set (#213).
    Set(SetSpec),
    /// Spreadsheet / calculator (#214).
    Spreadsheet(SpreadsheetSpec),
    /// Balance (#215).
    Balance(BalanceSpec),
    /// Recycle block (#216).
    RecycleBlock(RecycleBlockSpec),
    /// Control block (#217).
    ControlBlock(ControlBlockSpec),
    /// User / CAPE-OPEN unit (#218).
    UserUnit(UserUnitSpec),
}

/// Stream multiplier / duplicator specification (#173).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MultiplierSpec {}

/// Largest number of selector inlets.
pub const MAX_SELECTOR_INLETS: u16 = 12;

/// Stream selector specification (#174). Inlet `k` is the variable port
/// `in k`.
#[derive(Debug, Clone, PartialEq)]
pub struct SelectorSpec {
    /// Number of inlets, 2 to [`MAX_SELECTOR_INLETS`].
    pub inlets: u16,
}

impl Default for SelectorSpec {
    fn default() -> Self {
        Self { inlets: 2 }
    }
}

/// Stream cutter specification (#175).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StreamCutterSpec {}

/// Tank specification (#178).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TankSpec {}

/// Adjust specification (#212).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AdjustSpec {}

/// Set specification (#213).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SetSpec {}

/// Spreadsheet specification (#214).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SpreadsheetSpec {}

/// Balance specification (#215).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BalanceSpec {}

/// Recycle block specification (#216).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RecycleBlockSpec {}

/// Control block specification (#217).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ControlBlockSpec {}

/// User / CAPE-OPEN unit specification (#218).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UserUnitSpec {}

/// Default specification of a simple/logical kind.
pub(crate) fn default_spec(kind: UnitKind) -> SimpleLogicalSpec {
    use SimpleLogicalSpec as S;
    match kind {
        UnitKind::ComponentSplitter => S::ComponentSplitter(Default::default()),
        UnitKind::Multiplier => S::Multiplier(Default::default()),
        UnitKind::Selector => S::Selector(Default::default()),
        UnitKind::StreamCutter => S::StreamCutter(Default::default()),
        UnitKind::Tank => S::Tank(Default::default()),
        UnitKind::Adjust => S::Adjust(Default::default()),
        UnitKind::Set => S::Set(Default::default()),
        UnitKind::Spreadsheet => S::Spreadsheet(Default::default()),
        UnitKind::Balance => S::Balance(Default::default()),
        UnitKind::RecycleBlock => S::RecycleBlock(Default::default()),
        UnitKind::ControlBlock => S::ControlBlock(Default::default()),
        UnitKind::UserUnit => S::UserUnit(Default::default()),
        _ => unreachable!("{kind:?} has no SimpleLogicalSpec"),
    }
}

impl GroupSpec for SimpleLogicalSpec {
    fn kind(&self) -> UnitKind {
        use SimpleLogicalSpec as S;
        match self {
            S::ComponentSplitter(_) => UnitKind::ComponentSplitter,
            S::Multiplier(_) => UnitKind::Multiplier,
            S::Selector(_) => UnitKind::Selector,
            S::StreamCutter(_) => UnitKind::StreamCutter,
            S::Tank(_) => UnitKind::Tank,
            S::Adjust(_) => UnitKind::Adjust,
            S::Set(_) => UnitKind::Set,
            S::Spreadsheet(_) => UnitKind::Spreadsheet,
            S::Balance(_) => UnitKind::Balance,
            S::RecycleBlock(_) => UnitKind::RecycleBlock,
            S::ControlBlock(_) => UnitKind::ControlBlock,
            S::UserUnit(_) => UnitKind::UserUnit,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        match self {
            SimpleLogicalSpec::ComponentSplitter(s) => s.issues(),
            SimpleLogicalSpec::Selector(s) if !(2..=MAX_SELECTOR_INLETS).contains(&s.inlets) => {
                vec![SpecIssue {
                    field: "inlets",
                    problem: SpecProblem::CountOutOfRange,
                }]
            }
            _ => Vec::new(),
        }
    }

    fn variable_ports(&self) -> Vec<PortSpec> {
        match self {
            SimpleLogicalSpec::Selector(s) => numbered(PortSpec::inlet("in"), s.inlets).collect(),
            _ => Vec::new(),
        }
    }
}

/// Calculates a simple/logical unit.
pub(crate) fn calculate(
    kind: UnitKind,
    spec: &SimpleLogicalSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    match spec {
        SimpleLogicalSpec::ComponentSplitter(s) => {
            component_splitter::calculate(s, &inputs.material[0], ctx)
        }
        _ => Err(super::not_implemented(kind)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UnitGroup;

    #[test]
    fn simple_logical_kinds_follow_the_scaffolding_contract() {
        super::super::testing::check_drawn_only_contract(UnitGroup::SimpleLogical);
    }

    #[test]
    fn selector_inlets_are_numbered() {
        let spec = SimpleLogicalSpec::Selector(SelectorSpec { inlets: 3 });
        let labels: Vec<_> = spec.variable_ports().iter().map(|p| p.label()).collect();
        assert_eq!(labels, ["in 1", "in 2", "in 3"]);
    }
}
