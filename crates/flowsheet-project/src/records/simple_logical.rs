//! Project records of the simple/logical group (#172–#175, #178,
//! #212–#218); owned by the simple/logical work package. Add fields with
//! `#[serde(default)]` behind format version 4.

use flowsheet_core::units::SimpleLogicalSpec;
use flowsheet_core::units::simple_logical::{MAX_SELECTOR_INLETS, SelectorSpec};
use serde::{Deserialize, Serialize};

use super::count_in;

fn default_inlets() -> u16 {
    SelectorSpec::default().inlets
}

/// Simple/logical specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum SimpleLogicalRecord {
    ComponentSplitter {},
    Multiplier {},
    Selector {
        #[serde(default = "default_inlets")]
        inlets: u16,
    },
    StreamCutter {},
    Tank {},
    Adjust {},
    Set {},
    Spreadsheet {},
    Balance {},
    RecycleBlock {},
    ControlBlock {},
    UserUnit {},
}

impl SimpleLogicalRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&SimpleLogicalSpec> for SimpleLogicalRecord {
    fn from(s: &SimpleLogicalSpec) -> Self {
        use SimpleLogicalSpec as S;
        match s {
            S::ComponentSplitter(_) => Self::ComponentSplitter {},
            S::Multiplier(_) => Self::Multiplier {},
            S::Selector(x) => Self::Selector { inlets: x.inlets },
            S::StreamCutter(_) => Self::StreamCutter {},
            S::Tank(_) => Self::Tank {},
            S::Adjust(_) => Self::Adjust {},
            S::Set(_) => Self::Set {},
            S::Spreadsheet(_) => Self::Spreadsheet {},
            S::Balance(_) => Self::Balance {},
            S::RecycleBlock(_) => Self::RecycleBlock {},
            S::ControlBlock(_) => Self::ControlBlock {},
            S::UserUnit(_) => Self::UserUnit {},
        }
    }
}

impl TryFrom<SimpleLogicalRecord> for SimpleLogicalSpec {
    type Error = String;

    fn try_from(r: SimpleLogicalRecord) -> Result<Self, String> {
        use SimpleLogicalRecord as R;
        Ok(match r {
            R::ComponentSplitter {} => Self::ComponentSplitter(Default::default()),
            R::Multiplier {} => Self::Multiplier(Default::default()),
            R::Selector { inlets } => Self::Selector(SelectorSpec {
                inlets: count_in(inlets, "inlets", 2..=MAX_SELECTOR_INLETS)?,
            }),
            R::StreamCutter {} => Self::StreamCutter(Default::default()),
            R::Tank {} => Self::Tank(Default::default()),
            R::Adjust {} => Self::Adjust(Default::default()),
            R::Set {} => Self::Set(Default::default()),
            R::Spreadsheet {} => Self::Spreadsheet(Default::default()),
            R::Balance {} => Self::Balance(Default::default()),
            R::RecycleBlock {} => Self::RecycleBlock(Default::default()),
            R::ControlBlock {} => Self::ControlBlock(Default::default()),
            R::UserUnit {} => Self::UserUnit(Default::default()),
        })
    }
}
