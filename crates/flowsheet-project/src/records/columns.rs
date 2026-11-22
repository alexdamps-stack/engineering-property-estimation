//! Project records of the columns group (#185–#187, #190); owned by the
//! columns work package. Add fields with `#[serde(default)]` behind format
//! version 4.

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    DistillationColumnSpec, MAX_COLUMN_PORTS, PetroleumFractionatorSpec,
};
use serde::{Deserialize, Serialize};

use super::count_in;

fn one() -> u16 {
    1
}

/// Column specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ColumnRecord {
    ShortcutColumn {},
    DistillationColumn {
        #[serde(default = "one")]
        feeds: u16,
        #[serde(default)]
        side_draws: u16,
    },
    Absorber {},
    PetroleumFractionator {
        #[serde(default)]
        side_draws: u16,
    },
}

impl ColumnRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&ColumnSpec> for ColumnRecord {
    fn from(s: &ColumnSpec) -> Self {
        match s {
            ColumnSpec::Shortcut(_) => Self::ShortcutColumn {},
            ColumnSpec::Distillation(x) => Self::DistillationColumn {
                feeds: x.feeds,
                side_draws: x.side_draws,
            },
            ColumnSpec::Absorber(_) => Self::Absorber {},
            ColumnSpec::PetroleumFractionator(x) => Self::PetroleumFractionator {
                side_draws: x.side_draws,
            },
        }
    }
}

impl TryFrom<ColumnRecord> for ColumnSpec {
    type Error = String;

    fn try_from(r: ColumnRecord) -> Result<Self, String> {
        Ok(match r {
            ColumnRecord::ShortcutColumn {} => Self::Shortcut(Default::default()),
            ColumnRecord::DistillationColumn { feeds, side_draws } => {
                Self::Distillation(DistillationColumnSpec {
                    feeds: count_in(feeds, "feeds", 1..=MAX_COLUMN_PORTS)?,
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                })
            }
            ColumnRecord::Absorber {} => Self::Absorber(Default::default()),
            ColumnRecord::PetroleumFractionator { side_draws } => {
                Self::PetroleumFractionator(PetroleumFractionatorSpec {
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                })
            }
        })
    }
}
