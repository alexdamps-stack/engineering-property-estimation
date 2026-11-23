//! Project records of the columns group (#185–#187, #190); owned by the
//! columns work package. Add fields with `#[serde(default)]` behind format
//! version 4.

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    DistillationColumnSpec, MAX_COLUMN_PORTS, PetroleumFractionatorSpec, ShortcutColumnSpec,
    ShortcutReflux,
};
use serde::{Deserialize, Serialize};

use super::count_in;

fn one() -> u16 {
    1
}

fn unset_reflux() -> ShortcutRefluxRecord {
    ShortcutRefluxRecord::MinimumMultiple(None)
}

/// Shortcut column reflux choice (#185).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ShortcutRefluxRecord {
    Ratio(Option<f64>),
    MinimumMultiple(Option<f64>),
}

/// Column specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ColumnRecord {
    /// Shortcut fields were added behind format version 4; the earlier empty
    /// record loads with every value unset.
    ShortcutColumn {
        #[serde(default)]
        light_key: Option<String>,
        #[serde(default)]
        heavy_key: Option<String>,
        #[serde(default)]
        light_key_recovery: Option<f64>,
        #[serde(default)]
        heavy_key_recovery: Option<f64>,
        #[serde(default)]
        pressure_pa: Option<f64>,
        #[serde(default = "unset_reflux")]
        reflux: ShortcutRefluxRecord,
    },
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
        match self {
            Self::ShortcutColumn {
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux,
                ..
            } => {
                let reflux = match reflux {
                    ShortcutRefluxRecord::Ratio(v) | ShortcutRefluxRecord::MinimumMultiple(v) => *v,
                };
                [
                    *light_key_recovery,
                    *heavy_key_recovery,
                    *pressure_pa,
                    reflux,
                ]
                .into_iter()
                .flatten()
                .collect()
            }
            _ => Vec::new(),
        }
    }
}

impl From<&ColumnSpec> for ColumnRecord {
    fn from(s: &ColumnSpec) -> Self {
        match s {
            ColumnSpec::Shortcut(x) => Self::ShortcutColumn {
                light_key: x.light_key.clone(),
                heavy_key: x.heavy_key.clone(),
                light_key_recovery: x.light_key_recovery,
                heavy_key_recovery: x.heavy_key_recovery,
                pressure_pa: x.pressure_pa,
                reflux: match x.reflux {
                    ShortcutReflux::Ratio(v) => ShortcutRefluxRecord::Ratio(v),
                    ShortcutReflux::MinimumMultiple(v) => ShortcutRefluxRecord::MinimumMultiple(v),
                },
            },
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
            ColumnRecord::ShortcutColumn {
                light_key,
                heavy_key,
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux,
            } => Self::Shortcut(ShortcutColumnSpec {
                light_key,
                heavy_key,
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux: match reflux {
                    ShortcutRefluxRecord::Ratio(v) => ShortcutReflux::Ratio(v),
                    ShortcutRefluxRecord::MinimumMultiple(v) => ShortcutReflux::MinimumMultiple(v),
                },
            }),
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
