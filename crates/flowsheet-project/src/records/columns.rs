//! Project records of the columns group (#185–#187, #190); owned by the
//! columns work package. Add fields with `#[serde(default)]` behind format
//! version 4.

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::columns::{
    AbsorberSpec, ColumnInternals, DistillationColumnSpec, LiquidSideDraw, MAX_COLUMN_PORTS,
    PetroleumFractionatorSpec, ShortcutColumnSpec, ShortcutReflux,
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

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum InternalsRecord {
    Trays(u16),
    Packed {
        hetp_m: Option<f64>,
        basis: PackingRecord,
    },
}

impl InternalsRecord {
    fn numbers(self) -> impl Iterator<Item = f64> {
        let (hetp, height) = match self {
            Self::Trays(_) => (None, None),
            Self::Packed { hetp_m, basis } => (
                hetp_m,
                match basis {
                    PackingRecord::Height { height_m, .. } => height_m,
                    PackingRecord::TheoreticalStages(_) => None,
                },
            ),
        };
        [hetp, height].into_iter().flatten()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PackingRecord {
    Height {
        height_m: Option<f64>,
        rounding: RoundingRecord,
    },
    TheoreticalStages(u16),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum RoundingRecord {
    Nearest,
    Down,
    Up,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct LiquidDrawRecord {
    stage: u16,
    molar_flow_mol_s: Option<f64>,
}

impl From<ColumnInternals> for InternalsRecord {
    fn from(value: ColumnInternals) -> Self {
        match value {
            ColumnInternals::Trays(count) => Self::Trays(count),
            ColumnInternals::Packed(spec) => Self::Packed {
                hetp_m: spec.hetp_m,
                basis: match spec.basis {
                    PackedSectionBasis::TheoreticalStages(count) => {
                        PackingRecord::TheoreticalStages(count)
                    }
                    PackedSectionBasis::PackedHeight { height_m, rounding } => {
                        PackingRecord::Height {
                            height_m,
                            rounding: match rounding {
                                StageRounding::Nearest => RoundingRecord::Nearest,
                                StageRounding::Down => RoundingRecord::Down,
                                StageRounding::Up => RoundingRecord::Up,
                            },
                        }
                    }
                },
            },
        }
    }
}

impl From<InternalsRecord> for ColumnInternals {
    fn from(value: InternalsRecord) -> Self {
        match value {
            InternalsRecord::Trays(count) => Self::Trays(count),
            InternalsRecord::Packed { hetp_m, basis } => Self::Packed(PackedSectionSpec {
                hetp_m,
                basis: match basis {
                    PackingRecord::TheoreticalStages(count) => {
                        PackedSectionBasis::TheoreticalStages(count)
                    }
                    PackingRecord::Height { height_m, rounding } => {
                        PackedSectionBasis::PackedHeight {
                            height_m,
                            rounding: match rounding {
                                RoundingRecord::Nearest => StageRounding::Nearest,
                                RoundingRecord::Down => StageRounding::Down,
                                RoundingRecord::Up => StageRounding::Up,
                            },
                        }
                    }
                },
            }),
        }
    }
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
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        feed_stages: Vec<u16>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        liquid_draws: Vec<LiquidDrawRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pressure_pa: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reflux_ratio: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        distillate_fraction: Option<f64>,
    },
    Absorber {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        pressure_pa: Option<f64>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
    },
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
            Self::DistillationColumn {
                liquid_draws,
                internals,
                pressure_pa,
                reflux_ratio,
                distillate_fraction,
                ..
            } => {
                let mut values: Vec<f64> = [*pressure_pa, *reflux_ratio, *distillate_fraction]
                    .into_iter()
                    .flatten()
                    .collect();
                values.extend(liquid_draws.iter().filter_map(|draw| draw.molar_flow_mol_s));
                values.extend(internals.iter().flat_map(|i| i.numbers()));
                values
            }
            Self::Absorber {
                pressure_pa,
                internals,
            } => pressure_pa
                .iter()
                .copied()
                .chain(internals.iter().flat_map(|i| i.numbers()))
                .collect(),
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
                feed_stages: x.feed_stages.clone(),
                liquid_draws: x
                    .liquid_draws
                    .iter()
                    .map(|d| LiquidDrawRecord {
                        stage: d.stage,
                        molar_flow_mol_s: d.molar_flow_mol_s,
                    })
                    .collect(),
                internals: (x.internals != DistillationColumnSpec::default().internals)
                    .then(|| x.internals.into()),
                pressure_pa: x.pressure_pa,
                reflux_ratio: x.reflux_ratio,
                distillate_fraction: x.distillate_fraction,
            },
            ColumnSpec::Absorber(x) => Self::Absorber {
                pressure_pa: x.pressure_pa,
                internals: (x.internals != AbsorberSpec::default().internals)
                    .then(|| x.internals.into()),
            },
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
            ColumnRecord::DistillationColumn {
                feeds,
                side_draws,
                feed_stages,
                liquid_draws,
                internals,
                pressure_pa,
                reflux_ratio,
                distillate_fraction,
            } => {
                if feed_stages.len() > usize::from(feeds)
                    || liquid_draws.len() > usize::from(side_draws)
                {
                    return Err("column has more feed or draw specifications than ports".into());
                }
                Self::Distillation(DistillationColumnSpec {
                    feeds: count_in(feeds, "feeds", 1..=MAX_COLUMN_PORTS)?,
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                    feed_stages,
                    liquid_draws: liquid_draws
                        .into_iter()
                        .map(|d| LiquidSideDraw {
                            stage: d.stage,
                            molar_flow_mol_s: d.molar_flow_mol_s,
                        })
                        .collect(),
                    internals: internals
                        .map_or_else(|| DistillationColumnSpec::default().internals, Into::into),
                    pressure_pa,
                    reflux_ratio,
                    distillate_fraction,
                })
            }
            ColumnRecord::Absorber {
                pressure_pa,
                internals,
            } => Self::Absorber(AbsorberSpec {
                pressure_pa,
                internals: internals.map_or_else(|| AbsorberSpec::default().internals, Into::into),
            }),
            ColumnRecord::PetroleumFractionator { side_draws } => {
                Self::PetroleumFractionator(PetroleumFractionatorSpec {
                    side_draws: count_in(side_draws, "side_draws", 0..=MAX_COLUMN_PORTS)?,
                })
            }
        })
    }
}
