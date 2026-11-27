//! Project records of the columns group (#185–#187, #190); owned by the
//! columns work package. Add fields with `#[serde(default)]` behind format
//! version 4.

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::columns::{
    AbsorberSpec, BottomSpecification, ColumnInternals, ColumnPressureDrop, ComponentTarget,
    DistillationColumnSpec, DistillationCondenser, LiquidSideDraw, MAX_COLUMN_PORTS,
    PetroleumFractionatorSpec, ShortcutColumnSpec, ShortcutCondenser, ShortcutReflux,
    TopSpecification, VaporSideDraw,
};
use serde::{Deserialize, Serialize};

use super::count_in;

fn one() -> u16 {
    1
}

fn unset_reflux() -> ShortcutRefluxRecord {
    ShortcutRefluxRecord::MinimumMultiple(None)
}

/// Shortcut column reflux or stage-count choice (#185).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ShortcutRefluxRecord {
    Ratio(Option<f64>),
    MinimumMultiple(Option<f64>),
    /// Added behind format version 4.
    TheoreticalStages(Option<f64>),
}

impl From<ShortcutReflux> for ShortcutRefluxRecord {
    fn from(r: ShortcutReflux) -> Self {
        match r {
            ShortcutReflux::Ratio(v) => Self::Ratio(v),
            ShortcutReflux::MinimumMultiple(v) => Self::MinimumMultiple(v),
            ShortcutReflux::TheoreticalStages(v) => Self::TheoreticalStages(v),
        }
    }
}

impl From<ShortcutRefluxRecord> for ShortcutReflux {
    fn from(r: ShortcutRefluxRecord) -> Self {
        match r {
            ShortcutRefluxRecord::Ratio(v) => Self::Ratio(v),
            ShortcutRefluxRecord::MinimumMultiple(v) => Self::MinimumMultiple(v),
            ShortcutRefluxRecord::TheoreticalStages(v) => Self::TheoreticalStages(v),
        }
    }
}

/// Shortcut column condenser (#185), added behind format version 4. Records
/// without it load as a total condenser, and a total condenser is not written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ShortcutCondenserRecord {
    #[default]
    Total,
    Partial,
}

impl ShortcutCondenserRecord {
    fn is_total(&self) -> bool {
        *self == Self::Total
    }
}

impl From<ShortcutCondenser> for ShortcutCondenserRecord {
    fn from(c: ShortcutCondenser) -> Self {
        match c {
            ShortcutCondenser::Total => Self::Total,
            ShortcutCondenser::Partial => Self::Partial,
        }
    }
}

impl From<ShortcutCondenserRecord> for ShortcutCondenser {
    fn from(c: ShortcutCondenserRecord) -> Self {
        match c {
            ShortcutCondenserRecord::Total => Self::Total,
            ShortcutCondenserRecord::Partial => Self::Partial,
        }
    }
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

fn is_zero(v: &f64) -> bool {
    *v == 0.0
}

fn is_zero_count(v: &u16) -> bool {
    *v == 0
}

/// Rigorous column condenser (#186), added behind format version 4. Records
/// without it load as a total condenser, and a total condenser is not written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum DistillationCondenserRecord {
    #[default]
    Total,
    PartialVapor,
    PartialVaporLiquid {
        liquid_distillate_mol_s: Option<f64>,
    },
}

impl DistillationCondenserRecord {
    fn is_total(&self) -> bool {
        *self == Self::Total
    }
}

/// Stage pressure drop below the first tray (#186), added behind format
/// version 4. Records without it load as zero drop (uniform pressure).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum PressureDropRecord {
    PerStagePa(f64),
    ColumnPa(f64),
}

impl Default for PressureDropRecord {
    fn default() -> Self {
        Self::PerStagePa(0.0)
    }
}

impl PressureDropRecord {
    fn is_zero(&self) -> bool {
        *self == Self::PerStagePa(0.0)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ComponentTargetRecord {
    component: Option<String>,
    value: Option<f64>,
}

impl From<&ComponentTarget> for ComponentTargetRecord {
    fn from(t: &ComponentTarget) -> Self {
        Self {
            component: t.component.clone(),
            value: t.value,
        }
    }
}

impl From<ComponentTargetRecord> for ComponentTarget {
    fn from(t: ComponentTargetRecord) -> Self {
        Self {
            component: t.component,
            value: t.value,
        }
    }
}

/// First column specification other than the reflux ratio, which keeps the
/// original `reflux_ratio` field (#186, added behind format version 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum TopSpecificationRecord {
    RefluxRateMolS(Option<f64>),
    DistillateRateMolS(Option<f64>),
    CondenserDutyW(Option<f64>),
    DistillatePurity(ComponentTargetRecord),
    DistillateRecovery(ComponentTargetRecord),
}

/// Second column specification other than the distillate fraction, which
/// keeps the original `distillate_fraction` field (#186, format version 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum BottomSpecificationRecord {
    BottomsRateMolS(Option<f64>),
    ReboilerDutyW(Option<f64>),
    BoilupRatio(Option<f64>),
    BottomsPurity(ComponentTargetRecord),
    BottomsRecovery(ComponentTargetRecord),
}

impl TopSpecificationRecord {
    fn number(&self) -> Option<f64> {
        match self {
            Self::RefluxRateMolS(v) | Self::DistillateRateMolS(v) | Self::CondenserDutyW(v) => *v,
            Self::DistillatePurity(t) | Self::DistillateRecovery(t) => t.value,
        }
    }
}

impl BottomSpecificationRecord {
    fn number(&self) -> Option<f64> {
        match self {
            Self::BottomsRateMolS(v) | Self::ReboilerDutyW(v) | Self::BoilupRatio(v) => *v,
            Self::BottomsPurity(t) | Self::BottomsRecovery(t) => t.value,
        }
    }
}

/// Splits a top specification into the original reflux-ratio field and the
/// specification record.
fn top_record(top: &TopSpecification) -> (Option<f64>, Option<TopSpecificationRecord>) {
    use TopSpecificationRecord as R;
    match top {
        TopSpecification::RefluxRatio(v) => (*v, None),
        TopSpecification::RefluxRate(v) => (None, Some(R::RefluxRateMolS(*v))),
        TopSpecification::DistillateRate(v) => (None, Some(R::DistillateRateMolS(*v))),
        TopSpecification::CondenserDuty(v) => (None, Some(R::CondenserDutyW(*v))),
        TopSpecification::DistillatePurity(t) => (None, Some(R::DistillatePurity(t.into()))),
        TopSpecification::DistillateRecovery(t) => (None, Some(R::DistillateRecovery(t.into()))),
    }
}

fn bottom_record(bottom: &BottomSpecification) -> (Option<f64>, Option<BottomSpecificationRecord>) {
    use BottomSpecificationRecord as R;
    match bottom {
        BottomSpecification::DistillateFraction(v) => (*v, None),
        BottomSpecification::BottomsRate(v) => (None, Some(R::BottomsRateMolS(*v))),
        BottomSpecification::ReboilerDuty(v) => (None, Some(R::ReboilerDutyW(*v))),
        BottomSpecification::BoilupRatio(v) => (None, Some(R::BoilupRatio(*v))),
        BottomSpecification::BottomsPurity(t) => (None, Some(R::BottomsPurity(t.into()))),
        BottomSpecification::BottomsRecovery(t) => (None, Some(R::BottomsRecovery(t.into()))),
    }
}

fn top_spec(
    reflux_ratio: Option<f64>,
    record: Option<TopSpecificationRecord>,
) -> Result<TopSpecification, String> {
    use TopSpecificationRecord as R;
    Ok(match record {
        None => TopSpecification::RefluxRatio(reflux_ratio),
        Some(_) if reflux_ratio.is_some() => {
            return Err("column has both a reflux ratio and another top specification".into());
        }
        Some(R::RefluxRateMolS(v)) => TopSpecification::RefluxRate(v),
        Some(R::DistillateRateMolS(v)) => TopSpecification::DistillateRate(v),
        Some(R::CondenserDutyW(v)) => TopSpecification::CondenserDuty(v),
        Some(R::DistillatePurity(t)) => TopSpecification::DistillatePurity(t.into()),
        Some(R::DistillateRecovery(t)) => TopSpecification::DistillateRecovery(t.into()),
    })
}

fn bottom_spec(
    distillate_fraction: Option<f64>,
    record: Option<BottomSpecificationRecord>,
) -> Result<BottomSpecification, String> {
    use BottomSpecificationRecord as R;
    Ok(match record {
        None => BottomSpecification::DistillateFraction(distillate_fraction),
        Some(_) if distillate_fraction.is_some() => {
            return Err(
                "column has both a distillate fraction and another bottom specification".into(),
            );
        }
        Some(R::BottomsRateMolS(v)) => BottomSpecification::BottomsRate(v),
        Some(R::ReboilerDutyW(v)) => BottomSpecification::ReboilerDuty(v),
        Some(R::BoilupRatio(v)) => BottomSpecification::BoilupRatio(v),
        Some(R::BottomsPurity(t)) => BottomSpecification::BottomsPurity(t.into()),
        Some(R::BottomsRecovery(t)) => BottomSpecification::BottomsRecovery(t.into()),
    })
}

fn draw_records(draws: impl Iterator<Item = (u16, Option<f64>)>) -> Vec<LiquidDrawRecord> {
    draws
        .map(|(stage, molar_flow_mol_s)| LiquidDrawRecord {
            stage,
            molar_flow_mol_s,
        })
        .collect()
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
        #[serde(default, skip_serializing_if = "ShortcutCondenserRecord::is_total")]
        condenser: ShortcutCondenserRecord,
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
        // The fields below were added behind format version 4 (#186
        // specification set); records without them load as the original
        // total-condenser, uniform-pressure column.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        top_specification: Option<TopSpecificationRecord>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        bottom_specification: Option<BottomSpecificationRecord>,
        #[serde(default, skip_serializing_if = "DistillationCondenserRecord::is_total")]
        condenser: DistillationCondenserRecord,
        #[serde(default, skip_serializing_if = "is_zero")]
        condenser_pressure_drop_pa: f64,
        #[serde(default, skip_serializing_if = "PressureDropRecord::is_zero")]
        pressure_drop: PressureDropRecord,
        #[serde(default, skip_serializing_if = "is_zero_count")]
        vapor_side_draws: u16,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        vapor_draws: Vec<LiquidDrawRecord>,
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
                    ShortcutRefluxRecord::Ratio(v)
                    | ShortcutRefluxRecord::MinimumMultiple(v)
                    | ShortcutRefluxRecord::TheoreticalStages(v) => *v,
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
                top_specification,
                bottom_specification,
                condenser,
                condenser_pressure_drop_pa,
                pressure_drop,
                vapor_draws,
                ..
            } => {
                let liquid_distillate = match condenser {
                    DistillationCondenserRecord::PartialVaporLiquid {
                        liquid_distillate_mol_s,
                    } => *liquid_distillate_mol_s,
                    _ => None,
                };
                let (PressureDropRecord::PerStagePa(drop) | PressureDropRecord::ColumnPa(drop)) =
                    *pressure_drop;
                let mut values: Vec<f64> = [
                    *pressure_pa,
                    *reflux_ratio,
                    *distillate_fraction,
                    top_specification.as_ref().and_then(|t| t.number()),
                    bottom_specification.as_ref().and_then(|b| b.number()),
                    liquid_distillate,
                    Some(*condenser_pressure_drop_pa),
                    Some(drop),
                ]
                .into_iter()
                .flatten()
                .collect();
                values.extend(liquid_draws.iter().filter_map(|draw| draw.molar_flow_mol_s));
                values.extend(vapor_draws.iter().filter_map(|draw| draw.molar_flow_mol_s));
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
                reflux: x.reflux.into(),
                condenser: x.condenser.into(),
            },
            ColumnSpec::Distillation(x) => {
                let (reflux_ratio, top_specification) = top_record(&x.top);
                let (distillate_fraction, bottom_specification) = bottom_record(&x.bottom);
                Self::DistillationColumn {
                    feeds: x.feeds,
                    side_draws: x.side_draws,
                    feed_stages: x.feed_stages.clone(),
                    liquid_draws: draw_records(
                        x.liquid_draws.iter().map(|d| (d.stage, d.molar_flow_mol_s)),
                    ),
                    internals: (x.internals != DistillationColumnSpec::default().internals)
                        .then(|| x.internals.into()),
                    pressure_pa: x.pressure_pa,
                    reflux_ratio,
                    distillate_fraction,
                    top_specification,
                    bottom_specification,
                    condenser: match x.condenser {
                        DistillationCondenser::Total => DistillationCondenserRecord::Total,
                        DistillationCondenser::PartialVapor => {
                            DistillationCondenserRecord::PartialVapor
                        }
                        DistillationCondenser::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        } => DistillationCondenserRecord::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        },
                    },
                    condenser_pressure_drop_pa: x.condenser_pressure_drop_pa,
                    pressure_drop: match x.pressure_drop {
                        ColumnPressureDrop::PerStage(v) => PressureDropRecord::PerStagePa(v),
                        ColumnPressureDrop::Column(v) => PressureDropRecord::ColumnPa(v),
                    },
                    vapor_side_draws: x.vapor_side_draws,
                    vapor_draws: draw_records(
                        x.vapor_draws.iter().map(|d| (d.stage, d.molar_flow_mol_s)),
                    ),
                }
            }
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
                condenser,
            } => Self::Shortcut(ShortcutColumnSpec {
                light_key,
                heavy_key,
                light_key_recovery,
                heavy_key_recovery,
                pressure_pa,
                reflux: reflux.into(),
                condenser: condenser.into(),
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
                top_specification,
                bottom_specification,
                condenser,
                condenser_pressure_drop_pa,
                pressure_drop,
                vapor_side_draws,
                vapor_draws,
            } => {
                if feed_stages.len() > usize::from(feeds)
                    || liquid_draws.len() > usize::from(side_draws)
                    || vapor_draws.len() > usize::from(vapor_side_draws)
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
                    vapor_side_draws: count_in(
                        vapor_side_draws,
                        "vapor_side_draws",
                        0..=MAX_COLUMN_PORTS,
                    )?,
                    vapor_draws: vapor_draws
                        .into_iter()
                        .map(|d| VaporSideDraw {
                            stage: d.stage,
                            molar_flow_mol_s: d.molar_flow_mol_s,
                        })
                        .collect(),
                    pressure_pa,
                    condenser_pressure_drop_pa,
                    pressure_drop: match pressure_drop {
                        PressureDropRecord::PerStagePa(v) => ColumnPressureDrop::PerStage(v),
                        PressureDropRecord::ColumnPa(v) => ColumnPressureDrop::Column(v),
                    },
                    condenser: match condenser {
                        DistillationCondenserRecord::Total => DistillationCondenser::Total,
                        DistillationCondenserRecord::PartialVapor => {
                            DistillationCondenser::PartialVapor
                        }
                        DistillationCondenserRecord::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        } => DistillationCondenser::PartialVaporLiquid {
                            liquid_distillate_mol_s,
                        },
                    },
                    top: top_spec(reflux_ratio, top_specification)?,
                    bottom: bottom_spec(distillate_fraction, bottom_specification)?,
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
