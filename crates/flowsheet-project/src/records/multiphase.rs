//! Project records of the multiphase group (#176, #177, #188, #189); owned
//! by the multiphase work package. Add fields with `#[serde(default)]`
//! behind format version 4.
//!
//! Liquid–liquid extractor fields (#189) were added behind format version 4;
//! a version-4 record without them (`{"kind": "liquid_extractor"}`, written
//! while the kind was drawn only) loads with the default specification.
//!
//! Three-phase column fields (#188) were added behind format version 4 in the
//! same way: `{"kind": "three_phase_column"}` opens with five ideal trays,
//! feed at stage 2, light-liquid total reflux and unset values.

use flowsheet_core::units::MultiphaseSpec;
use flowsheet_core::units::columns::ColumnInternals;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::multiphase::{
    DecanterReflux, ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec,
    PhaseReflux, ThreePhaseColumnSpec,
};
use serde::{Deserialize, Serialize};

use super::columns::InternalsRecord;

/// Multiphase specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum MultiphaseRecord {
    ThreePhaseSeparator {},
    Decanter {},
    ThreePhaseColumn {
        /// `None` (a drawn-only record) opens with five ideal trays.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        internals: Option<InternalsRecord>,
        #[serde(default = "default_feed_stage")]
        feed_stage: u16,
        #[serde(default)]
        pressure_pa: Option<f64>,
        #[serde(default)]
        decanter_temperature_k: Option<f64>,
        #[serde(default)]
        reflux: DecanterRefluxRecord,
        #[serde(default)]
        bottoms_fraction: Option<f64>,
    },
    LiquidExtractor {
        #[serde(default)]
        internals: ExtractorInternalsRecord,
        #[serde(default)]
        pressure_pa: Option<f64>,
        #[serde(default)]
        thermal: ExtractorThermalRecord,
        #[serde(default)]
        feed_phase: FeedPhaseRecord,
    },
}

fn default_feed_stage() -> u16 {
    ThreePhaseColumnSpec::default().feed_stage
}

/// Reflux phase and amount of the three-phase column decanter (#188).
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum DecanterRefluxRecord {
    #[default]
    LightTotal,
    HeavyTotal,
    LightRatio {
        #[serde(default)]
        ratio: Option<f64>,
    },
    HeavyRatio {
        #[serde(default)]
        ratio: Option<f64>,
    },
    Fractions {
        #[serde(default)]
        light: Option<f64>,
        #[serde(default)]
        heavy: Option<f64>,
    },
}

impl DecanterRefluxRecord {
    fn numbers(self) -> Vec<f64> {
        match self {
            Self::LightTotal | Self::HeavyTotal => Vec::new(),
            Self::LightRatio { ratio } | Self::HeavyRatio { ratio } => ratio.into_iter().collect(),
            Self::Fractions { light, heavy } => light.into_iter().chain(heavy).collect(),
        }
    }
}

impl From<DecanterReflux> for DecanterRefluxRecord {
    fn from(value: DecanterReflux) -> Self {
        match value {
            DecanterReflux::Light(PhaseReflux::Total) => Self::LightTotal,
            DecanterReflux::Heavy(PhaseReflux::Total) => Self::HeavyTotal,
            DecanterReflux::Light(PhaseReflux::Ratio(ratio)) => Self::LightRatio { ratio },
            DecanterReflux::Heavy(PhaseReflux::Ratio(ratio)) => Self::HeavyRatio { ratio },
            DecanterReflux::Fractions { light, heavy } => Self::Fractions { light, heavy },
        }
    }
}

impl From<DecanterRefluxRecord> for DecanterReflux {
    fn from(value: DecanterRefluxRecord) -> Self {
        match value {
            DecanterRefluxRecord::LightTotal => Self::Light(PhaseReflux::Total),
            DecanterRefluxRecord::HeavyTotal => Self::Heavy(PhaseReflux::Total),
            DecanterRefluxRecord::LightRatio { ratio } => Self::Light(PhaseReflux::Ratio(ratio)),
            DecanterRefluxRecord::HeavyRatio { ratio } => Self::Heavy(PhaseReflux::Ratio(ratio)),
            DecanterRefluxRecord::Fractions { light, heavy } => Self::Fractions { light, heavy },
        }
    }
}

/// Every number of a column internals record.
fn internals_numbers(internals: Option<InternalsRecord>) -> Vec<f64> {
    match internals.map(ColumnInternals::from) {
        Some(ColumnInternals::Packed(spec)) => {
            let mut out: Vec<f64> = spec.hetp_m.into_iter().collect();
            if let PackedSectionBasis::PackedHeight {
                height_m: Some(h), ..
            } = spec.basis
            {
                out.push(h);
            }
            out
        }
        _ => Vec::new(),
    }
}

/// Ideal stages or one packed section of the extractor.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ExtractorInternalsRecord {
    Stages(u16),
    Packed {
        hetp_m: Option<f64>,
        basis: ExtractorPackingRecord,
    },
}

impl Default for ExtractorInternalsRecord {
    fn default() -> Self {
        (&ExtractorInternals::default()).into()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ExtractorPackingRecord {
    Height {
        height_m: Option<f64>,
        rounding: ExtractorRoundingRecord,
    },
    TheoreticalStages(u16),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ExtractorRoundingRecord {
    Nearest,
    Down,
    Up,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ExtractorThermalRecord {
    Isothermal {
        #[serde(default)]
        temperature_k: Option<f64>,
    },
    Adiabatic,
}

impl Default for ExtractorThermalRecord {
    fn default() -> Self {
        Self::Isothermal {
            temperature_k: None,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum FeedPhaseRecord {
    #[default]
    Detect,
    Heavy,
    Light,
}

impl From<&ExtractorInternals> for ExtractorInternalsRecord {
    fn from(value: &ExtractorInternals) -> Self {
        match *value {
            ExtractorInternals::Stages(count) => Self::Stages(count),
            ExtractorInternals::Packed(spec) => Self::Packed {
                hetp_m: spec.hetp_m,
                basis: match spec.basis {
                    PackedSectionBasis::TheoreticalStages(count) => {
                        ExtractorPackingRecord::TheoreticalStages(count)
                    }
                    PackedSectionBasis::PackedHeight { height_m, rounding } => {
                        ExtractorPackingRecord::Height {
                            height_m,
                            rounding: match rounding {
                                StageRounding::Nearest => ExtractorRoundingRecord::Nearest,
                                StageRounding::Down => ExtractorRoundingRecord::Down,
                                StageRounding::Up => ExtractorRoundingRecord::Up,
                            },
                        }
                    }
                },
            },
        }
    }
}

impl From<ExtractorInternalsRecord> for ExtractorInternals {
    fn from(value: ExtractorInternalsRecord) -> Self {
        match value {
            ExtractorInternalsRecord::Stages(count) => Self::Stages(count),
            ExtractorInternalsRecord::Packed { hetp_m, basis } => Self::Packed(PackedSectionSpec {
                hetp_m,
                basis: match basis {
                    ExtractorPackingRecord::TheoreticalStages(count) => {
                        PackedSectionBasis::TheoreticalStages(count)
                    }
                    ExtractorPackingRecord::Height { height_m, rounding } => {
                        PackedSectionBasis::PackedHeight {
                            height_m,
                            rounding: match rounding {
                                ExtractorRoundingRecord::Nearest => StageRounding::Nearest,
                                ExtractorRoundingRecord::Down => StageRounding::Down,
                                ExtractorRoundingRecord::Up => StageRounding::Up,
                            },
                        }
                    }
                },
            }),
        }
    }
}

impl MultiphaseRecord {
    /// Every number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        match self {
            Self::LiquidExtractor {
                internals,
                pressure_pa,
                thermal,
                ..
            } => {
                let mut out: Vec<f64> = pressure_pa.iter().copied().collect();
                if let ExtractorThermalRecord::Isothermal {
                    temperature_k: Some(t),
                } = thermal
                {
                    out.push(*t);
                }
                if let ExtractorInternalsRecord::Packed { hetp_m, basis } = internals {
                    out.extend(hetp_m);
                    if let ExtractorPackingRecord::Height {
                        height_m: Some(h), ..
                    } = basis
                    {
                        out.push(*h);
                    }
                }
                out
            }
            Self::ThreePhaseColumn {
                internals,
                pressure_pa,
                decanter_temperature_k,
                reflux,
                bottoms_fraction,
                ..
            } => internals_numbers(*internals)
                .into_iter()
                .chain(*pressure_pa)
                .chain(*decanter_temperature_k)
                .chain(reflux.numbers())
                .chain(*bottoms_fraction)
                .collect(),
            _ => Vec::new(),
        }
    }
}

impl From<&MultiphaseSpec> for MultiphaseRecord {
    fn from(s: &MultiphaseSpec) -> Self {
        match s {
            MultiphaseSpec::ThreePhaseSeparator(_) => Self::ThreePhaseSeparator {},
            MultiphaseSpec::Decanter(_) => Self::Decanter {},
            MultiphaseSpec::ThreePhaseColumn(c) => Self::ThreePhaseColumn {
                internals: Some(c.internals.into()),
                feed_stage: c.feed_stage,
                pressure_pa: c.pressure_pa,
                decanter_temperature_k: c.decanter_temperature_k,
                reflux: c.reflux.into(),
                bottoms_fraction: c.bottoms_fraction,
            },
            MultiphaseSpec::LiquidExtractor(e) => Self::LiquidExtractor {
                internals: (&e.internals).into(),
                pressure_pa: e.pressure_pa,
                thermal: match e.thermal {
                    ExtractorThermal::Isothermal { temperature_k } => {
                        ExtractorThermalRecord::Isothermal { temperature_k }
                    }
                    ExtractorThermal::Adiabatic => ExtractorThermalRecord::Adiabatic,
                },
                feed_phase: match e.feed_phase {
                    FeedPhase::Detect => FeedPhaseRecord::Detect,
                    FeedPhase::Heavy => FeedPhaseRecord::Heavy,
                    FeedPhase::Light => FeedPhaseRecord::Light,
                },
            },
        }
    }
}

impl TryFrom<MultiphaseRecord> for MultiphaseSpec {
    type Error = String;

    fn try_from(r: MultiphaseRecord) -> Result<Self, String> {
        Ok(match r {
            MultiphaseRecord::ThreePhaseSeparator {} => {
                Self::ThreePhaseSeparator(Default::default())
            }
            MultiphaseRecord::Decanter {} => Self::Decanter(Default::default()),
            MultiphaseRecord::ThreePhaseColumn {
                internals,
                feed_stage,
                pressure_pa,
                decanter_temperature_k,
                reflux,
                bottoms_fraction,
            } => Self::ThreePhaseColumn(ThreePhaseColumnSpec {
                internals: internals
                    .map(Into::into)
                    .unwrap_or(ThreePhaseColumnSpec::default().internals),
                feed_stage,
                pressure_pa,
                decanter_temperature_k,
                reflux: reflux.into(),
                bottoms_fraction,
            }),
            MultiphaseRecord::LiquidExtractor {
                internals,
                pressure_pa,
                thermal,
                feed_phase,
            } => Self::LiquidExtractor(LiquidExtractorSpec {
                internals: internals.into(),
                pressure_pa,
                thermal: match thermal {
                    ExtractorThermalRecord::Isothermal { temperature_k } => {
                        ExtractorThermal::Isothermal { temperature_k }
                    }
                    ExtractorThermalRecord::Adiabatic => ExtractorThermal::Adiabatic,
                },
                feed_phase: match feed_phase {
                    FeedPhaseRecord::Detect => FeedPhase::Detect,
                    FeedPhaseRecord::Heavy => FeedPhase::Heavy,
                    FeedPhaseRecord::Light => FeedPhase::Light,
                },
            }),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trip(spec: LiquidExtractorSpec) {
        let spec = MultiphaseSpec::LiquidExtractor(spec);
        let record = MultiphaseRecord::from(&spec);
        let json = serde_json::to_string(&record).unwrap();
        let back: MultiphaseRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(MultiphaseSpec::try_from(back).unwrap(), spec);
    }

    #[test]
    fn extractor_specifications_round_trip() {
        round_trip(LiquidExtractorSpec::default());
        round_trip(LiquidExtractorSpec {
            internals: ExtractorInternals::Stages(12),
            pressure_pa: Some(2e5),
            thermal: ExtractorThermal::Isothermal {
                temperature_k: Some(310.0),
            },
            feed_phase: FeedPhase::Heavy,
        });
        round_trip(LiquidExtractorSpec {
            internals: ExtractorInternals::Packed(PackedSectionSpec {
                hetp_m: Some(0.5),
                basis: PackedSectionBasis::PackedHeight {
                    height_m: Some(3.2),
                    rounding: StageRounding::Down,
                },
            }),
            pressure_pa: Some(1e5),
            thermal: ExtractorThermal::Adiabatic,
            feed_phase: FeedPhase::Light,
        });
    }

    #[test]
    fn drawn_only_extractor_record_loads_with_defaults() {
        let record: MultiphaseRecord =
            serde_json::from_str(r#"{"kind": "liquid_extractor"}"#).unwrap();
        assert_eq!(
            MultiphaseSpec::try_from(record).unwrap(),
            MultiphaseSpec::LiquidExtractor(LiquidExtractorSpec::default())
        );
        assert!(
            serde_json::from_str::<MultiphaseRecord>(
                r#"{"kind": "liquid_extractor", "solvent_ratio": 1.0}"#
            )
            .is_err()
        );
    }

    fn column_round_trip(spec: ThreePhaseColumnSpec) {
        let spec = MultiphaseSpec::ThreePhaseColumn(spec);
        let record = MultiphaseRecord::from(&spec);
        let json = serde_json::to_string(&record).unwrap();
        let back: MultiphaseRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(MultiphaseSpec::try_from(back).unwrap(), spec);
    }

    #[test]
    fn three_phase_column_specifications_round_trip() {
        column_round_trip(ThreePhaseColumnSpec::default());
        column_round_trip(ThreePhaseColumnSpec {
            internals: ColumnInternals::Trays(9),
            feed_stage: 4,
            pressure_pa: Some(101_325.0),
            decanter_temperature_k: Some(320.0),
            reflux: DecanterReflux::Heavy(PhaseReflux::Ratio(Some(3.0))),
            bottoms_fraction: Some(0.7),
        });
        let spec = ThreePhaseColumnSpec {
            internals: ColumnInternals::Packed(PackedSectionSpec {
                hetp_m: Some(0.4),
                basis: PackedSectionBasis::PackedHeight {
                    height_m: Some(2.0),
                    rounding: StageRounding::Up,
                },
            }),
            reflux: DecanterReflux::Fractions {
                light: Some(0.9),
                heavy: None,
            },
            ..Default::default()
        };
        column_round_trip(spec.clone());
        let record = MultiphaseRecord::from(&MultiphaseSpec::ThreePhaseColumn(spec));
        assert_eq!(record.numbers(), [0.4, 2.0, 0.9]);
        column_round_trip(ThreePhaseColumnSpec {
            reflux: DecanterReflux::Light(PhaseReflux::Ratio(None)),
            ..Default::default()
        });
    }

    #[test]
    fn drawn_only_three_phase_column_record_loads_with_defaults() {
        let record: MultiphaseRecord =
            serde_json::from_str(r#"{"kind": "three_phase_column"}"#).unwrap();
        assert_eq!(
            MultiphaseSpec::try_from(record).unwrap(),
            MultiphaseSpec::ThreePhaseColumn(ThreePhaseColumnSpec::default())
        );
        assert!(
            serde_json::from_str::<MultiphaseRecord>(
                r#"{"kind": "three_phase_column", "boilup_ratio": 1.0}"#
            )
            .is_err()
        );
    }
}
