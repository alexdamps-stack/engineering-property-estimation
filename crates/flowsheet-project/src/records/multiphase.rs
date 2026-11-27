//! Project records of the multiphase group (#176, #177, #188, #189); owned
//! by the multiphase work package. Add fields with `#[serde(default)]`
//! behind format version 4.
//!
//! Liquid–liquid extractor fields (#189) were added behind format version 4;
//! a version-4 record without them (`{"kind": "liquid_extractor"}`, written
//! while the kind was drawn only) loads with the default specification.

use flowsheet_core::units::MultiphaseSpec;
use flowsheet_core::units::columns::packing::{
    PackedSectionBasis, PackedSectionSpec, StageRounding,
};
use flowsheet_core::units::multiphase::{
    ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec,
};
use serde::{Deserialize, Serialize};

/// Multiphase specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum MultiphaseRecord {
    ThreePhaseSeparator {},
    Decanter {},
    ThreePhaseColumn {},
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
            _ => Vec::new(),
        }
    }
}

impl From<&MultiphaseSpec> for MultiphaseRecord {
    fn from(s: &MultiphaseSpec) -> Self {
        match s {
            MultiphaseSpec::ThreePhaseSeparator(_) => Self::ThreePhaseSeparator {},
            MultiphaseSpec::Decanter(_) => Self::Decanter {},
            MultiphaseSpec::ThreePhaseColumn(_) => Self::ThreePhaseColumn {},
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
            MultiphaseRecord::ThreePhaseColumn {} => Self::ThreePhaseColumn(Default::default()),
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
}
