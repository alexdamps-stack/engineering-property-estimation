//! Project records of the rotating group (#182, A14); owned by the rotating
//! work package. Pump, compressor and expander keep the core `pressurize` /
//! `expand` records in `records/mod.rs`. Add fields with `#[serde(default)]`
//! behind format version 4.

use flowsheet_core::DEFAULT_EFFICIENCY;
use flowsheet_core::units::RotatingSpec;
use flowsheet_core::units::rotating::{MAX_COMPRESSOR_STAGES, MultistageCompressorSpec};
use serde::{Deserialize, Serialize};

use super::count_in;

fn default_stages() -> u16 {
    MultistageCompressorSpec::default().stages
}

fn default_efficiency() -> f64 {
    DEFAULT_EFFICIENCY
}

/// Rotating-group specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum RotatingRecord {
    /// Fields added behind format version 4; an earlier empty record loads
    /// with two stages, the default efficiency and the pressures unset.
    MultistageCompressor {
        #[serde(default = "default_stages")]
        stages: u16,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        outlet_pressure_pa: Option<f64>,
        #[serde(default = "default_efficiency")]
        efficiency: f64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        intercooler_temperature_k: Option<f64>,
    },
}

impl RotatingRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        match self {
            Self::MultistageCompressor {
                outlet_pressure_pa,
                efficiency,
                intercooler_temperature_k,
                ..
            } => std::iter::once(*efficiency)
                .chain(*outlet_pressure_pa)
                .chain(*intercooler_temperature_k)
                .collect(),
        }
    }
}

impl From<&RotatingSpec> for RotatingRecord {
    fn from(s: &RotatingSpec) -> Self {
        match s {
            RotatingSpec::MultistageCompressor(s) => Self::MultistageCompressor {
                stages: s.stages,
                outlet_pressure_pa: s.outlet_pressure_pa,
                efficiency: s.efficiency,
                intercooler_temperature_k: s.intercooler_temperature_k,
            },
        }
    }
}

impl TryFrom<RotatingRecord> for RotatingSpec {
    type Error = String;

    fn try_from(r: RotatingRecord) -> Result<Self, String> {
        Ok(match r {
            RotatingRecord::MultistageCompressor {
                stages,
                outlet_pressure_pa,
                efficiency,
                intercooler_temperature_k,
            } => Self::MultistageCompressor(MultistageCompressorSpec {
                stages: count_in(stages, "stages", 1..=MAX_COMPRESSOR_STAGES)?,
                outlet_pressure_pa,
                efficiency,
                intercooler_temperature_k,
            }),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multistage_record_round_trips_and_loads_legacy_empty_records() {
        let spec = RotatingSpec::MultistageCompressor(MultistageCompressorSpec {
            stages: 3,
            outlet_pressure_pa: Some(6e6),
            efficiency: 0.8,
            intercooler_temperature_k: Some(310.0),
        });
        let json = serde_json::to_string(&RotatingRecord::from(&spec)).unwrap();
        let back: RotatingRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(RotatingSpec::try_from(back).unwrap(), spec);

        let legacy: RotatingRecord =
            serde_json::from_str(r#"{"kind": "multistage_compressor"}"#).unwrap();
        assert_eq!(
            RotatingSpec::try_from(legacy).unwrap(),
            RotatingSpec::MultistageCompressor(MultistageCompressorSpec::default())
        );
        let bad: RotatingRecord =
            serde_json::from_str(r#"{"kind": "multistage_compressor", "stages": 0}"#).unwrap();
        assert!(RotatingSpec::try_from(bad).is_err());
    }
}
