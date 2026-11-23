//! Project records of the rotating group (#182, A14); owned by the rotating
//! work package. Pump, compressor and expander keep the core `pressurize` /
//! `expand` records in `records/mod.rs`. Add fields with `#[serde(default)]`
//! behind format version 4.

use flowsheet_core::units::RotatingSpec;
use serde::{Deserialize, Serialize};

/// Rotating-group specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum RotatingRecord {
    MultistageCompressor {},
}

impl RotatingRecord {
    /// Every number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&RotatingSpec> for RotatingRecord {
    fn from(s: &RotatingSpec) -> Self {
        match s {
            RotatingSpec::MultistageCompressor(_) => Self::MultistageCompressor {},
        }
    }
}

impl TryFrom<RotatingRecord> for RotatingSpec {
    type Error = String;

    fn try_from(r: RotatingRecord) -> Result<Self, String> {
        Ok(match r {
            RotatingRecord::MultistageCompressor {} => {
                Self::MultistageCompressor(Default::default())
            }
        })
    }
}
