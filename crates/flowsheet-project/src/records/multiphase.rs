//! Project records of the multiphase group (#176, #177, #188, #189); owned
//! by the multiphase work package. Add fields with `#[serde(default)]`
//! behind format version 4.

use flowsheet_core::units::MultiphaseSpec;
use serde::{Deserialize, Serialize};

/// Multiphase specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum MultiphaseRecord {
    ThreePhaseSeparator {},
    Decanter {},
    ThreePhaseColumn {},
    LiquidExtractor {},
}

impl MultiphaseRecord {
    /// Every number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&MultiphaseSpec> for MultiphaseRecord {
    fn from(s: &MultiphaseSpec) -> Self {
        match s {
            MultiphaseSpec::ThreePhaseSeparator(_) => Self::ThreePhaseSeparator {},
            MultiphaseSpec::Decanter(_) => Self::Decanter {},
            MultiphaseSpec::ThreePhaseColumn(_) => Self::ThreePhaseColumn {},
            MultiphaseSpec::LiquidExtractor(_) => Self::LiquidExtractor {},
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
            MultiphaseRecord::LiquidExtractor {} => Self::LiquidExtractor(Default::default()),
        })
    }
}
