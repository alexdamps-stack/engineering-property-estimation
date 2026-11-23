//! Project records of the solids group (#198–#211); owned by the solids work
//! package. Add fields with `#[serde(default)]` behind format version 4.

use flowsheet_core::units::SolidsSpec;
use serde::{Deserialize, Serialize};

/// Generates the record enum and its conversions for kinds whose
/// specification has no fields yet. Move a kind out of the macro into
/// explicit code when it gains fields.
macro_rules! solids_records {
    ($($record:ident => $spec:ident),* $(,)?) => {
        /// Solids specification, tagged by unit kind.
        #[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
        #[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
        pub(crate) enum SolidsRecord {
            $($record {},)*
        }

        impl From<&SolidsSpec> for SolidsRecord {
            fn from(s: &SolidsSpec) -> Self {
                match s {
                    $(SolidsSpec::$spec(_) => Self::$record {},)*
                }
            }
        }

        impl TryFrom<SolidsRecord> for SolidsSpec {
            type Error = String;

            fn try_from(r: SolidsRecord) -> Result<Self, String> {
                Ok(match r {
                    $(SolidsRecord::$record {} => Self::$spec(Default::default()),)*
                })
            }
        }
    };
}

solids_records! {
    Cyclone => Cyclone,
    Hydrocyclone => Hydrocyclone,
    Baghouse => Baghouse,
    RotaryFilter => RotaryFilter,
    SolidsSeparator => SolidsSeparator,
    Crystallizer => Crystallizer,
    Centrifuge => Centrifuge,
    ElectrostaticPrecipitator => ElectrostaticPrecipitator,
    VenturiScrubber => VenturiScrubber,
    Crusher => Crusher,
    Screen => Screen,
    Dryer => Dryer,
    Granulator => Granulator,
    SolidsWasher => SolidsWasher,
}

impl SolidsRecord {
    /// Every number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}
