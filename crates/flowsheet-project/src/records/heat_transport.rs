//! Project records of the heat-transport group (#179–#181, #183, #184);
//! owned by the heat-transport work package. Add fields with
//! `#[serde(default)]` behind format version 4.

use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{MAX_PASSES, MultiStreamExchangerSpec};
use serde::{Deserialize, Serialize};

use super::count_in;

fn default_passes() -> u16 {
    MultiStreamExchangerSpec::default().passes
}

/// Heat-transport specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum HeatTransportRecord {
    AirCooler {},
    FiredHeater {},
    MultiStreamExchanger {
        #[serde(default = "default_passes")]
        passes: u16,
    },
    PipeSegment {},
    ReliefValve {},
}

impl HeatTransportRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&HeatTransportSpec> for HeatTransportRecord {
    fn from(s: &HeatTransportSpec) -> Self {
        match s {
            HeatTransportSpec::AirCooler(_) => Self::AirCooler {},
            HeatTransportSpec::FiredHeater(_) => Self::FiredHeater {},
            HeatTransportSpec::MultiStreamExchanger(x) => {
                Self::MultiStreamExchanger { passes: x.passes }
            }
            HeatTransportSpec::PipeSegment(_) => Self::PipeSegment {},
            HeatTransportSpec::ReliefValve(_) => Self::ReliefValve {},
        }
    }
}

impl TryFrom<HeatTransportRecord> for HeatTransportSpec {
    type Error = String;

    fn try_from(r: HeatTransportRecord) -> Result<Self, String> {
        Ok(match r {
            HeatTransportRecord::AirCooler {} => Self::AirCooler(Default::default()),
            HeatTransportRecord::FiredHeater {} => Self::FiredHeater(Default::default()),
            HeatTransportRecord::MultiStreamExchanger { passes } => {
                Self::MultiStreamExchanger(MultiStreamExchangerSpec {
                    passes: count_in(passes, "passes", 2..=MAX_PASSES)?,
                })
            }
            HeatTransportRecord::PipeSegment {} => Self::PipeSegment(Default::default()),
            HeatTransportRecord::ReliefValve {} => Self::ReliefValve(Default::default()),
        })
    }
}
