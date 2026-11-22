//! Project records of the reactors group (#191–#197); owned by the reactors
//! work package. Add fields with `#[serde(default)]` behind format version 4.

use flowsheet_core::units::ReactorSpec;
use serde::{Deserialize, Serialize};

/// Reactor specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ReactorRecord {
    ConversionReactor {},
    YieldReactor {},
    EquilibriumReactor {},
    GibbsReactor {},
    Cstr {},
    PlugFlowReactor {},
    BatchReactor {},
}

impl ReactorRecord {
    /// Every number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        Vec::new()
    }
}

impl From<&ReactorSpec> for ReactorRecord {
    fn from(s: &ReactorSpec) -> Self {
        match s {
            ReactorSpec::Conversion(_) => Self::ConversionReactor {},
            ReactorSpec::Yield(_) => Self::YieldReactor {},
            ReactorSpec::Equilibrium(_) => Self::EquilibriumReactor {},
            ReactorSpec::Gibbs(_) => Self::GibbsReactor {},
            ReactorSpec::Cstr(_) => Self::Cstr {},
            ReactorSpec::PlugFlow(_) => Self::PlugFlowReactor {},
            ReactorSpec::Batch(_) => Self::BatchReactor {},
        }
    }
}

impl TryFrom<ReactorRecord> for ReactorSpec {
    type Error = String;

    fn try_from(r: ReactorRecord) -> Result<Self, String> {
        Ok(match r {
            ReactorRecord::ConversionReactor {} => Self::Conversion(Default::default()),
            ReactorRecord::YieldReactor {} => Self::Yield(Default::default()),
            ReactorRecord::EquilibriumReactor {} => Self::Equilibrium(Default::default()),
            ReactorRecord::GibbsReactor {} => Self::Gibbs(Default::default()),
            ReactorRecord::Cstr {} => Self::Cstr(Default::default()),
            ReactorRecord::PlugFlowReactor {} => Self::PlugFlow(Default::default()),
            ReactorRecord::BatchReactor {} => Self::Batch(Default::default()),
        })
    }
}
