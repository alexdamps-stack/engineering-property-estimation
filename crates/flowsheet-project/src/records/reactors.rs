//! Legacy drawn-only reactor records from the A17 palette.

use flowsheet_core::{UnitKind, UnitSpec};
use serde::{Deserialize, Serialize};

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
    pub(crate) fn default_spec(self) -> UnitSpec {
        let kind = match self {
            Self::ConversionReactor {} => UnitKind::ConversionReactor,
            Self::YieldReactor {} => UnitKind::YieldReactor,
            Self::EquilibriumReactor {} => UnitKind::EquilibriumReactor,
            Self::GibbsReactor {} => UnitKind::GibbsReactor,
            Self::Cstr {} => UnitKind::Cstr,
            Self::PlugFlowReactor {} => UnitKind::PlugFlowReactor,
            Self::BatchReactor {} => UnitKind::BatchReactor,
        };
        kind.default_spec()
    }
}
