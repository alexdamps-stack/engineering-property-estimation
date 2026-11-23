//! Project records of the simple/logical group (#172–#175, #178,
//! #212–#218); owned by the simple/logical work package. Add fields with
//! `#[serde(default)]` behind format version 4.

use std::collections::BTreeMap;

use flowsheet_core::SeparatorThermal;
use flowsheet_core::units::SimpleLogicalSpec;
use flowsheet_core::units::simple_logical::{
    ComponentSplitterSpec, MAX_SELECTOR_INLETS, MultiplierSpec, SelectorSpec, TankSpec,
};
use serde::{Deserialize, Serialize};

use super::{ThermalRecord, count_in};

fn default_inlets() -> u16 {
    SelectorSpec::default().inlets
}

fn adiabatic() -> ThermalRecord {
    ThermalRecord::Adiabatic
}

/// Simple/logical specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum SimpleLogicalRecord {
    ComponentSplitter {
        #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
        overhead_recovery: BTreeMap<String, f64>,
    },
    Multiplier {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        factor: Option<f64>,
    },
    Selector {
        #[serde(default = "default_inlets")]
        inlets: u16,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        selected_inlet: Option<u16>,
    },
    StreamCutter {},
    Tank {
        #[serde(default = "adiabatic")]
        thermal: ThermalRecord,
        #[serde(default)]
        pressure_drop_pa: f64,
    },
    Adjust {},
    Set {},
    Spreadsheet {},
    Balance {},
    RecycleBlock {},
    ControlBlock {},
    UserUnit {},
}

impl SimpleLogicalRecord {
    /// Every number set in the record (counts are integers and always finite).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        match self {
            Self::ComponentSplitter { overhead_recovery } => {
                overhead_recovery.values().copied().collect()
            }
            Self::Multiplier { factor } => factor.iter().copied().collect(),
            Self::Tank {
                thermal,
                pressure_drop_pa,
            } => {
                let set_value = match *thermal {
                    ThermalRecord::OutletTemperatureK(t) | ThermalRecord::DutyW(t) => t,
                    _ => None,
                };
                std::iter::once(*pressure_drop_pa)
                    .chain(set_value)
                    .collect()
            }
            _ => Vec::new(),
        }
    }
}

impl From<&SimpleLogicalSpec> for SimpleLogicalRecord {
    fn from(s: &SimpleLogicalSpec) -> Self {
        use SimpleLogicalSpec as S;
        match s {
            S::ComponentSplitter(s) => Self::ComponentSplitter {
                overhead_recovery: s.overhead_recovery.clone(),
            },
            S::Multiplier(s) => Self::Multiplier { factor: s.factor },
            S::Selector(s) => Self::Selector {
                inlets: s.inlets,
                selected_inlet: s.selected_inlet,
            },
            S::StreamCutter(_) => Self::StreamCutter {},
            S::Tank(s) => Self::Tank {
                thermal: match s.thermal {
                    SeparatorThermal::Adiabatic => ThermalRecord::Adiabatic,
                    SeparatorThermal::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
                    SeparatorThermal::Duty(q) => ThermalRecord::DutyW(q),
                },
                pressure_drop_pa: s.pressure_drop_pa,
            },
            S::Adjust(_) => Self::Adjust {},
            S::Set(_) => Self::Set {},
            S::Spreadsheet(_) => Self::Spreadsheet {},
            S::Balance(_) => Self::Balance {},
            S::RecycleBlock(_) => Self::RecycleBlock {},
            S::ControlBlock(_) => Self::ControlBlock {},
            S::UserUnit(_) => Self::UserUnit {},
        }
    }
}

impl TryFrom<SimpleLogicalRecord> for SimpleLogicalSpec {
    type Error = String;

    fn try_from(r: SimpleLogicalRecord) -> Result<Self, String> {
        use SimpleLogicalRecord as R;
        Ok(match r {
            R::ComponentSplitter { overhead_recovery } => {
                Self::ComponentSplitter(ComponentSplitterSpec { overhead_recovery })
            }
            R::Multiplier { factor } => Self::Multiplier(MultiplierSpec { factor }),
            R::Selector {
                inlets,
                selected_inlet,
            } => Self::Selector(SelectorSpec {
                inlets: count_in(inlets, "inlets", 2..=MAX_SELECTOR_INLETS)?,
                selected_inlet,
            }),
            R::StreamCutter {} => Self::StreamCutter(Default::default()),
            R::Tank {
                thermal,
                pressure_drop_pa,
            } => Self::Tank(TankSpec {
                thermal: match thermal {
                    ThermalRecord::Adiabatic => SeparatorThermal::Adiabatic,
                    ThermalRecord::OutletTemperatureK(t) => SeparatorThermal::OutletTemperature(t),
                    ThermalRecord::DutyW(q) => SeparatorThermal::Duty(q),
                    ThermalRecord::EnergyStream => {
                        return Err(
                            "tank thermal mode cannot take its duty from an energy stream".into(),
                        );
                    }
                },
                pressure_drop_pa,
            }),
            R::Adjust {} => Self::Adjust(Default::default()),
            R::Set {} => Self::Set(Default::default()),
            R::Spreadsheet {} => Self::Spreadsheet(Default::default()),
            R::Balance {} => Self::Balance(Default::default()),
            R::RecycleBlock {} => Self::RecycleBlock(Default::default()),
            R::ControlBlock {} => Self::ControlBlock(Default::default()),
            R::UserUnit {} => Self::UserUnit(Default::default()),
        })
    }
}
