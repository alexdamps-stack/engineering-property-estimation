//! Project records of the heat-transport group (#179–#181, #183, #184);
//! owned by the heat-transport work package. Add fields with
//! `#[serde(default)]` behind format version 4.

use flowsheet_core::ThermalSpec;
use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{
    AirCoolerSpec, AirCoolerThermal, FiredHeaterSpec, MAX_PASSES, MultiStreamExchangerSpec,
};
use serde::{Deserialize, Serialize};

use super::{ThermalRecord, count_in};

fn default_passes() -> u16 {
    MultiStreamExchangerSpec::default().passes
}

fn unset_temperature() -> ThermalRecord {
    ThermalRecord::OutletTemperatureK(None)
}

/// Heat-transport specification, tagged by unit kind.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum HeatTransportRecord {
    AirCooler {
        #[serde(default = "unset_temperature")]
        thermal: ThermalRecord,
        #[serde(default)]
        pressure_drop_pa: f64,
        #[serde(default)]
        ambient_temperature_k: Option<f64>,
        #[serde(default)]
        air_mass_flow_kg_s: Option<f64>,
        #[serde(default)]
        fan_power_w: f64,
    },
    FiredHeater {
        #[serde(default = "unset_temperature")]
        thermal: ThermalRecord,
        #[serde(default)]
        pressure_drop_pa: f64,
    },
    MultiStreamExchanger {
        #[serde(default = "default_passes")]
        passes: u16,
    },
    PipeSegment {},
    ReliefValve {},
}

impl HeatTransportRecord {
    /// Every floating-point number set in the record.
    pub(crate) fn numbers(&self) -> Vec<f64> {
        let thermal_value = |thermal: &ThermalRecord| match thermal {
            ThermalRecord::OutletTemperatureK(t) | ThermalRecord::DutyW(t) => *t,
            ThermalRecord::Adiabatic | ThermalRecord::EnergyStream => None,
        };
        match self {
            Self::AirCooler {
                thermal,
                pressure_drop_pa,
                ambient_temperature_k,
                air_mass_flow_kg_s,
                fan_power_w,
            } => [
                thermal_value(thermal),
                Some(*pressure_drop_pa),
                *ambient_temperature_k,
                *air_mass_flow_kg_s,
                Some(*fan_power_w),
            ]
            .into_iter()
            .flatten()
            .collect(),
            Self::FiredHeater {
                thermal,
                pressure_drop_pa,
            } => [Some(*pressure_drop_pa), thermal_value(thermal)]
                .into_iter()
                .flatten()
                .collect(),
            _ => Vec::new(),
        }
    }
}

impl From<&HeatTransportSpec> for HeatTransportRecord {
    fn from(s: &HeatTransportSpec) -> Self {
        match s {
            HeatTransportSpec::AirCooler(s) => Self::AirCooler {
                thermal: match s.thermal {
                    AirCoolerThermal::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
                    AirCoolerThermal::Duty(q) => ThermalRecord::DutyW(q),
                },
                pressure_drop_pa: s.pressure_drop_pa,
                ambient_temperature_k: s.ambient_temperature_k,
                air_mass_flow_kg_s: s.air_mass_flow_kg_s,
                fan_power_w: s.fan_power_w,
            },
            HeatTransportSpec::FiredHeater(s) => Self::FiredHeater {
                thermal: match s.thermal {
                    ThermalSpec::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
                    ThermalSpec::Duty(q) => ThermalRecord::DutyW(q),
                    ThermalSpec::EnergyStream => ThermalRecord::EnergyStream,
                },
                pressure_drop_pa: s.pressure_drop_pa,
            },
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
            HeatTransportRecord::AirCooler {
                thermal,
                pressure_drop_pa,
                ambient_temperature_k,
                air_mass_flow_kg_s,
                fan_power_w,
            } => Self::AirCooler(AirCoolerSpec {
                thermal: match thermal {
                    ThermalRecord::OutletTemperatureK(t) => AirCoolerThermal::OutletTemperature(t),
                    ThermalRecord::DutyW(q) => AirCoolerThermal::Duty(q),
                    _ => return Err("unexpected thermal specification for an air cooler".into()),
                },
                pressure_drop_pa,
                ambient_temperature_k,
                air_mass_flow_kg_s,
                fan_power_w,
            }),
            HeatTransportRecord::FiredHeater {
                thermal,
                pressure_drop_pa,
            } => Self::FiredHeater(FiredHeaterSpec {
                thermal: match thermal {
                    ThermalRecord::OutletTemperatureK(t) => ThermalSpec::OutletTemperature(t),
                    ThermalRecord::DutyW(q) => ThermalSpec::Duty(q),
                    ThermalRecord::EnergyStream => ThermalSpec::EnergyStream,
                    ThermalRecord::Adiabatic => {
                        return Err("unexpected adiabatic specification for a fired heater".into());
                    }
                },
                pressure_drop_pa,
            }),
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
