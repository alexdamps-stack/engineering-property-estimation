//! Project records of the heat-transport group (#179–#181, #183, #184);
//! owned by the heat-transport work package. Add fields with
//! `#[serde(default)]` behind format version 4.

use flowsheet_core::ThermalSpec;
use flowsheet_core::units::HeatTransportSpec;
use flowsheet_core::units::heat_transport::{
    AirCoolerSpec, AirCoolerThermal, ExchangerPassDirection, ExchangerPassSpec,
    ExchangerPassThermal, FiredHeaterSpec, MAX_PASSES, MultiStreamExchangerSpec,
};
use serde::{Deserialize, Serialize};

use super::{ThermalRecord, count_in};

fn default_passes() -> u16 {
    MultiStreamExchangerSpec::default().passes.len() as u16
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
        #[serde(default)]
        pass_specs: Option<Vec<ExchangerPassRecord>>,
        #[serde(default)]
        minimum_approach_k: f64,
    },
    PipeSegment {},
    ReliefValve {},
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ExchangerPassRecord {
    direction: ExchangerDirectionRecord,
    thermal: ExchangerPassThermalRecord,
    #[serde(default)]
    pressure_drop_pa: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum ExchangerDirectionRecord {
    Hot,
    Cold,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
enum ExchangerPassThermalRecord {
    OutletTemperatureK(Option<f64>),
    DutyW(Option<f64>),
    TerminalApproach {
        reference_pass: u16,
        approach_k: Option<f64>,
    },
    Balance,
}

impl From<&ExchangerPassSpec> for ExchangerPassRecord {
    fn from(pass: &ExchangerPassSpec) -> Self {
        Self {
            direction: match pass.direction {
                ExchangerPassDirection::Hot => ExchangerDirectionRecord::Hot,
                ExchangerPassDirection::Cold => ExchangerDirectionRecord::Cold,
            },
            thermal: match pass.thermal {
                ExchangerPassThermal::OutletTemperature(t) => {
                    ExchangerPassThermalRecord::OutletTemperatureK(t)
                }
                ExchangerPassThermal::Duty(q) => ExchangerPassThermalRecord::DutyW(q),
                ExchangerPassThermal::TerminalApproach {
                    reference_pass,
                    approach_k,
                } => ExchangerPassThermalRecord::TerminalApproach {
                    reference_pass,
                    approach_k,
                },
                ExchangerPassThermal::Balance => ExchangerPassThermalRecord::Balance,
            },
            pressure_drop_pa: pass.pressure_drop_pa,
        }
    }
}

impl From<ExchangerPassRecord> for ExchangerPassSpec {
    fn from(pass: ExchangerPassRecord) -> Self {
        Self {
            direction: match pass.direction {
                ExchangerDirectionRecord::Hot => ExchangerPassDirection::Hot,
                ExchangerDirectionRecord::Cold => ExchangerPassDirection::Cold,
            },
            thermal: match pass.thermal {
                ExchangerPassThermalRecord::OutletTemperatureK(t) => {
                    ExchangerPassThermal::OutletTemperature(t)
                }
                ExchangerPassThermalRecord::DutyW(q) => ExchangerPassThermal::Duty(q),
                ExchangerPassThermalRecord::TerminalApproach {
                    reference_pass,
                    approach_k,
                } => ExchangerPassThermal::TerminalApproach {
                    reference_pass,
                    approach_k,
                },
                ExchangerPassThermalRecord::Balance => ExchangerPassThermal::Balance,
            },
            pressure_drop_pa: pass.pressure_drop_pa,
        }
    }
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
            Self::MultiStreamExchanger {
                pass_specs,
                minimum_approach_k,
                ..
            } => std::iter::once(*minimum_approach_k)
                .chain(pass_specs.iter().flatten().flat_map(|pass| {
                    let value = match pass.thermal {
                        ExchangerPassThermalRecord::OutletTemperatureK(t)
                        | ExchangerPassThermalRecord::DutyW(t) => t,
                        ExchangerPassThermalRecord::TerminalApproach { approach_k, .. } => {
                            approach_k
                        }
                        ExchangerPassThermalRecord::Balance => None,
                    };
                    [Some(pass.pressure_drop_pa), value].into_iter().flatten()
                }))
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
            HeatTransportSpec::MultiStreamExchanger(x) => Self::MultiStreamExchanger {
                passes: u16::try_from(x.passes.len()).unwrap_or(u16::MAX),
                pass_specs: Some(x.passes.iter().map(Into::into).collect()),
                minimum_approach_k: x.minimum_approach_k,
            },
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
            HeatTransportRecord::MultiStreamExchanger {
                passes,
                pass_specs,
                minimum_approach_k,
            } => {
                let count = count_in(passes, "passes", 2..=MAX_PASSES)?;
                let passes = match pass_specs {
                    Some(specs) => {
                        if specs.len() != usize::from(count) {
                            return Err(
                                "multi-stream exchanger pass specs do not match the pass count"
                                    .into(),
                            );
                        }
                        specs.into_iter().map(Into::into).collect()
                    }
                    None => MultiStreamExchangerSpec::with_pass_count(count).passes,
                };
                Self::MultiStreamExchanger(MultiStreamExchangerSpec {
                    passes,
                    minimum_approach_k,
                })
            }
            HeatTransportRecord::PipeSegment {} => Self::PipeSegment(Default::default()),
            HeatTransportRecord::ReliefValve {} => Self::ReliefValve(Default::default()),
        })
    }
}
