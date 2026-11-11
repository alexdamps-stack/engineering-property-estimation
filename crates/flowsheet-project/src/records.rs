//! Serializable mirrors of the `flowsheet-core` topology and specifications.
//!
//! The engine types stay free of serialization concerns; these records define
//! the stable JSON spelling shared by project format versions 1 and 2.

use flowsheet_core::{
    ExchangerSpec, MixerPressure, PressureRise, PressureSpec, SeparatorThermal, ThermalSpec,
    UnitKind, UnitSpec,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum KindRecord {
    Feed,
    Product,
    Mixer,
    Heater,
    Cooler,
    Valve,
    Separator,
    Splitter,
    HeatExchanger,
    Pump,
    Compressor,
    Expander,
}

impl From<UnitKind> for KindRecord {
    fn from(k: UnitKind) -> Self {
        match k {
            UnitKind::Feed => KindRecord::Feed,
            UnitKind::Product => KindRecord::Product,
            UnitKind::Mixer => KindRecord::Mixer,
            UnitKind::Heater => KindRecord::Heater,
            UnitKind::Cooler => KindRecord::Cooler,
            UnitKind::Valve => KindRecord::Valve,
            UnitKind::Separator => KindRecord::Separator,
            UnitKind::Splitter => KindRecord::Splitter,
            UnitKind::HeatExchanger => KindRecord::HeatExchanger,
            UnitKind::Pump => KindRecord::Pump,
            UnitKind::Compressor => KindRecord::Compressor,
            UnitKind::Expander => KindRecord::Expander,
        }
    }
}

impl From<KindRecord> for UnitKind {
    fn from(k: KindRecord) -> Self {
        match k {
            KindRecord::Feed => UnitKind::Feed,
            KindRecord::Product => UnitKind::Product,
            KindRecord::Mixer => UnitKind::Mixer,
            KindRecord::Heater => UnitKind::Heater,
            KindRecord::Cooler => UnitKind::Cooler,
            KindRecord::Valve => UnitKind::Valve,
            KindRecord::Separator => UnitKind::Separator,
            KindRecord::Splitter => UnitKind::Splitter,
            KindRecord::HeatExchanger => UnitKind::HeatExchanger,
            KindRecord::Pump => UnitKind::Pump,
            KindRecord::Compressor => UnitKind::Compressor,
            KindRecord::Expander => UnitKind::Expander,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum MixerPressureRecord {
    LowestInlet,
    SpecifiedPa(Option<f64>),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ThermalRecord {
    Adiabatic,
    OutletTemperatureK(Option<f64>),
    DutyW(Option<f64>),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
// The `_pa` suffix names the SI unit in the stored JSON spelling.
#[allow(clippy::enum_variant_names)]
pub(crate) enum PressureRecord {
    OutletPressurePa(Option<f64>),
    PressureDropPa(Option<f64>),
    PressureIncreasePa(Option<f64>),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum ExchangerRecord {
    HotOutletTemperatureK(Option<f64>),
    ColdOutletTemperatureK(Option<f64>),
    DutyW(Option<f64>),
}

/// Specification in SI units with absolute pressure; `null` means unset.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum SpecRecord {
    Feed {
        temperature_k: Option<f64>,
        pressure_pa: Option<f64>,
        molar_flow_mol_s: Option<f64>,
    },
    Product,
    Mixer {
        outlet_pressure: MixerPressureRecord,
    },
    HeatExchange {
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    Valve {
        pressure: PressureRecord,
    },
    Separator {
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    Splitter {
        fraction_to_first: Option<f64>,
    },
    HeatExchanger {
        thermal: ExchangerRecord,
        hot_pressure_drop_pa: f64,
        cold_pressure_drop_pa: f64,
    },
    Pressurize {
        pressure: PressureRecord,
        efficiency: f64,
    },
    Expand {
        pressure: PressureRecord,
        efficiency: f64,
    },
}

impl From<UnitSpec> for SpecRecord {
    fn from(s: UnitSpec) -> Self {
        let drop_or_outlet = |p: PressureSpec| match p {
            PressureSpec::OutletPressure(v) => PressureRecord::OutletPressurePa(v),
            PressureSpec::PressureDrop(v) => PressureRecord::PressureDropPa(v),
        };
        match s {
            UnitSpec::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            } => SpecRecord::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            },
            UnitSpec::Product => SpecRecord::Product,
            UnitSpec::Mixer { outlet_pressure } => SpecRecord::Mixer {
                outlet_pressure: match outlet_pressure {
                    MixerPressure::LowestInlet => MixerPressureRecord::LowestInlet,
                    MixerPressure::Specified(p) => MixerPressureRecord::SpecifiedPa(p),
                },
            },
            UnitSpec::HeatExchange {
                thermal,
                pressure_drop_pa,
            } => SpecRecord::HeatExchange {
                thermal: match thermal {
                    ThermalSpec::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
                    ThermalSpec::Duty(q) => ThermalRecord::DutyW(q),
                },
                pressure_drop_pa,
            },
            UnitSpec::Valve { pressure } => SpecRecord::Valve {
                pressure: drop_or_outlet(pressure),
            },
            UnitSpec::Separator {
                thermal,
                pressure_drop_pa,
            } => SpecRecord::Separator {
                thermal: match thermal {
                    SeparatorThermal::Adiabatic => ThermalRecord::Adiabatic,
                    SeparatorThermal::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
                    SeparatorThermal::Duty(q) => ThermalRecord::DutyW(q),
                },
                pressure_drop_pa,
            },
            UnitSpec::Splitter { fraction_to_first } => SpecRecord::Splitter { fraction_to_first },
            UnitSpec::HeatExchanger {
                thermal,
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
            } => SpecRecord::HeatExchanger {
                thermal: match thermal {
                    ExchangerSpec::HotOutletTemperature(t) => {
                        ExchangerRecord::HotOutletTemperatureK(t)
                    }
                    ExchangerSpec::ColdOutletTemperature(t) => {
                        ExchangerRecord::ColdOutletTemperatureK(t)
                    }
                    ExchangerSpec::Duty(q) => ExchangerRecord::DutyW(q),
                },
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
            },
            UnitSpec::Pressurize {
                pressure,
                efficiency,
            } => SpecRecord::Pressurize {
                pressure: match pressure {
                    PressureRise::OutletPressure(p) => PressureRecord::OutletPressurePa(p),
                    PressureRise::PressureIncrease(p) => PressureRecord::PressureIncreasePa(p),
                },
                efficiency,
            },
            UnitSpec::Expand {
                pressure,
                efficiency,
            } => SpecRecord::Expand {
                pressure: drop_or_outlet(pressure),
                efficiency,
            },
        }
    }
}

impl TryFrom<SpecRecord> for UnitSpec {
    type Error = String;

    fn try_from(s: SpecRecord) -> Result<Self, String> {
        let unexpected = |what: &str| Err(format!("{what} is not valid for this specification"));
        let pressure = |p: PressureRecord| match p {
            PressureRecord::OutletPressurePa(v) => Ok(PressureSpec::OutletPressure(v)),
            PressureRecord::PressureDropPa(v) => Ok(PressureSpec::PressureDrop(v)),
            PressureRecord::PressureIncreasePa(_) => {
                Err("pressure_increase_pa is only valid for pumps and compressors".to_string())
            }
        };
        Ok(match s {
            SpecRecord::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            } => UnitSpec::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            },
            SpecRecord::Product => UnitSpec::Product,
            SpecRecord::Mixer { outlet_pressure } => UnitSpec::Mixer {
                outlet_pressure: match outlet_pressure {
                    MixerPressureRecord::LowestInlet => MixerPressure::LowestInlet,
                    MixerPressureRecord::SpecifiedPa(p) => MixerPressure::Specified(p),
                },
            },
            SpecRecord::HeatExchange {
                thermal,
                pressure_drop_pa,
            } => UnitSpec::HeatExchange {
                thermal: match thermal {
                    ThermalRecord::OutletTemperatureK(t) => ThermalSpec::OutletTemperature(t),
                    ThermalRecord::DutyW(q) => ThermalSpec::Duty(q),
                    ThermalRecord::Adiabatic => return unexpected("adiabatic"),
                },
                pressure_drop_pa,
            },
            SpecRecord::Valve { pressure: p } => UnitSpec::Valve {
                pressure: pressure(p)?,
            },
            SpecRecord::Separator {
                thermal,
                pressure_drop_pa,
            } => UnitSpec::Separator {
                thermal: match thermal {
                    ThermalRecord::Adiabatic => SeparatorThermal::Adiabatic,
                    ThermalRecord::OutletTemperatureK(t) => SeparatorThermal::OutletTemperature(t),
                    ThermalRecord::DutyW(q) => SeparatorThermal::Duty(q),
                },
                pressure_drop_pa,
            },
            SpecRecord::Splitter { fraction_to_first } => UnitSpec::Splitter { fraction_to_first },
            SpecRecord::HeatExchanger {
                thermal,
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
            } => UnitSpec::HeatExchanger {
                thermal: match thermal {
                    ExchangerRecord::HotOutletTemperatureK(t) => {
                        ExchangerSpec::HotOutletTemperature(t)
                    }
                    ExchangerRecord::ColdOutletTemperatureK(t) => {
                        ExchangerSpec::ColdOutletTemperature(t)
                    }
                    ExchangerRecord::DutyW(q) => ExchangerSpec::Duty(q),
                },
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
            },
            SpecRecord::Pressurize {
                pressure,
                efficiency,
            } => UnitSpec::Pressurize {
                pressure: match pressure {
                    PressureRecord::OutletPressurePa(p) => PressureRise::OutletPressure(p),
                    PressureRecord::PressureIncreasePa(p) => PressureRise::PressureIncrease(p),
                    PressureRecord::PressureDropPa(_) => return unexpected("pressure_drop_pa"),
                },
                efficiency,
            },
            SpecRecord::Expand {
                pressure: p,
                efficiency,
            } => UnitSpec::Expand {
                pressure: pressure(p)?,
                efficiency,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_default_spec_round_trips() {
        for kind in UnitKind::ALL {
            let spec = kind.default_spec();
            let json = serde_json::to_string(&SpecRecord::from(spec)).unwrap();
            let back: SpecRecord = serde_json::from_str(&json).unwrap();
            assert_eq!(UnitSpec::try_from(back).unwrap(), spec, "{kind}");
            let k: KindRecord = kind.into();
            assert_eq!(UnitKind::from(k), kind);
        }
    }

    #[test]
    fn mismatched_variants_are_rejected() {
        let bad = SpecRecord::HeatExchange {
            thermal: ThermalRecord::Adiabatic,
            pressure_drop_pa: 0.0,
        };
        assert!(UnitSpec::try_from(bad).is_err());
        let bad = SpecRecord::Valve {
            pressure: PressureRecord::PressureIncreasePa(Some(1.0)),
        };
        assert!(UnitSpec::try_from(bad).is_err());
    }
}
