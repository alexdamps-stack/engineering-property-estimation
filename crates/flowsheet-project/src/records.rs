//! Serializable mirrors of the `flowsheet-core` topology and specifications.
//!
//! The engine types stay free of serialization concerns; these records define
//! the stable JSON spelling of project format versions 1 to 4. Version 4 adds
//! the reactor kinds and their specifications; earlier versions share the
//! remaining records unchanged.

use flowsheet_core::{
    ExchangerSpec, MixerPressure, PressureRise, PressureSpec, Reaction, ReactionExtent,
    SeparatorThermal, ThermalSpec, UnitKind, UnitSpec,
};
use serde::{Deserialize, Serialize};

use crate::reactor_records::{
    BatchThermalRecord, EquilibriumOptionsRecord, EquilibriumReactionRecord, KineticsRecord,
    PlugFlowThermalRecord, equilibrium_numbers, kinetics_numbers,
};

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
    ConversionReactor,
    YieldReactor,
    EquilibriumReactor,
    GibbsReactor,
    Cstr,
    PlugFlowReactor,
    BatchReactor,
}

impl KindRecord {
    /// First project format version that can store this kind.
    pub(crate) fn since_format(self) -> u64 {
        match self {
            KindRecord::ConversionReactor
            | KindRecord::YieldReactor
            | KindRecord::EquilibriumReactor
            | KindRecord::GibbsReactor
            | KindRecord::Cstr
            | KindRecord::PlugFlowReactor
            | KindRecord::BatchReactor => 4,
            _ => 1,
        }
    }
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
            UnitKind::ConversionReactor => KindRecord::ConversionReactor,
            UnitKind::YieldReactor => KindRecord::YieldReactor,
            UnitKind::EquilibriumReactor => KindRecord::EquilibriumReactor,
            UnitKind::GibbsReactor => KindRecord::GibbsReactor,
            UnitKind::Cstr => KindRecord::Cstr,
            UnitKind::PlugFlowReactor => KindRecord::PlugFlowReactor,
            UnitKind::BatchReactor => KindRecord::BatchReactor,
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
            KindRecord::ConversionReactor => UnitKind::ConversionReactor,
            KindRecord::YieldReactor => UnitKind::YieldReactor,
            KindRecord::EquilibriumReactor => UnitKind::EquilibriumReactor,
            KindRecord::GibbsReactor => UnitKind::GibbsReactor,
            KindRecord::Cstr => UnitKind::Cstr,
            KindRecord::PlugFlowReactor => UnitKind::PlugFlowReactor,
            KindRecord::BatchReactor => UnitKind::BatchReactor,
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

/// One stoichiometric coefficient (format 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CoefficientRecord {
    pub component: String,
    pub coefficient: f64,
}

/// Extent of one reaction (format 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub(crate) enum ExtentRecord {
    Conversion {
        key: Option<String>,
        fraction: Option<f64>,
    },
    MolarExtentMolS(Option<f64>),
}

/// One conversion-reactor reaction (format 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ReactionRecord {
    pub coefficients: Vec<CoefficientRecord>,
    pub extent: ExtentRecord,
}

/// One relative mass yield (format 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct YieldRecord {
    pub component: String,
    pub mass_yield: f64,
}

/// Specification in SI units with absolute pressure; `null` means unset.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    ConversionReactor {
        reactions: Vec<ReactionRecord>,
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    YieldReactor {
        mass_yields: Vec<YieldRecord>,
        inerts: Vec<String>,
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    EquilibriumReactor {
        reactions: Vec<EquilibriumReactionRecord>,
        equilibrium: EquilibriumOptionsRecord,
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    GibbsReactor {
        inerts: Vec<String>,
        equilibrium: EquilibriumOptionsRecord,
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    Cstr {
        kinetics: KineticsRecord,
        volume_m3: Option<f64>,
        thermal: ThermalRecord,
        pressure_drop_pa: f64,
    },
    PlugFlow {
        kinetics: KineticsRecord,
        length_m: Option<f64>,
        diameter_m: Option<f64>,
        thermal: PlugFlowThermalRecord,
        pressure_drop_pa: f64,
    },
    Batch {
        kinetics: KineticsRecord,
        reaction_time_s: Option<f64>,
        thermal: BatchThermalRecord,
    },
}

impl SpecRecord {
    /// Every number set in the specification (unset values are skipped).
    pub(crate) fn numbers(&self) -> Vec<f64> {
        let pressure = |p: &PressureRecord| match *p {
            PressureRecord::OutletPressurePa(v)
            | PressureRecord::PressureDropPa(v)
            | PressureRecord::PressureIncreasePa(v) => v,
        };
        let thermal = |t: &ThermalRecord| match *t {
            ThermalRecord::Adiabatic => None,
            ThermalRecord::OutletTemperatureK(v) | ThermalRecord::DutyW(v) => v,
        };
        let values: Vec<Option<f64>> = match self {
            SpecRecord::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            } => vec![*temperature_k, *pressure_pa, *molar_flow_mol_s],
            SpecRecord::Product => vec![],
            SpecRecord::Mixer { outlet_pressure } => match outlet_pressure {
                MixerPressureRecord::LowestInlet => vec![],
                MixerPressureRecord::SpecifiedPa(p) => vec![*p],
            },
            SpecRecord::HeatExchange {
                thermal: t,
                pressure_drop_pa,
            }
            | SpecRecord::Separator {
                thermal: t,
                pressure_drop_pa,
            } => vec![thermal(t), Some(*pressure_drop_pa)],
            SpecRecord::Valve { pressure: p } => vec![pressure(p)],
            SpecRecord::Splitter { fraction_to_first } => vec![*fraction_to_first],
            SpecRecord::HeatExchanger {
                thermal: t,
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
            } => vec![
                match *t {
                    ExchangerRecord::HotOutletTemperatureK(v)
                    | ExchangerRecord::ColdOutletTemperatureK(v)
                    | ExchangerRecord::DutyW(v) => v,
                },
                Some(*hot_pressure_drop_pa),
                Some(*cold_pressure_drop_pa),
            ],
            SpecRecord::Pressurize {
                pressure: p,
                efficiency,
            }
            | SpecRecord::Expand {
                pressure: p,
                efficiency,
            } => vec![pressure(p), Some(*efficiency)],
            SpecRecord::ConversionReactor {
                reactions,
                thermal: t,
                pressure_drop_pa,
            } => {
                let mut v = vec![thermal(t), Some(*pressure_drop_pa)];
                for r in reactions {
                    v.extend(r.coefficients.iter().map(|c| Some(c.coefficient)));
                    v.push(match &r.extent {
                        ExtentRecord::Conversion { fraction, .. } => *fraction,
                        ExtentRecord::MolarExtentMolS(e) => *e,
                    });
                }
                v
            }
            SpecRecord::YieldReactor {
                mass_yields,
                thermal: t,
                pressure_drop_pa,
                ..
            } => {
                let mut v = vec![thermal(t), Some(*pressure_drop_pa)];
                v.extend(mass_yields.iter().map(|y| Some(y.mass_yield)));
                v
            }
            SpecRecord::EquilibriumReactor {
                reactions,
                equilibrium,
                thermal: t,
                pressure_drop_pa,
            } => {
                let mut v = vec![thermal(t), Some(*pressure_drop_pa)];
                v.extend(equilibrium_numbers(reactions, equilibrium));
                v
            }
            SpecRecord::GibbsReactor {
                equilibrium,
                thermal: t,
                pressure_drop_pa,
                ..
            } => {
                let mut v = vec![thermal(t), Some(*pressure_drop_pa)];
                v.extend(equilibrium_numbers(&[], equilibrium));
                v
            }
            SpecRecord::Cstr {
                kinetics,
                volume_m3,
                thermal: t,
                pressure_drop_pa,
            } => {
                let mut v = vec![thermal(t), Some(*pressure_drop_pa), *volume_m3];
                v.extend(kinetics_numbers(kinetics));
                v
            }
            SpecRecord::PlugFlow {
                kinetics,
                length_m,
                diameter_m,
                thermal: t,
                pressure_drop_pa,
            } => {
                let mut v = vec![*length_m, *diameter_m, Some(*pressure_drop_pa)];
                match *t {
                    PlugFlowThermalRecord::IsothermalK(x) => v.push(x),
                    PlugFlowThermalRecord::Adiabatic => {}
                    PlugFlowThermalRecord::Coolant {
                        temperature_k,
                        u_w_m2_k,
                    } => v.extend([temperature_k, u_w_m2_k]),
                }
                v.extend(kinetics_numbers(kinetics));
                v
            }
            SpecRecord::Batch {
                kinetics,
                reaction_time_s,
                thermal: t,
            } => {
                let mut v = vec![*reaction_time_s];
                if let BatchThermalRecord::IsothermalK(x) = *t {
                    v.push(x);
                }
                v.extend(kinetics_numbers(kinetics));
                v
            }
        };
        values.into_iter().flatten().collect()
    }
}

fn thermal_record(t: SeparatorThermal) -> ThermalRecord {
    match t {
        SeparatorThermal::Adiabatic => ThermalRecord::Adiabatic,
        SeparatorThermal::OutletTemperature(t) => ThermalRecord::OutletTemperatureK(t),
        SeparatorThermal::Duty(q) => ThermalRecord::DutyW(q),
    }
}

fn separator_thermal(t: ThermalRecord) -> SeparatorThermal {
    match t {
        ThermalRecord::Adiabatic => SeparatorThermal::Adiabatic,
        ThermalRecord::OutletTemperatureK(t) => SeparatorThermal::OutletTemperature(t),
        ThermalRecord::DutyW(q) => SeparatorThermal::Duty(q),
    }
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
            UnitSpec::ConversionReactor {
                reactions,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::ConversionReactor {
                reactions: reactions
                    .into_iter()
                    .map(|r| ReactionRecord {
                        coefficients: r
                            .coefficients
                            .into_iter()
                            .map(|(component, coefficient)| CoefficientRecord {
                                component,
                                coefficient,
                            })
                            .collect(),
                        extent: match r.extent {
                            ReactionExtent::Conversion { key, fraction } => {
                                ExtentRecord::Conversion { key, fraction }
                            }
                            ReactionExtent::MolarExtent(e) => ExtentRecord::MolarExtentMolS(e),
                        },
                    })
                    .collect(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::YieldReactor {
                mass_yields,
                inerts,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::YieldReactor {
                mass_yields: mass_yields
                    .into_iter()
                    .map(|(component, mass_yield)| YieldRecord {
                        component,
                        mass_yield,
                    })
                    .collect(),
                inerts,
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::EquilibriumReactor {
                reactions,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::EquilibriumReactor {
                reactions: reactions.into_iter().map(Into::into).collect(),
                equilibrium: equilibrium.into(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::GibbsReactor {
                inerts,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::GibbsReactor {
                inerts,
                equilibrium: equilibrium.into(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::Cstr {
                kinetics,
                volume_m3,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::Cstr {
                kinetics: kinetics.into(),
                volume_m3,
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::PlugFlow {
                kinetics,
                length_m,
                diameter_m,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::PlugFlow {
                kinetics: kinetics.into(),
                length_m,
                diameter_m,
                thermal: thermal.into(),
                pressure_drop_pa,
            },
            UnitSpec::Batch {
                kinetics,
                reaction_time_s,
                thermal,
            } => SpecRecord::Batch {
                kinetics: kinetics.into(),
                reaction_time_s,
                thermal: thermal.into(),
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
            SpecRecord::ConversionReactor {
                reactions,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::ConversionReactor {
                reactions: reactions
                    .into_iter()
                    .map(|r| Reaction {
                        coefficients: r
                            .coefficients
                            .into_iter()
                            .map(|c| (c.component, c.coefficient))
                            .collect(),
                        extent: match r.extent {
                            ExtentRecord::Conversion { key, fraction } => {
                                ReactionExtent::Conversion { key, fraction }
                            }
                            ExtentRecord::MolarExtentMolS(e) => ReactionExtent::MolarExtent(e),
                        },
                    })
                    .collect(),
                thermal: separator_thermal(thermal),
                pressure_drop_pa,
            },
            SpecRecord::YieldReactor {
                mass_yields,
                inerts,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::YieldReactor {
                mass_yields: mass_yields
                    .into_iter()
                    .map(|y| (y.component, y.mass_yield))
                    .collect(),
                inerts,
                thermal: separator_thermal(thermal),
                pressure_drop_pa,
            },
            SpecRecord::EquilibriumReactor {
                reactions,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::EquilibriumReactor {
                reactions: reactions.into_iter().map(Into::into).collect(),
                equilibrium: equilibrium.into(),
                thermal: separator_thermal(thermal),
                pressure_drop_pa,
            },
            SpecRecord::GibbsReactor {
                inerts,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::GibbsReactor {
                inerts,
                equilibrium: equilibrium.into(),
                thermal: separator_thermal(thermal),
                pressure_drop_pa,
            },
            SpecRecord::Cstr {
                kinetics,
                volume_m3,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::Cstr {
                kinetics: kinetics.into(),
                volume_m3,
                thermal: separator_thermal(thermal),
                pressure_drop_pa,
            },
            SpecRecord::PlugFlow {
                kinetics,
                length_m,
                diameter_m,
                thermal,
                pressure_drop_pa,
            } => UnitSpec::PlugFlow {
                kinetics: kinetics.into(),
                length_m,
                diameter_m,
                thermal: thermal.into(),
                pressure_drop_pa,
            },
            SpecRecord::Batch {
                kinetics,
                reaction_time_s,
                thermal,
            } => UnitSpec::Batch {
                kinetics: kinetics.into(),
                reaction_time_s,
                thermal: thermal.into(),
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
            let json = serde_json::to_string(&SpecRecord::from(spec.clone())).unwrap();
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
