//! Serializable mirrors of the `flowsheet-core` topology and specifications.
//!
//! The engine types stay free of serialization concerns; these records define
//! the stable JSON spelling of the project format. Versions 1 to 3 share the
//! core records below; version 4 (A17 foundation) adds the energy records and
//! one nested record per unit group, `{"type": "<group>", "unit": {"kind":
//! "<kind>", ...}}`, owned by `records/<group>.rs`. Groups add fields with
//! `#[serde(default)]` behind version 4; only the integrator bumps
//! [`crate::FORMAT_VERSION`] (`docs/a17-architecture.md`).

use flowsheet_core::{
    EnergySourceDuty, ExchangerSpec, MixerPressure, PressureRise, PressureSpec, Reaction,
    ReactionExtent, SeparatorThermal, ThermalSpec, UnitKind, UnitSpec,
};
use serde::{Deserialize, Serialize};

use crate::reactor_records::{
    BatchThermalRecord, EquilibriumOptionsRecord, EquilibriumReactionRecord, KineticsRecord,
    PlugFlowThermalRecord, equilibrium_numbers, kinetics_numbers,
};

mod columns;
mod heat_transport;
mod multiphase;
mod reactors;
mod rotating;
mod simple_logical;
mod solids;

pub(crate) use columns::ColumnRecord;
pub(crate) use heat_transport::HeatTransportRecord;
pub(crate) use multiphase::MultiphaseRecord;
pub(crate) use reactors::ReactorRecord;
pub(crate) use rotating::RotatingRecord;
pub(crate) use simple_logical::SimpleLogicalRecord;
pub(crate) use solids::SolidsRecord;

/// Generates [`KindRecord`] with conversions from and to [`UnitKind`]; the
/// JSON name is the snake-case variant name.
macro_rules! kind_records {
    ($($kind:ident),* $(,)?) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
        #[serde(rename_all = "snake_case")]
        pub(crate) enum KindRecord {
            $($kind,)*
        }

        impl From<UnitKind> for KindRecord {
            fn from(k: UnitKind) -> Self {
                match k {
                    $(UnitKind::$kind => KindRecord::$kind,)*
                }
            }
        }

        impl From<KindRecord> for UnitKind {
            fn from(k: KindRecord) -> Self {
                match k {
                    $(KindRecord::$kind => UnitKind::$kind,)*
                }
            }
        }
    };
}

kind_records! {
    Feed, Product, Mixer, Heater, Cooler, Valve, Separator, Splitter, HeatExchanger, Pump,
    Compressor, Expander,
    // Format version 4 (A17).
    EnergySource, EnergySink, ComponentSplitter, Multiplier, Selector, StreamCutter,
    ThreePhaseSeparator, Decanter, Tank, AirCooler, FiredHeater, MultiStreamExchanger,
    MultistageCompressor, PipeSegment, ReliefValve, ShortcutColumn, DistillationColumn, Absorber,
    ThreePhaseColumn, LiquidExtractor, PetroleumFractionator, ConversionReactor, YieldReactor,
    EquilibriumReactor, GibbsReactor, Cstr, PlugFlowReactor, BatchReactor, Cyclone, Hydrocyclone,
    Baghouse, RotaryFilter, SolidsSeparator, Crystallizer, Centrifuge, ElectrostaticPrecipitator,
    VenturiScrubber, Crusher, Screen, Dryer, Granulator, SolidsWasher, Adjust, Set, Spreadsheet,
    Balance, RecycleBlock, ControlBlock, UserUnit,
}

impl KindRecord {
    pub(crate) fn since_format(self) -> u64 {
        match self {
            Self::Feed
            | Self::Product
            | Self::Mixer
            | Self::Heater
            | Self::Cooler
            | Self::Valve
            | Self::Separator
            | Self::Splitter
            | Self::HeatExchanger
            | Self::Pump
            | Self::Compressor
            | Self::Expander => 1,
            _ => 4,
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
    /// Heater duty taken from its energy stream (version 4).
    EnergyStream,
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
    HotEndApproachK(Option<f64>),
    ColdEndApproachK(Option<f64>),
}

/// Energy source duty (version 4).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub(crate) enum EnergyDutyRecord {
    Calculated,
    SpecifiedW(Option<f64>),
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
    // Format version 4 (A17).
    EnergySource {
        duty: EnergyDutyRecord,
    },
    EnergySink,
    Rotating {
        unit: RotatingRecord,
    },
    HeatTransport {
        unit: HeatTransportRecord,
    },
    SimpleLogical {
        unit: SimpleLogicalRecord,
    },
    Reactors {
        unit: ReactorRecord,
    },
    Columns {
        unit: ColumnRecord,
    },
    Multiphase {
        unit: MultiphaseRecord,
    },
    Solids {
        unit: SolidsRecord,
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
            ThermalRecord::Adiabatic | ThermalRecord::EnergyStream => None,
            ThermalRecord::OutletTemperatureK(v) | ThermalRecord::DutyW(v) => v,
        };
        let values: Vec<Option<f64>> = match self {
            SpecRecord::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            } => vec![*temperature_k, *pressure_pa, *molar_flow_mol_s],
            SpecRecord::Product | SpecRecord::EnergySink => vec![],
            SpecRecord::EnergySource { duty } => match duty {
                EnergyDutyRecord::Calculated => vec![],
                EnergyDutyRecord::SpecifiedW(q) => vec![*q],
            },
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
                    | ExchangerRecord::DutyW(v)
                    | ExchangerRecord::HotEndApproachK(v)
                    | ExchangerRecord::ColdEndApproachK(v) => v,
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
            SpecRecord::Rotating { unit } => return unit.numbers(),
            SpecRecord::HeatTransport { unit } => return unit.numbers(),
            SpecRecord::SimpleLogical { unit } => return unit.numbers(),
            SpecRecord::Reactors { .. } => return Vec::new(),
            SpecRecord::Columns { unit } => return unit.numbers(),
            SpecRecord::Multiphase { unit } => return unit.numbers(),
            SpecRecord::Solids { unit } => return unit.numbers(),
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

fn separator_thermal(t: ThermalRecord) -> Result<SeparatorThermal, String> {
    Ok(match t {
        ThermalRecord::Adiabatic => SeparatorThermal::Adiabatic,
        ThermalRecord::OutletTemperatureK(t) => SeparatorThermal::OutletTemperature(t),
        ThermalRecord::DutyW(q) => SeparatorThermal::Duty(q),
        ThermalRecord::EnergyStream => {
            return Err("energy_stream is not valid for a reactor".into());
        }
    })
}

impl From<&UnitSpec> for SpecRecord {
    fn from(s: &UnitSpec) -> Self {
        let drop_or_outlet = |p: PressureSpec| match p {
            PressureSpec::OutletPressure(v) => PressureRecord::OutletPressurePa(v),
            PressureSpec::PressureDrop(v) => PressureRecord::PressureDropPa(v),
        };
        match *s {
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
            UnitSpec::EnergySource { duty } => SpecRecord::EnergySource {
                duty: match duty {
                    EnergySourceDuty::Calculated => EnergyDutyRecord::Calculated,
                    EnergySourceDuty::Specified(q) => EnergyDutyRecord::SpecifiedW(q),
                },
            },
            UnitSpec::EnergySink => SpecRecord::EnergySink,
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
                    ThermalSpec::EnergyStream => ThermalRecord::EnergyStream,
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
                    ExchangerSpec::HotEndApproach(dt) => ExchangerRecord::HotEndApproachK(dt),
                    ExchangerSpec::ColdEndApproach(dt) => ExchangerRecord::ColdEndApproachK(dt),
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
                ref reactions,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::ConversionReactor {
                reactions: reactions
                    .iter()
                    .cloned()
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
                ref mass_yields,
                ref inerts,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::YieldReactor {
                mass_yields: mass_yields
                    .iter()
                    .cloned()
                    .map(|(component, mass_yield)| YieldRecord {
                        component,
                        mass_yield,
                    })
                    .collect(),
                inerts: inerts.clone(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::EquilibriumReactor {
                ref reactions,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::EquilibriumReactor {
                reactions: reactions.iter().cloned().map(Into::into).collect(),
                equilibrium: equilibrium.into(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::GibbsReactor {
                ref inerts,
                equilibrium,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::GibbsReactor {
                inerts: inerts.clone(),
                equilibrium: equilibrium.into(),
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::Cstr {
                ref kinetics,
                volume_m3,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::Cstr {
                kinetics: kinetics.clone().into(),
                volume_m3,
                thermal: thermal_record(thermal),
                pressure_drop_pa,
            },
            UnitSpec::PlugFlow {
                ref kinetics,
                length_m,
                diameter_m,
                thermal,
                pressure_drop_pa,
            } => SpecRecord::PlugFlow {
                kinetics: kinetics.clone().into(),
                length_m,
                diameter_m,
                thermal: thermal.into(),
                pressure_drop_pa,
            },
            UnitSpec::Batch {
                ref kinetics,
                reaction_time_s,
                thermal,
            } => SpecRecord::Batch {
                kinetics: kinetics.clone().into(),
                reaction_time_s,
                thermal: thermal.into(),
            },
            UnitSpec::Rotating(ref s) => SpecRecord::Rotating { unit: s.into() },
            UnitSpec::HeatTransport(ref s) => SpecRecord::HeatTransport { unit: s.into() },
            UnitSpec::SimpleLogical(ref s) => SpecRecord::SimpleLogical { unit: s.into() },
            UnitSpec::Columns(ref s) => SpecRecord::Columns { unit: s.into() },
            UnitSpec::Multiphase(ref s) => SpecRecord::Multiphase { unit: s.into() },
            UnitSpec::Solids(ref s) => SpecRecord::Solids { unit: s.into() },
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
            SpecRecord::EnergySource { duty } => UnitSpec::EnergySource {
                duty: match duty {
                    EnergyDutyRecord::Calculated => EnergySourceDuty::Calculated,
                    EnergyDutyRecord::SpecifiedW(q) => EnergySourceDuty::Specified(q),
                },
            },
            SpecRecord::EnergySink => UnitSpec::EnergySink,
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
                    ThermalRecord::EnergyStream => ThermalSpec::EnergyStream,
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
                    ThermalRecord::EnergyStream => return unexpected("energy_stream"),
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
                    ExchangerRecord::HotEndApproachK(dt) => ExchangerSpec::HotEndApproach(dt),
                    ExchangerRecord::ColdEndApproachK(dt) => ExchangerSpec::ColdEndApproach(dt),
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
                thermal: separator_thermal(thermal)?,
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
                thermal: separator_thermal(thermal)?,
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
                thermal: separator_thermal(thermal)?,
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
                thermal: separator_thermal(thermal)?,
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
                thermal: separator_thermal(thermal)?,
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
            SpecRecord::Rotating { unit } => UnitSpec::Rotating(unit.try_into()?),
            SpecRecord::HeatTransport { unit } => UnitSpec::HeatTransport(unit.try_into()?),
            SpecRecord::SimpleLogical { unit } => UnitSpec::SimpleLogical(unit.try_into()?),
            SpecRecord::Reactors { unit } => unit.default_spec(),
            SpecRecord::Columns { unit } => UnitSpec::Columns(unit.try_into()?),
            SpecRecord::Multiphase { unit } => UnitSpec::Multiphase(unit.try_into()?),
            SpecRecord::Solids { unit } => UnitSpec::Solids(unit.try_into()?),
        })
    }
}

/// A count stored in a group record, checked against the range the group
/// accepts (`u16` keeps the JSON value integral).
pub(crate) fn count_in(
    value: u16,
    field: &str,
    range: std::ops::RangeInclusive<u16>,
) -> Result<u16, String> {
    if range.contains(&value) {
        Ok(value)
    } else {
        Err(format!(
            "{field} {value} is outside {}..={}",
            range.start(),
            range.end()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_default_spec_round_trips() {
        for kind in UnitKind::ALL {
            let spec = kind.default_spec();
            let json = serde_json::to_string(&SpecRecord::from(&spec)).unwrap();
            let back: SpecRecord = serde_json::from_str(&json).unwrap();
            assert_eq!(UnitSpec::try_from(back).unwrap(), spec, "{kind}");
            let k: KindRecord = kind.into();
            assert_eq!(UnitKind::from(k), kind);
        }
    }

    #[test]
    fn group_records_nest_under_a_type_and_name_their_kind() {
        let legacy: SpecRecord =
            serde_json::from_str(r#"{"type":"reactors","unit":{"kind":"cstr"}}"#).unwrap();
        let spec = UnitSpec::try_from(legacy).unwrap();
        assert_eq!(spec, UnitKind::Cstr.default_spec());
        let json = serde_json::to_value(SpecRecord::from(&spec)).unwrap();
        assert_eq!(json["type"], "cstr");
        let json = serde_json::to_string(&SpecRecord::from(
            &UnitKind::DistillationColumn.default_spec(),
        ))
        .unwrap();
        assert_eq!(
            json,
            r#"{"type":"columns","unit":{"kind":"distillation_column","feeds":1,"side_draws":0}}"#
        );
        // Fields added later read with their default from older version-4 files.
        let older: SpecRecord =
            serde_json::from_str(r#"{"type":"columns","unit":{"kind":"distillation_column"}}"#)
                .unwrap();
        assert_eq!(
            UnitSpec::try_from(older).unwrap(),
            UnitKind::DistillationColumn.default_spec()
        );
        // Unknown fields are rejected rather than silently dropped.
        assert!(
            serde_json::from_str::<SpecRecord>(
                r#"{"type":"reactors","unit":{"kind":"cstr","volume_m3":1.0}}"#
            )
            .is_err()
        );
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
        let bad = SpecRecord::Separator {
            thermal: ThermalRecord::EnergyStream,
            pressure_drop_pa: 0.0,
        };
        assert!(UnitSpec::try_from(bad).is_err());
    }
}
