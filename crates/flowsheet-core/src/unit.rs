//! Unit operation kinds, ports, and user specifications.
//!
//! Specification sets follow section 2.3 of `docs/flowsheet-and-ui.md`.
//! Mutually exclusive choices (for example outlet temperature *or* duty) are
//! modelled as enums, so contradictory combinations cannot be represented.

use std::fmt;

/// Stable identifier of a unit operation within one [`crate::Flowsheet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnitId(pub u64);

impl fmt::Display for UnitId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "U{}", self.0)
    }
}

/// Whether material enters or leaves a unit through a port.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PortDirection {
    /// Material enters the unit.
    Inlet,
    /// Material leaves the unit.
    Outlet,
}

/// Static description of one material port of a unit kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PortSpec {
    /// Short, user-facing port name.
    pub name: &'static str,
    /// Flow direction relative to the unit.
    pub direction: PortDirection,
}

const fn inlet(name: &'static str) -> PortSpec {
    PortSpec {
        name,
        direction: PortDirection::Inlet,
    }
}

const fn outlet(name: &'static str) -> PortSpec {
    PortSpec {
        name,
        direction: PortDirection::Outlet,
    }
}

/// Unit operations of the demonstration process, plus flowsheet boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum UnitKind {
    /// Boundary source of a material stream.
    Feed,
    /// Boundary sink of a material stream.
    Product,
    /// Combines two inlet streams into one outlet stream.
    Mixer,
    /// Adds heat to a stream.
    Heater,
    /// Removes heat from a stream.
    Cooler,
    /// Isenthalpic pressure reduction.
    Valve,
    /// Two-phase vapor/liquid separator (drum).
    Separator,
    /// Divides one stream into two streams of identical state.
    Splitter,
    /// Two-stream heat exchanger (hot side and cold side).
    HeatExchanger,
    /// Raises the pressure of a liquid stream.
    Pump,
    /// Raises the pressure of a gas stream.
    Compressor,
    /// Lowers the pressure of a gas stream while recovering work.
    Expander,
}

impl UnitKind {
    /// All kinds, in palette order.
    pub const ALL: [UnitKind; 12] = [
        UnitKind::Feed,
        UnitKind::Mixer,
        UnitKind::Splitter,
        UnitKind::Heater,
        UnitKind::Cooler,
        UnitKind::HeatExchanger,
        UnitKind::Valve,
        UnitKind::Pump,
        UnitKind::Compressor,
        UnitKind::Expander,
        UnitKind::Separator,
        UnitKind::Product,
    ];

    /// Whether the steady-state engine is planned to calculate this kind in
    /// the first usable version. Pumps, compressors, and expanders need a PS
    /// flash and efficiency models and are deferred (proposal section 2.3);
    /// they can be drawn and specified but not yet solved.
    pub fn calculation_planned(self) -> bool {
        !matches!(
            self,
            UnitKind::Pump | UnitKind::Compressor | UnitKind::Expander
        )
    }

    /// User-facing name of the kind.
    pub fn label(self) -> &'static str {
        match self {
            UnitKind::Feed => "Feed",
            UnitKind::Product => "Product",
            UnitKind::Mixer => "Mixer",
            UnitKind::Heater => "Heater",
            UnitKind::Cooler => "Cooler",
            UnitKind::Valve => "Valve",
            UnitKind::Separator => "Separator",
            UnitKind::Splitter => "Splitter",
            UnitKind::HeatExchanger => "Heat exchanger",
            UnitKind::Pump => "Pump",
            UnitKind::Compressor => "Compressor",
            UnitKind::Expander => "Expander",
        }
    }

    /// Tag series shared by kinds that are numbered together; heaters,
    /// coolers, and heat exchangers are all `E-…`.
    pub(crate) fn tag_series(self) -> &'static str {
        match self {
            UnitKind::Feed => "Feed",
            UnitKind::Product => "Product",
            UnitKind::Mixer => "M",
            UnitKind::Heater | UnitKind::Cooler | UnitKind::HeatExchanger => "E",
            UnitKind::Valve => "V",
            UnitKind::Separator => "D",
            UnitKind::Splitter => "TEE",
            UnitKind::Pump => "P",
            UnitKind::Compressor => "K",
            UnitKind::Expander => "EX",
        }
    }

    /// Generated tag for the `n`-th unit (1-based) of this kind's series,
    /// following common PFD conventions: `Feed 1`, `M-100`, `E-100`,
    /// `V-100`, `D-100`, `P-100`, `K-100`, `Product 1`.
    pub(crate) fn tag(self, n: u64) -> String {
        match self {
            UnitKind::Feed | UnitKind::Product => format!("{} {n}", self.tag_series()),
            _ => format!("{}-{}", self.tag_series(), 99 + n),
        }
    }

    /// Material ports of this kind. Port indices in [`PortRef`] refer to
    /// positions in this slice. Energy streams are not modelled yet.
    pub fn ports(self) -> &'static [PortSpec] {
        const FEED: [PortSpec; 1] = [outlet("out")];
        const PRODUCT: [PortSpec; 1] = [inlet("in")];
        const MIXER: [PortSpec; 3] = [inlet("in 1"), inlet("in 2"), outlet("out")];
        const SINGLE: [PortSpec; 2] = [inlet("in"), outlet("out")];
        const SEPARATOR: [PortSpec; 3] = [inlet("in"), outlet("vapor"), outlet("liquid")];
        const SPLITTER: [PortSpec; 3] = [inlet("in"), outlet("out 1"), outlet("out 2")];
        const EXCHANGER: [PortSpec; 4] = [
            inlet("hot in"),
            inlet("cold in"),
            outlet("hot out"),
            outlet("cold out"),
        ];
        match self {
            UnitKind::Feed => &FEED,
            UnitKind::Product => &PRODUCT,
            UnitKind::Mixer => &MIXER,
            UnitKind::Heater
            | UnitKind::Cooler
            | UnitKind::Valve
            | UnitKind::Pump
            | UnitKind::Compressor
            | UnitKind::Expander => &SINGLE,
            UnitKind::Separator => &SEPARATOR,
            UnitKind::Splitter => &SPLITTER,
            UnitKind::HeatExchanger => &EXCHANGER,
        }
    }

    /// Default specification for a newly created unit of this kind
    /// (defaults column of the proposal's section 2.3).
    pub fn default_spec(self) -> UnitSpec {
        match self {
            UnitKind::Feed => UnitSpec::Feed {
                temperature_k: None,
                pressure_pa: None,
                molar_flow_mol_s: None,
            },
            UnitKind::Product => UnitSpec::Product,
            UnitKind::Mixer => UnitSpec::Mixer {
                outlet_pressure: MixerPressure::LowestInlet,
            },
            UnitKind::Heater | UnitKind::Cooler => UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::Valve => UnitSpec::Valve {
                pressure: PressureSpec::OutletPressure(None),
            },
            UnitKind::Separator => UnitSpec::Separator {
                thermal: SeparatorThermal::Adiabatic,
                pressure_drop_pa: 0.0,
            },
            UnitKind::Splitter => UnitSpec::Splitter {
                fraction_to_first: None,
            },
            UnitKind::HeatExchanger => UnitSpec::HeatExchanger {
                thermal: ExchangerSpec::HotOutletTemperature(None),
                hot_pressure_drop_pa: 0.0,
                cold_pressure_drop_pa: 0.0,
            },
            UnitKind::Pump | UnitKind::Compressor => UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(None),
                efficiency: DEFAULT_EFFICIENCY,
            },
            UnitKind::Expander => UnitSpec::Expand {
                pressure: PressureSpec::OutletPressure(None),
                efficiency: DEFAULT_EFFICIENCY,
            },
        }
    }
}

/// Default isentropic (compressor, expander) or hydraulic (pump) efficiency,
/// shown in the form and editable; 75 % is a common simulator default.
pub const DEFAULT_EFFICIENCY: f64 = 0.75;

impl fmt::Display for UnitKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.label())
    }
}

/// Mixer outlet pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MixerPressure {
    /// Lowest inlet pressure (default).
    LowestInlet,
    /// Specified absolute pressure in Pa; must not exceed the lowest inlet
    /// pressure (checked once inlet states are calculated).
    Specified(Option<f64>),
}

/// Heater/cooler thermal specification: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalSpec {
    /// Outlet temperature in K.
    OutletTemperature(Option<f64>),
    /// Duty in W, positive into the unit.
    Duty(Option<f64>),
}

/// Separator thermal specification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SeparatorThermal {
    /// No heat exchange (default).
    Adiabatic,
    /// Outlet temperature in K.
    OutletTemperature(Option<f64>),
    /// Duty in W, positive into the unit.
    Duty(Option<f64>),
}

/// Valve pressure specification: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureSpec {
    /// Absolute outlet pressure in Pa.
    OutletPressure(Option<f64>),
    /// Pressure drop in Pa.
    PressureDrop(Option<f64>),
}

/// Pump/compressor pressure specification: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureRise {
    /// Absolute outlet pressure in Pa.
    OutletPressure(Option<f64>),
    /// Pressure increase in Pa.
    PressureIncrease(Option<f64>),
}

/// Two-stream heat exchanger specification: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExchangerSpec {
    /// Hot-side outlet temperature in K.
    HotOutletTemperature(Option<f64>),
    /// Cold-side outlet temperature in K.
    ColdOutletTemperature(Option<f64>),
    /// Heat transferred from the hot to the cold side in W.
    Duty(Option<f64>),
}

/// User-entered specification of a unit, in SI units with absolute pressure.
///
/// `None` means "not yet specified". The feed composition is not part of
/// this model yet because component identities are defined with work
/// package A2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitSpec {
    /// Feed state and flow.
    Feed {
        /// Temperature in K.
        temperature_k: Option<f64>,
        /// Absolute pressure in Pa.
        pressure_pa: Option<f64>,
        /// Total molar flow in mol/s.
        molar_flow_mol_s: Option<f64>,
    },
    /// Product boundary; no specification.
    Product,
    /// Adiabatic mixer.
    Mixer {
        /// Outlet pressure choice.
        outlet_pressure: MixerPressure,
    },
    /// Heater or cooler.
    HeatExchange {
        /// Outlet temperature or duty.
        thermal: ThermalSpec,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Isenthalpic valve.
    Valve {
        /// Outlet pressure or pressure drop.
        pressure: PressureSpec,
    },
    /// Two-phase separator.
    Separator {
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Splitter with two outlets (n − 1 = 1 split fraction).
    Splitter {
        /// Fraction of the inlet flow leaving through `out 1`, 0 to 1.
        fraction_to_first: Option<f64>,
    },
    /// Two-stream heat exchanger.
    HeatExchanger {
        /// Hot outlet temperature, cold outlet temperature, or duty.
        thermal: ExchangerSpec,
        /// Hot-side pressure drop in Pa.
        hot_pressure_drop_pa: f64,
        /// Cold-side pressure drop in Pa.
        cold_pressure_drop_pa: f64,
    },
    /// Pump or compressor.
    Pressurize {
        /// Outlet pressure or pressure increase.
        pressure: PressureRise,
        /// Efficiency, greater than 0 and at most 1.
        efficiency: f64,
    },
    /// Expander.
    Expand {
        /// Outlet pressure or pressure drop.
        pressure: PressureSpec,
        /// Isentropic efficiency, greater than 0 and at most 1.
        efficiency: f64,
    },
}

/// What is wrong with one specification field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecProblem {
    /// The field is required but unset.
    Missing,
    /// The value is NaN or infinite.
    NotFinite,
    /// The value must be strictly positive (temperature, absolute pressure).
    NotPositive,
    /// The value must not be negative (flows, pressure drops).
    Negative,
    /// The value must lie between 0 and 1 (split fractions).
    NotAFraction,
    /// The value must be greater than 0 and at most 1 (efficiencies).
    NotAnEfficiency,
}

impl fmt::Display for SpecProblem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            SpecProblem::Missing => "is missing",
            SpecProblem::NotFinite => "is not a finite number",
            SpecProblem::NotPositive => "must be positive",
            SpecProblem::Negative => "must not be negative",
            SpecProblem::NotAFraction => "must be between 0 and 1",
            SpecProblem::NotAnEfficiency => "must be greater than 0 and at most 1",
        })
    }
}

/// One problem found by [`UnitSpec::issues`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpecIssue {
    /// User-facing field name.
    pub field: &'static str,
    /// The problem.
    pub problem: SpecProblem,
}

#[derive(Clone, Copy)]
enum Sign {
    Positive,
    NonNegative,
    Fraction,
    Efficiency,
    Any,
}

impl UnitSpec {
    /// Missing and invalid fields. Checks that need calculated inlet states
    /// (for example a valve outlet pressure above its inlet pressure) belong
    /// to the flowsheet engine (A8).
    pub fn issues(&self) -> Vec<SpecIssue> {
        let mut out = Vec::new();
        let mut opt = |field, value: Option<f64>, sign| match value {
            None => out.push(SpecIssue {
                field,
                problem: SpecProblem::Missing,
            }),
            Some(v) => check_value(&mut out, field, v, sign),
        };
        let mut fixed = Vec::new();
        match *self {
            UnitSpec::Product
            | UnitSpec::Mixer {
                outlet_pressure: MixerPressure::LowestInlet,
            }
            | UnitSpec::Separator {
                thermal: SeparatorThermal::Adiabatic,
                ..
            } => {}
            UnitSpec::Feed {
                temperature_k,
                pressure_pa,
                molar_flow_mol_s,
            } => {
                opt("temperature", temperature_k, Sign::Positive);
                opt("pressure", pressure_pa, Sign::Positive);
                opt("molar flow", molar_flow_mol_s, Sign::NonNegative);
            }
            UnitSpec::Mixer {
                outlet_pressure: MixerPressure::Specified(p),
            } => opt("outlet pressure", p, Sign::Positive),
            UnitSpec::HeatExchange { thermal, .. } => match thermal {
                ThermalSpec::OutletTemperature(t) => opt("outlet temperature", t, Sign::Positive),
                ThermalSpec::Duty(q) => opt("duty", q, Sign::Any),
            },
            UnitSpec::Valve { pressure } => match pressure {
                PressureSpec::OutletPressure(p) => opt("outlet pressure", p, Sign::Positive),
                PressureSpec::PressureDrop(dp) => opt("pressure drop", dp, Sign::NonNegative),
            },
            UnitSpec::Separator { thermal, .. } => match thermal {
                SeparatorThermal::Adiabatic => {}
                SeparatorThermal::OutletTemperature(t) => {
                    opt("outlet temperature", t, Sign::Positive)
                }
                SeparatorThermal::Duty(q) => opt("duty", q, Sign::Any),
            },
            UnitSpec::Splitter { fraction_to_first } => {
                opt("split fraction", fraction_to_first, Sign::Fraction)
            }
            UnitSpec::HeatExchanger { thermal, .. } => match thermal {
                ExchangerSpec::HotOutletTemperature(t) => {
                    opt("hot outlet temperature", t, Sign::Positive)
                }
                ExchangerSpec::ColdOutletTemperature(t) => {
                    opt("cold outlet temperature", t, Sign::Positive)
                }
                ExchangerSpec::Duty(q) => opt("duty", q, Sign::NonNegative),
            },
            UnitSpec::Pressurize { pressure, .. } => match pressure {
                PressureRise::OutletPressure(p) => opt("outlet pressure", p, Sign::Positive),
                PressureRise::PressureIncrease(dp) => {
                    opt("pressure increase", dp, Sign::NonNegative)
                }
            },
            UnitSpec::Expand { pressure, .. } => match pressure {
                PressureSpec::OutletPressure(p) => opt("outlet pressure", p, Sign::Positive),
                PressureSpec::PressureDrop(dp) => opt("pressure drop", dp, Sign::NonNegative),
            },
        }
        match *self {
            UnitSpec::HeatExchange {
                pressure_drop_pa, ..
            }
            | UnitSpec::Separator {
                pressure_drop_pa, ..
            } => check_value(
                &mut fixed,
                "pressure drop",
                pressure_drop_pa,
                Sign::NonNegative,
            ),
            UnitSpec::HeatExchanger {
                hot_pressure_drop_pa,
                cold_pressure_drop_pa,
                ..
            } => {
                check_value(
                    &mut fixed,
                    "hot-side pressure drop",
                    hot_pressure_drop_pa,
                    Sign::NonNegative,
                );
                check_value(
                    &mut fixed,
                    "cold-side pressure drop",
                    cold_pressure_drop_pa,
                    Sign::NonNegative,
                );
            }
            UnitSpec::Pressurize { efficiency, .. } | UnitSpec::Expand { efficiency, .. } => {
                check_value(&mut fixed, "efficiency", efficiency, Sign::Efficiency)
            }
            _ => {}
        }
        out.extend(fixed);
        out
    }
}

fn check_value(out: &mut Vec<SpecIssue>, field: &'static str, v: f64, sign: Sign) {
    let problem = if !v.is_finite() {
        Some(SpecProblem::NotFinite)
    } else {
        match sign {
            Sign::Positive if v <= 0.0 => Some(SpecProblem::NotPositive),
            Sign::NonNegative if v < 0.0 => Some(SpecProblem::Negative),
            Sign::Fraction if !(0.0..=1.0).contains(&v) => Some(SpecProblem::NotAFraction),
            Sign::Efficiency if v <= 0.0 || v > 1.0 => Some(SpecProblem::NotAnEfficiency),
            _ => None,
        }
    };
    if let Some(problem) = problem {
        out.push(SpecIssue { field, problem });
    }
}

/// A unit operation placed on the flowsheet.
#[derive(Debug, Clone, PartialEq)]
pub struct Unit {
    /// Identifier, unique within its flowsheet.
    pub id: UnitId,
    /// Kind of unit operation.
    pub kind: UnitKind,
    /// User-facing tag, e.g. `M-100`.
    pub name: String,
    /// User specification.
    pub spec: UnitSpec,
}

/// Reference to one material port of one unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PortRef {
    /// Owning unit.
    pub unit: UnitId,
    /// Index into [`UnitKind::ports`].
    pub port: usize,
}

impl PortRef {
    /// Creates a port reference.
    pub fn new(unit: UnitId, port: usize) -> Self {
        Self { unit, port }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields(spec: UnitSpec) -> Vec<(&'static str, SpecProblem)> {
        spec.issues()
            .into_iter()
            .map(|i| (i.field, i.problem))
            .collect()
    }

    #[test]
    fn every_kind_has_ports_consistent_with_its_role() {
        for kind in UnitKind::ALL {
            let ports = kind.ports();
            let inlets = ports
                .iter()
                .filter(|p| p.direction == PortDirection::Inlet)
                .count();
            let outlets = ports.len() - inlets;
            match kind {
                UnitKind::Feed => assert_eq!((inlets, outlets), (0, 1)),
                UnitKind::Product => assert_eq!((inlets, outlets), (1, 0)),
                UnitKind::Mixer => assert_eq!((inlets, outlets), (2, 1)),
                UnitKind::Separator | UnitKind::Splitter => {
                    assert_eq!((inlets, outlets), (1, 2))
                }
                UnitKind::HeatExchanger => assert_eq!((inlets, outlets), (2, 2)),
                _ => assert_eq!((inlets, outlets), (1, 1)),
            }
        }
    }

    #[test]
    fn tags_follow_pfd_conventions() {
        assert_eq!(UnitKind::Mixer.tag(1), "M-100");
        assert_eq!(UnitKind::Cooler.tag(2), "E-101");
        assert_eq!(UnitKind::Separator.tag(1), "D-100");
        assert_eq!(UnitKind::Feed.tag(2), "Feed 2");
        assert_eq!(UnitKind::Heater.tag_series(), UnitKind::Cooler.tag_series());
        assert_eq!(
            UnitKind::HeatExchanger.tag_series(),
            UnitKind::Heater.tag_series()
        );
        assert_eq!(UnitKind::Pump.tag(1), "P-100");
        assert_eq!(UnitKind::Compressor.tag(1), "K-100");
        assert_eq!(UnitKind::Expander.tag(1), "EX-100");
        assert_eq!(UnitKind::Splitter.tag(1), "TEE-100");
    }

    #[test]
    fn rotating_equipment_is_drawable_but_not_yet_calculated() {
        for kind in UnitKind::ALL {
            let deferred = matches!(
                kind,
                UnitKind::Pump | UnitKind::Compressor | UnitKind::Expander
            );
            assert_eq!(kind.calculation_planned(), !deferred, "{kind}");
        }
    }

    #[test]
    fn new_unit_defaults_and_ranges() {
        use SpecProblem::*;
        assert_eq!(
            fields(UnitKind::Splitter.default_spec()),
            [("split fraction", Missing)]
        );
        assert_eq!(
            fields(UnitKind::HeatExchanger.default_spec()),
            [("hot outlet temperature", Missing)]
        );
        // Efficiency has a visible default, so only the pressure is required.
        for kind in [UnitKind::Pump, UnitKind::Compressor, UnitKind::Expander] {
            assert_eq!(fields(kind.default_spec()), [("outlet pressure", Missing)]);
        }
        for (f, ok) in [
            (0.0, true),
            (1.0, true),
            (0.3, true),
            (-0.1, false),
            (1.2, false),
        ] {
            let issues = fields(UnitSpec::Splitter {
                fraction_to_first: Some(f),
            });
            assert_eq!(issues.is_empty(), ok, "fraction {f}");
        }
        for (eta, ok) in [(1.0, true), (0.5, true), (0.0, false), (1.01, false)] {
            let issues = fields(UnitSpec::Pressurize {
                pressure: PressureRise::PressureIncrease(Some(1.0e5)),
                efficiency: eta,
            });
            assert_eq!(issues.is_empty(), ok, "efficiency {eta}");
            if !ok {
                assert_eq!(issues, [("efficiency", NotAnEfficiency)]);
            }
        }
        let hx = UnitSpec::HeatExchanger {
            thermal: ExchangerSpec::Duty(Some(-1.0)),
            hot_pressure_drop_pa: 0.0,
            cold_pressure_drop_pa: -5.0,
        };
        assert_eq!(
            fields(hx),
            [("duty", Negative), ("cold-side pressure drop", Negative)]
        );
    }

    #[test]
    fn defaults_require_only_the_inputs_without_a_default() {
        use SpecProblem::Missing;
        assert_eq!(
            fields(UnitKind::Feed.default_spec()),
            [
                ("temperature", Missing),
                ("pressure", Missing),
                ("molar flow", Missing)
            ]
        );
        assert_eq!(
            fields(UnitKind::Cooler.default_spec()),
            [("outlet temperature", Missing)]
        );
        assert_eq!(
            fields(UnitKind::Valve.default_spec()),
            [("outlet pressure", Missing)]
        );
        // Lowest-inlet mixer pressure and an adiabatic separator need no input.
        assert!(UnitKind::Mixer.default_spec().issues().is_empty());
        assert!(UnitKind::Separator.default_spec().issues().is_empty());
        assert!(UnitKind::Product.default_spec().issues().is_empty());
    }

    #[test]
    fn invalid_values_are_reported_per_field() {
        let feed = UnitSpec::Feed {
            temperature_k: Some(0.0),
            pressure_pa: Some(f64::NAN),
            molar_flow_mol_s: Some(-1.0),
        };
        assert_eq!(
            fields(feed),
            [
                ("temperature", SpecProblem::NotPositive),
                ("pressure", SpecProblem::NotFinite),
                ("molar flow", SpecProblem::Negative)
            ]
        );
        let cooler = UnitSpec::HeatExchange {
            thermal: ThermalSpec::Duty(Some(-5.0e3)),
            pressure_drop_pa: -1.0,
        };
        // Negative duty is valid (heat removed); negative pressure drop is not.
        assert_eq!(fields(cooler), [("pressure drop", SpecProblem::Negative)]);
        let valve = UnitSpec::Valve {
            pressure: PressureSpec::PressureDrop(Some(0.0)),
        };
        assert!(valve.issues().is_empty());
    }
}
