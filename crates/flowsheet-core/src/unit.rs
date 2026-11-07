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
}

impl UnitKind {
    /// All kinds, in palette order.
    pub const ALL: [UnitKind; 7] = [
        UnitKind::Feed,
        UnitKind::Mixer,
        UnitKind::Heater,
        UnitKind::Cooler,
        UnitKind::Valve,
        UnitKind::Separator,
        UnitKind::Product,
    ];

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
        }
    }

    /// Tag series shared by kinds that are numbered together; heaters and
    /// coolers are both heat exchangers (`E-100`, `E-101`, …).
    pub(crate) fn tag_series(self) -> &'static str {
        match self {
            UnitKind::Feed => "Feed",
            UnitKind::Product => "Product",
            UnitKind::Mixer => "M",
            UnitKind::Heater | UnitKind::Cooler => "E",
            UnitKind::Valve => "V",
            UnitKind::Separator => "D",
        }
    }

    /// Generated tag for the `n`-th unit (1-based) of this kind's series,
    /// following common PFD conventions: `Feed 1`, `M-100`, `E-100`,
    /// `V-100`, `D-100`, `Product 1`.
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
        match self {
            UnitKind::Feed => &FEED,
            UnitKind::Product => &PRODUCT,
            UnitKind::Mixer => &MIXER,
            UnitKind::Heater | UnitKind::Cooler | UnitKind::Valve => &SINGLE,
            UnitKind::Separator => &SEPARATOR,
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
        }
    }
}

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
}

impl fmt::Display for SpecProblem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            SpecProblem::Missing => "is missing",
            SpecProblem::NotFinite => "is not a finite number",
            SpecProblem::NotPositive => "must be positive",
            SpecProblem::Negative => "must not be negative",
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
        }
        if let UnitSpec::HeatExchange {
            pressure_drop_pa, ..
        }
        | UnitSpec::Separator {
            pressure_drop_pa, ..
        } = *self
        {
            check_value(
                &mut fixed,
                "pressure drop",
                pressure_drop_pa,
                Sign::NonNegative,
            );
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
                UnitKind::Separator => assert_eq!((inlets, outlets), (1, 2)),
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
