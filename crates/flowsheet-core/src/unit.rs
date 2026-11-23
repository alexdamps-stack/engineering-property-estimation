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
    /// Stoichiometric reactor: reactions with specified conversion or extent
    /// (Aspen HYSYS Conversion Reactor, Aspen Plus RStoic).
    ConversionReactor,
    /// Yield reactor: specified product mass yields without reaction
    /// stoichiometry (Aspen HYSYS Yield Shift Reactor, Aspen Plus RYield).
    YieldReactor,
    /// Reactions at chemical equilibrium (Aspen HYSYS Equilibrium Reactor,
    /// Aspen Plus REquil, DWSIM Equilibrium Reactor).
    EquilibriumReactor,
    /// Minimum Gibbs energy subject to element balances (Aspen Plus RGibbs,
    /// HYSYS and DWSIM Gibbs Reactor).
    GibbsReactor,
    /// Continuous stirred-tank reactor with kinetic reactions (RCSTR).
    Cstr,
    /// Plug-flow (tubular) reactor with kinetic reactions (RPlug).
    PlugFlowReactor,
    /// Batch reactor at constant pressure, represented in the steady-state
    /// flowsheet by its time-averaged product (RBatch).
    BatchReactor,
}

impl UnitKind {
    /// All kinds, in palette order.
    pub const ALL: [UnitKind; 19] = [
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
        UnitKind::ConversionReactor,
        UnitKind::YieldReactor,
        UnitKind::EquilibriumReactor,
        UnitKind::GibbsReactor,
        UnitKind::Cstr,
        UnitKind::PlugFlowReactor,
        UnitKind::BatchReactor,
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
            UnitKind::ConversionReactor => "Conversion reactor",
            UnitKind::YieldReactor => "Yield reactor",
            UnitKind::EquilibriumReactor => "Equilibrium reactor",
            UnitKind::GibbsReactor => "Gibbs reactor",
            UnitKind::Cstr => "CSTR",
            UnitKind::PlugFlowReactor => "Plug-flow reactor",
            UnitKind::BatchReactor => "Batch reactor",
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
            UnitKind::ConversionReactor
            | UnitKind::YieldReactor
            | UnitKind::EquilibriumReactor
            | UnitKind::GibbsReactor
            | UnitKind::Cstr
            | UnitKind::PlugFlowReactor
            | UnitKind::BatchReactor => "R",
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

    /// Inverse of [`Self::tag`]: the series number of a generated tag, or
    /// `None` for a user-chosen name.
    pub(crate) fn tag_number(self, name: &str) -> Option<u64> {
        let n = match self {
            UnitKind::Feed | UnitKind::Product => name
                .strip_prefix(self.tag_series())?
                .strip_prefix(' ')?
                .parse::<u64>()
                .ok()?,
            _ => name
                .strip_prefix(self.tag_series())?
                .strip_prefix('-')?
                .parse::<u64>()
                .ok()?
                .checked_sub(99)?,
        };
        (n > 0 && self.tag(n) == name).then_some(n)
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
            | UnitKind::Expander
            | UnitKind::ConversionReactor
            | UnitKind::YieldReactor
            | UnitKind::EquilibriumReactor
            | UnitKind::GibbsReactor
            | UnitKind::Cstr
            | UnitKind::PlugFlowReactor
            | UnitKind::BatchReactor => &SINGLE,
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
            // Reactors have no adiabatic default: an adiabatic exothermic reaction
            // readily leaves the package temperature domain, so the user chooses.
            UnitKind::ConversionReactor => UnitSpec::ConversionReactor {
                reactions: Vec::new(),
                thermal: SeparatorThermal::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::YieldReactor => UnitSpec::YieldReactor {
                mass_yields: Vec::new(),
                inerts: Vec::new(),
                thermal: SeparatorThermal::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::EquilibriumReactor => UnitSpec::EquilibriumReactor {
                reactions: Vec::new(),
                equilibrium: EquilibriumOptions::default(),
                thermal: SeparatorThermal::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::GibbsReactor => UnitSpec::GibbsReactor {
                inerts: Vec::new(),
                equilibrium: EquilibriumOptions::default(),
                thermal: SeparatorThermal::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::Cstr => UnitSpec::Cstr {
                kinetics: Kinetics::default(),
                volume_m3: None,
                thermal: SeparatorThermal::OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::PlugFlowReactor => UnitSpec::PlugFlow {
                kinetics: Kinetics::default(),
                length_m: None,
                diameter_m: None,
                thermal: PlugFlowThermal::Isothermal(None),
                pressure_drop_pa: 0.0,
            },
            UnitKind::BatchReactor => UnitSpec::Batch {
                kinetics: Kinetics::default(),
                reaction_time_s: None,
                thermal: BatchThermal::Isothermal(None),
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

/// Separator and reactor thermal specification.
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
/// `None` means "not yet specified". Feed compositions are supplied separately
/// as [`crate::FeedCompositions`], with explicit identities and mole basis.
/// Reactor chemistry names components by package identifier; whether those
/// identifiers exist in the package and whether reactions conserve elements is
/// checked by the engine before solving.
#[derive(Debug, Clone, PartialEq)]
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
    /// Stoichiometric (conversion) reactor with one outlet.
    ConversionReactor {
        /// Reactions, applied in order; each step sees the outlet of the
        /// previous one.
        reactions: Vec<Reaction>,
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Yield reactor with one outlet.
    YieldReactor {
        /// Relative mass yields `(component, weight)` of the reacting part of
        /// the feed. The weights are normalized so that mass is conserved;
        /// unlisted, non-inert components leave with zero flow.
        mass_yields: Vec<(String, f64)>,
        /// Components that pass through unchanged and take no part in the yield.
        inerts: Vec<String>,
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Equilibrium reactor: the listed reactions reach chemical equilibrium.
    EquilibriumReactor {
        /// Independent reactions with their equilibrium constants.
        reactions: Vec<EquilibriumReaction>,
        /// Fugacity basis and temperature approach.
        equilibrium: EquilibriumOptions,
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Gibbs reactor: minimum Gibbs energy over all non-inert components
    /// subject to element conservation.
    GibbsReactor {
        /// Components that pass through unchanged.
        inerts: Vec<String>,
        /// Fugacity basis and temperature approach.
        equilibrium: EquilibriumOptions,
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Continuous stirred-tank reactor: perfectly mixed contents at outlet
    /// conditions.
    Cstr {
        /// Kinetic reactions and catalyst loading.
        kinetics: Kinetics,
        /// Reactor volume in m³.
        volume_m3: Option<f64>,
        /// Adiabatic, outlet temperature, or duty.
        thermal: SeparatorThermal,
        /// Pressure drop in Pa.
        pressure_drop_pa: f64,
    },
    /// Plug-flow reactor: kinetic reactions integrated along a tube.
    PlugFlow {
        /// Kinetic reactions and catalyst loading.
        kinetics: Kinetics,
        /// Tube length in m.
        length_m: Option<f64>,
        /// Inner tube diameter in m.
        diameter_m: Option<f64>,
        /// Isothermal, adiabatic, or exchange with a coolant.
        thermal: PlugFlowThermal,
        /// Pressure drop in Pa, linear along the length.
        pressure_drop_pa: f64,
    },
    /// Batch reactor at the inlet pressure. Each charge reacts for the
    /// reaction time; the outlet is the time-averaged product, so outlet flow
    /// equals the inlet mass flow.
    Batch {
        /// Kinetic reactions and catalyst loading.
        kinetics: Kinetics,
        /// Reaction time per batch in s.
        reaction_time_s: Option<f64>,
        /// Isothermal or adiabatic.
        thermal: BatchThermal,
    },
}

/// One reaction of an equilibrium reactor.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumReaction {
    /// Stoichiometric coefficients by package component identifier: negative
    /// for reactants, positive for products. Must conserve every element.
    pub coefficients: Vec<(String, f64)>,
    /// Source of the equilibrium constant.
    pub constant: EquilibriumConstant,
}

/// Equilibrium constant of one reaction (fugacity basis, standard state the
/// ideal gas at 1 bar): exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EquilibriumConstant {
    /// From the ideal-gas standard Gibbs energies of the package data,
    /// `ln K = -sum_i nu_i g0_i(T) / (R T)` (Aspen and DWSIM default).
    GibbsEnergy,
    /// User correlation `ln K = a + b/T + c ln T + d T`, T in K (the
    /// Aspen HYSYS and DWSIM form).
    Correlation(Correlation),
}

/// `ln K = a + b/T + c ln T + d T` with T in K. Used for equilibrium
/// constants and LHHW driving-force and adsorption constants.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Correlation {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl Correlation {
    /// `ln K` at `t` K.
    pub fn ln_value(&self, t: f64) -> f64 {
        self.a + self.b / t + self.c * t.ln() + self.d * t
    }
}

/// Fugacities used in the equilibrium conditions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FugacityBasis {
    /// Fugacity coefficients of the property package for the reacting fluid
    /// (single phase, minimum-Gibbs root).
    #[default]
    Package,
    /// Ideal gas: fugacity equals partial pressure.
    IdealGas,
}

/// Options shared by the equilibrium and Gibbs reactors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquilibriumOptions {
    /// Fugacity basis of the equilibrium conditions.
    pub fugacity: FugacityBasis,
    /// Temperature approach in K: equilibrium is evaluated at the outlet
    /// temperature plus this value (Aspen REquil/RGibbs). Zero by default.
    pub temperature_approach_k: f64,
}

impl Default for EquilibriumOptions {
    fn default() -> Self {
        Self {
            fugacity: FugacityBasis::Package,
            temperature_approach_k: 0.0,
        }
    }
}

/// Kinetic reactions of a CSTR, plug-flow or batch reactor.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Kinetics {
    /// The reactions; rates add.
    pub reactions: Vec<KineticReaction>,
    /// Catalyst bulk density in kg of catalyst per m³ of reactor, required
    /// when a rate is per catalyst mass.
    pub catalyst_bulk_density_kg_m3: Option<f64>,
}

/// One kinetic reaction.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticReaction {
    /// Stoichiometric coefficients by package component identifier: negative
    /// for reactants, positive for products. Must conserve every element.
    pub coefficients: Vec<(String, f64)>,
    /// Phase whose composition and volume the rate refers to.
    pub phase: ReactingPhase,
    /// Quantity the orders apply to.
    pub basis: ConcentrationBasis,
    /// Whether the rate is per reactor volume or per catalyst mass.
    pub rate_basis: RateBasis,
    /// The rate expression.
    pub rate: RateLaw,
}

/// Reacting phase of a kinetic reaction. A single supercritical phase
/// counts as either.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReactingPhase {
    #[default]
    Vapor,
    Liquid,
}

/// Concentration basis of the rate expression (Aspen and DWSIM options).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConcentrationBasis {
    /// mol/m³ of the reacting phase.
    #[default]
    MolarConcentration,
    /// kg/m³ of the reacting phase.
    MassConcentration,
    /// Mole fraction in the reacting phase.
    MoleFraction,
    /// Mass fraction in the reacting phase.
    MassFraction,
    /// Partial pressure `x_i p` in Pa.
    PartialPressure,
    /// Fugacity `x_i phi_i p` in Pa.
    Fugacity,
}

/// What the rate is expressed per.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RateBasis {
    /// mol/(m³ s) of reacting-phase volume.
    #[default]
    ReactorVolume,
    /// mol/(kg s) of catalyst; the catalyst mass is the bulk density times
    /// the reactor volume.
    CatalystMass,
}

/// `k = k0 (T/T_ref)^n exp(-E/R (1/T - 1/T_ref))`, or without a reference
/// temperature `k = k0 T^n exp(-E/(R T))` (Aspen power-law form).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Arrhenius {
    /// Pre-exponential factor in the units implied by the bases and orders.
    pub pre_exponential: Option<f64>,
    /// Activation energy in J/mol.
    pub activation_energy_j_mol: f64,
    /// Temperature exponent n.
    pub temperature_exponent: f64,
    /// Reference temperature in K, if any.
    pub reference_temperature_k: Option<f64>,
}

impl Arrhenius {
    /// Rate constant at `t` K; `None` while the pre-exponential is unset.
    pub fn value(&self, t: f64) -> Option<f64> {
        let k0 = self.pre_exponential?;
        let r = thermo_types::MOLAR_GAS_CONSTANT;
        Some(match self.reference_temperature_k {
            Some(t0) => {
                k0 * (t / t0).powf(self.temperature_exponent)
                    * (-self.activation_energy_j_mol / r * (1.0 / t - 1.0 / t0)).exp()
            }
            None => {
                k0 * t.powf(self.temperature_exponent)
                    * (-self.activation_energy_j_mol / (r * t)).exp()
            }
        })
    }
}

/// A rate constant with exponents by component.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct RateTerm {
    pub constant: Arrhenius,
    /// Exponents of the basis quantities; unlisted components have order 0.
    pub orders: Vec<(String, f64)>,
}

/// An LHHW constant `exp(a + b/T + c ln T + d T)` with exponents.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LhhwTerm {
    pub constant: Correlation,
    pub orders: Vec<(String, f64)>,
}

/// Rate expression of one kinetic reaction.
#[derive(Debug, Clone, PartialEq)]
pub enum RateLaw {
    /// `r = k_f prod C^a - k_r prod C^b` (reverse optional).
    PowerLaw {
        forward: RateTerm,
        reverse: Option<RateTerm>,
    },
    /// Langmuir–Hinshelwood–Hougen–Watson:
    /// `r = k (K1 prod C^a - K2 prod C^b) / (sum_i K_i prod C^c_i)^m`.
    Lhhw {
        /// Kinetic factor k.
        kinetic: Arrhenius,
        /// Forward driving-force term K1 prod C^a.
        driving_forward: LhhwTerm,
        /// Reverse driving-force term K2 prod C^b, if reversible.
        driving_reverse: Option<LhhwTerm>,
        /// Adsorption terms; an empty list means a denominator of 1.
        adsorption: Vec<LhhwTerm>,
        /// Adsorption exponent m.
        exponent: f64,
    },
}

/// Thermal specification of a plug-flow reactor: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlugFlowThermal {
    /// Constant temperature along the tube (the feed is brought to it at the
    /// inlet); the duty is reported.
    Isothermal(Option<f64>),
    /// No heat exchange.
    Adiabatic,
    /// Wall heat flux `U (T_c - T)` per tube surface with a constant coolant
    /// temperature.
    Coolant {
        /// Coolant temperature in K.
        temperature_k: Option<f64>,
        /// Overall heat-transfer coefficient in W/(m² K).
        u_w_m2_k: Option<f64>,
    },
}

/// Thermal specification of a batch reactor: exactly one of these.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BatchThermal {
    /// Constant temperature (the charge is brought to it at the start).
    Isothermal(Option<f64>),
    /// No heat exchange.
    Adiabatic,
}

/// One reaction of a conversion reactor.
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    /// Stoichiometric coefficients by package component identifier: negative
    /// for reactants, positive for products. Must conserve every element.
    pub coefficients: Vec<(String, f64)>,
    /// How far the reaction proceeds.
    pub extent: ReactionExtent,
}

/// Extent specification of one reaction: exactly one of these.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionExtent {
    /// Fraction, 0 to 1, of the key reactant entering this reaction step
    /// that is consumed.
    Conversion {
        /// Key reactant; must have a negative coefficient.
        key: Option<String>,
        /// Fractional conversion of the key reactant.
        fraction: Option<f64>,
    },
    /// Molar extent in mol/s, not negative; component `i` changes by
    /// `coefficient_i * extent`.
    MolarExtent(Option<f64>),
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
    /// The sign contradicts the unit kind (a heater duty below zero or a
    /// cooler duty above zero).
    ContradictsUnitKind,
    /// A stoichiometric coefficient is zero.
    Zero,
    /// The same component is listed twice.
    Duplicate,
    /// A reaction lacks a reactant (negative coefficient) or a product.
    NeedsReactantAndProduct,
    /// The conversion key is not a reactant of its reaction.
    KeyNotReactant,
    /// All mass yields are zero.
    AllZero,
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
            SpecProblem::ContradictsUnitKind => {
                "contradicts the unit kind (a heater adds heat, a cooler removes it)"
            }
            SpecProblem::Zero => "must not be zero",
            SpecProblem::Duplicate => "lists a component twice",
            SpecProblem::NeedsReactantAndProduct => {
                "needs at least one reactant (negative) and one product (positive)"
            }
            SpecProblem::KeyNotReactant => "must be a reactant of its reaction",
            SpecProblem::AllZero => "must not all be zero",
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
        let mut chemistry = Vec::new();
        match *self {
            UnitSpec::ConversionReactor { ref reactions, .. } => {
                reaction_issues(&mut chemistry, reactions);
            }
            UnitSpec::YieldReactor {
                ref mass_yields,
                ref inerts,
                ..
            } => yield_issues(&mut chemistry, mass_yields, inerts),
            UnitSpec::EquilibriumReactor { ref reactions, .. } => {
                if reactions.is_empty() {
                    chemistry.push(SpecIssue {
                        field: "reactions",
                        problem: SpecProblem::Missing,
                    });
                }
                for r in reactions {
                    coefficient_issues(&mut chemistry, &r.coefficients);
                    if let EquilibriumConstant::Correlation(c) = r.constant {
                        for v in [c.a, c.b, c.c, c.d] {
                            check_value(&mut chemistry, "equilibrium constant", v, Sign::Any);
                        }
                    }
                }
            }
            UnitSpec::GibbsReactor { ref inerts, .. } => {
                let mut seen = Vec::new();
                for id in inerts {
                    component_issue(&mut chemistry, &mut seen, "inert component", id);
                }
            }
            UnitSpec::Cstr { ref kinetics, .. }
            | UnitSpec::PlugFlow { ref kinetics, .. }
            | UnitSpec::Batch { ref kinetics, .. } => kinetics_issues(&mut chemistry, kinetics),
            _ => {}
        }
        match *self {
            UnitSpec::Cstr {
                thermal, volume_m3, ..
            } => {
                opt("volume", volume_m3, Sign::Positive);
                match thermal {
                    SeparatorThermal::Adiabatic => {}
                    SeparatorThermal::OutletTemperature(t) => {
                        opt("outlet temperature", t, Sign::Positive)
                    }
                    SeparatorThermal::Duty(q) => opt("duty", q, Sign::Any),
                }
            }
            UnitSpec::PlugFlow {
                length_m,
                diameter_m,
                thermal,
                ..
            } => {
                opt("length", length_m, Sign::Positive);
                opt("diameter", diameter_m, Sign::Positive);
                match thermal {
                    PlugFlowThermal::Adiabatic => {}
                    PlugFlowThermal::Isothermal(t) => opt("temperature", t, Sign::Positive),
                    PlugFlowThermal::Coolant {
                        temperature_k,
                        u_w_m2_k,
                    } => {
                        opt("coolant temperature", temperature_k, Sign::Positive);
                        opt("heat-transfer coefficient", u_w_m2_k, Sign::NonNegative);
                    }
                }
            }
            UnitSpec::Batch {
                reaction_time_s,
                thermal,
                ..
            } => {
                opt("reaction time", reaction_time_s, Sign::Positive);
                if let BatchThermal::Isothermal(t) = thermal {
                    opt("temperature", t, Sign::Positive);
                }
            }
            UnitSpec::ConversionReactor { thermal, .. }
            | UnitSpec::YieldReactor { thermal, .. }
            | UnitSpec::EquilibriumReactor { thermal, .. }
            | UnitSpec::GibbsReactor { thermal, .. } => match thermal {
                SeparatorThermal::Adiabatic => {}
                SeparatorThermal::OutletTemperature(t) => {
                    opt("outlet temperature", t, Sign::Positive)
                }
                SeparatorThermal::Duty(q) => opt("duty", q, Sign::Any),
            },
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
            }
            | UnitSpec::ConversionReactor {
                pressure_drop_pa, ..
            }
            | UnitSpec::YieldReactor {
                pressure_drop_pa, ..
            }
            | UnitSpec::Cstr {
                pressure_drop_pa, ..
            }
            | UnitSpec::PlugFlow {
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
            UnitSpec::EquilibriumReactor {
                equilibrium,
                pressure_drop_pa,
                ..
            }
            | UnitSpec::GibbsReactor {
                equilibrium,
                pressure_drop_pa,
                ..
            } => {
                check_value(
                    &mut fixed,
                    "temperature approach",
                    equilibrium.temperature_approach_k,
                    Sign::Any,
                );
                check_value(
                    &mut fixed,
                    "pressure drop",
                    pressure_drop_pa,
                    Sign::NonNegative,
                );
            }
            _ => {}
        }
        out.extend(chemistry);
        out.extend(fixed);
        out
    }
}

fn component_issue(
    out: &mut Vec<SpecIssue>,
    seen: &mut Vec<String>,
    field: &'static str,
    id: &str,
) {
    let problem = if id.trim().is_empty() {
        SpecProblem::Missing
    } else if seen.iter().any(|s| s == id) {
        SpecProblem::Duplicate
    } else {
        seen.push(id.to_owned());
        return;
    };
    out.push(SpecIssue { field, problem });
}

/// Components, coefficients and reactant/product presence of one reaction.
fn coefficient_issues(out: &mut Vec<SpecIssue>, coefficients: &[(String, f64)]) {
    let mut seen = Vec::new();
    for (id, nu) in coefficients {
        component_issue(out, &mut seen, "reaction component", id);
        if !nu.is_finite() {
            check_value(out, "stoichiometric coefficient", *nu, Sign::Any);
        } else if *nu == 0.0 {
            out.push(SpecIssue {
                field: "stoichiometric coefficient",
                problem: SpecProblem::Zero,
            });
        }
    }
    let has = |f: fn(f64) -> bool| coefficients.iter().any(|(_, nu)| nu.is_finite() && f(*nu));
    if !has(|nu| nu < 0.0) || !has(|nu| nu > 0.0) {
        out.push(SpecIssue {
            field: "reaction",
            problem: SpecProblem::NeedsReactantAndProduct,
        });
    }
}

fn order_issues(out: &mut Vec<SpecIssue>, orders: &[(String, f64)]) {
    let mut seen = Vec::new();
    for (id, a) in orders {
        component_issue(out, &mut seen, "order component", id);
        check_value(out, "order", *a, Sign::Any);
    }
}

fn arrhenius_issues(out: &mut Vec<SpecIssue>, k: &Arrhenius) {
    match k.pre_exponential {
        None => out.push(SpecIssue {
            field: "pre-exponential factor",
            problem: SpecProblem::Missing,
        }),
        Some(v) => check_value(out, "pre-exponential factor", v, Sign::NonNegative),
    }
    check_value(
        out,
        "activation energy",
        k.activation_energy_j_mol,
        Sign::Any,
    );
    check_value(
        out,
        "temperature exponent",
        k.temperature_exponent,
        Sign::Any,
    );
    if let Some(t) = k.reference_temperature_k {
        check_value(out, "reference temperature", t, Sign::Positive);
    }
}

fn lhhw_issues(out: &mut Vec<SpecIssue>, term: &LhhwTerm) {
    let c = term.constant;
    for v in [c.a, c.b, c.c, c.d] {
        check_value(out, "LHHW constant", v, Sign::Any);
    }
    order_issues(out, &term.orders);
}

fn kinetics_issues(out: &mut Vec<SpecIssue>, kinetics: &Kinetics) {
    if kinetics.reactions.is_empty() {
        out.push(SpecIssue {
            field: "reactions",
            problem: SpecProblem::Missing,
        });
    }
    for r in &kinetics.reactions {
        coefficient_issues(out, &r.coefficients);
        match &r.rate {
            RateLaw::PowerLaw { forward, reverse } => {
                for term in std::iter::once(forward).chain(reverse) {
                    arrhenius_issues(out, &term.constant);
                    order_issues(out, &term.orders);
                }
            }
            RateLaw::Lhhw {
                kinetic,
                driving_forward,
                driving_reverse,
                adsorption,
                exponent,
            } => {
                arrhenius_issues(out, kinetic);
                for term in std::iter::once(driving_forward)
                    .chain(driving_reverse)
                    .chain(adsorption)
                {
                    lhhw_issues(out, term);
                }
                check_value(out, "adsorption exponent", *exponent, Sign::NonNegative);
            }
        }
    }
    let catalytic = kinetics
        .reactions
        .iter()
        .any(|r| r.rate_basis == RateBasis::CatalystMass);
    match kinetics.catalyst_bulk_density_kg_m3 {
        None if catalytic => out.push(SpecIssue {
            field: "catalyst bulk density",
            problem: SpecProblem::Missing,
        }),
        Some(rho) => check_value(out, "catalyst bulk density", rho, Sign::Positive),
        None => {}
    }
}

fn reaction_issues(out: &mut Vec<SpecIssue>, reactions: &[Reaction]) {
    if reactions.is_empty() {
        out.push(SpecIssue {
            field: "reactions",
            problem: SpecProblem::Missing,
        });
    }
    for r in reactions {
        coefficient_issues(out, &r.coefficients);
        match &r.extent {
            ReactionExtent::Conversion { key, fraction } => {
                match key {
                    None => out.push(SpecIssue {
                        field: "key component",
                        problem: SpecProblem::Missing,
                    }),
                    Some(key) if !r.coefficients.iter().any(|(id, nu)| id == key && *nu < 0.0) => {
                        out.push(SpecIssue {
                            field: "key component",
                            problem: SpecProblem::KeyNotReactant,
                        })
                    }
                    Some(_) => {}
                }
                match fraction {
                    None => out.push(SpecIssue {
                        field: "conversion",
                        problem: SpecProblem::Missing,
                    }),
                    Some(x) => check_value(out, "conversion", *x, Sign::Fraction),
                }
            }
            ReactionExtent::MolarExtent(None) => out.push(SpecIssue {
                field: "extent",
                problem: SpecProblem::Missing,
            }),
            ReactionExtent::MolarExtent(Some(e)) => {
                check_value(out, "extent", *e, Sign::NonNegative)
            }
        }
    }
}

fn yield_issues(out: &mut Vec<SpecIssue>, yields: &[(String, f64)], inerts: &[String]) {
    if yields.is_empty() {
        out.push(SpecIssue {
            field: "mass yields",
            problem: SpecProblem::Missing,
        });
    }
    let mut seen = Vec::new();
    for (id, y) in yields {
        component_issue(out, &mut seen, "yield component", id);
        check_value(out, "mass yield", *y, Sign::NonNegative);
    }
    if !yields.is_empty() && yields.iter().all(|(_, y)| *y == 0.0) {
        out.push(SpecIssue {
            field: "mass yields",
            problem: SpecProblem::AllZero,
        });
    }
    for id in inerts {
        component_issue(out, &mut seen, "inert component", id);
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

impl Unit {
    /// [`UnitSpec::issues`] plus checks that depend on the unit kind: a heater
    /// with a negative duty or a cooler with a positive duty is contradictory
    /// (proposal section 2.4) and is reported before solving. Direction checks
    /// that need the inlet state (an outlet temperature on the wrong side of the
    /// inlet temperature) are reported as unit warnings by the engine.
    pub fn issues(&self) -> Vec<SpecIssue> {
        let mut out = self.spec.issues();
        if let UnitSpec::HeatExchange {
            thermal: ThermalSpec::Duty(Some(q)),
            ..
        } = self.spec
        {
            let contradictory = match self.kind {
                UnitKind::Heater => q < 0.0,
                UnitKind::Cooler => q > 0.0,
                _ => false,
            };
            if contradictory {
                out.push(SpecIssue {
                    field: "duty",
                    problem: SpecProblem::ContradictsUnitKind,
                });
            }
        }
        out
    }
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

    #[test]
    fn duty_sign_must_match_heater_or_cooler() {
        let unit = |kind, q| Unit {
            id: UnitId(1),
            kind,
            name: "E-100".into(),
            spec: UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(q)),
                pressure_drop_pa: 0.0,
            },
        };
        let contradiction = [SpecIssue {
            field: "duty",
            problem: SpecProblem::ContradictsUnitKind,
        }];
        assert_eq!(unit(UnitKind::Heater, -1.0).issues(), contradiction);
        assert_eq!(unit(UnitKind::Cooler, 1.0).issues(), contradiction);
        for (kind, q) in [
            (UnitKind::Heater, 1.0),
            (UnitKind::Heater, 0.0),
            (UnitKind::Cooler, -1.0),
            (UnitKind::Cooler, 0.0),
        ] {
            assert!(unit(kind, q).issues().is_empty(), "{kind} {q}");
        }
    }
}
