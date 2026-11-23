//! Unit operation kinds, ports, and user specifications.
//!
//! Specification sets follow section 2.3 of `docs/flowsheet-and-ui.md`.
//! Mutually exclusive choices (for example outlet temperature *or* duty) are
//! modelled as enums, so contradictory combinations cannot be represented.
//!
//! The kinds added by work package A17 belong to groups whose metadata,
//! specifications and calculations live in [`crate::units`]; see
//! `docs/a17-architecture.md`. This file holds the shared registry (the
//! [`UnitKind`] enum, the [`UnitSpec`] dispatch) and the core kinds.

use std::fmt;

use crate::units::{
    self, ColumnSpec, GroupSpec, HeatTransportSpec, MultiphaseSpec, ReactorSpec, RotatingSpec,
    SimpleLogicalSpec, SolidsSpec,
};

/// Stable identifier of a unit operation within one [`crate::Flowsheet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UnitId(pub u64);

impl fmt::Display for UnitId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "U{}", self.0)
    }
}

/// Whether material or energy enters or leaves a unit through a port.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PortDirection {
    /// Material (or energy) enters the unit.
    Inlet,
    /// Material (or energy) leaves the unit.
    Outlet,
}

/// What a port, and therefore the stream connected to it, carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PortMedium {
    /// A material stream (flow, composition, state).
    Material,
    /// An energy stream: a heat duty or shaft work in W (#171).
    Energy,
}

/// Static description of one port of a unit.
///
/// Ports are identified within a unit by `(name, number, direction, medium)`;
/// numbered ports (`number: Some(n)`) belong to a variable group such as the
/// feeds or side draws of a column, whose size is part of the specification
/// (see [`Unit::ports`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PortSpec {
    /// Short, user-facing base name, e.g. `in`, `feed`, `duty`.
    pub name: &'static str,
    /// Position within a variable port group (1-based), if numbered.
    pub number: Option<u16>,
    /// Flow direction relative to the unit.
    pub direction: PortDirection,
    /// Material or energy.
    pub medium: PortMedium,
    /// An optional port may stay unconnected without a diagnostic. Every
    /// energy port of the A17 scaffolding is optional.
    pub optional: bool,
}

impl PortSpec {
    /// Required material inlet.
    pub const fn inlet(name: &'static str) -> Self {
        Self {
            name,
            number: None,
            direction: PortDirection::Inlet,
            medium: PortMedium::Material,
            optional: false,
        }
    }

    /// Required material outlet.
    pub const fn outlet(name: &'static str) -> Self {
        Self {
            direction: PortDirection::Outlet,
            ..Self::inlet(name)
        }
    }

    /// Optional energy inlet (energy enters the unit).
    pub const fn energy_inlet(name: &'static str) -> Self {
        Self {
            medium: PortMedium::Energy,
            optional: true,
            ..Self::inlet(name)
        }
    }

    /// Optional energy outlet (energy leaves the unit).
    pub const fn energy_outlet(name: &'static str) -> Self {
        Self {
            direction: PortDirection::Outlet,
            ..Self::energy_inlet(name)
        }
    }

    /// The same port as member `n` (1-based) of a numbered group.
    pub const fn numbered(self, n: u16) -> Self {
        Self {
            number: Some(n),
            ..self
        }
    }

    /// The same port, optional.
    pub const fn optional(self) -> Self {
        Self {
            optional: true,
            ..self
        }
    }

    /// User-facing label: the name, followed by the number if numbered.
    pub fn label(&self) -> String {
        match self.number {
            Some(n) => format!("{} {n}", self.name),
            None => self.name.to_string(),
        }
    }

    /// Whether `other` denotes the same port (used to keep connections when
    /// a variable port group changes size).
    pub fn same_port(&self, other: &PortSpec) -> bool {
        self.name == other.name
            && self.number == other.number
            && self.direction == other.direction
            && self.medium == other.medium
    }
}

impl fmt::Display for PortSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.label())
    }
}

const fn inlet(name: &'static str) -> PortSpec {
    PortSpec::inlet(name)
}

const fn outlet(name: &'static str) -> PortSpec {
    PortSpec::outlet(name)
}

/// The A17 work-package group that owns a kind's specification, project
/// record, form and calculation (`docs/a17-architecture.md`). `Core` kinds
/// are maintained by the integrator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum UnitGroup {
    /// Boundaries, energy streams and the original A8 units.
    Core,
    /// Pumps, compressors, expanders, multistage compressor (A14, #182).
    Rotating,
    /// Air cooler, fired heater, multi-stream exchanger, pipe, relief valve.
    HeatTransport,
    /// Component splitter, manipulators, tank and logical operations.
    SimpleLogical,
    /// Reactors (#191–#197).
    Reactors,
    /// Shortcut and rigorous columns (#185–#187, #190).
    Columns,
    /// VLLE separators and three-phase / liquid-liquid columns.
    Multiphase,
    /// Solids handling (#198–#211).
    Solids,
}

/// Palette section of a kind (`docs/flowsheet-and-ui.md`, section 3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum PaletteCategory {
    /// Material and energy stream boundaries.
    Streams,
    /// Mixers, splitters and stream manipulators.
    MixingSplitting,
    /// Separators and vessels.
    Separators,
    /// Heat transfer equipment.
    HeatTransfer,
    /// Pressure change and transport.
    PressureChange,
    /// Columns.
    Columns,
    /// Reactors.
    Reactors,
    /// Solids handling.
    Solids,
    /// Logical operations and flowsheet control.
    Logical,
    /// User and CAPE-OPEN units.
    Extensibility,
}

impl PaletteCategory {
    /// All categories in palette order.
    pub const ALL: [PaletteCategory; 10] = [
        PaletteCategory::Streams,
        PaletteCategory::MixingSplitting,
        PaletteCategory::Separators,
        PaletteCategory::HeatTransfer,
        PaletteCategory::PressureChange,
        PaletteCategory::Columns,
        PaletteCategory::Reactors,
        PaletteCategory::Solids,
        PaletteCategory::Logical,
        PaletteCategory::Extensibility,
    ];

    /// User-facing heading.
    pub fn label(self) -> &'static str {
        match self {
            PaletteCategory::Streams => "Streams",
            PaletteCategory::MixingSplitting => "Mixing/splitting",
            PaletteCategory::Separators => "Separators",
            PaletteCategory::HeatTransfer => "Heat transfer",
            PaletteCategory::PressureChange => "Pressure change",
            PaletteCategory::Columns => "Columns",
            PaletteCategory::Reactors => "Reactors",
            PaletteCategory::Solids => "Solids",
            PaletteCategory::Logical => "Logical",
            PaletteCategory::Extensibility => "Extensibility",
        }
    }

    /// Kinds of this category in palette order.
    pub fn kinds(self) -> impl Iterator<Item = UnitKind> {
        UnitKind::ALL
            .into_iter()
            .filter(move |k| k.category() == self)
    }
}

/// Static metadata of a unit kind. Group kinds supply it from their module
/// in [`crate::units`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KindInfo {
    /// User-facing name.
    pub label: &'static str,
    /// Tag series shared by kinds numbered together (`E`, `R`, `T`, …).
    pub tag_series: &'static str,
    /// Palette section.
    pub category: PaletteCategory,
    /// Fixed material ports, in port-index order.
    pub ports: &'static [PortSpec],
    /// Fixed energy ports; they follow the material ports in index order.
    pub energy_ports: &'static [PortSpec],
    /// GitHub sub-issue of A17 that implements the kind, if any.
    pub issue: Option<u16>,
}

/// Whether the engine calculates a kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalculationStatus {
    /// Calculated by the engine.
    Calculated,
    /// Drawable and specifiable only; a solve containing it stops in
    /// preflight with this reason (completing "<kind> …").
    DrawnOnly(&'static str),
}

/// Unit operations, plus flowsheet boundaries for material and energy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum UnitKind {
    /// Boundary source of a material stream.
    Feed,
    /// Boundary sink of a material stream.
    Product,
    /// Boundary source of an energy stream: a specified duty, or the duty
    /// calculated by the unit it is connected to (#171).
    EnergySource,
    /// Boundary sink of an energy stream (#171).
    EnergySink,
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
    // --- A17 kinds; metadata in `crate::units::<group>` ---
    /// Component splitter (#172).
    ComponentSplitter,
    /// Stream multiplier / duplicator (#173).
    Multiplier,
    /// Stream selector (#174).
    Selector,
    /// Property-package / stream cutter (#175).
    StreamCutter,
    /// Three-phase separator (#176).
    ThreePhaseSeparator,
    /// Liquid–liquid decanter (#177).
    Decanter,
    /// Tank / storage vessel (#178).
    Tank,
    /// Air cooler (#179).
    AirCooler,
    /// Fired heater / furnace (#180).
    FiredHeater,
    /// Multi-stream (LNG / plate-fin) exchanger (#181).
    MultiStreamExchanger,
    /// Multistage compressor with intercooling (#182).
    MultistageCompressor,
    /// Pipe segment / pipeline (#183).
    PipeSegment,
    /// Relief valve (#184).
    ReliefValve,
    /// Shortcut distillation (#185).
    ShortcutColumn,
    /// Rigorous distillation column (#186).
    DistillationColumn,
    /// Absorber / stripper (#187).
    Absorber,
    /// Three-phase distillation (#188).
    ThreePhaseColumn,
    /// Liquid–liquid extractor (#189).
    LiquidExtractor,
    /// Petroleum / multi-column fractionator (#190).
    PetroleumFractionator,
    /// Conversion (stoichiometric) reactor (#191).
    ConversionReactor,
    /// Yield reactor (#192).
    YieldReactor,
    /// Equilibrium reactor (#193).
    EquilibriumReactor,
    /// Gibbs reactor (#194).
    GibbsReactor,
    /// Continuous stirred-tank reactor (#195).
    Cstr,
    /// Plug-flow reactor (#196).
    PlugFlowReactor,
    /// Batch reactor (#197).
    BatchReactor,
    /// Cyclone (#198).
    Cyclone,
    /// Hydrocyclone (#199).
    Hydrocyclone,
    /// Baghouse / fabric filter (#200).
    Baghouse,
    /// Rotary vacuum / cake filter (#201).
    RotaryFilter,
    /// Simple solids separator (#202).
    SolidsSeparator,
    /// Crystallizer (#203).
    Crystallizer,
    /// Centrifuge (#204).
    Centrifuge,
    /// Electrostatic precipitator (#205).
    ElectrostaticPrecipitator,
    /// Venturi scrubber (#206).
    VenturiScrubber,
    /// Crusher (#207).
    Crusher,
    /// Screen / classifier (#208).
    Screen,
    /// Dryer (#209).
    Dryer,
    /// Granulator / fluidised bed (#210).
    Granulator,
    /// Solids washer / counter-current decantation (#211).
    SolidsWasher,
    /// Adjust / design specification (#212).
    Adjust,
    /// Set (linked variable) (#213).
    Set,
    /// Spreadsheet / calculator (#214).
    Spreadsheet,
    /// Balance (#215).
    Balance,
    /// Explicit recycle block (#216).
    RecycleBlock,
    /// Control blocks: PID, ratio, surge, transfer function (#217).
    ControlBlock,
    /// User / CAPE-OPEN unit operation (#218).
    UserUnit,
}

impl UnitKind {
    /// All kinds, in palette order (grouped by [`PaletteCategory`]).
    pub const ALL: [UnitKind; 61] = [
        // Streams
        UnitKind::Feed,
        UnitKind::Product,
        UnitKind::EnergySource,
        UnitKind::EnergySink,
        // Mixing/splitting
        UnitKind::Mixer,
        UnitKind::Splitter,
        UnitKind::ComponentSplitter,
        UnitKind::Multiplier,
        UnitKind::Selector,
        UnitKind::StreamCutter,
        // Separators
        UnitKind::Separator,
        UnitKind::ThreePhaseSeparator,
        UnitKind::Decanter,
        UnitKind::Tank,
        // Heat transfer
        UnitKind::Heater,
        UnitKind::Cooler,
        UnitKind::HeatExchanger,
        UnitKind::AirCooler,
        UnitKind::FiredHeater,
        UnitKind::MultiStreamExchanger,
        // Pressure change
        UnitKind::Valve,
        UnitKind::Pump,
        UnitKind::Compressor,
        UnitKind::Expander,
        UnitKind::MultistageCompressor,
        UnitKind::PipeSegment,
        UnitKind::ReliefValve,
        // Columns
        UnitKind::ShortcutColumn,
        UnitKind::DistillationColumn,
        UnitKind::Absorber,
        UnitKind::ThreePhaseColumn,
        UnitKind::LiquidExtractor,
        UnitKind::PetroleumFractionator,
        // Reactors
        UnitKind::ConversionReactor,
        UnitKind::YieldReactor,
        UnitKind::EquilibriumReactor,
        UnitKind::GibbsReactor,
        UnitKind::Cstr,
        UnitKind::PlugFlowReactor,
        UnitKind::BatchReactor,
        // Solids
        UnitKind::Cyclone,
        UnitKind::Hydrocyclone,
        UnitKind::Baghouse,
        UnitKind::RotaryFilter,
        UnitKind::SolidsSeparator,
        UnitKind::Crystallizer,
        UnitKind::Centrifuge,
        UnitKind::ElectrostaticPrecipitator,
        UnitKind::VenturiScrubber,
        UnitKind::Crusher,
        UnitKind::Screen,
        UnitKind::Dryer,
        UnitKind::Granulator,
        UnitKind::SolidsWasher,
        // Logical
        UnitKind::Adjust,
        UnitKind::Set,
        UnitKind::Spreadsheet,
        UnitKind::Balance,
        UnitKind::RecycleBlock,
        UnitKind::ControlBlock,
        // Extensibility
        UnitKind::UserUnit,
    ];

    /// The work-package group that owns this kind.
    pub fn group(self) -> UnitGroup {
        use UnitKind::*;
        match self {
            Feed | Product | EnergySource | EnergySink | Mixer | Heater | Cooler | Valve
            | Separator | Splitter | HeatExchanger => UnitGroup::Core,
            Pump | Compressor | Expander | MultistageCompressor => UnitGroup::Rotating,
            AirCooler | FiredHeater | MultiStreamExchanger | PipeSegment | ReliefValve => {
                UnitGroup::HeatTransport
            }
            ComponentSplitter | Multiplier | Selector | StreamCutter | Tank | Adjust | Set
            | Spreadsheet | Balance | RecycleBlock | ControlBlock | UserUnit => {
                UnitGroup::SimpleLogical
            }
            ConversionReactor | YieldReactor | EquilibriumReactor | GibbsReactor | Cstr
            | PlugFlowReactor | BatchReactor => UnitGroup::Reactors,
            ShortcutColumn | DistillationColumn | Absorber | PetroleumFractionator => {
                UnitGroup::Columns
            }
            ThreePhaseSeparator | Decanter | ThreePhaseColumn | LiquidExtractor => {
                UnitGroup::Multiphase
            }
            Cyclone
            | Hydrocyclone
            | Baghouse
            | RotaryFilter
            | SolidsSeparator
            | Crystallizer
            | Centrifuge
            | ElectrostaticPrecipitator
            | VenturiScrubber
            | Crusher
            | Screen
            | Dryer
            | Granulator
            | SolidsWasher => UnitGroup::Solids,
        }
    }

    /// Static metadata (label, tags, palette category, fixed ports).
    pub fn info(self) -> KindInfo {
        match self.group() {
            UnitGroup::Core => core_info(self),
            UnitGroup::Rotating => units::rotating::info(self),
            UnitGroup::HeatTransport => units::heat_transport::info(self),
            UnitGroup::SimpleLogical => units::simple_logical::info(self),
            UnitGroup::Reactors => units::reactors::info(self),
            UnitGroup::Columns => units::columns::info(self),
            UnitGroup::Multiphase => units::multiphase::info(self),
            UnitGroup::Solids => units::solids::info(self),
        }
    }

    /// Whether the engine calculates this kind, or why not.
    pub fn calculation_status(self) -> CalculationStatus {
        match self.group() {
            UnitGroup::Core => CalculationStatus::Calculated,
            UnitGroup::Rotating => units::rotating::status(self),
            UnitGroup::HeatTransport => units::heat_transport::status(self),
            UnitGroup::SimpleLogical => units::simple_logical::status(self),
            UnitGroup::Reactors => units::reactors::status(self),
            UnitGroup::Columns => units::columns::status(self),
            UnitGroup::Multiphase => units::multiphase::status(self),
            UnitGroup::Solids => units::solids::status(self),
        }
    }

    /// Whether the steady-state engine calculates this kind. Pumps,
    /// compressors, and expanders need a PS flash and efficiency models and
    /// are deferred (A14); most A17 kinds are drawn only until their group
    /// implements them. They can be drawn and specified but not yet solved.
    pub fn calculation_planned(self) -> bool {
        self.calculation_status() == CalculationStatus::Calculated
    }

    /// User-facing name of the kind.
    pub fn label(self) -> &'static str {
        self.info().label
    }

    /// Palette section of the kind.
    pub fn category(self) -> PaletteCategory {
        self.info().category
    }

    /// A17 sub-issue implementing the kind, if any.
    pub fn issue(self) -> Option<u16> {
        self.info().issue
    }

    /// Whether the kind is a flowsheet boundary (no internal balance).
    pub fn is_boundary(self) -> bool {
        matches!(
            self,
            UnitKind::Feed | UnitKind::Product | UnitKind::EnergySource | UnitKind::EnergySink
        )
    }

    /// Tag series shared by kinds that are numbered together; heaters,
    /// coolers, and heat exchangers are all `E-…`.
    pub(crate) fn tag_series(self) -> &'static str {
        self.info().tag_series
    }

    /// Generated tag for the `n`-th unit (1-based) of this kind's series,
    /// following common PFD conventions: `Feed 1`, `M-100`, `E-100`,
    /// `V-100`, `D-100`, `P-100`, `K-100`, `Product 1`.
    pub(crate) fn tag(self, n: u64) -> String {
        if self.is_boundary() {
            format!("{} {n}", self.tag_series())
        } else {
            format!("{}-{}", self.tag_series(), 99 + n)
        }
    }

    /// Inverse of [`Self::tag`]: the series number of a generated tag, or
    /// `None` for a user-chosen name.
    pub(crate) fn tag_number(self, name: &str) -> Option<u64> {
        let rest = name.strip_prefix(self.tag_series())?;
        let n = if self.is_boundary() {
            rest.strip_prefix(' ')?.parse::<u64>().ok()?
        } else {
            rest.strip_prefix('-')?
                .parse::<u64>()
                .ok()?
                .checked_sub(99)?
        };
        (n > 0 && self.tag(n) == name).then_some(n)
    }

    /// Fixed **material** ports of this kind. Port indices in [`PortRef`]
    /// refer to positions in [`Unit::ports`], which starts with this slice,
    /// followed by [`Self::energy_ports`] and the variable ports of the
    /// specification.
    pub fn ports(self) -> &'static [PortSpec] {
        self.info().ports
    }

    /// Fixed (optional) energy ports, after the material ports.
    pub fn energy_ports(self) -> &'static [PortSpec] {
        self.info().energy_ports
    }

    /// Every port of a new unit of this kind (default specification).
    pub fn default_ports(self) -> Vec<PortSpec> {
        ports_for(self, &self.default_spec())
    }

    /// Default specification for a newly created unit of this kind
    /// (defaults column of the proposal's section 2.3).
    pub fn default_spec(self) -> UnitSpec {
        match self.group() {
            UnitGroup::Core => core_default_spec(self),
            UnitGroup::Rotating => match self {
                UnitKind::Pump | UnitKind::Compressor => UnitSpec::Pressurize {
                    pressure: PressureRise::OutletPressure(None),
                    efficiency: DEFAULT_EFFICIENCY,
                },
                UnitKind::Expander => UnitSpec::Expand {
                    pressure: PressureSpec::OutletPressure(None),
                    efficiency: DEFAULT_EFFICIENCY,
                },
                _ => UnitSpec::Rotating(units::rotating::default_spec(self)),
            },
            UnitGroup::HeatTransport => {
                UnitSpec::HeatTransport(units::heat_transport::default_spec(self))
            }
            UnitGroup::SimpleLogical => {
                UnitSpec::SimpleLogical(units::simple_logical::default_spec(self))
            }
            UnitGroup::Reactors => UnitSpec::Reactors(units::reactors::default_spec(self)),
            UnitGroup::Columns => UnitSpec::Columns(units::columns::default_spec(self)),
            UnitGroup::Multiphase => UnitSpec::Multiphase(units::multiphase::default_spec(self)),
            UnitGroup::Solids => UnitSpec::Solids(units::solids::default_spec(self)),
        }
    }
}

fn core_info(kind: UnitKind) -> KindInfo {
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
    // Energy boundaries carry one required energy port; heater and cooler
    // expose an optional one in the direction the heat flows.
    const ENERGY_OUT: [PortSpec; 1] = [PortSpec {
        optional: false,
        ..PortSpec::energy_outlet("Q")
    }];
    const ENERGY_IN: [PortSpec; 1] = [PortSpec {
        optional: false,
        ..PortSpec::energy_inlet("Q")
    }];
    const HEAT_IN: [PortSpec; 1] = [PortSpec::energy_inlet("Q")];
    const HEAT_OUT: [PortSpec; 1] = [PortSpec::energy_outlet("Q")];
    let info = |label, tag_series, category, ports, energy_ports| KindInfo {
        label,
        tag_series,
        category,
        ports,
        energy_ports,
        issue: None,
    };
    use PaletteCategory as C;
    match kind {
        UnitKind::Feed => info("Feed", "Feed", C::Streams, &FEED, &[]),
        UnitKind::Product => info("Product", "Product", C::Streams, &PRODUCT, &[]),
        UnitKind::EnergySource => KindInfo {
            issue: Some(171),
            ..info("Energy source", "Q-in", C::Streams, &[], &ENERGY_OUT)
        },
        UnitKind::EnergySink => KindInfo {
            issue: Some(171),
            ..info("Energy sink", "Q-out", C::Streams, &[], &ENERGY_IN)
        },
        UnitKind::Mixer => info("Mixer", "M", C::MixingSplitting, &MIXER, &[]),
        UnitKind::Splitter => info("Splitter", "TEE", C::MixingSplitting, &SPLITTER, &[]),
        UnitKind::Heater => info("Heater", "E", C::HeatTransfer, &SINGLE, &HEAT_IN),
        UnitKind::Cooler => info("Cooler", "E", C::HeatTransfer, &SINGLE, &HEAT_OUT),
        UnitKind::HeatExchanger => info("Heat exchanger", "E", C::HeatTransfer, &EXCHANGER, &[]),
        UnitKind::Valve => info("Valve", "V", C::PressureChange, &SINGLE, &[]),
        UnitKind::Separator => info("Separator", "D", C::Separators, &SEPARATOR, &[]),
        _ => unreachable!("{kind:?} is not a core kind"),
    }
}

fn core_default_spec(kind: UnitKind) -> UnitSpec {
    match kind {
        UnitKind::Feed => UnitSpec::Feed {
            temperature_k: None,
            pressure_pa: None,
            molar_flow_mol_s: None,
        },
        UnitKind::Product => UnitSpec::Product,
        UnitKind::EnergySource => UnitSpec::EnergySource {
            duty: EnergySourceDuty::Calculated,
        },
        UnitKind::EnergySink => UnitSpec::EnergySink,
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
        _ => unreachable!("{kind:?} is not a core kind"),
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
    /// Duty taken from the energy stream connected to the unit's energy
    /// inlet (heaters only; #171).
    EnergyStream,
}

/// Thermal specification for a separator or steady-state tank.
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
    /// Hot inlet temperature minus cold outlet temperature, in K.
    HotEndApproach(Option<f64>),
    /// Hot outlet temperature minus cold inlet temperature, in K.
    ColdEndApproach(Option<f64>),
}

/// Duty of an energy source boundary (#171).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnergySourceDuty {
    /// The connected unit calculates the duty from its own specification;
    /// the source reports it (default).
    Calculated,
    /// Duty in W delivered to the connected unit, which must take its duty
    /// from the energy stream.
    Specified(Option<f64>),
}

/// How a unit treats the energy stream at one of its energy ports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyRole {
    /// The unit determines the duty (from its own specification).
    Determines,
    /// The unit takes the duty from the stream.
    Consumes,
}

/// User-entered specification of a unit, in SI units with absolute pressure.
///
/// `None` means "not yet specified". Feed compositions are supplied separately
/// as [`crate::FeedCompositions`], with explicit identities and mole basis.
/// The group variants wrap a specification type owned by that group's module
/// in [`crate::units`].
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
    /// Energy stream source boundary.
    EnergySource {
        /// Calculated by the connected unit, or specified.
        duty: EnergySourceDuty,
    },
    /// Energy stream sink boundary; no specification.
    EnergySink,
    /// Adiabatic mixer.
    Mixer {
        /// Outlet pressure choice.
        outlet_pressure: MixerPressure,
    },
    /// Heater or cooler.
    HeatExchange {
        /// Outlet temperature, duty, or duty from the energy stream.
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
    /// Pump or compressor (calculation owned by the rotating group).
    Pressurize {
        /// Outlet pressure or pressure increase.
        pressure: PressureRise,
        /// Efficiency, greater than 0 and at most 1.
        efficiency: f64,
    },
    /// Expander (calculation owned by the rotating group).
    Expand {
        /// Outlet pressure or pressure drop.
        pressure: PressureSpec,
        /// Isentropic efficiency, greater than 0 and at most 1.
        efficiency: f64,
    },
    /// Rotating group kinds other than pump, compressor and expander.
    Rotating(RotatingSpec),
    /// Heat-transport group.
    HeatTransport(HeatTransportSpec),
    /// Simple and logical group.
    SimpleLogical(SimpleLogicalSpec),
    /// Reactors group.
    Reactors(ReactorSpec),
    /// Columns group.
    Columns(ColumnSpec),
    /// Multiphase group.
    Multiphase(MultiphaseSpec),
    /// Solids group.
    Solids(SolidsSpec),
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
    /// The choice is not available for this unit kind (for example a cooler
    /// taking its duty from an energy stream; its energy port is an outlet).
    NotAvailableForKind,
    /// A count is outside its admissible range (for example zero column feeds).
    CountOutOfRange,
    /// A value lies outside a model-specific admissible range (for example a
    /// key recovery of exactly 1 or a multiple of minimum reflux of 1 or less).
    OutOfRange,
    /// The selection repeats another selection (for example equal light and
    /// heavy keys).
    Duplicate,
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
            SpecProblem::NotAvailableForKind => "is not available for this unit kind",
            SpecProblem::CountOutOfRange => "is outside the admissible range",
            SpecProblem::OutOfRange => "is outside the admissible range",
            SpecProblem::Duplicate => "duplicates another selection",
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

/// Admissible range of a specification value, for [`check_value`].
#[derive(Clone, Copy)]
pub enum Sign {
    /// Strictly positive.
    Positive,
    /// Zero or positive.
    NonNegative,
    /// 0 to 1.
    Fraction,
    /// Greater than 0, at most 1.
    Efficiency,
    /// Any finite value.
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
            | UnitSpec::EnergySink
            | UnitSpec::EnergySource {
                duty: EnergySourceDuty::Calculated,
            }
            | UnitSpec::Mixer {
                outlet_pressure: MixerPressure::LowestInlet,
            }
            | UnitSpec::Separator {
                thermal: SeparatorThermal::Adiabatic,
                ..
            } => {}
            UnitSpec::EnergySource {
                duty: EnergySourceDuty::Specified(q),
            } => opt("duty", q, Sign::Any),
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
                ThermalSpec::EnergyStream => {}
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
                ExchangerSpec::HotEndApproach(t) => opt("hot-end approach", t, Sign::NonNegative),
                ExchangerSpec::ColdEndApproach(t) => opt("cold-end approach", t, Sign::NonNegative),
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
            UnitSpec::Rotating(ref s) => out.extend(s.issues()),
            UnitSpec::HeatTransport(ref s) => out.extend(s.issues()),
            UnitSpec::SimpleLogical(ref s) => out.extend(s.issues()),
            UnitSpec::Reactors(ref s) => out.extend(s.issues()),
            UnitSpec::Columns(ref s) => out.extend(s.issues()),
            UnitSpec::Multiphase(ref s) => out.extend(s.issues()),
            UnitSpec::Solids(ref s) => out.extend(s.issues()),
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

    /// Whether this specification belongs to `kind`. [`crate::Flowsheet`]
    /// rejects a specification of another kind.
    pub fn fits(&self, kind: UnitKind) -> bool {
        use UnitKind as K;
        match self {
            UnitSpec::Feed { .. } => kind == K::Feed,
            UnitSpec::Product => kind == K::Product,
            UnitSpec::EnergySource { .. } => kind == K::EnergySource,
            UnitSpec::EnergySink => kind == K::EnergySink,
            UnitSpec::Mixer { .. } => kind == K::Mixer,
            UnitSpec::HeatExchange { .. } => matches!(kind, K::Heater | K::Cooler),
            UnitSpec::Valve { .. } => kind == K::Valve,
            UnitSpec::Separator { .. } => kind == K::Separator,
            UnitSpec::Splitter { .. } => kind == K::Splitter,
            UnitSpec::HeatExchanger { .. } => kind == K::HeatExchanger,
            UnitSpec::Pressurize { .. } => matches!(kind, K::Pump | K::Compressor),
            UnitSpec::Expand { .. } => kind == K::Expander,
            UnitSpec::Rotating(s) => s.kind() == kind,
            UnitSpec::HeatTransport(s) => s.kind() == kind,
            UnitSpec::SimpleLogical(s) => s.kind() == kind,
            UnitSpec::Reactors(s) => s.kind() == kind,
            UnitSpec::Columns(s) => s.kind() == kind,
            UnitSpec::Multiphase(s) => s.kind() == kind,
            UnitSpec::Solids(s) => s.kind() == kind,
        }
    }

    /// Ports defined by the specification (numbered feeds, side draws,
    /// exchanger passes, …); they follow the kind's fixed ports.
    pub fn variable_ports(&self) -> Vec<PortSpec> {
        match self {
            UnitSpec::Rotating(s) => s.variable_ports(),
            UnitSpec::HeatTransport(s) => s.variable_ports(),
            UnitSpec::SimpleLogical(s) => s.variable_ports(),
            UnitSpec::Reactors(s) => s.variable_ports(),
            UnitSpec::Columns(s) => s.variable_ports(),
            UnitSpec::Multiphase(s) => s.variable_ports(),
            UnitSpec::Solids(s) => s.variable_ports(),
            _ => Vec::new(),
        }
    }

    /// Whether the unit determines or consumes the duty of the energy
    /// stream at `port` (`docs/a17-architecture.md`, energy streams).
    pub fn energy_role(&self, port: &PortSpec) -> EnergyRole {
        let consumes = match self {
            UnitSpec::EnergySource { duty } => *duty == EnergySourceDuty::Calculated,
            UnitSpec::EnergySink => true,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::EnergyStream,
                ..
            } => port.direction == PortDirection::Inlet,
            UnitSpec::Rotating(s) => s.consumes_energy(port),
            UnitSpec::HeatTransport(s) => s.consumes_energy(port),
            UnitSpec::SimpleLogical(s) => s.consumes_energy(port),
            UnitSpec::Reactors(s) => s.consumes_energy(port),
            UnitSpec::Columns(s) => s.consumes_energy(port),
            UnitSpec::Multiphase(s) => s.consumes_energy(port),
            UnitSpec::Solids(s) => s.consumes_energy(port),
            _ => false,
        };
        if consumes {
            EnergyRole::Consumes
        } else {
            EnergyRole::Determines
        }
    }
}

/// Appends an issue if `v` is not finite or outside `sign`. Group modules
/// use it for their own fields.
pub fn check_value(out: &mut Vec<SpecIssue>, field: &'static str, v: f64, sign: Sign) {
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

/// [`check_value`] for an optional value: `None` is reported as missing.
pub fn check_optional(out: &mut Vec<SpecIssue>, field: &'static str, v: Option<f64>, sign: Sign) {
    match v {
        None => out.push(SpecIssue {
            field,
            problem: SpecProblem::Missing,
        }),
        Some(v) => check_value(out, field, v, sign),
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
        if let UnitSpec::HeatExchange { thermal, .. } = self.spec {
            let problem = match (self.kind, thermal) {
                (UnitKind::Heater, ThermalSpec::Duty(Some(q))) if q < 0.0 => {
                    Some(("duty", SpecProblem::ContradictsUnitKind))
                }
                (UnitKind::Cooler, ThermalSpec::Duty(Some(q))) if q > 0.0 => {
                    Some(("duty", SpecProblem::ContradictsUnitKind))
                }
                // A cooler's energy port is an outlet: its duty is a result
                // the energy stream carries away, never an input (#171).
                (UnitKind::Cooler, ThermalSpec::EnergyStream) => {
                    Some(("duty from energy stream", SpecProblem::NotAvailableForKind))
                }
                _ => None,
            };
            if let Some((field, problem)) = problem {
                out.push(SpecIssue { field, problem });
            }
        }
        out
    }

    /// Every port of this unit in port-index order: the kind's fixed
    /// material ports, its fixed energy ports, then the variable ports of
    /// the specification.
    pub fn ports(&self) -> Vec<PortSpec> {
        ports_for(self.kind, &self.spec)
    }
}

/// Every port of a unit of `kind` with specification `spec`, in port-index
/// order (see [`Unit::ports`]).
pub fn ports_for(kind: UnitKind, spec: &UnitSpec) -> Vec<PortSpec> {
    let mut ports = kind.ports().to_vec();
    ports.extend_from_slice(kind.energy_ports());
    ports.extend(spec.variable_ports());
    ports
}

/// Reference to one port of one unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PortRef {
    /// Owning unit.
    pub unit: UnitId,
    /// Index into [`Unit::ports`].
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
            assert!(ports.iter().all(|p| p.medium == PortMedium::Material));
            assert!(
                kind.energy_ports()
                    .iter()
                    .all(|p| p.medium == PortMedium::Energy)
            );
            let inlets = ports
                .iter()
                .filter(|p| p.direction == PortDirection::Inlet)
                .count();
            let outlets = ports.len() - inlets;
            match kind {
                UnitKind::Feed => assert_eq!((inlets, outlets), (0, 1)),
                UnitKind::Product => assert_eq!((inlets, outlets), (1, 0)),
                UnitKind::EnergySource | UnitKind::EnergySink => {
                    assert_eq!((inlets, outlets), (0, 0))
                }
                UnitKind::Mixer => assert_eq!((inlets, outlets), (2, 1)),
                UnitKind::Separator | UnitKind::Splitter => {
                    assert_eq!((inlets, outlets), (1, 2))
                }
                UnitKind::HeatExchanger => assert_eq!((inlets, outlets), (2, 2)),
                UnitKind::Heater
                | UnitKind::Cooler
                | UnitKind::Valve
                | UnitKind::Pump
                | UnitKind::Compressor
                | UnitKind::Expander => assert_eq!((inlets, outlets), (1, 1)),
                _ => {}
            }
        }
    }

    #[test]
    fn port_labels_are_unique_within_every_default_unit() {
        for kind in UnitKind::ALL {
            let ports = kind.default_ports();
            for (i, a) in ports.iter().enumerate() {
                for b in &ports[i + 1..] {
                    assert!(!a.same_port(b), "{kind}: duplicate port {a}");
                    assert_ne!(a.label(), b.label(), "{kind}: duplicate label");
                }
            }
        }
    }

    #[test]
    fn every_kind_has_metadata_a_fitting_default_and_a_category() {
        let mut issues = Vec::new();
        for kind in UnitKind::ALL {
            let info = kind.info();
            assert!(!info.label.is_empty());
            assert!(!info.tag_series.is_empty());
            assert!(kind.default_spec().fits(kind), "{kind}");
            assert!(kind.category().kinds().any(|k| k == kind));
            let existing = matches!(
                kind,
                UnitKind::Pump | UnitKind::Compressor | UnitKind::Expander
            );
            if kind.group() != UnitGroup::Core && !existing {
                issues.push(kind.issue().expect("A17 kinds name their issue"));
            }
        }
        // #171 (energy streams) is core; one kind per remaining sub-issue.
        issues.sort();
        assert_eq!(issues, (172..=218).collect::<Vec<u16>>());
        // ALL is grouped by category in palette order.
        let order: Vec<_> = UnitKind::ALL.iter().map(|k| k.category()).collect();
        let mut sorted = order.clone();
        sorted.sort();
        assert_eq!(order, sorted);
        let mut unique = UnitKind::ALL.to_vec();
        unique.sort();
        unique.dedup();
        assert_eq!(unique.len(), UnitKind::ALL.len());
    }

    #[test]
    fn specifications_fit_only_their_kind() {
        for a in UnitKind::ALL {
            for b in UnitKind::ALL {
                let shared = matches!(
                    (a, b),
                    (UnitKind::Heater, UnitKind::Cooler)
                        | (UnitKind::Cooler, UnitKind::Heater)
                        | (UnitKind::Pump, UnitKind::Compressor)
                        | (UnitKind::Compressor, UnitKind::Pump)
                );
                assert_eq!(a.default_spec().fits(b), a == b || shared, "{a} / {b}");
            }
        }
    }

    #[test]
    fn tags_follow_pfd_conventions() {
        assert_eq!(UnitKind::Mixer.tag(1), "M-100");
        assert_eq!(UnitKind::Cooler.tag(2), "E-101");
        assert_eq!(UnitKind::Separator.tag(1), "D-100");
        assert_eq!(UnitKind::Feed.tag(2), "Feed 2");
        assert_eq!(UnitKind::EnergySource.tag(1), "Q-in 1");
        assert_eq!(UnitKind::EnergySource.tag_number("Q-in 3"), Some(3));
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
    fn core_kinds_are_calculated() {
        // Group kinds report their own status (`units::<group>::status`)
        // and are checked in their module's tests.
        for kind in UnitKind::ALL
            .into_iter()
            .filter(|k| k.group() == UnitGroup::Core)
        {
            assert!(kind.calculation_planned(), "{kind}");
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
        assert!(UnitKind::EnergySource.default_spec().issues().is_empty());
        assert_eq!(
            fields(UnitSpec::EnergySource {
                duty: EnergySourceDuty::Specified(None)
            }),
            [("duty", Missing)]
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
        let unit = |kind, thermal| Unit {
            id: UnitId(1),
            kind,
            name: "E-100".into(),
            spec: UnitSpec::HeatExchange {
                thermal,
                pressure_drop_pa: 0.0,
            },
        };
        let contradiction = [SpecIssue {
            field: "duty",
            problem: SpecProblem::ContradictsUnitKind,
        }];
        let duty = |q| ThermalSpec::Duty(Some(q));
        assert_eq!(unit(UnitKind::Heater, duty(-1.0)).issues(), contradiction);
        assert_eq!(unit(UnitKind::Cooler, duty(1.0)).issues(), contradiction);
        for (kind, q) in [
            (UnitKind::Heater, 1.0),
            (UnitKind::Heater, 0.0),
            (UnitKind::Cooler, -1.0),
            (UnitKind::Cooler, 0.0),
        ] {
            assert!(unit(kind, duty(q)).issues().is_empty(), "{kind} {q}");
        }
        assert!(
            unit(UnitKind::Heater, ThermalSpec::EnergyStream)
                .issues()
                .is_empty()
        );
        assert_eq!(
            unit(UnitKind::Cooler, ThermalSpec::EnergyStream).issues()[0].problem,
            SpecProblem::NotAvailableForKind
        );
    }

    #[test]
    fn energy_roles_follow_the_specification() {
        let q_in = PortSpec::energy_inlet("Q");
        let q_out = PortSpec::energy_outlet("Q");
        let heater = |thermal| UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa: 0.0,
        };
        assert_eq!(
            heater(ThermalSpec::EnergyStream).energy_role(&q_in),
            EnergyRole::Consumes
        );
        assert_eq!(
            heater(ThermalSpec::Duty(Some(1.0))).energy_role(&q_in),
            EnergyRole::Determines
        );
        assert_eq!(
            UnitKind::EnergySource.default_spec().energy_role(&q_out),
            EnergyRole::Consumes
        );
        assert_eq!(
            UnitSpec::EnergySource {
                duty: EnergySourceDuty::Specified(Some(1.0))
            }
            .energy_role(&q_out),
            EnergyRole::Determines
        );
        assert_eq!(
            UnitSpec::EnergySink.energy_role(&q_in),
            EnergyRole::Consumes
        );
    }
}
