//! Streams, unit operations, balances, and calculation order.
//!
//! This crate provides the **flowsheet topology and sequential engine**:
//! unit operations, their material and energy ports, the streams connecting
//! them, user specifications in SI units with absolute pressure, and
//! pre-calculation diagnostics, balances and flash-based unit operations
//! (A8). See the specification sets in `docs/flowsheet-and-ui.md`, section
//! 2.3. A10 adds [`Flowsheet::solve_network`] with deterministic automatic
//! tears, or [`Flowsheet::solve_recycles`] with explicit starting states and
//! limits. A17 adds energy streams (#171), the full unit palette and the
//! per-group modules in [`units`] (`docs/a17-architecture.md`).
//!
//! The model is independent of any GUI framework. Graphical layout (block
//! positions, zoom, selection) belongs to the consuming application.
//!
//! Every mutation increments [`Flowsheet::revision`]. [`Calculation::is_current`]
//! compares the solve-relevant topology/specification snapshot (not names or
//! the revision counter), the feed compositions and immutable package identity,
//! so results from an older or independently edited project are not current.

mod engine;
mod flowsheet;
mod numerics;
mod reaction;
mod reactors;
mod recycle;
mod unit;
pub mod units;

pub use engine::*;
pub use reaction::{ReactionBalance, ReactionTerms};
pub use recycle::*;

pub use flowsheet::{
    ConnectError, Diagnostic, DiagnosticKind, Flowsheet, MAX_IDENTIFIER, RestoreError,
    SpecChangeError, Stream, StreamId, UnitIdsExhausted, UnitStatus,
};
pub use unit::{
    Arrhenius, BatchThermal, CalculationStatus, ConcentrationBasis, Correlation,
    DEFAULT_EFFICIENCY, EnergyRole, EnergySourceDuty, EquilibriumConstant, EquilibriumOptions,
    EquilibriumReaction, ExchangerSpec, FugacityBasis, KindInfo, KineticReaction, Kinetics,
    LhhwTerm, MixerPressure, PaletteCategory, PlugFlowThermal, PortDirection, PortMedium, PortRef,
    PortSpec, PressureRise, PressureSpec, RateBasis, RateLaw, RateTerm, ReactingPhase, Reaction,
    ReactionExtent, SeparatorThermal, Sign, SpecIssue, SpecProblem, ThermalSpec, Unit, UnitGroup,
    UnitId, UnitKind, UnitSpec, check_optional, check_value, ports_for,
};
