//! Streams, unit operations, balances, and calculation order.
//!
//! This crate currently provides the **flowsheet topology model** only:
//! unit operations, their material ports, the streams connecting them, user
//! specifications in SI units with absolute pressure, and pre-calculation
//! diagnostics. It performs no thermodynamic calculation. Balances,
//! flash-based unit operations, and calculation order are introduced with
//! work package A8 (see `docs/development-plan.md` and the specification
//! sets in `docs/flowsheet-and-ui.md`, section 2.3).
//!
//! The model is independent of any GUI framework. Graphical layout (block
//! positions, zoom, selection) belongs to the consuming application.
//!
//! Every mutation increments [`Flowsheet::revision`]. Future calculation
//! results must record the revision they were computed from, so that results
//! derived from an older revision are never presented as current.

mod flowsheet;
mod unit;

pub use flowsheet::{
    ConnectError, Diagnostic, DiagnosticKind, Flowsheet, Stream, StreamId, UnitStatus,
};
pub use unit::{
    MixerPressure, PortDirection, PortRef, PortSpec, PressureSpec, SeparatorThermal, SpecIssue,
    SpecProblem, ThermalSpec, Unit, UnitId, UnitKind, UnitSpec,
};
