//! Streams, unit operations, balances, and calculation order.
//!
//! This crate provides the **flowsheet topology and sequential engine**:
//! unit operations, their material ports, the streams connecting them, user
//! specifications in SI units with absolute pressure, and pre-calculation
//! diagnostics, balances and flash-based unit operations (A8). See the
//! specification sets in `docs/flowsheet-and-ui.md`, section 2.3.
//!
//! The model is independent of any GUI framework. Graphical layout (block
//! positions, zoom, selection) belongs to the consuming application.
//!
//! Every mutation increments [`Flowsheet::revision`]. [`Calculation::is_current`]
//! compares the full topology/input snapshot and immutable package identity,
//! so results from an older or independently edited project are not current.

mod engine;
mod flowsheet;
mod unit;

pub use engine::*;

pub use flowsheet::{
    ConnectError, Diagnostic, DiagnosticKind, Flowsheet, RestoreError, Stream, StreamId, UnitStatus,
};
pub use unit::{
    DEFAULT_EFFICIENCY, ExchangerSpec, MixerPressure, PortDirection, PortRef, PortSpec,
    PressureRise, PressureSpec, SeparatorThermal, SpecIssue, SpecProblem, ThermalSpec, Unit,
    UnitId, UnitKind, UnitSpec,
};
