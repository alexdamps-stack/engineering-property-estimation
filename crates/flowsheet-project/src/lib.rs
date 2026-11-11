//! Versioned project files, property-package selection, result snapshots and
//! result invalidation for the flowsheet editor (A9).
//!
//! This crate is independent of any GUI framework. It turns the A8 engine
//! (`flowsheet-core`) into a persistent, recalculable project:
//!
//! - [`PackageSelection`] records the exact cubic package choice (variant,
//!   ordered components, explicit interaction matrix and its source, PH
//!   search interval) and builds the immutable [`CubicPackage`].
//! - [`Project`] holds the flowsheet, feed mole fractions, layout and display
//!   settings, and reads/writes the versioned JSON project file
//!   (`docs/flowsheet-and-ui.md`, section 5). Unknown format versions fail.
//! - [`Fingerprints`] identify the inputs each unit result depends on: the
//!   unit's own kind and specification, its feed composition, the package and
//!   result provenance, and recursively every upstream unit. An edit
//!   therefore invalidates exactly the affected unit and everything
//!   downstream of it.
//! - [`ResultSet`] is a serializable snapshot of a successful [`Calculation`].
//!   Every entry carries the fingerprint it was computed for, so stored or
//!   loaded results are shown as current only while their inputs and the
//!   engine/model/data provenance are unchanged.
//!
//! [`CubicPackage`]: property_packages::CubicPackage
//! [`Calculation`]: flowsheet_core::Calculation

mod acceptance;
mod demonstration;
mod error;
mod fingerprint;
mod package;
mod project;
mod records;
mod results;
mod time;

pub use acceptance::{AcceptanceReport, run_acceptance};
pub use demonstration::{
    DEMONSTRATION_COMPONENTS, DEMONSTRATION_INTERACTIONS, demonstration, recycle_demonstration,
};
pub use error::ProjectError;
pub use fingerprint::{Fingerprint, Fingerprints};
pub use package::{CubicVariantChoice, PackageSelection, Provenance};
pub use project::{
    DisplayUnitSet, FORMAT_VERSION, FeedFractions, InputIssue, Project, SolveFailure, solve_project,
};
pub use results::{
    Currency, PhaseResult, RecycleResult, ResultSet, StateResult, StreamResult, UnitResult,
};
pub use time::format_utc;

/// Revision of the flowsheet engine numerics recorded in result provenance.
/// Change it whenever `flowsheet-core` unit-operation equations, balances or
/// routing conventions change, so that stored results load as stale.
pub const ENGINE_REVISION: &str = "flowsheet-core a10-tears-v2";
