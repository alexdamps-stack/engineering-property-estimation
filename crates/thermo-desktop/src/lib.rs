//! Native flowsheet editor (egui).
//!
//! The editor follows `docs/flowsheet-and-ui.md`. It edits the
//! GUI-independent topology model from `flowsheet-core`, solves it in the
//! background through the A8 engine with the immutable PR/PR-BM package, and
//! reads and writes versioned project files through `flowsheet-project`.
//! Results carry input fingerprints, so every edit marks exactly the
//! affected units and streams stale until the next solve, and failed solves
//! never leave old values marked current.
//!
//! The application is a library so that headless UI tests in `tests/` can
//! drive it with `egui_kittest` without opening a window.

pub mod app;
pub mod canvas;
pub mod editor;
pub mod panels;
mod solve;
pub mod symbols;
pub mod units;
pub mod view;

pub use app::FlowsheetApp;
