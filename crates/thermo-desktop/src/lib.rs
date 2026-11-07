//! Native flowsheet editor (egui GUI prototype).
//!
//! This prototype evaluates the canvas, connections, and properties form
//! called for by development plan section 8 and laid out in
//! `docs/flowsheet-and-ui.md`. It edits the GUI-independent topology model
//! from `flowsheet-core`; it performs no thermodynamic calculation and
//! displays no results.
//!
//! The application is a library so that headless UI tests in `tests/` can
//! drive it with `egui_kittest` without opening a window.

pub mod app;
pub mod canvas;
pub mod editor;
pub mod symbols;
pub mod units;
pub mod view;

pub use app::FlowsheetApp;
