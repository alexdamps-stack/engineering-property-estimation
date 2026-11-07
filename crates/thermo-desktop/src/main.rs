//! Native flowsheet editor (egui GUI prototype). See the library crate docs.
//!
//! Pass `--demo` to start with the demonstration flowsheet loaded.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use thermo_desktop::FlowsheetApp;

fn main() -> eframe::Result {
    let demo = std::env::args().any(|a| a == "--demo");
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_title("Property Estimation — Flowsheet")
            .with_inner_size([1280.0, 800.0])
            .with_min_inner_size([800.0, 500.0]),
        ..Default::default()
    };
    eframe::run_native(
        "thermo-desktop",
        options,
        Box::new(move |_cc| {
            Ok(Box::new(if demo {
                FlowsheetApp::with_demonstration()
            } else {
                FlowsheetApp::default()
            }))
        }),
    )
}
