//! Native flowsheet editor. See the library crate docs.
//!
//! Usage:
//!
//! - `thermo-desktop` opens an empty editor.
//! - `thermo-desktop --demo` starts with the demonstration project.
//! - `thermo-desktop <file.flowsheet.json>` opens a project file.
//! - `thermo-desktop --self-test <report.json>` runs the headless A9
//!   acceptance workflow (solve, invalidate, save, reload, recalculate,
//!   failure and cancellation) without opening a window, writes a JSON
//!   report and exits with status 0 on success and 1 on failure. Release
//!   builds have no console on Windows, so the report file is the output.
//! - `thermo-desktop --version` writes the version to the report path given
//!   after it, or to standard output.

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::path::PathBuf;
use std::process::ExitCode;

use thermo_desktop::FlowsheetApp;

fn self_test(report: PathBuf) -> ExitCode {
    let dir = std::env::temp_dir().join(format!("thermo-desktop-self-test-{}", std::process::id()));
    let outcome = std::fs::create_dir_all(&dir)
        .map_err(|e| e.to_string())
        .and_then(|()| flowsheet_project::run_acceptance(&dir));
    let _ = std::fs::remove_dir_all(&dir);
    let (json, code) = match outcome {
        Ok(r) => (
            serde_json::json!({ "status": "passed", "version": env!("CARGO_PKG_VERSION"), "report": r }),
            ExitCode::SUCCESS,
        ),
        Err(e) => (
            serde_json::json!({ "status": "failed", "version": env!("CARGO_PKG_VERSION"), "error": e }),
            ExitCode::FAILURE,
        ),
    };
    let text = serde_json::to_string_pretty(&json).unwrap_or_default();
    if std::fs::write(&report, format!("{text}\n")).is_err() {
        println!("{text}");
    }
    code
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("--self-test") => {
            let report = args.get(1).map_or_else(
                || PathBuf::from("thermo-desktop-self-test.json"),
                PathBuf::from,
            );
            return self_test(report);
        }
        Some("--version") => {
            let text = format!("thermo-desktop {}\n", env!("CARGO_PKG_VERSION"));
            match args.get(1) {
                Some(path) => {
                    let _ = std::fs::write(path, text);
                }
                None => print!("{text}"),
            }
            return ExitCode::SUCCESS;
        }
        _ => {}
    }
    let demo = args.iter().any(|a| a == "--demo");
    let file = args
        .iter()
        .find(|a| !a.starts_with("--"))
        .map(PathBuf::from);
    let options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_title("Property Estimation Flowsheet")
            .with_inner_size([1280.0, 800.0])
            .with_min_inner_size([800.0, 500.0]),
        ..Default::default()
    };
    let result = eframe::run_native(
        "thermo-desktop",
        options,
        Box::new(move |_cc| {
            Ok(Box::new(match (file, demo) {
                (Some(path), _) => FlowsheetApp::with_file(path),
                (None, true) => FlowsheetApp::with_demonstration(),
                (None, false) => FlowsheetApp::default(),
            }))
        }),
    );
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("thermo-desktop: {e}");
            ExitCode::FAILURE
        }
    }
}
