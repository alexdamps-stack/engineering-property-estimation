//! Specification forms of the rotating group (pump, compressor, expander,
//! multistage compressor); owned by the rotating work package.

use eframe::egui::Ui;
use flowsheet_core::units::{GroupSpec, RotatingSpec};
use flowsheet_core::{PressureRise, PressureSpec, UnitKind, UnitSpec};

use super::{choice, drawn_only_note, note_row, opt_field, value_field};
use crate::units::{DisplayUnits, Quantity};

/// Pumps, compressors, and expanders need a PS flash (proposal section 2.3).
const ROTATING_DEFERRED: &str = "Calculation deferred: needs a PS flash and efficiency model.";

/// Form of the core pump/compressor ([`UnitSpec::Pressurize`]) and
/// expander ([`UnitSpec::Expand`]) specifications.
pub(super) fn core_form(ui: &mut Ui, spec: &mut UnitSpec, u: DisplayUnits) {
    match spec {
        UnitSpec::Pressurize {
            pressure,
            efficiency,
        } => {
            let current = matches!(pressure, PressureRise::PressureIncrease(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet p", "Pressure increase"]) {
                Some(0) => *pressure = PressureRise::OutletPressure(None),
                Some(_) => *pressure = PressureRise::PressureIncrease(None),
                None => {}
            }
            match pressure {
                PressureRise::OutletPressure(p) => {
                    opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u)
                }
                PressureRise::PressureIncrease(dp) => {
                    opt_field(ui, "Pressure increase", dp, Quantity::PressureDrop, u)
                }
            }
            value_field(ui, "Efficiency", efficiency, Quantity::Fraction, u);
            note_row(ui, ROTATING_DEFERRED);
        }
        UnitSpec::Expand {
            pressure,
            efficiency,
        } => {
            let current = matches!(pressure, PressureSpec::PressureDrop(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet p", "Pressure drop"]) {
                Some(0) => *pressure = PressureSpec::OutletPressure(None),
                Some(_) => *pressure = PressureSpec::PressureDrop(None),
                None => {}
            }
            match pressure {
                PressureSpec::OutletPressure(p) => {
                    opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u)
                }
                PressureSpec::PressureDrop(dp) => {
                    opt_field(ui, "Pressure drop", dp, Quantity::PressureDrop, u)
                }
            }
            value_field(
                ui,
                "Isentropic efficiency",
                efficiency,
                Quantity::Fraction,
                u,
            );
            note_row(ui, ROTATING_DEFERRED);
        }
        _ => unreachable!("not a pump, compressor or expander specification"),
    }
}

/// Form of a [`RotatingSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut RotatingSpec, _u: DisplayUnits) {
    let kind: UnitKind = spec.kind();
    match spec {
        RotatingSpec::MultistageCompressor(_) => {}
    }
    drawn_only_note(ui, kind);
}
