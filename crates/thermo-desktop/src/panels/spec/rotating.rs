//! Specification forms of the rotating group (pump, compressor, expander,
//! multistage compressor); owned by the rotating work package.

use eframe::egui::Ui;
use flowsheet_core::units::RotatingSpec;
use flowsheet_core::units::rotating::MAX_COMPRESSOR_STAGES;
use flowsheet_core::{PressureRise, PressureSpec, UnitSpec};

use super::{choice, count_field, note_row, opt_field, value_field};
use crate::units::{DisplayUnits, Quantity};

const COMPRESSION_NOTE: &str = "Isentropic efficiency: Δh = Δh_is / η (PS then PH flash).";
const EXPANSION_NOTE: &str = "Isentropic efficiency: Δh = η · Δh_is (PS then PH flash).";

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
            value_field(
                ui,
                "Isentropic efficiency",
                efficiency,
                Quantity::Fraction,
                u,
            );
            note_row(ui, COMPRESSION_NOTE);
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
            note_row(ui, EXPANSION_NOTE);
        }
        _ => unreachable!("not a pump, compressor or expander specification"),
    }
}

/// Form of a [`RotatingSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut RotatingSpec, u: DisplayUnits) {
    match spec {
        RotatingSpec::MultistageCompressor(s) => {
            count_field(ui, "Stages", &mut s.stages, 1..=MAX_COMPRESSOR_STAGES);
            opt_field(
                ui,
                "Outlet pressure",
                &mut s.outlet_pressure_pa,
                Quantity::Pressure,
                u,
            );
            value_field(
                ui,
                "Stage isentropic efficiency",
                &mut s.efficiency,
                Quantity::Fraction,
                u,
            );
            if s.stages > 1 {
                opt_field(
                    ui,
                    "Intercooler outlet T",
                    &mut s.intercooler_temperature_k,
                    Quantity::Temperature,
                    u,
                );
            }
            note_row(
                ui,
                "Equal stage pressure ratio; intercoolers between stages, no aftercooler.",
            );
        }
    }
}
