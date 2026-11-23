//! Specification forms of the heat-transport group; owned by the
//! heat-transport work package.

use eframe::egui::Ui;
use flowsheet_core::units::heat_transport::{AirCoolerThermal, MAX_PASSES};
use flowsheet_core::units::{GroupSpec, HeatTransportSpec};

use super::{
    choice, count_field, drawn_only_note, note_row, opt_field, process_thermal_fields, value_field,
};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`HeatTransportSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut HeatTransportSpec, u: DisplayUnits) {
    let kind = spec.kind();
    match spec {
        HeatTransportSpec::AirCooler(s) => {
            let current = matches!(s.thermal, AirCoolerThermal::Duty(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet T", "Duty"]) {
                Some(0) => s.thermal = AirCoolerThermal::OutletTemperature(None),
                Some(_) => s.thermal = AirCoolerThermal::Duty(None),
                None => {}
            }
            match &mut s.thermal {
                AirCoolerThermal::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u);
                }
                AirCoolerThermal::Duty(q) => {
                    opt_field(ui, "Process duty (≤ 0)", q, Quantity::Duty, u);
                }
            }
            value_field(
                ui,
                "Pressure drop",
                &mut s.pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
            opt_field(
                ui,
                "Ambient air temperature",
                &mut s.ambient_temperature_k,
                Quantity::Temperature,
                u,
            );
            opt_field(
                ui,
                "Dry-air mass flow",
                &mut s.air_mass_flow_kg_s,
                Quantity::MassFlow,
                u,
            );
            value_field(ui, "Fan power (≥ 0)", &mut s.fan_power_w, Quantity::Duty, u);
            note_row(
                ui,
                "Fan power warms the discharge air. A calculated energy source can report it.",
            );
        }
        HeatTransportSpec::FiredHeater(s) => {
            process_thermal_fields(
                ui,
                &mut s.thermal,
                &mut s.pressure_drop_pa,
                u,
                Some("The energy stream at the duty port supplies the heat."),
            );
            note_row(
                ui,
                "Process-side heat duty only. Combustion and fuel use are outside this model.",
            );
        }
        HeatTransportSpec::MultiStreamExchanger(s) => {
            count_field(ui, "Passes", &mut s.passes, 2..=MAX_PASSES);
            note_row(ui, "Each pass k has the ports in k and out k.");
        }
        HeatTransportSpec::PipeSegment(_) | HeatTransportSpec::ReliefValve(_) => {}
    }
    drawn_only_note(ui, kind);
}
