//! Specification forms of the heat-transport group; owned by the
//! heat-transport work package.

use eframe::egui::{ComboBox, Ui};
use flowsheet_core::units::heat_transport::{
    AirCoolerThermal, ExchangerPassDirection, ExchangerPassSpec, ExchangerPassThermal, MAX_PASSES,
    MultiStreamExchangerSpec,
};
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
            let current = match s.thermal {
                AirCoolerThermal::OutletTemperature(_) => 0,
                AirCoolerThermal::Duty(_) => 1,
                AirCoolerThermal::HotEndApproach(_) => 2,
                AirCoolerThermal::ColdEndApproach(_) => 3,
            };
            match choice(
                ui,
                "Specify",
                current,
                &[
                    "Outlet T",
                    "Duty",
                    "Fixed hot-end approach",
                    "Fixed cold-end approach",
                ],
            ) {
                Some(0) => s.thermal = AirCoolerThermal::OutletTemperature(None),
                Some(1) => s.thermal = AirCoolerThermal::Duty(None),
                Some(2) => s.thermal = AirCoolerThermal::HotEndApproach(None),
                Some(_) => s.thermal = AirCoolerThermal::ColdEndApproach(None),
                None => {}
            }
            match &mut s.thermal {
                AirCoolerThermal::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u);
                }
                AirCoolerThermal::Duty(q) => {
                    opt_field(ui, "Process duty (≤ 0)", q, Quantity::Duty, u);
                }
                AirCoolerThermal::HotEndApproach(dt) => {
                    opt_field(
                        ui,
                        "Process inlet − pre-fan air outlet",
                        dt,
                        Quantity::TemperatureDifference,
                        u,
                    );
                    note_row(
                        ui,
                        "Sets the air outlet before fan heat and solves duty and process outlet.",
                    );
                }
                AirCoolerThermal::ColdEndApproach(dt) => {
                    opt_field(
                        ui,
                        "Process outlet − ambient air inlet",
                        dt,
                        Quantity::TemperatureDifference,
                        u,
                    );
                    note_row(
                        ui,
                        "Sets the process outlet and solves duty and air discharge.",
                    );
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
        HeatTransportSpec::MultiStreamExchanger(s) => multi_stream_form(ui, s, u),
        HeatTransportSpec::PipeSegment(_) | HeatTransportSpec::ReliefValve(_) => {}
    }
    drawn_only_note(ui, kind);
}

fn multi_stream_form(ui: &mut Ui, spec: &mut MultiStreamExchangerSpec, u: DisplayUnits) {
    use ExchangerPassDirection::{Cold, Hot};
    use ExchangerPassThermal::{Balance, Duty, OutletTemperature, TerminalApproach};

    let mut count = spec.passes.len() as u16;
    if count_field(ui, "Passes", &mut count, 2..=MAX_PASSES) {
        spec.passes.resize(
            usize::from(count),
            ExchangerPassSpec {
                direction: Cold,
                thermal: OutletTemperature(None),
                pressure_drop_pa: 0.0,
            },
        );
    }
    note_row(
        ui,
        "Countercurrent flow. Each inlet goes to the outlet with the same number.",
    );
    value_field(
        ui,
        "Minimum terminal approach",
        &mut spec.minimum_approach_k,
        Quantity::TemperatureDifference,
        u,
    );
    note_row(
        ui,
        "Assign one pass to heat balance. Approach checks cover the two ends of the exchanger.",
    );

    for i in 0..spec.passes.len() {
        let number = i + 1;
        let cold = spec.passes[i].direction == Cold;
        if let Some(pick) = choice(
            ui,
            &format!("Pass {number}"),
            usize::from(cold),
            &["Hot", "Cold"],
        ) {
            spec.passes[i].direction = if pick == 0 { Hot } else { Cold };
        }
        let references: Vec<u16> = spec
            .passes
            .iter()
            .enumerate()
            .filter(|(_, p)| p.direction != spec.passes[i].direction)
            .map(|(j, _)| (j + 1) as u16)
            .collect();
        let mode = match spec.passes[i].thermal {
            OutletTemperature(_) => 0,
            Duty(_) => 1,
            TerminalApproach { .. } => 2,
            Balance => 3,
        };
        if let Some(pick) = choice(
            ui,
            "Specify",
            mode,
            &["Outlet T", "Duty", "Fixed approach", "Heat balance"],
        ) {
            let thermal = match pick {
                0 => OutletTemperature(None),
                1 => Duty(None),
                2 => TerminalApproach {
                    reference_pass: references.first().copied().unwrap_or(0),
                    approach_k: None,
                },
                _ => {
                    for pass in &mut spec.passes {
                        if pass.thermal == Balance {
                            pass.thermal = OutletTemperature(None);
                        }
                    }
                    Balance
                }
            };
            spec.passes[i].thermal = thermal;
        }
        let pass = &mut spec.passes[i];
        match &mut pass.thermal {
            OutletTemperature(t) => {
                opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
            }
            Duty(q) => {
                let label = match pass.direction {
                    Hot => "Duty (≤ 0)",
                    Cold => "Duty (≥ 0)",
                };
                opt_field(ui, label, q, Quantity::Duty, u);
            }
            TerminalApproach {
                reference_pass,
                approach_k,
            } => {
                let label = ui.label("Reference inlet");
                let selected = if references.contains(reference_pass) {
                    format!("Pass {reference_pass} inlet")
                } else {
                    "Select opposite-direction inlet".into()
                };
                ComboBox::from_id_salt(("exchanger-reference-inlet", number))
                    .selected_text(selected)
                    .show_ui(ui, |ui| {
                        for reference in &references {
                            ui.selectable_value(
                                reference_pass,
                                *reference,
                                format!("Pass {reference} inlet"),
                            );
                        }
                    })
                    .response
                    .labelled_by(label.id);
                ui.end_row();
                let label = if pass.direction == Hot {
                    "Hot outlet − cold inlet"
                } else {
                    "Hot inlet − cold outlet"
                };
                opt_field(ui, label, approach_k, Quantity::TemperatureDifference, u);
            }
            Balance => note_row(
                ui,
                "Calculates duty and outlet temperature from the other passes.",
            ),
        }
        value_field(
            ui,
            "Pressure drop",
            &mut pass.pressure_drop_pa,
            Quantity::PressureDrop,
            u,
        );
    }
}
