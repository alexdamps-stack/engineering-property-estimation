//! Unit specification forms.
//!
//! Core kinds are edited here; each A17 group owns its form in
//! `spec/<group>.rs` (`docs/a17-architecture.md`). Group forms use the
//! field helpers of this module (`choice`, `opt_field`, `value_field`,
//! `count_field`, `note_row`, `drawn_only_note`).

use eframe::egui::{DragValue, RichText, Ui};
use flowsheet_core::{
    CalculationStatus, EnergySourceDuty, ExchangerSpec, MixerPressure, PressureSpec,
    SeparatorThermal, ThermalSpec, UnitKind, UnitSpec,
};

use crate::symbols::Theme;
use crate::units::{DisplayUnits, DisplayUnitsExt, Quantity, starting_value};

mod columns;
mod heat_transport;
mod multiphase;
mod reactors;
mod rotating;
mod simple_logical;
mod solids;

/// Specification form; mutually exclusive choices use radio buttons
/// (proposal section 3).
pub(super) fn spec_form(
    ui: &mut Ui,
    kind: UnitKind,
    spec: &mut UnitSpec,
    u: DisplayUnits,
    components: &[String],
) {
    match spec {
        UnitSpec::Product | UnitSpec::EnergySink => {
            ui.label(RichText::new("No specification.").weak());
            ui.end_row();
        }
        UnitSpec::EnergySource { duty } => {
            let current = matches!(duty, EnergySourceDuty::Specified(_)) as usize;
            match choice(ui, "Duty", current, &["Calculated", "Specified"]) {
                Some(0) => *duty = EnergySourceDuty::Calculated,
                Some(_) => *duty = EnergySourceDuty::Specified(None),
                None => {}
            }
            match duty {
                EnergySourceDuty::Calculated => note_row(
                    ui,
                    "Reports the duty the connected unit calculates from its own specification.",
                ),
                EnergySourceDuty::Specified(q) => {
                    opt_field(ui, "Duty (+ to the unit)", q, Quantity::Duty, u);
                    note_row(ui, "The connected unit must take its duty from the stream.");
                }
            }
        }
        UnitSpec::Rotating(s) => rotating::form(ui, s, u),
        UnitSpec::HeatTransport(s) => heat_transport::form(ui, s, u),
        UnitSpec::SimpleLogical(s) => simple_logical::form(ui, s, u, components),
        UnitSpec::Reactors(s) => reactors::form(ui, s, u),
        UnitSpec::Columns(s) => columns::form(ui, s, u),
        UnitSpec::Multiphase(s) => multiphase::form(ui, s, u),
        UnitSpec::Solids(s) => solids::form(ui, s, u),
        UnitSpec::Feed {
            temperature_k,
            pressure_pa,
            molar_flow_mol_s,
        } => {
            opt_field(ui, "Temperature", temperature_k, Quantity::Temperature, u);
            opt_field(ui, "Pressure", pressure_pa, Quantity::Pressure, u);
            opt_field(ui, "Molar flow", molar_flow_mol_s, Quantity::MolarFlow, u);
        }
        UnitSpec::Mixer { outlet_pressure } => {
            let current = matches!(outlet_pressure, MixerPressure::Specified(_)) as usize;
            match choice(ui, "Outlet p", current, &["Lowest inlet", "Specified"]) {
                Some(0) => *outlet_pressure = MixerPressure::LowestInlet,
                Some(_) => *outlet_pressure = MixerPressure::Specified(None),
                None => {}
            }
            if let MixerPressure::Specified(p) = outlet_pressure {
                opt_field(ui, "Outlet pressure", p, Quantity::Pressure, u);
            }
            note_row(ui, "Adiabatic.");
        }
        UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa,
        } => {
            let current = match thermal {
                ThermalSpec::OutletTemperature(_) => 0,
                ThermalSpec::Duty(_) => 1,
                ThermalSpec::EnergyStream => 2,
            };
            // Only a heater can take its duty from an energy stream: a
            // cooler's energy port is an outlet (#171).
            let options: &[&str] = if kind == UnitKind::Heater {
                &["Outlet T", "Duty", "Energy stream"]
            } else {
                &["Outlet T", "Duty"]
            };
            match choice(ui, "Specify", current, options) {
                Some(0) => *thermal = ThermalSpec::OutletTemperature(None),
                Some(1) => *thermal = ThermalSpec::Duty(None),
                Some(_) => *thermal = ThermalSpec::EnergyStream,
                None => {}
            }
            match thermal {
                ThermalSpec::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
                }
                ThermalSpec::Duty(q) => opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u),
                ThermalSpec::EnergyStream => note_row(ui, "Duty from the energy stream at port Q."),
            }
            value_field(
                ui,
                "Pressure drop",
                pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
        UnitSpec::Valve { pressure } => {
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
            note_row(ui, "Isenthalpic.");
        }
        UnitSpec::Separator {
            thermal,
            pressure_drop_pa,
        } => vessel_thermal_fields(ui, thermal, pressure_drop_pa, u),
        UnitSpec::Splitter { fraction_to_first } => {
            opt_field(
                ui,
                "Fraction to out 1",
                fraction_to_first,
                Quantity::Fraction,
                u,
            );
            note_row(ui, "Out 2 receives the rest; same state in both.");
        }
        UnitSpec::HeatExchanger {
            thermal,
            hot_pressure_drop_pa,
            cold_pressure_drop_pa,
        } => {
            let current = match thermal {
                ExchangerSpec::HotOutletTemperature(_) => 0,
                ExchangerSpec::ColdOutletTemperature(_) => 1,
                ExchangerSpec::Duty(_) => 2,
            };
            match choice(ui, "Specify", current, &["Hot out T", "Cold out T", "Duty"]) {
                Some(0) => *thermal = ExchangerSpec::HotOutletTemperature(None),
                Some(1) => *thermal = ExchangerSpec::ColdOutletTemperature(None),
                Some(_) => *thermal = ExchangerSpec::Duty(None),
                None => {}
            }
            match thermal {
                ExchangerSpec::HotOutletTemperature(t) => {
                    opt_field(ui, "Hot outlet temperature", t, Quantity::Temperature, u)
                }
                ExchangerSpec::ColdOutletTemperature(t) => {
                    opt_field(ui, "Cold outlet temperature", t, Quantity::Temperature, u)
                }
                ExchangerSpec::Duty(q) => opt_field(ui, "Duty (hot to cold)", q, Quantity::Duty, u),
            }
            value_field(
                ui,
                "Hot-side pressure drop",
                hot_pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
            value_field(
                ui,
                "Cold-side pressure drop",
                cold_pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
        UnitSpec::Pressurize { .. } | UnitSpec::Expand { .. } => rotating::core_form(ui, spec, u),
    }
}

/// Specification controls shared by the separator and steady-state tank.
fn vessel_thermal_fields(
    ui: &mut Ui,
    thermal: &mut SeparatorThermal,
    pressure_drop_pa: &mut f64,
    u: DisplayUnits,
) {
    let current = match thermal {
        SeparatorThermal::Adiabatic => 0,
        SeparatorThermal::OutletTemperature(_) => 1,
        SeparatorThermal::Duty(_) => 2,
    };
    match choice(ui, "Thermal", current, &["Adiabatic", "Outlet T", "Duty"]) {
        Some(0) => *thermal = SeparatorThermal::Adiabatic,
        Some(1) => *thermal = SeparatorThermal::OutletTemperature(None),
        Some(_) => *thermal = SeparatorThermal::Duty(None),
        None => {}
    }
    match thermal {
        SeparatorThermal::Adiabatic => {}
        SeparatorThermal::OutletTemperature(t) => {
            opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
        }
        SeparatorThermal::Duty(q) => opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u),
    }
    value_field(
        ui,
        "Pressure drop",
        pressure_drop_pa,
        Quantity::PressureDrop,
        u,
    );
}

/// A row of radio buttons; returns the newly chosen index if it changed.
fn choice(ui: &mut Ui, label: &str, current: usize, options: &[&str]) -> Option<usize> {
    ui.label(label);
    let mut picked = current;
    ui.horizontal_wrapped(|ui| {
        for (i, text) in options.iter().enumerate() {
            ui.radio_value(&mut picked, i, *text);
        }
    });
    ui.end_row();
    (picked != current).then_some(picked)
}

fn note_row(ui: &mut Ui, text: &str) {
    ui.label("");
    ui.label(RichText::new(text).small().weak());
    ui.end_row();
}

fn drag_value<'a>(shown: &'a mut f64, q: Quantity, u: DisplayUnits) -> DragValue<'a> {
    DragValue::new(shown)
        .speed(u.drag_speed(q))
        .range(u.display_min(q)..=u.display_max(q))
        // Limits apply to user edits only: a stored out-of-range value (for
        // example from a file) is shown as is and never rewritten on display.
        .clamp_existing_to_range(false)
        .max_decimals(6)
        .suffix(format!(" {}", u.symbol(q)))
}

/// An optional SI value edited in display units. Unset values are
/// highlighted and offer an "Enter value" button.
fn opt_field(ui: &mut Ui, label: &str, value: &mut Option<f64>, q: Quantity, u: DisplayUnits) {
    ui.label(label);
    ui.horizontal(|ui| match value {
        Some(si) => {
            let mut shown = u.si_to_display(q, *si);
            if ui.add(drag_value(&mut shown, q, u)).changed() {
                *si = u.display_to_si(q, shown);
            }
            if ui.small_button("×").on_hover_text("Clear").clicked() {
                *value = None;
            }
        }
        None => {
            let warning = Theme::from_visuals(ui.visuals()).warning;
            ui.label(RichText::new("required").italics().color(warning));
            if ui
                .small_button("Enter value")
                .on_hover_text("Enter a value (starts from standard conditions)")
                .clicked()
            {
                *value = Some(starting_value(q));
            }
        }
    });
    ui.end_row();
}

/// A value that always has a default (for example a pressure drop of zero).
fn value_field(ui: &mut Ui, label: &str, si: &mut f64, q: Quantity, u: DisplayUnits) {
    ui.label(label);
    let mut shown = u.si_to_display(q, *si);
    if ui.add(drag_value(&mut shown, q, u)).changed() {
        *si = u.display_to_si(q, shown);
    }
    ui.end_row();
}

/// A count (number of feeds, passes, …) edited as an integer in `range`.
/// Returns whether it changed.
fn count_field(
    ui: &mut Ui,
    label: &str,
    n: &mut u16,
    range: std::ops::RangeInclusive<u16>,
) -> bool {
    let label = ui.label(label);
    let changed = ui
        .add(DragValue::new(n).range(range).speed(0.05))
        .labelled_by(label.id)
        .changed();
    ui.end_row();
    changed
}

/// Note shown in the form of a kind that is drawn only.
fn drawn_only_note(ui: &mut Ui, kind: UnitKind) {
    if let CalculationStatus::DrawnOnly(reason) = kind.calculation_status() {
        let issue = kind
            .issue()
            .map_or_else(String::new, |n| format!(" (#{n})"));
        let warning = Theme::from_visuals(ui.visuals()).warning;
        ui.label("");
        ui.label(
            RichText::new(format!("Not calculated: {kind} {reason}{issue}."))
                .small()
                .color(warning),
        );
        ui.end_row();
    }
}
