//! Unit specification forms.

use eframe::egui::{DragValue, RichText, Ui};
use flowsheet_core::{
    ExchangerSpec, MixerPressure, PressureRise, PressureSpec, SeparatorThermal, ThermalSpec,
    UnitSpec,
};

use crate::symbols::Theme;
use crate::units::{DisplayUnits, DisplayUnitsExt, Quantity, starting_value};

/// Specification form; mutually exclusive choices use radio buttons
/// (proposal section 3).
pub(super) fn spec_form(ui: &mut Ui, spec: &mut UnitSpec, u: DisplayUnits) {
    match spec {
        UnitSpec::Product => {
            ui.label(RichText::new("No specification.").weak());
            ui.end_row();
        }
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
            let current = matches!(thermal, ThermalSpec::Duty(_)) as usize;
            match choice(ui, "Specify", current, &["Outlet T", "Duty"]) {
                Some(0) => *thermal = ThermalSpec::OutletTemperature(None),
                Some(_) => *thermal = ThermalSpec::Duty(None),
                None => {}
            }
            match thermal {
                ThermalSpec::OutletTemperature(t) => {
                    opt_field(ui, "Outlet temperature", t, Quantity::Temperature, u)
                }
                ThermalSpec::Duty(q) => opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u),
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
        } => {
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
                SeparatorThermal::Duty(q) => {
                    opt_field(ui, "Duty (+ into unit)", q, Quantity::Duty, u)
                }
            }
            value_field(
                ui,
                "Pressure drop",
                pressure_drop_pa,
                Quantity::PressureDrop,
                u,
            );
        }
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
    }
}

/// Pumps, compressors, and expanders need a PS flash (proposal section 2.3).
const ROTATING_DEFERRED: &str = "Calculation deferred: needs a PS flash and efficiency model.";

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
