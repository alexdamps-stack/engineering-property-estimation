//! Specification forms of the columns group; owned by the columns work
//! package.

use eframe::egui::{ComboBox, Ui};
use flowsheet_core::units::columns::{MAX_COLUMN_PORTS, ShortcutColumnSpec, ShortcutReflux};
use flowsheet_core::units::{ColumnSpec, GroupSpec};

use super::{choice, count_field, drawn_only_note, note_row, opt_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`ColumnSpec`]. Feeds and side draws are numbered ports; a
/// connected port cannot be removed (disconnect it first).
pub(super) fn form(ui: &mut Ui, spec: &mut ColumnSpec, u: DisplayUnits, components: &[String]) {
    let kind = spec.kind();
    match spec {
        ColumnSpec::Shortcut(s) => shortcut_form(ui, s, u, components),
        ColumnSpec::Distillation(s) => {
            count_field(ui, "Feeds", &mut s.feeds, 1..=MAX_COLUMN_PORTS);
            count_field(ui, "Side draws", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
        }
        ColumnSpec::PetroleumFractionator(s) => {
            count_field(ui, "Side draws", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
        }
        ColumnSpec::Absorber(_) => {}
    }
    drawn_only_note(ui, kind);
}

fn key_field(ui: &mut Ui, label: &str, key: &mut Option<String>, components: &[String]) {
    let text = ui.label(label);
    let shown = match key.as_deref() {
        None => "Choose a component".to_owned(),
        Some(id) if components.iter().any(|c| c == id) => id.to_owned(),
        Some(id) => format!("{id} (outside package)"),
    };
    ComboBox::from_id_salt(label)
        .selected_text(shown)
        .show_ui(ui, |ui| {
            ui.selectable_value(key, None, "Choose a component");
            for id in components {
                ui.selectable_value(key, Some(id.clone()), id);
            }
        })
        .response
        .labelled_by(text.id);
    ui.end_row();
}

fn shortcut_form(ui: &mut Ui, s: &mut ShortcutColumnSpec, u: DisplayUnits, components: &[String]) {
    key_field(ui, "Light key", &mut s.light_key, components);
    key_field(ui, "Heavy key", &mut s.heavy_key, components);
    opt_field(
        ui,
        "Light key to distillate",
        &mut s.light_key_recovery,
        Quantity::Fraction,
        u,
    );
    opt_field(
        ui,
        "Heavy key to bottoms",
        &mut s.heavy_key_recovery,
        Quantity::Fraction,
        u,
    );
    opt_field(
        ui,
        "Column pressure",
        &mut s.pressure_pa,
        Quantity::Pressure,
        u,
    );
    let current = matches!(s.reflux, ShortcutReflux::Ratio(_)) as usize;
    match choice(
        ui,
        "Reflux",
        current,
        &["Multiple of minimum", "Reflux ratio"],
    ) {
        Some(0) => s.reflux = ShortcutReflux::MinimumMultiple(None),
        Some(_) => s.reflux = ShortcutReflux::Ratio(None),
        None => {}
    }
    match &mut s.reflux {
        ShortcutReflux::MinimumMultiple(m) => {
            opt_field(ui, "R / R_min (> 1)", m, Quantity::Factor, u)
        }
        ShortcutReflux::Ratio(r) => opt_field(ui, "Reflux ratio L/D", r, Quantity::Factor, u),
    }
    note_row(
        ui,
        "Fenske–Underwood–Gilliland estimate. Recoveries lie strictly between 0 and 1 and sum to more than 1.",
    );
    note_row(
        ui,
        "Volatilities come from the package. Keys must be adjacent in volatility.",
    );
    note_row(
        ui,
        "Total condenser and reboiler at column pressure. Both products leave as saturated liquid.",
    );
    note_row(
        ui,
        "Feed components are limited to methane, ethane, propane, nitrogen and CO2 between 200 and 350 K.",
    );
}
