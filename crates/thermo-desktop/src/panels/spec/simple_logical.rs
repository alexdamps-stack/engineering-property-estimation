//! Specification forms of the simple/logical group; owned by the
//! simple/logical work package.

use eframe::egui::{ComboBox, Ui};
use flowsheet_core::units::simple_logical::MAX_SELECTOR_INLETS;
use flowsheet_core::units::{GroupSpec, SimpleLogicalSpec};

use super::{count_field, drawn_only_note, note_row, opt_field};
use crate::units::{DisplayUnits, Quantity};

/// Form of a [`SimpleLogicalSpec`].
pub(super) fn form(
    ui: &mut Ui,
    spec: &mut SimpleLogicalSpec,
    u: DisplayUnits,
    components: &[String],
) {
    let kind = spec.kind();
    match spec {
        SimpleLogicalSpec::ComponentSplitter(s) => {
            note_row(
                ui,
                "Enter each component's fraction to overhead. Bottoms receives the rest.",
            );
            for id in components {
                let before = s.overhead_recovery.get(id).copied();
                let mut value = before;
                opt_field(ui, id, &mut value, Quantity::Fraction, u);
                match value {
                    Some(v) if Some(v) != before => {
                        s.overhead_recovery.insert(id.clone(), v);
                    }
                    None if before.is_some() => {
                        s.overhead_recovery.remove(id);
                    }
                    _ => {}
                }
            }
            let unknown: Vec<_> = s
                .overhead_recovery
                .keys()
                .filter(|id| !components.contains(id))
                .cloned()
                .collect();
            for id in unknown {
                ui.label(format!("{id} (outside package)"));
                if ui.small_button("Remove recovery").clicked() {
                    s.overhead_recovery.remove(&id);
                }
                ui.end_row();
            }
            note_row(
                ui,
                "Both outlets use the inlet temperature and pressure. The splitter calculates the duty.",
            );
        }
        SimpleLogicalSpec::Multiplier(s) => {
            opt_field(ui, "Flow factor", &mut s.factor, Quantity::Factor, u);
            note_row(
                ui,
                "Scales flow at the inlet state. Results show the material and enthalpy added or removed.",
            );
        }
        SimpleLogicalSpec::Selector(s) => {
            if count_field(ui, "Inlets", &mut s.inlets, 2..=MAX_SELECTOR_INLETS)
                && s.selected_inlet.is_some_and(|n| n > s.inlets)
            {
                s.selected_inlet = None;
            }
            let label = ui.label("Selected inlet");
            ComboBox::from_id_salt("selected_inlet")
                .selected_text(
                    s.selected_inlet
                        .map_or_else(|| "Choose an inlet".into(), |n| format!("Inlet {n}")),
                )
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut s.selected_inlet, None, "Choose an inlet");
                    for n in 1..=s.inlets.clamp(2, MAX_SELECTOR_INLETS) {
                        ui.selectable_value(&mut s.selected_inlet, Some(n), format!("Inlet {n}"));
                    }
                })
                .response
                .labelled_by(label.id);
            ui.end_row();
            note_row(
                ui,
                "Forwards the selected inlet. Results show the material and enthalpy removed with the other inlets.",
            );
        }
        _ => {}
    }
    drawn_only_note(ui, kind);
}
