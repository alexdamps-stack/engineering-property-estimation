//! Specification forms of the simple/logical group; owned by the
//! simple/logical work package.

use eframe::egui::Ui;
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
        SimpleLogicalSpec::Selector(s) => {
            count_field(ui, "Inlets", &mut s.inlets, 2..=MAX_SELECTOR_INLETS);
        }
        _ => {}
    }
    drawn_only_note(ui, kind);
}
