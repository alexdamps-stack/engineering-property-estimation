//! Specification forms of the simple/logical group; owned by the
//! simple/logical work package.

use eframe::egui::Ui;
use flowsheet_core::units::simple_logical::MAX_SELECTOR_INLETS;
use flowsheet_core::units::{GroupSpec, SimpleLogicalSpec};

use super::{count_field, drawn_only_note};
use crate::units::DisplayUnits;

/// Form of a [`SimpleLogicalSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut SimpleLogicalSpec, _u: DisplayUnits) {
    let kind = spec.kind();
    if let SimpleLogicalSpec::Selector(s) = spec {
        count_field(ui, "Inlets", &mut s.inlets, 2..=MAX_SELECTOR_INLETS);
    }
    drawn_only_note(ui, kind);
}
