//! Specification forms of the reactors group; owned by the reactors work
//! package.

use eframe::egui::Ui;
use flowsheet_core::units::{GroupSpec, ReactorSpec};

use super::drawn_only_note;
use crate::units::DisplayUnits;

/// Form of a [`ReactorSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut ReactorSpec, _u: DisplayUnits) {
    drawn_only_note(ui, spec.kind());
}
