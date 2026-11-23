//! Specification forms of the solids group; owned by the solids work
//! package.

use eframe::egui::Ui;
use flowsheet_core::units::{GroupSpec, SolidsSpec};

use super::drawn_only_note;
use crate::units::DisplayUnits;

/// Form of a [`SolidsSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut SolidsSpec, _u: DisplayUnits) {
    drawn_only_note(ui, spec.kind());
}
