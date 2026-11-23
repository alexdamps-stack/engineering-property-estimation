//! Specification forms of the multiphase group; owned by the multiphase
//! work package.

use eframe::egui::Ui;
use flowsheet_core::units::{GroupSpec, MultiphaseSpec};

use super::drawn_only_note;
use crate::units::DisplayUnits;

/// Form of a [`MultiphaseSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut MultiphaseSpec, _u: DisplayUnits) {
    drawn_only_note(ui, spec.kind());
}
