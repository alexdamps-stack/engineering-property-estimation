//! Specification forms of the columns group; owned by the columns work
//! package.

use eframe::egui::Ui;
use flowsheet_core::units::columns::MAX_COLUMN_PORTS;
use flowsheet_core::units::{ColumnSpec, GroupSpec};

use super::{count_field, drawn_only_note};
use crate::units::DisplayUnits;

/// Form of a [`ColumnSpec`]. Feeds and side draws are numbered ports; a
/// connected port cannot be removed (disconnect it first).
pub(super) fn form(ui: &mut Ui, spec: &mut ColumnSpec, _u: DisplayUnits) {
    let kind = spec.kind();
    match spec {
        ColumnSpec::Distillation(s) => {
            count_field(ui, "Feeds", &mut s.feeds, 1..=MAX_COLUMN_PORTS);
            count_field(ui, "Side draws", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
        }
        ColumnSpec::PetroleumFractionator(s) => {
            count_field(ui, "Side draws", &mut s.side_draws, 0..=MAX_COLUMN_PORTS);
        }
        ColumnSpec::Shortcut(_) | ColumnSpec::Absorber(_) => {}
    }
    drawn_only_note(ui, kind);
}
