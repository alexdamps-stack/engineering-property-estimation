//! Specification forms of the heat-transport group; owned by the
//! heat-transport work package.

use eframe::egui::Ui;
use flowsheet_core::units::heat_transport::MAX_PASSES;
use flowsheet_core::units::{GroupSpec, HeatTransportSpec};

use super::{count_field, drawn_only_note, note_row};
use crate::units::DisplayUnits;

/// Form of a [`HeatTransportSpec`].
pub(super) fn form(ui: &mut Ui, spec: &mut HeatTransportSpec, _u: DisplayUnits) {
    let kind = spec.kind();
    match spec {
        HeatTransportSpec::MultiStreamExchanger(s) => {
            count_field(ui, "Passes", &mut s.passes, 2..=MAX_PASSES);
            note_row(ui, "Each pass k has the ports in k and out k.");
        }
        HeatTransportSpec::AirCooler(_)
        | HeatTransportSpec::FiredHeater(_)
        | HeatTransportSpec::PipeSegment(_)
        | HeatTransportSpec::ReliefValve(_) => {}
    }
    drawn_only_note(ui, kind);
}
