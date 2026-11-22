//! Shared result presentation model: which stream properties are shown, how
//! each is read from a stored [`StreamResult`], and the per-stream columns of
//! the stream table. The stream properties panel, the stream table and its
//! tab-separated export all read results through this module.

use flowsheet_core::StreamId;
use flowsheet_project::StreamResult;

use crate::editor::{Editor, ItemStatus};
use crate::units::{DisplayUnits, DisplayUnitsExt, Quantity, format_quantity, format_value};

/// A scalar stream property shown in the properties panel and the stream table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum StreamProperty {
    Temperature,
    Pressure,
    MolarFlow,
    MassFlow,
    VaporFraction,
    MolarEnthalpy,
    MolarEntropy,
}

impl StreamProperty {
    /// Row order of the stream table and its export.
    pub(super) const TABLE: [StreamProperty; 7] = [
        StreamProperty::Temperature,
        StreamProperty::Pressure,
        StreamProperty::MolarFlow,
        StreamProperty::MassFlow,
        StreamProperty::VaporFraction,
        StreamProperty::MolarEnthalpy,
        StreamProperty::MolarEntropy,
    ];

    /// Flow rows of the stream properties panel, shown for every result.
    pub(super) const PANEL_FLOWS: [StreamProperty; 2] =
        [StreamProperty::MolarFlow, StreamProperty::MassFlow];

    /// State rows of the stream properties panel, shown when the stream has
    /// a thermodynamic state.
    pub(super) const PANEL_STATE: [StreamProperty; 5] = [
        StreamProperty::Temperature,
        StreamProperty::Pressure,
        StreamProperty::VaporFraction,
        StreamProperty::MolarEnthalpy,
        StreamProperty::MolarEntropy,
    ];

    /// Row label.
    pub(super) fn label(self) -> &'static str {
        match self {
            StreamProperty::Temperature => "Temperature",
            StreamProperty::Pressure => "Pressure",
            StreamProperty::MolarFlow => "Molar flow",
            StreamProperty::MassFlow => "Mass flow",
            StreamProperty::VaporFraction => "Vapor fraction",
            StreamProperty::MolarEnthalpy => "Molar enthalpy",
            StreamProperty::MolarEntropy => "Molar entropy",
        }
    }

    /// Quantity that sets the display unit; `None` for dimensionless values
    /// shown as stored.
    pub(super) fn quantity(self) -> Option<Quantity> {
        match self {
            StreamProperty::Temperature => Some(Quantity::Temperature),
            StreamProperty::Pressure => Some(Quantity::Pressure),
            StreamProperty::MolarFlow => Some(Quantity::MolarFlow),
            StreamProperty::MassFlow => Some(Quantity::MassFlow),
            StreamProperty::VaporFraction => None,
            StreamProperty::MolarEnthalpy => Some(Quantity::MolarEnthalpy),
            StreamProperty::MolarEntropy => Some(Quantity::MolarEntropy),
        }
    }

    /// Unit symbol in the display set, or `–` for dimensionless values.
    pub(super) fn unit_symbol(self, u: DisplayUnits) -> &'static str {
        self.quantity().map_or("–", |q| u.symbol(q))
    }

    /// Stored SI value, if the result has one.
    pub(super) fn si_value(self, r: &StreamResult) -> Option<f64> {
        let state = r.state.as_ref();
        match self {
            StreamProperty::Temperature => state.map(|s| s.temperature_k),
            StreamProperty::Pressure => state.map(|s| s.pressure_pa),
            StreamProperty::MolarFlow => Some(r.molar_flow_mol_s),
            StreamProperty::MassFlow => Some(r.mass_flow_kg_s),
            StreamProperty::VaporFraction => state.and_then(|s| s.vapor_fraction),
            StreamProperty::MolarEnthalpy => state.map(|s| s.enthalpy_j_mol),
            StreamProperty::MolarEntropy => state.map(|s| s.entropy_j_mol_k),
        }
    }

    /// Value converted to the display set, without a unit.
    pub(super) fn display_value(self, r: &StreamResult, u: DisplayUnits) -> Option<f64> {
        self.si_value(r)
            .map(|v| self.quantity().map_or(v, |q| u.si_to_display(q, v)))
    }

    /// Value with its unit as shown in the stream properties panel.
    pub(super) fn panel_text(self, r: &StreamResult, u: DisplayUnits) -> String {
        match (self.si_value(r), self.quantity()) {
            (Some(v), Some(q)) => format_quantity(u, q, v),
            (Some(v), None) => format_value(v),
            (None, _) if self == StreamProperty::VaporFraction => "— (supercritical)".into(),
            (None, _) => "—".into(),
        }
    }
}

/// One stream of the stream table.
pub(super) struct StreamColumn {
    pub(super) name: String,
    pub(super) status: Option<ItemStatus>,
    pub(super) result: Option<StreamResult>,
}

/// Every material stream in flowsheet order with its status and, when it has a
/// status, its stored result.
pub(super) fn stream_columns(ed: &mut Editor) -> Vec<(StreamId, StreamColumn)> {
    let fps = ed.fingerprints();
    // Material streams only; energy streams are shown in their own panel.
    let streams: Vec<_> = ed
        .fs
        .material_streams()
        .map(|s| (s.id, s.name.clone()))
        .collect();
    streams
        .into_iter()
        .map(|(id, name)| {
            let status = ed.stream_status(id, &fps);
            let result = status.and(
                ed.results
                    .as_ref()
                    .and_then(|r| r.streams.get(&id).cloned()),
            );
            (
                id,
                StreamColumn {
                    name,
                    status,
                    result,
                },
            )
        })
        .collect()
}

/// Component names of the current package, in package order.
pub(super) fn component_names(ed: &Editor) -> Vec<String> {
    ed.package
        .as_ref()
        .map(|p| p.components.clone())
        .unwrap_or_default()
}

/// Mole fraction `i` of a stored result as display text. Results that do not
/// match the package component list (for example from a hand-edited file)
/// show a dash instead of panicking.
pub(super) fn fraction_text(fractions: &[f64], i: usize) -> String {
    fractions
        .get(i)
        .map_or_else(|| "—".to_string(), |x| format_value(*x))
}

/// One-line summary of a stream result: temperature, pressure and molar flow.
pub(super) fn stream_summary(s: &StreamResult, u: DisplayUnits) -> String {
    match &s.state {
        Some(state) => format!(
            "{}, {}, {}",
            format_quantity(u, Quantity::Temperature, state.temperature_k),
            format_quantity(u, Quantity::Pressure, state.pressure_pa),
            format_quantity(u, Quantity::MolarFlow, s.molar_flow_mol_s)
        ),
        None => "no flow (absent phase)".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn out_of_range_mole_fractions_render_as_a_dash() {
        assert_eq!(fraction_text(&[0.25, 0.75], 1), format_value(0.75));
        assert_eq!(fraction_text(&[0.25], 3), "—");
    }

    #[test]
    fn panel_rows_cover_every_table_row_once() {
        let mut panel: Vec<_> = StreamProperty::PANEL_FLOWS
            .iter()
            .chain(&StreamProperty::PANEL_STATE)
            .map(|p| p.label())
            .collect();
        let mut table: Vec<_> = StreamProperty::TABLE.iter().map(|p| p.label()).collect();
        panel.sort_unstable();
        table.sort_unstable();
        assert_eq!(panel, table);
    }
}
