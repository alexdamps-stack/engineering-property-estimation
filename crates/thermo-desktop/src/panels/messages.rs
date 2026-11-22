//! Message list of the bottom panel.

use eframe::egui::{RichText, ScrollArea, Ui};
use flowsheet_core::{DiagnosticKind, UnitId};

use super::status_style;
use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;

/// One message in the list.
pub(super) struct Message {
    unit: Option<UnitId>,
    status: ItemStatus,
    text: String,
}

pub(super) fn messages_of(ed: &mut Editor) -> Vec<Message> {
    let mut out = Vec::new();
    if ed.results.as_ref().is_some_and(|r| r.recycle.is_some()) {
        let fps = ed.fingerprints();
        let current = ed.results.as_ref().is_some_and(|results| {
            // Results of deleted units have no inputs left to compare.
            results
                .units
                .keys()
                .filter(|id| ed.fs.unit(**id).is_some())
                .all(|id| results.unit_currency(*id, &fps) == flowsheet_project::Currency::Current)
        });
        let r = ed
            .results
            .as_ref()
            .and_then(|r| r.recycle.as_ref())
            .expect("recycle result");
        out.push(Message {
            unit: None,
            status: if current {
                ItemStatus::Solved
            } else {
                ItemStatus::Stale
            },
            text: format!(
                "Recycle {}: {} iterations, max scaled residual {:.2e}; tears {:?}",
                if current { "converged" } else { "result stale" },
                r.iterations,
                r.max_scaled_residual,
                r.tears
                    .iter()
                    .map(|(id, _)| format!("S{id}"))
                    .collect::<Vec<_>>()
            ),
        });
    }
    if let Some(f) = &ed.failure {
        out.push(Message {
            unit: f.unit,
            status: ItemStatus::Failed,
            text: format!("Solve failed: {}", f.message),
        });
    }
    for d in ed.fs.recycle_diagnostics() {
        let unit = match d.kind {
            DiagnosticKind::UnconnectedPort(p) => Some(p.unit),
            DiagnosticKind::MissingSpecification(u, _)
            | DiagnosticKind::InvalidSpecification(u, _) => Some(u),
            DiagnosticKind::RecycleUnsupported(ref c) => c.first().copied(),
            DiagnosticKind::NoFeed => None,
        };
        out.push(Message {
            unit,
            status: ItemStatus::Underspecified,
            text: d.message,
        });
    }
    for i in ed.input_issues() {
        out.push(Message {
            unit: Some(i.unit),
            status: ItemStatus::Underspecified,
            text: i.message,
        });
    }
    let fps = ed.fingerprints();
    if let Some(r) = &ed.results {
        for (id, u) in &r.units {
            if ed.fs.unit(*id).is_none() || u.warnings.is_empty() {
                continue;
            }
            let current = r.unit_currency(*id, &fps) == flowsheet_project::Currency::Current;
            for w in &u.warnings {
                out.push(Message {
                    unit: Some(*id),
                    status: if current {
                        ItemStatus::Warning
                    } else {
                        ItemStatus::Stale
                    },
                    text: format!("{}: {w}", ed.fs.unit(*id).map_or("?", |u| u.name.as_str())),
                });
            }
        }
    }
    out
}

pub(super) fn messages(ui: &mut Ui, ed: &mut Editor, items: Vec<Message>) {
    let theme = Theme::from_visuals(ui.visuals());
    ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        if items.is_empty() {
            ui.label(
                RichText::new("No open items: topology and specifications are complete.").weak(),
            );
        }
        for m in items {
            let (glyph, _, color) = status_style(m.status, &theme);
            let text = RichText::new(format!("{glyph} {}", m.text)).color(color);
            match m.unit {
                Some(unit) => {
                    if ui.link(text).clicked() {
                        ed.selection = Some(Selection::Unit(unit));
                    }
                }
                None => {
                    ui.label(text);
                }
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use flowsheet_core::UnitKind;
    use flowsheet_project::CubicVariantChoice;

    #[test]
    fn recycle_message_ignores_results_of_deleted_units() {
        let mut ed = Editor::default();
        ed.load_project(flowsheet_project::recycle_demonstration(
            CubicVariantChoice::PengRobinsonBostonMathias,
        ));
        ed.start_solve(None);
        ed.wait_for_solve();
        assert!(ed.failure.is_none(), "{:?}", ed.failure);
        let recycle = |ed: &mut Editor| {
            messages_of(ed)
                .into_iter()
                .find(|m| m.text.starts_with("Recycle"))
                .map(|m| m.status)
        };
        assert_eq!(recycle(&mut ed), Some(ItemStatus::Solved));
        // Delete a product sink: no remaining unit's inputs change.
        let product = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Product)
            .unwrap()
            .id;
        ed.selection = Some(Selection::Unit(product));
        ed.delete_selection();
        assert_eq!(recycle(&mut ed), Some(ItemStatus::Solved));
    }
}
