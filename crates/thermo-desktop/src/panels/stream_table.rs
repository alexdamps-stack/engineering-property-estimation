//! Stream table (proposal section 4.2) and its tab-separated export.

use eframe::egui::{self, RichText, ScrollArea, Ui};

use super::presentation::{
    StreamProperty, component_names, fraction_text, fraction_value, stream_columns,
};
use super::status_style;
use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;
use crate::units::format_value;

/// Stream table as tab-separated text with display units, a status row and
/// one row per property; empty cells where no value exists.
pub fn stream_table_tsv(ed: &mut Editor) -> String {
    let u = ed.units;
    let fps = ed.fingerprints();
    let cols = stream_columns(ed, &fps);
    let components = component_names(ed);
    let mut out = String::from("Stream\tUnit");
    for (_, c) in &cols {
        out.push('\t');
        out.push_str(&tsv_cell(c.name));
    }
    out.push_str("\nStatus\t");
    for (_, c) in &cols {
        out.push('\t');
        out.push_str(match c.status {
            Some(ItemStatus::Solved) => "solved",
            Some(ItemStatus::Warning) => "warning",
            Some(ItemStatus::Stale) => "stale",
            _ => "",
        });
    }
    for property in StreamProperty::TABLE {
        out.push('\n');
        out.push_str(property.label());
        out.push('\t');
        out.push_str(property.unit_symbol(u));
        for (_, c) in &cols {
            out.push('\t');
            if let Some(shown) = c.result.and_then(|r| property.display_value(r, u))
                && shown.is_finite()
            {
                out.push_str(&shown.to_string());
            }
        }
    }
    for (i, name) in components.iter().enumerate() {
        out.push_str(&format!("\nx {}\tmol/mol", tsv_cell(name)));
        for (_, c) in &cols {
            out.push('\t');
            if let Some(x) = c
                .result
                .and_then(|r| r.state.as_ref())
                .and_then(|s| fraction_value(&s.mole_fractions, i))
            {
                out.push_str(&x.to_string());
            }
        }
    }
    out.push('\n');
    out
}

/// Streams as columns, properties as rows (proposal section 4.2).
pub(super) fn stream_table(ui: &mut Ui, ed: &mut Editor) {
    let fps = ed.fingerprints();
    let cols = stream_columns(ed, &fps);
    if cols.is_empty() {
        ui.label(RichText::new("No streams yet.").weak());
        return;
    }
    let theme = Theme::from_visuals(ui.visuals());
    let u = ed.units;
    let components = component_names(ed);
    let mut selection = ed.selection;
    ScrollArea::both().auto_shrink(false).show(ui, |ui| {
        egui::Grid::new("stream_table")
            .striped(true)
            .min_col_width(70.0)
            .show(ui, |ui| {
                ui.label("");
                ui.label("");
                for (id, c) in &cols {
                    let selected = selection == Some(Selection::Stream(*id));
                    if ui
                        .selectable_label(selected, RichText::new(c.name).strong())
                        .clicked()
                    {
                        selection = Some(Selection::Stream(*id));
                    }
                }
                ui.end_row();
                for label in ["From", "To"] {
                    ui.label(label);
                    ui.label("");
                    for (id, _) in &cols {
                        let s = ed.fs.stream(*id).expect("listed stream");
                        let port = if label == "From" { s.from } else { s.to };
                        let name = ed
                            .fs
                            .unit(port.unit)
                            .map(|u| u.name.as_str())
                            .unwrap_or("?");
                        ui.label(RichText::new(name).weak());
                    }
                    ui.end_row();
                }
                ui.label("Status");
                ui.label("");
                for (_, c) in &cols {
                    match c.status {
                        Some(s) => {
                            let (glyph, label, color) = status_style(s, &theme);
                            let short = if s == ItemStatus::Stale {
                                "stale"
                            } else {
                                label
                            };
                            ui.label(RichText::new(format!("{glyph} {short}")).color(color));
                        }
                        None => {
                            ui.label(RichText::new("—").weak());
                        }
                    }
                }
                ui.end_row();
                for property in StreamProperty::TABLE {
                    ui.label(property.label());
                    ui.label(RichText::new(property.unit_symbol(u)).weak());
                    for (_, c) in &cols {
                        let stale = c.status == Some(ItemStatus::Stale);
                        match c.result.and_then(|r| property.display_value(r, u)) {
                            Some(v) => {
                                let text = format_value(v);
                                ui.label(if stale {
                                    RichText::new(text).italics().weak()
                                } else {
                                    RichText::new(text)
                                });
                            }
                            None => {
                                ui.label(RichText::new("—").weak());
                            }
                        }
                    }
                    ui.end_row();
                }
                for (i, name) in components.iter().enumerate() {
                    ui.label(format!("x {name}"));
                    ui.label(RichText::new("mol/mol").weak());
                    for (_, c) in &cols {
                        let stale = c.status == Some(ItemStatus::Stale);
                        match c.result.and_then(|r| r.state.as_ref()) {
                            Some(s) => {
                                let text = fraction_text(&s.mole_fractions, i);
                                ui.label(if stale {
                                    RichText::new(text).italics().weak()
                                } else {
                                    RichText::new(text)
                                });
                            }
                            None => {
                                ui.label(RichText::new("—").weak());
                            }
                        }
                    }
                    ui.end_row();
                }
            });
    });
    ed.selection = selection;
}

/// A user-entered name as one tab-separated cell: tabs and line breaks would
/// split the cell or the row, so they become spaces.
fn tsv_cell(name: &str) -> String {
    name.replace(['\t', '\n', '\r'], " ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_cells_keep_one_cell_per_name() {
        assert_eq!(tsv_cell("a\tb\nc\r\nd"), "a b c  d");
        assert_eq!(tsv_cell("plain"), "plain");
    }

    #[test]
    fn copied_stream_table_survives_names_with_tabs_and_newlines() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        let id = ed.fs.streams().next().unwrap().id;
        ed.fs.rename_stream(id, "bad\tname\nwith\rbreaks");
        let tsv = stream_table_tsv(&mut ed);
        let lines: Vec<&str> = tsv.lines().collect();
        let columns = lines[0].split('\t').count();
        assert_eq!(columns, 2 + ed.fs.streams().count());
        assert!(lines[0].contains("bad name with breaks"));
        assert!(lines.iter().all(|l| l.split('\t').count() == columns));
    }

    #[test]
    fn stream_table_with_results_shorter_than_the_component_list_does_not_panic() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        ed.start_solve(None);
        ed.wait_for_solve();
        // A package with one component more than the stored results, as in a
        // hand-edited file whose results were left in place.
        ed.package.as_mut().unwrap().components.push("extra".into());
        let tsv = stream_table_tsv(&mut ed);
        assert!(tsv.contains("\nx extra\tmol/mol"), "{tsv}");
    }
}
