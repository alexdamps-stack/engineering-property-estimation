//! Properties panel (specifications, feed composition, results), message
//! list, and stream table (proposal sections 3, 4.2 and 4.6).
//!
//! One submodule per view; `presentation` is the shared result model that
//! the stream views and the stream table export read results through.

use eframe::egui::{Color32, RichText, ScrollArea, Ui};

use crate::editor::{Editor, ItemStatus, Selection};
use crate::symbols::Theme;

mod messages;
mod presentation;
mod project;
mod spec;
mod stream;
mod stream_table;
mod unit;

use messages::{messages, messages_of};
use project::project_properties;
use stream::stream_properties;
use stream_table::stream_table;
pub use stream_table::stream_table_tsv;
use unit::unit_properties;

/// Tabs of the bottom panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BottomTab {
    #[default]
    Messages,
    StreamTable,
}

/// Glyph, label and colour of a status.
pub fn status_style(status: ItemStatus, theme: &Theme) -> (&'static str, &'static str, Color32) {
    match status {
        ItemStatus::Underspecified => ("!", "Underspecified", theme.warning),
        ItemStatus::Ready => ("-", "Ready (not solved)", theme.muted),
        ItemStatus::Solved => ("✔", "Solved", theme.ok),
        ItemStatus::Warning => ("⚠", "Warning", theme.warning),
        ItemStatus::Failed => ("✖", "Failed", theme.error),
        ItemStatus::Stale => (
            "~",
            "Stale: inputs changed since the last solve",
            theme.muted,
        ),
    }
}

fn status_line(ui: &mut Ui, status: ItemStatus) {
    let theme = Theme::from_visuals(ui.visuals());
    let (glyph, label, color) = status_style(status, &theme);
    ui.horizontal(|ui| {
        ui.label("Status:");
        let text = RichText::new(format!("{glyph} {label}")).color(color);
        ui.label(if status == ItemStatus::Stale {
            text.italics()
        } else {
            text
        });
    });
}

/// A result value; stale values are grey italic and marked.
fn value_text(text: String, stale: bool) -> RichText {
    if stale {
        RichText::new(format!("{text} (stale)")).italics().weak()
    } else {
        RichText::new(text)
    }
}

/// Properties panel. Returns `true` if the user asked for a new project.
pub fn properties(ui: &mut Ui, ed: &mut Editor) -> bool {
    let mut new_project = false;
    ScrollArea::vertical().show(ui, |ui| {
        ui.add_space(6.0);
        match ed.selection {
            Some(Selection::Unit(id)) if ed.fs.unit(id).is_some() => unit_properties(ui, ed, id),
            Some(Selection::Stream(id)) if ed.fs.stream(id).is_some() => {
                stream_properties(ui, ed, id)
            }
            _ => new_project = project_properties(ui, ed),
        }
    });
    new_project
}

/// Bottom panel with messages and the stream table.
pub fn bottom_panel(ui: &mut Ui, ed: &mut Editor, tab: &mut BottomTab) {
    let items = messages_of(ed);
    ui.horizontal(|ui| {
        ui.selectable_value(
            tab,
            BottomTab::Messages,
            format!("Messages ({})", items.len()),
        );
        ui.selectable_value(tab, BottomTab::StreamTable, "Stream table");
        if *tab == BottomTab::StreamTable
            && ui
                .small_button("Copy")
                .on_hover_text("Copy the stream table as tab-separated text for spreadsheets")
                .clicked()
        {
            let tsv = stream_table_tsv(ed);
            ui.ctx().copy_text(tsv);
            ed.status = Some("Stream table copied as tab-separated text.".into());
        }
        if let Some(status) = &ed.status {
            ui.separator();
            ui.colored_label(Theme::from_visuals(ui.visuals()).warning, status);
            if ui.small_button("×").clicked() {
                ed.status = None;
            }
        }
    });
    ui.separator();
    match tab {
        BottomTab::Messages => messages(ui, ed, items),
        BottomTab::StreamTable => stream_table(ui, ed),
    }
}
