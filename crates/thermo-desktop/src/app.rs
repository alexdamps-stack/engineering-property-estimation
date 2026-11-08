//! Application shell following the main-window layout of
//! `docs/flowsheet-and-ui.md`, section 3: menu and toolbar, palette, canvas,
//! properties panel, bottom panel (messages and stream table), status bar,
//! plus the new-project and unsaved-changes dialogs.

use std::path::PathBuf;

use eframe::egui::{
    self, Id, Key, KeyboardShortcut, Modifiers, Panel, RichText, ScrollArea, Sense, Stroke, Ui,
    ViewportCommand, vec2,
};
use flowsheet_core::UnitKind;
use flowsheet_project::{CubicVariantChoice, PackageSelection, format_utc};

use crate::canvas;
use crate::editor::{Editor, ItemStatus};
use crate::panels::{self, BottomTab};
use crate::symbols::{Theme, draw_unit, icon_rect};
use crate::units::DisplayUnits;

/// File-name suffix proposed for project files.
pub const PROJECT_SUFFIX: &str = ".flowsheet.json";

/// Interaction note recorded when a project is created in the dialog. The
/// dialog offers no interaction data, so the choice is explicit.
pub const DIALOG_INTERACTIONS: &str = "explicit synthetic zero-kij approximation selected in the desktop New project dialog; not measured mixture accuracy";

/// New-project dialog state (proposal section 4.1, step 1).
#[derive(Debug, Clone)]
pub struct NewProjectDialog {
    /// Available components and whether each is selected, in database order.
    pub components: Vec<(String, bool)>,
    pub variant: CubicVariantChoice,
    pub ph_minimum_temperature_k: f64,
    pub ph_maximum_temperature_k: f64,
    pub error: Option<String>,
}

impl NewProjectDialog {
    fn open() -> Self {
        let (components, error) = match PackageSelection::available_components() {
            Ok(ids) => (ids.into_iter().map(|id| (id, false)).collect(), None),
            Err(e) => (Vec::new(), Some(e.to_string())),
        };
        Self {
            components,
            variant: CubicVariantChoice::PengRobinsonBostonMathias,
            ph_minimum_temperature_k: 200.0,
            ph_maximum_temperature_k: 600.0,
            error,
        }
    }

    /// The selection the dialog would create.
    pub fn selection(&self) -> PackageSelection {
        let mut s = PackageSelection::zero_interactions(
            self.variant,
            self.components
                .iter()
                .filter(|(_, on)| *on)
                .map(|(id, _)| id.clone())
                .collect(),
            DIALOG_INTERACTIONS,
        );
        s.ph_minimum_temperature_k = self.ph_minimum_temperature_k;
        s.ph_maximum_temperature_k = self.ph_maximum_temperature_k;
        s
    }
}

/// An action that replaces or closes the project and therefore asks about
/// unsaved changes first.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Pending {
    NewProject,
    Open,
    Demonstration,
    RecycleDemonstration,
    Close,
}

/// The thermo-desktop application.
#[derive(Default)]
pub struct FlowsheetApp {
    editor: Editor,
    bottom_tab: BottomTab,
    new_project: Option<NewProjectDialog>,
    confirm: Option<Pending>,
    allow_close: bool,
    title: String,
}

impl FlowsheetApp {
    /// Starts with the demonstration flowsheet loaded.
    pub fn with_demonstration() -> Self {
        let mut app = Self::default();
        app.editor.load_demonstration();
        app
    }

    /// Starts with the A10 recycle demonstration loaded.
    pub fn with_recycle_demonstration() -> Self {
        let mut app = Self::default();
        app.editor
            .load_project(flowsheet_project::recycle_demonstration(
                CubicVariantChoice::PengRobinsonBostonMathias,
            ));
        app
    }

    /// Starts with a project file opened; failures are shown in the status line.
    pub fn with_file(path: PathBuf) -> Self {
        let mut app = Self::default();
        if let Err(e) = app.editor.open_path(&path) {
            app.editor.status = Some(e);
        }
        app
    }

    /// Editor state (for tests and embedding).
    pub fn editor(&self) -> &Editor {
        &self.editor
    }

    /// Mutable editor state (for tests and embedding).
    pub fn editor_mut(&mut self) -> &mut Editor {
        &mut self.editor
    }

    /// The open new-project dialog, if any (for tests).
    pub fn new_project_dialog(&mut self) -> Option<&mut NewProjectDialog> {
        self.new_project.as_mut()
    }

    fn request(&mut self, action: Pending) {
        if self.editor.is_modified() {
            self.confirm = Some(action);
        } else {
            self.perform(action);
        }
    }

    fn perform(&mut self, action: Pending) {
        match action {
            Pending::NewProject => self.new_project = Some(NewProjectDialog::open()),
            Pending::Demonstration => self.editor.load_demonstration(),
            Pending::RecycleDemonstration => {
                self.editor
                    .load_project(flowsheet_project::recycle_demonstration(
                        CubicVariantChoice::PengRobinsonBostonMathias,
                    ))
            }
            Pending::Open => {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("Flowsheet project", &["json"])
                    .pick_file()
                    && let Err(e) = self.editor.open_path(&path)
                {
                    self.editor.status = Some(e);
                }
            }
            Pending::Close => self.allow_close = true,
        }
    }

    /// Saves to the current file, or asks for one. Returns `true` if saved.
    fn save(&mut self, choose: bool) -> bool {
        let ed = &mut self.editor;
        if ed.package.is_none() {
            ed.status = Some("Create a project with a property package before saving.".into());
            return false;
        }
        let path = match (&ed.path, choose) {
            (Some(p), false) => Some(p.clone()),
            _ => rfd::FileDialog::new()
                .add_filter("Flowsheet project", &["json"])
                .set_file_name(format!("project{PROJECT_SUFFIX}"))
                .save_file(),
        };
        let Some(path) = path else {
            return false;
        };
        match ed.save_to(&path) {
            Ok(()) => true,
            Err(e) => {
                ed.status = Some(e);
                false
            }
        }
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        let shortcut = |m, k| ctx.input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(m, k)));
        if shortcut(Modifiers::COMMAND | Modifiers::SHIFT, Key::S) {
            self.save(true);
        } else if shortcut(Modifiers::COMMAND, Key::S) {
            self.save(false);
        }
        if shortcut(Modifiers::COMMAND, Key::O) {
            self.request(Pending::Open);
        }
        if shortcut(Modifiers::COMMAND, Key::N) {
            self.request(Pending::NewProject);
        }
        if !ctx.egui_wants_keyboard_input() {
            if ctx.input(|i| i.key_pressed(Key::F5)) {
                self.editor.start_solve(Some(ctx.clone()));
            }
            if self.editor.is_solving() && ctx.input(|i| i.key_pressed(Key::Escape)) {
                self.editor.cancel_solve();
            }
        }
    }

    fn update_title(&mut self, ctx: &egui::Context) {
        let name = self
            .editor
            .path
            .as_ref()
            .and_then(|p| p.file_name())
            .map_or_else(
                || "Untitled".to_string(),
                |n| n.to_string_lossy().into_owned(),
            );
        let marker = if self.editor.is_modified() {
            " •"
        } else {
            ""
        };
        let title = format!("{name}{marker} — Property Estimation Flowsheet");
        if title != self.title {
            ctx.send_viewport_cmd(ViewportCommand::Title(title.clone()));
            self.title = title;
        }
    }
}

impl eframe::App for FlowsheetApp {
    fn ui(&mut self, ui: &mut Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        self.editor.poll_solve();
        if ctx.input(|i| i.viewport().close_requested()) && !self.allow_close {
            if self.editor.is_modified() {
                ctx.send_viewport_cmd(ViewportCommand::CancelClose);
                self.confirm = Some(Pending::Close);
            } else {
                self.editor.cancel_solve();
            }
        }
        self.shortcuts(&ctx);
        self.update_title(&ctx);

        Panel::top("menu").show(ui, |ui| self.menu_bar(ui));
        Panel::top("toolbar").show(ui, |ui| self.toolbar(ui));
        Panel::bottom("status_bar").show(ui, |ui| status_bar(ui, &mut self.editor));
        let tab = &mut self.bottom_tab;
        let ed = &mut self.editor;
        Panel::bottom("bottom_panel")
            .resizable(true)
            .default_size(180.0)
            .show(ui, |ui| panels::bottom_panel(ui, ed, tab));
        Panel::left("palette")
            .resizable(false)
            .exact_size(175.0)
            .show(ui, |ui| palette(ui, ed));
        let mut open_dialog = false;
        Panel::right("properties")
            .resizable(true)
            .default_size(320.0)
            .show(ui, |ui| open_dialog = panels::properties(ui, ed));
        egui::CentralPanel::no_frame().show(ui, |ui| canvas::show(ui, ed));
        if open_dialog {
            self.request(Pending::NewProject);
        }

        self.new_project_window(&ctx);
        self.confirm_window(&ctx);

        if self.allow_close {
            ctx.send_viewport_cmd(ViewportCommand::Close);
        }
        // Automatic recalculation waits until no pointer drag is in progress,
        // so dragging a value does not start a solve per frame.
        if !ctx.input(|i| i.pointer.any_down()) && ctx.dragged_id().is_none() {
            self.editor.auto_solve_if_due(Some(ctx.clone()));
        }
        if self.editor.is_solving() {
            ctx.request_repaint_after(std::time::Duration::from_millis(100));
        }
    }
}

impl FlowsheetApp {
    fn menu_bar(&mut self, ui: &mut Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("New project…  (Ctrl+N)").clicked() {
                    self.request(Pending::NewProject);
                }
                if ui.button("Open…  (Ctrl+O)").clicked() {
                    self.request(Pending::Open);
                }
                if ui.button("Save  (Ctrl+S)").clicked() {
                    self.save(false);
                }
                if ui.button("Save as…  (Ctrl+Shift+S)").clicked() {
                    self.save(true);
                }
                ui.separator();
                if ui.button("Load demonstration").clicked() {
                    self.request(Pending::Demonstration);
                }
                if ui.button("Load recycle demonstration").clicked() {
                    self.request(Pending::RecycleDemonstration);
                }
                ui.separator();
                if ui.button("Exit").clicked() {
                    ui.ctx().send_viewport_cmd(ViewportCommand::Close);
                }
            });
            ui.menu_button("Flowsheet", |ui| {
                if ui.button("Solve  (F5)").clicked() {
                    self.editor.start_solve(Some(ui.ctx().clone()));
                }
                if ui
                    .add_enabled(self.editor.is_solving(), egui::Button::new("Cancel  (Esc)"))
                    .clicked()
                {
                    self.editor.cancel_solve();
                }
                ui.checkbox(&mut self.editor.auto_solve, "Auto-solve after edits");
            });
            ui.menu_button("View", |ui| {
                if ui.button("Fit view").clicked() {
                    self.editor.fit_requested = true;
                }
                if ui.button("Copy stream table (tab-separated)").clicked() {
                    let tsv = panels::stream_table_tsv(&mut self.editor);
                    ui.ctx().copy_text(tsv);
                    self.editor.status = Some("Stream table copied as tab-separated text.".into());
                }
            });
        });
    }

    fn toolbar(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            if ui
                .button("New")
                .on_hover_text("New project: choose components and package (Ctrl+N)")
                .clicked()
            {
                self.request(Pending::NewProject);
            }
            if ui
                .button("Open")
                .on_hover_text("Open a project file (Ctrl+O)")
                .clicked()
            {
                self.request(Pending::Open);
            }
            if ui
                .button("Save")
                .on_hover_text("Save the project (Ctrl+S)")
                .clicked()
            {
                self.save(false);
            }
            if ui
                .button("Demonstration")
                .on_hover_text("Mix two gas streams, cool, reduce pressure, separate phases")
                .clicked()
            {
                self.request(Pending::Demonstration);
            }
            ui.separator();
            let ed = &mut self.editor;
            let blocker = ed.solve_blocker();
            let solve = ui.add_enabled(
                blocker.is_none() && !ed.is_solving(),
                egui::Button::new("Solve (F5)"),
            );
            let solve = match &blocker {
                Some(reason) => solve.on_disabled_hover_text(reason),
                None => solve.on_hover_text("Calculate the flowsheet in the background"),
            };
            if solve.clicked() {
                ed.start_solve(Some(ui.ctx().clone()));
            }
            if ui
                .add_enabled(ed.is_solving(), egui::Button::new("Cancel"))
                .on_hover_text("Stop the calculation; results stay as they were (Esc)")
                .clicked()
            {
                ed.cancel_solve();
            }
            ui.checkbox(&mut ed.auto_solve, "Auto-solve")
                .on_hover_text("Recalculate after each completed edit");
            if let Some((done, total, iteration)) = ed.solve_progress() {
                ui.add(egui::Spinner::new());
                ui.label(format!("Solving… {done}/{total} units"));
                if iteration > 0 {
                    ui.label(format!("Recycle iteration {iteration}"));
                }
            }
            ui.separator();
            if ui.button("Fit view").clicked() {
                ed.fit_requested = true;
            }
            ui.label(RichText::new(format!("{:.0} %", ed.viewport.zoom * 100.0)).weak());
            ui.separator();
            let label = ui.label("Units:");
            egui::ComboBox::from_id_salt("display_units")
                .selected_text(ed.units.label())
                .show_ui(ui, |ui| {
                    for u in DisplayUnits::ALL {
                        ui.selectable_value(&mut ed.units, u, u.label());
                    }
                })
                .response
                .labelled_by(label.id);
        });
    }

    fn new_project_window(&mut self, ctx: &egui::Context) {
        let Some(dialog) = &mut self.new_project else {
            return;
        };
        let mut create = None;
        let mut close = false;
        let response = egui::Modal::new(Id::unique("new_project")).show(ctx, |ui| {
            ui.set_width(420.0);
            ui.heading("New project");
            ui.add_space(4.0);
            ui.strong("Components");
            ui.label(
                RichText::new("Order of selection is the database order shown here.")
                    .small()
                    .weak(),
            );
            egui::Grid::new("components").num_columns(2).show(ui, |ui| {
                for (i, (id, on)) in dialog.components.iter_mut().enumerate() {
                    ui.checkbox(on, id.as_str());
                    if i % 2 == 1 {
                        ui.end_row();
                    }
                }
            });
            ui.add_space(6.0);
            ui.strong("Property package");
            ui.horizontal(|ui| {
                for v in CubicVariantChoice::ALL {
                    ui.radio_value(&mut dialog.variant, v, v.label());
                }
            });
            ui.label(
                RichText::new(
                    "Binary interactions: all kij = 0, recorded in the project as an explicit \
                     synthetic approximation. No measured interaction data is applied.",
                )
                .small(),
            );
            ui.horizontal(|ui| {
                ui.label("PH search interval");
                ui.add(
                    egui::DragValue::new(&mut dialog.ph_minimum_temperature_k)
                        .range(200.0..=600.0)
                        .suffix(" K"),
                );
                ui.label("to");
                ui.add(
                    egui::DragValue::new(&mut dialog.ph_maximum_temperature_k)
                        .range(200.0..=600.0)
                        .suffix(" K"),
                );
            });
            if let Some(e) = &dialog.error {
                ui.colored_label(Theme::from_visuals(ui.visuals()).error, e);
            }
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                let any = dialog.components.iter().any(|(_, on)| *on);
                if ui.add_enabled(any, egui::Button::new("Create")).clicked() {
                    let selection = dialog.selection();
                    match selection.build() {
                        Ok(_) => create = Some(selection),
                        Err(e) => dialog.error = Some(e.to_string()),
                    }
                }
                if ui.button("Cancel").clicked() {
                    close = true;
                }
            });
        });
        if let Some(selection) = create {
            self.editor.new_project(selection);
            self.new_project = None;
        } else if close || response.should_close() {
            self.new_project = None;
        }
    }

    fn confirm_window(&mut self, ctx: &egui::Context) {
        let Some(action) = self.confirm.clone() else {
            return;
        };
        let mut choice = None;
        let response = egui::Modal::new(Id::unique("unsaved_changes")).show(ctx, |ui| {
            ui.set_width(360.0);
            ui.heading("Unsaved changes");
            ui.label("The project has changes that are not saved. Save them first?");
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Save").clicked() {
                    choice = Some(0);
                }
                if ui.button("Discard").clicked() {
                    choice = Some(1);
                }
                if ui.button("Cancel").clicked() {
                    choice = Some(2);
                }
            });
        });
        match choice {
            Some(0) => {
                self.confirm = None;
                if self.save(false) {
                    self.perform(action);
                }
            }
            Some(1) => {
                self.confirm = None;
                self.perform(action);
            }
            Some(_) => self.confirm = None,
            None if response.should_close() => self.confirm = None,
            None => {}
        }
    }
}

fn status_bar(ui: &mut Ui, ed: &mut Editor) {
    let fps = ed.fingerprints();
    let issues = ed.input_issues();
    let theme = Theme::from_visuals(ui.visuals());
    let counts = ed.fs.units().fold([0usize; 4], |mut c, u| {
        match ed.unit_status(u.id, &fps, &issues) {
            ItemStatus::Solved => c[0] += 1,
            ItemStatus::Warning => c[1] += 1,
            ItemStatus::Failed => c[2] += 1,
            ItemStatus::Stale => c[3] += 1,
            _ => {}
        }
        c
    });
    ui.horizontal(|ui| {
        match &ed.package {
            Some(p) => {
                let data = ed
                    .package_cache
                    .get(p)
                    .provenance
                    .as_ref()
                    .map_or_else(|| "data ?".to_string(), |v| format!("data {}", v.data_version));
                ui.label(format!(
                    "Package: {} · {} components · {data}",
                    p.variant.label(),
                    p.components.len()
                ))
                .on_hover_text(&p.interaction_source);
            }
            None => {
                ui.label(RichText::new("Package: none — File > New project").weak());
            }
        }
        ui.separator();
        if let Some((done, total, iteration)) = ed.solve_progress() {
            ui.label(format!("Solving {done}/{total}…"));
            if iteration > 0 { ui.label(format!("Recycle iteration {iteration}")); }
        } else if counts[2] > 0 || ed.failure.is_some() {
            ui.colored_label(theme.error, "✖ Last solve failed — see Messages");
        } else if let Some(r) = &ed.results {
            let when = format_utc(r.solved_at_unix_s);
            if counts[3] > 0 {
                ui.label(
                    RichText::new(format!("Stale: {} unit(s) changed since {when}", counts[3]))
                        .italics(),
                );
            } else {
                ui.label(format!("Solved {when}"));
            }
        } else {
            ui.label(RichText::new("Not solved").weak());
        }
        ui.separator();
        ui.label(
            RichText::new(format!(
                "{} units, {} streams · {} solved, {} warning, {} failed, {} stale · {} open item(s)",
                ed.fs.units().count(),
                ed.fs.streams().count(),
                counts[0],
                counts[1],
                counts[2],
                counts[3],
                ed.fs.recycle_diagnostics().len() + issues.len()
            ))
            .weak(),
        );
    });
}

fn palette(ui: &mut Ui, ed: &mut Editor) {
    ScrollArea::vertical().show(ui, |ui| palette_items(ui, ed));
}

fn palette_items(ui: &mut Ui, ed: &mut Editor) {
    ui.add_space(6.0);
    ui.strong("Unit operations");
    ui.add_space(4.0);
    for kind in UnitKind::ALL {
        let (rect, response) =
            ui.allocate_exact_size(vec2(ui.available_width(), 36.0), Sense::click());
        let visuals = ui.style().interact(&response);
        if response.hovered() {
            ui.painter().rect_filled(rect, 4, visuals.weak_bg_fill);
        }
        let icon_area = egui::Rect::from_min_size(rect.min + vec2(8.0, 5.0), vec2(36.0, 26.0));
        draw_unit(
            ui.painter(),
            kind,
            icon_rect(icon_area, kind),
            Stroke::new(1.2, visuals.text_color()),
            egui::Color32::TRANSPARENT,
        );
        ui.painter().text(
            rect.left_center() + vec2(54.0, 0.0),
            egui::Align2::LEFT_CENTER,
            kind.label(),
            egui::FontId::proportional(14.0),
            visuals.text_color(),
        );
        let label = kind.label();
        response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, label));
        if response
            .on_hover_text(format!("Place a {}", kind.label().to_lowercase()))
            .clicked()
        {
            let hint = ed
                .content_bounds()
                .map(|b| egui::pos2(b.right() + 100.0, b.center().y))
                .unwrap_or(egui::Pos2::ZERO);
            ed.add_unit_near(kind, hint);
        }
    }
    ui.add_space(12.0);
    ui.separator();
    ui.label(
        RichText::new(
            "Drag from an outlet to an inlet to connect.\n\
             Drag blocks to move them.\n\
             Drag the background or scroll to pan; Ctrl+scroll to zoom.\n\
             Del removes the selection. F5 solves.",
        )
        .small()
        .weak(),
    );
}
