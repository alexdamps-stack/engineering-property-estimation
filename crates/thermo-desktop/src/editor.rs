//! Editor state: the project inputs (flowsheet, package, feed compositions),
//! results and their currency, plus GUI-only layout, viewport and selection.

use std::cell::RefCell;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use eframe::egui::{Pos2, Rect, Vec2, pos2};
use flowsheet_core::{Flowsheet, PortRef, StreamId, UnitId, UnitKind, UnitStatus};
use flowsheet_project::{
    CubicVariantChoice, Currency, DisplayUnitSet, FeedFractions, Fingerprint, Fingerprints,
    InputIssue, PackageSelection, Project, ResultSet,
};

use crate::solve::{PackageCache, SolveJob, SolveOutcome};
use crate::units::DisplayUnits;
use crate::view::{Viewport, block_rect};

/// The selected canvas element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Selection {
    Unit(UnitId),
    Stream(StreamId),
}

/// An in-progress pointer drag on the canvas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Drag {
    /// Panning the view.
    Panning,
    /// Moving a unit; `grab` is the pointer offset from the unit centre (world).
    Moving { unit: UnitId, grab: Vec2 },
    /// Drawing a new stream from a free port.
    Connecting { from: PortRef },
}

/// Status of a unit or stream as shown to the user (proposal section 4.6).
/// Exactly one applies; values are never shown without one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ItemStatus {
    /// A specification, composition or connection is missing or invalid.
    Underspecified,
    /// Complete, but no result exists.
    Ready,
    /// Current result without warnings.
    Solved,
    /// Current result with at least one warning.
    Warning,
    /// The last solve failed at this unit and its inputs are unchanged.
    Failed,
    /// A result exists but inputs changed since it was calculated.
    Stale,
}

/// The most recent failed solve attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Failure {
    /// Unit the failure is attributed to, if any.
    pub unit: Option<UnitId>,
    /// Fingerprint of that unit's inputs when it failed.
    pub unit_fingerprint: Option<Fingerprint>,
    /// User-facing cause.
    pub message: String,
    /// Combined fingerprint of all inputs of the failed attempt.
    pub inputs: Fingerprint,
}

/// Inputs the modified state depends on besides the saved text, so the
/// project is serialized only when one of them changed rather than every
/// frame. The flowsheet is represented by [`Flowsheet::revision`], which
/// every mutation advances; a flowsheet replaced wholesale must go through
/// [`Editor::load_project`] or [`Editor::clear`], which reset the caches.
struct ModifiedCache {
    revision: u64,
    package: Option<PackageSelection>,
    feeds: FeedFractions,
    layout: HashMap<UnitId, Pos2>,
    units: DisplayUnits,
    results: Option<ResultSet>,
    modified: bool,
}

impl ModifiedCache {
    fn matches(&self, ed: &Editor) -> bool {
        self.revision == ed.fs.revision()
            && self.units == ed.units
            && self.package == ed.package
            && self.feeds == ed.feeds
            && self.layout == ed.layout
            && self.results == ed.results
    }
}

/// Fingerprints together with the inputs they were computed from.
struct FingerprintCache {
    revision: u64,
    package: Option<PackageSelection>,
    feeds: FeedFractions,
    fingerprints: Fingerprints,
}

/// Complete editor state.
pub struct Editor {
    pub fs: Flowsheet,
    /// Property package of the project; `None` until a project is created.
    pub package: Option<PackageSelection>,
    /// Feed mole fractions in package component order.
    pub feeds: FeedFractions,
    /// Most recent successful results (check currency before showing them).
    pub results: Option<ResultSet>,
    /// Most recent failed solve, if it happened after the last success.
    pub failure: Option<Failure>,
    /// Recalculate automatically after each completed edit.
    pub auto_solve: bool,
    /// File the project was last saved to or opened from.
    pub path: Option<PathBuf>,
    /// World-space centre of each unit block. Layout is GUI-only data.
    pub layout: HashMap<UnitId, Pos2>,
    pub viewport: Viewport,
    pub selection: Option<Selection>,
    pub drag: Option<Drag>,
    pub units: DisplayUnits,
    /// Fit the view to the content on the next canvas frame.
    pub fit_requested: bool,
    /// Last user-facing status message (e.g. a rejected connection).
    pub status: Option<String>,
    /// Screen rectangle of the canvas in the last frame (used by UI tests to
    /// convert drawing coordinates to pointer positions).
    pub canvas_rect: Rect,
    pub(crate) job: Option<SolveJob>,
    pub(crate) package_cache: PackageCache,
    /// Combined input fingerprint of the last solve attempt (any outcome).
    pub(crate) last_attempt: Option<Fingerprint>,
    /// Serialized project at the last save or open, to detect changes. Set it
    /// through [`Self::set_saved`] so the modified-state cache is dropped.
    pub(crate) saved: Option<String>,
    /// Cached [`Self::is_modified`] answer and the inputs it was computed for.
    modified_cache: RefCell<Option<ModifiedCache>>,
    /// Cached [`Self::fingerprints`] and the inputs they were computed for.
    fingerprint_cache: Option<FingerprintCache>,
}

impl Default for Editor {
    fn default() -> Self {
        Self {
            fs: Flowsheet::new(),
            package: None,
            feeds: FeedFractions::new(),
            results: None,
            failure: None,
            auto_solve: false,
            path: None,
            layout: HashMap::new(),
            viewport: Viewport::default(),
            selection: None,
            drag: None,
            units: DisplayUnits::default(),
            fit_requested: false,
            status: None,
            canvas_rect: Rect::NOTHING,
            job: None,
            package_cache: PackageCache::default(),
            last_attempt: None,
            saved: None,
            modified_cache: RefCell::new(None),
            fingerprint_cache: None,
        }
    }
}

impl Editor {
    /// Replaces the project with an empty one without a package, keeping
    /// display preferences.
    pub fn clear(&mut self) {
        *self = Self {
            units: self.units,
            auto_solve: self.auto_solve,
            ..Self::default()
        };
    }

    /// Starts an empty project with the given package.
    pub fn new_project(&mut self, package: PackageSelection) {
        self.clear();
        self.package = Some(package);
    }

    /// Adds a unit at `center` (world) and selects it.
    pub fn add_unit_at(&mut self, kind: UnitKind, center: Pos2) -> UnitId {
        let id = self.fs.add_unit(kind);
        self.layout.insert(id, center);
        self.selection = Some(Selection::Unit(id));
        id
    }

    /// Adds a unit near `hint` (world), offset so it does not cover an existing block.
    pub fn add_unit_near(&mut self, kind: UnitKind, hint: Pos2) -> UnitId {
        let mut center = pos2(
            (hint.x / 10.0).round() * 10.0,
            (hint.y / 10.0).round() * 10.0,
        );
        while self.layout.values().any(|p| p.distance(center) < 30.0) {
            center += Vec2::splat(30.0);
        }
        self.add_unit_at(kind, center)
    }

    /// Whether a new stream between the two ports would be accepted, in either order.
    pub fn can_connect(&self, a: PortRef, b: PortRef) -> bool {
        self.oriented(a, b).is_some_and(|(from, to)| {
            from.unit != to.unit
                && self.fs.stream_at(from).is_none()
                && self.fs.stream_at(to).is_none()
        })
    }

    /// Orders two ports as (outlet, inlet), if they have opposite directions.
    fn oriented(&self, a: PortRef, b: PortRef) -> Option<(PortRef, PortRef)> {
        use flowsheet_core::PortDirection::{Inlet, Outlet};
        match (
            self.fs.port_spec(a)?.direction,
            self.fs.port_spec(b)?.direction,
        ) {
            (Outlet, Inlet) => Some((a, b)),
            (Inlet, Outlet) => Some((b, a)),
            _ => None,
        }
    }

    /// Connects two ports dragged in either direction; reports rejections in
    /// the status line and selects the new stream on success.
    pub fn connect(&mut self, a: PortRef, b: PortRef) -> Option<StreamId> {
        if a == b {
            return None;
        }
        let Some((from, to)) = self.oriented(a, b) else {
            self.status = Some("Connect an outlet to an inlet.".into());
            return None;
        };
        match self.fs.connect(from, to) {
            Ok(id) => {
                self.status = None;
                self.selection = Some(Selection::Stream(id));
                Some(id)
            }
            Err(e) => {
                self.status = Some(format!("Cannot connect: {e}."));
                None
            }
        }
    }

    /// Deletes the selected unit (with its streams and composition) or stream.
    /// A failure attributed to a deleted unit is dropped with it, so the
    /// status bar does not keep reporting a failure nobody can inspect.
    pub fn delete_selection(&mut self) {
        match self.selection.take() {
            Some(Selection::Unit(id)) => {
                self.fs.remove_unit(id);
                self.layout.remove(&id);
                self.feeds.remove(&id);
                if self.failure.as_ref().is_some_and(|f| f.unit == Some(id)) {
                    self.failure = None;
                }
            }
            Some(Selection::Stream(id)) => {
                self.fs.disconnect(id);
            }
            None => {}
        }
    }

    /// World-space bounding box of all blocks.
    pub fn content_bounds(&self) -> Option<Rect> {
        self.fs
            .units()
            .filter_map(|u| Some(block_rect(u.kind, *self.layout.get(&u.id)?)))
            .reduce(|a, b| a.union(b))
    }

    /// Fits the viewport to the content (or resets it when empty).
    pub fn fit_to_content(&mut self, canvas: Rect) {
        match self.content_bounds() {
            Some(bounds) => self.viewport.fit(canvas, bounds.expand(30.0)),
            None => self.viewport = Viewport::default(),
        }
    }

    /// Replaces the project with the qualified A8 demonstration: two natural
    /// gas feeds are mixed, cooled, expanded and separated, with PR-BM and
    /// the explicit synthetic zero-kij approximation (proposal section 4.1).
    pub fn load_demonstration(&mut self) {
        let units = self.units;
        self.load_project(flowsheet_project::demonstration(
            CubicVariantChoice::PengRobinsonBostonMathias,
        ));
        self.units = units;
        self.set_saved(None);
    }

    /// Replaces the editor content with a project. The file path is kept by
    /// the caller ([`Self::open_path`]).
    pub fn load_project(&mut self, project: Project) {
        let auto_solve = self.auto_solve;
        *self = Self {
            auto_solve,
            ..Self::default()
        };
        self.layout = project
            .layout
            .iter()
            .map(|(id, [x, y])| (*id, pos2(*x, *y)))
            .collect();
        // Units without a stored position are placed in a row below.
        let mut spare = 0.0;
        for u in project.flowsheet.units() {
            self.layout.entry(u.id).or_insert_with(|| {
                spare += 120.0;
                pos2(spare, 300.0)
            });
        }
        self.fs = project.flowsheet;
        self.package = Some(project.package);
        self.feeds = project.feeds;
        self.results = project.results;
        self.units = match project.display_units {
            DisplayUnitSet::ProcessMetric => DisplayUnits::ProcessMetric,
            DisplayUnitSet::Si => DisplayUnits::Si,
        };
        self.fit_requested = true;
    }

    /// The current project, or `None` before a package is chosen.
    pub fn to_project(&self) -> Option<Project> {
        Some(Project {
            flowsheet: self.fs.clone(),
            package: self.package.clone()?,
            feeds: self.feeds.clone(),
            layout: self
                .layout
                .iter()
                .filter(|(id, _)| self.fs.unit(**id).is_some())
                .map(|(id, p)| (*id, [p.x, p.y]))
                .collect(),
            display_units: match self.units {
                DisplayUnits::ProcessMetric => DisplayUnitSet::ProcessMetric,
                DisplayUnits::Si => DisplayUnitSet::Si,
            },
            results: self.results.clone(),
        })
    }

    /// Writes the project to `path` and remembers it as the project file.
    pub fn save_to(&mut self, path: &Path) -> Result<(), String> {
        let project = self
            .to_project()
            .ok_or("Create a project with a property package before saving.")?;
        let text = project.to_json().map_err(|e| e.to_string())?;
        std::fs::write(path, &text).map_err(|e| format!("Saving failed: {e}"))?;
        self.path = Some(path.to_path_buf());
        self.set_saved(Some(text));
        self.status = Some(format!("Saved {}.", path.display()));
        Ok(())
    }

    /// Opens a project file. On failure the current project is unchanged.
    pub fn open_path(&mut self, path: &Path) -> Result<(), String> {
        let project = Project::load(path).map_err(|e| format!("Cannot open: {e}."))?;
        self.load_project(project);
        self.path = Some(path.to_path_buf());
        let saved = self.to_project().and_then(|p| p.to_json().ok());
        self.set_saved(saved);
        let stale = self.count_stale();
        self.status = Some(if stale > 0 {
            format!(
                "Opened {}. {stale} stored result(s) are stale; solve to update.",
                path.display()
            )
        } else {
            format!("Opened {}.", path.display())
        });
        Ok(())
    }

    /// Records the serialized project at the last save or open and drops the
    /// cached modified state.
    pub(crate) fn set_saved(&mut self, saved: Option<String>) {
        self.saved = saved;
        self.modified_cache.get_mut().take();
    }

    /// Whether the project changed since it was last saved or opened.
    ///
    /// The answer is cached and recomputed (by serializing the project) only
    /// when the flowsheet revision or another project input changed, so
    /// calling this every frame is cheap.
    pub fn is_modified(&self) -> bool {
        if let Some(cache) = self.modified_cache.borrow().as_ref()
            && cache.matches(self)
        {
            return cache.modified;
        }
        let modified = match (&self.saved, self.to_project()) {
            (_, None) => self.fs.units().next().is_some(),
            (None, Some(_)) => true,
            (Some(saved), Some(p)) => p.to_json().ok().as_ref() != Some(saved),
        };
        *self.modified_cache.borrow_mut() = Some(ModifiedCache {
            revision: self.fs.revision(),
            package: self.package.clone(),
            feeds: self.feeds.clone(),
            layout: self.layout.clone(),
            units: self.units,
            results: self.results.clone(),
            modified,
        });
        modified
    }

    /// Current input fingerprints (none without a package). Cached on the
    /// flowsheet revision, the feed compositions and the package.
    pub fn fingerprints(&mut self) -> Fingerprints {
        if let Some(c) = &self.fingerprint_cache
            && c.revision == self.fs.revision()
            && c.package == self.package
            && c.feeds == self.feeds
        {
            return c.fingerprints.clone();
        }
        let package_fp = match &self.package {
            Some(p) => self.package_cache.get(p).fingerprint,
            None => None,
        };
        let fingerprints = Fingerprints::compute_with(&self.fs, &self.feeds, package_fp);
        self.fingerprint_cache = Some(FingerprintCache {
            revision: self.fs.revision(),
            package: self.package.clone(),
            feeds: self.feeds.clone(),
            fingerprints: fingerprints.clone(),
        });
        fingerprints
    }

    /// Feed composition problems (missing package, missing or invalid fractions).
    pub fn input_issues(&self) -> Vec<InputIssue> {
        match &self.package {
            Some(package) => Project {
                flowsheet: self.fs.clone(),
                package: package.clone(),
                feeds: self.feeds.clone(),
                layout: Default::default(),
                display_units: DisplayUnitSet::default(),
                results: None,
            }
            .input_issues(),
            None => Vec::new(),
        }
    }

    /// Display status of a unit.
    pub fn unit_status(&self, id: UnitId, fps: &Fingerprints, issues: &[InputIssue]) -> ItemStatus {
        if self.fs.unit_status(id) != Some(UnitStatus::Ready)
            || issues.iter().any(|i| i.unit == id)
            || (self.package.is_none()
                && self.fs.unit(id).is_some_and(|u| u.kind == UnitKind::Feed))
        {
            return ItemStatus::Underspecified;
        }
        if let Some(f) = &self.failure
            && f.unit == Some(id)
            && f.unit_fingerprint.is_some()
            && f.unit_fingerprint == fps.units.get(&id).copied()
        {
            return ItemStatus::Failed;
        }
        let Some(results) = &self.results else {
            return ItemStatus::Ready;
        };
        match results.unit_currency(id, fps) {
            Currency::Current if results.units[&id].warnings.is_empty() => ItemStatus::Solved,
            Currency::Current => ItemStatus::Warning,
            Currency::Stale => ItemStatus::Stale,
            Currency::Missing => ItemStatus::Ready,
        }
    }

    /// Display status of a stream result, or `None` if it has none.
    pub fn stream_status(&self, id: StreamId, fps: &Fingerprints) -> Option<ItemStatus> {
        let results = self.results.as_ref()?;
        match results.stream_currency(&self.fs, id, fps) {
            Currency::Current if results.streams[&id].warnings.is_empty() => {
                Some(ItemStatus::Solved)
            }
            Currency::Current => Some(ItemStatus::Warning),
            Currency::Stale => Some(ItemStatus::Stale),
            Currency::Missing => None,
        }
    }

    /// Number of stored unit and stream results that are stale.
    pub fn count_stale(&mut self) -> usize {
        let fps = self.fingerprints();
        let Some(r) = &self.results else {
            return 0;
        };
        r.units
            .keys()
            .filter(|id| self.fs.unit(**id).is_some())
            .filter(|id| r.unit_currency(**id, &fps) == Currency::Stale)
            .count()
            + r.streams
                .keys()
                .filter(|id| r.stream_currency(&self.fs, **id, &fps) == Currency::Stale)
                .count()
    }

    /// Why a solve cannot start now, if it cannot.
    pub fn solve_blocker(&mut self) -> Option<String> {
        let Some(package) = &self.package else {
            return Some(
                "Create a project first: choose components and a package (File > New project)."
                    .into(),
            );
        };
        if let Err(e) = &self.package_cache.get(package).built {
            return Some(format!("The property package cannot be built: {e}"));
        }
        let open = self.fs.recycle_diagnostics().len() + self.input_issues().len();
        (open > 0)
            .then(|| format!("{open} open item(s) must be resolved before solving; see Messages."))
    }

    /// Starts a background solve of the current inputs. Without `repaint`
    /// the caller must poll [`Self::poll_solve`] (tests do this).
    pub fn start_solve(&mut self, repaint: Option<eframe::egui::Context>) {
        if let Some(reason) = self.solve_blocker() {
            self.status = Some(reason);
            return;
        }
        let Some(project) = self.to_project() else {
            return;
        };
        let package = match &self.package_cache.get(&project.package).built {
            Ok(p) => p.clone(),
            Err(_) => return,
        };
        if let Some(job) = self.job.take() {
            job.cancel();
        }
        let inputs = self.fingerprints();
        self.last_attempt = Some(inputs.combined());
        self.job = Some(SolveJob::spawn(project, package, inputs, repaint));
        self.status = None;
    }

    /// Cancels a running solve. The result, if it still arrives, is discarded
    /// and earlier results keep their currency.
    ///
    /// The cancelled inputs stay recorded as the last attempt on purpose:
    /// clearing it would make auto-solve restart the very calculation the
    /// user just stopped on the next frame. Auto-solve resumes with the next
    /// input change; F5 (or Solve) retries the unchanged inputs at once.
    pub fn cancel_solve(&mut self) {
        if let Some(job) = self.job.take() {
            job.cancel();
            self.status = Some(
                "Calculation cancelled; no results were changed. Press F5 to solve again.".into(),
            );
        }
    }

    /// Whether a background solve is running.
    pub fn is_solving(&self) -> bool {
        self.job.is_some()
    }

    /// `(completed, total, recycle_iteration)` of the running solve.
    pub fn solve_progress(&self) -> Option<(usize, usize, usize)> {
        self.job.as_ref().map(SolveJob::progress)
    }

    /// Collects a finished background solve. Returns `true` if one finished.
    pub fn poll_solve(&mut self) -> bool {
        let Some(outcome) = self.job.as_ref().and_then(SolveJob::try_take) else {
            return false;
        };
        let job = self.job.take().expect("job present");
        match outcome {
            SolveOutcome::Solved(results) => {
                self.results = Some(results);
                self.failure = None;
                self.status = None;
            }
            SolveOutcome::Failed(f) if f.cancelled => {
                self.status = Some(f.message);
            }
            SolveOutcome::Failed(f) => {
                self.status = None;
                self.failure = Some(Failure {
                    unit: f.unit,
                    unit_fingerprint: f.unit.and_then(|u| job.inputs.units.get(&u).copied()),
                    message: f.message,
                    inputs: job.inputs.combined(),
                });
            }
        }
        true
    }

    /// Blocks until a running solve finishes (tests and the self-test).
    pub fn wait_for_solve(&mut self) {
        while self.is_solving() {
            if !self.poll_solve() {
                std::thread::sleep(std::time::Duration::from_millis(2));
            }
        }
    }

    /// Starts an automatic solve if enabled, idle, complete and changed
    /// since the last attempt.
    pub fn auto_solve_if_due(&mut self, repaint: Option<eframe::egui::Context>) {
        if !self.auto_solve || self.is_solving() {
            return;
        }
        let combined = self.fingerprints().combined();
        if self.last_attempt == Some(combined) || self.solve_blocker().is_some() {
            return;
        }
        self.start_solve(repaint);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn connect_accepts_either_drag_direction() {
        let mut ed = Editor::default();
        let feed = ed.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let valve = ed.add_unit_at(UnitKind::Valve, pos2(200.0, 0.0));
        // Dragged from the inlet back to the outlet.
        let id = ed
            .connect(PortRef::new(valve, 0), PortRef::new(feed, 0))
            .expect("connected");
        let s = ed.fs.stream(id).unwrap();
        assert_eq!(s.from, PortRef::new(feed, 0));
        assert_eq!(s.to, PortRef::new(valve, 0));
        assert_eq!(ed.selection, Some(Selection::Stream(id)));
    }

    #[test]
    fn rejected_connection_sets_status() {
        let mut ed = Editor::default();
        let a = ed.add_unit_at(UnitKind::Feed, pos2(-100.0, 0.0));
        let b = ed.add_unit_at(UnitKind::Feed, pos2(200.0, 0.0));
        assert!(!ed.can_connect(PortRef::new(a, 0), PortRef::new(b, 0)));
        assert!(ed.connect(PortRef::new(a, 0), PortRef::new(b, 0)).is_none());
        assert!(ed.status.is_some());
        assert_eq!(ed.fs.streams().count(), 0);
    }

    #[test]
    fn deleting_a_unit_removes_layout_streams_and_composition() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        let feed = ed.fs.units().find(|u| u.kind == UnitKind::Feed).unwrap().id;
        let mixer = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Mixer)
            .unwrap()
            .id;
        ed.selection = Some(Selection::Unit(feed));
        ed.delete_selection();
        assert!(ed.fs.unit(feed).is_none());
        assert!(!ed.layout.contains_key(&feed));
        assert!(!ed.feeds.contains_key(&feed));
        assert_eq!(ed.fs.streams().count(), 6);
        assert_eq!(
            ed.fs.unit_status(mixer),
            Some(flowsheet_core::UnitStatus::Underspecified)
        );
    }

    #[test]
    fn add_unit_near_avoids_stacking() {
        let mut ed = Editor::default();
        let a = ed.add_unit_near(UnitKind::Heater, pos2(100.0, 100.0));
        let b = ed.add_unit_near(UnitKind::Heater, pos2(100.0, 100.0));
        assert_ne!(ed.layout[&a], ed.layout[&b]);
    }

    #[test]
    fn demonstration_is_complete_and_ready_to_solve() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        assert_eq!(ed.fs.units().count(), 8);
        assert_eq!(ed.fs.streams().count(), 7);
        assert!(ed.fs.recycle_diagnostics().is_empty());
        assert!(ed.input_issues().is_empty());
        assert_eq!(ed.solve_blocker(), None);
        let fps = ed.fingerprints();
        assert!(
            ed.fs
                .units()
                .all(|u| ed.unit_status(u.id, &fps, &[]) == ItemStatus::Ready)
        );
        assert!(ed.content_bounds().is_some());
        assert!(ed.is_modified());
    }

    #[test]
    fn solving_marks_results_current_and_edits_mark_them_stale() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        ed.start_solve(None);
        ed.wait_for_solve();
        assert!(ed.failure.is_none(), "{:?}", ed.failure);
        let fps = ed.fingerprints();
        for u in ed.fs.units() {
            assert_eq!(
                ed.unit_status(u.id, &fps, &[]),
                ItemStatus::Solved,
                "{}",
                u.name
            );
        }
        let valve = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Valve)
            .unwrap()
            .id;
        ed.fs.set_spec(
            valve,
            flowsheet_core::UnitSpec::Valve {
                pressure: flowsheet_core::PressureSpec::OutletPressure(Some(6.0e6)),
            },
        );
        let fps = ed.fingerprints();
        assert_eq!(ed.unit_status(valve, &fps, &[]), ItemStatus::Stale);
        ed.start_solve(None);
        ed.wait_for_solve();
        let fps = ed.fingerprints();
        assert_eq!(ed.unit_status(valve, &fps, &[]), ItemStatus::Failed);
        let failure = ed.failure.as_ref().unwrap();
        assert!(failure.message.contains("V-100"));
        let mixer = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Mixer)
            .unwrap()
            .id;
        assert_eq!(ed.unit_status(mixer, &fps, &[]), ItemStatus::Solved);
        // The separator keeps its old, now stale, values; none are current.
        let sep = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Separator)
            .unwrap()
            .id;
        assert_eq!(ed.unit_status(sep, &fps, &[]), ItemStatus::Stale);
    }

    #[test]
    fn auto_solve_runs_once_per_input_change() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        ed.auto_solve = true;
        ed.auto_solve_if_due(None);
        assert!(ed.is_solving());
        ed.wait_for_solve();
        ed.auto_solve_if_due(None);
        assert!(!ed.is_solving(), "unchanged inputs must not solve again");
        let feed = ed.fs.units().find(|u| u.kind == UnitKind::Feed).unwrap().id;
        ed.feeds.get_mut(&feed).unwrap().swap(0, 1);
        ed.auto_solve_if_due(None);
        assert!(ed.is_solving());
        ed.cancel_solve();
        assert!(!ed.is_solving());
        assert!(ed.status.as_deref().unwrap().contains("cancelled"));
    }

    #[test]
    fn cancelled_auto_solve_is_not_restarted_until_inputs_change() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        ed.auto_solve = true;
        ed.auto_solve_if_due(None);
        assert!(ed.is_solving());
        ed.cancel_solve();
        ed.auto_solve_if_due(None);
        assert!(!ed.is_solving(), "cancel must not restart at once");
        assert!(ed.status.as_deref().unwrap().contains("F5"));
        let valve = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Valve)
            .unwrap()
            .id;
        ed.fs.set_spec(
            valve,
            flowsheet_core::UnitSpec::Valve {
                pressure: flowsheet_core::PressureSpec::OutletPressure(Some(2.0e6)),
            },
        );
        ed.auto_solve_if_due(None);
        assert!(ed.is_solving(), "an input change solves again");
        ed.cancel_solve();
    }

    #[test]
    fn deleting_the_failed_unit_clears_the_failure() {
        let mut ed = Editor::default();
        ed.load_demonstration();
        let valve = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Valve)
            .unwrap()
            .id;
        ed.fs.set_spec(
            valve,
            flowsheet_core::UnitSpec::Valve {
                pressure: flowsheet_core::PressureSpec::OutletPressure(Some(6.0e6)),
            },
        );
        ed.start_solve(None);
        ed.wait_for_solve();
        assert_eq!(ed.failure.as_ref().and_then(|f| f.unit), Some(valve));
        // Deleting another unit keeps the failure.
        let product = ed
            .fs
            .units()
            .find(|u| u.kind == UnitKind::Product)
            .unwrap()
            .id;
        ed.selection = Some(Selection::Unit(product));
        ed.delete_selection();
        assert!(ed.failure.is_some());
        ed.selection = Some(Selection::Unit(valve));
        ed.delete_selection();
        assert!(ed.failure.is_none());
    }

    #[test]
    fn modified_state_and_fingerprints_follow_every_kind_of_edit() {
        let dir =
            std::env::temp_dir().join(format!("thermo-desktop-modified-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("m.flowsheet.json");
        let mut ed = Editor::default();
        ed.load_demonstration();
        assert!(ed.is_modified());
        ed.save_to(&path).unwrap();
        assert!(!ed.is_modified());
        assert!(!ed.is_modified(), "cached answer is stable");

        let feed = ed.fs.units().find(|u| u.kind == UnitKind::Feed).unwrap().id;
        // Flowsheet edit (advances the revision).
        ed.fs.rename_unit(feed, "F-1");
        assert!(ed.is_modified());
        ed.save_to(&path).unwrap();
        assert!(!ed.is_modified());
        // Layout edit (no revision change).
        *ed.layout.get_mut(&feed).unwrap() += Vec2::new(10.0, 0.0);
        assert!(ed.is_modified());
        ed.save_to(&path).unwrap();
        // Feed composition edit (no revision change) also changes fingerprints.
        let before = ed.fingerprints();
        assert_eq!(ed.fingerprints(), before);
        ed.feeds.get_mut(&feed).unwrap().swap(0, 1);
        assert!(ed.is_modified());
        assert_ne!(ed.fingerprints(), before);
        ed.save_to(&path).unwrap();
        // Display units and results.
        ed.units = DisplayUnits::Si;
        assert!(ed.is_modified());
        ed.save_to(&path).unwrap();
        ed.start_solve(None);
        ed.wait_for_solve();
        assert!(ed.is_modified());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn save_and_open_round_trip_through_a_file() {
        let dir =
            std::env::temp_dir().join(format!("thermo-desktop-editor-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("demo.flowsheet.json");
        let mut ed = Editor::default();
        ed.load_demonstration();
        ed.start_solve(None);
        ed.wait_for_solve();
        ed.units = DisplayUnits::Si;
        ed.save_to(&path).unwrap();
        assert!(!ed.is_modified());

        let mut other = Editor::default();
        other.open_path(&path).unwrap();
        assert!(!other.is_modified());
        assert_eq!(other.units, DisplayUnits::Si);
        assert_eq!(other.count_stale(), 0);
        let fps = other.fingerprints();
        assert!(
            other
                .fs
                .units()
                .all(|u| other.unit_status(u.id, &fps, &[]) == ItemStatus::Solved)
        );
        assert_eq!(other.layout, ed.layout);

        // A file that cannot be read leaves the open project untouched.
        std::fs::write(dir.join("bad.json"), "{\"format_version\": 99}").unwrap();
        let err = other.open_path(&dir.join("bad.json")).unwrap_err();
        assert!(err.contains("version 99"), "{err}");
        assert_eq!(other.path.as_deref(), Some(path.as_path()));
        assert_eq!(other.fs.units().count(), 8);
    }
}
