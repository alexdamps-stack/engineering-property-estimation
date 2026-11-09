//! Background calculation with progress and cancellation, and the cached
//! property package of the current project.

use std::ops::ControlFlow;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, TryRecvError, channel};

use eframe::egui::Context;
use flowsheet_project::{
    Fingerprint, Fingerprints, PackageSelection, PreparedPackage, Project, Provenance, ResultSet,
    SolveFailure, solve_project,
};

/// The built package for one selection, rebuilt only when the selection changes.
#[derive(Default)]
pub(crate) struct PackageCache {
    entry: Option<(PackageSelection, CachedPackage)>,
}

pub(crate) struct CachedPackage {
    pub built: Result<PreparedPackage, String>,
    pub fingerprint: Option<Fingerprint>,
    pub provenance: Option<Provenance>,
}

impl PackageCache {
    pub(crate) fn get(&mut self, selection: &PackageSelection) -> &CachedPackage {
        if self.entry.as_ref().is_none_or(|(s, _)| s != selection) {
            // One model build per selection; provenance and fingerprint
            // come from the prepared package.
            let built = selection.prepare().map_err(|e| e.to_string());
            let cached = CachedPackage {
                fingerprint: built.as_ref().ok().map(PreparedPackage::fingerprint),
                provenance: built.as_ref().ok().map(|p| p.provenance().clone()),
                built,
            };
            self.entry = Some((selection.clone(), cached));
        }
        &self.entry.as_ref().expect("filled above").1
    }
}

pub(crate) enum SolveOutcome {
    Solved(ResultSet),
    Failed(SolveFailure),
}

/// A solve running on a worker thread. Dropping or cancelling the job
/// discards its outcome; the worker stops at the next unit boundary.
pub(crate) struct SolveJob {
    receiver: Receiver<SolveOutcome>,
    cancel: Arc<AtomicBool>,
    completed: Arc<AtomicUsize>,
    total: Arc<AtomicUsize>,
    iteration: Arc<AtomicUsize>,
    /// Fingerprints of the solved inputs.
    pub inputs: Fingerprints,
}

impl SolveJob {
    pub(crate) fn spawn(
        project: Project,
        package: PreparedPackage,
        inputs: Fingerprints,
        repaint: Option<Context>,
    ) -> Self {
        let (sender, receiver) = channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let completed = Arc::new(AtomicUsize::new(0));
        let total = Arc::new(AtomicUsize::new(project.flowsheet.units().count()));
        let iteration = Arc::new(AtomicUsize::new(0));
        let pass = iteration.clone();
        let (c, done, all) = (cancel.clone(), completed.clone(), total.clone());
        std::thread::Builder::new()
            .name("flowsheet-solve".into())
            .spawn(move || {
                let progress = |p: flowsheet_core::SolveProgress| {
                    pass.store(p.iteration, Ordering::Relaxed);
                    done.store(p.completed, Ordering::Relaxed);
                    all.store(p.total, Ordering::Relaxed);
                    if let Some(ctx) = &repaint {
                        ctx.request_repaint();
                    }
                    if c.load(Ordering::Relaxed) {
                        ControlFlow::Break(())
                    } else {
                        ControlFlow::Continue(())
                    }
                };
                let now = std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .map_or(0, |d| d.as_secs());
                let outcome = match solve_project(&project, &package, progress, now) {
                    Ok(r) => SolveOutcome::Solved(r),
                    Err(f) => SolveOutcome::Failed(f),
                };
                // The receiver is gone if the job was cancelled; that is fine.
                let _ = sender.send(outcome);
                if let Some(ctx) = &repaint {
                    ctx.request_repaint();
                }
            })
            .expect("spawn solver thread");
        Self {
            receiver,
            cancel,
            completed,
            total,
            iteration,
            inputs,
        }
    }

    pub(crate) fn cancel(self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub(crate) fn progress(&self) -> (usize, usize, usize) {
        (
            self.completed.load(Ordering::Relaxed),
            self.total.load(Ordering::Relaxed),
            self.iteration.load(Ordering::Relaxed),
        )
    }

    pub(crate) fn try_take(&self) -> Option<SolveOutcome> {
        match self.receiver.try_recv() {
            Ok(outcome) => Some(outcome),
            Err(TryRecvError::Empty) => None,
            Err(TryRecvError::Disconnected) => Some(SolveOutcome::Failed(SolveFailure {
                unit: None,
                message:
                    "The calculation stopped unexpectedly (worker thread ended without a result)."
                        .into(),
                cancelled: false,
            })),
        }
    }
}

impl Drop for SolveJob {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
    }
}
