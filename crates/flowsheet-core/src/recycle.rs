//! Explicit tear streams and bounded direct substitution for cyclic networks.
use crate::engine::molar_thermal_scale;
use crate::{
    Calculation, Diagnostic, DiagnosticKind, FeedCompositions, Flowsheet, MaterialStream,
    PortDirection, PortRef, SolveError, SolveProgress, StreamId, UnitKind, UnitSpec,
};
use property_packages::{CubicPackage, EquilibriumState};
use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;

/// Bounded direct substitution. There is no convergence claim for arbitrary loops.
#[derive(Debug, Clone, Copy)]
pub struct RecycleOptions {
    pub maximum_iterations: usize,
    /// Maximum scaled mismatch in component flows, energy flow, T, p, h and x.
    /// Denominators are max(|old|, |new|, 1 SI unit); energy flow and h also
    /// use the larger thermal scale n·cp·T (cp·T) as a floor, so the target
    /// does not depend on the caloric reference point. x uses absolute error.
    pub tolerance: f64,
}

impl Default for RecycleOptions {
    fn default() -> Self {
        Self {
            maximum_iterations: 500,
            tolerance: 1e-10,
        }
    }
}

/// Closure evidence; final unit balances use the actual connected stream values.
#[derive(Debug, Clone)]
pub struct RecycleConvergence {
    pub iterations: usize,
    pub max_scaled_residual: f64,
    pub tear_residuals: BTreeMap<StreamId, f64>,
}

fn error(message: impl Into<String>) -> SolveError {
    SolveError {
        unit: None,
        message: message.into(),
        diagnostics: vec![],
        cancelled: false,
    }
}

impl Flowsheet {
    /// Recycle-capable solve using deterministic automatic tears and defaults.
    pub fn solve_network(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
    ) -> Result<Calculation, SolveError> {
        self.solve_network_with_progress(package, feeds, |_| ControlFlow::Continue(()))
    }
    /// Structural diagnostics for the recycle-capable solver. A loop itself is
    /// valid input; missing ports/specifications remain errors.
    pub fn recycle_diagnostics(&self) -> Vec<Diagnostic> {
        self.diagnostics()
            .into_iter()
            .filter(|d| !matches!(d.kind, DiagnosticKind::RecycleUnsupported(_)))
            .collect()
    }

    /// Deterministically cuts one stream of each remaining directed cycle.
    /// All other branches remain intact. This is not a minimum-tear algorithm.
    pub fn automatic_tears(&self) -> BTreeSet<StreamId> {
        let mut graph = self.clone();
        let mut tears = BTreeSet::new();
        while let Some(cycle) = graph.find_cycle() {
            // Prefer a mixer inlet: a fresh flowing inlet can start a recycle
            // from zero flow without sending an empty stream into a heater.
            let id = graph
                .streams()
                .filter(|s| {
                    cycle.iter().enumerate().any(|(i, from)| {
                        *from == s.from.unit && cycle[(i + 1) % cycle.len()] == s.to.unit
                    })
                })
                .min_by_key(|s| {
                    (
                        graph.unit(s.to.unit).expect("unit").kind != UnitKind::Mixer,
                        s.id,
                    )
                })
                .expect("cycle edge")
                .id;
            graph.disconnect(id);
            tears.insert(id);
        }
        tears
    }

    /// Solves acyclic networks normally, or automatically tears loops and starts
    /// mixer-inlet tears at zero flow with the first feed's flashed state; other
    /// tears start with the flow and flashed state of the nearest upstream feed.
    /// Each call starts afresh. For process-specific starting states use
    /// [`Self::solve_recycles`].
    pub fn solve_network_with_progress(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
        mut progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Result<Calculation, SolveError> {
        let ids = self.automatic_tears();
        if ids.is_empty() {
            return self.solve_with_progress(package, feeds, progress);
        }
        self.tear_order(&ids)?;
        self.validate_inputs(package, feeds)?;
        let feed = self
            .units()
            .find(|u| u.kind == UnitKind::Feed)
            .expect("preflight feed");
        if progress(SolveProgress {
            iteration: 0,
            completed: 0,
            total: self.units().count(),
            unit: feed.id,
        })
        .is_break()
        {
            let mut e = error("Calculation cancelled; no results were produced.");
            e.cancelled = true;
            return Err(e);
        }
        let mut flashed = BTreeMap::new();
        let mut feed_state = |id: crate::UnitId| -> Result<(f64, EquilibriumState), SolveError> {
            if let Some(found) = flashed.get(&id) {
                return Ok(Clone::clone(found));
            }
            let UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(p),
                molar_flow_mol_s: Some(n),
            } = self.unit(id).expect("feed unit").spec
            else {
                unreachable!("preflight specification")
            };
            let z: Vec<_> = feeds[&id].entries().iter().map(|(_, x)| *x).collect();
            let state = package
                .tp(&z, t, p)
                .map_err(|e| error(format!("Recycle initialization: {e}")))?;
            flashed.insert(id, (n, state.clone()));
            Ok((n, state))
        };
        let (_, first) = feed_state(feed.id)?;
        let mut guesses = BTreeMap::new();
        for id in ids {
            let stream = self.stream(id).expect("tear stream");
            // A mixer tear starts empty: its fresh inlet carries the first pass.
            // Other units would only pass an empty tear through, so they start
            // with the flow and flashed state of the nearest upstream feed.
            let upstream = (self.unit(stream.to.unit).expect("unit").kind != UnitKind::Mixer)
                .then(|| self.upstream_feed(stream.from.unit))
                .flatten();
            let guess = match upstream {
                Some(feed) => {
                    let (n, state) = feed_state(feed)?;
                    MaterialStream {
                        molar_flow_mol_s: n,
                        state: Some(state),
                    }
                }
                None => MaterialStream {
                    molar_flow_mol_s: 0.0,
                    state: Some(first.clone()),
                },
            };
            guesses.insert(id, guess);
        }
        self.solve_recycles(
            package,
            feeds,
            &guesses,
            RecycleOptions::default(),
            progress,
        )
    }

    /// Explicit tear states are generated by the caller's *same* package. They
    /// are re-flashed from composition, pressure and enthalpy before use, so
    /// cached phase properties and provenance cannot leak into the calculation.
    /// All cyclic dependencies must be cut. Each pass reads the previous tear
    /// values even when its source unit happens to run first (Jacobi iteration).
    /// Cancellation is checked before every unit, on every pass. Failure or
    /// iteration exhaustion returns no partially current calculation.
    pub fn solve_recycles(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
        guesses: &BTreeMap<StreamId, MaterialStream>,
        options: RecycleOptions,
        mut progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Result<Calculation, SolveError> {
        if options.maximum_iterations == 0
            || !options.tolerance.is_finite()
            || options.tolerance <= 0.0
            || options.tolerance > 1e-8
        {
            return Err(error(
                "Recycle options require positive iterations and a finite tolerance in (0, 1e-8].",
            ));
        }
        let ids = guesses.keys().copied().collect();
        let order = self.tear_order(&ids)?;
        if guesses.is_empty() {
            return self.solve_with_progress(package, feeds, progress);
        }
        self.validate_inputs(package, feeds)?;
        let mut tears = BTreeMap::new();
        for (&id, guess) in guesses {
            if progress(SolveProgress {
                iteration: 0,
                completed: 0,
                total: order.len(),
                unit: self.stream(id).expect("tear stream").to.unit,
            })
            .is_break()
            {
                let mut e = error("Calculation cancelled; no results were produced.");
                e.cancelled = true;
                return Err(e);
            }
            if !guess.molar_flow_mol_s.is_finite() || guess.molar_flow_mol_s < 0.0 {
                return Err(error(format!(
                    "{id}: tear flow must be finite and nonnegative"
                )));
            }
            let state = match &guess.state {
                Some(s) => Some(
                    package
                        .ph(&s.mole_fractions, s.pressure_pa, s.enthalpy())
                        .map_err(|e| error(format!("{id}: invalid tear starting state: {e}")))?,
                ),
                None if guess.molar_flow_mol_s == 0.0 => None,
                None => return Err(error(format!("{id}: positive tear flow requires a state"))),
            };
            tears.insert(
                id,
                MaterialStream {
                    molar_flow_mol_s: guess.molar_flow_mol_s,
                    state,
                },
            );
        }
        let mut last = f64::INFINITY;
        let mut balance_failure = None;
        for iteration in 1..=options.maximum_iterations {
            let mut report = self.solve_pass(package, feeds, order.clone(), &tears, |mut p| {
                p.iteration = iteration;
                progress(p)
            })?;
            let residuals: BTreeMap<_, _> = tears
                .iter()
                .map(|(&id, old)| {
                    (
                        id,
                        mismatch(old, &report.streams[&id], package.model().fractions().len()),
                    )
                })
                .collect();
            last = residuals.values().copied().fold(0.0_f64, f64::max);
            if last <= options.tolerance {
                match self.close_balances(&mut report, package) {
                    Ok(()) => {
                        report.recycle = Some(RecycleConvergence {
                            iterations: iteration,
                            max_scaled_residual: last,
                            tear_residuals: residuals,
                        });
                        return Ok(report);
                    }
                    Err(e) => balance_failure = Some(e),
                }
            }
            for (&id, value) in &mut tears {
                *value = report.streams[&id].clone();
            }
        }
        Err(error(format!(
            "Recycle did not converge after {} iterations; maximum scaled tear residual {last:.6e} (target {:.6e}); tears: {:?}. {}No results were produced.",
            options.maximum_iterations,
            options.tolerance,
            guesses.keys().collect::<Vec<_>>(),
            balance_failure.map_or_else(String::new, |e| format!(
                "Final connected-stream balance: {e}. "
            ))
        )))
    }

    /// Nearest feed upstream of `unit` (breadth-first against the stream
    /// direction, ties broken by stream identifier), if any.
    fn upstream_feed(&self, unit: crate::UnitId) -> Option<crate::UnitId> {
        let mut seen = BTreeSet::from([unit]);
        let mut queue = std::collections::VecDeque::from([unit]);
        while let Some(id) = queue.pop_front() {
            if self.unit(id)?.kind == UnitKind::Feed {
                return Some(id);
            }
            for s in self.streams().filter(|s| s.to.unit == id) {
                if seen.insert(s.from.unit) {
                    queue.push_back(s.from.unit);
                }
            }
        }
        None
    }

    fn tear_order(&self, tears: &BTreeSet<StreamId>) -> Result<Vec<crate::UnitId>, SolveError> {
        let diagnostics = self.recycle_diagnostics();
        if !diagnostics.is_empty() {
            return Err(SolveError {
                message: diagnostics
                    .iter()
                    .map(|d| d.message.as_str())
                    .collect::<Vec<_>>()
                    .join("\n"),
                diagnostics,
                ..error("")
            });
        }
        if let Some(id) = tears.iter().find(|id| self.stream(**id).is_none()) {
            return Err(error(format!("Unknown tear stream {id}")));
        }
        let mut incoming: BTreeMap<_, usize> = self.units().map(|u| (u.id, 0)).collect();
        for s in self.streams().filter(|s| !tears.contains(&s.id)) {
            *incoming.get_mut(&s.to.unit).expect("unit") += 1;
        }
        let mut ready: BTreeSet<_> = incoming
            .iter()
            .filter_map(|(&id, &n)| (n == 0).then_some(id))
            .collect();
        let mut order = vec![];
        while let Some(id) = ready.pop_first() {
            order.push(id);
            for s in self
                .streams()
                .filter(|s| s.from.unit == id && !tears.contains(&s.id))
            {
                let n = incoming.get_mut(&s.to.unit).expect("unit");
                *n -= 1;
                if *n == 0 {
                    ready.insert(s.to.unit);
                }
            }
        }
        if order.len() != incoming.len() {
            return Err(error("Tear streams do not break every recycle loop."));
        }
        Ok(order)
    }

    fn close_balances(
        &self,
        report: &mut Calculation,
        package: &CubicPackage,
    ) -> Result<(), String> {
        for u in self
            .units()
            .filter(|u| !matches!(u.kind, UnitKind::Feed | UnitKind::Product))
        {
            let mut inputs = vec![];
            let mut outputs = vec![];
            for (i, port) in u.kind.ports().iter().enumerate() {
                let s = self.stream_at(PortRef::new(u.id, i)).expect("checked port");
                let values = if port.direction == PortDirection::Inlet {
                    &mut inputs
                } else {
                    &mut outputs
                };
                values.push(report.streams[&s.id].clone());
            }
            let unit = report.units.get_mut(&u.id).expect("calculated unit");
            unit.balance = Some(
                crate::engine::balance(
                    &inputs,
                    &outputs,
                    unit.duty_w,
                    package.model().fractions().len(),
                )
                .map_err(|e| format!("{}: {e}", u.name))?,
            );
        }
        Ok(())
    }
}

fn mismatch(a: &MaterialStream, b: &MaterialStream, count: usize) -> f64 {
    fn scaled(a: f64, b: f64) -> f64 {
        scaled_by(a, b, 0.0)
    }
    /// Like `scaled`, with an additional reference-independent floor.
    fn scaled_by(a: f64, b: f64, floor: f64) -> f64 {
        (a - b).abs() / a.abs().max(b.abs()).max(floor).max(1.0)
    }
    // Enthalpy depends on the caloric reference point, so |h| can be near zero
    // for an ordinary stream. The thermal scale n cp T (cp T per mole) keeps
    // the enthalpy targets physically relative; flows, T, p and x are unchanged.
    let mut r = scaled(a.molar_flow_mol_s, b.molar_flow_mol_s).max(scaled_by(
        a.enthalpy_flow_w(),
        b.enthalpy_flow_w(),
        a.thermal_scale_w().max(b.thermal_scale_w()),
    ));
    for (x, y) in a
        .component_flows_mol_s(count)
        .into_iter()
        .zip(b.component_flows_mol_s(count))
    {
        r = r.max(scaled(x, y));
    }
    if a.molar_flow_mol_s == 0.0 && b.molar_flow_mol_s == 0.0 {
        return r;
    }
    match (&a.state, &b.state) {
        (Some(a), Some(b)) => {
            r = r
                .max(scaled(a.temperature_k, b.temperature_k))
                .max(scaled(a.pressure_pa, b.pressure_pa))
                .max(scaled_by(
                    a.enthalpy(),
                    b.enthalpy(),
                    molar_thermal_scale(a).max(molar_thermal_scale(b)),
                ));
            for (x, y) in a.mole_fractions.iter().zip(&b.mole_fractions) {
                r = r.max((x - y).abs());
            }
            r
        }
        _ => f64::INFINITY,
    }
}
