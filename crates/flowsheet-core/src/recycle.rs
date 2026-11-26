//! Explicit tear streams and bounded direct substitution for cyclic networks.
use crate::engine::{empty, molar_thermal_scale, validate_empty_thermal_spec};
use crate::{
    Calculation, Diagnostic, DiagnosticKind, ExchangerSpec, FeedCompositions, Flowsheet,
    MaterialStream, PortMedium, PortRef, SeparatorThermal, SolveError, SolveProgress, StreamId,
    ThermalSpec, UnitKind, UnitSpec,
};
use property_packages::FlowsheetPackage;
use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};
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
        package: &dyn FlowsheetPackage,
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
            // Prefer a material stream (energy streams cannot be torn) and
            // a mixer inlet: a fresh flowing inlet can start a recycle from
            // zero flow without sending an empty stream into a heater.
            let id = graph
                .streams()
                .filter(|s| {
                    cycle.iter().enumerate().any(|(i, from)| {
                        *from == s.from.unit && cycle[(i + 1) % cycle.len()] == s.to.unit
                    })
                })
                .min_by_key(|s| {
                    (
                        graph.stream_medium(s.id) != Some(PortMedium::Material),
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
    /// mixer-inlet tears at zero flow. Other tears start from the nearest feed
    /// along the same material path and use the source outlet temperature
    /// when specified. For disconnected loops or process-specific starts use
    /// [`Self::solve_recycles`].
    pub fn solve_network_with_progress(
        &self,
        package: &dyn FlowsheetPackage,
        feeds: &FeedCompositions,
        mut progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Result<Calculation, SolveError> {
        let ids = self.automatic_tears();
        if ids.is_empty() {
            return self.solve_with_progress(package, feeds, progress);
        }
        self.tear_order(&ids)?;
        self.validate_inputs(package, feeds)?;
        let mut guesses = BTreeMap::new();
        for id in ids {
            let stream = self.stream(id).expect("tear stream");
            if self.unit(stream.to.unit).expect("tear inlet unit").kind == UnitKind::Mixer {
                guesses.insert(id, empty());
                continue;
            }
            if progress(SolveProgress {
                iteration: 0,
                completed: 0,
                total: self.units().count(),
                unit: stream.to.unit,
            })
            .is_break()
            {
                return Err(SolveError::cancelled());
            }
            let feed = self.upstream_feed(stream.from).ok_or_else(|| {
                error(format!("{id}: recycle initialization requires a materially reachable feed; supply explicit tear states"))
            })?;
            let UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(p),
                molar_flow_mol_s: Some(n),
            } = self.unit(feed).expect("feed unit").spec
            else {
                unreachable!("preflight specification")
            };
            // A specified equipment outlet can be much warmer than the fresh feed.
            let source = self.unit(stream.from.unit).expect("tear outlet unit");
            let t = match source.spec {
                UnitSpec::HeatExchange {
                    thermal: ThermalSpec::OutletTemperature(Some(t)),
                    ..
                }
                | UnitSpec::Separator {
                    thermal: SeparatorThermal::OutletTemperature(Some(t)),
                    ..
                }
                | UnitSpec::ConversionReactor {
                    thermal: SeparatorThermal::OutletTemperature(Some(t)),
                    ..
                }
                | UnitSpec::YieldReactor {
                    thermal: SeparatorThermal::OutletTemperature(Some(t)),
                    ..
                } => t,
                UnitSpec::HeatExchanger {
                    thermal: ExchangerSpec::HotOutletTemperature(Some(t)),
                    ..
                } if stream.from.port == 2 => t,
                UnitSpec::HeatExchanger {
                    thermal: ExchangerSpec::ColdOutletTemperature(Some(t)),
                    ..
                } if stream.from.port == 3 => t,
                _ => t,
            };
            let z: Vec<_> = feeds[&feed].entries().iter().map(|(_, x)| *x).collect();
            let state = package
                .tp(&z, t, p)
                .map_err(|e| error(format!("{id}: recycle initialization: {e}")))?;
            guesses.insert(
                id,
                MaterialStream {
                    molar_flow_mol_s: n,
                    state: Some(state),
                },
            );
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
        package: &dyn FlowsheetPackage,
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
        if let Some(id) = guesses
            .keys()
            .find(|id| self.stream_medium(**id) == Some(PortMedium::Energy))
        {
            return Err(error(format!(
                "{id} is an energy stream; tear a material stream of the loop instead."
            )));
        }
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
                return Err(SolveError::cancelled());
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
                        mismatch(old, &report.streams[&id], package.component_ids().len()),
                    )
                })
                .collect();
            last = residuals.values().copied().fold(0.0_f64, f64::max);
            if last <= options.tolerance {
                match self.close_balances(&mut report, package)? {
                    None => {
                        report.recycle = Some(RecycleConvergence {
                            iterations: iteration,
                            max_scaled_residual: last,
                            tear_residuals: residuals,
                        });
                        return Ok(report);
                    }
                    Some(e) => balance_failure = Some(e),
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

    /// Nearest feed upstream of an outlet along material connections.
    fn upstream_feed(&self, outlet: PortRef) -> Option<crate::UnitId> {
        let mut seen = HashSet::from([outlet]);
        let mut queue = VecDeque::from([outlet]);
        while let Some(outlet) = queue.pop_front() {
            let unit = self.unit(outlet.unit).expect("material path unit");
            if unit.kind == UnitKind::Feed {
                return Some(unit.id);
            }
            for s in self.material_streams().filter(|s| {
                s.to.unit == unit.id
                    && (unit.kind != UnitKind::HeatExchanger || s.to.port + 2 == outlet.port)
            }) {
                if seen.insert(s.from) {
                    queue.push_back(s.from);
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
        let order = self.kahn_order(|s| !tears.contains(&s.id));
        if order.len() != self.units().count() {
            return Err(error("Tear streams do not break every recycle loop."));
        }
        Ok(order)
    }

    // A balance mismatch can improve on the next pass. An infeasible duty at
    // the converged inlet state must be reported to the caller immediately.
    fn close_balances(
        &self,
        report: &mut Calculation,
        package: &dyn FlowsheetPackage,
    ) -> Result<Option<String>, SolveError> {
        for u in self.units().filter(|u| !u.kind.is_boundary()) {
            let ports = self.unit_ports(u);
            let stream = |s: &Option<StreamId>| {
                s.map_or_else(crate::engine::empty, |s| report.streams[&s].clone())
            };
            // An unconnected optional outlet is not part of the network and
            // its value is not kept, so the unit keeps its pass balance.
            if ports.material_outlets.iter().any(Option::is_none) {
                continue;
            }
            let inputs: Vec<_> = ports.material_inlets.iter().map(stream).collect();
            let outputs: Vec<_> = ports.material_outlets.iter().map(stream).collect();
            let unit = report.units.get_mut(&u.id).expect("calculated unit");
            let terms = unit
                .balance
                .as_ref()
                .and_then(|b| b.reaction.as_ref())
                .map(|r| r.terms.clone());
            let chemistry = terms
                .as_ref()
                .map(|_| crate::reaction::Chemistry::of(package))
                .transpose()
                .map_err(|e| error(format!("{}: {e}", u.name)))?;
            let balance = crate::engine::balance(
                &inputs,
                &outputs,
                unit.duty_w,
                unit.material_exchange.as_ref(),
                package.component_ids().len(),
                terms.as_ref().zip(chemistry.as_ref()),
            );
            match balance {
                Ok(balance) => unit.balance = Some(balance),
                Err(e) => return Ok(Some(format!("{}: {e}", u.name))),
            }
            validate_empty_thermal_spec(&u.spec, &inputs, &outputs, package).map_err(|e| {
                SolveError {
                    unit: Some(u.id),
                    ..error(format!("{}: {e}", u.name))
                }
            })?;
        }
        Ok(None)
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
