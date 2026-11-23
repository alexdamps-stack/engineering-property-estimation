//! Deterministic sequential-modular execution. No GUI dependencies.
use crate::units::{self, CalcContext, UnitInputs, UnitOutput};
use crate::{
    CalculationStatus, Diagnostic, EnergySourceDuty, ExchangerSpec, Flowsheet, MixerPressure,
    PortDirection, PortMedium, PortRef, PressureSpec, SeparatorThermal, Stream, StreamId,
    ThermalSpec, Unit, UnitId, UnitKind, UnitSpec,
};
use property_packages::{CubicPackage, EquilibriumState};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::ControlFlow;
use thermo_types::{Composition, CompositionBasis, Phase};

/// Solved energy stream (#171): a heat duty or shaft work in W, positive in
/// the stream direction, that is, into the destination unit and out of the
/// source unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyStream {
    /// Duty in W, positive into the destination unit.
    pub duty_w: f64,
}

/// Feed fractions include identities and basis; they must match package order.
pub type FeedCompositions = BTreeMap<UnitId, Composition>;

/// Material result. An absent separator phase has zero flow and no state.
/// A splitter may carry a known state at zero flow. No thermodynamic property
/// is inferred for an absent phase.
#[derive(Debug, Clone)]
pub struct MaterialStream {
    pub molar_flow_mol_s: f64,
    pub state: Option<EquilibriumState>,
}

impl MaterialStream {
    /// Mass flow from the phase molar masses (rho * molar volume), kg/s.
    /// The empty outlet carries zero mass without a fabricated molar mass.
    pub fn mass_flow_kg_s(&self) -> f64 {
        self.state.as_ref().map_or(0.0, |s| {
            self.molar_flow_mol_s
                * s.phases
                    .iter()
                    .map(|p| p.phase_fraction * p.state.mass_density * p.state.molar_volume)
                    .sum::<f64>()
        })
    }

    pub fn enthalpy_flow_w(&self) -> f64 {
        self.state
            .as_ref()
            .map_or(0.0, |s| self.molar_flow_mol_s * s.enthalpy())
    }

    /// Thermal energy-flow scale n·cp·T in W, independent of the caloric
    /// reference point (phase-weighted cp, latent heat not included). Used to
    /// scale enthalpy comparisons where |n·h| alone can be close to zero.
    pub fn thermal_scale_w(&self) -> f64 {
        self.state
            .as_ref()
            .map_or(0.0, |s| self.molar_flow_mol_s * molar_thermal_scale(s))
    }

    pub fn component_flows_mol_s(&self, count: usize) -> Vec<f64> {
        self.state.as_ref().map_or_else(
            || vec![0.0; count],
            |s| {
                s.mole_fractions
                    .iter()
                    .map(|x| x * self.molar_flow_mol_s)
                    .collect()
            },
        )
    }
}

/// Molar thermal scale cp·T in J/mol (phase-weighted cp). Nonfinite or
/// nonpositive values return 0, so callers keep their absolute floors.
pub(crate) fn molar_thermal_scale(s: &EquilibriumState) -> f64 {
    let cp: f64 = s.phases.iter().map(|p| p.phase_fraction * p.state.cp).sum();
    let scale = cp * s.temperature_k;
    if scale.is_finite() && scale > 0.0 {
        scale
    } else {
        0.0
    }
}

/// Material crossing the modeled flowsheet boundary at a logical unit.
/// Flow is positive into the unit. Enthalpy flow uses that same direction
/// and the package's signed molar enthalpy.
#[derive(Debug, Clone)]
pub struct MaterialExchange {
    pub molar_flow_mol_s: f64,
    /// Signed molar flows in package component order.
    pub component_flows_mol_s: Vec<f64>,
    pub enthalpy_flow_w: f64,
}

/// Outlet minus inlet residuals after accounting for external material and
/// duty. Boundary feeds/products have no internal balance.
#[derive(Debug, Clone)]
pub struct Balance {
    pub component_residuals_mol_s: Vec<f64>,
    pub max_scaled_component_residual: f64,
    pub energy_residual_w: f64,
    pub scaled_energy_residual: f64,
}

#[derive(Debug, Clone)]
pub struct UnitCalculation {
    /// External duty, positive into the unit; it enters the energy balance.
    pub duty_w: f64,
    /// Heat transferred internally from the hot to the cold side of a
    /// heat exchanger, W. `None` for other units. It is not an external duty
    /// and does not enter the energy balance.
    pub transferred_duty_w: Option<f64>,
    /// Additional flowsheet boundary contribution from material and its enthalpy.
    pub material_exchange: Option<MaterialExchange>,
    pub balance: Option<Balance>,
    pub warnings: Vec<String>,
}

/// Atomic successful calculation. Failure returns no partially solved report.
/// Inputs and package identity are retained to detect obsolete results.
#[derive(Debug, Clone)]
pub struct Calculation {
    topology: Flowsheet,
    feeds: FeedCompositions,
    package: CubicPackage,
    pub order: Vec<UnitId>,
    pub streams: BTreeMap<StreamId, MaterialStream>,
    /// Solved energy streams; material streams are in `streams`.
    pub energy_streams: BTreeMap<StreamId, EnergyStream>,
    pub units: BTreeMap<UnitId, UnitCalculation>,
    /// Tear-stream closure for a successful recycle calculation.
    pub recycle: Option<crate::RecycleConvergence>,
}

impl Calculation {
    /// Whether this result still answers the given inputs. Only solve-relevant
    /// content is compared ([`Flowsheet::same_calculation_inputs`]): renaming a
    /// unit or stream, or a revision bump alone, does not make results stale.
    pub fn is_current(
        &self,
        flowsheet: &Flowsheet,
        feeds: &FeedCompositions,
        package: &CubicPackage,
    ) -> bool {
        self.topology.same_calculation_inputs(flowsheet)
            && self.feeds == *feeds
            && self.package.same_configuration(package)
    }
}

/// A preflight problem or a named unit failure. Numerical errors remain visible.
#[derive(Debug, Clone)]
pub struct SolveError {
    pub unit: Option<UnitId>,
    pub message: String,
    pub diagnostics: Vec<Diagnostic>,
    /// The caller stopped the solve through [`Flowsheet::solve_with_progress`];
    /// no numerical failure occurred and no result is returned.
    pub cancelled: bool,
}

/// Progress of a running solve, reported before each unit is calculated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SolveProgress {
    /// Recycle pass (1-based), or zero for an acyclic solve/initialization.
    pub iteration: usize,
    /// Units already calculated.
    pub completed: usize,
    /// Units in the calculation order.
    pub total: usize,
    /// The unit about to be calculated.
    pub unit: UnitId,
}

impl SolveError {
    /// The error returned when the progress callback stops a solve.
    pub(crate) fn cancelled() -> Self {
        Self {
            unit: None,
            message: "Calculation cancelled; no results were produced.".into(),
            diagnostics: Vec::new(),
            cancelled: true,
        }
    }
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}
impl std::error::Error for SolveError {}

impl Flowsheet {
    /// Validates the full topology/specification set before returning a stable
    /// topological order. Multiple connections between two units are counted
    /// individually; branched/rejoined graphs are supported, recycles rejected.
    pub fn calculation_order(&self) -> Result<Vec<UnitId>, SolveError> {
        let diagnostics = self.diagnostics();
        if !diagnostics.is_empty() {
            return Err(SolveError {
                unit: None,
                message: diagnostics
                    .iter()
                    .map(|d| d.message.as_str())
                    .collect::<Vec<_>>()
                    .join("\n"),
                diagnostics,
                cancelled: false,
            });
        }
        Ok(self.kahn_order(|_| true))
    }

    /// Kahn topological order over the streams accepted by `keep`, ties broken by
    /// unit identifier. Each kept connection counts individually. Units on a cycle
    /// of kept streams are omitted, so a short result means the kept graph is cyclic.
    pub(crate) fn kahn_order(&self, keep: impl Fn(&Stream) -> bool) -> Vec<UnitId> {
        let mut incoming: BTreeMap<_, usize> = self.units().map(|u| (u.id, 0)).collect();
        for s in self.streams().filter(|s| keep(s)) {
            *incoming.get_mut(&s.to.unit).expect("checked topology") += 1;
        }
        let mut ready: BTreeSet<_> = incoming
            .iter()
            .filter_map(|(&id, &n)| (n == 0).then_some(id))
            .collect();
        let mut order = Vec::new();
        while let Some(id) = ready.pop_first() {
            order.push(id);
            for s in self.streams().filter(|s| s.from.unit == id && keep(s)) {
                let n = incoming.get_mut(&s.to.unit).expect("checked topology");
                *n -= 1;
                if *n == 0 {
                    ready.insert(s.to.unit);
                }
            }
        }
        order
    }

    /// Solves all units with a single immutable package. Each call computes a
    /// fresh report, so a failed solve cannot reuse old downstream results.
    pub fn solve(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
    ) -> Result<Calculation, SolveError> {
        self.solve_with_progress(package, feeds, |_| ControlFlow::Continue(()))
    }

    /// [`Self::solve`] with a callback before each unit calculation. Returning
    /// [`ControlFlow::Break`] stops the solve between units and returns a
    /// [`SolveError`] with `cancelled` set; no partial result is returned.
    /// The callback cannot interrupt a property calculation in progress.
    pub fn solve_with_progress(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
        progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Result<Calculation, SolveError> {
        let order = self.calculation_order()?;
        self.validate_inputs(package, feeds)?;
        self.solve_pass(package, feeds, order, &BTreeMap::new(), progress)
    }

    pub(crate) fn validate_inputs(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
    ) -> Result<(), SolveError> {
        let fail = |id: UnitId, message: String| SolveError {
            unit: Some(id),
            message: format!("{}: {message}", self.unit(id).expect("unit exists").name),
            diagnostics: Vec::new(),
            cancelled: false,
        };
        // Complete preflight before invoking any property calculation.
        for u in self.units() {
            if let CalculationStatus::DrawnOnly(reason) = u.kind.calculation_status() {
                let issue = u
                    .kind
                    .issue()
                    .map_or_else(String::new, |n| format!(" (#{n})"));
                return Err(fail(u.id, format!("{} {reason}{issue}", u.kind)));
            }
            if u.kind == UnitKind::Feed {
                let c = feeds
                    .get(&u.id)
                    .ok_or_else(|| fail(u.id, "missing feed composition".into()))?;
                if c.basis() != CompositionBasis::Mole
                    || c.entries()
                        .iter()
                        .map(|(id, _)| id.as_str())
                        .collect::<Vec<_>>()
                        != package.model().component_ids()
                {
                    return Err(fail(
                        u.id,
                        "feed component identities/order and mole basis must match the package"
                            .into(),
                    ));
                }
                package
                    .validate_composition(&c.entries().iter().map(|(_, x)| *x).collect::<Vec<_>>())
                    .map_err(|e| fail(u.id, e.to_string()))?;
            }
        }
        if feeds
            .keys()
            .any(|id| self.unit(*id).is_none_or(|u| u.kind != UnitKind::Feed))
        {
            return Err(SolveError {
                unit: None,
                message: "composition supplied for an unknown or non-feed unit".into(),
                diagnostics: Vec::new(),
                cancelled: false,
            });
        }
        Ok(())
    }

    /// One sequential pass over `order`. Callers must have run
    /// [`Self::validate_inputs`] once for these inputs; recycle iterations
    /// reuse the same inputs and do not repeat it.
    pub(crate) fn solve_pass(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
        order: Vec<UnitId>,
        tears: &BTreeMap<StreamId, MaterialStream>,
        mut progress: impl FnMut(SolveProgress) -> ControlFlow<()>,
    ) -> Result<Calculation, SolveError> {
        let fail = |id: UnitId, message: String| SolveError {
            unit: Some(id),
            message: format!("{}: {message}", self.unit(id).expect("unit exists").name),
            diagnostics: Vec::new(),
            cancelled: false,
        };
        let mut report = Calculation {
            topology: self.clone(),
            feeds: feeds.clone(),
            package: package.clone(),
            order: order.clone(),
            streams: BTreeMap::new(),
            energy_streams: BTreeMap::new(),
            units: BTreeMap::new(),
            recycle: None,
        };
        let total = order.len();
        for (completed, id) in order.into_iter().enumerate() {
            let step = SolveProgress {
                iteration: 0,
                completed,
                total,
                unit: id,
            };
            if progress(step).is_break() {
                return Err(SolveError::cancelled());
            }
            let u = self.unit(id).expect("checked topology");
            let ports = self.unit_ports(u);
            let inlets: Vec<MaterialStream> = ports
                .material_inlets
                .iter()
                .map(|s| match s {
                    Some(s) => tears.get(s).unwrap_or_else(|| &report.streams[s]).clone(),
                    None => empty(),
                })
                .collect();
            // `None`: unconnected, or a source that takes its duty from this unit.
            let energy_in: Vec<Option<f64>> = ports
                .energy
                .iter()
                .filter(|p| p.direction == PortDirection::Inlet)
                .map(|p| {
                    p.stream
                        .and_then(|s| report.energy_streams.get(&s).map(|e| e.duty_w))
                })
                .collect();
            let inputs = UnitInputs {
                material: &inlets,
                energy: &energy_in,
            };
            let ctx = CalcContext::new(package, feeds.get(&id));
            let UnitOutput {
                outlets,
                duty_w,
                transferred_duty_w,
                material_exchange,
                warnings,
                energy,
            } = calculate_unit(u.kind, &u.spec, &inputs, &ctx).map_err(|e| fail(id, e))?;
            if !duty_w.is_finite()
                || transferred_duty_w.is_some_and(|q| !q.is_finite())
                || energy.iter().flatten().any(|q| !q.is_finite())
                || outlets.iter().any(|s| {
                    !s.molar_flow_mol_s.is_finite()
                        || !s.enthalpy_flow_w().is_finite()
                        || !s.mass_flow_kg_s().is_finite()
                })
            {
                return Err(fail(
                    id,
                    "nonfinite extensive stream property or duty".into(),
                ));
            }
            if outlets.len() != ports.material_outlets.len() {
                return Err(fail(
                    id,
                    format!(
                        "calculation returned {} outlets for {} outlet ports",
                        outlets.len(),
                        ports.material_outlets.len()
                    ),
                ));
            }
            let balance = if u.kind.is_boundary() {
                None
            } else {
                Some(
                    balance(
                        &inlets,
                        &outlets,
                        duty_w,
                        material_exchange.as_ref(),
                        package.model().fractions().len(),
                    )
                    .map_err(|e| fail(id, e))?,
                )
            };
            transfer_energy(&ports.energy, duty_w, energy, &mut report.energy_streams)
                .map_err(|e| fail(id, e))?;
            for (stream, result) in ports.material_outlets.iter().zip(outlets) {
                if let Some(s) = stream {
                    report.streams.insert(*s, result);
                }
            }
            report.units.insert(
                id,
                UnitCalculation {
                    duty_w,
                    transferred_duty_w,
                    material_exchange,
                    balance,
                    warnings,
                },
            );
        }
        Ok(report)
    }
}

/// One energy port of a unit with the stream attached to it, if any.
#[derive(Debug, Clone, Copy)]
pub(crate) struct EnergyPort {
    pub(crate) direction: PortDirection,
    pub(crate) stream: Option<StreamId>,
}

/// Streams at a unit's ports, by medium and direction, in port order.
/// Unconnected (optional) ports are `None`.
#[derive(Debug, Clone, Default)]
pub(crate) struct UnitPorts {
    pub(crate) material_inlets: Vec<Option<StreamId>>,
    pub(crate) material_outlets: Vec<Option<StreamId>>,
    pub(crate) energy: Vec<EnergyPort>,
}

impl Flowsheet {
    /// Streams attached to each port of `unit`.
    pub(crate) fn unit_ports(&self, unit: &Unit) -> UnitPorts {
        let mut out = UnitPorts::default();
        for (i, p) in unit.ports().iter().enumerate() {
            let stream = self.stream_at(PortRef::new(unit.id, i)).map(|s| s.id);
            match (p.medium, p.direction) {
                (PortMedium::Material, PortDirection::Inlet) => out.material_inlets.push(stream),
                (PortMedium::Material, PortDirection::Outlet) => out.material_outlets.push(stream),
                (PortMedium::Energy, direction) => {
                    out.energy.push(EnergyPort { direction, stream })
                }
            }
        }
        out
    }
}

/// Relative tolerance between the duty an energy stream delivers and the
/// duty the receiving unit used.
const ENERGY_TOLERANCE: f64 = 1e-9;

/// Records the energy-port values of a calculated unit (#171).
///
/// `values` holds one entry per energy port (positive in the port's
/// direction); empty means the default: a single energy port carries the
/// unit's `duty_w` (negated for an outlet). An outlet value becomes the
/// stream's duty. An inlet whose stream already has a duty must match the
/// value the unit used; an inlet whose stream has none (its source is an
/// energy source set to "calculated") receives the unit's value.
fn transfer_energy(
    ports: &[EnergyPort],
    duty_w: f64,
    values: Vec<Option<f64>>,
    streams: &mut BTreeMap<StreamId, EnergyStream>,
) -> Result<(), String> {
    let values = if !values.is_empty() {
        values
    } else if let [port] = ports {
        vec![Some(match port.direction {
            PortDirection::Inlet => duty_w,
            PortDirection::Outlet => -duty_w,
        })]
    } else {
        vec![None; ports.len()]
    };
    if values.len() != ports.len() {
        return Err(format!(
            "calculation returned {} energy values for {} energy ports",
            values.len(),
            ports.len()
        ));
    }
    for (port, value) in ports.iter().zip(values) {
        let Some(id) = port.stream else { continue };
        match (port.direction, streams.get(&id), value) {
            // An energy source left to the connected unit leaves no value.
            (PortDirection::Outlet, _, None) => {}
            (PortDirection::Outlet, _, Some(q)) => {
                streams.insert(id, EnergyStream { duty_w: q });
            }
            (PortDirection::Inlet, Some(given), Some(q)) => {
                if (given.duty_w - q).abs() > ENERGY_TOLERANCE * given.duty_w.abs().max(1.0) {
                    return Err(format!(
                        "energy stream {id} delivers {:.6e} W but the unit used {q:.6e} W",
                        given.duty_w
                    ));
                }
            }
            (PortDirection::Inlet, None, Some(q)) => {
                streams.insert(id, EnergyStream { duty_w: q });
            }
            (PortDirection::Inlet, _, None) => {
                return Err(format!("the unit did not determine energy stream {id}"));
            }
        }
    }
    Ok(())
}

/// Zero-flow outlet without a thermodynamic state.
pub(crate) fn empty() -> MaterialStream {
    MaterialStream {
        molar_flow_mol_s: 0.0,
        state: None,
    }
}

const ZERO_FLOW_WARNING: &str = "Zero inlet flow: passed through as a zero-flow outlet without a thermodynamic state; no duty applied.";

/// Calculates one unit. Core kinds are calculated here; group kinds are
/// routed to their module in [`crate::units`].
fn calculate_unit(
    kind: UnitKind,
    spec: &UnitSpec,
    unit_inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    match spec {
        UnitSpec::Pressurize { .. } | UnitSpec::Expand { .. } => {
            return units::rotating::calculate_core(kind, spec, unit_inputs, ctx);
        }
        UnitSpec::Rotating(s) => return units::rotating::calculate(kind, s, unit_inputs, ctx),
        UnitSpec::HeatTransport(s) => {
            return units::heat_transport::calculate(kind, s, unit_inputs, ctx);
        }
        UnitSpec::SimpleLogical(s) => {
            return units::simple_logical::calculate(kind, s, unit_inputs, ctx);
        }
        UnitSpec::Reactors(s) => return units::reactors::calculate(kind, s, unit_inputs, ctx),
        UnitSpec::Columns(s) => return units::columns::calculate(kind, s, unit_inputs, ctx),
        UnitSpec::Multiphase(s) => return units::multiphase::calculate(kind, s, unit_inputs, ctx),
        UnitSpec::Solids(s) => return units::solids::calculate(kind, s, unit_inputs, ctx),
        _ => {}
    }
    let inputs = unit_inputs.material;
    let feed: Option<&Composition> = ctx.feed();
    let pkg: &CubicPackage = ctx.package();
    let spec = spec.clone();
    fn state(s: &MaterialStream) -> Result<&EquilibriumState, String> {
        s.state.as_ref().ok_or_else(|| {
            "inlet has no material phase; cannot apply a thermodynamic specification".to_string()
        })
    }
    let positive_flow = |s: &MaterialStream| {
        if s.molar_flow_mol_s > 0.0 {
            Ok(s.molar_flow_mol_s)
        } else {
            Err("a thermal operation requires positive inlet flow".to_string())
        }
    };
    // An empty separator outlet or a splitter endpoint is a valid zero-flow
    // stream; it passes through single-inlet units unchanged (section 2.3).
    let zero = |s: &MaterialStream| s.molar_flow_mol_s == 0.0;
    let pass_through =
        |count: usize| UnitOutput::new(vec![empty(); count], 0.0, vec![ZERO_FLOW_WARNING.into()]);
    let pressure = |p: f64, drop: f64| {
        if (p - drop).is_finite() && p - drop > 0.0 {
            Ok(p - drop)
        } else {
            Err("pressure drop leaves nonpositive absolute outlet pressure".to_string())
        }
    };
    let ph =
        |s: &EquilibriumState, p, h| pkg.ph(&s.mole_fractions, p, h).map_err(|e| e.to_string());
    let tp =
        |s: &EquilibriumState, p, t| pkg.tp(&s.mole_fractions, t, p).map_err(|e| e.to_string());
    let stream = |flow, state| MaterialStream {
        molar_flow_mol_s: flow,
        state: Some(state),
    };
    match spec {
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(p),
            molar_flow_mol_s: Some(n),
        } => {
            let z: Vec<_> = feed
                .expect("preflight feed")
                .entries()
                .iter()
                .map(|(_, x)| *x)
                .collect();
            Ok(UnitOutput::new(
                vec![stream(n, pkg.tp(&z, t, p).map_err(|e| e.to_string())?)],
                0.0,
                vec![],
            ))
        }
        UnitSpec::Product => Ok(UnitOutput::new(vec![], 0.0, vec![])),
        // Energy boundaries (#171): a specified source sets its stream; a
        // calculated one leaves it to the connected unit; a sink records
        // what it receives.
        UnitSpec::EnergySource { duty } => Ok(UnitOutput {
            energy: vec![match duty {
                EnergySourceDuty::Specified(Some(q)) => Some(q),
                EnergySourceDuty::Specified(None) => {
                    return Err("missing energy source duty".into());
                }
                EnergySourceDuty::Calculated => None,
            }],
            ..UnitOutput::new(vec![], 0.0, vec![])
        }),
        UnitSpec::EnergySink => match unit_inputs.energy {
            [Some(q)] => Ok(UnitOutput {
                energy: vec![Some(*q)],
                ..UnitOutput::new(vec![], 0.0, vec![])
            }),
            _ => Err("energy sink receives no determined duty".into()),
        },
        UnitSpec::Mixer { outlet_pressure } => {
            let active: Vec<_> = inputs.iter().filter(|s| s.molar_flow_mol_s > 0.0).collect();
            let n: f64 = active.iter().map(|s| s.molar_flow_mol_s).sum();
            if n <= 0.0 || !n.is_finite() {
                return Err("mixer needs a finite positive total inlet flow".into());
            }
            let mut pmin = f64::INFINITY;
            let mut z = vec![0.0; pkg.model().fractions().len()];
            let mut h = 0.0;
            for s in active {
                let a = state(s)?;
                pmin = pmin.min(a.pressure_pa);
                for (out, x) in z.iter_mut().zip(&a.mole_fractions) {
                    *out += s.molar_flow_mol_s / n * x;
                }
                h += s.molar_flow_mol_s / n * a.enthalpy();
            }
            let p = match outlet_pressure {
                MixerPressure::LowestInlet => pmin,
                MixerPressure::Specified(Some(p)) if p <= pmin => p,
                _ => {
                    return Err(
                        "mixer outlet pressure exceeds lowest flowing inlet pressure".into(),
                    );
                }
            };
            Ok(UnitOutput::new(
                vec![stream(n, pkg.ph(&z, p, h).map_err(|e| e.to_string())?)],
                0.0,
                vec![],
            ))
        }
        UnitSpec::HeatExchange { .. } if zero(&inputs[0]) => Ok(pass_through(1)),
        UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa,
        } => {
            let a = state(&inputs[0])?;
            let n = positive_flow(&inputs[0])?;
            let p = pressure(a.pressure_pa, pressure_drop_pa)?;
            // A heater may take its duty from the energy stream at its
            // energy inlet (#171); otherwise the duty is specified or follows
            // from the outlet temperature.
            let specified = match thermal {
                ThermalSpec::Duty(Some(q)) => Some(q),
                ThermalSpec::EnergyStream => match unit_inputs.energy {
                    [Some(q)] => Some(*q),
                    _ => return Err("no duty is supplied by an energy stream".into()),
                },
                _ => None,
            };
            let out = match (thermal, specified) {
                (_, Some(q)) => ph(a, p, a.enthalpy() + q / n)?,
                (ThermalSpec::OutletTemperature(Some(t)), _) => tp(a, p, t)?,
                _ => return Err("missing thermal specification".into()),
            };
            let q = specified.unwrap_or_else(|| n * (out.enthalpy() - a.enthalpy()));
            // A duty sign that contradicts the kind is rejected before solving
            // (`Unit::issues`); an outlet temperature on the wrong side of the
            // inlet is only known here and depends on the iterate inside a
            // recycle, so it is reported on the result instead of aborting.
            let (t_in, t_out) = (a.temperature_k, out.temperature_k);
            let warnings = match kind {
                UnitKind::Heater if q < 0.0 || t_out < t_in => vec![format!(
                    "Heater removes heat: outlet temperature {t_out:.2} K is below inlet {t_in:.2} K (duty {q:.6e} W); use a cooler or check the specification."
                )],
                UnitKind::Cooler if q > 0.0 || t_out > t_in => vec![format!(
                    "Cooler adds heat: outlet temperature {t_out:.2} K is above inlet {t_in:.2} K (duty {q:.6e} W); use a heater or check the specification."
                )],
                _ => vec![],
            };
            Ok(UnitOutput::new(vec![stream(n, out)], q, warnings))
        }
        UnitSpec::Valve { .. } if zero(&inputs[0]) => Ok(pass_through(1)),
        UnitSpec::Valve { pressure: spec } => {
            let a = state(&inputs[0])?;
            let p = match spec {
                PressureSpec::OutletPressure(Some(p)) if p <= a.pressure_pa => p,
                PressureSpec::PressureDrop(Some(drop)) => pressure(a.pressure_pa, drop)?,
                _ => return Err("valve outlet pressure exceeds inlet pressure".into()),
            };
            Ok(UnitOutput::new(
                vec![stream(inputs[0].molar_flow_mol_s, ph(a, p, a.enthalpy())?)],
                0.0,
                vec![],
            ))
        }
        UnitSpec::Separator { .. } if zero(&inputs[0]) => Ok(pass_through(2)),
        UnitSpec::Separator {
            thermal,
            pressure_drop_pa,
        } => {
            let a = state(&inputs[0])?;
            let n = positive_flow(&inputs[0])?;
            let p = pressure(a.pressure_pa, pressure_drop_pa)?;
            // Preserve the actual inlet equilibrium phases at unchanged conditions;
            // reflashing a saturation endpoint can select an arbitrary root.
            let out = match thermal {
                SeparatorThermal::Adiabatic if pressure_drop_pa == 0.0 => a.clone(),
                SeparatorThermal::Adiabatic => ph(a, p, a.enthalpy())?,
                SeparatorThermal::OutletTemperature(Some(t)) => tp(a, p, t)?,
                SeparatorThermal::Duty(Some(q)) => ph(a, p, a.enthalpy() + q / n)?,
                _ => return Err("missing separator thermal specification".into()),
            };
            let q = match thermal {
                SeparatorThermal::Adiabatic => 0.0,
                SeparatorThermal::Duty(Some(q)) => q,
                _ => n * (out.enthalpy() - a.enthalpy()),
            };
            let mut outlets = vec![empty(); 2];
            for phase in &out.phases {
                if phase.phase_fraction == 0.0 {
                    continue;
                }
                let index = match phase.phase {
                    Phase::Vapor => 0,
                    Phase::Liquid => 1,
                    Phase::Supercritical => 0,
                };
                if outlets[index].state.is_some() {
                    return Err("separator requires at most one vapor and one liquid phase".into());
                }
                let mut single = phase.clone();
                single.phase_fraction = 1.0;
                outlets[index] = stream(
                    n * phase.phase_fraction,
                    EquilibriumState {
                        mole_fractions: phase.composition.clone(),
                        phases: vec![single],
                        ..out.clone()
                    },
                );
            }
            let warnings = if out.phases.iter().any(|p| p.phase == Phase::Supercritical) {
                vec!["Single supercritical fluid routed to the vapor port by convention; phase label retained. Absent liquid outlet has zero flow and no thermodynamic state.".into()]
            } else if outlets.iter().any(|s| s.state.is_none()) {
                vec!["Single phase: absent separator outlet has zero flow and no thermodynamic state.".into()]
            } else {
                vec![]
            };
            Ok(UnitOutput::new(outlets, q, warnings))
        }
        UnitSpec::Splitter {
            fraction_to_first: Some(f),
        } => {
            let mut a = inputs[0].clone();
            let mut b = a.clone();
            a.molar_flow_mol_s *= f;
            b.molar_flow_mol_s *= 1.0 - f;
            Ok(UnitOutput::new(vec![a, b], 0.0, vec![]))
        }
        UnitSpec::HeatExchanger {
            thermal,
            hot_pressure_drop_pa,
            cold_pressure_drop_pa,
        } if zero(&inputs[0]) || zero(&inputs[1]) => {
            if matches!(
                thermal,
                ExchangerSpec::HotEndApproach(_) | ExchangerSpec::ColdEndApproach(_)
            ) {
                return Err(
                    "fixed exchanger approach requires positive flow on both inlet sides".into(),
                );
            }
            // Without flow on one side no heat can be transferred: the empty
            // side passes through empty, the flowing side keeps its enthalpy
            // and only takes its pressure drop.
            let mut outlets = Vec::new();
            for (s, drop) in inputs
                .iter()
                .zip([hot_pressure_drop_pa, cold_pressure_drop_pa])
            {
                outlets.push(if zero(s) {
                    empty()
                } else {
                    let a = state(s)?;
                    let n = positive_flow(s)?;
                    if drop == 0.0 {
                        stream(n, a.clone())
                    } else {
                        stream(n, ph(a, pressure(a.pressure_pa, drop)?, a.enthalpy())?)
                    }
                });
            }
            Ok(UnitOutput {
                outlets,
                duty_w: 0.0,
                transferred_duty_w: Some(0.0),
                material_exchange: None,
                energy: Vec::new(),
                warnings: vec![
                    "Zero flow on one or both exchanger sides: no heat transferred; empty sides pass through as zero-flow outlets without a thermodynamic state.".into(),
                ],
            })
        }
        UnitSpec::HeatExchanger {
            thermal,
            hot_pressure_drop_pa,
            cold_pressure_drop_pa,
        } => {
            let hot = state(&inputs[0])?;
            let cold = state(&inputs[1])?;
            let nh = positive_flow(&inputs[0])?;
            let nc = positive_flow(&inputs[1])?;
            let hp = pressure(hot.pressure_pa, hot_pressure_drop_pa)?;
            let cp = pressure(cold.pressure_pa, cold_pressure_drop_pa)?;
            let approach_duty = |flow, h_high: f64, h_low: f64| {
                let delta_h = h_high - h_low;
                let q: f64 = flow * delta_h;
                // At zero transfer, allow only roundoff in the two TP enthalpies.
                let roundoff = 32.0 * f64::EPSILON * h_high.abs().max(h_low.abs()).max(1.0);
                if !q.is_finite() || delta_h < -roundoff {
                    return Err(
                        "fixed exchanger approach requires finite, nonnegative hot-to-cold duty"
                            .to_string(),
                    );
                }
                Ok(q.max(0.0))
            };
            let (q, hot_target, cold_target) = match thermal {
                ExchangerSpec::Duty(Some(q)) => (q, None, None),
                ExchangerSpec::HotOutletTemperature(Some(t)) => (
                    nh * (hot.enthalpy() - tp(hot, hp, t)?.enthalpy()),
                    None,
                    None,
                ),
                ExchangerSpec::ColdOutletTemperature(Some(t)) => (
                    nc * (tp(cold, cp, t)?.enthalpy() - cold.enthalpy()),
                    None,
                    None,
                ),
                ExchangerSpec::HotEndApproach(Some(dt)) => {
                    let out = tp(cold, cp, hot.temperature_k - dt)?;
                    let q = approach_duty(nc, out.enthalpy(), cold.enthalpy())?;
                    (q, None, Some(out))
                }
                ExchangerSpec::ColdEndApproach(Some(dt)) => {
                    let out = tp(hot, hp, cold.temperature_k + dt)?;
                    let q = approach_duty(nh, hot.enthalpy(), out.enthalpy())?;
                    (q, Some(out), None)
                }
                _ => return Err("missing exchanger thermal specification".into()),
            };
            let fixed_approach = hot_target.is_some() || cold_target.is_some();
            let ho = match hot_target {
                Some(out) => out,
                None => ph(hot, hp, hot.enthalpy() - q / nh)?,
            };
            let co = match cold_target {
                Some(out) => out,
                None => ph(cold, cp, cold.enthalpy() + q / nc)?,
            };
            if fixed_approach {
                // Allow numerical noise from the PH inversion at an exact pinch.
                const TERMINAL_TOLERANCE_K: f64 = 1e-7;
                if hot.temperature_k - co.temperature_k < -TERMINAL_TOLERANCE_K
                    || ho.temperature_k - cold.temperature_k < -TERMINAL_TOLERANCE_K
                {
                    return Err(
                        "fixed exchanger approach requires no terminal temperature crossover"
                            .into(),
                    );
                }
            } else if q < 0.0
                || hot.temperature_k < cold.temperature_k
                || ho.temperature_k < cold.temperature_k
                || co.temperature_k > hot.temperature_k
            {
                return Err("heat exchanger requires hot-to-cold duty without terminal temperature crossover".into());
            }
            // The transferred heat is internal: the external duty stays zero.
            Ok(UnitOutput {
                outlets: vec![stream(nh, ho), stream(nc, co)],
                duty_w: 0.0,
                transferred_duty_w: Some(q),
                material_exchange: None,
                energy: Vec::new(),
                warnings: vec![],
            })
        }
        _ => Err("unsupported or incomplete unit specification".into()),
    }
}

pub(crate) fn balance(
    inputs: &[MaterialStream],
    outputs: &[MaterialStream],
    duty: f64,
    exchange: Option<&MaterialExchange>,
    count: usize,
) -> Result<Balance, String> {
    let mut residual = vec![0.0; count];
    let mut energy = -duty;
    let mut energy_scale = duty.abs();
    let mut thermal_scale = 0.0;
    let mut external_inflow = 0.0;
    if let Some(x) = exchange {
        if x.component_flows_mol_s.len() != count
            || !x.molar_flow_mol_s.is_finite()
            || !x.enthalpy_flow_w.is_finite()
            || x.component_flows_mol_s.iter().any(|n| !n.is_finite())
        {
            return Err("invalid external material exchange".into());
        }
        for (r, n) in residual.iter_mut().zip(&x.component_flows_mol_s) {
            *r -= n;
        }
        energy -= x.enthalpy_flow_w;
        energy_scale += x.enthalpy_flow_w.abs();
        external_inflow = x.molar_flow_mol_s.max(0.0);
    }
    for (sign, streams) in [(-1.0, inputs), (1.0, outputs)] {
        for s in streams {
            if !s.molar_flow_mol_s.is_finite() || s.molar_flow_mol_s < 0.0 {
                return Err("nonfinite or negative outlet flow".into());
            }
            for (r, n) in residual.iter_mut().zip(s.component_flows_mol_s(count)) {
                *r += sign * n;
            }
            let h = s.enthalpy_flow_w();
            energy += sign * h;
            energy_scale += h.abs();
            thermal_scale += s.thermal_scale_w();
        }
    }
    let flow_scale = (inputs.iter().map(|s| s.molar_flow_mol_s).sum::<f64>() + external_inflow)
        .max(f64::MIN_POSITIVE);
    let max_scaled = residual
        .iter()
        .map(|r| r.abs() / flow_scale)
        .fold(0.0_f64, f64::max);
    // |Q| + sum |n h| depends on the caloric reference point and can collapse
    // when h is near zero; sum n cp T keeps a reference-independent scale.
    let scaled_energy = energy / energy_scale.max(thermal_scale).max(1.0);
    if residual.iter().any(|r| !r.is_finite())
        || !energy.is_finite()
        || !energy_scale.is_finite()
        || !thermal_scale.is_finite()
        || !flow_scale.is_finite()
        || !scaled_energy.is_finite()
        || max_scaled > 1e-8
        || scaled_energy.abs() > 2e-9
    {
        return Err("component or energy balance failed (1e-8 relative component; 2e-9 scaled energy budget)".into());
    }
    Ok(Balance {
        component_residuals_mol_s: residual,
        max_scaled_component_residual: max_scaled,
        energy_residual_w: energy,
        scaled_energy_residual: scaled_energy,
    })
}
