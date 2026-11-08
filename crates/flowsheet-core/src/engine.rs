//! Deterministic, acyclic, sequential-modular execution. No GUI dependencies.
use crate::{
    Diagnostic, ExchangerSpec, Flowsheet, MixerPressure, PortDirection, PortRef, PressureSpec,
    SeparatorThermal, StreamId, ThermalSpec, UnitId, UnitKind, UnitSpec,
};
use property_packages::{CubicPackage, EquilibriumState};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use thermo_types::{Composition, CompositionBasis, Phase};

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

/// Outlet minus inlet component and energy residuals; duty is positive into
/// the unit. Boundary feeds/products have no internal balance.
#[derive(Debug, Clone)]
pub struct Balance {
    pub component_residuals_mol_s: Vec<f64>,
    pub max_scaled_component_residual: f64,
    pub energy_residual_w: f64,
    pub scaled_energy_residual: f64,
}

#[derive(Debug, Clone)]
pub struct UnitCalculation {
    pub duty_w: f64,
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
    pub units: BTreeMap<UnitId, UnitCalculation>,
}

impl Calculation {
    pub fn is_current(
        &self,
        flowsheet: &Flowsheet,
        feeds: &FeedCompositions,
        package: &CubicPackage,
    ) -> bool {
        self.topology == *flowsheet
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
            });
        }
        let mut incoming: BTreeMap<_, usize> = self.units().map(|u| (u.id, 0)).collect();
        for s in self.streams() {
            *incoming.get_mut(&s.to.unit).expect("checked topology") += 1;
        }
        let mut ready: BTreeSet<_> = incoming
            .iter()
            .filter_map(|(&id, &n)| (n == 0).then_some(id))
            .collect();
        let mut order = Vec::new();
        while let Some(id) = ready.pop_first() {
            order.push(id);
            for s in self.streams().filter(|s| s.from.unit == id) {
                let n = incoming.get_mut(&s.to.unit).expect("checked topology");
                *n -= 1;
                if *n == 0 {
                    ready.insert(s.to.unit);
                }
            }
        }
        Ok(order)
    }

    /// Solves all units with a single immutable package. Each call computes a
    /// fresh report, so a failed solve cannot reuse old downstream results.
    pub fn solve(
        &self,
        package: &CubicPackage,
        feeds: &FeedCompositions,
    ) -> Result<Calculation, SolveError> {
        let order = self.calculation_order()?;
        let fail = |id: UnitId, message: String| SolveError {
            unit: Some(id),
            message: format!("{}: {message}", self.unit(id).expect("unit exists").name),
            diagnostics: Vec::new(),
        };
        // Complete preflight before invoking any property calculation.
        for u in self.units() {
            if !u.kind.calculation_planned() {
                return Err(fail(
                    u.id,
                    format!(
                        "{} requires a PS flash and efficiency model (deferred)",
                        u.kind
                    ),
                ));
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
            });
        }
        let mut report = Calculation {
            topology: self.clone(),
            feeds: feeds.clone(),
            package: package.clone(),
            order: order.clone(),
            streams: BTreeMap::new(),
            units: BTreeMap::new(),
        };
        for id in order {
            let u = self.unit(id).expect("checked topology");
            let mut inlets = Vec::new();
            for (i, p) in u.kind.ports().iter().enumerate() {
                if p.direction == PortDirection::Inlet {
                    let s = self
                        .stream_at(PortRef::new(id, i))
                        .expect("checked topology");
                    inlets.push(report.streams[&s.id].clone());
                }
            }
            let (outlets, duty_w, warnings) =
                calculate_unit(u.spec, &inlets, feeds.get(&id), package)
                    .map_err(|e| fail(id, e))?;
            if !duty_w.is_finite()
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
            let balance = if matches!(u.kind, UnitKind::Feed | UnitKind::Product) {
                None
            } else {
                Some(
                    balance(&inlets, &outlets, duty_w, package.model().fractions().len())
                        .map_err(|e| fail(id, e))?,
                )
            };
            let ports = u
                .kind
                .ports()
                .iter()
                .enumerate()
                .filter(|(_, p)| p.direction == PortDirection::Outlet);
            for ((i, _), result) in ports.zip(outlets) {
                let s = self
                    .stream_at(PortRef::new(id, i))
                    .expect("checked topology");
                report.streams.insert(s.id, result);
            }
            report.units.insert(
                id,
                UnitCalculation {
                    duty_w,
                    balance,
                    warnings,
                },
            );
        }
        Ok(report)
    }
}

type UnitOutput = (Vec<MaterialStream>, f64, Vec<String>);

fn calculate_unit(
    spec: UnitSpec,
    inputs: &[MaterialStream],
    feed: Option<&Composition>,
    pkg: &CubicPackage,
) -> Result<UnitOutput, String> {
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
            Ok((
                vec![stream(n, pkg.tp(&z, t, p).map_err(|e| e.to_string())?)],
                0.0,
                vec![],
            ))
        }
        UnitSpec::Product => Ok((vec![], 0.0, vec![])),
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
            Ok((
                vec![stream(n, pkg.ph(&z, p, h).map_err(|e| e.to_string())?)],
                0.0,
                vec![],
            ))
        }
        UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa,
        } => {
            let a = state(&inputs[0])?;
            let n = positive_flow(&inputs[0])?;
            let p = pressure(a.pressure_pa, pressure_drop_pa)?;
            let out = match thermal {
                ThermalSpec::OutletTemperature(Some(t)) => tp(a, p, t)?,
                ThermalSpec::Duty(Some(q)) => ph(a, p, a.enthalpy() + q / n)?,
                _ => return Err("missing thermal specification".into()),
            };
            let q = match thermal {
                ThermalSpec::Duty(Some(q)) => q,
                _ => n * (out.enthalpy() - a.enthalpy()),
            };
            Ok((vec![stream(n, out)], q, vec![]))
        }
        UnitSpec::Valve { pressure: spec } => {
            let a = state(&inputs[0])?;
            let p = match spec {
                PressureSpec::OutletPressure(Some(p)) if p <= a.pressure_pa => p,
                PressureSpec::PressureDrop(Some(drop)) => pressure(a.pressure_pa, drop)?,
                _ => return Err("valve outlet pressure exceeds inlet pressure".into()),
            };
            Ok((
                vec![stream(inputs[0].molar_flow_mol_s, ph(a, p, a.enthalpy())?)],
                0.0,
                vec![],
            ))
        }
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
            let mut outlets = vec![
                MaterialStream {
                    molar_flow_mol_s: 0.0,
                    state: None
                };
                2
            ];
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
            Ok((outlets, q, warnings))
        }
        UnitSpec::Splitter {
            fraction_to_first: Some(f),
        } => {
            let mut a = inputs[0].clone();
            let mut b = a.clone();
            a.molar_flow_mol_s *= f;
            b.molar_flow_mol_s *= 1.0 - f;
            Ok((vec![a, b], 0.0, vec![]))
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
            let q = match thermal {
                ExchangerSpec::Duty(Some(q)) => q,
                ExchangerSpec::HotOutletTemperature(Some(t)) => {
                    nh * (hot.enthalpy() - tp(hot, hp, t)?.enthalpy())
                }
                ExchangerSpec::ColdOutletTemperature(Some(t)) => {
                    nc * (tp(cold, cp, t)?.enthalpy() - cold.enthalpy())
                }
                _ => return Err("missing exchanger thermal specification".into()),
            };
            let ho = ph(hot, hp, hot.enthalpy() - q / nh)?;
            let co = ph(cold, cp, cold.enthalpy() + q / nc)?;
            if q < 0.0
                || hot.temperature_k < cold.temperature_k
                || ho.temperature_k < cold.temperature_k
                || co.temperature_k > hot.temperature_k
            {
                return Err("heat exchanger requires hot-to-cold duty without terminal temperature crossover".into());
            }
            Ok((vec![stream(nh, ho), stream(nc, co)], 0.0, vec![]))
        }
        _ => Err("unsupported or incomplete unit specification".into()),
    }
}

fn balance(
    inputs: &[MaterialStream],
    outputs: &[MaterialStream],
    duty: f64,
    count: usize,
) -> Result<Balance, String> {
    let mut residual = vec![0.0; count];
    let mut energy = -duty;
    let mut energy_scale = duty.abs();
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
        }
    }
    let flow_scale = inputs
        .iter()
        .map(|s| s.molar_flow_mol_s)
        .sum::<f64>()
        .max(f64::MIN_POSITIVE);
    let max_scaled = residual
        .iter()
        .map(|r| r.abs() / flow_scale)
        .fold(0.0_f64, f64::max);
    let scaled_energy = energy / energy_scale.max(1.0);
    if residual.iter().any(|r| !r.is_finite())
        || !energy.is_finite()
        || !energy_scale.is_finite()
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
