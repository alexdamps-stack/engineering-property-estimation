//! Serializable snapshots of successful calculations.

use std::collections::BTreeMap;

use flowsheet_core::{Calculation, Flowsheet, MaterialStream, StreamId, UnitId};
use serde::{Deserialize, Serialize};
use thermo_types::{DiagnosticSeverity, Phase, ValidityStatus};

use crate::{Fingerprint, Fingerprints, Provenance};

/// Whether a stored result belongs to the current inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Currency {
    /// Computed for exactly the current inputs and provenance.
    Current,
    /// A result exists but its inputs, package or provenance changed.
    Stale,
    /// No result exists for this item.
    Missing,
}

/// Result of one phase of an equilibrium state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhaseResult {
    /// `vapor`, `liquid` or `supercritical` (a label, not a stability claim).
    pub label: String,
    /// Molar phase fraction.
    pub fraction: f64,
    /// Phase mole fractions in package component order.
    pub mole_fractions: Vec<f64>,
    /// Mass density, kg/m³.
    pub density_kg_m3: f64,
    /// Compressibility factor.
    pub compressibility: f64,
    /// Molar enthalpy, J/mol, on the package caloric reference.
    pub enthalpy_j_mol: f64,
    /// Molar entropy, J/(mol K).
    pub entropy_j_mol_k: f64,
    /// Computational validity of the phase state.
    pub validity: String,
}

/// Equilibrium state of a stream.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StateResult {
    /// Temperature, K.
    pub temperature_k: f64,
    /// Absolute pressure, Pa.
    pub pressure_pa: f64,
    /// Overall mole fractions in package component order.
    pub mole_fractions: Vec<f64>,
    /// Molar vapor fraction; `None` for a single supercritical fluid.
    pub vapor_fraction: Option<f64>,
    /// Molar enthalpy, J/mol.
    pub enthalpy_j_mol: f64,
    /// Molar entropy, J/(mol K).
    pub entropy_j_mol_k: f64,
    /// Phases with nonzero fraction.
    pub phases: Vec<PhaseResult>,
    /// Equilibrium solver iterations.
    pub iterations: usize,
    /// Largest log-fugacity residual of the equilibrium solve.
    pub max_ln_fugacity_residual: f64,
    /// Solver diagnostics, `code: message`.
    pub diagnostics: Vec<String>,
}

/// Result of one material stream.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StreamResult {
    /// Fingerprint of the source unit and outlet port when calculated.
    pub fingerprint: Fingerprint,
    /// Total molar flow, mol/s.
    pub molar_flow_mol_s: f64,
    /// Total mass flow, kg/s.
    pub mass_flow_kg_s: f64,
    /// Equilibrium state; `None` for an absent separator phase (zero flow).
    pub state: Option<StateResult>,
    /// Warnings derived from state validity and solver diagnostics.
    pub warnings: Vec<String>,
}

/// Result of one unit operation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnitResult {
    /// Fingerprint of the unit inputs when calculated.
    pub fingerprint: Fingerprint,
    /// External duty, W, positive into the unit. It enters the energy
    /// balance and is zero for heat exchangers.
    pub duty_w: f64,
    /// Heat exchangers only: heat transferred internally from the hot side
    /// (inlet port 0) to the cold side (inlet port 1), W. A solved exchanger
    /// never reports a negative value. It is not an external duty and is not
    /// part of `duty_w`. `None` for other units, and for exchanger results
    /// read from format version 1 or 2 files, which did not store it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transferred_duty_w: Option<f64>,
    /// Material and enthalpy crossing the modeled flowsheet boundary at this unit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub material_exchange: Option<MaterialExchangeResult>,
    /// Shortcut column design estimates (#185, added behind format version 4).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub shortcut_column: Option<ShortcutColumnEstimate>,
    /// Engine warnings (for example a single-phase separator).
    pub warnings: Vec<String>,
    /// Largest scaled component balance residual (internal units only).
    pub max_scaled_component_residual: Option<f64>,
    /// Scaled energy balance residual (internal units only).
    pub scaled_energy_residual: Option<f64>,
}

/// Fenske–Underwood–Gilliland estimates of a shortcut column (#185). Stage
/// counts include the reboiler and exclude the total condenser.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShortcutColumnEstimate {
    pub relative_volatility: f64,
    pub feed_condition_q: f64,
    pub minimum_stages: f64,
    pub underwood_root: f64,
    pub minimum_reflux_ratio: f64,
    pub reflux_ratio: f64,
    pub gilliland_x: f64,
    pub theoretical_stages: f64,
    pub rectifying_stages: f64,
    /// W, positive out of the column.
    pub condenser_duty_w: f64,
    /// W, positive into the column.
    pub reboiler_duty_w: f64,
    pub top_vapor_temperature_k: f64,
    pub volatility_iterations: usize,
}

impl ShortcutColumnEstimate {
    fn numbers(&self) -> [f64; 12] {
        [
            self.relative_volatility,
            self.feed_condition_q,
            self.minimum_stages,
            self.underwood_root,
            self.minimum_reflux_ratio,
            self.reflux_ratio,
            self.gilliland_x,
            self.theoretical_stages,
            self.rectifying_stages,
            self.condenser_duty_w,
            self.reboiler_duty_w,
            self.top_vapor_temperature_k,
        ]
    }
}

impl From<&flowsheet_core::units::columns::ShortcutColumnResult> for ShortcutColumnEstimate {
    fn from(r: &flowsheet_core::units::columns::ShortcutColumnResult) -> Self {
        Self {
            relative_volatility: r.relative_volatility,
            feed_condition_q: r.feed_condition_q,
            minimum_stages: r.minimum_stages,
            underwood_root: r.underwood_root,
            minimum_reflux_ratio: r.minimum_reflux_ratio,
            reflux_ratio: r.reflux_ratio,
            gilliland_x: r.gilliland_x,
            theoretical_stages: r.theoretical_stages,
            rectifying_stages: r.rectifying_stages,
            condenser_duty_w: r.condenser_duty_w,
            reboiler_duty_w: r.reboiler_duty_w,
            top_vapor_temperature_k: r.top_vapor_temperature_k,
            volatility_iterations: r.volatility_iterations,
        }
    }
}

/// External material contribution, positive flow into the unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialExchangeResult {
    pub molar_flow_mol_s: f64,
    /// Signed molar flow of each component in package order.
    pub component_flows_mol_s: Vec<f64>,
    /// Enthalpy carried by the signed flow, on the package reference.
    pub enthalpy_flow_w: f64,
}

/// Result of one energy stream (#171, format version 4).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergyStreamResult {
    /// Fingerprint of both connected units when calculated.
    pub fingerprint: Fingerprint,
    /// Duty in W, positive in the stream direction (into the destination
    /// unit, out of the source unit).
    pub duty_w: f64,
}

/// Closure evidence for a successful cyclic network calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecycleResult {
    pub iterations: usize,
    pub max_scaled_residual: f64,
    /// Stream identifier and scaled closure residual, in identifier order.
    pub tears: Vec<(u64, f64)>,
}

/// Snapshot of one successful flowsheet calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResultSet {
    /// Closure evidence for cyclic networks; absent for older or acyclic results.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recycle: Option<RecycleResult>,
    /// Engine, model and data identity that produced the values.
    pub provenance: Provenance,
    /// Completion time, seconds since the Unix epoch (UTC).
    pub solved_at_unix_s: u64,
    /// Unit identifiers in calculation order.
    #[serde(with = "ids::list")]
    pub order: Vec<UnitId>,
    /// Per-unit results.
    #[serde(with = "ids::map")]
    pub units: BTreeMap<UnitId, UnitResult>,
    /// Per-stream results.
    #[serde(with = "ids::map")]
    pub streams: BTreeMap<StreamId, StreamResult>,
    /// Per-energy-stream results (format version 4); absent in older files.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty", with = "ids::map")]
    pub energy_streams: BTreeMap<StreamId, EnergyStreamResult>,
}

impl ResultSet {
    /// Converts a successful calculation. `fingerprints` must describe the
    /// inputs that were solved, not later edits.
    pub fn from_calculation(
        calc: &Calculation,
        fs: &Flowsheet,
        fingerprints: &Fingerprints,
        provenance: Provenance,
        solved_at_unix_s: u64,
    ) -> Self {
        let units = calc
            .units
            .iter()
            .filter_map(|(id, u)| {
                Some((
                    *id,
                    UnitResult {
                        fingerprint: *fingerprints.units.get(id)?,
                        duty_w: u.duty_w,
                        transferred_duty_w: u.transferred_duty_w,
                        material_exchange: u.material_exchange.as_ref().map(|x| {
                            MaterialExchangeResult {
                                molar_flow_mol_s: x.molar_flow_mol_s,
                                component_flows_mol_s: x.component_flows_mol_s.clone(),
                                enthalpy_flow_w: x.enthalpy_flow_w,
                            }
                        }),
                        shortcut_column: u.shortcut_column.as_ref().map(Into::into),
                        warnings: u.warnings.clone(),
                        max_scaled_component_residual: u
                            .balance
                            .as_ref()
                            .map(|b| b.max_scaled_component_residual),
                        scaled_energy_residual: u
                            .balance
                            .as_ref()
                            .map(|b| b.scaled_energy_residual),
                    },
                ))
            })
            .collect();
        let streams = calc
            .streams
            .iter()
            .filter_map(|(id, s)| Some((*id, stream_result(s, fingerprints.stream(fs, *id)?))))
            .collect();
        let energy_streams = calc
            .energy_streams
            .iter()
            .filter_map(|(id, e)| {
                Some((
                    *id,
                    EnergyStreamResult {
                        fingerprint: fingerprints.stream(fs, *id)?,
                        duty_w: e.duty_w,
                    },
                ))
            })
            .collect();
        Self {
            recycle: calc.recycle.as_ref().map(|r| RecycleResult {
                iterations: r.iterations,
                max_scaled_residual: r.max_scaled_residual,
                tears: r.tear_residuals.iter().map(|(id, r)| (id.0, *r)).collect(),
            }),
            provenance,
            solved_at_unix_s,
            order: calc.order.clone(),
            units,
            streams,
            energy_streams,
        }
    }

    /// Currency of a unit result against the current fingerprints.
    pub fn unit_currency(&self, id: UnitId, current: &Fingerprints) -> Currency {
        match self.units.get(&id) {
            None => Currency::Missing,
            Some(r) if current.units.get(&id) == Some(&r.fingerprint) => Currency::Current,
            Some(_) => Currency::Stale,
        }
    }

    /// Currency of a stream result against the current fingerprints. A
    /// stream that no longer exists has no result.
    pub fn stream_currency(
        &self,
        fs: &Flowsheet,
        id: StreamId,
        current: &Fingerprints,
    ) -> Currency {
        let stored = self
            .streams
            .get(&id)
            .map(|r| r.fingerprint)
            .or_else(|| self.energy_streams.get(&id).map(|r| r.fingerprint));
        match stored {
            _ if fs.stream(id).is_none() => Currency::Missing,
            None => Currency::Missing,
            Some(fp) if current.stream(fs, id) == Some(fp) => Currency::Current,
            Some(_) => Currency::Stale,
        }
    }

    /// Describes the first stored number that is NaN or infinite, if any.
    /// JSON cannot represent such numbers, so they must not be saved.
    pub(crate) fn first_non_finite(&self) -> Option<String> {
        let bad = |v: f64| !v.is_finite();
        let opt = |v: Option<f64>| v.is_some_and(bad);
        if let Some(r) = &self.recycle
            && (bad(r.max_scaled_residual) || r.tears.iter().any(|(_, v)| bad(*v)))
        {
            return Some("recycle closure result".into());
        }
        for (id, u) in &self.units {
            if bad(u.duty_w)
                || opt(u.transferred_duty_w)
                || u.material_exchange.as_ref().is_some_and(|x| {
                    bad(x.molar_flow_mol_s)
                        || bad(x.enthalpy_flow_w)
                        || x.component_flows_mol_s.iter().copied().any(bad)
                })
                || u.shortcut_column
                    .as_ref()
                    .is_some_and(|r| r.numbers().into_iter().any(bad))
                || opt(u.max_scaled_component_residual)
                || opt(u.scaled_energy_residual)
            {
                return Some(format!("result of unit {id}"));
            }
        }
        for (id, e) in &self.energy_streams {
            if bad(e.duty_w) {
                return Some(format!("result of energy stream {id}"));
            }
        }
        for (id, s) in &self.streams {
            let state_bad = s.state.as_ref().is_some_and(|t| {
                [
                    t.temperature_k,
                    t.pressure_pa,
                    t.enthalpy_j_mol,
                    t.entropy_j_mol_k,
                    t.max_ln_fugacity_residual,
                ]
                .into_iter()
                .any(bad)
                    || opt(t.vapor_fraction)
                    || t.mole_fractions.iter().any(|v| bad(*v))
                    || t.phases.iter().any(|p| {
                        [
                            p.fraction,
                            p.density_kg_m3,
                            p.compressibility,
                            p.enthalpy_j_mol,
                            p.entropy_j_mol_k,
                        ]
                        .into_iter()
                        .any(bad)
                            || p.mole_fractions.iter().any(|v| bad(*v))
                    })
            });
            if bad(s.molar_flow_mol_s) || bad(s.mass_flow_kg_s) || state_bad {
                return Some(format!("result of stream {id}"));
            }
        }
        None
    }

    /// Describes the first stored composition or material exchange with a
    /// component count different from the package, if any.
    pub(crate) fn first_component_mismatch(&self, components: usize) -> Option<String> {
        for (id, u) in &self.units {
            if let Some(x) = &u.material_exchange
                && x.component_flows_mol_s.len() != components
            {
                return Some(format!(
                    "stored material exchange of unit {id} has {} flows for {components} components",
                    x.component_flows_mol_s.len()
                ));
            }
        }
        for (id, s) in &self.streams {
            let Some(state) = &s.state else { continue };
            if state.mole_fractions.len() != components {
                return Some(format!(
                    "stored result of stream {id} has {} mole fractions for {components} components",
                    state.mole_fractions.len()
                ));
            }
            if let Some(p) = state
                .phases
                .iter()
                .find(|p| p.mole_fractions.len() != components)
            {
                return Some(format!(
                    "stored {} phase result of stream {id} has {} mole fractions for {components} components",
                    p.label,
                    p.mole_fractions.len()
                ));
            }
        }
        None
    }

    /// Removes entries for units and streams that no longer exist.
    pub fn retain_existing(&mut self, fs: &Flowsheet) {
        self.units.retain(|id, _| fs.unit(*id).is_some());
        self.streams.retain(|id, _| fs.stream(*id).is_some());
        self.energy_streams.retain(|id, _| fs.stream(*id).is_some());
        self.order.retain(|id| fs.unit(*id).is_some());
    }
}

fn stream_result(s: &MaterialStream, fingerprint: Fingerprint) -> StreamResult {
    let mut warnings = Vec::new();
    let state = s.state.as_ref().map(|a| {
        let label = |p: Phase| match p {
            Phase::Vapor => "vapor",
            Phase::Liquid => "liquid",
            Phase::Supercritical => "supercritical",
        };
        let phases: Vec<PhaseResult> = a
            .phases
            .iter()
            .filter(|p| p.phase_fraction > 0.0)
            .map(|p| {
                let validity = match p.state.validity {
                    ValidityStatus::WithinDeclaredDomain => "within declared domain",
                    ValidityStatus::Extrapolated => "extrapolated",
                    ValidityStatus::Unassessed => "unassessed",
                };
                if p.state.validity == ValidityStatus::Extrapolated {
                    warnings.push(format!(
                        "{} phase state is outside the declared model domain (extrapolated).",
                        label(p.phase)
                    ));
                }
                PhaseResult {
                    label: label(p.phase).into(),
                    fraction: p.phase_fraction,
                    mole_fractions: p.composition.clone(),
                    density_kg_m3: p.state.mass_density,
                    compressibility: p.state.compressibility,
                    enthalpy_j_mol: p.state.enthalpy,
                    entropy_j_mol_k: p.state.entropy,
                    validity: validity.into(),
                }
            })
            .collect();
        let vapor_fraction =
            (!a.phases.iter().any(|p| p.phase == Phase::Supercritical)).then(|| {
                a.phases
                    .iter()
                    .filter(|p| p.phase == Phase::Vapor)
                    .fold(0.0, |sum, p| sum + p.phase_fraction)
            });
        for d in &a.diagnostics {
            if d.severity == DiagnosticSeverity::Warning {
                warnings.push(format!("{}: {}", d.code, d.message));
            }
        }
        StateResult {
            temperature_k: a.temperature_k,
            pressure_pa: a.pressure_pa,
            mole_fractions: a.mole_fractions.clone(),
            vapor_fraction,
            enthalpy_j_mol: a.enthalpy(),
            entropy_j_mol_k: a.entropy(),
            phases,
            iterations: a.iterations,
            max_ln_fugacity_residual: a.convergence.max_ln_fugacity_residual,
            diagnostics: a
                .diagnostics
                .iter()
                .map(|d| format!("{}: {}", d.code, d.message))
                .collect(),
        }
    });
    StreamResult {
        fingerprint,
        molar_flow_mol_s: s.molar_flow_mol_s,
        mass_flow_kg_s: s.mass_flow_kg_s(),
        state,
        warnings,
    }
}

/// Serialization of engine identifiers, which carry no serde support.
pub(crate) mod ids {
    use std::collections::BTreeMap;

    use flowsheet_core::{StreamId, UnitId};
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub(crate) trait Id: Copy + Ord {
        fn raw(self) -> u64;
        fn from_raw(raw: u64) -> Self;
    }
    impl Id for UnitId {
        fn raw(self) -> u64 {
            self.0
        }
        fn from_raw(raw: u64) -> Self {
            UnitId(raw)
        }
    }
    impl Id for StreamId {
        fn raw(self) -> u64 {
            self.0
        }
        fn from_raw(raw: u64) -> Self {
            StreamId(raw)
        }
    }

    pub(crate) mod map {
        use super::*;

        pub(crate) fn serialize<S: Serializer, K: Id, V: Serialize>(
            m: &BTreeMap<K, V>,
            s: S,
        ) -> Result<S::Ok, S::Error> {
            let raw: BTreeMap<u64, &V> = m.iter().map(|(k, v)| (k.raw(), v)).collect();
            raw.serialize(s)
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>, K: Id, V: Deserialize<'de>>(
            d: D,
        ) -> Result<BTreeMap<K, V>, D::Error> {
            let raw = BTreeMap::<u64, V>::deserialize(d)?;
            Ok(raw.into_iter().map(|(k, v)| (K::from_raw(k), v)).collect())
        }
    }

    pub(crate) mod list {
        use super::*;

        pub(crate) fn serialize<S: Serializer, K: Id>(v: &[K], s: S) -> Result<S::Ok, S::Error> {
            v.iter().map(|k| k.raw()).collect::<Vec<_>>().serialize(s)
        }

        pub(crate) fn deserialize<'de, D: Deserializer<'de>, K: Id>(
            d: D,
        ) -> Result<Vec<K>, D::Error> {
            Ok(Vec::<u64>::deserialize(d)?
                .into_iter()
                .map(K::from_raw)
                .collect())
        }
    }
}
