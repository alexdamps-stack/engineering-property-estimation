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
    /// Duty, W, positive into the unit.
    pub duty_w: f64,
    /// Engine warnings (for example a single-phase separator).
    pub warnings: Vec<String>,
    /// Largest scaled component balance residual (internal units only).
    pub max_scaled_component_residual: Option<f64>,
    /// Scaled energy balance residual (internal units only).
    pub scaled_energy_residual: Option<f64>,
}

/// Snapshot of one successful flowsheet calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ResultSet {
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
        Self {
            provenance,
            solved_at_unix_s,
            order: calc.order.clone(),
            units,
            streams,
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
        match self.streams.get(&id) {
            _ if fs.stream(id).is_none() => Currency::Missing,
            None => Currency::Missing,
            Some(r) if current.stream(fs, id) == Some(r.fingerprint) => Currency::Current,
            Some(_) => Currency::Stale,
        }
    }

    /// Removes entries for units and streams that no longer exist.
    pub fn retain_existing(&mut self, fs: &Flowsheet) {
        self.units.retain(|id, _| fs.unit(*id).is_some());
        self.streams.retain(|id, _| fs.stream(*id).is_some());
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
