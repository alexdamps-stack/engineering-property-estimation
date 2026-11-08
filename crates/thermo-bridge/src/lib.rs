//! A12 external-integration boundary.
//!
//! A deliberately narrow, versioned interface for foreign hosts (the DWSIM
//! .NET client in `integrations/dwsim`): UTF-8 JSON in, UTF-8 JSON out, over a
//! four-function C ABI. The scope is the qualified cubic package (`PR`,
//! `PR-BM`) TP and PH flash through [`property_packages::CubicPackage`]; no
//! solver or model code lives here.
//!
//! Contract (schema `thermocore.bridge.*.v1`, ABI version [`ABI_VERSION`]):
//! - SI only, with units in every field name (`temperature_K`, `pressure_Pa`
//!   absolute, `enthalpy_J_per_mol` on the reported caloric reference).
//! - Unknown fields are rejected; there are no defaults for physics inputs
//!   (component order, complete `kij` and its source are explicit).
//! - Every call returns a response object with `status` `ok` or `error`.
//!   Errors carry a stable `kind` and a message; a failed calculation never
//!   becomes a partial success. Panics are caught at the ABI and reported as
//!   `internal-panic`.
//! - Successful results retain model, data and caloric-reference provenance,
//!   validity status, convergence residuals and solver diagnostics.

use property_packages::{CubicPackage, EquilibriumState};
use serde::Deserialize;
use serde_json::{Value, json};
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, DiagnosticSeverity,
    NormalizationPolicy, Phase, ValidityStatus,
};

pub mod ffi;

/// Incremented only for incompatible changes of the C functions or the JSON schema.
pub const ABI_VERSION: u32 = 1;
pub const REQUEST_SCHEMA: &str = "thermocore.bridge.request.v1";
pub const RESPONSE_SCHEMA: &str = "thermocore.bridge.response.v1";
pub const INFO_SCHEMA: &str = "thermocore.bridge.info.v1";
pub const BRIDGE_VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    schema: String,
    #[serde(default)]
    request_id: Option<String>,
    calculation: String,
    package: PackageSpec,
    state: StateSpec,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PackageSpec {
    model: String,
    components: Vec<String>,
    kij: Vec<Vec<f64>>,
    interaction_source: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StateSpec {
    mole_fractions: Vec<f64>,
    #[serde(default, rename = "temperature_K")]
    temperature_k: Option<f64>,
    #[serde(default, rename = "pressure_Pa")]
    pressure_pa: Option<f64>,
    #[serde(default)]
    enthalpy_j_per_mol: Option<f64>,
}

/// Stable error categories of the v1 response.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeError {
    pub kind: &'static str,
    pub message: String,
}

impl BridgeError {
    fn request(message: impl Into<String>) -> Self {
        Self {
            kind: "invalid-request",
            message: message.into(),
        }
    }
}

impl From<CalculationError> for BridgeError {
    fn from(e: CalculationError) -> Self {
        let (kind, message) = match e {
            CalculationError::InvalidInput(i) => ("invalid-input", format!("{i:?}")),
            CalculationError::MissingParameter(m) => ("missing-parameter", m),
            CalculationError::UnsupportedCapability(m) => ("unsupported-capability", m),
            CalculationError::OutsideValidityDomain(m) => ("outside-validity-domain", m),
            CalculationError::NonConvergence(m) => ("non-convergence", m),
            CalculationError::NumericalFailure(m) => ("numerical-failure", m),
        };
        Self { kind, message }
    }
}

impl From<thermo_types::InputError> for BridgeError {
    fn from(e: thermo_types::InputError) -> Self {
        CalculationError::from(e).into()
    }
}

fn envelope(request_id: Option<&str>) -> Value {
    json!({
        "schema": RESPONSE_SCHEMA,
        "abi_version": ABI_VERSION,
        "bridge_version": BRIDGE_VERSION,
        "request_id": request_id,
    })
}

/// Error response in the v1 envelope.
pub fn error_response(request_id: Option<&str>, error: &BridgeError) -> Value {
    let mut v = envelope(request_id);
    v["status"] = json!("error");
    v["error"] = json!({"kind": error.kind, "message": error.message});
    v
}

/// Bridge identity, capabilities and the embedded component catalogue (with
/// CAS numbers so hosts can map their own compound names explicitly).
pub fn info() -> Value {
    let catalogue = ComponentDatabase::embedded().map(|db| {
        let components: Vec<Value> = db
            .components()
            .iter()
            .map(|c| {
                json!({
                    "id": c.id(),
                    "cas": c.cas(),
                    "molar_mass_kg_per_mol": c.molar_mass().value(),
                })
            })
            .collect();
        let m = db.metadata();
        (
            components,
            json!({"id": m.identifier, "version": m.version, "source": m.source}),
        )
    });
    let (components, data) = match catalogue {
        Ok(c) => c,
        Err(e) => {
            let e = BridgeError::from(e);
            (
                vec![],
                json!({"error": {"kind": e.kind, "message": e.message}}),
            )
        }
    };
    json!({
        "schema": INFO_SCHEMA,
        "abi_version": ABI_VERSION,
        "bridge_version": BRIDGE_VERSION,
        "request_schema": REQUEST_SCHEMA,
        "response_schema": RESPONSE_SCHEMA,
        "models": ["PR", "PR-BM"],
        "calculations": ["tp-flash", "ph-flash"],
        "units": {"temperature": "K", "pressure": "Pa (absolute)", "molar_enthalpy": "J/mol",
                  "molar_entropy": "J/(mol K)", "composition": "mole fraction"},
        "temperature_envelope_K": [200.0, 600.0],
        "max_phases": 2,
        "data": data,
        "components": components,
    })
}

/// Evaluates one v1 request. Always returns a complete v1 response.
pub fn calculate(request: &str) -> Value {
    let parsed: Result<Request, _> = serde_json::from_str(request);
    let request = match parsed {
        Ok(r) => r,
        Err(e) => {
            // Echo request_id when the document is at least valid JSON.
            let id = serde_json::from_str::<Value>(request)
                .ok()
                .and_then(|v| v["request_id"].as_str().map(str::to_owned));
            return error_response(
                id.as_deref(),
                &BridgeError::request(format!("malformed v1 request: {e}")),
            );
        }
    };
    let id = request.request_id.clone();
    match evaluate(request) {
        Ok(result) => {
            let mut v = envelope(id.as_deref());
            v["status"] = json!("ok");
            v["result"] = result;
            v
        }
        Err(e) => error_response(id.as_deref(), &e),
    }
}

fn evaluate(r: Request) -> Result<Value, BridgeError> {
    if r.schema != REQUEST_SCHEMA {
        return Err(BridgeError::request(format!(
            "unsupported request schema {:?}; this bridge accepts {REQUEST_SCHEMA}",
            r.schema
        )));
    }
    let variant = match r.package.model.as_str() {
        "PR" => CubicVariant::PengRobinson,
        "PR-BM" => CubicVariant::PengRobinsonBostonMathias,
        other => {
            return Err(BridgeError {
                kind: "unsupported-capability",
                message: format!("model {other:?} is not exposed by bridge v1; use PR or PR-BM"),
            });
        }
    };
    let n = r.package.components.len();
    if r.state.mole_fractions.len() != n {
        return Err(BridgeError::request(format!(
            "{} mole fractions for {n} components",
            r.state.mole_fractions.len()
        )));
    }
    let entries = r
        .package
        .components
        .iter()
        .zip(&r.state.mole_fractions)
        .map(|(id, x)| Ok((ComponentId::new(id)?, *x)))
        .collect::<Result<Vec<_>, thermo_types::InputError>>()?;
    let composition = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::RequireUnitSum,
    )?;
    let db = ComponentDatabase::embedded()?;
    let model = CubicModel::new(
        &db,
        &composition,
        variant,
        r.package.kij,
        r.package.interaction_source,
    )?;
    let molar_masses: Vec<f64> = r
        .package
        .components
        .iter()
        .map(|id| db.get(id).map(|c| c.molar_mass().value()))
        .collect::<Result<_, _>>()?;
    let fractions: Vec<f64> = composition.entries().iter().map(|(_, x)| *x).collect();
    let package = CubicPackage::new(model, PhFlashOptions::default());
    let s = &r.state;
    let (state, target_h) = match r.calculation.as_str() {
        "tp-flash" => {
            if s.enthalpy_j_per_mol.is_some() {
                return Err(BridgeError::request(
                    "enthalpy_j_per_mol is an output of tp-flash; omit it",
                ));
            }
            let t = s
                .temperature_k
                .ok_or_else(|| BridgeError::request("tp-flash requires state.temperature_K"))?;
            let p = s
                .pressure_pa
                .ok_or_else(|| BridgeError::request("tp-flash requires state.pressure_Pa"))?;
            (package.tp(&fractions, t, p)?, None)
        }
        "ph-flash" => {
            if s.temperature_k.is_some() {
                return Err(BridgeError::request(
                    "temperature_K is an output of ph-flash; omit it",
                ));
            }
            let p = s
                .pressure_pa
                .ok_or_else(|| BridgeError::request("ph-flash requires state.pressure_Pa"))?;
            let h = s.enthalpy_j_per_mol.ok_or_else(|| {
                BridgeError::request("ph-flash requires state.enthalpy_j_per_mol")
            })?;
            (package.ph(&fractions, p, h)?, Some(h))
        }
        other => {
            return Err(BridgeError {
                kind: "unsupported-capability",
                message: format!("calculation {other:?}; bridge v1 supports tp-flash, ph-flash"),
            });
        }
    };
    Ok(result_json(
        &r.calculation,
        &r.package.components,
        &molar_masses,
        &state,
        target_h,
    ))
}

fn label(p: Phase) -> &'static str {
    match p {
        Phase::Vapor => "vapor",
        Phase::Liquid => "liquid",
        Phase::Supercritical => "supercritical",
    }
}

fn validity(v: ValidityStatus) -> &'static str {
    match v {
        ValidityStatus::WithinDeclaredDomain => "within-declared-domain",
        ValidityStatus::Extrapolated => "extrapolated",
        ValidityStatus::Unassessed => "unassessed",
    }
}

fn vapor_fraction(state: &EquilibriumState) -> Option<f64> {
    match state.phases.as_slice() {
        [single] => match single.phase {
            Phase::Vapor => Some(1.0),
            Phase::Liquid => Some(0.0),
            Phase::Supercritical => None,
        },
        [vapor, _] => Some(vapor.phase_fraction),
        _ => None,
    }
}

fn result_json(
    calculation: &str,
    components: &[String],
    molar_masses: &[f64],
    state: &EquilibriumState,
    target_h: Option<f64>,
) -> Value {
    let phases: Vec<Value> = state
        .phases
        .iter()
        .map(|ph| {
            let s = &ph.state;
            let mw: f64 = ph
                .composition
                .iter()
                .zip(molar_masses)
                .map(|(x, m)| x * m)
                .sum();
            json!({
                "label": label(ph.phase),
                "molar_phase_fraction": ph.phase_fraction,
                "mole_fractions": ph.composition,
                "molar_mass_kg_per_mol": mw,
                "Z": s.compressibility,
                "molar_volume_m3_per_mol": s.molar_volume,
                "mass_density_kg_per_m3": s.mass_density,
                "enthalpy_J_per_mol": s.enthalpy,
                "entropy_J_per_mol_K": s.entropy,
                "cp_J_per_mol_K": s.cp,
                "cv_J_per_mol_K": s.cv,
                "speed_of_sound_m_per_s": s.speed_of_sound,
                "ln_fugacity_coefficients": s.ln_fugacity_coefficients,
                "validity": validity(s.validity),
            })
        })
        .collect();
    let first = &state.phases[0].state;
    json!({
        "calculation": calculation,
        "components": components,
        "feed_mole_fractions": state.mole_fractions,
        "temperature_K": state.temperature_k,
        "pressure_Pa": state.pressure_pa,
        "target_enthalpy_J_per_mol": target_h,
        "enthalpy_J_per_mol": state.enthalpy(),
        "entropy_J_per_mol_K": state.entropy(),
        "phase_count": state.phases.len(),
        "molar_vapor_fraction": vapor_fraction(state),
        "phases": phases,
        "convergence": {
            "iterations": state.iterations,
            "max_abs_ln_fugacity_residual": state.convergence.max_ln_fugacity_residual,
            "max_abs_component_balance_residual": state.convergence.max_component_balance_residual,
            "scaled_energy_residual": state.convergence.scaled_energy_residual,
        },
        "provenance": {
            "model": {"id": first.model.identifier, "implementation": first.model.implementation_version,
                      "formulation": first.model.formulation},
            "data": {"id": first.data.identifier, "version": first.data.version, "source": first.data.source},
            "caloric_reference": {"id": first.caloric_reference.identifier,
                                  "definition": first.caloric_reference.definition},
        },
        "validity_note": "Two-phase fluid flash inside the cubic computational envelope; physical domain unassessed unless the phase validity says otherwise; no solids, no three-phase split.",
        "diagnostics": state.diagnostics.iter().map(|d| json!({
            "code": d.code,
            "severity": match d.severity { DiagnosticSeverity::Information => "information", DiagnosticSeverity::Warning => "warning" },
            "message": d.message,
        })).collect::<Vec<_>>(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn natural_gas_tp() -> Value {
        json!({
            "schema": REQUEST_SCHEMA, "request_id": "ng-sep", "calculation": "tp-flash",
            "package": {"model": "PR-BM",
                "components": ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"],
                "kij": vec![vec![0.0; 5]; 5],
                "interaction_source": "Explicit synthetic zero-kij approximation; not fitted VLE data"},
            "state": {"mole_fractions": [0.85, 0.075, 0.04666666666666667, 0.018333333333333333, 0.01],
                      "temperature_K": 220.0, "pressure_Pa": 1_500_000.0}
        })
    }

    #[test]
    fn tp_result_carries_si_provenance_and_convergence() {
        let r = calculate(&natural_gas_tp().to_string());
        assert_eq!(r["status"], "ok", "{r}");
        assert_eq!(r["schema"], RESPONSE_SCHEMA);
        assert_eq!(r["request_id"], "ng-sep");
        let res = &r["result"];
        assert_eq!(res["phase_count"], 2);
        assert_eq!(res["temperature_K"], 220.0);
        assert_eq!(res["pressure_Pa"], 1_500_000.0);
        let beta = res["molar_vapor_fraction"].as_f64().unwrap();
        assert!(beta > 0.0 && beta < 1.0);
        assert!(
            res["convergence"]["max_abs_ln_fugacity_residual"]
                .as_f64()
                .unwrap()
                <= 1e-10
        );
        assert!(
            res["provenance"]["model"]["formulation"]
                .as_str()
                .unwrap()
                .contains("zero-kij")
        );
        assert!(res["provenance"]["caloric_reference"]["id"].is_string());
        assert!(res["provenance"]["data"]["version"].is_string());
    }

    #[test]
    fn tp_ph_round_trip_through_bridge() {
        let tp = calculate(&natural_gas_tp().to_string());
        let h = tp["result"]["enthalpy_J_per_mol"].as_f64().unwrap();
        let mut ph = natural_gas_tp();
        ph["calculation"] = json!("ph-flash");
        ph["state"].as_object_mut().unwrap().remove("temperature_K");
        ph["state"]["enthalpy_j_per_mol"] = json!(h);
        let r = calculate(&ph.to_string());
        assert_eq!(r["status"], "ok", "{r}");
        assert!((r["result"]["temperature_K"].as_f64().unwrap() - 220.0).abs() < 1e-7);
        assert!(
            r["result"]["convergence"]["scaled_energy_residual"]
                .as_f64()
                .unwrap()
                .abs()
                <= 1e-10
        );
    }

    #[test]
    fn failures_are_typed_and_never_partial() {
        type Mutation = Box<dyn Fn(&mut Value)>;
        let cases: Vec<(Mutation, &str)> = vec![
            (
                Box::new(|v| v["schema"] = json!("thermocore.bridge.request.v0")),
                "invalid-request",
            ),
            (
                Box::new(|v| v["state"]["gauge_pressure"] = json!(1.0)),
                "invalid-request",
            ),
            (
                Box::new(|v| v["state"]["pressure_Pa"] = json!(-1.0)),
                "invalid-input",
            ),
            (
                Box::new(|v| v["state"]["temperature_K"] = json!(150.0)),
                "unsupported-capability",
            ),
            (
                Box::new(|v| v["package"]["components"][4] = json!("hydrogen")),
                "missing-parameter",
            ),
            (
                Box::new(|v| v["package"]["model"] = json!("SRK")),
                "unsupported-capability",
            ),
            (
                Box::new(|v| v["package"]["kij"][0][1] = json!(0.1)),
                "unsupported-capability",
            ),
            (
                Box::new(|v| v["package"]["interaction_source"] = json!(" ")),
                "missing-parameter",
            ),
            (
                Box::new(|v| v["state"]["mole_fractions"][0] = json!(0.5)),
                "invalid-input",
            ),
            (
                Box::new(|v| v["state"]["enthalpy_j_per_mol"] = json!(0.0)),
                "invalid-request",
            ),
            (
                Box::new(|v| v["calculation"] = json!("ps-flash")),
                "unsupported-capability",
            ),
        ];
        for (mutate, kind) in cases {
            let mut v = natural_gas_tp();
            mutate(&mut v);
            let r = calculate(&v.to_string());
            assert_eq!(r["status"], "error", "{v}");
            assert_eq!(r["error"]["kind"], kind, "{r}");
            assert!(r.get("result").is_none());
            assert_eq!(r["request_id"], "ng-sep");
        }
        let r = calculate("{not json");
        assert_eq!(r["error"]["kind"], "invalid-request");
        assert!(r["request_id"].is_null());
    }

    #[test]
    fn info_lists_catalogue_with_cas() {
        let i = info();
        assert_eq!(i["abi_version"], ABI_VERSION);
        let methane = i["components"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["id"] == "methane")
            .unwrap();
        assert_eq!(methane["cas"], "74-82-8");
        assert!(i["data"]["version"].is_string());
    }
}
