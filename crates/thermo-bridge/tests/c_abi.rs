//! Exercises the exported C functions exactly as a foreign host does, and
//! checks that the bridge adds no numerical layer over the engine.
#![allow(unsafe_code)]

use serde_json::{Value, json};
use std::ffi::{CStr, CString};
use thermocore_bridge::ffi::{
    tc_bridge_abi_version, tc_bridge_calculate, tc_bridge_free, tc_bridge_info,
};

fn take(ptr: *mut std::ffi::c_char) -> Value {
    assert!(!ptr.is_null());
    // SAFETY: returned by the bridge, NUL-terminated, freed once below.
    let text = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap().to_owned();
    unsafe { tc_bridge_free(ptr) };
    serde_json::from_str(&text).unwrap()
}

fn call(request: &[u8]) -> Value {
    let c = CString::new(request).unwrap();
    take(unsafe { tc_bridge_calculate(c.as_ptr()) })
}

fn request() -> Value {
    serde_json::from_str(include_str!(
        "../../../integrations/dwsim/requests/natural-gas-separator-tp.json"
    ))
    .unwrap()
}

#[test]
fn abi_version_info_and_free_contract() {
    assert_eq!(tc_bridge_abi_version(), 1);
    let info = take(tc_bridge_info());
    assert_eq!(info["schema"], "thermocore.bridge.info.v1");
    assert_eq!(info["abi_version"], 1);
    unsafe { tc_bridge_free(std::ptr::null_mut()) };
}

#[test]
fn null_and_non_utf8_requests_return_typed_errors() {
    let r = take(unsafe { tc_bridge_calculate(std::ptr::null()) });
    assert_eq!(r["status"], "error");
    assert_eq!(r["error"]["kind"], "invalid-request");
    let r = call(&[0xff, 0xfe, b'{']);
    assert_eq!(r["error"]["kind"], "invalid-request");
}

#[test]
fn c_abi_result_is_identical_to_direct_engine_call() {
    use thermo_data::ComponentDatabase;
    use thermo_models::{CubicModel, CubicVariant};
    use thermo_solvers::{FlashOptions, tp_flash};
    use thermo_types::{
        ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure, Temperature,
    };
    let req = request();
    let r = call(req.to_string().as_bytes());
    assert_eq!(r["status"], "ok", "{r}");
    let pkg = &req["package"];
    let ids: Vec<String> = serde_json::from_value(pkg["components"].clone()).unwrap();
    let z: Vec<f64> = serde_json::from_value(req["state"]["mole_fractions"].clone()).unwrap();
    let kij: Vec<Vec<f64>> = serde_json::from_value(pkg["kij"].clone()).unwrap();
    let comp = Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(&z)
            .map(|(i, x)| (ComponentId::new(i).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &comp,
        CubicVariant::PengRobinsonBostonMathias,
        kij,
        pkg["interaction_source"].as_str().unwrap(),
    )
    .unwrap();
    let direct = tp_flash(
        &model,
        Temperature::new(req["state"]["temperature_K"].as_f64().unwrap()).unwrap(),
        Pressure::new(req["state"]["pressure_Pa"].as_f64().unwrap()).unwrap(),
        FlashOptions::default(),
    )
    .unwrap();
    let res = &r["result"];
    assert_eq!(res["phase_count"], json!(direct.phases.len()));
    assert_eq!(
        res["molar_vapor_fraction"].as_f64(),
        direct.vapor_fraction(),
        "bit-identical through JSON (shortest round-trip formatting)"
    );
    for (p, d) in res["phases"].as_array().unwrap().iter().zip(&direct.phases) {
        let x: Vec<f64> = serde_json::from_value(p["mole_fractions"].clone()).unwrap();
        assert_eq!(x, d.composition);
        assert_eq!(p["Z"].as_f64().unwrap(), d.state.compressibility);
    }
}
