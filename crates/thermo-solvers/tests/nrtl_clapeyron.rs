//! Mandatory offline A7 NRTL/gamma–phi checks against independent Clapeyron.jl
//! references. The live check supplies a fresh candidate through
//! THERMO_NRTL_REFERENCE without replacing the committed fixture.
use serde_json::Value;
use thermo_data::{AlphaTerms, ComponentDatabase, LiquidDatabase, TauTerms};
use thermo_models::{ActivityModel, NrtlModel, NrtlParameters, VaporModel};
use thermo_solvers::{
    ActivityOptions, activity_bubble_pressure, activity_bubble_temperature, activity_dew_pressure,
    activity_dew_temperature, activity_ph_flash, activity_tp_flash,
};
use thermo_types::{
    ComponentId, Composition, CompositionBasis, NormalizationPolicy, Phase, Pressure, Temperature,
};

fn fixture() -> Value {
    let input = match std::env::var_os("THERMO_NRTL_REFERENCE") {
        Some(path) => std::fs::read_to_string(path).expect("read live NRTL candidate"),
        None => include_str!("../../../validation/fixtures/nrtl-v1.json").into(),
    };
    let f: Value = serde_json::from_str(&input).expect("valid NRTL reference JSON");
    assert_eq!(f["schema"], "clapeyron-nrtl-reference-v1");
    assert_eq!(f["clapeyron_version"], "0.6.30");
    assert_eq!(f["julia_version"], "1.11.7");
    for (key, count) in [
        ("activity", 30),
        ("ternary", 9),
        ("pure", 10),
        ("saturation", 36),
        ("tp", 8),
        ("ph", 16),
    ] {
        assert_eq!(f[key].as_array().unwrap().len(), count, "{key} count");
    }
    f
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

thread_local! {
    static WORST: std::cell::RefCell<std::collections::BTreeMap<String, (f64, f64)>> =
        const { std::cell::RefCell::new(std::collections::BTreeMap::new()) };
}

/// Largest absolute deviation and budget fraction per quantity (printed with --nocapture).
fn report(test: &str) {
    WORST.with(|w| {
        for (key, (dev, used)) in w.borrow_mut().split_off("") {
            println!("{test}: {key}: max |deviation| {dev:e}, max budget fraction {used:.3e}");
        }
    });
}

fn close(f: &Value, key: &str, actual: f64, expected: f64, context: &str) {
    let b = numbers(&f["tolerances"][key]);
    let limit = b[0] + b[1] * expected.abs();
    let dev = (actual - expected).abs();
    WORST.with(|w| {
        let mut w = w.borrow_mut();
        let e = w.entry(key.into()).or_insert((0.0, 0.0));
        e.0 = e.0.max(dev);
        e.1 = e.1.max(if limit > 0.0 { dev / limit } else { 0.0 });
    });
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= limit,
        "{context}: {key}: Rust {actual:e}, Clapeyron {expected:e}, deviation {:e}, budget {limit:e}",
        (actual - expected).abs()
    );
}

fn close_all(f: &Value, key: &str, actual: &[f64], expected: &[f64], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: {key} length");
    for (a, e) in actual.iter().zip(expected) {
        close(f, key, *a, *e, context);
    }
}

fn model(z: &[f64]) -> ActivityModel {
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![
            (ComponentId::new("methanol").unwrap(), z[0]),
            (ComponentId::new("water").unwrap(), z[1]),
        ],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition,
        VaporModel::IdealGas,
    )
    .unwrap()
}

fn check_excess(f: &Value, nrtl: &NrtlModel, case: &Value, context: &str) {
    let t = case["T"].as_f64().unwrap();
    let ex = nrtl.excess(t, &numbers(&case["x"])).unwrap();
    let c = format!("{context} T={t} x={}", case["x"]);
    close_all(
        f,
        "ln_gamma",
        &ex.ln_activity_coefficients,
        &numbers(&case["ln_gamma"]),
        &c,
    );
    close_all(
        f,
        "d_ln_gamma_dT",
        &ex.d_ln_activity_coefficients_dt,
        &numbers(&case["d_ln_gamma_dT"]),
        &c,
    );
    for (key, actual) in [
        ("g_excess", ex.excess_gibbs),
        ("h_excess", ex.excess_enthalpy),
        ("s_excess", ex.excess_entropy),
        ("cp_excess", ex.excess_heat_capacity),
    ] {
        close(f, key, actual, case[key].as_f64().unwrap(), &c);
    }
}

#[test]
fn binary_methanol_water_activity_matches_clapeyron_nrtl() {
    let f = fixture();
    let m = model(&[0.5, 0.5]);
    for case in f["activity"].as_array().unwrap() {
        check_excess(&f, m.nrtl(), case, "methanol/water");
    }
    report("binary_methanol_water_activity_matches_clapeyron_nrtl");
}

#[test]
fn directed_ternary_with_all_temperature_terms_matches_clapeyron_aspen_nrtl() {
    let f = fixture();
    let p = &f["ternary_parameters"];
    let m =
        |key: &str| -> Vec<Vec<f64>> { p[key].as_array().unwrap().iter().map(numbers).collect() };
    let (t0, t1, t2, t3, a0, a1) = (m("t0"), m("t1"), m("t2"), m("t3"), m("a0"), m("a1"));
    let ids: Vec<ComponentId> = f["ternary_ids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| ComponentId::new(v.as_str().unwrap()).unwrap())
        .collect();
    let tau = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| TauTerms {
                    t0: t0[i][j],
                    t1: t1[i][j],
                    t2: t2[i][j],
                    t3: t3[i][j],
                })
                .collect()
        })
        .collect();
    let alpha = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| AlphaTerms {
                    a0: a0[i][j],
                    a1: a1[i][j],
                })
                .collect()
        })
        .collect();
    let nrtl = NrtlModel::new(
        NrtlParameters::new(
            ids,
            tau,
            alpha,
            "synthetic directed ternary for formulation checks",
        )
        .unwrap(),
    );
    for case in f["ternary"].as_array().unwrap() {
        check_excess(&f, &nrtl, case, "ternary");
    }
    report("directed_ternary_with_all_temperature_terms_matches_clapeyron_aspen_nrtl");
}

#[test]
fn pure_liquid_correlations_and_caloric_standard_states_match() {
    let f = fixture();
    for case in f["pure"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let index = usize::from(id == "water");
        let z = if index == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let m = model(&z);
        let t = Temperature::new(case["T"].as_f64().unwrap()).unwrap();
        let s = m.pure_liquid(index, t).unwrap();
        let c = format!("{id} T={}", t.value());
        close(
            &f,
            "psat",
            s.vapor_pressure,
            case["psat"].as_f64().unwrap(),
            &c,
        );
        close(
            &f,
            "dlnpsat_dT",
            s.d_ln_vapor_pressure_dt,
            case["dlnpsat_dT"].as_f64().unwrap(),
            &c,
        );
        close(
            &f,
            "d2lnpsat_dT2",
            s.d2_ln_vapor_pressure_dt2,
            case["d2lnpsat_dT2"].as_f64().unwrap(),
            &c,
        );
        close(&f, "h", s.enthalpy, case["h_liquid"].as_f64().unwrap(), &c);
        close(&f, "s", s.entropy, case["s_liquid"].as_f64().unwrap(), &c);
        close(
            &f,
            "cp",
            s.heat_capacity,
            case["cp_liquid"].as_f64().unwrap(),
            &c,
        );
        let rho = LiquidDatabase::embedded()
            .unwrap()
            .pure(id)
            .unwrap()
            .saturated_liquid_molar_density(t)
            .unwrap();
        close(
            &f,
            "liquid_molar_density",
            rho,
            case["liquid_molar_density"].as_f64().unwrap(),
            &c,
        );
    }
    report("pure_liquid_correlations_and_caloric_standard_states_match");
}

#[test]
fn bubble_and_dew_points_match_independent_modified_raoult_solves() {
    let f = fixture();
    let options = ActivityOptions::default();
    for case in f["saturation"].as_array().unwrap() {
        let kind = case["kind"].as_str().unwrap();
        let (liquid, vapor) = (numbers(&case["liquid"]), numbers(&case["vapor"]));
        let feed = if kind.starts_with("bubble") {
            &liquid
        } else {
            &vapor
        };
        let m = model(feed);
        let t = Temperature::new(case["T"].as_f64().unwrap()).unwrap();
        let p = Pressure::new(case["p"].as_f64().unwrap()).unwrap();
        let result = match kind {
            "bubble-pressure" => activity_bubble_pressure(&m, t),
            "dew-pressure" => activity_dew_pressure(&m, t, options),
            "bubble-temperature" => activity_bubble_temperature(&m, p, options),
            "dew-temperature" => activity_dew_temperature(&m, p, options),
            other => panic!("unknown saturation kind {other}"),
        }
        .unwrap();
        let c = format!("{kind} T={} p={}", t.value(), p.value());
        close(&f, "temperature", result.temperature_k, t.value(), &c);
        close(&f, "pressure", result.pressure_pa, p.value(), &c);
        close_all(&f, "composition", &result.liquid, &liquid, &c);
        close_all(&f, "composition", &result.vapor, &vapor, &c);
        assert!(result.max_ln_fugacity_residual <= 1e-10, "{c}: residual");
    }
    report("bubble_and_dew_points_match_independent_modified_raoult_solves");
}

fn phase_name(phase: Phase) -> &'static str {
    match phase {
        Phase::Vapor => "vapor",
        Phase::Liquid => "liquid",
        Phase::Supercritical => "supercritical",
    }
}

#[test]
fn tp_flash_matches_independent_split_and_caloric_assembly() {
    let f = fixture();
    for case in f["tp"].as_array().unwrap() {
        let z = numbers(&case["z"]);
        let (t, p) = (case["T"].as_f64().unwrap(), case["p"].as_f64().unwrap());
        let c = format!("TP T={t} p={p} z={z:?}");
        let flash = activity_tp_flash(
            &model(&z),
            Temperature::new(t).unwrap(),
            Pressure::new(p).unwrap(),
            ActivityOptions::default(),
        )
        .unwrap();
        close(
            &f,
            "pressure",
            flash.bubble_pressure_pa,
            case["bubble_pressure"].as_f64().unwrap(),
            &c,
        );
        close(
            &f,
            "pressure",
            flash.dew_pressure_pa,
            case["dew_pressure"].as_f64().unwrap(),
            &c,
        );
        let expected = case["phases"].as_array().unwrap();
        assert_eq!(flash.phases.len(), expected.len(), "{c}: phase count");
        for (actual, e) in flash.phases.iter().zip(expected) {
            assert_eq!(
                phase_name(actual.phase),
                e["phase"].as_str().unwrap(),
                "{c}"
            );
            close(
                &f,
                "fraction",
                actual.phase_fraction,
                e["fraction"].as_f64().unwrap(),
                &c,
            );
            close_all(
                &f,
                "composition",
                &actual.composition,
                &numbers(&e["composition"]),
                &c,
            );
            let s = &e["state"];
            close(&f, "h", actual.state.enthalpy, s["h"].as_f64().unwrap(), &c);
            close(&f, "s", actual.state.entropy, s["s"].as_f64().unwrap(), &c);
            close(&f, "cp", actual.state.cp, s["cp"].as_f64().unwrap(), &c);
            close(
                &f,
                "molar_volume",
                actual.state.molar_volume,
                s["molar_volume"].as_f64().unwrap(),
                &c,
            );
        }
        close(&f, "h", flash.enthalpy(), case["h"].as_f64().unwrap(), &c);
        assert!(flash.max_ln_fugacity_residual <= 1e-10);
        assert!(flash.max_component_balance_residual <= 1e-12);
    }
    report("tp_flash_matches_independent_split_and_caloric_assembly");
}

#[test]
fn ph_flash_recovers_independent_states_including_pure_coexistence() {
    let f = fixture();
    for case in f["ph"].as_array().unwrap() {
        let z = numbers(&case["z"]);
        let p = case["p"].as_f64().unwrap();
        let target = case["target_h"].as_f64().unwrap();
        let c = format!("PH {}", case["id"]);
        let ph = activity_ph_flash(
            &model(&z),
            Pressure::new(p).unwrap(),
            target,
            ActivityOptions::default(),
        )
        .unwrap_or_else(|e| panic!("{c}: {e}"));
        close(
            &f,
            "temperature",
            ph.temperature_k,
            case["T"].as_f64().unwrap(),
            &c,
        );
        let vapor: f64 = ph
            .phases
            .iter()
            .filter(|x| x.phase == Phase::Vapor)
            .map(|x| x.phase_fraction)
            .sum();
        close(
            &f,
            "fraction",
            vapor,
            case["vapor_fraction"].as_f64().unwrap(),
            &c,
        );
        assert!(ph.scaled_energy_residual.abs() <= 1e-10, "{c}");
    }
    report("ph_flash_recovers_independent_states_including_pure_coexistence");
}
