//! Characterisation correlations against the published worked examples of
//! Fahim, Al-Sahhaf and Elkilani (2010), transcribed in
//! `validation/pseudo-components/fahim-2010-worked-examples.json`.
use serde_json::Value;
use thermo_data::{
    PseudoComponentInput, d86_mean_average_boiling_point, d86_to_tbp,
    d86_volume_average_boiling_point, specific_gravity_from_api, specific_gravity_from_watson_k,
    watson_k,
};
use thermo_types::MOLAR_GAS_CONSTANT as R;

const EXAMPLES: &str =
    include_str!("../../../validation/pseudo-components/fahim-2010-worked-examples.json");

fn example(id: &str) -> Value {
    let all: Value = serde_json::from_str(EXAMPLES).unwrap();
    all["examples"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["id"] == id)
        .unwrap_or_else(|| panic!("example {id}"))
        .clone()
}

fn expected(example: &Value, quantity: &str) -> (f64, f64, f64) {
    let row = example["expected"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["quantity"] == quantity)
        .unwrap_or_else(|| panic!("{quantity}"));
    (
        row["value"].as_f64().unwrap(),
        row["atol"].as_f64().unwrap(),
        row["rtol"].as_f64().unwrap(),
    )
}

#[track_caller]
fn check(example: &Value, quantity: &str, actual: f64) {
    let (value, atol, rtol) = expected(example, quantity);
    assert!(actual.is_finite(), "{} {quantity}", example["id"]);
    assert!(
        (actual - value).abs() <= atol + rtol * value.abs(),
        "{} {quantity}: {actual} vs printed {value} (atol {atol}, rtol {rtol})",
        example["id"]
    );
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn d86_points(e: &Value) -> Vec<(f64, f64)> {
    numbers(&e["inputs"]["d86_volume_percent"])
        .into_iter()
        .zip(numbers(&e["inputs"]["d86_temperature_degC"]))
        .map(|(v, c)| (v, c + 273.15))
        .collect()
}

fn fahrenheit(k: f64) -> f64 {
    k * 1.8 - 459.67
}

#[test]
fn e3_1_d86_to_tbp_api_method() {
    let e = example("E3.1");
    let tbp = d86_to_tbp(&d86_points(&e)).unwrap();
    let row = &e["expected"][0];
    let atol = row["atol"].as_f64().unwrap();
    for ((_, t), printed) in tbp.iter().zip(numbers(&row["values"])) {
        let c = t - 273.15;
        assert!((c - printed).abs() <= atol, "{c} vs {printed}");
    }
}

#[test]
fn e3_2_average_boiling_points_and_watson_k() {
    let e = example("E3.2");
    let d86 = d86_points(&e);
    let vabp = d86_volume_average_boiling_point(&d86).unwrap();
    let meabp = d86_mean_average_boiling_point(&d86).unwrap();
    check(&e, "vabp_degF", fahrenheit(vabp));
    check(
        &e,
        "meabp_correction_D_degF",
        fahrenheit(vabp) - fahrenheit(meabp),
    );
    check(&e, "meabp_degF", fahrenheit(meabp));
    let sg = specific_gravity_from_api(e["inputs"]["api_gravity"].as_f64().unwrap());
    check(&e, "specific_gravity", sg);
    check(&e, "watson_k", watson_k(meabp, sg));
}

#[test]
fn e3_6_gravity_at_constant_watson_k_and_molar_mass() {
    let e = example("E3.6");
    let tb = e["inputs"]["normal_boiling_point_degC"].as_f64().unwrap() + 273.15;
    let sg = specific_gravity_from_watson_k(tb, e["inputs"]["watson_k"].as_f64().unwrap());
    check(&e, "specific_gravity", sg);
    let c = PseudoComponentInput::new("pseudo:e3-6", tb, sg)
        .characterise()
        .unwrap();
    check(&e, "molar_mass_g_per_mol", c.molar_mass_kg_per_mol * 1000.0);
}

fn e3_7_oil(e: &Value) -> PseudoComponentInput {
    let tb = e["inputs"]["normal_boiling_point_degC"].as_f64().unwrap() + 273.15;
    let sg = specific_gravity_from_api(e["inputs"]["api_gravity"].as_f64().unwrap());
    PseudoComponentInput::new("pseudo:e3-7", tb, sg)
}

#[test]
fn e3_7_gravity_watson_k_and_molar_mass() {
    let e = example("E3.7");
    let c = e3_7_oil(&e).characterise().unwrap();
    check(&e, "specific_gravity", c.input.specific_gravity);
    check(&e, "watson_k", c.watson_k);
    check(&e, "molar_mass_g_per_mol", c.molar_mass_kg_per_mol * 1000.0);
}

#[test]
fn e3_9_critical_pressure_and_acentric_factor() {
    let e = example("E3.9");
    let c = e3_7_oil(&e).characterise().unwrap();
    check(&e, "critical_pressure_bar", c.critical_pressure_pa / 1e5);
    check(&e, "acentric_factor", c.acentric_factor);
    // The printed Tc is recorded as excluded, with its reason.
    assert!(e["excluded"][0]["quantity"] == "critical_temperature_K");
}

#[test]
fn e3_11_critical_constants_and_ideal_gas_heat_capacity() {
    let e = example("E3.11");
    let c = e3_7_oil(&e).characterise().unwrap();
    check(&e, "critical_temperature_K", c.critical_temperature_k);
    check(&e, "critical_pressure_bar", c.critical_pressure_pa / 1e5);
    check(&e, "acentric_factor", c.acentric_factor);
    check(&e, "molar_mass_g_per_mol", c.molar_mass_kg_per_mol * 1000.0);
    // cp/R coefficients times R in J/(mol K) are kJ/(kmol K) per power of T.
    let [a, b, cc] = c.ideal_gas_cp_over_r.map(|x| x * R);
    check(&e, "cp_A_kJ_per_kmol_K", a);
    check(&e, "cp_B_kJ_per_kmol_K2", b);
    check(&e, "cp_C_kJ_per_kmol_K3", cc);
}

/// The excluded E3.11 enthalpy total is inconsistent with the book's own
/// printed A, B and C; this guards the recorded reason.
#[test]
fn e3_11_printed_enthalpy_total_does_not_follow_from_printed_constants() {
    let e = example("E3.11");
    let (a, _, _) = expected(&e, "cp_A_kJ_per_kmol_K");
    let (b, _, _) = expected(&e, "cp_B_kJ_per_kmol_K2");
    let (c, _, _) = expected(&e, "cp_C_kJ_per_kmol_K3");
    let (t1, t2): (f64, f64) = (298.0, 723.15);
    let dh = a * (t2 - t1) + b / 2.0 * (t2 * t2 - t1 * t1) + c / 3.0 * (t2.powi(3) - t1.powi(3));
    assert!((dh - 250_200.0).abs() < 100.0, "{dh}");
    let printed = e["excluded"][0]["printed_kJ_per_kmol"].as_f64().unwrap();
    assert!((dh - printed).abs() > 50_000.0);
}
