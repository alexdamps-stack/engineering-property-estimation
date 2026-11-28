//! Product TBP and ASTM D86 curves (#190) against the independent
//! standard-library references of
//! `validation/columns/generate-petroleum-curve-references.py` and the
//! printed example E3.1 of Fahim et al. (2010) read backwards.

use serde_json::Value;
use thermo_data::{CurveComponent, D86_TO_TBP_API_CONSTANTS, TbpCurve, tbp_to_d86};

const CURVES: &str = include_str!("../../../validation/columns/petroleum-curves-v1.json");
const COLUMNS: &str = include_str!("../../../validation/columns/petroleum-v1.json");

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn points(v: &Value) -> Vec<(f64, f64)> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|p| (number(&p[0]), number(&p[1])))
        .collect()
}

fn tolerance(fixture: &Value, field: &str) -> f64 {
    number(&fixture["tolerances"][field][0])
}

#[test]
fn api_constants_match_the_reference_transcription() {
    let fixture: Value = serde_json::from_str(CURVES).unwrap();
    for (rust, python) in D86_TO_TBP_API_CONSTANTS
        .iter()
        .zip(fixture["api_d86_tbp_constants"].as_array().unwrap())
    {
        assert_eq!(rust.0, number(&python[0]));
        assert_eq!(rust.1, number(&python[1]));
        assert_eq!(rust.2, number(&python[2]));
    }
}

#[test]
fn published_e3_1_d86_is_recovered_from_its_tbp() {
    let fixture: Value = serde_json::from_str(CURVES).unwrap();
    let e = &fixture["published_example"];
    let atol = number(&e["atol_degC"]);
    let tbp: Vec<(f64, f64)> = D86_TO_TBP_API_CONSTANTS
        .iter()
        .zip(e["tbp_degC"].as_array().unwrap())
        .map(|((v, _, _), t)| (*v, number(t) + 273.15))
        .collect();
    let d86 = tbp_to_d86(&tbp).unwrap();
    for ((_, actual), printed) in d86.iter().zip(e["printed_d86_degC"].as_array().unwrap()) {
        let printed = number(printed) + 273.15;
        assert!(
            (actual - printed).abs() <= atol,
            "{actual} K differs from the printed {printed} K"
        );
    }
}

#[test]
fn synthetic_curves_match_the_independent_construction() {
    let fixture: Value = serde_json::from_str(CURVES).unwrap();
    let tol = tolerance(&fixture, "curve_temperature_K");
    for case in fixture["synthetic_curves"].as_array().unwrap() {
        let components: Vec<CurveComponent> = points(&case["components"])
            .into_iter()
            .map(|(tb, v)| CurveComponent {
                normal_boiling_point_k: tb,
                liquid_volume: v,
            })
            .collect();
        let curve = TbpCurve::from_components(&components).unwrap();
        for (v, t) in points(&case["tbp_K"]) {
            let actual = curve.temperature_at(v).unwrap();
            assert!(
                (actual - t).abs() <= tol,
                "{} TBP {v}: {actual} {t}",
                case["id"]
            );
        }
        for (v, t) in points(&case["d86_K"]) {
            let actual = curve.d86_temperature_at(v).unwrap();
            assert!(
                (actual - t).abs() <= tol,
                "{} D86 {v}: {actual} {t}",
                case["id"]
            );
        }
    }
}

#[test]
fn reference_products_match_at_their_reference_compositions() {
    let fixture: Value = serde_json::from_str(CURVES).unwrap();
    let columns: Value = serde_json::from_str(COLUMNS).unwrap();
    let tol = tolerance(&fixture, "curve_temperature_K");
    let products = fixture["products"].as_array().unwrap();
    assert_eq!(products.len(), 13);
    for product in products {
        let case = columns["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["id"] == product["case"])
            .unwrap();
        let label = product["product"].as_str().unwrap();
        let index = match label {
            "overhead" => 0,
            "bottoms" => 1,
            draw => 1 + draw["side draw ".len()..].parse::<usize>().unwrap(),
        };
        let z = case["outlets"][index]["z"].as_array().unwrap();
        let mut excluded = 0.0;
        let components: Vec<CurveComponent> = case["components"]
            .as_array()
            .unwrap()
            .iter()
            .zip(z)
            .filter_map(|(c, x)| match c.get("standard_liquid_volume_m3_per_mol") {
                Some(v) => Some(CurveComponent {
                    normal_boiling_point_k: number(&c["normal_boiling_point_K"]),
                    liquid_volume: number(x).max(0.0) * number(v),
                }),
                None => {
                    excluded += number(x);
                    None
                }
            })
            .collect();
        assert!(
            (excluded - number(&product["excluded_mole_fraction"])).abs()
                <= tolerance(&fixture, "excluded_mole_fraction")
        );
        let curve = TbpCurve::from_components(&components).unwrap();
        for (v, t) in points(&product["tbp_K"]) {
            assert!(
                (curve.temperature_at(v).unwrap() - t).abs() <= tol,
                "{label} TBP {v}"
            );
        }
        for (v, t) in points(&product["d86_K"]) {
            assert!(
                (curve.d86_temperature_at(v).unwrap() - t).abs() <= tol,
                "{label} D86 {v}"
            );
        }
    }
}
