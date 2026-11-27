//! #226 package operations against the independent Clapeyron references of
//! `validation/fixtures/package-operations-v1.json`
//! (`validation/julia/generate_packages.jl`, no Rust calls): GERG-2008
//! bubble/dew points, and pumps, compressors and expanders under GERG-2008,
//! NRTL-IG and EOS-CG-2021 solved through the flowsheet. The live Clapeyron
//! check (`validation/julia/check.jl`) points `THERMO_PACKAGE_REFERENCE` at a
//! freshly generated candidate; otherwise the committed fixture is used.

use flowsheet_core::*;
use property_packages::{EosCgPackage, FlowsheetPackage, GergPackage, NrtlPackage, SaturationKind};
use serde_json::Value;
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, GergModel, VaporModel};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

const POWER: usize = 2;

fn fixture() -> Value {
    let path = std::env::var("THERMO_PACKAGE_REFERENCE").unwrap_or_else(|_| {
        format!(
            "{}/../../validation/fixtures/package-operations-v1.json",
            env!("CARGO_MANIFEST_DIR")
        )
    });
    let fixture: Value = serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
    assert_eq!(fixture["schema"], "clapeyron-package-operations-v1");
    fixture
}

fn number(v: &Value) -> f64 {
    v.as_f64().unwrap_or_else(|| panic!("not a number: {v}"))
}

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(number).collect()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|s| s.as_str().unwrap().to_owned())
        .collect()
}

fn check(fixture: &Value, what: &str, key: &str, actual: f64, expected: f64) {
    let t = numbers(&fixture["tolerances"][key]);
    let allowed = t[0] + t[1] * expected.abs();
    assert!(
        (actual - expected).abs() <= allowed,
        "{what}: {actual} vs reference {expected} (allowed {allowed})"
    );
}

fn composition(ids: &[String], z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn gerg(ids: &[String], z: &[f64]) -> GergPackage {
    GergPackage::new(
        GergModel::new(&composition(ids, z)).unwrap(),
        GergPackage::TEMPERATURE_RANGE_K,
    )
    .unwrap()
}

#[test]
fn gerg_bubble_and_dew_points_match_clapeyron() {
    let fixture = fixture();
    for case in fixture["gerg_saturation"].as_array().unwrap() {
        let ids = strings(&case["ids"]);
        let z = numbers(&case["z"]);
        let kind = case["kind"].as_str().unwrap();
        let (t, p) = (number(&case["T"]), number(&case["p"]));
        let package = gerg(&ids, &z);
        let saturation = match kind {
            "bubble-pressure" => package.saturation_pressure(&z, t, SaturationKind::Bubble),
            "dew-pressure" => package.saturation_pressure(&z, t, SaturationKind::Dew),
            "bubble-temperature" => package.saturation_temperature(&z, p, SaturationKind::Bubble),
            "dew-temperature" => package.saturation_temperature(&z, p, SaturationKind::Dew),
            other => panic!("unknown kind {other}"),
        };
        let what = format!("{ids:?} {kind} T={t} p={p}");
        let s = saturation.unwrap_or_else(|e| panic!("{what}: {e}"));
        check(
            &fixture,
            &format!("{what} T"),
            "temperature",
            s.temperature_k,
            t,
        );
        check(&fixture, &format!("{what} p"), "pressure", s.pressure_pa, p);
        for (a, e) in s.liquid_composition.iter().zip(numbers(&case["liquid"])) {
            check(&fixture, &format!("{what} x"), "composition", *a, e);
        }
        for (a, e) in s.vapor_composition.iter().zip(numbers(&case["vapor"])) {
            check(&fixture, &format!("{what} y"), "composition", *a, e);
        }
        check(
            &fixture,
            &format!("{what} vl"),
            "volume",
            s.liquid.molar_volume,
            number(&case["vl"]),
        );
        check(
            &fixture,
            &format!("{what} vv"),
            "volume",
            s.vapor.molar_volume,
            number(&case["vv"]),
        );
    }
}

/// Feed → machine → product with a calculated power boundary.
fn machine(
    case: &Value,
    ids: &[String],
    z: &[f64],
) -> (Flowsheet, FeedCompositions, UnitId, StreamId) {
    let inlet = &case["inlet"];
    let p_out = number(&case["outlet"]["p"]);
    let efficiency = number(&case["efficiency"]);
    let (kind, spec) = match case["kind"].as_str().unwrap() {
        "pump" => (
            UnitKind::Pump,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        "compressor" => (
            UnitKind::Compressor,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        "expander" => (
            UnitKind::Expander,
            UnitSpec::Expand {
                pressure: PressureSpec::OutletPressure(Some(p_out)),
                efficiency,
            },
        ),
        other => panic!("unknown machine {other}"),
    };
    let mut flowsheet = Flowsheet::new();
    let source = flowsheet.add_unit(UnitKind::Feed);
    let unit = flowsheet.add_unit(kind);
    let product = flowsheet.add_unit(UnitKind::Product);
    assert!(flowsheet.set_spec(
        source,
        UnitSpec::Feed {
            temperature_k: Some(number(&inlet["T"])),
            pressure_pa: Some(number(&inlet["p"])),
            molar_flow_mol_s: Some(number(&case["molar_flow_mol_s"])),
        }
    ));
    flowsheet.try_set_spec(unit, spec).unwrap();
    flowsheet
        .connect(PortRef::new(source, 0), PortRef::new(unit, 0))
        .unwrap();
    let outlet = flowsheet
        .connect(PortRef::new(unit, 1), PortRef::new(product, 0))
        .unwrap();
    if kind == UnitKind::Expander {
        let e = flowsheet.add_unit(UnitKind::EnergySink);
        flowsheet
            .connect(PortRef::new(unit, POWER), PortRef::new(e, 0))
            .unwrap();
    } else {
        let e = flowsheet.add_unit(UnitKind::EnergySource);
        assert!(flowsheet.set_spec(
            e,
            UnitSpec::EnergySource {
                duty: EnergySourceDuty::Calculated
            }
        ));
        flowsheet
            .connect(PortRef::new(e, 0), PortRef::new(unit, POWER))
            .unwrap();
    }
    let mut feeds = FeedCompositions::new();
    feeds.insert(source, composition(ids, z));
    (flowsheet, feeds, unit, outlet)
}

/// Solves one reference machine and compares temperatures, power and either
/// absolute calorics (same reference) or differences from the inlet.
fn compare_machine(fixture: &Value, case: &Value, package: &dyn FlowsheetPackage, absolute: bool) {
    let id = case["id"].as_str().unwrap();
    let ids = strings(&case["ids"]);
    let z = numbers(&case["z"]);
    let (fs, feeds, unit, outlet) = machine(case, &ids, &z);
    let calc = fs
        .solve(package, &feeds)
        .unwrap_or_else(|e| panic!("{id}: {e:?}"));
    let report = &calc.units[&unit];
    let balance = report.balance.as_ref().unwrap();
    assert!(balance.max_scaled_component_residual <= 1e-8, "{id}");
    assert!(balance.scaled_energy_residual <= 2e-9, "{id}");
    let stage = &report.rotating.as_ref().unwrap().stages[0];
    let (ideal, reference) = (&case["isentropic_outlet"], &case["outlet"]);
    check(
        fixture,
        &format!("{id} isentropic T"),
        "temperature",
        stage.isentropic_outlet_temperature_k,
        number(&ideal["T"]),
    );
    check(
        fixture,
        &format!("{id} outlet T"),
        "temperature",
        stage.outlet_temperature_k,
        number(&reference["T"]),
    );
    check(
        fixture,
        &format!("{id} isentropic power"),
        "power",
        stage.isentropic_power_w,
        number(&case["isentropic_power_W"]),
    );
    check(
        fixture,
        &format!("{id} power"),
        "power",
        stage.power_w,
        number(&case["power_W"]),
    );
    let out = calc.streams[&outlet].state.as_ref().unwrap();
    assert_eq!(
        out.phases.len() as u64,
        reference["phases"].as_u64().unwrap(),
        "{id}"
    );
    if let Some(beta) = reference["vapor_fraction"].as_f64() {
        check(
            fixture,
            &format!("{id} vapor fraction"),
            "fraction",
            out.vapor_fraction().unwrap(),
            beta,
        );
    }
    if absolute {
        check(
            fixture,
            &format!("{id} outlet h"),
            "h",
            out.enthalpy(),
            number(&reference["h"]),
        );
        check(
            fixture,
            &format!("{id} outlet s"),
            "s",
            out.entropy(),
            number(&reference["s"]),
        );
    } else {
        let inlet = package
            .tp(&z, number(&case["inlet"]["T"]), number(&case["inlet"]["p"]))
            .unwrap();
        check(
            fixture,
            &format!("{id} outlet Δh"),
            "h",
            out.enthalpy() - inlet.enthalpy(),
            number(&reference["h_minus_inlet"]),
        );
        check(
            fixture,
            &format!("{id} outlet Δs"),
            "s",
            out.entropy() - inlet.entropy(),
            number(&reference["s_minus_inlet"]),
        );
    }
}

#[test]
fn gerg_machines_match_clapeyron() {
    let fixture = fixture();
    for case in fixture["gerg_units"].as_array().unwrap() {
        let ids = strings(&case["ids"]);
        let package = gerg(&ids, &numbers(&case["z"]));
        compare_machine(&fixture, case, &package, true);
    }
}

#[test]
fn nrtl_machines_match_the_clapeyron_gamma_phi_reference() {
    let fixture = fixture();
    for case in fixture["nrtl_units"].as_array().unwrap() {
        let ids = strings(&case["ids"]);
        let model = ActivityModel::from_databases(
            &ComponentDatabase::embedded().unwrap(),
            &LiquidDatabase::embedded().unwrap(),
            &composition(&ids, &numbers(&case["z"])),
            VaporModel::IdealGas,
        )
        .unwrap();
        let package = NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap();
        compare_machine(&fixture, case, &package, true);
    }
}

#[test]
fn eos_cg_gas_machines_match_clapeyron() {
    let fixture = fixture();
    for case in fixture["eoscg_units"].as_array().unwrap() {
        let package =
            EosCgPackage::new(strings(&case["ids"]), EosCgPackage::TEMPERATURE_RANGE_K).unwrap();
        compare_machine(&fixture, case, &package, false);
    }
}
