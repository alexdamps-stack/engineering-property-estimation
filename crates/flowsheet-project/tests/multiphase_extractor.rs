//! Liquid–liquid extractor (#189) through the project solve entry point:
//! independent reference, save/reload with results, legacy record.

use std::ops::ControlFlow;

use flowsheet_core::units::MultiphaseSpec;
use flowsheet_core::units::multiphase::{
    ExtractorInternals, ExtractorThermal, FeedPhase, LiquidExtractorSpec,
};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{PackageChoice, PackageSelection, Project, solve_project};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../../validation/multiphase/extractor-v1.json");

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn project(case: &Value) -> (Project, UnitId, [StreamId; 2]) {
    let ids = case["components"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_owned())
        .collect();
    let mut project = Project::new(PackageSelection::default_for(
        PackageChoice::NrtlIdealGas,
        ids,
        "",
    ));
    let (t, p) = (
        case["temperature_K"].as_f64().unwrap(),
        case["pressure_Pa"].as_f64().unwrap(),
    );
    let fs = &mut project.flowsheet;
    let unit = fs.add_unit(UnitKind::LiquidExtractor);
    fs.set_spec(
        unit,
        UnitSpec::Multiphase(MultiphaseSpec::LiquidExtractor(LiquidExtractorSpec {
            internals: ExtractorInternals::Stages(case["stages"].as_u64().unwrap() as u16),
            pressure_pa: Some(p),
            thermal: ExtractorThermal::Isothermal {
                temperature_k: Some(t),
            },
            feed_phase: FeedPhase::Heavy,
        })),
    );
    for (port, inlet) in ["feed", "solvent"].into_iter().enumerate() {
        let source = fs.add_unit(UnitKind::Feed);
        fs.set_spec(
            source,
            UnitSpec::Feed {
                temperature_k: Some(t),
                pressure_pa: Some(p),
                molar_flow_mol_s: case[inlet]["molar_flow_mol_s"].as_f64(),
            },
        );
        project.feeds.insert(source, numbers(&case[inlet]["z"]));
        fs.connect(PortRef::new(source, 0), PortRef::new(unit, port))
            .unwrap();
    }
    let outlets = [2, 3].map(|port| {
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(PortRef::new(unit, port), PortRef::new(product, 0))
            .unwrap()
    });
    (project, unit, outlets)
}

#[test]
fn extractor_solves_saves_and_reloads() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let mut case = fixture["cases"][1].clone();
    case["components"] = fixture["components"].clone();
    let (mut project, unit, outlets) = project(&case);
    let result = solve_project(
        &project,
        &project.package.prepare().unwrap(),
        |_| ControlFlow::Continue(()),
        0,
    )
    .unwrap_or_else(|e| panic!("{}", e.message));
    let saved_unit = &result.units[&unit];
    assert!(saved_unit.max_scaled_component_residual.unwrap() < 1e-8);
    assert!(saved_unit.scaled_energy_residual.unwrap().abs() < 2e-9);
    let profile = saved_unit.liquid_extractor.as_ref().unwrap();
    let reference = case["profile"].as_array().unwrap();
    let last = reference.last().unwrap();
    let raffinate = &result.streams[&outlets[1]];
    let tol = &fixture["tolerances"];
    for (a, b) in raffinate
        .state
        .as_ref()
        .unwrap()
        .mole_fractions
        .iter()
        .zip(numbers(&last["x_raffinate"]))
    {
        assert!((a - b).abs() <= tol["composition"][0].as_f64().unwrap());
    }
    assert_eq!(profile.stage_temperatures_k.len(), reference.len());
    assert!(profile.feed_is_heavy);
    assert_eq!(profile.duty_w, saved_unit.duty_w);
    project.results = Some(result);
    let saved = project.to_json().unwrap();
    let reopened = Project::from_json(&saved).unwrap();
    assert_eq!(reopened.to_json().unwrap(), saved);
    assert_eq!(
        reopened.results.as_ref().unwrap().units[&unit].liquid_extractor,
        project.results.as_ref().unwrap().units[&unit].liquid_extractor
    );
    assert_eq!(
        reopened.flowsheet.unit(unit).unwrap().spec,
        project.flowsheet.unit(unit).unwrap().spec
    );
}

#[test]
fn drawn_only_extractor_record_opens_with_its_connections() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let mut case = fixture["cases"][0].clone();
    case["components"] = fixture["components"].clone();
    let (project, unit, _) = project(&case);
    let mut old: Value = serde_json::from_str(&project.to_json().unwrap()).unwrap();
    let record = old["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "liquid_extractor")
        .unwrap();
    record["spec"]["unit"] = json!({"kind": "liquid_extractor"});
    let reopened = Project::from_json(&old.to_string()).unwrap();
    let extractor = reopened.flowsheet.unit(unit).unwrap();
    assert_eq!(
        extractor.spec,
        UnitSpec::Multiphase(MultiphaseSpec::LiquidExtractor(Default::default()))
    );
    for stream in project.flowsheet.streams() {
        assert_eq!(reopened.flowsheet.stream(stream.id).unwrap(), stream);
    }
    assert!(
        extractor
            .issues()
            .iter()
            .any(|issue| issue.field == "column pressure")
    );
}
