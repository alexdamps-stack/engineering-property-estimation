//! Three-phase distillation column (#188) through the project solve entry
//! point: independent reference, save/reload with results, legacy record.

use std::ops::ControlFlow;

use flowsheet_core::units::MultiphaseSpec;
use flowsheet_core::units::columns::ColumnInternals;
use flowsheet_core::units::multiphase::{DecanterReflux, PhaseReflux, ThreePhaseColumnSpec};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{PackageChoice, PackageSelection, Project, solve_project};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../../validation/multiphase/three-phase-column-v1.json");

fn numbers(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

/// Project of the first reference case: light-liquid total reflux.
fn project(case: &Value, components: &Value) -> (Project, UnitId, [StreamId; 3], StreamId) {
    let ids = components
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
    assert_eq!(case["reflux"]["kind"], "light_total");
    let fs = &mut project.flowsheet;
    let unit = fs.add_unit(UnitKind::ThreePhaseColumn);
    fs.set_spec(
        unit,
        UnitSpec::Multiphase(MultiphaseSpec::ThreePhaseColumn(ThreePhaseColumnSpec {
            internals: ColumnInternals::Trays(case["internal_stages"].as_u64().unwrap() as u16),
            feed_stage: case["feed_stage"].as_u64().unwrap() as u16,
            pressure_pa: case["pressure_Pa"].as_f64(),
            decanter_temperature_k: case["decanter_temperature_K"].as_f64(),
            reflux: DecanterReflux::Light(PhaseReflux::Total),
            bottoms_fraction: case["bottoms_fraction"].as_f64(),
        })),
    );
    let feed = &case["feed"];
    let source = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        source,
        UnitSpec::Feed {
            temperature_k: feed["temperature_K"].as_f64(),
            pressure_pa: feed["pressure_Pa"].as_f64(),
            molar_flow_mol_s: feed["molar_flow_mol_s"].as_f64(),
        },
    );
    project.feeds.insert(source, numbers(&feed["z"]));
    fs.connect(PortRef::new(source, 0), PortRef::new(unit, 0))
        .unwrap();
    let outlets = [1, 2, 3].map(|port| {
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(PortRef::new(unit, port), PortRef::new(product, 0))
            .unwrap()
    });
    let duty = fs.add_unit(UnitKind::EnergySource);
    fs.set_spec(
        duty,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Calculated,
        },
    );
    let reboiler = fs
        .connect(PortRef::new(duty, 0), PortRef::new(unit, 5))
        .unwrap();
    (project, unit, outlets, reboiler)
}

#[test]
fn three_phase_column_solves_saves_and_reloads() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let case = &fixture["cases"][0];
    let (mut project, unit, outlets, reboiler) = project(case, &fixture["components"]);
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
    let profile = saved_unit.three_phase_column.as_ref().unwrap();
    let tol = &fixture["tolerances"];
    let qr = case["reboiler_duty_W"].as_f64().unwrap();
    let [atol, rtol] = [
        tol["duty_W"][0].as_f64().unwrap(),
        tol["duty_W"][1].as_f64().unwrap(),
    ];
    assert!((profile.reboiler_duty_w - qr).abs() <= atol + rtol * qr);
    assert!((result.energy_streams[&reboiler].duty_w - qr).abs() <= atol + rtol * qr);
    for (stream, name) in outlets.iter().zip(["distillate", "decant", "bottoms"]) {
        let product = &case["products"][name];
        let s = &result.streams[stream];
        assert!((s.molar_flow_mol_s - product["molar_flow_mol_s"].as_f64().unwrap()).abs() <= 1e-8);
        for (a, b) in s
            .state
            .as_ref()
            .unwrap()
            .mole_fractions
            .iter()
            .zip(numbers(&product["x"]))
        {
            assert!(
                (a - b).abs() <= tol["composition"][0].as_f64().unwrap(),
                "{name}"
            );
        }
    }
    project.results = Some(result);
    let saved = project.to_json().unwrap();
    let reopened = Project::from_json(&saved).unwrap();
    assert_eq!(reopened.to_json().unwrap(), saved);
    assert_eq!(
        reopened.results.as_ref().unwrap().units[&unit].three_phase_column,
        project.results.as_ref().unwrap().units[&unit].three_phase_column
    );
    assert_eq!(
        reopened.flowsheet.unit(unit).unwrap().spec,
        project.flowsheet.unit(unit).unwrap().spec
    );
}

#[test]
fn drawn_only_three_phase_column_record_opens_with_its_connections() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let (project, unit, _, _) = project(&fixture["cases"][0], &fixture["components"]);
    let mut old: Value = serde_json::from_str(&project.to_json().unwrap()).unwrap();
    let record = old["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "three_phase_column")
        .unwrap();
    record["spec"]["unit"] = json!({"kind": "three_phase_column"});
    let reopened = Project::from_json(&old.to_string()).unwrap();
    let column = reopened.flowsheet.unit(unit).unwrap();
    assert_eq!(
        column.spec,
        UnitSpec::Multiphase(MultiphaseSpec::ThreePhaseColumn(Default::default()))
    );
    for stream in project.flowsheet.streams() {
        assert_eq!(reopened.flowsheet.stream(stream.id).unwrap(), stream);
    }
    assert!(
        column
            .issues()
            .iter()
            .any(|issue| issue.field == "decanter temperature")
    );
}
