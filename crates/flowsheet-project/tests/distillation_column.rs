//! Normal flowsheet execution, independent references, and project persistence.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    ColumnInternals, DistillationColumnSpec, LiquidSideDraw,
    packing::{PackedSectionBasis, PackedSectionSpec, StageRounding},
};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../../validation/columns/distillation-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn count(value: &Value) -> u16 {
    u16::try_from(value.as_u64().unwrap()).unwrap()
}

struct ColumnCase {
    project: Project,
    column: UnitId,
    products: Vec<StreamId>,
    energy_ports: [StreamId; 2],
}

impl ColumnCase {
    fn new(case: &Value) -> Self {
        let variant = match case["variant"].as_str().unwrap() {
            "PR" => CubicVariantChoice::PengRobinson,
            "PR-BM" => CubicVariantChoice::PengRobinsonBostonMathias,
            other => panic!("unexpected fixture variant {other}"),
        };
        let ids: Vec<_> = case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|id| id.as_str().unwrap().to_owned())
            .collect();
        let mut project = Project::new(PackageSelection::zero_interactions(
            variant,
            ids,
            "zero interactions in the independent column references",
        ));
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::DistillationColumn);
        let feeds = case["feeds"].as_array().unwrap();
        let draws = case["draws"].as_array().unwrap();
        let spec = DistillationColumnSpec {
            feeds: u16::try_from(feeds.len()).unwrap(),
            side_draws: u16::try_from(draws.len()).unwrap(),
            feed_stages: feeds.iter().map(|feed| count(&feed["stage"])).collect(),
            liquid_draws: draws
                .iter()
                .map(|draw| LiquidSideDraw {
                    stage: count(&draw["stage"]),
                    molar_flow_mol_s: Some(number(&draw["molar_flow_mol_s"])),
                })
                .collect(),
            internals: match case["packed_height_m"].as_f64() {
                None => ColumnInternals::Trays(count(&case["internal_stages"])),
                Some(height_m) => ColumnInternals::Packed(PackedSectionSpec {
                    hetp_m: Some(number(&case["hetp_m"])),
                    basis: PackedSectionBasis::PackedHeight {
                        height_m: Some(height_m),
                        rounding: StageRounding::Nearest,
                    },
                }),
            },
            pressure_pa: Some(number(&case["pressure_Pa"])),
            reflux_ratio: Some(number(&case["reflux_ratio"])),
            distillate_fraction: Some(number(&case["distillate_fraction"])),
        };
        fs.set_spec(column, UnitSpec::Columns(ColumnSpec::Distillation(spec)));
        for (i, f) in feeds.iter().enumerate() {
            let feed = fs.add_unit(UnitKind::Feed);
            fs.set_spec(
                feed,
                UnitSpec::Feed {
                    temperature_k: Some(number(&f["temperature_K"])),
                    pressure_pa: Some(number(&f["pressure_Pa"])),
                    molar_flow_mol_s: Some(number(&f["molar_flow_mol_s"])),
                },
            );
            project.feeds.insert(
                feed,
                f["z"].as_array().unwrap().iter().map(number).collect(),
            );
            fs.connect(PortRef::new(feed, 0), PortRef::new(column, 4 + i))
                .unwrap();
        }
        let product_ports = [0, 1]
            .into_iter()
            .chain((0..draws.len()).map(|i| 4 + feeds.len() + i));
        let products = product_ports
            .map(|i| {
                let product = fs.add_unit(UnitKind::Product);
                fs.connect(PortRef::new(column, i), PortRef::new(product, 0))
                    .unwrap()
            })
            .collect();
        let sink = fs.add_unit(UnitKind::EnergySink);
        let source = fs.add_unit(UnitKind::EnergySource);
        fs.set_spec(
            source,
            UnitSpec::EnergySource {
                duty: EnergySourceDuty::Calculated,
            },
        );
        let condenser = fs
            .connect(PortRef::new(column, 2), PortRef::new(sink, 0))
            .unwrap();
        let reboiler = fs
            .connect(PortRef::new(source, 0), PortRef::new(column, 3))
            .unwrap();
        Self {
            project,
            column,
            products,
            energy_ports: [condenser, reboiler],
        }
    }

    fn solve(&self) -> Result<ResultSet, SolveFailure> {
        solve_project(
            &self.project,
            &self.project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        )
    }
}

#[test]
fn single_feed_trays_and_two_feed_packing_match_independent_mesh_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let mut rig = ColumnCase::new(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{}: {}", case["id"], e.message));
        let unit = &result.units[&rig.column];
        let column = unit.distillation_column.as_ref().unwrap();
        let compare = |field: &str, actual: f64, expected: f64| {
            let range = &fixture["tolerances"][field];
            let tolerance = number(&range[0]) + number(&range[1]) * expected.abs();
            assert!(
                (actual - expected).abs() <= tolerance,
                "{} {field}: {actual} differs from {expected} by more than {tolerance}",
                case["id"]
            );
        };
        let temperatures = case["stage_temperatures_K"].as_array().unwrap();
        assert_eq!(column.stage_temperatures_k.len(), temperatures.len());
        for (actual, expected) in column.stage_temperatures_k.iter().zip(temperatures) {
            compare("temperature_K", *actual, number(expected));
        }
        for (port, duty, field) in [
            (
                rig.energy_ports[0],
                column.condenser_duty_w,
                "condenser_duty_W",
            ),
            (
                rig.energy_ports[1],
                column.reboiler_duty_w,
                "reboiler_duty_W",
            ),
        ] {
            compare("duty_W", duty, number(&case[field]));
            assert_eq!(result.energy_streams[&port].duty_w, duty);
        }
        assert_eq!(column.packed_height_m, case["packed_height_m"].as_f64());
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
        for (stream, expected) in rig.products.iter().zip(case["outlets"].as_array().unwrap()) {
            let actual = &result.streams[stream];
            compare(
                "molar_flow_mol_s",
                actual.molar_flow_mol_s,
                number(&expected["molar_flow_mol_s"]),
            );
            let phase = actual.state.as_ref().unwrap();
            assert_eq!(phase.vapor_fraction, Some(0.0));
            for (actual, expected) in phase
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                compare("composition", *actual, number(expected));
            }
        }
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved);
        assert_eq!(
            reopened.flowsheet.unit(rig.column).unwrap().spec,
            rig.project.flowsheet.unit(rig.column).unwrap().spec
        );
    }
}

#[test]
fn ternary_column_closes_with_an_inactive_package_component() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let mut case = fixture["cases"][0].clone();
    case["ids"] = json!(["methane", "nitrogen", "ethane", "propane"]);
    case["feeds"][0]["z"] = json!([0.06, 0.0, 0.39, 0.55]);
    let rig = ColumnCase::new(&case);
    let result = rig.solve().unwrap();
    let unit = &result.units[&rig.column];
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
    assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
    let product = |stream| {
        &result.streams[&stream]
            .state
            .as_ref()
            .unwrap()
            .mole_fractions
    };
    let top = product(rig.products[0]);
    let bottom = product(rig.products[1]);
    assert_eq!((top[1], bottom[1]), (0.0, 0.0));
    assert!(top[0] > 0.06 && bottom[0] < 0.06);
}

#[test]
fn legacy_column_keeps_its_feed_and_draw_connections() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let rig = ColumnCase::new(&fixture["cases"][1]);
    let mut old: Value = serde_json::from_str(&rig.project.to_json().unwrap()).unwrap();
    let record = old["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "distillation_column")
        .unwrap();
    record["spec"]["unit"] = json!({"kind": "distillation_column", "feeds": 2, "side_draws": 1});
    let reopened = Project::from_json(&old.to_string()).unwrap();
    let fs = &reopened.flowsheet;
    let column = fs.unit(rig.column).unwrap();
    assert_eq!(
        column.ports(),
        rig.project.flowsheet.unit(rig.column).unwrap().ports()
    );
    for stream in rig.project.flowsheet.streams() {
        assert_eq!(fs.stream(stream.id).unwrap(), stream);
    }
    assert!(
        column
            .issues()
            .iter()
            .any(|issue| issue.field == "column pressure")
    );
    assert!(
        column
            .issues()
            .iter()
            .any(|issue| issue.field == "side-draw molar flow")
    );
}

#[test]
fn impossible_draw_or_feed_pressure_fails_at_the_column() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let mut rig = ColumnCase::new(&fixture["cases"][1]);
    let original = rig.project.flowsheet.unit(rig.column).unwrap().spec.clone();
    let UnitSpec::Columns(ColumnSpec::Distillation(mut spec)) = original else {
        unreachable!();
    };
    spec.liquid_draws[0].molar_flow_mol_s = Some(80.0);
    rig.project.flowsheet.set_spec(
        rig.column,
        UnitSpec::Columns(ColumnSpec::Distillation(spec.clone())),
    );
    let error = rig.solve().unwrap_err();
    assert_eq!(error.unit, Some(rig.column));
    assert!(
        error.message.contains("positive bottoms flow"),
        "{}",
        error.message
    );

    spec.liquid_draws[0].molar_flow_mol_s = Some(6.0);
    spec.pressure_pa = Some(2.1e6);
    rig.project.flowsheet.set_spec(
        rig.column,
        UnitSpec::Columns(ColumnSpec::Distillation(spec)),
    );
    let error = rig.solve().unwrap_err();
    assert_eq!(error.unit, Some(rig.column));
    assert!(
        error.message.contains("exceeds feed 2 pressure"),
        "{}",
        error.message
    );
}
