//! Countercurrent column calculations through the project solve entry point.

use std::ops::ControlFlow;

use flowsheet_core::units::{
    ColumnSpec,
    columns::{
        AbsorberSpec, ColumnInternals,
        packing::{PackedSectionBasis, PackedSectionSpec, StageRounding},
    },
};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../../validation/columns/absorber-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

struct AbsorberCase {
    project: Project,
    column: UnitId,
    gas: UnitId,
    outlets: [StreamId; 2],
}

impl AbsorberCase {
    fn new(case: &Value) -> Self {
        let variant = match case["variant"].as_str().unwrap() {
            "PR" => CubicVariantChoice::PengRobinson,
            "PR-BM" => CubicVariantChoice::PengRobinsonBostonMathias,
            other => panic!("unexpected fixture variant {other}"),
        };
        let ids = case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|id| id.as_str().unwrap().to_owned())
            .collect();
        let mut project = Project::new(PackageSelection::zero_interactions(
            variant,
            ids,
            "zero interactions in the independent absorber references",
        ));
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::Absorber);
        fs.set_spec(
            column,
            UnitSpec::Columns(ColumnSpec::Absorber(AbsorberSpec {
                pressure_pa: Some(number(&case["pressure_Pa"])),
                internals: match case["packed_height_m"].as_f64() {
                    None => ColumnInternals::Trays(
                        u16::try_from(case["stages"].as_u64().unwrap()).unwrap(),
                    ),
                    Some(height_m) => ColumnInternals::Packed(PackedSectionSpec {
                        hetp_m: Some(number(&case["hetp_m"])),
                        basis: PackedSectionBasis::PackedHeight {
                            height_m: Some(height_m),
                            rounding: StageRounding::Nearest,
                        },
                    }),
                },
                ..Default::default()
            })),
        );
        let [gas, _liquid] = ["gas", "liquid"].map(|port| {
            let feed = fs.add_unit(UnitKind::Feed);
            let inlet = &case[port];
            fs.set_spec(
                feed,
                UnitSpec::Feed {
                    temperature_k: Some(number(&inlet["temperature_K"])),
                    pressure_pa: Some(number(&inlet["pressure_Pa"])),
                    molar_flow_mol_s: Some(number(&inlet["molar_flow_mol_s"])),
                },
            );
            project.feeds.insert(
                feed,
                inlet["z"].as_array().unwrap().iter().map(number).collect(),
            );
            let index = usize::from(port == "liquid");
            fs.connect(PortRef::new(feed, 0), PortRef::new(column, index))
                .unwrap();
            feed
        });
        let outlets = [2, 3].map(|port| {
            let product = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(column, port), PortRef::new(product, 0))
                .unwrap()
        });
        Self {
            project,
            column,
            gas,
            outlets,
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
fn adiabatic_trays_and_packing_match_independent_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let mut rig = AbsorberCase::new(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{}: {}", case["id"], e.message));
        let unit = &result.units[&rig.column];
        let column = unit.absorber.as_ref().unwrap();
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
        assert_eq!(column.packed_height_m, case["packed_height_m"].as_f64());
        assert_eq!(unit.duty_w, 0.0);
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
        for (i, expected) in case["outlets"].as_array().unwrap().iter().enumerate() {
            let actual = &result.streams[&rig.outlets[i]];
            compare(
                "molar_flow_mol_s",
                actual.molar_flow_mol_s,
                number(&expected["molar_flow_mol_s"]),
            );
            let state = actual.state.as_ref().unwrap();
            assert_eq!(state.vapor_fraction, Some(if i == 0 { 1.0 } else { 0.0 }));
            assert_eq!(state.pressure_pa, number(&case["pressure_Pa"]));
            compare(
                "temperature_K",
                state.temperature_k,
                number(if i == 0 {
                    &temperatures[0]
                } else {
                    temperatures.last().unwrap()
                }),
            );
            for (actual, expected) in state
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                compare("composition", *actual, number(expected));
            }
        }
        let top = &result.streams[&rig.outlets[0]];
        let z = &top.state.as_ref().unwrap().mole_fractions;
        let gas_flow = number(&case["gas"]["molar_flow_mol_s"]);
        if case["id"] == "propane-absorption-PR" {
            assert!(top.molar_flow_mol_s * z[1] < gas_flow * number(&case["gas"]["z"][1]));
        }
        if case["id"] == "ethane-stripping-packed-PR-BM" {
            assert!(top.molar_flow_mol_s * z[0] > gas_flow * number(&case["gas"]["z"][0]));
        }
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved);
    }
}

#[test]
fn legacy_absorber_keeps_the_four_connected_ports() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let rig = AbsorberCase::new(&fixture["cases"][1]);
    let mut old: Value = serde_json::from_str(&rig.project.to_json().unwrap()).unwrap();
    let record = old["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "absorber")
        .unwrap();
    record["spec"]["unit"] = json!({"kind": "absorber"});
    let reopened = Project::from_json(&old.to_string()).unwrap();
    let fs = &reopened.flowsheet;
    let column = fs.unit(rig.column).unwrap();
    assert_eq!(
        column.spec,
        UnitSpec::Columns(ColumnSpec::Absorber(AbsorberSpec::default()))
    );
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
}

#[test]
fn gas_feed_rejects_liquid_and_insufficient_pressure() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let mut rig = AbsorberCase::new(&fixture["cases"][0]);
    for (temperature, pressure, message) in [
        (220.0, 2.2e6, "gas inlet requires a single Vapor phase"),
        (330.0, 1.9e6, "absorber pressure exceeds gas inlet pressure"),
    ] {
        rig.project.flowsheet.set_spec(
            rig.gas,
            UnitSpec::Feed {
                temperature_k: Some(temperature),
                pressure_pa: Some(pressure),
                molar_flow_mol_s: Some(70.0),
            },
        );
        let error = rig.solve().unwrap_err();
        assert_eq!(error.unit, Some(rig.column));
        assert!(error.message.contains(message), "{}", error.message);
    }
}
