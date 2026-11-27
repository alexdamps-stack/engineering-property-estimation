//! Equilibrium-stage columns on the package domain (#186, #187): independent
//! references outside the original dry-gas domain, gamma–phi (NRTL-IG)
//! columns, the stage liquid–liquid check and unsupported packages.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{AbsorberSpec, ColumnInternals, DistillationColumnSpec};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    PackageChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::Value;

const FIXTURE: &str = include_str!("../../../validation/columns/extended-v1.json");
const LIQUID_SPLIT: &str = "liquid splits into two liquids; use the three-phase column";
const VLLE_UNCHECKED: &str = "stage liquids were not tested for a second liquid phase";

#[derive(Clone, Copy)]
struct Inlet<'a> {
    temperature_k: f64,
    pressure_pa: f64,
    molar_flow_mol_s: f64,
    z: &'a [f64],
}

impl Inlet<'_> {
    fn from_json(value: &Value) -> (Vec<f64>, [f64; 3]) {
        (
            value["z"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect(),
            [
                value["temperature_K"].as_f64().unwrap(),
                value["pressure_Pa"].as_f64().unwrap(),
                value["molar_flow_mol_s"].as_f64().unwrap(),
            ],
        )
    }
}

struct Rig {
    project: Project,
    column: UnitId,
    products: Vec<StreamId>,
}

impl Rig {
    fn add_feed(project: &mut Project, inlet: Inlet<'_>, to: PortRef) {
        let fs = &mut project.flowsheet;
        let feed = fs.add_unit(UnitKind::Feed);
        fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(inlet.temperature_k),
                pressure_pa: Some(inlet.pressure_pa),
                molar_flow_mol_s: Some(inlet.molar_flow_mol_s),
            },
        );
        project.feeds.insert(feed, inlet.z.to_vec());
        fs.connect(PortRef::new(feed, 0), to).unwrap();
    }

    fn products(project: &mut Project, column: UnitId, ports: [usize; 2]) -> Vec<StreamId> {
        let fs = &mut project.flowsheet;
        ports
            .into_iter()
            .map(|port| {
                let product = fs.add_unit(UnitKind::Product);
                fs.connect(PortRef::new(column, port), PortRef::new(product, 0))
                    .unwrap()
            })
            .collect()
    }

    fn distillation(
        selection: PackageSelection,
        spec: DistillationColumnSpec,
        feed: Inlet<'_>,
    ) -> Self {
        let mut project = Project::new(selection);
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::DistillationColumn);
        fs.set_spec(column, UnitSpec::Columns(ColumnSpec::Distillation(spec)));
        let sink = fs.add_unit(UnitKind::EnergySink);
        let source = fs.add_unit(UnitKind::EnergySource);
        fs.set_spec(
            source,
            UnitSpec::EnergySource {
                duty: EnergySourceDuty::Calculated,
            },
        );
        fs.connect(PortRef::new(column, 2), PortRef::new(sink, 0))
            .unwrap();
        fs.connect(PortRef::new(source, 0), PortRef::new(column, 3))
            .unwrap();
        Self::add_feed(&mut project, feed, PortRef::new(column, 4));
        let products = Self::products(&mut project, column, [0, 1]);
        Self {
            project,
            column,
            products,
        }
    }

    fn absorber(
        selection: PackageSelection,
        spec: AbsorberSpec,
        gas: Inlet<'_>,
        liquid: Inlet<'_>,
    ) -> Self {
        let mut project = Project::new(selection);
        let column = project.flowsheet.add_unit(UnitKind::Absorber);
        project
            .flowsheet
            .set_spec(column, UnitSpec::Columns(ColumnSpec::Absorber(spec)));
        Self::add_feed(&mut project, gas, PortRef::new(column, 0));
        Self::add_feed(&mut project, liquid, PortRef::new(column, 1));
        let products = Self::products(&mut project, column, [2, 3]);
        Self {
            project,
            column,
            products,
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

fn ids(components: &[&str]) -> Vec<String> {
    components.iter().map(|id| id.to_string()).collect()
}

fn cubic(case: &Value) -> PackageSelection {
    assert_eq!(case["variant"], "PR");
    let mut selection = PackageSelection::zero_interactions(
        PackageChoice::PengRobinson,
        case["ids"]
            .as_array()
            .unwrap()
            .iter()
            .map(|id| id.as_str().unwrap().to_owned())
            .collect(),
        case["kij_note"].as_str().unwrap(),
    );
    selection.kij = case["kij"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect()
        })
        .collect();
    selection
}

fn embedded(variant: PackageChoice, components: &[&str]) -> PackageSelection {
    PackageSelection::default_for(variant, ids(components), "no interaction matrix")
}

fn distillation_spec(
    stages: u16,
    feed_stage: u16,
    pressure: f64,
    reflux: f64,
    d: f64,
) -> DistillationColumnSpec {
    DistillationColumnSpec {
        feed_stages: vec![feed_stage],
        internals: ColumnInternals::Trays(stages),
        pressure_pa: Some(pressure),
        reflux_ratio: Some(reflux),
        distillate_fraction: Some(d),
        ..Default::default()
    }
}

fn rig_for(case: &Value) -> Rig {
    let number = |v: &Value| v.as_f64().unwrap();
    let count = |v: &Value| u16::try_from(v.as_u64().unwrap()).unwrap();
    match case["kind"].as_str().unwrap() {
        "distillation" => {
            let feed = &case["feeds"][0];
            assert_eq!(case["feeds"].as_array().unwrap().len(), 1);
            let (z, [t, p, f]) = Inlet::from_json(feed);
            Rig::distillation(
                cubic(case),
                distillation_spec(
                    count(&case["internal_stages"]),
                    count(&feed["stage"]),
                    number(&case["pressure_Pa"]),
                    number(&case["reflux_ratio"]),
                    number(&case["distillate_fraction"]),
                ),
                Inlet {
                    temperature_k: t,
                    pressure_pa: p,
                    molar_flow_mol_s: f,
                    z: &z,
                },
            )
        }
        "absorber" => {
            let (gz, [gt, gp, gf]) = Inlet::from_json(&case["gas"]);
            let (lz, [lt, lp, lf]) = Inlet::from_json(&case["liquid"]);
            Rig::absorber(
                cubic(case),
                AbsorberSpec {
                    pressure_pa: Some(number(&case["pressure_Pa"])),
                    internals: ColumnInternals::Trays(count(&case["stages"])),
                },
                Inlet {
                    temperature_k: gt,
                    pressure_pa: gp,
                    molar_flow_mol_s: gf,
                    z: &gz,
                },
                Inlet {
                    temperature_k: lt,
                    pressure_pa: lp,
                    molar_flow_mol_s: lf,
                    z: &lz,
                },
            )
        }
        other => panic!("unexpected case kind {other}"),
    }
}

#[test]
fn columns_with_water_methanol_and_hot_reboilers_match_independent_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let id = case["id"].as_str().unwrap();
        let mut rig = rig_for(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{id}: {}", e.message));
        let compare = |field: &str, actual: f64, expected: &Value| {
            let expected = expected.as_f64().unwrap();
            let range = &fixture["tolerances"][field];
            let tolerance =
                range[0].as_f64().unwrap() + range[1].as_f64().unwrap() * expected.abs();
            assert!(
                (actual - expected).abs() <= tolerance,
                "{id} {field}: {actual} differs from {expected} by more than {tolerance}"
            );
        };
        let unit = &result.units[&rig.column];
        assert!(unit.warnings.is_empty(), "{id}: {:?}", unit.warnings);
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
        let temperatures = match &unit.distillation_column {
            Some(column) => {
                compare("duty_W", column.condenser_duty_w, &case["condenser_duty_W"]);
                compare("duty_W", column.reboiler_duty_w, &case["reboiler_duty_W"]);
                &column.stage_temperatures_k
            }
            None => &unit.absorber.as_ref().unwrap().stage_temperatures_k,
        };
        let expected = case["stage_temperatures_K"].as_array().unwrap();
        assert_eq!(temperatures.len(), expected.len(), "{id}");
        for (actual, expected) in temperatures.iter().zip(expected) {
            compare("temperature_K", *actual, expected);
        }
        let outlets = case["outlets"].as_array().unwrap();
        for (stream, expected) in rig.products.iter().zip(outlets) {
            let actual = &result.streams[stream];
            compare(
                "molar_flow_mol_s",
                actual.molar_flow_mol_s,
                &expected["molar_flow_mol_s"],
            );
            let state = actual.state.as_ref().unwrap();
            for (actual, expected) in state
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                compare("composition", *actual, expected);
            }
        }
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        assert_eq!(
            Project::from_json(&saved).unwrap().to_json().unwrap(),
            saved
        );
    }
}

fn methanol_water_feed(z: &[f64]) -> Inlet<'_> {
    Inlet {
        temperature_k: 340.0,
        pressure_pa: 1.2e5,
        molar_flow_mol_s: 100.0,
        z,
    }
}

#[test]
fn nrtl_methanol_water_column_solves_through_the_package_interface() {
    let rig = Rig::distillation(
        embedded(PackageChoice::NrtlIdealGas, &["methanol", "water"]),
        distillation_spec(6, 5, 101_325.0, 1.5, 0.5),
        methanol_water_feed(&[0.5, 0.5]),
    );
    let result = rig.solve().unwrap_or_else(|e| panic!("{}", e.message));
    let unit = &result.units[&rig.column];
    assert!(unit.warnings.is_empty(), "{:?}", unit.warnings);
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
    assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
    let column = unit.distillation_column.as_ref().unwrap();
    // Methanol boils at 337.7 K and water at 373.1 K at 1 atm.
    let temperatures = &column.stage_temperatures_k;
    assert!(temperatures.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(temperatures[0] > 337.0 && temperatures[7] < 373.2);
    let methanol = |stream: StreamId| {
        result.streams[&stream]
            .state
            .as_ref()
            .unwrap()
            .mole_fractions[0]
    };
    assert!(methanol(rig.products[0]) > 0.9 && methanol(rig.products[1]) < 0.1);
    assert!(column.condenser_duty_w > 0.0 && column.reboiler_duty_w > 0.0);
}

#[test]
fn nrtl_absorber_solves_through_the_package_interface() {
    // Methanol vapor stripped into, and partly absorbed by, cold water.
    let rig = Rig::absorber(
        embedded(PackageChoice::NrtlIdealGas, &["methanol", "water"]),
        AbsorberSpec {
            pressure_pa: Some(101_325.0),
            internals: ColumnInternals::Trays(4),
        },
        Inlet {
            temperature_k: 360.0,
            pressure_pa: 1.1e5,
            molar_flow_mol_s: 20.0,
            z: &[0.8, 0.2],
        },
        Inlet {
            temperature_k: 300.0,
            pressure_pa: 1.1e5,
            molar_flow_mol_s: 60.0,
            z: &[0.0, 1.0],
        },
    );
    let result = rig.solve().unwrap_or_else(|e| panic!("{}", e.message));
    let unit = &result.units[&rig.column];
    assert!(unit.warnings.is_empty(), "{:?}", unit.warnings);
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
    assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
}

#[test]
fn stage_liquid_that_splits_into_two_liquids_fails_the_column() {
    // Water/1-butanol (ChemSep/DECHEMA NRTL): the overhead approaches the
    // heterogeneous azeotrope, whose condensate separates into two liquids.
    let rig = Rig::distillation(
        embedded(PackageChoice::NrtlIdealGas, &["water", "1-butanol"]),
        distillation_spec(6, 5, 101_325.0, 1.0, 0.05),
        Inlet {
            temperature_k: 350.0,
            pressure_pa: 1.2e5,
            molar_flow_mol_s: 100.0,
            z: &[0.996, 0.004],
        },
    );
    let error = rig.solve().unwrap_err();
    assert_eq!(error.unit, Some(rig.column));
    assert!(error.message.contains(LIQUID_SPLIT), "{}", error.message);

    // PR with kij = 0 splits methanol/water into two liquids near the reboiler.
    let mut selection = embedded(PackageChoice::PengRobinson, &["methanol", "water"]);
    selection.interaction_source = "zero interactions; an approximation".into();
    let rig = Rig::distillation(
        selection,
        distillation_spec(6, 5, 101_325.0, 1.5, 0.5),
        methanol_water_feed(&[0.5, 0.5]),
    );
    let error = rig.solve().unwrap_err();
    assert!(error.message.contains(LIQUID_SPLIT), "{}", error.message);
}

#[test]
fn package_without_three_phase_flash_warns_that_liquids_are_unchecked() {
    let rig = Rig::distillation(
        embedded(PackageChoice::Gerg2008, &["ethane", "propane"]),
        distillation_spec(4, 4, 2e6, 2.5, 0.45),
        Inlet {
            temperature_k: 295.0,
            pressure_pa: 2.3e6,
            molar_flow_mol_s: 100.0,
            z: &[0.45, 0.55],
        },
    );
    let result = rig.solve().unwrap_or_else(|e| panic!("{}", e.message));
    let unit = &result.units[&rig.column];
    assert_eq!(unit.warnings.len(), 1);
    assert!(
        unit.warnings[0].contains(VLLE_UNCHECKED),
        "{:?}",
        unit.warnings
    );
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
}

#[test]
fn package_without_column_operations_fails_explicitly() {
    let selection = embedded(PackageChoice::EosCg2021, &["carbon-dioxide", "hydrogen"]);
    let gas = Inlet {
        temperature_k: 340.0,
        pressure_pa: 1.2e6,
        molar_flow_mol_s: 100.0,
        z: &[0.5, 0.5],
    };
    let distillation = Rig::distillation(
        selection.clone(),
        distillation_spec(4, 3, 1e6, 2.5, 0.45),
        gas,
    );
    let absorber = Rig::absorber(
        selection,
        AbsorberSpec {
            pressure_pa: Some(1e6),
            internals: ColumnInternals::Trays(3),
        },
        gas,
        Inlet {
            z: &[0.9, 0.1],
            ..gas
        },
    );
    for rig in [distillation, absorber] {
        let error = rig.solve().unwrap_err();
        assert_eq!(error.unit, Some(rig.column));
        assert!(error.message.contains("unsupported"), "{}", error.message);
    }
}
