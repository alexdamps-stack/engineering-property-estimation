//! Rigorous distillation (#186) specification set: partial condensers,
//! pressure profiles, vapor side draws and alternative column
//! specifications against independent references, explicit failures of
//! infeasible specifications, and project persistence.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{
    BottomSpecification, ColumnInternals, ColumnPressureDrop, ComponentTarget,
    DistillationColumnSpec, DistillationCondenser, LiquidSideDraw, TopSpecification, VaporSideDraw,
};
use flowsheet_core::{EnergySourceDuty, PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    PackageChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

const FIXTURE: &str = include_str!("../../../validation/columns/specs-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn count(value: &Value) -> u16 {
    u16::try_from(value.as_u64().unwrap()).unwrap()
}

fn fixture() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}

fn case(id: &str) -> Value {
    fixture()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == id)
        .unwrap_or_else(|| panic!("fixture case {id}"))
        .clone()
}

fn selection(case: &Value) -> PackageSelection {
    let variant = match case["variant"].as_str().unwrap() {
        "PR" => PackageChoice::PengRobinson,
        "PR-BM" => PackageChoice::PengRobinsonBostonMathias,
        other => panic!("unexpected fixture variant {other}"),
    };
    let mut selection =
        PackageSelection::zero_interactions(variant, ids(case), case["kij_note"].as_str().unwrap());
    selection.kij = case["kij"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| row.as_array().unwrap().iter().map(number).collect())
        .collect();
    selection
}

fn ids(case: &Value) -> Vec<String> {
    case["ids"]
        .as_array()
        .unwrap()
        .iter()
        .map(|id| id.as_str().unwrap().to_owned())
        .collect()
}

fn target(case: &Value, value: &Value) -> ComponentTarget {
    ComponentTarget {
        component: Some(ids(case)[value["component"].as_u64().unwrap() as usize].clone()),
        value: Some(number(&value["value"])),
    }
}

/// The column specification of a fixture case.
fn spec_for(case: &Value) -> DistillationColumnSpec {
    let draws = |key: &str| case[key].as_array().unwrap().clone();
    let (top, bottom) = (&case["top"], &case["bottom"]);
    DistillationColumnSpec {
        feeds: u16::try_from(case["feeds"].as_array().unwrap().len()).unwrap(),
        feed_stages: draws("feeds").iter().map(|f| count(&f["stage"])).collect(),
        side_draws: u16::try_from(draws("draws").len()).unwrap(),
        liquid_draws: draws("draws")
            .iter()
            .map(|d| LiquidSideDraw {
                stage: count(&d["stage"]),
                molar_flow_mol_s: Some(number(&d["molar_flow_mol_s"])),
            })
            .collect(),
        vapor_side_draws: u16::try_from(draws("vapor_draws").len()).unwrap(),
        vapor_draws: draws("vapor_draws")
            .iter()
            .map(|d| VaporSideDraw {
                stage: count(&d["stage"]),
                molar_flow_mol_s: Some(number(&d["molar_flow_mol_s"])),
            })
            .collect(),
        internals: ColumnInternals::Trays(count(&case["internal_stages"])),
        pressure_pa: Some(number(&case["pressure_Pa"])),
        condenser_pressure_drop_pa: number(&case["condenser_pressure_drop_Pa"]),
        pressure_drop: ColumnPressureDrop::PerStage(number(&case["stage_pressure_drop_Pa"])),
        condenser: match case["condenser"].as_str().unwrap() {
            "total" => DistillationCondenser::Total,
            "partial_vapor" => DistillationCondenser::PartialVapor,
            "partial_vapor_liquid" => DistillationCondenser::PartialVaporLiquid {
                liquid_distillate_mol_s: Some(number(&case["liquid_distillate_mol_s"])),
            },
            other => panic!("unexpected condenser {other}"),
        },
        top: if let Some(r) = top.get("reflux_ratio") {
            TopSpecification::RefluxRatio(Some(number(r)))
        } else if let Some(q) = top.get("condenser_duty_W") {
            TopSpecification::CondenserDuty(Some(number(q)))
        } else if let Some(t) = top.get("distillate_mole_fraction") {
            TopSpecification::DistillatePurity(target(case, t))
        } else {
            panic!("unexpected top specification {top}")
        },
        bottom: if let Some(f) = bottom.get("distillate_fraction") {
            BottomSpecification::DistillateFraction(Some(number(f)))
        } else if let Some(q) = bottom.get("reboiler_duty_W") {
            BottomSpecification::ReboilerDuty(Some(number(q)))
        } else if let Some(r) = bottom.get("boilup_ratio") {
            BottomSpecification::BoilupRatio(Some(number(r)))
        } else if let Some(t) = bottom.get("bottoms_recovery") {
            BottomSpecification::BottomsRecovery(target(case, t))
        } else {
            panic!("unexpected bottom specification {bottom}")
        },
    }
}

struct Rig {
    project: Project,
    column: UnitId,
    /// Material products in output order: distillate, bottoms, liquid draws,
    /// vapor draws, liquid distillate.
    products: Vec<StreamId>,
    energy: [StreamId; 2],
}

impl Rig {
    fn new(selection: PackageSelection, spec: DistillationColumnSpec, feeds: &[Value]) -> Self {
        let mut project = Project::new(selection);
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::DistillationColumn);
        fs.set_spec(
            column,
            UnitSpec::Columns(ColumnSpec::Distillation(spec.clone())),
        );
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
        let ports = fs.unit(column).unwrap().ports();
        let outlets: Vec<usize> = [0, 1]
            .into_iter()
            .chain(
                (4..ports.len())
                    .filter(|&i| !ports[i].label().starts_with("feed"))
                    .collect::<Vec<_>>(),
            )
            .collect();
        let products = outlets
            .into_iter()
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
            energy: [condenser, reboiler],
        }
    }

    fn from_case(case: &Value) -> Self {
        Self::new(
            selection(case),
            spec_for(case),
            case["feeds"].as_array().unwrap(),
        )
    }

    fn set(&mut self, spec: DistillationColumnSpec) {
        self.project.flowsheet.set_spec(
            self.column,
            UnitSpec::Columns(ColumnSpec::Distillation(spec)),
        );
    }

    fn solve(&self) -> Result<ResultSet, SolveFailure> {
        solve_project(
            &self.project,
            &self.project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        )
    }

    fn error(&self) -> String {
        let error = self.solve().unwrap_err();
        assert_eq!(error.unit, Some(self.column), "{}", error.message);
        error.message
    }
}

#[test]
fn partial_condensers_pressure_profiles_and_specifications_match_independent_references() {
    let fixture = fixture();
    let compare = |id: &Value, field: &str, actual: f64, expected: f64| {
        let range = &fixture["tolerances"][field];
        let tolerance = number(&range[0]) + number(&range[1]) * expected.abs();
        assert!(
            (actual - expected).abs() <= tolerance,
            "{id} {field}: {actual} differs from {expected} by more than {tolerance}"
        );
    };
    for case in fixture["cases"].as_array().unwrap() {
        let id = &case["id"];
        let mut rig = Rig::from_case(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{id}: {}", e.message));
        let unit = &result.units[&rig.column];
        let column = unit.distillation_column.as_ref().unwrap();
        let temperatures = case["stage_temperatures_K"].as_array().unwrap();
        assert_eq!(column.stage_temperatures_k.len(), temperatures.len());
        for (actual, expected) in column.stage_temperatures_k.iter().zip(temperatures) {
            compare(id, "temperature_K", *actual, number(expected));
        }
        let pressures: Vec<f64> = case["stage_pressures_Pa"]
            .as_array()
            .unwrap()
            .iter()
            .map(number)
            .collect();
        for (actual, expected) in column.stage_pressures_pa.iter().zip(&pressures) {
            assert!(
                (actual - expected).abs() <= 1e-9 * expected,
                "{id} pressure"
            );
        }
        for (port, duty, field) in [
            (rig.energy[0], column.condenser_duty_w, "condenser_duty_W"),
            (rig.energy[1], column.reboiler_duty_w, "reboiler_duty_W"),
        ] {
            compare(id, "duty_W", duty, number(&case[field]));
            assert_eq!(result.energy_streams[&port].duty_w, duty);
        }
        for (actual, field) in [
            (column.distillate_mol_s, "distillate_mol_s"),
            (column.reflux_mol_s, "reflux_mol_s"),
            (column.bottoms_mol_s, "bottoms_mol_s"),
            (column.boilup_mol_s, "boilup_mol_s"),
        ] {
            compare(
                id,
                "molar_flow_mol_s",
                actual.unwrap(),
                number(&case[field]),
            );
        }
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9, "{id}");
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9, "{id}");
        let outlets = case["outlets"].as_array().unwrap();
        assert_eq!(rig.products.len(), outlets.len(), "{id}");
        for (stream, expected) in rig.products.iter().zip(outlets) {
            let actual = &result.streams[stream];
            compare(
                id,
                "molar_flow_mol_s",
                actual.molar_flow_mol_s,
                number(&expected["molar_flow_mol_s"]),
            );
            let state = actual.state.as_ref().unwrap();
            let vapor = if expected["phase"] == "vapor" {
                1.0
            } else {
                0.0
            };
            assert_eq!(state.vapor_fraction, Some(vapor), "{id}");
            for (actual, expected) in state
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                compare(id, "composition", *actual, number(expected));
            }
        }
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved, "{id}");
        assert_eq!(
            reopened.flowsheet.unit(rig.column).unwrap().spec,
            rig.project.flowsheet.unit(rig.column).unwrap().spec
        );
    }
}

#[test]
fn whole_column_pressure_drop_equals_the_same_drop_per_stage() {
    let case = case("pressure-drop-methanol-water-PR");
    let rig = Rig::from_case(&case);
    let per_stage = rig.solve().unwrap();
    let mut spec = spec_for(&case);
    let stages = f64::from(count(&case["internal_stages"]));
    spec.pressure_drop =
        ColumnPressureDrop::Column(number(&case["stage_pressure_drop_Pa"]) * stages);
    let mut rig = rig;
    rig.set(spec);
    let whole = rig.solve().unwrap();
    let profile = |r: &ResultSet| r.units[&rig.column].distillation_column.clone().unwrap();
    let (a, b) = (profile(&per_stage), profile(&whole));
    for (x, y) in a.stage_pressures_pa.iter().zip(&b.stage_pressures_pa) {
        assert!((x - y).abs() <= 1e-9 * x);
    }
    for (x, y) in a.stage_temperatures_k.iter().zip(&b.stage_temperatures_k) {
        assert!((x - y).abs() <= 1e-8, "{x} {y}");
    }
}

#[test]
fn equivalent_flow_specifications_reproduce_the_reflux_ratio_column() {
    // Reflux rate with bottoms rate, and distillate rate with the reboiler
    // duty of the reference, reach the same column as reflux ratio with
    // distillate fraction.
    let case = case("reboiler-duty-PR");
    let reference = number(&case["reboiler_duty_W"]);
    let mut rig = Rig::from_case(&case);
    let base = rig.solve().unwrap().units[&rig.column]
        .distillation_column
        .clone()
        .unwrap();
    let (d, l, b) = (
        base.distillate_mol_s.unwrap(),
        base.reflux_mol_s.unwrap(),
        base.bottoms_mol_s.unwrap(),
    );
    for (top, bottom) in [
        (
            TopSpecification::RefluxRate(Some(l)),
            BottomSpecification::BottomsRate(Some(b)),
        ),
        (
            TopSpecification::DistillateRate(Some(d)),
            BottomSpecification::ReboilerDuty(Some(reference)),
        ),
        (
            TopSpecification::RefluxRatio(Some(l / d)),
            BottomSpecification::BoilupRatio(Some(base.boilup_mol_s.unwrap() / b)),
        ),
    ] {
        let mut spec = spec_for(&case);
        spec.top = top.clone();
        spec.bottom = bottom.clone();
        rig.set(spec);
        let other = rig
            .solve()
            .unwrap_or_else(|e| panic!("{top:?} {bottom:?}: {}", e.message));
        let other = other.units[&rig.column]
            .distillation_column
            .clone()
            .unwrap();
        for (x, y) in base
            .stage_temperatures_k
            .iter()
            .zip(&other.stage_temperatures_k)
        {
            assert!((x - y).abs() <= 1e-6, "{top:?} {bottom:?}: {x} {y}");
        }
        assert!((other.reboiler_duty_w - reference).abs() <= 1e-7 * reference);
    }
}

#[test]
fn infeasible_specifications_fail_explicitly() {
    let case = case("purity-boilup-PR-BM");
    let mut rig = Rig::from_case(&case);
    let base = spec_for(&case);

    // A distillate purity two stages cannot reach.
    let mut spec = base.clone();
    spec.internals = ColumnInternals::Trays(1);
    spec.feed_stages = vec![2];
    spec.top = TopSpecification::DistillatePurity(ComponentTarget {
        component: Some("ethane".into()),
        value: Some(0.999_999),
    });
    rig.set(spec);
    let message = rig.error();
    assert!(
        message.contains("specifications cannot be met"),
        "{message}"
    );

    // A reboiler duty beyond what boiling up all but a vanishing bottoms
    // flow at this reflux ratio can absorb.
    let mut spec = base.clone();
    spec.top = TopSpecification::RefluxRatio(Some(2.5));
    spec.bottom = BottomSpecification::ReboilerDuty(Some(1.0e9));
    rig.set(spec);
    let message = rig.error();
    assert!(
        message.contains("specifications cannot be met"),
        "{message}"
    );

    // Two specifications that both fix the distillate.
    let mut spec = base.clone();
    spec.top = TopSpecification::DistillateRate(Some(40.0));
    spec.bottom = BottomSpecification::BottomsRate(Some(60.0));
    rig.set(spec.clone());
    assert!(
        rig.project
            .flowsheet
            .unit(rig.column)
            .unwrap()
            .issues()
            .iter()
            .any(|issue| issue.field == "bottom specification")
    );
    // The flowsheet rejects it before calculation.
    let message = rig.solve().unwrap_err().message;
    assert!(
        message.contains("bottom specification duplicates another selection"),
        "{message}"
    );

    // A purity component outside the package, and one outside the feed.
    let methane = TopSpecification::DistillatePurity(ComponentTarget {
        component: Some("methane".into()),
        value: Some(0.9),
    });
    let mut spec = base.clone();
    spec.top = methane;
    rig.set(spec.clone());
    assert!(rig.error().contains("is not a package component"));
    let mut ternary = case.clone();
    ternary["ids"] = json!(["methane", "ethane", "propane"]);
    ternary["kij"] = json!([[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
    ternary["feeds"][0]["z"] = json!([0.0, 0.45, 0.55]);
    let absent = Rig::new(
        selection(&ternary),
        spec,
        ternary["feeds"].as_array().unwrap(),
    );
    assert!(absent.error().contains("is not in the column feed"));

    // A liquid distillate that takes the whole distillate.
    let mut spec = base.clone();
    spec.top = TopSpecification::RefluxRatio(Some(2.5));
    spec.bottom = BottomSpecification::DistillateFraction(Some(0.3));
    spec.condenser = DistillationCondenser::PartialVaporLiquid {
        liquid_distillate_mol_s: Some(30.0),
    };
    let partial = Rig::new(selection(&case), spec, case["feeds"].as_array().unwrap());
    assert!(
        partial
            .error()
            .contains("smaller than the total distillate")
    );

    // A feed below its stage pressure after the pressure drop.
    let mut spec = base;
    spec.top = TopSpecification::RefluxRatio(Some(2.5));
    spec.bottom = BottomSpecification::DistillateFraction(Some(0.45));
    spec.condenser_pressure_drop_pa = 4.0e5;
    rig.set(spec);
    assert!(rig.error().contains("exceeds feed 1 pressure"));
}

#[test]
fn new_fields_round_trip_and_legacy_records_load_unchanged() {
    let case = case("partial-vapor-liquid-vapor-draw-PR-BM");
    let mut rig = Rig::from_case(&case);
    let mut spec = spec_for(&case);
    spec.condenser_pressure_drop_pa = 1500.0;
    spec.pressure_drop = ColumnPressureDrop::Column(4000.0);
    spec.top = TopSpecification::DistillateRecovery(ComponentTarget {
        component: Some("ethane".into()),
        value: Some(0.9),
    });
    spec.bottom = BottomSpecification::BottomsRate(Some(52.0));
    rig.set(spec.clone());
    let saved: Value = serde_json::from_str(&rig.project.to_json().unwrap()).unwrap();
    let record = saved["units"]
        .as_array()
        .unwrap()
        .iter()
        .find(|u| u["kind"] == "distillation_column")
        .unwrap()["spec"]["unit"]
        .clone();
    assert_eq!(
        record,
        json!({
            "kind": "distillation_column",
            "feeds": 1,
            "side_draws": 0,
            "feed_stages": [4],
            "pressure_pa": 2.0e6,
            "top_specification": {"distillate_recovery": {"component": "ethane", "value": 0.9}},
            "bottom_specification": {"bottoms_rate_mol_s": 52.0},
            "condenser": {"partial_vapor_liquid": {"liquid_distillate_mol_s": 20.0}},
            "condenser_pressure_drop_pa": 1500.0,
            "pressure_drop": {"column_pa": 4000.0},
            "vapor_side_draws": 1,
            "vapor_draws": [{"stage": 5, "molar_flow_mol_s": 8.0}],
        })
    );
    let reopened = Project::from_json(&saved.to_string()).unwrap();
    assert_eq!(
        reopened.flowsheet.unit(rig.column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::Distillation(spec))
    );
    assert_eq!(
        reopened.flowsheet.unit(rig.column).unwrap().ports(),
        rig.project.flowsheet.unit(rig.column).unwrap().ports()
    );

    // A version-4 record written before the specification set loads as a
    // total condenser at uniform pressure with reflux ratio and distillate
    // fraction, and writes back the same fields.
    let legacy = json!({
        "kind": "distillation_column", "feeds": 1, "side_draws": 0,
        "pressure_pa": 2.0e6, "reflux_ratio": 2.5, "distillate_fraction": 0.45,
    });
    let mut plain = Project::new(selection(&case));
    let column = plain.flowsheet.add_unit(UnitKind::DistillationColumn);
    let mut plain_json: Value = serde_json::from_str(&plain.to_json().unwrap()).unwrap();
    plain_json["units"][0]["spec"]["unit"] = legacy.clone();
    plain = Project::from_json(&plain_json.to_string()).unwrap();
    let UnitSpec::Columns(ColumnSpec::Distillation(loaded)) =
        &plain.flowsheet.unit(column).unwrap().spec
    else {
        unreachable!()
    };
    assert_eq!(
        *loaded,
        DistillationColumnSpec {
            pressure_pa: Some(2.0e6),
            top: TopSpecification::RefluxRatio(Some(2.5)),
            bottom: BottomSpecification::DistillateFraction(Some(0.45)),
            ..Default::default()
        }
    );
    let rewritten: Value = serde_json::from_str(&plain.to_json().unwrap()).unwrap();
    assert_eq!(rewritten["units"][0]["spec"]["unit"], legacy);

    // A record with both the original field and a new specification is rejected.
    let mut both = plain_json.clone();
    both["units"][0]["spec"]["unit"]["top_specification"] = json!({"reflux_rate_mol_s": 10.0});
    assert!(Project::from_json(&both.to_string()).is_err());
}
