//! Shortcut column (#185) through normal project solving and v4 persistence.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{ShortcutColumnSpec, ShortcutCondenser, ShortcutReflux};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::{Value, json};

struct ColumnProject {
    project: Project,
    column: UnitId,
    distillate: StreamId,
    bottoms: StreamId,
}

fn spec() -> ShortcutColumnSpec {
    ShortcutColumnSpec {
        light_key: Some("ethane".into()),
        heavy_key: Some("propane".into()),
        light_key_recovery: Some(0.98),
        heavy_key_recovery: Some(0.97),
        pressure_pa: Some(2.0e6),
        reflux: ShortcutReflux::MinimumMultiple(Some(1.3)),
        condenser: ShortcutCondenser::Total,
    }
}

fn shortcut_record(value: &mut Value) -> &mut Value {
    &mut value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "shortcut_column")
        .unwrap()["spec"]["unit"]
}

impl ColumnProject {
    fn new() -> Self {
        let mut project = Project::new(PackageSelection::zero_interactions(
            CubicVariantChoice::PengRobinsonBostonMathias,
            vec!["ethane".into(), "propane".into()],
            "explicit zero interactions for the shortcut column project test",
        ));
        let fs = &mut project.flowsheet;
        let feed = fs.add_unit(UnitKind::Feed);
        let column = fs.add_unit(UnitKind::ShortcutColumn);
        let top = fs.add_unit(UnitKind::Product);
        let bottom = fs.add_unit(UnitKind::Product);
        fs.set_spec(
            feed,
            UnitSpec::Feed {
                temperature_k: Some(295.0),
                pressure_pa: Some(2.0e6),
                molar_flow_mol_s: Some(100.0),
            },
        );
        fs.try_set_spec(column, UnitSpec::Columns(ColumnSpec::Shortcut(spec())))
            .unwrap();
        project.feeds.insert(feed, vec![0.4, 0.6]);
        fs.connect(PortRef::new(feed, 0), PortRef::new(column, 0))
            .unwrap();
        let distillate = fs
            .connect(PortRef::new(column, 1), PortRef::new(top, 0))
            .unwrap();
        let bottoms = fs
            .connect(PortRef::new(column, 2), PortRef::new(bottom, 0))
            .unwrap();
        Self {
            project,
            column,
            distillate,
            bottoms,
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
fn shortcut_specification_and_estimates_survive_a_v4_round_trip() {
    let mut line = ColumnProject::new();
    let results = line.solve().unwrap();
    let estimate = results.units[&line.column].shortcut_column.clone().unwrap();
    // Reference case c2c3-binary-PR-BM of validation/columns/shortcut-v1.json.
    assert!((estimate.minimum_stages - 7.694_873_089_732_829).abs() < 1e-5);
    assert!((estimate.theoretical_stages - 16.673_548_499_478_848).abs() < 2e-5);
    assert!((estimate.condenser_duty_w - 1_147_736.584_977_614).abs() < 1.2);
    let d = results.streams[&line.distillate].molar_flow_mol_s;
    let b = results.streams[&line.bottoms].molar_flow_mol_s;
    assert!((d - 41.0).abs() < 1e-9 && (b - 59.0).abs() < 1e-9);
    let recoveries: Vec<f64> = estimate
        .distillate_recoveries
        .iter()
        .flatten()
        .copied()
        .collect();
    assert!((recoveries[0] - 0.98).abs() < 1e-12 && (recoveries[1] - 0.03).abs() < 1e-12);

    line.project.results = Some(results.clone());
    let saved = line.project.to_json().unwrap();
    assert!(saved.contains("\"shortcut_column\""));
    // A total condenser is the default and is not written.
    assert!(!saved.contains("\"condenser\""));
    let reopened = Project::from_json(&saved).unwrap();
    assert_eq!(reopened.to_json().unwrap(), saved);
    assert_eq!(
        reopened.flowsheet.unit(line.column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::Shortcut(spec()))
    );
    assert_eq!(
        reopened.results.as_ref().unwrap().units[&line.column].shortcut_column,
        Some(estimate)
    );
    // Recalculating the reopened project reproduces the stored estimates.
    line.project = reopened;
    let again = line.solve().unwrap();
    assert_eq!(again.units[&line.column], results.units[&line.column]);
}

#[test]
fn partial_condenser_and_stage_count_survive_a_v4_round_trip() {
    let mut line = ColumnProject::new();
    // Reference case c2c3-partial-PR-BM of validation/columns/shortcut-v1.json.
    let mut partial = spec();
    partial.condenser = ShortcutCondenser::Partial;
    line.project
        .flowsheet
        .try_set_spec(
            line.column,
            UnitSpec::Columns(ColumnSpec::Shortcut(partial.clone())),
        )
        .unwrap();
    let results = line.solve().unwrap();
    let estimate = results.units[&line.column].shortcut_column.clone().unwrap();
    assert!((estimate.condenser_duty_w - 774_150.458_2).abs() < 1.0);
    let top = results.streams[&line.distillate].state.as_ref().unwrap();
    assert!((top.temperature_k - 271.1137).abs() < 1e-3);

    // Stage count instead of reflux: the estimate's own N reproduces R.
    let mut staged = partial.clone();
    staged.reflux = ShortcutReflux::TheoreticalStages(Some(estimate.theoretical_stages));
    line.project
        .flowsheet
        .try_set_spec(
            line.column,
            UnitSpec::Columns(ColumnSpec::Shortcut(staged.clone())),
        )
        .unwrap();
    let results = line.solve().unwrap();
    let by_stages = results.units[&line.column].shortcut_column.clone().unwrap();
    assert!((by_stages.reflux_ratio / estimate.reflux_ratio - 1.0).abs() < 1e-9);

    line.project.results = Some(results.clone());
    let saved = line.project.to_json().unwrap();
    let mut value: Value = serde_json::from_str(&saved).unwrap();
    assert_eq!(value["format_version"], 4);
    let record = shortcut_record(&mut value);
    assert_eq!(record["condenser"], "partial");
    assert_eq!(
        record["reflux"],
        json!({"theoretical_stages": estimate.theoretical_stages})
    );
    let reopened = Project::from_json(&saved).unwrap();
    assert_eq!(reopened.to_json().unwrap(), saved);
    assert_eq!(
        reopened.flowsheet.unit(line.column).unwrap().spec,
        UnitSpec::Columns(ColumnSpec::Shortcut(staged))
    );
    assert_eq!(
        reopened.results.as_ref().unwrap().units[&line.column].shortcut_column,
        Some(by_stages)
    );

    // An unknown condenser or operating choice is refused, not defaulted.
    for (field, bad) in [
        ("condenser", json!("refluxed")),
        ("reflux", json!({"stage_count": 12.0})),
    ] {
        let mut value: Value = serde_json::from_str(&saved).unwrap();
        shortcut_record(&mut value)[field] = bad;
        assert!(Project::from_json(&value.to_string()).is_err(), "{field}");
    }
}

#[test]
fn shortcut_estimate_saved_before_the_distribution_loads() {
    let line = ColumnProject::new();
    let mut project = line.project.clone();
    project.results = Some(line.solve().unwrap());
    let mut value: Value = serde_json::from_str(&project.to_json().unwrap()).unwrap();
    fn strip(v: &mut Value) -> usize {
        match v {
            Value::Object(map) => {
                let here = map.remove("distillate_recoveries").is_some() as usize;
                here + map.values_mut().map(strip).sum::<usize>()
            }
            Value::Array(items) => items.iter_mut().map(strip).sum(),
            _ => 0,
        }
    }
    assert_eq!(strip(&mut value["results"]), 1);
    let reopened = Project::from_json(&value.to_string()).unwrap();
    let estimate = reopened.results.unwrap().units[&line.column]
        .shortcut_column
        .clone()
        .unwrap();
    assert!(estimate.distillate_recoveries.is_empty());
}

#[test]
fn legacy_empty_shortcut_record_loads_unset_and_reports_the_missing_inputs() {
    let mut line = ColumnProject::new();
    let mut value: Value = serde_json::from_str(&line.project.to_json().unwrap()).unwrap();
    assert_eq!(value["format_version"], 4);
    let record = &mut value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "shortcut_column")
        .unwrap()["spec"]["unit"];
    assert_eq!(record["reflux"], json!({"minimum_multiple": 1.3}));
    *record = json!({"kind": "shortcut_column"});
    line.project = Project::from_json(&value.to_string()).unwrap();
    assert_eq!(
        line.project.flowsheet.unit(line.column).unwrap().spec,
        UnitKind::ShortcutColumn.default_spec()
    );
    assert!(line.solve().is_err());

    // Unknown fields are refused rather than dropped.
    let mut value: Value = serde_json::from_str(&line.project.to_json().unwrap()).unwrap();
    value["units"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|u| u["kind"] == "shortcut_column")
        .unwrap()["spec"]["unit"]["relative_volatility"] = json!(2.0);
    assert!(Project::from_json(&value.to_string()).is_err());
}
