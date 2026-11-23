//! Shortcut column (#185) through normal project solving and v4 persistence.

use std::ops::ControlFlow;

use flowsheet_core::units::ColumnSpec;
use flowsheet_core::units::columns::{ShortcutColumnSpec, ShortcutReflux};
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
    }
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

    line.project.results = Some(results.clone());
    let saved = line.project.to_json().unwrap();
    assert!(saved.contains("\"shortcut_column\""));
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
