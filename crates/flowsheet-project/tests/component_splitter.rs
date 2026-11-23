//! Project-format and flowsheet-solve coverage for #172.
use std::ops::ControlFlow;

use flowsheet_core::units::SimpleLogicalSpec;
use flowsheet_core::units::simple_logical::ComponentSplitterSpec;
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{CubicVariantChoice, PackageSelection, Project, solve_project};

struct SplitterProject {
    project: Project,
    splitter: UnitId,
    inlet: StreamId,
    overhead: StreamId,
    bottoms: StreamId,
    energy: StreamId,
}

fn splitter_project() -> SplitterProject {
    let mut project = Project::new(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        vec!["propane".into(), "methane".into()],
        "explicit zero-interaction approximation for the binary splitter test",
    ));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(250.0),
            pressure_pa: Some(3e6),
            molar_flow_mol_s: Some(10.0),
        },
    );
    project.feeds.insert(feed, vec![0.5, 0.5]);
    let splitter = fs.add_unit(UnitKind::ComponentSplitter);
    let overhead_product = fs.add_unit(UnitKind::Product);
    let bottoms_product = fs.add_unit(UnitKind::Product);
    let energy_source = fs.add_unit(UnitKind::EnergySource);
    let inlet = fs
        .connect(PortRef::new(feed, 0), PortRef::new(splitter, 0))
        .unwrap();
    let overhead = fs
        .connect(PortRef::new(splitter, 1), PortRef::new(overhead_product, 0))
        .unwrap();
    let bottoms = fs
        .connect(PortRef::new(splitter, 2), PortRef::new(bottoms_product, 0))
        .unwrap();
    let energy = fs
        .connect(PortRef::new(energy_source, 0), PortRef::new(splitter, 3))
        .unwrap();
    SplitterProject {
        project,
        splitter,
        inlet,
        overhead,
        bottoms,
        energy,
    }
}

fn spec<const N: usize>(recoveries: [(&str, f64); N]) -> UnitSpec {
    UnitSpec::SimpleLogical(SimpleLogicalSpec::ComponentSplitter(
        ComponentSplitterSpec {
            overhead_recovery: recoveries
                .into_iter()
                .map(|(id, recovery)| (id.into(), recovery))
                .collect(),
        },
    ))
}

#[test]
fn loaded_component_splitter_closes_component_and_energy_balances() {
    let mut case = splitter_project();
    case.project
        .flowsheet
        .set_spec(case.splitter, spec([("propane", 0.2), ("methane", 0.8)]));
    let json = case.project.to_json().unwrap();
    let project = Project::from_json(&json).unwrap();
    assert_eq!(project.to_json().unwrap(), json);
    let package = project.package.prepare().unwrap();
    let result = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    let inlet = &result.streams[&case.inlet];
    let overhead = &result.streams[&case.overhead];
    let bottoms = &result.streams[&case.bottoms];
    let inlet_state = inlet.state.as_ref().unwrap();
    for (stream, expected) in [(overhead, [0.2, 0.8]), (bottoms, [0.8, 0.2])] {
        assert!((stream.molar_flow_mol_s - 5.0).abs() < 1e-12);
        let state = stream.state.as_ref().unwrap();
        assert_eq!(state.temperature_k, inlet_state.temperature_k);
        assert_eq!(state.pressure_pa, inlet_state.pressure_pa);
        for (&actual, expected) in state.mole_fractions.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }
    let unit = &result.units[&case.splitter];
    assert!(unit.max_scaled_component_residual.unwrap() < 1e-12);
    assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-12);
    assert!(unit.duty_w.abs() > 1.0);
    assert_eq!(result.energy_streams[&case.energy].duty_w, unit.duty_w);
}

#[test]
fn complete_recovery_leaves_the_other_outlet_empty() {
    let mut case = splitter_project();
    case.project
        .flowsheet
        .set_spec(case.splitter, spec([("propane", 1.0), ("methane", 1.0)]));
    let package = case.project.package.prepare().unwrap();
    let r = solve_project(&case.project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    assert_eq!(
        r.streams[&case.overhead].state,
        r.streams[&case.inlet].state
    );
    assert_eq!(r.streams[&case.bottoms].molar_flow_mol_s, 0.0);
    assert_eq!(r.streams[&case.bottoms].state, None);
    assert_eq!(r.units[&case.splitter].duty_w, 0.0);

    let feed = case.project.flowsheet.stream(case.inlet).unwrap().from.unit;
    case.project.flowsheet.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(250.0),
            pressure_pa: Some(3e6),
            molar_flow_mol_s: Some(0.0),
        },
    );
    let r = solve_project(&case.project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    for outlet in [case.overhead, case.bottoms] {
        assert_eq!(r.streams[&outlet].molar_flow_mol_s, 0.0);
        assert_eq!(r.streams[&outlet].state, None);
    }
    assert_eq!(r.energy_streams[&case.energy].duty_w, 0.0);
}

#[test]
fn original_v4_record_and_invalid_recoveries_are_reported_at_the_splitter() {
    let mut case = splitter_project();
    // The foundation's v4 record had no recovery field.
    let json = case.project.to_json().unwrap();
    assert!(!json.contains("overhead_recovery"));
    case.project = Project::from_json(&json).unwrap();
    let package = case.project.package.prepare().unwrap();
    let solve =
        |project: &Project| solve_project(project, &package, |_| ControlFlow::Continue(()), 0);
    let err = solve(&case.project).unwrap_err();
    assert_eq!(err.unit, Some(case.splitter));
    assert!(
        err.message
            .contains("missing overhead recovery for propane")
    );

    case.project.flowsheet.set_spec(
        case.splitter,
        spec([("propane", 0.5), ("methane", 0.5), ("nitrogen", 0.5)]),
    );
    let err = solve(&case.project).unwrap_err();
    assert_eq!(err.unit, Some(case.splitter));
    assert!(
        err.message
            .contains("nitrogen is outside the property package")
    );

    case.project
        .flowsheet
        .set_spec(case.splitter, spec([("propane", -0.2), ("methane", 0.5)]));
    let err = solve(&case.project).unwrap_err();
    let tag = &case.project.flowsheet.unit(case.splitter).unwrap().name;
    assert!(err.message.contains(tag), "{}", err.message);
    assert!(err.message.contains("overhead recovery"));
}
