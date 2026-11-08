use flowsheet_core::{PortRef, ThermalSpec, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, Currency, Project, recycle_demonstration, solve_project,
};
use std::ops::ControlFlow;

fn solve(project: &Project) -> flowsheet_project::ResultSet {
    solve_project(
        project,
        &project.package.build().unwrap(),
        |_| ControlFlow::Continue(()),
        0,
    )
    .unwrap()
}

#[test]
fn cyclic_results_survive_reload_and_edits_invalidate_the_loop_and_downstream() {
    let mut project = recycle_demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    project.results = Some(solve(&project));
    assert_eq!(project.results.as_ref().unwrap().units.len(), 5);
    assert_eq!(project.results.as_ref().unwrap().streams.len(), 5);
    let mut loaded = Project::from_json(&project.to_json().unwrap()).unwrap();
    let fps = loaded.fingerprints();
    let results = loaded.results.as_ref().unwrap();
    assert!(results.recycle.as_ref().unwrap().max_scaled_residual <= 1e-10);
    for u in loaded.flowsheet.units() {
        assert_eq!(results.unit_currency(u.id, &fps), Currency::Current);
    }
    for s in loaded.flowsheet.streams() {
        assert_eq!(
            results.stream_currency(&loaded.flowsheet, s.id, &fps),
            Currency::Current
        );
    }
    let heater = loaded
        .flowsheet
        .units()
        .find(|u| u.kind == UnitKind::Heater)
        .unwrap()
        .id;
    loaded.flowsheet.set_spec(
        heater,
        UnitSpec::HeatExchange {
            thermal: ThermalSpec::OutletTemperature(Some(360.0)),
            pressure_drop_pa: 0.0,
        },
    );
    let fps = loaded.fingerprints();
    for u in loaded.flowsheet.units() {
        assert_eq!(
            results.unit_currency(u.id, &fps),
            if u.kind == UnitKind::Feed {
                Currency::Current
            } else {
                Currency::Stale
            }
        );
    }
    let recalculated = solve(&loaded);
    assert!(
        recalculated
            .units
            .keys()
            .all(|id| recalculated.unit_currency(*id, &fps) == Currency::Current)
    );
    let stream = loaded
        .flowsheet
        .streams()
        .find(|s| s.name == "recycle")
        .unwrap();
    assert!(
        (recalculated.streams[&stream.id]
            .state
            .as_ref()
            .unwrap()
            .temperature_k
            - 360.0)
            .abs()
            < 1e-8
    );
}

#[test]
fn cyclic_fingerprints_cover_feed_composition_connections_and_keep_other_branches_current() {
    let mut project = recycle_demonstration(CubicVariantChoice::PengRobinson);
    let feed = project.flowsheet.add_unit(UnitKind::Feed);
    let product = project.flowsheet.add_unit(UnitKind::Product);
    project.flowsheet.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(10.0),
        },
    );
    project
        .flowsheet
        .connect(PortRef::new(feed, 0), PortRef::new(product, 0))
        .unwrap();
    project
        .feeds
        .insert(feed, project.feeds.values().next().unwrap().clone());
    let before = project.fingerprints();
    let first_feed = *project.feeds.keys().next().unwrap();
    project
        .feeds
        .insert(first_feed, vec![0.70, 0.15, 0.13, 0.01, 0.01]);
    let after = project.fingerprints();
    for u in project.flowsheet.units() {
        assert_eq!(
            before.units[&u.id] == after.units[&u.id],
            u.id == feed || u.id == product
        );
    }
    let recycle = project
        .flowsheet
        .streams()
        .find(|s| s.name == "recycle")
        .unwrap()
        .clone();
    let fresh = project
        .flowsheet
        .streams()
        .find(|s| s.name == "fresh")
        .unwrap()
        .clone();
    project.flowsheet.disconnect(recycle.id);
    project.flowsheet.disconnect(fresh.id);
    project.flowsheet.connect(recycle.from, fresh.to).unwrap();
    project.flowsheet.connect(fresh.from, recycle.to).unwrap();
    let changed = project.fingerprints();
    assert_ne!(after.combined(), changed.combined());
    let mixer = project
        .flowsheet
        .units()
        .find(|u| u.kind == UnitKind::Mixer)
        .unwrap()
        .id;
    assert_ne!(after.units[&mixer], changed.units[&mixer]);
    assert_eq!(after.units[&product], changed.units[&product]);
}

#[test]
fn project_cancellation_and_failed_loop_return_no_new_results() {
    let project = recycle_demonstration(CubicVariantChoice::PengRobinson);
    let pkg = project.package.build().unwrap();
    let error = solve_project(
        &project,
        &pkg,
        |p| {
            if p.iteration >= 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        0,
    )
    .unwrap_err();
    assert!(error.cancelled);
    let mut project = project;
    project.results = Some(solve(&project));
    let old_results = project.results.clone();
    let splitter = project
        .flowsheet
        .units()
        .find(|u| u.kind == UnitKind::Splitter)
        .unwrap()
        .id;
    project.flowsheet.set_spec(
        splitter,
        UnitSpec::Splitter {
            fraction_to_first: Some(1.0),
        },
    );
    let error = solve_project(&project, &pkg, |_| ControlFlow::Continue(()), 0).unwrap_err();
    assert!(!error.cancelled);
    assert!(error.message.contains("did not converge after 500"));
    assert_eq!(project.results, old_results);
    let fps = project.fingerprints();
    assert_eq!(
        project
            .results
            .as_ref()
            .unwrap()
            .unit_currency(splitter, &fps),
        Currency::Stale
    );
}
