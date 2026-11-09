//! Persistence, invalidation and recalculation of projects (A9, #54), and the
//! demonstration project against the independent A2 process fixture.

use std::ops::ControlFlow;
use std::path::PathBuf;

use flowsheet_core::{PortRef, PressureSpec, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, Currency, Fingerprint, Fingerprints, Project, ProjectError, ResultSet,
    demonstration, run_acceptance, solve_project,
};

fn solve(project: &Project) -> ResultSet {
    let package = project.package.prepare().unwrap();
    solve_project(project, &package, |_| ControlFlow::Continue(()), 0).unwrap()
}

fn unit(project: &Project, kind: UnitKind) -> UnitId {
    project
        .flowsheet
        .units()
        .find(|u| u.kind == kind)
        .unwrap()
        .id
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("flowsheet-project-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn demonstration_project_reproduces_the_independent_process_fixture() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let variant = if case["variant"] == "PR" {
            CubicVariantChoice::PengRobinson
        } else {
            CubicVariantChoice::PengRobinsonBostonMathias
        };
        let close = |key: &str, actual: f64, expected: f64| {
            let budget = &fixture["tolerances"][key];
            assert!(
                (actual - expected).abs()
                    <= budget[0].as_f64().unwrap() + budget[1].as_f64().unwrap() * expected.abs(),
                "{} {key}: {actual} != {expected}",
                case["variant"]
            );
        };
        let project = demonstration(variant);
        let r = solve(&project);
        for expected in case["streams"].as_array().unwrap() {
            let name = expected["id"].as_str().unwrap();
            let s = project
                .flowsheet
                .streams()
                .find(|s| s.name == name)
                .unwrap();
            let result = &r.streams[&s.id];
            let state = result.state.as_ref().unwrap();
            close(
                "temperature_K",
                state.temperature_k,
                expected["T"].as_f64().unwrap(),
            );
            assert_eq!(state.pressure_pa, expected["p"].as_f64().unwrap());
            close(
                "enthalpy_J_per_mol",
                state.enthalpy_j_mol,
                expected["target_h"].as_f64().unwrap(),
            );
            close(
                "phase_fraction",
                result.molar_flow_mol_s / 120.0,
                expected["flow_mol_per_s"].as_f64().unwrap() / 120.0,
            );
            for (x, z) in state
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                close("composition", *x, z.as_f64().unwrap());
            }
            if state.phases.len() == 2 {
                close(
                    "phase_fraction",
                    state.vapor_fraction.unwrap(),
                    expected["vapor_fraction"].as_f64().unwrap(),
                );
            }
        }
        let cooler = unit(&project, UnitKind::Cooler);
        close(
            "enthalpy_J_per_mol",
            r.units[&cooler].duty_w / 120.0,
            case["cooler_duty_W"].as_f64().unwrap() / 120.0,
        );
    }
}

#[test]
fn acceptance_workflow_passes() {
    let dir = scratch("acceptance");
    let report = run_acceptance(&dir).unwrap();
    assert_eq!(report.steps.len(), 7, "{:#?}", report.steps);
    assert!(report.max_reload_relative_difference <= 1e-9);
    assert!(report.separator_temperature_k[1] > report.separator_temperature_k[0]);
}

#[test]
fn save_reload_preserves_inputs_selection_layout_and_results() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.results = Some(solve(&project));
    project.display_units = flowsheet_project::DisplayUnitSet::Si;
    let text = project.to_json().unwrap();
    let loaded = Project::from_json(&text).unwrap();
    assert_eq!(loaded.package, project.package);
    assert_eq!(loaded.feeds, project.feeds);
    assert_eq!(loaded.layout, project.layout);
    assert_eq!(loaded.display_units, project.display_units);
    assert_eq!(loaded.results, project.results);
    let a: Vec<_> = project.flowsheet.units().cloned().collect();
    let b: Vec<_> = loaded.flowsheet.units().cloned().collect();
    assert_eq!(a, b);
    let a: Vec<_> = project.flowsheet.streams().cloned().collect();
    let b: Vec<_> = loaded.flowsheet.streams().cloned().collect();
    assert_eq!(a, b);
    assert_eq!(loaded.fingerprints(), project.fingerprints());
    assert_eq!(loaded.to_json().unwrap(), text);

    // Through the file system as well.
    let path = scratch("roundtrip").join("p.flowsheet.json");
    project.save(&path).unwrap();
    assert_eq!(Project::load(&path).unwrap().to_json().unwrap(), text);
}

#[test]
fn edits_invalidate_exactly_the_affected_results() {
    let mut project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let r = solve(&project);
    let fs = &project.flowsheet;
    let feed2 = fs
        .units()
        .filter(|u| u.kind == UnitKind::Feed)
        .nth(1)
        .unwrap()
        .id;
    let feed1 = unit(&project, UnitKind::Feed);

    // Names are not inputs.
    let stream = project.flowsheet.streams().next().unwrap().id;
    project.flowsheet.rename_stream(stream, "renamed");
    project.flowsheet.rename_unit(feed1, "Natural gas");
    let fps = project.fingerprints();
    assert!(
        project
            .flowsheet
            .units()
            .all(|u| r.unit_currency(u.id, &fps) == Currency::Current)
    );

    // A feed composition edit invalidates that feed and everything downstream.
    project.feeds.get_mut(&feed2).unwrap()[0] = 0.69;
    project.feeds.get_mut(&feed2).unwrap()[1] = 0.16;
    let fps = project.fingerprints();
    assert_eq!(r.unit_currency(feed1, &fps), Currency::Current);
    for u in project.flowsheet.units().filter(|u| u.id != feed1) {
        assert_eq!(r.unit_currency(u.id, &fps), Currency::Stale, "{}", u.name);
    }
    let a = project
        .flowsheet
        .stream_at(PortRef::new(feed1, 0))
        .unwrap()
        .id;
    assert_eq!(
        r.stream_currency(&project.flowsheet, a, &fps),
        Currency::Current
    );

    // A package change invalidates everything.
    let mut other = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let r = solve(&other);
    other.package.variant = CubicVariantChoice::PengRobinson;
    let fps = other.fingerprints();
    assert!(
        other
            .flowsheet
            .units()
            .all(|u| r.unit_currency(u.id, &fps) == Currency::Stale)
    );

    // A disconnected stream leaves no result for the downstream unit.
    let mut cut = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let r = solve(&cut);
    let sep = unit(&cut, UnitKind::Separator);
    let into_sep = cut.flowsheet.stream_at(PortRef::new(sep, 0)).unwrap().id;
    cut.flowsheet.disconnect(into_sep);
    let fps = cut.fingerprints();
    assert_eq!(r.unit_currency(sep, &fps), Currency::Stale);
    assert_eq!(
        r.stream_currency(&cut.flowsheet, into_sep, &fps),
        Currency::Missing
    );
}

#[test]
fn results_of_another_package_or_provenance_are_stale() {
    let project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let r = solve(&project);
    let own = project.fingerprints();
    assert_eq!(own.package, Fingerprints::package(&project.package));
    // The package fingerprint covers the selection and the provenance
    // (engine revision, model and data versions) of the running program, so
    // results produced under any other package fingerprint are stale.
    let foreign =
        Fingerprints::compute_with(&project.flowsheet, &project.feeds, Some(Fingerprint(1)));
    let unbuildable = Fingerprints::compute_with(&project.flowsheet, &project.feeds, None);
    for fps in [&foreign, &unbuildable] {
        assert!(
            project
                .flowsheet
                .units()
                .all(|u| r.unit_currency(u.id, fps) == Currency::Stale)
        );
        assert!(project.flowsheet.streams().all(|s| r.stream_currency(
            &project.flowsheet,
            s.id,
            fps
        ) == Currency::Stale));
    }
    assert_ne!(own.combined(), foreign.combined());
}

#[test]
fn failed_and_cancelled_solves_return_no_results() {
    let mut project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let package = project.package.prepare().unwrap();
    let valve = unit(&project, UnitKind::Valve);
    project.flowsheet.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(9.0e6)),
        },
    );
    let f = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap_err();
    assert_eq!(f.unit, Some(valve));
    assert!(!f.cancelled);
    assert!(f.message.contains("V-100"), "{}", f.message);

    let mut calls = 0;
    let f = solve_project(
        &project,
        &package,
        |p| {
            calls += 1;
            assert_eq!(p.total, 8);
            ControlFlow::Break(())
        },
        0,
    )
    .unwrap_err();
    assert!(f.cancelled);
    assert_eq!(calls, 1);

    // Missing composition is reported before any calculation.
    let feed = unit(&project, UnitKind::Feed);
    project.feeds.remove(&feed);
    let f = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap_err();
    assert_eq!(f.unit, Some(feed));
    assert!(f.message.contains("missing composition"));
    assert_eq!(project.input_issues().len(), 1);
    project.feeds.insert(feed, vec![0.5, 0.5, 0.0, 0.0, 0.1]);
    assert!(project.input_issues()[0].message.contains("normalize"));
}

#[test]
fn malformed_and_foreign_files_are_rejected() {
    let project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let text = project.to_json().unwrap();
    let edit = |from: &str, to: &str| {
        assert!(text.contains(from), "{from}");
        Project::from_json(&text.replacen(from, to, 1))
    };
    assert_eq!(
        edit("\"format_version\": 2", "\"format_version\": 3").unwrap_err(),
        ProjectError::UnsupportedVersion(Some(3))
    );
    assert_eq!(
        edit("\"format_version\": 2", "\"format_version\": 0").unwrap_err(),
        ProjectError::UnsupportedVersion(Some(0))
    );
    assert_eq!(
        edit("\"format_version\": 2", "\"format\": 2").unwrap_err(),
        ProjectError::UnsupportedVersion(None)
    );
    assert!(matches!(
        edit("\"generator\"", "\"unexpected\": 0, \"generator\""),
        Err(ProjectError::Format(_))
    ));
    assert!(matches!(
        edit("0.88,", "0.88, 0.0,"),
        Err(ProjectError::Invalid(_))
    ));
    assert!(matches!(
        Project::from_json("not json"),
        Err(ProjectError::Format(_))
    ));
    assert!(matches!(
        Project::load(std::path::Path::new("does/not/exist.json")),
        Err(ProjectError::Io(_))
    ));
}

#[test]
fn reloaded_projects_continue_identifiers_and_tags() {
    let mut project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    project.results = Some(solve(&project));
    // Delete the last stream and unit, then save: no identifier may be reused
    // for a result that belonged to the deleted items.
    let liquid = project
        .flowsheet
        .units()
        .find(|u| u.name == "Liquid")
        .unwrap()
        .id;
    project.flowsheet.remove_unit(liquid);
    let text = project.to_json().unwrap();
    let mut loaded = Project::from_json(&text).unwrap();
    assert!(!loaded.results.as_ref().unwrap().units.contains_key(&liquid));
    let p = loaded.flowsheet.add_unit(UnitKind::Product);
    assert_eq!(loaded.flowsheet.unit(p).unwrap().name, "Product 2");
    let m = loaded.flowsheet.add_unit(UnitKind::Mixer);
    assert_eq!(loaded.flowsheet.unit(m).unwrap().name, "M-101");
    let fps: Fingerprints = loaded.fingerprints();
    assert_eq!(
        loaded.results.as_ref().unwrap().unit_currency(p, &fps),
        Currency::Missing
    );
}

#[test]
fn version_2_is_written_and_version_1_files_still_load() {
    assert_eq!(flowsheet_project::FORMAT_VERSION, 2);
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.results = Some(solve(&project));
    let text = project.to_json().unwrap();
    assert!(text.contains("\"format_version\": 2"));

    // A version-1 file (same schema) is read and written back as version 2.
    let v1 = text.replacen("\"format_version\": 2", "\"format_version\": 1", 1);
    let loaded = Project::from_json(&v1).unwrap();
    assert_eq!(loaded.results, project.results);
    assert_eq!(loaded.to_json().unwrap(), text);

    // The committed A10 example is a version-1 input-only project.
    let example = include_str!("../../../examples/recycle.flowsheet.json");
    assert!(example.contains("\"format_version\": 1"));
    let loaded = Project::from_json(example).unwrap();
    assert!(loaded.to_json().unwrap().contains("\"format_version\": 2"));
}

fn directory_entries(dir: &std::path::Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn saving_replaces_the_file_atomically_and_failures_keep_the_original() {
    let dir = scratch("atomic");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("p.flowsheet.json");
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.save(&path).unwrap();
    let original = std::fs::read_to_string(&path).unwrap();

    // Overwriting an existing file works and leaves no temporary files.
    project.display_units = flowsheet_project::DisplayUnitSet::Si;
    project.save(&path).unwrap();
    let replaced = std::fs::read_to_string(&path).unwrap();
    assert_ne!(replaced, original);
    assert_eq!(
        directory_entries(&dir),
        vec!["p.flowsheet.json".to_string()]
    );

    // A save that fails before writing leaves the file untouched.
    let mut broken = project.clone();
    let feed = unit(&broken, UnitKind::Feed);
    broken.feeds.get_mut(&feed).unwrap()[0] = f64::NAN;
    assert!(matches!(broken.save(&path), Err(ProjectError::Invalid(_))));
    assert_eq!(std::fs::read_to_string(&path).unwrap(), replaced);

    // A save whose final rename fails cleans up its temporary file.
    let blocked = dir.join("blocked");
    std::fs::create_dir_all(blocked.join("occupied")).unwrap();
    assert!(matches!(project.save(&blocked), Err(ProjectError::Io(_))));
    assert!(blocked.join("occupied").is_dir());
    assert_eq!(
        directory_entries(&dir),
        vec!["blocked".to_string(), "p.flowsheet.json".to_string()]
    );

    // Windows refuses to replace a read-only file: the original survives.
    #[cfg(windows)]
    {
        let mut perms = std::fs::metadata(&path).unwrap().permissions();
        perms.set_readonly(true);
        std::fs::set_permissions(&path, perms.clone()).unwrap();
        project.display_units = flowsheet_project::DisplayUnitSet::ProcessMetric;
        let result = project.save(&path);
        #[allow(clippy::permissions_set_readonly_false)]
        perms.set_readonly(false);
        std::fs::set_permissions(&path, perms).unwrap();
        assert!(matches!(result, Err(ProjectError::Io(_))), "{result:?}");
        assert_eq!(std::fs::read_to_string(&path).unwrap(), replaced);
        assert_eq!(
            directory_entries(&dir),
            vec!["blocked".to_string(), "p.flowsheet.json".to_string()]
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn non_finite_numbers_are_rejected_when_saving() {
    let base = demonstration(CubicVariantChoice::PengRobinson);
    let invalid = |p: &Project| matches!(p.to_json(), Err(ProjectError::Invalid(_)));
    assert!(base.to_json().is_ok());

    let mut p = base.clone();
    let cooler = unit(&p, UnitKind::Cooler);
    let UnitSpec::HeatExchange { thermal, .. } = p.flowsheet.unit(cooler).unwrap().spec else {
        panic!("cooler spec");
    };
    p.flowsheet.set_spec(
        cooler,
        UnitSpec::HeatExchange {
            thermal,
            pressure_drop_pa: f64::INFINITY,
        },
    );
    assert!(invalid(&p));

    let mut p = base.clone();
    let valve = unit(&p, UnitKind::Valve);
    p.flowsheet.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(f64::NAN)),
        },
    );
    assert!(invalid(&p));

    let mut p = base.clone();
    let feed = unit(&p, UnitKind::Feed);
    p.feeds.get_mut(&feed).unwrap()[1] = f64::NEG_INFINITY;
    assert!(invalid(&p));

    let mut p = base.clone();
    p.layout.insert(feed, [f32::NAN, 0.0]);
    assert!(invalid(&p));

    let mut p = base.clone();
    p.package.kij[0][1] = f64::NAN;
    assert!(invalid(&p));

    let mut p = base.clone();
    let mut r = solve(&p);
    r.units.values_mut().next().unwrap().duty_w = f64::NAN;
    p.results = Some(r);
    assert!(invalid(&p));

    let mut p = base.clone();
    let mut r = solve(&p);
    let s = r.streams.values_mut().next().unwrap();
    s.state.as_mut().unwrap().mole_fractions[0] = f64::INFINITY;
    p.results = Some(r);
    assert!(invalid(&p));
}

#[test]
fn solving_with_a_package_of_another_selection_fails() {
    let project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    let mut other = project.package.clone();
    other.variant = CubicVariantChoice::PengRobinson;
    let stale = other.prepare().unwrap();
    let f = solve_project(&project, &stale, |_| ControlFlow::Continue(()), 0).unwrap_err();
    assert!(!f.cancelled);
    assert!(
        f.message.contains("different package selection"),
        "{}",
        f.message
    );

    let mut options = project.package.clone();
    options.ph_maximum_temperature_k = 550.0;
    let f = solve_project(
        &project,
        &options.prepare().unwrap(),
        |_| ControlFlow::Continue(()),
        0,
    )
    .unwrap_err();
    assert!(f.message.contains("different package selection"));

    // The prepared package carries the same identity as the selection.
    let own = project.package.prepare().unwrap();
    assert_eq!(
        Some(own.fingerprint()),
        Fingerprints::package(&project.package)
    );
    assert_eq!(own.provenance(), &project.package.provenance().unwrap());
    let r = solve_project(&project, &own, |_| ControlFlow::Continue(()), 0).unwrap();
    let fps = project.fingerprints();
    assert!(
        project
            .flowsheet
            .units()
            .all(|u| r.unit_currency(u.id, &fps) == Currency::Current)
    );
}

#[test]
fn stored_results_must_match_the_component_count() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.results = Some(solve(&project));
    let text = project.to_json().unwrap();
    let mut value: serde_json::Value = serde_json::from_str(&text).unwrap();
    let streams = value["results"]["streams"].as_object_mut().unwrap();
    let state = streams.values_mut().next().unwrap()["state"]
        .as_object_mut()
        .unwrap();
    state["mole_fractions"].as_array_mut().unwrap().pop();
    let err = Project::from_json(&value.to_string()).unwrap_err();
    assert!(
        matches!(&err, ProjectError::Invalid(m) if m.contains("mole fractions")),
        "{err}"
    );

    let mut value: serde_json::Value = serde_json::from_str(&text).unwrap();
    let streams = value["results"]["streams"].as_object_mut().unwrap();
    let phase = &mut streams.values_mut().next().unwrap()["state"]["phases"][0];
    phase["mole_fractions"]
        .as_array_mut()
        .unwrap()
        .push(0.0.into());
    assert!(matches!(
        Project::from_json(&value.to_string()),
        Err(ProjectError::Invalid(_))
    ));
}

#[test]
fn packages_without_components_are_rejected_on_load() {
    let project = demonstration(CubicVariantChoice::PengRobinson);
    let mut value: serde_json::Value = serde_json::from_str(&project.to_json().unwrap()).unwrap();
    value["package"]["components"] = serde_json::json!([]);
    value["package"]["kij"] = serde_json::json!([]);
    value["feeds"] = serde_json::json!([]);
    let err = Project::from_json(&value.to_string()).unwrap_err();
    assert!(
        matches!(&err, ProjectError::Package(m) if m.contains("at least one component")),
        "{err}"
    );
}
