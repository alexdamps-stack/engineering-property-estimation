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
    let package = project.package.build().unwrap();
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
    let package = project.package.build().unwrap();
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
fn reads_v1_writes_v2_and_rejects_unknown_versions_before_content() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.results = Some(solve(&project));
    let mut value: serde_json::Value = serde_json::from_str(&project.to_json().unwrap()).unwrap();
    assert_eq!(value["format_version"], 2);
    value["format_version"] = 1.into();
    // Pre-recycle v1 files omitted this optional field.
    value["results"].as_object_mut().unwrap().remove("recycle");
    let loaded = Project::from_json(&value.to_string()).unwrap();
    assert!(loaded.to_json().unwrap().contains("\"format_version\": 2"));
    assert!(loaded.results.unwrap().recycle.is_none());
    assert_eq!(
        Project::from_json(r#"{"format_version":99,"units":"invalid"}"#).unwrap_err(),
        ProjectError::UnsupportedVersion(Some(99))
    );
}

#[test]
fn rejects_inconsistent_overall_and_phase_result_dimensions() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    project.results = Some(solve(&project));
    for phase in [false, true] {
        for too_long in [false, true] {
            let mut edited = project.clone();
            let state = edited
                .results
                .as_mut()
                .unwrap()
                .streams
                .values_mut()
                .find_map(|s| s.state.as_mut())
                .unwrap();
            let fractions = if phase {
                &mut state.phases[0].mole_fractions
            } else {
                &mut state.mole_fractions
            };
            if too_long {
                fractions.push(0.0);
            } else {
                fractions.pop();
            }
            assert!(matches!(
                Project::from_json(&edited.to_json().unwrap()),
                Err(ProjectError::Invalid(_))
            ));
        }
    }
}

#[test]
fn atomic_save_replaces_file_and_cleans_up_on_failed_replace() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("project.json");
    let project = demonstration(CubicVariantChoice::PengRobinson);
    project.save(&path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        project.to_json().unwrap()
    );
    std::fs::write(&path, "old content").unwrap();
    project.save(&path).unwrap();
    assert_eq!(
        std::fs::read_to_string(&path).unwrap(),
        project.to_json().unwrap()
    );
    // A directory cannot be replaced by a project file. Its existing content
    // and the successful save remain intact, and no temporary file leaks.
    let blocked = directory.path().join("blocked");
    std::fs::create_dir(&blocked).unwrap();
    std::fs::write(blocked.join("original"), "preserve me").unwrap();
    assert!(matches!(project.save(&blocked), Err(ProjectError::Io(_))));
    assert_eq!(
        std::fs::read_to_string(blocked.join("original")).unwrap(),
        "preserve me"
    );
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 2);
}

#[cfg(windows)]
#[test]
fn failed_atomic_save_preserves_locked_original() {
    use std::os::windows::fs::OpenOptionsExt;
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("locked.json");
    std::fs::write(&path, "original content").unwrap();
    let _lock = std::fs::OpenOptions::new()
        .read(true)
        .share_mode(1)
        .open(&path)
        .unwrap();
    let project = demonstration(CubicVariantChoice::PengRobinson);
    assert!(matches!(project.save(&path), Err(ProjectError::Io(_))));
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "original content");
    assert_eq!(std::fs::read_dir(directory.path()).unwrap().count(), 1);
}

#[test]
fn atomic_save_accepts_a_bare_relative_filename() {
    // Do not change cwd: other tests may be running concurrently.
    let cleanup = tempfile::NamedTempFile::new_in(".")
        .unwrap()
        .into_temp_path();
    let relative = PathBuf::from(cleanup.file_name().unwrap());
    let project = demonstration(CubicVariantChoice::PengRobinson);
    project.save(&relative).unwrap();
    assert_eq!(
        Project::load(&relative).unwrap().to_json().unwrap(),
        project.to_json().unwrap()
    );
    project.save(&relative).unwrap();
    drop(cleanup);
}

#[test]
fn previous_engine_results_load_as_stale() {
    let mut project = demonstration(CubicVariantChoice::PengRobinson);
    let mut results = solve(&project);
    results.provenance.engine = "flowsheet-core a10-tears-v1".into();
    // Recreate the historical package checksum with the previous engine
    // identity, while keeping every physical input unchanged.
    let text = serde_json::to_string(&(&project.package, &results.provenance)).unwrap();
    let mut checksum = 0xcbf2_9ce4_8422_2325_u64;
    for field in ["package", &text] {
        for byte in (field.len() as u64)
            .to_le_bytes()
            .iter()
            .chain(field.as_bytes())
        {
            checksum ^= u64::from(*byte);
            checksum = checksum.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    let old = Fingerprints::compute_with(
        &project.flowsheet,
        &project.feeds,
        Some(Fingerprint(checksum)),
    );
    for (id, result) in &mut results.units {
        result.fingerprint = old.units[id];
    }
    for (id, result) in &mut results.streams {
        result.fingerprint = old.stream(&project.flowsheet, *id).unwrap();
    }
    project.results = Some(results);
    let project = Project::from_json(&project.to_json().unwrap()).unwrap();
    let results = project.results.as_ref().unwrap();
    assert_eq!(results.provenance.engine, "flowsheet-core a10-tears-v1");
    assert_ne!(
        results.provenance.engine,
        flowsheet_project::ENGINE_REVISION
    );
    let current = project.fingerprints();
    assert!(
        project
            .flowsheet
            .units()
            .all(|u| results.unit_currency(u.id, &current) == Currency::Stale)
    );
    assert!(project.flowsheet.streams().all(|s| results.stream_currency(
        &project.flowsheet,
        s.id,
        &current
    ) == Currency::Stale));
}
