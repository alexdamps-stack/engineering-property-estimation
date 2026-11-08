//! Headless run of the A9 acceptance workflow (`docs/flowsheet-and-ui.md`,
//! section 4.1) through the same project API the desktop uses. The desktop
//! executable exposes it as `--self-test`, so an installed Windows build can
//! verify calculation, invalidation and persistence without a display.

use std::ops::ControlFlow;
use std::path::Path;

use flowsheet_core::{PressureSpec, UnitId, UnitKind, UnitSpec};
use serde::Serialize;

use crate::time::now_unix;
use crate::{CubicVariantChoice, Currency, ProjectError, ResultSet, demonstration, solve_project};

/// Outcome of [`run_acceptance`].
#[derive(Debug, Clone, Serialize)]
pub struct AcceptanceReport {
    /// Program and engine identity.
    pub engine: String,
    /// Completed checks, in order.
    pub steps: Vec<String>,
    /// Separator temperature, K, at the original and changed valve pressure.
    pub separator_temperature_k: [f64; 2],
    /// Separator vapor fraction at the original and changed valve pressure.
    pub separator_vapor_fraction: [f64; 2],
    /// Largest relative difference between results before saving and after
    /// reloading and recalculating.
    pub max_reload_relative_difference: f64,
}

fn check(ok: bool, what: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(format!("acceptance check failed: {what}"))
    }
}

fn separator_state(r: &ResultSet, project: &crate::Project) -> Result<(f64, f64), String> {
    let fs = &project.flowsheet;
    let stream = fs
        .streams()
        .find(|s| s.name == "expand")
        .ok_or("expand stream missing")?;
    let state = r
        .streams
        .get(&stream.id)
        .and_then(|s| s.state.as_ref())
        .ok_or("separator inlet has no state")?;
    Ok((
        state.temperature_k,
        state
            .vapor_fraction
            .ok_or("separator inlet is supercritical")?,
    ))
}

fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// Largest relative difference of stream temperatures, pressures, flows,
/// enthalpies, compositions and unit duties between two result sets.
fn max_difference(a: &ResultSet, b: &ResultSet) -> Result<f64, String> {
    check(a.streams.len() == b.streams.len(), "same stream results")?;
    let mut worst = 0.0_f64;
    for (id, x) in &a.streams {
        let y = b.streams.get(id).ok_or("stream missing after reload")?;
        worst = worst.max(relative(x.molar_flow_mol_s, y.molar_flow_mol_s));
        match (&x.state, &y.state) {
            (Some(s), Some(t)) => {
                worst = worst
                    .max(relative(s.temperature_k, t.temperature_k))
                    .max(relative(s.pressure_pa, t.pressure_pa))
                    .max(relative(s.enthalpy_j_mol, t.enthalpy_j_mol));
                for (u, v) in s.mole_fractions.iter().zip(&t.mole_fractions) {
                    worst = worst.max((u - v).abs());
                }
            }
            (None, None) => {}
            _ => return Err("phase presence changed after reload".into()),
        }
    }
    for (id, x) in &a.units {
        let y = b.units.get(id).ok_or("unit missing after reload")?;
        if x.duty_w != 0.0 || y.duty_w != 0.0 {
            worst = worst.max(relative(x.duty_w, y.duty_w));
        }
    }
    Ok(worst)
}

/// Runs the acceptance workflow, writing temporary project files to `dir`.
pub fn run_acceptance(dir: &Path) -> Result<AcceptanceReport, String> {
    let mut steps = Vec::new();
    let e = |e: ProjectError| e.to_string();

    // 1. Create the demonstration and solve it.
    let mut project = demonstration(CubicVariantChoice::PengRobinsonBostonMathias);
    check(
        project.flowsheet.diagnostics().is_empty(),
        "complete topology",
    )?;
    check(
        project.input_issues().is_empty(),
        "complete feed compositions",
    )?;
    let package = project.package.build().map_err(e)?;
    let first = solve_project(
        &project,
        &package,
        |_| ControlFlow::Continue(()),
        now_unix(),
    )
    .map_err(|f| f.message)?;
    project.results = Some(first.clone());
    let fps = project.fingerprints();
    check(
        project
            .flowsheet
            .units()
            .all(|u| first.unit_currency(u.id, &fps) == Currency::Current)
            && project.flowsheet.streams().all(|s| {
                first.stream_currency(&project.flowsheet, s.id, &fps) == Currency::Current
            }),
        "all results current after solving",
    )?;
    for u in first.units.values() {
        if let (Some(c), Some(en)) = (u.max_scaled_component_residual, u.scaled_energy_residual) {
            check(
                c <= 1e-8 && en.abs() <= 2e-9,
                "unit balances within A8 limits",
            )?;
        }
    }
    let (t1, v1) = separator_state(&first, &project)?;
    check(v1 > 0.0 && v1 < 1.0, "two phases enter the separator")?;
    steps.push(format!(
        "solved demonstration: {} units, {} streams, separator inlet {t1:.3} K, vapor fraction {v1:.4}",
        first.units.len(),
        first.streams.len()
    ));

    // 2. Change the valve pressure: only the valve and downstream go stale.
    let valve = unit(&project, UnitKind::Valve)?;
    project.flowsheet.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(2.0e6)),
        },
    );
    let fps = project.fingerprints();
    let fs = &project.flowsheet;
    for u in fs.units() {
        let downstream = matches!(
            u.kind,
            UnitKind::Valve | UnitKind::Separator | UnitKind::Product
        );
        let expected = if downstream {
            Currency::Stale
        } else {
            Currency::Current
        };
        check(
            first.unit_currency(u.id, &fps) == expected,
            &format!("{} currency after valve edit", u.name),
        )?;
    }
    for s in fs.streams() {
        let downstream = s.from.unit == valve
            || fs.unit(s.from.unit).map(|u| u.kind) == Some(UnitKind::Separator);
        let expected = if downstream {
            Currency::Stale
        } else {
            Currency::Current
        };
        check(
            first.stream_currency(fs, s.id, &fps) == expected,
            &format!("stream {} currency after valve edit", s.name),
        )?;
    }
    steps.push("valve pressure 15 -> 20 bar(a): valve, separator, products and their streams stale; upstream current".into());

    // 3. Recalculate.
    let second = solve_project(
        &project,
        &package,
        |_| ControlFlow::Continue(()),
        now_unix(),
    )
    .map_err(|f| f.message)?;
    project.results = Some(second.clone());
    let (t2, v2) = separator_state(&second, &project)?;
    check(t2 != t1, "changed separator temperature")?;
    steps.push(format!(
        "recalculated: separator inlet {t2:.3} K, vapor fraction {v2:.4}"
    ));

    // 4. Save, reload and compare.
    let path = dir.join("acceptance.flowsheet.json");
    project.save(&path).map_err(e)?;
    let saved = std::fs::read_to_string(&path).map_err(|e| e.to_string())?;
    let mut loaded = crate::Project::load(&path).map_err(e)?;
    check(
        loaded.to_json().map_err(e)? == saved,
        "lossless save/reload",
    )?;
    let loaded_results = loaded.results.clone().ok_or("results not reloaded")?;
    let fps = loaded.fingerprints();
    check(
        loaded
            .flowsheet
            .units()
            .all(|u| loaded_results.unit_currency(u.id, &fps) == Currency::Current),
        "reloaded results are current",
    )?;
    let package_loaded = loaded.package.build().map_err(e)?;
    let third = solve_project(
        &loaded,
        &package_loaded,
        |_| ControlFlow::Continue(()),
        now_unix(),
    )
    .map_err(|f| f.message)?;
    let difference = max_difference(&second, &third)?;
    check(
        difference <= 1e-9,
        "recalculated reload matches saved results",
    )?;
    steps.push(format!(
        "saved, reloaded and recalculated: largest relative difference {difference:.3e}"
    ));

    // 5. A failure produces no result and is attributed to its unit.
    let valve = unit(&loaded, UnitKind::Valve)?;
    loaded.flowsheet.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(6.0e6)),
        },
    );
    match solve_project(
        &loaded,
        &package_loaded,
        |_| ControlFlow::Continue(()),
        now_unix(),
    ) {
        Ok(_) => return Err("a valve pressure rise must fail".into()),
        Err(f) => check(
            f.unit == Some(valve) && !f.cancelled,
            "failure names the valve",
        )?,
    }
    let fps = loaded.fingerprints();
    check(
        loaded_results.unit_currency(valve, &fps) == Currency::Stale,
        "old valve result is not current after failed edit",
    )?;
    steps.push(
        "valve pressure above inlet: solve fails at the valve; old valve result stale".into(),
    );

    // 6. Cancellation returns no result.
    let cancelled = solve_project(
        &project,
        &package,
        |p| {
            if p.completed >= 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        now_unix(),
    );
    check(
        matches!(cancelled, Err(ref f) if f.cancelled),
        "cancelled solve returns no result",
    )?;
    steps.push("cancellation after two units returns no result".into());

    // 7. Unknown format versions are rejected.
    let newer = saved.replacen("\"format_version\": 1", "\"format_version\": 2", 1);
    check(
        matches!(
            crate::Project::from_json(&newer),
            Err(ProjectError::UnsupportedVersion(Some(2)))
        ),
        "format version 2 is rejected",
    )?;
    steps.push("format version 2 rejected".into());
    let _ = std::fs::remove_file(&path);

    Ok(AcceptanceReport {
        engine: crate::ENGINE_REVISION.into(),
        steps,
        separator_temperature_k: [t1, t2],
        separator_vapor_fraction: [v1, v2],
        max_reload_relative_difference: difference,
    })
}

fn unit(project: &crate::Project, kind: UnitKind) -> Result<UnitId, String> {
    project
        .flowsheet
        .units()
        .find(|u| u.kind == kind)
        .map(|u| u.id)
        .ok_or_else(|| format!("no {kind} in demonstration"))
}
