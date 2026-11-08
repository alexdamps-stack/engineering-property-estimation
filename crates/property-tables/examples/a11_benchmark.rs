//! A11 reproducible profile, accuracy and speed benchmark.
//!
//! `cargo run --release -p property-tables --example a11_benchmark`
//! (`-- --quick` uses fewer points and repetitions for a smoke run).
//! Timings are wall-clock on the current machine, single thread; they are
//! evidence for docs/a11-status.md, never a pull-request gate.
use property_tables::experiment::*;
use property_tables::{OutsidePolicy, QualificationCriteria, TableDomain, TabulatedMixture};
use std::hint::black_box;
use std::time::Instant;
use thermo_models::RootSelection;
use thermo_solvers::{
    FlashOptions, PhFlashOptions, StabilityOptions, assess_stability, ph_flash, tp_flash,
};
use thermo_types::{Pressure, Temperature};

/// Run `f` over all points `reps` times; report (min, median, max) ns per call.
fn time(points: &[(f64, f64)], reps: usize, mut f: impl FnMut(f64, f64)) -> (f64, f64, f64) {
    let mut per_call: Vec<f64> = (0..reps)
        .map(|_| {
            let start = Instant::now();
            for &(t, p) in points {
                f(black_box(t), black_box(p));
            }
            start.elapsed().as_nanos() as f64 / points.len() as f64
        })
        .collect();
    per_call.sort_by(f64::total_cmp);
    (per_call[0], per_call[reps / 2], per_call[reps - 1])
}

fn tp(t: f64, p: f64) -> (Temperature, Pressure) {
    (Temperature::new(t).unwrap(), Pressure::new(p).unwrap())
}

fn row(name: &str, (lo, med, hi): (f64, f64, f64)) -> f64 {
    println!("| {name} | {lo:.0} | {med:.0} | {hi:.0} |");
    med
}

fn main() {
    let quick = std::env::args().any(|a| a == "--quick");
    let (n_fast, n_slow, reps, n_err) = if quick {
        (200, 20, 3, 2_000)
    } else {
        (5_000, 500, 11, 200_000)
    };
    let model = model();
    let points = sample_points(2026, n_fast);
    let slow = &points[..n_slow];

    println!("# A11 benchmark (quick={quick})");
    println!("mixture {COMPONENTS:?} x={FRACTIONS:?}, exact PR, zero kij");
    println!("domain {DOMAIN:?}");
    println!(
        "timing: {n_fast} points (flash/stability/PH: {n_slow}), {reps} repetitions, ns/call min|median|max\n"
    );

    println!("## Direct-calculation profile (#60)\n");
    println!("| operation | min ns | median ns | max ns |\n|---|---|---|---|");
    row(
        "CubicModel::compressibility_roots",
        time(&points, reps, |t, p| {
            let (t, p) = tp(t, p);
            black_box(model.compressibility_roots(t, p).unwrap());
        }),
    );
    row(
        "CubicModel::ln_fugacity",
        time(&points, reps, |t, p| {
            let (t, p) = tp(t, p);
            black_box(
                model
                    .ln_fugacity(t, p, RootSelection::MinimumGibbs)
                    .unwrap(),
            );
        }),
    );
    let direct = row(
        "CubicModel::evaluate (all properties, no phase check)",
        time(&points, reps, |t, p| {
            let (t, p) = tp(t, p);
            black_box(model.evaluate(t, p, RootSelection::MinimumGibbs).unwrap());
        }),
    );
    row(
        "assess_stability (default TPD trials)",
        time(slow, reps, |t, p| {
            let (t, p) = tp(t, p);
            black_box(assess_stability(&model, t, p, StabilityOptions::default()).unwrap());
        }),
    );
    let flash = row(
        "tp_flash (stability + single-phase properties)",
        time(slow, reps, |t, p| {
            let (t, p) = tp(t, p);
            black_box(tp_flash(&model, t, p, FlashOptions::default()).unwrap());
        }),
    );
    let targets: Vec<(f64, f64)> = slow
        .iter()
        .map(|&(t, p)| {
            let (tt, pp) = tp(t, p);
            (
                model
                    .evaluate(tt, pp, RootSelection::MinimumGibbs)
                    .unwrap()
                    .enthalpy,
                p,
            )
        })
        .collect();
    let ph = row(
        "ph_flash (single-phase target h)",
        time(&targets, reps, |h, p| {
            black_box(
                ph_flash(
                    &model,
                    Pressure::new(p).unwrap(),
                    h,
                    PhFlashOptions::default(),
                )
                .unwrap(),
            );
        }),
    );

    println!("\n## Table build\n");
    let start = Instant::now();
    let mixture = mixture(OutsidePolicy::Error).unwrap();
    let build = start.elapsed();
    let q = mixture.table().qualification();
    println!(
        "nodes {}x{} = {}, node data {} bytes, build incl. qualification {:.2} s, qualification points {}, min reduced bulk modulus {:.4}, min nontrivial TPD {:?}",
        DOMAIN.temperature_nodes,
        DOMAIN.pressure_nodes,
        DOMAIN.temperature_nodes * DOMAIN.pressure_nodes,
        mixture.table().memory_bytes(),
        build.as_secs_f64(),
        q.points_checked,
        q.minimum_reduced_bulk_modulus,
        q.minimum_nontrivial_tpd
    );

    println!("\n## Table evaluation speed (#61)\n");
    println!("| operation | min ns | median ns | max ns |\n|---|---|---|---|");
    let table = row(
        "GibbsTable::evaluate (same property set)",
        time(&points, reps, |t, p| {
            black_box(mixture.table().evaluate(t, p).unwrap());
        }),
    );
    row(
        "TabulatedMixture::evaluate (domain check + table)",
        time(&points, reps, |t, p| {
            black_box(mixture.evaluate(t, p).unwrap());
        }),
    );
    let table_ph = row(
        "TabulatedMixture::ph (table PH inversion)",
        time(&targets, reps, |h, p| {
            black_box(mixture.ph(p, h).unwrap());
        }),
    );
    println!(
        "\nmedian speedup vs CubicModel::evaluate: {:.1}x; vs tp_flash: {:.0}x; table PH vs ph_flash: {:.0}x",
        direct / table,
        flash / table,
        ph / table_ph
    );

    // PH accuracy: table temperature against the direct PH flash temperature for the same target.
    let (mut dt_max, mut dh_max, mut iterations) = (0.0_f64, 0.0_f64, 0usize);
    for (&(t, _), &(h, p)) in slow.iter().zip(&targets) {
        let tab = mixture.ph(p, h).unwrap().properties;
        let dir = ph_flash(
            &model,
            Pressure::new(p).unwrap(),
            h,
            PhFlashOptions::default(),
        )
        .unwrap();
        iterations += dir.iterations;
        dt_max = dt_max
            .max((tab.temperature_k - dir.temperature_k).abs())
            .max((tab.temperature_k - t).abs());
        let (tt, pp) = tp(tab.temperature_k, p);
        dh_max = dh_max.max(
            (model
                .evaluate(tt, pp, RootSelection::MinimumGibbs)
                .unwrap()
                .enthalpy
                - h)
                .abs(),
        );
    }
    println!(
        "PH over {} targets: max |T_table - T_direct| = {dt_max:.2e} K; max |h_direct(T_table) - h_target| = {dh_max:.2e} J/mol; mean ph_flash outer iterations {:.1}",
        slow.len(),
        iterations as f64 / slow.len() as f64
    );

    println!("\n## Accuracy against the direct model (#61)\n");
    println!("values: g,h abs J/mol; s abs J/(mol K); v, rho relative. derivatives: relative.\n");
    println!(
        "| grid | points | g | h | s | v | rho | cp | cv | dv/dT_p | dv/dp_T | w |\n|---|---|---|---|---|---|---|---|---|---|---|---|"
    );
    let print = |label: &str, e: ErrorSummary| {
        println!(
            "| {label} | {} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} |",
            e.points,
            e.gibbs_abs,
            e.enthalpy_abs,
            e.entropy_abs,
            e.volume_rel,
            e.density_rel,
            e.cp_rel,
            e.cv_rel,
            e.dv_dt_rel,
            e.dv_dp_rel,
            e.speed_of_sound_rel
        )
    };
    let approved = format!(
        "{}x{} random",
        DOMAIN.temperature_nodes, DOMAIN.pressure_nodes
    );
    print(
        &approved,
        compare(mixture.table(), &model, &sample_points(7, n_err)).unwrap(),
    );
    print(
        "approved, cell centres",
        compare(mixture.table(), &model, &cell_centres(&DOMAIN)).unwrap(),
    );
    if !quick {
        for (nt, np) in [(29, 31), (29, 61), (57, 61), (113, 181)] {
            let d = TableDomain {
                temperature_nodes: nt,
                pressure_nodes: np,
                ..DOMAIN
            };
            let c = QualificationCriteria {
                samples_per_cell: 1,
                ..Default::default()
            };
            let m = TabulatedMixture::new(model.clone(), d, c, OutsidePolicy::Error).unwrap();
            print(
                &format!("{nt}x{np} random (study)"),
                compare(m.table(), &model, &sample_points(7, 20_000)).unwrap(),
            );
        }
    }

    println!("\n## Outside the approved domain (#62)\n");
    let fallback = mixture.with_policy(OutsidePolicy::DirectModel);
    for (t, p) in [
        (350.0, 4.0e6),
        (450.0, 2.0e6),
        (300.0, 1.2e7),
        (200.0, 3.0e6),
    ] {
        let strict = mixture.evaluate(t, p).map(|e| e.source);
        let direct = fallback.evaluate(t, p).map(|e| e.source);
        println!(
            "T={t} K p={p} Pa: Error policy -> {}; DirectModel policy -> {}",
            strict.map_or_else(|e| e.to_string(), |s| format!("{s:?}")),
            direct.map_or_else(|e| e.to_string(), |s| format!("{s:?}"))
        );
    }
}
