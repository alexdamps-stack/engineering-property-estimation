//! Equilibrium-stage foundation of the columns group against the independent
//! exact-rational fixture `validation/columns/stages-v1.json`, and a stage
//! built from a cubic TP flash.

use flowsheet_core::MaterialStream;
use flowsheet_core::units::CalcContext;
use flowsheet_core::units::columns::stages::{
    ColumnEnds, Stage, StageCascade, StageFeed, StageProperties,
};
use property_packages::CubicPackage;
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

fn fixture() -> Value {
    serde_json::from_str(include_str!("../../../validation/columns/stages-v1.json")).unwrap()
}

/// Parses a decimal or `a/b` string.
fn num(v: &Value) -> f64 {
    let s = v.as_str().unwrap();
    match s.split_once('/') {
        Some((a, b)) => a.parse::<f64>().unwrap() / b.parse::<f64>().unwrap(),
        None => s.parse().unwrap(),
    }
}

fn row(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(num).collect()
}

fn rows(v: &Value) -> Vec<Vec<f64>> {
    v.as_array().unwrap().iter().map(row).collect()
}

struct Case {
    cascade: StageCascade,
    props: Vec<StageProperties>,
    k: Vec<Vec<f64>>,
    x: Vec<Vec<f64>>,
    liquid_summation: Vec<f64>,
    vapor_summation: Vec<f64>,
}

fn case(name: &str) -> Case {
    let doc = fixture();
    let c = &doc["cases"][name];
    let (vapor, liquid) = (row(&c["vapor"]), row(&c["liquid"]));
    let (u, w, q) = (
        row(&c["liquid_draw"]),
        row(&c["vapor_draw"]),
        row(&c["duty_w"]),
    );
    let (k, x, y) = (rows(&c["k"]), rows(&c["x"]), rows(&c["y"]));
    let (hl, hv) = (row(&c["h_liquid"]), row(&c["h_vapor"]));
    let stages = (0..vapor.len())
        .map(|j| Stage {
            liquid_draw_mol_s: u[j],
            vapor_draw_mol_s: w[j],
            duty_w: q[j],
            ..Stage::equilibrium(350.0, 1e5, liquid[j], vapor[j], x[j].clone(), y[j].clone())
        })
        .collect();
    let feeds = c["feeds"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| StageFeed {
            stage: f["stage"].as_u64().unwrap() as usize,
            molar_flow_mol_s: num(&f["flow"]),
            z: row(&f["z"]),
            molar_enthalpy_j_mol: num(&f["h"]),
        })
        .collect();
    let props = (0..vapor.len())
        .map(|j| StageProperties {
            k: k[j].clone(),
            liquid_enthalpy_j_mol: hl[j],
            vapor_enthalpy_j_mol: hv[j],
        })
        .collect();
    Case {
        cascade: StageCascade { stages, feeds },
        props,
        k,
        x,
        liquid_summation: row(&c["liquid_summation"]),
        vapor_summation: row(&c["vapor_summation"]),
    }
}

fn close(a: f64, b: f64, rel: f64, abs: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs()) + abs
}

#[test]
fn tridiagonal_component_balance_matches_the_exact_reference() {
    for name in ["distillation", "absorber"] {
        let c = case(name);
        let solved = c.cascade.solve_liquid_fractions(&c.k).unwrap();
        for (j, (s, r)) in solved.iter().zip(&c.x).enumerate() {
            for (i, (a, b)) in s.iter().zip(r).enumerate() {
                assert!(
                    close(*a, *b, 1e-12, 1e-15),
                    "{name} stage {} x{i}: {a} vs {b}",
                    j + 1
                );
            }
        }
    }
}

#[test]
fn reference_profiles_close_every_mesh_residual() {
    for name in ["distillation", "absorber"] {
        let c = case(name);
        let flow_scale: f64 = c.cascade.feeds.iter().map(|f| f.molar_flow_mol_s).sum();
        for r in c.cascade.total_flow_residuals() {
            assert!(r.abs() <= 1e-12 * flow_scale, "{name}: total flow {r}");
        }
        let res = c.cascade.residuals(&c.props).unwrap();
        let energy_scale = flow_scale * 5e4;
        for (j, r) in res.iter().enumerate() {
            for m in &r.component_mol_s {
                assert!(
                    m.abs() <= 1e-12 * flow_scale,
                    "{name} stage {}: M {m}",
                    j + 1
                );
            }
            for e in &r.equilibrium {
                assert!(e.abs() <= 1e-14, "{name} stage {}: E {e}", j + 1);
            }
            assert!(r.energy_w.abs() <= 1e-12 * energy_scale, "{name}: H {r:?}");
            // The fixed-K profile is not normalized: summations match the reference.
            assert!(close(
                r.liquid_summation,
                c.liquid_summation[j],
                1e-11,
                1e-14
            ));
            assert!(close(r.vapor_summation, c.vapor_summation[j], 1e-11, 1e-14));
        }
        let b = c.cascade.balance(&c.props).unwrap();
        for v in &b.component_residuals_mol_s {
            assert!(v.abs() <= 1e-12 * flow_scale, "{name}: {b:?}");
        }
        assert!(
            b.energy_residual_w.abs() <= 1e-12 * energy_scale,
            "{name}: {b:?}"
        );
    }
}

#[test]
fn stage_residuals_sum_to_the_overall_balance_and_end_duties() {
    let mut c = case("distillation");
    // Perturb compositions and the reboiler duty: stage residuals no longer
    // vanish, but their sums still equal the overall cascade residuals.
    c.cascade.stages[2].x[0] += 0.01;
    c.cascade.stages[1].y[2] -= 0.003;
    c.cascade.stages[4].duty_w += 250.0;
    let res = c.cascade.residuals(&c.props).unwrap();
    let b = c.cascade.balance(&c.props).unwrap();
    for i in 0..3 {
        let sum: f64 = res.iter().map(|r| r.component_mol_s[i]).sum();
        assert!(close(sum, b.component_residuals_mol_s[i], 1e-10, 1e-10));
    }
    let energy: f64 = res.iter().map(|r| r.energy_w).sum();
    assert!(close(energy, b.energy_residual_w, 1e-10, 1e-6));
    assert!((res[4].energy_w - 250.0).abs() < 1e-6);

    // Total condenser and reboiler: only the end stages carry duty.
    let duties: Vec<f64> = c.cascade.stages.iter().map(|s| s.duty_w).collect();
    let ends = ColumnEnds::DISTILLATION;
    assert!(
        ends.check(&c.cascade)
            .unwrap_err()
            .contains("no energy port")
    );
    for s in &mut c.cascade.stages[1..4] {
        s.duty_w = 0.0;
    }
    let ports = ends.port_duties(&c.cascade).unwrap();
    assert_eq!(ports, [Some(-duties[0]), Some(duties[4])]);
    assert!(duties[0] < 0.0 && duties[4] > 0.0);
}

fn package() -> CubicPackage {
    let composition = Composition::new(
        CompositionBasis::Mole,
        vec![
            (ComponentId::new("ethane").unwrap(), 0.5),
            (ComponentId::new("propane").unwrap(), 0.5),
        ],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition,
        CubicVariant::PengRobinson,
        vec![vec![0.0; 2]; 2],
        "synthetic ethane/propane pair with zero binary interaction for stage tests",
    )
    .unwrap();
    CubicPackage::new(model, Default::default())
}

/// A single adiabatic stage fed by a two-phase TP flash state is that flash:
/// cubic stage properties at the flash phases close every MESH residual.
#[test]
fn cubic_stage_properties_reproduce_a_tp_flash() {
    let pkg = package();
    let ctx = CalcContext::new(&pkg, None);
    let (t, p, n) = (260.0, 0.8e6, 10.0);
    let z = [0.5, 0.5];
    let state = ctx.tp(&z, t, p).unwrap();
    assert_eq!(state.phases.len(), 2, "the test point must be two-phase");
    let vapor = state
        .phases
        .iter()
        .max_by(|a, b| a.state.molar_volume.total_cmp(&b.state.molar_volume))
        .unwrap();
    let liquid = state
        .phases
        .iter()
        .min_by(|a, b| a.state.molar_volume.total_cmp(&b.state.molar_volume))
        .unwrap();
    let feed = StageFeed::from_stream(
        1,
        &MaterialStream {
            molar_flow_mol_s: n,
            state: Some(state.clone()),
        },
        2,
    )
    .unwrap();
    let cascade = StageCascade {
        stages: vec![Stage::equilibrium(
            t,
            p,
            n * liquid.phase_fraction,
            n * vapor.phase_fraction,
            liquid.composition.clone(),
            vapor.composition.clone(),
        )],
        feeds: vec![feed],
    };
    let props = ctx.stage_properties(&cascade.stages[0]).unwrap();
    assert!(props.k[0] > 1.0 && props.k[1] < 1.0, "{props:?}");
    let r = &cascade.residuals(&[props]).unwrap()[0];
    for m in &r.component_mol_s {
        assert!(m.abs() < 1e-9 * n, "{r:?}");
    }
    for e in &r.equilibrium {
        assert!(e.abs() < 1e-7, "{r:?}");
    }
    let scale = n * state.enthalpy().abs().max(1e3);
    assert!(r.energy_w.abs() < 1e-9 * scale, "{r:?}");
}
