//! A11 correctness, domain and fallback tests. Accuracy targets were fixed before
//! the grid was chosen (docs/a11-status.md); they are not relaxed here.
use property_tables::experiment::*;
use property_tables::*;
use std::sync::OnceLock;
use thermo_models::{CubicModel, CubicVariant, RootSelection};
use thermo_solvers::{PhFlashOptions, ph_flash};
use thermo_types::{
    CalculationError, ComponentId, Composition, CompositionBasis, NormalizationPolicy, Pressure,
    Temperature,
};

// Value targets.
const ENTHALPY_ABS: f64 = 0.1; // J/mol
const ENTROPY_ABS: f64 = 1e-3; // J/(mol K)
const DENSITY_REL: f64 = 1e-6;
// Derivative targets (relative).
const DERIVATIVE_REL: f64 = 1e-3;

fn strict() -> &'static TabulatedMixture {
    static M: OnceLock<TabulatedMixture> = OnceLock::new();
    M.get_or_init(|| mixture(OutsidePolicy::Error).expect("approved domain qualifies"))
}

/// The same qualified table with the direct fallback policy.
fn fallback() -> &'static TabulatedMixture {
    static M: OnceLock<TabulatedMixture> = OnceLock::new();
    M.get_or_init(|| strict().with_policy(OutsidePolicy::DirectModel))
}

fn direct(model: &CubicModel, t: f64, p: f64) -> Properties {
    Properties::from_direct(
        &model
            .evaluate(
                Temperature::new(t).unwrap(),
                Pressure::new(p).unwrap(),
                RootSelection::MinimumGibbs,
            )
            .unwrap(),
    )
}

#[test]
fn approved_domain_qualification_is_recorded() {
    let q = strict().table().qualification();
    let d = DOMAIN;
    let expected = ((d.temperature_nodes - 1) * 4 + 1) * ((d.pressure_nodes - 1) * 4 + 1);
    assert_eq!(q.points_checked, expected);
    assert!(q.minimum_reduced_bulk_modulus > q.criteria.minimum_reduced_bulk_modulus);
    assert_eq!(
        strict().table().memory_bytes(),
        d.temperature_nodes * d.pressure_nodes * 32
    );
}

#[test]
fn value_and_derivative_errors_meet_targets_on_fixed_seed_sample() {
    let m = strict();
    for points in [sample_points(61, 1500), cell_centres(&DOMAIN)] {
        let e = compare(m.table(), m.model(), &points).unwrap();
        assert!(e.enthalpy_abs <= ENTHALPY_ABS, "{e:?}");
        assert!(e.entropy_abs <= ENTROPY_ABS, "{e:?}");
        assert!(
            e.density_rel <= DENSITY_REL && e.volume_rel <= DENSITY_REL,
            "{e:?}"
        );
        for d in [
            e.cp_rel,
            e.cv_rel,
            e.dv_dt_rel,
            e.dv_dp_rel,
            e.speed_of_sound_rel,
        ] {
            assert!(d <= DERIVATIVE_REL, "{e:?}");
        }
    }
}

#[test]
fn nodes_reproduce_the_direct_model() {
    let m = strict();
    let d = DOMAIN;
    for (t, p) in [
        (d.minimum_temperature_k, d.minimum_pressure_pa),
        (d.maximum_temperature_k, d.maximum_pressure_pa),
        (
            d.minimum_temperature_k + 2.5 * 20.0,
            d.minimum_pressure_pa
                * (d.maximum_pressure_pa / d.minimum_pressure_pa).powf(40.0 / 90.0),
        ),
    ] {
        let a = m.table().evaluate(t, p).unwrap();
        let b = direct(m.model(), t, p);
        assert!((a.gibbs_energy - b.gibbs_energy).abs() < 1e-8);
        assert!((a.entropy - b.entropy).abs() < 1e-10);
        assert!(((a.molar_volume - b.molar_volume) / b.molar_volume).abs() < 1e-11);
        assert!(((a.dv_dt_p - b.dv_dt_p) / b.dv_dt_p).abs() < 1e-10);
    }
}

#[test]
fn table_properties_are_thermodynamically_consistent() {
    let m = strict();
    for (t, p) in sample_points(5, 50) {
        let s = m.table().evaluate(t, p).unwrap();
        assert!((s.gibbs_energy + t * s.entropy - s.enthalpy).abs() < 1e-9);
        assert!((s.enthalpy - p * s.molar_volume - s.internal_energy).abs() < 1e-9);
        // Central differences on the table itself, within one cell.
        let (dt, dp) = (1e-4, p * 1e-6);
        let (tl, th) = (t - dt, t + dt);
        if !DOMAIN.contains(tl, p)
            || !DOMAIN.contains(th, p)
            || !DOMAIN.contains(t, p + dp)
            || !DOMAIN.contains(t, p - dp)
        {
            continue;
        }
        let cp = (m.table().evaluate(th, p).unwrap().enthalpy
            - m.table().evaluate(tl, p).unwrap().enthalpy)
            / (2.0 * dt);
        assert!(((cp - s.cp) / s.cp).abs() < 1e-5, "cp {cp} vs {}", s.cp);
        // Maxwell: (∂s/∂p)_T = -(∂v/∂T)_p.
        let ds_dp = (m.table().evaluate(t, p + dp).unwrap().entropy
            - m.table().evaluate(t, p - dp).unwrap().entropy)
            / (2.0 * dp);
        assert!(
            ((ds_dp + s.dv_dt_p) / s.dv_dt_p).abs() < 1e-4,
            "{ds_dp} vs {}",
            s.dv_dt_p
        );
        let dv_dp = (m.table().evaluate(t, p + dp).unwrap().molar_volume
            - m.table().evaluate(t, p - dp).unwrap().molar_volume)
            / (2.0 * dp);
        assert!(((dv_dp - s.dv_dp_t) / s.dv_dp_t).abs() < 1e-5);
    }
}

#[test]
fn domain_edges_are_inclusive_and_outside_is_a_clear_error() {
    let m = strict();
    let d = DOMAIN;
    for (t, p) in [
        (d.minimum_temperature_k, d.minimum_pressure_pa),
        (d.maximum_temperature_k, d.minimum_pressure_pa),
        (d.minimum_temperature_k, d.maximum_pressure_pa),
        (d.maximum_temperature_k, d.maximum_pressure_pa),
    ] {
        assert_eq!(m.evaluate(t, p).unwrap().source, Source::Table);
    }
    for (t, p) in [
        (d.minimum_temperature_k - 1e-6, 1e6),
        (d.maximum_temperature_k + 1e-6, 1e6),
        (350.0, d.minimum_pressure_pa * (1.0 - 1e-9)),
        (350.0, d.maximum_pressure_pa * (1.0 + 1e-9)),
    ] {
        let e = m.evaluate(t, p).unwrap_err();
        assert!(
            matches!(e, TableError::OutsideApprovedDomain { .. }),
            "{e:?}"
        );
        assert!(e.to_string().contains("outside the approved table domain"));
        assert!(m.table().evaluate(t, p).is_err());
    }
    for (t, p) in [(f64::NAN, 1e6), (350.0, f64::INFINITY)] {
        for mixture in [strict(), fallback()] {
            assert!(matches!(
                mixture.evaluate(t, p),
                Err(TableError::Calculation(CalculationError::InvalidInput(_)))
            ));
        }
    }
}

#[test]
fn direct_fallback_matches_direct_model_and_reports_its_source() {
    let m = fallback();
    assert_eq!(m.evaluate(350.0, 4e6).unwrap().source, Source::Table);
    for (t, p) in [(450.0, 2e6), (300.0, 1.2e7), (279.0, 5e5)] {
        let e = m.evaluate(t, p).unwrap();
        assert_eq!(e.source, Source::DirectModel);
        let b = direct(m.model(), t, p);
        assert!((e.properties.enthalpy - b.enthalpy).abs() < 1e-9 * b.enthalpy.abs().max(1.0));
        assert!(((e.properties.mass_density - b.mass_density) / b.mass_density).abs() < 1e-12);
    }
    // Inside the two-phase region the direct flash splits: a clear error, no table value.
    assert_eq!(
        m.evaluate(200.0, 3e6).unwrap_err(),
        TableError::DirectResultNotSinglePhase { phases: 2 }
    );
}

#[test]
fn composition_outside_the_fixed_mixture_is_not_tabulated() {
    let other = [0.80, 0.10, 0.05, 0.03, 0.02];
    assert_eq!(
        strict().evaluate_at(&other, 350.0, 4e6).unwrap_err(),
        TableError::CompositionMismatch
    );
    assert_eq!(
        strict().evaluate_at(&FRACTIONS, 350.0, 4e6).unwrap().source,
        Source::Table
    );
    let e = fallback().evaluate_at(&other, 350.0, 4e6).unwrap();
    assert_eq!(e.source, Source::DirectModel);
    let b = direct(
        &fallback().model().with_fractions(&other).unwrap(),
        350.0,
        4e6,
    );
    assert!(((e.properties.mass_density - b.mass_density) / b.mass_density).abs() < 1e-12);
    assert!(strict().evaluate_at(&other[..4], 350.0, 4e6).is_err());
}

#[test]
fn table_ph_agrees_with_direct_ph_flash_and_falls_back_outside() {
    let m = strict();
    for (t, p) in sample_points(9, 4) {
        let h = direct(m.model(), t, p).enthalpy;
        let e = m.ph(p, h).unwrap();
        assert_eq!(e.source, Source::Table);
        assert!(
            (e.properties.temperature_k - t).abs() < 1e-4,
            "{} vs {t}",
            e.properties.temperature_k
        );
        let flash = ph_flash(
            m.model(),
            Pressure::new(p).unwrap(),
            h,
            PhFlashOptions::default(),
        )
        .unwrap();
        assert!((e.properties.temperature_k - flash.temperature_k).abs() < 1e-4);
    }
    let h_hot = direct(m.model(), 450.0, 2e6).enthalpy;
    assert!(matches!(
        m.ph(2e6, h_hot),
        Err(TableError::EnthalpyOutsideApprovedDomain { .. })
    ));
    assert!(matches!(
        m.ph(1.2e7, 0.0),
        Err(TableError::PressureOutsideApprovedDomain { .. })
    ));
    let e = fallback().ph(2e6, h_hot).unwrap();
    assert_eq!(e.source, Source::DirectModel);
    assert!((e.properties.temperature_k - 450.0).abs() < 1e-6);
}

#[test]
fn qualification_rejects_phase_boundaries_and_the_critical_region() {
    // Crosses the dew curve of this gas (two-phase at 200 K, 30 bar).
    let crossing = TableDomain {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 300.0,
        minimum_pressure_pa: 1e6,
        maximum_pressure_pa: 5e6,
        temperature_nodes: 11,
        pressure_nodes: 5,
    };
    let e = GibbsTable::build(&model(), crossing, QualificationCriteria::default()).unwrap_err();
    assert!(matches!(e, TableError::QualificationFailed { .. }), "{e:?}");
    // The approved domain's closest approach to criticality is ≈0.62 in reduced
    // bulk modulus; a stricter bound rejects it with the critical-region reason.
    let coarse = TableDomain {
        temperature_nodes: 8,
        pressure_nodes: 8,
        ..DOMAIN
    };
    let strict_bound = QualificationCriteria {
        minimum_reduced_bulk_modulus: 0.7,
        samples_per_cell: 1,
        ..Default::default()
    };
    match GibbsTable::build(&model(), coarse, strict_bound).unwrap_err() {
        TableError::QualificationFailed { reason, .. } => {
            assert!(reason.contains("critical region"))
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn invalid_configurations_are_rejected() {
    let c = QualificationCriteria {
        samples_per_cell: 1,
        ..Default::default()
    };
    for d in [
        TableDomain {
            temperature_nodes: 1,
            ..DOMAIN
        },
        TableDomain {
            minimum_temperature_k: 430.0,
            ..DOMAIN
        },
        TableDomain {
            minimum_pressure_pa: 0.0,
            ..DOMAIN
        },
        TableDomain {
            maximum_pressure_pa: f64::NAN,
            ..DOMAIN
        },
    ] {
        assert!(matches!(
            GibbsTable::build(&model(), d, c),
            Err(TableError::InvalidDomain(_))
        ));
    }
    let composition = Composition::new(
        CompositionBasis::Mole,
        COMPONENTS
            .iter()
            .zip(FRACTIONS)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let bm = CubicModel::new(
        &thermo_data::ComponentDatabase::embedded().unwrap(),
        &composition,
        CubicVariant::PengRobinsonBostonMathias,
        vec![vec![0.0; 5]; 5],
        "synthetic zero kij",
    )
    .unwrap();
    let small = TableDomain {
        temperature_nodes: 3,
        pressure_nodes: 3,
        ..DOMAIN
    };
    assert!(matches!(
        GibbsTable::build(&bm, small, c),
        Err(TableError::Calculation(
            CalculationError::UnsupportedCapability(_)
        ))
    ));
}
