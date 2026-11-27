//! Format 4 reactors: specifications and reacting results survive a save and
//! load, and a reactor in a file declaring an older format is rejected.

use std::ops::ControlFlow;

use flowsheet_core::{PortRef, Reaction, ReactionExtent, SeparatorThermal, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, Currency, PackageSelection, Project, ProjectError, solve_project,
};

const COMPONENTS: [&str; 5] = ["methane", "oxygen", "nitrogen", "carbon-dioxide", "water"];

fn reactor_project(spec: UnitSpec, kind: UnitKind) -> Project {
    let mut project = Project::new(PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        COMPONENTS.iter().map(|c| c.to_string()).collect(),
        "explicit synthetic zero-kij approximation; A17 reactor persistence test",
    ));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(450.0),
            pressure_pa: Some(1.5e5),
            molar_flow_mol_s: Some(2.0),
        },
    );
    let reactor = fs.add_unit(kind);
    fs.set_spec(reactor, spec);
    let product = fs.add_unit(UnitKind::Product);
    fs.connect(PortRef::new(feed, 0), PortRef::new(reactor, 0))
        .unwrap();
    fs.connect(PortRef::new(reactor, 1), PortRef::new(product, 0))
        .unwrap();
    project.feeds.insert(feed, vec![0.08, 0.2, 0.72, 0.0, 0.0]);
    project
}

fn combustion() -> UnitSpec {
    UnitSpec::ConversionReactor {
        reactions: vec![Reaction {
            coefficients: vec![
                ("methane".into(), -1.0),
                ("oxygen".into(), -2.0),
                ("carbon-dioxide".into(), 1.0),
                ("water".into(), 2.0),
            ],
            extent: ReactionExtent::Conversion {
                key: Some("methane".into()),
                fraction: Some(0.75),
            },
        }],
        thermal: SeparatorThermal::OutletTemperature(Some(560.0)),
        pressure_drop_pa: 2.0e4,
    }
}

fn yields() -> UnitSpec {
    UnitSpec::YieldReactor {
        mass_yields: vec![("carbon-dioxide".into(), 2.75), ("water".into(), 2.25)],
        inerts: vec!["nitrogen".into()],
        thermal: SeparatorThermal::OutletTemperature(Some(500.0)),
        pressure_drop_pa: 0.0,
    }
}

#[test]
fn reactor_specifications_and_results_round_trip() {
    for (spec, kind) in [
        (combustion(), UnitKind::ConversionReactor),
        (yields(), UnitKind::YieldReactor),
    ] {
        let mut project = reactor_project(spec.clone(), kind);
        let package = project.package.prepare().unwrap();
        let results = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
        let reactor = project
            .flowsheet
            .units()
            .find(|u| u.kind == kind)
            .unwrap()
            .id;
        let stored = results.units[&reactor].reaction.clone().unwrap();
        assert_eq!(stored.generation_mol_s.len(), COMPONENTS.len());
        assert!(stored.standard_reaction_enthalpy_w < 0.0, "{kind}");
        assert_eq!(
            stored.extents_mol_s.len(),
            usize::from(kind == UnitKind::ConversionReactor)
        );
        project.results = Some(results);
        let text = project.to_json().unwrap();
        assert!(text.contains("\"format_version\": 4"));
        let loaded = Project::from_json(&text).unwrap();
        assert_eq!(loaded.flowsheet.unit(reactor).unwrap().spec, spec);
        assert_eq!(loaded.results, project.results);
        assert_eq!(
            loaded
                .results
                .as_ref()
                .unwrap()
                .unit_currency(reactor, &loaded.fingerprints()),
            Currency::Current
        );
        assert_eq!(loaded.to_json().unwrap(), text);

        // Editing the chemistry makes the stored result stale.
        let mut edited = loaded.clone();
        let mut spec = edited.flowsheet.unit(reactor).unwrap().spec.clone();
        match &mut spec {
            UnitSpec::ConversionReactor { reactions, .. } => {
                reactions[0].extent = ReactionExtent::MolarExtent(Some(0.01));
            }
            UnitSpec::YieldReactor { mass_yields, .. } => mass_yields[0].1 = 3.0,
            _ => unreachable!(),
        }
        edited.flowsheet.set_spec(reactor, spec);
        assert_eq!(
            edited
                .results
                .as_ref()
                .unwrap()
                .unit_currency(reactor, &edited.fingerprints()),
            Currency::Stale
        );

        // A reactor cannot appear in a file that declares an older format.
        for old in [1, 2, 3] {
            let older = text.replacen(
                "\"format_version\": 4",
                &format!("\"format_version\": {old}"),
                1,
            );
            let err = Project::from_json(&older).unwrap_err();
            assert!(
                matches!(&err, ProjectError::Invalid(m) if m.contains("requires format version 4")),
                "{err}"
            );
        }
    }
}

#[test]
fn malformed_reactor_records_and_values_are_rejected() {
    let project = reactor_project(combustion(), UnitKind::ConversionReactor);
    let text = project.to_json().unwrap();
    let edit = |from: &str, to: &str| {
        assert!(text.contains(from), "{from}");
        Project::from_json(&text.replacen(from, to, 1))
    };
    assert!(matches!(
        edit("\"coefficient\": -1.0", "\"coefficient\": \"x\""),
        Err(ProjectError::Format(_))
    ));
    assert!(matches!(
        edit("\"component\": \"methane\"", "\"species\": \"methane\""),
        Err(ProjectError::Format(_))
    ));
    // Nonfinite chemistry cannot be written.
    let mut bad = project.clone();
    let reactor = bad
        .flowsheet
        .units()
        .find(|u| u.kind == UnitKind::ConversionReactor)
        .unwrap()
        .id;
    let mut spec = combustion();
    if let UnitSpec::ConversionReactor { reactions, .. } = &mut spec {
        reactions[0].coefficients[0].1 = f64::NAN;
    }
    bad.flowsheet.set_spec(reactor, spec);
    assert!(matches!(bad.to_json(), Err(ProjectError::Invalid(m)) if m.contains("not finite")));
}

/// Every option of the equilibrium, Gibbs and kinetic reactors survives a
/// save and load unchanged.
#[test]
fn equilibrium_and_kinetic_specifications_round_trip() {
    use flowsheet_core::{
        Arrhenius, BatchThermal, ConcentrationBasis, Correlation, EquilibriumConstant,
        EquilibriumOptions, EquilibriumReaction, FugacityBasis, KineticReaction, Kinetics,
        LhhwTerm, PlugFlowThermal, RateBasis, RateLaw, RateTerm, ReactingPhase,
    };
    let shift = vec![
        ("carbon-dioxide".to_string(), -1.0),
        ("water".to_string(), 1.0),
    ];
    let arrhenius = Arrhenius {
        pre_exponential: Some(2.5e-3),
        activation_energy_j_mol: 6.0e4,
        temperature_exponent: 0.5,
        reference_temperature_k: Some(500.0),
    };
    let lhhw = |a: f64| LhhwTerm {
        constant: Correlation {
            a,
            b: 1000.0,
            c: -0.5,
            d: 1e-3,
        },
        orders: vec![("methane".into(), 0.5)],
    };
    let kinetics = Kinetics {
        reactions: vec![
            KineticReaction {
                coefficients: shift.clone(),
                phase: ReactingPhase::Liquid,
                basis: ConcentrationBasis::Fugacity,
                rate_basis: RateBasis::CatalystMass,
                rate: RateLaw::PowerLaw {
                    forward: RateTerm {
                        constant: arrhenius,
                        orders: vec![("water".into(), 1.5)],
                    },
                    reverse: Some(RateTerm {
                        constant: Arrhenius {
                            pre_exponential: None,
                            ..arrhenius
                        },
                        orders: vec![],
                    }),
                },
            },
            KineticReaction {
                coefficients: shift.clone(),
                phase: ReactingPhase::Vapor,
                basis: ConcentrationBasis::MassFraction,
                rate_basis: RateBasis::ReactorVolume,
                rate: RateLaw::Lhhw {
                    kinetic: arrhenius,
                    driving_forward: lhhw(1.0),
                    driving_reverse: Some(lhhw(2.0)),
                    adsorption: vec![lhhw(3.0), lhhw(4.0)],
                    exponent: 2.0,
                },
            },
        ],
        catalyst_bulk_density_kg_m3: Some(900.0),
        liquid_volume: Default::default(),
    };
    let equilibrium = EquilibriumOptions {
        fugacity: FugacityBasis::IdealGas,
        temperature_approach_k: -15.0,
    };
    let specs = [
        (
            UnitKind::EquilibriumReactor,
            UnitSpec::EquilibriumReactor {
                reactions: vec![
                    EquilibriumReaction {
                        coefficients: shift.clone(),
                        constant: EquilibriumConstant::GibbsEnergy,
                    },
                    EquilibriumReaction {
                        coefficients: shift.clone(),
                        constant: EquilibriumConstant::Correlation(Correlation {
                            a: 1.0,
                            b: 2.0,
                            c: 3.0,
                            d: 4.0,
                        }),
                    },
                ],
                equilibrium,
                thermal: SeparatorThermal::Duty(Some(-1.0e3)),
                pressure_drop_pa: 5.0e3,
            },
        ),
        (
            UnitKind::GibbsReactor,
            UnitSpec::GibbsReactor {
                inerts: vec!["nitrogen".into()],
                equilibrium,
                thermal: SeparatorThermal::Adiabatic,
                pressure_drop_pa: 0.0,
            },
        ),
        (
            UnitKind::Cstr,
            UnitSpec::Cstr {
                kinetics: kinetics.clone(),
                volume_m3: Some(2.0),
                thermal: SeparatorThermal::OutletTemperature(Some(450.0)),
                pressure_drop_pa: 1.0e3,
            },
        ),
        (
            UnitKind::PlugFlowReactor,
            UnitSpec::PlugFlow {
                kinetics: kinetics.clone(),
                length_m: Some(3.0),
                diameter_m: None,
                thermal: PlugFlowThermal::Coolant {
                    temperature_k: Some(400.0),
                    u_w_m2_k: Some(150.0),
                },
                pressure_drop_pa: 0.0,
            },
        ),
        (
            UnitKind::BatchReactor,
            UnitSpec::Batch {
                kinetics,
                reaction_time_s: Some(3600.0),
                thermal: BatchThermal::Adiabatic,
            },
        ),
    ];
    for (kind, spec) in specs {
        let project = reactor_project(spec.clone(), kind);
        let text = project.to_json().unwrap();
        let loaded = Project::from_json(&text).unwrap();
        let unit = loaded.flowsheet.units().find(|u| u.kind == kind).unwrap();
        assert_eq!(unit.spec, spec, "{kind}");
        assert_eq!(loaded.to_json().unwrap(), text, "{kind}");
    }
}

#[test]
fn a_solved_cstr_stores_its_extents() {
    use flowsheet_core::{
        Arrhenius, ConcentrationBasis, KineticReaction, Kinetics, RateBasis, RateLaw, RateTerm,
        ReactingPhase,
    };
    let spec = UnitSpec::Cstr {
        kinetics: Kinetics {
            reactions: vec![KineticReaction {
                coefficients: vec![
                    ("methane".into(), -1.0),
                    ("oxygen".into(), -2.0),
                    ("carbon-dioxide".into(), 1.0),
                    ("water".into(), 2.0),
                ],
                phase: ReactingPhase::Vapor,
                basis: ConcentrationBasis::MolarConcentration,
                rate_basis: RateBasis::ReactorVolume,
                rate: RateLaw::PowerLaw {
                    forward: RateTerm {
                        constant: Arrhenius {
                            pre_exponential: Some(0.05),
                            ..Default::default()
                        },
                        orders: vec![("methane".into(), 1.0)],
                    },
                    reverse: None,
                },
            }],
            catalyst_bulk_density_kg_m3: None,
            liquid_volume: Default::default(),
        },
        volume_m3: Some(1.0),
        thermal: SeparatorThermal::OutletTemperature(Some(500.0)),
        pressure_drop_pa: 0.0,
    };
    let mut project = reactor_project(spec, UnitKind::Cstr);
    let package = project.package.prepare().unwrap();
    let results = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    let reactor = project
        .flowsheet
        .units()
        .find(|u| u.kind == UnitKind::Cstr)
        .unwrap()
        .id;
    let stored = results.units[&reactor].reaction.clone().unwrap();
    assert_eq!(stored.extents_mol_s.len(), 1);
    assert!(stored.extents_mol_s[0] > 0.0 && stored.extents_mol_s[0] < 0.16);
    project.results = Some(results);
    let loaded = Project::from_json(&project.to_json().unwrap()).unwrap();
    assert_eq!(loaded.results, project.results);
}
