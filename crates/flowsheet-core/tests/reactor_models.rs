//! A17 equilibrium, Gibbs and kinetic reactors against the independent
//! Cantera/SciPy reference (`validation/fixtures/reactors-v1.json`), plus
//! consistency and failure tests.

use flowsheet_core::*;
use property_packages::CubicPackage;
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy, Temperature};

const FIXTURE: &str = include_str!("../../../validation/fixtures/reactors-v1.json");

fn fixture() -> Value {
    serde_json::from_str(FIXTURE).unwrap()
}

fn components() -> Vec<String> {
    fixture()["components"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c.as_str().unwrap().to_owned())
        .collect()
}

fn composition(fractions: &Value) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        components()
            .iter()
            .map(|id| {
                (
                    ComponentId::new(id).unwrap(),
                    fractions.get(id).and_then(Value::as_f64).unwrap_or(0.0),
                )
            })
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

fn package(variant: CubicVariant) -> CubicPackage {
    let n = components().len();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &composition(&serde_json::json!({"nitrogen": 1.0})),
        variant,
        vec![vec![0.0; n]; n],
        "explicit synthetic zero-kij approximation; A17 reactor models, not measured accuracy",
    )
    .unwrap();
    CubicPackage::new(
        model,
        PhFlashOptions {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 600.0,
            ..Default::default()
        },
    )
}

fn line(spec: UnitSpec, kind: UnitKind, feed: &Value) -> (Flowsheet, FeedCompositions, UnitId) {
    let mut fs = Flowsheet::new();
    let f = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        f,
        UnitSpec::Feed {
            temperature_k: feed["temperature_K"].as_f64(),
            pressure_pa: feed["pressure_Pa"].as_f64(),
            molar_flow_mol_s: feed["flow_mol_per_s"].as_f64(),
        },
    );
    let r = fs.add_unit(kind);
    assert!(fs.set_spec(r, spec));
    let p = fs.add_unit(UnitKind::Product);
    fs.connect(PortRef::new(f, 0), PortRef::new(r, 0)).unwrap();
    fs.connect(PortRef::new(r, 1), PortRef::new(p, 0)).unwrap();
    (fs, [(f, composition(&feed["mole_fractions"]))].into(), r)
}

fn pairs(v: &Value) -> Vec<(String, f64)> {
    v.as_object()
        .unwrap()
        .iter()
        .map(|(k, v)| (k.clone(), v.as_f64().unwrap()))
        .collect()
}

fn ideal(fugacity: FugacityBasis, case: &Value) -> EquilibriumOptions {
    EquilibriumOptions {
        fugacity,
        temperature_approach_k: case["temperature_approach_K"].as_f64().unwrap(),
    }
}

fn equilibrium_spec(
    case: &Value,
    thermal: SeparatorThermal,
    fugacity: FugacityBasis,
) -> (UnitKind, UnitSpec) {
    let pressure_drop_pa = case["pressure_drop_Pa"].as_f64().unwrap();
    if case["kind"] == "gibbs" {
        (
            UnitKind::GibbsReactor,
            UnitSpec::GibbsReactor {
                inerts: case["inerts"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_str().unwrap().into())
                    .collect(),
                equilibrium: ideal(fugacity, case),
                thermal,
                pressure_drop_pa,
            },
        )
    } else {
        (
            UnitKind::EquilibriumReactor,
            UnitSpec::EquilibriumReactor {
                reactions: case["reactions"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|r| EquilibriumReaction {
                        coefficients: pairs(&r["coefficients"]),
                        constant: match &r["constant"] {
                            Value::String(_) => EquilibriumConstant::GibbsEnergy,
                            c => EquilibriumConstant::Correlation(Correlation {
                                a: c["a"].as_f64().unwrap(),
                                b: c["b"].as_f64().unwrap(),
                                c: c["c"].as_f64().unwrap(),
                                d: c["d"].as_f64().unwrap(),
                            }),
                        },
                    })
                    .collect(),
                equilibrium: ideal(fugacity, case),
                thermal,
                pressure_drop_pa,
            },
        )
    }
}

fn arrhenius(k: &Value) -> Arrhenius {
    Arrhenius {
        pre_exponential: k["k0"].as_f64(),
        activation_energy_j_mol: k["E"].as_f64().unwrap(),
        temperature_exponent: k["n"].as_f64().unwrap(),
        reference_temperature_k: k["Tref"].as_f64(),
    }
}

fn lhhw(t: &Value) -> LhhwTerm {
    LhhwTerm {
        constant: Correlation {
            a: t["a"].as_f64().unwrap(),
            b: t["b"].as_f64().unwrap(),
            c: t["c"].as_f64().unwrap(),
            d: t["d"].as_f64().unwrap(),
        },
        orders: pairs(&t["orders"]),
    }
}

fn kinetics(case: &Value) -> Kinetics {
    Kinetics {
        reactions: case["reactions"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| KineticReaction {
                coefficients: pairs(&r["coefficients"]),
                phase: ReactingPhase::Vapor,
                basis: match r["basis"].as_str().unwrap() {
                    "partial_pressure" => ConcentrationBasis::PartialPressure,
                    "mole_fraction" => ConcentrationBasis::MoleFraction,
                    "molar_concentration" => ConcentrationBasis::MolarConcentration,
                    other => panic!("{other}"),
                },
                rate_basis: if r["rate_basis"] == "catalyst" {
                    RateBasis::CatalystMass
                } else {
                    RateBasis::ReactorVolume
                },
                rate: if r["law"] == "power" {
                    RateLaw::PowerLaw {
                        forward: RateTerm {
                            constant: arrhenius(&r["forward"]),
                            orders: pairs(&r["forward"]["orders"]),
                        },
                        reverse: (!r["reverse"].is_null()).then(|| RateTerm {
                            constant: arrhenius(&r["reverse"]),
                            orders: pairs(&r["reverse"]["orders"]),
                        }),
                    }
                } else {
                    RateLaw::Lhhw {
                        kinetic: arrhenius(&r["kinetic"]),
                        driving_forward: lhhw(&r["driving_forward"]),
                        driving_reverse: None,
                        adsorption: r["adsorption"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .map(lhhw)
                            .collect(),
                        exponent: r["exponent"].as_f64().unwrap(),
                    }
                },
            })
            .collect(),
        catalyst_bulk_density_kg_m3: case["catalyst_bulk_density_kg_m3"].as_f64(),
        liquid_volume: Default::default(),
    }
}

fn kinetic_spec(case: &Value) -> (UnitKind, UnitSpec) {
    let k = kinetics(case);
    match case["kind"].as_str().unwrap() {
        "cstr" => (
            UnitKind::Cstr,
            UnitSpec::Cstr {
                kinetics: k,
                volume_m3: case["volume_m3"].as_f64(),
                thermal: SeparatorThermal::OutletTemperature(case["outlet_temperature_K"].as_f64()),
                pressure_drop_pa: case["pressure_drop_Pa"].as_f64().unwrap(),
            },
        ),
        "pfr" => (
            UnitKind::PlugFlowReactor,
            UnitSpec::PlugFlow {
                kinetics: k,
                length_m: case["length_m"].as_f64(),
                diameter_m: case["diameter_m"].as_f64(),
                thermal: PlugFlowThermal::Isothermal(case["temperature_K"].as_f64()),
                pressure_drop_pa: case["pressure_drop_Pa"].as_f64().unwrap(),
            },
        ),
        _ => (
            UnitKind::BatchReactor,
            UnitSpec::Batch {
                kinetics: k,
                reaction_time_s: case["reaction_time_s"].as_f64(),
                thermal: BatchThermal::Isothermal(case["temperature_K"].as_f64()),
            },
        ),
    }
}

fn outlet<'a>(fs: &Flowsheet, r: &'a Calculation, unit: UnitId) -> &'a MaterialStream {
    &r.streams[&fs.stream_at(PortRef::new(unit, 1)).unwrap().id]
}

fn inlet<'a>(fs: &Flowsheet, r: &'a Calculation, unit: UnitId) -> &'a MaterialStream {
    &r.streams[&fs.stream_at(PortRef::new(unit, 0)).unwrap().id]
}

fn residual_flow(s: &MaterialStream) -> f64 {
    let state = s.state.as_ref().unwrap();
    s.molar_flow_mol_s
        * state
            .phases
            .iter()
            .map(|p| p.phase_fraction * p.state.residual_enthalpy)
            .sum::<f64>()
}

/// Solves one fixture case and compares flows, extents, standard heat of
/// reaction and the ideal-gas duty.
fn check(case: &Value, kind: UnitKind, spec: UnitSpec, pkg: &CubicPackage, flow_key: &str) {
    let fx = fixture();
    let flow_tol = fx["tolerances"][flow_key].as_f64().unwrap();
    let duty_tol = fx["tolerances"]["duty_relative"].as_f64().unwrap();
    let name = case["id"].as_str().unwrap();
    let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
    let r = fs
        .solve(pkg, &feeds)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let out = outlet(&fs, &r, unit);
    let feed_flow = case["feed"]["flow_mol_per_s"].as_f64().unwrap();
    let ids = components();
    for (id, n) in ids.iter().zip(out.component_flows_mol_s(ids.len())) {
        let expected = case["outlet_flows_mol_per_s"][id].as_f64().unwrap();
        assert!(
            (n - expected).abs() <= flow_tol * feed_flow,
            "{name} {id}: {n} vs {expected}"
        );
    }
    let u = &r.units[&unit];
    let reacting = u.balance.as_ref().unwrap().reaction.as_ref().unwrap();
    for (ours, theirs) in reacting
        .terms
        .extents_mol_s
        .iter()
        .zip(case["extents_mol_per_s"].as_array().unwrap())
    {
        let theirs = theirs.as_f64().unwrap();
        assert!(
            (ours - theirs).abs() <= flow_tol * feed_flow,
            "{name}: extent {ours} vs {theirs}"
        );
    }
    let scale = case["duty_scale_W"].as_f64().unwrap();
    let standard = case["standard_reaction_enthalpy_W"].as_f64().unwrap();
    assert!(
        (reacting.terms.standard_reaction_enthalpy_w - standard).abs() <= duty_tol * scale,
        "{name}: standard heat"
    );
    let ideal_duty = u.duty_w - (residual_flow(out) - residual_flow(inlet(&fs, &r, unit)));
    let expected = case["ideal_gas_duty_W"].as_f64().unwrap();
    assert!(
        (ideal_duty - expected).abs() <= duty_tol * scale,
        "{name}: ideal duty {ideal_duty} vs {expected}"
    );
}

#[test]
fn standard_gibbs_energies_match_cantera() {
    let fx = fixture();
    let tol = fx["tolerances"]["gibbs_energy_J_per_mol"].as_f64().unwrap();
    let db = ComponentDatabase::embedded().unwrap();
    for (id, values) in fx["standard_gibbs_energy_J_per_mol"].as_object().unwrap() {
        for (t, g) in values.as_object().unwrap() {
            let t: f64 = t.parse().unwrap();
            let ours = db
                .get(id)
                .unwrap()
                .standard_gibbs_energy(Temperature::new(t).unwrap())
                .unwrap();
            assert!((ours - g.as_f64().unwrap()).abs() < tol, "{id} {t}: {ours}");
        }
    }
}

#[test]
fn equilibrium_and_gibbs_reactors_match_the_independent_reference() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        let pkg = package(variant);
        for case in fixture()["equilibrium_cases"].as_array().unwrap() {
            let t = case["outlet_temperature_K"].as_f64();
            let (kind, spec) = equilibrium_spec(
                case,
                SeparatorThermal::OutletTemperature(t),
                FugacityBasis::IdealGas,
            );
            check(case, kind, spec, &pkg, "equilibrium_flow_relative");
        }
    }
}

#[test]
fn kinetic_reactors_match_the_independent_reference() {
    // The batch case uses the PR molar volume of the reference.
    let pkg = package(CubicVariant::PengRobinson);
    for case in fixture()["kinetic_cases"].as_array().unwrap() {
        let (kind, spec) = kinetic_spec(case);
        check(case, kind, spec, &pkg, "kinetic_flow_relative");
    }
}

/// Package fugacities satisfy the equilibrium condition at the outlet, and
/// approach the ideal-gas result at low pressure.
#[test]
fn package_fugacity_equilibrium_satisfies_the_equilibrium_condition() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = &fixture()["equilibrium_cases"][2];
    let t = case["outlet_temperature_K"].as_f64().unwrap();
    let (kind, spec) = equilibrium_spec(
        case,
        SeparatorThermal::OutletTemperature(Some(t)),
        FugacityBasis::Package,
    );
    let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let out = outlet(&fs, &r, unit).state.as_ref().unwrap();
    let x = &out.mole_fractions;
    let ln_phi = pkg.ln_fugacity_coefficients(x, t, out.pressure_pa).unwrap();
    let ids = components();
    let idx = |id: &str| ids.iter().position(|c| c == id).unwrap();
    let ln_f = |id: &str| (x[idx(id)] * out.pressure_pa / 1e5).ln() + ln_phi[idx(id)];
    let c = &case["reactions"][0]["constant"];
    let ln_k = c["a"].as_f64().unwrap() + c["b"].as_f64().unwrap() / t;
    let lhs = ln_f("methanol") - ln_f("carbon-monoxide") - 2.0 * ln_f("hydrogen");
    assert!((lhs - ln_k).abs() < 1e-9, "{lhs} vs {ln_k}");
    // Nonideality shifts methanol synthesis at 19 bar measurably.
    let ideal_methanol = case["outlet_flows_mol_per_s"]["methanol"].as_f64().unwrap();
    let methanol = outlet(&fs, &r, unit).component_flows_mol_s(ids.len())[idx("methanol")];
    assert!((methanol - ideal_methanol).abs() > 1e-4, "{methanol}");
}

/// Adiabatic and duty specifications are consistent with the
/// outlet-temperature result for every thermal-solve reactor.
#[test]
fn adiabatic_and_duty_specifications_close_the_energy_balance() {
    let pkg = package(CubicVariant::PengRobinson);
    let fx = fixture();
    type SpecFor = Box<dyn Fn(SeparatorThermal) -> UnitSpec>;
    let cases: Vec<(UnitKind, SpecFor, Value)> = vec![
        {
            let c = fx["equilibrium_cases"][1].clone();
            let cc = c.clone();
            (
                UnitKind::EquilibriumReactor,
                Box::new(move |th| equilibrium_spec(&cc, th, FugacityBasis::Package).1),
                c,
            )
        },
        {
            let c = fx["equilibrium_cases"][0].clone();
            let cc = c.clone();
            (
                UnitKind::GibbsReactor,
                Box::new(move |th| equilibrium_spec(&cc, th, FugacityBasis::Package).1),
                c,
            )
        },
        {
            let c = fx["kinetic_cases"][0].clone();
            let cc = c.clone();
            (
                UnitKind::Cstr,
                Box::new(move |th| match kinetic_spec(&cc).1 {
                    UnitSpec::Cstr {
                        kinetics,
                        volume_m3,
                        pressure_drop_pa,
                        ..
                    } => UnitSpec::Cstr {
                        kinetics,
                        volume_m3,
                        thermal: th,
                        pressure_drop_pa,
                    },
                    _ => unreachable!(),
                }),
                c,
            )
        },
    ];
    for (kind, spec, mut case) in cases {
        // Dilute with nitrogen so the adiabatic outlet stays below 600 K.
        let x = case["feed"]["mole_fractions"].as_object_mut().unwrap();
        let mut others = 0.0;
        for (id, v) in x.iter_mut() {
            if id != "nitrogen" {
                *v = (v.as_f64().unwrap() * 0.05).into();
                others += v.as_f64().unwrap();
            }
        }
        x.insert("nitrogen".into(), (1.0 - others).into());
        let t = case["outlet_temperature_K"].as_f64().unwrap();
        let (fs, feeds, unit) = line(
            spec(SeparatorThermal::OutletTemperature(Some(t))),
            kind,
            &case["feed"],
        );
        let r = fs
            .solve(&pkg, &feeds)
            .unwrap_or_else(|e| panic!("{kind}: {e}"));
        let q = r.units[&unit].duty_w;
        let flows = outlet(&fs, &r, unit).component_flows_mol_s(7);
        // The duty found at T reproduces T.
        let (fs, feeds, unit) = line(spec(SeparatorThermal::Duty(Some(q))), kind, &case["feed"]);
        let r = fs
            .solve(&pkg, &feeds)
            .unwrap_or_else(|e| panic!("{kind} duty: {e}"));
        let out = outlet(&fs, &r, unit);
        let t2 = out.state.as_ref().unwrap().temperature_k;
        assert!((t2 - t).abs() < 1e-6, "{kind}: {t2} vs {t}");
        for (a, b) in out.component_flows_mol_s(7).iter().zip(&flows) {
            assert!((a - b).abs() < 1e-8, "{kind}: {a} vs {b}");
        }
        // Adiabatic: zero duty and a closed balance.
        let (fs, feeds, unit) = line(spec(SeparatorThermal::Adiabatic), kind, &case["feed"]);
        let r = fs
            .solve(&pkg, &feeds)
            .unwrap_or_else(|e| panic!("{kind} adiabatic: {e}"));
        assert_eq!(r.units[&unit].duty_w, 0.0);
        let b = r.units[&unit].balance.as_ref().unwrap();
        assert!(b.scaled_energy_residual.abs() < 2e-9, "{kind}");
    }
}

#[test]
fn plug_flow_thermal_modes_are_consistent() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = fixture()["kinetic_cases"][1].clone();
    let with = |thermal| match kinetic_spec(&case).1 {
        UnitSpec::PlugFlow {
            kinetics,
            length_m,
            diameter_m,
            pressure_drop_pa,
            ..
        } => UnitSpec::PlugFlow {
            kinetics,
            length_m,
            diameter_m,
            thermal,
            pressure_drop_pa,
        },
        _ => unreachable!(),
    };
    let solve = |thermal| {
        let (fs, feeds, unit) = line(with(thermal), UnitKind::PlugFlowReactor, &case["feed"]);
        let r = fs.solve(&pkg, &feeds).unwrap();
        let out = outlet(&fs, &r, unit).clone();
        (out, r.units[&unit].duty_w)
    };
    // No wall heat transfer is the adiabatic reactor.
    let (adiabatic, q) = solve(PlugFlowThermal::Adiabatic);
    assert_eq!(q, 0.0);
    let (coolant, q0) = solve(PlugFlowThermal::Coolant {
        temperature_k: Some(300.0),
        u_w_m2_k: Some(0.0),
    });
    assert_eq!(q0, 0.0);
    let t = |s: &MaterialStream| s.state.as_ref().unwrap().temperature_k;
    assert!((t(&adiabatic) - t(&coolant)).abs() < 1e-9);
    // A coolant above the feed temperature heats the gas: positive duty and a
    // hotter outlet than adiabatic (the exothermic reactions may overshoot
    // the coolant temperature).
    let (cooled, q_cooled) = solve(PlugFlowThermal::Coolant {
        temperature_k: Some(560.0),
        u_w_m2_k: Some(200.0),
    });
    assert!(q_cooled > 0.0, "{q_cooled}");
    assert!(t(&cooled) > t(&adiabatic), "{}", t(&cooled));
    // The hotter gas methanates more than the adiabatic reactor.
    let methane = |s: &MaterialStream| s.component_flows_mol_s(7)[0];
    assert!(methane(&cooled) > methane(&adiabatic));
}

#[test]
fn reactor_model_failures_are_reported() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = fixture()["equilibrium_cases"][1].clone();
    // Linearly dependent reactions.
    let (kind, mut spec) = equilibrium_spec(
        &case,
        SeparatorThermal::OutletTemperature(Some(520.0)),
        FugacityBasis::IdealGas,
    );
    if let UnitSpec::EquilibriumReactor { reactions, .. } = &mut spec {
        let mut double = reactions[0].clone();
        for c in &mut double.coefficients {
            c.1 *= 2.0;
        }
        reactions.push(double);
    }
    let (fs, feeds, _) = line(spec, kind, &case["feed"]);
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert!(e.message.contains("linearly dependent"), "{}", e.message);
    // Catalyst-mass rates need a bulk density.
    let k = fixture()["kinetic_cases"][1].clone();
    let (_, mut spec) = kinetic_spec(&k);
    if let UnitSpec::PlugFlow { kinetics, .. } = &mut spec {
        kinetics.catalyst_bulk_density_kg_m3 = None;
    }
    assert!(
        spec.issues()
            .iter()
            .any(|i| i.field == "catalyst bulk density" && i.problem == SpecProblem::Missing)
    );
    // An order on a component that is not in the package.
    let (kind, mut spec) = kinetic_spec(&fixture()["kinetic_cases"][0]);
    if let UnitSpec::Cstr { kinetics, .. } = &mut spec
        && let RateLaw::PowerLaw { forward, .. } = &mut kinetics.reactions[0].rate
    {
        forward.orders.push(("ammonia".into(), 1.0));
    }
    let (fs, feeds, _) = line(spec, kind, &fixture()["kinetic_cases"][0]["feed"]);
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert!(
        e.message
            .contains("\"ammonia\" is not in the property package"),
        "{}",
        e.message
    );
}

#[test]
fn defaults_list_the_required_inputs() {
    use SpecProblem::Missing;
    let fields = |kind: UnitKind| {
        kind.default_spec()
            .issues()
            .into_iter()
            .map(|i| (i.field, i.problem))
            .collect::<Vec<_>>()
    };
    assert_eq!(
        fields(UnitKind::EquilibriumReactor),
        [("outlet temperature", Missing), ("reactions", Missing)]
    );
    assert_eq!(
        fields(UnitKind::GibbsReactor),
        [("outlet temperature", Missing)]
    );
    assert_eq!(
        fields(UnitKind::Cstr),
        [
            ("volume", Missing),
            ("outlet temperature", Missing),
            ("reactions", Missing)
        ]
    );
    assert_eq!(
        fields(UnitKind::PlugFlowReactor),
        [
            ("length", Missing),
            ("diameter", Missing),
            ("temperature", Missing),
            ("reactions", Missing)
        ]
    );
    assert_eq!(
        fields(UnitKind::BatchReactor),
        [
            ("reaction time", Missing),
            ("temperature", Missing),
            ("reactions", Missing)
        ]
    );
}
