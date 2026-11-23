//! A17 conversion and yield reactors against the independent Cantera
//! reference (`validation/fixtures/reactions-v1.json`), plus failure paths.

use flowsheet_core::*;
use property_packages::CubicPackage;
use serde_json::Value;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

const FIXTURE: &str = include_str!("../../../validation/fixtures/reactions-v1.json");

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
    let any = composition(&serde_json::json!({"nitrogen": 1.0}));
    let n = components().len();
    let model = CubicModel::new(
        &ComponentDatabase::embedded().unwrap(),
        &any,
        variant,
        vec![vec![0.0; n]; n],
        "explicit synthetic zero-kij approximation; A17 reactor balances, not measured accuracy",
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

/// Feed -> unit -> product, with the given feed state.
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
    let feeds = [(f, composition(&feed["mole_fractions"]))].into();
    (fs, feeds, r)
}

fn reaction(r: &Value) -> Reaction {
    Reaction {
        coefficients: r["coefficients"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(id, v)| (id.clone(), v.as_f64().unwrap()))
            .collect(),
        extent: match r.get("conversion") {
            Some(c) => ReactionExtent::Conversion {
                key: Some(c["key"].as_str().unwrap().into()),
                fraction: c["fraction"].as_f64(),
            },
            None => ReactionExtent::MolarExtent(r["extent_mol_per_s"].as_f64()),
        },
    }
}

fn case_spec(case: &Value, thermal: SeparatorThermal) -> (UnitKind, UnitSpec) {
    let pressure_drop_pa = case["pressure_drop_Pa"].as_f64().unwrap();
    if case["kind"] == "conversion" {
        (
            UnitKind::ConversionReactor,
            UnitSpec::ConversionReactor {
                reactions: case["reactions"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(reaction)
                    .collect(),
                thermal,
                pressure_drop_pa,
            },
        )
    } else {
        (
            UnitKind::YieldReactor,
            UnitSpec::YieldReactor {
                mass_yields: case["mass_yields"]
                    .as_object()
                    .unwrap()
                    .iter()
                    .map(|(id, v)| (id.clone(), v.as_f64().unwrap()))
                    .collect(),
                inerts: case["inerts"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_str().unwrap().to_owned())
                    .collect(),
                thermal,
                pressure_drop_pa,
            },
        )
    }
}

fn outlet<'a>(fs: &Flowsheet, r: &'a Calculation, unit: UnitId) -> &'a MaterialStream {
    &r.streams[&fs.stream_at(PortRef::new(unit, 1)).unwrap().id]
}

fn inlet<'a>(fs: &Flowsheet, r: &'a Calculation, unit: UnitId) -> &'a MaterialStream {
    &r.streams[&fs.stream_at(PortRef::new(unit, 0)).unwrap().id]
}

/// Phase-weighted residual enthalpy flow n * h_res, W.
fn residual_flow(s: &MaterialStream) -> f64 {
    let state = s.state.as_ref().unwrap();
    s.molar_flow_mol_s
        * state
            .phases
            .iter()
            .map(|p| p.phase_fraction * p.state.residual_enthalpy)
            .sum::<f64>()
}

#[test]
fn formation_enthalpies_match_cantera() {
    let fx = fixture();
    let tol = fx["tolerances"]["formation_enthalpy_J_per_mol"]
        .as_f64()
        .unwrap();
    let db = ComponentDatabase::embedded().unwrap();
    for (id, h) in fx["formation_enthalpy_J_per_mol"].as_object().unwrap() {
        let ours = db.get(id).unwrap().standard_formation_enthalpy();
        assert!((ours - h.as_f64().unwrap()).abs() < tol, "{id}: {ours}");
    }
}

#[test]
fn reactor_cases_match_the_independent_reference() {
    let fx = fixture();
    let flow_tol = fx["tolerances"]["component_flow_relative"]
        .as_f64()
        .unwrap();
    let duty_tol = fx["tolerances"]["duty_relative"].as_f64().unwrap();
    let ids = components();
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        let pkg = package(variant);
        for case in fx["cases"].as_array().unwrap() {
            let name = case["id"].as_str().unwrap();
            let t_out = case["outlet_temperature_K"].as_f64().unwrap();
            let (kind, spec) = case_spec(case, SeparatorThermal::OutletTemperature(Some(t_out)));
            let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
            let r = fs
                .solve(&pkg, &feeds)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            let out = outlet(&fs, &r, unit);
            let inp = inlet(&fs, &r, unit);
            let feed_flow = case["feed"]["flow_mol_per_s"].as_f64().unwrap();
            for (id, n) in ids.iter().zip(out.component_flows_mol_s(ids.len())) {
                let expected = case["outlet_flows_mol_per_s"][id].as_f64().unwrap();
                assert!(
                    (n - expected).abs() <= flow_tol * feed_flow,
                    "{name} {id}: {n} vs {expected}"
                );
            }
            let state = out.state.as_ref().unwrap();
            assert_eq!(state.temperature_k, t_out);
            let p_in = case["feed"]["pressure_Pa"].as_f64().unwrap();
            assert_eq!(
                state.pressure_pa,
                p_in - case["pressure_drop_Pa"].as_f64().unwrap()
            );
            let unit_result = &r.units[&unit];
            let reacting = unit_result
                .balance
                .as_ref()
                .unwrap()
                .reaction
                .as_ref()
                .unwrap();
            let scale = case["duty_scale_W"].as_f64().unwrap();
            let standard = case["standard_reaction_enthalpy_W"].as_f64().unwrap();
            assert!(
                (reacting.terms.standard_reaction_enthalpy_w - standard).abs() <= duty_tol * scale,
                "{name}: standard heat {} vs {standard}",
                reacting.terms.standard_reaction_enthalpy_w
            );
            for (ours, theirs) in reacting
                .terms
                .extents_mol_s
                .iter()
                .zip(case["extents_mol_per_s"].as_array().unwrap())
            {
                assert!((ours - theirs.as_f64().unwrap()).abs() <= flow_tol * feed_flow);
            }
            // Remove the cubic departure; the remainder is the ideal-gas duty.
            let ideal = unit_result.duty_w - (residual_flow(out) - residual_flow(inp));
            let expected = case["ideal_gas_duty_W"].as_f64().unwrap();
            assert!(
                (ideal - expected).abs() <= duty_tol * scale,
                "{name} {variant:?}: ideal duty {ideal} vs {expected}"
            );
            for (element, net) in reacting
                .elements
                .iter()
                .zip(&reacting.element_residuals_mol_s)
            {
                let expected = case["element_generation_mol_per_s"][element]
                    .as_f64()
                    .unwrap();
                assert!(
                    (net - expected).abs() <= flow_tol * 10.0 * feed_flow,
                    "{name} {element}: {net} vs {expected}"
                );
            }
            let unbalanced = case["kind"] == "yield";
            if unbalanced {
                assert!(reacting.mass_residual_kg_s.abs() < 1e-12);
            } else {
                // Element-balanced, but the tabulated molar masses differ from
                // atom sums at the 1e-6 level.
                assert!(reacting.mass_residual_kg_s.abs() < 1e-6 * feed_flow * 0.05);
            }
            assert_eq!(
                unit_result
                    .warnings
                    .iter()
                    .any(|w| w.contains("not elements")),
                unbalanced,
                "{name}: {:?}",
                unit_result.warnings
            );
        }
    }
}

fn methane_case() -> Value {
    fixture()["cases"][0].clone()
}

#[test]
fn duty_and_adiabatic_specifications_are_consistent_with_outlet_temperature() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = methane_case();
    let (kind, spec) = case_spec(&case, SeparatorThermal::OutletTemperature(Some(560.0)));
    let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
    let q = fs.solve(&pkg, &feeds).unwrap().units[&unit].duty_w;
    // Duty round trip recovers the outlet temperature.
    let (kind, spec) = case_spec(&case, SeparatorThermal::Duty(Some(q)));
    let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let t = outlet(&fs, &r, unit).state.as_ref().unwrap().temperature_k;
    assert!((t - 560.0).abs() < 2e-6, "{t}");
    // A small adiabatic conversion stays in the domain; holding its outlet
    // temperature reproduces zero duty.
    let mut small = case.clone();
    small["reactions"][0]["conversion"]["fraction"] = 0.04.into();
    let (kind, spec) = case_spec(&small, SeparatorThermal::Adiabatic);
    let (fs, feeds, unit) = line(spec, kind, &small["feed"]);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let t = outlet(&fs, &r, unit).state.as_ref().unwrap().temperature_k;
    assert!(t > 450.0 && t < 600.0, "exothermic adiabatic rise: {t}");
    assert_eq!(r.units[&unit].duty_w, 0.0);
    let (kind, spec) = case_spec(&small, SeparatorThermal::OutletTemperature(Some(t)));
    let (fs, feeds, unit) = line(spec, kind, &small["feed"]);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let scale = r.units[&unit]
        .balance
        .as_ref()
        .unwrap()
        .reaction
        .as_ref()
        .unwrap()
        .terms
        .standard_reaction_enthalpy_w
        .abs();
    assert!(r.units[&unit].duty_w.abs() < 1e-6 * scale);
}

#[test]
fn adiabatic_combustion_beyond_the_package_domain_fails_visibly() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = methane_case();
    let (kind, spec) = case_spec(&case, SeparatorThermal::Adiabatic);
    let (fs, feeds, unit) = line(spec, kind, &case["feed"]);
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert_eq!(e.unit, Some(unit));
    assert!(e.message.contains("reactor outlet"), "{}", e.message);
}

fn conversion(coefficients: &[(&str, f64)], extent: ReactionExtent) -> UnitSpec {
    UnitSpec::ConversionReactor {
        reactions: vec![Reaction {
            coefficients: coefficients
                .iter()
                .map(|(id, v)| ((*id).to_owned(), *v))
                .collect(),
            extent,
        }],
        thermal: SeparatorThermal::OutletTemperature(Some(500.0)),
        pressure_drop_pa: 0.0,
    }
}

fn key(id: &str, x: f64) -> ReactionExtent {
    ReactionExtent::Conversion {
        key: Some(id.into()),
        fraction: Some(x),
    }
}

const COMBUSTION: [(&str, f64); 4] = [
    ("methane", -1.0),
    ("oxygen", -2.0),
    ("carbon-dioxide", 1.0),
    ("water", 2.0),
];

#[test]
fn chemistry_is_checked_against_the_package_before_solving() {
    let pkg = package(CubicVariant::PengRobinson);
    let feed = methane_case()["feed"].clone();
    for (spec, expected) in [
        (
            conversion(
                &[
                    ("methane", -1.0),
                    ("oxygen", -1.0),
                    ("carbon-dioxide", 1.0),
                    ("water", 2.0),
                ],
                key("methane", 0.5),
            ),
            "reaction 1 does not conserve O",
        ),
        (
            conversion(
                &[("methane", -1.0), ("hydrogen", 2.0), ("carbon", 1.0)],
                key("methane", 0.5),
            ),
            "\"hydrogen\" is not in the property package",
        ),
        (
            UnitSpec::YieldReactor {
                mass_yields: vec![("propane".into(), 1.0)],
                inerts: vec![],
                thermal: SeparatorThermal::Adiabatic,
                pressure_drop_pa: 0.0,
            },
            "\"propane\" is not in the property package",
        ),
    ] {
        let kind = match spec {
            UnitSpec::YieldReactor { .. } => UnitKind::YieldReactor,
            _ => UnitKind::ConversionReactor,
        };
        let (fs, feeds, unit) = line(spec, kind, &feed);
        let e = fs.solve(&pkg, &feeds).unwrap_err();
        assert_eq!(e.unit, Some(unit));
        assert!(e.message.contains(expected), "{}", e.message);
    }
}

#[test]
fn a_reaction_cannot_consume_more_than_is_present() {
    let pkg = package(CubicVariant::PengRobinson);
    let feed = methane_case()["feed"].clone();
    // 0.16 mol/s methane needs 0.32 mol/s oxygen; only 0.4 is fed, so an
    // extent of 0.25 mol/s would consume 0.5 mol/s oxygen.
    let spec = conversion(&COMBUSTION, ReactionExtent::MolarExtent(Some(0.25)));
    let (fs, feeds, _) = line(spec, UnitKind::ConversionReactor, &feed);
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert!(
        e.message
            .contains("reaction 1 consumes more methane than is available")
            || e.message.contains("consumes more oxygen"),
        "{}",
        e.message
    );
    // Complete conversion of the key is exact, not a shortage.
    let spec = conversion(&COMBUSTION, key("methane", 1.0));
    let (fs, feeds, unit) = line(spec, UnitKind::ConversionReactor, &feed);
    let r = fs.solve(&pkg, &feeds).unwrap();
    assert_eq!(outlet(&fs, &r, unit).component_flows_mol_s(7)[0], 0.0);
}

#[test]
fn zero_inlet_flow_passes_through_unless_an_extent_or_duty_requires_material() {
    let pkg = package(CubicVariant::PengRobinson);
    let mut feed = methane_case()["feed"].clone();
    feed["flow_mol_per_s"] = 0.0.into();
    let (fs, feeds, unit) = line(
        conversion(&COMBUSTION, key("methane", 0.5)),
        UnitKind::ConversionReactor,
        &feed,
    );
    let r = fs.solve(&pkg, &feeds).unwrap();
    assert_eq!(outlet(&fs, &r, unit).molar_flow_mol_s, 0.0);
    assert!(outlet(&fs, &r, unit).state.is_none());
    let (fs, feeds, _) = line(
        conversion(&COMBUSTION, ReactionExtent::MolarExtent(Some(0.1))),
        UnitKind::ConversionReactor,
        &feed,
    );
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert!(e.message.contains("positive molar extent"), "{}", e.message);
}

#[test]
fn reactor_specifications_report_each_problem() {
    use SpecProblem::*;
    let fields = |spec: UnitSpec| {
        spec.issues()
            .into_iter()
            .map(|i| (i.field, i.problem))
            .collect::<Vec<_>>()
    };
    assert_eq!(
        fields(UnitKind::ConversionReactor.default_spec()),
        [("outlet temperature", Missing), ("reactions", Missing)]
    );
    assert_eq!(
        fields(UnitKind::YieldReactor.default_spec()),
        [("outlet temperature", Missing), ("mass yields", Missing)]
    );
    assert_eq!(
        fields(conversion(
            &[("methane", -1.0), ("methane", 0.0), ("", 1.0)],
            key("water", 2.0)
        )),
        [
            ("reaction component", Duplicate),
            ("stoichiometric coefficient", Zero),
            ("reaction component", Missing),
            ("key component", KeyNotReactant),
            ("conversion", NotAFraction),
        ]
    );
    assert_eq!(
        fields(conversion(
            &[("methane", -1.0), ("oxygen", -2.0)],
            ReactionExtent::MolarExtent(Some(-1.0))
        )),
        [("reaction", NeedsReactantAndProduct), ("extent", Negative)]
    );
    assert_eq!(
        fields(UnitSpec::YieldReactor {
            mass_yields: vec![("methane".into(), 0.0), ("water".into(), -0.0)],
            inerts: vec!["water".into()],
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: -1.0,
        }),
        [
            ("mass yields", AllZero),
            ("inert component", Duplicate),
            ("pressure drop", Negative),
        ]
    );
}

/// A reactor inside a recycle: its reacting terms must close the balance
/// that recycle convergence re-evaluates.
#[test]
fn reactor_inside_a_recycle_closes_its_reacting_balance() {
    let pkg = package(CubicVariant::PengRobinson);
    let case = methane_case();
    let mut fs = Flowsheet::new();
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(450.0),
            pressure_pa: Some(1.5e5),
            molar_flow_mol_s: Some(2.0),
        },
    );
    let mixer = fs.add_unit(UnitKind::Mixer);
    let reactor = fs.add_unit(UnitKind::ConversionReactor);
    let (_, mut spec) = case_spec(&case, SeparatorThermal::OutletTemperature(Some(450.0)));
    if let UnitSpec::ConversionReactor {
        reactions,
        pressure_drop_pa,
        ..
    } = &mut spec
    {
        reactions[0].extent = key("methane", 0.3);
        *pressure_drop_pa = 0.0;
    }
    fs.set_spec(reactor, spec);
    let split = fs.add_unit(UnitKind::Splitter);
    fs.set_spec(
        split,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.6),
        },
    );
    let product = fs.add_unit(UnitKind::Product);
    for (from, fp, to, tp) in [
        (feed, 0, mixer, 0),
        (mixer, 2, reactor, 0),
        (reactor, 1, split, 0),
        (split, 1, product, 0),
        (split, 2, mixer, 1),
    ] {
        fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
            .unwrap();
    }
    let feeds = [(feed, composition(&case["feed"]["mole_fractions"]))].into();
    let r = fs.solve_network(&pkg, &feeds).unwrap();
    let reacting = r.units[&reactor]
        .balance
        .as_ref()
        .unwrap()
        .reaction
        .clone()
        .unwrap();
    // Overall: the product carries feed plus generation.
    let product_in = inlet(&fs, &r, product).component_flows_mol_s(7);
    let feed_out =
        r.streams[&fs.stream_at(PortRef::new(feed, 0)).unwrap().id].component_flows_mol_s(7);
    for i in 0..7 {
        let net = product_in[i] - feed_out[i] - reacting.terms.generation_mol_s[i];
        assert!(net.abs() < 1e-7, "component {i}: {net}");
    }
    assert!(reacting.terms.generation_mol_s[0] < 0.0);
}
