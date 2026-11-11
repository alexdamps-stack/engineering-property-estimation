#[path = "../examples/support/mod.rs"]
mod support;

use flowsheet_core::*;
use thermo_models::CubicVariant;
use thermo_types::{Composition, CompositionBasis, NormalizationPolicy};

fn single(
    kind: UnitKind,
    spec: UnitSpec,
    t: f64,
    flow: f64,
) -> (
    Flowsheet,
    property_packages::CubicPackage,
    FeedCompositions,
    UnitId,
) {
    let (_, pkg, feeds) = support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
    let composition = feeds.values().next().unwrap().clone();
    let mut fs = Flowsheet::new();
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(flow),
        },
    );
    let u = fs.add_unit(kind);
    fs.set_spec(u, spec);
    fs.connect(PortRef::new(feed, 0), PortRef::new(u, 0))
        .unwrap();
    for (i, p) in kind.ports().iter().enumerate() {
        if p.direction == PortDirection::Outlet {
            let product = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(u, i), PortRef::new(product, 0))
                .unwrap();
        }
    }
    (fs, pkg, [(feed, composition)].into(), u)
}

#[test]
fn specified_duty_runs_the_independently_qualified_demonstration() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    for (case, variant) in fixture["cases"].as_array().unwrap().iter().zip([
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ]) {
        let (mut fs, pkg, feeds) = support::demonstration(variant, 1.0);
        let cooler = unit(&fs, UnitKind::Cooler);
        let q = case["cooler_duty_W"].as_f64().unwrap();
        fs.set_spec(
            cooler,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(q)),
                pressure_drop_pa: 0.0,
            },
        );
        let valve = unit(&fs, UnitKind::Valve);
        fs.set_spec(
            valve,
            UnitSpec::Valve {
                pressure: PressureSpec::PressureDrop(Some(3.5e6)),
            },
        );
        let mixer = unit(&fs, UnitKind::Mixer);
        fs.set_spec(
            mixer,
            UnitSpec::Mixer {
                outlet_pressure: MixerPressure::Specified(Some(5e6)),
            },
        );
        let r = fs.solve(&pkg, &feeds).unwrap();
        assert!(
            (stream(&fs, &r, "cool")
                .state
                .as_ref()
                .unwrap()
                .temperature_k
                - 240.0)
                .abs()
                < 2e-6
        );
        assert!(
            (stream(&fs, &r, "expand")
                .state
                .as_ref()
                .unwrap()
                .temperature_k
                - case["streams"][4]["T"].as_f64().unwrap())
            .abs()
                < 2e-6
        );
        assert_eq!(r.units[&cooler].duty_w, q);
    }
}

#[test]
fn heater_temperature_duty_and_pressure_drop_obey_energy_balance() {
    for thermal in [
        ThermalSpec::OutletTemperature(Some(350.0)),
        ThermalSpec::Duty(Some(20_000.0)),
    ] {
        let (fs, pkg, feeds, u) = single(
            UnitKind::Heater,
            UnitSpec::HeatExchange {
                thermal,
                pressure_drop_pa: 1e5,
            },
            330.0,
            10.0,
        );
        let r = fs.solve(&pkg, &feeds).unwrap();
        let out = &r.streams[&fs.stream_at(PortRef::new(u, 1)).unwrap().id];
        assert_eq!(out.state.as_ref().unwrap().pressure_pa, 4.9e6);
        assert!(r.units[&u].duty_w > 0.0);
        assert!(
            r.units[&u]
                .balance
                .as_ref()
                .unwrap()
                .scaled_energy_residual
                .abs()
                < 1e-9
        );
        match thermal {
            ThermalSpec::OutletTemperature(Some(t)) => {
                assert_eq!(out.state.as_ref().unwrap().temperature_k, t)
            }
            ThermalSpec::Duty(Some(q)) => assert_eq!(r.units[&u].duty_w, q),
            _ => unreachable!(),
        }
    }
}

#[test]
fn separator_thermal_modes_and_empty_outlet_have_closed_balances() {
    for thermal in [
        SeparatorThermal::Adiabatic,
        SeparatorThermal::OutletTemperature(Some(240.0)),
        SeparatorThermal::Duty(Some(-10_000.0)),
    ] {
        let (fs, pkg, feeds, u) = single(
            UnitKind::Separator,
            UnitSpec::Separator {
                thermal,
                pressure_drop_pa: 1e5,
            },
            330.0,
            10.0,
        );
        let r = fs.solve(&pkg, &feeds).unwrap();
        assert!(
            r.units[&u]
                .balance
                .as_ref()
                .unwrap()
                .max_scaled_component_residual
                < 1e-9
        );
        assert!(
            r.units[&u]
                .balance
                .as_ref()
                .unwrap()
                .scaled_energy_residual
                .abs()
                < 1e-9
        );
        for s in fs.streams().filter(|s| s.from.unit == u) {
            let out = &r.streams[&s.id];
            if let Some(a) = &out.state {
                assert_eq!(a.pressure_pa, 4.9e6);
                assert_eq!(a.phases.len(), 1);
                assert_eq!(a.phases[0].phase_fraction, 1.0);
            } else {
                assert_eq!(out.molar_flow_mol_s, 0.0);
                assert_eq!(out.enthalpy_flow_w(), 0.0);
                assert_eq!(r.units[&u].warnings.len(), 1);
            }
        }
    }
    let (fs, pkg, feeds, u) = single(
        UnitKind::Separator,
        UnitSpec::Separator {
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
        330.0,
        10.0,
    );
    let r = fs.solve(&pkg, &feeds).unwrap();
    let absent = &r.streams[&fs.stream_at(PortRef::new(u, 2)).unwrap().id];
    assert_eq!(absent.molar_flow_mol_s, 0.0);
    assert!(absent.state.is_none());
    assert_eq!(r.units[&u].warnings.len(), 1);
}

#[test]
fn splitter_endpoints_and_branch_rejoin_preserve_states_and_balances() {
    for fraction in [0.0, 0.3, 1.0] {
        let (mut fs, pkg, feeds, u) = single(
            UnitKind::Splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(fraction),
            },
            330.0,
            10.0,
        );
        // Rejoin both branches into a mixer, including an empty branch.
        let ends: Vec<_> = fs
            .units()
            .filter(|u| u.kind == UnitKind::Product)
            .map(|u| u.id)
            .collect();
        for id in ends {
            fs.remove_unit(id);
        }
        let mix = fs.add_unit(UnitKind::Mixer);
        for i in 0..2 {
            fs.connect(PortRef::new(u, i + 1), PortRef::new(mix, i))
                .unwrap();
        }
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(PortRef::new(mix, 2), PortRef::new(product, 0))
            .unwrap();
        let r = fs.solve(&pkg, &feeds).unwrap();
        let a = &r.streams[&fs.stream_at(PortRef::new(u, 1)).unwrap().id];
        let b = &r.streams[&fs.stream_at(PortRef::new(u, 2)).unwrap().id];
        assert_eq!(a.molar_flow_mol_s, 10.0 * fraction);
        assert_eq!(b.molar_flow_mol_s, 10.0 * (1.0 - fraction));
        assert_eq!(
            a.state.as_ref().unwrap().enthalpy(),
            b.state.as_ref().unwrap().enthalpy()
        );
        let out = &r.streams[&fs.stream_at(PortRef::new(mix, 2)).unwrap().id];
        assert!((out.state.as_ref().unwrap().temperature_k - 330.0).abs() < 2e-6);
    }
}

#[test]
fn exchanger_conserves_energy_for_all_specifications() {
    for thermal in [
        ExchangerSpec::Duty(Some(1000.0)),
        ExchangerSpec::HotOutletTemperature(Some(340.0)),
        ExchangerSpec::ColdOutletTemperature(Some(332.0)),
    ] {
        let (mut fs, pkg, feeds) =
            support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
        let feed_ids: Vec<_> = feeds.keys().copied().collect();
        let remove: Vec<_> = fs
            .units()
            .filter(|u| u.kind != UnitKind::Feed)
            .map(|u| u.id)
            .collect();
        for id in remove {
            fs.remove_unit(id);
        }
        let ex = fs.add_unit(UnitKind::HeatExchanger);
        fs.set_spec(
            ex,
            UnitSpec::HeatExchanger {
                thermal,
                hot_pressure_drop_pa: 1e5,
                cold_pressure_drop_pa: 2e5,
            },
        );
        fs.connect(PortRef::new(feed_ids[1], 0), PortRef::new(ex, 0))
            .unwrap();
        fs.connect(PortRef::new(feed_ids[0], 0), PortRef::new(ex, 1))
            .unwrap();
        for i in 2..4 {
            let product = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(ex, i), PortRef::new(product, 0))
                .unwrap();
        }
        let r = fs.solve(&pkg, &feeds).unwrap();
        // The transferred heat is internal and recorded separately (#119).
        assert_eq!(r.units[&ex].duty_w, 0.0);
        let transferred = r.units[&ex].transferred_duty_w.unwrap();
        assert!(transferred > 0.0);
        let hot_in = &r.streams[&fs.stream_at(PortRef::new(ex, 0)).unwrap().id];
        let hot_out = &r.streams[&fs.stream_at(PortRef::new(ex, 2)).unwrap().id];
        let released = hot_in.enthalpy_flow_w() - hot_out.enthalpy_flow_w();
        assert!((released / transferred - 1.0).abs() < 1e-8);
        if let ExchangerSpec::Duty(Some(q)) = thermal {
            assert_eq!(transferred, q);
        }
        let balance = r.units[&ex].balance.as_ref().unwrap();
        assert!(
            balance.max_scaled_component_residual < 1e-9
                && balance.scaled_energy_residual.abs() < 1e-9
        );
    }
}

#[test]
fn invalid_specs_zero_flow_domain_and_unsupported_units_fail_explicitly() {
    // A zero-flow heater inlet is passed through (#89); see
    // `zero_flow_streams_pass_through_downstream_units`.
    let cases = [
        (
            UnitKind::Heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(450.0)),
                pressure_drop_pa: 0.0,
            },
            330.0,
            1.0,
            "search interval",
        ),
        (
            UnitKind::Heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(f64::NAN)),
                pressure_drop_pa: 0.0,
            },
            330.0,
            1.0,
            "finite",
        ),
        (
            UnitKind::Valve,
            UnitSpec::Valve {
                pressure: PressureSpec::PressureDrop(Some(5e6)),
            },
            330.0,
            1.0,
            "nonpositive",
        ),
        (
            UnitKind::Valve,
            UnitSpec::Valve {
                pressure: PressureSpec::OutletPressure(None),
            },
            330.0,
            1.0,
            "missing",
        ),
        (
            UnitKind::Pump,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(6e6)),
                efficiency: 0.75,
            },
            330.0,
            1.0,
            "PS flash",
        ),
        (
            UnitKind::Splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(1.1),
            },
            330.0,
            1.0,
            "between 0 and 1",
        ),
    ];
    for (kind, spec, t, n, message) in cases {
        let (fs, pkg, feeds, _) = single(kind, spec, t, n);
        let e = fs.solve(&pkg, &feeds).unwrap_err();
        assert!(e.message.contains(message), "{}", e.message);
    }
    let (mut fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let mixer = unit(&fs, UnitKind::Mixer);
    fs.set_spec(
        mixer,
        UnitSpec::Mixer {
            outlet_pressure: MixerPressure::Specified(Some(6e6)),
        },
    );
    assert_eq!(fs.solve(&pkg, &feeds).unwrap_err().unit, Some(mixer));
    let (_, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    assert!(
        Flowsheet::new()
            .solve(&pkg, &feeds)
            .unwrap_err()
            .message
            .contains("no feed")
    );
}

#[test]
fn ph_nonconvergence_is_attached_to_its_unit() {
    let (fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let pkg = property_packages::CubicPackage::new(
        pkg.model().clone(),
        thermo_solvers::PhFlashOptions {
            max_iterations: 0,
            ..Default::default()
        },
    );
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert_eq!(e.unit, Some(unit(&fs, UnitKind::Mixer)));
    assert!(e.message.contains("iteration limit"));
}

#[test]
fn cycles_and_disconnected_ports_are_rejected_before_execution() {
    let (mut fs, _, _) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let connection = fs.streams().next().unwrap().id;
    fs.disconnect(connection);
    assert!(!fs.calculation_order().unwrap_err().diagnostics.is_empty());
    let mut fs = Flowsheet::new();
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(1.0),
        },
    );
    let mix = fs.add_unit(UnitKind::Mixer);
    let split = fs.add_unit(UnitKind::Splitter);
    fs.set_spec(
        split,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.5),
        },
    );
    let product = fs.add_unit(UnitKind::Product);
    fs.connect(PortRef::new(feed, 0), PortRef::new(mix, 0))
        .unwrap();
    fs.connect(PortRef::new(mix, 2), PortRef::new(split, 0))
        .unwrap();
    fs.connect(PortRef::new(split, 1), PortRef::new(mix, 1))
        .unwrap();
    fs.connect(PortRef::new(split, 2), PortRef::new(product, 0))
        .unwrap();
    assert!(
        fs.calculation_order()
            .unwrap_err()
            .diagnostics
            .iter()
            .any(|d| matches!(d.kind, DiagnosticKind::RecycleUnsupported(_)))
    );
}

fn unit(fs: &Flowsheet, kind: UnitKind) -> UnitId {
    fs.units().find(|u| u.kind == kind).unwrap().id
}

fn stream<'a>(fs: &Flowsheet, r: &'a Calculation, name: &str) -> &'a MaterialStream {
    &r.streams[&fs.streams().find(|s| s.name == name).unwrap().id]
}

#[test]
fn connected_demonstration_matches_independent_references_and_scales() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../../validation/fixtures/process-v1.json")).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let variant = if case["variant"] == "PR" {
            CubicVariant::PengRobinson
        } else {
            CubicVariant::PengRobinsonBostonMathias
        };
        let close = |key: &str, actual: f64, expected: f64| {
            let budget = &fixture["tolerances"][key];
            assert!(
                actual.is_finite()
                    && (actual - expected).abs()
                        <= budget[0].as_f64().unwrap()
                            + budget[1].as_f64().unwrap() * expected.abs(),
                "{} {key}: {actual} != {expected}",
                case["variant"]
            );
        };
        for scale in [0.01, 1.0, 100.0] {
            let (fs, pkg, feeds) = support::demonstration(variant, scale);
            let r = fs.solve(&pkg, &feeds).unwrap();
            assert!(r.is_current(&fs, &feeds, &pkg));
            assert_eq!(r.order.len(), 8);
            for connection in fs.streams() {
                assert!(
                    r.order.iter().position(|id| *id == connection.from.unit)
                        < r.order.iter().position(|id| *id == connection.to.unit)
                );
            }
            for expected in case["streams"].as_array().unwrap() {
                let s = stream(&fs, &r, expected["id"].as_str().unwrap());
                let a = s.state.as_ref().unwrap();
                close(
                    "temperature_K",
                    a.temperature_k,
                    expected["T"].as_f64().unwrap(),
                );
                assert_eq!(a.pressure_pa, expected["p"].as_f64().unwrap());
                close(
                    "phase_fraction",
                    s.molar_flow_mol_s / (120.0 * scale),
                    expected["flow_mol_per_s"].as_f64().unwrap() / 120.0,
                );
                close(
                    "enthalpy_J_per_mol",
                    a.enthalpy(),
                    expected["target_h"].as_f64().unwrap(),
                );
                assert_eq!(
                    a.phases.len(),
                    expected["phases"].as_u64().unwrap() as usize
                );
                for (x, z) in a
                    .mole_fractions
                    .iter()
                    .zip(expected["z"].as_array().unwrap())
                {
                    close("composition", *x, z.as_f64().unwrap());
                }
                if a.phases.len() == 2 {
                    close(
                        "phase_fraction",
                        a.phases[0].phase_fraction,
                        expected["vapor_fraction"].as_f64().unwrap(),
                    );
                    for (phase, key) in a.phases.iter().zip(["vapor", "liquid"]) {
                        for (x, z) in phase
                            .composition
                            .iter()
                            .zip(expected[key].as_array().unwrap())
                        {
                            close("composition", *x, z.as_f64().unwrap());
                        }
                    }
                }
                for p in &a.phases {
                    assert!(
                        !p.state.model.identifier.is_empty()
                            && !p.state.data.version.is_empty()
                            && !p.state.caloric_reference.identifier.is_empty()
                    );
                }
            }
            let cooler = unit(&fs, UnitKind::Cooler);
            close(
                "enthalpy_J_per_mol",
                r.units[&cooler].duty_w / (120.0 * scale),
                case["cooler_duty_W"].as_f64().unwrap() / 120.0,
            );
            close(
                "enthalpy_J_per_mol",
                stream(&fs, &r, "expand").state.as_ref().unwrap().enthalpy(),
                stream(&fs, &r, "cool").state.as_ref().unwrap().enthalpy(),
            );
            for b in r.units.values().filter_map(|u| u.balance.as_ref()) {
                assert!(
                    b.max_scaled_component_residual < 1e-9 && b.scaled_energy_residual.abs() < 1e-9
                );
            }
            // Whole-process closure, independently of individual-unit checks.
            let source = [stream(&fs, &r, "A"), stream(&fs, &r, "B")];
            let products = [
                stream(&fs, &r, "separator_vapor"),
                stream(&fs, &r, "separator_liquid"),
            ];
            let energy_in: f64 = source.iter().map(|s| s.enthalpy_flow_w()).sum();
            let energy_out: f64 = products.iter().map(|s| s.enthalpy_flow_w()).sum();
            assert!((energy_out - energy_in - r.units[&cooler].duty_w).abs() / scale < 2e-4);
            let mass_in: f64 = source.iter().map(|s| s.mass_flow_kg_s()).sum();
            let mass_out: f64 = products.iter().map(|s| s.mass_flow_kg_s()).sum();
            assert!((mass_out - mass_in).abs() / scale < 1e-9);
            for i in 0..support::COMPONENTS.len() {
                let n_in: f64 = source.iter().map(|s| s.component_flows_mol_s(5)[i]).sum();
                let n_out: f64 = products.iter().map(|s| s.component_flows_mol_s(5)[i]).sum();
                assert!((n_out - n_in).abs() / scale < 1e-7);
            }
        }
    }
}

#[test]
fn preflight_failures_and_input_changes_cannot_leave_current_results() {
    let (mut fs, pkg, mut feeds) =
        support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let valve = unit(&fs, UnitKind::Valve);
    fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(6e6)),
        },
    );
    assert!(!r.is_current(&fs, &feeds, &pkg));
    let error = fs.solve(&pkg, &feeds).unwrap_err();
    assert_eq!(error.unit, Some(valve));
    assert!(error.message.contains("V-100") && error.message.contains("pressure"));
    let feed = *feeds.keys().next().unwrap();
    feeds.remove(&feed);
    assert!(
        fs.solve(&pkg, &feeds)
            .unwrap_err()
            .message
            .contains("composition")
    );
    let (fs, pkg, mut feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let feed = *feeds.keys().next().unwrap();
    feeds.insert(feed, support::composition(&[0.7, 0.15, 0.1, 0.04, 0.01]));
    assert!(!r.is_current(&fs, &feeds, &pkg));
    let entries = feeds[&feed].entries().iter().rev().cloned().collect();
    feeds.insert(
        feed,
        Composition::new(
            CompositionBasis::Mole,
            entries,
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap(),
    );
    assert!(
        fs.solve(&pkg, &feeds)
            .unwrap_err()
            .message
            .contains("identities/order")
    );
    let (_, other, _) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    assert!(!r.is_current(&fs, &feeds, &other));
}

#[test]
fn liquid_separator_and_zero_feed_network_behave_explicitly() {
    let (fs, pkg, mut feeds, u) = single(
        UnitKind::Separator,
        UnitSpec::Separator {
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
        240.0,
        10.0,
    );
    *feeds.values_mut().next().unwrap() = support::composition(&[0.0, 0.0, 1.0, 0.0, 0.0]);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let vapor = &r.streams[&fs.stream_at(PortRef::new(u, 1)).unwrap().id];
    let liquid = &r.streams[&fs.stream_at(PortRef::new(u, 2)).unwrap().id];
    assert!(vapor.state.is_none());
    assert_eq!(vapor.molar_flow_mol_s, 0.0);
    assert_eq!(liquid.molar_flow_mol_s, 10.0);
    assert_eq!(
        liquid.state.as_ref().unwrap().phases[0].phase,
        thermo_types::Phase::Liquid
    );
    let (mut fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    for &id in feeds.keys() {
        fs.set_spec(
            id,
            UnitSpec::Feed {
                temperature_k: Some(330.0),
                pressure_pa: Some(5e6),
                molar_flow_mol_s: Some(0.0),
            },
        );
    }
    assert!(
        fs.solve(&pkg, &feeds)
            .unwrap_err()
            .message
            .contains("positive total inlet flow")
    );
}

#[test]
fn currency_checks_distinguish_forks_with_the_same_revision() {
    let (mut fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let mut fork = fs.clone();
    let valve = unit(&fs, UnitKind::Valve);
    fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(1.6e6)),
        },
    );
    fork.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(1.7e6)),
        },
    );
    assert_eq!(fs.revision(), fork.revision());
    let r = fs.solve(&pkg, &feeds).unwrap();
    assert!(r.is_current(&fs, &feeds, &pkg.clone()));
    assert!(!r.is_current(&fork, &feeds, &pkg));
    let mut extra = feeds.clone();
    extra.insert(valve, feeds.values().next().unwrap().clone());
    assert!(
        fs.solve(&pkg, &extra)
            .unwrap_err()
            .message
            .contains("non-feed")
    );
    let mut mass = feeds.clone();
    *mass.values_mut().next().unwrap() = Composition::new(
        CompositionBasis::Mass,
        feeds.values().next().unwrap().entries().to_vec(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    assert!(
        fs.solve(&pkg, &mass)
            .unwrap_err()
            .message
            .contains("mole basis")
    );
}

#[test]
fn invalid_package_bounds_and_nonfinite_extensive_values_are_rejected() {
    let (fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let bad = property_packages::CubicPackage::new(
        pkg.model().clone(),
        thermo_solvers::PhFlashOptions {
            minimum_temperature_k: f64::NAN,
            ..Default::default()
        },
    );
    assert!(
        fs.solve(&bad, &feeds)
            .unwrap_err()
            .message
            .contains("finite temperature bounds")
    );
    let (fs, pkg, feeds, _) = single(
        UnitKind::Splitter,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.5),
        },
        330.0,
        f64::MAX,
    );
    assert!(
        fs.solve(&pkg, &feeds)
            .unwrap_err()
            .message
            .contains("nonfinite extensive")
    );
}

#[test]
fn progress_reports_every_unit_in_order_and_cancellation_returns_nothing() {
    let (fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
    let mut seen = Vec::new();
    let r = fs
        .solve_with_progress(&pkg, &feeds, |p| {
            assert_eq!((p.completed, p.total), (seen.len(), 8));
            seen.push(p.unit);
            std::ops::ControlFlow::Continue(())
        })
        .unwrap();
    assert_eq!(seen, r.order);
    assert!(r.is_current(&fs, &feeds, &pkg));

    for stop_at in [0, 3, 7] {
        let error = fs
            .solve_with_progress(&pkg, &feeds, |p| {
                if p.completed == stop_at {
                    std::ops::ControlFlow::Break(())
                } else {
                    std::ops::ControlFlow::Continue(())
                }
            })
            .unwrap_err();
        assert!(error.cancelled && error.unit.is_none());
        assert!(error.message.contains("cancelled"));
    }
    // Preflight failures are reported as failures, not cancellations, and
    // happen before the first progress call.
    let mut empty = Flowsheet::new();
    empty.add_unit(UnitKind::Valve);
    let error = empty
        .solve_with_progress(&pkg, &feeds, |_| panic!("no unit may start"))
        .unwrap_err();
    assert!(!error.cancelled);
}

/// Replaces the product on `port` of `u` with `units` in series, ending in a
/// new product. Returns the inserted unit identifiers.
fn chain(
    fs: &mut Flowsheet,
    u: UnitId,
    port: usize,
    units: &[(UnitKind, UnitSpec)],
) -> Vec<UnitId> {
    let old = fs.stream_at(PortRef::new(u, port)).unwrap().to.unit;
    fs.remove_unit(old);
    let mut from = PortRef::new(u, port);
    let mut ids = vec![];
    for &(kind, spec) in units {
        let id = fs.add_unit(kind);
        assert!(fs.set_spec(id, spec));
        fs.connect(from, PortRef::new(id, 0)).unwrap();
        from = PortRef::new(id, 1);
        ids.push(id);
    }
    let product = fs.add_unit(UnitKind::Product);
    fs.connect(from, PortRef::new(product, 0)).unwrap();
    ids
}

fn valve_and_heater() -> [(UnitKind, UnitSpec); 2] {
    [
        (
            UnitKind::Valve,
            UnitSpec::Valve {
                pressure: PressureSpec::PressureDrop(Some(1e5)),
            },
        ),
        (
            UnitKind::Heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(5e3)),
                pressure_drop_pa: 0.0,
            },
        ),
    ]
}

fn assert_empty_pass_through(fs: &Flowsheet, r: &Calculation, id: UnitId) {
    for (i, p) in fs.unit(id).unwrap().kind.ports().iter().enumerate() {
        if p.direction == PortDirection::Outlet {
            let s = &r.streams[&fs.stream_at(PortRef::new(id, i)).unwrap().id];
            assert_eq!(s.molar_flow_mol_s, 0.0);
            assert!(s.state.is_none());
        }
    }
    let u = &r.units[&id];
    assert_eq!(u.duty_w, 0.0);
    assert_eq!(u.warnings.len(), 1, "{:?}", u.warnings);
    assert!(u.warnings[0].contains("Zero"), "{}", u.warnings[0]);
    let b = u.balance.as_ref().unwrap();
    assert_eq!(b.max_scaled_component_residual, 0.0);
    assert_eq!(b.energy_residual_w, 0.0);
}

#[test]
fn zero_flow_streams_pass_through_downstream_units() {
    // All-vapour feed: the separator liquid line is empty; the valve, heater
    // and a second separator downstream report it instead of failing (#89).
    let (mut fs, pkg, feeds, sep) = single(
        UnitKind::Separator,
        UnitSpec::Separator {
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
        330.0,
        10.0,
    );
    let mut line = valve_and_heater().to_vec();
    line.push((UnitKind::Separator, UnitKind::Separator.default_spec()));
    let ids = chain(&mut fs, sep, 2, &line);
    // `chain` connected the second separator's vapour port; add its liquid.
    let liquid_product = fs.add_unit(UnitKind::Product);
    fs.connect(PortRef::new(ids[2], 2), PortRef::new(liquid_product, 0))
        .unwrap();
    let r = fs.solve(&pkg, &feeds).unwrap();
    let liquid = &r.streams[&fs.stream_at(PortRef::new(sep, 2)).unwrap().id];
    assert!(liquid.state.is_none() && liquid.molar_flow_mol_s == 0.0);
    for id in ids {
        assert_empty_pass_through(&fs, &r, id);
    }

    // Splitter endpoints: the empty branch passes through a valve and heater,
    // the full branch is calculated normally.
    for fraction in [0.0, 1.0] {
        let (mut fs, pkg, feeds, split) = single(
            UnitKind::Splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(fraction),
            },
            330.0,
            10.0,
        );
        let first = chain(&mut fs, split, 1, &valve_and_heater());
        let second = chain(&mut fs, split, 2, &valve_and_heater());
        let r = fs.solve(&pkg, &feeds).unwrap();
        let (empty, full) = if fraction == 0.0 {
            (first, second)
        } else {
            (second, first)
        };
        for id in empty {
            assert_empty_pass_through(&fs, &r, id);
        }
        let heater = full[1];
        assert!(r.units[&heater].warnings.is_empty());
        assert_eq!(r.units[&heater].duty_w, 5e3);
        let out = &r.streams[&fs.stream_at(PortRef::new(heater, 1)).unwrap().id];
        assert_eq!(out.molar_flow_mol_s, 10.0);
    }

    // An exchanger with an empty side transfers no heat; the flowing side
    // keeps its enthalpy and takes its own pressure drop.
    for fraction in [0.0, 1.0] {
        let (mut fs, pkg, feeds, split) = single(
            UnitKind::Splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(fraction),
            },
            330.0,
            10.0,
        );
        let ends: Vec<_> = fs
            .units()
            .filter(|u| u.kind == UnitKind::Product)
            .map(|u| u.id)
            .collect();
        for id in ends {
            fs.remove_unit(id);
        }
        let ex = fs.add_unit(UnitKind::HeatExchanger);
        fs.set_spec(
            ex,
            UnitSpec::HeatExchanger {
                thermal: ExchangerSpec::Duty(Some(1e3)),
                hot_pressure_drop_pa: 1e5,
                cold_pressure_drop_pa: 2e5,
            },
        );
        for i in 0..2 {
            fs.connect(PortRef::new(split, i + 1), PortRef::new(ex, i))
                .unwrap();
            let product = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(ex, i + 2), PortRef::new(product, 0))
                .unwrap();
        }
        let r = fs.solve(&pkg, &feeds).unwrap();
        let u = &r.units[&ex];
        assert_eq!((u.duty_w, u.transferred_duty_w), (0.0, Some(0.0)));
        assert_eq!(u.warnings.len(), 1);
        // fraction 0: the hot side (first outlet) is empty.
        let (empty, full, drop) = if fraction == 0.0 {
            (2, 3, 2e5)
        } else {
            (3, 2, 1e5)
        };
        let empty = &r.streams[&fs.stream_at(PortRef::new(ex, empty)).unwrap().id];
        assert!(empty.state.is_none() && empty.molar_flow_mol_s == 0.0);
        let full = &r.streams[&fs.stream_at(PortRef::new(ex, full)).unwrap().id];
        let inlet = r.streams[&fs.stream_at(PortRef::new(split, 0)).unwrap().id]
            .state
            .clone()
            .unwrap();
        let state = full.state.as_ref().unwrap();
        assert_eq!(full.molar_flow_mol_s, 10.0);
        assert_eq!(state.pressure_pa, 5e6 - drop);
        assert!((state.enthalpy() - inlet.enthalpy()).abs() <= 1e-9 * inlet.enthalpy().abs());
        assert!(u.balance.as_ref().unwrap().scaled_energy_residual.abs() < 1e-9);
    }
}

#[test]
fn heater_and_cooler_direction_contradictions_are_reported() {
    // Contradictory duty signs are rejected before any property calculation.
    for (kind, q) in [(UnitKind::Heater, -1e3), (UnitKind::Cooler, 1e3)] {
        let (fs, pkg, feeds, u) = single(
            kind,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(q)),
                pressure_drop_pa: 0.0,
            },
            330.0,
            10.0,
        );
        assert_eq!(fs.unit_status(u), Some(UnitStatus::Underspecified));
        let e = fs
            .solve_with_progress(&pkg, &feeds, |_| panic!("no unit may start"))
            .unwrap_err();
        assert!(e.diagnostics.iter().any(|d| d.kind
            == DiagnosticKind::InvalidSpecification(
                u,
                SpecIssue {
                    field: "duty",
                    problem: SpecProblem::ContradictsUnitKind,
                }
            )));
        assert!(
            e.message.contains("contradicts the unit kind"),
            "{}",
            e.message
        );
    }
    // An outlet temperature on the wrong side of the inlet temperature needs
    // the inlet state; the unit is solved and carries a warning.
    for (kind, t, warning) in [
        (UnitKind::Heater, 300.0, Some("Heater removes heat")),
        (UnitKind::Cooler, 350.0, Some("Cooler adds heat")),
        (UnitKind::Heater, 350.0, None),
        (UnitKind::Cooler, 300.0, None),
    ] {
        let (fs, pkg, feeds, u) = single(
            kind,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(t)),
                pressure_drop_pa: 0.0,
            },
            330.0,
            10.0,
        );
        let r = fs.solve(&pkg, &feeds).unwrap();
        let warnings = &r.units[&u].warnings;
        match warning {
            Some(w) => {
                assert_eq!(warnings.len(), 1);
                assert!(warnings[0].contains(w), "{}", warnings[0]);
            }
            None => assert!(warnings.is_empty(), "{warnings:?}"),
        }
    }
}

#[test]
fn renaming_or_revision_changes_alone_keep_results_current() {
    let (mut fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let r = fs.solve(&pkg, &feeds).unwrap();
    let valve = unit(&fs, UnitKind::Valve);
    let stream = fs.streams().next().unwrap().id;
    let revision = fs.revision();
    assert!(fs.rename_unit(valve, "PV-1"));
    assert!(fs.rename_stream(stream, "renamed"));
    assert!(fs.revision() > revision);
    assert!(r.is_current(&fs, &feeds, &pkg));
    fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(1.4e6)),
        },
    );
    assert!(!r.is_current(&fs, &feeds, &pkg));
}

#[test]
fn energy_balance_is_scaled_by_the_larger_of_enthalpy_flow_and_n_cp_t() {
    let (fs, pkg, feeds) = support::demonstration(CubicVariant::PengRobinson, 1.0);
    let r = fs.solve(&pkg, &feeds).unwrap();
    for u in fs.units() {
        let Some(b) = &r.units[&u.id].balance else {
            continue;
        };
        let mut enthalpy = r.units[&u.id].duty_w.abs();
        let mut thermal = 0.0;
        for i in 0..u.kind.ports().len() {
            let s = &r.streams[&fs.stream_at(PortRef::new(u.id, i)).unwrap().id];
            enthalpy += s.enthalpy_flow_w().abs();
            thermal += s.thermal_scale_w();
            if let Some(state) = &s.state {
                let cp: f64 = state
                    .phases
                    .iter()
                    .map(|p| p.phase_fraction * p.state.cp)
                    .sum();
                let expected = s.molar_flow_mol_s * cp * state.temperature_k;
                assert!((s.thermal_scale_w() - expected).abs() <= 1e-12 * expected.abs());
                assert!(expected > 0.0 || s.molar_flow_mol_s == 0.0);
            }
        }
        let scale = enthalpy.max(thermal).max(1.0);
        let expected = b.energy_residual_w / scale;
        assert!((b.scaled_energy_residual - expected).abs() <= 1e-12 * expected.abs());
    }
}
