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
fn absent_separator_outlet_propagates_through_thermal_units_and_splitter() {
    for liquid_feed in [false, true] {
        let (mut fs, pkg, mut feeds, separator) = single(
            UnitKind::Separator,
            UnitSpec::Separator {
                thermal: SeparatorThermal::Adiabatic,
                pressure_drop_pa: 0.0,
            },
            if liquid_feed { 240.0 } else { 330.0 },
            10.0,
        );
        if liquid_feed {
            *feeds.values_mut().next().unwrap() = support::composition(&[0.0, 0.0, 1.0, 0.0, 0.0]);
        }
        let absent_port = if liquid_feed { 1 } else { 2 };
        let end = fs
            .stream_at(PortRef::new(separator, absent_port))
            .unwrap()
            .to
            .unit;
        fs.remove_unit(end);
        let valve = fs.add_unit(UnitKind::Valve);
        fs.set_spec(
            valve,
            UnitSpec::Valve {
                pressure: PressureSpec::OutletPressure(Some(1e6)),
            },
        );
        let heater = fs.add_unit(UnitKind::Heater);
        fs.set_spec(
            heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(350.0)),
                pressure_drop_pa: 0.0,
            },
        );
        let secondary = fs.add_unit(UnitKind::Separator);
        fs.set_spec(
            secondary,
            UnitSpec::Separator {
                thermal: SeparatorThermal::Duty(Some(0.0)),
                pressure_drop_pa: 0.0,
            },
        );
        let splitter = fs.add_unit(UnitKind::Splitter);
        fs.set_spec(
            splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(0.0),
            },
        );
        for (from, fp, to, tp) in [
            (separator, absent_port, valve, 0),
            (valve, 1, heater, 0),
            (heater, 1, secondary, 0),
            (secondary, 1, splitter, 0),
        ] {
            fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
                .unwrap();
        }
        for (from, port) in [(secondary, 2), (splitter, 1), (splitter, 2)] {
            let p = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(from, port), PortRef::new(p, 0))
                .unwrap();
        }
        let report = fs.solve(&pkg, &feeds).unwrap();
        for u in [valve, heater, secondary] {
            assert_eq!(report.units[&u].duty_w, 0.0);
            assert!(!report.units[&u].warnings.is_empty());
            let b = report.units[&u].balance.as_ref().unwrap();
            assert_eq!(b.energy_residual_w, 0.0);
            assert_eq!(b.max_scaled_component_residual, 0.0);
        }
        for s in fs
            .streams()
            .filter(|s| [valve, heater, secondary, splitter].contains(&s.from.unit))
        {
            assert_eq!(report.streams[&s.id].molar_flow_mol_s, 0.0);
            assert!(report.streams[&s.id].state.is_none());
        }
        fs.set_spec(
            heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(1.0)),
                pressure_drop_pa: 0.0,
            },
        );
        assert!(
            fs.solve(&pkg, &feeds)
                .unwrap_err()
                .message
                .contains("nonzero specified duty")
        );
    }
}

#[test]
fn splitter_endpoint_empty_branch_can_feed_a_heater() {
    for fraction in [0.0, 1.0] {
        let (mut fs, pkg, feeds, splitter) = single(
            UnitKind::Splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(fraction),
            },
            330.0,
            10.0,
        );
        let port = if fraction == 0.0 { 1 } else { 2 };
        let product = fs.stream_at(PortRef::new(splitter, port)).unwrap().to.unit;
        fs.remove_unit(product);
        let heater = fs.add_unit(UnitKind::Heater);
        let product = fs.add_unit(UnitKind::Product);
        fs.connect(PortRef::new(splitter, port), PortRef::new(heater, 0))
            .unwrap();
        fs.connect(PortRef::new(heater, 1), PortRef::new(product, 0))
            .unwrap();
        for thermal in [
            ThermalSpec::OutletTemperature(Some(350.0)),
            ThermalSpec::Duty(Some(0.0)),
        ] {
            fs.set_spec(
                heater,
                UnitSpec::HeatExchange {
                    thermal,
                    pressure_drop_pa: 0.0,
                },
            );
            let report = fs.solve(&pkg, &feeds).unwrap();
            let outlet = &report.streams[&fs.stream_at(PortRef::new(heater, 1)).unwrap().id];
            assert_eq!(outlet.molar_flow_mol_s, 0.0);
            assert!(outlet.state.is_none());
            assert_eq!(report.units[&heater].duty_w, 0.0);
            assert!(!report.units[&heater].warnings.is_empty());
        }
    }
}

#[test]
fn exchanger_with_empty_sides_transfers_no_heat_and_preserves_active_enthalpy() {
    for empty_hot in [false, true] {
        for both_empty in [false, true] {
            let (_, pkg, initial_feeds) =
                support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
            let mut fs = Flowsheet::new();
            let hx = fs.add_unit(UnitKind::HeatExchanger);
            fs.set_spec(
                hx,
                UnitSpec::HeatExchanger {
                    thermal: ExchangerSpec::Duty(Some(0.0)),
                    hot_pressure_drop_pa: 1e5,
                    cold_pressure_drop_pa: 2e5,
                },
            );
            let mut feeds = FeedCompositions::new();
            for side in 0..2 {
                let feed = fs.add_unit(UnitKind::Feed);
                let absent = both_empty || (side == 0) == empty_hot;
                fs.set_spec(
                    feed,
                    UnitSpec::Feed {
                        temperature_k: Some(if side == 0 { 350.0 } else { 330.0 }),
                        pressure_pa: Some(5e6),
                        molar_flow_mol_s: Some(if absent { 0.0 } else { 10.0 }),
                    },
                );
                feeds.insert(feed, initial_feeds.values().next().unwrap().clone());
                fs.connect(PortRef::new(feed, 0), PortRef::new(hx, side))
                    .unwrap();
                let product = fs.add_unit(UnitKind::Product);
                fs.connect(PortRef::new(hx, side + 2), PortRef::new(product, 0))
                    .unwrap();
            }
            let report = fs.solve(&pkg, &feeds).unwrap();
            assert!(!report.units[&hx].warnings.is_empty());
            assert_eq!(report.units[&hx].duty_w, 0.0);
            for side in 0..2 {
                let inlet = &report.streams[&fs.stream_at(PortRef::new(hx, side)).unwrap().id];
                let outlet = &report.streams[&fs.stream_at(PortRef::new(hx, side + 2)).unwrap().id];
                assert_eq!(inlet.molar_flow_mol_s, outlet.molar_flow_mol_s);
                if inlet.molar_flow_mol_s == 0.0 {
                    assert!(outlet.state.is_none());
                } else {
                    let a = inlet.state.as_ref().unwrap();
                    let b = outlet.state.as_ref().unwrap();
                    assert!((a.enthalpy() - b.enthalpy()).abs() < 1e-5);
                    assert_eq!(b.pressure_pa, 5e6 - if side == 0 { 1e5 } else { 2e5 });
                }
            }
            // An outlet-temperature specification on absent material is
            // inactive; on the flowing side it cannot demand heat transfer.
            fs.set_spec(
                hx,
                UnitSpec::HeatExchanger {
                    thermal: if empty_hot {
                        ExchangerSpec::HotOutletTemperature(Some(350.0))
                    } else {
                        ExchangerSpec::ColdOutletTemperature(Some(330.0))
                    },
                    hot_pressure_drop_pa: 0.0,
                    cold_pressure_drop_pa: 0.0,
                },
            );
            fs.solve(&pkg, &feeds).unwrap();
            if !both_empty {
                fs.set_spec(
                    hx,
                    UnitSpec::HeatExchanger {
                        thermal: if empty_hot {
                            ExchangerSpec::ColdOutletTemperature(Some(340.0))
                        } else {
                            ExchangerSpec::HotOutletTemperature(Some(340.0))
                        },
                        hot_pressure_drop_pa: 0.0,
                        cold_pressure_drop_pa: 0.0,
                    },
                );
                assert!(
                    fs.solve(&pkg, &feeds)
                        .unwrap_err()
                        .message
                        .contains("both exchanger sides")
                );
            }
            fs.set_spec(
                hx,
                UnitSpec::HeatExchanger {
                    thermal: ExchangerSpec::Duty(Some(1.0)),
                    hot_pressure_drop_pa: 0.0,
                    cold_pressure_drop_pa: 0.0,
                },
            );
            assert!(
                fs.solve(&pkg, &feeds)
                    .unwrap_err()
                    .message
                    .contains("both exchanger sides")
            );
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
        assert_eq!(r.units[&ex].duty_w, 0.0);
        let balance = r.units[&ex].balance.as_ref().unwrap();
        assert!(
            balance.max_scaled_component_residual < 1e-9
                && balance.scaled_energy_residual.abs() < 1e-9
        );
    }
}

#[test]
fn invalid_specs_zero_flow_domain_and_unsupported_units_fail_explicitly() {
    let cases = [
        (
            UnitKind::Heater,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::Duty(Some(1.0)),
                pressure_drop_pa: 0.0,
            },
            330.0,
            0.0,
            "positive inlet flow",
        ),
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
