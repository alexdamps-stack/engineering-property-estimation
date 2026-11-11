#[path = "../examples/support/mod.rs"]
mod support;

use flowsheet_core::*;
use property_packages::CubicPackage;
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use thermo_models::CubicVariant;

fn network(
    variant: CubicVariant,
    fraction: f64,
    scale: f64,
) -> (Flowsheet, CubicPackage, FeedCompositions, StreamId) {
    let (_, pkg, feeds) = support::demonstration(variant, scale);
    let mut fs = Flowsheet::new();
    let feed = fs.add_unit(UnitKind::Feed);
    let mix = fs.add_unit(UnitKind::Mixer);
    let heat = fs.add_unit(UnitKind::Heater);
    let split = fs.add_unit(UnitKind::Splitter);
    let product = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(100.0 * scale),
        },
    );
    fs.set_spec(
        heat,
        UnitSpec::HeatExchange {
            thermal: ThermalSpec::OutletTemperature(Some(350.0)),
            pressure_drop_pa: 0.0,
        },
    );
    fs.set_spec(
        split,
        UnitSpec::Splitter {
            fraction_to_first: Some(fraction),
        },
    );
    for (from, fp, to, tp) in [
        (feed, 0, mix, 0),
        (mix, 2, heat, 0),
        (heat, 1, split, 0),
        (split, 2, product, 0),
    ] {
        fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
            .unwrap();
    }
    let tear = fs
        .connect(PortRef::new(split, 1), PortRef::new(mix, 1))
        .unwrap();
    (
        fs,
        pkg,
        [(feed, feeds.values().next().unwrap().clone())].into(),
        tear,
    )
}

fn start(pkg: &CubicPackage, flow: f64, temperature: f64) -> MaterialStream {
    MaterialStream {
        molar_flow_mol_s: flow,
        state: Some(
            pkg.tp(&[0.88, 0.06, 0.03, 0.02, 0.01], temperature, 5e6)
                .unwrap(),
        ),
    }
}

#[test]
fn recycle_matches_analytic_flow_component_and_energy_solution_reproducibly() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        for scale in [0.01, 1.0, 100.0] {
            let (fs, pkg, feeds, tear) = network(variant, 0.5, scale);
            let run = || {
                fs.solve_network_with_progress(&pkg, &feeds, |_| ControlFlow::Continue(()))
                    .unwrap()
            };
            let report = run();
            let repeat = run();
            let closure = report.recycle.as_ref().unwrap();
            assert_eq!(closure.iterations, repeat.recycle.unwrap().iterations);
            assert!(closure.max_scaled_residual <= 1e-10);
            assert_eq!(fs.automatic_tears(), [tear].into());
            assert_eq!(
                report.streams[&tear].molar_flow_mol_s,
                repeat.streams[&tear].molar_flow_mol_s
            );
            assert!(report.is_current(&fs, &feeds, &pkg));
            let n = 100.0 * scale;
            assert!((report.streams[&tear].molar_flow_mol_s / n - 1.0).abs() < 2e-10);
            let product = fs.units().find(|u| u.kind == UnitKind::Product).unwrap();
            let stream = &report.streams[&fs.stream_at(PortRef::new(product.id, 0)).unwrap().id];
            assert!((stream.molar_flow_mol_s / n - 1.0).abs() < 2e-10);
            let inlet = start(&pkg, n, 330.0);
            let outlet = start(&pkg, n, 350.0);
            let expected_duty = outlet.enthalpy_flow_w() - inlet.enthalpy_flow_w();
            let duty: f64 = report.units.values().map(|u| u.duty_w).sum();
            assert!((duty / expected_duty - 1.0).abs() < 2e-8);
            for (actual, expected) in stream
                .component_flows_mol_s(5)
                .iter()
                .zip(inlet.component_flows_mol_s(5))
            {
                assert!((actual - expected).abs() <= n * 2e-10);
            }
            for u in report.units.values() {
                if let Some(b) = &u.balance {
                    assert!(b.max_scaled_component_residual <= 1e-8);
                    assert!(b.scaled_energy_residual.abs() <= 2e-9);
                }
            }
        }
    }
}

#[test]
fn alternate_initial_flow_temperature_and_composition_reach_the_same_fixed_point() {
    let (fs, pkg, feeds, tear) = network(CubicVariant::PengRobinsonBostonMathias, 0.5, 1.0);
    for flow in [0.0, 5.0, 500.0] {
        let mut guess = start(&pkg, flow, 300.0);
        guess.state = Some(pkg.tp(&[0.70, 0.15, 0.13, 0.01, 0.01], 300.0, 5e6).unwrap());
        let report = fs
            .solve_recycles(
                &pkg,
                &feeds,
                &[(tear, guess)].into(),
                RecycleOptions::default(),
                |_| ControlFlow::Continue(()),
            )
            .unwrap();
        assert!((report.streams[&tear].molar_flow_mol_s - 100.0).abs() < 2e-8);
        for (x, y) in report.streams[&tear]
            .state
            .as_ref()
            .unwrap()
            .mole_fractions
            .iter()
            .zip([0.88, 0.06, 0.03, 0.02, 0.01])
        {
            assert!((x - y).abs() < 2e-10);
        }
    }
}

#[test]
fn cancellation_on_a_later_pass_and_iteration_exhaustion_return_no_report() {
    let (fs, pkg, feeds, tear) = network(CubicVariant::PengRobinson, 0.5, 1.0);
    let guesses = [(tear, start(&pkg, 0.0, 330.0))].into();
    let failure = fs
        .solve_recycles(&pkg, &feeds, &guesses, RecycleOptions::default(), |p| {
            if p.iteration == 3 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .unwrap_err();
    assert!(failure.cancelled);
    let failure = fs
        .solve_recycles(
            &pkg,
            &feeds,
            &guesses,
            RecycleOptions {
                maximum_iterations: 1,
                ..Default::default()
            },
            |_| ControlFlow::Continue(()),
        )
        .unwrap_err();
    assert!(!failure.cancelled);
    assert!(failure.message.contains("did not converge after 1"));
    assert!(failure.message.contains("residual"));
    // There is no finite solution with positive fresh feed and 100% recycle.
    let (fs, pkg, feeds, tear) = network(CubicVariant::PengRobinson, 1.0, 1.0);
    let failure = fs
        .solve_recycles(
            &pkg,
            &feeds,
            &[(tear, start(&pkg, 0.0, 330.0))].into(),
            RecycleOptions {
                maximum_iterations: 8,
                ..Default::default()
            },
            |_| ControlFlow::Continue(()),
        )
        .unwrap_err();
    assert!(failure.message.contains("did not converge after 8"));
}

#[test]
fn incomplete_tears_bad_guesses_invalid_options_and_missing_inputs_fail() {
    let (mut fs, pkg, feeds, tear) = network(CubicVariant::PengRobinson, 0.5, 1.0);
    let run = |guesses: &BTreeMap<_, _>, options| {
        fs.solve_recycles(
            &pkg,
            &feeds,
            guesses,
            options,
            |_| ControlFlow::Continue(()),
        )
    };
    assert!(
        run(&BTreeMap::new(), RecycleOptions::default())
            .unwrap_err()
            .message
            .contains("break every")
    );
    assert!(
        run(
            &[(StreamId(999), start(&pkg, 0.0, 330.0))].into(),
            RecycleOptions::default()
        )
        .unwrap_err()
        .message
        .contains("Unknown tear")
    );
    for flow in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            run(
                &[(tear, start(&pkg, flow, 330.0))].into(),
                RecycleOptions::default()
            )
            .is_err()
        );
    }
    for tolerance in [0.0, -1.0, f64::NAN, 1.0] {
        assert!(
            run(
                &[(tear, start(&pkg, 0.0, 330.0))].into(),
                RecycleOptions {
                    tolerance,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }
    assert!(
        run(
            &[(
                tear,
                MaterialStream {
                    molar_flow_mol_s: 1.0,
                    state: None
                }
            )]
            .into(),
            RecycleOptions::default()
        )
        .is_err()
    );
    fs.disconnect(tear);
    let mut called = false;
    assert!(
        fs.solve_network_with_progress(&pkg, &feeds, |_| {
            called = true;
            ControlFlow::Continue(())
        })
        .is_err()
    );
    assert!(!called);
}

#[test]
fn two_recycles_and_a_large_downstream_chain_converge() {
    let (mut fs, pkg, mut feeds, _) = network(CubicVariant::PengRobinson, 0.5, 1.0);
    let a = fs.add_unit(UnitKind::Feed);
    let m = fs.add_unit(UnitKind::Mixer);
    let s = fs.add_unit(UnitKind::Splitter);
    let p = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        a,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(20.0),
        },
    );
    fs.set_spec(
        s,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.25),
        },
    );
    feeds.insert(a, feeds.values().next().unwrap().clone());
    for (from, fp, to, tp) in [(a, 0, m, 0), (m, 2, s, 0), (s, 1, m, 1)] {
        fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
            .unwrap();
    }
    let mut from = s;
    let mut port = 2;
    for _ in 0..100 {
        let next = fs.add_unit(UnitKind::Heater);
        fs.set_spec(
            next,
            UnitSpec::HeatExchange {
                thermal: ThermalSpec::OutletTemperature(Some(330.0)),
                pressure_drop_pa: 0.0,
            },
        );
        fs.connect(PortRef::new(from, port), PortRef::new(next, 0))
            .unwrap();
        from = next;
        port = 1;
    }
    fs.connect(PortRef::new(from, port), PortRef::new(p, 0))
        .unwrap();
    assert_eq!(fs.automatic_tears().len(), 2);
    let report = fs
        .solve_network_with_progress(&pkg, &feeds, |_| ControlFlow::Continue(()))
        .unwrap();
    assert_eq!(report.units.len(), 109);
    assert_eq!(report.recycle.unwrap().tear_residuals.len(), 2);
    let product = &report.streams[&fs.stream_at(PortRef::new(p, 0)).unwrap().id];
    assert!((product.molar_flow_mol_s - 20.0).abs() < 1e-8);
}

#[test]
fn interacting_loops_share_units_and_match_the_analytic_solution() {
    let (mut fs, pkg, feeds, _) = network(CubicVariant::PengRobinsonBostonMathias, 0.5, 1.0);
    let mixer = fs.units().find(|u| u.kind == UnitKind::Mixer).unwrap().id;
    let heater = fs.units().find(|u| u.kind == UnitKind::Heater).unwrap().id;
    let splitter = fs
        .units()
        .find(|u| u.kind == UnitKind::Splitter)
        .unwrap()
        .id;
    let product = fs.units().find(|u| u.kind == UnitKind::Product).unwrap().id;
    let a = fs.stream_at(PortRef::new(mixer, 2)).unwrap().id;
    let b = fs.stream_at(PortRef::new(product, 0)).unwrap().id;
    fs.disconnect(a);
    fs.disconnect(b);
    let second_mixer = fs.add_unit(UnitKind::Mixer);
    let second_splitter = fs.add_unit(UnitKind::Splitter);
    fs.set_spec(
        second_splitter,
        UnitSpec::Splitter {
            fraction_to_first: Some(0.25),
        },
    );
    for (from, fp, to, tp) in [
        (mixer, 2, second_mixer, 0),
        (second_mixer, 2, heater, 0),
        (splitter, 2, second_splitter, 0),
        (second_splitter, 1, second_mixer, 1),
        (second_splitter, 2, product, 0),
    ] {
        fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
            .unwrap();
    }
    assert_eq!(fs.automatic_tears().len(), 2);
    let report = fs.solve_network(&pkg, &feeds).unwrap();
    let n = report.streams[&fs.stream_at(PortRef::new(heater, 1)).unwrap().id].molar_flow_mol_s;
    assert!((n / (100.0 / 0.375) - 1.0).abs() < 1e-9);
    assert!(
        (report.streams[&fs.stream_at(PortRef::new(product, 0)).unwrap().id].molar_flow_mol_s
            - 100.0)
            .abs()
            < 1e-7
    );
}

#[test]
fn automatic_tears_start_with_reordered_units_and_zero_recycle() {
    let (fs, pkg, feeds, _) = network(CubicVariant::PengRobinson, 0.0, 1.0);
    let reorder = |id: UnitId| UnitId(20 - id.0);
    let units = fs.units().cloned().map(|mut u| {
        u.id = reorder(u.id);
        u
    });
    let streams = fs.streams().cloned().map(|mut s| {
        s.from.unit = reorder(s.from.unit);
        s.to.unit = reorder(s.to.unit);
        s
    });
    let mut fs = Flowsheet::restore(units, streams).unwrap();
    let feeds = feeds.into_iter().map(|(id, c)| (reorder(id), c)).collect();
    for fraction in [0.0, 0.5, 0.9] {
        let splitter = fs
            .units()
            .find(|u| u.kind == UnitKind::Splitter)
            .unwrap()
            .id;
        fs.set_spec(
            splitter,
            UnitSpec::Splitter {
                fraction_to_first: Some(fraction),
            },
        );
        let report = fs.solve_network(&pkg, &feeds).unwrap();
        let tear = *fs.automatic_tears().first().unwrap();
        assert_eq!(
            fs.unit(fs.stream(tear).unwrap().to.unit).unwrap().kind,
            UnitKind::Mixer
        );
        assert!(
            (report.streams[&tear].molar_flow_mol_s - 100.0 * fraction / (1.0 - fraction)).abs()
                < 2e-6
        );
    }
}

#[test]
fn small_absolute_tear_changes_cannot_bypass_connected_component_balances() {
    let (fs, pkg, feeds, tear) = network(CubicVariant::PengRobinson, 0.5, 1e-16);
    let guesses = [(tear, start(&pkg, 0.0, 350.0))].into();
    // At this scale the tear differences are below the absolute SI floors
    // almost immediately, while the relative component imbalance is large.
    let error = fs
        .solve_recycles(
            &pkg,
            &feeds,
            &guesses,
            RecycleOptions {
                maximum_iterations: 2,
                ..Default::default()
            },
            |_| ControlFlow::Continue(()),
        )
        .unwrap_err();
    assert!(error.message.contains("Final connected-stream balance"));
    let report = fs
        .solve_recycles(&pkg, &feeds, &guesses, RecycleOptions::default(), |_| {
            ControlFlow::Continue(())
        })
        .unwrap();
    assert!(report.recycle.as_ref().unwrap().iterations >= 20);
    assert!((report.streams[&tear].molar_flow_mol_s / 1e-14 - 1.0).abs() < 3e-8);
}

#[test]
fn feed_effluent_exchanger_loop_without_a_mixer_converges_for_either_tear() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        for reverse_connections in [false, true] {
            let (_, pkg, initial_feeds) = support::demonstration(variant, 1.0);
            let mut fs = Flowsheet::new();
            // The lowest-ID feed belongs to a different process branch.
            let unrelated = fs.add_unit(UnitKind::Feed);
            let unrelated_product = fs.add_unit(UnitKind::Product);
            fs.set_spec(
                unrelated,
                UnitSpec::Feed {
                    temperature_k: Some(250.0),
                    pressure_pa: Some(5e6),
                    molar_flow_mol_s: Some(10.0),
                },
            );
            fs.connect(
                PortRef::new(unrelated, 0),
                PortRef::new(unrelated_product, 0),
            )
            .unwrap();
            let feed = fs.add_unit(UnitKind::Feed);
            let hx = fs.add_unit(UnitKind::HeatExchanger);
            let heater = fs.add_unit(UnitKind::Heater);
            let product = fs.add_unit(UnitKind::Product);
            fs.set_spec(
                feed,
                UnitSpec::Feed {
                    temperature_k: Some(330.0),
                    pressure_pa: Some(5e6),
                    molar_flow_mol_s: Some(10.0),
                },
            );
            fs.set_spec(
                hx,
                UnitSpec::HeatExchanger {
                    thermal: ExchangerSpec::Duty(Some(1000.0)),
                    hot_pressure_drop_pa: 0.0,
                    cold_pressure_drop_pa: 0.0,
                },
            );
            fs.set_spec(
                heater,
                UnitSpec::HeatExchange {
                    thermal: ThermalSpec::OutletTemperature(Some(350.0)),
                    pressure_drop_pa: 0.0,
                },
            );
            let mut edges = vec![
                (feed, 0, hx, 1),
                (hx, 3, heater, 0),
                (heater, 1, hx, 0),
                (hx, 2, product, 0),
            ];
            if reverse_connections {
                edges.reverse();
            }
            for (from, fp, to, tp) in edges {
                fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
                    .unwrap();
            }
            let feeds = [
                (feed, initial_feeds.values().next().unwrap().clone()),
                (
                    unrelated,
                    support::composition(&[0.70, 0.15, 0.13, 0.01, 0.01]),
                ),
            ]
            .into();
            let report = fs.solve_network(&pkg, &feeds).unwrap();
            assert!(report.recycle.as_ref().unwrap().max_scaled_residual <= 1e-10);
            for s in report.streams.values() {
                assert!((s.molar_flow_mol_s - 10.0).abs() < 1e-9);
            }
            for u in report.units.values() {
                if let Some(b) = &u.balance {
                    assert!(b.max_scaled_component_residual <= 1e-8);
                    assert!(b.scaled_energy_residual.abs() <= 2e-9);
                }
            }
            let outlet = &report.streams[&fs.stream_at(PortRef::new(product, 0)).unwrap().id];
            let inlet = pkg.tp(&[0.88, 0.06, 0.03, 0.02, 0.01], 330.0, 5e6).unwrap();
            assert!(
                (10.0 * (outlet.state.as_ref().unwrap().enthalpy() - inlet.enthalpy())
                    - report.units[&heater].duty_w)
                    .abs()
                    < 1e-5
            );
        }
    }
}

#[test]
fn exchanger_sides_do_not_create_a_material_feed_for_a_closed_loop() {
    for variant in [
        CubicVariant::PengRobinson,
        CubicVariant::PengRobinsonBostonMathias,
    ] {
        for reverse_connections in [false, true] {
            for loop_side in [0, 1] {
                let (_, pkg, initial_feeds) = support::demonstration(variant, 1.0);
                let mut fs = Flowsheet::new();
                let feed = fs.add_unit(UnitKind::Feed);
                let hx = fs.add_unit(UnitKind::HeatExchanger);
                let heater = fs.add_unit(UnitKind::Heater);
                let product = fs.add_unit(UnitKind::Product);
                fs.set_spec(
                    feed,
                    UnitSpec::Feed {
                        temperature_k: Some(if loop_side == 0 { 330.0 } else { 360.0 }),
                        pressure_pa: Some(5e6),
                        molar_flow_mol_s: Some(10.0),
                    },
                );
                fs.set_spec(
                    hx,
                    UnitSpec::HeatExchanger {
                        thermal: ExchangerSpec::Duty(Some(0.0)),
                        hot_pressure_drop_pa: 0.0,
                        cold_pressure_drop_pa: 0.0,
                    },
                );
                fs.set_spec(
                    heater,
                    UnitSpec::HeatExchange {
                        thermal: ThermalSpec::OutletTemperature(Some(350.0)),
                        pressure_drop_pa: 0.0,
                    },
                );
                let flowing_side = 1 - loop_side;
                let mut edges = vec![
                    (feed, 0, hx, flowing_side),
                    (hx, flowing_side + 2, product, 0),
                    (hx, loop_side + 2, heater, 0),
                    (heater, 1, hx, loop_side),
                ];
                if reverse_connections {
                    edges.reverse();
                }
                for (from, fp, to, tp) in edges {
                    fs.connect(PortRef::new(from, fp), PortRef::new(to, tp))
                        .unwrap();
                }
                let feeds = [(feed, initial_feeds.values().next().unwrap().clone())].into();
                let error = fs.solve_network(&pkg, &feeds).unwrap_err();
                assert!(
                    error.message.contains("materially reachable feed"),
                    "{}",
                    error.message
                );
                assert!(error.message.contains("explicit tear states"));
                // A closed material loop can be solved only when the caller
                // supplies its independent circulation flow and starting state.
                let tear = *fs.automatic_tears().first().unwrap();
                let explicit = [(tear, start(&pkg, 7.0, 350.0))].into();
                let report = fs
                    .solve_recycles(&pkg, &feeds, &explicit, RecycleOptions::default(), |_| {
                        ControlFlow::Continue(())
                    })
                    .unwrap();
                let loop_stream = fs.stream_at(PortRef::new(hx, loop_side + 2)).unwrap().id;
                assert!((report.streams[&loop_stream].molar_flow_mol_s - 7.0).abs() < 1e-10);
            }
        }
    }
}
