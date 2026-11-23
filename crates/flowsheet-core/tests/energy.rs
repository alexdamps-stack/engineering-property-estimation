//! Energy streams (#171), optional and variable ports, and drawn-only kinds
//! (A17 foundation).
#[path = "../examples/support/mod.rs"]
mod support;

use flowsheet_core::units::{ColumnSpec, columns::DistillationColumnSpec};
use flowsheet_core::*;
use thermo_models::CubicVariant;

/// Energy port of heaters and coolers (after `in` and `out`).
const Q: usize = 2;

fn feed_spec(t: f64, n: f64) -> UnitSpec {
    UnitSpec::Feed {
        temperature_k: Some(t),
        pressure_pa: Some(5e6),
        molar_flow_mol_s: Some(n),
    }
}

fn heater_spec(thermal: ThermalSpec) -> UnitSpec {
    UnitSpec::HeatExchange {
        thermal,
        pressure_drop_pa: 0.0,
    }
}

/// Feed → heater → product, with the given heater specification. Returns
/// the flowsheet, package, feeds and the heater.
fn heater_line(
    thermal: ThermalSpec,
) -> (
    Flowsheet,
    property_packages::CubicPackage,
    FeedCompositions,
    UnitId,
) {
    let (_, pkg, feeds) = support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
    let mut fs = Flowsheet::new();
    let feed = fs.add_unit(UnitKind::Feed);
    let heater = fs.add_unit(UnitKind::Heater);
    let product = fs.add_unit(UnitKind::Product);
    assert!(fs.set_spec(feed, feed_spec(300.0, 10.0)));
    assert!(fs.set_spec(heater, heater_spec(thermal)));
    fs.connect(PortRef::new(feed, 0), PortRef::new(heater, 0))
        .unwrap();
    fs.connect(PortRef::new(heater, 1), PortRef::new(product, 0))
        .unwrap();
    let composition = feeds.values().next().unwrap().clone();
    (fs, pkg, [(feed, composition)].into(), heater)
}

fn energy_source(fs: &mut Flowsheet, duty: EnergySourceDuty, to: UnitId) -> (UnitId, StreamId) {
    let source = fs.add_unit(UnitKind::EnergySource);
    assert!(fs.set_spec(source, UnitSpec::EnergySource { duty }));
    let s = fs
        .connect(PortRef::new(source, 0), PortRef::new(to, Q))
        .unwrap();
    (source, s)
}

fn outlet_temperature(fs: &Flowsheet, r: &Calculation, unit: UnitId) -> f64 {
    let s = fs.stream_at(PortRef::new(unit, 1)).unwrap().id;
    r.streams[&s].state.as_ref().unwrap().temperature_k
}

#[test]
fn heater_and_energy_ports_are_optional_and_typed() {
    let (mut fs, _, _, heater) = heater_line(ThermalSpec::OutletTemperature(Some(320.0)));
    assert!(fs.diagnostics().is_empty(), "{:?}", fs.diagnostics());
    let ports = fs.unit(heater).unwrap().ports();
    assert_eq!(ports.len(), 3);
    assert_eq!(ports[Q].medium, PortMedium::Energy);
    assert!(ports[Q].optional);
    // Energy ports accept only energy streams.
    let feed = fs.add_unit(UnitKind::Feed);
    assert_eq!(
        fs.connect(PortRef::new(feed, 0), PortRef::new(heater, Q)),
        Err(ConnectError::MediumMismatch(
            PortRef::new(feed, 0),
            PortRef::new(heater, Q)
        ))
    );
    let sink = fs.add_unit(UnitKind::EnergySink);
    let cooler = fs.add_unit(UnitKind::Cooler);
    let id = fs
        .connect(PortRef::new(cooler, Q), PortRef::new(sink, 0))
        .unwrap();
    assert_eq!(fs.stream(id).unwrap().name, format!("Q{}", id.0));
    assert_eq!(fs.stream_medium(id), Some(PortMedium::Energy));
    assert_eq!(fs.energy_streams().count(), 1);
    assert_eq!(fs.material_streams().count(), 2);
}

#[test]
fn specified_energy_source_drives_a_heater_like_a_specified_duty() {
    let q = 5.0e3;
    let (fs, pkg, feeds, heater) = heater_line(ThermalSpec::Duty(Some(q)));
    let direct = fs.solve(&pkg, &feeds).unwrap();

    let (mut fs, pkg, feeds, heater2) = heater_line(ThermalSpec::EnergyStream);
    let kinds: Vec<_> = fs.diagnostics().into_iter().map(|d| d.kind).collect();
    assert_eq!(
        kinds,
        [DiagnosticKind::EnergyStreamMissing(PortRef::new(
            heater2, Q
        ))]
    );
    let (_, s) = energy_source(&mut fs, EnergySourceDuty::Specified(Some(q)), heater2);
    assert!(fs.diagnostics().is_empty(), "{:?}", fs.diagnostics());
    let r = fs.solve(&pkg, &feeds).unwrap();
    assert_eq!(r.energy_streams[&s].duty_w, q);
    assert_eq!(r.units[&heater2].duty_w, q);
    assert!(r.units[&heater2].balance.is_some());
    let (a, b) = (
        outlet_temperature(&fs, &r, heater2),
        outlet_temperature(&fs, &direct, heater),
    );
    assert!((a - b).abs() < 1e-9, "{a} vs {b}");
}

#[test]
fn calculated_energy_source_reports_the_duty_of_its_unit() {
    let (mut fs, pkg, feeds, heater) = heater_line(ThermalSpec::OutletTemperature(Some(330.0)));
    let (source, s) = energy_source(&mut fs, EnergySourceDuty::Calculated, heater);
    assert!(fs.diagnostics().is_empty(), "{:?}", fs.diagnostics());
    let r = fs.solve(&pkg, &feeds).unwrap();
    let q = r.units[&heater].duty_w;
    assert!(q > 0.0);
    assert_eq!(r.energy_streams[&s].duty_w, q);
    // Boundaries carry no balance.
    assert!(r.units[&source].balance.is_none());
}

#[test]
fn cooler_heat_is_transferred_to_a_heater() {
    // The demonstration cooler (240 K outlet) heats a cold side stream.
    let (mut fs, pkg, mut feeds) =
        support::demonstration(CubicVariant::PengRobinsonBostonMathias, 1.0);
    let cooler = fs.units().find(|u| u.kind == UnitKind::Cooler).unwrap().id;
    let cold = fs.add_unit(UnitKind::Feed);
    let heater = fs.add_unit(UnitKind::Heater);
    let product = fs.add_unit(UnitKind::Product);
    fs.set_spec(cold, feed_spec(260.0, 400.0));
    fs.set_spec(heater, heater_spec(ThermalSpec::EnergyStream));
    fs.connect(PortRef::new(cold, 0), PortRef::new(heater, 0))
        .unwrap();
    fs.connect(PortRef::new(heater, 1), PortRef::new(product, 0))
        .unwrap();
    let q = fs
        .connect(PortRef::new(cooler, Q), PortRef::new(heater, Q))
        .unwrap();
    feeds.insert(cold, support::composition(&[0.88, 0.06, 0.03, 0.02, 0.01]));
    assert!(fs.diagnostics().is_empty(), "{:?}", fs.diagnostics());
    let r = fs.solve(&pkg, &feeds).unwrap();
    let removed = r.units[&cooler].duty_w;
    assert!(removed < 0.0);
    // Positive in the stream direction: out of the cooler, into the heater.
    assert_eq!(r.energy_streams[&q].duty_w, -removed);
    assert_eq!(r.units[&heater].duty_w, -removed);
    assert!(outlet_temperature(&fs, &r, heater) > 260.0);
    for u in [cooler, heater] {
        let b = r.units[&u].balance.as_ref().unwrap();
        assert!(b.scaled_energy_residual.abs() <= 2e-9);
    }
    // The heater depends on the cooler through the energy stream only.
    let pos = |id| r.order.iter().position(|x| *x == id).unwrap();
    assert!(pos(cooler) < pos(heater));
}

#[test]
fn energy_streams_must_be_determined_by_exactly_one_end() {
    // Specified source and a heater with its own outlet temperature.
    let (mut fs, pkg, feeds, heater) = heater_line(ThermalSpec::OutletTemperature(Some(330.0)));
    let (source, s) = energy_source(&mut fs, EnergySourceDuty::Specified(Some(1e3)), heater);
    let kinds: Vec<_> = fs.diagnostics().into_iter().map(|d| d.kind).collect();
    assert_eq!(kinds, [DiagnosticKind::EnergyOverspecified(s)]);
    assert!(fs.solve(&pkg, &feeds).is_err());
    // Calculated source and a heater taking its duty from the stream.
    fs.set_spec(source, UnitKind::EnergySource.default_spec());
    fs.set_spec(heater, heater_spec(ThermalSpec::EnergyStream));
    let kinds: Vec<_> = fs.diagnostics().into_iter().map(|d| d.kind).collect();
    assert_eq!(kinds, [DiagnosticKind::EnergyUnderspecified(s)]);
    // A cooler cannot take its duty from its outgoing energy stream.
    let cooler = fs.add_unit(UnitKind::Cooler);
    fs.set_spec(cooler, heater_spec(ThermalSpec::EnergyStream));
    assert!(fs.diagnostics().iter().any(|d| d.kind
        == DiagnosticKind::InvalidSpecification(
            cooler,
            SpecIssue {
                field: "duty from energy stream",
                problem: SpecProblem::NotAvailableForKind,
            }
        )));
}

#[test]
fn a_zero_flow_heater_rejects_a_delivered_duty() {
    let (mut fs, pkg, feeds, heater) = heater_line(ThermalSpec::EnergyStream);
    energy_source(&mut fs, EnergySourceDuty::Specified(Some(1e3)), heater);
    let feed = fs.units().find(|u| u.kind == UnitKind::Feed).unwrap().id;
    fs.set_spec(feed, feed_spec(300.0, 0.0));
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert_eq!(e.unit, Some(heater));
    assert!(e.message.contains("delivers"), "{}", e.message);
}

#[test]
fn energy_streams_are_never_torn() {
    // Heater → cooler with the cooler's heat returned to the heater: the
    // only cycle runs through one material and one energy stream.
    let (mut fs, pkg, feeds, heater) = heater_line(ThermalSpec::EnergyStream);
    let product = fs.units().find(|u| u.kind == UnitKind::Product).unwrap().id;
    let link = fs.stream_at(PortRef::new(heater, 1)).unwrap().id;
    fs.disconnect(link);
    let cooler = fs.add_unit(UnitKind::Cooler);
    fs.set_spec(
        cooler,
        heater_spec(ThermalSpec::OutletTemperature(Some(290.0))),
    );
    let material = fs
        .connect(PortRef::new(heater, 1), PortRef::new(cooler, 0))
        .unwrap();
    fs.connect(PortRef::new(cooler, 1), PortRef::new(product, 0))
        .unwrap();
    let energy = fs
        .connect(PortRef::new(cooler, Q), PortRef::new(heater, Q))
        .unwrap();
    assert_eq!(fs.automatic_tears(), [material].into());
    let guess = MaterialStream {
        molar_flow_mol_s: 1.0,
        state: None,
    };
    let e = fs
        .solve_recycles(
            &pkg,
            &feeds,
            &[(energy, guess)].into(),
            RecycleOptions::default(),
            |_| std::ops::ControlFlow::Continue(()),
        )
        .unwrap_err();
    assert!(e.message.contains("energy stream"), "{}", e.message);
}

#[test]
fn drawn_only_kinds_stop_in_preflight_with_their_issue() {
    // Any drawn-only kind with one inlet and one outlet whose specification
    // is complete by default; the test passes vacuously once none is left.
    let Some(kind) = UnitKind::ALL.into_iter().find(|k| {
        !k.calculation_planned()
            && k.issue().is_some()
            && k.ports().len() == 2
            && k.ports()[0].direction == PortDirection::Inlet
            && k.default_spec().issues().is_empty()
    }) else {
        return;
    };
    let (mut fs, pkg, feeds, heater) = heater_line(ThermalSpec::OutletTemperature(Some(320.0)));
    let product = fs.stream_at(PortRef::new(heater, 1)).unwrap().to.unit;
    let link = fs.stream_at(PortRef::new(heater, 1)).unwrap().id;
    fs.disconnect(link);
    let unit = fs.add_unit(kind);
    fs.connect(PortRef::new(heater, 1), PortRef::new(unit, 0))
        .unwrap();
    fs.connect(PortRef::new(unit, 1), PortRef::new(product, 0))
        .unwrap();
    // Energy ports are optional.
    assert!(
        fs.diagnostics().is_empty(),
        "{kind}: {:?}",
        fs.diagnostics()
    );
    let e = fs.solve(&pkg, &feeds).unwrap_err();
    assert_eq!(e.unit, Some(unit));
    assert!(e.message.contains("drawn only"), "{}", e.message);
    let issue = format!("(#{})", kind.issue().unwrap());
    assert!(e.message.ends_with(&issue), "{}", e.message);
}

#[test]
fn variable_ports_keep_their_streams_when_the_group_changes() {
    let mut fs = Flowsheet::new();
    let column = fs.add_unit(UnitKind::DistillationColumn);
    let feed = fs.add_unit(UnitKind::Feed);
    let draw = fs.add_unit(UnitKind::Product);
    // distillate, bottoms, condenser duty, reboiler duty, feed 1
    let ports = fs.unit(column).unwrap().ports();
    assert_eq!(ports[4].label(), "feed 1");
    let f = fs
        .connect(PortRef::new(feed, 0), PortRef::new(column, 4))
        .unwrap();
    #[allow(clippy::needless_update)]
    let spec = |feeds, side_draws| {
        UnitSpec::Columns(ColumnSpec::Distillation(DistillationColumnSpec {
            feeds,
            side_draws,
            // Further fields the columns group adds keep their defaults.
            ..Default::default()
        }))
    };
    fs.try_set_spec(column, spec(2, 1)).unwrap();
    let ports = fs.unit(column).unwrap().ports();
    assert_eq!(ports[6].label(), "side draw 1");
    let d = fs
        .connect(PortRef::new(column, 6), PortRef::new(draw, 0))
        .unwrap();
    // A third feed moves the side draw from index 6 to 7.
    fs.try_set_spec(column, spec(3, 1)).unwrap();
    assert_eq!(fs.stream(d).unwrap().from, PortRef::new(column, 7));
    assert_eq!(fs.stream(f).unwrap().to, PortRef::new(column, 4));
    // Removing the connected side draw is rejected and changes nothing.
    let revision = fs.revision();
    assert_eq!(
        fs.try_set_spec(column, spec(3, 0)),
        Err(SpecChangeError::ConnectedPortRemoved(PortRef::new(
            column, 7
        )))
    );
    assert_eq!(fs.revision(), revision);
    assert_eq!(
        fs.try_set_spec(column, UnitKind::Cstr.default_spec()),
        Err(SpecChangeError::WrongKind(column))
    );
    // Restoring re-validates connections against the variable ports.
    let restored = Flowsheet::restore(fs.units().cloned(), fs.streams().cloned()).expect("restore");
    assert!(restored.same_calculation_inputs(&fs));
}
