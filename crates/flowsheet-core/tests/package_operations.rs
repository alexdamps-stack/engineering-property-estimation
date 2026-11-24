//! Units that used cubic-only operations (#226) solved under the non-cubic
//! flowsheet packages: rotating machines through the package PS flash, a
//! shortcut column through package bubble/dew points, and reactors rejected
//! explicitly where the package has no formation basis. Component and energy
//! balances must close (relative 1e-8 and 2e-9 scaled energy).

use flowsheet_core::*;
use property_packages::{FlowsheetPackage, GergPackage, NrtlPackage};
use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_models::{ActivityModel, GergModel, VaporModel};
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

/// Energy ports follow the two material ports.
const POWER: usize = 2;

fn composition(ids: &[&str], z: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        ids.iter()
            .zip(z)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

const NG: [&str; 5] = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"];
const NG_Z: [f64; 5] = [0.88, 0.06, 0.03, 0.02, 0.01];

fn gerg(ids: &[&str], z: &[f64]) -> GergPackage {
    GergPackage::new(
        GergModel::new(&composition(ids, z)).unwrap(),
        GergPackage::TEMPERATURE_RANGE_K,
    )
    .unwrap()
}

fn nrtl(z: &[f64]) -> NrtlPackage {
    let model = ActivityModel::from_databases(
        &ComponentDatabase::embedded().unwrap(),
        &LiquidDatabase::embedded().unwrap(),
        &composition(&["methanol", "water"], z),
        VaporModel::IdealGas,
    )
    .unwrap();
    NrtlPackage::new(model, NrtlPackage::TEMPERATURE_RANGE_K).unwrap()
}

/// Feed → unit → product; a rotating unit gets a calculated power boundary.
fn line(
    ids: &[&str],
    z: &[f64],
    inlet: (f64, f64, f64),
    kind: UnitKind,
    spec: UnitSpec,
) -> (Flowsheet, FeedCompositions, UnitId, StreamId) {
    let mut flowsheet = Flowsheet::new();
    let source = flowsheet.add_unit(UnitKind::Feed);
    let unit = flowsheet.add_unit(kind);
    let product = flowsheet.add_unit(UnitKind::Product);
    assert!(flowsheet.set_spec(
        source,
        UnitSpec::Feed {
            temperature_k: Some(inlet.0),
            pressure_pa: Some(inlet.1),
            molar_flow_mol_s: Some(inlet.2),
        }
    ));
    flowsheet.try_set_spec(unit, spec).unwrap();
    flowsheet
        .connect(PortRef::new(source, 0), PortRef::new(unit, 0))
        .unwrap();
    let outlet = flowsheet
        .connect(PortRef::new(unit, 1), PortRef::new(product, 0))
        .unwrap();
    match kind {
        UnitKind::Pump | UnitKind::Compressor => {
            let e = flowsheet.add_unit(UnitKind::EnergySource);
            assert!(flowsheet.set_spec(
                e,
                UnitSpec::EnergySource {
                    duty: EnergySourceDuty::Calculated
                }
            ));
            flowsheet
                .connect(PortRef::new(e, 0), PortRef::new(unit, POWER))
                .unwrap();
        }
        UnitKind::Expander => {
            let e = flowsheet.add_unit(UnitKind::EnergySink);
            flowsheet
                .connect(PortRef::new(unit, POWER), PortRef::new(e, 0))
                .unwrap();
        }
        _ => {}
    }
    let mut feeds = FeedCompositions::new();
    feeds.insert(source, composition(ids, z));
    (flowsheet, feeds, unit, outlet)
}

fn closed(calc: &Calculation, unit: UnitId) {
    let balance = calc.units[&unit].balance.as_ref().unwrap();
    assert!(balance.max_scaled_component_residual <= 1e-8);
    assert!(balance.scaled_energy_residual <= 2e-9);
}

#[test]
fn gerg_compressor_and_expander_use_the_gerg_ps_flash() {
    let package = gerg(&NG, &NG_Z);
    for (kind, spec, inlet, p_out) in [
        (
            UnitKind::Compressor,
            UnitSpec::Pressurize {
                pressure: PressureRise::OutletPressure(Some(4e6)),
                efficiency: 0.75,
            },
            (260.0, 2e6, 10.0),
            4e6,
        ),
        (
            UnitKind::Expander,
            UnitSpec::Expand {
                pressure: PressureSpec::OutletPressure(Some(3e6)),
                efficiency: 0.8,
            },
            (340.0, 6e6, 10.0),
            3e6,
        ),
    ] {
        let (fs, feeds, unit, outlet) = line(&NG, &NG_Z, inlet, kind, spec);
        let calc = fs
            .solve(&package, &feeds)
            .unwrap_or_else(|e| panic!("{kind}: {e:?}"));
        closed(&calc, unit);
        let result = calc.units[&unit].rotating.as_ref().unwrap();
        let stage = &result.stages[0];
        // The isentropic state reproduces the inlet entropy on the GERG reference.
        let inlet_state = package.tp(&NG_Z, inlet.0, inlet.1).unwrap();
        let ideal = package
            .tp(&NG_Z, stage.isentropic_outlet_temperature_k, p_out)
            .unwrap();
        assert!(
            (ideal.entropy() - inlet_state.entropy()).abs()
                <= 1e-8 * inlet_state.entropy().abs().max(1.0),
            "{kind}"
        );
        let out = calc.streams[&outlet].state.as_ref().unwrap();
        // Irreversibility: the outlet entropy exceeds the inlet entropy.
        assert!(out.entropy() > inlet_state.entropy(), "{kind}");
        if kind == UnitKind::Compressor {
            assert!(stage.power_w > stage.isentropic_power_w && stage.isentropic_power_w > 0.0);
        } else {
            assert!(stage.power_w > stage.isentropic_power_w && stage.power_w < 0.0);
        }
    }
}

#[test]
fn nrtl_pump_reports_the_pressure_independent_liquid_explicitly() {
    let z = [0.5, 0.5];
    let package = nrtl(&z);
    let (fs, feeds, unit, outlet) = line(
        &["methanol", "water"],
        &z,
        (320.0, 101_325.0, 5.0),
        UnitKind::Pump,
        UnitSpec::Pressurize {
            pressure: PressureRise::OutletPressure(Some(500_000.0)),
            efficiency: 0.7,
        },
    );
    let calc = fs.solve(&package, &feeds).unwrap();
    closed(&calc, unit);
    let report = &calc.units[&unit];
    assert!(
        report
            .warnings
            .iter()
            .any(|w| w.contains("no Poynting term")),
        "{:?}",
        report.warnings
    );
    let out = calc.streams[&outlet].state.as_ref().unwrap();
    assert!((out.temperature_k - 320.0).abs() < 1e-6);
    assert_eq!(out.pressure_pa, 500_000.0);
}

#[test]
fn nrtl_expander_of_vapor_partially_condenses_on_the_gamma_phi_isentrope() {
    let z = [0.5, 0.5];
    let package = nrtl(&z);
    let (fs, feeds, unit, outlet) = line(
        &["methanol", "water"],
        &z,
        (420.0, 300_000.0, 2.0),
        UnitKind::Expander,
        UnitSpec::Expand {
            pressure: PressureSpec::OutletPressure(Some(101_325.0)),
            efficiency: 0.8,
        },
    );
    let calc = fs.solve(&package, &feeds).unwrap();
    closed(&calc, unit);
    let stage = &calc.units[&unit].rotating.as_ref().unwrap().stages[0];
    let inlet = package.tp(&z, 420.0, 300_000.0).unwrap();
    let ideal = package.ps(&z, 101_325.0, inlet.entropy()).unwrap();
    assert_eq!(stage.isentropic_outlet_temperature_k, ideal.temperature_k);
    assert!(stage.power_w < 0.0);
    assert!(calc.streams[&outlet].state.is_some());
}

/// Ethane/propane shortcut column at `p` under GERG-2008, fed at `feed_t`.
fn gerg_shortcut(
    ids: &[&str],
    z: &[f64],
    p: f64,
    feed_t: f64,
) -> (Result<Calculation, SolveError>, UnitId) {
    let package = gerg(ids, z);
    let mut flowsheet = Flowsheet::new();
    let source = flowsheet.add_unit(UnitKind::Feed);
    let column = flowsheet.add_unit(UnitKind::ShortcutColumn);
    let top = flowsheet.add_unit(UnitKind::Product);
    let bottom = flowsheet.add_unit(UnitKind::Product);
    assert!(flowsheet.set_spec(
        source,
        UnitSpec::Feed {
            temperature_k: Some(feed_t),
            pressure_pa: Some(p),
            molar_flow_mol_s: Some(100.0),
        }
    ));
    flowsheet
        .try_set_spec(
            column,
            UnitSpec::Columns(units::ColumnSpec::Shortcut(
                units::columns::ShortcutColumnSpec {
                    light_key: Some("ethane".into()),
                    heavy_key: Some("propane".into()),
                    light_key_recovery: Some(0.95),
                    heavy_key_recovery: Some(0.95),
                    pressure_pa: Some(p),
                    reflux: units::columns::ShortcutReflux::Ratio(Some(1.1)),
                },
            )),
        )
        .unwrap();
    flowsheet
        .connect(PortRef::new(source, 0), PortRef::new(column, 0))
        .unwrap();
    flowsheet
        .connect(PortRef::new(column, 1), PortRef::new(top, 0))
        .unwrap();
    flowsheet
        .connect(PortRef::new(column, 2), PortRef::new(bottom, 0))
        .unwrap();
    let mut feeds = FeedCompositions::new();
    feeds.insert(source, composition(ids, z));
    (flowsheet.solve(&package, &feeds), column)
}

#[test]
fn shortcut_column_under_gerg_uses_gerg_bubble_and_dew_points() {
    let (calc, column) = gerg_shortcut(&["ethane", "propane"], &[0.5, 0.5], 1.5e6, 280.0);
    let calc = calc.unwrap_or_else(|e| panic!("{e:?}"));
    closed(&calc, column);
    let result = calc.units[&column].shortcut_column.as_ref().unwrap();
    // Model-to-model plausibility only (no GERG column reference): ethane/propane
    // volatility at 1.5 MPa lies between 1 and 4 and Fenske stages are finite.
    assert!(result.relative_volatility > 1.0 && result.relative_volatility < 4.0);
    assert!(result.minimum_stages.is_finite() && result.minimum_stages > 1.0);
    // The ternary of validation/columns/shortcut-v1.json has a methane-rich
    // distillate whose bubble point lies below the 230 K GERG bound: the
    // column fails at that operation instead of extrapolating.
    let (failed, column) = gerg_shortcut(
        &["methane", "ethane", "propane"],
        &[0.2, 0.45, 0.35],
        2.5e6,
        260.0,
    );
    let e = failed.unwrap_err();
    assert_eq!(e.unit, Some(column));
    assert!(
        e.message.contains("distillate bubble point"),
        "{}",
        e.message
    );
}

#[test]
fn reactors_fail_explicitly_without_a_formation_basis() {
    let package = gerg(&NG, &NG_Z);
    let (fs, feeds, unit, _) = line(
        &NG,
        &NG_Z,
        (300.0, 2e6, 1.0),
        UnitKind::ConversionReactor,
        UnitSpec::ConversionReactor {
            reactions: vec![Reaction {
                coefficients: vec![("ethane".into(), -1.0), ("methane".into(), 2.0)],
                extent: ReactionExtent::MolarExtent(Some(0.0)),
            }],
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
    );
    let e = fs.solve(&package, &feeds).unwrap_err();
    assert_eq!(e.unit, Some(unit));
    assert!(e.message.contains("formation basis"), "{}", e.message);
    let dynamic: &dyn FlowsheetPackage = &package;
    assert!(dynamic.formation_components().is_err());
}
