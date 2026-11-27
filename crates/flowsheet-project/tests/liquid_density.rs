//! #228 engineering liquid density through the project solver: separate
//! reporting beside the EOS density, volumetric liquid flow, pump hydraulics,
//! reactor liquid-volume basis, explicit unavailability, and persistence.
use std::ops::ControlFlow;

use flowsheet_core::{
    Arrhenius, ConcentrationBasis, EnergySourceDuty, KineticReaction, Kinetics, LiquidVolumeBasis,
    PortRef, PressureRise, RateBasis, RateLaw, RateTerm, ReactingPhase, SeparatorThermal, StreamId,
    UnitId, UnitKind, UnitSpec,
};
use flowsheet_project::{CubicVariantChoice, PackageSelection, Project, ResultSet, solve_project};
use property_packages::SaturationKind;
use thermo_models::CostaldLiquidDensity;

fn selection(components: &[&str]) -> PackageSelection {
    PackageSelection::zero_interactions(
        CubicVariantChoice::PengRobinson,
        components.iter().map(|c| c.to_string()).collect(),
        "explicit zero-interaction approximation for the #228 tests",
    )
}

fn solve(project: &Project) -> ResultSet {
    solve_project(
        project,
        &project.package.prepare().unwrap(),
        |_| ControlFlow::Continue(()),
        0,
    )
    .unwrap()
}

struct PumpLine {
    project: Project,
    pump: UnitId,
    inlet: StreamId,
}

/// Feed (T, p, 2 mol/s) → pump to `p_out` → product.
fn pump_line(components: &[&str], z: Vec<f64>, t: f64, p: f64, p_out: f64) -> PumpLine {
    let mut project = Project::new(selection(components));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(p),
            molar_flow_mol_s: Some(2.0),
        },
    );
    project.feeds.insert(feed, z);
    let pump = fs.add_unit(UnitKind::Pump);
    fs.set_spec(
        pump,
        UnitSpec::Pressurize {
            pressure: PressureRise::OutletPressure(Some(p_out)),
            efficiency: 0.75,
        },
    );
    let product = fs.add_unit(UnitKind::Product);
    let source = fs.add_unit(UnitKind::EnergySource);
    fs.set_spec(
        source,
        UnitSpec::EnergySource {
            duty: EnergySourceDuty::Calculated,
        },
    );
    let inlet = fs
        .connect(PortRef::new(feed, 0), PortRef::new(pump, 0))
        .unwrap();
    fs.connect(PortRef::new(pump, 1), PortRef::new(product, 0))
        .unwrap();
    fs.connect(PortRef::new(source, 0), PortRef::new(pump, 2))
        .unwrap();
    PumpLine {
        project,
        pump,
        inlet,
    }
}

#[test]
fn liquid_pump_reports_engineering_density_volume_flow_and_hydraulic_power() {
    let ids = ["ethane", "propane"];
    let (t, p_in, p_out) = (260.0, 3e6, 6e6);
    let mut line = pump_line(&ids, vec![0.3, 0.7], t, p_in, p_out);
    let results = solve(&line.project);

    let inlet = &results.streams[&line.inlet];
    let state = inlet.state.as_ref().unwrap();
    assert_eq!(state.phases.len(), 1);
    let phase = &state.phases[0];
    assert_eq!(phase.label, "liquid");
    let engineering = phase.engineering_liquid_density.as_ref().unwrap();
    let rho = engineering.density_kg_m3.unwrap();
    assert_eq!(
        engineering.validity.as_deref(),
        Some("within declared domain")
    );
    assert!(
        engineering.basis.starts_with("COSTALD"),
        "{}",
        engineering.basis
    );
    // Reported beside, not instead of, the EOS density.
    assert!((rho / phase.density_kg_m3 - 1.0).abs() > 1e-3);
    let direct = CostaldLiquidDensity::embedded(&ids)
        .unwrap()
        .evaluate(&phase.mole_fractions, t, p_in)
        .unwrap();
    assert!((rho / direct.mass_density - 1.0).abs() < 1e-14);
    let flow = inlet.engineering_liquid_volume_flow_m3_s.unwrap();
    assert!((flow / (2.0 * direct.molar_volume) - 1.0).abs() < 1e-14);
    assert!((flow * rho / inlet.mass_flow_kg_s - 1.0).abs() < 1e-12);

    let unit = &results.units[&line.pump];
    let rotating = unit.rotating.as_ref().unwrap();
    let hydraulics = rotating.pump_hydraulics.as_ref().unwrap();
    let power = hydraulics.hydraulic_power_w.unwrap();
    assert!((power / (flow * (p_out - p_in)) - 1.0).abs() < 1e-14);
    let head = hydraulics.differential_head_m.unwrap();
    assert!((head / ((p_out - p_in) / (rho * 9.80665)) - 1.0).abs() < 1e-14);
    // The shaft power still comes from the enthalpy change (EOS) and closes
    // the energy balance; the hydraulic power is a separate labelled figure.
    assert!(rotating.power_w() > 0.0 && power > 0.0);
    assert!((unit.duty_w - rotating.power_w()).abs() <= 1e-9 * unit.duty_w.abs());
    assert!(unit.scaled_energy_residual.unwrap().abs() <= 2e-9);
    // NPSH available on the package's own bubble pressure at the inlet.
    let bubble = line
        .project
        .package
        .prepare()
        .unwrap()
        .package()
        .saturation_pressure(&state.mole_fractions, t, SaturationKind::Bubble)
        .unwrap();
    let npsh = hydraulics.npsh_available_m.unwrap();
    assert!((npsh / ((p_in - bubble.pressure_pa) / (rho * 9.80665)) - 1.0).abs() < 1e-12);
    assert!(npsh > 0.0 && hydraulics.npsh_note.is_none());

    // Saved and reloaded unchanged.
    line.project.results = Some(results);
    let loaded = Project::from_json(&line.project.to_json().unwrap()).unwrap();
    assert_eq!(loaded.results, line.project.results);
}

#[test]
fn water_and_methanol_pumps_use_the_chemsep_liquid_density() {
    // Pure water: ChemSep DIPPR 106, qualified (about 996.5 kg/m³ at 300 K).
    let line = pump_line(&["water"], vec![1.0], 300.0, 2e5, 5e5);
    let results = solve(&line.project);
    let phase = &results.streams[&line.inlet].state.as_ref().unwrap().phases[0];
    let engineering = phase.engineering_liquid_density.as_ref().unwrap();
    let rho = engineering.density_kg_m3.unwrap();
    assert!((rho / 996.5 - 1.0).abs() < 0.005, "{rho}");
    assert_eq!(
        engineering.validity.as_deref(),
        Some("within declared domain")
    );
    assert!(
        engineering.basis.contains("chemsep-8.3"),
        "{}",
        engineering.basis
    );
    let hydraulics = results.units[&line.pump]
        .rotating
        .as_ref()
        .unwrap()
        .pump_hydraulics
        .as_ref()
        .unwrap();
    let power = hydraulics.hydraulic_power_w.unwrap();
    let flow = results.streams[&line.inlet]
        .engineering_liquid_volume_flow_m3_s
        .unwrap();
    assert!((power / (flow * 3e5) - 1.0).abs() < 1e-14);
    // The NPSH needs the package bubble pressure; outside the cubic A13
    // domain it is reported as unavailable, not substituted.
    assert!(hydraulics.npsh_available_m.is_some() != hydraulics.npsh_note.is_some());

    // Methanol/water liquid: ideal mixing, computed but unassessed.
    let line = pump_line(&["methanol", "water"], vec![0.4, 0.6], 300.0, 2e5, 5e5);
    let results = solve(&line.project);
    let phase = &results.streams[&line.inlet].state.as_ref().unwrap().phases[0];
    let engineering = phase.engineering_liquid_density.as_ref().unwrap();
    assert!(engineering.density_kg_m3.is_some());
    assert_eq!(engineering.validity.as_deref(), Some("unassessed"));
}

#[test]
fn pump_hydraulics_are_unavailable_outside_the_density_domain() {
    // Water above 10 MPa: outside the qualified pressure range of the
    // ChemSep/Tait route.
    let line = pump_line(&["water"], vec![1.0], 300.0, 1.2e7, 1.5e7);
    let results = solve(&line.project);
    let phase = &results.streams[&line.inlet].state.as_ref().unwrap().phases[0];
    let engineering = phase.engineering_liquid_density.as_ref().unwrap();
    assert_eq!(engineering.density_kg_m3, None);
    assert!(
        engineering.basis.contains("qualified up to"),
        "{}",
        engineering.basis
    );
    assert_eq!(
        results.streams[&line.inlet].engineering_liquid_volume_flow_m3_s,
        None
    );
    let hydraulics = results.units[&line.pump]
        .rotating
        .as_ref()
        .unwrap()
        .pump_hydraulics
        .as_ref()
        .unwrap();
    assert_eq!(hydraulics.hydraulic_power_w, None);
    assert!(hydraulics.inlet_density.basis.starts_with("unavailable"));
    // The pump itself still solves on the EOS (shaft power, outlet state).
    assert!(results.units[&line.pump].duty_w > 0.0);
}

#[test]
fn vapor_streams_and_non_cubic_packages_have_no_engineering_density() {
    let line = pump_line(&["ethane", "propane"], vec![0.3, 0.7], 260.0, 3e6, 6e6);
    let results = solve(&line.project);
    for stream in results.streams.values() {
        for phase in &stream.state.as_ref().unwrap().phases {
            assert_eq!(
                phase.engineering_liquid_density.is_some(),
                phase.label == "liquid"
            );
        }
    }
    let nrtl = PackageSelection::default_for(
        flowsheet_project::PackageChoice::NrtlIdealGas,
        vec!["methanol".into(), "water".into()],
        "",
    );
    let package = nrtl.prepare().unwrap();
    assert!(matches!(
        package
            .package()
            .engineering_liquid_density(&[0.5, 0.5], 300.0, 1e5),
        Err(thermo_types::CalculationError::UnsupportedCapability(_))
    ));
}

/// Liquid-phase CSTR, 2 ethane → methane + propane, first order in ethane
/// concentration, isothermal at 250 K and 5 MPa.
fn liquid_cstr(basis: LiquidVolumeBasis, t: f64, z: Vec<f64>) -> (Project, UnitId) {
    let mut project = Project::new(selection(&["methane", "ethane", "propane"]));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(t),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(1.0),
        },
    );
    project.feeds.insert(feed, z);
    let reactor = fs.add_unit(UnitKind::Cstr);
    fs.set_spec(
        reactor,
        UnitSpec::Cstr {
            kinetics: Kinetics {
                reactions: vec![KineticReaction {
                    coefficients: vec![
                        ("ethane".into(), -2.0),
                        ("methane".into(), 1.0),
                        ("propane".into(), 1.0),
                    ],
                    phase: ReactingPhase::Liquid,
                    basis: ConcentrationBasis::MolarConcentration,
                    rate_basis: RateBasis::ReactorVolume,
                    rate: RateLaw::PowerLaw {
                        forward: RateTerm {
                            constant: Arrhenius {
                                pre_exponential: Some(1e-4),
                                ..Default::default()
                            },
                            orders: vec![("ethane".into(), 1.0)],
                        },
                        reverse: None,
                    },
                }],
                catalyst_bulk_density_kg_m3: None,
                liquid_volume: basis,
            },
            volume_m3: Some(0.05),
            thermal: SeparatorThermal::OutletTemperature(Some(t)),
            pressure_drop_pa: 0.0,
        },
    );
    let product = fs.add_unit(UnitKind::Product);
    fs.connect(PortRef::new(feed, 0), PortRef::new(reactor, 0))
        .unwrap();
    fs.connect(PortRef::new(reactor, 1), PortRef::new(product, 0))
        .unwrap();
    (project, reactor)
}

#[test]
fn reactor_liquid_volume_basis_is_explicit_and_persisted() {
    let extent = |basis| {
        let (mut project, reactor) = liquid_cstr(basis, 250.0, vec![0.0, 0.4, 0.6]);
        let results = solve(&project);
        let unit = &results.units[&reactor];
        assert!(unit.max_scaled_component_residual.unwrap() <= 1e-8);
        let xi = unit.reaction.as_ref().unwrap().extents_mol_s[0];
        // Outlet liquid (the only stream containing methane): engineering
        // over EOS density.
        let outlet = results
            .streams
            .values()
            .find(|s| s.state.as_ref().unwrap().mole_fractions[0] > 0.0)
            .unwrap();
        let phase = &outlet.state.as_ref().unwrap().phases[0];
        assert_eq!(phase.label, "liquid");
        let ratio = phase
            .engineering_liquid_density
            .as_ref()
            .unwrap()
            .density_kg_m3
            .unwrap()
            / phase.density_kg_m3;
        project.results = Some(results);
        let text = project.to_json().unwrap();
        assert_eq!(
            text.contains("engineering_density"),
            basis == LiquidVolumeBasis::EngineeringDensity
        );
        let loaded = Project::from_json(&text).unwrap();
        assert_eq!(loaded.results, project.results);
        let spec = &loaded.flowsheet.unit(reactor).unwrap().spec;
        assert!(matches!(spec, UnitSpec::Cstr { kinetics, .. } if kinetics.liquid_volume == basis));
        (xi, ratio)
    };
    let (eos, ratio) = extent(LiquidVolumeBasis::EquationOfState);
    let (engineering, _) = extent(LiquidVolumeBasis::EngineeringDensity);
    assert!(eos > 0.0 && engineering > 0.0);
    assert!((ratio - 1.0).abs() > 0.01, "{ratio}");
    // Liquid concentrations scale with the density used, so the extent of
    // the first-order reaction moves in the direction of the density ratio.
    assert_eq!(
        engineering > eos,
        ratio > 1.0,
        "{engineering} vs {eos}, {ratio}"
    );
    assert!((engineering / eos - 1.0).abs() > 0.01);
}

#[test]
fn reactor_engineering_basis_fails_explicitly_outside_the_density_domain() {
    // Compressed liquid ethane at 295 K (T/Tc = 0.966): the EOS basis solves,
    // the engineering basis refuses to extrapolate beyond T/Tcm = 0.95.
    for basis in [
        LiquidVolumeBasis::EquationOfState,
        LiquidVolumeBasis::EngineeringDensity,
    ] {
        let (project, _) = liquid_cstr(basis, 295.0, vec![0.0, 1.0, 0.0]);
        let result = solve_project(
            &project,
            &project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        );
        match basis {
            LiquidVolumeBasis::EquationOfState => assert!(result.is_ok()),
            LiquidVolumeBasis::EngineeringDensity => {
                let e = result.unwrap_err();
                assert!(
                    e.message.contains("engineering density")
                        && e.message.contains("outside validity domain"),
                    "{}",
                    e.message
                );
            }
        }
    }
}
