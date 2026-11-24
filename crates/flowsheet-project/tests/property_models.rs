use std::ops::ControlFlow;

use flowsheet_core::{PortRef, PressureSpec, UnitKind, UnitSpec};
use flowsheet_project::{PackageChoice, PackageSelection, Project, solve_project};

#[test]
fn eos_cg_2021_accepts_all_components_in_one_saved_project() {
    let components = PackageSelection::available_components_for(PackageChoice::EosCg2021).unwrap();
    assert_eq!(components.len(), 16);
    let hydrogen = components.iter().position(|id| id == "hydrogen").unwrap();
    let project = Project::new(PackageSelection::default_for(
        PackageChoice::EosCg2021,
        components,
        "",
    ));
    let project = Project::from_json(&project.to_json().unwrap()).unwrap();
    let prepared = project.package.prepare().unwrap();
    let mut fractions = vec![1e-5; 16];
    fractions[hydrogen] = 1.0 - 15e-5;
    let state = prepared.package().tp(&fractions, 345.0, 120_000.0).unwrap();
    assert_eq!(state.phases[0].state.ln_fugacity_coefficients.len(), 16);
    assert!(
        state.phases[0]
            .state
            .speed_of_sound
            .is_some_and(|speed| speed > 0.0)
    );
}

#[test]
fn eos_cg_2021_flowsheet_recovers_reference_gas_density_and_solves_a_valve() {
    let mut project = Project::new(PackageSelection::default_for(
        PackageChoice::EosCg2021,
        vec!["hydrogen".into(), "carbon-dioxide".into()],
        "",
    ));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    let valve = fs.add_unit(UnitKind::Valve);
    let product = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(350.0),
            // x_CO2=0.9 and rho=250 mol/m³ in the independent teqp fixture.
            pressure_pa: Some(715_094.847_477_150_2),
            molar_flow_mol_s: Some(10.0),
        },
    );
    fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(500_000.0)),
        },
    );
    project.feeds.insert(feed, vec![0.1, 0.9]);
    let inlet = fs
        .connect(PortRef::new(feed, 0), PortRef::new(valve, 0))
        .unwrap();
    let outlet = fs
        .connect(PortRef::new(valve, 1), PortRef::new(product, 0))
        .unwrap();
    let project = Project::from_json(&project.to_json().unwrap()).unwrap();
    let package = project.package.prepare().unwrap();
    let result = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    assert_eq!(result.provenance.model, "EOS-CG-2021");
    let inlet = result.streams[&inlet].state.as_ref().unwrap();
    let outlet = result.streams[&outlet].state.as_ref().unwrap();
    let implied_density = inlet.pressure_pa
        / (thermo_types::MOLAR_GAS_CONSTANT
            * inlet.temperature_k
            * inlet.phases[0].compressibility);
    assert!((implied_density - 250.0).abs() < 1e-7);
    assert!((inlet.enthalpy_j_mol - outlet.enthalpy_j_mol).abs() < 1e-5);
    for (fraction, expected) in outlet.phases[0].mole_fractions.iter().zip([0.1, 0.9]) {
        assert!((fraction - expected).abs() < 1e-14);
    }
    assert!(
        outlet.phases[0]
            .speed_of_sound_m_s
            .is_some_and(|speed| speed > 0.0)
    );
    assert_eq!(outlet.phases[0].dynamic_viscosity_pa_s, None);
    assert!(package.package().tp(&[0.1, 0.9], 300.0, 500_000.0).is_err());
}

#[test]
fn gerg_valve_reproduces_the_independent_reference_after_project_reload() {
    let mut package = PackageSelection::default_for(
        PackageChoice::Gerg2008,
        ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"]
            .map(String::from)
            .to_vec(),
        "",
    );
    package.ph_minimum_temperature_k = 300.0;
    let mut project = Project::new(package);
    let flowsheet = &mut project.flowsheet;
    let feed = flowsheet.add_unit(UnitKind::Feed);
    let valve = flowsheet.add_unit(UnitKind::Valve);
    let product = flowsheet.add_unit(UnitKind::Product);
    flowsheet.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(10.0),
        },
    );
    flowsheet.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(3e6)),
        },
    );
    project
        .feeds
        .insert(feed, vec![0.88, 0.06, 0.03, 0.02, 0.01]);
    flowsheet
        .connect(PortRef::new(feed, 0), PortRef::new(valve, 0))
        .unwrap();
    let outlet = flowsheet
        .connect(PortRef::new(valve, 1), PortRef::new(product, 0))
        .unwrap();

    let project = Project::from_json(&project.to_json().unwrap()).unwrap();
    let package = project.package.prepare().unwrap();
    let result = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    assert_eq!(result.provenance.model, "GERG-2008");
    assert_eq!(
        result.provenance.caloric_reference.as_deref(),
        Some("gerg-native-ideal")
    );
    let outlet = result.streams[&outlet].state.as_ref().unwrap();
    // ng-isenthalpic-valve, validation/fixtures/gerg-equilibrium-v1.json.
    assert!((outlet.temperature_k - 322.30329325632096).abs() < 2e-6);
    assert!((outlet.enthalpy_j_mol - 424.5677333317572).abs() < 2e-5);
    assert!(outlet.phases[0].speed_of_sound_m_s.is_some_and(|v| v > 0.0));
}

#[test]
fn nrtl_separator_keeps_inferred_regression_provenance_and_the_reference_split() {
    let mut project = Project::new(PackageSelection::default_for(
        PackageChoice::NrtlIdealGas,
        vec!["methanol".into(), "water".into()],
        "",
    ));
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    let separator = fs.add_unit(UnitKind::Separator);
    let vapor = fs.add_unit(UnitKind::Product);
    let liquid = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(355.0),
            pressure_pa: Some(101_325.0),
            molar_flow_mol_s: Some(10.0),
        },
    );
    project.feeds.insert(feed, vec![0.3, 0.7]);
    let inlet = fs
        .connect(PortRef::new(feed, 0), PortRef::new(separator, 0))
        .unwrap();
    let top = fs
        .connect(PortRef::new(separator, 1), PortRef::new(vapor, 0))
        .unwrap();
    let bottom = fs
        .connect(PortRef::new(separator, 2), PortRef::new(liquid, 0))
        .unwrap();

    let project = Project::from_json(&project.to_json().unwrap()).unwrap();
    let package = project.package.prepare().unwrap();
    let result = solve_project(&project, &package, |_| ControlFlow::Continue(()), 0).unwrap();
    assert_eq!(result.provenance.model, "NRTL-IG");
    let inlet = &result.streams[&inlet];
    let liquid = &result.streams[&bottom];
    let vapor = &result.streams[&top];
    // The 355 K methanol/water TP case in validation/fixtures/nrtl-v1.json.
    assert!((vapor.molar_flow_mol_s / 10.0 - 0.27684751579190725).abs() < 5e-8);
    assert!((vapor.molar_flow_mol_s + liquid.molar_flow_mol_s - 10.0).abs() < 1e-10);
    assert!(
        inlet
            .warnings
            .iter()
            .any(|w| w.contains("A7-NRTL-INFERRED-REGRESSION-INTERVAL"))
    );
    assert_eq!(
        liquid.state.as_ref().unwrap().phases[0].nrtl_regression_interval_inferred,
        Some(true)
    );
    let liquid_phase = &liquid.state.as_ref().unwrap().phases[0];
    assert!(liquid_phase.cp_j_mol_k.is_some());
    assert_eq!(liquid_phase.cv_j_mol_k, None);
    assert_eq!(liquid_phase.speed_of_sound_m_s, None);
    let vapor_phase = &vapor.state.as_ref().unwrap().phases[0];
    assert!(vapor_phase.cv_j_mol_k.is_some_and(|v| v > 0.0));
    assert!(vapor_phase.speed_of_sound_m_s.is_some_and(|v| v > 0.0));
    assert_eq!(vapor_phase.dynamic_viscosity_pa_s, None);
    assert_eq!(vapor_phase.thermal_conductivity_w_m_k, None);
}
