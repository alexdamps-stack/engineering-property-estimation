//! Petroleum pseudo-components in a saved project (#229): an assay TBP curve
//! is cut into pseudo-components, stored with the package selection, solved
//! through a heater and a flash separator with closed balances, and reloaded.
use std::ops::ControlFlow;

use flowsheet_core::{PortRef, SeparatorThermal, ThermalSpec, UnitKind, UnitSpec};
use flowsheet_project::{PackageChoice, PackageSelection, Project, solve_project};
use thermo_data::{AssayCut, cut_mole_fractions, cut_tbp_curve};

const KIJ_SOURCE: &str = "explicit approximation: kij = 0 for methane and every hydrocarbon pseudo-component pair (Fahim et al. 2010, sec. 3.6.1); no fitted values";

/// A kerosene-range TBP curve (vol%, K) cut into five pseudo-components of
/// constant Watson K = 11.8.
fn kerosene_cuts() -> Vec<AssayCut> {
    let tbp = [(0.0, 430.0), (30.0, 465.0), (70.0, 505.0), (100.0, 540.0)];
    cut_tbp_curve(&tbp, 5, 11.8, "pseudo:kero-").unwrap()
}

fn selection(variant: PackageChoice) -> (PackageSelection, Vec<f64>) {
    let cuts = kerosene_cuts();
    let x_cuts = cut_mole_fractions(&cuts).unwrap();
    let n = cuts.len() + 1;
    let selection = PackageSelection::with_pseudo_components(
        variant,
        vec!["methane".into()],
        cuts.into_iter().map(|c| c.input).collect(),
        vec![vec![0.0; n]; n],
        KIJ_SOURCE,
    );
    let mut x = vec![0.02];
    x.extend(x_cuts.iter().map(|v| v * 0.98));
    (selection, x)
}

#[test]
fn pseudo_component_flash_closes_balances_and_survives_reload() {
    let (package, x) = selection(PackageChoice::PengRobinson);
    let mut project = Project::new(package);
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    let heater = fs.add_unit(UnitKind::Heater);
    let separator = fs.add_unit(UnitKind::Separator);
    let vapor = fs.add_unit(UnitKind::Product);
    let liquid = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(320.0),
            pressure_pa: Some(150_000.0),
            molar_flow_mol_s: Some(5.0),
        },
    );
    fs.set_spec(
        heater,
        UnitSpec::HeatExchange {
            thermal: ThermalSpec::OutletTemperature(Some(490.0)),
            pressure_drop_pa: 0.0,
        },
    );
    fs.set_spec(
        separator,
        UnitSpec::Separator {
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
    );
    project.feeds.insert(feed, x.clone());
    fs.connect(PortRef::new(feed, 0), PortRef::new(heater, 0))
        .unwrap();
    fs.connect(PortRef::new(heater, 1), PortRef::new(separator, 0))
        .unwrap();
    let top = fs
        .connect(PortRef::new(separator, 1), PortRef::new(vapor, 0))
        .unwrap();
    let bottom = fs
        .connect(PortRef::new(separator, 2), PortRef::new(liquid, 0))
        .unwrap();

    let text = project.to_json().unwrap();
    assert!(text.contains("\"pseudo_components\""));
    assert!(text.contains("pseudo:kero-3"));
    let project = Project::from_json(&text).unwrap();
    assert_eq!(project.package.pseudo_components.len(), 5);
    let prepared = project.package.prepare().unwrap();
    assert!(
        prepared.provenance().data_version.contains("+pseudo-"),
        "{}",
        prepared.provenance().data_version
    );
    let result = solve_project(&project, &prepared, |_| ControlFlow::Continue(()), 0).unwrap();
    for id in [heater, separator] {
        let unit = &result.units[&id];
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-8);
        assert!(unit.scaled_energy_residual.unwrap().abs() < 2e-9);
    }
    let top = &result.streams[&top];
    let bottom = &result.streams[&bottom];
    assert!(top.molar_flow_mol_s > 0.0 && bottom.molar_flow_mol_s > 0.0);
    assert!((top.molar_flow_mol_s + bottom.molar_flow_mol_s - 5.0).abs() < 1e-9);

    // Results survive a save/reload and stay current.
    let mut solved = project.clone();
    solved.results = Some(result);
    let reloaded = Project::from_json(&solved.to_json().unwrap()).unwrap();
    let results = reloaded.results.as_ref().unwrap();
    let fingerprints = reloaded.fingerprints();
    for unit in reloaded.flowsheet.units() {
        assert_eq!(
            results.unit_currency(unit.id, &fingerprints),
            flowsheet_project::Currency::Current
        );
    }
    // Editing a characterisation input makes the results stale.
    let mut edited = reloaded.clone();
    edited.package.pseudo_components[0].specific_gravity += 0.01;
    let fingerprints = edited.fingerprints();
    assert_eq!(
        results.unit_currency(separator, &fingerprints),
        flowsheet_project::Currency::Stale
    );
}

#[test]
fn pseudo_component_selections_fail_explicitly() {
    let (ok, _) = selection(PackageChoice::PengRobinsonBostonMathias);
    assert!(ok.build().is_ok());
    // PR <-> SRK keeps the characterisation and the PH interval.
    let srk = ok.with_model(PackageChoice::SoaveRedlichKwong, KIJ_SOURCE);
    assert_eq!(srk.pseudo_components, ok.pseudo_components);
    assert!(srk.build().is_ok());
    // Non-cubic packages refuse pseudo-components.
    let mut gerg = ok.clone();
    gerg.variant = PackageChoice::Gerg2008;
    gerg.kij.clear();
    gerg.interaction_source.clear();
    assert!(gerg.build().unwrap_err().to_string().contains("pseudo"));
    // A missing definition, a partial interaction matrix, a PH interval
    // outside the ideal-gas domain.
    let mut missing = ok.clone();
    missing.pseudo_components.pop();
    assert!(missing.build().is_err());
    let mut partial = ok.clone();
    partial.kij.pop();
    assert!(partial.build().is_err());
    let mut interval = ok.clone();
    interval.ph_minimum_temperature_k = 200.0;
    assert!(interval.build().is_err());
    // Unknown fields in a stored characterisation are refused.
    let text = Project::new(ok).to_json().unwrap();
    let bad = text.replacen(
        "\"specific_gravity\"",
        "\"critical_volume\": 1.0, \"specific_gravity\"",
        1,
    );
    assert!(Project::from_json(&bad).is_err());
}

#[test]
fn reactors_reject_pseudo_components_before_solving() {
    let (package, x) = selection(PackageChoice::PengRobinson);
    let mut project = Project::new(package);
    let fs = &mut project.flowsheet;
    let feed = fs.add_unit(UnitKind::Feed);
    let reactor = fs.add_unit(UnitKind::GibbsReactor);
    let product = fs.add_unit(UnitKind::Product);
    fs.set_spec(
        feed,
        UnitSpec::Feed {
            temperature_k: Some(500.0),
            pressure_pa: Some(101_325.0),
            molar_flow_mol_s: Some(1.0),
        },
    );
    fs.set_spec(
        reactor,
        UnitSpec::GibbsReactor {
            inerts: Vec::new(),
            equilibrium: Default::default(),
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        },
    );
    project.feeds.insert(feed, x);
    fs.connect(PortRef::new(feed, 0), PortRef::new(reactor, 0))
        .unwrap();
    fs.connect(PortRef::new(reactor, 1), PortRef::new(product, 0))
        .unwrap();
    let prepared = project.package.prepare().unwrap();
    let failure = solve_project(&project, &prepared, |_| ControlFlow::Continue(()), 0).unwrap_err();
    assert!(
        failure.message.contains("pseudo-components"),
        "{}",
        failure.message
    );
}
