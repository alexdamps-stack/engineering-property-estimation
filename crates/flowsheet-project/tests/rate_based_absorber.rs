//! Rate-based packed absorber (#241) through the project solve entry point.

use std::ops::ControlFlow;

use flowsheet_core::units::{
    ColumnSpec,
    columns::{
        AbsorberSpec, ColumnInternals,
        packing::{
            PackedCalculation, PackedSectionBasis, PackedSectionSpec, RateBasedSpec, StageRounding,
        },
        packing_data::Packing,
    },
};
use flowsheet_core::{PortRef, StreamId, UnitId, UnitKind, UnitSpec};
use flowsheet_project::{
    CubicVariantChoice, PackageSelection, Project, ResultSet, SolveFailure, solve_project,
};
use serde_json::Value;

const FIXTURE: &str = include_str!("../../../validation/columns/rate-based-v1.json");

fn number(value: &Value) -> f64 {
    value.as_f64().unwrap()
}

fn strings(value: &Value) -> Vec<String> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().to_owned())
        .collect()
}

#[derive(Clone, Copy)]
struct Inlet<'a> {
    temperature_k: f64,
    pressure_pa: f64,
    molar_flow_mol_s: f64,
    z: &'a [f64],
}

struct Rig {
    project: Project,
    column: UnitId,
    gas: UnitId,
    outlets: [StreamId; 2],
}

impl Rig {
    fn new(
        variant: CubicVariantChoice,
        ids: &[&str],
        spec: AbsorberSpec,
        gas: Inlet,
        liquid: Inlet,
    ) -> Self {
        let mut project = Project::new(PackageSelection::zero_interactions(
            variant,
            ids.iter().map(|s| s.to_string()).collect(),
            "zero interactions in the independent rate-based references",
        ));
        let fs = &mut project.flowsheet;
        let column = fs.add_unit(UnitKind::Absorber);
        fs.set_spec(column, UnitSpec::Columns(ColumnSpec::Absorber(spec)));
        let mut feeds = Vec::new();
        for (index, inlet) in [gas, liquid].into_iter().enumerate() {
            let fs = &mut project.flowsheet;
            let feed = fs.add_unit(UnitKind::Feed);
            fs.set_spec(
                feed,
                UnitSpec::Feed {
                    temperature_k: Some(inlet.temperature_k),
                    pressure_pa: Some(inlet.pressure_pa),
                    molar_flow_mol_s: Some(inlet.molar_flow_mol_s),
                },
            );
            fs.connect(PortRef::new(feed, 0), PortRef::new(column, index))
                .unwrap();
            project.feeds.insert(feed, inlet.z.to_vec());
            feeds.push(feed);
        }
        let fs = &mut project.flowsheet;
        let outlets = [2, 3].map(|port| {
            let product = fs.add_unit(UnitKind::Product);
            fs.connect(PortRef::new(column, port), PortRef::new(product, 0))
                .unwrap()
        });
        Self {
            project,
            column,
            gas: feeds[0],
            outlets,
        }
    }

    fn solve(&self) -> Result<ResultSet, SolveFailure> {
        solve_project(
            &self.project,
            &self.project.package.prepare().unwrap(),
            |_| ControlFlow::Continue(()),
            0,
        )
    }

    fn spec(&self) -> AbsorberSpec {
        match &self.project.flowsheet.unit(self.column).unwrap().spec {
            UnitSpec::Columns(ColumnSpec::Absorber(spec)) => spec.clone(),
            _ => unreachable!(),
        }
    }

    fn edit(&mut self, change: impl FnOnce(&mut AbsorberSpec)) {
        let mut spec = self.spec();
        change(&mut spec);
        self.project
            .flowsheet
            .set_spec(self.column, UnitSpec::Columns(ColumnSpec::Absorber(spec)));
    }

    fn edit_rate(&mut self, change: impl FnOnce(&mut RateBasedSpec)) {
        self.edit(|spec| match &mut spec.calculation {
            PackedCalculation::RateBased(rate) => change(rate),
            PackedCalculation::Equilibrium => unreachable!(),
        });
    }

    fn issue_fields(&self) -> Vec<&'static str> {
        self.project
            .flowsheet
            .unit(self.column)
            .unwrap()
            .issues()
            .iter()
            .map(|i| i.field)
            .collect()
    }
}

fn packed(height_m: f64, rate: RateBasedSpec) -> AbsorberSpec {
    AbsorberSpec {
        pressure_pa: None,
        internals: ColumnInternals::Packed(PackedSectionSpec {
            hetp_m: None,
            basis: PackedSectionBasis::PackedHeight {
                height_m: Some(height_m),
                rounding: StageRounding::Nearest,
            },
        }),
        calculation: PackedCalculation::RateBased(rate),
    }
}

fn inlet(value: &Value) -> (Vec<f64>, [f64; 3]) {
    (
        value["z"].as_array().unwrap().iter().map(number).collect(),
        [
            number(&value["temperature_K"]),
            number(&value["pressure_Pa"]),
            number(&value["molar_flow_mol_s"]),
        ],
    )
}

fn fixture_rig(case: &Value) -> Rig {
    let variant = match case["variant"].as_str().unwrap() {
        "PR" => CubicVariantChoice::PengRobinson,
        "PR-BM" => CubicVariantChoice::PengRobinsonBostonMathias,
        other => panic!("unexpected variant {other}"),
    };
    let ids = strings(&case["ids"]);
    let ids: Vec<&str> = ids.iter().map(String::as_str).collect();
    let mut spec = packed(
        number(&case["packed_height_m"]),
        RateBasedSpec {
            diameter_m: Some(number(&case["diameter_m"])),
            segments: case["segments"].as_u64().unwrap() as u16,
            packing: match case["packing"].as_str().unwrap() {
                "metal_pall_ring_25" => Packing::MetalPallRing25,
                "metal_pall_ring_50" => Packing::MetalPallRing50,
                other => panic!("unexpected packing {other}"),
            },
            interfacial_area_factor: number(&case["area_factor"]),
            critical_surface_tension_n_m: Some(number(&case["critical_surface_tension_N_m"])),
            gas_transport_neglected: strings(&case["gas_neglected"]),
            liquid_transport_neglected: strings(&case["liquid_neglected"]),
            ..Default::default()
        },
    );
    spec.pressure_pa = Some(number(&case["pressure_Pa"]));
    let (gz, [gt, gp, gf]) = inlet(&case["gas"]);
    let (lz, [lt, lp, lf]) = inlet(&case["liquid"]);
    Rig::new(
        variant,
        &ids,
        spec,
        Inlet {
            temperature_k: gt,
            pressure_pa: gp,
            molar_flow_mol_s: gf,
            z: &gz,
        },
        Inlet {
            temperature_k: lt,
            pressure_pa: lp,
            molar_flow_mol_s: lf,
            z: &lz,
        },
    )
}

#[test]
fn segment_profiles_match_independent_references() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    for case in fixture["cases"].as_array().unwrap() {
        let mut rig = fixture_rig(case);
        let result = rig
            .solve()
            .unwrap_or_else(|e| panic!("{}: {}", case["id"], e.message));
        let id = case["id"].as_str().unwrap();
        let compare = |field: &str, what: &str, actual: f64, expected: f64| {
            let range = &fixture["tolerances"][field];
            let tolerance = number(&range[0]) + number(&range[1]) * expected.abs();
            assert!(
                (actual - expected).abs() <= tolerance,
                "{id} {what}: {actual} differs from {expected} by more than {tolerance}"
            );
        };
        let unit = &result.units[&rig.column];
        assert_eq!(unit.duty_w, 0.0);
        assert!(unit.max_scaled_component_residual.unwrap() < 1e-9);
        assert!(unit.scaled_energy_residual.unwrap().abs() < 1e-9);
        let absorber = unit.absorber.as_ref().unwrap();
        let rate = absorber.rate_based.as_ref().unwrap();
        let profiles = &case["profiles"];
        let scalar = |key: &str, actual: &[f64], field: &str| {
            let expected = profiles[key].as_array().unwrap();
            assert_eq!(actual.len(), expected.len(), "{id} {key}");
            for (k, (a, e)) in actual.iter().zip(expected).enumerate() {
                compare(field, &format!("{key}[{k}]"), *a, number(e));
            }
        };
        let vector = |key: &str, actual: &[Vec<f64>], field: &str| {
            let expected = profiles[key].as_array().unwrap();
            for (k, (a, e)) in actual.iter().zip(expected).enumerate() {
                for (i, (a, e)) in a.iter().zip(e.as_array().unwrap()).enumerate() {
                    compare(field, &format!("{key}[{k}][{i}]"), *a, number(e));
                }
            }
        };
        scalar(
            "vapor_temperatures_K",
            &rate.vapor_temperatures_k,
            "temperature_K",
        );
        scalar(
            "liquid_temperatures_K",
            &rate.liquid_temperatures_k,
            "temperature_K",
        );
        scalar(
            "interface_temperatures_K",
            &rate.interface_temperatures_k,
            "temperature_K",
        );
        scalar(
            "vapor_flows_mol_s",
            &rate.vapor_flows_mol_s,
            "molar_flow_mol_s",
        );
        scalar(
            "liquid_flows_mol_s",
            &rate.liquid_flows_mol_s,
            "molar_flow_mol_s",
        );
        vector("vapor_z", &rate.vapor_mole_fractions, "composition");
        vector("liquid_z", &rate.liquid_mole_fractions, "composition");
        vector("fluxes_mol_s", &rate.fluxes_mol_s, "flux_mol_s");
        scalar(
            "interfacial_area_m2_m3",
            &rate.interfacial_area_m2_m3,
            "hydraulics",
        );
        scalar("flooding_fractions", &rate.flooding_fractions, "hydraulics");
        scalar(
            "pressure_drops_Pa_m",
            &rate.pressure_drops_pa_m,
            "hydraulics",
        );
        compare(
            "hydraulics",
            "pressure drop",
            rate.pressure_drop_pa,
            number(&case["pressure_drop_Pa"]),
        );
        compare(
            "hydraulics",
            "max flooding",
            rate.max_flooding_fraction,
            number(&case["max_flooding_fraction"]),
        );
        assert_eq!(absorber.stage_temperatures_k, rate.liquid_temperatures_k);
        assert_eq!(
            absorber.packed_height_m,
            Some(number(&case["packed_height_m"]))
        );
        assert!(
            unit.warnings
                .iter()
                .any(|w| w.contains("neglect trace components"))
        );
        for (k, expected) in case["outlets"].as_array().unwrap().iter().enumerate() {
            let stream = &result.streams[&rig.outlets[k]];
            compare(
                "molar_flow_mol_s",
                "outlet flow",
                stream.molar_flow_mol_s,
                number(&expected["molar_flow_mol_s"]),
            );
            let state = stream.state.as_ref().unwrap();
            assert_eq!(state.pressure_pa, number(&case["pressure_Pa"]));
            compare(
                "temperature_K",
                "outlet T",
                state.temperature_k,
                number(&expected["temperature_K"]),
            );
            for (a, e) in state
                .mole_fractions
                .iter()
                .zip(expected["z"].as_array().unwrap())
            {
                compare("composition", "outlet z", *a, number(e));
            }
        }
        // Methanol is absorbed: less leaves with the gas than enters.
        let top = &result.streams[&rig.outlets[0]];
        let gas_in = number(&case["gas"]["molar_flow_mol_s"]) * number(&case["gas"]["z"][1]);
        assert!(top.molar_flow_mol_s * top.state.as_ref().unwrap().mole_fractions[1] < gas_in);

        // Project round trip with the rate-based record and its results.
        rig.project.results = Some(result);
        let saved = rig.project.to_json().unwrap();
        let reopened = Project::from_json(&saved).unwrap();
        assert_eq!(reopened.to_json().unwrap(), saved);
        assert_eq!(
            reopened.flowsheet.unit(rig.column).unwrap().spec,
            rig.project.flowsheet.unit(rig.column).unwrap().spec
        );
    }
}

/// Nitrogen/carbon dioxide gas washed with water: the second system, used
/// for the limiting behaviour (the package's PH flash of this feed blend
/// succeeds, so the equilibrium mode solves the same column).
fn co2_wash(rate_based: bool) -> Rig {
    let mut spec = packed(
        2.0,
        RateBasedSpec {
            diameter_m: Some(0.8),
            segments: 4,
            packing: Packing::MetalPallRing50,
            // Case input of the Onda wetted area (no library value).
            critical_surface_tension_n_m: Some(0.075),
            gas_transport_neglected: vec!["water".into()],
            liquid_transport_neglected: vec!["nitrogen".into(), "carbon-dioxide".into()],
            ..Default::default()
        },
    );
    spec.pressure_pa = Some(3e5);
    if !rate_based {
        spec.calculation = PackedCalculation::Equilibrium;
    }
    Rig::new(
        CubicVariantChoice::PengRobinson,
        &["nitrogen", "carbon-dioxide", "water"],
        spec,
        Inlet {
            temperature_k: 300.0,
            pressure_pa: 3e5,
            molar_flow_mol_s: 10.0,
            z: &[0.8, 0.2, 0.0],
        },
        Inlet {
            temperature_k: 288.0,
            pressure_pa: 3.1e5,
            molar_flow_mol_s: 200.0,
            z: &[0.0, 0.0, 1.0],
        },
    )
}

#[test]
fn large_transfer_rates_approach_the_equilibrium_stages() {
    let mut rig = co2_wash(true);
    let finite = rig.solve().unwrap();
    rig.edit_rate(|rate| rate.interfacial_area_factor = 1e5);
    let limit = rig.solve().unwrap();
    // The same column with one ideal stage per segment.
    let mut ideal = co2_wash(false);
    ideal.edit(|spec| spec.internals = ColumnInternals::Trays(4));
    let ideal_result = ideal.solve().unwrap();
    let stages = &ideal_result.units[&ideal.column]
        .absorber
        .as_ref()
        .unwrap()
        .stage_temperatures_k;
    let rate = limit.units[&rig.column].absorber.as_ref().unwrap();
    let profile = rate.rate_based.as_ref().unwrap();
    for (k, &t) in stages.iter().enumerate() {
        for actual in [
            profile.vapor_temperatures_k[k],
            profile.liquid_temperatures_k[k],
            profile.interface_temperatures_k[k],
        ] {
            assert!((actual - t).abs() < 2e-3, "segment {k}: {actual} vs {t}");
        }
    }
    let co2_out = |result: &ResultSet, rig: &Rig| {
        let s = &result.streams[&rig.outlets[0]];
        s.molar_flow_mol_s * s.state.as_ref().unwrap().mole_fractions[1]
    };
    let (ideal_co2, limit_co2, finite_co2) = (
        co2_out(&ideal_result, &ideal),
        co2_out(&limit, &rig),
        co2_out(&finite, &rig),
    );
    assert!(
        (limit_co2 / ideal_co2 - 1.0).abs() < 1e-4,
        "{limit_co2} vs {ideal_co2}"
    );
    // Finite transfer rates absorb less than ideal stages.
    assert!(
        finite_co2 > limit_co2 * (1.0 + 1e-6),
        "{finite_co2} {limit_co2}"
    );
}

#[test]
fn hydraulic_and_transport_limits_fail_explicitly() {
    let fixture: Value = serde_json::from_str(FIXTURE).unwrap();
    let case = &fixture["cases"][0];
    // Narrow column: the gas load exceeds flooding.
    let mut rig = fixture_rig(case);
    rig.edit_rate(|rate| rate.diameter_m = Some(0.15));
    let error = rig.solve().unwrap_err();
    assert!(
        error.message.contains("flooding exceeded"),
        "{}",
        error.message
    );
    // Wide column: the liquid no longer wets the packing.
    let mut rig = fixture_rig(case);
    rig.edit_rate(|rate| rate.diameter_m = Some(1.2));
    let error = rig.solve().unwrap_err();
    assert!(
        error.message.contains("below the minimum"),
        "{}",
        error.message
    );
    // Without the trace basis the liquid transport properties are outside
    // their domain (dissolved nitrogen above 0.75 Tc).
    let mut rig = fixture_rig(case);
    rig.edit_rate(|rate| rate.liquid_transport_neglected.clear());
    let error = rig.solve().unwrap_err();
    assert!(
        error.message.contains("outside validity domain"),
        "{}",
        error.message
    );
    // A neglected component above the allowed trace fraction.
    let mut rig = fixture_rig(case);
    rig.project.flowsheet.set_spec(
        rig.gas,
        UnitSpec::Feed {
            temperature_k: Some(300.0),
            pressure_pa: Some(2.6e5),
            molar_flow_mol_s: Some(15.0),
        },
    );
    rig.project.feeds.insert(rig.gas, vec![0.975, 0.025, 0.0]);
    let error = rig.solve().unwrap_err();
    assert!(
        error.message.contains("allowed for a component neglected"),
        "{}",
        error.message
    );
    // An unknown component in the transport basis.
    let mut rig = fixture_rig(case);
    rig.edit_rate(|rate| rate.gas_transport_neglected.push("argon".into()));
    let error = rig.solve().unwrap_err();
    assert!(error.message.contains("argon"), "{}", error.message);
    // Above 3 bar absolute the Robbins correlation is not applicable.
    let mut rig = co2_wash(true);
    rig.edit(|spec| spec.pressure_pa = Some(3.05e5));
    rig.project.flowsheet.set_spec(
        rig.gas,
        UnitSpec::Feed {
            temperature_k: Some(300.0),
            pressure_pa: Some(3.1e5),
            molar_flow_mol_s: Some(10.0),
        },
    );
    let error = rig.solve().unwrap_err();
    assert!(
        error.message.contains("not suitable above 3 bar"),
        "{}",
        error.message
    );
}

#[test]
fn rate_based_mode_needs_a_packed_height_and_valid_inputs() {
    let mut rig = co2_wash(true);
    assert!(rig.issue_fields().is_empty(), "{:?}", rig.issue_fields());
    rig.edit(|spec| spec.internals = ColumnInternals::Trays(4));
    assert_eq!(rig.issue_fields(), ["packed height"]);
    rig.edit_rate(|rate| {
        rate.diameter_m = None;
        rate.critical_surface_tension_n_m = None;
        rate.segments = 0;
        rate.interfacial_area_factor = 0.0;
        rate.gas_transport_neglected = vec!["water".into(), "water".into()];
    });
    assert_eq!(
        rig.issue_fields(),
        [
            "packed height",
            "column diameter",
            "packing critical surface tension",
            "segments",
            "interfacial area factor",
            "gas transport neglected components"
        ]
    );
}

#[test]
fn legacy_absorber_records_load_as_equilibrium_stages() {
    let rig = co2_wash(false);
    let saved = rig.project.to_json().unwrap();
    assert!(!saved.contains("rate_based"));
    let reopened = Project::from_json(&saved).unwrap();
    match &reopened.flowsheet.unit(rig.column).unwrap().spec {
        UnitSpec::Columns(ColumnSpec::Absorber(spec)) => {
            assert_eq!(spec.calculation, PackedCalculation::Equilibrium)
        }
        other => panic!("{other:?}"),
    }
}
