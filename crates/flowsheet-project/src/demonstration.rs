//! The A8 demonstration process as an editable project.

use std::collections::BTreeMap;

use flowsheet_core::{PortRef, PressureSpec, ThermalSpec, UnitKind, UnitSpec};

use crate::{CubicVariantChoice, DisplayUnitSet, FeedFractions, PackageSelection, Project};

/// Component order of the demonstration.
pub const DEMONSTRATION_COMPONENTS: [&str; 5] =
    ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"];

/// Interaction source of the demonstration package, identical to the
/// headless A8 example and the A2 process qualification.
pub const DEMONSTRATION_INTERACTIONS: &str = "explicit synthetic zero-kij approximation; A2 fixed demonstration, not measured natural-gas accuracy";

/// The qualified A8 demonstration: two natural-gas feeds (330 K and 350 K at
/// 5 MPa, 100 and 20 mol/s) are mixed, cooled to 240 K, expanded to 1.5 MPa
/// and separated. The package is PR or PR-BM with the explicit synthetic
/// zero-kij approximation and a 200..400 K PH interval, exactly as in
/// `crates/flowsheet-core/examples/support/mod.rs`. Stream names match the
/// independent process fixture `validation/fixtures/process-v1.json`.
pub fn demonstration(variant: CubicVariantChoice) -> Project {
    let mut package = PackageSelection::zero_interactions(
        variant,
        DEMONSTRATION_COMPONENTS.map(String::from).to_vec(),
        DEMONSTRATION_INTERACTIONS,
    );
    package.ph_minimum_temperature_k = 200.0;
    package.ph_maximum_temperature_k = 400.0;
    let mut project = Project::new(package);
    let fs = &mut project.flowsheet;
    let mut layout = BTreeMap::new();
    let mut place = |fs: &mut flowsheet_core::Flowsheet, kind, x: f32, y: f32| {
        let id = fs.add_unit(kind);
        layout.insert(id, [x, y]);
        id
    };
    let a = place(fs, UnitKind::Feed, -100.0, 0.0);
    let b = place(fs, UnitKind::Feed, -100.0, 120.0);
    let mix = place(fs, UnitKind::Mixer, 160.0, 60.0);
    let cool = place(fs, UnitKind::Cooler, 300.0, 60.0);
    let valve = place(fs, UnitKind::Valve, 440.0, 60.0);
    let sep = place(fs, UnitKind::Separator, 580.0, 60.0);
    let vapor = place(fs, UnitKind::Product, 740.0, -60.0);
    let liquid = place(fs, UnitKind::Product, 740.0, 180.0);
    fs.rename_unit(vapor, "Vapor");
    fs.rename_unit(liquid, "Liquid");
    fs.set_spec(
        a,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(100.0),
        },
    );
    fs.set_spec(
        b,
        UnitSpec::Feed {
            temperature_k: Some(350.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(20.0),
        },
    );
    fs.set_spec(
        cool,
        UnitSpec::HeatExchange {
            thermal: ThermalSpec::OutletTemperature(Some(240.0)),
            pressure_drop_pa: 0.0,
        },
    );
    fs.set_spec(
        valve,
        UnitSpec::Valve {
            pressure: PressureSpec::OutletPressure(Some(1.5e6)),
        },
    );
    for (name, from, fp, to, tp) in [
        ("A", a, 0, mix, 0),
        ("B", b, 0, mix, 1),
        ("mix", mix, 2, cool, 0),
        ("cool", cool, 1, valve, 0),
        ("expand", valve, 1, sep, 0),
        ("separator_vapor", sep, 1, vapor, 0),
        ("separator_liquid", sep, 2, liquid, 0),
    ] {
        let id = fs
            .connect(PortRef::new(from, fp), PortRef::new(to, tp))
            .expect("fixed demonstration topology");
        fs.rename_stream(id, name);
    }
    project.layout = layout;
    project.feeds = FeedFractions::from([
        (a, vec![0.88, 0.06, 0.03, 0.02, 0.01]),
        (b, vec![0.70, 0.15, 0.13, 0.01, 0.01]),
    ]);
    project.display_units = DisplayUnitSet::ProcessMetric;
    project
}
