//! Canonical A8 demonstration inputs; no reference answers are read here.
use flowsheet_core::{
    FeedCompositions, Flowsheet, PortRef, PressureSpec, ThermalSpec, UnitKind, UnitSpec,
};
use property_packages::CubicPackage;
use thermo_data::ComponentDatabase;
use thermo_models::{CubicModel, CubicVariant};
use thermo_solvers::PhFlashOptions;
use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

pub const COMPONENTS: [&str; 5] = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"];

pub fn composition(fractions: &[f64]) -> Composition {
    Composition::new(
        CompositionBasis::Mole,
        COMPONENTS
            .iter()
            .zip(fractions)
            .map(|(id, x)| (ComponentId::new(*id).unwrap(), *x))
            .collect(),
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap()
}

pub fn demonstration(
    variant: CubicVariant,
    scale: f64,
) -> (Flowsheet, CubicPackage, FeedCompositions) {
    assert!(scale.is_finite() && scale > 0.0);
    let mut fs = Flowsheet::new();
    // Create consumers before their sources to exercise dependency ordering.
    let vapor = fs.add_unit(UnitKind::Product);
    let liquid = fs.add_unit(UnitKind::Product);
    let sep = fs.add_unit(UnitKind::Separator);
    let valve = fs.add_unit(UnitKind::Valve);
    let cool = fs.add_unit(UnitKind::Cooler);
    let mix = fs.add_unit(UnitKind::Mixer);
    let a = fs.add_unit(UnitKind::Feed);
    let b = fs.add_unit(UnitKind::Feed);
    let ca = composition(&[0.88, 0.06, 0.03, 0.02, 0.01]);
    let cb = composition(&[0.70, 0.15, 0.13, 0.01, 0.01]);
    fs.set_spec(
        a,
        UnitSpec::Feed {
            temperature_k: Some(330.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(100.0 * scale),
        },
    );
    fs.set_spec(
        b,
        UnitSpec::Feed {
            temperature_k: Some(350.0),
            pressure_pa: Some(5e6),
            molar_flow_mol_s: Some(20.0 * scale),
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
            .unwrap();
        fs.rename_stream(id, name);
    }
    let model = CubicModel::new(&ComponentDatabase::embedded().unwrap(), &ca, variant, vec![vec![0.0; 5]; 5], "explicit synthetic zero-kij approximation; A2 fixed demonstration, not measured natural-gas accuracy").unwrap();
    let pkg = CubicPackage::new(
        model,
        PhFlashOptions {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 400.0,
            ..Default::default()
        },
    );
    (fs, pkg, [(a, ca), (b, cb)].into())
}
