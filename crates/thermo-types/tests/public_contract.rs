use thermo_types::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-14 + 1e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}
fn id(name: &str) -> ComponentId {
    ComponentId::new(name).unwrap()
}
fn masses() -> Vec<(ComponentId, MolarMass)> {
    vec![
        (id("a"), MolarMass::new(0.010).unwrap()),
        (id("b"), MolarMass::new(0.030).unwrap()),
    ]
}
fn composition(
    basis: CompositionBasis,
    a: f64,
    b: f64,
    policy: NormalizationPolicy,
) -> Result<Composition, InputError> {
    Composition::new(basis, vec![(id("a"), a), (id("b"), b)], policy)
}
fn reference() -> ReferenceConditions {
    ReferenceConditions::new(
        "named test: 300 K, 100 kPa absolute",
        Temperature::new(300.0).unwrap(),
        Pressure::new(100_000.0).unwrap(),
    )
    .unwrap()
}

#[test]
fn finite_and_absolute_boundaries_are_enforced() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            FiniteSi::new(value),
            Err(InputError::NonFinite(_))
        ));
        assert!(Temperature::new(value).is_err());
        assert!(Pressure::new(value).is_err());
        assert!(MolarMass::new(value).is_err());
        assert!(MolarFlow::new(value).is_err());
        assert!(MassFlow::new(value).is_err());
        assert!(PositiveDimensionless::new(value).is_err());
    }
    for value in [0.0, -0.0, -1.0] {
        assert!(matches!(
            Temperature::new(value),
            Err(InputError::NotPositive(_))
        ));
        assert!(Pressure::new(value).is_err());
        assert!(MolarMass::new(value).is_err());
        assert!(PositiveDimensionless::new(value).is_err());
    }
    assert!(Temperature::new(f64::MIN_POSITIVE).is_ok());
    assert!(MolarFlow::new(-1.0).is_err());
    assert!(MassFlow::new(-1.0).is_err());
    close(MolarFlow::new(0.0).unwrap().value(), 0.0);
    close(FiniteSi::new(-123.0).unwrap().value(), -123.0);
}

#[test]
fn explicit_si_conversions_have_known_values() {
    close(Temperature::from_celsius(25.0).unwrap().value(), 298.15);
    assert!(Temperature::from_celsius(-273.15).is_err());
    assert!(Temperature::from_celsius(f64::NAN).is_err());
    close(Pressure::from_bar_absolute(2.5).unwrap().value(), 250_000.0);
    assert!(Pressure::from_bar_absolute(f64::MAX).is_err());
    let mass = MolarFlow::new(10.0)
        .unwrap()
        .to_mass(MolarMass::new(0.020).unwrap())
        .unwrap();
    close(mass.value(), 0.2);
    close(
        mass.to_molar(MolarMass::new(0.020).unwrap())
            .unwrap()
            .value(),
        10.0,
    );
}

#[test]
fn component_identity_is_stable_and_case_sensitive() {
    assert_eq!(id("cas:74-82-8").as_str(), "cas:74-82-8");
    assert_ne!(id("Methane"), id("methane"));
    for invalid in ["", " methane", "methane\n", "two names", "µ"] {
        assert_eq!(
            ComponentId::new(invalid),
            Err(InputError::InvalidIdentifier)
        );
    }
}

#[test]
fn normalization_is_an_explicit_policy() {
    assert_eq!(
        composition(
            CompositionBasis::Mole,
            2.0,
            3.0,
            NormalizationPolicy::RequireUnitSum
        ),
        Err(InputError::NotNormalized)
    );
    let c = composition(
        CompositionBasis::Mass,
        2.0,
        3.0,
        NormalizationPolicy::Normalize,
    )
    .unwrap();
    assert_eq!(c.basis(), CompositionBasis::Mass);
    close(c.entries()[0].1, 0.4);
    close(c.entries()[1].1, 0.6);
    // Accepted near-unit inputs are also explicitly rescaled.
    let c = composition(
        CompositionBasis::Mole,
        0.5,
        0.5 + 5e-11,
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    close(c.entries().iter().map(|(_, f)| f).sum(), 1.0);
    assert!(
        composition(
            CompositionBasis::Mole,
            0.5,
            0.5 + 2e-10,
            NormalizationPolicy::RequireUnitSum
        )
        .is_err()
    );
    let huge = composition(
        CompositionBasis::Mole,
        f64::MAX,
        f64::MAX,
        NormalizationPolicy::Normalize,
    )
    .unwrap();
    close(huge.entries()[0].1, 0.5);
}

#[test]
fn invalid_compositions_have_distinct_failures() {
    assert_eq!(
        Composition::new(
            CompositionBasis::Mole,
            vec![],
            NormalizationPolicy::Normalize
        ),
        Err(InputError::EmptyComposition)
    );
    assert_eq!(
        composition(
            CompositionBasis::Mole,
            0.0,
            0.0,
            NormalizationPolicy::Normalize
        ),
        Err(InputError::ZeroTotal)
    );
    assert!(matches!(
        composition(
            CompositionBasis::Mole,
            -0.1,
            1.1,
            NormalizationPolicy::Normalize
        ),
        Err(InputError::Negative(_))
    ));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            composition(
                CompositionBasis::Mole,
                bad,
                1.0,
                NormalizationPolicy::Normalize
            ),
            Err(InputError::NonFinite(_))
        ));
    }
    assert_eq!(
        Composition::new(
            CompositionBasis::Mass,
            vec![(id("a"), 1.0), (id("a"), 0.0)],
            NormalizationPolicy::Normalize
        ),
        Err(InputError::DuplicateComponent)
    );
}

#[test]
fn mass_mole_conversion_matches_hand_calculated_values() {
    let mole = composition(
        CompositionBasis::Mole,
        0.5,
        0.5,
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let mass = mole.convert(CompositionBasis::Mass, &masses()).unwrap();
    close(mass.entries()[0].1, 0.25);
    close(mass.entries()[1].1, 0.75);
    close(mole.molar_mass(&masses()).unwrap().value(), 0.020);
    close(mass.molar_mass(&masses()).unwrap().value(), 0.020);
    let recovered = mass.convert(CompositionBasis::Mole, &masses()).unwrap();
    close(recovered.entries()[0].1, 0.5);
}

#[test]
fn zero_and_trace_components_preserve_identity_and_order() {
    for trace in [0.0, 1e-20] {
        let c = composition(
            CompositionBasis::Mole,
            trace,
            1.0,
            NormalizationPolicy::Normalize,
        )
        .unwrap();
        let mass = c.convert(CompositionBasis::Mass, &masses()).unwrap();
        assert_eq!(mass.entries().len(), 2);
        assert_eq!(mass.entries()[0].0, id("a"));
        close(mass.entries()[0].1, trace / 3.0);
        if trace > 0.0 {
            assert!(mass.entries()[0].1 > 0.0);
        }
        close(mass.molar_mass(&masses()).unwrap().value(), 0.030);
    }
    let pure = Composition::new(
        CompositionBasis::Mass,
        vec![(id("a"), 1.0)],
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    close(pure.molar_mass(&masses()[..1]).unwrap().value(), 0.010);
}

#[test]
fn molar_mass_order_is_checked_and_joint_permutation_is_invariant() {
    let c = composition(
        CompositionBasis::Mole,
        0.2,
        0.8,
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let mut reversed_masses = masses();
    reversed_masses.reverse();
    assert_eq!(
        c.convert(CompositionBasis::Mass, &reversed_masses),
        Err(InputError::ComponentMismatch)
    );
    assert_eq!(
        c.molar_mass(&masses()[..1]),
        Err(InputError::ComponentMismatch)
    );
    let mut entries = c.entries().to_vec();
    entries.reverse();
    let reversed = Composition::new(
        CompositionBasis::Mole,
        entries,
        NormalizationPolicy::RequireUnitSum,
    )
    .unwrap();
    let expected = c.convert(CompositionBasis::Mass, &masses()).unwrap();
    let actual = reversed
        .convert(CompositionBasis::Mass, &reversed_masses)
        .unwrap();
    close(expected.entries()[0].1, actual.entries()[1].1);
    close(
        c.molar_mass(&masses()).unwrap().value(),
        reversed.molar_mass(&reversed_masses).unwrap().value(),
    );
}

#[test]
fn standard_volume_requires_conditions_and_matches_known_gas_law_values() {
    let v = MolarFlow::new(2.0)
        .unwrap()
        .to_standard_volume(reference(), GasVolumeAssumption::IdealGas)
        .unwrap();
    close(v.value(), 0.049_886_775_708_919_44);
    assert_eq!(v.reference().name(), reference().name());
    close(v.to_molar().unwrap().value(), 2.0);
    let z = GasVolumeAssumption::CompressibilityFactor(PositiveDimensionless::new(0.8).unwrap());
    let real = MolarFlow::new(2.0)
        .unwrap()
        .to_standard_volume(reference(), z)
        .unwrap();
    close(real.value(), 0.039_909_420_567_135_55);
    close(real.to_molar().unwrap().value(), 2.0);
    assert_eq!(real.assumption(), z);
    let hourly = StandardVolumeFlow::from_cubic_metres_per_hour(
        3600.0,
        reference(),
        GasVolumeAssumption::IdealGas,
    )
    .unwrap();
    close(hourly.value(), 1.0);
    close(hourly.to_molar().unwrap().value(), 40.090_785_014_242_016);
    let other = ReferenceConditions::new(
        "600 K",
        Temperature::new(600.0).unwrap(),
        reference().pressure(),
    )
    .unwrap();
    close(
        MolarFlow::new(2.0)
            .unwrap()
            .to_standard_volume(other, GasVolumeAssumption::IdealGas)
            .unwrap()
            .value(),
        2.0 * v.value(),
    );
}

#[test]
fn standard_volume_rejects_invalid_values_and_names() {
    for name in ["", "  ", "line\nbreak"] {
        assert!(
            ReferenceConditions::new(name, reference().temperature(), reference().pressure())
                .is_err()
        );
    }
    for value in [
        -1.0,
        -f64::from_bits(1),
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(
            StandardVolumeFlow::new(value, reference(), GasVolumeAssumption::IdealGas).is_err()
        );
        assert!(
            StandardVolumeFlow::from_cubic_metres_per_hour(
                value,
                reference(),
                GasVolumeAssumption::IdealGas
            )
            .is_err()
        );
    }
    close(
        StandardVolumeFlow::new(0.0, reference(), GasVolumeAssumption::IdealGas)
            .unwrap()
            .to_molar()
            .unwrap()
            .value(),
        0.0,
    );
}

#[test]
fn arithmetic_range_failures_do_not_become_successful_nan_or_zero() {
    assert!(matches!(
        MolarFlow::new(f64::MAX)
            .unwrap()
            .to_mass(MolarMass::new(2.0).unwrap()),
        Err(InputError::NonFinite(_))
    ));
    assert_eq!(
        MolarFlow::new(f64::from_bits(1))
            .unwrap()
            .to_mass(MolarMass::new(0.1).unwrap()),
        Err(InputError::NumericalRange)
    );
    assert!(
        MassFlow::new(f64::MAX)
            .unwrap()
            .to_molar(MolarMass::new(0.1).unwrap())
            .is_err()
    );
    assert_eq!(
        composition(
            CompositionBasis::Mole,
            f64::from_bits(1),
            f64::MAX,
            NormalizationPolicy::Normalize
        ),
        Err(InputError::NumericalRange)
    );
    assert!(
        MolarFlow::new(f64::MAX)
            .unwrap()
            .to_standard_volume(reference(), GasVolumeAssumption::IdealGas)
            .unwrap()
            .to_molar()
            .is_ok()
    );
}
