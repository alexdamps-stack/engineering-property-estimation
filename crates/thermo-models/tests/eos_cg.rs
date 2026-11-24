use thermo_models::EosCg2021;

#[test]
fn every_pure_component_and_binary_pair_evaluates_in_the_dilute_region() {
    for (i, &first) in EosCg2021::COMPONENTS.iter().enumerate() {
        let pure = EosCg2021::new(&[first], &[1.0]).unwrap();
        let state = pure
            .evaluate(345.0, 1.0)
            .unwrap_or_else(|error| panic!("{first}: {error}"));
        assert!(state.cp > state.cv && state.speed_of_sound > 0.0, "{first}");
        for &second in &EosCg2021::COMPONENTS[i + 1..] {
            let binary = EosCg2021::new(&[first, second], &[0.3, 0.7]).unwrap();
            let reverse = EosCg2021::new(&[second, first], &[0.7, 0.3]).unwrap();
            let state = binary
                .evaluate(345.0, 1.0)
                .unwrap_or_else(|error| panic!("{first}/{second}: {error}"));
            let reversed = reverse.evaluate(345.0, 1.0).unwrap();
            assert!(state.pressure.is_finite() && state.cp > state.cv);
            assert_eq!(state.pressure, reversed.pressure);
            assert!(
                (state.ln_fugacity_coefficients[0] - reversed.ln_fugacity_coefficients[1]).abs()
                    < 1e-12
            );
        }
    }
}

#[test]
fn binary_pressure_and_residual_helmholtz_match_teqp() {
    let reference = include_str!("../../../validation/fixtures/eoscg2021-co2-h2-teqp.csv");
    assert_eq!(reference.lines().count(), 37);
    for line in reference.lines().skip(1) {
        let values: Vec<f64> = line
            .split(',')
            .take(5)
            .map(|value| value.parse().unwrap())
            .collect();
        let model = EosCg2021::new(
            &["carbon-dioxide", "hydrogen"],
            &[values[0], 1.0 - values[0]],
        )
        .unwrap();
        let (pressure, slope, alpha) = model.isotherm(values[1], values[2]).unwrap();
        assert!((alpha - values[3]).abs() < 1e-10, "{line}: alpha = {alpha}",);
        assert!(
            (pressure / values[4] - 1.0).abs() < 1e-9,
            "{line}: pressure = {pressure}",
        );
        if slope <= 0.0 || pressure <= 0.0 {
            assert!(model.evaluate(values[1], values[2]).is_err(), "{line}");
        }
    }
}
