//! The #224 ChemSep extension (1-butanol and its NRTL pairs) against the source
//! correlations transcribed in `validation/prepare-chemsep-nrtl-source-data.py`.

use thermo_data::{ComponentDatabase, LiquidDatabase};
use thermo_types::{CalculationError, MOLAR_GAS_CONSTANT as R, Temperature};

fn t(k: f64) -> Temperature {
    Temperature::new(k).unwrap()
}

#[test]
fn butanol_ideal_gas_is_the_chemsep_rpp_polynomial_with_its_formation_basis() {
    let db = ComponentDatabase::embedded().unwrap();
    let b = db.get("1-butanol").unwrap();
    assert_eq!(b.cas(), "71-36-3");
    assert_eq!(b.critical_temperature(), 563.05);
    assert_eq!(b.critical_pressure(), 4.423e6);
    assert_eq!(b.acentric_factor(), 0.591);
    assert_eq!(b.molar_mass().value(), 0.0741216);
    // ChemSep RPPHeatCapacityCp, equation 100, J/(kmol K).
    let rpp = [37141.0, 136.32, 0.55607, -0.00080567, 3.2127e-07];
    for temperature in [200.0, 298.15, 450.0, 999.0, 1000.0] {
        let cp = rpp
            .iter()
            .enumerate()
            .map(|(k, c)| c * f64::powi(temperature, k as i32))
            .sum::<f64>()
            / 1000.0;
        let ig = b.ideal_gas(t(temperature)).unwrap();
        assert!((ig.heat_capacity - cp).abs() < 1e-9 * cp, "{temperature} K");
    }
    // Ideal-gas formation enthalpy and absolute entropy of ChemSep.
    assert!((b.standard_formation_enthalpy() + 274_600.0).abs() < 1e-6);
    assert!((b.standard_entropy() - 361.48).abs() < 1e-9);
    // The polynomial is not extrapolated beyond its 1000 K source range.
    assert!(matches!(
        b.ideal_gas(t(1000.5)),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    let mixed = db.metadata_for(["water", "1-butanol"]);
    assert!(mixed.version.contains("+chemsep-") && mixed.source.contains("components-chemsep"));
    assert_eq!(db.metadata_for(["water"]), db.metadata());
    let _ = R;
}

#[test]
fn butanol_liquid_correlations_and_pairs_follow_the_chemsep_extract() {
    let liquids = LiquidDatabase::embedded().unwrap();
    let b = liquids.pure("1-butanol").unwrap();
    // ChemSep equation 101 (Pa) at 350 K and its analytical derivatives.
    let (a, bb, c, d, e) = (94.3424, -9140.942, -10.0038, 1.706963e-06, 2.0);
    let ln_p = |x: f64| a + bb / x + c * x.ln() + d * x.powf(e);
    let vp = b.vapor_pressure(t(350.0)).unwrap();
    assert!((vp.ln_pressure - ln_p(350.0)).abs() < 1e-12);
    let h = 1e-3;
    let numeric = (ln_p(350.0 + h) - ln_p(350.0 - h)) / (2.0 * h);
    assert!((vp.d_ln_pressure_dt - numeric).abs() < 1e-8);
    let numeric2 = (ln_p(350.0 + h) - 2.0 * ln_p(350.0) + ln_p(350.0 - h)) / (h * h);
    assert!((vp.d2_ln_pressure_dt2 - numeric2).abs() < 1e-5);
    // Normal boiling point of the correlation, about 391.3 K.
    let p_nbp = b.vapor_pressure(t(391.31)).unwrap().pressure;
    assert!((p_nbp / 101_325.0 - 1.0).abs() < 1e-3);
    // ChemSep equation 105 (kmol/m3).
    let rho = 1000.0 * 0.7931 / 0.24201f64.powf(1.0 + (1.0 - 298.15 / 563.05f64).powf(0.22263));
    assert!((b.saturated_liquid_molar_density(t(298.15)).unwrap() - rho).abs() < 1e-9 * rho);
    assert!(matches!(
        b.saturated_liquid_molar_density(t(550.0)),
        Err(CalculationError::OutsideValidityDomain(_))
    ));
    // Directed DECHEMA NRTL pair, A_ij in cal/mol: tau_ij = A_ij 4.184/(R T).
    let pair = liquids.nrtl_binary("water", "1-butanol").unwrap();
    assert!((pair.tau_ij.t1 - 2633.6951 * 4.184 / R).abs() < 1e-9);
    assert!((pair.tau_ji.t1 - 504.0381 * 4.184 / R).abs() < 1e-9);
    assert_eq!(pair.alpha.a0, 0.4447);
    assert!(pair.regression_interval_inferred);
    let reverse = liquids.nrtl_binary("1-butanol", "water").unwrap();
    assert_eq!(reverse.tau_ij, pair.tau_ji);
    assert!(liquids.nrtl_binary("methanol", "1-butanol").is_ok());
    let mixed = liquids.metadata_for(["water", "1-butanol"]);
    assert!(mixed.version.contains("+chemsep-"));
    assert_eq!(
        liquids.metadata_for(["methanol", "water"]),
        liquids.metadata()
    );
}
