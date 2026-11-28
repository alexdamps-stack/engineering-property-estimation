//! Molecular diffusivities for rate-based mass transfer (#241).
//!
//! Package-independent side properties, like the transport properties of
//! #225: nothing here feeds back into equilibrium or caloric results. Every
//! constant was checked against a readable reproduction of the original;
//! see `docs/a17-columns-rate-based.md` § "Verification of constants".
//!
//! * **Gas binaries**: Fuller, Schettler and Giddings (1966), Ind. Eng. Chem.
//!   58(5), 18–27, in the original form
//!   \(D_{AB} = 10^{-3}\,T^{1.75}(1/M_A + 1/M_B)^{1/2} / (P[\Sigma_{vA}^{1/3} + \Sigma_{vB}^{1/3}]^2)\),
//!   cm²/s with T in K, P in atm and M in g/mol, and the 1966 atomic and
//!   molecular diffusion volumes (reproduced in Cussler, *Diffusion*, 3rd ed.
//!   2009, Eq. 5.1-9 and Table 5.1-4).
//! * **Gas mixtures**: Wilke (1950) for a component diffusing through the
//!   other (stagnant) components, \(D_{i,m} = (1 - y_i) / \sum_{j \ne i} y_j / D_{ij}\)
//!   (as given, citing Wilke, in eq. 13 of Getzinger, LBNL report
//!   "Steady-state diffusion in ternary gas mixtures", 1962,
//!   <https://escholarship.org/uc/item/0w40t47j>).
//! * **Liquids**: Wilke and Chang (1955), AIChE J. 1, 264–270,
//!   \(D = 7.4\times10^{-8} (\phi_B M_B)^{1/2} T / (\mu_B V_A^{0.6})\) in cm²/s
//!   with μ in cP and \(V_A\) the solute molar volume at its normal boiling
//!   point in cm³/mol, for a solute A dilute in a solvent B; association
//!   factors φ: water 2.6, methanol 1.9, unassociated solvents 1.0.
//! * **Normal-boiling molar volume**: Tyn and Calus (1975),
//!   \(V_b = 0.285\,V_c^{1.048}\) (cm³/mol), from the ChemSep critical
//!   volume of the transport extract.
//!
//! Declared domain: components with a diffusion volume in [`fuller_volume`]
//! (the eight transport components), gas 200–600 K and at most
//! [`FULLER_MAX_PRESSURE_PA`] (a declared limit: the Fuller method is a
//! low-pressure correlation and no dense-gas correction is applied), liquid
//! 250–500 K. Outside it the result is `OutsideValidityDomain`.
use thermo_data::{ComponentDatabase, TransportDatabase};
use thermo_types::{CalculationError, InputError};

/// Implementation version of the diffusivity formulations; bump on any change.
pub const DIFFUSIVITY_MODEL_VERSION: &str = "diffusivity-v2";
/// Largest pressure of the Fuller gas diffusivity, Pa (declared, not a
/// literature constant).
pub const FULLER_MAX_PRESSURE_PA: f64 = 2.0e6;
/// Temperature domain of the gas diffusivity, K.
pub const GAS_DIFFUSIVITY_TEMPERATURE_RANGE_K: [f64; 2] = [200.0, 600.0];
/// Temperature domain of the liquid diffusivity, K.
pub const LIQUID_DIFFUSIVITY_TEMPERATURE_RANGE_K: [f64; 2] = [250.0, 500.0];
/// Standard atmosphere, Pa.
const ATMOSPHERE_PA: f64 = 101_325.0;

/// Fuller–Schettler–Giddings (1966) diffusion volume Σv of a component,
/// cm³/mol. Simple molecules use their tabulated molecular volume (N2 17.9,
/// O2 16.6, CO2 26.9, H2O 12.7); hydrocarbons and methanol are sums of the
/// atomic increments C 16.5, H 1.98, O 5.48.
pub fn fuller_volume(id: &str) -> Option<f64> {
    const C: f64 = 16.5;
    const H: f64 = 1.98;
    const O: f64 = 5.48;
    Some(match id {
        "nitrogen" => 17.9,
        "oxygen" => 16.6,
        "carbon-dioxide" => 26.9,
        "water" => 12.7,
        "methane" => C + 4.0 * H,
        "ethane" => 2.0 * C + 6.0 * H,
        "propane" => 3.0 * C + 8.0 * H,
        "methanol" => C + 4.0 * H + O,
        _ => return None,
    })
}

/// Wilke–Chang association factor φ of a solvent (Wilke and Chang 1955):
/// water 2.6, methanol 1.9, unassociated solvents 1.0.
pub fn wilke_chang_association(id: &str) -> f64 {
    match id {
        "water" => 2.6,
        "methanol" => 1.9,
        _ => 1.0,
    }
}

/// Fuller–Schettler–Giddings binary gas diffusivity, m²/s, at `t_k` and
/// `p_pa` for molar masses in g/mol and diffusion volumes in cm³/mol.
pub fn fuller_binary(t_k: f64, p_pa: f64, m_a: f64, m_b: f64, v_a: f64, v_b: f64) -> f64 {
    let sum = v_a.cbrt() + v_b.cbrt();
    1e-3 * t_k.powf(1.75) * (1.0 / m_a + 1.0 / m_b).sqrt() / (p_pa / ATMOSPHERE_PA * sum * sum)
        * 1e-4
}

/// Wilke–Chang liquid diffusivity, m²/s, of a dilute solute (normal-boiling
/// molar volume in cm³/mol) in a solvent with association factor `phi`,
/// molar mass in g/mol and viscosity in Pa s.
pub fn wilke_chang(
    t_k: f64,
    phi: f64,
    solvent_molar_mass: f64,
    viscosity_pa_s: f64,
    solute_volume_cm3: f64,
) -> f64 {
    7.4e-8 * (phi * solvent_molar_mass).sqrt() * t_k
        / (viscosity_pa_s * 1e3 * solute_volume_cm3.powf(0.6))
        * 1e-4
}

/// Tyn–Calus molar volume at the normal boiling point, cm³/mol, from the
/// critical volume in cm³/mol.
pub fn tyn_calus_volume(critical_volume_cm3: f64) -> f64 {
    0.285 * critical_volume_cm3.powf(1.048)
}

#[derive(Debug, Clone)]
struct Species {
    id: String,
    molar_mass_g_mol: f64,
    fuller_volume: f64,
    association: f64,
    boiling_volume_cm3: f64,
}

/// Diffusivity model for an ordered component selection.
#[derive(Debug, Clone)]
pub struct DiffusivityModel {
    species: Vec<Species>,
}

fn outside(message: String) -> CalculationError {
    CalculationError::OutsideValidityDomain(message)
}

impl DiffusivityModel {
    /// Model for `ids` in package order. Every component needs a Fuller
    /// volume, a component record (molar mass) and a transport record
    /// (critical volume); otherwise `MissingParameter`.
    pub fn new(ids: &[&str]) -> Result<Self, CalculationError> {
        if ids.is_empty() {
            return Err(CalculationError::InvalidInput(InputError::EmptyComposition));
        }
        let components = ComponentDatabase::embedded()?;
        let transport = TransportDatabase::embedded()?;
        let species = ids
            .iter()
            .map(|&id| {
                let fuller = fuller_volume(id).ok_or_else(|| {
                    CalculationError::MissingParameter(format!(
                        "no Fuller diffusion volume for {id}"
                    ))
                })?;
                Ok(Species {
                    id: id.to_owned(),
                    molar_mass_g_mol: components.get(id)?.molar_mass().value() * 1e3,
                    fuller_volume: fuller,
                    association: wilke_chang_association(id),
                    boiling_volume_cm3: tyn_calus_volume(
                        transport.get(id)?.critical_volume() * 1e6,
                    ),
                })
            })
            .collect::<Result<Vec<_>, CalculationError>>()?;
        Ok(Self { species })
    }

    /// Number of components.
    pub fn len(&self) -> usize {
        self.species.len()
    }

    /// Whether the model has no components (never true for a built model).
    pub fn is_empty(&self) -> bool {
        self.species.is_empty()
    }

    /// Fuller binary diffusivity of components `i` and `j`, m²/s.
    pub fn gas_binary(&self, i: usize, j: usize, t_k: f64, p_pa: f64) -> f64 {
        let (a, b) = (&self.species[i], &self.species[j]);
        fuller_binary(
            t_k,
            p_pa,
            a.molar_mass_g_mol,
            b.molar_mass_g_mol,
            a.fuller_volume,
            b.fuller_volume,
        )
    }

    fn check(&self, x: &[f64]) -> Result<(), CalculationError> {
        if x.len() != self.species.len() {
            return Err(CalculationError::InvalidInput(
                InputError::ComponentMismatch,
            ));
        }
        if x.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(CalculationError::InvalidInput(InputError::NonFinite(
                "diffusivity mole fractions",
            )));
        }
        Ok(())
    }

    /// Effective (Wilke 1950) gas diffusivity of every component in the
    /// mixture `y`, m²/s. A component without another component present fails.
    pub fn gas_effective(
        &self,
        y: &[f64],
        t_k: f64,
        p_pa: f64,
    ) -> Result<Vec<f64>, CalculationError> {
        self.check(y)?;
        let [lo, hi] = GAS_DIFFUSIVITY_TEMPERATURE_RANGE_K;
        if !(lo..=hi).contains(&t_k) {
            return Err(outside(format!(
                "gas diffusivity: T = {t_k} K outside {lo}..={hi} K"
            )));
        }
        if !(p_pa > 0.0 && p_pa <= FULLER_MAX_PRESSURE_PA) {
            return Err(outside(format!(
                "gas diffusivity: p = {p_pa} Pa outside the Fuller low-pressure domain (at most {FULLER_MAX_PRESSURE_PA} Pa)"
            )));
        }
        let total: f64 = y.iter().sum();
        (0..y.len())
            .map(|i| {
                let others = total - y[i];
                let resistance: f64 = (0..y.len())
                    .filter(|&j| j != i)
                    .map(|j| y[j] / self.gas_binary(i, j, t_k, p_pa))
                    .sum();
                if others <= 0.0 || resistance <= 0.0 {
                    return Err(CalculationError::InvalidInput(InputError::NotPositive(
                        "gas diffusivity needs a second component",
                    )));
                }
                Ok(others / resistance)
            })
            .collect()
    }

    /// Wilke–Chang diffusivity of every component other than `solvent`,
    /// treated as a dilute solute in that solvent, m²/s, with the liquid's
    /// viscosity (Pa s) as the solvent viscosity. The solvent's own entry is
    /// zero.
    pub fn liquid_in_solvent(
        &self,
        solvent: usize,
        t_k: f64,
        viscosity_pa_s: f64,
    ) -> Result<Vec<f64>, CalculationError> {
        let [lo, hi] = LIQUID_DIFFUSIVITY_TEMPERATURE_RANGE_K;
        if !(lo..=hi).contains(&t_k) {
            return Err(outside(format!(
                "liquid diffusivity: T = {t_k} K outside {lo}..={hi} K"
            )));
        }
        if !(viscosity_pa_s.is_finite() && viscosity_pa_s > 0.0) {
            return Err(CalculationError::InvalidInput(InputError::NotPositive(
                "liquid diffusivity viscosity",
            )));
        }
        let b = self
            .species
            .get(solvent)
            .ok_or(CalculationError::InvalidInput(
                InputError::ComponentMismatch,
            ))?;
        Ok(self
            .species
            .iter()
            .enumerate()
            .map(|(i, a)| {
                if i == solvent {
                    0.0
                } else {
                    wilke_chang(
                        t_k,
                        b.association,
                        b.molar_mass_g_mol,
                        viscosity_pa_s,
                        a.boiling_volume_cm3,
                    )
                }
            })
            .collect())
    }

    /// Component identifiers in model order.
    pub fn component_ids(&self) -> Vec<&str> {
        self.species.iter().map(|s| s.id.as_str()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fuller_reproduces_published_worked_examples() {
        // Cussler, Diffusion, 3rd ed. (2009), Example 5.1-2: H2 in N2 at
        // 21 °C (294 K) and 2 atm with V(H2) = 7.07, V(N2) = 17.9 and
        // M = 2.02, 28.0: D = 0.37 cm²/s (measured 0.38).
        let d = fuller_binary(294.0, 2.0 * ATMOSPHERE_PA, 2.02, 28.0, 7.07, 17.9) * 1e4;
        assert!((d - 0.37).abs() < 0.005, "{d}");
        // n-butanol in air at 273 K and 1 atm, Σv = 4(16.5) + 10(1.98) +
        // 5.48 = 91.28 and 20.1, M = 74.1 and 29: D = 7.73e-6 m²/s (worked
        // example of the University of Jordan Transport Phenomena II notes,
        // "03- Diffusion coefficients", after Geankoplis).
        let v_butanol: f64 = 4.0 * 16.5 + 10.0 * 1.98 + 5.48;
        assert!((v_butanol - 91.28).abs() < 1e-12);
        // The printed result carries the example's intermediate rounding:
        // the same inputs give 7.701e-6 m²/s (0.4 % below the printed value).
        let d = fuller_binary(273.0, ATMOSPHERE_PA, 74.1, 29.0, v_butanol, 20.1);
        assert!((d / 7.73e-6 - 1.0).abs() < 0.005, "{d}");
        // D ∝ 1/p.
        let half = fuller_binary(273.0, 2.0 * ATMOSPHERE_PA, 74.1, 29.0, v_butanol, 20.1);
        assert!((half * 2.0 / d - 1.0).abs() < 1e-14);
    }

    #[test]
    fn wilke_chang_reproduces_a_published_worked_example() {
        // Ethylbenzene in water at 293 K: φ = 2.6, M = 18.0, μ = 1.0 cP,
        // V_A = 139.5 cm³/mol gives 0.77e-5 cm²/s (University of Jordan
        // Transport Phenomena II notes, after Wilke and Chang 1955).
        let d = wilke_chang(293.0, 2.6, 18.0, 1.0e-3, 139.5) * 1e4;
        assert!((d - 0.77e-5).abs() < 0.005e-5, "{d}");
        assert!((tyn_calus_volume(94.0) - 0.285 * 94.0_f64.powf(1.048)).abs() < 1e-12);
        let model = DiffusivityModel::new(&["carbon-dioxide", "water"]).unwrap();
        let d = model.liquid_in_solvent(1, 298.15, 0.89e-3).unwrap();
        assert_eq!(d[1], 0.0);
        assert!(d[0] > 1e-9 && d[0] < 3e-9, "{d:?}");
        assert!(model.liquid_in_solvent(1, 600.0, 0.89e-3).is_err());
    }

    #[test]
    fn wilke_mixture_rule_reduces_to_the_binary_and_domain_is_explicit() {
        let model = DiffusivityModel::new(&["methanol", "nitrogen", "water"]).unwrap();
        let d = model.gas_effective(&[0.2, 0.8, 0.0], 300.0, 1e5).unwrap();
        let binary = model.gas_binary(0, 1, 300.0, 1e5);
        assert!((d[0] / binary - 1.0).abs() < 1e-14);
        assert!((d[1] / binary - 1.0).abs() < 1e-14);
        assert!(matches!(
            model.gas_effective(&[0.2, 0.8, 0.0], 300.0, 3e6),
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(model.gas_effective(&[1.0, 0.0, 0.0], 300.0, 1e5).is_err());
        assert!(matches!(
            DiffusivityModel::new(&["hydrogen", "water"]),
            Err(CalculationError::MissingParameter(_))
        ));
    }
}
