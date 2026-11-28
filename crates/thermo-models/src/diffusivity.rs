//! Molecular diffusivities for rate-based mass transfer (#241).
//!
//! Package-independent side properties, like the transport properties of
//! #225: nothing here feeds back into equilibrium or caloric results.
//!
//! * **Gas binaries**: Fuller, Schettler and Giddings (1966), Ind. Eng. Chem.
//!   58(5), 18–27, with the atomic and molecular diffusion volumes of Fuller,
//!   Ensley and Giddings (1969), J. Phys. Chem. 73, 3679–3685, in the form of
//!   Poling, Prausnitz and O'Connell, *The Properties of Gases and Liquids*,
//!   5th ed. (2001), eq. 11-4.4 and Table 11-1:
//!   \(D_{AB} = 0.00143\,T^{1.75} / (P\,M_{AB}^{1/2}[\Sigma_v{}_A^{1/3} + \Sigma_v{}_B^{1/3}]^2)\)
//!   in cm²/s with T in K, P in bar and \(M_{AB} = 2/(1/M_A + 1/M_B)\) in g/mol.
//! * **Gas mixtures**: Blanc's law for a component diffusing through the
//!   rest of the mixture (Poling et al. 2001, eq. 11-7.2),
//!   \(D_{i,m} = (1 - y_i) / \sum_{j \ne i} y_j / D_{ij}\).
//! * **Liquids**: Wilke and Chang (1955), AIChE J. 1, 264–270 (Poling et al.
//!   2001, eq. 11-9.1), \(D = 7.4\times10^{-8} (\phi M)^{1/2} T / (\eta V_A^{0.6})\)
//!   in cm²/s with η in cP and \(V_A\) the solute molar volume at its normal
//!   boiling point in cm³/mol. For a mixed solvent \(\phi M\) is the
//!   solvent-fraction-weighted \(\sum_{j\ne A} x_j \phi_j M_j / \sum_{j \ne A} x_j\)
//!   (Poling et al. 2001, §11-10) and η is the liquid-mixture viscosity.
//!   Association factors φ: water 2.6, methanol 1.9, unassociated 1.0
//!   (Wilke and Chang 1955).
//! * **Normal-boiling molar volume**: Tyn and Calus (1975), Processing 21(4),
//!   16–17 (recommended with Wilke–Chang in Poling et al. 2001, §11-9),
//!   \(V_b = 0.285\,V_c^{1.048}\), cm³/mol, from the ChemSep critical volume
//!   of the transport extract.
//!
//! Declared domain: components with a Fuller volume in [`fuller_volume`]
//! (the eight transport components), gas 200–600 K and at most
//! [`FULLER_MAX_PRESSURE_PA`] (the Fuller method is a low-pressure
//! correlation; no Takahashi dense-gas correction is applied), liquid
//! 250–500 K. Outside it the result is `OutsideValidityDomain`.
use thermo_data::{ComponentDatabase, TransportDatabase};
use thermo_types::{CalculationError, InputError};

/// Implementation version of the diffusivity formulations; bump on any change.
pub const DIFFUSIVITY_MODEL_VERSION: &str = "diffusivity-v1";
/// Largest pressure of the Fuller gas diffusivity, Pa.
pub const FULLER_MAX_PRESSURE_PA: f64 = 2.0e6;
/// Temperature domain of the gas diffusivity, K.
pub const GAS_DIFFUSIVITY_TEMPERATURE_RANGE_K: [f64; 2] = [200.0, 600.0];
/// Temperature domain of the liquid diffusivity, K.
pub const LIQUID_DIFFUSIVITY_TEMPERATURE_RANGE_K: [f64; 2] = [250.0, 500.0];

/// Fuller diffusion volume Σv of a component, cm³/mol (Poling et al. 2001,
/// Table 11-1). Simple molecules use their tabulated molecular volume; the
/// hydrocarbons and methanol are sums of the atomic increments
/// C 15.9, H 2.31 and O 6.11.
pub fn fuller_volume(id: &str) -> Option<f64> {
    const C: f64 = 15.9;
    const H: f64 = 2.31;
    const O: f64 = 6.11;
    Some(match id {
        "nitrogen" => 18.5,
        "oxygen" => 16.3,
        "carbon-dioxide" => 26.7,
        "water" => 13.1,
        "methane" => C + 4.0 * H,
        "ethane" => 2.0 * C + 6.0 * H,
        "propane" => 3.0 * C + 8.0 * H,
        "methanol" => C + 4.0 * H + O,
        _ => return None,
    })
}

/// Wilke–Chang association factor φ of a solvent component (Wilke and
/// Chang 1955): water 2.6, methanol 1.9, unassociated solvents 1.0.
pub fn wilke_chang_association(id: &str) -> f64 {
    match id {
        "water" => 2.6,
        "methanol" => 1.9,
        _ => 1.0,
    }
}

/// Fuller binary gas diffusivity, m²/s, at `t_k` and `p_pa` for molar masses
/// in g/mol and diffusion volumes in cm³/mol.
pub fn fuller_binary(t_k: f64, p_pa: f64, m_a: f64, m_b: f64, v_a: f64, v_b: f64) -> f64 {
    let m_ab = 2.0 / (1.0 / m_a + 1.0 / m_b);
    let sum = v_a.cbrt() + v_b.cbrt();
    0.00143 * t_k.powf(1.75) / (p_pa * 1e-5 * m_ab.sqrt() * sum * sum) * 1e-4
}

/// Wilke–Chang liquid diffusivity, m²/s, for the solvent group `phi_m`
/// (φM in g/mol), solvent viscosity in Pa s and solute normal-boiling
/// molar volume in cm³/mol.
pub fn wilke_chang(t_k: f64, phi_m: f64, viscosity_pa_s: f64, solute_volume_cm3: f64) -> f64 {
    7.4e-8 * phi_m.sqrt() * t_k / (viscosity_pa_s * 1e3 * solute_volume_cm3.powf(0.6)) * 1e-4
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

    /// Effective (Blanc) gas diffusivity of every component in the mixture
    /// `y`, m²/s. A component without another component present fails.
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

    /// Wilke–Chang diffusivity of every component of the liquid `x` through
    /// the other components, m²/s, with the liquid-mixture viscosity in Pa s.
    pub fn liquid_effective(
        &self,
        x: &[f64],
        t_k: f64,
        viscosity_pa_s: f64,
    ) -> Result<Vec<f64>, CalculationError> {
        self.check(x)?;
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
        (0..x.len())
            .map(|i| {
                let (mut weight, mut phi_m) = (0.0, 0.0);
                for (j, s) in self.species.iter().enumerate() {
                    if j != i {
                        weight += x[j];
                        phi_m += x[j] * s.association * s.molar_mass_g_mol;
                    }
                }
                if weight <= 0.0 {
                    return Err(CalculationError::InvalidInput(InputError::NotPositive(
                        "liquid diffusivity needs a solvent",
                    )));
                }
                Ok(wilke_chang(
                    t_k,
                    phi_m / weight,
                    viscosity_pa_s,
                    self.species[i].boiling_volume_cm3,
                ))
            })
            .collect()
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
    fn fuller_co2_nitrogen_agrees_with_measurement() {
        // CO2–N2 at 298.2 K and 1 atm: measured 0.167 cm²/s (Cussler,
        // Diffusion, 3rd ed. 2009, Table 5.1-1). Fuller's method has a
        // typical error of a few per cent (Poling et al. 2001, Table 11-2).
        let d = fuller_binary(298.2, 101_325.0, 44.0095, 28.0134, 26.7, 18.5);
        assert!((d / 0.167e-4 - 1.0).abs() < 0.05, "{d}");
        // Hand evaluation of eq. 11-4.4 for the same inputs.
        let m_ab: f64 = 2.0 / (1.0 / 44.0095 + 1.0 / 28.0134);
        let s: f64 = 26.7_f64.cbrt() + 18.5_f64.cbrt();
        let hand = 0.00143 * 298.2_f64.powf(1.75) / (1.01325 * m_ab.sqrt() * s * s);
        assert!((d * 1e4 / hand - 1.0).abs() < 1e-14);
        // D ∝ 1/p.
        assert!(
            (fuller_binary(298.2, 2.0 * 101_325.0, 44.0095, 28.0134, 26.7, 18.5) * 2.0 / d - 1.0)
                .abs()
                < 1e-14
        );
    }

    #[test]
    fn wilke_chang_co2_in_water_is_near_measurement() {
        // CO2 in water at 298 K: measured about 1.9e-9 m²/s (Cussler 2009,
        // Table 5.2-1); Wilke–Chang's typical error is 10–20 %.
        let model = DiffusivityModel::new(&["carbon-dioxide", "water"]).unwrap();
        let d = model
            .liquid_effective(&[1e-4, 1.0 - 1e-4], 298.15, 0.89e-3)
            .unwrap()[0];
        assert!((d / 1.92e-9 - 1.0).abs() < 0.2, "{d}");
        // A pure liquid has no solvent for its only component.
        assert!(
            model
                .liquid_effective(&[0.0, 1.0], 298.15, 0.89e-3)
                .is_err()
        );
        assert!((tyn_calus_volume(94.0) - 0.285 * 94.0_f64.powf(1.048)).abs() < 1e-12);
    }

    #[test]
    fn blanc_reduces_to_the_binary_and_domain_is_explicit() {
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
