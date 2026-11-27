//! Phase transport properties (#225): dynamic viscosity, thermal conductivity
//! and surface tension. Side properties only: nothing here feeds back into
//! equilibrium or caloric results.
//!
//! * **Gas and dense fluid** (vapor and supercritical phases): Chung, Ajlan,
//!   Lee and Starling (1988), Ind. Eng. Chem. Res. 27, 671–679, in the form of
//!   Poling, Prausnitz and O'Connell (2001), eqs. 9-4.10, 9-5.24–9-5.44,
//!   9-6.18–9-6.21, 10-3.14 and 10-5.5–10-5.7. The phase's own EOS molar
//!   volume is the density input and NASA7 supplies the ideal-gas cv.
//! * **Liquid**: ChemSep 8.3 saturated-liquid correlations of each component
//!   with explicitly labelled mixing rules. The pressure effect on the liquid
//!   is not modelled; the declared pressure bound limits it.
//!
//! Outside the declared domain a property is `OutsideValidityDomain`, never an
//! extrapolated value. See `docs/transport-properties.md`.
use thermo_data::{Component, ComponentDatabase, TransportComponent, TransportDatabase};
use thermo_types::{
    CalculationError, DataMetadata, InputError, MOLAR_GAS_CONSTANT as R, ModelMetadata,
    Temperature, ValidityStatus,
};

/// Implementation version of the transport formulations; bump on any change.
pub const TRANSPORT_MODEL_VERSION: &str = "transport-v1";
/// Declared temperature domain of the gas method, K. The lower bound is the
/// NASA7 ideal-gas cv range; the upper bound is the verified range.
pub const GAS_TRANSPORT_TEMPERATURE_RANGE_K: [f64; 2] = [200.0, 600.0];
/// Largest reduced density ρ·Vc,m of the gas method for nonpolar mixtures.
pub const GAS_TRANSPORT_MAX_REDUCED_DENSITY: f64 = 1.5;
/// Largest reduced density of the gas viscosity when a polar component is present.
pub const POLAR_GAS_TRANSPORT_MAX_REDUCED_DENSITY: f64 = 0.1;
/// Components with a dipole of at least this many debye, or a nonzero Chung
/// association factor, are polar for the domain rules.
pub const POLAR_DIPOLE_THRESHOLD_DEBYE: f64 = 0.5;
/// Largest T/Tc of every present component for the liquid correlations.
pub const LIQUID_TRANSPORT_MAX_REDUCED_TEMPERATURE: f64 = 0.75;
/// Largest pressure of the liquid viscosity and thermal conductivity, Pa
/// (saturated-liquid correlations; compression is not modelled).
pub const LIQUID_TRANSPORT_MAX_PRESSURE_PA: f64 = 5.0e6;
/// Reduced density separating the reported dilute-gas and dense-fluid regimes.
pub const DENSE_FLUID_REDUCED_DENSITY: f64 = 0.1;

/// Chung viscosity coefficients a, b, c, d of E_i = a + b ω + c μr⁴ + d κ
/// (Chung et al. 1988, Table II; Poling et al. 2001, Table 9-6).
const VISCOSITY_E: [[f64; 4]; 10] = [
    [6.324, 50.412, -51.680, 1189.0],
    [1.210e-3, -1.154e-3, -6.257e-3, 0.03728],
    [5.283, 254.209, -168.48, 3898.0],
    [6.623, 38.096, -8.464, 31.42],
    [19.745, 7.630, -14.354, 31.53],
    [-1.900, -12.537, 4.985, -18.15],
    [24.275, 3.450, -11.291, 69.35],
    [0.7972, 1.117, 0.01235, -4.117],
    [-0.2382, 0.06770, -0.8163, 4.025],
    [0.06863, 0.3479, 0.5926, -0.727],
];
/// Chung thermal-conductivity coefficients of B_i (Poling et al. 2001, Table 10-5).
const CONDUCTIVITY_B: [[f64; 4]; 7] = [
    [2.4166, 0.74824, -0.91858, 121.72],
    [-0.50924, -1.5094, -49.991, 69.983],
    [6.6107, 5.6207, 64.760, 27.039],
    [14.543, -8.9139, -5.6379, 74.344],
    [0.79274, 0.82019, -0.69369, 6.3173],
    [-5.8634, 12.801, 9.5893, 65.529],
    [91.089, 128.11, -54.217, 523.81],
];

/// Which method a phase result used.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransportRegime {
    /// Chung method, reduced density at most [`DENSE_FLUID_REDUCED_DENSITY`].
    DiluteGas,
    /// Chung method above that density (dense gas or supercritical fluid).
    DenseFluid,
    /// ChemSep liquid correlations.
    Liquid,
}

impl TransportRegime {
    pub const fn label(self) -> &'static str {
        match self {
            Self::DiluteGas => "gas",
            Self::DenseFluid => "dense fluid",
            Self::Liquid => "liquid",
        }
    }
}

/// Transport properties of one phase. Each property is a value or the reason
/// it is unavailable (`OutsideValidityDomain`, `MissingParameter`,
/// `UnsupportedCapability`). Units: Pa s, W/(m K), N/m.
#[derive(Debug, Clone, PartialEq)]
pub struct TransportProperties {
    pub regime: TransportRegime,
    pub viscosity: Result<f64, CalculationError>,
    pub thermal_conductivity: Result<f64, CalculationError>,
    /// Liquid phases only (against the coexisting vapor).
    pub surface_tension: Result<f64, CalculationError>,
    /// `WithinDeclaredDomain` for one present component; `Unassessed` for
    /// mixtures, whose mixing rules are implementation-verified only.
    pub validity: ValidityStatus,
    pub model: ModelMetadata,
    pub data: DataMetadata,
}

#[derive(Debug, Clone)]
struct Pure {
    component: Component,
    transport: TransportComponent,
}

impl Pure {
    fn polar(&self) -> bool {
        self.transport.dipole_moment_debye() >= POLAR_DIPOLE_THRESHOLD_DEBYE
            || self.transport.association_factor() > 0.0
    }
}

/// Transport model for an ordered component selection.
#[derive(Debug, Clone)]
pub struct TransportModel {
    components: Vec<Pure>,
    data: DataMetadata,
}

/// Chung mixture parameters (Chung units: K, cm³/mol, g/mol, debye).
#[derive(Debug, Clone, Copy)]
struct ChungMixture {
    tc: f64,
    vc: f64,
    omega: f64,
    molar_mass: f64,
    reduced_dipole: f64,
    association: f64,
}

impl ChungMixture {
    fn coefficient(&self, row: &[f64; 4]) -> f64 {
        row[0]
            + row[1] * self.omega
            + row[2] * self.reduced_dipole.powi(4)
            + row[3] * self.association
    }
    fn fc(&self) -> f64 {
        1.0 - 0.2756 * self.omega + 0.059035 * self.reduced_dipole.powi(4) + self.association
    }
}

/// Neufeld et al. (1972) collision integral, Poling eq. 9-4.3 (six parameters).
fn collision_integral(t_star: f64) -> f64 {
    1.16145 * t_star.powf(-0.14874)
        + 0.52487 * (-0.77320 * t_star).exp()
        + 2.16178 * (-2.43787 * t_star).exp()
}

/// Chung's G2 group for coefficients `c` = (E or B)_1..5 and y = ρVc/6.
fn chung_g2(c: [f64; 5], y: f64) -> f64 {
    let g1 = (1.0 - 0.5 * y) / (1.0 - y).powi(3);
    (c[0] * (-(-c[3] * y).exp_m1()) / y + c[1] * g1 * (c[4] * y).exp() + c[2] * g1)
        / (c[0] * c[3] + c[1] + c[2])
}

impl TransportModel {
    /// Model for `ids` in package order. Every component needs a transport
    /// record and a component-database record, otherwise `MissingParameter`.
    pub fn new(ids: &[&str]) -> Result<Self, CalculationError> {
        if ids.is_empty() {
            return Err(CalculationError::InvalidInput(InputError::EmptyComposition));
        }
        let database = ComponentDatabase::embedded()?;
        let transport = TransportDatabase::embedded()?;
        let components = ids
            .iter()
            .map(|id| {
                Ok(Pure {
                    component: database.get(id)?.clone(),
                    transport: transport.get(id)?.clone(),
                })
            })
            .collect::<Result<Vec<_>, CalculationError>>()?;
        let base = database.metadata_for(ids.iter().copied());
        let t = transport.metadata();
        let data = DataMetadata {
            // The identifier names the transport extract; the version also
            // records the component-constant dataset version it combines with.
            identifier: t.identifier,
            version: format!("{}+components-{}", t.version, base.version),
            source: format!(
                "{}; Tc, pc, acentric factor, molar mass and NASA7 cv from {}",
                t.source, base.source
            ),
        };
        Ok(Self { components, data })
    }

    pub fn component_ids(&self) -> Vec<&str> {
        self.components.iter().map(|c| c.component.id()).collect()
    }

    pub fn data(&self) -> &DataMetadata {
        &self.data
    }

    /// Method identity of a regime.
    pub fn metadata(regime: TransportRegime) -> ModelMetadata {
        match regime {
            TransportRegime::DiluteGas | TransportRegime::DenseFluid => ModelMetadata {
                identifier: "chung-1988".into(),
                implementation_version: TRANSPORT_MODEL_VERSION.into(),
                formulation: "Chung et al. (1988) dense-fluid viscosity and thermal conductivity; six-parameter Neufeld collision integral; Chung mixing rules with xi_ij = zeta_ij = 1 and kappa_ij = sqrt(kappa_i kappa_j); phase EOS molar volume; NASA7 ideal-gas cv; no critical enhancement".into(),
            },
            TransportRegime::Liquid => ModelMetadata {
                identifier: "chemsep-liquid-correlations".into(),
                implementation_version: TRANSPORT_MODEL_VERSION.into(),
                formulation: "ChemSep 8.3 saturated-liquid correlations at T (pressure effect not modelled); viscosity ln(mu) = sum x_i ln(mu_i) (Grunberg-Nissan with G_ij = 0, labelled approximation); thermal conductivity DIPPR 9H (sum w_i k_i^-2)^-1/2 in mass fractions; surface tension sum x_i sigma_i (labelled approximation)".into(),
            },
        }
    }

    fn check_fractions(&self, x: &[f64]) -> Result<(), CalculationError> {
        if x.len() != self.components.len() {
            return Err(CalculationError::InvalidInput(
                InputError::ComponentMismatch,
            ));
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(CalculationError::InvalidInput(InputError::NonFinite(
                "transport mole fraction",
            )));
        }
        if x.iter().any(|v| *v < 0.0) {
            return Err(CalculationError::InvalidInput(InputError::Negative(
                "transport mole fraction",
            )));
        }
        if (x.iter().sum::<f64>() - 1.0).abs() > 1e-8 {
            return Err(CalculationError::InvalidInput(InputError::NotNormalized));
        }
        Ok(())
    }

    fn present<'a>(&'a self, x: &'a [f64]) -> impl Iterator<Item = (&'a Pure, f64)> + 'a {
        self.components
            .iter()
            .zip(x.iter().copied())
            .filter(|(_, xi)| *xi > 0.0)
    }

    fn validity(x: &[f64]) -> ValidityStatus {
        if x.iter().filter(|v| **v > 0.0).count() == 1 {
            ValidityStatus::WithinDeclaredDomain
        } else {
            ValidityStatus::Unassessed
        }
    }

    /// Chung mixing rules (Poling eqs. 9-5.25–9-5.44).
    fn chung_mixture(&self, x: &[f64]) -> ChungMixture {
        let p: Vec<_> = self
            .components
            .iter()
            .map(|c| {
                let vc = c.transport.critical_volume() * 1e6;
                let m = c.component.molar_mass().value() * 1e3;
                (
                    0.809 * vc.cbrt(),
                    c.component.critical_temperature() / 1.2593,
                    c.component.acentric_factor(),
                    m,
                    c.transport.dipole_moment_debye(),
                    c.transport.association_factor(),
                )
            })
            .collect();
        let (mut s3, mut e, mut w, mut msum, mut mu4, mut kappa) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for (i, a) in p.iter().enumerate() {
            for (j, b) in p.iter().enumerate() {
                let xx = x[i] * x[j];
                if xx == 0.0 {
                    continue;
                }
                let sij = (a.0 * b.0).sqrt();
                let sij3 = sij.powi(3);
                let eij = (a.1 * b.1).sqrt();
                let mij = 2.0 * a.3 * b.3 / (a.3 + b.3);
                s3 += xx * sij3;
                e += xx * eij * sij3;
                w += xx * 0.5 * (a.2 + b.2) * sij3;
                msum += xx * eij * sij * sij * mij.sqrt();
                mu4 += xx * a.4 * a.4 * b.4 * b.4 / sij3;
                kappa += xx * (a.5 * b.5).sqrt();
            }
        }
        let e = e / s3;
        let molar_mass = (msum / (e * s3.powf(2.0 / 3.0))).powi(2);
        let vc = s3 / 0.809f64.powi(3);
        let tc = 1.2593 * e;
        let dipole = (s3 * mu4).powf(0.25);
        ChungMixture {
            tc,
            vc,
            omega: w / s3,
            molar_mass,
            reduced_dipole: 131.3 * dipole / (vc * tc).sqrt(),
            association: kappa,
        }
    }

    /// Gas or supercritical phase at `temperature` with the phase's molar
    /// volume (m³/mol) and mole fractions `x` in component order.
    pub fn gas(
        &self,
        x: &[f64],
        temperature: Temperature,
        molar_volume: f64,
    ) -> Result<TransportProperties, CalculationError> {
        self.check_fractions(x)?;
        if !molar_volume.is_finite() || molar_volume <= 0.0 {
            return Err(CalculationError::InvalidInput(InputError::NotPositive(
                "transport molar volume",
            )));
        }
        let t = temperature.value();
        let mix = self.chung_mixture(x);
        // ρ in mol/cm³ times Vc in cm³/mol.
        let reduced_density = mix.vc * 1e-6 / molar_volume;
        let regime = if reduced_density > DENSE_FLUID_REDUCED_DENSITY {
            TransportRegime::DenseFluid
        } else {
            TransportRegime::DiluteGas
        };
        let polar: Vec<&str> = self
            .present(x)
            .filter(|(c, _)| c.polar())
            .map(|(c, _)| c.component.id())
            .collect();
        let [lo, hi] = GAS_TRANSPORT_TEMPERATURE_RANGE_K;
        let t_star = 1.2593 * t / mix.tc;
        let common = if !(lo..=hi).contains(&t) {
            Err(CalculationError::OutsideValidityDomain(format!(
                "Chung gas transport: temperature {t} K outside {lo}..={hi} K"
            )))
        } else if reduced_density > GAS_TRANSPORT_MAX_REDUCED_DENSITY {
            Err(CalculationError::OutsideValidityDomain(format!(
                "Chung gas transport: reduced density {reduced_density:.4} above {GAS_TRANSPORT_MAX_REDUCED_DENSITY}"
            )))
        } else if !(0.3..=100.0).contains(&t_star) {
            Err(CalculationError::OutsideValidityDomain(format!(
                "Chung gas transport: T* = {t_star:.4} outside the collision-integral range 0.3..=100"
            )))
        } else {
            Ok(())
        };
        let viscosity = common.clone().and_then(|()| {
            if !polar.is_empty() && reduced_density > POLAR_GAS_TRANSPORT_MAX_REDUCED_DENSITY {
                return Err(CalculationError::OutsideValidityDomain(format!(
                    "Chung gas viscosity with polar {}: reduced density {reduced_density:.4} above {POLAR_GAS_TRANSPORT_MAX_REDUCED_DENSITY}",
                    polar.join(", ")
                )));
            }
            finite_positive(chung_viscosity(&mix, t, reduced_density), "Chung viscosity")
        });
        let thermal_conductivity = common.and_then(|()| {
            if !polar.is_empty() {
                return Err(CalculationError::OutsideValidityDomain(format!(
                    "Chung gas thermal conductivity is not qualified for polar components ({})",
                    polar.join(", ")
                )));
            }
            let mut cv = 0.0;
            for (c, xi) in self.present(x) {
                cv += xi * (c.component.ideal_gas(temperature)?.heat_capacity - R);
            }
            finite_positive(
                chung_conductivity(&mix, t, reduced_density, cv),
                "Chung thermal conductivity",
            )
        });
        Ok(TransportProperties {
            regime,
            viscosity,
            thermal_conductivity,
            surface_tension: Err(CalculationError::UnsupportedCapability(
                "surface tension is reported for liquid phases only".into(),
            )),
            validity: Self::validity(x),
            model: Self::metadata(regime),
            data: self.data.clone(),
        })
    }

    /// Liquid phase at `temperature` and `pressure_pa` with mole fractions `x`.
    pub fn liquid(
        &self,
        x: &[f64],
        temperature: Temperature,
        pressure_pa: f64,
    ) -> Result<TransportProperties, CalculationError> {
        self.check_fractions(x)?;
        if !pressure_pa.is_finite() || pressure_pa <= 0.0 {
            return Err(CalculationError::InvalidInput(InputError::NotPositive(
                "transport pressure",
            )));
        }
        let t = temperature.value();
        let reduced = self.present(x).try_for_each(|(c, _)| {
            let tc = c.component.critical_temperature();
            if t > LIQUID_TRANSPORT_MAX_REDUCED_TEMPERATURE * tc {
                Err(CalculationError::OutsideValidityDomain(format!(
                    "liquid transport: T = {t} K above {LIQUID_TRANSPORT_MAX_REDUCED_TEMPERATURE} Tc of {} ({tc} K)",
                    c.component.id()
                )))
            } else {
                Ok(())
            }
        });
        let pressure = reduced.clone().and_then(|()| {
            if pressure_pa > LIQUID_TRANSPORT_MAX_PRESSURE_PA {
                Err(CalculationError::OutsideValidityDomain(format!(
                    "liquid transport: pressure {pressure_pa} Pa above {LIQUID_TRANSPORT_MAX_PRESSURE_PA} Pa (compression not modelled)"
                )))
            } else {
                Ok(())
            }
        });
        let viscosity = pressure.clone().and_then(|()| {
            let mut ln_mu = 0.0;
            for (c, xi) in self.present(x) {
                ln_mu += xi * c.transport.liquid_viscosity_at(temperature)?.ln();
            }
            finite_positive(ln_mu.exp(), "liquid viscosity")
        });
        let thermal_conductivity = pressure.and_then(|()| {
            let mass: f64 = self
                .present(x)
                .map(|(c, xi)| xi * c.component.molar_mass().value())
                .sum();
            let mut sum = 0.0;
            for (c, xi) in self.present(x) {
                let k = c.transport.liquid_thermal_conductivity_at(temperature)?;
                sum += xi * c.component.molar_mass().value() / mass / (k * k);
            }
            finite_positive(sum.sqrt().recip(), "liquid thermal conductivity")
        });
        let surface_tension = reduced.and_then(|()| {
            let mut sigma = 0.0;
            for (c, xi) in self.present(x) {
                sigma += xi * c.transport.surface_tension_at(temperature)?;
            }
            finite_positive(sigma, "surface tension")
        });
        Ok(TransportProperties {
            regime: TransportRegime::Liquid,
            viscosity,
            thermal_conductivity,
            surface_tension,
            validity: Self::validity(x),
            model: Self::metadata(TransportRegime::Liquid),
            data: self.data.clone(),
        })
    }
}

fn finite_positive(v: f64, what: &str) -> Result<f64, CalculationError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(CalculationError::NumericalFailure(format!(
            "{what} is not finite and positive"
        )))
    }
}

/// Chung dense-fluid viscosity, Pa s (Poling eqs. 9-6.18–9-6.21).
fn chung_viscosity(mix: &ChungMixture, t: f64, reduced_density: f64) -> f64 {
    let e: Vec<f64> = VISCOSITY_E.iter().map(|r| mix.coefficient(r)).collect();
    let t_star = 1.2593 * t / mix.tc;
    let y = reduced_density / 6.0;
    let g2 = chung_g2([e[0], e[1], e[2], e[3], e[4]], y);
    let eta_kk = e[6] * y * y * g2 * (e[7] + e[8] / t_star + e[9] / (t_star * t_star)).exp();
    let eta_k = t_star.sqrt() / collision_integral(t_star) * mix.fc() * (1.0 / g2 + e[5] * y);
    // μP to Pa s.
    (eta_k + eta_kk) * 36.344 * (mix.molar_mass * mix.tc).sqrt() / mix.vc.powf(2.0 / 3.0) * 1e-7
}

/// Chung low-pressure viscosity, Pa s (Poling eq. 9-4.10).
fn chung_dilute_viscosity(mix: &ChungMixture, t: f64) -> f64 {
    let t_star = 1.2593 * t / mix.tc;
    40.785 * mix.fc() * (mix.molar_mass * t).sqrt()
        / (mix.vc.powf(2.0 / 3.0) * collision_integral(t_star))
        * 1e-7
}

/// Chung dense-fluid thermal conductivity, W/(m K) (Poling eqs. 10-3.14,
/// 10-5.5–10-5.7); `cv` is the ideal-gas molar isochoric heat capacity.
// 0.6366 is Chung's published coefficient of Ψ, not an approximation of 2/π.
#[allow(clippy::approx_constant)]
fn chung_conductivity(mix: &ChungMixture, t: f64, reduced_density: f64, cv: f64) -> f64 {
    let b: Vec<f64> = CONDUCTIVITY_B.iter().map(|r| mix.coefficient(r)).collect();
    let tr = t / mix.tc;
    let alpha = cv / R - 1.5;
    let beta = 0.7862 - 0.7109 * mix.omega + 1.3168 * mix.omega * mix.omega;
    let z = 2.0 + 10.5 * tr * tr;
    let psi = 1.0
        + alpha * (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z)
            / (0.6366 + beta * z + 1.061 * alpha * beta);
    let y = reduced_density / 6.0;
    let g2 = chung_g2([b[0], b[1], b[2], b[3], b[4]], y);
    let m = mix.molar_mass * 1e-3;
    let q = 3.586e-3 * (mix.tc / m).sqrt() / mix.vc.powf(2.0 / 3.0);
    31.2 * chung_dilute_viscosity(mix, t) * psi / m * (1.0 / g2 + b[5] * y)
        + q * b[6] * y * y * tr.sqrt() * g2
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(k: f64) -> Temperature {
        Temperature::new(k).unwrap()
    }

    #[test]
    fn dilute_limit_of_the_dense_forms_is_the_low_pressure_method() {
        let model = TransportModel::new(&["methane", "nitrogen"]).unwrap();
        let mix = model.chung_mixture(&[0.7, 0.3]);
        let dilute = chung_dilute_viscosity(&mix, 300.0);
        let dense = chung_viscosity(&mix, 300.0, 1e-9);
        // The published constants agree to their rounding:
        // 36.344 sqrt(1.2593) = 40.7853 against 40.785 (7e-6 relative).
        assert!((dense / dilute - 1.0).abs() < 1e-5, "{dense} {dilute}");
    }

    #[test]
    fn pure_mixing_rules_return_the_pure_constants() {
        let model = TransportModel::new(&["water", "methane"]).unwrap();
        let mix = model.chung_mixture(&[1.0, 0.0]);
        let db = TransportDatabase::embedded().unwrap();
        let water = db.get("water").unwrap();
        assert!((mix.vc / (water.critical_volume() * 1e6) - 1.0).abs() < 1e-14);
        assert!((mix.tc / 647.096 - 1.0).abs() < 1e-3);
        assert!((mix.molar_mass / 18.015 - 1.0).abs() < 1e-3);
        assert!((mix.association - 0.075908).abs() < 1e-15);
    }

    #[test]
    fn domain_violations_are_reported_not_extrapolated() {
        let model = TransportModel::new(&["methane", "water"]).unwrap();
        let gas = model.gas(&[1.0, 0.0], t(650.0), 0.02).unwrap();
        assert!(matches!(
            gas.viscosity,
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        let wet = model.gas(&[0.9, 0.1], t(400.0), 0.03).unwrap();
        assert!(wet.viscosity.is_ok());
        assert!(matches!(
            wet.thermal_conductivity,
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert_eq!(wet.validity, ValidityStatus::Unassessed);
        let liquid = model.liquid(&[0.0, 1.0], t(300.0), 2.0e7).unwrap();
        assert!(matches!(
            liquid.viscosity,
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(liquid.surface_tension.is_ok());
        assert_eq!(liquid.validity, ValidityStatus::WithinDeclaredDomain);
        // A trace of methane above 0.75 Tc makes the liquid mixture unavailable.
        let traced = model.liquid(&[1e-6, 1.0 - 1e-6], t(300.0), 1e5).unwrap();
        assert!(matches!(
            traced.viscosity,
            Err(CalculationError::OutsideValidityDomain(_))
        ));
        assert!(matches!(
            TransportModel::new(&["methane", "hydrogen"]),
            Err(CalculationError::MissingParameter(_))
        ));
        assert!(matches!(
            model.gas(&[0.5, 0.4], t(300.0), 0.02),
            Err(CalculationError::InvalidInput(_))
        ));
    }
}
