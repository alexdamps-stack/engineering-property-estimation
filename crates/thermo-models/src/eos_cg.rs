//! EOS-CG-2021 multiparameter Helmholtz mixture model.

use std::ops::{Add, Mul, Sub};
use std::sync::OnceLock;

use serde::Deserialize;
use thermo_types::{CalculationError as Error, MOLAR_GAS_CONSTANT as R};

#[derive(Debug, Clone, Copy, Default)]
struct Jet {
    value: f64,
    d: f64,
    t: f64,
    dd: f64,
    dt: f64,
    tt: f64,
}

impl Jet {
    fn constant(value: f64) -> Self {
        Self {
            value,
            ..Self::default()
        }
    }

    fn map(self, value: f64, slope: f64, curvature: f64) -> Self {
        Self {
            value,
            d: slope * self.d,
            t: slope * self.t,
            dd: slope * self.dd + curvature * self.d * self.d,
            dt: slope * self.dt + curvature * self.d * self.t,
            tt: slope * self.tt + curvature * self.t * self.t,
        }
    }

    fn exp(self) -> Self {
        let value = self.value.exp();
        self.map(value, value, value)
    }

    fn ln(self) -> Self {
        let inverse = self.value.recip();
        self.map(self.value.ln(), inverse, -inverse * inverse)
    }

    fn pow(self, power: f64) -> Self {
        if power == 0.0 {
            return Self::constant(1.0);
        }
        if power == 1.0 {
            return self;
        }
        let value = self.value.powf(power);
        let slope = power * self.value.powf(power - 1.0);
        self.map(
            value,
            slope,
            power * (power - 1.0) * self.value.powf(power - 2.0),
        )
    }

    fn abs_pow(self, power: f64) -> Self {
        if self.value >= 0.0 {
            self.pow(power)
        } else {
            (self * -1.0).pow(power)
        }
    }
}

impl Add for Jet {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            d: self.d + rhs.d,
            t: self.t + rhs.t,
            dd: self.dd + rhs.dd,
            dt: self.dt + rhs.dt,
            tt: self.tt + rhs.tt,
        }
    }
}

impl Sub for Jet {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + rhs * -1.0
    }
}

impl Mul<f64> for Jet {
    type Output = Self;

    fn mul(self, factor: f64) -> Self {
        self.map(self.value * factor, factor, 0.0)
    }
}

impl Mul for Jet {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            d: self.d * rhs.value + self.value * rhs.d,
            t: self.t * rhs.value + self.value * rhs.t,
            dd: self.dd * rhs.value + 2.0 * self.d * rhs.d + self.value * rhs.dd,
            dt: self.dt * rhs.value + self.d * rhs.t + self.t * rhs.d + self.value * rhs.dt,
            tt: self.tt * rhs.value + 2.0 * self.t * rhs.t + self.value * rhs.tt,
        }
    }
}

#[derive(Deserialize)]
#[serde(tag = "type")]
enum Ideal {
    #[serde(rename = "IdealGasHelmholtzLead")]
    Lead { a1: f64, a2: f64 },
    #[serde(rename = "IdealGasHelmholtzEnthalpyEntropyOffset")]
    Offset { a1: f64, a2: f64 },
    #[serde(rename = "IdealGasHelmholtzLogTau")]
    LogTau { a: f64 },
    #[serde(rename = "IdealGasHelmholtzPower")]
    Power { n: Vec<f64>, t: Vec<f64> },
    #[serde(rename = "IdealGasHelmholtzCP0PolyT")]
    HeatCapacity {
        c: Vec<f64>,
        t: Vec<f64>,
        #[serde(rename = "T0")]
        reference_temperature: f64,
        #[serde(rename = "Tc")]
        reducing_temperature: f64,
    },
    #[serde(rename = "IdealGasHelmholtzPlanckEinstein")]
    PlanckEinstein { n: Vec<f64>, t: Vec<f64> },
    #[serde(rename = "IdealGasHelmholtzPlanckEinsteinFunctionT")]
    PlanckEinsteinKelvin {
        n: Vec<f64>,
        v: Vec<f64>,
        #[serde(rename = "Tcrit")]
        critical_temperature: f64,
    },
}

impl Ideal {
    fn evaluate(&self, delta: Jet, tau: Jet) -> Jet {
        let planck_einstein =
            |n: f64, theta: f64| (Jet::constant(1.0) - (tau * -theta).exp()).ln() * n;
        match self {
            Self::Lead { a1, a2 } => delta.ln() + Jet::constant(*a1) + tau * *a2,
            Self::Offset { a1, a2 } => Jet::constant(*a1) + tau * *a2,
            Self::LogTau { a } => tau.ln() * *a,
            Self::Power { n, t } => n
                .iter()
                .zip(t)
                .fold(Jet::default(), |sum, (&n, &t)| sum + tau.pow(t) * n),
            Self::HeatCapacity {
                c,
                t,
                reference_temperature,
                reducing_temperature,
            } => {
                let tau0 = reducing_temperature / reference_temperature;
                let reference = Jet::constant(tau0);
                c.iter()
                    .zip(t)
                    .fold(Jet::default(), |sum, (&coefficient, &power)| {
                        let integral = if power == 0.0 {
                            (tau * tau0.recip()).ln() - (tau - reference) * tau0.recip()
                        } else if power == -1.0 {
                            tau - reference - tau * (tau * tau0.recip()).ln()
                        } else {
                            let exponent = -power;
                            (tau.pow(exponent)
                                - Jet::constant(tau0.powf(exponent))
                                - (tau - reference) * (exponent * tau0.powf(exponent - 1.0)))
                                * (-reducing_temperature.powf(power) / (power * (power + 1.0)))
                        };
                        sum + integral
                            * if power == -1.0 {
                                coefficient / reducing_temperature
                            } else {
                                coefficient
                            }
                    })
            }
            Self::PlanckEinstein { n, t } => n
                .iter()
                .zip(t)
                .fold(Jet::default(), |sum, (&n, &t)| sum + planck_einstein(n, t)),
            Self::PlanckEinsteinKelvin {
                n,
                v,
                critical_temperature,
            } => n.iter().zip(v).fold(Jet::default(), |sum, (&n, &v)| {
                sum + planck_einstein(n, v / critical_temperature)
            }),
        }
    }
}

#[derive(Deserialize)]
#[serde(tag = "type")]
enum Residual {
    #[serde(rename = "ResidualHelmholtzPower")]
    Power {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        l: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzGaussian")]
    Gaussian {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        eta: Vec<f64>,
        epsilon: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzGERG2008")]
    Gerg {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        eta: Vec<f64>,
        epsilon: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
    },
    #[serde(rename = "Gaussian+Exponential")]
    Hybrid {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        l: Vec<f64>,
        eta: Vec<f64>,
        epsilon: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzGaoB")]
    GaoB {
        n: Vec<f64>,
        d: Vec<f64>,
        t: Vec<f64>,
        eta: Vec<f64>,
        epsilon: Vec<f64>,
        beta: Vec<f64>,
        gamma: Vec<f64>,
        b: Vec<f64>,
    },
    #[serde(rename = "ResidualHelmholtzNonAnalytic")]
    NonAnalytic {
        n: Vec<f64>,
        a: Vec<f64>,
        b: Vec<f64>,
        beta: Vec<f64>,
        #[serde(rename = "A")]
        theta: Vec<f64>,
        #[serde(rename = "B")]
        density: Vec<f64>,
        #[serde(rename = "C")]
        density_decay: Vec<f64>,
        #[serde(rename = "D")]
        temperature_decay: Vec<f64>,
    },
}

impl Residual {
    fn evaluate(&self, delta: Jet, tau: Jet) -> Jet {
        let mut sum = Jet::default();
        match self {
            Self::Power { n, d, t, l } => {
                for i in 0..n.len() {
                    let term = delta.pow(d[i]) * tau.pow(t[i]) * n[i];
                    sum = sum
                        + if l[i] == 0.0 {
                            term
                        } else {
                            term * (delta.pow(l[i]) * -1.0).exp()
                        };
                }
            }
            Self::Gaussian {
                n,
                d,
                t,
                eta,
                epsilon,
                beta,
                gamma,
            } => {
                for i in 0..n.len() {
                    let exponent = (delta - Jet::constant(epsilon[i])).pow(2.0) * -eta[i]
                        + (tau - Jet::constant(gamma[i])).pow(2.0) * -beta[i];
                    sum = sum + delta.pow(d[i]) * tau.pow(t[i]) * exponent.exp() * n[i];
                }
            }
            Self::NonAnalytic {
                n,
                a,
                b,
                beta,
                theta,
                density,
                density_decay,
                temperature_decay,
            } => {
                let dr = delta - Jet::constant(1.0);
                let tr = tau - Jet::constant(1.0);
                for i in 0..n.len() {
                    let theta = dr.abs_pow(beta[i].recip()) * theta[i] - tr;
                    let shape = theta.pow(2.0) + dr.abs_pow(2.0 * a[i]) * density[i];
                    let decay = (dr.pow(2.0) * -density_decay[i]
                        + tr.pow(2.0) * -temperature_decay[i])
                        .exp();
                    sum = sum + delta * shape.pow(b[i]) * decay * n[i];
                }
            }
            Self::Gerg {
                n,
                d,
                t,
                eta,
                epsilon,
                beta,
                gamma,
            } => {
                for i in 0..n.len() {
                    let exponent = (delta - Jet::constant(epsilon[i])).pow(2.0) * -eta[i]
                        + (delta - Jet::constant(gamma[i])) * -beta[i];
                    sum = sum + delta.pow(d[i]) * tau.pow(t[i]) * exponent.exp() * n[i];
                }
            }
            Self::Hybrid {
                n,
                d,
                t,
                l,
                eta,
                epsilon,
                beta,
                gamma,
            } => {
                for i in 0..n.len() {
                    let exponent = if l[i] > 0.0 {
                        delta.pow(l[i]) * -1.0
                    } else {
                        (delta - Jet::constant(epsilon[i])).pow(2.0) * -eta[i]
                            + (tau - Jet::constant(gamma[i])).pow(2.0) * -beta[i]
                    };
                    sum = sum + delta.pow(d[i]) * tau.pow(t[i]) * exponent.exp() * n[i];
                }
            }
            Self::GaoB {
                n,
                d,
                t,
                eta,
                epsilon,
                beta,
                gamma,
                b,
            } => {
                for i in 0..n.len() {
                    let exponent = (delta - Jet::constant(epsilon[i])).pow(2.0) * eta[i]
                        + ((tau - Jet::constant(gamma[i])).pow(2.0) * beta[i]
                            + Jet::constant(b[i]))
                        .pow(-1.0);
                    sum = sum + delta.pow(d[i]) * tau.pow(t[i]) * exponent.exp() * n[i];
                }
            }
        }
        sum
    }
}

#[derive(Deserialize)]
struct Component {
    critical_density: f64,
    critical_temperature: f64,
    molar_mass: f64,
    source_gas_constant: f64,
    ideal: Vec<Ideal>,
    residual: Vec<Residual>,
}

#[derive(Deserialize)]
struct Pair {
    first: usize,
    second: usize,
    beta_temperature: f64,
    gamma_temperature: f64,
    beta_volume: f64,
    gamma_volume: f64,
    departure_scale: f64,
    #[serde(default)]
    departure: Vec<Residual>,
}

#[derive(Deserialize)]
struct Parameters {
    components: [Component; 16],
    pairs: Vec<Pair>,
}

fn parameters() -> &'static Parameters {
    static PARAMETERS: OnceLock<Parameters> = OnceLock::new();
    PARAMETERS.get_or_init(|| {
        serde_json::from_str(include_str!("../../../data/eos-cg-2021.json"))
            .expect("embedded EOS-CG-2021 mixture parameters")
    })
}

/// A homogeneous point before any phase stability assessment. SI molar units.
#[derive(Debug, Clone)]
pub struct EosCg2021State {
    pub pressure: f64,
    pub compressibility: f64,
    pub alpha_residual: f64,
    pub mass_density: f64,
    pub enthalpy: f64,
    pub residual_enthalpy: f64,
    pub entropy: f64,
    pub cv: f64,
    pub cp: f64,
    pub speed_of_sound: f64,
    pub ln_fugacity_coefficients: Vec<f64>,
}

/// Native EOS-CG-2021 with fugacity coefficients returned in input order.
#[derive(Debug, Clone)]
pub struct EosCg2021 {
    fractions: [Jet; 16],
    indices: Vec<usize>,
}

impl EosCg2021 {
    pub const COMPONENTS: [&str; 16] = [
        "carbon-dioxide",
        "water",
        "nitrogen",
        "oxygen",
        "argon",
        "carbon-monoxide",
        "hydrogen",
        "methane",
        "hydrogen-sulfide",
        "sulfur-dioxide",
        "monoethanolamine",
        "diethanolamine",
        "hydrogen-chloride",
        "chlorine",
        "ammonia",
        "methyldiethanolamine",
    ];

    pub fn new(components: &[&str], fractions: &[f64]) -> Result<Self, Error> {
        if components.is_empty()
            || components.len() != fractions.len()
            || fractions.iter().any(|value| !(0.0..=1.0).contains(value))
            || (fractions.iter().sum::<f64>() - 1.0).abs() > 1e-10
        {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 requires normalized mole fractions in component order".into(),
            ));
        }
        let mut indices = Vec::with_capacity(components.len());
        let mut composition = [Jet::default(); 16];
        for (&name, &fraction) in components.iter().zip(fractions) {
            let index = Self::COMPONENTS
                .iter()
                .position(|id| *id == name)
                .ok_or_else(|| {
                    Error::UnsupportedCapability(format!(
                        "unsupported EOS-CG-2021 component {name}"
                    ))
                })?;
            if indices.contains(&index) {
                return Err(Error::OutsideValidityDomain(format!(
                    "duplicate EOS-CG-2021 component {name}"
                )));
            }
            indices.push(index);
            composition[index] = Jet::constant(fraction);
        }
        Ok(Self {
            fractions: composition,
            indices,
        })
    }

    fn residual(x: &[Jet; 16], density: Jet, inverse_temperature: Jet) -> Jet {
        let table = parameters();
        let mut reducing_temperature = Jet::default();
        let mut reducing_volume = Jet::default();
        for (fraction, component) in x.iter().zip(&table.components) {
            let square = *fraction * *fraction;
            reducing_temperature = reducing_temperature + square * component.critical_temperature;
            reducing_volume = reducing_volume + square * component.critical_density.recip();
        }
        for pair in &table.pairs {
            let first = x[pair.first];
            let second = x[pair.second];
            if first.value + second.value == 0.0 {
                continue;
            }
            let a = &table.components[pair.first];
            let b = &table.components[pair.second];
            let cross = |beta: f64, gamma: f64, combining: f64| {
                first
                    * second
                    * (first + second)
                    * (first * beta.powi(2) + second).pow(-1.0)
                    * (2.0 * beta * gamma * combining)
            };
            reducing_temperature = reducing_temperature
                + cross(
                    pair.beta_temperature,
                    pair.gamma_temperature,
                    (a.critical_temperature * b.critical_temperature).sqrt(),
                );
            reducing_volume = reducing_volume
                + cross(
                    pair.beta_volume,
                    pair.gamma_volume,
                    (0.5 * (a.critical_density.cbrt().recip() + b.critical_density.cbrt().recip()))
                        .powi(3),
                );
        }
        let delta = density * reducing_volume;
        let tau = inverse_temperature * reducing_temperature;
        let evaluate = |blocks: &[Residual]| {
            blocks.iter().fold(Jet::default(), |sum, block| {
                sum + block.evaluate(delta, tau)
            })
        };
        let present = |fraction: Jet| fraction.value != 0.0 || fraction.d != 0.0;
        let mut residual = Jet::default();
        for (&fraction, component) in x.iter().zip(&table.components) {
            if present(fraction) {
                residual = residual + fraction * evaluate(&component.residual);
            }
        }
        for pair in &table.pairs {
            let weight = x[pair.first] * x[pair.second] * pair.departure_scale;
            if present(weight) {
                residual = residual + weight * evaluate(&pair.departure);
            }
        }
        residual
    }

    fn residual_derivatives(
        &self,
        temperature: f64,
        density: f64,
    ) -> Result<(Jet, Jet, Jet), Error> {
        if !(250.0..=500.0).contains(&temperature) || !(0.001..=30_000.0).contains(&density) {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 homogeneous evaluation requires 250–500 K and 0.001–30000 mol/m³"
                    .into(),
            ));
        }
        let rho = Jet {
            value: density,
            d: density,
            dd: density,
            ..Jet::default()
        };
        let inv_t = temperature.recip();
        let tau = Jet {
            value: inv_t,
            t: inv_t,
            tt: inv_t,
            ..Jet::default()
        };
        let ar = Self::residual(&self.fractions, rho, tau);
        Ok((ar, rho, tau))
    }

    /// Mathematical isotherm, including its mechanically unstable branches.
    /// Returns (pressure in Pa, dp/dρ at fixed T, dimensionless residual Helmholtz energy).
    pub fn isotherm(&self, temperature: f64, density: f64) -> Result<(f64, f64, f64), Error> {
        let (ar, _, _) = self.residual_derivatives(temperature, density)?;
        let rt = R * temperature;
        let pressure = density * rt * (1.0 + ar.d);
        let slope = rt * (1.0 + ar.d + ar.dd);
        if [pressure, slope, ar.value].into_iter().all(f64::is_finite) {
            Ok((pressure, slope, ar.value))
        } else {
            Err(Error::NumericalFailure(
                "nonfinite EOS-CG-2021 isotherm".into(),
            ))
        }
    }

    pub fn evaluate(&self, temperature: f64, density: f64) -> Result<EosCg2021State, Error> {
        let (ar, rho, tau) = self.residual_derivatives(temperature, density)?;
        let inv_t = tau.value;
        let mut ideal = Jet::default();
        let mut molar_mass = 0.0;
        for (fraction, component) in self.fractions.iter().zip(&parameters().components) {
            let x = fraction.value;
            molar_mass += x * component.molar_mass;
            if x == 0.0 {
                continue;
            }
            let delta = rho * component.critical_density.recip();
            let tau = tau * component.critical_temperature;
            // Preserve the published ideal density coefficient of one and
            // scale only the temperature-dependent and reference terms.
            let density_term = delta.ln();
            let pure = component
                .ideal
                .iter()
                .fold(Jet::default(), |sum, term| sum + term.evaluate(delta, tau));
            ideal = ideal
                + (density_term
                    + (pure - density_term) * (component.source_gas_constant / R)
                    + Jet::constant(x.ln()))
                    * x;
        }
        let alpha = ideal + ar;
        let z = 1.0 + ar.d;
        let stiffness = 1.0 + ar.d + ar.dd;
        let coupling = z - ar.dt;
        let cv = R * (alpha.t - alpha.tt);
        if !(z > 0.0 && stiffness > 0.0 && cv > 0.0) {
            return Err(Error::OutsideValidityDomain(
                "EOS-CG-2021 homogeneous state is mechanically or calorically inadmissible".into(),
            ));
        }
        let cp = cv + R * coupling * coupling / stiffness;
        let derivatives: Vec<_> = self
            .indices
            .iter()
            .map(|&index| {
                let mut composition = self.fractions;
                composition[index].d = 1.0;
                Self::residual(&composition, Jet::constant(density), Jet::constant(inv_t)).d
            })
            .collect();
        let weighted_derivative: f64 = derivatives
            .iter()
            .zip(&self.indices)
            .map(|(derivative, &index)| derivative * self.fractions[index].value)
            .sum();
        let common = ar.value + ar.d - z.ln() - weighted_derivative;
        let state = EosCg2021State {
            pressure: density * R * temperature * z,
            compressibility: z,
            alpha_residual: ar.value,
            mass_density: density * molar_mass,
            enthalpy: R * temperature * (alpha.t + z),
            residual_enthalpy: R * temperature * (ar.t + ar.d),
            entropy: R * (alpha.t - alpha.value),
            cv,
            cp,
            speed_of_sound: (R * temperature * stiffness * cp / (molar_mass * cv)).sqrt(),
            ln_fugacity_coefficients: derivatives
                .into_iter()
                .map(|derivative| common + derivative)
                .collect(),
        };
        if [
            state.pressure,
            state.enthalpy,
            state.entropy,
            state.cp,
            state.cv,
            state.speed_of_sound,
        ]
        .into_iter()
        .chain(state.ln_fugacity_coefficients.iter().copied())
        .all(f64::is_finite)
        {
            Ok(state)
        } else {
            Err(Error::NumericalFailure(
                "nonfinite EOS-CG-2021 homogeneous result".into(),
            ))
        }
    }
}
