//! Validated SI inputs and explicit errors for homogeneous property calculations.

use serde::{Deserialize, Serialize};
use std::{error::Error, fmt};

/// Molar gas constant in J/(mol K), used by every implementation layer.
pub const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

#[derive(Debug, Clone, PartialEq)]
pub enum ThermoError {
    InvalidInput(String),
    MissingData(String),
    Unsupported(String),
    NumericalFailure(String),
}

impl fmt::Display for ThermoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let (kind, message) = match self {
            Self::InvalidInput(m) => ("invalid input", m),
            Self::MissingData(m) => ("missing data", m),
            Self::Unsupported(m) => ("unsupported", m),
            Self::NumericalFailure(m) => ("numerical failure", m),
        };
        write!(f, "{kind}: {message}")
    }
}

impl Error for ThermoError {}

macro_rules! positive_quantity {
    ($name:ident, $unit:literal) => {
        #[doc = concat!("Finite, strictly positive quantity in ", $unit, ".")]
        #[derive(Debug, Clone, Copy, PartialEq)]
        pub struct $name(f64);

        impl $name {
            pub fn new(value: f64) -> Result<Self, ThermoError> {
                if !value.is_finite() || value <= 0.0 {
                    return Err(ThermoError::InvalidInput(
                        concat!(stringify!($name), " must be finite and positive").into(),
                    ));
                }
                Ok(Self(value))
            }

            pub fn value(self) -> f64 {
                self.0
            }
        }
    };
}

positive_quantity!(Temperature, "K");
positive_quantity!(Pressure, "Pa absolute");
positive_quantity!(MolarMass, "kg/mol");

#[derive(Debug, Clone)]
pub struct Component {
    id: String,
    critical_temperature: Temperature,
    critical_pressure: Pressure,
    molar_mass: MolarMass,
    acentric_factor: f64,
}

impl Component {
    pub fn new(
        id: impl Into<String>,
        critical_temperature_k: f64,
        critical_pressure_pa: f64,
        molar_mass_kg_mol: f64,
        acentric_factor: f64,
    ) -> Result<Self, ThermoError> {
        let id = id.into();
        if id.trim().is_empty() || !acentric_factor.is_finite() {
            return Err(ThermoError::InvalidInput(
                "component requires an ID and a finite acentric factor".into(),
            ));
        }
        Ok(Self {
            id,
            critical_temperature: Temperature::new(critical_temperature_k)?,
            critical_pressure: Pressure::new(critical_pressure_pa)?,
            molar_mass: MolarMass::new(molar_mass_kg_mol)?,
            acentric_factor,
        })
    }

    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn critical_temperature(&self) -> Temperature {
        self.critical_temperature
    }
    pub fn critical_pressure(&self) -> Pressure {
        self.critical_pressure
    }
    pub fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }
    pub fn acentric_factor(&self) -> f64 {
        self.acentric_factor
    }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum CompositionBasis {
    Mole,
    Mass,
}

/// A normalized composition with distinct, ordered component IDs.
#[derive(Debug, Clone)]
pub struct Mixture {
    components: Vec<Component>,
    mole_fractions: Vec<f64>,
    binary_interactions: Vec<Vec<f64>>,
}

impl Mixture {
    /// Fractions sum to one within 1e-10; only roundoff is normalized.
    /// Requires a complete symmetric, constant kij matrix with a zero diagonal.
    pub fn new(
        components: Vec<Component>,
        fractions: Vec<f64>,
        basis: CompositionBasis,
        binary_interactions: Vec<Vec<f64>>,
    ) -> Result<Self, ThermoError> {
        let n = components.len();
        if n == 0 || fractions.len() != n || n > 64 {
            return Err(ThermoError::InvalidInput(
                "require 1–64 components and matching fractions".into(),
            ));
        }
        for (i, component) in components.iter().enumerate() {
            if components[..i].iter().any(|c| c.id() == component.id()) {
                return Err(ThermoError::InvalidInput("duplicate component IDs".into()));
            }
        }
        if fractions.iter().any(|v| !v.is_finite() || *v < 0.0)
            || (fractions.iter().sum::<f64>() - 1.0).abs() > 1e-10
        {
            return Err(ThermoError::InvalidInput(
                "fractions must be finite, nonnegative, and sum to one within 1e-10".into(),
            ));
        }
        if binary_interactions.len() != n || binary_interactions.iter().any(|row| row.len() != n) {
            return Err(ThermoError::MissingData(
                "a complete square binary_interactions matrix is required".into(),
            ));
        }
        for (i, row) in binary_interactions.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                if !value.is_finite()
                    || value >= 1.0
                    || (i == j && value != 0.0)
                    || value != binary_interactions[j][i]
                {
                    return Err(ThermoError::InvalidInput(
                        "kij must be finite, symmetric, less than one, with a zero diagonal".into(),
                    ));
                }
            }
        }
        let amounts: Vec<_> = fractions
            .iter()
            .zip(&components)
            .map(|(v, c)| match basis {
                CompositionBasis::Mole => *v,
                CompositionBasis::Mass => v / c.molar_mass().value(),
            })
            .collect();
        let sum: f64 = amounts.iter().sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(ThermoError::InvalidInput(
                "composition cannot be normalized".into(),
            ));
        }
        Ok(Self {
            components,
            mole_fractions: amounts.iter().map(|v| v / sum).collect(),
            binary_interactions,
        })
    }

    pub fn components(&self) -> &[Component] {
        &self.components
    }
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }
    pub fn binary_interactions(&self) -> &[Vec<f64>] {
        &self.binary_interactions
    }
    pub fn molar_mass_kg_mol(&self) -> f64 {
        self.components
            .iter()
            .zip(&self.mole_fractions)
            .map(|(c, x)| x * c.molar_mass().value())
            .sum()
    }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq)]
pub enum CubicModel {
    #[serde(rename = "PR")]
    PengRobinson,
    #[serde(rename = "PR-BM")]
    PengRobinsonBostonMathias,
}

/// Selects an algebraic homogeneous root, not an equilibrium phase.
#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum RootSelection {
    Smallest,
    Largest,
}

#[cfg(test)]
mod tests {
    use super::*;
    fn components() -> Vec<Component> {
        vec![
            Component::new("a", 200.0, 4e6, 0.016, 0.01).unwrap(),
            Component::new("b", 300.0, 5e6, 0.032, 0.1).unwrap(),
        ]
    }
    fn kij() -> Vec<Vec<f64>> {
        vec![vec![0.0; 2]; 2]
    }

    #[test]
    fn rejects_invalid_si_quantities() {
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Temperature::new(value).is_err());
            assert!(Pressure::new(value).is_err());
            assert!(MolarMass::new(value).is_err());
        }
    }

    #[test]
    fn converts_mass_to_mole_and_preserves_zero_components() {
        let m = Mixture::new(components(), vec![0.5, 0.5], CompositionBasis::Mass, kij()).unwrap();
        assert!((m.mole_fractions()[0] - 2.0 / 3.0).abs() < 1e-14);
        let pure =
            Mixture::new(components(), vec![1.0, 0.0], CompositionBasis::Mole, kij()).unwrap();
        assert_eq!(pure.mole_fractions(), &[1.0, 0.0]);
    }

    #[test]
    fn rejects_bad_compositions_and_parameter_matrices() {
        for fractions in [vec![], vec![0.2, 0.2], vec![-0.1, 1.1], vec![f64::NAN, 0.0]] {
            assert!(Mixture::new(components(), fractions, CompositionBasis::Mole, kij()).is_err());
        }
        for matrix in [
            vec![],
            vec![vec![0.0]],
            vec![vec![0.0, 0.1], vec![0.0, 0.0]],
            vec![vec![0.1, 0.0], vec![0.0, 0.0]],
            vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        ] {
            assert!(
                Mixture::new(components(), vec![0.5, 0.5], CompositionBasis::Mole, matrix).is_err()
            );
        }
        let c = components()[0].clone();
        assert!(
            Mixture::new(
                vec![c.clone(), c],
                vec![0.5, 0.5],
                CompositionBasis::Mole,
                kij()
            )
            .is_err()
        );
    }
}
