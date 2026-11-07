use crate::{InputError, MolarMass, quantity::checked_product};

/// Stable, case-sensitive identity, independent of display names and array indices.
/// ASCII letters, digits, `_`, `-`, `.`, and `:` are allowed; aliases are not resolved.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ComponentId(String);

impl ComponentId {
    pub fn new(value: impl Into<String>) -> Result<Self, InputError> {
        let value = value.into();
        if value.is_empty()
            || !value
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"_-.:".contains(&c))
        {
            return Err(InputError::InvalidIdentifier);
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Composition is always explicitly molar or mass based.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompositionBasis {
    Mole,
    Mass,
}

/// Input policy is mandatory; accepted fractions are always rescaled to unit sum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizationPolicy {
    /// Reject a sum differing from one by more than `COMPOSITION_SUM_TOLERANCE`.
    RequireUnitSum,
    /// Treat nonnegative entries as relative amounts and explicitly normalize.
    Normalize,
}

/// Absolute dimensionless acceptance tolerance, not a property accuracy claim.
pub const COMPOSITION_SUM_TOLERANCE: f64 = 1e-10;

/// Ordered component identities and fractions; zero entries remain present.
#[derive(Debug, Clone, PartialEq)]
pub struct Composition {
    basis: CompositionBasis,
    entries: Vec<(ComponentId, f64)>,
}

impl Composition {
    pub fn new(
        basis: CompositionBasis,
        mut entries: Vec<(ComponentId, f64)>,
        policy: NormalizationPolicy,
    ) -> Result<Self, InputError> {
        if entries.is_empty() {
            return Err(InputError::EmptyComposition);
        }
        let mut maximum = 0.0_f64;
        for (i, (id, value)) in entries.iter().enumerate() {
            if !value.is_finite() {
                return Err(InputError::NonFinite("composition"));
            }
            if *value < 0.0 {
                return Err(InputError::Negative("composition"));
            }
            if entries[..i].iter().any(|(other, _)| other == id) {
                return Err(InputError::DuplicateComponent);
            }
            maximum = maximum.max(*value);
        }
        if maximum == 0.0 {
            return Err(InputError::ZeroTotal);
        }
        if policy == NormalizationPolicy::RequireUnitSum {
            let total: f64 = entries.iter().map(|(_, value)| value).sum();
            if !total.is_finite() || (total - 1.0).abs() > COMPOSITION_SUM_TOLERANCE {
                return Err(InputError::NotNormalized);
            }
        }
        // Scaling avoids overflow when normalizing large relative amounts.
        let total: f64 = entries.iter().map(|(_, value)| value / maximum).sum();
        if !total.is_finite() {
            return Err(InputError::NonFinite("composition sum"));
        }
        for (_, value) in &mut entries {
            let normalized = (*value / maximum) / total;
            if *value > 0.0 && normalized == 0.0 {
                return Err(InputError::NumericalRange);
            }
            *value = normalized;
        }
        Ok(Self { basis, entries })
    }

    pub fn basis(&self) -> CompositionBasis {
        self.basis
    }
    pub fn entries(&self) -> &[(ComponentId, f64)] {
        &self.entries
    }

    fn validate_masses(&self, masses: &[(ComponentId, MolarMass)]) -> Result<(), InputError> {
        if masses.len() != self.entries.len()
            || self
                .entries
                .iter()
                .zip(masses)
                .any(|((id, _), (other, _))| id != other)
        {
            return Err(InputError::ComponentMismatch);
        }
        Ok(())
    }

    /// Convert basis with molar masses explicitly keyed in the same order.
    /// Masses are required even for zero components and same-basis conversions.
    pub fn convert(
        &self,
        basis: CompositionBasis,
        masses: &[(ComponentId, MolarMass)],
    ) -> Result<Self, InputError> {
        self.validate_masses(masses)?;
        if basis == self.basis {
            return Ok(self.clone());
        }
        let scale = match self.basis {
            CompositionBasis::Mole => masses
                .iter()
                .map(|(_, m)| m.value())
                .fold(0.0_f64, f64::max),
            CompositionBasis::Mass => masses
                .iter()
                .map(|(_, m)| m.value())
                .fold(f64::INFINITY, f64::min),
        };
        let entries = self
            .entries
            .iter()
            .zip(masses)
            .map(|((id, fraction), (_, mass))| {
                let factor = match self.basis {
                    CompositionBasis::Mole => mass.value() / scale,
                    CompositionBasis::Mass => scale / mass.value(),
                };
                if *fraction > 0.0 && factor == 0.0 {
                    return Err(InputError::NumericalRange);
                }
                Ok((id.clone(), checked_product(*fraction, factor)?))
            })
            .collect::<Result<Vec<_>, InputError>>()?;
        Self::new(basis, entries, NormalizationPolicy::Normalize)
    }

    /// Mixture molar mass in kg/mol: sum(x_i M_i), or its mass-basis equivalent.
    pub fn molar_mass(&self, masses: &[(ComponentId, MolarMass)]) -> Result<MolarMass, InputError> {
        self.validate_masses(masses)?;
        let mole = self.convert(CompositionBasis::Mole, masses)?;
        let mut total = 0.0;
        for ((_, fraction), (_, mass)) in mole.entries.iter().zip(masses) {
            total += checked_product(*fraction, mass.value())?;
        }
        MolarMass::new(total)
    }
}
