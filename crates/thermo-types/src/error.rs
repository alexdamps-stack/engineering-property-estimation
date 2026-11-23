use std::fmt;

/// Invalid inputs and arithmetic range failures are distinguishable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InputError {
    /// NaN or either infinity, including arithmetic overflow.
    NonFinite(&'static str),
    /// A quantity required to be strictly positive was zero or negative.
    NotPositive(&'static str),
    /// A nonnegative quantity was negative.
    Negative(&'static str),
    /// Empty or malformed identifier/name.
    InvalidIdentifier,
    /// An empty composition is not a mixture.
    EmptyComposition,
    /// Every component amount was zero.
    ZeroTotal,
    /// Component identifiers must be unique.
    DuplicateComponent,
    /// Fractions did not sum to one within the declared tolerance.
    NotNormalized,
    /// Parameter identities or ordering did not match the composition.
    ComponentMismatch,
    /// A positive contribution underflowed to zero during conversion.
    NumericalRange,
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(name) => write!(f, "{name} must be finite"),
            Self::NotPositive(name) => write!(f, "{name} must be strictly positive"),
            Self::Negative(name) => write!(f, "{name} must be nonnegative"),
            Self::InvalidIdentifier => f.write_str("identifier must be nonempty and well formed"),
            Self::EmptyComposition => f.write_str("composition must contain at least one component"),
            Self::ZeroTotal => f.write_str("composition must have a positive total amount"),
            Self::DuplicateComponent => f.write_str("composition contains duplicate component identifiers"),
            Self::NotNormalized => f.write_str("composition fractions must sum to one within the declared tolerance"),
            Self::ComponentMismatch => f.write_str("component identities or ordering do not match the composition"),
            Self::NumericalRange => f.write_str("a positive contribution underflowed during conversion"),
        }
    }
}

impl std::error::Error for InputError {}

/// Stable failure categories for model evaluation, equilibrium solves and host interfaces.
/// Invalid solver configuration uses `UnsupportedCapability`; physical states outside
/// a model's supported envelope use `OutsideValidityDomain`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CalculationError {
    InvalidInput(InputError),
    MissingParameter(String),
    UnsupportedCapability(String),
    OutsideValidityDomain(String),
    NonConvergence(String),
    NumericalFailure(String),
}

impl fmt::Display for CalculationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(error) => write!(f, "invalid calculation input: {error}"),
            Self::MissingParameter(message) => write!(f, "missing parameter: {message}"),
            Self::UnsupportedCapability(message) => write!(f, "unsupported calculation: {message}"),
            Self::OutsideValidityDomain(message) => write!(f, "outside validity domain: {message}"),
            Self::NonConvergence(message) => write!(f, "calculation did not converge: {message}"),
            Self::NumericalFailure(message) => write!(f, "numerical failure: {message}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_preserves_plain_messages_without_debug_escaping() {
        let error = CalculationError::NonConvergence("trial \"liquid\" exhausted\nretry".into());
        assert_eq!(error.to_string(), "calculation did not converge: trial \"liquid\" exhausted\nretry");
        assert_eq!(CalculationError::from(InputError::NonFinite("temperature")).to_string(), "invalid calculation input: temperature must be finite");
        for error in [InputError::NotPositive("pressure"), InputError::Negative("amount"), InputError::InvalidIdentifier, InputError::EmptyComposition, InputError::ZeroTotal, InputError::DuplicateComponent, InputError::NotNormalized, InputError::ComponentMismatch, InputError::NumericalRange] {
            let message = error.to_string();
            assert!(!message.is_empty());
            assert!(!message.contains("InputError"));
        }
    }
}

impl std::error::Error for CalculationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInput(error) => Some(error),
            _ => None,
        }
    }
}

impl From<InputError> for CalculationError {
    fn from(value: InputError) -> Self {
        Self::InvalidInput(value)
    }
}
