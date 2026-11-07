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
        write!(f, "invalid calculation input: {self:?}")
    }
}

impl std::error::Error for InputError {}

/// Failure vocabulary for future calculations; no model or solver is implemented.
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
        write!(f, "calculation failed: {self:?}")
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
