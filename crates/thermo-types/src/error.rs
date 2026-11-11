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
        f.write_str("invalid calculation input: ")?;
        match self {
            Self::NonFinite(what) => write!(f, "{what} is not finite"),
            Self::NotPositive(what) => write!(f, "{what} must be strictly positive"),
            Self::Negative(what) => write!(f, "{what} must not be negative"),
            Self::InvalidIdentifier => f.write_str("empty or malformed identifier"),
            Self::EmptyComposition => f.write_str("the composition is empty"),
            Self::ZeroTotal => f.write_str("every component amount is zero"),
            Self::DuplicateComponent => f.write_str("component identifiers are not unique"),
            Self::NotNormalized => f.write_str("fractions do not sum to one within tolerance"),
            Self::ComponentMismatch => {
                f.write_str("parameter identities or order do not match the composition")
            }
            Self::NumericalRange => {
                f.write_str("a positive contribution underflowed to zero during conversion")
            }
        }
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
