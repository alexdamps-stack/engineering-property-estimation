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

/// Failure vocabulary shared by models, solvers, packages and their hosts.
///
/// Kinds are chosen by what is wrong, so equivalent mistakes map to the same kind
/// across entry points:
/// - invalid solver or package options (search bounds, tolerances, thresholds,
///   oversized iteration budgets) are
///   [`UnsupportedCapability`](Self::UnsupportedCapability);
/// - a zero iteration limit is [`NonConvergence`](Self::NonConvergence);
/// - a requested state (temperature, pressure, density) outside the declared
///   range of a model, solver or package is
///   [`OutsideValidityDomain`](Self::OutsideValidityDomain).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CalculationError {
    /// A quantity or composition failed its input invariants.
    InvalidInput(InputError),
    /// A required parameter (component, interaction, coefficient) is unavailable.
    MissingParameter(String),
    /// The request, or an option it carries, is not supported by this implementation.
    UnsupportedCapability(String),
    /// The requested state is outside the declared validity range.
    OutsideValidityDomain(String),
    /// An iterative method did not converge within its budget.
    NonConvergence(String),
    /// Arithmetic or numerical breakdown, such as a singular system.
    NumericalFailure(String),
}

impl fmt::Display for CalculationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(error) => write!(f, "{error}"),
            Self::MissingParameter(m) => write!(f, "missing parameter: {m}"),
            Self::UnsupportedCapability(m) => write!(f, "unsupported: {m}"),
            Self::OutsideValidityDomain(m) => write!(f, "outside validity domain: {m}"),
            Self::NonConvergence(m) => write!(f, "did not converge: {m}"),
            Self::NumericalFailure(m) => write!(f, "numerical failure: {m}"),
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
