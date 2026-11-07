use crate::InputError;

macro_rules! quantity {
    ($name:ident, $doc:literal, $predicate:expr, $error:ident) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
        pub struct $name(f64);
        impl $name {
            /// Construct from the documented SI unit, checking finiteness and range.
            pub fn new(value: f64) -> Result<Self, InputError> {
                if !value.is_finite() {
                    return Err(InputError::NonFinite(stringify!($name)));
                }
                if !($predicate)(value) {
                    return Err(InputError::$error(stringify!($name)));
                }
                Ok(Self(value))
            }
            /// Return the value in the documented SI unit.
            pub fn value(self) -> f64 {
                self.0
            }
        }
    };
}

quantity!(
    FiniteSi,
    "Finite scalar in a caller-documented SI unit; signed values allowed.",
    |_: f64| true,
    NonFinite
);
quantity!(
    Temperature,
    "Strictly positive absolute temperature in kelvin (K).",
    |v| v > 0.0,
    NotPositive
);
quantity!(
    Pressure,
    "Strictly positive absolute pressure in pascals (Pa).",
    |v| v > 0.0,
    NotPositive
);
quantity!(
    MolarMass,
    "Strictly positive molar mass in kg/mol.",
    |v| v > 0.0,
    NotPositive
);
quantity!(
    MolarFlow,
    "Nonnegative molar flow in mol/s; zero is allowed.",
    |v| v >= 0.0,
    Negative
);
quantity!(
    MassFlow,
    "Nonnegative mass flow in kg/s; zero is allowed.",
    |v| v >= 0.0,
    Negative
);
quantity!(
    PositiveDimensionless,
    "Strictly positive finite dimensionless scalar (for example Z).",
    |v| v > 0.0,
    NotPositive
);

impl Temperature {
    /// Convert degrees Celsius to absolute kelvin; absolute zero is rejected.
    pub fn from_celsius(value: f64) -> Result<Self, InputError> {
        if !value.is_finite() {
            return Err(InputError::NonFinite("Celsius temperature"));
        }
        Self::new(value + 273.15)
    }
}

impl Pressure {
    /// Convert absolute bar to Pa. Gauge pressure requires an explicit upstream conversion.
    pub fn from_bar_absolute(value: f64) -> Result<Self, InputError> {
        Self::new(value * 100_000.0)
    }
}

pub(crate) fn checked_product(a: f64, b: f64) -> Result<f64, InputError> {
    let value = a * b;
    if !value.is_finite() {
        Err(InputError::NonFinite("conversion result"))
    } else if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(InputError::NumericalRange)
    } else {
        Ok(value)
    }
}

impl MolarFlow {
    /// Convert using a supplied mixture molar mass in kg/mol.
    pub fn to_mass(self, molar_mass: MolarMass) -> Result<MassFlow, InputError> {
        MassFlow::new(checked_product(self.value(), molar_mass.value())?)
    }
}

impl MassFlow {
    /// Convert using a supplied mixture molar mass in kg/mol.
    pub fn to_molar(self, molar_mass: MolarMass) -> Result<MolarFlow, InputError> {
        let value = self.value() / molar_mass.value();
        if self.value() > 0.0 && value == 0.0 {
            return Err(InputError::NumericalRange);
        }
        MolarFlow::new(value)
    }
}
