use crate::{
    InputError, MolarFlow, PositiveDimensionless, Pressure, Temperature, quantity::checked_product,
};

/// SI molar gas constant, J/(mol K), from exact SI defining constants.
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Named reference conditions; there is deliberately no default "standard".
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceConditions {
    name: String,
    temperature: Temperature,
    pressure: Pressure,
}

impl ReferenceConditions {
    pub fn new(
        name: impl Into<String>,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<Self, InputError> {
        let name = name.into();
        if name.trim().is_empty() || name.chars().any(char::is_control) {
            return Err(InputError::InvalidIdentifier);
        }
        Ok(Self {
            name,
            temperature,
            pressure,
        })
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub fn pressure(&self) -> Pressure {
        self.pressure
    }
}

/// Explicit gas-volume assumption at the named reference conditions.
/// The caller must establish applicability and supply Z for the same mixture/state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GasVolumeAssumption {
    IdealGas,
    CompressibilityFactor(PositiveDimensionless),
}

impl GasVolumeAssumption {
    fn molar_volume(self, reference: &ReferenceConditions) -> Result<f64, InputError> {
        let z = match self {
            Self::IdealGas => 1.0,
            Self::CompressibilityFactor(z) => z.value(),
        };
        let rt = checked_product(MOLAR_GAS_CONSTANT, reference.temperature.value())?;
        let volume = checked_product(z, rt)? / reference.pressure.value();
        if !volume.is_finite() {
            return Err(InputError::NonFinite("reference molar volume"));
        }
        if volume == 0.0 {
            return Err(InputError::NumericalRange);
        }
        Ok(volume)
    }
}

/// Nonnegative gas reference-volume flow in m^3/s, carrying conditions and assumption.
/// This is neither actual process volume nor a liquid-volume conversion.
#[derive(Debug, Clone, PartialEq)]
pub struct StandardVolumeFlow {
    value: f64,
    reference: ReferenceConditions,
    assumption: GasVolumeAssumption,
}

impl StandardVolumeFlow {
    pub fn new(
        value: f64,
        reference: ReferenceConditions,
        assumption: GasVolumeAssumption,
    ) -> Result<Self, InputError> {
        if !value.is_finite() {
            return Err(InputError::NonFinite("standard volume flow"));
        }
        if value < 0.0 {
            return Err(InputError::Negative("standard volume flow"));
        }
        assumption.molar_volume(&reference)?;
        Ok(Self {
            value,
            reference,
            assumption,
        })
    }

    /// Convert m^3/h explicitly to the stored m^3/s.
    pub fn from_cubic_metres_per_hour(
        value: f64,
        reference: ReferenceConditions,
        assumption: GasVolumeAssumption,
    ) -> Result<Self, InputError> {
        if !value.is_finite() {
            return Err(InputError::NonFinite("standard volume flow"));
        }
        if value < 0.0 {
            return Err(InputError::Negative("standard volume flow"));
        }
        if value > 0.0 && value / 3600.0 == 0.0 {
            return Err(InputError::NumericalRange);
        }
        Self::new(value / 3600.0, reference, assumption)
    }
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn reference(&self) -> &ReferenceConditions {
        &self.reference
    }
    pub fn assumption(&self) -> GasVolumeAssumption {
        self.assumption
    }

    /// n_dot = p_ref V_dot / (Z_ref R T_ref).
    pub fn to_molar(&self) -> Result<MolarFlow, InputError> {
        let value = self.value / self.assumption.molar_volume(&self.reference)?;
        if self.value > 0.0 && value == 0.0 {
            return Err(InputError::NumericalRange);
        }
        MolarFlow::new(value)
    }
}

impl MolarFlow {
    /// V_dot = n_dot Z_ref R T_ref / p_ref; all reference choices are explicit.
    pub fn to_standard_volume(
        self,
        reference: ReferenceConditions,
        assumption: GasVolumeAssumption,
    ) -> Result<StandardVolumeFlow, InputError> {
        let value = checked_product(self.value(), assumption.molar_volume(&reference)?)?;
        StandardVolumeFlow::new(value, reference, assumption)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_molar_volume_matches_prescribed_ideal_gas_state() {
        let reference = ReferenceConditions::new(
            "test",
            Temperature::new(300.0).unwrap(),
            Pressure::new(100_000.0).unwrap(),
        )
        .unwrap();
        let volume = GasVolumeAssumption::IdealGas
            .molar_volume(&reference)
            .unwrap();
        assert!((volume - 0.024_943_387_854_459_72).abs() < 1e-16);
    }

    #[test]
    fn unrepresentable_reference_volume_is_an_error() {
        let reference = ReferenceConditions::new(
            "range test",
            Temperature::new(f64::MAX).unwrap(),
            Pressure::new(f64::MAX).unwrap(),
        )
        .unwrap();
        assert!(matches!(
            GasVolumeAssumption::IdealGas.molar_volume(&reference),
            Err(InputError::NonFinite(_))
        ));
    }
}
