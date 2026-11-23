//! Display unit sets (proposal section 4.3). The flowsheet model stores SI
//! values with absolute pressure; conversions here affect presentation and
//! editing only.

/// Selectable display unit set: the unit set stored with the project
/// (SI, or process metric — the default for process engineers). The
/// presentation behaviour lives in [`DisplayUnitsExt`].
pub use flowsheet_project::DisplayUnitSet as DisplayUnits;

/// Physical quantity of an editable field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantity {
    /// Temperature, stored in K.
    Temperature,
    /// Absolute pressure, stored in Pa.
    Pressure,
    /// Pressure difference, stored in Pa.
    PressureDrop,
    /// Molar flow, stored in mol/s.
    MolarFlow,
    /// Heat duty, stored in W, positive into the unit.
    Duty,
    /// Dimensionless fraction from 0 to 1 (split fraction, efficiency).
    Fraction,
    /// Mass flow, stored in kg/s.
    MassFlow,
    /// Molar enthalpy, stored in J/mol.
    MolarEnthalpy,
    /// Molar entropy, stored in J/(mol K).
    MolarEntropy,
    /// Mass density, stored in kg/m³.
    Density,
    /// Volume, stored in m³.
    Volume,
    /// Length, stored in m.
    Length,
    /// Time, stored in s.
    Time,
    /// Temperature difference, stored in K (may be negative).
    TemperatureDifference,
    /// Heat-transfer coefficient, stored in W/(m² K).
    HeatTransferCoefficient,
}

/// Labels, unit symbols and conversions of a [`DisplayUnits`] set.
pub trait DisplayUnitsExt: Copy {
    /// All unit sets, in menu order.
    const ALL: [DisplayUnits; 2];

    /// Menu label.
    fn label(self) -> &'static str;

    /// Unit symbol for a quantity. Absolute pressures are always marked `(a)`.
    fn symbol(self, q: Quantity) -> &'static str;

    /// Converts a stored SI value to this display set.
    fn si_to_display(self, q: Quantity, si: f64) -> f64;

    /// Converts a displayed value back to SI.
    fn display_to_si(self, q: Quantity, shown: f64) -> f64;

    /// Drag speed for editing a displayed value.
    fn drag_speed(self, q: Quantity) -> f64;

    /// Lowest displayed value the editor allows. Duty may be negative (heat
    /// removed); all other quantities are non-negative in SI. Strict
    /// positivity is reported by the flowsheet diagnostics, not clamped away.
    fn display_min(self, q: Quantity) -> f64;

    /// Highest displayed value the editor allows (fractions are capped at 1).
    fn display_max(self, q: Quantity) -> f64;
}

impl DisplayUnitsExt for DisplayUnits {
    const ALL: [DisplayUnits; 2] = [DisplayUnits::ProcessMetric, DisplayUnits::Si];

    fn label(self) -> &'static str {
        match self {
            DisplayUnits::Si => "SI",
            DisplayUnits::ProcessMetric => "Process metric",
        }
    }

    fn symbol(self, q: Quantity) -> &'static str {
        match (self, q) {
            (DisplayUnits::Si, Quantity::Temperature) => "K",
            (DisplayUnits::Si, Quantity::Pressure) => "Pa(a)",
            (DisplayUnits::Si, Quantity::PressureDrop) => "Pa",
            (DisplayUnits::Si, Quantity::MolarFlow) => "mol/s",
            (DisplayUnits::Si, Quantity::Duty) => "W",
            (DisplayUnits::Si, Quantity::Fraction) => "–",
            (DisplayUnits::Si, Quantity::MassFlow) => "kg/s",
            (DisplayUnits::Si, Quantity::MolarEnthalpy) => "J/mol",
            (DisplayUnits::Si, Quantity::MolarEntropy) => "J/(mol K)",
            (_, Quantity::Density) => "kg/m³",
            (_, Quantity::Volume) => "m³",
            (_, Quantity::Length) => "m",
            (_, Quantity::Time) => "s",
            (_, Quantity::TemperatureDifference) => "K",
            (_, Quantity::HeatTransferCoefficient) => "W/(m² K)",
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => "°C",
            (DisplayUnits::ProcessMetric, Quantity::Pressure) => "bar(a)",
            (DisplayUnits::ProcessMetric, Quantity::PressureDrop) => "bar",
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => "kmol/h",
            (DisplayUnits::ProcessMetric, Quantity::Duty) => "kW",
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => "%",
            (DisplayUnits::ProcessMetric, Quantity::MassFlow) => "kg/h",
            (DisplayUnits::ProcessMetric, Quantity::MolarEnthalpy) => "kJ/kmol",
            (DisplayUnits::ProcessMetric, Quantity::MolarEntropy) => "kJ/(kmol K)",
        }
    }

    fn si_to_display(self, q: Quantity, si: f64) -> f64 {
        match (self, q) {
            (DisplayUnits::Si, _) => si,
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => si - 273.15,
            (DisplayUnits::ProcessMetric, Quantity::Pressure | Quantity::PressureDrop) => {
                si / 1.0e5
            }
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => si * 3.6,
            (DisplayUnits::ProcessMetric, Quantity::Duty) => si / 1.0e3,
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => si * 100.0,
            (DisplayUnits::ProcessMetric, Quantity::MassFlow) => si * 3600.0,
            // kJ/kmol and kJ/(kmol K) equal J/mol and J/(mol K) numerically.
            (
                DisplayUnits::ProcessMetric,
                Quantity::MolarEnthalpy | Quantity::MolarEntropy | Quantity::Density,
            ) => si,
            (
                DisplayUnits::ProcessMetric,
                Quantity::Volume
                | Quantity::Length
                | Quantity::Time
                | Quantity::TemperatureDifference
                | Quantity::HeatTransferCoefficient,
            ) => si,
        }
    }

    fn display_to_si(self, q: Quantity, shown: f64) -> f64 {
        match (self, q) {
            (DisplayUnits::Si, _) => shown,
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => shown + 273.15,
            (DisplayUnits::ProcessMetric, Quantity::Pressure | Quantity::PressureDrop) => {
                shown * 1.0e5
            }
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => shown / 3.6,
            (DisplayUnits::ProcessMetric, Quantity::Duty) => shown * 1.0e3,
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => shown / 100.0,
            (DisplayUnits::ProcessMetric, Quantity::MassFlow) => shown / 3600.0,
            (
                DisplayUnits::ProcessMetric,
                Quantity::MolarEnthalpy | Quantity::MolarEntropy | Quantity::Density,
            ) => shown,
            (
                DisplayUnits::ProcessMetric,
                Quantity::Volume
                | Quantity::Length
                | Quantity::Time
                | Quantity::TemperatureDifference
                | Quantity::HeatTransferCoefficient,
            ) => shown,
        }
    }

    fn drag_speed(self, q: Quantity) -> f64 {
        match (self, q) {
            (DisplayUnits::Si, Quantity::Fraction) => 0.001,
            (DisplayUnits::Si, Quantity::Pressure | Quantity::PressureDrop | Quantity::Duty) => {
                100.0
            }
            (_, Quantity::MolarFlow) | (DisplayUnits::ProcessMetric, Quantity::PressureDrop) => {
                0.01
            }
            _ => 0.1,
        }
    }

    fn display_min(self, q: Quantity) -> f64 {
        match q {
            Quantity::Duty | Quantity::TemperatureDifference => f64::NEG_INFINITY,
            _ => self.si_to_display(q, 0.0),
        }
    }

    fn display_max(self, q: Quantity) -> f64 {
        match q {
            Quantity::Fraction => self.si_to_display(q, 1.0),
            _ => f64::INFINITY,
        }
    }
}

/// Initial SI value when the user chooses to specify a previously unset
/// field: 298.15 K, 101 325 Pa, 1 mol/s, and zero for differences and duty.
/// These are editing starting points, not defaults applied to calculations.
pub fn starting_value(q: Quantity) -> f64 {
    match q {
        Quantity::Temperature => 298.15,
        Quantity::Pressure => 101_325.0,
        Quantity::MolarFlow => 1.0,
        Quantity::PressureDrop | Quantity::Duty => 0.0,
        Quantity::Fraction => 0.5,
        Quantity::MassFlow => 1.0,
        Quantity::MolarEnthalpy | Quantity::MolarEntropy | Quantity::Density => 0.0,
        Quantity::Volume | Quantity::Length => 1.0,
        Quantity::Time => 60.0,
        Quantity::TemperatureDifference => 0.0,
        Quantity::HeatTransferCoefficient => 100.0,
    }
}

/// Formats a calculated value for display with about six significant
/// digits. Nonfinite values are never shown as numbers.
pub fn format_value(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    if v == 0.0 {
        return "0".into();
    }
    let magnitude = v.abs().log10().floor() as i32;
    if !(-3..6).contains(&magnitude) {
        return format!("{v:.4e}");
    }
    let decimals = (5 - magnitude).clamp(0, 6) as usize;
    format!("{v:.decimals$}")
}

/// Scaled residuals at or below this size are round-off and shown as a bound.
pub const RESIDUAL_DISPLAY_FLOOR: f64 = 1e-12;

/// Formats a dimensionless scaled balance residual with two significant
/// digits. Values at or below [`RESIDUAL_DISPLAY_FLOOR`] are shown as
/// `< 1e-12` regardless of sign: their digits are floating-point round-off
/// that changes with solver numerics and carries no information.
pub fn format_residual(v: f64) -> String {
    if !v.is_finite() {
        return "—".into();
    }
    if v.abs() <= RESIDUAL_DISPLAY_FLOOR {
        return format!("< {RESIDUAL_DISPLAY_FLOOR:.0e}");
    }
    format!("{v:.1e}")
}

/// Formats a value of quantity `q` in the display set, with its unit.
pub fn format_quantity(u: DisplayUnits, q: Quantity, si: f64) -> String {
    format!("{} {}", format_value(u.si_to_display(q, si)), u.symbol(q))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_Q: [Quantity; 10] = [
        Quantity::MassFlow,
        Quantity::MolarEnthalpy,
        Quantity::MolarEntropy,
        Quantity::Density,
        Quantity::Fraction,
        Quantity::Temperature,
        Quantity::Pressure,
        Quantity::PressureDrop,
        Quantity::MolarFlow,
        Quantity::Duty,
    ];

    #[test]
    fn conversions_round_trip() {
        for units in DisplayUnits::ALL {
            for q in ALL_Q {
                let si = 1234.5;
                let back = units.display_to_si(q, units.si_to_display(q, si));
                assert!((back - si).abs() <= 1e-9 * si, "{units:?} {q:?}");
            }
        }
    }

    #[test]
    fn process_metric_matches_definitions() {
        let m = DisplayUnits::ProcessMetric;
        assert_eq!(DisplayUnits::default(), m);
        assert!(m.si_to_display(Quantity::Temperature, 273.15).abs() < 1e-12);
        assert!((m.si_to_display(Quantity::Pressure, 101_325.0) - 1.01325).abs() < 1e-12);
        assert!((m.si_to_display(Quantity::MolarFlow, 1.0) - 3.6).abs() < 1e-12);
        assert!((m.si_to_display(Quantity::Duty, -2500.0) + 2.5).abs() < 1e-12);
        // A pressure difference has no offset and no absolute marker.
        assert_eq!(m.si_to_display(Quantity::PressureDrop, 0.0), 0.0);
        assert_eq!(m.symbol(Quantity::PressureDrop), "bar");
        assert_eq!(m.symbol(Quantity::Pressure), "bar(a)");
        assert!((m.display_min(Quantity::Temperature) + 273.15).abs() < 1e-12);
        assert_eq!(m.display_min(Quantity::Duty), f64::NEG_INFINITY);
        assert!((m.si_to_display(Quantity::MassFlow, 1.0) - 3600.0).abs() < 1e-9);
    }

    #[test]
    fn values_are_formatted_without_nan() {
        assert_eq!(format_value(f64::NAN), "—");
        assert_eq!(format_value(f64::INFINITY), "—");
        assert_eq!(format_value(0.0), "0");
        assert_eq!(format_value(-58.9655), "-58.9655");
        assert_eq!(format_value(15.0), "15.0000");
        assert_eq!(format_value(0.964993), "0.964993");
        assert_eq!(format_value(-555838.418), "-555838");
        assert_eq!(format_value(2.5e6), "2.5000e6");
        assert_eq!(format_value(1.2e-7), "1.2000e-7");
        assert_eq!(
            format_quantity(DisplayUnits::ProcessMetric, Quantity::Pressure, 1.5e6),
            "15.0000 bar(a)"
        );
    }

    #[test]
    fn round_off_residuals_are_shown_as_a_bound() {
        for v in [0.0, 1.8e-16, 9.9e-17, -2.9e-17, -1e-12, 1e-12] {
            assert_eq!(format_residual(v), "< 1e-12", "{v:e}");
        }
        assert_eq!(format_residual(3.4e-9), "3.4e-9");
        assert_eq!(format_residual(-2.5e-6), "-2.5e-6");
        assert_eq!(format_residual(f64::NAN), "—");
    }
}
