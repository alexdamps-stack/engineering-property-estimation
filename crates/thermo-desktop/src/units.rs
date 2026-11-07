//! Display unit sets (proposal section 4.3). The flowsheet model stores SI
//! values with absolute pressure; conversions here affect presentation and
//! editing only.

/// Selectable display unit set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DisplayUnits {
    /// K, Pa, mol/s, W.
    Si,
    /// Â°C, bar, kmol/h, kW â the default for process engineers.
    #[default]
    ProcessMetric,
}

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
}

impl DisplayUnits {
    /// All unit sets, in menu order.
    pub const ALL: [DisplayUnits; 2] = [DisplayUnits::ProcessMetric, DisplayUnits::Si];

    /// Menu label.
    pub fn label(self) -> &'static str {
        match self {
            DisplayUnits::Si => "SI",
            DisplayUnits::ProcessMetric => "Process metric",
        }
    }

    /// Unit symbol for a quantity. Absolute pressures are always marked `(a)`.
    pub fn symbol(self, q: Quantity) -> &'static str {
        match (self, q) {
            (DisplayUnits::Si, Quantity::Temperature) => "K",
            (DisplayUnits::Si, Quantity::Pressure) => "Pa(a)",
            (DisplayUnits::Si, Quantity::PressureDrop) => "Pa",
            (DisplayUnits::Si, Quantity::MolarFlow) => "mol/s",
            (DisplayUnits::Si, Quantity::Duty) => "W",
            (DisplayUnits::Si, Quantity::Fraction) => "–",
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => "Â°C",
            (DisplayUnits::ProcessMetric, Quantity::Pressure) => "bar(a)",
            (DisplayUnits::ProcessMetric, Quantity::PressureDrop) => "bar",
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => "kmol/h",
            (DisplayUnits::ProcessMetric, Quantity::Duty) => "kW",
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => "%",
        }
    }

    /// Converts a stored SI value to this display set.
    pub fn si_to_display(self, q: Quantity, si: f64) -> f64 {
        match (self, q) {
            (DisplayUnits::Si, _) => si,
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => si - 273.15,
            (DisplayUnits::ProcessMetric, Quantity::Pressure | Quantity::PressureDrop) => {
                si / 1.0e5
            }
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => si * 3.6,
            (DisplayUnits::ProcessMetric, Quantity::Duty) => si / 1.0e3,
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => si * 100.0,
        }
    }

    /// Converts a displayed value back to SI.
    pub fn display_to_si(self, q: Quantity, shown: f64) -> f64 {
        match (self, q) {
            (DisplayUnits::Si, _) => shown,
            (DisplayUnits::ProcessMetric, Quantity::Temperature) => shown + 273.15,
            (DisplayUnits::ProcessMetric, Quantity::Pressure | Quantity::PressureDrop) => {
                shown * 1.0e5
            }
            (DisplayUnits::ProcessMetric, Quantity::MolarFlow) => shown / 3.6,
            (DisplayUnits::ProcessMetric, Quantity::Duty) => shown * 1.0e3,
            (DisplayUnits::ProcessMetric, Quantity::Fraction) => shown / 100.0,
        }
    }

    /// Drag speed for editing a displayed value.
    pub fn drag_speed(self, q: Quantity) -> f64 {
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

    /// Lowest displayed value the editor allows. Duty may be negative (heat
    /// removed); all other quantities are non-negative in SI. Strict
    /// positivity is reported by the flowsheet diagnostics, not clamped away.
    pub fn display_min(self, q: Quantity) -> f64 {
        match q {
            Quantity::Duty => f64::NEG_INFINITY,
            _ => self.si_to_display(q, 0.0),
        }
    }

    /// Highest displayed value the editor allows (fractions are capped at 1).
    pub fn display_max(self, q: Quantity) -> f64 {
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_Q: [Quantity; 6] = [
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
    }
}
