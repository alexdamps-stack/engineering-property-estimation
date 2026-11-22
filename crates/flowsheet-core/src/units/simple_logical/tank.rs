//! Steady-state tank with one bulk material outlet (#178).

use crate::engine::{MaterialStream, empty};
use crate::unit::{SeparatorThermal, Sign, SpecIssue, check_optional, check_value};
use crate::units::{CalcContext, UnitOutput};

/// No-accumulation tank with an assigned pressure drop and thermal mode.
#[derive(Debug, Clone, PartialEq)]
pub struct TankSpec {
    pub thermal: SeparatorThermal,
    /// Pressure loss in Pa, subtracted from the inlet absolute pressure.
    pub pressure_drop_pa: f64,
}

impl Default for TankSpec {
    fn default() -> Self {
        Self {
            thermal: SeparatorThermal::Adiabatic,
            pressure_drop_pa: 0.0,
        }
    }
}

impl TankSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        check_value(
            &mut issues,
            "pressure drop",
            self.pressure_drop_pa,
            Sign::NonNegative,
        );
        match self.thermal {
            SeparatorThermal::Adiabatic => {}
            SeparatorThermal::OutletTemperature(t) => {
                check_optional(&mut issues, "outlet temperature", t, Sign::Positive);
            }
            SeparatorThermal::Duty(q) => check_optional(&mut issues, "duty", q, Sign::Any),
        }
        issues
    }
}

pub(super) fn calculate(
    spec: &TankSpec,
    inlet: &MaterialStream,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    let n = inlet.molar_flow_mol_s;
    if n == 0.0 {
        return Ok(UnitOutput::new(
            vec![empty()],
            0.0,
            vec!["Zero inlet flow: the tank outlet is empty and no duty is applied.".into()],
        ));
    }
    let state = inlet
        .state
        .as_ref()
        .ok_or("tank inlet has no thermodynamic state")?;
    let p = state.pressure_pa - spec.pressure_drop_pa;
    if !p.is_finite() || p <= 0.0 {
        return Err("tank pressure drop leaves nonpositive absolute outlet pressure".into());
    }
    let z = &state.mole_fractions;
    let h = state.enthalpy();
    let (outlet, duty) = match spec.thermal {
        // Preserve the inlet's phase fractions when the conditions are unchanged.
        SeparatorThermal::Adiabatic if spec.pressure_drop_pa == 0.0 => (state.clone(), 0.0),
        SeparatorThermal::Adiabatic => (ctx.ph(z, p, h)?, 0.0),
        SeparatorThermal::Duty(Some(q)) => (ctx.ph(z, p, h + q / n)?, q),
        SeparatorThermal::OutletTemperature(Some(t)) => {
            let out = ctx.tp(z, t, p)?;
            let q = n * (out.enthalpy() - h);
            (out, q)
        }
        _ => return Err("missing tank thermal specification".into()),
    };
    Ok(UnitOutput::new(
        vec![MaterialStream {
            molar_flow_mol_s: n,
            state: Some(outlet),
        }],
        duty,
        vec![],
    ))
}
