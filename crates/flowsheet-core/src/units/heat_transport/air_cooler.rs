//! Specified-duty or specified-outlet-temperature air cooler (#179).

use crate::engine::{MaterialStream, empty};
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional, check_value};
use crate::units::{CalcContext, UnitOutput};

const DRY_AIR_CP_J_KG_K: f64 = 1005.0;

/// The process-side thermal specification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AirCoolerThermal {
    /// Process outlet temperature in K.
    OutletTemperature(Option<f64>),
    /// Signed process duty in W, zero or negative for cooling.
    Duty(Option<f64>),
}

impl Default for AirCoolerThermal {
    fn default() -> Self {
        Self::OutletTemperature(None)
    }
}

/// Steady-state air cooler with prescribed dry-air flow and a constant air
/// heat capacity of 1005 J/(kg K).
///
/// With heat transfer, the process outlet must be at least as warm as ambient
/// and the exchanger air outlet at most as warm as the process inlet. These
/// are terminal bounds with zero minimum approach. Equipment sizing, UA and
/// performance correlations are outside the model.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AirCoolerSpec {
    /// Process-side thermal mode.
    pub thermal: AirCoolerThermal,
    /// Process pressure loss in Pa.
    pub pressure_drop_pa: f64,
    /// Ambient air inlet temperature in K.
    pub ambient_temperature_k: Option<f64>,
    /// Dry-air mass flow in kg/s.
    pub air_mass_flow_kg_s: Option<f64>,
    /// Assigned parasitic fan power in W, dissipated into the air after the exchanger.
    pub fan_power_w: f64,
}

impl AirCoolerSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        match self.thermal {
            AirCoolerThermal::OutletTemperature(t) => {
                check_optional(&mut issues, "outlet temperature", t, Sign::Positive);
            }
            AirCoolerThermal::Duty(Some(q)) if q.is_finite() && q > 0.0 => {
                issues.push(SpecIssue {
                    field: "duty",
                    problem: SpecProblem::ContradictsUnitKind,
                });
            }
            AirCoolerThermal::Duty(q) => check_optional(&mut issues, "duty", q, Sign::Any),
        }
        check_value(
            &mut issues,
            "pressure drop",
            self.pressure_drop_pa,
            Sign::NonNegative,
        );
        check_optional(
            &mut issues,
            "ambient temperature",
            self.ambient_temperature_k,
            Sign::Positive,
        );
        check_optional(
            &mut issues,
            "air mass flow",
            self.air_mass_flow_kg_s,
            Sign::Positive,
        );
        check_value(
            &mut issues,
            "fan power",
            self.fan_power_w,
            Sign::NonNegative,
        );
        issues
    }

    /// Discharge air temperature in K after fan heat. Pass the solved, signed
    /// process duty in W. Missing or invalid air inputs, invalid duty or a
    /// nonfinite discharge temperature yield `None`.
    pub fn air_outlet_temperature_k(&self, process_duty_w: f64) -> Option<f64> {
        self.air_temperature_k(process_duty_w, self.fan_power_w)
    }

    fn air_temperature_k(&self, process_duty_w: f64, fan_power_w: f64) -> Option<f64> {
        let ambient = self.ambient_temperature_k?;
        let air_mass_flow = self.air_mass_flow_kg_s?;
        if !ambient.is_finite()
            || ambient <= 0.0
            || !air_mass_flow.is_finite()
            || air_mass_flow <= 0.0
            || !fan_power_w.is_finite()
            || fan_power_w < 0.0
            || !process_duty_w.is_finite()
            || process_duty_w > 0.0
        {
            return None;
        }

        let capacity = air_mass_flow * DRY_AIR_CP_J_KG_K;
        let rise = |power: f64| {
            if capacity.is_finite() {
                power / capacity
            } else {
                (power / air_mass_flow) / DRY_AIR_CP_J_KG_K
            }
        };
        let temperature = ambient + rise(-process_duty_w) + rise(fan_power_w);
        temperature.is_finite().then_some(temperature)
    }
}

pub(super) fn calculate(
    spec: &AirCoolerSpec,
    inlet: &MaterialStream,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("air cooler {} {}", issue.field, issue.problem));
    }

    let flow = inlet.molar_flow_mol_s;
    let mut output = if flow == 0.0 {
        UnitOutput::new(
            vec![empty()],
            0.0,
            vec![
                "Zero process inlet flow. The process outlet is empty and process duty is zero."
                    .into(),
            ],
        )
    } else {
        if !flow.is_finite() || flow < 0.0 {
            return Err("air cooler needs a finite nonnegative process inlet flow".into());
        }
        let state = inlet
            .state
            .as_ref()
            .ok_or("air cooler process inlet has no thermodynamic state")?;
        let pressure = state.pressure_pa - spec.pressure_drop_pa;
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(
                "air cooler pressure drop leaves nonpositive absolute outlet pressure".into(),
            );
        }
        let z = &state.mole_fractions;
        let h_in = state.enthalpy();
        let (outlet, duty) = match spec.thermal {
            AirCoolerThermal::OutletTemperature(Some(t)) => {
                let outlet = ctx.tp(z, t, pressure)?;
                let enthalpy_change = outlet.enthalpy() - h_in;
                if enthalpy_change > 0.0 {
                    return Err(
                        "air cooler thermal specification would add heat to the process".into(),
                    );
                }
                let duty = flow * enthalpy_change;
                (outlet, duty)
            }
            AirCoolerThermal::Duty(Some(q)) if q == 0.0 && spec.pressure_drop_pa == 0.0 => {
                (state.clone(), 0.0)
            }
            AirCoolerThermal::Duty(Some(q)) => (ctx.ph(z, pressure, h_in + q / flow)?, q),
            _ => return Err("missing air cooler thermal specification".into()),
        };
        if !duty.is_finite() {
            return Err("air cooler process duty is nonfinite".into());
        }
        if duty < 0.0 {
            let ambient = spec
                .ambient_temperature_k
                .ok_or("missing ambient temperature")?;
            if outlet.temperature_k < ambient {
                return Err(
                    "air cooler process outlet temperature falls below ambient air temperature"
                        .into(),
                );
            }
            if state.temperature_k <= ambient {
                return Err(
                    "air cooler process inlet must be warmer than ambient air to transfer heat"
                        .into(),
                );
            }
            let air_outlet = spec
                .air_temperature_k(duty, 0.0)
                .ok_or("air cooler exchanger air outlet temperature is nonfinite")?;
            if air_outlet > state.temperature_k {
                return Err(
                    "air cooler exchanger air outlet temperature exceeds the process inlet limit"
                        .into(),
                );
            }
        }

        UnitOutput::new(
            vec![MaterialStream {
                molar_flow_mol_s: flow,
                state: Some(outlet),
            }],
            duty,
            vec![],
        )
    };

    spec.air_outlet_temperature_k(output.duty_w)
        .ok_or("air cooler discharge air temperature is nonfinite")?;
    output.energy = vec![Some(spec.fan_power_w)];
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::balance;
    use property_packages::CubicPackage;
    use thermo_data::ComponentDatabase;
    use thermo_models::{CubicModel, CubicVariant};
    use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

    fn package() -> CubicPackage {
        let composition = Composition::new(
            CompositionBasis::Mole,
            vec![
                (ComponentId::new("methane").unwrap(), 0.8),
                (ComponentId::new("nitrogen").unwrap(), 0.2),
            ],
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            CubicVariant::PengRobinson,
            vec![vec![0.0; 2]; 2],
            "synthetic stream with zero binary interaction for air-cooler tests",
        )
        .unwrap();
        CubicPackage::new(model, Default::default())
    }

    fn process_inlet(ctx: &CalcContext) -> MaterialStream {
        MaterialStream {
            molar_flow_mol_s: 40.0,
            state: Some(ctx.tp(&[0.8, 0.2], 420.0, 2.0e6).unwrap()),
        }
    }

    fn spec() -> AirCoolerSpec {
        AirCoolerSpec {
            thermal: AirCoolerThermal::Duty(Some(-50_250.0)),
            pressure_drop_pa: 50_000.0,
            ambient_temperature_k: Some(300.0),
            air_mass_flow_kg_s: Some(2.0),
            fan_power_w: 10_050.0,
        }
    }

    #[test]
    fn sensible_air_reference_and_process_balance_in_both_modes() {
        let pkg = package();
        let ctx = CalcContext::new(&pkg, None);
        let inlet = process_inlet(&ctx);
        let mut spec = spec();
        let ph = calculate(&spec, &inlet, &ctx).unwrap();
        let state = ph.outlets[0].state.as_ref().unwrap();
        spec.thermal = AirCoolerThermal::OutletTemperature(Some(state.temperature_k));
        let tp = calculate(&spec, &inlet, &ctx).unwrap();

        for output in [ph, tp] {
            // With the assigned constant cp, 50,250 W into 2 kg/s of dry air
            // raises its temperature 25 K by hand calculation. The 10,050 W
            // downstream fan adds 5 K to the air that entered at 300 K.
            assert!((output.duty_w + 50_250.0).abs() < 1e-4);
            assert!((spec.air_outlet_temperature_k(output.duty_w).unwrap() - 330.0).abs() < 1e-7);
            assert_eq!(output.energy, [Some(10_050.0)]);
            let state = output.outlets[0].state.as_ref().unwrap();
            assert!((state.pressure_pa - 1.95e6).abs() < 1e-6);
            assert!(state.temperature_k < 420.0);
            balance(
                std::slice::from_ref(&inlet),
                &output.outlets,
                output.duty_w,
                None,
                2,
            )
            .unwrap();
        }
    }

    #[test]
    fn zero_process_flow_preserves_fan_reporting() {
        let pkg = package();
        let ctx = CalcContext::new(&pkg, None);
        let spec = spec();
        let output = calculate(&spec, &empty(), &ctx).unwrap();

        assert!(output.outlets[0].molar_flow_mol_s.abs() < f64::EPSILON);
        assert!(output.outlets[0].state.is_none());
        assert!(output.duty_w.abs() < f64::EPSILON);
        assert_eq!(output.energy, [Some(10_050.0)]);
        assert!((spec.air_outlet_temperature_k(output.duty_w).unwrap() - 305.0).abs() < 1e-12);
    }

    #[test]
    fn zero_heat_transfer_allows_isenthalpic_pressure_drop_below_ambient() {
        let pkg = package();
        let ctx = CalcContext::new(&pkg, None);
        let inlet = process_inlet(&ctx);
        let spec = AirCoolerSpec {
            thermal: AirCoolerThermal::Duty(Some(0.0)),
            ambient_temperature_k: Some(450.0),
            ..spec()
        };
        let output = calculate(&spec, &inlet, &ctx).unwrap();

        let state = output.outlets[0].state.as_ref().unwrap();
        assert!((state.pressure_pa - 1.95e6).abs() < 1e-6);
        assert!(state.temperature_k < 450.0);
        assert!((output.outlets[0].enthalpy_flow_w() - inlet.enthalpy_flow_w()).abs() < 1e-4);
        assert!(output.duty_w.abs() < f64::EPSILON);
    }

    #[test]
    fn rejects_process_heating_and_nonpositive_outlet_pressure() {
        let pkg = package();
        let ctx = CalcContext::new(&pkg, None);
        let inlet = process_inlet(&ctx);
        let mut spec = spec();
        spec.thermal = AirCoolerThermal::OutletTemperature(Some(440.0));
        assert!(
            calculate(&spec, &inlet, &ctx)
                .unwrap_err()
                .contains("add heat")
        );
        spec.pressure_drop_pa = 2.0e6;
        assert!(
            calculate(&spec, &inlet, &ctx)
                .unwrap_err()
                .contains("outlet pressure")
        );
    }

    #[test]
    fn enforces_terminal_bounds_before_fan_heat() {
        let pkg = package();
        let ctx = CalcContext::new(&pkg, None);
        let inlet = process_inlet(&ctx);
        let mut spec = spec();
        spec.thermal = AirCoolerThermal::OutletTemperature(Some(299.0));
        assert!(
            calculate(&spec, &inlet, &ctx)
                .unwrap_err()
                .contains("below ambient")
        );
        spec.thermal = AirCoolerThermal::OutletTemperature(Some(300.0));
        assert!(calculate(&spec, &inlet, &ctx).is_ok());

        spec.thermal = AirCoolerThermal::Duty(Some(-50_250.0));
        spec.ambient_temperature_k = Some(410.0);
        assert!(
            calculate(&spec, &inlet, &ctx)
                .unwrap_err()
                .contains("below ambient")
        );
        spec.ambient_temperature_k = Some(300.0);
        spec.air_mass_flow_kg_s = Some(0.4);
        assert!(
            calculate(&spec, &inlet, &ctx)
                .unwrap_err()
                .contains("process inlet limit")
        );
        spec.thermal = AirCoolerThermal::Duty(Some(-60_300.0));
        spec.air_mass_flow_kg_s = Some(0.5);
        assert!(calculate(&spec, &inlet, &ctx).is_ok());

        spec.thermal = AirCoolerThermal::Duty(Some(-50_250.0));
        spec.air_mass_flow_kg_s = Some(2.0);
        spec.fan_power_w = 241_200.0;
        let output = calculate(&spec, &inlet, &ctx).unwrap();
        assert!((spec.air_outlet_temperature_k(output.duty_w).unwrap() - 445.0).abs() < 1e-12);
    }

    #[test]
    fn validates_required_values_and_input_ranges() {
        let default = AirCoolerSpec::default();
        assert_eq!(default.thermal, AirCoolerThermal::OutletTemperature(None));
        assert!(default.fan_power_w.abs() < f64::EPSILON);
        assert_eq!(
            default
                .issues()
                .iter()
                .map(|issue| issue.problem)
                .collect::<Vec<_>>(),
            [SpecProblem::Missing; 3]
        );

        let mut spec = spec();
        spec.thermal = AirCoolerThermal::Duty(Some(1.0));
        spec.pressure_drop_pa = f64::INFINITY;
        spec.ambient_temperature_k = Some(0.0);
        spec.air_mass_flow_kg_s = Some(f64::NAN);
        spec.fan_power_w = -1.0;
        assert_eq!(
            spec.issues()
                .iter()
                .map(|issue| issue.problem)
                .collect::<Vec<_>>(),
            [
                SpecProblem::ContradictsUnitKind,
                SpecProblem::NotFinite,
                SpecProblem::NotPositive,
                SpecProblem::NotFinite,
                SpecProblem::Negative,
            ]
        );
        assert!(spec.air_outlet_temperature_k(-50_250.0).is_none());
    }
}
