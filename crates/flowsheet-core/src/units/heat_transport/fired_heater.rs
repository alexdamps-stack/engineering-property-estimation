//! Duty-only fired heater with a single process stream (#180).

use crate::engine::{MaterialStream, empty};
use crate::unit::{Sign, SpecIssue, SpecProblem, ThermalSpec, check_optional, check_value};
use crate::units::{CalcContext, UnitOutput};

/// Steady-state fired heater specified by its net heat input to the process.
/// The model excludes combustion, fuel and flue gas, efficiency, and furnace design.
#[derive(Debug, Clone, PartialEq)]
pub struct FiredHeaterSpec {
    /// Outlet temperature (K), nonnegative duty (W), or duty from the energy inlet.
    pub thermal: ThermalSpec,
    /// Pressure loss in Pa, subtracted from the inlet absolute pressure.
    pub pressure_drop_pa: f64,
}

impl Default for FiredHeaterSpec {
    fn default() -> Self {
        Self {
            thermal: ThermalSpec::OutletTemperature(None),
            pressure_drop_pa: 0.0,
        }
    }
}

impl FiredHeaterSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        match self.thermal {
            ThermalSpec::OutletTemperature(t) => {
                check_optional(&mut issues, "outlet temperature", t, Sign::Positive);
            }
            ThermalSpec::Duty(Some(q)) if q.is_finite() && q < 0.0 => issues.push(SpecIssue {
                field: "duty",
                problem: SpecProblem::ContradictsUnitKind,
            }),
            ThermalSpec::Duty(q) => check_optional(&mut issues, "duty", q, Sign::Any),
            ThermalSpec::EnergyStream => {}
        }
        check_value(
            &mut issues,
            "pressure drop",
            self.pressure_drop_pa,
            Sign::NonNegative,
        );
        issues
    }
}

pub(super) fn calculate(
    spec: &FiredHeaterSpec,
    inlet: &MaterialStream,
    energy: &[Option<f64>],
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("{} {}", issue.field, issue.problem));
    }
    let flow = inlet.molar_flow_mol_s;
    if !flow.is_finite() || flow < 0.0 {
        return Err("fired heater requires finite, nonnegative inlet flow".into());
    }

    let delivered_duty = energy.first().copied().flatten();
    let specified_duty = match spec.thermal {
        ThermalSpec::OutletTemperature(_) => None,
        ThermalSpec::Duty(q) => q,
        ThermalSpec::EnergyStream => {
            let q = delivered_duty.ok_or("no duty is supplied by an energy stream")?;
            if !q.is_finite() || q < 0.0 {
                return Err("fired heater requires finite, nonnegative delivered duty".into());
            }
            Some(q)
        }
    };
    let inlet_and_outlet_pressure = inlet
        .state
        .as_ref()
        .map(|state| {
            let p = state.pressure_pa - spec.pressure_drop_pa;
            if !p.is_finite() || p <= 0.0 {
                return Err("fired heater requires a finite, positive absolute outlet pressure");
            }
            Ok((state, p))
        })
        .transpose()?;

    let (outlet, duty, warnings) = if flow == 0.0 {
        if specified_duty.is_some_and(|q| q != 0.0) || delivered_duty.is_some_and(|q| q != 0.0) {
            return Err("fired heater cannot absorb nonzero duty at zero inlet flow".into());
        }
        (
            empty(),
            0.0,
            vec![
                "Zero inlet flow: the fired heater outlet is empty and no duty is applied.".into(),
            ],
        )
    } else {
        let (state, p) =
            inlet_and_outlet_pressure.ok_or("fired heater inlet has no thermodynamic state")?;
        let (out, q) = match (spec.thermal, specified_duty) {
            // Preserve the inlet phase fractions when neither pressure nor enthalpy changes.
            (_, Some(0.0)) if spec.pressure_drop_pa == 0.0 => (state.clone(), 0.0),
            (_, Some(q)) => (
                ctx.ph(&state.mole_fractions, p, state.enthalpy() + q / flow)?,
                q,
            ),
            (ThermalSpec::OutletTemperature(Some(t)), None) => {
                let out = ctx.tp(&state.mole_fractions, t, p)?;
                let q = flow * (out.enthalpy() - state.enthalpy());
                if !q.is_finite() {
                    return Err("fired heater outlet temperature implies nonfinite duty".into());
                }
                if q < 0.0 {
                    return Err("fired heater outlet temperature requires heat removal".into());
                }
                (out, q)
            }
            _ => return Err("missing fired heater thermal specification".into()),
        };
        (
            MaterialStream {
                molar_flow_mol_s: flow,
                state: Some(out),
            },
            q,
            vec![],
        )
    };
    let mut result = UnitOutput::new(vec![outlet], duty, warnings);
    if spec.thermal == ThermalSpec::EnergyStream {
        result.energy = vec![Some(duty)];
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use property_packages::CubicPackage;
    use thermo_data::ComponentDatabase;
    use thermo_models::{CubicModel, CubicVariant};
    use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

    const Z: [f64; 2] = [0.4, 0.6];

    fn package() -> CubicPackage {
        let composition = Composition::new(
            CompositionBasis::Mole,
            ["methane", "nitrogen"]
                .into_iter()
                .zip(Z)
                .map(|(id, x)| (ComponentId::new(id).unwrap(), x))
                .collect(),
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            CubicVariant::PengRobinson,
            vec![vec![0.0; 2]; 2],
            "zero interactions for the dilute-gas heater test",
        )
        .unwrap();
        CubicPackage::new(model, Default::default())
    }

    fn inlet(ctx: &CalcContext) -> MaterialStream {
        MaterialStream {
            molar_flow_mol_s: 4.0,
            state: Some(ctx.tp(&Z, 300.0, 1.0).unwrap()),
        }
    }

    #[test]
    fn thermal_modes_close_analytic_balances() {
        let package = package();
        let ctx = CalcContext::new(&package, None);
        let inlet = inlet(&ctx);
        // Independent ideal-gas limit from 4 mol/s times the mole-weighted
        // integral of the NASA7 cp polynomials from 300 K to 400 K, using
        // data/source-extracts/components.provisional.json (2026-11-07.1).
        // The 0.01 W and 0.001 K budgets allow PR departures over 1 to 0.5 Pa.
        let ideal_duty = 13_084.355_399_572_9;
        for (thermal, energy) in [
            (ThermalSpec::OutletTemperature(Some(400.0)), None),
            (ThermalSpec::Duty(Some(ideal_duty)), None),
            (ThermalSpec::EnergyStream, Some(ideal_duty)),
        ] {
            let spec = FiredHeaterSpec {
                thermal,
                pressure_drop_pa: 0.5,
            };
            let result = calculate(&spec, &inlet, &[energy], &ctx).unwrap();
            let [outlet] = result.outlets.as_slice() else {
                panic!("expected one process outlet");
            };
            let state = outlet.state.as_ref().unwrap();
            assert_eq!(state.pressure_pa, 0.5);
            assert!((state.temperature_k - 400.0).abs() < 0.001);
            assert_eq!(outlet.molar_flow_mol_s, 4.0);
            for (actual, expected) in outlet.component_flows_mol_s(2).into_iter().zip([1.6, 2.4]) {
                assert!((actual - expected).abs() < 1e-12);
            }
            assert!((result.duty_w - ideal_duty).abs() < 0.01);
            assert!(
                (outlet.enthalpy_flow_w() - inlet.enthalpy_flow_w() - result.duty_w).abs() < 1e-5
            );
            if let Some(q) = energy {
                assert_eq!(result.energy, [Some(q)]);
            } else {
                assert!(result.energy.is_empty());
            }
        }
    }

    #[test]
    fn zero_duty_with_pressure_drop_is_isenthalpic() {
        let package = package();
        let ctx = CalcContext::new(&package, None);
        let inlet = MaterialStream {
            molar_flow_mol_s: 4.0,
            state: Some(ctx.tp(&Z, 300.0, 1e6).unwrap()),
        };
        let spec = FiredHeaterSpec {
            thermal: ThermalSpec::Duty(Some(0.0)),
            pressure_drop_pa: 1e5,
        };
        let result = calculate(&spec, &inlet, &[], &ctx).unwrap();
        assert_eq!(result.duty_w, 0.0);
        assert_eq!(result.outlets[0].state.as_ref().unwrap().pressure_pa, 9e5);
        assert!((result.outlets[0].enthalpy_flow_w() - inlet.enthalpy_flow_w()).abs() < 1e-5);
    }

    #[test]
    fn invalid_specs_identify_the_field_and_problem() {
        let default = FiredHeaterSpec::default();
        assert_eq!(default.pressure_drop_pa, 0.0);
        assert_eq!(default.thermal, ThermalSpec::OutletTemperature(None));
        for (thermal, field, problem) in [
            (default.thermal, "outlet temperature", SpecProblem::Missing),
            (
                ThermalSpec::OutletTemperature(Some(0.0)),
                "outlet temperature",
                SpecProblem::NotPositive,
            ),
            (ThermalSpec::Duty(None), "duty", SpecProblem::Missing),
            (
                ThermalSpec::Duty(Some(-1.0)),
                "duty",
                SpecProblem::ContradictsUnitKind,
            ),
            (
                ThermalSpec::Duty(Some(f64::NAN)),
                "duty",
                SpecProblem::NotFinite,
            ),
        ] {
            assert_eq!(
                FiredHeaterSpec { thermal, ..default }.issues(),
                [SpecIssue { field, problem }]
            );
        }
        for (pressure_drop_pa, problem) in [
            (-1.0, SpecProblem::Negative),
            (f64::INFINITY, SpecProblem::NotFinite),
        ] {
            let spec = FiredHeaterSpec {
                thermal: ThermalSpec::EnergyStream,
                pressure_drop_pa,
            };
            assert_eq!(
                spec.issues(),
                [SpecIssue {
                    field: "pressure drop",
                    problem
                }]
            );
        }
    }

    #[test]
    fn rejects_heat_removal_and_invalid_outlet_pressure() {
        let package = package();
        let ctx = CalcContext::new(&package, None);
        let inlet = inlet(&ctx);
        let temperature = FiredHeaterSpec {
            thermal: ThermalSpec::OutletTemperature(Some(290.0)),
            ..Default::default()
        };
        let error = calculate(&temperature, &inlet, &[], &ctx).unwrap_err();
        assert!(error.contains("requires heat removal"), "{error}");
        for pressure_drop_pa in [1.0, 2.0] {
            let spec = FiredHeaterSpec {
                pressure_drop_pa,
                ..temperature
            };
            let error = calculate(&spec, &inlet, &[], &ctx).unwrap_err();
            assert!(error.contains("absolute outlet pressure"), "{error}");
        }
    }

    #[test]
    fn rejects_missing_or_invalid_delivered_duty() {
        let package = package();
        let ctx = CalcContext::new(&package, None);
        let inlet = inlet(&ctx);
        let spec = FiredHeaterSpec {
            thermal: ThermalSpec::EnergyStream,
            ..Default::default()
        };
        let error = calculate(&spec, &inlet, &[None], &ctx).unwrap_err();
        assert!(error.contains("no duty is supplied"), "{error}");
        for q in [-1.0, f64::INFINITY, f64::NAN] {
            let error = calculate(&spec, &inlet, &[Some(q)], &ctx).unwrap_err();
            assert!(
                error.contains("finite, nonnegative delivered duty"),
                "{error}"
            );
        }
    }

    #[test]
    fn zero_flow_rejects_nonzero_heat_and_accepts_zero_duty() {
        let package = package();
        let ctx = CalcContext::new(&package, None);
        let inlet = empty();
        for (thermal, energy) in [
            (ThermalSpec::Duty(Some(1.0)), None),
            (ThermalSpec::EnergyStream, Some(1.0)),
        ] {
            let spec = FiredHeaterSpec {
                thermal,
                ..Default::default()
            };
            let error = calculate(&spec, &inlet, &[energy], &ctx).unwrap_err();
            assert!(error.contains("nonzero duty at zero inlet flow"), "{error}");
        }
        for (thermal, energy) in [
            (ThermalSpec::OutletTemperature(Some(400.0)), None),
            (ThermalSpec::Duty(Some(0.0)), None),
            (ThermalSpec::EnergyStream, Some(0.0)),
        ] {
            let spec = FiredHeaterSpec {
                thermal,
                ..Default::default()
            };
            let result = calculate(&spec, &inlet, &[energy], &ctx).unwrap();
            assert_eq!(result.duty_w, 0.0);
            assert_eq!(result.outlets[0].molar_flow_mol_s, 0.0);
            assert!(result.outlets[0].state.is_none());
            assert_eq!(result.warnings.len(), 1);
        }
    }
}
