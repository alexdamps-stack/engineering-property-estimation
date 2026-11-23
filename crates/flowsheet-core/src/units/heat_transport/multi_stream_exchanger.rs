//! Steady-state multi-stream exchanger with a shared countercurrent core (#181).

use crate::engine::{MaterialStream, empty};
use crate::unit::{Sign, SpecIssue, SpecProblem, check_optional, check_value};
use crate::units::{CalcContext, UnitOutput};

/// Largest number of passes of a multi-stream exchanger.
pub const MAX_PASSES: u16 = 12;

/// Heat transfer direction of one material pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExchangerPassDirection {
    /// Duty into the pass is zero or negative.
    Hot,
    /// Duty into the pass is zero or positive.
    Cold,
}

impl ExchangerPassDirection {
    fn allows_duty(self, duty: f64) -> bool {
        match self {
            Self::Hot => duty <= 0.0,
            Self::Cold => duty >= 0.0,
        }
    }
}

/// Thermal specification of one material pass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExchangerPassThermal {
    /// Specified outlet temperature in K.
    OutletTemperature(Option<f64>),
    /// Specified duty in W, positive into the pass.
    Duty(Option<f64>),
    /// Determine the duty from the other passes by closing the heat balance.
    Balance,
}

/// Specifications for one independent inlet/outlet pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExchangerPassSpec {
    /// Whether the pass supplies or receives heat.
    pub direction: ExchangerPassDirection,
    /// Temperature, duty, or an unknown duty solved by the energy balance.
    pub thermal: ExchangerPassThermal,
    /// Pressure loss in Pa, subtracted from the inlet absolute pressure.
    pub pressure_drop_pa: f64,
}

/// Adiabatic exchanger with 2 to [`MAX_PASSES`] independent material passes.
/// Each numbered inlet feeds its matching outlet without mixing.
/// The minimum approach applies at both terminals of the shared countercurrent
/// core. Interior temperature profiles, pinch analysis, and rating are outside
/// the model.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiStreamExchangerSpec {
    /// Pass specifications in port order.
    pub passes: Vec<ExchangerPassSpec>,
    /// Minimum permitted hot-minus-cold terminal temperature difference in K.
    pub minimum_approach_k: f64,
}

impl Default for MultiStreamExchangerSpec {
    fn default() -> Self {
        Self::with_pass_count(3)
    }
}

impl MultiStreamExchangerSpec {
    /// Initialize numbered passes with one hot pass and the last pass balancing.
    /// The pass count is checked during specification validation.
    pub fn with_pass_count(count: u16) -> Self {
        let passes = (0..count)
            .map(|i| ExchangerPassSpec {
                direction: if i == 0 {
                    ExchangerPassDirection::Hot
                } else {
                    ExchangerPassDirection::Cold
                },
                thermal: if i + 1 == count {
                    ExchangerPassThermal::Balance
                } else {
                    ExchangerPassThermal::OutletTemperature(None)
                },
                pressure_drop_pa: 0.0,
            })
            .collect();
        Self {
            passes,
            minimum_approach_k: 0.0,
        }
    }

    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        use ExchangerPassDirection::{Cold, Hot};
        use ExchangerPassThermal::{Balance, Duty, OutletTemperature};

        if !(2..=usize::from(MAX_PASSES)).contains(&self.passes.len()) {
            return vec![SpecIssue {
                field: "passes",
                problem: SpecProblem::CountOutOfRange,
            }];
        }

        let mut issues = Vec::new();
        check_value(
            &mut issues,
            "minimum approach",
            self.minimum_approach_k,
            Sign::NonNegative,
        );
        if self.passes.iter().filter(|p| p.thermal == Balance).count() != 1 {
            issues.push(SpecIssue {
                field: "balancing pass count (exactly one required)",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        for (direction, field) in [(Hot, "hot pass"), (Cold, "cold pass")] {
            if !self.passes.iter().any(|p| p.direction == direction) {
                issues.push(SpecIssue {
                    field,
                    problem: SpecProblem::Missing,
                });
            }
        }
        for (pass, &(temperature, duty, drop)) in self.passes.iter().zip(&PASS_FIELDS) {
            match pass.thermal {
                OutletTemperature(t) => check_optional(&mut issues, temperature, t, Sign::Positive),
                Duty(Some(q)) if q.is_finite() && !pass.direction.allows_duty(q) => {
                    issues.push(SpecIssue {
                        field: duty,
                        problem: SpecProblem::ContradictsUnitKind,
                    });
                }
                Duty(q) => check_optional(&mut issues, duty, q, Sign::Any),
                Balance => {}
            }
            check_value(&mut issues, drop, pass.pressure_drop_pa, Sign::NonNegative);
        }
        issues
    }
}

// SpecIssue stores static field names; include the same one-based indices as the ports.
macro_rules! pass_fields {
    ($($number:literal),+) => {
        [$(
            (
                concat!("pass ", $number, " outlet temperature"),
                concat!("pass ", $number, " duty"),
                concat!("pass ", $number, " pressure drop"),
            ),
        )+]
    };
}
const PASS_FIELDS: [(&str, &str, &str); MAX_PASSES as usize] =
    pass_fields!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);

pub(super) fn calculate(
    spec: &MultiStreamExchangerSpec,
    inlets: &[MaterialStream],
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    if let Some(issue) = spec.issues().first() {
        return Err(format!("{} {}", issue.field, issue.problem));
    }
    if inlets.len() != spec.passes.len() {
        return Err("multi-stream exchanger requires one inlet per pass".into());
    }

    let balance = spec
        .passes
        .iter()
        .position(|p| p.thermal == ExchangerPassThermal::Balance)
        .ok_or("multi-stream exchanger requires a balancing pass")?;
    // Solve the balancing pass last while keeping each outlet in port order.
    let order = (0..inlets.len())
        .filter(|&i| i != balance)
        .chain(std::iter::once(balance));
    let mut outlets = vec![empty(); inlets.len()];
    let mut warnings = Vec::new();
    let mut net_duty = 0.0;
    let mut transferred_duty = 0.0;
    for i in order {
        let pass = spec.passes[i];
        let (outlet, duty) = calculate_pass(pass, &inlets[i], -net_duty, ctx)
            .map_err(|error| format!("pass {} {error}", i + 1))?;
        if outlet.state.is_none() {
            warnings.push(format!(
                "Zero inlet flow on pass {}. The pass was bypassed with an empty outlet.",
                i + 1
            ));
        }
        net_duty += duty;
        if pass.direction == ExchangerPassDirection::Hot {
            transferred_duty -= duty;
        }
        if !net_duty.is_finite() || !transferred_duty.is_finite() {
            return Err(format!(
                "pass {} causes a nonfinite total exchanger duty",
                i + 1
            ));
        }
        outlets[i] = outlet;
    }

    check_terminals(spec, inlets, &outlets)?;
    let mut result = UnitOutput::new(outlets, 0.0, warnings);
    result.transferred_duty_w = Some(transferred_duty);
    Ok(result)
}

fn calculate_pass(
    pass: ExchangerPassSpec,
    inlet: &MaterialStream,
    balance_duty: f64,
    ctx: &CalcContext,
) -> Result<(MaterialStream, f64), String> {
    use ExchangerPassThermal::{Balance, Duty, OutletTemperature};

    let flow = inlet.molar_flow_mol_s;
    if !flow.is_finite() || flow < 0.0 {
        return Err("requires finite, nonnegative inlet flow".into());
    }
    if flow == 0.0 {
        return if matches!(pass.thermal, Duty(Some(0.0))) {
            Ok((empty(), 0.0))
        } else {
            Err("requires an explicitly specified zero duty at zero inlet flow".into())
        };
    }

    let state = inlet
        .state
        .as_ref()
        .ok_or("inlet has no thermodynamic state")?;
    let pressure = state.pressure_pa - pass.pressure_drop_pa;
    if !pressure.is_finite() || pressure <= 0.0 {
        return Err(
            "requires a finite, positive absolute outlet pressure after pressure drop".into(),
        );
    }
    let h_in = state.enthalpy();
    if !state.temperature_k.is_finite() || !h_in.is_finite() {
        return Err("requires finite inlet temperature and enthalpy".into());
    }
    let z = &state.mole_fractions;
    let thermal = match pass.thermal {
        Balance => Duty(Some(balance_duty)),
        specified => specified,
    };
    let (outlet, duty) = match thermal {
        OutletTemperature(Some(t)) => {
            let outlet = ctx
                .tp(z, t, pressure)
                .map_err(|e| format!("TP flash failed ({e})"))?;
            let duty = flow * (outlet.enthalpy() - h_in);
            check_duty(pass, duty)?;
            (outlet, duty)
        }
        Duty(Some(duty)) => {
            check_duty(pass, duty)?;
            // Preserve the equilibrium state when both pressure and enthalpy are unchanged.
            let outlet = if duty == 0.0 && pass.pressure_drop_pa == 0.0 {
                state.clone()
            } else {
                ctx.ph(z, pressure, h_in + duty / flow)
                    .map_err(|e| format!("PH flash failed ({e})"))?
            };
            (outlet, duty)
        }
        _ => return Err("is missing a thermal specification".into()),
    };
    Ok((
        MaterialStream {
            molar_flow_mol_s: flow,
            state: Some(outlet),
        },
        duty,
    ))
}

fn check_duty(pass: ExchangerPassSpec, duty: f64) -> Result<(), String> {
    let label = if pass.thermal == ExchangerPassThermal::Balance {
        "balancing duty"
    } else {
        "duty"
    };
    if !duty.is_finite() {
        return Err(format!("{label} is not finite"));
    }
    if !pass.direction.allows_duty(duty) {
        return Err(match pass.direction {
            ExchangerPassDirection::Hot => format!("{label} must be nonpositive for a hot pass"),
            ExchangerPassDirection::Cold => format!("{label} must be nonnegative for a cold pass"),
        });
    }
    Ok(())
}

fn check_terminals(
    spec: &MultiStreamExchangerSpec,
    inlets: &[MaterialStream],
    outlets: &[MaterialStream],
) -> Result<(), String> {
    use ExchangerPassDirection::{Cold, Hot};

    let flowing = spec
        .passes
        .iter()
        .zip(inlets)
        .zip(outlets)
        .enumerate()
        .filter_map(|(i, ((pass, inlet), outlet))| {
            Some((
                i + 1,
                pass.direction,
                inlet.state.as_ref()?.temperature_k,
                outlet.state.as_ref()?.temperature_k,
            ))
        });
    let min = spec.minimum_approach_k;
    for (hot, _, hot_in, hot_out) in flowing.clone().filter(|p| p.1 == Hot) {
        for (cold, _, cold_in, cold_out) in flowing.clone().filter(|p| p.1 == Cold) {
            if hot_in - cold_out < min || hot_out - cold_in < min {
                return Err(format!(
                    "hot pass {hot} and cold pass {cold} violate the minimum terminal approach of {min} K"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use property_packages::CubicPackage;
    use thermo_data::ComponentDatabase;
    use thermo_models::{CubicModel, CubicVariant};
    use thermo_types::{ComponentId, Composition, CompositionBasis, NormalizationPolicy};

    #[test]
    fn defaults_report_missing_temperature_by_pass_number() {
        let spec = MultiStreamExchangerSpec::default();
        assert_eq!(spec.passes.len(), 3);
        assert_eq!(spec.passes[0].direction, ExchangerPassDirection::Hot);
        assert_eq!(spec.passes[2].direction, ExchangerPassDirection::Cold);
        assert_eq!(spec.passes[2].thermal, ExchangerPassThermal::Balance);
        let issues = spec.issues();
        assert_eq!(
            issues.iter().map(|issue| issue.field).collect::<Vec<_>>(),
            ["pass 1 outlet temperature", "pass 2 outlet temperature"]
        );
        assert!(
            issues
                .iter()
                .all(|issue| issue.problem == SpecProblem::Missing)
        );

        let mut largest = MultiStreamExchangerSpec::with_pass_count(MAX_PASSES);
        largest.passes.last_mut().unwrap().pressure_drop_pa = f64::INFINITY;
        assert_eq!(
            largest.issues().last(),
            Some(&SpecIssue {
                field: "pass 12 pressure drop",
                problem: SpecProblem::NotFinite,
            })
        );
    }

    #[test]
    fn checks_all_flowing_terminal_pairs_at_zero_duty() {
        use ExchangerPassDirection::{Cold, Hot};
        use ExchangerPassThermal::{Balance, Duty};

        let composition = Composition::new(
            CompositionBasis::Mole,
            vec![(ComponentId::new("methane").unwrap(), 1.0)],
            NormalizationPolicy::RequireUnitSum,
        )
        .unwrap();
        let model = CubicModel::new(
            &ComponentDatabase::embedded().unwrap(),
            &composition,
            CubicVariant::PengRobinson,
            vec![vec![0.0]],
            "pure methane",
        )
        .unwrap();
        let package = CubicPackage::new(model, Default::default());
        let ctx = CalcContext::new(&package, None);
        let mut inlets = [300.0, 350.0, 305.0, 325.0].map(|t| MaterialStream {
            molar_flow_mol_s: 1.0,
            state: Some(ctx.tp(&[1.0], t, 100_000.0).unwrap()),
        });
        let mut spec = MultiStreamExchangerSpec {
            passes: [
                (Cold, Balance),
                (Hot, Duty(Some(0.0))),
                (Cold, Duty(Some(0.0))),
                (Hot, Duty(Some(0.0))),
            ]
            .map(|(direction, thermal)| ExchangerPassSpec {
                direction,
                thermal,
                pressure_drop_pa: 0.0,
            })
            .to_vec(),
            minimum_approach_k: 20.0,
        };
        let result = calculate(&spec, &inlets, &ctx).unwrap();
        assert_eq!(result.duty_w, 0.0);
        assert_eq!(result.transferred_duty_w, Some(0.0));

        spec.minimum_approach_k = 21.0;
        let error = calculate(&spec, &inlets, &ctx).unwrap_err();
        assert!(error.contains("hot pass 4 and cold pass 3"), "{error}");

        // A splitter can supply a state at zero flow. The exchanger still bypasses it.
        inlets[3].molar_flow_mol_s = 0.0;
        let result = calculate(&spec, &inlets, &ctx).unwrap();
        assert!(result.outlets[3].state.is_none());
        assert!(result.warnings[0].contains("pass 4"));

        inlets[1].molar_flow_mol_s = 0.0;
        assert_eq!(
            calculate(&spec, &inlets, &ctx).unwrap().transferred_duty_w,
            Some(0.0)
        );
    }
}
