//! Isothermal component recovery splitter (#172).

use std::collections::BTreeMap;

use crate::engine::{MaterialStream, empty};
use crate::unit::{Sign, SpecIssue, SpecProblem, check_value};
use crate::units::{CalcContext, UnitOutput};

/// Component splitter at inlet temperature and pressure. The energy duty is
/// calculated from the inlet and outlet enthalpy flows.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ComponentSplitterSpec {
    /// Fraction of each inlet component sent to overhead, keyed by the
    /// property package's component identifier. Bottoms receives the rest.
    pub overhead_recovery: BTreeMap<String, f64>,
}

impl ComponentSplitterSpec {
    pub(super) fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        if self.overhead_recovery.is_empty() {
            issues.push(SpecIssue {
                field: "overhead recoveries",
                problem: SpecProblem::Missing,
            });
        }
        for &recovery in self.overhead_recovery.values() {
            check_value(&mut issues, "overhead recovery", recovery, Sign::Fraction);
        }
        issues
    }

    /// Check that every recovery belongs to the selected package and every
    /// package component has a recovery, including components absent at the inlet.
    pub fn check_components(&self, components: &[&str]) -> Result<(), String> {
        if let Some(id) = components
            .iter()
            .find(|&&id| !self.overhead_recovery.contains_key(id))
        {
            return Err(format!("missing overhead recovery for {id}"));
        }
        if let Some(id) = self
            .overhead_recovery
            .keys()
            .find(|id| !components.contains(&id.as_str()))
        {
            return Err(format!(
                "overhead recovery for {id} is outside the property package"
            ));
        }
        Ok(())
    }
}

pub(super) fn calculate(
    spec: &ComponentSplitterSpec,
    inlet: &MaterialStream,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    let components = ctx.component_ids();
    spec.check_components(&components)?;
    let inlet_flow = inlet.molar_flow_mol_s;
    if inlet_flow == 0.0 {
        return Ok(UnitOutput::new(vec![empty(), empty()], 0.0, vec![]));
    }
    let state = inlet
        .state
        .as_ref()
        .ok_or("component splitter inlet has no state")?;
    let recoveries: Vec<_> = components
        .iter()
        .map(|&id| spec.overhead_recovery[id])
        .collect();
    let mut outlets = Vec::with_capacity(2);
    for overhead in [true, false] {
        let mut amounts: Vec<_> = state
            .mole_fractions
            .iter()
            .zip(&recoveries)
            .map(|(&z, &r)| z * if overhead { r } else { 1.0 - r })
            .collect();
        let share: f64 = amounts.iter().sum();
        let flow = inlet_flow * share;
        let outlet = if flow == 0.0 {
            empty()
        } else {
            for x in &mut amounts {
                *x /= share;
            }
            let outlet_state = if amounts == state.mole_fractions {
                state.clone()
            } else {
                ctx.tp(&amounts, state.temperature_k, state.pressure_pa)
                    .map_err(|error| {
                        let port = if overhead { "overhead" } else { "bottoms" };
                        format!("{port} outlet TP flash: {error}")
                    })?
            };
            MaterialStream {
                molar_flow_mol_s: flow,
                state: Some(outlet_state),
            }
        };
        outlets.push(outlet);
    }
    let duty = outlets
        .iter()
        .map(MaterialStream::enthalpy_flow_w)
        .sum::<f64>()
        - inlet.enthalpy_flow_w();
    Ok(UnitOutput::new(outlets, duty, vec![]))
}
