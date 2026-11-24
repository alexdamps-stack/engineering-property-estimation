//! Rotating group (A14, #182): pumps, compressors, expanders and the
//! multistage compressor with intercooling.
//!
//! Owned by the rotating work package (`docs/a17-architecture.md`). The
//! pump/compressor/expander specifications are the shared
//! [`UnitSpec::Pressurize`] and [`UnitSpec::Expand`] variants in `unit.rs`;
//! their calculation is routed here through `calculate_core`.
//!
//! Every machine uses the isentropic-efficiency model: a PS flash at the
//! outlet pressure and the inlet entropy gives the ideal enthalpy change
//! \(\Delta h_s\). Pumps and compressors take \(\Delta h = \Delta h_s/\eta\),
//! expanders \(\Delta h = \eta\,\Delta h_s\), and a PH flash at the outlet
//! pressure gives the actual outlet state. The shaft power
//! \(W = \dot n\,\Delta h\) is the unit's external duty (positive into the
//! fluid) and closes the energy balance.

use super::{CalcContext, GroupSpec, UnitInputs, UnitOutput, kind_info};
use crate::UnitSpec;
use crate::engine::{MaterialStream, empty};
use crate::unit::{
    CalculationStatus, DEFAULT_EFFICIENCY, KindInfo, PaletteCategory, PortSpec, PressureRise,
    PressureSpec, Sign, SpecIssue, SpecProblem, UnitKind, check_optional, check_value,
};
use property_packages::EquilibriumState;
use thermo_types::Phase;

const SINGLE: [PortSpec; 2] = [PortSpec::inlet("in"), PortSpec::outlet("out")];
const POWER_IN: [PortSpec; 1] = [PortSpec::energy_inlet("power")];
const POWER_OUT: [PortSpec; 1] = [PortSpec::energy_outlet("power")];
const MULTISTAGE_ENERGY: [PortSpec; 2] = [
    PortSpec::energy_inlet("power"),
    PortSpec::energy_outlet("intercooler duty"),
];

/// Largest number of compression stages in a multistage compressor.
pub const MAX_COMPRESSOR_STAGES: u16 = 8;

/// Metadata of a rotating kind.
pub(crate) fn info(kind: UnitKind) -> KindInfo {
    use PaletteCategory::PressureChange as P;
    let core = |label, tag_series, energy_ports| KindInfo {
        issue: None,
        ..kind_info(label, tag_series, P, &SINGLE, energy_ports, 0)
    };
    match kind {
        UnitKind::Pump => core("Pump", "P", &POWER_IN),
        UnitKind::Compressor => core("Compressor", "K", &POWER_IN),
        UnitKind::Expander => core("Expander", "EX", &POWER_OUT),
        UnitKind::MultistageCompressor => kind_info(
            "Multistage compressor",
            "K",
            P,
            &SINGLE,
            &MULTISTAGE_ENERGY,
            182,
        ),
        _ => unreachable!("{kind:?} is not a rotating kind"),
    }
}

/// Calculation status of a rotating kind: all four use the isentropic
/// efficiency model.
pub(crate) fn status(_kind: UnitKind) -> CalculationStatus {
    CalculationStatus::Calculated
}

/// Specifications of the rotating kinds other than pump, compressor and
/// expander.
#[derive(Debug, Clone, PartialEq)]
pub enum RotatingSpec {
    /// Multistage compressor with intercooling (#182).
    MultistageCompressor(MultistageCompressorSpec),
}

/// Multistage compressor with intercoolers between stages (#182).
///
/// The stages share one isentropic efficiency and an equal pressure ratio
/// \((p_{out}/p_{in})^{1/N}\). Each of the \(N-1\) intercoolers cools the
/// discharge to the same outlet temperature at constant pressure. There is no
/// aftercooler: the last stage discharges directly.
#[derive(Debug, Clone, PartialEq)]
pub struct MultistageCompressorSpec {
    /// Number of compression stages, 1 to [`MAX_COMPRESSOR_STAGES`].
    pub stages: u16,
    /// Absolute discharge pressure of the last stage, Pa.
    pub outlet_pressure_pa: Option<f64>,
    /// Isentropic efficiency of every stage, greater than 0 and at most 1.
    pub efficiency: f64,
    /// Intercooler outlet temperature, K. Required for two or more stages.
    pub intercooler_temperature_k: Option<f64>,
}

impl Default for MultistageCompressorSpec {
    fn default() -> Self {
        Self {
            stages: 2,
            outlet_pressure_pa: None,
            efficiency: DEFAULT_EFFICIENCY,
            intercooler_temperature_k: None,
        }
    }
}

impl MultistageCompressorSpec {
    fn issues(&self) -> Vec<SpecIssue> {
        let mut issues = Vec::new();
        if !(1..=MAX_COMPRESSOR_STAGES).contains(&self.stages) {
            issues.push(SpecIssue {
                field: "stages",
                problem: SpecProblem::CountOutOfRange,
            });
        }
        check_optional(
            &mut issues,
            "outlet pressure",
            self.outlet_pressure_pa,
            Sign::Positive,
        );
        check_value(&mut issues, "efficiency", self.efficiency, Sign::Efficiency);
        if self.stages > 1 {
            check_optional(
                &mut issues,
                "intercooler temperature",
                self.intercooler_temperature_k,
                Sign::Positive,
            );
        } else if let Some(t) = self.intercooler_temperature_k {
            check_value(&mut issues, "intercooler temperature", t, Sign::Positive);
        }
        issues
    }
}

/// Default specification of a rotating group kind.
pub(crate) fn default_spec(kind: UnitKind) -> RotatingSpec {
    match kind {
        UnitKind::MultistageCompressor => {
            RotatingSpec::MultistageCompressor(MultistageCompressorSpec::default())
        }
        _ => unreachable!("{kind:?} has no RotatingSpec"),
    }
}

impl GroupSpec for RotatingSpec {
    fn kind(&self) -> UnitKind {
        match self {
            RotatingSpec::MultistageCompressor(_) => UnitKind::MultistageCompressor,
        }
    }

    fn issues(&self) -> Vec<SpecIssue> {
        match self {
            RotatingSpec::MultistageCompressor(s) => s.issues(),
        }
    }
}

impl CalcContext<'_> {
    /// PS equilibrium at the specified absolute pressure and molar entropy.
    pub fn ps(
        &self,
        z: &[f64],
        p_pa: f64,
        s_j_mol_k: f64,
    ) -> Result<property_packages::EquilibriumState, String> {
        self.cubic_package()?
            .ps(z, p_pa, s_j_mol_k)
            .map_err(|e| e.to_string())
    }
}

/// One compression or expansion step of a solved machine.
#[derive(Debug, Clone, PartialEq)]
pub struct RotatingStage {
    /// Absolute inlet pressure, Pa.
    pub inlet_pressure_pa: f64,
    /// Absolute outlet pressure, Pa.
    pub outlet_pressure_pa: f64,
    /// Temperature of the isentropic outlet state, K.
    pub isentropic_outlet_temperature_k: f64,
    /// Actual outlet temperature, K.
    pub outlet_temperature_k: f64,
    /// Isentropic power \(\dot n\,\Delta h_s\), W, positive into the fluid.
    pub isentropic_power_w: f64,
    /// Actual shaft power \(\dot n\,\Delta h\), W, positive into the fluid.
    pub power_w: f64,
}

/// Solved stages of a pump, compressor, expander or multistage compressor.
#[derive(Debug, Clone, PartialEq)]
pub struct RotatingResult {
    /// Isentropic efficiency used for every stage.
    pub efficiency: f64,
    /// One entry per stage, in flow order. Empty at zero flow.
    pub stages: Vec<RotatingStage>,
    /// Heat removed by all intercoolers, W (nonnegative; zero for single
    /// machines).
    pub intercooler_duty_w: f64,
    /// Pumps with flow only: hydraulic quantities on the engineering liquid
    /// density (#228), or why they are unavailable.
    pub pump_hydraulics: Option<PumpHydraulics>,
}

/// Hydraulic quantities of a pump on the package's engineering liquid
/// density (#228). They are reported beside the shaft power, which comes
/// from the enthalpy change and closes the energy balance; they never
/// change the outlet state.
#[derive(Debug, Clone, PartialEq)]
pub enum PumpHydraulics {
    Calculated(PumpHydraulicValues),
    /// The package has no engineering liquid density, or the inlet is
    /// outside its domain or not liquid. Nothing is substituted.
    Unavailable(String),
}

/// Pump hydraulics at the inlet conditions.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpHydraulicValues {
    /// Engineering liquid density of the inlet, kg/m³.
    pub inlet_density_kg_m3: f64,
    /// Inlet volumetric flow on that density, m³/s.
    pub inlet_volumetric_flow_m3_s: f64,
    /// Differential head Δp / (ρ g), m, with standard gravity.
    pub differential_head_m: f64,
    /// Hydraulic power V̇ Δp, W.
    pub hydraulic_power_w: f64,
    /// Model and data identity of the density, for labelling.
    pub density_basis: String,
    /// Whether the density lies in its qualified domain or is unassessed.
    pub density_validity: thermo_types::ValidityStatus,
}

/// Standard acceleration of gravity, m/s² (CGPM 1901).
pub const STANDARD_GRAVITY_M_S2: f64 = 9.80665;

/// Hydraulic quantities of a pump raising `inlet` (molar flow `n`) by `dp`.
fn pump_hydraulics(inlet: &EquilibriumState, n: f64, dp: f64, ctx: &CalcContext) -> PumpHydraulics {
    if inlet.phases.is_empty() || inlet.phases.iter().any(|p| p.phase != Phase::Liquid) {
        return PumpHydraulics::Unavailable(
            "the engineering liquid density applies to liquid inlets only".into(),
        );
    }
    let mut volume = 0.0;
    let mut mass = 0.0;
    let mut validity = thermo_types::ValidityStatus::WithinDeclaredDomain;
    let mut basis = String::new();
    for p in &inlet.phases {
        let d = match ctx.package().engineering_liquid_density(
            &p.composition,
            p.state.temperature_k,
            p.state.pressure_pa,
        ) {
            Ok(d) => d,
            Err(e) => return PumpHydraulics::Unavailable(e.to_string()),
        };
        volume += p.phase_fraction * d.molar_volume;
        mass += p.phase_fraction * d.molar_volume * d.mass_density;
        if d.validity != thermo_types::ValidityStatus::WithinDeclaredDomain {
            validity = d.validity;
        }
        basis = format!(
            "{} v{} ({} {})",
            d.model.identifier, d.model.implementation_version, d.data.identifier, d.data.version
        );
    }
    let density = mass / volume;
    let flow = n * volume;
    PumpHydraulics::Calculated(PumpHydraulicValues {
        inlet_density_kg_m3: density,
        inlet_volumetric_flow_m3_s: flow,
        differential_head_m: dp / (density * STANDARD_GRAVITY_M_S2),
        hydraulic_power_w: flow * dp,
        density_basis: basis,
        density_validity: validity,
    })
}

impl RotatingResult {
    /// Total shaft power, W, positive into the fluid (negative for an
    /// expander).
    pub fn power_w(&self) -> f64 {
        self.stages.iter().map(|s| s.power_w).sum()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Machine {
    Pump,
    Compressor,
    Expander,
}

impl Machine {
    fn label(self) -> &'static str {
        match self {
            Machine::Pump => "pump",
            Machine::Compressor => "compressor",
            Machine::Expander => "expander",
        }
    }
}

fn phase_fraction(state: &EquilibriumState, phase: Phase) -> f64 {
    state
        .phases
        .iter()
        .filter(|p| p.phase == phase)
        .map(|p| p.phase_fraction)
        .sum()
}

/// Inlet phase rules: a pump takes no vapor, a compressor no liquid. A single
/// supercritical phase is accepted by both.
fn check_inlet_phase(machine: Machine, state: &EquilibriumState) -> Result<(), String> {
    match machine {
        Machine::Pump => {
            let vapor = phase_fraction(state, Phase::Vapor);
            if vapor > 0.0 {
                return Err(format!(
                    "pump inlet contains vapor (molar vapor fraction {vapor:.6}); \
                     use a compressor or subcool the feed"
                ));
            }
        }
        Machine::Compressor => {
            let liquid = phase_fraction(state, Phase::Liquid);
            if liquid > 0.0 {
                return Err(format!(
                    "compressor inlet contains liquid (molar liquid fraction {liquid:.6}); \
                     separate the liquid or use a pump"
                ));
            }
        }
        Machine::Expander => {}
    }
    Ok(())
}

/// One isentropic-efficiency step from `inlet` to `p_out` at molar flow `n`.
fn stage(
    machine: Machine,
    inlet: &EquilibriumState,
    p_out: f64,
    efficiency: f64,
    n: f64,
    ctx: &CalcContext,
) -> Result<(EquilibriumState, RotatingStage, Vec<String>), String> {
    let name = machine.label();
    check_inlet_phase(machine, inlet)?;
    let p_in = inlet.pressure_pa;
    let compress = machine != Machine::Expander;
    if !p_out.is_finite() || p_out <= 0.0 {
        return Err(format!(
            "{name} requires a finite, positive outlet pressure"
        ));
    }
    if compress && p_out < p_in {
        return Err(format!(
            "{name} outlet pressure {p_out} Pa is below the inlet pressure {p_in} Pa"
        ));
    }
    if !compress && p_out > p_in {
        return Err(format!(
            "{name} outlet pressure {p_out} Pa exceeds the inlet pressure {p_in} Pa"
        ));
    }
    let z = &inlet.mole_fractions;
    let h_in = inlet.enthalpy();
    let (ideal, out) = if p_out == p_in {
        // No pressure change: keep the inlet equilibrium phases unchanged.
        (inlet.clone(), inlet.clone())
    } else {
        let ideal = ctx
            .ps(z, p_out, inlet.entropy())
            .map_err(|e| format!("{name} isentropic outlet: {e}"))?;
        let dh_s = ideal.enthalpy() - h_in;
        if !dh_s.is_finite() || (compress && dh_s < 0.0) || (!compress && dh_s > 0.0) {
            return Err(format!(
                "{name} isentropic enthalpy change {dh_s} J/mol has the wrong sign"
            ));
        }
        let out = if efficiency == 1.0 {
            ideal.clone()
        } else {
            let dh = if compress {
                dh_s / efficiency
            } else {
                dh_s * efficiency
            };
            ctx.ph(z, p_out, h_in + dh)
                .map_err(|e| format!("{name} outlet: {e}"))?
        };
        (ideal, out)
    };
    let mut warnings = Vec::new();
    match machine {
        Machine::Pump if phase_fraction(&out, Phase::Vapor) > 0.0 => warnings
            .push("Pump outlet contains vapor; check the specification and inlet state.".into()),
        Machine::Compressor if phase_fraction(&out, Phase::Liquid) > 0.0 => {
            warnings.push("Compressor discharge contains liquid (retrograde condensation).".into())
        }
        Machine::Expander if phase_fraction(&out, Phase::Liquid) > 0.0 => {
            warnings.push(format!(
                "Expander outlet is two-phase or liquid: molar liquid fraction {:.6}.",
                phase_fraction(&out, Phase::Liquid)
            ));
        }
        _ => {}
    }
    let result = RotatingStage {
        inlet_pressure_pa: p_in,
        outlet_pressure_pa: p_out,
        isentropic_outlet_temperature_k: ideal.temperature_k,
        outlet_temperature_k: out.temperature_k,
        isentropic_power_w: n * (ideal.enthalpy() - h_in),
        power_w: n * (out.enthalpy() - h_in),
    };
    Ok((out, result, warnings))
}

fn zero_flow_output(kind: UnitKind, efficiency: f64, energy_ports: usize) -> UnitOutput {
    let mut out = UnitOutput::new(
        vec![empty()],
        0.0,
        vec![format!(
            "Zero inlet flow: the {} outlet is empty and no power is transferred.",
            kind.label().to_lowercase()
        )],
    );
    out.energy = vec![Some(0.0); energy_ports];
    out.rotating = Some(RotatingResult {
        efficiency,
        stages: Vec::new(),
        intercooler_duty_w: 0.0,
        pump_hydraulics: None,
    });
    out
}

fn inlet_state<'a>(
    kind: UnitKind,
    inputs: &'a UnitInputs,
) -> Result<Option<(&'a EquilibriumState, f64)>, String> {
    let inlet = inputs
        .material
        .first()
        .ok_or_else(|| format!("{kind} has no inlet"))?;
    let n = inlet.molar_flow_mol_s;
    if !n.is_finite() || n < 0.0 {
        return Err(format!("{kind} requires a finite, nonnegative inlet flow"));
    }
    if n == 0.0 {
        return Ok(None);
    }
    let state = inlet
        .state
        .as_ref()
        .ok_or_else(|| format!("{kind} inlet has no thermodynamic state"))?;
    Ok(Some((state, n)))
}

fn first_issue(issues: Vec<SpecIssue>) -> Result<(), String> {
    match issues.first() {
        Some(issue) => Err(format!("{} {}", issue.field, issue.problem)),
        None => Ok(()),
    }
}

/// Pump, compressor and expander ([`UnitSpec::Pressurize`] /
/// [`UnitSpec::Expand`]).
pub(crate) fn calculate_core(
    kind: UnitKind,
    spec: &UnitSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    first_issue(spec.issues())?;
    let (machine, efficiency) = match (kind, spec) {
        (UnitKind::Pump, UnitSpec::Pressurize { efficiency, .. }) => (Machine::Pump, *efficiency),
        (UnitKind::Compressor, UnitSpec::Pressurize { efficiency, .. }) => {
            (Machine::Compressor, *efficiency)
        }
        (UnitKind::Expander, UnitSpec::Expand { efficiency, .. }) => {
            (Machine::Expander, *efficiency)
        }
        _ => return Err(format!("{kind} has a specification of another kind")),
    };
    let Some((state, n)) = inlet_state(kind, inputs)? else {
        return Ok(zero_flow_output(kind, efficiency, 1));
    };
    let p_in = state.pressure_pa;
    let p_out = match *spec {
        UnitSpec::Pressurize {
            pressure: PressureRise::OutletPressure(Some(p)),
            ..
        }
        | UnitSpec::Expand {
            pressure: PressureSpec::OutletPressure(Some(p)),
            ..
        } => p,
        UnitSpec::Pressurize {
            pressure: PressureRise::PressureIncrease(Some(dp)),
            ..
        } => p_in + dp,
        UnitSpec::Expand {
            pressure: PressureSpec::PressureDrop(Some(dp)),
            ..
        } => p_in - dp,
        _ => {
            return Err(format!(
                "missing {} pressure specification",
                machine.label()
            ));
        }
    };
    let (out, result, warnings) = stage(machine, state, p_out, efficiency, n, ctx)?;
    let pump_hydraulics =
        (machine == Machine::Pump).then(|| pump_hydraulics(state, n, p_out - p_in, ctx));
    let power = result.power_w;
    let mut output = UnitOutput::new(
        vec![MaterialStream {
            molar_flow_mol_s: n,
            state: Some(out),
        }],
        power,
        warnings,
    );
    output.rotating = Some(RotatingResult {
        efficiency,
        stages: vec![result],
        intercooler_duty_w: 0.0,
        pump_hydraulics,
    });
    Ok(output)
}

/// Calculates a [`RotatingSpec`] unit.
pub(crate) fn calculate(
    kind: UnitKind,
    spec: &RotatingSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    match spec {
        RotatingSpec::MultistageCompressor(s) => multistage(kind, s, inputs, ctx),
    }
}

fn multistage(
    kind: UnitKind,
    spec: &MultistageCompressorSpec,
    inputs: &UnitInputs,
    ctx: &CalcContext,
) -> Result<UnitOutput, String> {
    first_issue(spec.issues())?;
    let Some((state, n)) = inlet_state(kind, inputs)? else {
        return Ok(zero_flow_output(kind, spec.efficiency, 2));
    };
    let p_in = state.pressure_pa;
    let p_out = spec.outlet_pressure_pa.ok_or("missing outlet pressure")?;
    if p_out < p_in {
        return Err(format!(
            "multistage compressor outlet pressure {p_out} Pa is below the inlet pressure {p_in} Pa"
        ));
    }
    let count = spec.stages;
    let ratio = (p_out / p_in).powf(1.0 / f64::from(count));
    let mut current = state.clone();
    let mut stages = Vec::with_capacity(usize::from(count));
    let mut warnings = Vec::new();
    let mut cooled_w = 0.0;
    for k in 1..=count {
        let p = if k == count {
            p_out
        } else {
            p_in * ratio.powi(i32::from(k))
        };
        let (discharge, result, stage_warnings) =
            stage(Machine::Compressor, &current, p, spec.efficiency, n, ctx)
                .map_err(|e| format!("stage {k}: {e}"))?;
        warnings.extend(
            stage_warnings
                .into_iter()
                .map(|w| format!("Stage {k}: {w}")),
        );
        stages.push(result);
        current = discharge;
        if k < count {
            let t = spec
                .intercooler_temperature_k
                .ok_or("missing intercooler temperature")?;
            if t > current.temperature_k {
                return Err(format!(
                    "intercooler {k} outlet temperature {t} K exceeds the stage discharge \
                     temperature {} K",
                    current.temperature_k
                ));
            }
            let cooled = ctx
                .tp(&current.mole_fractions, t, p)
                .map_err(|e| format!("intercooler {k}: {e}"))?;
            if phase_fraction(&cooled, Phase::Liquid) > 0.0 {
                return Err(format!(
                    "intercooler {k} condenses liquid at {t} K; raise the intercooler \
                     temperature or separate the condensate"
                ));
            }
            cooled_w += n * (current.enthalpy() - cooled.enthalpy());
            current = cooled;
        }
    }
    let result = RotatingResult {
        efficiency: spec.efficiency,
        stages,
        intercooler_duty_w: cooled_w,
        pump_hydraulics: None,
    };
    let power = result.power_w();
    let mut output = UnitOutput::new(
        vec![MaterialStream {
            molar_flow_mol_s: n,
            state: Some(current),
        }],
        power - cooled_w,
        warnings,
    );
    output.energy = vec![Some(power), Some(cooled_w)];
    output.rotating = Some(result);
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotating_kinds_are_calculated() {
        for kind in [
            UnitKind::Pump,
            UnitKind::Compressor,
            UnitKind::Expander,
            UnitKind::MultistageCompressor,
        ] {
            assert_eq!(kind.calculation_status(), CalculationStatus::Calculated);
        }
    }

    #[test]
    fn multistage_specification_issues() {
        let mut spec = MultistageCompressorSpec::default();
        let fields = |s: &MultistageCompressorSpec| {
            s.issues().into_iter().map(|i| i.field).collect::<Vec<_>>()
        };
        assert_eq!(
            fields(&spec),
            ["outlet pressure", "intercooler temperature"]
        );
        spec.stages = 1;
        assert_eq!(fields(&spec), ["outlet pressure"]);
        spec.stages = MAX_COMPRESSOR_STAGES + 1;
        spec.outlet_pressure_pa = Some(6e6);
        spec.intercooler_temperature_k = Some(310.0);
        spec.efficiency = 1.5;
        assert_eq!(fields(&spec), ["stages", "efficiency"]);
    }
}
