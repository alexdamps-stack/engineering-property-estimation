//! Bounded PR/PR-BM and GERG-2008 bubble/dew boundaries; no critical-point or
//! retrograde continuation.
//!
//! #227 extends the cubic domain beyond the dry C1–C3/N2/CO2 set to every
//! component of the cubic package (including water, methanol and
//! pseudo-components). A boundary whose active set is not dry (an *extended*
//! set) additionally gets aqueous and hydrocarbon incipient-liquid seeds when
//! water is present, multiphase (VLLE) tangent-plane evidence for its phases,
//! and, for a bubble point whose specified liquid splits into two liquids, the
//! three-phase (heterogeneous) bubble point bracketed with the VLLE flash.
use crate::equilibrium_model::EquilibriumModel;
use crate::gerg_flash::{GergEquilibrium, validate};
use crate::gerg_ph::{gerg_pure_saturation_pressure, gerg_pure_saturation_temperature};
use crate::stability::assess_stability_impl;
use crate::{
    MultiphaseOptions, MultiphaseStability, MultiphaseTpFlash, PhaseRole,
    SATURATION_MAX_REDUCED_TEMPERATURE, StabilityAssessment, StabilityOptions, assess_stability,
    gerg_flash::QUALIFIED_DRY_COMPONENTS, pure_saturation_pressure,
    saturation::single_active_component, vlle_stability, vlle_tp_flash,
};
use thermo_models::{
    CUBIC_PRESSURE_RANGE_PA, CUBIC_TEMPERATURE_RANGE_K, CubicModel, GergModel, HomogeneousState,
    RootSelection,
};
use thermo_types::{CalculationError as Error, Pressure, Temperature};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicSaturationKind {
    Bubble,
    Dew,
}

/// Bounds are part of the requested branch definition. Multiple physical roots
/// inside them are rejected rather than silently choosing a retrograde branch.
#[derive(Debug, Clone, Copy)]
pub struct SaturationOptions {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub minimum_pressure_pa: f64,
    pub maximum_pressure_pa: f64,
    pub max_iterations: usize,
}
impl Default for SaturationOptions {
    fn default() -> Self {
        Self {
            minimum_temperature_k: CUBIC_TEMPERATURE_RANGE_K[0],
            maximum_temperature_k: 350.0,
            minimum_pressure_pa: 10_000.0,
            maximum_pressure_pa: 10_000_000.0,
            max_iterations: 100,
        }
    }
}
impl SaturationOptions {
    pub fn validate(&self) -> Result<(), Error> {
        let values = [
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.minimum_pressure_pa,
            self.maximum_pressure_pa,
        ];
        let [t_lo, t_hi] = CUBIC_TEMPERATURE_RANGE_K;
        let [p_lo, p_hi] = CUBIC_PRESSURE_RANGE_PA;
        // Invalid options use the same kinds as the PH and flash validators.
        if values.iter().any(|v| !v.is_finite())
            || self.minimum_temperature_k < t_lo
            || self.maximum_temperature_k > t_hi
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || self.minimum_pressure_pa < p_lo
            || self.maximum_pressure_pa > p_hi
            || self.minimum_pressure_pa >= self.maximum_pressure_pa
            || self.max_iterations > 2000
        {
            return Err(Error::UnsupportedCapability(format!(
                "A13 requires ordered finite bounds inside {t_lo}..{t_hi} K, {p_lo}..{p_hi:e} Pa and iteration budget 1..2000"
            )));
        }
        if self.max_iterations == 0 {
            return Err(Error::NonConvergence(
                "A13 saturation iteration limit is zero".into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone)]
pub struct SaturationPoint {
    pub kind: CubicSaturationKind,
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub liquid_composition: Vec<f64>,
    pub vapor_composition: Vec<f64>,
    pub liquid: HomogeneousState,
    pub vapor: HomogeneousState,
    pub max_ln_fugacity_residual: f64,
    /// Iterations in the accepted Newton path (other seeds are not included).
    pub iterations: usize,
    /// Finite-trial TPD evidence for liquid and vapor, in that order.
    pub phase_stability: [StabilityAssessment; 2],
    /// Multiphase (VLLE) TPD evidence for an extended (non-dry) component set:
    /// liquid then vapor for a two-phase point, the vapor for a three-phase
    /// bubble point. Empty for the dry A13 set and for GERG-2008.
    pub multiphase_stability: Vec<MultiphaseStability>,
    /// The second liquid of a three-phase (heterogeneous) bubble point. The
    /// specified liquid has then split: `liquid` and `liquid_composition` are the
    /// lighter liquid (lower mass density), and this is the heavier one. `None`
    /// for an ordinary two-phase boundary.
    pub second_liquid: Option<SecondLiquid>,
}

/// Heavier liquid of a three-phase bubble point, in equilibrium with the lighter
/// liquid and the incipient vapor.
#[derive(Debug, Clone)]
pub struct SecondLiquid {
    pub composition: Vec<f64>,
    pub state: HomogeneousState,
    /// Moles of this liquid per mole of the specified (overall) liquid.
    pub material_fraction: f64,
}

/// Temperatures of the solid–liquid–vapor triple points of components that
/// freeze inside the cubic temperature interval, K. Solids are not modelled: a
/// boundary below one of these with a liquid that is at least
/// [`FREEZING_LIQUID_FRACTION`] of that component is outside the domain.
pub const TRIPLE_POINT_TEMPERATURES_K: [(&str, f64); 2] =
    [("water", 273.16), ("carbon-dioxide", 216.592)];
/// Liquid mole fraction of a freezing component at or above which a boundary
/// below its triple-point temperature is rejected.
pub const FREEZING_LIQUID_FRACTION: f64 = 0.9;
/// Minimum max_i |x_i^L1 − x_i^L2| of a three-phase bubble point; closer liquids
/// are treated as a plait (liquid–liquid critical) point.
pub const PLAIT_LIQUID_SEPARATION: f64 = 0.03;

/// Composition is the specified liquid for bubble, specified vapor for dew.
pub fn saturation_pressure(
    model: &CubicModel,
    temperature: Temperature,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    solve(model, temperature.value(), true, kind, options).map_err(|r| r.error)
}
pub fn saturation_temperature(
    model: &CubicModel,
    pressure: Pressure,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    solve(model, pressure.value(), false, kind, options).map_err(|r| r.error)
}

/// Bubble/dew pressure of a native GERG-2008 fluid with the bounded A13 algorithm
/// (the same seeds, fugacity re-check, near-critical, metastability and TPD
/// rejections as the cubic solver), restricted to the GERG equilibrium domain
/// (qualified dry components, 230..=350 K, 0.1..=10 MPa).
pub fn gerg_saturation_pressure(
    model: &GergModel,
    temperature: Temperature,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    let adapter = gerg_saturation_model(model, options)?;
    solve(&adapter, temperature.value(), true, kind, options).map_err(|r| r.error)
}

/// Bubble/dew temperature of a native GERG-2008 fluid; see [`gerg_saturation_pressure`].
pub fn gerg_saturation_temperature(
    model: &GergModel,
    pressure: Pressure,
    kind: CubicSaturationKind,
    options: SaturationOptions,
) -> Result<SaturationPoint, Error> {
    let adapter = gerg_saturation_model(model, options)?;
    solve(&adapter, pressure.value(), false, kind, options).map_err(|r| r.error)
}

fn gerg_saturation_model(
    model: &GergModel,
    options: SaturationOptions,
) -> Result<GergEquilibrium, Error> {
    options.validate()?;
    if options.minimum_temperature_k < 230.0
        || options.maximum_temperature_k > 350.0
        || options.minimum_pressure_pa < 1e5
        || options.maximum_pressure_pa > 1e7
    {
        return Err(Error::OutsideValidityDomain(
            "GERG saturation bounds must lie in the equilibrium envelope 230..=350 K, 0.1..=10 MPa"
                .into(),
        ));
    }
    // Rejects components outside the qualified dry set before any solve.
    validate(
        model,
        Temperature::new(options.minimum_temperature_k)?,
        Pressure::new(1e5)?,
    )?;
    GergEquilibrium::new(model)
}

/// Model operations of the bounded saturation algorithm beyond the shared
/// equilibrium contract: pure-fluid saturation and finite-trial stability.
trait SaturationModel: EquilibriumModel {
    fn component_ids(&self) -> Vec<&str>;
    fn critical_temperature(&self, index: usize) -> f64;
    /// Saturation pressure and iterations of the canonical pure model at T.
    fn pure_saturation_pressure(
        &self,
        t: f64,
        max_iterations: usize,
    ) -> Result<(f64, usize), Error>;
    /// Saturation temperature and iterations of the canonical pure model at p,
    /// searched in `bracket` (already limited to 0.95 Tc).
    fn pure_saturation_temperature(
        &self,
        p: f64,
        bracket: [f64; 2],
        max_iterations: usize,
    ) -> Result<(f64, usize), Rejection>;
    fn assess_stability(&self, t: Temperature, p: Pressure) -> Result<StabilityAssessment, Error>;
    /// True when an active component lies outside the dry A13 set, so that the
    /// extended (multiphase-aware) treatment applies. GERG-2008 is never extended.
    fn extended(&self) -> bool {
        false
    }
    /// Multiphase TPD evidence of this composition; only for extended models.
    fn multiphase_stability(
        &self,
        _t: Temperature,
        _p: Pressure,
    ) -> Result<MultiphaseStability, Error> {
        Err(Error::UnsupportedCapability(
            "multiphase stability unavailable for this saturation model".into(),
        ))
    }
    /// Multiphase TP flash of this composition; only for extended models.
    fn multiphase_flash(
        &self,
        _t: Temperature,
        _p: Pressure,
    ) -> Result<MultiphaseTpFlash<HomogeneousState>, Error> {
        Err(Error::UnsupportedCapability(
            "multiphase flash unavailable for this saturation model".into(),
        ))
    }
}

impl SaturationModel for CubicModel {
    fn component_ids(&self) -> Vec<&str> {
        self.component_ids()
    }
    fn extended(&self) -> bool {
        self.component_ids()
            .iter()
            .zip(self.fractions())
            .any(|(id, x)| *x > 0.0 && !QUALIFIED_DRY_COMPONENTS.contains(id))
    }
    fn multiphase_stability(
        &self,
        t: Temperature,
        p: Pressure,
    ) -> Result<MultiphaseStability, Error> {
        vlle_stability(self, t, p, MultiphaseOptions::default())
    }
    fn multiphase_flash(
        &self,
        t: Temperature,
        p: Pressure,
    ) -> Result<MultiphaseTpFlash<HomogeneousState>, Error> {
        vlle_tp_flash(self, t, p, MultiphaseOptions::default())
    }
    fn critical_temperature(&self, index: usize) -> f64 {
        self.components()[index].critical_temperature()
    }
    fn pure_saturation_pressure(
        &self,
        t: f64,
        max_iterations: usize,
    ) -> Result<(f64, usize), Error> {
        let r = pure_saturation_pressure(self, temperature(t)?, max_iterations)?;
        Ok((r.pressure_pa, r.iterations))
    }
    fn pure_saturation_temperature(
        &self,
        fixed: f64,
        [mut a, mut b]: [f64; 2],
        max_iterations: usize,
    ) -> Result<(f64, usize), Rejection> {
        let p_a = pure_saturation_pressure(self, temperature(a)?, 200)?.pressure_pa;
        let p_b = pure_saturation_pressure(self, temperature(b)?, 200)?.pressure_pa;
        if fixed < p_a || fixed > p_b {
            return Err(Rejection::near_critical(
                "A13 pure saturation not bracketed; near-critical or outside requested range",
            ));
        }
        let mut count = 0;
        for i in 1..=max_iterations {
            count = i;
            let mid = (a + b) / 2.0;
            let r = pure_saturation_pressure(self, temperature(mid)?, 200)?;
            if r.pressure_pa < fixed {
                a = mid;
            } else {
                b = mid;
            }
            if (b - a) < 1e-9 {
                break;
            }
        }
        if b - a >= 1e-9 {
            return Err(
                Error::NonConvergence("A13 pure temperature iteration limit".into()).into(),
            );
        }
        Ok(((a + b) / 2.0, count))
    }
    fn assess_stability(&self, t: Temperature, p: Pressure) -> Result<StabilityAssessment, Error> {
        assess_stability(self, t, p, StabilityOptions::default())
    }
}

impl SaturationModel for GergEquilibrium {
    fn component_ids(&self) -> Vec<&str> {
        self.model.component_ids()
    }
    fn critical_temperature(&self, index: usize) -> f64 {
        self.critical_temperature(index)
    }
    fn pure_saturation_pressure(
        &self,
        t: f64,
        max_iterations: usize,
    ) -> Result<(f64, usize), Error> {
        gerg_pure_saturation_pressure(&self.model, t, max_iterations)
    }
    fn pure_saturation_temperature(
        &self,
        p: f64,
        [a, b]: [f64; 2],
        max_iterations: usize,
    ) -> Result<(f64, usize), Rejection> {
        match gerg_pure_saturation_temperature(&self.model, pressure(p)?, a, b, max_iterations)? {
            Some(sat) => Ok((sat.temperature_k, sat.iterations)),
            None => Err(Rejection::near_critical(
                "A13 pure saturation not bracketed; near-critical or outside requested range",
            )),
        }
    }
    fn assess_stability(&self, t: Temperature, p: Pressure) -> Result<StabilityAssessment, Error> {
        let options = StabilityOptions::default();
        options.validate()?;
        assess_stability_impl(self, t, p, options)
    }
}

/// Why a boundary was not returned, ordered from least to most specific. When no
/// verified root exists, the most specific rejection over all seeds is reported.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum RejectionKind {
    /// Input, bound or model-evaluation error, or an inner solver failure.
    Other,
    /// Newton converged but the 1e-9 fugacity re-check failed.
    FugacityRecheck,
    /// Trivial/coalescing roots, singular boundary equations or T/Tc above the
    /// supported pure-saturation limit.
    NearCritical,
    /// A selected phase root is metastable relative to another root.
    MetastablePhase,
    /// Boundary phases unstable in finite TPD trials (third phase or LLE).
    UnstableBoundaryPhases,
    /// Below a triple point with an essentially pure freezing liquid.
    SolidPhase,
}

#[derive(Debug, Clone)]
struct Rejection {
    kind: RejectionKind,
    error: Error,
}

impl Rejection {
    fn new(kind: RejectionKind, error: Error) -> Self {
        Self { kind, error }
    }
    fn near_critical(message: &str) -> Self {
        Self::new(
            RejectionKind::NearCritical,
            Error::OutsideValidityDomain(message.into()),
        )
    }
}

impl From<Error> for Rejection {
    fn from(error: Error) -> Self {
        Self::new(RejectionKind::Other, error)
    }
}

/// Keep the most specific rejection; the first seen wins among equals.
fn keep_most_specific(best: &mut Option<Rejection>, candidate: Rejection) {
    if best.as_ref().is_none_or(|b| candidate.kind > b.kind) {
        *best = Some(candidate);
    }
}

/// max_i |ln K_i| below this marks trivial/coalescing phases. Unlike an absolute
/// |x_i − y_i| test, it does not reject dilute mixtures whose phase compositions
/// differ little only because the feed is dilute.
const NEAR_CRITICAL_LN_K: f64 = 0.03;

fn solve<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    kind: CubicSaturationKind,
    o: SaturationOptions,
) -> Result<SaturationPoint, Rejection> {
    o.validate()?;
    let (lo, hi) = if pressure_unknown {
        (o.minimum_pressure_pa.ln(), o.maximum_pressure_pa.ln())
    } else {
        (o.minimum_temperature_k.ln(), o.maximum_temperature_k.ln())
    };
    let (flo, fhi) = if pressure_unknown {
        (o.minimum_temperature_k, o.maximum_temperature_k)
    } else {
        (o.minimum_pressure_pa, o.maximum_pressure_pa)
    };
    if !fixed.is_finite() || fixed < flo || fixed > fhi {
        return Err(Error::OutsideValidityDomain(
            "A13 fixed coordinate outside requested supported range".into(),
        )
        .into());
    }
    if let Some(index) = single_active_component(model.fractions()) {
        // with_fractions accepts a small unit-sum error without normalizing.
        // Canonicalize the single active component for both pure saturation
        // paths and their returned phase states, retaining inactive slots.
        let mut fractions = vec![0.0; model.fractions().len()];
        fractions[index] = 1.0;
        let pure_model = model.with_fractions(&fractions)?;
        let model = &pure_model;
        let tc = model.critical_temperature(index);
        let (t, p, it) = if pressure_unknown {
            if fixed / tc > SATURATION_MAX_REDUCED_TEMPERATURE {
                return Err(Rejection::near_critical(&format!(
                    "A13 near-critical diagnostic: pure saturation supported for T/Tc <= {SATURATION_MAX_REDUCED_TEMPERATURE}"
                )));
            }
            let (p, iterations) = model.pure_saturation_pressure(fixed, o.max_iterations)?;
            (fixed, p, iterations)
        } else {
            let a = o.minimum_temperature_k;
            let mut b = o
                .maximum_temperature_k
                .min(SATURATION_MAX_REDUCED_TEMPERATURE * tc);
            // 0.95 Tc can round to a temperature whose T/Tc exceeds 0.95, which the
            // pure-pressure solver rejects; step down to the largest admitted value.
            while b / tc > SATURATION_MAX_REDUCED_TEMPERATURE {
                b = b.next_down();
            }
            if b <= a {
                return Err(Rejection::near_critical(
                    "A13 near-critical diagnostic: no supported pure-temperature interval",
                ));
            }
            let (t, iterations) =
                model.pure_saturation_temperature(fixed, [a, b], o.max_iterations)?;
            (t, fixed, iterations)
        };
        if p < o.minimum_pressure_pa || p > o.maximum_pressure_pa {
            return Err(Error::OutsideValidityDomain(
                "A13 saturation pressure outside requested range".into(),
            )
            .into());
        }
        return point(
            model,
            t,
            p,
            kind,
            model.fractions().to_vec(),
            model.fractions().to_vec(),
            it,
        );
    }
    let active: Vec<_> = model
        .fractions()
        .iter()
        .enumerate()
        .filter_map(|(i, x)| (*x > 0.0).then_some(i))
        .collect();
    let extended = model.extended();
    let mut solutions: Vec<SaturationPoint> = Vec::new();
    let mut rejection: Option<Rejection> = None;
    for seed in 0..=12 {
        let u = lo + (hi - lo) * seed as f64 / 12.0;
        let (t, p) = if pressure_unknown {
            (fixed, u.exp())
        } else {
            (u.exp(), fixed)
        };
        let k = model.wilson_values(temperature(t)?, pressure(p)?);
        let mut starts = vec![active.iter().map(|i| k[*i].ln()).collect::<Vec<_>>()];
        if extended && kind == CubicSaturationKind::Dew {
            starts.extend(incipient_liquid_seeds(model, &active));
        }
        for ln_k in starts {
            let mut q = vec![u];
            q.extend(ln_k);
            newton_from_seed(
                model,
                fixed,
                pressure_unknown,
                kind,
                &active,
                q,
                (lo, hi),
                o.max_iterations,
                &mut solutions,
                &mut rejection,
            );
        }
    }
    match solutions.len() {
        1 => Ok(solutions.remove(0)),
        0 if extended && kind == CubicSaturationKind::Bubble => {
            // A liquid that splits into two liquids has no two-phase bubble point;
            // try the three-phase (heterogeneous) bubble point instead.
            match three_phase_bubble(model, fixed, pressure_unknown, &active, (lo, hi), o) {
                Ok(point) => Ok(point),
                Err(ThreePhase::NotBracketed) => Err(no_boundary(rejection)),
                // A specific three-phase rejection (plait point, fugacity re-check,
                // metastable or unstable phases) is reported; a generic flash or
                // iteration failure does not mask a two-phase rejection.
                Err(ThreePhase::Rejected(r)) => match rejection {
                    Some(two_phase) if r.kind == RejectionKind::Other => Err(two_phase),
                    _ => Err(r),
                },
            }
        }
        // Report the most specific reason a converged root was rejected.
        0 => Err(no_boundary(rejection)),
        _ => Err(Error::UnsupportedCapability("A13 multiple saturation branches in requested interval; narrow bounds; retrograde branch selection unsupported".into()).into()),
    }
}

/// Outcome of the three-phase bubble attempt that did not return a point.
enum ThreePhase {
    /// The VLLE scan found no two-liquid state next to a vapor-bearing one: there
    /// is no heterogeneous bubble point to report, and the two-phase rejection stands.
    NotBracketed,
    Rejected(Rejection),
}

impl From<Error> for ThreePhase {
    fn from(error: Error) -> Self {
        Self::Rejected(error.into())
    }
}

/// Phase set of the specified liquid at one coordinate of the VLLE scan.
enum ScanState {
    TwoLiquids(MultiphaseTpFlash<HomogeneousState>),
    Vapor(MultiphaseTpFlash<HomogeneousState>),
    Other,
}

/// Number of scan intervals of the three-phase bracket search.
const THREE_PHASE_SCAN_INTERVALS: usize = 24;
/// Bisection stops once the bracket is this narrow in ln T or ln p.
const THREE_PHASE_BRACKET_WIDTH: f64 = 1e-6;

fn coordinates(fixed: f64, pressure_unknown: bool, u: f64) -> (f64, f64) {
    if pressure_unknown {
        (fixed, u.exp())
    } else {
        (u.exp(), fixed)
    }
}

fn scan_state<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    u: f64,
) -> Result<ScanState, Error> {
    let (t, p) = coordinates(fixed, pressure_unknown, u);
    let flash = model.multiphase_flash(temperature(t)?, pressure(p)?)?;
    let vapor = flash.phases.iter().any(|ph| ph.role == PhaseRole::Vapor);
    Ok(if vapor {
        ScanState::Vapor(flash)
    } else if flash.phases.len() == 2 {
        ScanState::TwoLiquids(flash)
    } else {
        ScanState::Other
    })
}

/// Three-phase (heterogeneous) bubble point of a specified liquid that splits
/// into two liquids. The VLLE flash scans from the liquid side (high pressure or
/// low temperature) for the first vapor-bearing state after a two-liquid one,
/// bisects that bracket, and Newton's method then solves
/// ln K2_i + ln φ_i^L2 − ln φ_i^L1 = 0, ln Kv_i + ln φ_i^V − ln φ_i^L1 = 0,
/// Σ (x1_i − x2_i) = 0 and ln Σ y_i = 0 for ln T or ln p, ln K2, the heavy-liquid
/// amount β and ln Kv, with x1 = z/(1 + β(K2 − 1)), x2 = K2 x1, y = Kv x1.
fn three_phase_bubble<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    active: &[usize],
    (lo, hi): (f64, f64),
    o: SaturationOptions,
) -> Result<SaturationPoint, ThreePhase> {
    // exp(ln bound) can round just outside the bound (and the EOS envelope).
    let (lo, hi) = (lo + 1e-12, hi - 1e-12);
    let at = |k: usize| {
        let s = k as f64 / THREE_PHASE_SCAN_INTERVALS as f64;
        if pressure_unknown {
            hi - (hi - lo) * s
        } else {
            lo + (hi - lo) * s
        }
    };
    let mut previous: Option<(f64, MultiphaseTpFlash<HomogeneousState>)> = None;
    let mut bracket = None;
    for k in 0..=THREE_PHASE_SCAN_INTERVALS {
        let u = at(k);
        match scan_state(model, fixed, pressure_unknown, u)? {
            ScanState::TwoLiquids(flash) => previous = Some((u, flash)),
            ScanState::Vapor(flash) => {
                if let Some((ua, liquids)) = previous.take() {
                    bracket = Some((ua, liquids, u, flash));
                }
                break;
            }
            ScanState::Other => previous = None,
        }
    }
    let Some((mut ua, mut liquids, mut ub, mut vapor)) = bracket else {
        return Err(ThreePhase::NotBracketed);
    };
    let mut bisections = 0;
    while (ub - ua).abs() > THREE_PHASE_BRACKET_WIDTH {
        if bisections == o.max_iterations {
            return Err(ThreePhase::Rejected(
                Error::NonConvergence("A13 three-phase bubble bracket iteration limit".into())
                    .into(),
            ));
        }
        bisections += 1;
        let mid = 0.5 * (ua + ub);
        match scan_state(model, fixed, pressure_unknown, mid)? {
            ScanState::TwoLiquids(flash) => (ua, liquids) = (mid, flash),
            ScanState::Vapor(flash) => (ub, vapor) = (mid, flash),
            ScanState::Other => {
                return Err(ThreePhase::Rejected(
                    Error::NonConvergence(
                        "A13 three-phase bubble: the liquid-liquid split vanishes inside the bracket; no verified point returned".into(),
                    )
                    .into(),
                ));
            }
        }
    }
    let z = model.fractions();
    let (light, heavy) = (&liquids.phases[0], &liquids.phases[1]);
    let Some(incipient) = vapor.phases.iter().find(|ph| ph.role == PhaseRole::Vapor) else {
        return Err(ThreePhase::NotBracketed);
    };
    let m = active.len();
    let mut q = vec![0.5 * (ua + ub)];
    q.extend(
        active
            .iter()
            .map(|i| (heavy.composition[*i] / light.composition[*i]).ln()),
    );
    q.push(heavy.phase_fraction);
    q.extend(
        active
            .iter()
            .map(|i| (incipient.composition[*i] / light.composition[*i]).ln()),
    );
    if q.iter().any(|v| !v.is_finite()) || z.len() != light.composition.len() {
        return Err(ThreePhase::Rejected(
            Error::NumericalFailure("A13 three-phase bubble: invalid initial phases".into()).into(),
        ));
    }
    let mut rejection: Option<Rejection> = None;
    for iteration in 1..=o.max_iterations {
        let Ok((r, phases)) = residual3(model, fixed, pressure_unknown, active, &q) else {
            break;
        };
        let norm = r.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        if norm < 2e-11 {
            let (t, p) = coordinates(fixed, pressure_unknown, q[0]);
            let [x1, x2, y] = phases;
            return three_phase_point(model, t, p, x1, x2, y, q[1 + m], iteration)
                .map_err(ThreePhase::Rejected);
        }
        let n = q.len();
        let mut jac = vec![vec![0.0; n]; n];
        for j in 0..n {
            let mut shifted = q.clone();
            let step = 1e-7_f64.max(1e-7 * q[j].abs());
            // Difference inward at the coordinate bound and the β = 1 limit.
            let inward = (j == 0 && q[0] + step > hi) || (j == 1 + m && q[j] + step >= 1.0);
            shifted[j] += if inward { -step } else { step };
            let h = shifted[j] - q[j];
            let Ok((v, _)) = residual3(model, fixed, pressure_unknown, active, &shifted) else {
                return Err(ThreePhase::Rejected(
                    Error::NonConvergence(
                        "A13 three-phase bubble: Jacobian evaluation failed".into(),
                    )
                    .into(),
                ));
            };
            for i in 0..n {
                jac[i][j] = (v[i] - r[i]) / h;
            }
        }
        let Some(step) = linear(jac, r.iter().map(|v| -v).collect()) else {
            rejection = Some(Rejection::near_critical(
                "A13 near-critical diagnostic: singular three-phase bubble equations (plait or critical point); no verified point returned",
            ));
            break;
        };
        let mut accepted = false;
        let mut scale = 1.0;
        for _ in 0..30 {
            let next: Vec<_> = q.iter().zip(&step).map(|(a, b)| a + scale * b).collect();
            let beta = next[1 + m];
            if next[0] >= lo
                && next[0] <= hi
                && beta > 0.0
                && beta < 1.0
                && next.iter().all(|v| v.is_finite() && v.abs() < 100.0)
                && let Ok((v, _)) = residual3(model, fixed, pressure_unknown, active, &next)
                && v.iter().fold(0.0_f64, |a, v| a.max(v.abs())) < norm
            {
                q = next;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    Err(ThreePhase::Rejected(rejection.unwrap_or_else(|| {
        Error::NonConvergence(
            "A13 three-phase bubble point did not converge inside the requested bounds".into(),
        )
        .into()
    })))
}

type ThreePhaseResidual = (Vec<f64>, [Vec<f64>; 3]);

fn residual3<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    active: &[usize],
    q: &[f64],
) -> Result<ThreePhaseResidual, Error> {
    let m = active.len();
    let (t, p) = coordinates(fixed, pressure_unknown, q[0]);
    let beta = q[1 + m];
    let z = model.fractions();
    let mut x1 = vec![0.0; z.len()];
    let mut x2 = x1.clone();
    let mut y = x1.clone();
    for (j, &i) in active.iter().enumerate() {
        let k2 = q[1 + j].exp();
        let kv = q[2 + m + j].exp();
        let denominator = 1.0 + beta * (k2 - 1.0);
        if denominator.is_nan() || denominator <= 0.0 {
            return Err(Error::NumericalFailure(
                "A13 three-phase bubble: nonpositive liquid amount".into(),
            ));
        }
        x1[i] = z[i] / denominator;
        x2[i] = k2 * x1[i];
        y[i] = kv * x1[i];
    }
    let sums = [&x1, &x2, &y].map(|x| x.iter().sum::<f64>());
    if sums.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(Error::NumericalFailure(
            "A13 three-phase bubble: invalid phase composition".into(),
        ));
    }
    let normalized = [&x1, &x2, &y]
        .iter()
        .zip(sums)
        .map(|(x, s)| x.iter().map(|v| v / s).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let t = Temperature::new(t)?;
    let p = Pressure::new(p)?;
    let l1 =
        model
            .with_fractions(&normalized[0])?
            .ln_fugacity(t, p, RootSelection::SmallestVolume)?;
    let l2 =
        model
            .with_fractions(&normalized[1])?
            .ln_fugacity(t, p, RootSelection::SmallestVolume)?;
    let v =
        model
            .with_fractions(&normalized[2])?
            .ln_fugacity(t, p, RootSelection::LargestVolume)?;
    let mut r = vec![sums[2].ln(), sums[0] - sums[1]];
    r.extend(active.iter().enumerate().map(|(j, i)| {
        q[1 + j] + l2.ln_fugacity_coefficients[*i] - l1.ln_fugacity_coefficients[*i]
    }));
    r.extend(active.iter().enumerate().map(|(j, i)| {
        q[2 + m + j] + v.ln_fugacity_coefficients[*i] - l1.ln_fugacity_coefficients[*i]
    }));
    let [a, b, c]: [Vec<f64>; 3] = normalized
        .try_into()
        .map_err(|_| Error::NumericalFailure("A13 three-phase bubble phases".into()))?;
    Ok((r, [a, b, c]))
}

/// Verifies a converged three-phase bubble point with the checks of [`point`]
/// applied to all three phases, plus a plait-point guard on the two liquids.
#[allow(clippy::too_many_arguments)]
fn three_phase_point<M: SaturationModel>(
    model: &M,
    t: f64,
    p: f64,
    x1: Vec<f64>,
    x2: Vec<f64>,
    y: Vec<f64>,
    beta: f64,
    iterations: usize,
) -> Result<SaturationPoint, Rejection> {
    check_freezing(model, t, &[&x1, &x2])?;
    let tm = temperature(t)?;
    let pr = pressure(p)?;
    let models = [
        model.with_fractions(&x1)?,
        model.with_fractions(&x2)?,
        model.with_fractions(&y)?,
    ];
    let roots = [
        RootSelection::SmallestVolume,
        RootSelection::SmallestVolume,
        RootSelection::LargestVolume,
    ];
    let mut states = Vec::new();
    for (phase, root) in models.iter().zip(roots) {
        states.push(phase.evaluate(tm, pr, root)?);
    }
    let separation = x1
        .iter()
        .zip(&x2)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    if separation < PLAIT_LIQUID_SEPARATION {
        return Err(Rejection::near_critical(&format!(
            "A13 near-critical diagnostic: the two liquids of the three-phase bubble point differ by less than {PLAIT_LIQUID_SEPARATION} in mole fraction (plait point)"
        )));
    }
    if states[2].molar_volume / states[0].molar_volume < 1.5
        || states[2].molar_volume / states[1].molar_volume < 1.5
    {
        return Err(Rejection::near_critical(
            "A13 near-critical diagnostic: vapor/liquid volume ratio <1.5 at the three-phase bubble point",
        ));
    }
    if beta.is_nan() || beta <= 1e-9 || beta >= 1.0 - 1e-9 {
        return Err(Error::NonConvergence(
            "A13 three-phase bubble: a liquid amount vanishes; no verified point returned".into(),
        )
        .into());
    }
    let compositions = [&x1, &x2, &y];
    let mut residual: f64 = 0.0;
    for i in 0..x1.len() {
        if x1[i] > 0.0 && x2[i] > 0.0 && y[i] > 0.0 {
            let reference = x1[i].ln() + states[0].ln_fugacity_coefficients[i];
            for k in 1..3 {
                let value = compositions[k][i].ln() + states[k].ln_fugacity_coefficients[i];
                residual = residual.max((value - reference).abs());
            }
        }
    }
    if residual > 1e-9 {
        return Err(Rejection::new(
            RejectionKind::FugacityRecheck,
            Error::NonConvergence("A13 three-phase fugacity equality failed".into()),
        ));
    }
    for ((phase, composition), selected) in models.iter().zip(compositions).zip(&states) {
        let stable = phase.ln_fugacity(tm, pr, RootSelection::MinimumGibbs)?;
        let excess: f64 = composition
            .iter()
            .zip(&selected.ln_fugacity_coefficients)
            .zip(&stable.ln_fugacity_coefficients)
            .map(|((w, a), b)| w * (a - b))
            .sum();
        if excess > 1e-8 {
            return Err(Rejection::new(
                RejectionKind::MetastablePhase,
                Error::UnsupportedCapability(
                    "A13 selected three-phase root is metastable relative to another root".into(),
                ),
            ));
        }
    }
    let phase_stability = [
        models[0].assess_stability(tm, pr)?,
        models[2].assess_stability(tm, pr)?,
    ];
    let multiphase = models[2].multiphase_stability(tm, pr)?;
    if phase_stability.iter().any(|s| !s.stable_in_trials) || !multiphase.stable_in_trials {
        return Err(Rejection::new(
            RejectionKind::UnstableBoundaryPhases,
            Error::UnsupportedCapability(
                "A13 three-phase bubble phases unstable in TPD trials (a fourth phase or a lower tangent plane)".into(),
            ),
        ));
    }
    let [l1, l2, v]: [HomogeneousState; 3] = states
        .try_into()
        .map_err(|_| Error::NumericalFailure("A13 three-phase states".into()))?;
    // The lighter liquid is reported as `liquid`, the heavier as the second liquid.
    let (light, light_x, heavy, heavy_x, heavy_fraction) = if l1.mass_density <= l2.mass_density {
        (l1, x1, l2, x2, beta)
    } else {
        (l2, x2, l1, x1, 1.0 - beta)
    };
    Ok(SaturationPoint {
        kind: CubicSaturationKind::Bubble,
        temperature_k: t,
        pressure_pa: p,
        liquid_composition: light_x,
        vapor_composition: y,
        liquid: light,
        vapor: v,
        max_ln_fugacity_residual: residual,
        iterations,
        phase_stability,
        multiphase_stability: vec![multiphase],
        second_liquid: Some(SecondLiquid {
            composition: heavy_x,
            state: heavy,
            material_fraction: heavy_fraction,
        }),
    })
}

fn no_boundary(rejection: Option<Rejection>) -> Rejection {
    rejection.unwrap_or_else(|| {
        Error::NonConvergence("A13 no verified boundary inside requested bounds; iteration limit, absent bracket or unsupported branch".into()).into()
    })
}

/// Incipient-liquid starts for a dew point of a water-bearing vapor, in ln K of
/// the active components: an aqueous liquid (0.999 water) and a hydrocarbon
/// liquid (the water-free vapor with 0.001 water). These mirror the aqueous and
/// hydrocarbon trials of the multiphase stability test.
fn incipient_liquid_seeds<M: SaturationModel>(model: &M, active: &[usize]) -> Vec<Vec<f64>> {
    let ids = model.component_ids();
    let y = model.fractions();
    let Some(water) = active.iter().copied().find(|i| ids[*i] == "water") else {
        return Vec::new();
    };
    let dry: f64 = active.iter().filter(|i| **i != water).map(|i| y[*i]).sum();
    if dry <= 0.0 {
        return Vec::new();
    }
    [0.999, 0.001]
        .iter()
        .map(|water_fraction| {
            active
                .iter()
                .map(|i| {
                    let x = if *i == water {
                        *water_fraction
                    } else {
                        (1.0 - water_fraction) * y[*i] / dry
                    };
                    (y[*i] / x).ln()
                })
                .collect()
        })
        .collect()
}

/// Bounded finite-difference Newton from one seed; a converged root is verified
/// by [`point`] and either added to `solutions` (deduplicated) or recorded as a
/// rejection.
#[allow(clippy::too_many_arguments)]
fn newton_from_seed<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    kind: CubicSaturationKind,
    active: &[usize],
    mut q: Vec<f64>,
    (lo, hi): (f64, f64),
    max_iterations: usize,
    solutions: &mut Vec<SaturationPoint>,
    rejection: &mut Option<Rejection>,
) {
    let active = active.to_vec();
    for iteration in 1..=max_iterations {
        let Ok((r, x, y)) = residual(model, fixed, pressure_unknown, kind, &active, &q) else {
            break;
        };
        let norm = r.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        if norm < 2e-11 {
            let (t, p) = if pressure_unknown {
                (fixed, q[0].exp())
            } else {
                (q[0].exp(), fixed)
            };
            match point(model, t, p, kind, x, y, iteration) {
                Ok(s) => {
                    if !solutions.iter().any(|old| {
                        (old.temperature_k - t).abs() < 1e-5
                            && (old.pressure_pa / p - 1.0).abs() < 1e-7
                    }) {
                        solutions.push(s);
                    }
                }
                Err(r) => keep_most_specific(rejection, r),
            }
            break;
        }
        let n = q.len();
        let mut jac = vec![vec![0.0; n]; n];
        let mut valid = true;
        for j in 0..n {
            let mut shifted = q.clone();
            shifted[j] += 1e-5;
            if let Ok((v, _, _)) = residual(model, fixed, pressure_unknown, kind, &active, &shifted)
            {
                for i in 0..n {
                    jac[i][j] = (v[i] - r[i]) / 1e-5;
                }
            } else {
                valid = false;
                break;
            }
        }
        if !valid {
            break;
        }
        let Some(step) = linear(jac, r.iter().map(|v| -v).collect()) else {
            keep_most_specific(
                rejection,
                Rejection::near_critical(
                    "A13 near-critical diagnostic: singular boundary equations; no verified point returned",
                ),
            );
            break;
        };
        let mut accepted = false;
        let mut scale = 1.0;
        for _ in 0..20 {
            let next: Vec<_> = q.iter().zip(&step).map(|(a, b)| a + scale * b).collect();
            if next[0] >= lo
                && next[0] <= hi
                && next.iter().all(|v| v.is_finite() && v.abs() < 100.0)
                && let Ok((v, _, _)) =
                    residual(model, fixed, pressure_unknown, kind, &active, &next)
                && v.iter().fold(0.0_f64, |a, v| a.max(v.abs())) < norm
            {
                q = next;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            break;
        }
    }
}

fn temperature(value: f64) -> Result<Temperature, Error> {
    Ok(Temperature::new(value)?)
}

fn pressure(value: f64) -> Result<Pressure, Error> {
    Ok(Pressure::new(value)?)
}
type ResidualEvaluation = (Vec<f64>, Vec<f64>, Vec<f64>);

fn residual<M: SaturationModel>(
    model: &M,
    fixed: f64,
    pressure_unknown: bool,
    kind: CubicSaturationKind,
    active: &[usize],
    q: &[f64],
) -> Result<ResidualEvaluation, Error> {
    let (t, p) = if pressure_unknown {
        (fixed, q[0].exp())
    } else {
        (q[0].exp(), fixed)
    };
    let z = model.fractions();
    let mut x = z.to_vec();
    let mut y = z.to_vec();
    for (j, i) in active.iter().enumerate() {
        match kind {
            CubicSaturationKind::Bubble => y[*i] = z[*i] * q[j + 1].exp(),
            CubicSaturationKind::Dew => x[*i] = z[*i] * (-q[j + 1]).exp(),
        }
    }
    let sum = match kind {
        CubicSaturationKind::Bubble => y.iter().sum::<f64>(),
        CubicSaturationKind::Dew => x.iter().sum::<f64>(),
    };
    if !sum.is_finite() || sum <= 0.0 {
        return Err(Error::NumericalFailure(
            "A13 invalid incipient composition".into(),
        ));
    }
    match kind {
        CubicSaturationKind::Bubble => y.iter_mut().for_each(|v| *v /= sum),
        CubicSaturationKind::Dew => x.iter_mut().for_each(|v| *v /= sum),
    };
    let t = Temperature::new(t)?;
    let p = Pressure::new(p)?;
    let l = model
        .with_fractions(&x)?
        .ln_fugacity(t, p, RootSelection::SmallestVolume)?;
    let v = model
        .with_fractions(&y)?
        .ln_fugacity(t, p, RootSelection::LargestVolume)?;
    let mut r = vec![sum.ln()];
    r.extend(
        active.iter().enumerate().map(|(j, i)| {
            q[j + 1] + v.ln_fugacity_coefficients[*i] - l.ln_fugacity_coefficients[*i]
        }),
    );
    Ok((r, x, y))
}
fn point<M: SaturationModel>(
    model: &M,
    t: f64,
    p: f64,
    kind: CubicSaturationKind,
    x: Vec<f64>,
    y: Vec<f64>,
    iterations: usize,
) -> Result<SaturationPoint, Rejection> {
    let tm = temperature(t)?;
    let pr = pressure(p)?;
    let lm = model.with_fractions(&x)?;
    let vm = model.with_fractions(&y)?;
    check_freezing(model, t, &[&x])?;
    let l = lm.evaluate(tm, pr, RootSelection::SmallestVolume)?;
    let v = vm.evaluate(tm, pr, RootSelection::LargestVolume)?;
    // K → 1 test on the active components (a pure fluid has K = 1 by definition).
    let max_ln_k = x
        .iter()
        .zip(&y)
        .filter(|(a, b)| **a > 0.0 && **b > 0.0)
        .map(|(a, b)| (b / a).ln().abs())
        .fold(0.0, f64::max);
    if v.molar_volume / l.molar_volume < 1.5
        || (max_ln_k < NEAR_CRITICAL_LN_K && single_active_component(model.fractions()).is_none())
    {
        return Err(Rejection::near_critical(&format!(
            "A13 near-critical diagnostic: phase volume ratio <1.5 or max |ln K| <{NEAR_CRITICAL_LN_K}"
        )));
    }
    let residual = x
        .iter()
        .zip(&y)
        .enumerate()
        .filter(|(_, (a, b))| **a > 0.0 && **b > 0.0)
        .map(|(i, (a, b))| {
            (a.ln() + l.ln_fugacity_coefficients[i] - b.ln() - v.ln_fugacity_coefficients[i]).abs()
        })
        .fold(0.0, f64::max);
    if residual > 1e-9 {
        return Err(Rejection::new(
            RejectionKind::FugacityRecheck,
            Error::NonConvergence("A13 fugacity equality failed".into()),
        ));
    }
    // TPD uses the minimum-Gibbs root. Ensure it assesses the selected phase,
    // rather than silently assessing a different root at the same composition.
    for (phase, composition, selected) in [(&lm, &x, &l), (&vm, &y, &v)] {
        let stable = phase.ln_fugacity(tm, pr, RootSelection::MinimumGibbs)?;
        let excess: f64 = composition
            .iter()
            .zip(&selected.ln_fugacity_coefficients)
            .zip(&stable.ln_fugacity_coefficients)
            .map(|((w, a), b)| w * (a - b))
            .sum();
        if excess > 1e-8 {
            return Err(Rejection::new(
                RejectionKind::MetastablePhase,
                Error::UnsupportedCapability(
                    "A13 selected phase is metastable relative to another root".into(),
                ),
            ));
        }
    }
    let phase_stability = [lm.assess_stability(tm, pr)?, vm.assess_stability(tm, pr)?];
    for phase in &phase_stability {
        if !phase.stable_in_trials {
            return Err(Rejection::new(
                RejectionKind::UnstableBoundaryPhases,
                Error::UnsupportedCapability(
                    "A13 boundary phases unstable in finite TPD trials".into(),
                ),
            ));
        }
    }
    // Extended sets: the multiphase trials (aqueous, hydrocarbon, every near-pure
    // component) must not find a third phase below either phase's tangent plane.
    let mut multiphase_stability = Vec::new();
    if model.extended() && single_active_component(model.fractions()).is_none() {
        for phase in [&lm, &vm] {
            let evidence = phase.multiphase_stability(tm, pr)?;
            if !evidence.stable_in_trials {
                return Err(Rejection::new(
                    RejectionKind::UnstableBoundaryPhases,
                    Error::UnsupportedCapability(
                        "A13 boundary phases unstable in multiphase TPD trials (a third phase, or a liquid-liquid split)".into(),
                    ),
                ));
            }
            multiphase_stability.push(evidence);
        }
    }
    Ok(SaturationPoint {
        kind,
        temperature_k: t,
        pressure_pa: p,
        liquid_composition: x,
        vapor_composition: y,
        liquid: l,
        vapor: v,
        max_ln_fugacity_residual: residual,
        iterations,
        phase_stability,
        multiphase_stability,
        second_liquid: None,
    })
}

/// Rejects a boundary below a component's triple-point temperature whose liquid
/// is essentially that component: the solid (ice, dry ice) is not modelled.
fn check_freezing<M: SaturationModel>(
    model: &M,
    t: f64,
    liquids: &[&[f64]],
) -> Result<(), Rejection> {
    let ids = model.component_ids();
    for (id, triple) in TRIPLE_POINT_TEMPERATURES_K {
        let Some(i) = ids.iter().position(|c| *c == id) else {
            continue;
        };
        if t < triple && liquids.iter().any(|x| x[i] >= FREEZING_LIQUID_FRACTION) {
            return Err(Rejection::new(
                RejectionKind::SolidPhase,
                Error::OutsideValidityDomain(format!(
                    "A13 boundary at {t} K is below the {id} triple point ({triple} K) with a liquid of at least {FREEZING_LIQUID_FRACTION} {id}; solid phases are not modelled"
                )),
            ));
        }
    }
    Ok(())
}
fn linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for k in 0..n {
        let pivot = (k..n).max_by(|i, j| a[*i][k].abs().total_cmp(&a[*j][k].abs()))?;
        if a[pivot][k].abs() < 1e-12 {
            return None;
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        let pivot_row = a[k].clone();
        for i in k + 1..n {
            let f = a[i][k] / pivot_row[k];
            for (entry, pivot) in a[i][k..].iter_mut().zip(&pivot_row[k..]) {
                *entry -= f * pivot;
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        x[i] = (b[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    Some(x)
}

#[derive(Debug, Clone)]
pub enum EnvelopeStop {
    RequestedRangeComplete,
    NearCritical {
        temperature_k: f64,
        diagnostic: String,
    },
    UnverifiedBoundary {
        temperature_k: f64,
        diagnostic: String,
    },
}
#[derive(Debug, Clone)]
pub struct PhaseEnvelope {
    pub kind: CubicSaturationKind,
    pub requested_range: SaturationOptions,
    /// Actual extrema of verified samples; None when the first sample failed.
    pub supported_temperature_range_k: Option<[f64; 2]>,
    pub supported_pressure_range_pa: Option<[f64; 2]>,
    pub points: Vec<SaturationPoint>,
    pub stop: EnvelopeStop,
}
/// Temperature-parameterized bounded sampling of one boundary. Stops at the
/// first failed/near-critical point; never bridges gaps or extrapolates branches.
pub fn trace_phase_envelope(
    model: &CubicModel,
    kind: CubicSaturationKind,
    step_k: f64,
    maximum_points: usize,
    o: SaturationOptions,
) -> Result<PhaseEnvelope, Error> {
    o.validate()?;
    if !step_k.is_finite() || step_k <= 0.0 || maximum_points == 0 || maximum_points > 1000 {
        return Err(Error::UnsupportedCapability(
            "A13 trace requires positive finite step and 1..1000 points".into(),
        ));
    }
    let mut result = PhaseEnvelope {
        kind,
        requested_range: o,
        supported_temperature_range_k: None,
        supported_pressure_range_pa: None,
        points: Vec::new(),
        stop: EnvelopeStop::RequestedRangeComplete,
    };
    let intervals = ((o.maximum_temperature_k - o.minimum_temperature_k) / step_k).ceil();
    if !intervals.is_finite() || intervals >= maximum_points as f64 {
        return Err(Error::UnsupportedCapability(
            "A13 trace point budget cannot cover requested range".into(),
        ));
    }
    // Interior samples strictly below the maximum, then the maximum exactly once.
    // `ceil` of a float quotient can otherwise yield a last interior sample equal to
    // (or within roundoff of) the maximum, which would be computed twice.
    let mut temperatures: Vec<f64> = (0..intervals as usize)
        .map(|i| o.minimum_temperature_k + i as f64 * step_k)
        .filter(|t| o.maximum_temperature_k - t > 1e-9 * step_k.max(1.0))
        .collect();
    temperatures.push(o.maximum_temperature_k);
    for t in temperatures {
        match solve(model, t, true, kind, o) {
            Ok(p) => result.points.push(p),
            Err(rejection) => {
                let message = rejection.error.to_string();
                result.stop = if rejection.kind == RejectionKind::NearCritical {
                    EnvelopeStop::NearCritical {
                        temperature_k: t,
                        diagnostic: message,
                    }
                } else {
                    EnvelopeStop::UnverifiedBoundary {
                        temperature_k: t,
                        diagnostic: message,
                    }
                };
                break;
            }
        }
    }
    if let (Some(first), Some(last)) = (result.points.first(), result.points.last()) {
        result.supported_temperature_range_k = Some([first.temperature_k, last.temperature_k]);
        result.supported_pressure_range_pa = Some([
            result
                .points
                .iter()
                .map(|p| p.pressure_pa)
                .fold(f64::INFINITY, f64::min),
            result
                .points
                .iter()
                .map(|p| p.pressure_pa)
                .fold(0.0, f64::max),
        ]);
    }
    Ok(result)
}
