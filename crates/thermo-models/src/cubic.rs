use std::f64::consts::SQRT_2;
use thermo_data::{Component, ComponentDatabase, REFERENCE_PRESSURE};
use thermo_types::{
    CalculationError as Error, CaloricReference, Composition, CompositionBasis, DataMetadata,
    InputError, MOLAR_GAS_CONSTANT as R, ModelMetadata, Pressure, Temperature, ValidityStatus,
};

pub const OMEGA_A: f64 = 0.457_235_528_921_382_2;
pub const OMEGA_B: f64 = 0.077_796_073_903_888_46;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CubicVariant {
    PengRobinson,
    PengRobinsonBostonMathias,
}

/// Root request by volume, deliberately without a stable-phase assertion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootSelection {
    SmallestVolume,
    LargestVolume,
    MinimumGibbs,
}

/// All quantities are molar SI unless named otherwise. Vector order matches input.
#[derive(Debug, Clone)]
pub struct HomogeneousState {
    pub temperature_k: f64,
    pub pressure_pa: f64,
    pub compressibility: f64,
    pub molar_volume: f64,
    pub mass_density: f64,
    pub ln_fugacity_coefficients: Vec<f64>,
    /// J/mol relative to ideal gas at the same T,p,x.
    pub residual_enthalpy: f64,
    /// J/(mol K) relative to ideal gas at the same T,p,x.
    pub residual_entropy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub internal_energy: f64,
    pub cp: f64,
    pub cv: f64,
    pub speed_of_sound: f64,
    /// Pa/K at fixed molar volume.
    pub dp_dt: f64,
    /// Pa/(m³/mol) at fixed T.
    pub dp_dv: f64,
    pub model: ModelMetadata,
    pub data: DataMetadata,
    pub caloric_reference: CaloricReference,
    /// Computational envelope only; physical phase validity is unassessed.
    pub validity: ValidityStatus,
    pub selection: RootSelection,
}

/// Volumetric root and fugacity coefficients only; no caloric or stability claim.
#[derive(Debug, Clone)]
pub struct FugacityState {
    pub compressibility: f64,
    /// m³/mol
    pub molar_volume: f64,
    pub ln_fugacity_coefficients: Vec<f64>,
    /// sum x_i ln(phi_i) = g_res/(RT) at this T,p,x.
    pub reduced_residual_gibbs: f64,
    /// Number of admissible real roots found at this T,p,x (1 or 3, rarely 2).
    pub root_count: usize,
}

/// Private immutable data and fully specified symmetric, constant interactions.
#[derive(Debug, Clone)]
pub struct CubicModel {
    variant: CubicVariant,
    components: Vec<Component>,
    fractions: Vec<f64>,
    kij: Vec<Vec<f64>>,
    interaction_source: String,
    data: DataMetadata,
}

#[derive(Debug)]
struct Parameters {
    a: f64,
    da: f64,
    d2a: f64,
    b: f64,
    bi: Vec<f64>,
    cross: Vec<f64>,
}

impl CubicModel {
    /// Matrix order MUST match composition. No missing-pair defaults are provided.
    pub fn new(
        database: &ComponentDatabase,
        composition: &Composition,
        variant: CubicVariant,
        kij: Vec<Vec<f64>>,
        interaction_source: impl Into<String>,
    ) -> Result<Self, Error> {
        if composition.basis() != CompositionBasis::Mole {
            return Err(Error::UnsupportedCapability(
                "convert mass composition explicitly before EOS evaluation".into(),
            ));
        }
        let components: Vec<_> = composition
            .entries()
            .iter()
            .map(|(id, _)| database.get(id.as_str()).cloned())
            .collect::<Result<_, _>>()?;
        let n = components.len();
        if kij.len() != n || kij.iter().any(|row| row.len() != n) {
            return Err(Error::MissingParameter(
                "complete interaction matrix required".into(),
            ));
        }
        let interaction_source = interaction_source.into();
        if interaction_source.trim().is_empty() {
            return Err(Error::MissingParameter(
                "interaction source or explicit approximation description required".into(),
            ));
        }
        for (i, row) in kij.iter().enumerate() {
            for (j, _) in row.iter().enumerate() {
                if !kij[i][j].is_finite()
                    || kij[i][j] != kij[j][i]
                    || (i == j && kij[i][j] != 0.0)
                    || kij[i][j] >= 1.0
                {
                    return Err(Error::UnsupportedCapability(
                        "kij must be finite, symmetric, <1, with zero diagonal".into(),
                    ));
                }
            }
        }
        Ok(Self {
            variant,
            components,
            fractions: composition.entries().iter().map(|(_, x)| *x).collect(),
            kij,
            interaction_source,
            data: database.metadata(),
        })
    }

    pub fn component_ids(&self) -> Vec<&str> {
        self.components.iter().map(Component::id).collect()
    }

    pub fn metadata(&self) -> ModelMetadata {
        ModelMetadata {
            identifier: match self.variant {
                CubicVariant::PengRobinson => "PR-exact-1976",
                CubicVariant::PengRobinsonBostonMathias => "PR-exact-BM-Clapeyron-0.6.30",
            }
            .into(),
            implementation_version: "a3-v1".into(),
            formulation: format!(
                "OmegaA={OMEGA_A}; OmegaB={OMEGA_B}; R={R}; quadratic attraction; linear b; no translation; constant symmetric kij: {}; ordered components={:?}; kij={:?}",
                self.interaction_source,
                self.component_ids(),
                self.kij
            ),
        }
    }

    fn parameters(&self, t: f64) -> Result<Parameters, Error> {
        self.parameters_at(t, &self.fractions)
    }

    fn parameters_at(&self, t: f64, fractions: &[f64]) -> Result<Parameters, Error> {
        let mut bi = Vec::new();
        let mut g = Vec::new();
        let mut dg = Vec::new();
        let mut d2g = Vec::new();
        for c in &self.components {
            let (alpha, d, d2) = alpha_derivatives(
                self.variant,
                t / c.critical_temperature(),
                c.acentric_factor(),
            )?;
            let ac = OMEGA_A * R * R * c.critical_temperature().powi(2) / c.critical_pressure();
            let tc = c.critical_temperature();
            let root = (ac * alpha).sqrt();
            g.push(root);
            dg.push(ac * d / tc / (2.0 * root));
            d2g.push(
                ac * d2 / (tc * tc) / (2.0 * root)
                    - ac * ac * d * d / (tc * tc) / (4.0 * root.powi(3)),
            );
            bi.push(OMEGA_B * R * tc / c.critical_pressure());
        }
        let mut a = 0.0;
        let mut da = 0.0;
        let mut d2a = 0.0;
        let mut cross = vec![0.0; g.len()];
        for i in 0..g.len() {
            for j in 0..g.len() {
                let weight = fractions[i] * fractions[j] * (1.0 - self.kij[i][j]);
                a += weight * g[i] * g[j];
                da += weight * (dg[i] * g[j] + g[i] * dg[j]);
                d2a += weight * (d2g[i] * g[j] + 2.0 * dg[i] * dg[j] + g[i] * d2g[j]);
                cross[i] += fractions[j] * (1.0 - self.kij[i][j]) * g[i] * g[j];
            }
        }
        let b = bi.iter().zip(fractions).map(|(b, x)| b * x).sum();
        if ![a, da, d2a, b].iter().all(|v| v.is_finite()) || a <= 0.0 || b <= 0.0 {
            return Err(Error::NumericalFailure(
                "invalid mixed cubic parameters".into(),
            ));
        }
        Ok(Parameters {
            a,
            da,
            d2a,
            b,
            bi,
            cross,
        })
    }

    fn check_domain(t: f64, p: f64) -> Result<(), Error> {
        if !(200.0..=600.0).contains(&t) || !(0.01..=100_000_000.0).contains(&p) {
            return Err(Error::OutsideValidityDomain("homogeneous EOS envelope: 200..=600 K, 0.01..=1e8 Pa; no phase/stability qualification".into()));
        }
        Ok(())
    }

    /// All admissible real Z roots, in increasing volume order, including unstable roots.
    /// Near-coalescent roots return a diagnostic error instead of an accuracy claim.
    pub fn compressibility_roots(
        &self,
        temperature: Temperature,
        pressure: Pressure,
    ) -> Result<Vec<f64>, Error> {
        let (t, p) = (temperature.value(), pressure.value());
        Self::check_domain(t, p)?;
        let par = self.parameters(t)?;
        let a = par.a * p / (R * t).powi(2);
        let b = par.b * p / (R * t);
        let mut roots = real_roots(
            b - 1.0,
            a - 3.0 * b * b - 2.0 * b,
            -a * b + b * b + b.powi(3),
            b,
            128,
        )?;
        roots.retain(|z| *z > b && z.is_finite());
        if roots.is_empty() {
            return Err(Error::NumericalFailure("no admissible Z > B root".into()));
        }
        Ok(roots)
    }

    /// Evaluate a mechanically stable homogeneous root. MinimumGibbs compares
    /// roots at this composition only, and does not establish mixture stability.
    pub fn evaluate(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        selection: RootSelection,
    ) -> Result<HomogeneousState, Error> {
        let (t, p) = (temperature.value(), pressure.value());
        let roots = self.compressibility_roots(temperature, pressure)?;
        if self.variant == CubicVariant::PengRobinsonBostonMathias
            && self
                .components
                .iter()
                .zip(&self.fractions)
                .any(|(c, x)| *x > 0.0 && t == c.critical_temperature())
        {
            return Err(Error::UnsupportedCapability("BM second temperature derivative has one-sided values at T=Tc; cp/cv not defined there".into()));
        }
        let par = self.parameters(t)?;
        let log_phi = |z: f64| ln_phi(&par, t, p, z);
        let z = select_root(&roots, selection, |z| {
            self.reduced_residual_gibbs(&log_phi(z))
        });
        let v = z * R * t / p;
        let denom = v * v + 2.0 * par.b * v - par.b * par.b;
        let pt = R / (v - par.b) - par.da / denom;
        let pv = -R * t / (v - par.b).powi(2) + par.a * (2.0 * v + 2.0 * par.b) / denom.powi(2);
        if pv >= 0.0 || -pv * v / (R * t / v) < 1e-8 {
            return Err(Error::NumericalFailure(
                "mechanically unstable or spinodal root".into(),
            ));
        }
        let q = attraction_integral(v, par.b);
        let hr = R * t * (z - 1.0) + (t * par.da - par.a) * q;
        let sr = R * (z - par.b * p / (R * t)).ln() + par.da * q;
        let mut h0 = 0.0;
        let mut s0 = -R * (p / REFERENCE_PRESSURE).ln();
        let mut cp0 = 0.0;
        let mut mass = 0.0;
        for (c, x) in self.components.iter().zip(&self.fractions) {
            let ideal = c.ideal_gas(temperature)?;
            h0 += x * ideal.enthalpy;
            s0 += x * ideal.entropy;
            if *x > 0.0 {
                s0 -= R * x * x.ln();
            }
            cp0 += x * ideal.heat_capacity;
            mass += x * c.molar_mass().value();
        }
        let cv = cp0 - R + t * par.d2a * q;
        let cp = cv - t * pt * pt / pv;
        let sound = (-v * v * pv * cp / (mass * cv)).sqrt();
        let ln_phi = log_phi(z);
        if cv <= 0.0
            || cp <= 0.0
            || ![hr, sr, h0, s0, cp, cv, sound, pt, pv, mass / v]
                .iter()
                .chain(ln_phi.iter())
                .all(|v| v.is_finite())
        {
            return Err(Error::NumericalFailure(
                "nonfinite or calorically unstable root".into(),
            ));
        }
        Ok(HomogeneousState { temperature_k:t, pressure_pa:p, compressibility:z, molar_volume:v, mass_density:mass/v,
            ln_fugacity_coefficients:ln_phi, residual_enthalpy:hr, residual_entropy:sr, enthalpy:h0+hr, entropy:s0+sr,
            internal_energy:h0+hr-p*v, cp,cv,speed_of_sound:sound,dp_dt:pt,dp_dv:pv,model:self.metadata(),data:self.data.clone(),
            caloric_reference: CaloricReference { identifier: "ideal-pure-zero-298.15K-1bar-cp-integrated-v1".into(), definition: "NASA7 cp integrals; pure ideal h=s=0 at 298.15 K, 100000 Pa; ideal mixing entropy; EOS departures at same T,p,x; J/mol and J/(mol K)".into() },
            validity: ValidityStatus::Unassessed, selection })
    }

    /// The same model and interactions at another mole-fraction vector in model order.
    /// Fractions must be finite, nonnegative and sum to one within 1e-10; zero entries
    /// keep the component present at infinite dilution for fugacity evaluation.
    pub fn with_fractions(&self, fractions: &[f64]) -> Result<Self, Error> {
        if fractions.len() != self.components.len() {
            return Err(InputError::ComponentMismatch.into());
        }
        if fractions.iter().any(|x| !x.is_finite()) {
            return Err(InputError::NonFinite("mole fraction").into());
        }
        if fractions.iter().any(|x| *x < 0.0) {
            return Err(InputError::Negative("mole fraction").into());
        }
        if (fractions.iter().sum::<f64>() - 1.0).abs() > 1e-10 {
            return Err(InputError::NotNormalized.into());
        }
        Ok(Self {
            fractions: fractions.to_vec(),
            ..self.clone()
        })
    }

    pub fn components(&self) -> &[Component] {
        &self.components
    }

    /// Mole fractions in model order.
    pub fn fractions(&self) -> &[f64] {
        &self.fractions
    }

    pub fn variant(&self) -> CubicVariant {
        self.variant
    }

    /// Mixture covolume b in m³/mol at the model fractions (temperature independent).
    pub fn covolume(&self) -> f64 {
        self.components
            .iter()
            .zip(&self.fractions)
            .map(|(c, x)| x * OMEGA_B * R * c.critical_temperature() / c.critical_pressure())
            .sum()
    }

    /// Z and ln(phi) of the requested root, without caloric evaluation. Intended for
    /// stability and flash iterations. MinimumGibbs compares roots at this composition
    /// only. A selected root with nonnegative (dp/dv)_T is rejected as mechanically
    /// unstable rather than returned.
    pub fn ln_fugacity(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        selection: RootSelection,
    ) -> Result<FugacityState, Error> {
        let (t, p) = (temperature.value(), pressure.value());
        let roots = self.compressibility_roots(temperature, pressure)?;
        let par = self.parameters(t)?;
        let z = select_root(&roots, selection, |z| {
            self.reduced_residual_gibbs(&ln_phi(&par, t, p, z))
        });
        let v = z * R * t / p;
        let denom = v * v + 2.0 * par.b * v - par.b * par.b;
        let pv = -R * t / (v - par.b).powi(2) + par.a * (2.0 * v + 2.0 * par.b) / denom.powi(2);
        if pv >= 0.0 || !pv.is_finite() {
            return Err(Error::NumericalFailure(
                "selected root is mechanically unstable".into(),
            ));
        }
        let ln_fugacity_coefficients = ln_phi(&par, t, p, z);
        if !ln_fugacity_coefficients.iter().all(|v| v.is_finite()) {
            return Err(Error::NumericalFailure("nonfinite ln(phi)".into()));
        }
        Ok(FugacityState {
            compressibility: z,
            molar_volume: v,
            root_count: roots.len(),
            reduced_residual_gibbs: self.reduced_residual_gibbs(&ln_fugacity_coefficients),
            ln_fugacity_coefficients,
        })
    }

    fn reduced_residual_gibbs(&self, ln_phi: &[f64]) -> f64 {
        ln_phi.iter().zip(&self.fractions).map(|(l, x)| l * x).sum()
    }

    /// EOS pressure and first/second volume derivatives; useful for critical constraints.
    pub fn pressure_derivatives(
        &self,
        temperature: Temperature,
        molar_volume: f64,
    ) -> Result<(f64, f64, f64), Error> {
        if !(100.0..=1000.0).contains(&temperature.value()) {
            return Err(Error::OutsideValidityDomain(
                "mechanical derivatives: 100..=1000 K".into(),
            ));
        }
        let par = self.parameters(temperature.value())?;
        let v = molar_volume;
        if !v.is_finite() || v <= par.b {
            return Err(InputError::NotPositive("molar volume minus b").into());
        }
        let d = v * v + 2.0 * par.b * v - par.b * par.b;
        let dd = 2.0 * (v + par.b);
        let rt = R * temperature.value();
        let result = (
            rt / (v - par.b) - par.a / d,
            -rt / (v - par.b).powi(2) + par.a * dd / d.powi(2),
            2.0 * rt / (v - par.b).powi(3) + 2.0 * par.a / d.powi(2)
                - 2.0 * par.a * dd * dd / d.powi(3),
        );
        if ![result.0, result.1, result.2].iter().all(|v| v.is_finite()) {
            return Err(Error::NumericalFailure(
                "pressure derivative overflow".into(),
            ));
        }
        Ok(result)
    }

    /// Extensive residual Helmholtz energy in J for volume in m³ and amount in mol.
    /// Accepts 200–600 K and finite V/n > b, independently of the TP pressure bounds.
    pub fn residual_helmholtz(
        &self,
        temperature: Temperature,
        volume: f64,
        moles: f64,
    ) -> Result<f64, Error> {
        Self::check_domain(temperature.value(), 1.0)?;
        if !volume.is_finite() || volume <= 0.0 || !moles.is_finite() || moles <= 0.0 {
            return Err(InputError::NotPositive("volume and moles").into());
        }
        let par = self.parameters(temperature.value())?;
        let v = volume / moles;
        if !v.is_finite() {
            return Err(Error::NumericalFailure("molar volume overflow".into()));
        }
        if v <= par.b {
            return Err(Error::OutsideValidityDomain("V/n must exceed b".into()));
        }
        let q = attraction_integral(v, par.b);
        let result = moles * (-R * temperature.value() * (-par.b / v).ln_1p() - par.a * q);
        if !result.is_finite() {
            return Err(Error::NumericalFailure("Helmholtz overflow".into()));
        }
        Ok(result)
    }
}

fn ln_phi(par: &Parameters, t: f64, p: f64, z: f64) -> Vec<f64> {
    let b_red = par.b * p / (R * t);
    let q = attraction_integral(z * R * t / p, par.b);
    (0..par.bi.len())
        .map(|i| {
            par.bi[i] / par.b * (z - 1.0)
                - (z - b_red).ln()
                - (2.0 * par.cross[i] - par.a * par.bi[i] / par.b) * q / (R * t)
        })
        .collect()
}

fn select_root(roots: &[f64], selection: RootSelection, gibbs: impl Fn(f64) -> f64) -> f64 {
    match selection {
        RootSelection::SmallestVolume => roots[0],
        RootSelection::LargestVolume => roots[roots.len() - 1],
        RootSelection::MinimumGibbs => roots
            .iter()
            .copied()
            .min_by(|a, b| gibbs(*a).total_cmp(&gibbs(*b)))
            .unwrap_or(roots[0]),
    }
}

// Equivalent to ln((v+(1+sqrt(2))*b)/(v+(1-sqrt(2))*b))/(2*sqrt(2)*b).
// Form the small increment directly: the ratio can round to 1 in dilute states,
// destroying the attractive contribution before the logarithm is evaluated.
fn attraction_integral(v: f64, b: f64) -> f64 {
    let scale = 2.0 * SQRT_2 * b;
    (scale / (v + (1.0 - SQRT_2) * b)).ln_1p() / scale
}

/// Alpha, d(alpha)/dTr and d²(alpha)/dTr². Low branch owns Tr=1;
/// BM's second derivative there is a one-sided value, not a two-sided derivative.
pub fn alpha_derivatives(
    variant: CubicVariant,
    tr: f64,
    omega: f64,
) -> Result<(f64, f64, f64), Error> {
    if !tr.is_finite() || !omega.is_finite() {
        return Err(InputError::NonFinite("alpha input").into());
    }
    if tr <= 0.0 {
        return Err(InputError::NotPositive("reduced temperature").into());
    }
    if !(0.01..=100.0).contains(&tr) || !(-0.5..=1.5).contains(&omega) {
        return Err(Error::OutsideValidityDomain(
            "alpha: 0.01<=Tr<=100, -0.5<=omega<=1.5".into(),
        ));
    }
    let m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    let result = if variant == CubicVariant::PengRobinsonBostonMathias && tr > 1.0 {
        let d = 1.0 + m / 2.0;
        let c = 2.0 * (1.0 - 1.0 / d);
        let alpha = (c * (1.0 - tr.powf(d))).exp();
        let log_d = -c * d * tr.powf(d - 1.0);
        (
            alpha,
            alpha * log_d,
            alpha * (log_d * log_d - c * d * (d - 1.0) * tr.powf(d - 2.0)),
        )
    } else {
        let f = 1.0 + m * (1.0 - tr.sqrt());
        (
            f * f,
            -m * f / tr.sqrt(),
            m * m / (2.0 * tr) + m * f / (2.0 * tr.powf(1.5)),
        )
    };
    if ![result.0, result.1, result.2].iter().all(|v| v.is_finite()) || result.0 <= 0.0 {
        return Err(Error::NumericalFailure(
            "alpha outside representable positive range".into(),
        ));
    }
    Ok(result)
}

// The derivative partitions the cubic into monotone intervals. Isolate every
// root above B, avoiding Cardano's cancellation near the ideal-gas limit.
fn real_roots(c2: f64, c1: f64, c0: f64, lower: f64, limit: usize) -> Result<Vec<f64>, Error> {
    let poly = |z: f64| ((z + c2) * z + c1) * z + c0;
    let mut boundaries = vec![lower];
    let derivative_disc = c2 * c2 - 3.0 * c1;
    // At the triple root, roundoff can make the derivative discriminant
    // slightly negative and otherwise hide the stationary inflection.
    if derivative_disc.abs() <= 1e-14 * (c2 * c2 + 3.0 * c1.abs()) {
        let z = -c2 / 3.0;
        let scale = z.abs().powi(3) + c2.abs() * z * z + c1.abs() * z + c0.abs();
        if z > lower && poly(z).abs() <= 1e-12 * scale {
            return Err(Error::NumericalFailure(
                "coalescent triple cubic root: critical diagnostic".into(),
            ));
        }
    }
    if derivative_disc >= 0.0 {
        let q = -c2 - derivative_disc.sqrt().copysign(c2);
        let mut stationary = if q == 0.0 {
            vec![-c2 / 3.0]
        } else {
            vec![q / 3.0, c1 / q]
        };
        stationary.sort_by(f64::total_cmp);
        for z in stationary {
            if z <= lower {
                continue;
            }
            let scale = z.abs().powi(3) + c2.abs() * z * z + c1.abs() * z + c0.abs();
            if poly(z).abs() <= 1e-12 * scale {
                return Err(Error::NumericalFailure(
                    "coalescent admissible cubic roots: critical/spinodal diagnostic".into(),
                ));
            }
            if boundaries.last().is_none_or(|previous| *previous != z) {
                boundaries.push(z);
            }
        }
    }
    boundaries.push(1.0 + c2.abs().max(c1.abs()).max(c0.abs()));
    let mut roots = Vec::new();
    for pair in boundaries.windows(2) {
        let (mut lo, mut hi) = (pair[0], pair[1]);
        let mut flo = poly(lo);
        let fhi = poly(hi);
        if flo.signum() == fhi.signum() {
            continue;
        }
        let mut converged = false;
        for _ in 0..limit {
            let mid = lo + (hi - lo) / 2.0;
            let fm = poly(mid);
            if fm == 0.0 || hi - lo <= 2e-14 * mid.abs().max(lower) || mid == lo || mid == hi {
                lo = mid;
                hi = mid;
                converged = true;
                break;
            }
            if fm.signum() == flo.signum() {
                lo = mid;
                flo = fm;
            } else {
                hi = mid;
            }
        }
        if !converged {
            return Err(Error::NonConvergence(format!(
                "cubic bisection exhausted {limit} iterations"
            )));
        }
        let root = lo + (hi - lo) / 2.0;
        let scale = root.abs().powi(3) + c2.abs() * root * root + c1.abs() * root + c0.abs();
        if !root.is_finite() || poly(root).abs() > 1e-11 * scale {
            return Err(Error::NumericalFailure("cubic root residual failed".into()));
        }
        roots.push(root);
    }
    Ok(roots)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bounded_root_solver_reports_exhaustion() {
        assert!(matches!(
            real_roots(-1.0, 0.1, -0.01, 0.0, 0),
            Err(Error::NonConvergence(_))
        ));
    }
}
