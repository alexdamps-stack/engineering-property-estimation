//! GERG pressure inversion inside an explicitly supplied density bracket.
use thermo_models::{GERG_GAS_CONSTANT, GergModel, GergState};
use thermo_types::{
    CalculationError as Error, Diagnostic, DiagnosticSeverity, Pressure, Temperature,
};

/// The caller selects a density interval, not an automatically assigned phase.
/// Bounds are mol/m³. No default bracket or phase selection is provided.
#[derive(Debug, Clone, Copy)]
pub struct GergDensityOptions {
    pub minimum_molar_density: f64,
    pub maximum_molar_density: f64,
    pub max_iterations: usize,
    /// |p_eos-p_target|/max(p_target,1 Pa), positive and <=1e-7.
    pub pressure_tolerance: f64,
}

impl GergDensityOptions {
    pub const fn new(minimum_molar_density: f64, maximum_molar_density: f64) -> Self {
        Self {
            minimum_molar_density,
            maximum_molar_density,
            max_iterations: 100,
            pressure_tolerance: 1e-10,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GergDensityResult {
    pub state: GergState,
    pub target_pressure_pa: f64,
    /// Signed pressure error scaled by max(target pressure, 1 Pa).
    pub scaled_pressure_residual: f64,
    /// Interior iterations; an accepted endpoint reports zero.
    pub iterations: usize,
    pub initial_density_bracket: [f64; 2],
    pub final_density_bracket: [f64; 2],
    pub diagnostics: Vec<Diagnostic>,
}

/// Invert GERG pressure with safeguarded Newton/bisection in the caller's bracket.
/// Every visited point must have positive isothermal stiffness. The returned
/// state must also pass the model's local thermal and pressure checks.
/// This does not enumerate roots, establish interval monotonicity or determine
/// equilibrium stability. A bracket can intentionally select a metastable root.
pub fn gerg_density(
    model: &GergModel,
    temperature: Temperature,
    pressure: Pressure,
    options: GergDensityOptions,
) -> Result<GergDensityResult, Error> {
    let mut lo = options.minimum_molar_density;
    let mut hi = options.maximum_molar_density;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo >= hi
        || lo < 1e-3
        || hi > 30000.0
        || options.max_iterations == 0
        || !options.pressure_tolerance.is_finite()
        || options.pressure_tolerance <= 0.0
        || options.pressure_tolerance > 1e-7
    {
        return Err(Error::UnsupportedCapability("GERG density requires an ordered bracket within 1e-3..=30000 mol/m³, positive iteration limit and pressure tolerance <=1e-7".into()));
    }
    let target = pressure.value();
    if target > 35e6 {
        return Err(Error::OutsideValidityDomain(
            "GERG target pressure exceeds 35 MPa".into(),
        ));
    }
    let scale = target.max(1.0);
    let initial = [lo, hi];
    let residual = |rho| -> Result<(f64, f64), Error> {
        // Trial pressure may exceed the final 35 MPa envelope. The low-level
        // potential remains bounded by the model's T-density envelope.
        let a = model.helmholtz(temperature, rho)?;
        let rt = GERG_GAS_CONSTANT * temperature.value();
        let p = rho * rt * (1.0 + a.ar_d);
        let derivative = rt * (1.0 + 2.0 * a.ar_d + a.ar_dd);
        if !p.is_finite() || !derivative.is_finite() {
            return Err(Error::NumericalFailure(
                "nonfinite GERG density trial".into(),
            ));
        }
        if derivative <= 0.0 {
            return Err(Error::OutsideValidityDomain("GERG density trial has nonpositive isothermal stiffness; choose a regular root bracket".into()));
        }
        Ok(((p - target) / scale, derivative / scale))
    };
    let finish = |rho, iterations, bracket| -> Result<GergDensityResult, Error> {
        let state = model.evaluate(temperature, rho)?;
        let scaled_pressure_residual = (state.pressure_pa - target) / scale;
        if scaled_pressure_residual.abs() > options.pressure_tolerance {
            return Err(Error::NonConvergence(
                "GERG density failed final pressure check".into(),
            ));
        }
        Ok(GergDensityResult { state,target_pressure_pa:target,scaled_pressure_residual,
            iterations,initial_density_bracket:initial,final_density_bracket:bracket,
            diagnostics:vec![Diagnostic {code:"A6-GERG-BRACKETED-DENSITY".into(),severity:DiagnosticSeverity::Information,
                message:"Density root within the caller's explicit bracket; positive local stiffness and cv. Root completeness, interval monotonicity, phase identity and equilibrium stability are unassessed.".into()}] })
    };
    let (flo, _) = residual(lo)?;
    let (fhi, _) = residual(hi)?;
    if flo.abs() <= options.pressure_tolerance {
        return finish(lo, 0, [lo, hi]);
    }
    if fhi.abs() <= options.pressure_tolerance {
        return finish(hi, 0, [lo, hi]);
    }
    if flo > 0.0 || fhi < 0.0 {
        return Err(Error::OutsideValidityDomain(
            "target pressure is not bracketed by the GERG density endpoints".into(),
        ));
    }
    let mut rho = lo + 0.5 * (hi - lo);
    for iteration in 1..=options.max_iterations {
        let (f, derivative) = residual(rho)?;
        if f.abs() <= options.pressure_tolerance {
            return finish(rho, iteration, [lo, hi]);
        }
        if f < 0.0 {
            lo = rho;
        } else {
            hi = rho;
        }
        let newton = rho - f / derivative;
        let margin = 0.1 * (hi - lo);
        rho = if newton > lo + margin && newton < hi - margin {
            newton
        } else {
            lo + 0.5 * (hi - lo)
        };
        if rho <= lo || rho >= hi {
            return Err(Error::NonConvergence(
                "GERG density bracket collapsed without meeting pressure tolerance".into(),
            ));
        }
    }
    Err(Error::NonConvergence(
        "GERG density iteration limit exhausted".into(),
    ))
}
