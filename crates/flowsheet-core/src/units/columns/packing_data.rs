//! Packing library, mass-transfer and hydraulic correlations of rate-based
//! packed sections (#241). Every constant here was checked against a
//! readable source; see `docs/a17-columns-rate-based.md` § "Verification
//! of constants".
//!
//! * **Packing data**: Perry's Chemical Engineers' Handbook, 8th ed. (2008),
//!   Table 14-13 "Characteristics of random packings" (p. 14-60): specific
//!   area, void fraction, packing factor F_p (from Kister, *Distillation
//!   Design*, 1992) and dry packing factor F_pd (from Robbins 1991).
//! * **Mass transfer**: Onda, Takeuchi and Okumoto (1968), J. Chem. Eng.
//!   Japan 1(1), 56–62, eqs. (1)–(3): wetted area, liquid- and gas-side
//!   coefficients; constant 5.23, or 2.00 for packings smaller than 15 mm.
//!   The critical surface tension σc of the packing material is a user
//!   input: no σc table could be verified (the paper takes it from Onda,
//!   Takeuchi and Koyama, Kagaku Kōgaku 31 (1967) 126, not available here).
//! * **Pressure drop**: Robbins (1991), Chem. Eng. Prog. 87(5), 87–90, as
//!   given in Perry's 8th ed., eqs. 14-145 to 14-151 (p. 14-59), with
//!   C3 = 7.4e-8 and C4 = 2.7e-5; "not suitable above 3 bar absolute"
//!   (p. 14-61), which is enforced.
//! * **Flooding**: Kister and Gill flood pressure drop, Perry's 8th ed.
//!   eq. 14-142, \(\Delta P_\text{flood} = 0.12 F_p^{0.7}\) inH2O/ft with
//!   F_p in ft⁻¹, used with the Robbins pressure drop as Perry's recommends
//!   (p. 14-61); not for F_p > 60 ft⁻¹ (p. 14-58), which is enforced.
//! * **Minimum wetting**: Perry's 8th ed. eq. 14-166 with the Glitsch rules
//!   of thumb of Table 14-16 (p. 14-79): Q_MW = Q_table (200/a_p)^0.5,
//!   doubled for liquids above 25 mN/m as recommended there for water.
//!
//! Billet–Schultes pressure drop uses its own Pall-ring constants from Table
//! 2a of Billet and Schultes (1999), with the pre-loading holdup equation.
//! The reported pressure-drop gas-density range is 0.06–28 kg/m³; the
//! separate loading/flooding data reach only 1.37 kg/m³. Capacity therefore
//! remains the explicit Kister–Gill method unless a vendor rating is supplied.
//! Bravo–Rocha–Fair (structured packings) is not implemented.

/// Standard acceleration of gravity, m/s².
pub const GRAVITY_M_S2: f64 = 9.806_65;
/// Largest absolute pressure of the Robbins correlation, Pa (Perry's 8th ed.
/// p. 14-61: "not suitable above 3 bar absolute").
pub const ROBBINS_MAX_PRESSURE_PA: f64 = 3.0e5;
/// Gas-density limit of the reported Billet–Schultes pressure-drop data,
/// kg/m³ (1999 paper, Table 1).
pub const BILLET_SCHULTES_MAX_GAS_DENSITY_KG_M3: f64 = 28.0;
/// Largest packing factor of the Kister–Gill flood equation, 1/ft (Perry's
/// 8th ed. p. 14-58).
pub const KISTER_GILL_MAX_PACKING_FACTOR_FT: f64 = 60.0;
/// Surface tension above which the minimum wetting rate of Table 14-16 is
/// doubled (Perry's 8th ed. p. 14-79), N/m.
pub const MINIMUM_WETTING_LOW_TENSION_N_M: f64 = 0.025;
const FT: f64 = 0.3048;
const LB: f64 = 0.453_592_37;
const ATMOSPHERE_PA: f64 = 101_325.0;

/// Packing surface of the Glitsch minimum-wetting rules (Perry's 8th ed.
/// Table 14-16).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum WettingSurface {
    /// Unglazed ceramic (chemical stoneware), 0.5 m³/(m² h).
    UnglazedCeramic,
    /// Oxidized metal (carbon steel, copper), 0.7 m³/(m² h).
    OxidizedMetal,
    /// Surface-treated metal (etched stainless steel), 1.0 m³/(m² h).
    SurfaceTreatedMetal,
    /// Glazed ceramic, 2.0 m³/(m² h).
    GlazedCeramic,
    /// Glass, 2.5 m³/(m² h).
    Glass,
    /// Bright metal (stainless steel, tantalum, other alloys), 3.0 m³/(m² h).
    #[default]
    BrightMetal,
    /// PVC-CPVC, 3.5 m³/(m² h).
    PvcCpvc,
    /// Polypropylene, 4.0 m³/(m² h).
    Polypropylene,
    /// Fluoropolymers (PTFE type), 5.0 m³/(m² h).
    Fluoropolymer,
}

impl WettingSurface {
    /// Every surface, in table order.
    pub const ALL: [WettingSurface; 9] = [
        Self::UnglazedCeramic,
        Self::OxidizedMetal,
        Self::SurfaceTreatedMetal,
        Self::GlazedCeramic,
        Self::Glass,
        Self::BrightMetal,
        Self::PvcCpvc,
        Self::Polypropylene,
        Self::Fluoropolymer,
    ];

    /// Table 14-16 minimum wetting rate for CMR packings with a > 150
    /// m²/m³, m³/(m² h).
    pub fn table_rate_m_h(self) -> f64 {
        match self {
            Self::UnglazedCeramic => 0.5,
            Self::OxidizedMetal => 0.7,
            Self::SurfaceTreatedMetal => 1.0,
            Self::GlazedCeramic => 2.0,
            Self::Glass => 2.5,
            Self::BrightMetal => 3.0,
            Self::PvcCpvc => 3.5,
            Self::Polypropylene => 4.0,
            Self::Fluoropolymer => 5.0,
        }
    }

    /// Display name.
    pub fn name(self) -> &'static str {
        match self {
            Self::UnglazedCeramic => "Unglazed ceramic",
            Self::OxidizedMetal => "Oxidized metal",
            Self::SurfaceTreatedMetal => "Surface-treated metal",
            Self::GlazedCeramic => "Glazed ceramic",
            Self::Glass => "Glass",
            Self::BrightMetal => "Bright metal (stainless steel)",
            Self::PvcCpvc => "PVC-CPVC",
            Self::Polypropylene => "Polypropylene",
            Self::Fluoropolymer => "Fluoropolymer",
        }
    }
}

/// Packings of the library.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Packing {
    /// Metal Pall rings, 25 mm.
    #[default]
    MetalPallRing25,
    /// Metal Pall rings, 50 mm.
    MetalPallRing50,
}

/// Cited data of one packing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackingData {
    /// Display name.
    pub name: &'static str,
    /// Nominal size, m (the Onda characteristic length d_p).
    pub nominal_size_m: f64,
    /// Specific surface area a, m²/m³.
    pub specific_area_m2_m3: f64,
    /// Void fraction ε.
    pub void_fraction: f64,
    /// Packing factor F_p, 1/m (Kister).
    pub packing_factor_m: f64,
    /// Dry packing factor F_pd, 1/m (Robbins).
    pub dry_packing_factor_m: f64,
    /// Source of the values.
    pub source: &'static str,
}

/// Billet–Schultes pressure-drop parameters for a specific random packing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BilletSchultesData {
    /// Geometric surface area, m²/m³.
    pub specific_area_m2_m3: f64,
    /// Void fraction.
    pub void_fraction: f64,
    /// Dry-bed resistance constant C_P,0.
    pub dry_resistance: f64,
}

const PERRY_14_13: &str = "Perry's Chemical Engineers' Handbook, 8th ed. (2008), Table 14-13, p. 14-60 (metal Pall rings)";

impl Packing {
    /// Every packing of the library.
    pub const ALL: [Packing; 2] = [Packing::MetalPallRing25, Packing::MetalPallRing50];

    /// Cited data of the packing.
    pub fn data(self) -> PackingData {
        match self {
            Self::MetalPallRing25 => PackingData {
                name: "Pall rings, metal, 25 mm",
                nominal_size_m: 0.025,
                specific_area_m2_m3: 205.0,
                void_fraction: 0.94,
                packing_factor_m: 183.0,
                dry_packing_factor_m: 174.0,
                source: PERRY_14_13,
            },
            Self::MetalPallRing50 => PackingData {
                name: "Pall rings, metal, 50 mm",
                nominal_size_m: 0.050,
                specific_area_m2_m3: 105.0,
                void_fraction: 0.96,
                packing_factor_m: 89.0,
                dry_packing_factor_m: 79.0,
                source: PERRY_14_13,
            },
        }
    }

    /// Billet–Schultes (1999) Table 2a constants. These are the source's
    /// geometry values and intentionally remain separate from Perry's
    /// Onda/Robbins packing data above.
    pub fn billet_schultes_data(self) -> BilletSchultesData {
        match self {
            Self::MetalPallRing25 => BilletSchultesData {
                specific_area_m2_m3: 223.5,
                void_fraction: 0.954,
                dry_resistance: 0.957,
            },
            Self::MetalPallRing50 => BilletSchultesData {
                specific_area_m2_m3: 112.6,
                void_fraction: 0.951,
                dry_resistance: 0.763,
            },
        }
    }
}

/// Mass-transfer correlation of a rate-based packed section.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum MassTransferCorrelation {
    /// Onda, Takeuchi and Okumoto (1968), random packings.
    #[default]
    Onda1968,
}

/// Pressure-drop model used independently of the flooding-capacity model.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum PressureDropCorrelation {
    /// Robbins (1991), limited to 3 bar absolute.
    #[default]
    Robbins1991,
    /// Billet and Schultes (1999), random packings with tabulated constants.
    BilletSchultes1999,
    /// Measured or vendor-rated total-bed pressure drop.
    VendorRating,
}

/// Capacity model used independently of the pressure-drop model.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum FloodingCorrelation {
    /// Kister and Gill flood pressure-drop criterion.
    #[default]
    KisterGill,
    /// Measured or vendor-rated gas flow at flooding.
    VendorRating,
}

/// Phase properties the correlations need, SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FilmProperties {
    /// Superficial mass flux, kg/(m² s).
    pub mass_flux_kg_m2_s: f64,
    /// Mass density, kg/m³.
    pub density_kg_m3: f64,
    /// Dynamic viscosity, Pa s.
    pub viscosity_pa_s: f64,
}

/// Onda eq. (1) wetted area a_w, m²/m³, for the liquid properties, its
/// surface tension σ_L and the packing's critical surface tension σc (N/m).
pub fn onda_wetted_area(
    packing: &PackingData,
    liquid: &FilmProperties,
    surface_tension_n_m: f64,
    critical_surface_tension_n_m: f64,
) -> f64 {
    let a = packing.specific_area_m2_m3;
    let l = liquid.mass_flux_kg_m2_s;
    let rho = liquid.density_kg_m3;
    let re = l / (a * liquid.viscosity_pa_s);
    let fr = l * l * a / (rho * rho * GRAVITY_M_S2);
    let we = l * l / (rho * surface_tension_n_m * a);
    let ratio = critical_surface_tension_n_m / surface_tension_n_m;
    a * (1.0 - (-1.45 * ratio.powf(0.75) * re.powf(0.1) * fr.powf(-0.05) * we.powf(0.2)).exp())
}

/// Onda eq. (2) liquid-side coefficient k_L, m/s, for diffusivity `d_m2_s`
/// and the wetted area a_w.
pub fn onda_liquid_coefficient(
    packing: &PackingData,
    liquid: &FilmProperties,
    wetted_area_m2_m3: f64,
    d_m2_s: f64,
) -> f64 {
    let mu = liquid.viscosity_pa_s;
    let rho = liquid.density_kg_m3;
    let re = liquid.mass_flux_kg_m2_s / (wetted_area_m2_m3 * mu);
    let sc = mu / (rho * d_m2_s);
    0.0051
        * re.powf(2.0 / 3.0)
        * sc.powf(-0.5)
        * (packing.specific_area_m2_m3 * packing.nominal_size_m).powf(0.4)
        * (mu * GRAVITY_M_S2 / rho).cbrt()
}

/// Onda eq. (3) gas-side coefficient, m/s (k_G·RT of the original pressure
/// basis), for diffusivity `d_m2_s`: constant 5.23, or 2.00 for packings
/// smaller than 15 mm.
pub fn onda_gas_coefficient(packing: &PackingData, gas: &FilmProperties, d_m2_s: f64) -> f64 {
    let a = packing.specific_area_m2_m3;
    let constant = if packing.nominal_size_m < 0.015 {
        2.00
    } else {
        5.23
    };
    let re = gas.mass_flux_kg_m2_s / (a * gas.viscosity_pa_s);
    let sc = gas.viscosity_pa_s / (gas.density_kg_m3 * d_m2_s);
    constant * a * d_m2_s * re.powf(0.7) * sc.cbrt() * (a * packing.nominal_size_m).powi(-2)
}

/// lb/(ft² h) per kg/(m² s).
const LB_FT2_H_PER_KG_M2_S: f64 = 3600.0 / LB * FT * FT;
/// lb/ft³ per kg/m³.
const LB_FT3_PER_KG_M3: f64 = FT * FT * FT / LB;
/// Pa/m per inH2O/ft.
const PA_M_PER_INH2O_FT: f64 = 249.088_908_333 / FT;

/// Robbins loading factors (G_f, L_f) in lb/(ft² h), Perry's eqs. 14-148 to
/// 14-151.
pub fn robbins_loading_factors(
    packing: &PackingData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
    pressure_pa: f64,
) -> (f64, f64) {
    let fpd = packing.dry_packing_factor_m * FT;
    let rho_g = gas.density_kg_m3 * LB_FT3_PER_KG_M3;
    let rho_l = liquid.density_kg_m3 * LB_FT3_PER_KG_M3;
    let mu_cp = liquid.viscosity_pa_s * 1e3;
    // F_s = U ρ_G^0.5 in ft/s (lb/ft³)^0.5.
    let fs = gas.mass_flux_kg_m2_s * LB_FT2_H_PER_KG_M2_S / 3600.0 / rho_g.sqrt();
    let mut gf = 986.0 * fs * (fpd / 20.0).sqrt();
    if pressure_pa > ATMOSPHERE_PA {
        // Eq. 14-150, operating pressures above atmospheric.
        gf *= 10f64.powf(0.3 * rho_g);
    }
    let l = liquid.mass_flux_kg_m2_s * LB_FT2_H_PER_KG_M2_S;
    let lf = if fpd > 200.0 {
        l * (62.4 / rho_l) * (fpd / 20.0).sqrt() * mu_cp.powf(0.2)
    } else if fpd < 15.0 {
        l * (62.4 / rho_l) * (20.0 / fpd).sqrt() * mu_cp.powf(0.1)
    } else {
        l * (62.4 / rho_l) * (fpd / 20.0).sqrt() * mu_cp.powf(0.1)
    };
    (gf, lf)
}

fn robbins_inh2o_ft(gf: f64, lf: f64) -> f64 {
    let dry = 7.4e-8 * gf * gf * 10f64.powf(2.7e-5 * lf);
    dry + 0.4 * (lf / 20_000.0).powf(0.1) * dry.powi(4)
}

fn check_robbins(pressure_pa: f64) -> Result<(), String> {
    if pressure_pa > ROBBINS_MAX_PRESSURE_PA {
        Err(format!(
            "Robbins pressure-drop correlation is not suitable above 3 bar absolute (column pressure {pressure_pa} Pa)"
        ))
    } else {
        Ok(())
    }
}

/// Robbins (1991) pressure drop per height of packing, Pa/m.
pub fn robbins_pressure_drop_pa_m(
    packing: &PackingData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
    pressure_pa: f64,
) -> Result<f64, String> {
    check_robbins(pressure_pa)?;
    let (gf, lf) = robbins_loading_factors(packing, gas, liquid, pressure_pa);
    Ok(robbins_inh2o_ft(gf, lf) * PA_M_PER_INH2O_FT)
}

/// Billet–Schultes (1999) irrigated pressure drop per packed height.
///
/// Uses the published pre-loading holdup equation (11), dry resistance
/// equation (25), wall correction (24), and irrigated pressure-drop
/// equations (28)–(30). The packing constants are kept distinct from the
/// Perry data used by the Onda and Robbins methods.
pub fn billet_schultes_pressure_drop_pa_m(
    packing: BilletSchultesData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
    column_diameter_m: f64,
) -> Result<f64, String> {
    if gas.density_kg_m3 > BILLET_SCHULTES_MAX_GAS_DENSITY_KG_M3 {
        return Err(format!(
            "Billet–Schultes pressure-drop gas density {:.3} kg/m³ exceeds the reported data limit of {:.1} kg/m³",
            gas.density_kg_m3, BILLET_SCHULTES_MAX_GAS_DENSITY_KG_M3
        ));
    }
    let a = packing.specific_area_m2_m3;
    let eps = packing.void_fraction;
    let u_g = gas.mass_flux_kg_m2_s / gas.density_kg_m3;
    let u_l = liquid.mass_flux_kg_m2_s / liquid.density_kg_m3;
    let d_p = 6.0 * (1.0 - eps) / a;
    let wall_inverse = 1.0 + (2.0 / 3.0) / (1.0 - eps) * d_p / column_diameter_m;
    let nu_g = gas.viscosity_pa_s / gas.density_kg_m3;
    let re_g = u_g * d_p / ((1.0 - eps) * nu_g) / wall_inverse;
    let h_l =
        (12.0 * liquid.viscosity_pa_s * u_l * a * a / (GRAVITY_M_S2 * liquid.density_kg_m3)).cbrt();
    if !(u_g.is_finite()
        && u_g > 0.0
        && u_l.is_finite()
        && u_l > 0.0
        && re_g.is_finite()
        && re_g > 0.0
        && h_l.is_finite()
        && h_l < eps)
    {
        return Err(
            "Billet–Schultes pressure drop is outside its pre-loading holdup domain".into(),
        );
    }
    let f_v = u_g * gas.density_kg_m3.sqrt();
    let fr_l = u_l * u_l * a / GRAVITY_M_S2;
    let c1 = 13_300.0 / a.powf(1.5);
    let psi_l = packing.dry_resistance
        * (64.0 / re_g + 1.8 / re_g.powf(0.08))
        * ((eps - h_l) / eps).powf(1.5)
        * (c1 * fr_l.sqrt()).exp();
    let drop = psi_l * a / (eps - h_l).powi(3) * f_v * f_v * wall_inverse / 2.0;
    if drop.is_finite() && drop >= 0.0 {
        Ok(drop)
    } else {
        Err("Billet–Schultes pressure drop is nonfinite".into())
    }
}

/// Kister–Gill flood pressure drop, Pa/m (Perry's eq. 14-142).
pub fn flood_pressure_drop_pa_m(packing: &PackingData) -> Result<f64, String> {
    let fp = packing.packing_factor_m * FT;
    if fp > KISTER_GILL_MAX_PACKING_FACTOR_FT {
        return Err(format!(
            "Kister–Gill flood pressure drop is not suitable for packing factors above 60 1/ft ({} has {fp:.0} 1/ft)",
            packing.name
        ));
    }
    Ok(0.12 * fp.powf(0.7) * PA_M_PER_INH2O_FT)
}

/// Fraction of flooding: the gas mass flux divided by the gas mass flux at
/// which the Robbins pressure drop equals the Kister–Gill flood pressure
/// drop at the same liquid load and phase properties.
pub fn flooding_fraction(
    packing: &PackingData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
    pressure_pa: f64,
) -> Result<f64, String> {
    check_robbins(pressure_pa)?;
    let target = flood_pressure_drop_pa_m(packing)?;
    let at = |flux: f64| {
        let trial = FilmProperties {
            mass_flux_kg_m2_s: flux,
            ..*gas
        };
        let (gf, lf) = robbins_loading_factors(packing, &trial, liquid, pressure_pa);
        robbins_inh2o_ft(gf, lf) * PA_M_PER_INH2O_FT
    };
    let (mut lo, mut hi) = (0.0, gas.mass_flux_kg_m2_s.max(1e-6));
    while at(hi) < target {
        lo = hi;
        hi *= 2.0;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if at(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= 1e-14 * hi {
            break;
        }
    }
    Ok(gas.mass_flux_kg_m2_s / (0.5 * (lo + hi)))
}

/// Minimum wetting rate, m³/(m² h) of superficial liquid velocity (Perry's
/// eq. 14-166 and Table 14-16), doubled above 25 mN/m.
pub fn minimum_wetting_rate_m_h(
    packing: &PackingData,
    surface: WettingSurface,
    surface_tension_n_m: f64,
) -> f64 {
    let base = surface.table_rate_m_h() * (200.0 / packing.specific_area_m2_m3).sqrt();
    if surface_tension_n_m > MINIMUM_WETTING_LOW_TENSION_N_M {
        2.0 * base
    } else {
        base
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water() -> FilmProperties {
        FilmProperties {
            mass_flux_kg_m2_s: 5.0,
            density_kg_m3: 998.0,
            viscosity_pa_s: 1.0e-3,
        }
    }

    fn air(flux: f64) -> FilmProperties {
        FilmProperties {
            mass_flux_kg_m2_s: flux,
            density_kg_m3: 1.2,
            viscosity_pa_s: 1.8e-5,
        }
    }

    #[test]
    fn billet_schultes_pressure_drop_uses_source_constants_and_density() {
        let constants = Packing::MetalPallRing25.billet_schultes_data();
        assert_eq!(
            (
                constants.specific_area_m2_m3,
                constants.void_fraction,
                constants.dry_resistance,
            ),
            (223.5, 0.954, 0.957)
        );
        let dense_gas = FilmProperties {
            density_kg_m3: 8.0,
            viscosity_pa_s: 2.0e-5,
            ..air(1.0)
        };
        let drop =
            billet_schultes_pressure_drop_pa_m(constants, &dense_gas, &water(), 1.0).unwrap();
        let higher_load = billet_schultes_pressure_drop_pa_m(
            constants,
            &FilmProperties {
                mass_flux_kg_m2_s: 1.2,
                ..dense_gas
            },
            &water(),
            1.0,
        )
        .unwrap();
        assert!(drop.is_finite() && drop > 0.0);
        assert!(higher_load > drop);
    }

    #[test]
    fn library_entries_match_perrys_table() {
        let p25 = Packing::MetalPallRing25.data();
        assert_eq!(
            (
                p25.specific_area_m2_m3,
                p25.void_fraction,
                p25.packing_factor_m
            ),
            (205.0, 0.94, 183.0)
        );
        assert_eq!(p25.dry_packing_factor_m, 174.0);
        let p50 = Packing::MetalPallRing50.data();
        assert_eq!(
            (
                p50.specific_area_m2_m3,
                p50.void_fraction,
                p50.packing_factor_m
            ),
            (105.0, 0.96, 89.0)
        );
        assert_eq!(p50.dry_packing_factor_m, 79.0);
        for packing in Packing::ALL {
            assert!(packing.data().source.contains("Table 14-13"));
        }
    }

    #[test]
    fn robbins_reproduces_perrys_example_13() {
        // Perry's 8th ed. Example 13 (p. 14-62/63): air/water on 2-in metal
        // Pall rings, G = 1500 and L = 9000 lb/(ft² h), ρG = 0.074 and
        // ρL = 62.4 lb/ft³, μL = 1.0 cP, Fpd = 24 1/ft: Gf = 1653,
        // Lf = 9859 and ΔP = 0.40 inH2O/ft read from the Robbins chart.
        let packing = PackingData {
            dry_packing_factor_m: 24.0 / FT,
            ..Packing::MetalPallRing50.data()
        };
        let gas = FilmProperties {
            mass_flux_kg_m2_s: 1500.0 / LB_FT2_H_PER_KG_M2_S,
            density_kg_m3: 0.074 / LB_FT3_PER_KG_M3,
            viscosity_pa_s: 1.8e-5,
        };
        let liquid = FilmProperties {
            mass_flux_kg_m2_s: 9000.0 / LB_FT2_H_PER_KG_M2_S,
            density_kg_m3: 62.4 / LB_FT3_PER_KG_M3,
            viscosity_pa_s: 1.0e-3,
        };
        let (gf, lf) = robbins_loading_factors(&packing, &gas, &liquid, ATMOSPHERE_PA);
        // The example rounds Fs to 1.53 (exact 1.5317): Gf = 1654.4 here.
        assert!((gf - 1653.0).abs() < 2.0, "{gf}");
        assert!((lf - 9859.0).abs() < 2.0, "{lf}");
        let dp = robbins_pressure_drop_pa_m(&packing, &gas, &liquid, ATMOSPHERE_PA).unwrap()
            / PA_M_PER_INH2O_FT;
        // The equations give 0.38; the example reads 0.40 from the chart.
        assert!((dp - 0.40).abs() < 0.03, "{dp}");
        assert!(robbins_pressure_drop_pa_m(&packing, &gas, &liquid, 3.5e5).is_err());
    }

    #[test]
    fn kister_gill_flood_and_flooding_fraction() {
        let d = Packing::MetalPallRing50.data();
        let flood = flood_pressure_drop_pa_m(&d).unwrap() / PA_M_PER_INH2O_FT;
        assert!((flood - 0.12 * (89.0 * FT).powf(0.7)).abs() < 1e-12);
        let f = flooding_fraction(&d, &air(1.0), &water(), ATMOSPHERE_PA).unwrap();
        assert!(f > 0.0 && f < 1.0, "{f}");
        let at_flood =
            robbins_pressure_drop_pa_m(&d, &air(1.0 / f), &water(), ATMOSPHERE_PA).unwrap();
        assert!((at_flood / flood_pressure_drop_pa_m(&d).unwrap() - 1.0).abs() < 1e-9);
        let coarse = PackingData {
            packing_factor_m: 61.0 / FT,
            ..d
        };
        assert!(flood_pressure_drop_pa_m(&coarse).is_err());
    }

    #[test]
    fn onda_area_and_coefficients() {
        let d = Packing::MetalPallRing25.data();
        let aw = onda_wetted_area(&d, &water(), 0.072, 0.075);
        assert!(aw > 0.0 && aw < d.specific_area_m2_m3);
        // A better-wetting packing (higher σc) wets more area.
        assert!(onda_wetted_area(&d, &water(), 0.072, 0.09) > aw);
        let kl = onda_liquid_coefficient(&d, &water(), aw, 2e-9);
        let kg = onda_gas_coefficient(&d, &air(1.0), 2e-5);
        assert!((1e-5..1e-3).contains(&kl), "{kl}");
        assert!((1e-3..1e-1).contains(&kg), "{kg}");
    }

    #[test]
    fn minimum_wetting_follows_the_glitsch_rule() {
        let d = Packing::MetalPallRing50.data();
        let low = minimum_wetting_rate_m_h(&d, WettingSurface::BrightMetal, 0.02);
        assert!((low - 3.0 * (200.0_f64 / 105.0).sqrt()).abs() < 1e-12);
        assert_eq!(
            minimum_wetting_rate_m_h(&d, WettingSurface::BrightMetal, 0.07),
            2.0 * low
        );
    }
}
