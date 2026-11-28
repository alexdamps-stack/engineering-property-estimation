//! Packing library, mass-transfer and hydraulic correlations of rate-based
//! packed sections (#241). See `docs/a17-columns-rate-based.md`.
//!
//! * **Packing data**: geometry (specific area, void fraction) and the
//!   tower packing factor of metal Pall rings from the manufacturer's
//!   published data, Raschig GmbH, *RASCHIG PALL-Ring* brochure 350
//!   (edition v2025010), "Operational data for PALL-Rings made of metal"
//!   and "Tower packing factor", Table 1.
//! * **Mass transfer**: Onda, Takeuchi and Okumoto (1968), J. Chem. Eng.
//!   Japan 1, 56–62, for random packings (wetted area, liquid- and gas-side
//!   coefficients), with the critical surface tensions of the packing
//!   materials given there.
//! * **Pressure drop**: Robbins (1991), Chem. Eng. Prog. 87(5), 87–90, the
//!   analytic form of the generalized pressure-drop correlation (GPDC).
//! * **Flooding**: the flood pressure drop of Kister and Gill (1991), Chem.
//!   Eng. Prog. 87(2), 32–42, \(\Delta P_\text{flood} = 0.115 F_p^{0.7}\)
//!   inH2O/ft with \(F_p\) in ft⁻¹; the flooding fraction is the gas load
//!   divided by the gas load at which the Robbins pressure drop reaches it
//!   at the same liquid load.
//! * **Minimum wetting**: Morris and Jackson (1953), *Absorption Towers*
//!   (Butterworths), minimum wetting rate 0.08 m³/(m h) for rings up to
//!   75 mm, as recommended in Sinnott, *Coulson & Richardson's Chemical
//!   Engineering* Vol. 6, 4th ed. (2005), §11.14.
//!
//! Billet–Schultes (1999) and Bravo–Rocha–Fair (structured packings) are not
//! implemented: their packing-specific constants have not been taken from a
//! verified source in this batch.

/// Standard acceleration of gravity, m/s².
pub const GRAVITY_M_S2: f64 = 9.806_65;
/// Morris–Jackson minimum wetting rate for rings up to 75 mm, m³/(m s)
/// (0.08 m³/(m h)).
pub const MINIMUM_WETTING_RATE_M2_S: f64 = 0.08 / 3600.0;

/// Packing material with its Onda (1968) critical surface tension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackingMaterial {
    /// Steel, σc = 0.075 N/m (Onda et al. 1968).
    Steel,
}

impl PackingMaterial {
    /// Critical surface tension σc of the material, N/m.
    pub fn critical_surface_tension_n_m(self) -> f64 {
        match self {
            Self::Steel => 0.075,
        }
    }
}

/// Packings of the library.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Packing {
    /// Metal Pall rings, 25 mm (25 × 0.5/0.6 mm wall).
    #[default]
    MetalPallRing25,
    /// Metal Pall rings, 50 mm (50 × 0.8 mm wall).
    MetalPallRing50,
}

/// Cited data of one packing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PackingData {
    /// Display name.
    pub name: &'static str,
    /// Random (dumped) packing.
    pub random: bool,
    /// Nominal size, m (the Onda characteristic length d_p).
    pub nominal_size_m: f64,
    /// Specific surface area a, m²/m³.
    pub specific_area_m2_m3: f64,
    /// Void fraction ε.
    pub void_fraction: f64,
    /// Tower packing factor F_p, 1/m.
    pub packing_factor_m: f64,
    /// Material of the packing surface.
    pub material: PackingMaterial,
    /// Source of the values.
    pub source: &'static str,
}

const RASCHIG_PALL: &str = "Raschig GmbH, RASCHIG PALL-Ring brochure 350 (v2025010): operational data for metal PALL-Rings and tower packing factor Table 1";

impl Packing {
    /// Every packing of the library.
    pub const ALL: [Packing; 2] = [Packing::MetalPallRing25, Packing::MetalPallRing50];

    /// Cited data of the packing.
    pub fn data(self) -> PackingData {
        match self {
            Self::MetalPallRing25 => PackingData {
                name: "Pall rings, metal, 25 mm",
                random: true,
                nominal_size_m: 0.025,
                specific_area_m2_m3: 215.0,
                void_fraction: 0.95,
                packing_factor_m: 157.0,
                material: PackingMaterial::Steel,
                source: RASCHIG_PALL,
            },
            Self::MetalPallRing50 => PackingData {
                name: "Pall rings, metal, 50 mm",
                random: true,
                nominal_size_m: 0.050,
                specific_area_m2_m3: 105.0,
                void_fraction: 0.96,
                packing_factor_m: 66.0,
                material: PackingMaterial::Steel,
                source: RASCHIG_PALL,
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

/// Onda wetted (effective interfacial) area a_w, m²/m³, for the liquid
/// properties and surface tension σ_L in N/m.
pub fn onda_wetted_area(
    packing: &PackingData,
    liquid: &FilmProperties,
    surface_tension_n_m: f64,
) -> f64 {
    let a = packing.specific_area_m2_m3;
    let l = liquid.mass_flux_kg_m2_s;
    let rho = liquid.density_kg_m3;
    let re = l / (a * liquid.viscosity_pa_s);
    let fr = l * l * a / (rho * rho * GRAVITY_M_S2);
    let we = l * l / (rho * surface_tension_n_m * a);
    let ratio = packing.material.critical_surface_tension_n_m() / surface_tension_n_m;
    a * (1.0 - (-1.45 * ratio.powf(0.75) * re.powf(0.1) * fr.powf(-0.05) * we.powf(0.2)).exp())
}

/// Onda liquid-side coefficient k_L, m/s, for diffusivity `d_m2_s` and the
/// wetted area a_w.
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

/// Onda gas-side coefficient k_G, m/s (k_G·RT of the original pressure
/// basis), for diffusivity `d_m2_s`. The constant is 5.23 for packings
/// larger than 15 mm and 2.00 otherwise.
pub fn onda_gas_coefficient(packing: &PackingData, gas: &FilmProperties, d_m2_s: f64) -> f64 {
    let a = packing.specific_area_m2_m3;
    let constant = if packing.nominal_size_m > 0.015 {
        5.23
    } else {
        2.00
    };
    let re = gas.mass_flux_kg_m2_s / (a * gas.viscosity_pa_s);
    let sc = gas.viscosity_pa_s / (gas.density_kg_m3 * d_m2_s);
    constant * a * d_m2_s * re.powf(0.7) * sc.cbrt() * (a * packing.nominal_size_m).powi(-2)
}

// Robbins (1991) works in US customary units.
const LB_FT2_H_PER_KG_M2_S: f64 = 3600.0 / 0.453_592_37 * 0.3048 * 0.3048;
const LB_FT3_PER_KG_M3: f64 = 0.3048 * 0.3048 * 0.3048 / 0.453_592_37;
/// Pa/m per inH2O/ft.
const PA_M_PER_INH2O_FT: f64 = 249.088_908_333 / 0.3048;

fn robbins_inh2o_ft(packing: &PackingData, gas: &FilmProperties, liquid: &FilmProperties) -> f64 {
    let fpd = packing.packing_factor_m * 0.3048;
    let g = gas.mass_flux_kg_m2_s * LB_FT2_H_PER_KG_M2_S;
    let l = liquid.mass_flux_kg_m2_s * LB_FT2_H_PER_KG_M2_S;
    let rho_g = gas.density_kg_m3 * LB_FT3_PER_KG_M3;
    let rho_l = liquid.density_kg_m3 * LB_FT3_PER_KG_M3;
    let mu_cp = liquid.viscosity_pa_s * 1e3;
    let gf = g * (0.075 / rho_g).sqrt() * (fpd / 20.0).sqrt();
    let size = if fpd > 15.0 {
        (fpd / 20.0).sqrt()
    } else {
        (20.0 / fpd).sqrt()
    };
    let lf = l * (62.4 / rho_l) * size * mu_cp.powf(0.1);
    let dry = 7.4e-8 * gf * gf * 10f64.powf(2.7e-5 * lf);
    dry + 0.4 * (lf / 20_000.0).powf(0.1) * dry.powi(4)
}

/// Robbins (1991) GPDC pressure drop per height of packing, Pa/m.
pub fn robbins_pressure_drop_pa_m(
    packing: &PackingData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
) -> f64 {
    robbins_inh2o_ft(packing, gas, liquid) * PA_M_PER_INH2O_FT
}

/// Kister–Gill flood pressure drop, Pa/m.
pub fn flood_pressure_drop_pa_m(packing: &PackingData) -> f64 {
    0.115 * (packing.packing_factor_m * 0.3048).powf(0.7) * PA_M_PER_INH2O_FT
}

/// Fraction of flooding: the gas mass flux divided by the gas mass flux at
/// which the Robbins pressure drop equals the Kister–Gill flood pressure
/// drop at the same liquid load and phase properties.
pub fn flooding_fraction(
    packing: &PackingData,
    gas: &FilmProperties,
    liquid: &FilmProperties,
) -> f64 {
    let target = flood_pressure_drop_pa_m(packing);
    let at = |flux: f64| {
        let trial = FilmProperties {
            mass_flux_kg_m2_s: flux,
            ..*gas
        };
        robbins_pressure_drop_pa_m(packing, &trial, liquid)
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
    gas.mass_flux_kg_m2_s / (0.5 * (lo + hi))
}

/// Liquid wetting rate u_L / a, m³/(m s), for superficial liquid velocity
/// `u_l_m_s`.
pub fn wetting_rate_m2_s(packing: &PackingData, u_l_m_s: f64) -> f64 {
    u_l_m_s / packing.specific_area_m2_m3
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
    fn library_entries_are_physical_and_cited() {
        for packing in Packing::ALL {
            let d = packing.data();
            assert!(d.specific_area_m2_m3 > 0.0 && d.void_fraction < 1.0);
            assert!(d.packing_factor_m > 0.0 && d.nominal_size_m > 0.015);
            assert!(d.source.contains("RASCHIG"));
        }
    }

    #[test]
    fn onda_area_is_a_fraction_of_the_packing_area() {
        let d = Packing::MetalPallRing25.data();
        let aw = onda_wetted_area(&d, &water(), 0.072);
        assert!(aw > 0.0 && aw < d.specific_area_m2_m3);
        // Lower surface tension relative to σc wets more of the packing.
        assert!(onda_wetted_area(&d, &water(), 0.03) > aw);
        let kl = onda_liquid_coefficient(&d, &water(), aw, 2e-9);
        let kg = onda_gas_coefficient(&d, &air(1.0), 2e-5);
        // Typical magnitudes: k_L ~ 1e-4 m/s, k_G ~ 1e-2 m/s.
        assert!((1e-5..1e-3).contains(&kl), "{kl}");
        assert!((1e-3..1e-1).contains(&kg), "{kg}");
    }

    #[test]
    fn robbins_dry_pressure_drop_and_flooding() {
        let d = Packing::MetalPallRing50.data();
        let dry = FilmProperties {
            mass_flux_kg_m2_s: 0.0,
            ..water()
        };
        // Dry bed: ΔP = 7.4e-8 Gf² inH2O/ft.
        let g = 1.0 * LB_FT2_H_PER_KG_M2_S * (0.075 / (1.2 * LB_FT3_PER_KG_M3)).sqrt();
        let gf2 = g * g * (66.0 * 0.3048 / 20.0);
        let expected = 7.4e-8 * gf2 * (1.0 + 0.0) * PA_M_PER_INH2O_FT;
        let got = robbins_pressure_drop_pa_m(&d, &air(1.0), &dry);
        assert!((got / expected - 1.0).abs() < 1e-9, "{got} {expected}");
        let wet = robbins_pressure_drop_pa_m(&d, &air(1.0), &water());
        assert!(wet > got);
        let f = flooding_fraction(&d, &air(1.0), &water());
        assert!(f > 0.0 && f < 1.0, "{f}");
        let flood_flux = 1.0 / f;
        let at_flood = robbins_pressure_drop_pa_m(&d, &air(flood_flux), &water());
        assert!((at_flood / flood_pressure_drop_pa_m(&d) - 1.0).abs() < 1e-9);
    }
}
