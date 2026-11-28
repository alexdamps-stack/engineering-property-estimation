"""Independent rate-based packed absorber references (#241). Never calls Rust.

The documented segment model of docs/a17-columns-rate-based.md (mixed-flow
segments, two-film effective-diffusivity model with interface equilibrium,
Chilton-Colburn heat transfer in both phases, Onda 1968 coefficients,
Robbins 1991 pressure drop and Kister-Gill flood pressure drop as given in
Perry's 8th ed., Glitsch minimum wetting) is solved here in a separate
formulation:

* unknowns per segment are T_V, T_L, T_I, bulk x and y, interface x_I and
  y_I, total flows L and V and the component fluxes N, closed with the
  summation equations (the engine uses log component flows and log ratios);
* partial molar enthalpies are H_i = h_i^ig(T) - R T^2 (d ln phi_i/dT)_{p,x}
  from teqp fugacity coefficients by central differences in T and Cantera
  NASA7 ideal-gas enthalpies (the engine differentiates the phase enthalpy
  numerically in mole numbers);
* Newton's method on the whole column with a dense forward-difference
  Jacobian (relative steps) and a sum-of-squares line search, started from
  the feed compositions (with traces) at the flow-weighted inlet
  temperature and continued in the
  interfacial-area factor from 1e3 times the specified value.

Shared with the engine (not independent): the component constants, NASA7
data, the transport and liquid-density extracts, the PR/PR-BM formulation
(teqp adapters of validation/generate-references.py and
generate-flash-references.py), the transport methods (evaluated by
validation/transport/generate-transport-references.py), the engineering
liquid density (validation/generate-liquid-density-references.py, chemicals
1.5.2), the packing data and the correlation set itself. The fixture
therefore checks implementation agreement of the segment model, not the
accuracy of the correlations against experiment. The correlation functions
transcribed here are checked against published worked examples first
(`published_examples`), recorded in the fixture.

Usage (Python 3.14 with teqp 0.23.2, cantera 3.2.0, numpy 2.3.4,
scipy 1.16.3, CoolProp 7.2.0, chemicals 1.5.2; on Windows in WSL2):

    python -I validation/columns/generate-rate-based-references.py --output <candidate.json>

Review the difference before replacing rate-based-v1.json; never regenerate
to make a failing test pass.
"""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

PATH = Path(__file__).resolve()
ROOT = PATH.parents[2]
R = 8.31446261815324
G = 9.80665


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / "validation" / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


caloric = load("caloric", "generate-references.py")
equilibrium = load("equilibrium", "generate-flash-references.py")
transport = load("transport_reference", "transport/generate-transport-references.py")
density = load("density_reference", "generate-liquid-density-references.py")
TRANSPORT, BASE = transport.load()
DENSITY = density.load()

# Budgets fixed before any Rust comparison ([absolute, relative]; a value
# passes when |rust - ref| <= abs + rel*|ref|). Both sides converge their
# equations to ~1e-10 scaled residuals; the remaining difference is the
# partial molar enthalpy evaluation (engine: forward difference with a 1e-6
# mole-number step, truncation ~5e-3 J/mol; here central differences of
# ln phi in T, ~1e-6 J/mol). Propagated through fluxes of <~1 mol/s and
# liquid heat capacities of ~1.5 kW/K this is <~1e-5 K; budgets carry a
# factor of 10 above that.
TOLERANCES = {
    "temperature_K": [1e-4, 0.0],
    "composition": [1e-9, 1e-5],
    "molar_flow_mol_s": [1e-7, 1e-7],
    "flux_mol_s": [1e-8, 1e-4],
    "hydraulics": [0.0, 1e-5],
}

# Perry's Chemical Engineers' Handbook, 8th ed., Table 14-13 (p. 14-60),
# metal Pall rings: area m2/m3, voids, Fp (Kister) and Fpd (Robbins), 1/m.
PACKINGS = {
    "metal_pall_ring_25": dict(size=0.025, a=205.0, eps=0.94, fp=183.0, fpd=174.0),
    "metal_pall_ring_50": dict(size=0.050, a=105.0, eps=0.96, fp=89.0, fpd=79.0),
}
# Fuller, Schettler and Giddings (1966) diffusion volumes (Cussler,
# Diffusion, Table 5.1-4): molecules and atomic increments C, H, O.
FULLER = {"nitrogen": 17.9, "oxygen": 16.6, "carbon-dioxide": 26.9, "water": 12.7,
          "methane": 16.5 + 4 * 1.98, "ethane": 2 * 16.5 + 6 * 1.98,
          "propane": 3 * 16.5 + 8 * 1.98, "methanol": 16.5 + 4 * 1.98 + 5.48}
# Wilke and Chang (1955) association factors.
PHI = {"water": 2.6, "methanol": 1.9}
# Perry's 8th ed. Table 14-16 (Glitsch), bright metal, m3/(m2 h).
GLITSCH_BRIGHT_METAL = 3.0
MAX_NEGLECTED = 0.01
ATM = 101325.0


def inlet(t, p, flow, z):
    return dict(temperature_K=t, pressure_Pa=p, molar_flow_mol_s=flow, z=z)


IDS = ["nitrogen", "methanol", "water"]
# critical_surface_tension_N_m is a case input of the Onda wetted area (the
# engine has no sigma_c library value; 0.075 N/m is the value commonly quoted
# for steel and is not verified here).
CASES = [
    dict(id="methanol-scrubbing-pall50-PR", variant="PR", pressure_Pa=2.5e5,
         packed_height_m=2.0, diameter_m=0.4, segments=6, packing="metal_pall_ring_50",
         critical_surface_tension_N_m=0.075, area_factor=1.0,
         gas_neglected=["methanol", "water"], liquid_neglected=["nitrogen"],
         gas=inlet(298.0, 2.6e5, 15.0, [0.995, 0.005, 0.0]),
         liquid=inlet(288.0, 2.7e5, 30.0, [0.0, 0.0, 1.0])),
    dict(id="methanol-scrubbing-pall25-PR-BM", variant="PR-BM", pressure_Pa=2e5,
         packed_height_m=1.5, diameter_m=0.3, segments=5, packing="metal_pall_ring_25",
         critical_surface_tension_N_m=0.075, area_factor=1.0,
         gas_neglected=["methanol", "water"], liquid_neglected=["nitrogen"],
         gas=inlet(295.0, 2.1e5, 6.0, [0.996, 0.004, 0.0]),
         liquid=inlet(288.0, 2.2e5, 25.0, [0.0, 0.0, 1.0])),
]


def zero_kij():
    return np.zeros((len(IDS), len(IDS))).tolist()


def root_state(x, t, p, variant, phase):
    case = caloric.case(IDS, [float(v) for v in x], t, p, variant, zero_kij())
    stable = [r for r in case["roots"] if r["mechanically_stable"]]
    r = stable[0] if phase == "liquid" else stable[-1]
    return dict(h=r["h"], cp=r["cp"], rho=r["mass_density"], v=r["Z"] * R * t / p,
                ln_phi=np.array(r["ln_phi"]), roots=len(stable))


def partial_h(x, t, p, variant, phase):
    dt = 1e-3
    hi = root_state(x, t + dt, p, variant, phase)["ln_phi"]
    lo = root_state(x, t - dt, p, variant, phase)["ln_phi"]
    ideal = np.array([caloric.nasa(caloric.DATA[i], t)[0] for i in IDS])
    return ideal - R * t * t * (hi - lo) / (2 * dt)


def basis(x, neglected, strict=True):
    w = np.array(x, dtype=float)
    dropped = 0.0
    for name in neglected:
        i = IDS.index(name)
        assert not strict or w[i] <= MAX_NEGLECTED, (name, w[i])
        dropped = max(dropped, w[i])
        w[i] = 0.0
    return w / w.sum(), dropped


def molar_mass(i):
    return BASE[IDS[i]]["constants"]["molar_mass_kg_per_mol"] * 1e3


def fuller_cm2_s(t, p_atm, m_a, m_b, v_a, v_b):
    """Fuller, Schettler and Giddings (1966), cm2/s, P in atm."""
    return 1e-3 * t ** 1.75 * math.sqrt(1 / m_a + 1 / m_b) / (p_atm * (v_a ** (1 / 3) + v_b ** (1 / 3)) ** 2)


def fuller(i, j, t, p):
    return fuller_cm2_s(t, p / ATM, molar_mass(i), molar_mass(j), FULLER[IDS[i]], FULLER[IDS[j]]) * 1e-4


def wilke_chang_cm2_s(t, phi, m_solvent, mu_cp, v_solute):
    return 7.4e-8 * math.sqrt(phi * m_solvent) * t / (mu_cp * v_solute ** 0.6)


def gas_diffusivities(y, t, p):
    out = []
    for i in range(len(IDS)):
        resistance = sum(y[j] / fuller(i, j, t, p) for j in range(len(IDS)) if j != i)
        out.append((sum(y) - y[i]) / resistance)
    return out


def liquid_diffusivities(solvent, t, mu):
    """Wilke-Chang, each component dilute in the solvent; Tyn-Calus V_b."""
    out = []
    for i in range(len(IDS)):
        if i == solvent:
            out.append(0.0)
            continue
        vb = 0.285 * (TRANSPORT[IDS[i]]["critical_volume_m3_per_mol"] * 1e6) ** 1.048
        out.append(wilke_chang_cm2_s(t, PHI.get(IDS[solvent], 1.0), molar_mass(solvent),
                                     mu * 1e3, vb) * 1e-4)
    return out


def onda(pk, sigma_c, lflux, rho_l, mu_l, sigma, gflux, rho_g, mu_g):
    a, dp = pk["a"], pk["size"]
    aw = a * (1 - math.exp(-1.45 * (sigma_c / sigma) ** 0.75 * (lflux / (a * mu_l)) ** 0.1
                           * (lflux ** 2 * a / (rho_l ** 2 * G)) ** -0.05
                           * (lflux ** 2 / (rho_l * sigma * a)) ** 0.2))

    def k_liquid(d):
        return (0.0051 * (lflux / (aw * mu_l)) ** (2 / 3) * (mu_l / (rho_l * d)) ** -0.5
                * (a * dp) ** 0.4 / (rho_l / (mu_l * G)) ** (1 / 3))

    def k_gas(d):
        c = 2.0 if dp < 0.015 else 5.23
        return c * a * d * (gflux / (a * mu_g)) ** 0.7 * (mu_g / (rho_g * d)) ** (1 / 3) * (a * dp) ** -2

    return aw, k_liquid, k_gas


LB = 0.45359237
FT = 0.3048
PA_PER_INH2O_FT = 249.088908333 / FT


def robbins_factors(fpd_ft, g_lb, rho_g_lb, l_lb, rho_l_lb, mu_cp, p):
    """Robbins Gf, Lf (lb/ft2 h), Perry's 8th ed. eqs. 14-148 to 14-151."""
    fs = g_lb / 3600 / math.sqrt(rho_g_lb)
    gf = 986 * fs * math.sqrt(fpd_ft / 20)
    if p > ATM:
        gf *= 10 ** (0.3 * rho_g_lb)
    if fpd_ft > 200:
        lf = l_lb * (62.4 / rho_l_lb) * math.sqrt(fpd_ft / 20) * mu_cp ** 0.2
    elif fpd_ft < 15:
        lf = l_lb * (62.4 / rho_l_lb) * math.sqrt(20 / fpd_ft) * mu_cp ** 0.1
    else:
        lf = l_lb * (62.4 / rho_l_lb) * math.sqrt(fpd_ft / 20) * mu_cp ** 0.1
    return gf, lf


def robbins_inh2o_ft(gf, lf):
    """Perry's eqs. 14-145 to 14-147, C3 = 7.4e-8, C4 = 2.7e-5."""
    dry = 7.4e-8 * gf * gf * 10 ** (2.7e-5 * lf)
    return dry + 0.4 * (lf / 20000) ** 0.1 * dry ** 4


def robbins(pk, gflux, rho_g, lflux, rho_l, mu_l, p):
    """Robbins (1991) pressure drop, inH2O/ft, from SI inputs."""
    assert p <= 3e5, "Robbins is not suitable above 3 bar absolute"
    lb_flux = 3600 / LB * FT * FT
    gf, lf = robbins_factors(pk["fpd"] * FT, gflux * lb_flux, rho_g * FT ** 3 / LB,
                             lflux * lb_flux, rho_l * FT ** 3 / LB, mu_l * 1e3, p)
    return robbins_inh2o_ft(gf, lf)


def hydraulics(pk, gflux, rho_g, lflux, rho_l, mu_l, p):
    drop = robbins(pk, gflux, rho_g, lflux, rho_l, mu_l, p)
    fp_ft = pk["fp"] * FT
    assert fp_ft <= 60, "Kister-Gill flood equation not for Fp > 60 1/ft"
    flood = 0.12 * fp_ft ** 0.7
    g_flood = brentq(lambda g: robbins(pk, g, rho_g, lflux, rho_l, mu_l, p) - flood,
                     1e-9, 1e3, xtol=1e-15, rtol=1e-15)
    return drop * PA_PER_INH2O_FT, gflux / g_flood


def minimum_wetting_m_h(pk, sigma):
    """Perry's eq. 14-166 with Table 14-16, doubled above 25 mN/m."""
    base = GLITSCH_BRIGHT_METAL * math.sqrt(200 / pk["a"])
    return 2 * base if sigma > 0.025 else base


def published_examples():
    """Correlation functions above against published worked examples."""
    out = []
    # Cussler, Diffusion, 3rd ed., Example 5.1-2: H2-N2, 294 K, 2 atm -> 0.37 cm2/s.
    d = fuller_cm2_s(294.0, 2.0, 2.02, 28.0, 7.07, 17.9)
    assert abs(d - 0.37) < 0.005, d
    out.append(dict(example="Cussler Example 5.1-2 (Fuller, H2-N2)", published=0.37, computed=d,
                    unit="cm2/s"))
    # n-butanol in air, 273 K, 1 atm -> 7.73e-6 m2/s (U. Jordan notes).
    d = fuller_cm2_s(273.0, 1.0, 74.1, 29.0, 4 * 16.5 + 10 * 1.98 + 5.48, 20.1) * 1e-4
    assert abs(d / 7.73e-6 - 1) < 0.005, d
    out.append(dict(example="Fuller n-butanol-air (U. Jordan notes)", published=7.73e-6,
                    computed=d, unit="m2/s"))
    # Wilke-Chang ethylbenzene in water, 293 K, 1.0 cP, V 139.5 -> 0.77e-5 cm2/s.
    d = wilke_chang_cm2_s(293.0, 2.6, 18.0, 1.0, 139.5)
    assert abs(d - 0.77e-5) < 0.005e-5, d
    out.append(dict(example="Wilke-Chang ethylbenzene-water (U. Jordan notes)", published=0.77e-5,
                    computed=d, unit="cm2/s"))
    # Perry's 8th ed. Example 13: Gf = 1653, Lf = 9859, dP = 0.40 inH2O/ft (chart).
    gf, lf = robbins_factors(24.0, 1500.0, 0.074, 9000.0, 62.4, 1.0, ATM)
    dp = robbins_inh2o_ft(gf, lf)
    # The example rounds Fs to 1.53 (exact 1.5317), so Gf is 1654.4 here.
    assert abs(gf - 1653) < 2 and abs(lf - 9859) < 2 and abs(dp - 0.40) < 0.03, (gf, lf, dp)
    out.append(dict(example="Perry's 8th ed. Example 13 (Robbins)", published=[1653, 9859, 0.40],
                    computed=[gf, lf, dp], unit="lb/(ft2 h), lb/(ft2 h), inH2O/ft"))
    return out


class Column:
    def __init__(self, case):
        self.case = case
        # The 1 mol% limit of neglected components is checked on the solution,
        # not on Newton iterates (which may pass through larger fractions).
        self.strict = False
        self.n = case["segments"]
        self.c = len(IDS)
        self.p = case["pressure_Pa"]
        self.variant = case["variant"]
        self.pk = PACKINGS[case["packing"]]
        self.area = math.pi * case["diameter_m"] ** 2 / 4
        self.dz = case["packed_height_m"] / self.n
        self.feeds = {}
        for port, phase in (("gas", "vapor"), ("liquid", "liquid")):
            f = case[port]
            flash = equilibrium.flash_case(port, IDS, f["z"], f["temperature_K"], f["pressure_Pa"],
                                           self.variant)
            assert flash["phases"] == 1, (port, flash)
            h = root_state(f["z"], f["temperature_K"], f["pressure_Pa"], self.variant, phase)["h"]
            self.feeds[port] = (f["molar_flow_mol_s"], np.array(f["z"]), h)
        self.flow = self.feeds["gas"][0] + self.feeds["liquid"][0]
        self.ref_v = int(np.argmax(self.feeds["gas"][1]))
        self.ref_l = int(np.argmax(self.feeds["liquid"][1]))

    def size(self):
        return 5 * self.c + 5

    def unpack(self, values):
        c = self.c
        out = []
        for b in values.reshape(self.n, self.size()):
            out.append(dict(tv=b[0] * 100, tl=b[1] * 100, ti=b[2] * 100,
                            x=b[3:3 + c], y=b[3 + c:3 + 2 * c], xi=b[3 + 2 * c:3 + 3 * c],
                            yi=b[3 + 3 * c:3 + 4 * c], L=b[3 + 4 * c] * self.flow,
                            V=b[4 + 4 * c] * self.flow, N=b[5 + 4 * c:5 + 5 * c] * self.flow))
        return out

    def segment_props(self, s, factor):
        p, v = self.p, self.variant
        bv = root_state(s["y"], s["tv"], p, v, "vapor")
        bl = root_state(s["x"], s["tl"], p, v, "liquid")
        iv = root_state(s["yi"], s["ti"], p, v, "vapor")
        il = root_state(s["xi"], s["ti"], p, v, "liquid")
        gas_x, gas_drop = basis(s["y"], self.case["gas_neglected"], self.strict)
        liq_x, liq_drop = basis(s["x"], self.case["liquid_neglected"], self.strict)
        gv = root_state(gas_x, s["tv"], p, v, "vapor")["v"]
        gas_t, _ = transport.gas(IDS, list(gas_x), s["tv"], gv, TRANSPORT, BASE)
        liq_t = transport.liquid(IDS, list(liq_x), s["tl"], p, TRANSPORT, BASE)
        mu_v, lam_v = gas_t["viscosity"], gas_t["thermal_conductivity"]
        mu_l, lam_l, sigma = (liq_t["viscosity"], liq_t["thermal_conductivity"],
                              liq_t["surface_tension"])
        assert None not in (mu_v, lam_v, mu_l, lam_l, sigma)
        vol = density.expected_volume(DENSITY, IDS, list(s["x"]), s["tl"], p)[1]
        m_l = sum(xi * molar_mass(i) for i, xi in enumerate(s["x"])) * 1e-3
        rho_l = m_l / vol
        m_v = bv["rho"] * bv["v"]
        gflux = s["V"] * m_v / self.area
        lflux = s["L"] * m_l / self.area
        aw, k_liq, k_gas = onda(self.pk, self.case["critical_surface_tension_N_m"], lflux, rho_l,
                                mu_l, sigma, gflux, bv["rho"], mu_v)
        d_v = gas_diffusivities(s["y"], s["tv"], p)
        d_l = liquid_diffusivities(self.ref_l, s["tl"], mu_l)
        a_int = factor * aw * self.area * self.dz
        others_v = [i for i in range(self.c) if i != self.ref_v]
        others_l = [i for i in range(self.c) if i != self.ref_l]
        dv = sum(d_v[i] for i in others_v) / len(others_v)
        dl = sum(d_l[i] for i in others_l) / len(others_l)
        c_v, c_l = 1 / bv["v"], 1 / vol
        le_v = lam_v / (c_v * bv["cp"] * dv)
        le_l = lam_l / (c_l * bl["cp"] * dl)
        return dict(
            bv=bv, bl=bl, ln_k=il["ln_phi"] - iv["ln_phi"],
            hbar_v=partial_h(s["yi"], s["ti"], p, v, "vapor"),
            hbar_l=partial_h(s["xi"], s["ti"], p, v, "liquid"),
            c_v=c_v, c_l=c_l, a_int=a_int, area_m2_m3=factor * aw,
            kv=[k_gas(d_v[i]) for i in range(self.c)],
            kl=[k_liq(d_l[i]) for i in range(self.c)],
            ha_v=k_gas(dv) * c_v * bv["cp"] * le_v ** (2 / 3) * a_int,
            ha_l=k_liq(dl) * c_l * bl["cp"] * le_l ** (2 / 3) * a_int,
            gflux=gflux, lflux=lflux, rho_l=rho_l, mu_l=mu_l, u_l=lflux / rho_l, sigma=sigma,
            neglected=max(gas_drop, liq_drop))

    def residual(self, values, factor):
        segs = self.unpack(values)
        props = [self.segment_props(s, factor) for s in segs]
        f, scale = self.flow, self.flow * 1e4
        out = []
        for j, (s, pr) in enumerate(zip(segs, props)):
            if j > 0:
                up = segs[j - 1]
                l_in, hl_in = up["L"] * up["x"], up["L"] * props[j - 1]["bl"]["h"]
            else:
                fl, z, h = self.feeds["liquid"]
                l_in, hl_in = fl * z, fl * h
            if j + 1 < self.n:
                dn = segs[j + 1]
                v_in, hv_in = dn["V"] * dn["y"], dn["V"] * props[j + 1]["bv"]["h"]
            else:
                fg, z, h = self.feeds["gas"]
                v_in, hv_in = fg * z, fg * h
            n_t = s["N"].sum()
            out += list((v_in - s["V"] * s["y"] - s["N"]) / f)
            out += list((l_in - s["L"] * s["x"] + s["N"]) / f)
            out += [s["x"].sum() - 1, s["y"].sum() - 1]
            e_v = pr["ha_v"] * (s["tv"] - s["ti"]) + float(np.dot(s["N"], pr["hbar_v"]))
            e_l = pr["ha_l"] * (s["ti"] - s["tl"]) + float(np.dot(s["N"], pr["hbar_l"]))
            # Equations scaled by their film conductances as well (row scaling
            # only; it keeps the large-area continuation steps well scaled).
            sv, sl = scale + 100 * pr["ha_v"], scale + 100 * pr["ha_l"]
            out += [(hv_in - s["V"] * pr["bv"]["h"] - e_v) / sv,
                    (hl_in - s["L"] * pr["bl"]["h"] + e_l) / sl, (e_v - e_l) / (sv + sl - scale)]
            out += list(np.log(s["yi"]) - np.log(s["xi"]) - pr["ln_k"])
            out += [s["xi"].sum() - 1, s["yi"].sum() - 1]
            for i in range(self.c):
                if i != self.ref_v:
                    mean = 0.5 * (s["y"][i] + s["yi"][i])
                    g_v = pr["c_v"] * pr["kv"][i] * pr["a_int"]
                    out.append((s["N"][i] - g_v * (s["y"][i] - s["yi"][i]) - mean * n_t) / (f + g_v))
                if i != self.ref_l:
                    mean = 0.5 * (s["x"][i] + s["xi"][i])
                    g_l = pr["c_l"] * pr["kl"][i] * pr["a_int"]
                    out.append((s["N"][i] - g_l * (s["xi"][i] - s["x"][i]) - mean * n_t) / (f + g_l))
        return np.array(out), segs, props

    def initial(self):
        fg, zg, _ = self.feeds["gas"]
        fl, zl, _ = self.feeds["liquid"]
        z = (fg * zg + fl * zl) / self.flow
        t = (fg * self.case["gas"]["temperature_K"] + fl * self.case["liquid"]["temperature_K"]) / self.flow
        # Starting point only: each phase at its feed composition with traces
        # of the absent components, at the flow-weighted inlet temperature.
        del z
        x = (zl + 1e-6) / (zl + 1e-6).sum()
        y = (zg + 1e-6) / (zg + 1e-6).sum()
        beta = fg / self.flow
        block = np.r_[t / 100, t / 100, t / 100, x, y, x, y, (1 - beta), beta,
                      (fg * zg - beta * self.flow * y) / self.flow / self.n]
        return np.tile(block, self.n)

    def newton(self, values, factor):
        r, _, _ = self.residual(values, factor)
        for iteration in range(80):
            err = np.max(np.abs(r))
            print(f"{self.case['id']} factor {factor:g} iteration {iteration}: {err:.3e}",
                  file=sys.stderr, flush=True)
            # The rounding floor of the scaled residual is ~2e-11 (density
            # refinement and the temperature differences of ln phi).
            if err < 5e-11:
                return values
            jac = np.empty((len(r), len(values)))
            for k in range(len(values)):
                # Relative steps: trace mole fractions are ~1e-7.
                h = 1e-7 * max(abs(values[k]), 1e-9)
                trial = values.copy()
                trial[k] += h
                jac[:, k] = (self.residual(trial, factor)[0] - r) / h
            step = np.linalg.solve(jac, -r)
            # Backtracking on the sum of squares (the engine uses the max norm).
            merit = float(np.dot(r, r))
            for back in range(30):
                trial = values + 0.5 ** back * step
                try:
                    rt = self.residual(trial, factor)[0]
                except (AssertionError, ValueError, ZeroDivisionError, TypeError):
                    continue
                if np.all(np.isfinite(rt)) and float(np.dot(rt, rt)) < merit:
                    values, r = trial, rt
                    break
            else:
                # Below 1e-10 the scaled residual is at the rounding floor of
                # the temperature differences in the partial molar enthalpies.
                if err < 1e-10:
                    return values
                raise RuntimeError(f"line search stalled at {err:.3e}")
        raise RuntimeError("Newton did not converge")

    def solve(self):
        values = self.initial()
        target = self.case["area_factor"]
        # Newton along the area-factor continuation from the blend-flash guess.
        for multiple in (1e3, 1e2, 1e1, 1.0):
            values = self.newton(values, target * multiple)
        return values


def reference(case):
    col = Column(case)
    values = col.solve()
    col.strict = True
    r, segs, props = col.residual(values, case["area_factor"])
    assert np.max(np.abs(r)) < 1e-10
    profiles = dict(vapor_temperatures_K=[], liquid_temperatures_K=[], interface_temperatures_K=[],
                    vapor_flows_mol_s=[], liquid_flows_mol_s=[], vapor_z=[], liquid_z=[], fluxes_mol_s=[],
                    interfacial_area_m2_m3=[], flooding_fractions=[], pressure_drops_Pa_m=[])
    neglected = 0.0
    for s, pr in zip(segs, props):
        drop, flood = hydraulics(col.pk, pr["gflux"], pr["bv"]["rho"], pr["lflux"], pr["rho_l"],
                                 pr["mu_l"], col.p)
        assert flood < 1 and pr["u_l"] * 3600 >= minimum_wetting_m_h(col.pk, pr["sigma"])
        # The interface phases satisfy equal fugacities (residual above) and
        # are each stable by the A4 multi-start TPD search, so the interface
        # is a stable vapor-liquid split.
        model, _ = equilibrium.model_at(IDS, s["ti"], col.variant, zero_kij())
        for phase in (s["xi"], s["yi"]):
            tm = equilibrium.tpd_search(model, IDS, s["ti"], col.p, np.array(phase))[0]
            assert tm >= equilibrium.TOLERANCES["stable_tm_floor"], (case["id"], tm)
        neglected = max(neglected, pr["neglected"])
        for key, value in [("vapor_temperatures_K", s["tv"]), ("liquid_temperatures_K", s["tl"]),
                           ("interface_temperatures_K", s["ti"]), ("vapor_flows_mol_s", s["V"]),
                           ("liquid_flows_mol_s", s["L"]), ("vapor_z", list(s["y"])),
                           ("liquid_z", list(s["x"])), ("fluxes_mol_s", list(s["N"])),
                           ("interfacial_area_m2_m3", pr["area_m2_m3"]), ("flooding_fractions", flood),
                           ("pressure_drops_Pa_m", drop)]:
            profiles[key].append(value)
    top, bottom = segs[0], segs[-1]
    to_float = lambda v: [float(x) for x in v] if isinstance(v, list) else float(v)
    return dict(case, ids=IDS,
                profiles={k: [to_float(v) for v in vs] for k, vs in profiles.items()},
                pressure_drop_Pa=float(sum(profiles["pressure_drops_Pa_m"]) * col.dz),
                max_flooding_fraction=float(max(profiles["flooding_fractions"])),
                max_neglected_transport_fraction=float(neglected),
                outlets=[dict(molar_flow_mol_s=float(top["V"]), temperature_K=float(top["tv"]),
                              z=to_float(list(top["y"]))),
                         dict(molar_flow_mol_s=float(bottom["L"]), temperature_K=float(bottom["tl"]),
                              z=to_float(list(bottom["x"])))])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    watched = [PATH, ROOT / "validation/generate-references.py",
               ROOT / "validation/generate-flash-references.py",
               ROOT / "validation/transport/generate-transport-references.py",
               ROOT / "validation/generate-liquid-density-references.py",
               ROOT / "data/source-extracts/components.provisional.json",
               ROOT / "data/source-extracts/transport.provisional.json"]
    result = dict(
        schema="a17-rate-based-absorber-reference-v1",
        published_examples=published_examples(),
        packages={name: importlib.metadata.version(name)
                  for name in ("teqp", "cantera", "numpy", "scipy", "CoolProp", "chemicals")},
        sha256={str(path.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        tolerances=TOLERANCES, cases=[reference(case) for case in CASES],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
