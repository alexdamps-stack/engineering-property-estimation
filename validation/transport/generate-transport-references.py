"""Independent transport-property references for #225. Never calls Rust.

Two separate reference sets:

1. ``implementation`` - the documented formulations (Chung et al. 1988 dense
   gas method with Chung mixing rules; ChemSep saturated-liquid correlations
   with the labelled mixing rules) evaluated by this separate Python
   implementation from the committed source extracts, at prescribed T, p,
   molar volume and composition. Pure-component Chung conductivities are
   cross-checked at generation time against the independent ``chemicals``
   implementation (``Chung_dense``). Expected-unavailable states follow the
   declared validity domain of docs/transport-properties.md.
2. ``coolprop`` - model comparison of pure components against CoolProp 7.2.0
   reference transport correlations (REFPROP-quality fluid-specific models),
   at states inside the declared domain, away from saturation. Rust evaluates
   these with its own Peng-Robinson gas volume, so the budgets include the
   method and EOS-density differences; they are declared method accuracies,
   not implementation budgets.

Usage (pinned environment: Python 3.14, CoolProp 7.2.0, chemicals 1.5.2):

    python validation/transport/generate-transport-references.py --output <candidate.json>

Compare candidates with the committed fixture before adopting them; never
regenerate to make a failing test pass.
"""
import argparse
import hashlib
import json
import math
import pathlib
import sys

import CoolProp
import CoolProp.CoolProp as CP
import chemicals
from chemicals.thermal_conductivity import Chung_dense

ROOT = pathlib.Path(__file__).resolve().parents[2]
TRANSPORT = ROOT / "data/source-extracts/transport.provisional.json"
COMPONENTS = ROOT / "data/source-extracts/components.provisional.json"
R = 8.31446261815324

# Declared domain (docs/transport-properties.md); mirrors the engine constants.
GAS_T = (200.0, 600.0)
GAS_MAX_RHO_R = 1.5
POLAR_GAS_MAX_RHO_R = 0.1
POLAR_DIPOLE = 0.5
LIQUID_MAX_TR = 0.75
LIQUID_MAX_P = 5.0e6

# Budgets, declared per quantity: abs(actual-expected) <= atol + rtol*abs(expected).
TOLERANCES = {
    "implementation_viscosity": [1e-15, 1e-10],
    "implementation_thermal_conductivity": [1e-13, 1e-10],
    "implementation_surface_tension": [1e-15, 1e-10],
    "coolprop_dilute_gas_viscosity_nonpolar": [0.0, 0.06],
    "coolprop_dense_fluid_viscosity_nonpolar": [0.0, 0.12],
    "coolprop_gas_viscosity_polar": [0.0, 0.15],
    "coolprop_gas_thermal_conductivity_nonpolar": [0.0, 0.12],
    "coolprop_liquid_viscosity": [0.0, 0.15],
    "coolprop_liquid_thermal_conductivity": [0.0, 0.10],
    "coolprop_surface_tension": [0.0, 0.08],
}

COOLPROP_NAMES = {
    "methane": "Methane", "ethane": "Ethane", "propane": "n-Propane", "nitrogen": "Nitrogen",
    "carbon-dioxide": "CarbonDioxide", "oxygen": "Oxygen", "water": "Water", "methanol": "Methanol",
}

E_VISC = [(6.324, 50.412, -51.680, 1189.0), (1.210e-3, -1.154e-3, -6.257e-3, 0.03728),
          (5.283, 254.209, -168.48, 3898.0), (6.623, 38.096, -8.464, 31.42),
          (19.745, 7.630, -14.354, 31.53), (-1.900, -12.537, 4.985, -18.15),
          (24.275, 3.450, -11.291, 69.35), (0.7972, 1.117, 0.01235, -4.117),
          (-0.2382, 0.06770, -0.8163, 4.025), (0.06863, 0.3479, 0.5926, -0.727)]
B_COND = [(2.4166, 0.74824, -0.91858, 121.72), (-0.50924, -1.5094, -49.991, 69.983),
          (6.6107, 5.6207, 64.760, 27.039), (14.543, -8.9139, -5.6379, 74.344),
          (0.79274, 0.82019, -0.69369, 6.3173), (-5.8634, 12.801, 9.5893, 65.529),
          (91.089, 128.11, -54.217, 523.81)]


def load():
    transport = {c["id"]: c for c in json.loads(TRANSPORT.read_text(encoding="utf-8"))["components"]}
    base = {c["id"]: c for c in json.loads(COMPONENTS.read_text(encoding="utf-8"))["components"]}
    return transport, base


def cv_ideal(rec, t):
    lo, join, hi = rec["ideal_gas"]["temperature_breaks_K"]
    assert lo <= t <= hi
    a = rec["ideal_gas"]["coefficients"][0 if t < join else 1]
    return R * (a[0] + a[1] * t + a[2] * t**2 + a[3] * t**3 + a[4] * t**4) - R


def neufeld(ts):
    return 1.16145 * ts**-0.14874 + 0.52487 * math.exp(-0.77320 * ts) + 2.16178 * math.exp(-2.43787 * ts)


class Mixture:
    """Chung mixing rules, Poling et al. (2001) eqs. 9-5.25 to 9-5.44."""

    def __init__(self, ids, x, transport, base):
        n = len(ids)
        par = []
        for cid in ids:
            c, b = transport[cid], base[cid]["constants"]
            vc = c["critical_volume_m3_per_mol"] * 1e6
            par.append((0.809 * vc ** (1 / 3), b["critical_temperature_K"] / 1.2593, b["acentric_factor"],
                        b["molar_mass_kg_per_mol"] * 1e3, c["dipole_moment_debye"], c["association_factor"]))
        s3 = e = w = msum = mu4 = kap = 0.0
        for i in range(n):
            for j in range(n):
                xx = x[i] * x[j]
                if xx == 0:
                    continue
                a, b = par[i], par[j]
                sij = math.sqrt(a[0] * b[0])
                eij = math.sqrt(a[1] * b[1])
                mij = 2 * a[3] * b[3] / (a[3] + b[3])
                s3 += xx * sij**3
                e += xx * eij * sij**3
                w += xx * (a[2] + b[2]) / 2 * sij**3
                msum += xx * eij * sij**2 * math.sqrt(mij)
                mu4 += xx * a[4] ** 2 * b[4] ** 2 / sij**3
                kap += xx * math.sqrt(a[5] * b[5])
        e /= s3
        self.m = (msum / (e * s3 ** (2 / 3))) ** 2
        self.vc = s3 / 0.809**3
        self.tc = 1.2593 * e
        self.w = w / s3
        self.mur = 131.3 * (s3 * mu4) ** 0.25 / math.sqrt(self.vc * self.tc)
        self.kap = kap
        self.fc = 1 - 0.2756 * self.w + 0.059035 * self.mur**4 + self.kap

    def coef(self, row):
        return row[0] + row[1] * self.w + row[2] * self.mur**4 + row[3] * self.kap


def g2(c, y):
    g1 = (1 - 0.5 * y) / (1 - y) ** 3
    return (c[0] * (-math.expm1(-c[3] * y)) / y + c[1] * g1 * math.exp(c[4] * y) + c[2] * g1) / (c[0] * c[3] + c[1] + c[2])


def chung_mu(m, t, rho_r):
    e = [m.coef(r) for r in E_VISC]
    ts = 1.2593 * t / m.tc
    y = rho_r / 6
    g = g2(e, y)
    ekk = e[6] * y * y * g * math.exp(e[7] + e[8] / ts + e[9] / ts**2)
    ek = math.sqrt(ts) / neufeld(ts) * m.fc * (1 / g + e[5] * y)
    return (ek + ekk) * 36.344 * math.sqrt(m.m * m.tc) / m.vc ** (2 / 3) * 1e-7


def chung_mu0(m, t):
    return 40.785 * m.fc * math.sqrt(m.m * t) / (m.vc ** (2 / 3) * neufeld(1.2593 * t / m.tc)) * 1e-7


def chung_k(m, t, rho_r, cv):
    b = [m.coef(r) for r in B_COND]
    tr = t / m.tc
    alpha = cv / R - 1.5
    beta = 0.7862 - 0.7109 * m.w + 1.3168 * m.w**2
    z = 2 + 10.5 * tr * tr
    psi = 1 + alpha * (0.215 + 0.28288 * alpha - 1.061 * beta + 0.26665 * z) / (0.6366 + beta * z + 1.061 * alpha * beta)
    y = rho_r / 6
    g = g2(b, y)
    mk = m.m / 1000
    q = 3.586e-3 * math.sqrt(m.tc / mk) / m.vc ** (2 / 3)
    return 31.2 * chung_mu0(m, t) * psi / mk * (1 / g + b[5] * y) + q * b[6] * y * y * math.sqrt(tr) * g


def correlation(rec, t):
    lo, hi = rec["temperature_range_K"]
    if not lo <= t <= hi:
        return None
    a, b, c, d, e = rec["coefficients"]
    if rec["equation"] == 101:
        return math.exp(a + b / t + c * math.log(t) + d * t**e)
    assert rec["equation"] == 16
    return a + math.exp(b / t + c + d * t + e * t * t)


def polar(c):
    return c["dipole_moment_debye"] >= POLAR_DIPOLE or c["association_factor"] > 0


def gas(ids, x, t, v, transport, base):
    m = Mixture(ids, x, transport, base)
    rho_r = m.vc * 1e-6 / v
    present = [cid for cid, xi in zip(ids, x) if xi > 0]
    is_polar = any(polar(transport[c]) for c in present)
    ok = GAS_T[0] <= t <= GAS_T[1] and rho_r <= GAS_MAX_RHO_R
    mu = chung_mu(m, t, rho_r) if ok and (not is_polar or rho_r <= POLAR_GAS_MAX_RHO_R) else None
    k = None
    if ok and not is_polar:
        cv = sum(xi * cv_ideal(base[c], t) for c, xi in zip(ids, x) if xi > 0)
        k = chung_k(m, t, rho_r, cv)
        if len(present) == 1:
            c, b = transport[present[0]], base[present[0]]["constants"]
            ref = Chung_dense(t, m.m, m.tc, m.vc * 1e-6, m.w, cv, v, chung_mu0(m, t), m.mur * math.sqrt(m.vc * m.tc) / 131.3, m.kap)
            assert abs(ref / k - 1) < 1e-12, (ids, t, ref, k)
    return {"viscosity": mu, "thermal_conductivity": k, "surface_tension": None}, rho_r


def liquid(ids, x, t, p, transport, base):
    present = [(transport[c], base[c], xi) for c, xi in zip(ids, x) if xi > 0]
    reduced_ok = all(t <= LIQUID_MAX_TR * b["constants"]["critical_temperature_K"] for _, b, _ in present)
    p_ok = reduced_ok and p <= LIQUID_MAX_P
    out = {"viscosity": None, "thermal_conductivity": None, "surface_tension": None}
    if p_ok:
        mus = [correlation(c["liquid_viscosity"], t) for c, _, _ in present]
        if all(v is not None for v in mus):
            out["viscosity"] = math.exp(sum(xi * math.log(mu) for (_, _, xi), mu in zip(present, mus)))
        ks = [correlation(c["liquid_thermal_conductivity"], t) for c, _, _ in present]
        if all(v is not None for v in ks):
            mass = sum(xi * b["constants"]["molar_mass_kg_per_mol"] for _, b, xi in present)
            out["thermal_conductivity"] = sum(
                xi * b["constants"]["molar_mass_kg_per_mol"] / mass / k**2 for (_, b, xi), k in zip(present, ks)) ** -0.5
    if reduced_ok:
        ss = [correlation(c["surface_tension"], t) for c, _, _ in present]
        if all(v is not None for v in ss):
            out["surface_tension"] = sum(xi * s for (_, _, xi), s in zip(present, ss))
    return out


NG = (["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"], [0.85, 0.07, 0.03, 0.03, 0.02])
IMPLEMENTATION_GAS = [
    ("methane-dilute", ["methane"], [1.0], 300.0, 0.0249),
    ("methane-dense", ["methane"], [1.0], 300.0, 1.0e-4),
    ("natural-gas-dense", *NG, 280.0, 2.0e-4),
    ("natural-gas-dilute", *NG, 400.0, 1.0e-2),
    ("air-like", ["nitrogen", "oxygen"], [0.79, 0.21], 250.0, 0.02),
    ("co2-supercritical", ["carbon-dioxide"], [1.0], 320.0, 1.5e-4),
    ("zero-component-retained", ["methane", "water"], [1.0, 0.0], 350.0, 0.01),
    ("wet-methane", ["methane", "water"], [0.95, 0.05], 450.0, 0.03),
    ("water-vapor", ["water"], [1.0], 500.0, 0.04),
    ("methanol-vapor-dilute", ["methanol"], [1.0], 500.0, 2.0e-3),
    ("methanol-vapor-too-dense", ["methanol"], [1.0], 500.0, 5.0e-4),
    ("methane-too-hot", ["methane"], [1.0], 650.0, 0.05),
    ("propane-too-dense", ["propane"], [1.0], 400.0, 8.0e-5),
]
IMPLEMENTATION_LIQUID = [
    ("water-ambient", ["water"], [1.0], 300.0, 1e5),
    ("water-hot", ["water"], [1.0], 450.0, 2e6),
    ("methanol-ambient", ["methanol"], [1.0], 300.0, 1e5),
    ("propane-cold", ["propane"], [1.0], 250.0, 2e6),
    ("ethane-cold", ["ethane"], [1.0], 200.0, 1e6),
    ("methanol-water", ["water", "methanol"], [0.6, 0.4], 320.0, 1e5),
    ("ethane-propane", ["ethane", "propane"], [0.5, 0.5], 210.0, 1e6),
    ("water-compressed", ["water"], [1.0], 300.0, 2e7),
    ("propane-near-critical", ["propane"], [1.0], 300.0, 2e6),
    ("co2-below-correlation-range", ["carbon-dioxide"], [1.0], 210.0, 3e6),
    ("methane-trace-in-propane", ["methane", "propane"], [1e-6, 1 - 1e-6], 200.0, 1e6),
]


def coolprop_cases(transport, base):
    cases = []
    for cid, name in COOLPROP_NAMES.items():
        c = transport[cid]
        tc = base[cid]["constants"]["critical_temperature_K"]
        is_polar = polar(c)
        vc = c["critical_volume_m3_per_mol"]
        for t in (250.0, 350.0, 450.0, 550.0):
            for p in (1e4, 1e5, 2e6, 1e7):
                if t < CP.PropsSI("Tmin", name):
                    continue
                if t < CP.PropsSI("Tcrit", name):
                    if p > 0.5 * CP.PropsSI("P", "T", t, "Q", 1, name):
                        continue
                rho_r = CP.PropsSI("Dmolar", "T", t, "P", p, name) * vc
                # Margin to the engine's own PR-density domain checks.
                if rho_r > (0.08 if is_polar else 1.3):
                    continue
                if p == 1e4 and not is_polar:
                    continue
                cases.append({
                    "id": f"{cid}-gas-{t:g}K-{p:g}Pa", "component": cid, "phase": "gas",
                    "temperature_K": t, "pressure_Pa": p, "polar": is_polar, "coolprop_reduced_density": rho_r,
                    "expected": {
                        "viscosity": CP.PropsSI("V", "T", t, "P", p, name),
                        "thermal_conductivity": None if is_polar else CP.PropsSI("L", "T", t, "P", p, name),
                    },
                })
        lo = max(c[k]["temperature_range_K"][0] for k in ("liquid_viscosity", "liquid_thermal_conductivity", "surface_tension"))
        hi = min(c[k]["temperature_range_K"][1] for k in ("liquid_viscosity", "liquid_thermal_conductivity", "surface_tension"))
        for tr in (0.5, 0.6, 0.7, 0.74):
            t = round(tr * tc, 2)
            if not lo <= t <= hi or t < CP.PropsSI("Tmin", name) + 1 or t > LIQUID_MAX_TR * tc:
                continue
            psat = CP.PropsSI("P", "T", t, "Q", 0, name)
            for p in sorted({max(1.5 * psat, 1e5), 5e6}):
                if p < 1.5 * psat:
                    continue
                cases.append({
                    "id": f"{cid}-liquid-{t:g}K-{p:.6g}Pa", "component": cid, "phase": "liquid",
                    "temperature_K": t, "pressure_Pa": p, "polar": is_polar,
                    "expected": {
                        "viscosity": CP.PropsSI("V", "T", t, "P", p, name),
                        "thermal_conductivity": CP.PropsSI("L", "T", t, "P", p, name),
                        "surface_tension": CP.PropsSI("I", "T", t, "Q", 0, name),
                    },
                })
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    transport, base = load()
    implementation = []
    for cid, ids, x, t, v in IMPLEMENTATION_GAS:
        expected, rho_r = gas(ids, x, t, v, transport, base)
        implementation.append({"id": cid, "phase": "gas", "components": ids, "mole_fractions": x,
                               "temperature_K": t, "molar_volume_m3_per_mol": v, "reduced_density": rho_r,
                               "expected": expected})
    for cid, ids, x, t, p in IMPLEMENTATION_LIQUID:
        implementation.append({"id": cid, "phase": "liquid", "components": ids, "mole_fractions": x,
                               "temperature_K": t, "pressure_Pa": p, "expected": liquid(ids, x, t, p, transport, base)})
    fixture = {
        "schema": "transport-references-v1",
        "purpose": "implementation verification (matched formulation, independent Python) and model comparison against CoolProp 7.2.0 pure-fluid transport correlations",
        "generator": "validation/transport/generate-transport-references.py",
        "environment": {"python": sys.version.split()[0], "CoolProp": CoolProp.__version__, "chemicals": chemicals.__version__},
        "input_sha256": {str(p.relative_to(ROOT).as_posix()): hashlib.sha256(p.read_bytes()).hexdigest() for p in (TRANSPORT, COMPONENTS)},
        "comparison": "abs(actual-expected) <= atol + rtol*abs(expected); null expected means OutsideValidityDomain",
        "tolerances": TOLERANCES,
        "coolprop_gas_eos": "Rust evaluates gas states on the PR-exact-1976 largest-volume root; liquid states need no volume",
        "implementation": implementation,
        "coolprop": coolprop_cases(transport, base),
    }
    pathlib.Path(args.output).write_text(json.dumps(fixture, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
