"""Independent A17 equilibrium, Gibbs and kinetic reactor references.

Never invokes the Rust engine. Sources of independence:

* Equilibrium and Gibbs cases: Cantera 3.2.0 `equilibrate("TP")` of an
  ideal-gas mixture built from the pinned NASA7 records (including the
  formation constants a6, a7), or, for a user equilibrium-constant
  correlation, a SciPy bracketed solve of the single extent.
* Kinetic cases: rate laws written out here and solved with SciPy
  (`brentq` for the CSTR material balance, `solve_ivp` DOP853 for the plug-flow
  and batch profiles). The plug-flow and CSTR cases use partial-pressure and
  mole-fraction bases, so they need no equation of state. The batch case needs
  the reacting volume and uses a separate Peng–Robinson molar-volume
  implementation (exact critical constants derived below, zero kij).
* Duties: the ideal-gas part `H_out(T_out) - H_in(T_in)` on the formation basis
  from Cantera, compared by the Rust test after removing the PR residual
  enthalpy of the solved states.

Generate to an explicit candidate path and review before adoption.
"""
import argparse
import hashlib
import importlib.metadata
import json
import math
import platform
import re
from pathlib import Path

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
SOURCES = [ROOT / "data/source-extracts/components.provisional.json",
           ROOT / "data/source-extracts/components-reacting.provisional.json"]
DATA = {c["id"]: c for s in SOURCES for c in json.loads(s.read_text())["components"]}
COMPONENTS = ["methane", "water", "carbon-monoxide", "carbon-dioxide", "hydrogen", "nitrogen", "methanol"]
R = 8.31446261815324
P0 = 1e5
TOLERANCES = {"equilibrium_flow_relative": 1e-8, "kinetic_flow_relative": 1e-8,
              "duty_relative": 1e-7, "gibbs_energy_J_per_mol": 1e-6}

WGS = {"carbon-monoxide": -1, "water": -1, "carbon-dioxide": 1, "hydrogen": 1}
METHANATION = {"carbon-monoxide": -1, "hydrogen": -3, "methane": 1, "water": 1}

EQUILIBRIUM = [
    {"id": "gibbs-shift-methanation", "kind": "gibbs", "inerts": ["nitrogen"],
     "feed": {"temperature_K": 550.0, "pressure_Pa": 1e5, "flow_mol_per_s": 1.0,
              "mole_fractions": {"carbon-monoxide": 0.2, "water": 0.3, "carbon-dioxide": 0.05,
                                 "hydrogen": 0.05, "nitrogen": 0.4}},
     "outlet_temperature_K": 560.0, "temperature_approach_K": 0.0, "pressure_drop_Pa": 0.0},
    {"id": "equilibrium-water-gas-shift-approach", "kind": "equilibrium",
     "reactions": [{"coefficients": WGS, "constant": "gibbs"}],
     "feed": {"temperature_K": 500.0, "pressure_Pa": 2e5, "flow_mol_per_s": 1.0,
              "mole_fractions": {"carbon-monoxide": 0.25, "water": 0.35, "carbon-dioxide": 0.05,
                                 "hydrogen": 0.05, "nitrogen": 0.3}},
     "outlet_temperature_K": 520.0, "temperature_approach_K": 20.0, "pressure_drop_Pa": 0.0},
    {"id": "equilibrium-methanol-correlation", "kind": "equilibrium",
     "reactions": [{"coefficients": {"carbon-monoxide": -1, "hydrogen": -2, "methanol": 1},
                    "constant": {"a": -24.4, "b": 11000.0, "c": 0.0, "d": 0.0}}],
     "feed": {"temperature_K": 480.0, "pressure_Pa": 20e5, "flow_mol_per_s": 1.0,
              "mole_fractions": {"carbon-monoxide": 0.3, "hydrogen": 0.6, "nitrogen": 0.1}},
     "outlet_temperature_K": 500.0, "temperature_approach_K": 0.0, "pressure_drop_Pa": 1e5},
]

# Rate constants: k = k0 (T/Tref)^n exp(-E/R (1/T - 1/Tref)).
KINETIC = [
    {"id": "cstr-water-gas-shift-partial-pressure", "kind": "cstr",
     "feed": {"temperature_K": 500.0, "pressure_Pa": 5e5, "flow_mol_per_s": 1.0,
              "mole_fractions": {"carbon-monoxide": 0.2, "water": 0.3, "carbon-dioxide": 0.05,
                                 "hydrogen": 0.05, "nitrogen": 0.4}},
     "volume_m3": 0.5, "outlet_temperature_K": 550.0, "pressure_drop_Pa": 1e4,
     "catalyst_bulk_density_kg_m3": None,
     "reactions": [{"coefficients": WGS, "basis": "partial_pressure", "rate_basis": "volume",
                    "law": "power", "forward": {"k0": 2e-11, "E": 50000.0, "n": 0.0, "Tref": 500.0,
                                                "orders": {"carbon-monoxide": 1, "water": 1}},
                    "reverse": {"k0": 2e-12, "E": 80000.0, "n": 0.0, "Tref": 500.0,
                                "orders": {"carbon-dioxide": 1, "hydrogen": 1}}}]},
    {"id": "pfr-shift-and-lhhw-methanation-catalytic", "kind": "pfr",
     "feed": {"temperature_K": 500.0, "pressure_Pa": 6e5, "flow_mol_per_s": 0.5,
              "mole_fractions": {"carbon-monoxide": 0.15, "water": 0.25, "carbon-dioxide": 0.05,
                                 "hydrogen": 0.35, "nitrogen": 0.2}},
     "length_m": 2.0, "diameter_m": 0.1, "temperature_K": 560.0, "pressure_drop_Pa": 2e4,
     "catalyst_bulk_density_kg_m3": 1200.0,
     "reactions": [
         {"coefficients": WGS, "basis": "mole_fraction", "rate_basis": "catalyst",
          "law": "power", "forward": {"k0": 0.02, "E": 60000.0, "n": 0.0, "Tref": 560.0,
                                      "orders": {"carbon-monoxide": 1, "water": 1}},
          "reverse": {"k0": 0.004, "E": 90000.0, "n": 0.0, "Tref": 560.0,
                      "orders": {"carbon-dioxide": 1, "hydrogen": 1}}},
         {"coefficients": METHANATION, "basis": "partial_pressure", "rate_basis": "catalyst",
          "law": "lhhw", "kinetic": {"k0": 1e-3, "E": 70000.0, "n": 0.0, "Tref": 560.0},
          "driving_forward": {"a": -23.0, "b": 0.0, "c": 0.0, "d": 0.0,
                              "orders": {"carbon-monoxide": 1, "hydrogen": 1}},
          "adsorption": [{"a": 0.0, "b": 0.0, "c": 0.0, "d": 0.0, "orders": {}},
                         {"a": -14.0, "b": 2000.0, "c": 0.0, "d": 0.0, "orders": {"carbon-monoxide": 1}},
                         {"a": -6.9, "b": 0.0, "c": 0.0, "d": 0.0, "orders": {"hydrogen": 0.5}}],
          "exponent": 2.0}]},
    {"id": "batch-water-gas-shift-concentration", "kind": "batch",
     "feed": {"temperature_K": 480.0, "pressure_Pa": 3e5, "flow_mol_per_s": 2.0,
              "mole_fractions": {"carbon-monoxide": 0.2, "water": 0.3, "nitrogen": 0.5}},
     "reaction_time_s": 600.0, "temperature_K": 520.0, "catalyst_bulk_density_kg_m3": None,
     "reactions": [{"coefficients": WGS, "basis": "molar_concentration", "rate_basis": "volume",
                    "law": "power", "forward": {"k0": 5e-5, "E": 0.0, "n": 0.0, "Tref": 520.0,
                                                "orders": {"carbon-monoxide": 1, "water": 1}},
                    "reverse": None}]},
]


def elements(formula):
    counts = {}
    for symbol, digits in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
    return counts


def gas(ids):
    species = []
    for i in ids:
        low, high = DATA[i]["ideal_gas"]["coefficients"]
        s = ct.Species(i, elements(DATA[i]["formula"]))
        s.thermo = ct.NasaPoly2(200, 6000, 100000, [1000, *high, *low])
        species.append(s)
    return ct.Solution(thermo="ideal-gas", species=species)


GAS = None


def g0_over_rt(i, t):
    s = GAS.species(i).thermo
    return (s.h(t) - t * s.s(t)) / 1000.0 / (R * t)


def enthalpy_flow(flows, t):
    total = sum(flows.values())
    GAS.TPX = t, P0, {k: v / total for k, v in flows.items() if v > 0}
    return total * GAS.enthalpy_mole / 1000.0


def duty_scale(flows_in, t_in, flows_out, t_out):
    def s(flows, t):
        return sum(abs(n * GAS.species(k).thermo.h(t) / 1000.0) for k, n in flows.items())
    return s(flows_in, t_in) + s(flows_out, t_out)


def feed_flows(case):
    f = case["feed"]
    return {i: f["flow_mol_per_s"] * f["mole_fractions"].get(i, 0.0) for i in COMPONENTS}


def finish(case, inlet, outlet, t_out, extents=None):
    generation = {i: outlet[i] - inlet[i] for i in COMPONENTS}
    return {**case, "outlet_flows_mol_per_s": outlet, "extents_mol_per_s": extents or [],
            "standard_reaction_enthalpy_W": sum(n * GAS.species(i).thermo.h(298.15) / 1000.0
                                                for i, n in generation.items()),
            "ideal_gas_duty_W": enthalpy_flow(outlet, t_out) - enthalpy_flow(inlet, case["feed"]["temperature_K"]),
            "duty_scale_W": duty_scale(inlet, case["feed"]["temperature_K"], outlet, t_out)}


def equilibrium_case(case):
    inlet = feed_flows(case)
    t_eq = case["outlet_temperature_K"] + case["temperature_approach_K"]
    p = case["feed"]["pressure_Pa"] - case["pressure_drop_Pa"]
    if case["kind"] == "gibbs" or case["reactions"][0]["constant"] == "gibbs":
        # Cantera over exactly the species allowed to change, plus inerts as
        # spectators (they dilute but carry no shared element here).
        if case["kind"] == "gibbs":
            ids = COMPONENTS
        else:
            ids = sorted({*case["reactions"][0]["coefficients"],
                          *[i for i in COMPONENTS if inlet[i] > 0]})
        sub = gas(ids)
        total = sum(inlet.values())
        sub.TPX = t_eq, p, {i: inlet[i] / total for i in ids if inlet[i] > 0}
        sub.equilibrate("TP", rtol=1e-14, max_steps=10000)
        outlet = {i: 0.0 for i in COMPONENTS}
        # Cantera reports mole fractions; recover amounts from nitrogen, which
        # is inert and present in every equilibrium case.
        scale = inlet["nitrogen"] / sub.X[sub.species_index("nitrogen")]
        for i in ids:
            outlet[i] = scale * sub.X[sub.species_index(i)]
        extents = []
        if case["kind"] == "equilibrium":
            extents = [outlet["carbon-dioxide"] - inlet["carbon-dioxide"]]
    else:
        r = case["reactions"][0]
        nu = r["coefficients"]
        c = r["constant"]
        ln_k = c["a"] + c["b"] / t_eq + c["c"] * math.log(t_eq) + c["d"] * t_eq

        def residual(xi):
            n = {i: inlet[i] + nu.get(i, 0) * xi for i in COMPONENTS}
            total = sum(n.values())
            return sum(v * math.log(n[i] / total * p / P0) for i, v in nu.items()) - ln_k

        upper = min(inlet[i] / -v for i, v in nu.items() if v < 0)
        xi = brentq(residual, 1e-12, upper * (1 - 1e-12), xtol=1e-15, rtol=1e-15, maxiter=500)
        outlet = {i: inlet[i] + nu.get(i, 0) * xi for i in COMPONENTS}
        extents = [xi]
    return finish(case, inlet, outlet, case["outlet_temperature_K"], extents)


def arrhenius(k, t):
    return k["k0"] * (t / k["Tref"]) ** k["n"] * math.exp(-k["E"] / R * (1 / t - 1 / k["Tref"]))


def basis(reaction, x, p, molar_volume=None):
    b = reaction["basis"]
    if b == "partial_pressure":
        return {i: x[i] * p for i in COMPONENTS}
    if b == "mole_fraction":
        return dict(x)
    if b == "molar_concentration":
        return {i: x[i] / molar_volume for i in COMPONENTS}
    raise ValueError(b)


def product(q, orders):
    v = 1.0
    for i, a in orders.items():
        v *= q[i] ** a
    return v


def rate(reaction, t, x, p, molar_volume=None, bulk_density=None):
    q = basis(reaction, x, p, molar_volume)
    if reaction["law"] == "power":
        v = arrhenius(reaction["forward"], t) * product(q, reaction["forward"]["orders"])
        if reaction["reverse"]:
            v -= arrhenius(reaction["reverse"], t) * product(q, reaction["reverse"]["orders"])
    else:
        def lhhw(term):
            ln_k = term["a"] + term["b"] / t + term["c"] * math.log(t) + term["d"] * t
            return math.exp(ln_k) * product(q, term["orders"])
        v = arrhenius(reaction["kinetic"], t) * lhhw(reaction["driving_forward"])
        v /= sum(lhhw(a) for a in reaction["adsorption"]) ** reaction["exponent"]
    return v * (bulk_density if reaction["rate_basis"] == "catalyst" else 1.0)


def pr_omegas():
    # Critical conditions of Peng–Robinson: X^3 + 3X^2 + 3X ... closed form.
    x = (-1 + (6 * math.sqrt(2) + 8) ** (1 / 3) - (6 * math.sqrt(2) - 8) ** (1 / 3)) / 3
    return 8 * (5 * x + 1) / (49 - 37 * x), x / (x + 3)


def pr_molar_volume(x, t, p):
    """Largest-root (vapor) PR molar volume, PR 1976 alpha, zero kij."""
    omega_a, omega_b = pr_omegas()
    a = []
    b = []
    for i in COMPONENTS:
        c = DATA[i]["constants"]
        tc, pc, w = c["critical_temperature_K"], c["critical_pressure_Pa"], c["acentric_factor"]
        m = 0.37464 + 1.54226 * w - 0.26992 * w * w
        alpha = (1 + m * (1 - math.sqrt(t / tc))) ** 2
        a.append(omega_a * R * R * tc * tc / pc * alpha)
        b.append(omega_b * R * tc / pc)
    xs = [x[i] for i in COMPONENTS]
    am = sum(xi * xj * math.sqrt(ai * aj) for xi, ai in zip(xs, a) for xj, aj in zip(xs, a))
    bm = sum(xi * bi for xi, bi in zip(xs, b))
    big_a, big_b = am * p / (R * t) ** 2, bm * p / (R * t)
    roots = np.roots([1, -(1 - big_b), big_a - 3 * big_b ** 2 - 2 * big_b,
                      -(big_a * big_b - big_b ** 2 - big_b ** 3)])
    z = max(r.real for r in roots if abs(r.imag) < 1e-12)
    return z * R * t / p


def kinetic_case(case):
    inlet = feed_flows(case)
    reactions = case["reactions"]
    nu = [[r["coefficients"].get(i, 0) for i in COMPONENTS] for r in reactions]
    s, nr = len(COMPONENTS), len(reactions)
    rho = case["catalyst_bulk_density_kg_m3"]
    n_in = [inlet[i] for i in COMPONENTS]

    def x_of(n):
        total = sum(n)
        return {i: max(v, 0.0) / total for i, v in zip(COMPONENTS, n)}

    if case["kind"] == "cstr":
        t = case["outlet_temperature_K"]
        p = case["feed"]["pressure_Pa"] - case["pressure_drop_Pa"]
        assert nr == 1

        def residual(xi):
            n = [n_in[i] + nu[0][i] * xi for i in range(s)]
            return xi - case["volume_m3"] * rate(reactions[0], t, x_of(n), p, bulk_density=rho)

        upper = min(n_in[i] / -nu[0][i] for i in range(s) if nu[0][i] < 0)
        xi = brentq(residual, 0.0, upper * (1 - 1e-12), xtol=1e-16, rtol=1e-15, maxiter=500)
        extents = [xi]
        t_out = t
    elif case["kind"] == "pfr":
        t = case["temperature_K"]
        area = math.pi * case["diameter_m"] ** 2 / 4
        length = case["length_m"]
        p_in = case["feed"]["pressure_Pa"]

        def rhs(z, y):
            p = p_in - case["pressure_drop_Pa"] * z / length
            x = x_of(y[:s])
            r = [rate(rx, t, x, p, bulk_density=rho) for rx in reactions]
            return [area * sum(nu[j][i] * r[j] for j in range(nr)) for i in range(s)] + \
                   [area * rj for rj in r]

        sol = solve_ivp(rhs, (0, length), n_in + [0.0] * nr, method="DOP853",
                        rtol=1e-13, atol=1e-16)
        assert sol.success
        extents = list(sol.y[s:, -1])
        n_out = list(sol.y[:s, -1])
        t_out = t
    else:
        t = case["temperature_K"]
        p = case["feed"]["pressure_Pa"]
        z0 = [inlet[i] / sum(n_in) for i in COMPONENTS]

        def rhs(_, y):
            n = y[:s]
            x = x_of(n)
            v = pr_molar_volume(x, t, p)
            volume = sum(n) * v
            r = [volume * rate(rx, t, x, p, molar_volume=v, bulk_density=rho) for rx in reactions]
            return [sum(nu[j][i] * r[j] for j in range(nr)) for i in range(s)] + r

        sol = solve_ivp(rhs, (0, case["reaction_time_s"]), z0 + [0.0] * nr, method="DOP853",
                        rtol=1e-13, atol=1e-16)
        assert sol.success
        f = case["feed"]["flow_mol_per_s"]
        n_out = [f * v for v in sol.y[:s, -1]]
        extents = [f * v for v in sol.y[s:, -1]]
        t_out = t
    if case["kind"] == "cstr":
        n_out = [n_in[i] + nu[0][i] * extents[0] for i in range(s)]
    outlet = dict(zip(COMPONENTS, n_out))
    return finish(case, inlet, outlet, t_out, extents)


def main():
    global GAS
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    GAS = gas(COMPONENTS)
    omega_a, omega_b = pr_omegas()
    assert abs(omega_a - 0.4572355289213822) < 1e-15 and abs(omega_b - 0.07779607390388846) < 1e-16
    gibbs = {i: {str(t): g0_over_rt(i, t) * R * t for t in [300.0, 450.0, 600.0]} for i in COMPONENTS}
    output = {"schema": "a17-reactor-reference-v2", "generator": "validation/generate-reactor-references.py",
              "python": platform.python_version(),
              "packages": {n: importlib.metadata.version(n) for n in ["cantera", "numpy", "scipy"]},
              "input_sha256": {s.relative_to(ROOT).as_posix(): hashlib.sha256(s.read_bytes()).hexdigest()
                               for s in SOURCES},
              "independence": "Cantera NASA7 ideal-gas equilibrium and enthalpies; SciPy brentq/DOP853 solutions of rate laws written here; separate Python PR molar volume; no Rust calls",
              "scope": "Ideal-gas equilibrium (fugacity basis IdealGas) and isothermal kinetic reactors at fixed cases; the PR/PR-BM departure of the duty is removed by the Rust test",
              "components": COMPONENTS, "tolerances": TOLERANCES,
              "standard_gibbs_energy_J_per_mol": gibbs,
              "equilibrium_cases": [equilibrium_case(c) for c in EQUILIBRIUM],
              "kinetic_cases": [kinetic_case(c) for c in KINETIC]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {len(EQUILIBRIUM)} equilibrium and {len(KINETIC)} kinetic cases to {args.output}")


if __name__ == "__main__":
    main()
