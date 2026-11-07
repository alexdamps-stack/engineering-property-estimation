"""Explicit A4 reference regeneration, never called by Cargo tests.

Original adapter for this repository; no Rust code or recovered candidate code is used.
teqp 0.23.2 supplies residual Helmholtz derivatives and fugacity coefficients by
automatic differentiation; NumPy roots plus SciPy pressure refinement select
densities; SciPy solves phase equilibrium with a *different* algorithm from the
Rust successive substitution:

- pure saturation: teqp's cubic superancillary, refined by teqp ``pure_VLE_T``;
- stability: global multi-start L-BFGS-B minimization of the modified TPD over
  ln W (Wilson vapor/liquid starts plus 48 seeded Dirichlet starts);
- two-phase split: SciPy ``root`` (hybrid Powell) on ln K with a brentq
  Rachford-Rice inner solve, followed by the same multi-start TPD search on each
  phase. Phase stability from finite starts is evidence, not proof.

The model formulation matches validation/generate-references.py (exact PR factors,
Mathias-Copeman [m,0,0] for the PR alpha and Twu for the BM branch above Tc).
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
from scipy.optimize import brentq, minimize, root
import teqp

ROOT = Path(__file__).resolve().parents[1]
R = 8.31446261815324
OA, OB = 0.4572355289213822, 0.07779607390388846
SOURCE = ROOT / "data/source-extracts/components.provisional.json"
DATA = {c["id"]: c for c in json.loads(SOURCE.read_text())["components"]}

# Budgets fixed before any Rust comparison: implementation agreement for converged
# regular states (fugacity residual 1e-10 in Rust, ~1e-12 here), not accuracy claims.
TOLERANCES = {
    "vapor_fraction": [2e-8, 0.0],
    "composition": [2e-8, 1e-7],
    "Z": [2e-9, 2e-8],
    "saturation_pressure": [1e-6, 2e-9],
    "saturation_Z": [2e-10, 2e-9],
    "stable_tm_floor": -1e-8,
}
SEED = 0xA4A4
RANDOM_STARTS = 48


def model_at(ids, t, variant, kij):
    constants = [DATA[i]["constants"] for i in ids]
    alphas = []
    for c in constants:
        w = c["acentric_factor"]
        m = .37464 + 1.54226*w - .26992*w*w
        if variant == "PR-BM" and t > c["critical_temperature_K"]:
            d = 1 + m/2
            alphas.append({"type": "Twu", "c": [2*(1-1/d), 1, d]})
        else:
            alphas.append({"type": "Mathias-Copeman", "c": [m, 0, 0]})
    spec = {"type": "PR", "Tcrit / K": [c["critical_temperature_K"] for c in constants],
            "pcrit / Pa": [c["critical_pressure_Pa"] for c in constants],
            "acentric": [c["acentric_factor"] for c in constants],
            "alpha": alphas, "kmat": kij}
    model = teqp.make_model({"kind": "cubic", "model": spec})
    meta = model.get_meta()
    assert abs(meta["OmegaA"] - OA) < 1e-16 and abs(meta["OmegaB"] - OB) < 1e-16
    assert model.get_R(np.ones(len(ids))/len(ids)) == R
    return model, spec


def roots(model, t, p, x):
    """Mechanically stable roots as (Z, rho, ln_phi), refined against teqp pressure."""
    a, b = model.get_a(t, x), model.get_b(t, x)
    aa, bb = a*p/(R*t)**2, b*p/(R*t)
    found = sorted(z.real for z in np.roots([1, bb-1, aa-3*bb*bb-2*bb, -aa*bb+bb*bb+bb**3])
                   if abs(z.imag) < 1e-10 and z.real > bb)
    out = []
    for z in found:
        rho = p/(z*R*t)

        def residual(r):
            return r*R*t*(1+model.get_Ar01(t, r, x))-p
        rho = brentq(residual, rho*(1-1e-6), rho*(1+1e-6), xtol=min(1e-10, rho*1e-14), rtol=2e-14)
        if 1+2*model.get_Ar01(t, rho, x)+model.get_Ar02(t, rho, x) <= 0:
            continue
        ln_phi = np.log(model.get_fugacity_coefficients(t, rho*x))
        out.append((p/(rho*R*t), rho, ln_phi))
    assert out, "no mechanically stable root"
    return out


def min_gibbs(model, t, p, x):
    return min(roots(model, t, p, x), key=lambda r: float(np.dot(x, r[2])))


def wilson(ids, t, p):
    c = [DATA[i]["constants"] for i in ids]
    return np.array([k["critical_pressure_Pa"]/p*np.exp(5.373*(1+k["acentric_factor"])*(1-k["critical_temperature_K"]/t)) for k in c])


def tpd_search(model, ids, t, p, z):
    """Global multi-start minimization of tm over ln W; returns (min tm, w, Z)."""
    active = z > 0
    if active.sum() < 2:
        return 0.0, z.tolist(), min_gibbs(model, t, p, z)[0]
    zf = min_gibbs(model, t, p, z)
    d = np.where(active, np.log(np.where(active, z, 1)) + zf[2], 0.0)
    k = wilson(ids, t, p)
    rng = np.random.default_rng(SEED)
    starts = [z*k, z/k] + [rng.dirichlet(np.ones(len(z)))*np.where(active, 1, 0) for _ in range(RANDOM_STARTS)]

    def tm(ln_w_active):
        w_big = np.zeros(len(z))
        w_big[active] = np.exp(ln_w_active)
        w = w_big/w_big.sum()
        ln_phi = min_gibbs(model, t, p, w)[2]
        g = np.log(w_big[active]) + ln_phi[active] - d[active]
        value = 1 + np.sum(w_big[active]*(g - 1))
        return value, w_big[active]*g
    best = (np.inf, None)
    for s in starts:
        s = s[active]
        s = np.maximum(s, 1e-12)
        result = minimize(tm, np.log(s), jac=True, method="L-BFGS-B", bounds=[(-60.0, 5.0)]*len(s),
                          options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 2000})
        if result.fun < best[0]:
            best = (float(result.fun), result.x)
    w = np.zeros(len(z))
    w[active] = np.exp(best[1])
    w /= w.sum()
    return best[0], w.tolist(), min_gibbs(model, t, p, w)[0]


def rachford_rice(z, k):
    """Active components only (z > 0); inactive entries stay zero in both phases."""
    a = z > 0
    za, ka = z[a], k[a]

    def g(beta):
        return np.sum(za*(ka-1)/(1+beta*(ka-1)))
    lo, hi = 1/(1-ka.max()), 1/(1-ka.min())
    span = hi - lo
    beta = brentq(g, lo+1e-14*span, hi-1e-14*span, xtol=1e-16, rtol=1e-15, maxiter=500)
    x = np.zeros(len(z))
    x[a] = za/(1+beta*(ka-1))
    y = k*x
    return beta, x/x.sum(), y/y.sum()


def split(model, t, p, z, ln_k0):
    a = z > 0

    def residual(ln_k_active):
        ln_k = np.zeros(len(z))
        ln_k[a] = ln_k_active
        beta, x, y = rachford_rice(z, np.exp(ln_k))
        lv, ll = min_gibbs(model, t, p, y)[2], min_gibbs(model, t, p, x)[2]
        return ln_k_active - (ll - lv)[a]
    sol = root(residual, ln_k0[a], method="hybr", options={"xtol": 1e-13, "maxfev": 20000})
    # MINPACK may report "xtol too small" at machine precision; the explicit
    # fugacity residual below is the acceptance criterion.
    assert sol.success or "xtol" in sol.message, sol.message
    k = np.ones(len(z))
    k[a] = np.exp(sol.x)
    beta, x, y = rachford_rice(z, k)
    fv, fl = min_gibbs(model, t, p, y), min_gibbs(model, t, p, x)
    eq = np.max(np.abs(np.log(y[a])+fv[2][a]-np.log(x[a])-fl[2][a]))
    assert eq < 1e-11, f"fugacity residual {eq}"
    return beta, x, y, fv[0], fl[0], eq


def flash_case(case_id, ids, fractions, t, p, variant, kij=None):
    z = np.array(fractions, dtype=float)
    kij = kij or np.zeros((len(ids), len(ids))).tolist()
    model, spec = model_at(ids, t, variant, kij)
    tm, w, zw = tpd_search(model, ids, t, p, z)
    zf = min_gibbs(model, t, p, z)
    base = {"id": case_id, "ids": ids, "z": fractions, "T": t, "p": p, "variant": variant,
            "kij": kij, "teqp_spec": spec, "feed_min_tm": tm}
    if tm >= TOLERANCES["stable_tm_floor"]:
        return {**base, "phases": 1, "Z": zf[0]}
    vapor_trial = zw > zf[0]
    ratio = np.zeros(len(z))
    ratio[z > 0] = np.log(np.array(w)[z > 0]/z[z > 0])
    beta, x, y, z_v, z_l, eq = split(model, t, p, z, ratio if vapor_trial else -ratio)
    assert 0 < beta < 1
    for phase in (x, y):
        phase_tm = tpd_search(model, ids, t, p, phase)[0]
        assert phase_tm >= TOLERANCES["stable_tm_floor"], f"{case_id}: split phase unstable {phase_tm}"
    def g(w, r):
        a = w > 0
        return float(np.sum(w[a]*(np.log(w[a])+r[a])))
    dg = beta*g(y, min_gibbs(model, t, p, y)[2]) + (1-beta)*g(x, min_gibbs(model, t, p, x)[2]) - g(z, zf[2])
    assert dg < 0
    return {**base, "phases": 2, "vapor_fraction": beta, "vapor": y.tolist(), "liquid": x.tolist(),
            "Z_vapor": z_v, "Z_liquid": z_l, "fugacity_residual": eq, "reduced_gibbs_change": dg}


def saturation(component, tr, variant):
    tc = DATA[component]["constants"]["critical_temperature_K"]
    t = tr*tc
    model, _ = model_at([component], t, variant, [[0.0]])
    rl, rv = model.superanc_rhoLV(t)
    rl, rv = model.pure_VLE_T(t, rl, rv, 20)
    x = np.array([1.0])
    p = rv*R*t*(1+model.get_Ar01(t, rv, x))
    pl = rl*R*t*(1+model.get_Ar01(t, rl, x))
    assert abs(p-pl) <= 1e-9*p
    return {"component": component, "variant": variant, "T": t, "Tr": tr, "p": p,
            "Z_liquid": p/(rl*R*t), "Z_vapor": p/(rv*R*t)}


NG = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"]
# Mixed demonstration feed: 100 mol/s of A plus 20 mol/s of B (process.provisional.json).
MIX = [(.88*100+.70*20)/120, (.06*100+.15*20)/120, (.03*100+.13*20)/120,
       (.02*100+.01*20)/120, (.01*100+.01*20)/120]
CASES = [
    ("c1c3-250K-3MPa", ["methane", "propane"], [.5, .5], 250., 3e6, None),
    ("c1c3-330K-2MPa", ["methane", "propane"], [.5, .5], 330., 2e6, None),
    ("c1c3-lean-300K-1.5MPa", ["methane", "propane"], [.2, .8], 300., 1.5e6, None),
    ("c1c3-liquid-250K-5MPa", ["methane", "propane"], [.05, .95], 250., 5e6, None),
    ("c1c3-trace-n2-250K-3MPa", ["methane", "propane", "nitrogen"], [.5, .5, 0.], 250., 3e6, None),
    ("c1co2-kij-230K-4MPa", ["methane", "carbon-dioxide"], [.5, .5], 230., 4e6, [[0., .0919], [.0919, 0.]]),
    ("ng-A-330K-5MPa", NG, [.88, .06, .03, .02, .01], 330., 5e6, None),
    ("ng-mix-240K-5MPa", NG, MIX, 240., 5e6, None),
    ("ng-mix-220K-1.5MPa", NG, MIX, 220., 1.5e6, None),
    ("ng-B-250K-3MPa", NG, [.70, .15, .13, .01, .01], 250., 3e6, None),
    ("ng-B-230K-2MPa", NG, [.70, .15, .13, .01, .01], 230., 2e6, None),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True,
                        help="Explicit destination; review diffs before replacing fixtures")
    args = parser.parse_args()
    flashes = [flash_case(cid, ids, x, t, p, v, kij)
               for v in ["PR", "PR-BM"] for cid, ids, x, t, p, kij in CASES]
    sats = [saturation(c, tr, v) for v in ["PR", "PR-BM"]
            for c in ["ethane", "propane", "carbon-dioxide", "water", "methanol"]
            for tr in [.7, .8, .9] if 200 <= tr*DATA[c]["constants"]["critical_temperature_K"] <= 600]
    output = {"schema": "a4-flash-reference-v1", "generator": "validation/generate-flash-references.py",
              "python": platform.python_version(),
              "packages": {n: importlib.metadata.version(n) for n in ["teqp", "numpy", "scipy"]},
              "input_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "seed": SEED, "random_tpd_starts": RANDOM_STARTS,
              "independence": "teqp AD fugacities, NumPy/SciPy densities, multi-start L-BFGS-B TPD minimization, hybrid-Powell ln K solve, teqp superancillary saturation; shared source parameters; no Rust calls",
              "scope": "PR/PR-BM implementation agreement for stability classification, two-phase splits and pure saturation; synthetic kij; finite-start stability evidence; not experimental accuracy",
              "tolerances": TOLERANCES, "flash": flashes, "saturation": sats}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(f"Wrote {len(flashes)} flash and {len(sats)} saturation states to {args.output}")
    for f in flashes:
        print(f["variant"], f["id"], f["phases"], f.get("vapor_fraction"), round(f["feed_min_tm"], 6))


if __name__ == "__main__":
    main()
