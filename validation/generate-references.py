"""Explicit reference regeneration, never called by Cargo tests.

Original adapter for this repository; no recovered candidate code is used.
teqp evaluates residual Helmholtz derivatives by automatic differentiation.
NumPy solves the cubic; Cantera evaluates the original NASA7 correlations.
BM's above-Tc branch is represented algebraically by Twu [2*(1-1/d),1,d].
The branch is fixed for each requested temperature: derivative evidence is local,
and branch continuity is checked separately in Rust.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import urllib.request

import cantera as ct
import numpy as np
from scipy.optimize import brentq
import teqp

ROOT = Path(__file__).resolve().parents[1]
R = 8.31446261815324
OA, OB = 0.4572355289213822, 0.07779607390388846
SOURCE = ROOT / "data/source-extracts/components.provisional.json"
DATA = {c["id"]: c for c in json.loads(SOURCE.read_text())["components"]}

# Budgets selected before Rust comparisons, for regular double precision roots
# away from critical/spinodal singularities. They are implementation tolerances.
TOLERANCES = {
    "Z": [2e-10, 2e-9], "mass_density": [1e-7, 2e-9],
    "ln_phi": [2e-9, 2e-9], "h_res": [2e-5, 2e-9],
    "s_res": [2e-7, 2e-9], "h": [2e-5, 2e-9], "s": [2e-7, 2e-9],
    "cp": [2e-5, 2e-8], "cv": [2e-5, 2e-8], "sound": [2e-4, 2e-8],
}


def nasa(c, t):
    raw = c["ideal_gas"]
    low, high = raw["coefficients"]
    poly = ct.NasaPoly2(200, 6000, 100000, [1000, *high, *low])
    h = (poly.h(t) - poly.h(298.15)) / 1000
    s = (poly.s(t) - poly.s(298.15)) / 1000
    # Remove source fit offset at the join to reproduce cp-integrated convention.
    if t >= 1000:
        lp = ct.NasaPoly2(200, 6000, 100000, [1000, *low, *low])
        hp = ct.NasaPoly2(200, 6000, 100000, [1000, *high, *high])
        h = (hp.h(t) - hp.h(1000) + lp.h(1000) - lp.h(298.15)) / 1000
        s = (hp.s(t) - hp.s(1000) + lp.s(1000) - lp.s(298.15)) / 1000
        cp = hp.cp(t) / 1000
    else:
        cp = poly.cp(t) / 1000
    return h, s, cp


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


def case(ids, fractions, t, p, variant, kij=None):
    x = np.array(fractions)
    if kij is None:
        kij = np.zeros((len(ids), len(ids))).tolist()
    model, spec = model_at(ids, t, variant, kij)
    a, b = model.get_a(t, x), model.get_b(t, x)
    aa, bb = a*p/(R*t)**2, b*p/(R*t)
    roots = sorted(z.real for z in np.roots([1, bb-1, aa-3*bb*bb-2*bb, -aa*bb+bb*bb+bb**3])
                   if abs(z.imag) < 1e-10 and z.real > bb)
    result = []
    for z in roots:
        rho = p/(z*R*t)
        # Refine density using teqp pressure, independently of the polynomial.
        def residual(r):
            return r*R*t*(1+model.get_Ar01(t, r, x))-p
        rho = brentq(residual, rho*(1-1e-6), rho*(1+1e-6), xtol=min(1e-10,rho*1e-14), rtol=2e-14)
        assert abs(residual(rho)) <= 1e-10*p, "reference pressure refinement did not converge"
        z = p/(rho*R*t)
        ar = {key: getattr(model, "get_Ar"+key)(t, rho, x) for key in ["00", "01", "02", "10", "11", "20"]}
        stiffness = 1+2*ar["01"]+ar["02"]
        if stiffness <= 0:
            result.append({"Z": z, "mechanically_stable": False})
            continue
        assert stiffness > 1e-4
        pure = [nasa(DATA[i], t) for i in ids]
        h0 = sum(xx*v[0] for xx, v in zip(x, pure))
        s0 = sum(xx*v[1] for xx, v in zip(x, pure)) - R*np.log(p/100000) - R*sum(xx*np.log(xx) for xx in x if xx > 0)
        cp0 = sum(xx*v[2] for xx, v in zip(x, pure))
        hr = R*t*(ar["10"]+ar["01"])
        sr = R*(-ar["00"]+ar["10"]+np.log(z))
        cv = cp0-R-R*ar["20"]
        cp = cv+R*(1+ar["01"]-ar["11"])**2/stiffness
        mass = sum(xx*DATA[i]["constants"]["molar_mass_kg_per_mol"] for i, xx in zip(ids, x))
        result.append({"Z": z, "mechanically_stable": True, "mass_density": rho*mass,
                       "ln_phi": np.log(model.get_fugacity_coefficients(t, rho*x)).tolist(),
                       "h_res": hr, "s_res": sr, "h": h0+hr, "s": s0+sr,
                       "cp": cp, "cv": cv, "sound": float(np.sqrt(cp/cv*R*t*stiffness/mass)),
                       "stiffness": stiffness})
    return {"ids": ids, "x": fractions, "T": t, "p": p, "variant": variant,
            "kij": kij, "teqp_spec": spec, "roots": result}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True,
                        help="Explicit destination; review diffs before replacing fixtures")
    args = parser.parse_args()
    cases = []
    for variant in ["PR", "PR-BM"]:
        for ids, x in [(["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"], [.88, .06, .03, .02, .01]),
                       (["methane", "propane"], [.5, .5])]:
            for t in [320., 360., 400.]:
                for p in [1e5, 1e6, 5e6]:
                    cases.append(case(ids, x, t, p, variant))
        for component in DATA:
            cases.append(case([component], [1.], 400., 1e6, variant))
        cases.append(case(["propane"], [1.], 250., 1e5, variant))
        cases.append(case(["methane", "propane"], [.5, .5], 360., 2e6, variant, [[0., .035], [.035, 0.]]))
        cases.append(case(["methane", "propane", "nitrogen"], [.7, .3, 0.], 360., 2e6, variant))
        cases.append(case(["methane", "propane"], [.5, .5], 400., .01, variant))
    ideal = [{"id": i, "T": t, "h": nasa(c,t)[0], "s": nasa(c,t)[1], "cp": nasa(c,t)[2]}
             for i,c in DATA.items() for t in [200., 298.15, 600., 999.9, 1000., 1000.1, 3000., 6000.]]
    urls = ["https://raw.githubusercontent.com/usnistgov/teqp/v0.23.2/include/teqp/models/cubics/simple_cubics.hpp",
            "https://raw.githubusercontent.com/ClapeyronThermo/Clapeyron.jl/3440a15ba224664cac5b24441ac8adacbaf326f9/src/models/cubic/alphas/BM.jl"]
    source_hashes = {u: hashlib.sha256(urllib.request.urlopen(u).read()).hexdigest() for u in urls}
    output = {"schema": "a3-homogeneous-reference-v1", "generator": "validation/generate-references.py",
              "python": platform.python_version(), "packages": {n: importlib.metadata.version(n) for n in ["teqp", "numpy", "scipy", "cantera"]},
              "input_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(), "reference_source_sha256": source_hashes,
              "independence": "teqp residual AD, NumPy cubic roots, SciPy pressure refinement, Cantera NASA7; shared source parameters; no Rust calls or recovered code",
              "scope": "homogeneous fixed-composition roots, not global stability or flash, synthetic kij; implementation agreement, not experimental accuracy",
              "tolerances": TOLERANCES, "cases": cases, "ideal_gas": ideal}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(f"Wrote {len(cases)} cubic states, {len(ideal)} ideal states to {args.output}")


if __name__ == "__main__":
    main()
