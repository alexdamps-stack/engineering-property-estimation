"""Explicit A5 independent caloric/PH reference generation; no Rust calls.

Reuses the reviewed A3 teqp AD + Cantera NASA7 caloric adapter and A4
SciPy/teqp equilibrium adapter. PH targets come from independently calculated
TP states, pure superancillary coexistence and an isenthalpic valve solved
with SciPy brentq. Run manually; review the resulting fixture before replacing.
"""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "validation" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


caloric = load("caloric", "generate-references.py")
equilibrium = load("equilibrium", "generate-flash-references.py")
# Fixed before comparison: regular double-precision implementation budgets.
TOLERANCES = {"temperature_K": [2e-6, 0.0], "enthalpy_J_per_mol": [2e-5, 2e-9],
              "phase_fraction": [2e-8, 0.0], "composition": [2e-8, 1e-7]}


def enthalpy(ids, x, t, p, variant, root="minimum"):
    c = caloric.case(ids, list(x), t, p, variant)
    roots = [r for r in c["roots"] if r["mechanically_stable"]]
    if root == "liquid":
        return roots[0]["h"]
    if root == "vapor":
        return roots[-1]["h"]
    def gibbs(r):
        return np.dot(x, r["ln_phi"])
    return min(roots, key=gibbs)["h"]


def tp_target(case_id, ids, x, t, p, variant):
    c = equilibrium.flash_case(case_id, ids, x, t, p, variant)
    if c["phases"] == 1:
        c["target_h"] = enthalpy(ids, x, t, p, variant)
    else:
        beta = c["vapor_fraction"]
        hv = enthalpy(ids, c["vapor"], t, p, variant)
        hl = enthalpy(ids, c["liquid"], t, p, variant)
        c.update(target_h=beta*hv+(1-beta)*hl, vapor_h=hv, liquid_h=hl)
    return c


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = []
    ng = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"]
    xng = [.70, .15, .13, .01, .01]
    for v in ["PR", "PR-BM"]:
        for cid, ids, x, t, p in [
            ("binary-split", ["methane", "propane"], [.5, .5], 250., 3e6),
            ("binary-vapor", ["methane", "propane"], [.5, .5], 330., 2e6),
            ("binary-liquid", ["methane", "propane"], [.05, .95], 250., 5e6),
            ("ng-split", ng, xng, 230., 2e6),
        ]:
            cases.append(tp_target(cid, ids, x, t, p, v))
        for component in ["propane", "water"]:
            tr = 280. / caloric.DATA[component]["constants"]["critical_temperature_K"] if component == "propane" else .7
            sat = equilibrium.saturation(component, tr, v)
            t, p = sat["T"], sat["p"]
            hl = enthalpy([component], [1.], t, p, v, "liquid")
            hv = enthalpy([component], [1.], t, p, v, "vapor")
            for beta in [0., .25, .75, 1.]:
                cases.append({"id": f"{component}-saturation-{beta}", "ids": [component], "z": [1.], "kij": [[0.]],
                              "variant": v, "T": t, "p": p, "phases": 2, "vapor_fraction": beta,
                              "vapor": [1.], "liquid": [1.], "target_h": beta*hv+(1-beta)*hl,
                              "liquid_h": hl, "vapor_h": hv})
        # Narrow coexistence: trace ethane in propane; choose pressure so an
        # interior split exists at 280 K, verified by the independent solver.
        ps = equilibrium.saturation("propane", 280./caloric.DATA["propane"]["constants"]["critical_temperature_K"], v)["p"]
        cases.append(tp_target("narrow-boiling", ["ethane", "propane"], [.001, .999], 280., ps*1.002, v))
        assert cases[-1]["phases"] == 2
        inlet = tp_target("valve-inlet", ng, xng, 250., 3e6, v)
        h = inlet["target_h"]
        def residual(t):
            return tp_target("valve-trial", ng, xng, t, 1.5e6, v)["target_h"]-h
        tout = brentq(residual, 230., 250., xtol=1e-9)
        outlet = tp_target("isenthalpic-valve-outlet", ng, xng, tout, 1.5e6, v)
        assert abs(outlet["target_h"]-h) < 1e-6
        outlet.update(target_h=h, inlet_T=250., inlet_p=3e6)
        cases.append(outlet)
    output = {"schema": "a5-ph-reference-v1", "generator": "validation/generate-ph-references.py",
              "python": platform.python_version(),
              "packages": {n: importlib.metadata.version(n) for n in ["teqp", "numpy", "scipy", "cantera"]},
              "input_sha256": hashlib.sha256(caloric.SOURCE.read_bytes()).hexdigest(),
              "adapter_sha256": {name: hashlib.sha256((ROOT / "validation" / name).read_bytes()).hexdigest()
                                 for name in ["generate-references.py", "generate-flash-references.py", "generate-ph-references.py"]},
              "independence": "teqp residual AD, Cantera NASA7 offsets, SciPy equilibrium/valve and teqp pure VLE; shared parameters; no Rust calls",
              "scope": "Implementation agreement for selected regular PR/PR-BM PH states; finite-trial stability evidence; no experimental accuracy or critical qualification",
              "tolerances": TOLERANCES, "cases": cases}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(f"Wrote {len(cases)} PH reference states to {args.output}")


if __name__ == "__main__":
    main()
