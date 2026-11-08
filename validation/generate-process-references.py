"""Independent A2 process qualification; never invokes the Rust engine.

Retains the A3 caloric and A4 finite-start stability algorithms and budgets.
Generate to an explicit candidate path and review before adoption.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

from scipy.optimize import brentq

# The PH adapter loads the independently implemented caloric/equilibrium adapters.
import importlib.util

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ph", ROOT / "validation/generate-ph-references.py")
ph = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ph)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    process = json.loads((ROOT / "examples/process.provisional.json").read_text())
    ids = process["components"]
    a, b = process["feeds"]
    flow = a["flow_mol_per_s"] + b["flow_mol_per_s"]
    z = [(a["flow_mol_per_s"]*x+b["flow_mol_per_s"]*y)/flow
         for x, y in zip(a["mole_fractions"], b["mole_fractions"])]
    cases = []
    for variant in ["PR", "PR-BM"]:
        def state(name, x, t, p):
            print(variant, name, t, flush=True)
            return ph.tp_target(name, ids, x, t, p, variant)
        feeds = [state(f["id"], f["mole_fractions"], f["temperature_K"], f["pressure_Pa"])
                 for f in [a, b]]
        hmix = sum(f["flow_mol_per_s"]*s["target_h"] for f, s in zip([a, b], feeds))/flow
        # Mixer is on the vapor branch; qualify the solved state with full TPD below.
        tmix = brentq(lambda t: ph.enthalpy(ids, z, t, 5e6, variant)-hmix, 330., 350., xtol=1e-10)
        mix = state("mix", z, tmix, 5e6)
        assert mix["phases"] == 1 and abs(mix["target_h"]-hmix) < 1e-6
        cool = state("cool", z, 240., 5e6)
        tout = brentq(lambda t: ph.tp_target("trial", ids, z, t, 1.5e6, variant)["target_h"]-cool["target_h"],
                      200., 240., xtol=1e-9)
        expand = state("expand", z, tout, 1.5e6)
        assert abs(expand["target_h"]-cool["target_h"]) < 1e-6
        expand["target_h"] = cool["target_h"]
        streams = feeds + [mix, cool, expand]
        if expand["phases"] == 2:
            for name, key, fraction in [("separator_vapor", "vapor", expand["vapor_fraction"]),
                                        ("separator_liquid", "liquid", 1-expand["vapor_fraction"])]:
                s = state(name, expand[key], tout, 1.5e6)
                assert s["phases"] == 1
                s["flow_mol_per_s"] = flow*fraction
                streams.append(s)
        else:
            raise AssertionError("Review demonstration: no separator split")
        for s in streams:
            s.setdefault("flow_mol_per_s", a["flow_mol_per_s"] if s["id"] == "A" else
                         b["flow_mol_per_s"] if s["id"] == "B" else flow)
        cases.append({"variant": variant, "streams": streams,
                      "cooler_duty_W": flow*(cool["target_h"]-mix["target_h"])})
    files = ["examples/process.provisional.json", "data/source-extracts/components.provisional.json"]
    files += ["validation/"+n for n in ["generate-process-references.py", "generate-ph-references.py",
                                       "generate-references.py", "generate-flash-references.py"]]
    output = {"schema": "a2-process-reference-v1", "packages": {n: importlib.metadata.version(n)
              for n in ["teqp", "cantera", "numpy", "scipy"]},
              "input_sha256": {n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
              "independence": "teqp AD, Cantera NASA7, SciPy brentq energy solves and A4 multi-start TPD/split; no Rust calls",
              "scope": "Seven fixed streams per variant; synthetic zero kij; finite-start fluid stability, not experimental accuracy or a continuous state-domain qualification",
              "tolerances": ph.TOLERANCES, "cases": cases}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
