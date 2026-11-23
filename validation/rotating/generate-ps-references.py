"""Explicit A14 independent PS flash and rotating-equipment references; no Rust calls.

Reuses the reviewed A3 teqp AD + Cantera NASA7 caloric adapter
(validation/generate-references.py) and the A4 SciPy/teqp equilibrium adapter
(validation/generate-flash-references.py), as the A5 PH generator does.

- PS targets come from independently calculated TP states (TP -> PS round
  trips), pure teqp superancillary coexistence (lever rule in entropy) and a
  narrow-boiling ethane/propane split.
- Isentropic compression, pumping and expansion outlet states are solved with
  SciPy brentq on temperature (entropy, then the efficiency-corrected
  enthalpy) over independent TP flashes; the pure propane expansion into the
  two-phase region uses the entropy/enthalpy lever rule at the teqp saturation
  temperature. A two-stage intercooled compressor chains the same steps.

Run manually in the pinned reference environment:

    python validation/rotating/generate-ps-references.py --output <candidate.json>

Review the candidate before replacing validation/rotating/ps-v1.json.
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

ROOT = Path(__file__).resolve().parents[2]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "validation" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


caloric = load("caloric", "generate-references.py")
equilibrium = load("equilibrium", "generate-flash-references.py")
# Fixed before any Rust comparison: regular double-precision implementation budgets
# (the same magnitudes as the A5 PH fixture; entropy as the A3 caloric "s").
TOLERANCES = {"temperature_K": [2e-6, 0.0], "enthalpy_J_per_mol": [2e-5, 2e-9],
              "entropy_J_per_mol_K": [2e-7, 2e-9], "phase_fraction": [2e-8, 0.0],
              "composition": [2e-8, 1e-7], "power_W": [2e-5, 2e-8]}
NG = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"]
XNG = [.70, .15, .13, .01, .01]


def root_props(ids, x, t, p, variant, root="minimum"):
    c = caloric.case(ids, list(x), t, p, variant)
    roots = [r for r in c["roots"] if r["mechanically_stable"]]
    if root == "liquid":
        r = roots[0]
    elif root == "vapor":
        r = roots[-1]
    else:
        r = min(roots, key=lambda r: np.dot(x, r["ln_phi"]))
    return r["h"], r["s"]


def tp_state(case_id, ids, x, t, p, variant):
    """Independent TP equilibrium with phase-weighted enthalpy and entropy."""
    c = equilibrium.flash_case(case_id, ids, x, t, p, variant)
    if c["phases"] == 1:
        c["target_h"], c["target_s"] = root_props(ids, x, t, p, variant)
    else:
        beta = c["vapor_fraction"]
        hv, sv = root_props(ids, c["vapor"], t, p, variant)
        hl, sl = root_props(ids, c["liquid"], t, p, variant)
        c.update(target_h=beta*hv+(1-beta)*hl, target_s=beta*sv+(1-beta)*sl,
                 vapor_h=hv, liquid_h=hl, vapor_s=sv, liquid_s=sl)
    return c


def solve_t(f, lo, hi):
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-15, maxiter=200)


def mixture_step(ids, x, t_in, p_in, p_out, efficiency, variant, compress, bracket):
    """Isentropic then efficiency-corrected outlet of a mixture by brentq over TP flashes."""
    inlet = tp_state("inlet", ids, x, t_in, p_in, variant)
    h_in, s_in = inlet["target_h"], inlet["target_s"]
    t_s = solve_t(lambda t: tp_state("trial", ids, x, t, p_out, variant)["target_s"]-s_in, *bracket)
    ideal = tp_state("isentropic", ids, x, t_s, p_out, variant)
    dh_s = ideal["target_h"]-h_in
    h_out = h_in + (dh_s/efficiency if compress else dh_s*efficiency)
    t_out = solve_t(lambda t: tp_state("trial", ids, x, t, p_out, variant)["target_h"]-h_out, *bracket)
    outlet = tp_state("outlet", ids, x, t_out, p_out, variant)
    assert abs(outlet["target_h"]-h_out) < 1e-6 and abs(ideal["target_s"]-s_in) < 1e-8
    return inlet, ideal, outlet, h_out


def summary(c):
    keep = ["T", "p", "phases", "vapor_fraction", "vapor", "liquid", "target_h", "target_s"]
    return {k: c[k] for k in keep if k in c}


def compressor(variant, case_id, ids, x, t_in, p_in, p_out, efficiency, flow, bracket, compress=True, kind="compressor"):
    inlet, ideal, outlet, h_out = mixture_step(ids, x, t_in, p_in, p_out, efficiency, variant, compress, bracket)
    return {"id": case_id, "kind": kind, "ids": ids, "z": x, "variant": variant, "molar_flow_mol_s": flow,
            "efficiency": efficiency, "inlet": summary(inlet), "isentropic_outlet": summary(ideal),
            "outlet": summary(outlet),
            "power_W": flow*(h_out-inlet["target_h"])}


def pure_saturation_at(component, t, variant):
    tc = caloric.DATA[component]["constants"]["critical_temperature_K"]
    return equilibrium.saturation(component, t/tc, variant)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ps_cases, units = [], []
    for v in ["PR", "PR-BM"]:
        # TP -> PS round trips (uniquely recoverable regular states).
        for cid, ids, x, t, p in [
            ("binary-split", ["methane", "propane"], [.5, .5], 250., 3e6),
            ("binary-vapor", ["methane", "propane"], [.5, .5], 330., 2e6),
            ("binary-liquid", ["methane", "propane"], [.05, .95], 250., 5e6),
            ("ng-split", NG, XNG, 230., 2e6),
            ("ng-vapor", NG, XNG, 300., 2e6),
        ]:
            ps_cases.append(tp_state(cid, ids, x, t, p, v))
        # Pure coexistence: lever rule in entropy at the saturation temperature.
        for component in ["propane", "water"]:
            t = 280. if component == "propane" else .7*caloric.DATA[component]["constants"]["critical_temperature_K"]
            sat = pure_saturation_at(component, t, v)
            t, p = sat["T"], sat["p"]
            hl, sl = root_props([component], [1.], t, p, v, "liquid")
            hv, sv = root_props([component], [1.], t, p, v, "vapor")
            for beta in [0., .25, .75, 1.]:
                ps_cases.append({"id": f"{component}-saturation-{beta}", "ids": [component], "z": [1.],
                                 "kij": [[0.]], "variant": v, "T": t, "p": p, "phases": 2,
                                 "vapor_fraction": beta, "vapor": [1.], "liquid": [1.],
                                 "target_s": beta*sv+(1-beta)*sl, "target_h": beta*hv+(1-beta)*hl,
                                 "liquid_s": sl, "vapor_s": sv, "liquid_h": hl, "vapor_h": hv})
        # Narrow-boiling mixture, as the PH fixture.
        p_sat = pure_saturation_at("propane", 280., v)["p"]
        ps_cases.append(tp_state("narrow-boiling", ["ethane", "propane"], [.001, .999], 280., p_sat*1.002, v))
        assert ps_cases[-1]["phases"] == 2

        # Rotating equipment.
        units.append(compressor(v, "ng-compressor", NG, XNG, 300., 2e6, 6e6, .75, 10., (300., 450.)))
        units.append(compressor(v, "ng-expander", NG, XNG, 250., 5e6, 1.5e6, .8, 10., (200., 250.),
                                compress=False, kind="expander"))
        units.append(compressor(v, "water-pump", ["water"], [1.], 300., 1e5, 5e6, .7, 5., (299., 305.),
                                kind="pump"))
        # Pure propane: compressed liquid expanded into coexistence (lever rule).
        sat_in = pure_saturation_at("propane", 320., v)
        p_in = 1.2*sat_in["p"]
        h_in, s_in = root_props(["propane"], [1.], 320., p_in, v)
        sat_out = pure_saturation_at("propane", 270., v)
        t_sat, p_out = sat_out["T"], sat_out["p"]
        hl, sl = root_props(["propane"], [1.], t_sat, p_out, v, "liquid")
        hv, sv = root_props(["propane"], [1.], t_sat, p_out, v, "vapor")
        beta_s = (s_in-sl)/(sv-sl)
        h_s = beta_s*hv+(1-beta_s)*hl
        h_out = h_in-.8*(h_in-h_s)
        beta = (h_out-hl)/(hv-hl)
        assert 0 < beta_s < beta < 1
        units.append({"id": "propane-two-phase-expander", "kind": "expander", "ids": ["propane"], "z": [1.],
                      "variant": v, "molar_flow_mol_s": 2., "efficiency": .8,
                      "inlet": {"T": 320., "p": p_in, "phases": 1, "target_h": h_in, "target_s": s_in},
                      "isentropic_outlet": {"T": t_sat, "p": p_out, "phases": 2, "vapor_fraction": beta_s,
                                            "target_h": h_s, "target_s": s_in},
                      "outlet": {"T": t_sat, "p": p_out, "phases": 2, "vapor_fraction": beta,
                                 "target_h": h_out, "target_s": beta*sv+(1-beta)*sl},
                      "power_W": 2.*(h_out-h_in)})
        # Two-stage intercooled compressor, equal pressure ratio, intercooler to 310 K.
        p_mid = 2e6*np.sqrt(3.)
        first = compressor(v, "stage-1", NG, XNG, 300., 2e6, p_mid, .75, 10., (300., 400.))
        cooled = tp_state("intercooled", NG, XNG, 310., p_mid, v)
        assert cooled["phases"] == 1
        second = compressor(v, "stage-2", NG, XNG, 310., p_mid, 6e6, .75, 10., (310., 420.))
        duty = 10.*(cooled["target_h"]-first["outlet"]["target_h"])
        units.append({"id": "ng-two-stage-intercooled", "kind": "multistage_compressor", "ids": NG, "z": XNG,
                      "variant": v, "molar_flow_mol_s": 10., "efficiency": .75, "stages": 2,
                      "intercooler_T": 310., "inlet": first["inlet"], "stage_outlets": [first["outlet"], second["outlet"]],
                      "intercooler_outlet": summary(cooled),
                      "stage_power_W": [first["power_W"], second["power_W"]],
                      "power_W": first["power_W"]+second["power_W"], "intercooler_duty_W": duty,
                      "outlet": second["outlet"]})
    output = {"schema": "a14-ps-reference-v1", "generator": "validation/rotating/generate-ps-references.py",
              "python": platform.python_version(),
              "packages": {n: importlib.metadata.version(n) for n in ["teqp", "numpy", "scipy", "cantera"]},
              "input_sha256": hashlib.sha256(caloric.SOURCE.read_bytes()).hexdigest(),
              "adapter_sha256": {name: hashlib.sha256((ROOT / "validation" / name).read_bytes()).hexdigest()
                                 for name in ["generate-references.py", "generate-flash-references.py",
                                              "rotating/generate-ps-references.py"]},
              "independence": "teqp residual AD, Cantera NASA7 offsets, SciPy equilibrium and brentq outlet solves, teqp pure VLE; shared parameters; no Rust calls",
              "scope": "Implementation agreement for selected regular PR/PR-BM PS states and isentropic-efficiency pump/compressor/expander outlets; finite-trial stability evidence; no experimental accuracy or critical qualification",
              "tolerances": TOLERANCES, "ps": ps_cases, "units": units}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(f"Wrote {len(ps_cases)} PS states and {len(units)} rotating-unit references to {args.output}")
    for u in units:
        print(u["variant"], u["id"], u["outlet"]["T"], u["outlet"].get("phases"), u.get("power_W"))


if __name__ == "__main__":
    main()
