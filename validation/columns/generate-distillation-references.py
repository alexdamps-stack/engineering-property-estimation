"""Independent binary MESH references for A17 rigorous distillation (#186).

SciPy solves directly for temperatures, binary phase fractions and vapor
flows. teqp gives the fugacity coefficients and residual properties, and
Cantera supplies the ideal-gas enthalpy through the reviewed reference adapters.
No Rust binaries or outputs are used.
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

PATH = Path(__file__).resolve()
ADAPTER = PATH.with_name("generate-shortcut-references.py")
spec = importlib.util.spec_from_file_location("shortcut_reference", ADAPTER)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

IDS = ["ethane", "propane"]
KIJ = [[0, 0], [0, 0]]
# Allow propagation from the existing phase-reference budgets through a short
# nonlinear cascade. These budgets were selected before comparing to Rust.
TOLERANCES = {"temperature_K": [2e-5, 0], "composition": [2e-7, 0],
              "molar_flow_mol_s": [1e-8, 0], "duty_W": [0.02, 2e-6]}


def feed(temperature, pressure, molar_flow, ethane, stage):
    return dict(temperature_K=temperature, pressure_Pa=pressure,
                molar_flow_mol_s=molar_flow, z=[ethane, 1-ethane], stage=stage)


CASES = [
    dict(id="binary-trays-PR", variant="PR", pressure_Pa=2e6, internal_stages=4,
         reflux_ratio=2.5, distillate_fraction=0.45,
         feeds=[feed(295, 2.3e6, 100, 0.45, 4)], draws=[]),
    dict(id="binary-packed-two-feeds-PR-BM", variant="PR-BM", pressure_Pa=2e6,
         internal_stages=4, reflux_ratio=2.5, distillate_fraction=0.45,
         feeds=[feed(307, 2.3e6, 45, 0.58, 2), feed(292, 2e6, 55, 0.35, 5)],
         draws=[dict(stage=4, molar_flow_mol_s=6)], packed_height_m=1.6, hetp_m=0.4),
]


def reference(case):
    n = case["internal_stages"] + 2
    p, variant = case["pressure_Pa"], case["variant"]
    feed_totals, feed_ethane, feed_enthalpy = np.zeros((3, n))
    for inlet in case["feeds"]:
        j, f = inlet["stage"] - 1, inlet["molar_flow_mol_s"]
        flashed = ref.ps.tp_state(
            "feed", IDS, inlet["z"], inlet["temperature_K"], inlet["pressure_Pa"], variant
        )
        feed_totals[j] += f
        feed_ethane[j] += f*inlet["z"][0]
        feed_enthalpy[j] += f*flashed["target_h"]
    flow = feed_totals.sum()
    distillate = flow*case["distillate_fraction"]
    reflux = distillate*case["reflux_ratio"]
    draws = np.zeros(n)
    draws[0] = distillate
    for draw in case["draws"]:
        draws[draw["stage"] - 1] += draw["molar_flow_mol_s"]
    net = np.cumsum(feed_totals-draws)

    def properties(values):
        temperatures = values[:n]
        y = values[2*n-1:3*n-1]
        x = np.r_[y[1], values[n:2*n-1]]
        vapor = np.r_[0, reflux+distillate, values[3*n-1:]]
        liquid = np.r_[reflux, vapor[2:]+net[1:-1], net[-1]]
        equilibrium = np.empty((n, 2))
        h_l = np.empty(n)
        h_v = np.empty(n)
        for j, t in enumerate(temperatures):
            xl, yv = np.array([x[j], 1-x[j]]), np.array([y[j], 1-y[j]])
            model, _ = ref.equilibrium.model_at(IDS, t, variant, KIJ)
            ln_k = ref.ln_phi(model, t, p, xl, "liquid") - ref.ln_phi(model, t, p, yv, "vapor")
            equilibrium[j] = xl*np.exp(ln_k)-yv
            h_l[j] = ref.phase_h(IDS, xl, t, p, variant, KIJ, "liquid")
            h_v[j] = ref.phase_h(IDS, yv, t, p, variant, KIJ, "vapor")
        component = (np.r_[0, (liquid*x)[:-1]] + np.r_[(vapor*y)[1:], 0]
                     + feed_ethane - (liquid+draws)*x - vapor*y)
        energy = (np.r_[0, (liquid*h_l)[:-1]] + np.r_[(vapor*h_v)[1:], 0]
                  + feed_enthalpy - (liquid+draws)*h_l - vapor*h_v)
        return x, liquid, equilibrium, component, energy

    def residual(values):
        _, _, eq, component, energy = properties(values)
        return np.r_[eq.ravel(), component[1:]/flow, energy[1:-1]/(flow*10_000)]

    z = feed_ethane.sum()/flow
    bubble = ref.sat_point(IDS, [z, 1-z], p, variant, KIJ, "bubble")
    floor = np.maximum(0, -net[1:-1]) + 1e-8
    initial = np.r_[np.full(n, bubble["T"]), np.full(n-1, z),
                    np.full(n, bubble["vapor"][0]), floor+reflux+distillate]
    lo = np.r_[np.full(n, 200), np.full(2*n-1, 1e-9), floor]
    hi = np.r_[np.full(n, 350), np.full(2*n-1, 1-1e-9), np.full(n-2, np.inf)]
    fit = least_squares(residual, initial, bounds=(lo, hi), x_scale="jac",
                        ftol=1e-13, gtol=1e-13, xtol=1e-13)
    assert fit.success and np.max(np.abs(fit.fun)) < 1e-10, (case["id"], fit.message, fit.fun)
    x, liquid, _, _, energy = properties(fit.x)
    stage_numbers = [0, n-1] + [d["stage"]-1 for d in case["draws"]]
    flows = [distillate, float(liquid[-1])] + [d["molar_flow_mol_s"] for d in case["draws"]]
    assert energy[0] > 0 and energy[-1] < 0
    return dict(
        case, ids=IDS, stage_temperatures_K=fit.x[:n].tolist(),
        condenser_duty_W=float(energy[0]), reboiler_duty_W=float(-energy[-1]),
        outlets=[dict(molar_flow_mol_s=f, z=[float(x[j]), float(1-x[j])])
                 for j, f in zip(stage_numbers, flows)],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    watched = [PATH, ADAPTER, Path(ref.caloric.__file__), Path(ref.equilibrium.__file__),
               Path(ref.ps.__file__), ref.caloric.SOURCE]
    result = dict(
        schema="a17-rigorous-distillation-reference-v1",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "cantera", "numpy", "scipy")},
        sha256={str(path.relative_to(ref.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        tolerances=TOLERANCES, cases=[reference(case) for case in CASES],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
