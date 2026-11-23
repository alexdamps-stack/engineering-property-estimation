"""Independent binary countercurrent absorber and stripper references (#187).

SciPy solves the balances in temperature, composition and flow coordinates.
teqp and Cantera supply the cubic fugacities and phase enthalpies through
the same reference adapters used for the rigorous distillation cases.
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
# Allow propagation from the phase-property reference budgets through the
# calculated outlet flows and stage temperatures.
TOLERANCES = {"temperature_K": [2e-5, 0], "composition": [2e-7, 0],
              "molar_flow_mol_s": [2e-5, 0]}


def feed(temperature, pressure, flow, ethane):
    return dict(temperature_K=temperature, pressure_Pa=pressure,
                molar_flow_mol_s=flow, z=[ethane, 1-ethane])


CASES = [
    dict(id="propane-absorption-PR", variant="PR", pressure_Pa=2e6, stages=4,
         gas=feed(315, 2.2e6, 70, 0.65), liquid=feed(260, 2e6, 55, 0.18)),
    dict(id="ethane-stripping-packed-PR-BM", variant="PR-BM", pressure_Pa=2e6, stages=3,
         gas=feed(345, 2e6, 60, 0.15), liquid=feed(245, 2.3e6, 55, 0.65),
         packed_height_m=1.2, hetp_m=0.4),
    dict(id="single-contact-PR", variant="PR", pressure_Pa=2e6, stages=1,
         gas=feed(315, 2.2e6, 70, 0.65), liquid=feed(260, 2e6, 55, 0.18)),
]


def reference(case):
    n, p, variant = case["stages"], case["pressure_Pa"], case["variant"]
    feed_ethane, feed_enthalpy = np.zeros((2, n))
    for port, stage, limit in [("liquid", 0, "bubble"), ("gas", n-1, "dew")]:
        inlet = case[port]
        f, z = inlet["molar_flow_mol_s"], inlet["z"]
        t, pin = inlet["temperature_K"], inlet["pressure_Pa"]
        flashed = ref.ps.tp_state(port, IDS, z, t, pin, variant)
        saturation = ref.sat_point(IDS, z, pin, variant, KIJ, limit)["T"]
        assert flashed["phases"] == 1, (case["id"], port, flashed)
        assert (t < saturation if port == "liquid" else t > saturation), (case["id"], port, saturation)
        feed_ethane[stage] += f*z[0]
        feed_enthalpy[stage] += f*flashed["target_h"]
    liquid_feed = case["liquid"]["molar_flow_mol_s"]
    gas_feed = case["gas"]["molar_flow_mol_s"]
    flow = liquid_feed + gas_feed

    def residual(values):
        temperature, x, y, vapor = values.reshape(4, n)
        liquid = np.r_[vapor[1:] + liquid_feed - vapor[0], flow - vapor[0]]
        equilibrium = np.empty((n, 2))
        h_l, h_v = np.empty((2, n))
        for j, t in enumerate(temperature):
            xl, yv = np.array([x[j], 1-x[j]]), np.array([y[j], 1-y[j]])
            model, _ = ref.equilibrium.model_at(IDS, t, variant, KIJ)
            ln_k = ref.ln_phi(model, t, p, xl, "liquid") - ref.ln_phi(model, t, p, yv, "vapor")
            equilibrium[j] = xl*np.exp(ln_k)-yv
            h_l[j] = ref.phase_h(IDS, xl, t, p, variant, KIJ, "liquid")
            h_v[j] = ref.phase_h(IDS, yv, t, p, variant, KIJ, "vapor")
        component = (np.r_[0, (liquid*x)[:-1]] + np.r_[(vapor*y)[1:], 0]
                     + feed_ethane - liquid*x - vapor*y)
        energy = (np.r_[0, (liquid*h_l)[:-1]] + np.r_[(vapor*h_v)[1:], 0]
                  + feed_enthalpy - liquid*h_l - vapor*h_v)
        return np.r_[equilibrium.ravel(), component/flow, energy/(flow*10_000)]

    z = feed_ethane.sum()/flow
    bubble = ref.sat_point(IDS, [z, 1-z], p, variant, KIJ, "bubble")
    initial = np.r_[np.full(n, bubble["T"]), np.full(n, z),
                    np.full(n, bubble["vapor"][0]), np.full(n, gas_feed)]
    lo = np.r_[np.full(n, 200), np.full(2*n, 1e-9), np.full(n, 1e-8)]
    hi = np.r_[np.full(n, 350), np.full(2*n, 1-1e-9), flow-1e-8, np.full(n-1, np.inf)]
    fit = least_squares(residual, initial, bounds=(lo, hi), x_scale="jac",
                        ftol=1e-13, gtol=1e-13, xtol=1e-13)
    assert fit.success and np.max(np.abs(fit.fun)) < 1e-10, (case["id"], fit.message, fit.fun)
    temperature, x, y, vapor = fit.x.reshape(4, n)
    assert np.all(vapor[1:] + liquid_feed - vapor[0] > 0)
    for t, xl, yv in zip(temperature, x, y):
        midpoint = (xl+yv)/2
        tie = ref.equilibrium.flash_case("contact", IDS, [midpoint, 1-midpoint], t, p, variant)
        assert tie["phases"] == 2 and abs(tie["liquid"][0]-xl) < 1e-8
        assert abs(tie["vapor"][0]-yv) < 1e-8
    top_flow = float(vapor[0])
    return dict(case, ids=IDS, stage_temperatures_K=temperature.tolist(),
                outlets=[dict(molar_flow_mol_s=f, z=[float(z), float(1-z)])
                         for f, z in [(top_flow, y[0]), (flow-top_flow, x[-1])]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    watched = [PATH, ADAPTER, Path(ref.caloric.__file__), Path(ref.equilibrium.__file__),
               Path(ref.ps.__file__), ref.caloric.SOURCE]
    result = dict(
        schema="a17-countercurrent-absorber-reference-v1",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "cantera", "numpy", "scipy")},
        sha256={str(path.relative_to(ref.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        tolerances=TOLERANCES, cases=[reference(case) for case in CASES],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
