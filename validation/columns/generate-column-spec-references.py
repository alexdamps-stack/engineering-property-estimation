"""Independent MESH references for the extended #186 specification set.

Covers what the first distillation references (`distillation-v1.json`,
`extended-v1.json`) do not: a partial condenser (vapor distillate, and vapor
plus a fixed liquid distillate), a pressure profile (condenser drop and a
per-stage drop), a vapor side draw, and specifications other than the reflux
ratio and distillate fraction (reboiler duty, condenser duty, boilup ratio,
distillate mole fraction, bottoms recovery).

Formulation (Seader, Henley and Roper, Separation Process Principles, 3rd ed.,
Eqs. 10-1 to 10-4): every stage is an equilibrium stage at its own pressure;
the condenser is stage 1 and the partial reboiler stage N. A total condenser
returns liquid of the composition of the vapor from stage 2; a partial
condenser is an equilibrium stage whose vapor is the vapor distillate. The
unknowns are temperatures, log mole fractions with explicit summation
equations, log vapor flows above their positivity floor and, where a
specification does not fix them, the log distillate and log reflux flows.
Each specification that does not fix a flow adds one residual equation. The
damped Newton / Levenberg-Marquardt solver, the stage property evaluation
(teqp fugacities and residual properties, Cantera NASA7 ideal-gas enthalpy),
the TPD feed classification and the stage tie-line stability checks are
imported from `generate-extended-column-references.py`, which in turn uses
the reviewed reference adapters. Initial flows for free distillate or reflux
are stated in each case by hand. No Rust code or output is used.
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path

import numpy as np

PATH = Path(__file__).resolve()
EXTENDED = PATH.with_name("generate-extended-column-references.py")
spec = importlib.util.spec_from_file_location("extended_reference", EXTENDED)
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)
ref = ext.ref

# Budgets fixed before any Rust comparison ([absolute, relative]); the same
# propagation argument and values as extended-v1.json.
TOLERANCES = ext.TOLERANCES

C2C3 = ["ethane", "propane"]
ZERO = [[0, 0], [0, 0]]
MEOH_WATER = [[0, ext.MEOH_WATER_KIJ], [ext.MEOH_WATER_KIJ, 0]]


def column(**kwargs):
    base = dict(kind="distillation", variant="PR", ids=C2C3, kij=ZERO,
                kij_note=ext.ZERO_NOTE, condenser="total", liquid_distillate_mol_s=0.0,
                condenser_pressure_drop_Pa=0.0, stage_pressure_drop_Pa=0.0,
                draws=[], vapor_draws=[])
    base.update(kwargs)
    return base


FEED = ext.inlet(295.0, 2.3e6, 100.0, [0.45, 0.55], 4)

CASES = [
    column(id="partial-vapor-PR", pressure_Pa=2e6, internal_stages=4, feeds=[FEED],
           condenser="partial_vapor",
           top=dict(reflux_ratio=2.5), bottom=dict(distillate_fraction=0.45)),
    column(id="partial-vapor-liquid-vapor-draw-PR-BM", variant="PR-BM", pressure_Pa=2e6,
           internal_stages=5, condenser="partial_vapor_liquid", liquid_distillate_mol_s=20.0,
           feeds=[ext.inlet(290.0, 2.3e6, 100.0, [0.5, 0.5], 4)],
           vapor_draws=[dict(stage=5, molar_flow_mol_s=8.0)],
           top=dict(reflux_ratio=3.0), bottom=dict(distillate_fraction=0.45)),
    column(id="pressure-drop-methanol-water-PR", ids=["methanol", "water"], kij=MEOH_WATER,
           kij_note=ext.MEOH_WATER_NOTE, pressure_Pa=101325.0, condenser_pressure_drop_Pa=5000.0,
           stage_pressure_drop_Pa=700.0, internal_stages=6,
           feeds=[ext.inlet(340.0, 1.2e5, 100.0, [0.5, 0.5], 5)],
           top=dict(reflux_ratio=1.5), bottom=dict(distillate_fraction=0.5)),
    column(id="reboiler-duty-PR", pressure_Pa=2e6, internal_stages=4, feeds=[FEED],
           top=dict(reflux_ratio=2.5), bottom=dict(reboiler_duty_W=1.9e6),
           guess=dict(distillate=45.0)),
    column(id="purity-boilup-PR-BM", variant="PR-BM", pressure_Pa=2e6, internal_stages=6,
           feeds=[FEED],
           top=dict(distillate_mole_fraction=dict(component=0, value=0.95)),
           bottom=dict(boilup_ratio=2.0), guess=dict(distillate=42.0, reflux=110.0)),
    column(id="condenser-duty-recovery-PR", pressure_Pa=2e6, internal_stages=5,
           feeds=[FEED], top=dict(condenser_duty_W=1.6e6),
           bottom=dict(bottoms_recovery=dict(component=1, value=0.95)),
           guess=dict(distillate=45.0, reflux=100.0)),
]


def distillation(case):
    n, c = case["internal_stages"] + 2, len(case["ids"])
    first = case["pressure_Pa"] + case["condenser_pressure_drop_Pa"]
    pressures = np.r_[case["pressure_Pa"], first + case["stage_pressure_drop_Pa"]*np.arange(n-1)]
    partial = case["condenser"] != "total"
    liquid_distillate = case["liquid_distillate_mol_s"]
    feed_totals, feed_h = np.zeros((2, n))
    feed_components = np.zeros((n, c))
    for stream in case["feeds"]:
        j, f = stream["stage"] - 1, stream["molar_flow_mol_s"]
        assert stream["pressure_Pa"] >= pressures[j], (case["id"], "feed below stage pressure")
        feed_totals[j] += f
        feed_components[j] += f*np.asarray(stream["z"])
        feed_h[j] += f*ext.feed_enthalpy(case, stream, "liquid")
    flow = feed_totals.sum()
    liquid_draws, vapor_draws = np.zeros((2, n))
    for draw in case["draws"]:
        liquid_draws[draw["stage"] - 1] += draw["molar_flow_mol_s"]
    for draw in case["vapor_draws"]:
        vapor_draws[draw["stage"] - 1] += draw["molar_flow_mol_s"]
    available = flow - liquid_draws.sum() - vapor_draws.sum()

    top, bottom = case["top"], case["bottom"]
    fixed_distillate = None
    if "distillate_fraction" in bottom:
        fixed_distillate = flow*bottom["distillate_fraction"]
    elif "bottoms_rate_mol_s" in bottom:
        fixed_distillate = available - bottom["bottoms_rate_mol_s"]
    if "distillate_rate_mol_s" in top:
        assert fixed_distillate is None, case["id"]
        fixed_distillate = top["distillate_rate_mol_s"]
    free_reflux = not ({"reflux_ratio", "reflux_rate_mol_s"} & set(top))
    free = (fixed_distillate is None) + free_reflux
    guess = case.get("guess", {})

    def flows(extra):
        k = 0
        if fixed_distillate is None:
            distillate = extra[k]
            k += 1
        else:
            distillate = fixed_distillate
        if "reflux_ratio" in top:
            reflux = top["reflux_ratio"]*distillate
        elif "reflux_rate_mol_s" in top:
            reflux = top["reflux_rate_mol_s"]
        else:
            reflux = extra[k]
        return distillate, reflux

    rows_x = n if partial else n - 1

    def unpack(values):
        t = values[:n]
        k = n
        y = values[k:k + n*c].reshape(n, c)
        k += n*c
        lower = values[k:k + rows_x*c].reshape(rows_x, c)
        k += rows_x*c
        x = lower if partial else np.vstack([y[1], lower])
        inner = values[k:k + n - 2]
        distillate, reflux = flows(values[k + n - 2:])
        if distillate >= available or (partial and distillate <= liquid_distillate):
            raise ValueError("distillate outside its admissible range")
        draws = liquid_draws + vapor_draws
        draws[0] = distillate
        net = np.cumsum(feed_totals - draws)
        vapor = np.r_[distillate - liquid_distillate if partial else 0.0,
                      reflux + distillate, inner]
        liquid = np.r_[reflux, vapor[2:] + net[1:-1], net[-1]]
        if np.any(liquid <= 0) or np.any(vapor[1:] <= 0):
            raise ValueError("nonpositive internal flow")
        liquid_out = liquid + liquid_draws
        liquid_out[0] = reflux + (liquid_distillate if partial else distillate)
        vapor_out = vapor + vapor_draws
        return t, x, y, vapor, liquid, liquid_out, vapor_out, distillate, reflux

    def balances(values):
        t, x, y, vapor, liquid, liquid_out, vapor_out, distillate, reflux = unpack(values)
        eq = np.empty((n, c))
        h_l, h_v = np.empty((2, n))
        for j in range(n):
            eq[j], h_l[j], h_v[j] = ext.stage_properties(case, t[j], pressures[j], x[j], y[j])
        lx, vy = liquid[:, None]*x, vapor[:, None]*y
        component = (np.vstack([np.zeros(c), lx[:-1]]) + np.vstack([vy[1:], np.zeros(c)])
                     + feed_components - liquid_out[:, None]*x - vapor_out[:, None]*y)
        energy = (np.r_[0, (liquid*h_l)[:-1]] + np.r_[(vapor*h_v)[1:], 0]
                  + feed_h - liquid_out*h_l - vapor_out*h_v)
        top_flow = (vapor[0]*y[0] + (liquid_distillate*x[0] if partial else distillate*x[0]))
        x_distillate = top_flow/distillate
        return dict(t=t, x=x, y=y, vapor=vapor, liquid=liquid, eq=eq, component=component,
                    energy=energy, distillate=distillate, reflux=reflux,
                    x_distillate=x_distillate)

    def spec_residuals(s):
        out = []
        bottoms = s["liquid"][-1]
        scale = flow*10_000
        if "condenser_duty_W" in top:
            out.append((s["energy"][0] - top["condenser_duty_W"])/scale)
        if "distillate_mole_fraction" in top:
            target = top["distillate_mole_fraction"]
            out.append(s["x_distillate"][target["component"]] - target["value"])
        if "distillate_recovery" in top:
            target = top["distillate_recovery"]
            i = target["component"]
            out.append((s["distillate"]*s["x_distillate"][i]
                        - target["value"]*feed_components[:, i].sum())/flow)
        if "reboiler_duty_W" in bottom:
            out.append((-s["energy"][-1] - bottom["reboiler_duty_W"])/scale)
        if "boilup_ratio" in bottom:
            out.append((s["vapor"][-1] - bottom["boilup_ratio"]*bottoms)/flow)
        if "bottoms_mole_fraction" in bottom:
            target = bottom["bottoms_mole_fraction"]
            out.append(s["x"][-1][target["component"]] - target["value"])
        if "bottoms_recovery" in bottom:
            target = bottom["bottoms_recovery"]
            i = target["component"]
            out.append((bottoms*s["x"][-1][i] - target["value"]*feed_components[:, i].sum())/flow)
        assert len(out) == free, (case["id"], "specification count")
        return out

    def residual(values):
        s = balances(values)
        component = s["component"] if partial else s["component"][1:]
        x_rows = s["x"] if partial else s["x"][1:]
        return np.r_[s["eq"].ravel(), s["y"].sum(1) - 1, x_rows.sum(1) - 1,
                     component[:, :-1].ravel()/flow, s["energy"][1:-1]/(flow*10_000),
                     spec_residuals(s)]

    # Initial guess as in the extended references: a smoothed sharp split
    # at the starting distillate and every stage at its own bubble point.
    d0 = fixed_distillate if fixed_distillate is not None else guess["distillate"]
    if "reflux_ratio" in top:
        l0 = top["reflux_ratio"]*d0
    elif "reflux_rate_mol_s" in top:
        l0 = top["reflux_rate_mol_s"]
    else:
        l0 = guess["reflux"]
    z = feed_components.sum(0)/flow
    bubble = ref.sat_point(case["ids"], z, pressures[1], case["variant"], case["kij"], "bubble")
    sharp, remaining = np.zeros(c), d0
    for i in np.argsort(-bubble["K"]):
        sharp[i] = min(remaining, flow*z[i])
        remaining -= sharp[i]
    top_z = 0.9*sharp/d0 + 0.1*z
    bottom_z = (flow*z - d0*top_z)/(flow - d0)

    def stage_bubble(j, s, lowest):
        ref.T_RANGE = (lowest, ext.T_BOUNDS[1])
        try:
            return ref.sat_point(case["ids"], (1-s)*top_z + s*bottom_z, pressures[j],
                                 case["variant"], case["kij"], "bubble")
        finally:
            ref.T_RANGE = ext.T_BOUNDS

    condenser = stage_bubble(0, 0.0, ext.T_BOUNDS[0])
    guesses = [condenser] + [stage_bubble(j, s, condenser["T"])
                             for j, s in zip(range(1, n), np.linspace(0, 1, n)[1:])]
    liquids = guesses if partial else guesses[1:]
    draws0 = liquid_draws + vapor_draws
    draws0[0] = d0
    net0 = np.cumsum(feed_totals - draws0)
    floor = np.maximum(0, -net0[1:-1])
    initial = np.r_[[g["T"] for g in guesses], np.log([g["vapor"] for g in guesses]).ravel(),
                    np.log([g["liquid"] for g in liquids]).ravel(),
                    np.log(floor + l0 + d0)]
    if fixed_distillate is None:
        initial = np.r_[initial, np.log(d0)]
    if free_reflux:
        initial = np.r_[initial, np.log(l0)]
    def direct(u):
        return np.r_[u[:n], np.exp(u[n:])]

    solved = direct(ext.solve(case, lambda u: residual(direct(u)), initial))
    s = balances(solved)
    for j in range(n):
        ext.check_tie_lines(case, s["t"][j:j+1], pressures[j], s["x"][j:j+1], s["y"][j:j+1])
    # The overall component balance closes independently of the equations used.
    assert np.max(np.abs(s["component"].sum(0)))/flow < 1e-9, case["id"]
    energy = s["energy"]
    assert energy[0] > 0 and energy[-1] < 0, case["id"]
    bottoms = float(s["liquid"][-1])
    outlets = [dict(phase="vapor" if partial else "liquid",
                    molar_flow_mol_s=float(s["vapor"][0] if partial else s["distillate"]),
                    z=(s["y"][0] if partial else s["x"][0]).tolist()),
               dict(phase="liquid", molar_flow_mol_s=bottoms, z=s["x"][-1].tolist())]
    outlets += [dict(phase="liquid", molar_flow_mol_s=d["molar_flow_mol_s"],
                     z=s["x"][d["stage"]-1].tolist()) for d in case["draws"]]
    outlets += [dict(phase="vapor", molar_flow_mol_s=d["molar_flow_mol_s"],
                     z=s["y"][d["stage"]-1].tolist()) for d in case["vapor_draws"]]
    if partial and liquid_distillate > 0:
        outlets.append(dict(phase="liquid", molar_flow_mol_s=liquid_distillate,
                            z=s["x"][0].tolist()))
    return dict(case, stage_temperatures_K=s["t"].tolist(), stage_pressures_Pa=pressures.tolist(),
                condenser_duty_W=float(energy[0]), reboiler_duty_W=float(-energy[-1]),
                distillate_mol_s=float(s["distillate"]), reflux_mol_s=float(s["reflux"]),
                bottoms_mol_s=bottoms, boilup_mol_s=float(s["vapor"][-1]), outlets=outlets)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    watched = [PATH, EXTENDED, ext.ADAPTER, Path(ref.caloric.__file__),
               Path(ref.equilibrium.__file__), Path(ref.ps.__file__), ref.caloric.SOURCE]
    result = dict(
        schema="a17-column-specification-reference-v1",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "cantera", "numpy", "scipy")},
        sha256={str(path.relative_to(ref.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        tolerances=TOLERANCES,
        cases=[distillation(case) for case in CASES],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
