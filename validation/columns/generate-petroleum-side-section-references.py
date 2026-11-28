"""Independent MESH references for petroleum-fractionator pump-arounds and
side strippers (#190), written to `petroleum-v2.json`.

The crude-like feed, its characterisation, the property model (teqp 0.23.2
Peng-Robinson with the A4 adapter's root selection, Kesler-Lee / NASA7 ideal
gas) and the feed flash are those of `generate-petroleum-references.py`,
imported unchanged. This script adds its own MESH formulation for a network
of cascades (Seader, Henley and Roper, *Separation Process Principles*, 3rd
ed., Eqs. 10-1 to 10-4 for every stage):

- the main column of `generate-petroleum-references.py` (total or partial
  condenser, partial reboiler or an adiabatic steam-stripped bottom tray);
- **pump-arounds**: a fixed molar flow P of liquid leaves main stage i at the
  stage liquid composition and enthalpy and re-enters stage j <= i as a feed
  of the same composition, with the enthalpy of the liquid root at the return
  temperature and the return-stage pressure, or h_L,i - Q/P for a duty Q;
- **side strippers**: equilibrium cascades at the pressure of their draw
  stage. The liquid draw S (an unknown) leaves main stage i and enters the
  stripper's top stage with the main stage's liquid enthalpy; the stripper's
  top vapor enters main stage r <= i as a feed with its own enthalpy; its
  bottom liquid is the side product (molar flow, or standard liquid volume
  flow with the product rate as an unknown); it is stripped by steam to its
  adiabatic bottom stage or by a partial reboiler with a duty or boilup-ratio
  specification.

Unknowns are, for every stage of every cascade, the temperature, the log
liquid and vapor mole fractions with explicit summation equations and the
log liquid and vapor flows, plus the main condenser/reboiler duties, the free
distillate, each stripper's log draw rate and (boilup specification) its
reboiler duty, and the log rate of every volume-specified product. Total flows
are not eliminated and the specifications are residual equations; this is a
different arrangement from the Rust engine, which eliminates the summations
and total flows and substitutes the volume specifications in an outer loop.

Start: the Wang-Henke bubble-point iteration of the reference script for the
main column with each stripper replaced by its product drawn at the draw
stage and its steam fed to the return stage (no pump-arounds), stripper
stages at the draw-stage state and hand-stated flow guesses of the case
list; then the damped Newton / Levenberg-Marquardt solve of the reference
script. Every fitted stage phase of every cascade and each returned
pump-around liquid at its return temperature pass the A4 adapter's
multi-start TPD search, and every stage has distinct liquid and vapor roots.
No Rust code or output is used.

    python validation/columns/generate-petroleum-side-section-references.py --output /tmp/petroleum-v2-candidate.json
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

PATH = Path(__file__).resolve()
ROOT = PATH.parents[2]
BASE_PATH = ROOT / "validation/columns/generate-petroleum-references.py"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = load("petroleum_reference", BASE_PATH)
# Budgets fixed before any Rust comparison: those of petroleum-v1.json.
TOLERANCES = base.TOLERANCES
# Duty unknowns are carried in units of 100 kW, so that the reference
# Newton step limit (10 per unknown) does not throttle the duty changes a
# pump-around or side stripper needs.
DUTY_SCALE = 1.0e5


def fractionator(case):
    print(case["id"], file=sys.stderr, flush=True)
    components, kij, x_cuts = base.system(case["water"])
    case = dict(case, ids=[c["id"] for c in components], kij=kij)
    props = base.Properties(components, kij)
    n, c = case["internal_stages"] + 2, len(components)
    p = base.pressures(case)
    steam = case["steam"] is not None
    partial = case["condenser"] == "partial_vapor"
    feeds = [dict(f) for f in case["feeds"]]
    if steam:
        feeds.append(dict(case["steam"], stage=n))
    flows_in = np.zeros((n, c))
    h_in = np.zeros(n)
    for f in feeds:
        z = np.asarray(f["z"], float)
        stage = f["stage"] - 1
        assert f["pressure_Pa"] >= p[stage]
        if f is feeds[-1] and steam:
            h = props.single_phase_h(f["temperature_K"], f["pressure_Pa"], z, "vapor")
            f["vapor_fraction"] = 1.0
        else:
            state = props.flash(f["temperature_K"], f["pressure_Pa"], z)
            h = state["enthalpy"]
            f["vapor_fraction"] = state["vapor_fraction"]
        f["enthalpy_J_per_mol"] = h
        flows_in[stage] += f["molar_flow_mol_s"]*z
        h_in[stage] += f["molar_flow_mol_s"]*h
    total = flows_in.sum()
    volume = np.array([comp.get("standard_liquid_volume_m3_per_mol", np.nan) for comp in components])
    draws = case["draws"]
    pumps = case["pump_arounds"]
    reflux_ratio = case["top"]["reflux_ratio"]
    guess = case["guess"]
    is_volume = ["standard_liquid_volume_flow_m3_s" in d for d in draws]
    plain = [k for k, d in enumerate(draws) if "stripper" not in d]
    stripped = [k for k, d in enumerate(draws) if "stripper" in d]
    free_plain = [k for k in plain if is_volume[k]]

    strippers = []
    for k in stripped:
        s = dict(draws[k]["stripper"], product=k, draw=draws[k]["stage"] - 1)
        s["ret"] = s["return_stage"] - 1
        assert 1 <= s["ret"] <= s["draw"]
        s["pressure"] = p[s["draw"]]
        if "steam" in s:
            st = s["steam"]
            assert st["pressure_Pa"] >= s["pressure"]
            z = np.asarray(st["z"], float)
            st["enthalpy_J_per_mol"] = props.single_phase_h(st["temperature_K"], st["pressure_Pa"],
                                                            z, "vapor")
            st["vapor_fraction"] = 1.0
            s["steam_flow"] = st["molar_flow_mol_s"]
            s["steam_z"] = z
        else:
            s["steam_flow"] = 0.0
        strippers.append(s)

    # Variable layout.
    v_stages = list(range(0 if partial else 1, n))
    sizes = dict(t=n, u=n*c, w=n*c, ln_l=n, ln_v=len(v_stages), q_c=1,
                 q_r=0 if steam else 1, ln_d=1 if steam and not partial else 0,
                 ln_u=len(free_plain))
    for i, s in enumerate(strippers):
        ns = s["stages"]
        sizes.update({f"s{i}_t": ns, f"s{i}_u": ns*c, f"s{i}_w": ns*c, f"s{i}_l": ns,
                      f"s{i}_v": ns, f"s{i}_ln_s": 1,
                      f"s{i}_q": 1 if "boilup_ratio" in s else 0,
                      f"s{i}_ln_b": 1 if is_volume[s["product"]] else 0})
    offsets, at = {}, 0
    for name, size in sizes.items():
        offsets[name] = slice(at, at + size)
        at += size
    # Continuation factor on the pump-around circulation (see the solve).
    ramp = [1.0]

    def unpack(v):
        t = v[offsets["t"]]
        x = np.exp(v[offsets["u"]].reshape(n, c))
        y = np.exp(v[offsets["w"]].reshape(n, c))
        liquid = np.exp(v[offsets["ln_l"]])
        vapor = np.zeros(n)
        vapor[v_stages] = np.exp(v[offsets["ln_v"]])
        q = np.zeros(n)
        q[0] = DUTY_SCALE*v[offsets["q_c"]][0]
        if not steam:
            q[-1] = DUTY_SCALE*v[offsets["q_r"]][0]
        liquid_draw = np.zeros(n)
        for k in plain:
            if not is_volume[k]:
                liquid_draw[draws[k]["stage"] - 1] += draws[k]["molar_flow_mol_s"]
        for j, k in enumerate(free_plain):
            liquid_draw[draws[k]["stage"] - 1] += np.exp(v[offsets["ln_u"]][j])
        for pump in pumps:
            liquid_draw[pump["draw_stage"] - 1] += ramp[0]*pump["molar_flow_mol_s"]
        side = []
        for i, s in enumerate(strippers):
            ns = s["stages"]
            st = dict(t=v[offsets[f"s{i}_t"]], x=np.exp(v[offsets[f"s{i}_u"]].reshape(ns, c)),
                      y=np.exp(v[offsets[f"s{i}_w"]].reshape(ns, c)),
                      liquid=np.exp(v[offsets[f"s{i}_l"]]), vapor=np.exp(v[offsets[f"s{i}_v"]]),
                      draw=np.exp(v[offsets[f"s{i}_ln_s"]][0]))
            if "boilup_ratio" in s:
                st["q"] = DUTY_SCALE*v[offsets[f"s{i}_q"]][0]
            elif "reboiler_duty_W" in s:
                st["q"] = s["reboiler_duty_W"]
            else:
                st["q"] = 0.0
            st["bottoms"] = (np.exp(v[offsets[f"s{i}_ln_b"]][0]) if is_volume[s["product"]]
                             else draws[s["product"]]["molar_flow_mol_s"])
            liquid_draw[s["draw"]] += st["draw"]
            side.append(st)
        if partial:
            distillate = vapor[0]
        else:
            distillate = (np.exp(v[offsets["ln_d"]][0]) if steam
                          else case["bottom"]["distillate_fraction"]*total)
            liquid_draw[0] += distillate
        return t, x, y, liquid, vapor, q, liquid_draw, distillate, side

    def pump_return_h(pump, t, x):
        i, j = pump["draw_stage"] - 1, pump["return_stage"] - 1
        if "return_temperature_K" in pump:
            return props.phase(pump["return_temperature_K"], p[j], x[i], "liquid")[1]
        return props.phase(t[i], p[i], x[i], "liquid")[1] - pump["duty_W"]/pump["molar_flow_mol_s"]

    def residual(v):
        t, x, y, liquid, vapor, q, liquid_draw, distillate, side = unpack(v)
        comp_extra = np.zeros((n, c))
        h_extra = np.zeros(n)
        for pump in pumps:
            i, j = pump["draw_stage"] - 1, pump["return_stage"] - 1
            comp_extra[j] += ramp[0]*pump["molar_flow_mol_s"]*x[i]
            h_extra[j] += ramp[0]*pump["molar_flow_mol_s"]*pump_return_h(pump, t, x)
        for s, st in zip(strippers, side):
            comp_extra[s["ret"]] += st["vapor"][0]*st["y"][0]
            h_extra[s["ret"]] += st["vapor"][0]*props.phase(st["t"][0], s["pressure"], st["y"][0],
                                                            "vapor")[1]
        out = []
        for j in range(n):
            ln_l, h_l, _ = props.phase(t[j], p[j], x[j], "liquid")
            ln_v, h_v, _ = props.phase(t[j], p[j], y[j], "vapor")
            comp_in = flows_in[j] + comp_extra[j]
            energy = h_in[j] + h_extra[j] + q[j]
            if j > 0:
                comp_in = comp_in + liquid[j - 1]*x[j - 1]
                energy += liquid[j - 1]*props.phase(t[j - 1], p[j - 1], x[j - 1], "liquid")[1]
            if j + 1 < n:
                comp_in = comp_in + vapor[j + 1]*y[j + 1]
                energy += vapor[j + 1]*props.phase(t[j + 1], p[j + 1], y[j + 1], "vapor")[1]
            comp_out = (liquid[j] + liquid_draw[j])*x[j] + vapor[j]*y[j]
            energy -= (liquid[j] + liquid_draw[j])*h_l + vapor[j]*h_v
            out.extend((comp_in - comp_out)/total)
            out.extend(np.log(y[j]) - np.log(x[j]) - (ln_l - ln_v))
            out.append(x[j].sum() - 1)
            out.append(y[j].sum() - 1)
            out.append(energy/total/1e4)
        out.append((liquid[0] - reflux_ratio*distillate)/total)
        if partial and not steam:
            out.append((vapor[0] - case["bottom"]["distillate_fraction"]*total)/total)
        for k in free_plain:
            j = draws[k]["stage"] - 1
            out.append(np.exp(v[offsets["ln_u"]][free_plain.index(k)])*np.dot(x[j]/x[j].sum(), volume)
                       / draws[k]["standard_liquid_volume_flow_m3_s"] - 1)
        for s, st in zip(strippers, side):
            i, ns, ps = s["draw"], s["stages"], s["pressure"]
            ts, xs, ys, ls, vs = st["t"], st["x"], st["y"], st["liquid"], st["vapor"]
            h_draw = props.phase(t[i], p[i], x[i], "liquid")[1]
            for j in range(ns):
                ln_l, h_l, _ = props.phase(ts[j], ps, xs[j], "liquid")
                ln_v, h_v, _ = props.phase(ts[j], ps, ys[j], "vapor")
                if j == 0:
                    comp_in = st["draw"]*x[i]
                    energy = st["draw"]*h_draw
                else:
                    comp_in = ls[j - 1]*xs[j - 1]
                    energy = ls[j - 1]*props.phase(ts[j - 1], ps, xs[j - 1], "liquid")[1]
                if j + 1 < ns:
                    comp_in = comp_in + vs[j + 1]*ys[j + 1]
                    energy += vs[j + 1]*props.phase(ts[j + 1], ps, ys[j + 1], "vapor")[1]
                else:
                    if "steam" in s:
                        comp_in = comp_in + s["steam_flow"]*s["steam_z"]
                        energy += s["steam_flow"]*s["steam"]["enthalpy_J_per_mol"]
                    energy += st["q"]
                comp_out = ls[j]*xs[j] + vs[j]*ys[j]
                energy -= ls[j]*h_l + vs[j]*h_v
                out.extend((comp_in - comp_out)/total)
                out.extend(np.log(ys[j]) - np.log(xs[j]) - (ln_l - ln_v))
                out.append(xs[j].sum() - 1)
                out.append(ys[j].sum() - 1)
                out.append(energy/total/1e4)
            out.append((ls[-1] - st["bottoms"])/total)
            if "boilup_ratio" in s:
                out.append((vs[-1] - s["boilup_ratio"]*st["bottoms"])/total)
            if is_volume[s["product"]]:
                xb = xs[-1]/xs[-1].sum()
                out.append(st["bottoms"]*np.dot(xb, volume)
                           / draws[s["product"]]["standard_liquid_volume_flow_m3_s"] - 1)
        return np.array(out)

    # Start: the Wang-Henke iteration for the main column with every product
    # drawn at its stage (stripped products at their guessed rates), stripper
    # steam fed to the return stage, and no pump-arounds.
    d0 = guess["distillate_mol_s"]
    feed_stage = case["feeds"][0]["stage"] - 1
    rates = [d.get("molar_flow_mol_s") or guess["draws_mol_s"][k] for k, d in enumerate(draws)]
    u0 = np.zeros(n)
    for k, d in enumerate(draws):
        u0[d["stage"] - 1] += rates[k]
    start_flows, start_h = flows_in.copy(), h_in.copy()
    for s in strippers:
        if "steam" in s:
            start_flows[s["ret"]] += s["steam_flow"]*s["steam_z"]
            start_h[s["ret"]] += s["steam_flow"]*s["steam"]["enthalpy_J_per_mol"]
    t0, x0, y0, l0, vv0, q0 = base.bubble_point_method(
        props, p, start_flows, start_h, u0, d0, reflux_ratio, partial,
        np.linspace(guess["top_K"], guess["bottom_K"], n),
        (feed_stage, case["steam"]["molar_flow_mol_s"]) if steam else None)
    v0 = list(t0) + list(np.log(x0).ravel()) + list(np.log(y0).ravel()) + list(np.log(l0))
    v0 += list(np.log(vv0[v_stages]))
    v0.append(q0[0]/DUTY_SCALE)
    if not steam:
        v0.append(q0[-1]/DUTY_SCALE)
    if steam and not partial:
        v0.append(np.log(d0))
    for k in free_plain:
        v0.append(np.log(rates[k]))
    for s in strippers:
        i, ns, b = s["draw"], s["stages"], rates[s["product"]]
        vg = s["vapor_guess_mol_s"]
        v0 += [t0[i]]*ns + list(np.log(x0[i]))*ns + list(np.log(y0[i]))*ns
        v0 += [np.log(b - s["steam_flow"] + vg)]*(ns - 1) + [np.log(b)]
        v0 += [np.log(vg)]*ns
        v0.append(np.log(b + vg - s["steam_flow"]))
        if "boilup_ratio" in s:
            v0.append(3.0e4*s["boilup_ratio"]*b/DUTY_SCALE)
        if is_volume[s["product"]]:
            v0.append(np.log(b))
    assert len(v0) == at
    # Pump-arounds enter by continuation in their circulation (at a fixed
    # specific cooling), each step a converged solve from the previous one.
    v = np.array(v0, float)
    for step in ([0.0, 0.25, 0.5, 0.75, 1.0] if pumps else [1.0]):
        ramp[0] = step
        v = base.solve(residual, v)
    assert ramp[0] == 1.0
    t, x, y, liquid, vapor, q, liquid_draw, distillate, side = unpack(v)
    x = x/x.sum(axis=1, keepdims=True)
    y = y/y.sum(axis=1, keepdims=True)

    def check_stage(t_, p_, xl, yv, label):
        base.check_stable(props, t_, p_, xl)
        base.check_stable(props, t_, p_, yv)
        z_l = props.phase(t_, p_, xl, "liquid")[2]
        z_v = props.phase(t_, p_, yv, "vapor")[2]
        assert z_v > 1.05*z_l, (case["id"], label, z_l, z_v)

    for j in range(n):
        if j > 0 or partial:
            check_stage(t[j], p[j], x[j], y[j], j)
    pump_duties = []
    for pump in pumps:
        i = pump["draw_stage"] - 1
        h_ret = pump_return_h(pump, t, x)
        if "return_temperature_K" in pump:
            assert pump["return_temperature_K"] < t[i]
            base.check_stable(props, pump["return_temperature_K"], p[pump["return_stage"] - 1], x[i])
        pump_duties.append(float(pump["molar_flow_mol_s"]
                                 * (props.phase(t[i], p[i], x[i], "liquid")[1] - h_ret)))
    stripper_results = []
    products = {}
    for s, st in zip(strippers, side):
        xs = st["x"]/st["x"].sum(axis=1, keepdims=True)
        ys = st["y"]/st["y"].sum(axis=1, keepdims=True)
        for j in range(s["stages"]):
            check_stage(st["t"][j], s["pressure"], xs[j], ys[j], ("stripper", s["product"], j))
        stripper_results.append(dict(
            product=s["product"] + 1, stage_temperatures_K=st["t"].tolist(),
            draw_mol_s=float(st["draw"]), vapor_return_mol_s=float(st["vapor"][0]),
            bottoms_mol_s=float(st["bottoms"]), reboiler_duty_W=float(st["q"])))
        products[s["product"]] = (float(st["bottoms"]), xs[-1])
    outlets = [dict(phase="vapor" if partial else "liquid", molar_flow_mol_s=float(distillate),
                    z=(y[0] if partial else x[0]).tolist()),
               dict(phase="liquid", molar_flow_mol_s=float(liquid[-1]), z=x[-1].tolist())]
    for k, d in enumerate(draws):
        if k in products:
            flow, z = products[k]
        else:
            j = d["stage"] - 1
            flow = (d["molar_flow_mol_s"] if not is_volume[k]
                    else float(np.exp(v[offsets["ln_u"]][free_plain.index(k)])))
            z = x[j]
        outlets.append(dict(phase="liquid", molar_flow_mol_s=float(flow), z=z.tolist(),
                            standard_liquid_volume_flow_m3_s=(
                                float(flow*np.dot(z, volume)) if not case["water"] else None)))
    for s in strippers:
        if "steam" in s:
            s["steam"]["z"] = list(s["steam"]["z"])
    return dict(case, components=components, cut_mole_fractions=x_cuts,
                feeds=feeds[:len(case["feeds"])], steam=feeds[-1] if steam else None,
                stage_temperatures_K=t.tolist(), stage_pressures_Pa=p.tolist(),
                condenser_duty_W=float(-q[0]), reboiler_duty_W=float(q[-1]),
                distillate_mol_s=float(distillate), reflux_mol_s=float(liquid[0]),
                bottoms_mol_s=float(liquid[-1]), boilup_mol_s=float(vapor[-1]),
                pump_around_duties_W=pump_duties, side_strippers=stripper_results,
                outlets=outlets)


def cases():
    z = base.hydrocarbon_feed()
    steam = dict(temperature_K=520.0, pressure_Pa=3.0e5, z=[1.0] + [0.0]*base.CUTS)
    return [
        base.case(id="crude-pump-arounds-reboiled-PR",
                  feeds=[base.inlet(480.0, 2.0e5, 100.0, z, 9)],
                  draws=[dict(stage=4, molar_flow_mol_s=15.0), dict(stage=7, molar_flow_mol_s=15.0)],
                  pump_arounds=[dict(draw_stage=7, return_stage=5, molar_flow_mol_s=30.0,
                                     return_temperature_K=400.0),
                                dict(draw_stage=4, return_stage=2, molar_flow_mol_s=20.0,
                                     duty_W=3.0e5)],
                  top=dict(reflux_ratio=2.0), bottom=dict(distillate_fraction=0.3),
                  guess=dict(top_K=390.0, bottom_K=480.0, distillate_mol_s=30.0)),
        base.case(id="crude-reboiled-side-strippers-volume-PR",
                  feeds=[base.inlet(470.0, 2.0e5, 100.0, z, 9)],
                  draws=[dict(stage=4, molar_flow_mol_s=15.0,
                              stripper=dict(stages=2, return_stage=3, reboiler_duty_W=1.5e5,
                                            vapor_guess_mol_s=4.0)),
                         dict(stage=7, standard_liquid_volume_flow_m3_s=3.0e-3,
                              stripper=dict(stages=3, return_stage=6, boilup_ratio=0.5,
                                            vapor_guess_mol_s=8.0))],
                  pump_arounds=[],
                  top=dict(reflux_ratio=3.0), bottom=dict(distillate_fraction=0.3),
                  guess=dict(top_K=390.0, bottom_K=500.0, distillate_mol_s=30.0,
                             draws_mol_s=[15.0, 16.0])),
        base.case(id="crude-steam-side-stripper-pump-around-PR", water=True,
                  kij_note=base.WATER_NOTE, condenser="partial_vapor",
                  feeds=[base.inlet(480.0, 2.0e5, 100.0, base.with_water(z), 9)],
                  steam=dict(steam, molar_flow_mol_s=5.0),
                  draws=[dict(stage=4, molar_flow_mol_s=15.0,
                              stripper=dict(stages=3, return_stage=3,
                                            steam=dict(steam, molar_flow_mol_s=2.0),
                                            vapor_guess_mol_s=5.0)),
                         dict(stage=7, molar_flow_mol_s=15.0)],
                  pump_arounds=[dict(draw_stage=7, return_stage=6, molar_flow_mol_s=25.0,
                                     duty_W=2.0e5)],
                  top=dict(reflux_ratio=2.0), bottom=None,
                  guess=dict(top_K=395.0, bottom_K=455.0, distillate_mol_s=34.0)),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[],
                        help="generate only these case ids (review runs)")
    args = parser.parse_args()
    np.seterr(all="raise")
    watched = [PATH, BASE_PATH, ROOT / "validation/generate-flash-references.py",
               base.PSEUDO_PATH, ROOT / "data/source-extracts/components.provisional.json"]
    result = dict(
        schema="a17-petroleum-fractionator-side-section-reference-v2",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "numpy", "scipy")},
        sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        water_density_60F_kg_m3=base.WATER_DENSITY_60F,
        tolerances=TOLERANCES,
        cases=[fractionator(c) for c in cases() if c["id"] in args.case or not args.case],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
