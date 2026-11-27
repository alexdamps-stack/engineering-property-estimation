"""Independent MESH references for the petroleum fractionator (#190).

A crude-like feed is characterised exactly as documented in
`docs/pseudo-components-status.md`: a TBP curve is cut into equal-temperature
cuts at constant Watson K (Fahim, Al-Sahhaf and Elkilani, *Fundamentals of
Petroleum Refining*, 2010, sections 3.3.1-3.3.2), each cut's SG follows from
K, and Tc, Pc, omega, M and the ideal-gas cp come from the independent Python
correlations of `validation/generate-pseudo-component-references.py`. Cut
mole fractions follow n_i ~ V_i SG_i / M_i. The standard liquid volume of a
pseudo-component is M / (SG rho_w) with the water density at 60 degF and
1 atm stated below (IAPWS-95).

Formulation (Seader, Henley and Roper, *Separation Process Principles*, 3rd
ed., Eqs. 10-1 to 10-4): every stage is an equilibrium stage; stage 1 is the
condenser (total: bubble-point liquid of the vapor from stage 2; partial:
equilibrium stage whose vapor is the distillate) and stage N the partial
reboiler or, with stripping steam, an adiabatic bottom tray fed by the steam.
Unknowns are the stage temperatures, log liquid and vapor mole fractions with
explicit summation equations, log liquid and vapor flows, the condenser and
reboiler duties and, where a specification does not fix them, the distillate
and the molar flow of each side product specified by standard liquid volume.
This is a different arrangement from the Rust engine (normalized log ratios,
eliminated summations and total flows, successive substitution on the volume
specifications).

Properties: teqp 0.23.2 Peng-Robinson (exact Omega factors, Mathias-Copeman
[m, 0, 0] alpha, explicit kij) for fugacity coefficients and residual
enthalpy, with roots selected and refined by the reviewed A4 flash adapter
(`validation/generate-flash-references.py`). Ideal-gas enthalpy integrates
the Kesler-Lee cp of each cut from 298.15 K, and for water the lower NASA7
interval of `data/source-extracts/components.provisional.json`. The feed's
two-phase TP state comes from a Rachford-Rice / successive-substitution flash
written here; every fitted stage phase and both feed phases must pass the
adapter's multi-start TPD search (no second liquid), and every stage must have
distinct liquid and vapor roots. No Rust code or output is used.

Generate a candidate in the pinned reference environment and compare before
replacing `petroleum-v1.json`:

    python validation/columns/generate-petroleum-references.py --output /tmp/petroleum-candidate.json
"""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares
import teqp

PATH = Path(__file__).resolve()
ROOT = PATH.parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


equilibrium = load("equilibrium", ROOT / "validation/generate-flash-references.py")
# The pseudo-component reference imports CoolProp at module level for its
# surrogate data; only its plain-Python correlations are used here, so a stub
# module satisfies that import when CoolProp is not installed.
PSEUDO_PATH = ROOT / "validation/generate-pseudo-component-references.py"
try:
    import CoolProp  # noqa: F401
except ImportError:
    stub = types.ModuleType("CoolProp")
    stub.CoolProp = types.ModuleType("CoolProp.CoolProp")
    sys.modules["CoolProp"] = stub
    sys.modules["CoolProp.CoolProp"] = stub.CoolProp
pseudo = load("pseudo_reference", PSEUDO_PATH)

R = equilibrium.R
T0 = 298.15
# Water at 60 degF (288.7056 K) and 101.325 kPa, IAPWS-95 (CoolProp 7.2.0):
# the standard liquid volume basis M / (SG rho_w) of a pseudo-component.
WATER_DENSITY_60F = 999.0170824078306

# Budgets fixed before any Rust comparison ([absolute, relative]): the values
# and propagation argument of extended-v1.json / specs-v1.json. The standard
# liquid volume flow of a side product follows from its molar flow and
# composition budgets (2e-5 mol/s x 2e-4 m3/mol plus composition terms).
TOLERANCES = {"temperature_K": [2e-5, 0], "composition": [2e-7, 0],
              "molar_flow_mol_s": [2e-5, 0], "duty_W": [0.02, 2e-6],
              "standard_liquid_volume_flow_m3_s": [1e-8, 0]}

CRUDE_TBP = [(0.0, 340.0), (20.0, 400.0), (50.0, 460.0), (80.0, 510.0), (100.0, 550.0)]
WATSON_K = 11.8
CUTS = 8
ZERO_NOTE = ("kij = 0 between pseudo-components: an explicit approximation "
             "(Fahim et al. 2010, sec. 3.6.1), not fitted values")
WATER_KIJ = 0.5
WATER_NOTE = ("kij = 0 between pseudo-components and water/pseudo-component kij = 0.5: "
              "explicit approximations of the magnitude used for PR water/hydrocarbon pairs, "
              "for implementation agreement, not fitted values")
WATER = json.loads((ROOT / "data/source-extracts/components.provisional.json").read_text())
WATER = next(c for c in WATER["components"] if c["id"] == "water")


def cut_curve(tbp, count, watson_k, prefix):
    """Equal-temperature cuts of a TBP curve (vol%, K) at constant Watson K."""
    def volume_at(t):
        for (v0, t0), (v1, t1) in zip(tbp, tbp[1:]):
            if t <= t1:
                return v0 + (v1 - v0)*(t - t0)/(t1 - t0)
        (v0, t0), (v1, t1) = tbp[-2], tbp[-1]
        return v0 + (v1 - v0)*(t - t0)/(t1 - t0)
    first, last = tbp[0][1], tbp[-1][1]
    width = (last - first)/count
    span = tbp[-1][0] - tbp[0][0]
    cuts = []
    for i in range(count):
        lo = first + width*i
        hi = last if i + 1 == count else first + width*(i + 1)
        tb = 0.5*(lo + hi)
        sg = (1.8*tb)**(1.0/3.0)/watson_k
        cuts.append(dict(id=f"{prefix}{i + 1}", normal_boiling_point_K=tb, specific_gravity=sg,
                         liquid_volume_fraction=(volume_at(hi) - volume_at(lo))/span))
    return cuts


def characterise(cut):
    tb, sg = cut["normal_boiling_point_K"], cut["specific_gravity"]
    k = (1.8*tb)**(1.0/3.0)/sg
    tc = pseudo.kesler_lee_tc(tb, sg)
    pc = pseudo.kesler_lee_pc(tb, sg)
    omega = pseudo.lee_kesler_omega(tb, tc, pc, k)
    m = pseudo.riazi_daubert_m(tb, sg)
    return dict(cut, critical_temperature_K=tc, critical_pressure_Pa=pc, acentric_factor=omega,
                molar_mass_kg_per_mol=m, ideal_gas_cp_over_R=pseudo.kesler_lee_cp_over_r(k, omega, m),
                standard_liquid_volume_m3_per_mol=m/(sg*WATER_DENSITY_60F))


def system(water):
    """Components (water first when present), constants, kij and cut fractions."""
    cuts = [characterise(c) for c in cut_curve(CRUDE_TBP, CUTS, WATSON_K, "pseudo:crude-")]
    moles = np.array([c["liquid_volume_fraction"]*c["specific_gravity"]/c["molar_mass_kg_per_mol"]
                      for c in cuts])
    x_cuts = (moles/moles.sum()).tolist()
    components = ([dict(id="water", **WATER["constants"])] if water else []) + cuts
    n = len(components)
    kij = np.zeros((n, n))
    if water:
        kij[0, 1:] = kij[1:, 0] = WATER_KIJ
    return components, kij.tolist(), x_cuts


class Properties:
    """teqp PR with the case constants; ideal gas from Kesler-Lee / NASA7."""

    def __init__(self, components, kij):
        self.components = components
        alphas = []
        for c in components:
            w = c["acentric_factor"]
            alphas.append({"type": "Mathias-Copeman", "c": [.37464 + 1.54226*w - .26992*w*w, 0, 0]})
        spec = {"type": "PR", "Tcrit / K": [c["critical_temperature_K"] for c in components],
                "pcrit / Pa": [c["critical_pressure_Pa"] for c in components],
                "acentric": [c["acentric_factor"] for c in components],
                "alpha": alphas, "kmat": kij}
        self.model = teqp.make_model({"kind": "cubic", "model": spec})
        meta = self.model.get_meta()
        assert abs(meta["OmegaA"] - equilibrium.OA) < 1e-16
        assert abs(meta["OmegaB"] - equilibrium.OB) < 1e-16
        self.cache = {}

    def ideal_h(self, t):
        out = []
        for c in self.components:
            if "ideal_gas_cp_over_R" in c:
                a = c["ideal_gas_cp_over_R"]
                prim = lambda x: a[0]*x + a[1]*x**2/2 + a[2]*x**3/3
            else:
                a = WATER["ideal_gas"]["coefficients"][0]
                assert t < WATER["ideal_gas"]["temperature_breaks_K"][1]
                prim = lambda x: a[0]*x + a[1]*x**2/2 + a[2]*x**3/3 + a[3]*x**4/4 + a[4]*x**5/5
            out.append(R*(prim(t) - prim(T0)))
        return np.array(out)

    def wilson(self, t, p):
        return np.array([c["critical_pressure_Pa"]/p*np.exp(5.373*(1 + c["acentric_factor"])
                                                          * (1 - c["critical_temperature_K"]/t))
                         for c in self.components])

    def phase(self, t, p, x, kind):
        """(ln phi, molar enthalpy, Z) on the smallest (liquid) or largest (vapor) root."""
        x = np.asarray(x, float)
        x = x/x.sum()
        key = (t, p, x.tobytes(), kind)
        if key not in self.cache:
            roots = equilibrium.roots(self.model, t, p, x)
            z, rho, ln_phi = roots[0] if kind == "liquid" else roots[-1]
            h_res = R*t*(self.model.get_Ar10(t, rho, x) + self.model.get_Ar01(t, rho, x))
            if len(self.cache) > 200000:
                self.cache.clear()
            self.cache[key] = (ln_phi, float(np.dot(x, self.ideal_h(t)) + h_res), z)
        return self.cache[key]

    def flash(self, t, p, z):
        """Two-phase TP split by successive substitution and Rachford-Rice."""
        z = np.asarray(z, float)
        k = self.wilson(t, p)
        for _ in range(500):
            def rr(beta):
                return float(np.sum(z*(k - 1)/(1 + beta*(k - 1))))
            if rr(0.0) <= 0 or rr(1.0) >= 0:
                raise ValueError("feed is not two-phase")
            beta = brentq(rr, 0.0, 1.0, xtol=1e-15)
            x = z/(1 + beta*(k - 1))
            y = k*x
            x, y = x/x.sum(), y/y.sum()
            ln_k = self.phase(t, p, x, "liquid")[0] - self.phase(t, p, y, "vapor")[0]
            if np.max(np.abs(ln_k - np.log(k))) < 1e-13:
                break
            k = np.exp(ln_k)
        else:
            raise RuntimeError("feed flash did not converge")
        for phase in (x, y):
            check_stable(self, t, p, phase)
        h = beta*self.phase(t, p, y, "vapor")[1] + (1 - beta)*self.phase(t, p, x, "liquid")[1]
        return dict(vapor_fraction=beta, enthalpy=h)

    def single_phase_h(self, t, p, z, kind):
        z = np.asarray(z, float)
        check_stable(self, t, p, z)
        return self.phase(t, p, z, kind)[1]


def check_stable(props, t, p, x):
    """Adapter TPD search with this system's Wilson starts; no trial phase lowers G."""
    equilibrium.wilson = lambda ids, t_, p_: props.wilson(t_, p_)
    tm = equilibrium.tpd_search(props.model, None, t, p, np.asarray(x, float))[0]
    assert tm >= equilibrium.TOLERANCES["stable_tm_floor"], (t, p, tm)


def inlet(temperature, pressure, flow, z, stage):
    return dict(temperature_K=temperature, pressure_Pa=pressure, molar_flow_mol_s=flow, z=z,
                stage=stage)


def pressures(case):
    n = case["internal_stages"] + 2
    top = case["pressure_Pa"]
    first = top + case["condenser_pressure_drop_Pa"]
    return np.array([top] + [first + (j - 1)*case["stage_pressure_drop_Pa"] for j in range(1, n)])


def fractionator(case):
    print(case["id"], file=sys.stderr, flush=True)
    components, kij, x_cuts = system(case["water"])
    case = dict(case, ids=[c["id"] for c in components], kij=kij)
    props = Properties(components, kij)
    n, c = case["internal_stages"] + 2, len(components)
    p = pressures(case)
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
    draw_stage = [d["stage"] - 1 for d in draws]
    free_draws = [k for k, d in enumerate(draws) if "standard_liquid_volume_flow_m3_s" in d]
    reflux_ratio = case["top"]["reflux_ratio"]

    # Variable layout.
    v_stages = list(range(0 if partial else 1, n))
    sizes = dict(t=n, u=n*c, w=n*c, ln_l=n, ln_v=len(v_stages), q_c=1,
                 q_r=0 if steam else 1, ln_d=1 if steam and not partial else 0,
                 ln_u=len(free_draws))
    offsets, at = {}, 0
    for name, size in sizes.items():
        offsets[name] = slice(at, at + size)
        at += size

    def unpack(v):
        t = v[offsets["t"]]
        x = np.exp(v[offsets["u"]].reshape(n, c))
        y = np.exp(v[offsets["w"]].reshape(n, c))
        liquid = np.exp(v[offsets["ln_l"]])
        vapor = np.zeros(n)
        vapor[v_stages] = np.exp(v[offsets["ln_v"]])
        q = np.zeros(n)
        q[0] = v[offsets["q_c"]][0]
        if not steam:
            q[-1] = v[offsets["q_r"]][0]
        liquid_draw = np.zeros(n)
        for k, d in enumerate(draws):
            liquid_draw[draw_stage[k]] += d.get("molar_flow_mol_s", 0.0)
        for j, k in enumerate(free_draws):
            liquid_draw[draw_stage[k]] += np.exp(v[offsets["ln_u"]][j])
        if partial:
            distillate = vapor[0]
        else:
            distillate = (np.exp(v[offsets["ln_d"]][0]) if steam
                          else case["bottom"]["distillate_fraction"]*total)
            liquid_draw[0] += distillate
        return t, x, y, liquid, vapor, q, liquid_draw, distillate

    def residual(v):
        t, x, y, liquid, vapor, q, liquid_draw, distillate = unpack(v)
        out = []
        for j in range(n):
            ln_l, h_l, _ = props.phase(t[j], p[j], x[j], "liquid")
            ln_v, h_v, _ = props.phase(t[j], p[j], y[j], "vapor")
            comp_in = flows_in[j].copy()
            energy = h_in[j] + q[j]
            if j > 0:
                comp_in += liquid[j - 1]*x[j - 1]
                energy += liquid[j - 1]*props.phase(t[j - 1], p[j - 1], x[j - 1], "liquid")[1]
            if j + 1 < n:
                comp_in += vapor[j + 1]*y[j + 1]
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
        for k in free_draws:
            xk = x[draw_stage[k]]/x[draw_stage[k]].sum()
            out.append(liquid_draw[draw_stage[k]]*np.dot(xk, volume)
                       / draws[k]["standard_liquid_volume_flow_m3_s"] - 1)
        return np.array(out)

    # Start: the bubble-point (Wang-Henke) method of Seader, Henley and Roper
    # (sec. 10.3) at the hand-stated distillate and side-product guesses of the
    # case, then the simultaneous Newton solve above.
    guess = case["guess"]
    d0 = guess["distillate_mol_s"]
    feed_stage = case["feeds"][0]["stage"] - 1
    u0 = np.zeros(n)
    for k, d in enumerate(draws):
        u0[draw_stage[k]] += d.get("molar_flow_mol_s", 0.0)
    for k in free_draws:
        u0[draw_stage[k]] += guess["draws_mol_s"][k]
    t0, x0, y0, l0, vv0, q0 = bubble_point_method(
        props, p, flows_in, h_in, u0, d0, reflux_ratio, partial,
        np.linspace(guess["top_K"], guess["bottom_K"], n),
        (feed_stage, case["steam"]["molar_flow_mol_s"]) if steam else None)
    v0 = list(t0) + list(np.log(x0).ravel()) + list(np.log(y0).ravel()) + list(np.log(l0))
    v0 += list(np.log(vv0[v_stages]))
    v0.append(q0[0])
    if not steam:
        v0.append(q0[-1])
    if steam and not partial:
        v0.append(np.log(d0))
    for k in free_draws:
        v0.append(np.log(guess["draws_mol_s"][k]))
    v = solve(residual, np.array(v0, float))
    t, x, y, liquid, vapor, q, liquid_draw, distillate = unpack(v)
    x = x/x.sum(axis=1, keepdims=True)
    y = y/y.sum(axis=1, keepdims=True)
    for j in range(n):
        if j > 0 or partial:
            check_stable(props, t[j], p[j], x[j])
            check_stable(props, t[j], p[j], y[j])
        z_l = props.phase(t[j], p[j], x[j], "liquid")[2]
        z_v = props.phase(t[j], p[j], y[j], "vapor")[2]
        assert z_v > 1.05*z_l, (case["id"], j, z_l, z_v)
    outlets = [dict(phase="vapor" if partial else "liquid", molar_flow_mol_s=float(distillate),
                    z=(y[0] if partial else x[0]).tolist()),
               dict(phase="liquid", molar_flow_mol_s=float(liquid[-1]), z=x[-1].tolist())]
    for k, d in enumerate(draws):
        flow = float(liquid_draw[draw_stage[k]])
        outlets.append(dict(phase="liquid", molar_flow_mol_s=flow, z=x[draw_stage[k]].tolist(),
                            standard_liquid_volume_flow_m3_s=(
                                float(flow*np.dot(x[draw_stage[k]], volume))
                                if not case["water"] else None)))
    return dict(case, components=components, cut_mole_fractions=x_cuts,
                feeds=feeds[:len(case["feeds"])], steam=feeds[-1] if steam else None,
                stage_temperatures_K=t.tolist(), stage_pressures_Pa=p.tolist(),
                condenser_duty_W=float(-q[0]), reboiler_duty_W=float(q[-1]),
                distillate_mol_s=float(distillate), reflux_mol_s=float(liquid[0]),
                bottoms_mol_s=float(liquid[-1]), boilup_mol_s=float(vapor[-1]), outlets=outlets)


def bubble_t(props, p, x, t):
    """Bubble temperature of liquid x at p (Newton on ln sum K x), and its vapor."""
    y = x.copy()
    for _ in range(100):
        def f(tt):
            nonlocal y
            for _ in range(5):
                k = np.exp(props.phase(tt, p, x, "liquid")[0] - props.phase(tt, p, y, "vapor")[0])
                y = k*x/np.sum(k*x)
            return np.log(np.sum(k*x)), k
        g, k = f(t)
        if abs(g) < 1e-12:
            return t, k*x/np.sum(k*x)
        dg = (f(t + 1e-4)[0] - g)/1e-4
        t = t - float(np.clip(g/dg, -15.0, 15.0))
    raise RuntimeError("bubble point did not converge")


def bubble_point_method(props, p, feed, h_feed, draws, distillate, reflux, partial, t, stripping):
    """Wang-Henke iteration at fixed distillate and draws (start values only).

    With stripping steam (`stripping` = (feed stage index, steam flow)) the
    vapor below the feed is held at 1.6 times the steam flow (constant molar
    overflow); top-down energy updates are ill-conditioned there."""
    n, c = feed.shape
    u = draws.copy()
    v1 = distillate if partial else 0.0
    if not partial:
        u[0] += distillate
    net = np.cumsum(feed.sum(axis=1) - u) - v1
    vapor = np.full(n, (reflux + 1)*distillate)
    vapor[0] = v1
    if stripping:
        vapor[stripping[0] + 1:] = 1.6*stripping[1]
    z = feed.sum(axis=0)/feed.sum()
    x = np.tile(z, (n, 1))
    y = x.copy()
    for _ in range(200):
        liquid = np.append(vapor[1:], 0.0) + net
        liquid[0] = reflux*distillate
        k = np.array([np.exp(props.phase(t[j], p[j], x[j], "liquid")[0]
                             - props.phase(t[j], p[j], y[j], "vapor")[0]) for j in range(n)])
        new_x = np.empty((n, c))
        for i in range(c):
            a = np.zeros((n, n))
            for j in range(n):
                a[j, j] = -(liquid[j] + u[j] + vapor[j]*k[j, i])
                if j > 0:
                    a[j, j - 1] = liquid[j - 1]
                if j + 1 < n:
                    a[j, j + 1] = vapor[j + 1]*k[j + 1, i]
            new_x[:, i] = np.linalg.solve(a, -feed[:, i])
        new_x = np.maximum(new_x, 1e-30)
        x = new_x/new_x.sum(axis=1, keepdims=True)
        new_t = np.empty(n)
        for j in range(n):
            new_t[j], y[j] = bubble_t(props, p[j], x[j], t[j])
        h_l = np.array([props.phase(new_t[j], p[j], x[j], "liquid")[1] for j in range(n)])
        h_v = np.array([props.phase(new_t[j], p[j], y[j], "vapor")[1] for j in range(n)])
        new_v = vapor.copy()
        new_v[1] = liquid[0] + u[0] + v1
        last = stripping[0] if stripping else n - 1
        for j in range(1, last):
            above = liquid[0] if j == 1 else new_v[j] + net[j - 1]
            known = (net[j] + u[j])*h_l[j] + new_v[j]*h_v[j] - h_feed[j] - above*h_l[j - 1]
            # Keep the start's flows positive while the profile is far off.
            new_v[j + 1] = max(known/(h_v[j + 1] - h_l[j]), 1e-3*feed.sum(),
                               1e-3*feed.sum() - net[j + 1] if j + 1 < n - 1 else 0.0)
        change = max(np.max(np.abs(new_t - t)), np.max(np.abs(new_v - vapor))/feed.sum())
        t, vapor = new_t, 0.5*(vapor + new_v)
        vapor[0] = v1
        if change < 1e-6:
            break
    else:
        print(f"  bubble-point start stopped at change {change:.3e}", file=sys.stderr, flush=True)
    liquid = np.append(vapor[1:], 0.0) + net
    liquid[0] = reflux*distillate
    h_l = np.array([props.phase(t[j], p[j], x[j], "liquid")[1] for j in range(n)])
    h_v = np.array([props.phase(t[j], p[j], y[j], "vapor")[1] for j in range(n)])
    q = np.zeros(n)
    q[0] = (liquid[0] + u[0])*h_l[0] + v1*h_v[0] - vapor[1]*h_v[1]
    q[-1] = (liquid[-1]*h_l[-1] + vapor[-1]*h_v[-1] - h_feed[-1]
             - liquid[-2]*h_l[-2])
    return t, x, y, liquid, vapor, q


def solve(residual, v):
    """Damped Newton (forward-difference Jacobian, backtracking on the max norm),
    with SciPy Levenberg-Marquardt as a fallback start."""
    def jacobian(v, r):
        jac = np.empty((len(r), len(v)))
        for k in range(len(v)):
            h = 1e-7*max(abs(v[k]), 1.0)
            shifted = v.copy()
            shifted[k] += h
            jac[:, k] = (residual(shifted) - r)/h
        return jac
    for attempt in range(2):
        r = residual(v)
        for _ in range(80):
            error = np.max(np.abs(r))
            print(f"  newton residual {error:.3e} at {int(np.argmax(np.abs(r)))} of {len(r)}",
                  file=sys.stderr, flush=True)
            if error < 1e-12:
                return v
            step = np.linalg.solve(jacobian(v, r), -r)
            scale = min(1.0, 0.2/max(np.max(np.abs(step)), 1e-300)*50)
            for _ in range(40):
                trial = v + scale*step
                try:
                    rt = residual(trial)
                except (AssertionError, ValueError, FloatingPointError):
                    scale *= 0.5
                    continue
                if np.all(np.isfinite(rt)) and np.linalg.norm(rt) < np.linalg.norm(r):
                    v, r = trial, rt
                    break
                scale *= 0.5
            else:
                break
        if attempt == 0:
            v = least_squares(residual, v, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                              max_nfev=200*len(v)).x
    raise RuntimeError(f"column did not converge (residual {np.max(np.abs(residual(v))):.3e})")


def hydrocarbon_feed():
    _, _, x = system(False)
    return x


def with_water(z):
    return [0.0] + list(z)


def case(**kwargs):
    base = dict(variant="PR", water=False, steam=None, condenser="total", pressure_Pa=1.5e5,
                condenser_pressure_drop_Pa=0.0, stage_pressure_drop_Pa=0.0, internal_stages=10,
                tbp_curve_K=CRUDE_TBP, watson_k=WATSON_K, cuts=CUTS, kij_note=ZERO_NOTE)
    base.update(kwargs)
    return base


def cases():
    z = hydrocarbon_feed()
    feed = inlet(480.0, 2.0e5, 100.0, z, 9)
    return [
        case(id="crude-two-draws-reboiled-PR", feeds=[feed],
             draws=[dict(stage=4, molar_flow_mol_s=15.0), dict(stage=7, molar_flow_mol_s=15.0)],
             top=dict(reflux_ratio=2.0), bottom=dict(distillate_fraction=0.3),
             guess=dict(top_K=390.0, bottom_K=480.0, distillate_mol_s=30.0)),
        case(id="crude-volume-draws-partial-pressure-drop-PR", condenser="partial_vapor",
             condenser_pressure_drop_Pa=1.0e4, stage_pressure_drop_Pa=700.0,
             feeds=[inlet(490.0, 2.5e5, 100.0, z, 8)],
             draws=[dict(stage=4, standard_liquid_volume_flow_m3_s=2.5e-3),
                    dict(stage=6, molar_flow_mol_s=10.0),
                    dict(stage=10, standard_liquid_volume_flow_m3_s=2.0e-3)],
             top=dict(reflux_ratio=3.0), bottom=dict(distillate_fraction=0.25),
             guess=dict(top_K=400.0, bottom_K=490.0, distillate_mol_s=25.0, draws_mol_s=[17.0, 10.0, 11.0])),
        case(id="crude-steam-stripped-PR", water=True, kij_note=WATER_NOTE, condenser="partial_vapor",
             feeds=[inlet(480.0, 2.0e5, 100.0, with_water(z), 9)],
             steam=dict(temperature_K=520.0, pressure_Pa=3.0e5, molar_flow_mol_s=5.0,
                        z=[1.0] + [0.0]*CUTS),
             draws=[dict(stage=4, molar_flow_mol_s=15.0), dict(stage=7, molar_flow_mol_s=15.0)],
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
    watched = [PATH, ROOT / "validation/generate-flash-references.py", PSEUDO_PATH,
               ROOT / "data/source-extracts/components.provisional.json"]
    result = dict(
        schema="a17-petroleum-fractionator-reference-v1",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "numpy", "scipy")},
        sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        water_density_60F_kg_m3=WATER_DENSITY_60F,
        tolerances=TOLERANCES,
        cases=[fractionator(c) for c in cases() if c["id"] in args.case or not args.case],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
