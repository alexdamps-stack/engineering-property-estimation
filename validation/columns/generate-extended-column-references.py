"""Independent MESH references outside the original column domain (#186, #187).

The first column references (`distillation-v1.json`, `absorber-v1.json`) were
restricted to dry ethane/propane at 2 MPa and 245-345 K. These cases use
components and states that the columns accept since they take their domain
from the property package (A13 extended bubble/dew domain, #227):

- PR methanol/water at 101.325 kPa with an explicit kij (stage temperatures
  338-368 K, water and methanol active);
- PR propane/methanol at 2 MPa, with the reboiler near 430 K;
- a PR methane/CO2 absorber with a methanol solvent at 3 MPa.

A damped Newton iteration (forward-difference Jacobian, NumPy least-squares
steps, backtracking on the residual norm; SciPy Levenberg-Marquardt as a
fallback) solves the full multicomponent component, enthalpy, logarithmic
equilibrium and explicit summation equations in temperature, log mole
fraction and log vapor-flow coordinates. Distillation starts from a
smoothed sharp split with every stage at its own bubble point; the absorber
from a linear temperature profile. teqp gives the cubic fugacity
coefficients and residual properties and Cantera NASA7 the ideal-gas
enthalpy, through the reviewed reference adapters used by the other column
references. The independent TPD flash adapter classifies every inlet as
single-phase. At every fitted stage both phases must pass the adapter's
finite-start TPD search (which excludes a second liquid at the reference
solution), and the liquid and vapor roots must be the distinct
minimum-Gibbs roots of their compositions. No Rust code or output is used.
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
# The adapter's bubble/dew scan is limited to the old 200-350 K column domain.
# These cases use the cubic package interval instead.
T_BOUNDS = (200.0, 600.0)
ref.T_RANGE = T_BOUNDS

# Budgets fixed before any Rust comparison ([absolute, relative]); the same
# propagation argument and values as distillation-v1.json / absorber-v1.json.
TOLERANCES = {"temperature_K": [2e-5, 0], "composition": [2e-7, 0],
              "molar_flow_mol_s": [2e-5, 0], "duty_W": [0.02, 2e-6]}

# An explicit approximation, not a fitted parameter: PR with kij = 0 predicts a
# liquid-liquid split of methanol/water at the column bottom, so these cases use
# a constant negative kij of the magnitude reported for PR methanol/water VLE.
MEOH_WATER_KIJ = -0.075
MEOH_WATER_NOTE = ("methanol/water kij = -0.075 is an explicit approximation for "
                   "implementation agreement, not a fitted parameter")
ZERO_NOTE = "zero kij: an explicit approximation for implementation agreement"


def inlet(temperature, pressure, flow, z, stage=None):
    out = dict(temperature_K=temperature, pressure_Pa=pressure, molar_flow_mol_s=flow, z=z)
    if stage is not None:
        out["stage"] = stage
    return out


DISTILLATION = [
    dict(id="methanol-water-trays-PR", kind="distillation", variant="PR",
         ids=["methanol", "water"], kij=[[0, MEOH_WATER_KIJ], [MEOH_WATER_KIJ, 0]],
         kij_note=MEOH_WATER_NOTE, pressure_Pa=101325.0, internal_stages=6, reflux_ratio=1.5, distillate_fraction=0.5,
         feeds=[inlet(340.0, 1.2e5, 100.0, [0.5, 0.5], 5)], draws=[]),
    dict(id="propane-methanol-trays-PR", kind="distillation", variant="PR",
         ids=["propane", "methanol"], kij=[[0, 0], [0, 0]], kij_note=ZERO_NOTE,
         pressure_Pa=2e6, internal_stages=4, reflux_ratio=1.5, distillate_fraction=0.5,
         feeds=[inlet(300.0, 2.2e6, 100.0, [0.5, 0.5], 4)], draws=[]),
]

ABSORBER = [
    dict(id="methanol-solvent-absorber-PR", kind="absorber", variant="PR",
         ids=["methane", "carbon-dioxide", "methanol"], kij=np.zeros((3, 3)).tolist(),
         kij_note=ZERO_NOTE, pressure_Pa=3e6, stages=5,
         gas=inlet(290.0, 3.1e6, 50.0, [0.85, 0.15, 0.0]),
         liquid=inlet(250.0, 3.2e6, 40.0, [0.0, 0.0, 1.0])),
]


def stage_properties(case, t, p, x, y):
    ids, variant, kij = case["ids"], case["variant"], case["kij"]
    model, _ = ref.equilibrium.model_at(ids, t, variant, kij)
    ln_k = ref.ln_phi(model, t, p, x, "liquid") - ref.ln_phi(model, t, p, y, "vapor")
    return (np.log(x) + ln_k - np.log(y),
            ref.phase_h(ids, x, t, p, variant, kij, "liquid"),
            ref.phase_h(ids, y, t, p, variant, kij, "vapor"))


def feed_enthalpy(case, stream, phase):
    """Enthalpy of a single-phase inlet, classified by the independent TPD flash."""
    ids, variant, kij = case["ids"], case["variant"], case["kij"]
    t, p, z = stream["temperature_K"], stream["pressure_Pa"], stream["z"]
    flashed = ref.equilibrium.flash_case("inlet", ids, z, t, p, variant, kij)
    assert flashed["phases"] == 1, (case["id"], stream, flashed)
    if phase == "liquid":
        active = [i for i, v in enumerate(z) if v > 0]
        sub = [z[i] for i in active]
        if len(active) > 1:
            bubble = ref.sat_point([ids[i] for i in active], sub, p, variant,
                                   [[kij[i][k] for k in active] for i in active], "bubble")
            assert t < bubble["T"], (case["id"], "liquid inlet is not subcooled")
    return ref.phase_h(ids, np.asarray(z, float), t, p, variant, kij, phase)


def check_tie_lines(case, temperatures, p, x, y):
    """Each fitted stage is a stable vapor-liquid equilibrium.

    The MESH residuals already equate the fugacities. In addition, as the
    flash adapter requires of its own split phases, each phase must pass the
    finite-start TPD search (no trial phase, including a second liquid, lowers
    the Gibbs energy), the liquid and vapor roots used in the fit must be the
    minimum-Gibbs roots of their compositions, and they must be distinct."""
    equilibrium = ref.equilibrium
    floor = equilibrium.TOLERANCES["stable_tm_floor"]
    for t, xl, yv in zip(temperatures, x, y):
        model, _ = equilibrium.model_at(case["ids"], t, case["variant"], case["kij"])
        for phase in (xl, yv):
            tm = equilibrium.tpd_search(model, case["ids"], t, p, phase)[0]
            assert tm >= floor, (case["id"], t, "stage phase unstable", tm)
        z_liquid = equilibrium.roots(model, t, p, xl)[0][0]
        z_vapor = equilibrium.roots(model, t, p, yv)[-1][0]
        assert equilibrium.min_gibbs(model, t, p, xl)[0] == z_liquid, (case["id"], t)
        assert equilibrium.min_gibbs(model, t, p, yv)[0] == z_vapor, (case["id"], t)
        assert z_vapor > 1.05*z_liquid, (case["id"], t, z_liquid, z_vapor)


def solve(case, residual, initial):
    """Damped Newton on the square system in unbounded coordinates.

    Mole fractions are exp(u) with explicit summation residuals and vapor
    flows exp(w) above their positivity floor, so no bound is active. The
    Jacobian is a forward difference; each step is halved until the residual
    norm decreases. A trial point at which the reference adapters cannot
    evaluate a phase gets a large residual and is rejected. If Newton stalls,
    SciPy's Levenberg-Marquardt continues from the best point."""
    def guarded(u):
        try:
            r = residual(u)
        except (AssertionError, ValueError, RuntimeError, IndexError):
            return np.full(len(u), 1e3)
        return r if np.all(np.isfinite(r)) else np.full(len(u), 1e3)

    u, r = np.array(initial, float), guarded(initial)
    for _ in range(60):
        if np.max(np.abs(r)) < 1e-13:
            break
        jacobian = np.empty((len(r), len(u)))
        for k in range(len(u)):
            h = 1e-7*max(1.0, abs(u[k]))
            shifted = u.copy()
            shifted[k] += h
            jacobian[:, k] = (guarded(shifted) - r)/h
        step = np.linalg.lstsq(jacobian, -r, rcond=None)[0]
        for halving in range(30):
            trial = u + 0.5**halving*step
            rt = guarded(trial)
            if np.linalg.norm(rt) < np.linalg.norm(r):
                u, r = trial, rt
                break
        else:
            break
    if np.max(np.abs(r)) >= 1e-11:
        fit = least_squares(guarded, u, method="lm", x_scale="jac",
                            ftol=1e-15, xtol=1e-15, gtol=1e-15, max_nfev=200*len(u))
        u, r = fit.x, fit.fun
    worst = np.max(np.abs(r))
    assert worst < 1e-10, (case["id"], worst)
    return u


def distillation(case):
    n, c = case["internal_stages"] + 2, len(case["ids"])
    p = case["pressure_Pa"]
    feed_totals, feed_h = np.zeros((2, n))
    feed_components = np.zeros((n, c))
    for stream in case["feeds"]:
        j, f = stream["stage"] - 1, stream["molar_flow_mol_s"]
        feed_totals[j] += f
        feed_components[j] += f*np.asarray(stream["z"])
        feed_h[j] += f*feed_enthalpy(case, stream, "liquid")
    flow = feed_totals.sum()
    distillate = flow*case["distillate_fraction"]
    reflux = distillate*case["reflux_ratio"]
    draws = np.zeros(n)
    draws[0] = distillate
    for draw in case["draws"]:
        draws[draw["stage"] - 1] += draw["molar_flow_mol_s"]
    net = np.cumsum(feed_totals - draws)

    def unpack(values):
        t = values[:n]
        y = values[n:n + n*c].reshape(n, c)
        lower = values[n + n*c:n + n*c + (n-1)*c].reshape(n-1, c)
        x = np.vstack([y[1], lower])
        vapor = np.r_[0, reflux + distillate, values[n + n*c + (n-1)*c:]]
        liquid = np.r_[reflux, vapor[2:] + net[1:-1], net[-1]]
        return t, x, y, vapor, liquid

    def balances(values):
        t, x, y, vapor, liquid = unpack(values)
        eq = np.empty((n, c))
        h_l, h_v = np.empty((2, n))
        for j in range(n):
            eq[j], h_l[j], h_v[j] = stage_properties(case, t[j], p, x[j], y[j])
        lx, vy = liquid[:, None]*x, vapor[:, None]*y
        component = (np.vstack([np.zeros(c), lx[:-1]]) + np.vstack([vy[1:], np.zeros(c)])
                     + feed_components - (liquid + draws)[:, None]*x - vy)
        energy = (np.r_[0, (liquid*h_l)[:-1]] + np.r_[(vapor*h_v)[1:], 0]
                  + feed_h - (liquid + draws)*h_l - vapor*h_v)
        return t, x, y, liquid, eq, component, energy

    def residual(values):
        _, x, y, _, eq, component, energy = balances(values)
        return np.r_[eq.ravel(), y.sum(1) - 1, x[1:].sum(1) - 1,
                     component[1:, :-1].ravel()/flow, energy[1:-1]/(flow*10_000)]

    z = feed_components.sum(0)/flow
    bubble = ref.sat_point(case["ids"], z, p, case["variant"], case["kij"], "bubble")
    floor = np.maximum(0, -net[1:-1])
    fractions = (2*n-1)*c

    def direct(u):
        return np.r_[u[:n], np.exp(u[n:n+fractions]), floor + np.exp(u[n+fractions:])]

    # Initial guess: a smoothed sharp split by feed K-values gives distillate
    # and bottoms estimates; the liquid profile is linear between them and
    # every stage starts at the bubble point of its liquid.
    sharp, remaining = np.zeros(c), distillate
    for i in np.argsort(-bubble["K"]):
        sharp[i] = min(remaining, flow*z[i])
        remaining -= sharp[i]
    top = 0.9*sharp/distillate + 0.1*z
    bottom = (flow*z - distillate*top)/(flow - distillate)

    def stage_bubble(s, lowest):
        # The adapter accepts the first bracket of its upward scan; the column
        # condenser is its coldest stage, so the other stages scan above it.
        ref.T_RANGE = (lowest, T_BOUNDS[1])
        try:
            return ref.sat_point(case["ids"], (1-s)*top + s*bottom, p, case["variant"],
                                 case["kij"], "bubble")
        finally:
            ref.T_RANGE = T_BOUNDS

    condenser = stage_bubble(0.0, T_BOUNDS[0])
    guesses = [condenser] + [stage_bubble(s, condenser["T"]) for s in np.linspace(0, 1, n)[1:]]
    initial = np.r_[[g["T"] for g in guesses], np.log([g["vapor"] for g in guesses]).ravel(),
                    np.log([g["liquid"] for g in guesses[1:]]).ravel(),
                    np.full(n-2, np.log(reflux + distillate))]
    solved = direct(solve(case, lambda u: residual(direct(u)), initial))
    t, x, y, liquid, _, component, energy = balances(solved)
    check_tie_lines(case, t, p, x, y)
    # The condenser composition and the overall balance close independently.
    assert np.max(np.abs(component[0]))/flow < 1e-9
    stage_numbers = [0, n-1] + [d["stage"]-1 for d in case["draws"]]
    flows = [distillate, float(liquid[-1])] + [d["molar_flow_mol_s"] for d in case["draws"]]
    assert energy[0] > 0 and energy[-1] < 0
    return dict(case, stage_temperatures_K=t.tolist(),
                condenser_duty_W=float(energy[0]), reboiler_duty_W=float(-energy[-1]),
                outlets=[dict(molar_flow_mol_s=f, z=x[j].tolist())
                         for j, f in zip(stage_numbers, flows)])


def absorber(case):
    n, c = case["stages"], len(case["ids"])
    p = case["pressure_Pa"]
    feed_components = np.zeros((n, c))
    feed_h = np.zeros(n)
    for port, stage, phase in [("liquid", 0, "liquid"), ("gas", n-1, "vapor")]:
        stream = case[port]
        f = stream["molar_flow_mol_s"]
        feed_components[stage] += f*np.asarray(stream["z"])
        feed_h[stage] += f*feed_enthalpy(case, stream, phase)
    liquid_feed = case["liquid"]["molar_flow_mol_s"]
    gas_feed = case["gas"]["molar_flow_mol_s"]
    flow = liquid_feed + gas_feed

    def balances(values):
        t = values[:n]
        x = values[n:n + n*c].reshape(n, c)
        y = values[n + n*c:n + 2*n*c].reshape(n, c)
        vapor = values[n + 2*n*c:]
        liquid = np.r_[vapor[1:] + liquid_feed - vapor[0], flow - vapor[0]]
        eq = np.empty((n, c))
        h_l, h_v = np.empty((2, n))
        for j in range(n):
            eq[j], h_l[j], h_v[j] = stage_properties(case, t[j], p, x[j], y[j])
        lx, vy = liquid[:, None]*x, vapor[:, None]*y
        component = (np.vstack([np.zeros(c), lx[:-1]]) + np.vstack([vy[1:], np.zeros(c)])
                     + feed_components - lx - vy)
        energy = (np.r_[0, (liquid*h_l)[:-1]] + np.r_[(vapor*h_v)[1:], 0]
                  + feed_h - liquid*h_l - vapor*h_v)
        return t, x, y, vapor, liquid, eq, component, energy

    def residual(values):
        _, x, y, _, _, eq, component, energy = balances(values)
        return np.r_[eq.ravel(), x.sum(1) - 1, y.sum(1) - 1,
                     component[:, :-1].ravel()/flow, energy/(flow*10_000)]

    # Seed: temperatures linear between the inlets, each phase close to its
    # inlet composition with a tenth of the blend, the inlet gas flow.
    z = feed_components.sum(0)/flow
    t_liquid, t_gas = case["liquid"]["temperature_K"], case["gas"]["temperature_K"]
    x0 = 0.9*np.asarray(case["liquid"]["z"]) + 0.1*z
    y0 = 0.9*np.asarray(case["gas"]["z"]) + 0.1*z

    def direct(u):
        return np.r_[u[:n], np.exp(u[n:])]

    initial = np.r_[np.linspace(t_liquid, t_gas, n), np.log(np.tile(x0, n)),
                    np.log(np.tile(y0, n)), np.full(n, np.log(gas_feed))]
    solved = direct(solve(case, lambda u: residual(direct(u)), initial))
    t, x, y, vapor, liquid, _, _, _ = balances(solved)
    assert np.all(liquid > 0) and np.all(vapor > 0)
    check_tie_lines(case, t, p, x, y)
    top = float(vapor[0])
    return dict(case, stage_temperatures_K=t.tolist(),
                outlets=[dict(molar_flow_mol_s=top, z=y[0].tolist()),
                         dict(molar_flow_mol_s=flow - top, z=x[-1].tolist())])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    watched = [PATH, ADAPTER, Path(ref.caloric.__file__), Path(ref.equilibrium.__file__),
               Path(ref.ps.__file__), ref.caloric.SOURCE]
    result = dict(
        schema="a17-extended-column-reference-v1",
        packages={name: importlib.metadata.version(name) for name in ("teqp", "cantera", "numpy", "scipy")},
        sha256={str(path.relative_to(ref.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in watched},
        tolerances=TOLERANCES,
        cases=[distillation(case) for case in DISTILLATION] + [absorber(case) for case in ABSORBER],
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
