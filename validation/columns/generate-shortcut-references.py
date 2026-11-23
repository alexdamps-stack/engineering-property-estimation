"""Explicit A17 (#185) independent shortcut distillation references; no Rust calls.

Method (implemented exactly as the A17 #185 shortcut contract; see
validation/columns/README.md):

- Fenske, M. R. (1932). Fractionation of straight-run Pennsylvania gasoline.
  Ind. Eng. Chem. 24, 482-485. doi:10.1021/ie50269a003
  (minimum stages; Hengstebeck-Geddes form of the Fenske distribution of
  non-key components).
- Underwood, A. J. V. (1948). Fractional distillation of multicomponent
  mixtures. Chem. Eng. Prog. 44, 603-614 (minimum reflux, root theta between
  the key volatilities).
- Gilliland, E. R. (1940). Multicomponent rectification: estimation of the
  number of theoretical plates as a function of the reflux ratio.
  Ind. Eng. Chem. 32, 1220-1223.
- Molokanov, Yu. K., Korablina, T. P., Mazurina, N. I., Nikiforov, G. A. (1972).
  An approximation method for calculating the basic parameters of
  multicomponent fractionation. Int. Chem. Eng. 12, 209-212 (fit of the
  Gilliland correlation used here).
- Kirkbride, C. G. (1944). Process design procedure for multicomponent
  fractionators. Petrol. Refiner 23, 321-336 (feed-stage location).

Thermodynamics reuse the reviewed adapters, as the A14 PS generator does:
the A3 teqp AD + Cantera NASA7 caloric adapter (validation/generate-references.py)
for saturated-phase and feed enthalpies, the A4 SciPy/teqp equilibrium adapter
(validation/generate-flash-references.py) for root selection, fugacity
coefficients and the feed TP flash, and the A14 PS generator's phase-weighted TP
state (validation/rotating/generate-ps-references.py). Bubble and dew
temperatures at fixed pressure are solved here: an inner successive
substitution on ln K at fixed T (liquid = smallest, vapor = largest mechanically
stable volume root, as the contract specifies) inside an outer SciPy brentq on
T, bracketed by a 200..350 K scan with continuation. This is a different
algorithm from the Rust engine; only component parameters and the explicit
(zero) interaction matrix are shared.

Run manually in the pinned reference environment (validation/reference-requirements.txt;
in WSL2 Ubuntu on Windows):

    python validation/columns/generate-shortcut-references.py --output <candidate.json>

Review the candidate before replacing validation/columns/shortcut-v1.json.
"""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import math
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
ps = load("ps", "rotating/generate-ps-references.py")

# Budgets fixed before any Rust comparison ([absolute, relative]; a value passes when
# |rust - ref| <= abs + rel*|ref|). Regular double-precision implementation budgets of
# the same magnitude as ps-v1.json, propagated through the shortcut equations:
# - temperature_K, enthalpy_J_per_mol (feed TP state), composition: as ps-v1.json.
# - saturated_enthalpy_J_per_mol: the 2e-6 K temperature budget times a liquid cp of
#   up to ~250 J/mol/K gives 5e-4 J/mol on enthalpies evaluated at saturation points.
# - relative_volatility: ln phi budget 2e-9 per phase plus 2e-6 K times
#   |d ln alpha/dT| <~ 1e-2 /K -> 5e-8 relative.
# - stages (N_min, N, N_R, feed stage): ln(alpha_LK) >~ 0.3 amplifies the alpha
#   budget by <~ 3 in N_min; Gilliland/Kirkbride are smooth -> 1e-6 relative.
# - shortcut_scalar (q, theta, R_min, R, X, Y, Kirkbride ratio): q carries
#   5e-4 / |h_dew - h_bub| (>~ 5e3 J/mol) ~ 1e-7; Underwood adds the alpha budget
#   -> 1e-7 absolute plus 1e-6 relative.
# - molar_flow_mol_s: Fenske splits carry N_min ln(alpha_i) budgets (<~ 1e-6 rel).
# - duty_W: flows x saturated enthalpy differences (~1e4 J/mol) -> 1e-6 relative.
# - alpha_iterations: informational; a different but equivalent arrangement may
#   differ by one update.
TOLERANCES = {"temperature_K": [2e-6, 0.0], "enthalpy_J_per_mol": [2e-5, 2e-9],
              "saturated_enthalpy_J_per_mol": [5e-4, 2e-9], "composition": [2e-8, 1e-7],
              "relative_volatility": [0.0, 5e-8], "stages": [0.0, 1e-6],
              "shortcut_scalar": [1e-7, 1e-6], "molar_flow_mol_s": [1e-6, 1e-6],
              "duty_W": [1e-2, 1e-6], "alpha_iterations": [1, 0]}
T_RANGE = (200.0, 350.0)
ALPHA_TOL, ALPHA_MAX_ITER = 1e-10, 100


def ln_phi(model, t, p, x, phase):
    """ln phi of the smallest (liquid) or largest (vapor) mechanically stable root."""
    r = equilibrium.roots(model, t, p, np.asarray(x, dtype=float))
    return np.asarray(r[0][2] if phase == "liquid" else r[-1][2])


class Saturation:
    """Bubble (given liquid) or dew (given vapor) temperature at fixed pressure."""

    def __init__(self, ids, given, p, variant, kind, kij):
        self.ids, self.given, self.p, self.variant, self.kind, self.kij = ids, np.asarray(given, float), p, variant, kind, kij
        self.ln_k = None

    def inner(self, t):
        """Successive substitution on ln K at fixed T; returns the outer residual."""
        model, _ = equilibrium.model_at(self.ids, t, self.variant, self.kij)
        ln_k = np.log(equilibrium.wilson(self.ids, t, self.p)) if self.ln_k is None else self.ln_k.copy()
        g = self.given
        for _ in range(2000):
            k = np.exp(ln_k)
            if self.kind == "bubble":
                w = g*k/np.sum(g*k)
                new = ln_phi(model, t, self.p, g, "liquid") - ln_phi(model, t, self.p, w, "vapor")
            else:
                w = g/k/np.sum(g/k)
                new = ln_phi(model, t, self.p, w, "liquid") - ln_phi(model, t, self.p, g, "vapor")
            step = np.max(np.abs(new-ln_k))
            ln_k = new
            if step < 1e-13:
                break
        else:
            # Outside the bracket the incipient-phase substitution need not converge;
            # the scan skips such temperatures and brentq never accepts them.
            raise ValueError(f"{self.kind} inner substitution not converged at {t} K")
        if np.max(np.abs(ln_k)) < 1e-4:
            raise ValueError("trivial solution")
        self.ln_k = ln_k
        k = np.exp(ln_k)
        return math.log(np.sum(g*k)) if self.kind == "bubble" else math.log(np.sum(g/k))

    def solve(self):
        grid = np.arange(T_RANGE[0], T_RANGE[1]+1e-9, 2.5)
        prev = None
        for t in grid:
            try:
                f = self.inner(t)
            except ValueError:
                self.ln_k, prev = None, None
                continue
            if prev is not None and np.sign(f) != np.sign(prev[1]):
                lo = prev[0]
                break
            prev = (t, f, self.ln_k.copy())
        else:
            raise AssertionError(f"{self.kind} point not bracketed in {T_RANGE}")
        self.ln_k = prev[2]
        t = brentq(self.inner, lo, t, xtol=1e-12, rtol=1e-15, maxiter=200)
        self.inner(t)
        k = np.exp(self.ln_k)
        g = self.given
        if self.kind == "bubble":
            liquid, vapor = g, g*k/np.sum(g*k)
        else:
            liquid, vapor = g/k/np.sum(g/k), g
        # Final residual check of the converged equilibrium.
        model, _ = equilibrium.model_at(self.ids, t, self.variant, self.kij)
        resid = np.max(np.abs(np.log(vapor)+ln_phi(model, t, self.p, vapor, "vapor")
                              - np.log(liquid)-ln_phi(model, t, self.p, liquid, "liquid")))
        assert resid < 1e-11, resid
        assert abs(np.sum(liquid)-1) < 1e-12 and abs(np.sum(vapor)-1) < 1e-12
        return {"T": t, "liquid": liquid, "vapor": vapor, "K": k, "fugacity_residual": float(resid)}


def sat_point(ids, x, p, variant, kij, kind):
    return Saturation(ids, x, p, variant, kind, kij).solve()


def phase_h(ids, x, t, p, variant, kij, phase):
    c = caloric.case(ids, list(map(float, x)), t, p, variant, kij)
    stable = [r for r in c["roots"] if r["mechanically_stable"]]
    return (stable[0] if phase == "liquid" else stable[-1])["h"]


def alphas_from(sat, hk):
    return sat["K"]/sat["K"][hk]


def fenske(alpha, z, flow, lk, hk, r_lk, r_hk):
    n_min = math.log(r_lk/(1-r_lk)*r_hk/(1-r_hk))/math.log(alpha[lk])
    ln_db = math.log((1-r_hk)/r_hk) + n_min*np.log(alpha)
    rec = 1/(1+np.exp(-ln_db))
    d = flow*z*rec
    b = flow*z-d
    return n_min, d, b


def shortcut(case):
    ids, z, flow, variant = case["ids"], np.array(case["z"]), case["F_mol_s"], case["variant"]
    kij = np.zeros((len(ids), len(ids))).tolist()
    lk, hk = ids.index(case["light_key"]), ids.index(case["heavy_key"])
    r_lk, r_hk, p = case["r_LK"], case["r_HK"], case["p"]
    assert all(v > 0 for v in z) and p <= case["feed_p"] and r_lk+r_hk > 1
    # Feed: its own equilibrium state (T_F, p_F) -> h_F (isenthalpic let-down if p_F > p).
    feed = ps.tp_state("feed", ids, list(z), case["feed_T"], case["feed_p"], variant)
    h_f = feed["target_h"]
    # Thermal condition at column pressure.
    bub_z, dew_z = sat_point(ids, z, p, variant, kij, "bubble"), sat_point(ids, z, p, variant, kij, "dew")
    h_bub = phase_h(ids, z, bub_z["T"], p, variant, kij, "liquid")
    h_dew = phase_h(ids, z, dew_z["T"], p, variant, kij, "vapor")
    q = (h_dew-h_f)/(h_dew-h_bub)
    # Successive substitution on alpha.
    alpha = alphas_from(bub_z, hk)
    alpha0 = alpha.copy()
    iterations = 0
    for _ in range(ALPHA_MAX_ITER):
        assert alpha[lk] > 1
        n_min, d, b = fenske(alpha, z, flow, lk, hk, r_lk, r_hk)
        x_d, x_b = d/d.sum(), b/b.sum()
        top = sat_point(ids, x_d, p, variant, kij, "dew")
        bot = sat_point(ids, x_b, p, variant, kij, "bubble")
        new = np.sqrt(alphas_from(top, hk)*alphas_from(bot, hk))
        iterations += 1
        change = np.max(np.abs(np.log(new)-np.log(alpha)))
        alpha = new
        if change < ALPHA_TOL:
            break
    else:
        raise AssertionError("alpha iteration not converged")
    n_min, d, b = fenske(alpha, z, flow, lk, hk, r_lk, r_hk)
    assert abs(d[lk]/(flow*z[lk])-r_lk) < 1e-12 and abs(b[hk]/(flow*z[hk])-r_hk) < 1e-12
    dd, bb = d.sum(), b.sum()
    x_d, x_b = d/dd, b/bb
    # Adjacent keys.
    assert not any(1 < a < alpha[lk] for i, a in enumerate(alpha) if i not in (lk, hk))
    # Underwood.
    def uw(theta):
        return np.sum(alpha*z/(alpha-theta))-(1-q)
    span = alpha[lk]-1
    theta = brentq(uw, 1+1e-14*span, alpha[lk]-1e-14*span, xtol=1e-15, rtol=1e-15, maxiter=500)
    r_min = float(np.sum(alpha*x_d/(alpha-theta))-1)
    assert r_min > 0
    spec = case["reflux"]
    reflux = spec["value"]*r_min if spec["kind"] == "multiple_of_minimum" else spec["value"]
    assert reflux > r_min
    x_g = (reflux-r_min)/(reflux+1)
    y_g = 1-math.exp((1+54.4*x_g)/(11+117.2*x_g)*(x_g-1)/math.sqrt(x_g))
    n = (n_min+y_g)/(1-y_g)
    ratio = (bb/dd*z[hk]/z[lk]*(x_b[lk]/x_d[hk])**2)**0.206
    n_r = n*ratio/(1+ratio)
    # Products and duties at the final splits.
    top = sat_point(ids, x_d, p, variant, kij, "dew")
    dist = sat_point(ids, x_d, p, variant, kij, "bubble")
    bot = sat_point(ids, x_b, p, variant, kij, "bubble")
    assert bot["T"] > top["T"]
    h_vtop = phase_h(ids, x_d, top["T"], p, variant, kij, "vapor")
    h_d = phase_h(ids, x_d, dist["T"], p, variant, kij, "liquid")
    h_b = phase_h(ids, x_b, bot["T"], p, variant, kij, "liquid")
    q_c = (reflux+1)*dd*(h_vtop-h_d)
    q_r = dd*h_d+bb*h_b+q_c-flow*h_f
    assert q_c > 0 and q_r > 0
    out = {
        **case, "kij": kij,
        "feed_state": {"phases": feed["phases"], "vapor_fraction": feed.get("vapor_fraction"),
                       "vapor": feed.get("vapor"), "liquid": feed.get("liquid"),
                       "feed_min_tm": feed["feed_min_tm"]},
        "feed_phase_count": feed["phases"],
        "feed_isenthalpic_letdown": case["feed_p"] > p,
        "feed_two_phase_at_column_p_by_q": bool(0 < q < 1),
        "feed_bubble_T": bub_z["T"], "feed_dew_T": dew_z["T"],
        "h_bub_z": h_bub, "h_dew_z": h_dew, "h_F": h_f, "q": q,
        "alpha_initial": alpha0.tolist(), "alpha": alpha.tolist(), "alpha_iterations": iterations,
        "alpha_top": alphas_from(top, hk).tolist(), "alpha_bottom": alphas_from(bot, hk).tolist(),
        "N_min": n_min, "theta": theta, "R_min": r_min, "R": reflux, "X": x_g, "Y": y_g, "N": n,
        "kirkbride_NR_over_NS": ratio, "N_R": n_r, "N_S": n-n_r, "feed_stage_from_top": n_r+1,
        "D": dd, "B": bb, "d": d.tolist(), "b": b.tolist(), "x_D": x_d.tolist(), "x_B": x_b.tolist(),
        "T_top": top["T"], "top_liquid": top["liquid"].tolist(), "T_D": dist["T"],
        "T_bot": bot["T"], "bottom_vapor": bot["vapor"].tolist(),
        "h_D": h_d, "h_B": h_b, "h_V_top": h_vtop, "Q_C": q_c, "Q_R": q_r, "net_duty_W": q_r-q_c,
        "max_fugacity_residual": max(s["fugacity_residual"] for s in (bub_z, dew_z, top, dist, bot)),
    }
    return out


F = 100.0
CASES = [
    # a. Binary ethane/propane, slightly subcooled feed at column pressure, R = 1.3 R_min.
    {"id": "c2c3-binary-PR-BM", "ids": ["ethane", "propane"], "z": [0.4, 0.6], "F_mol_s": F,
     "feed_T": 295.0, "feed_p": 2.0e6, "variant": "PR-BM", "light_key": "ethane", "heavy_key": "propane",
     "r_LK": 0.98, "r_HK": 0.97, "p": 2.0e6, "reflux": {"kind": "multiple_of_minimum", "value": 1.3}},
    {"id": "c2c3-binary-PR", "ids": ["ethane", "propane"], "z": [0.4, 0.6], "F_mol_s": F,
     "feed_T": 295.0, "feed_p": 2.0e6, "variant": "PR", "light_key": "ethane", "heavy_key": "propane",
     "r_LK": 0.98, "r_HK": 0.97, "p": 2.0e6, "reflux": {"kind": "multiple_of_minimum", "value": 1.3}},
    # b. Ternary methane/ethane/propane, two-phase feed, reflux ratio given
    #    (R = 1.1, about 1.54 R_min = 0.7133 from an exploratory run of this script).
    {"id": "c1c2c3-ternary-PR-BM", "ids": ["methane", "ethane", "propane"], "z": [0.2, 0.45, 0.35],
     "F_mol_s": F, "feed_T": 260.0, "feed_p": 2.5e6, "variant": "PR-BM", "light_key": "ethane",
     "heavy_key": "propane", "r_LK": 0.95, "r_HK": 0.95, "p": 2.5e6, "reflux": {"kind": "ratio", "value": 1.1}},
    # c. Liquid feed at 4 MPa let down isenthalpically into a 2.5 MPa column.
    {"id": "c1c2c3-letdown-PR-BM", "ids": ["methane", "ethane", "propane"], "z": [0.1, 0.5, 0.4],
     "F_mol_s": F, "feed_T": 270.0, "feed_p": 4.0e6, "variant": "PR-BM", "light_key": "ethane",
     "heavy_key": "propane", "r_LK": 0.97, "r_HK": 0.96, "p": 2.5e6,
     "reflux": {"kind": "multiple_of_minimum", "value": 1.4}},
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = [shortcut(c) for c in CASES]
    output = {"schema": "a17-shortcut-reference-v1", "generator": "validation/columns/generate-shortcut-references.py",
              "python": platform.python_version(),
              "packages": {n: importlib.metadata.version(n) for n in ["teqp", "numpy", "scipy", "cantera"]},
              "input_sha256": hashlib.sha256(caloric.SOURCE.read_bytes()).hexdigest(),
              "adapter_sha256": {name: hashlib.sha256((ROOT / "validation" / name).read_bytes()).hexdigest()
                                 for name in ["generate-references.py", "generate-flash-references.py",
                                              "rotating/generate-ps-references.py",
                                              "columns/generate-shortcut-references.py"]},
              "method": "Fenske-Underwood-Gilliland(Molokanov)-Kirkbride, A17 #185 contract: alpha by successive "
                        "substitution (geometric mean of top dew and bottom bubble), Hengstebeck-Geddes Fenske "
                        "distribution, total condenser, partial reboiler counted in N, Q_R by overall balance",
              "units": {"T": "K", "p": "Pa", "h": "J/mol", "flow": "mol/s", "duty": "W (Q_C leaves, Q_R enters)"},
              "alpha_iterations_counting": "number of top-dew/bottom-bubble alpha updates until max |d ln alpha| < 1e-10, starting from alpha at the feed bubble point",
              "independence": "teqp residual AD fugacities and Cantera NASA7 caloric adapter, SciPy brentq bubble/dew "
                              "temperatures over inner ln K substitution, multi-start TPD feed flash; shared component "
                              "parameters and explicit zero kij; no Rust calls",
              "scope": "Implementation agreement for the shortcut column contract on regular subcritical PR/PR-BM "
                       "states (200..350 K); not experimental or rigorous-column accuracy",
              "tolerances": TOLERANCES, "cases": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(f"Wrote {len(results)} shortcut column references to {args.output}")
    for r in results:
        print(r["id"], f"phases={r['feed_phase_count']} q={r['q']:.6f} it={r['alpha_iterations']} "
              f"Nmin={r['N_min']:.6f} Rmin={r['R_min']:.6f} R={r['R']:.6f} N={r['N']:.6f} NR={r['N_R']:.6f} "
              f"QC={r['Q_C']:.3f} QR={r['Q_R']:.3f} Ttop={r['T_top']:.4f} TD={r['T_D']:.4f} Tbot={r['T_bot']:.4f} "
              f"Tbub(z)={r['feed_bubble_T']:.3f} Tdew(z)={r['feed_dew_T']:.3f}")


if __name__ == "__main__":
    main()
