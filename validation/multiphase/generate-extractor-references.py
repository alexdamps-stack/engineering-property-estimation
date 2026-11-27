"""Independent countercurrent liquid-liquid extractor references (#189).

Python standard library only. No Rust code or Rust output is used.

Model (matched formulation, implementation verification only):

- N ideal equilibrium stages at one temperature T and pressure p, numbered
  from the feed end. The feed (raffinate phase R) enters stage 1 and the
  solvent (extract phase E) enters stage N; the extract leaves stage 1 and
  the raffinate stage N (Seader, Henley and Roper, *Separation Process
  Principles*, 3rd ed., 2011, ch. 8 and section 10.5; Treybal, *Mass-Transfer
  Operations*, 3rd ed., 1980, ch. 10).
- Liquid-liquid equilibrium x_i^E gamma_i^E = x_i^R gamma_i^R. In the
  gamma-phi package both liquids share the pure-liquid vapor-pressure
  standard state, which cancels, so only the NRTL activity coefficients
  enter (Renon and Prausnitz, AIChE J. 14 (1968) 135):

      ln gamma_i = sum_j tau_ji G_ji x_j / sum_k G_ki x_k
                   + sum_j x_j G_ij / sum_k G_kj x_k
                     * (tau_ij - sum_m x_m tau_mj G_mj / sum_k G_kj x_k),
      G_ij = exp(-alpha_ij tau_ij),
      tau_ij = t0 + t1/T + t2 ln T + t3 T,  alpha_ij = a0 + a1 T.

- Parameters are read directly from the source transcriptions
  `data/source-extracts/liquid-nrtl.provisional.json` (methanol/water) and
  `liquid-chemsep-nrtl.provisional.json` (water/1-butanol and
  methanol/1-butanol; ChemSep/DECHEMA 1 atm VLE regressions as distributed
  with DWSIM, Artistic License 2.0). Shared parameters, no shared code.

Solution: Newton's method on all stage unknowns (both phase compositions and
both phase flows per stage) with a central-difference Jacobian, dense Gaussian
elimination with partial pivoting and a backtracking line search, started from
a Rachford-Rice LLE flash of the combined feed and solvent. Every stage tie
line is then checked against a composition-grid tangent-plane scan (step
1/400 in mole fraction): no grid point may have a tangent-plane distance below
-1e-7, which is evidence (not a proof) that the two liquids are the stable
split and no third liquid forms.

These parameters were regressed from 1 atm VLE data and do not describe
ambient mutual solubilities accurately (see docs/vlle-status.md); the
references verify the implementation, not experimental accuracy.
"""

import argparse
import hashlib
import json
import math
import platform
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTRACTS = [ROOT / "data/source-extracts/liquid-nrtl.provisional.json",
            ROOT / "data/source-extracts/liquid-chemsep-nrtl.provisional.json"]

IDS = ["water", "methanol", "1-butanol"]

# Declared before any comparison with the Rust engine ([absolute, relative]).
# The reference converges to a scaled residual below 1e-13 and the engine to
# 1e-11; the budgets allow for that and for the conditioning of the stage
# equations near the solvent end.
TOLERANCES = {
    "composition": [1e-8, 0.0],
    "molar_flow_mol_s": [1e-10, 1e-8],
}


def stream(flow, z):
    return dict(molar_flow_mol_s=flow, z=z)


CASES = [
    dict(id="methanol-from-water-3-stage-298K", stages=3, temperature_K=298.15,
         pressure_Pa=101325.0, feed=stream(1.0, [0.92, 0.08, 0.0]),
         solvent=stream(0.15, [0.0, 0.0, 1.0])),
    dict(id="wet-solvent-8-stage-310K", stages=8, temperature_K=310.0,
         pressure_Pa=200000.0, feed=stream(1.0, [0.95, 0.05, 0.0]),
         solvent=stream(0.2, [0.3, 0.0, 0.7])),
    dict(id="single-stage-298K", stages=1, temperature_K=298.15,
         pressure_Pa=101325.0, feed=stream(1.0, [0.9, 0.1, 0.0]),
         solvent=stream(0.15, [0.0, 0.0, 1.0])),
]


def load_parameters():
    pairs = {}
    for path in EXTRACTS:
        doc = json.loads(path.read_text(encoding="utf-8"))
        for b in doc["nrtl_binary"]:
            pairs[(b["i"], b["j"])] = (b["tau_ij"], b["tau_ji"], b["alpha_ij"])
    return pairs


PAIRS = load_parameters()


def tau_value(t, T):
    return t["t0"] + t["t1_K"] / T + t["t2"] * math.log(T) + t["t3_per_K"] * T


@lru_cache(maxsize=None)
def matrices(T):
    n = len(IDS)
    tau = [[0.0] * n for _ in range(n)]
    alpha = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (IDS[i], IDS[j]) in PAIRS:
                tij, tji, a = PAIRS[(IDS[i], IDS[j])]
            else:
                tji, tij, a = PAIRS[(IDS[j], IDS[i])]
            tau[i][j] = tau_value(tij, T)
            alpha[i][j] = a["a0"] + a["a1_per_K"] * T
    G = [[math.exp(-alpha[i][j] * tau[i][j]) for j in range(n)] for i in range(n)]
    return tau, G


def ln_gamma(x, T):
    tau, G = matrices(T)
    n = len(x)
    den = [sum(G[k][j] * x[k] for k in range(n)) for j in range(n)]
    num = [sum(x[m] * tau[m][j] * G[m][j] for m in range(n)) for j in range(n)]
    return [num[i] / den[i]
            + sum(x[j] * G[i][j] / den[j] * (tau[i][j] - num[j] / den[j]) for j in range(n))
            for i in range(n)]


def solve_linear(a, b):
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if abs(m[p][c]) == 0.0:
            raise ArithmeticError("singular Jacobian")
        m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            if f:
                for k in range(c, n + 1):
                    m[r][k] -= f * m[c][k]
    out = [0.0] * n
    for r in range(n - 1, -1, -1):
        out[r] = (m[r][n] - sum(m[r][k] * out[k] for k in range(r + 1, n))) / m[r][r]
    return out


def rachford_rice_lle(z, T):
    """Two-liquid flash: successive substitution on K = gamma^R / gamma^E."""
    n = len(z)
    # Trial phases: butanol-rich and water-rich.
    xe = [0.3, 0.1, 0.6]
    xr = [0.95, 0.04, 0.01]
    beta = 0.5
    for _ in range(5000):
        ge, gr = ln_gamma(xe, T), ln_gamma(xr, T)
        K = [math.exp(gr[i] - ge[i]) for i in range(n)]
        lo, hi = max(-1 / (k - 1) for k in K if k > 1), min(-1 / (k - 1) for k in K if k < 1)
        f = lambda b: sum(z[i] * (K[i] - 1) / (1 + b * (K[i] - 1)) for i in range(n))
        lo, hi = lo + 1e-15, hi - 1e-15
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        beta = 0.5 * (lo + hi)
        nxr = [z[i] / (1 + beta * (K[i] - 1)) for i in range(n)]
        nxe = [K[i] * nxr[i] for i in range(n)]
        sr, se = sum(nxr), sum(nxe)
        nxr = [v / sr for v in nxr]
        nxe = [v / se for v in nxe]
        change = max(abs(a - b) for a, b in zip(nxr + nxe, xr + xe))
        xr, xe = nxr, nxe
        if change < 1e-12:
            assert 0 < beta < 1, ("single liquid: negative flash", beta)
            return xe, xr, beta
    raise ArithmeticError("LLE flash did not converge")


def residuals(v, case):
    n, c = case["stages"], len(IDS)
    T = case["temperature_K"]
    F = [case["feed"]["molar_flow_mol_s"] * zi for zi in case["feed"]["z"]]
    S = [case["solvent"]["molar_flow_mol_s"] * zi for zi in case["solvent"]["z"]]
    scale = [F[i] + S[i] for i in range(c)]
    w = 2 * c + 2
    xe = [v[j * w: j * w + c] for j in range(n)]
    xr = [v[j * w + c: j * w + 2 * c] for j in range(n)]
    E = [v[j * w + 2 * c] for j in range(n)]
    R = [v[j * w + 2 * c + 1] for j in range(n)]
    out = []
    for j in range(n):
        rin = F if j == 0 else [R[j - 1] * x for x in xr[j - 1]]
        ein = S if j == n - 1 else [E[j + 1] * x for x in xe[j + 1]]
        for i in range(c):
            out.append((rin[i] + ein[i] - R[j] * xr[j][i] - E[j] * xe[j][i]) / scale[i])
        ge, gr = ln_gamma(xe[j], T), ln_gamma(xr[j], T)
        for i in range(c):
            out.append(xe[j][i] * math.exp(ge[i]) - xr[j][i] * math.exp(gr[i]))
        out.append(sum(xe[j]) - 1)
        out.append(sum(xr[j]) - 1)
    return out


def norm2(r):
    return math.sqrt(sum(x * x for x in r))


def solve(case):
    n, c = case["stages"], len(IDS)
    T = case["temperature_K"]
    fF, fS = case["feed"]["molar_flow_mol_s"], case["solvent"]["molar_flow_mol_s"]
    z = [(fF * a + fS * b) / (fF + fS)
         for a, b in zip(case["feed"]["z"], case["solvent"]["z"])]
    xe, xr, beta = rachford_rice_lle(z, T)
    v = []
    for _ in range(n):
        v += xe + xr + [beta * (fF + fS), (1 - beta) * (fF + fS)]
    r = residuals(v, case)
    for it in range(100):
        if max(abs(x) for x in r) < 1e-14:
            return v, it
        jac = [[0.0] * len(v) for _ in r]
        for k in range(len(v)):
            h = 1e-7 * max(abs(v[k]), 1e-3)
            vp, vm = v[:], v[:]
            vp[k] += h
            vm[k] -= h
            rp, rm = residuals(vp, case), residuals(vm, case)
            for row in range(len(r)):
                jac[row][k] = (rp[row] - rm[row]) / (2 * h)
        step = solve_linear(jac, [-x for x in r])
        base = norm2(r)
        lam = 1.0
        for _ in range(60):
            trial = [a + lam * b for a, b in zip(v, step)]
            if all(t > 0 for t in trial):
                rt = residuals(trial, case)
                if norm2(rt) < base:
                    break
            lam *= 0.5
        else:
            raise ArithmeticError(f"{case['id']}: line search failed at {base:.3e}")
        v, r = trial, rt
    raise ArithmeticError(f"{case['id']}: cascade Newton did not converge")


def tpd_minimum(x, T, step=400):
    lx = [math.log(a) + g for a, g in zip(x, ln_gamma(x, T))]
    best = math.inf
    for a in range(1, step):
        for b in range(1, step - a):
            w = [a / step, b / step, (step - a - b) / step]
            g = ln_gamma(w, T)
            d = sum(w[i] * (math.log(w[i]) + g[i] - lx[i]) for i in range(3))
            best = min(best, d)
    return best


def reference(case):
    v, iterations = solve(case)
    n, c = case["stages"], len(IDS)
    w = 2 * c + 2
    stages = []
    for j in range(n):
        xe, xr = v[j * w: j * w + c], v[j * w + c: j * w + 2 * c]
        E, R = v[j * w + 2 * c], v[j * w + 2 * c + 1]
        assert E > 0 and R > 0, case["id"]
        separation = max(abs(a - b) for a, b in zip(xe, xr))
        assert separation > 0.05, (case["id"], j, separation)
        tpd = tpd_minimum(xr, case["temperature_K"])
        assert tpd > -1e-7, (case["id"], j, tpd)
        stages.append(dict(stage=j + 1, extract_mol_s=E, raffinate_mol_s=R,
                           x_extract=xe, x_raffinate=xr, grid_tpd_minimum=tpd))
    # Overall balance closure of the reference itself.
    F = [case["feed"]["molar_flow_mol_s"] * zi for zi in case["feed"]["z"]]
    S = [case["solvent"]["molar_flow_mol_s"] * zi for zi in case["solvent"]["z"]]
    first, last = stages[0], stages[-1]
    for i in range(c):
        out = first["extract_mol_s"] * first["x_extract"][i] + \
            last["raffinate_mol_s"] * last["x_raffinate"][i]
        assert abs(out - F[i] - S[i]) < 1e-12, (case["id"], i)
    return dict(case, newton_iterations=iterations, profile=stages)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = dict(
        schema="a17-liquid-extractor-reference-v1",
        issue=189,
        components=IDS,
        method=("Newton on stage compositions and phase flows with NRTL activity "
                "coefficients (pure-liquid standard states cancel); Rachford-Rice "
                "LLE start; composition-grid TPD acceptance per stage. Python "
                "standard library, no Rust calls or Rust outputs."),
        scope=("Implementation agreement for the gamma-phi NRTL-IG package on "
               "the ChemSep/DECHEMA water/methanol/1-butanol pairs. Not experimental "
               "accuracy: the pairs are 1 atm VLE regressions."),
        generator=dict(python=platform.python_version(),
                       sha256=sha256(Path(__file__).resolve())),
        inputs={p.name: sha256(p) for p in EXTRACTS},
        tolerances=TOLERANCES,
        cases=[reference(case) for case in CASES],
    )
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    for case in out["cases"]:
        last = case["profile"][-1]
        print(case["id"], case["newton_iterations"], "iterations; raffinate methanol",
              last["x_raffinate"][1])


if __name__ == "__main__":
    main()
