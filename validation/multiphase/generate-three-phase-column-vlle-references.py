"""Independent three-phase distillation references with VLLE stages (#188).

Python standard library only. No Rust code or Rust output is used.

Extends the decanter-only references (`generate-three-phase-column-references.py`,
fixture `three-phase-column-v1.json`, unchanged) to columns whose trays and
reboiler carry two liquids. Model (matched formulation, implementation
verification only):

- One column pressure p. Stage 1 is a total condenser with a decanter at the
  specified temperature T_D (liquid-liquid equilibrium of the condensate);
  fractions phi_L and phi_H of the light and heavy liquid return to stage 2.
  Stages 2..N-1 are ideal trays, stage N a partial reboiler. Each stage
  carries a vapor and one liquid, or a vapor and two liquids in
  vapor-liquid-liquid equilibrium (Seader, Henley and Roper, *Separation
  Process Principles*, 3rd ed., 2011, section 10.1; three-phase stages:
  Cairns and Furzer, Ind. Eng. Chem. Res. 29 (1990) 1349). The liquid leaving
  a stage downward is the total of its liquids.
- Gamma-phi NRTL for water/1-butanol, exactly as in the decanter-only
  generator: ChemSep/DECHEMA pair (Renon and Prausnitz, AIChE J. 14 (1968)
  135), water Wagner-Pruss and ChemSep 1-butanol vapor pressures, ideal-gas
  vapor, NASA7 ideal-gas enthalpies anchored at 298.15 K, Clausius-Clapeyron
  pure-liquid enthalpies and the NRTL excess enthalpy (complex-step
  temperature derivatives, Squire and Trapp, SIAM Rev. 40 (1998) 110).
- Peng-Robinson (Peng and Robinson, Ind. Eng. Chem. Fundam. 15 (1976) 59) for
  water/propane with the exact critical factors Omega_a, Omega_b, the 1976
  m(omega) for every omega, van der Waals one-fluid mixing with
  k_ij(water, propane) = 0.5 and no volume translation. The k_ij is the
  explicit water-hydrocarbon approximation of the #224 VLLE references
  (docs/vlle-status.md), not a regressed value. Liquid: smallest-volume real
  root, vapor: largest; enthalpy h = h_ig(T) - h_ig(298.15 K) + h_res with the
  closed-form PR departure function; NASA7 low-range ideal-gas polynomials.
- Specifications: decanter temperature, reflux fractions, bottoms flow
  B = bottoms_fraction * F. The reboiler energy balance yields its duty.

Data are read directly from the source transcriptions
(`data/source-extracts/*.provisional.json`). Shared parameters, no shared code.

Solution: Newton's method with a central-difference Jacobian, dense Gaussian
elimination with partial pivoting and a backtracking line search, in explicit
stage temperatures, ln mole fractions with summation equations and explicit
total phase flows (the Rust engine uses ln component flows without summation
equations). The phase amounts of three-phase stages are unconstrained in sign
(a negative-flash formulation): the phase set is found by an active-set loop
that starts with one liquid on every stage, adds a second liquid where a
1/4000 (NRTL) or log-dense (PR) composition-grid tangent-plane scan of the
stage liquid is negative, and removes it where a converged liquid amount is
not positive. Accepted solutions pass the grid scan on every stage liquid
(both liquids on three-phase stages), the decanter tie-line scan (for PR
with vapor-like trial roots), a decanter bubble-point check (NRTL-IG: no
vapor), and component and energy balance closure. When one liquid on every
stage does not converge, the loop restarts with two liquids where the initial
profile liquid is unstable, or on every stage at the binary three-phase point
(PR: initial bubble points from Wilson K-values, since below its bubble point
a PR mixture has only its liquid root).
"""

import argparse
import cmath
import hashlib
import json
import math
import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXTRACT_DIR = ROOT / "data/source-extracts"
EXTRACTS = [EXTRACT_DIR / name for name in (
    "liquid-nrtl.provisional.json",
    "liquid-chemsep-nrtl.provisional.json",
    "components.provisional.json",
    "components-chemsep.provisional.json",
)]
R = 8.31446261815324
T_REF = 298.15
OMEGA_A = 0.4572355289213822
OMEGA_B = 0.07779607390388846
SQRT2 = math.sqrt(2.0)
WATER_HYDROCARBON_KIJ = 0.5

# Declared before any comparison with the Rust engine ([absolute, relative]);
# the same budgets as the decanter-only references.
TOLERANCES = {
    "temperature_K": [1e-6, 0.0],
    "composition": [1e-8, 0.0],
    "molar_flow_mol_s": [1e-9, 1e-7],
    "duty_W": [1e-3, 1e-7],
}


def nrtl_case(**kw):
    base = dict(package="NRTL-IG", components=["water", "1-butanol"],
                pressure_Pa=101325.0, decanter_temperature_K=320.0,
                feed=dict(temperature_K=330.0, pressure_Pa=2e5, molar_flow_mol_s=1.0,
                          z=[0.3, 0.7]),
                internal_stages=6, feed_stage=2, bottoms_fraction=0.71,
                initial_condensate=[0.72, 0.28])
    base.update(kw)
    return base


def pr_case(**kw):
    base = dict(package="PR", components=["water", "propane"],
                kij=[[0.0, WATER_HYDROCARBON_KIJ], [WATER_HYDROCARBON_KIJ, 0.0]],
                pressure_Pa=20e5, decanter_temperature_K=300.0,
                feed=dict(temperature_K=400.0, pressure_Pa=25e5, molar_flow_mol_s=1.0,
                          z=[0.02, 0.98]),
                internal_stages=6, feed_stage=3, reflux=dict(kind="light_total"),
                bottoms_fraction=0.9802, initial_condensate=[0.007, 0.993])
    base.update(kw)
    return base


CASES = [
    # Most of both liquids refluxed: the top section carries two liquids.
    nrtl_case(id="nrtl-butanol-fractions-0.95-0.9-split-trays",
              reflux=dict(kind="fractions", light=0.95, heavy=0.9)),
    # All light, 80 % heavy liquid refluxed: two liquids on the upper trays only.
    nrtl_case(id="nrtl-butanol-fractions-1-0.8-split-upper-trays",
              reflux=dict(kind="fractions", light=1.0, heavy=0.8)),
    # Wet propane drying with free-water decant: free water on the upper
    # trays, a dry stripping section and reboiler.
    pr_case(id="pr-propane-drying-free-water-upper-trays"),
    pr_case(id="pr-propane-drying-dry-reboiler", feed=dict(
        temperature_K=400.0, pressure_Pa=25e5, molar_flow_mol_s=1.0, z=[0.05, 0.95]),
        bottoms_fraction=0.951, initial_condensate=[0.007, 0.993]),
]


def reflux_fractions(r):
    kind = r["kind"]
    if kind == "light_total":
        return 1.0, 0.0
    if kind == "heavy_total":
        return 0.0, 1.0
    if kind == "light_ratio":
        return r["ratio"] / (1 + r["ratio"]), 0.0
    if kind == "heavy_ratio":
        return 0.0, r["ratio"] / (1 + r["ratio"])
    return r["light"], r["heavy"]


# ---------------------------------------------------------------- data


DOCS = [json.loads(p.read_text(encoding="utf-8")) for p in EXTRACTS]


def component(cid):
    for doc in DOCS:
        for c in doc.get("components", []):
            if c["id"] == cid:
                return c
    raise KeyError(cid)


def pure_vapor_pressure(cid):
    for doc in DOCS:
        for p in doc.get("pure_liquid", []):
            if p["id"] == cid:
                return p["vapor_pressure"]
    raise KeyError(cid)


def nrtl_pair(i, j):
    for doc in DOCS:
        for b in doc.get("nrtl_binary", []):
            if (b["i"], b["j"]) == (i, j):
                return b
    raise KeyError((i, j))


def nasa_h(coefficients, T):
    a1, a2, a3, a4, a5, a6, _ = coefficients
    return R * (a1 * T + a2 * T ** 2 / 2 + a3 * T ** 3 / 3 + a4 * T ** 4 / 4
                + a5 * T ** 5 / 5 + a6)


def ideal_gas_h(cid):
    ig = component(cid)["ideal_gas"]

    def h(T):
        real = T.real if isinstance(T, complex) else T
        c = ig["coefficients"][0 if real < ig["temperature_breaks_K"][1] else 1]
        return nasa_h(c, T) - nasa_h(ig["coefficients"][0], T_REF)
    return h


H_STEP = 1e-30


def d_dt(f, T):
    return f(complex(T, H_STEP)).imag / H_STEP


# ---------------------------------------------------------------- NRTL-IG


class Nrtl:
    """Binary gamma-phi model of the decanter-only generator."""

    def __init__(self, ids):
        assert ids == ["water", "1-butanol"]
        self.ids = ids
        self.vp = [pure_vapor_pressure(i) for i in ids]
        self.hig = [ideal_gas_h(i) for i in ids]
        self.pair = nrtl_pair("water", "1-butanol")

    def ln_psat(self, i, T):
        vp = self.vp[i]
        if vp["form"] == "pV":
            assert vp["using_tau_r"]
            tc = vp["reducing_temperature_K"]
            theta = 1 - T / tc
            s = sum(n * cmath.exp(t * cmath.log(theta)) for n, t in zip(vp["n"], vp["t"]))
            return math.log(vp["reducing_value"]) + tc / T * s
        assert vp["form"] == "chemsep-101"
        a, b, c, d, e = vp["coefficients"]
        return a + b / T + c * cmath.log(T) + d * T ** e

    def ln_gamma_c(self, x, T):
        def tau(t):
            return t["t0"] + t["t1_K"] / T + t["t2"] * cmath.log(T) + t["t3_per_K"] * T
        t12, t21 = tau(self.pair["tau_ij"]), tau(self.pair["tau_ji"])
        a = self.pair["alpha_ij"]["a0"] + self.pair["alpha_ij"]["a1_per_K"] * T
        g12, g21 = cmath.exp(-a * t12), cmath.exp(-a * t21)
        x1, x2 = x
        d1, d2 = x1 + x2 * g21, x2 + x1 * g12
        l1 = x2 * x2 * (t21 * (g21 / d1) ** 2 + t12 * g12 / d2 ** 2)
        l2 = x1 * x1 * (t12 * (g12 / d2) ** 2 + t21 * g21 / d1 ** 2)
        return [l1, l2]

    def ln_phi_liquid(self, x, T, p):
        g = self.ln_gamma_c(x, T)
        return [g[i].real + self.ln_psat(i, T).real - math.log(p) for i in range(2)]

    def ln_phi_vapor(self, y, T, p):
        return [0.0, 0.0]

    def h_liquid(self, x, T, p):
        def pure(i):
            return self.hig[i](T) - R * T * T * d_dt(lambda t: self.ln_psat(i, t), T)
        hE = -R * T * T * sum(x[i] * d_dt(lambda t, i=i: self.ln_gamma_c(x, t)[i], T)
                              for i in range(2))
        return sum(x[i] * pure(i) for i in range(2)) + hE

    def h_vapor(self, y, T, p):
        return sum(y[i] * self.hig[i](T) for i in range(2))

    def heavier(self, xa, xb):
        """Labels: the water-richer liquid is the heavy one."""
        return xa[0] > xb[0]

    def tpd_minimum(self, x, T, p):
        d = [math.log(a) + g for a, g in zip(x, self.ln_phi_liquid(x, T, p))]
        best = math.inf
        step = 4000
        for a in range(1, step):
            w = [a / step, 1 - a / step]
            g = self.ln_phi_liquid(w, T, p)
            best = min(best, sum(w[i] * (math.log(w[i]) + g[i] - d[i]) for i in range(2)))
        return best

    bubble_interval = (280.0, 440.0)


# ---------------------------------------------------------------- PR


def grid_values():
    v = [10.0 ** -k for k in range(14, 3, -1)]
    v += [k / 500 for k in range(1, 500)]
    v += [1 - 10.0 ** -k for k in range(3, 15)]
    return sorted(set(v))


class PengRobinson:
    def __init__(self, ids, kij):
        self.ids = ids
        self.c = [component(i)["constants"] for i in ids]
        self.hig = [ideal_gas_h(i) for i in ids]
        self.kij = kij
        self.tc = [c["critical_temperature_K"] for c in self.c]
        self.pc = [c["critical_pressure_Pa"] for c in self.c]
        self.m = [0.37464 + 1.54226 * c["acentric_factor"]
                  - 0.26992 * c["acentric_factor"] ** 2 for c in self.c]
        self.mw = [c["molar_mass_kg_per_mol"] for c in self.c]
        self.b = [OMEGA_B * R * tc / pc for tc, pc in zip(self.tc, self.pc)]

    def a_pure(self, i, T):
        """a_i and da_i/dT."""
        ac = OMEGA_A * (R * self.tc[i]) ** 2 / self.pc[i]
        sq = math.sqrt(T / self.tc[i])
        f = 1 + self.m[i] * (1 - sq)
        return ac * f * f, -ac * f * self.m[i] * sq / T

    def mixture(self, x, T):
        n = len(x)
        ap = [self.a_pure(i, T) for i in range(n)]
        aij = [[0.0] * n for _ in range(n)]
        daij = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                ai, dai = ap[i]
                aj, daj = ap[j]
                s = math.sqrt(ai * aj)
                aij[i][j] = (1 - self.kij[i][j]) * s
                daij[i][j] = (1 - self.kij[i][j]) * (dai * aj + ai * daj) / (2 * s)
        a = sum(x[i] * x[j] * aij[i][j] for i in range(n) for j in range(n))
        da = sum(x[i] * x[j] * daij[i][j] for i in range(n) for j in range(n))
        b = sum(x[i] * self.b[i] for i in range(n))
        return a, da, b, aij

    def roots(self, A, B):
        """Real roots Z > B of the PR cubic (Cardano, refined by Newton)."""
        c2, c1, c0 = -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)
        q = (3 * c1 - c2 * c2) / 9
        r = (9 * c2 * c1 - 27 * c0 - 2 * c2 ** 3) / 54
        disc = q ** 3 + r * r
        if disc > 0:
            s = math.copysign(abs(r + math.sqrt(disc)) ** (1 / 3), r + math.sqrt(disc))
            t = math.copysign(abs(r - math.sqrt(disc)) ** (1 / 3), r - math.sqrt(disc))
            zs = [s + t - c2 / 3]
        else:
            theta = math.acos(max(-1.0, min(1.0, r / math.sqrt(-q ** 3))))
            zs = [2 * math.sqrt(-q) * math.cos((theta + 2 * math.pi * k) / 3) - c2 / 3
                  for k in range(3)]
        out = []
        for z in zs:
            for _ in range(5):
                f = ((z + c2) * z + c1) * z + c0
                df = (3 * z + 2 * c2) * z + c1
                if df == 0:
                    break
                z -= f / df
            if z > B:
                out.append(z)
        return sorted(out)

    def state(self, x, T, p, root):
        a, da, b, aij = self.mixture(x, T)
        A, B = a * p / (R * T) ** 2, b * p / (R * T)
        zs = self.roots(A, B)
        if root == "liquid":
            candidates = [zs[0]]
        elif root == "vapor":
            candidates = [zs[-1]]
        else:
            candidates = zs
        best = None
        for Z in candidates:
            log_term = math.log((Z + (1 + SQRT2) * B) / (Z + (1 - SQRT2) * B))
            lnphi = []
            for i in range(len(x)):
                s = sum(x[j] * aij[i][j] for j in range(len(x)))
                lnphi.append(self.b[i] / b * (Z - 1) - math.log(Z - B)
                             - A / (2 * SQRT2 * B) * (2 * s / a - self.b[i] / b) * log_term)
            g = sum(xi * (math.log(xi) + f) for xi, f in zip(x, lnphi) if xi > 0)
            if best is None or g < best[0]:
                hres = R * T * (Z - 1) + (T * da - a) / (2 * SQRT2 * b) * log_term
                best = (g, Z, lnphi, hres)
        _, Z, lnphi, hres = best
        return dict(Z=Z, ln_phi=lnphi, h_res=hres,
                    density=sum(xi * m for xi, m in zip(x, self.mw)) / (Z * R * T / p))

    def ln_phi_liquid(self, x, T, p):
        return self.state(x, T, p, "liquid")["ln_phi"]

    def ln_phi_vapor(self, y, T, p):
        return self.state(y, T, p, "vapor")["ln_phi"]

    def h_ideal(self, x, T):
        return sum(x[i] * self.hig[i](T) for i in range(len(x)))

    def h_liquid(self, x, T, p):
        return self.h_ideal(x, T) + self.state(x, T, p, "liquid")["h_res"]

    def h_vapor(self, y, T, p):
        return self.h_ideal(y, T) + self.state(y, T, p, "vapor")["h_res"]

    def heavier(self, xa, xb, T, p):
        return (self.state(xa, T, p, "liquid")["density"]
                > self.state(xb, T, p, "liquid")["density"])

    def tpd_minimum(self, x, T, p):
        """Grid tangent-plane minimum from the liquid root of x; trial roots
        of minimum Gibbs energy (vapor-like and liquid-like trials)."""
        d = [math.log(a) + g for a, g in zip(x, self.ln_phi_liquid(x, T, p))]
        best = math.inf
        for a in grid_values():
            w = [a, 1 - a]
            g = self.state(w, T, p, "min")["ln_phi"]
            best = min(best, sum(w[i] * (math.log(w[i]) + g[i] - d[i]) for i in range(2)))
        return best

    bubble_interval = (300.0, 350.0)


def model_of(case):
    if case["package"] == "PR":
        return PengRobinson(case["components"], case["kij"])
    return Nrtl(case["components"])


def heavier(model, xa, xb, T, p):
    if isinstance(model, PengRobinson):
        return model.heavier(xa, xb, T, p)
    return model.heavier(xa, xb)


# ---------------------------------------------------------------- helpers


def solve_linear(a, b):
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if m[p][c] == 0.0:
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


def newton(residuals, v, kinds, tolerance, what, maxit=200):
    """Damped Newton; `kinds` marks each variable: 'ln' (ln mole fraction),
    'T', 'pos' (positive flow), 'free' (duty or a three-phase phase amount)."""
    r = residuals(v)
    for it in range(maxit):
        if max(map(abs, r)) < tolerance:
            return v, it
        jac = [[0.0] * len(v) for _ in r]
        for k in range(len(v)):
            h = 1e-7 if kinds[k] == "ln" else 1e-7 * max(abs(v[k]), 1e-2)
            vp, vm = v[:], v[:]
            vp[k] += h
            vm[k] -= h
            rp, rm = residuals(vp), residuals(vm)
            for row in range(len(r)):
                jac[row][k] = (rp[row] - rm[row]) / (2 * h)
        step = solve_linear(jac, [-x for x in r])
        limit = 1.0
        for k, d in enumerate(step):
            cap = {"ln": 2.0, "T": 10.0}.get(kinds[k])
            if cap and abs(d) > cap:
                limit = min(limit, cap / abs(d))
        base = math.sqrt(sum(x * x for x in r))
        lam = limit
        for _ in range(60):
            trial = [max(a + lam * b, 0.1 * a) if kinds[k] == "pos" else a + lam * b
                     for k, (a, b) in enumerate(zip(v, step))]
            try:
                rt = residuals(trial)
                if math.sqrt(sum(x * x for x in rt)) < base:
                    break
            except (ValueError, OverflowError, ZeroDivisionError):
                pass
            lam *= 0.5
        else:
            raise ArithmeticError(f"{what}: line search failed at {base:.3e}")
        v, r = trial, rt
    raise ArithmeticError(f"{what}: Newton did not converge")


def bubble_temperature(model, x, p):
    if isinstance(model, PengRobinson):
        return wilson_bubble(model, x, p)

    def f(T):
        y = x[:]
        for _ in range(100):
            lk = [a - b for a, b in zip(model.ln_phi_liquid(x, T, p), model.ln_phi_vapor(y, T, p))]
            ky = [math.exp(a) * b for a, b in zip(lk, x)]
            s = sum(ky)
            new = [a / s for a in ky]
            if max(abs(a - b) for a, b in zip(new, y)) < 1e-13:
                break
            y = new
        return s - 1, new

    lo, hi = model.bubble_interval
    assert f(lo)[0] < 0 < f(hi)[0], (x, f(lo)[0], f(hi)[0])
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid)[0] < 0:
            lo = mid
        else:
            hi = mid
    T = 0.5 * (lo + hi)
    return T, f(T)[1]


def wilson_bubble(model, x, p):
    """Initial-guess bubble point with Wilson K-values (Wilson, MIT 1968):
    below the bubble point a PR mixture has only its liquid root, so the
    model's own bubble function has the trivial solution there."""
    def k(T):
        return [model.pc[i] / p * math.exp(5.373 * (1 + model.c[i]["acentric_factor"])
                                            * (1 - model.tc[i] / T)) for i in range(len(x))]

    lo, hi = 150.0, 700.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if sum(a * b for a, b in zip(k(mid), x)) < 1:
            lo = mid
        else:
            hi = mid
    T = 0.5 * (lo + hi)
    y = [a * b for a, b in zip(k(T), x)]
    return T, [a / sum(y) for a in y]


def lle(model, T, p, guess):
    """Binary LLE by Newton on ln x (iso-fugacity plus summations)."""
    xa, xb = guess

    def res(v):
        a, b = [math.exp(v[0]), math.exp(v[1])], [math.exp(v[2]), math.exp(v[3])]
        fa, fb = model.ln_phi_liquid(a, T, p), model.ln_phi_liquid(b, T, p)
        return [v[i] + fa[i] - v[2 + i] - fb[i] for i in range(2)] + [sum(a) - 1, sum(b) - 1]

    v = [math.log(c) for c in xa + xb]
    v, _ = newton(res, v, ["ln"] * 4, 1e-14, "LLE")
    a, b = [math.exp(v[0]), math.exp(v[1])], [math.exp(v[2]), math.exp(v[3])]
    assert max(abs(i - j) for i, j in zip(a, b)) > 0.05
    return a, b


def three_phase_point(model, p, guess):
    """Binary vapor-liquid-liquid point at p: T, both liquids, vapor (Newton
    on T and ln mole fractions; iso-fugacity of all three phases)."""
    xa, xb = guess
    T0, _ = bubble_temperature(model, xa, p)

    def res(v):
        T = v[0]
        a, b, y = ([math.exp(c) for c in v[1:3]], [math.exp(c) for c in v[3:5]],
                   [math.exp(c) for c in v[5:7]])
        fa, fb = model.ln_phi_liquid(a, T, p), model.ln_phi_liquid(b, T, p)
        fv = model.ln_phi_vapor(y, T, p)
        return ([v[1 + i] + fa[i] - v[3 + i] - fb[i] for i in range(2)]
                + [v[5 + i] + fv[i] - v[1 + i] - fa[i] for i in range(2)]
                + [sum(a) - 1, sum(b) - 1, sum(y) - 1])

    # Start: the LLE at T0 and a vapor in equilibrium with its first liquid.
    a, b = lle(model, T0, p, guess)
    fa = model.ln_phi_liquid(a, T0, p)
    y = [a[i] * math.exp(fa[i]) for i in range(2)]
    y = [c / sum(y) for c in y]
    v = [T0] + [math.log(c) for c in a + b + y]
    v, _ = newton(res, v, ["T"] + ["ln"] * 6, 1e-13, "three-phase point")
    return v[0], [math.exp(c) for c in v[1:3]], [math.exp(c) for c in v[3:5]], \
        [math.exp(c) for c in v[5:7]]


# ---------------------------------------------------------------- column


class Column:
    def __init__(self, case, model):
        self.case, self.model = case, model
        self.C = len(case["components"])
        self.n = case["internal_stages"] + 1  # stages 2..N
        self.p = case["pressure_Pa"]
        self.td = case["decanter_temperature_K"]
        f = case["feed"]
        self.F = f["molar_flow_mol_s"]
        self.z = f["z"]
        self.feed_index = case["feed_stage"] - 2
        self.B = case["bottoms_fraction"] * self.F
        self.phi = reflux_fractions(case["reflux"])
        self.escale = self.F * 1e4
        self.hF = self.F * feed_enthalpy(model, f)

    # decanter: [ln xl (C), ln xh (C), A, Bh, Qc/escale]
    # two-phase stage: [T, ln x (C), ln y (C), L, V]
    # three-phase stage: [T, ln x1 (C), ln x2 (C), ln y (C), L1, L2, V]
    def layout(self, phases):
        C = self.C
        kinds = ["ln"] * (2 * C) + ["pos", "pos", "free"]
        for ph in phases:
            if ph == 2:
                kinds += ["T"] + ["ln"] * (2 * C) + ["pos", "pos"]
            else:
                kinds += ["T"] + ["ln"] * (3 * C) + ["free", "free", "pos"]
        return kinds

    def unpack(self, v, phases):
        C = self.C
        e = [math.exp(a) for a in v]
        d = dict(xl=e[0:C], xh=e[C:2 * C], A=v[2 * C], Bh=v[2 * C + 1],
                 qc=v[2 * C + 2] * self.escale)
        i = 2 * C + 3
        stages = []
        for ph in phases:
            if ph == 2:
                s = dict(T=v[i], x1=e[i + 1:i + 1 + C], y=e[i + 1 + C:i + 1 + 2 * C],
                         L1=v[i + 1 + 2 * C], V=v[i + 2 + 2 * C], x2=None, L2=0.0)
                i += 2 * C + 3
            else:
                s = dict(T=v[i], x1=e[i + 1:i + 1 + C], x2=e[i + 1 + C:i + 1 + 2 * C],
                         y=e[i + 1 + 2 * C:i + 1 + 3 * C], L1=v[i + 1 + 3 * C],
                         L2=v[i + 2 + 3 * C], V=v[i + 3 + 3 * C])
                i += 3 * C + 4
            stages.append(s)
        assert i == len(v)
        return d, stages

    def liquid_out(self, s):
        m, p = self.model, self.p
        flows = [s["L1"] * s["x1"][i] + (s["L2"] * s["x2"][i] if s["x2"] else 0.0)
                 for i in range(self.C)]
        h = s["L1"] * m.h_liquid(s["x1"], s["T"], p)
        if s["x2"]:
            h += s["L2"] * m.h_liquid(s["x2"], s["T"], p)
        return flows, h

    def residuals(self, v, phases):
        m, p, C = self.model, self.p, self.C
        d, st = self.unpack(v, phases)
        xl, xh, A, Bh, qc = d["xl"], d["xh"], d["A"], d["Bh"], d["qc"]
        out = []
        top = st[0]
        for i in range(C):
            out.append((top["V"] * top["y"][i] - A * xl[i] - Bh * xh[i]) / self.F)
        fl, fh = m.ln_phi_liquid(xl, self.td, p), m.ln_phi_liquid(xh, self.td, p)
        for i in range(C):
            out.append(math.log(xl[i]) + fl[i] - math.log(xh[i]) - fh[i])
        out.append(sum(xl) - 1)
        out.append(sum(xh) - 1)
        hl, hh = m.h_liquid(xl, self.td, p), m.h_liquid(xh, self.td, p)
        out.append((top["V"] * m.h_vapor(top["y"], top["T"], p) - qc - A * hl - Bh * hh)
                   / self.escale)
        a, b = self.phi
        reflux = [a * A * xl[i] + b * Bh * xh[i] for i in range(C)]
        reflux_h = a * A * hl + b * Bh * hh
        for k, s in enumerate(st):
            T = s["T"]
            if k == 0:
                lin, hlin = reflux, reflux_h
            else:
                lin, hlin = self.liquid_out(st[k - 1])
            if k + 1 < self.n:
                nxt = st[k + 1]
                vin = [nxt["V"] * c for c in nxt["y"]]
                hvin = nxt["V"] * m.h_vapor(nxt["y"], nxt["T"], p)
            else:
                vin, hvin = [0.0] * C, 0.0
            feed = [self.F * c for c in self.z] if k == self.feed_index else [0.0] * C
            hfeed = self.hF if k == self.feed_index else 0.0
            lout, hlout = self.liquid_out(s)
            for i in range(C):
                out.append((lin[i] + vin[i] + feed[i] - lout[i] - s["V"] * s["y"][i]) / self.F)
            f1 = m.ln_phi_liquid(s["x1"], T, p)
            fv = m.ln_phi_vapor(s["y"], T, p)
            for i in range(C):
                out.append(math.log(s["y"][i]) + fv[i] - math.log(s["x1"][i]) - f1[i])
            out.append(sum(s["x1"]) - 1)
            out.append(sum(s["y"]) - 1)
            if s["x2"]:
                f2 = m.ln_phi_liquid(s["x2"], T, p)
                for i in range(C):
                    out.append(math.log(s["x1"][i]) + f1[i] - math.log(s["x2"][i]) - f2[i])
                out.append(sum(s["x2"]) - 1)
            if k + 1 < self.n:
                out.append((hlin + hvin + hfeed - hlout - s["V"] * m.h_vapor(s["y"], T, p))
                           / self.escale)
            else:
                out.append((s["L1"] + s["L2"] - self.B) / self.F)
        return out

    def reboiler_duty(self, v, phases):
        d, st = self.unpack(v, phases)
        last = st[-1]
        if len(st) > 1:
            _, h_in = self.liquid_out(st[-2])
        else:
            a, b = self.phi
            h_in = (a * d["A"] * self.model.h_liquid(d["xl"], self.td, self.p)
                    + b * d["Bh"] * self.model.h_liquid(d["xh"], self.td, self.p))
        hfeed = self.hF if self.feed_index == self.n - 1 else 0.0
        _, h_out = self.liquid_out(last)
        return h_out + last["V"] * self.model.h_vapor(last["y"], last["T"], self.p) - h_in - hfeed

    def initial(self):
        """Linear liquid profile between an estimated top liquid and the
        bottoms, stage bubble temperatures, constant molar overflow."""
        m, C = self.model, self.C
        xl, xh = lle(m, self.td, self.p, guesses(self.case))
        if heavier(m, xl, xh, self.td, self.p):
            xl, xh = xh, xl
        yc = self.case["initial_condensate"]
        beta_h = (yc[0] - xl[0]) / (xh[0] - xl[0])
        fl, fh = self.phi
        share = (1 - fl) * (1 - beta_h) + (1 - fh) * beta_h
        D = self.F - self.B
        V = D / share
        A, Bh = (1 - beta_h) * V, beta_h * V
        prod = [(1 - fl) * A * xl[i] + (1 - fh) * Bh * xh[i] for i in range(C)]
        xb = [max(self.F * self.z[i] - prod[i], 1e-6 * self.B) for i in range(C)]
        xb = [a / sum(xb) for a in xb]
        reflux = V - D
        top = [(fl * A * xl[i] + fh * Bh * xh[i]) / reflux for i in range(C)]
        v = [math.log(c) for c in xl + xh] + [A, Bh, 0.0]
        for k in range(self.n):
            w = k / max(self.n - 1, 1)
            x = [(1 - w) * top[i] + w * xb[i] for i in range(C)]
            T, y = bubble_temperature(m, x, self.p)
            L = self.B if k == self.n - 1 else reflux + (self.F if k >= self.feed_index else 0)
            v += [T] + [math.log(c) for c in x + y] + [L, V]
        phases = [2] * self.n
        r = self.residuals(v, phases)
        v[2 * C + 2] = r[2 * C + 2]
        return v, phases

    def repack(self, v, phases, new):
        """Moves a solution to another phase set: a new second liquid starts
        from the stage LLE at the stage temperature (lever rule of the
        single liquid, 5..95 %); a removed one joins the first liquid."""
        m, C = self.model, self.C
        d, st = self.unpack(v, phases)
        out = [math.log(c) for c in d["xl"] + d["xh"]] + [d["A"], d["Bh"], d["qc"] / self.escale]
        for s, old, ph in zip(st, phases, new):
            T = s["T"]
            if ph == 2:
                flows, _ = self.liquid_out(s)
                L = sum(flows)
                x = [f / L for f in flows]
                out += [T] + [math.log(c) for c in x + s["y"]] + [L, s["V"]]
                continue
            if old == 3:
                out += ([T] + [math.log(c) for c in s["x1"] + s["x2"] + s["y"]]
                        + [s["L1"], s["L2"], s["V"]])
                continue
            x = s["x1"]
            xa, xb = lle(m, T, self.p, guesses(self.case))
            beta = min(max((x[0] - xa[0]) / (xb[0] - xa[0]), 0.05), 0.95)
            L = s["L1"]
            # The first liquid stays the one closer to the single liquid.
            if beta > 0.5:
                xa, xb, beta = xb, xa, 1 - beta
            out += ([T] + [math.log(c) for c in xa + xb + s["y"]]
                    + [(1 - beta) * L, beta * L, s["V"]])
        return out

    def at_three_phase_point(self, v, phases):
        """Binary: puts every three-phase stage at the invariant three-phase
        temperature and phase compositions of the column pressure."""
        T3, x1, x2, y = three_phase_point(self.model, self.p, guesses(self.case))
        d, st = self.unpack(v, phases)
        out = v[:2 * self.C + 3]
        for s, ph in zip(st, phases):
            if ph == 2:
                out += [s["T"]] + [math.log(c) for c in s["x1"] + s["y"]] + [s["L1"], s["V"]]
                continue
            L = s["L1"] + s["L2"]
            out += ([T3] + [math.log(c) for c in x1 + x2 + y]
                    + [0.95 * L, 0.05 * L, s["V"]])
        return out

    def solve(self):
        v, phases = self.initial()
        v0 = v
        seen = [list(phases)]
        total = 0
        for round in range(12):
            try:
                v, it = newton(lambda u: self.residuals(u, phases), v, self.layout(phases),
                               1e-12, self.case["id"], maxit=60 if round == 0 else 200)
            except ArithmeticError:
                if round > 0:
                    raise
                # One liquid everywhere does not converge: start again with
                # two liquids where the initial profile liquid is unstable,
                # or on every stage when none is (signed phase amounts then
                # remove the second liquid where it does not belong).
                _, st = self.unpack(v0, phases)
                new = [3 if self.model.tpd_minimum(s["x1"], s["T"], self.p) < -1e-9 else 2
                       for s in st]
                if 3 not in new:
                    new = [3] * self.n
                seen.append(new)
                v, phases = self.repack(v0, phases, new), new
                if self.C == 2:
                    v = self.at_three_phase_point(v, phases)
                continue
            total += it
            _, st = self.unpack(v, phases)
            new = list(phases)
            for k, s in enumerate(st):
                if phases[k] == 3:
                    if min(s["L1"], s["L2"]) <= 0:
                        new[k] = 2
                elif self.model.tpd_minimum(s["x1"], s["T"], self.p) < -1e-9:
                    new[k] = 3
            if new == phases:
                return v, phases, total
            assert new not in seen, (self.case["id"], "phase-set oscillation", new)
            seen.append(new)
            v = self.repack(v, phases, new)
            phases = new
        raise ArithmeticError(f"{self.case['id']}: phase set did not settle")


def guesses(case):
    if case["package"] == "PR":
        return [1e-3, 1 - 1e-3], [1 - 1e-10, 1e-10]
    return [0.6, 0.4], [0.99, 0.01]


def feed_enthalpy(model, f):
    """Single-phase feed enthalpy; the feed must be stable as one phase."""
    T, p, z = f["temperature_K"], f["pressure_Pa"], f["z"]
    if isinstance(model, PengRobinson):
        s = model.state(z, T, p, "min")
        # The feed is a superheated vapor: one real root, or the vapor root of
        # minimum Gibbs energy, and no negative grid tangent-plane distance.
        d = [math.log(a) + g for a, g in zip(z, s["ln_phi"])]
        for a in grid_values():
            w = [a, 1 - a]
            g = model.state(w, T, p, "min")["ln_phi"]
            assert sum(w[i] * (math.log(w[i]) + g[i] - d[i]) for i in range(2)) > -1e-9
        return model.h_ideal(z, T) + s["h_res"]
    # NRTL: a subcooled single liquid.
    lk = model.ln_phi_liquid(z, T, p)
    assert sum(z[i] * math.exp(lk[i]) for i in range(2)) < 1
    assert model.tpd_minimum(z, T, p) > -1e-9
    return model.h_liquid(z, T, p)


def reference(case):
    model = model_of(case)
    col = Column(case, model)
    v, phases, iterations = col.solve()
    d, st = col.unpack(v, phases)
    xl, xh, A, Bh, qc = d["xl"], d["xh"], d["A"], d["Bh"], d["qc"]
    qr = col.reboiler_duty(v, phases)
    assert qc > 0 and qr > 0, case["id"]
    td, p = col.td, col.p
    assert not heavier(model, xl, xh, td, p), case["id"]
    assert max(abs(a - b) for a, b in zip(xl, xh)) > 0.05
    decanter_tpd = model.tpd_minimum(xl, td, p)
    assert decanter_tpd > -1e-7, (case["id"], decanter_tpd)
    # No vapor at the decanter: for NRTL-IG the bubble pressure is below the
    # column pressure; for PR the grid scan above already includes
    # vapor-like (minimum-Gibbs root) trial phases.
    if isinstance(model, Nrtl):
        lk = model.ln_phi_liquid(xl, td, p)
        assert sum(xl[i] * math.exp(lk[i]) for i in range(col.C)) < 1, case["id"]
    stages = []
    for k, s in enumerate(st):
        T = s["T"]
        flows, _ = col.liquid_out(s)
        L = sum(flows)
        assert L > 0 and s["V"] > 0 and min(s["x1"] + s["y"]) > 0, (case["id"], k)
        record = dict(stage=k + 2, phase_count=phases[k], temperature_K=T, liquid_mol_s=L,
                      vapor_mol_s=s["V"], x=[f / L for f in flows], y=s["y"])
        tpd = model.tpd_minimum(s["x1"], T, p)
        assert tpd > -1e-9, (case["id"], k + 2, tpd)
        if phases[k] == 3:
            assert min(s["L1"], s["L2"]) > 0, (case["id"], k + 2)
            tpd = min(tpd, model.tpd_minimum(s["x2"], T, p))
            assert tpd > -1e-9, (case["id"], k + 2, tpd)
            a, b = (s["x1"], s["L1"]), (s["x2"], s["L2"])
            if heavier(model, a[0], b[0], T, p):
                a, b = b, a
            record.update(light_mol_s=a[1], x_light=a[0], heavy_mol_s=b[1], x_heavy=b[0])
        record["grid_tpd_minimum"] = tpd
        stages.append(record)
    fl, fh = col.phi
    kind = case["reflux"]["kind"]
    if kind.startswith("heavy"):
        distillate, decant = ((1 - fh) * Bh, xh), (A, xl)
    elif kind.startswith("light"):
        distillate, decant = ((1 - fl) * A, xl), (Bh, xh)
    else:
        distillate, decant = ((1 - fl) * A, xl), ((1 - fh) * Bh, xh)
    bflows, hb = col.liquid_out(st[-1])
    bottoms = (sum(bflows), [f / sum(bflows) for f in bflows])
    for i in range(col.C):
        out = distillate[0] * distillate[1][i] + decant[0] * decant[1][i] + bflows[i]
        assert abs(out - col.F * col.z[i]) < 1e-11, (case["id"], i)
    h_out = (distillate[0] * model.h_liquid(distillate[1], td, p)
             + decant[0] * model.h_liquid(decant[1], td, p) + hb)
    assert abs(h_out - col.hF - (qr - qc)) < 1e-6 * col.escale, case["id"]
    return dict(case, newton_iterations=iterations,
                decanter=dict(light_mol_s=A, x_light=xl, heavy_mol_s=Bh, x_heavy=xh,
                              grid_tpd_minimum=decanter_tpd),
                condenser_duty_W=qc, reboiler_duty_W=qr,
                products=dict(distillate=dict(molar_flow_mol_s=distillate[0], x=distillate[1]),
                              decant=dict(molar_flow_mol_s=decant[0], x=decant[1]),
                              bottoms=dict(molar_flow_mol_s=bottoms[0], x=bottoms[1])),
                profile=stages)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = dict(
        schema="a17-three-phase-column-vlle-reference-v1",
        issue=188,
        method=("Newton on explicit stage temperatures, ln mole fractions with summation "
                "equations and signed phase amounts; active-set phase-set loop with "
                "composition-grid TPD; decanter LLE at the specified temperature. NRTL-IG "
                "water/1-butanol (gamma-phi, complex-step caloric derivatives) and PR "
                "water/propane (closed-form departure functions). Python standard library, "
                "no Rust calls or Rust outputs."),
        scope=("Implementation agreement for the NRTL-IG and PR packages. Not experimental "
               "accuracy: the water/1-butanol pair is a 1 atm VLE regression and the PR "
               "water-propane kij = 0.5 is an explicit approximation (docs/vlle-status.md)."),
        generator=dict(python=platform.python_version(),
                       sha256=sha256(Path(__file__).resolve())),
        inputs={p.name: sha256(p) for p in EXTRACTS},
        tolerances=TOLERANCES,
        cases=[reference(case) for case in CASES],
    )
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    for case in out["cases"]:
        print(case["id"], case["newton_iterations"], "iterations; phases",
              [s["phase_count"] for s in case["profile"]], "T",
              [round(s["temperature_K"], 4) for s in case["profile"]],
              "Qc", round(case["condenser_duty_W"], 3), "Qr", round(case["reboiler_duty_W"], 3))


if __name__ == "__main__":
    main()
