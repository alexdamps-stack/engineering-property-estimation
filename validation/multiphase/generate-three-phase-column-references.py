"""Independent three-phase distillation column references (#188).

Python standard library only. No Rust code or Rust output is used.

Model (matched formulation, implementation verification only):

- One column pressure p. Stage 1 is a total condenser with a decanter at a
  specified temperature T_D: the vapor from stage 2 condenses completely and
  splits into a light and a heavy liquid (liquid-liquid equilibrium
  x^L_i gamma^L_i = x^H_i gamma^H_i). Fractions phi_L and phi_H of the two
  liquids return to stage 2 as reflux. Stages 2..N-1 are ideal trays, stage N a
  partial reboiler; each carries one liquid and one vapor
  (Seader, Henley and Roper, *Separation Process Principles*, 3rd ed., 2011,
  section 10.1; heterogeneous azeotropic distillation: Doherty and Malone,
  *Conceptual Design of Distillation Systems*, 2001, ch. 8).
- Gamma-phi with an ideal-gas vapor: K_i = gamma_i p_i^sat(T) / p. NRTL activity
  coefficients (Renon and Prausnitz, AIChE J. 14 (1968) 135).
- Caloric model: ideal-gas enthalpy from the NASA7 polynomials anchored at
  298.15 K; pure-liquid enthalpy h_i^L = h_i^ig - R T^2 d ln p_i^sat/dT
  (Clausius-Clapeyron with an ideal-gas vapor and negligible liquid volume);
  liquid mixture h^L = sum x_i h_i^L + h^E with h^E = -R T^2 sum x_i
  d ln gamma_i/dT. Temperature derivatives use the complex-step method
  (Squire and Trapp, SIAM Rev. 40 (1998) 110), exact to rounding.
- Specifications: decanter temperature, reflux fractions, bottoms flow
  B = bottoms_fraction * F. The reboiler energy balance yields its duty.

Data are read directly from the source transcriptions
(`data/source-extracts/*.provisional.json`): ChemSep/DECHEMA water/1-butanol
NRTL pair (1 atm VLE regression, as distributed with DWSIM, Artistic License
2.0), the water vapor-pressure ancillary (Wagner and Pruss form), the ChemSep
1-butanol vapor-pressure correlation, NASA7 ideal-gas coefficients. Shared
parameters, no shared code.

Solution: Newton's method with a central-difference Jacobian, dense Gaussian
elimination with partial pivoting and a backtracking line search, in explicit
stage temperatures, mole fractions and total phase flows with explicit
summation equations (the Rust engine uses ln component flows without summation
equations). Start: linear composition profile between an estimated top liquid
and the bottoms from the overall balance, stage bubble temperatures, constant
molar overflow. Accepted solutions pass: a 1/4000 composition-grid
tangent-plane scan of every stage liquid (no second liquid), of the decanter
light liquid (stable tie line), a bubble-pressure check of the decanter (no
vapor), and component and energy balance closure.
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
IDS = ["water", "1-butanol"]
R = 8.31446261815324
T_REF = 298.15

# Declared before any comparison with the Rust engine ([absolute, relative]).
# The reference converges to a scaled residual below 1e-12 and the engine to
# 1e-11 in its own variables; the budgets allow for the conditioning of the
# stripping section near pure 1-butanol (or pure water).
TOLERANCES = {
    "temperature_K": [1e-6, 0.0],
    "composition": [1e-8, 0.0],
    "molar_flow_mol_s": [1e-9, 1e-7],
    "duty_W": [1e-3, 1e-7],
}


def spec(**kw):
    base = dict(pressure_Pa=101325.0, decanter_temperature_K=320.0,
                feed=dict(temperature_K=330.0, pressure_Pa=2e5, molar_flow_mol_s=1.0,
                          z=[0.3, 0.7]),
                internal_stages=6, feed_stage=2, reflux=dict(kind="light_total"),
                bottoms_fraction=0.71)
    base.update(kw)
    return base


CASES = [
    spec(id="butanol-side-light-total-6-trays"),
    spec(id="butanol-side-light-ratio-8-trays", internal_stages=8, feed_stage=4,
         reflux=dict(kind="light_ratio", ratio=4.0), bottoms_fraction=0.65),
    spec(id="butanol-side-fractions-6-trays",
         reflux=dict(kind="fractions", light=0.9, heavy=0.2), bottoms_fraction=0.68),
    spec(id="water-side-heavy-total-6-trays",
         feed=dict(temperature_K=330.0, pressure_Pa=2e5, molar_flow_mol_s=1.0,
                   z=[0.996, 0.004]),
         reflux=dict(kind="heavy_total"), bottoms_fraction=0.992),
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


def load():
    docs = [json.loads(p.read_text(encoding="utf-8")) for p in EXTRACTS]
    pure, nasa, pair = {}, {}, None
    for doc in docs:
        for p in doc.get("pure_liquid", []):
            if p["id"] in IDS:
                pure[p["id"]] = p["vapor_pressure"]
        for b in doc.get("nrtl_binary", []):
            if (b["i"], b["j"]) == ("water", "1-butanol"):
                pair = b
        for c in doc.get("components", []):
            if c["id"] in IDS:
                nasa[c["id"]] = c["ideal_gas"]
    assert set(pure) == set(IDS) and set(nasa) == set(IDS) and pair
    return pure, nasa, pair


PURE, NASA, PAIR = load()


def ln_psat(i, T):
    vp = PURE[IDS[i]]
    if vp["form"] == "pV":
        assert vp["using_tau_r"]
        tc = vp["reducing_temperature_K"]
        theta = 1 - T / tc
        s = sum(n * cmath.exp(t * cmath.log(theta)) for n, t in zip(vp["n"], vp["t"]))
        return math.log(vp["reducing_value"]) + tc / T * s
    assert vp["form"] == "chemsep-101"
    a, b, c, d, e = vp["coefficients"]
    return a + b / T + c * cmath.log(T) + d * T ** e


def nasa_h(i, T):
    ig = NASA[IDS[i]]
    coeffs = ig["coefficients"][0 if (T.real if isinstance(T, complex) else T)
                                 < ig["temperature_breaks_K"][1] else 1]
    a1, a2, a3, a4, a5, a6, _ = coeffs
    return R * (a1 * T + a2 * T ** 2 / 2 + a3 * T ** 3 / 3 + a4 * T ** 4 / 4
                + a5 * T ** 5 / 5 + a6)


def h_ig(i, T):
    return nasa_h(i, T) - nasa_h(i, T_REF)


def tau_value(t, T):
    return t["t0"] + t["t1_K"] / T + t["t2"] * cmath.log(T) + t["t3_per_K"] * T


def ln_gamma_c(x, T):
    """Binary NRTL ln gamma (works with a complex temperature)."""
    t12, t21 = tau_value(PAIR["tau_ij"], T), tau_value(PAIR["tau_ji"], T)
    a = PAIR["alpha_ij"]["a0"] + PAIR["alpha_ij"]["a1_per_K"] * T
    g12, g21 = cmath.exp(-a * t12), cmath.exp(-a * t21)
    x1, x2 = x
    d1, d2 = x1 + x2 * g21, x2 + x1 * g12
    l1 = x2 * x2 * (t21 * (g21 / d1) ** 2 + t12 * g12 / d2 ** 2)
    l2 = x1 * x1 * (t12 * (g12 / d2) ** 2 + t21 * g21 / d1 ** 2)
    return [l1, l2]


def ln_gamma(x, T):
    return [v.real for v in ln_gamma_c(x, T)]


H = 1e-30


def d_dt(f, T):
    return f(complex(T, H)).imag / H


def psat(i, T):
    return math.exp(ln_psat(i, T).real)


def h_pure_liquid(i, T):
    return h_ig(i, T).real - R * T * T * d_dt(lambda t: ln_psat(i, t), T)


def h_liquid(x, T):
    hE = -R * T * T * sum(x[i] * d_dt(lambda t, i=i: ln_gamma_c(x, t)[i], T)
                          for i in range(2))
    return sum(x[i] * h_pure_liquid(i, T) for i in range(2)) + hE


def h_vapor(y, T):
    return sum(y[i] * h_ig(i, T).real for i in range(2))


def ln_k(x, T, p):
    g = ln_gamma(x, T)
    return [g[i] + ln_psat(i, T).real - math.log(p) for i in range(2)]


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


def bubble_temperature(x, p, lo=280.0, hi=440.0):
    f = lambda T: sum(x[i] * math.exp(ln_k(x, T, p)[i]) for i in range(2)) - 1
    assert f(lo) < 0 < f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lle(T, guess=(0.6, 0.99)):
    """Binary LLE by Newton on the two iso-activity equations."""
    v = list(guess)

    def res(v):
        xl, xh = [v[0], 1 - v[0]], [v[1], 1 - v[1]]
        gl, gh = ln_gamma(xl, T), ln_gamma(xh, T)
        return [math.log(xl[i]) + gl[i] - math.log(xh[i]) - gh[i] for i in range(2)]

    for _ in range(100):
        r = res(v)
        if max(map(abs, r)) < 1e-14:
            assert abs(v[0] - v[1]) > 0.05
            return v
        j = [[0.0, 0.0], [0.0, 0.0]]
        for k in range(2):
            h = 1e-7
            vp, vm = v[:], v[:]
            vp[k] += h
            vm[k] -= h
            rp, rm = res(vp), res(vm)
            for row in range(2):
                j[row][k] = (rp[row] - rm[row]) / (2 * h)
        s = solve_linear(j, [-a for a in r])
        v = [min(max(a + b, 1e-12), 1 - 1e-12) for a, b in zip(v, s)]
    raise ArithmeticError("LLE did not converge")


def tpd_minimum(x, T, step=4000):
    lx = [math.log(a) + g for a, g in zip(x, ln_gamma(x, T))]
    best = math.inf
    for a in range(1, step):
        w = [a / step, 1 - a / step]
        g = ln_gamma(w, T)
        best = min(best, sum(w[i] * (math.log(w[i]) + g[i] - lx[i]) for i in range(2)))
    return best


# ---------------------------------------------------------------- column


class Column:
    def __init__(self, case):
        self.case = case
        self.n = case["internal_stages"] + 1  # stages 2..N
        self.p = case["pressure_Pa"]
        self.td = case["decanter_temperature_K"]
        f = case["feed"]
        self.F = f["molar_flow_mol_s"]
        self.z = f["z"]
        self.hF = self.F * h_liquid(self.z, f["temperature_K"])
        self.feed_index = case["feed_stage"] - 2
        self.B = case["bottoms_fraction"] * self.F
        self.phi = reflux_fractions(case["reflux"])
        self.escale = self.F * 1e4

    # layout: decanter [xl0, xl1, xh0, xh1, A, Bh, Qc/escale],
    # stage k: [T, x0, x1, y0, y1, L, V]
    def unpack(self, v):
        d = v[:7]
        stages = [v[7 + 7 * k: 14 + 7 * k] for k in range(self.n)]
        return d, stages

    def residuals(self, v):
        d, st = self.unpack(v)
        xl, xh, A, Bh, qc = d[0:2], d[2:4], d[4], d[5], d[6] * self.escale
        out = []
        T2, x2, y2, L2, V2 = st[0][0], st[0][1:3], st[0][3:5], st[0][5], st[0][6]
        for i in range(2):
            out.append((V2 * y2[i] - A * xl[i] - Bh * xh[i]) / self.F)
        gl, gh = ln_gamma(xl, self.td), ln_gamma(xh, self.td)
        for i in range(2):
            out.append(math.log(xl[i]) + gl[i] - math.log(xh[i]) - gh[i])
        out.append(sum(xl) - 1)
        out.append(sum(xh) - 1)
        hl, hh = h_liquid(xl, self.td), h_liquid(xh, self.td)
        out.append((V2 * h_vapor(y2, T2) - qc - A * hl - Bh * hh) / self.escale)
        fl, fh = self.phi
        reflux = [fl * A * xl[i] + fh * Bh * xh[i] for i in range(2)]
        reflux_h = fl * A * hl + fh * Bh * hh
        for k in range(self.n):
            T, x, y, L, V = st[k][0], st[k][1:3], st[k][3:5], st[k][5], st[k][6]
            if k == 0:
                lin, hlin = reflux, reflux_h
            else:
                Tp, xp, Lp = st[k - 1][0], st[k - 1][1:3], st[k - 1][5]
                lin, hlin = [Lp * a for a in xp], Lp * h_liquid(xp, Tp)
            if k + 1 < self.n:
                Tn, yn, Vn = st[k + 1][0], st[k + 1][3:5], st[k + 1][6]
                vin, hvin = [Vn * a for a in yn], Vn * h_vapor(yn, Tn)
            else:
                vin, hvin = [0.0, 0.0], 0.0
            feed = [self.F * a for a in self.z] if k == self.feed_index else [0.0, 0.0]
            hfeed = self.hF if k == self.feed_index else 0.0
            for i in range(2):
                out.append((lin[i] + vin[i] + feed[i] - L * x[i] - V * y[i]) / self.F)
            lk = ln_k(x, T, self.p)
            for i in range(2):
                out.append(y[i] - math.exp(lk[i]) * x[i])
            out.append(sum(x) - 1)
            out.append(sum(y) - 1)
            if k + 1 < self.n:
                out.append((hlin + hvin + hfeed - L * h_liquid(x, T) - V * h_vapor(y, T))
                           / self.escale)
            else:
                out.append((L - self.B) / self.F)
        return out

    def reboiler_duty(self, v):
        _, st = self.unpack(v)
        T, x, y, L, V = st[-1][0], st[-1][1:3], st[-1][3:5], st[-1][5], st[-1][6]
        Tp, xp, Lp = st[-2][0], st[-2][1:3], st[-2][5]
        hfeed = self.hF if self.feed_index == self.n - 1 else 0.0
        return L * h_liquid(x, T) + V * h_vapor(y, T) - Lp * h_liquid(xp, Tp) - hfeed

    def initial(self):
        # Decanter phases at T_D and a heteroazeotrope-like condensate.
        xlw, xhw = lle(self.td)
        xl, xh = [xlw, 1 - xlw], [xhw, 1 - xhw]
        yc = 0.72
        beta_h = (yc - xlw) / (xhw - xlw)
        fl, fh = self.phi
        share = (1 - fl) * (1 - beta_h) + (1 - fh) * beta_h
        D = self.F - self.B
        V = D / share
        A, Bh = (1 - beta_h) * V, beta_h * V
        prod = [(1 - fl) * A * xl[i] + (1 - fh) * Bh * xh[i] for i in range(2)]
        xb = [max(self.F * self.z[i] - prod[i], 1e-6 * self.B) for i in range(2)]
        xb = [a / sum(xb) for a in xb]
        reflux = V - D
        top = [(fl * A * xl[i] + fh * Bh * xh[i]) / reflux for i in range(2)]
        v = xl + xh + [A, Bh, 0.0]
        for k in range(self.n):
            w = k / max(self.n - 1, 1)
            x = [(1 - w) * top[i] + w * xb[i] for i in range(2)]
            T = bubble_temperature(x, self.p)
            y = [math.exp(a) * b for a, b in zip(ln_k(x, T, self.p), x)]
            y = [a / sum(y) for a in y]
            L = self.B if k == self.n - 1 else reflux + (self.F if k >= self.feed_index else 0)
            v += [T] + x + y + [L, V]
        r = self.residuals(v)
        v[6] = r[6]  # decanter energy balance gives the starting duty
        return v

    def solve(self):
        v = self.initial()
        r = self.residuals(v)
        for it in range(200):
            if max(map(abs, r)) < 1e-12:
                return v, it
            jac = [[0.0] * len(v) for _ in r]
            for k in range(len(v)):
                h = 1e-7 * max(abs(v[k]), 1e-2)
                vp, vm = v[:], v[:]
                vp[k] += h
                vm[k] -= h
                rp, rm = self.residuals(vp), self.residuals(vm)
                for row in range(len(r)):
                    jac[row][k] = (rp[row] - rm[row]) / (2 * h)
            step = solve_linear(jac, [-x for x in r])
            base = math.sqrt(sum(x * x for x in r))
            lam = 1.0
            for _ in range(60):
                # Projected step: fractions, flows and temperatures never fall
                # below a tenth of their current value (index 6 is the duty).
                trial = [a + lam * b if k == 6 else max(a + lam * b, 0.1 * a)
                         for k, (a, b) in enumerate(zip(v, step))]
                try:
                    rt = self.residuals(trial)
                    if math.sqrt(sum(x * x for x in rt)) < base:
                        break
                except (ValueError, OverflowError):
                    pass
                lam *= 0.5
            else:
                raise ArithmeticError(f"{self.case['id']}: line search failed at {base:.3e}")
            v, r = trial, rt
        raise ArithmeticError(f"{self.case['id']}: Newton did not converge")


def reference(case):
    col = Column(case)
    feed = case["feed"]
    # The feed is a single subcooled liquid.
    assert sum(feed["z"][i] * math.exp(ln_gamma(feed["z"], feed["temperature_K"])[i])
               * psat(i, feed["temperature_K"]) for i in range(2)) < feed["pressure_Pa"]
    assert tpd_minimum(feed["z"], feed["temperature_K"]) > -1e-9, case["id"]
    v, iterations = col.solve()
    d, st = col.unpack(v)
    xl, xh, A, Bh, qc = d[0:2], d[2:4], d[4], d[5], d[6] * col.escale
    qr = col.reboiler_duty(v)
    assert qc > 0 and qr > 0, case["id"]
    td, p = col.td, col.p
    # Decanter: distinct liquids, light less dense is checked by the engine;
    # stable tie line and no vapor at T_D.
    assert abs(xl[0] - xh[0]) > 0.05
    decanter_tpd = tpd_minimum(xl, td)
    assert decanter_tpd > -1e-7, (case["id"], decanter_tpd)
    bubble_p = sum(xl[i] * math.exp(ln_gamma(xl, td)[i]) * psat(i, td) for i in range(2))
    assert bubble_p < p, case["id"]
    stages = []
    for k, s in enumerate(st):
        T, x, y, L, V = s[0], s[1:3], s[3:5], s[5], s[6]
        assert L > 0 and V > 0 and min(x + y) > 0, (case["id"], k)
        tpd = tpd_minimum(x, T)
        assert tpd > -1e-9, (case["id"], k + 2, tpd)
        stages.append(dict(stage=k + 2, temperature_K=T, liquid_mol_s=L, vapor_mol_s=V,
                           x=x, y=y, grid_tpd_minimum=tpd))
    fl, fh = col.phi
    kind = case["reflux"]["kind"]
    if kind.startswith("heavy"):
        distillate = ((1 - fh) * Bh, xh)
        decant = (A, xl)
    elif kind.startswith("light"):
        distillate = ((1 - fl) * A, xl)
        decant = (Bh, xh)
    else:
        distillate = ((1 - fl) * A, xl)
        decant = ((1 - fh) * Bh, xh)
    bottoms = (st[-1][5], st[-1][1:3])
    # Overall component and energy balance of the reference.
    for i in range(2):
        out = distillate[0] * distillate[1][i] + decant[0] * decant[1][i] + \
            bottoms[0] * bottoms[1][i]
        assert abs(out - col.F * col.z[i]) < 1e-11, (case["id"], i)
    h_out = (distillate[0] * h_liquid(distillate[1], td) + decant[0] * h_liquid(decant[1], td)
             + bottoms[0] * h_liquid(bottoms[1], st[-1][0]))
    assert abs(h_out - col.hF - (qr - qc)) < 1e-6 * col.escale, case["id"]
    return dict(case, newton_iterations=iterations,
                decanter=dict(light_mol_s=A, x_light=xl, heavy_mol_s=Bh, x_heavy=xh,
                              grid_tpd_minimum=decanter_tpd, bubble_pressure_Pa=bubble_p),
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
        schema="a17-three-phase-column-reference-v1",
        issue=188,
        components=IDS,
        method=("Newton on explicit stage temperatures, mole fractions and phase flows "
                "with summation equations; decanter LLE at the specified temperature; "
                "gamma-phi NRTL with ideal-gas vapor and Clausius-Clapeyron liquid "
                "enthalpies (complex-step derivatives); composition-grid TPD and "
                "bubble-pressure acceptance. Python standard library, no Rust calls "
                "or Rust outputs."),
        scope=("Implementation agreement for the NRTL-IG package on the ChemSep/DECHEMA "
               "water/1-butanol pair. Not experimental accuracy: the pair is a 1 atm "
               "VLE regression (see docs/vlle-status.md)."),
        generator=dict(python=platform.python_version(),
                       sha256=sha256(Path(__file__).resolve())),
        inputs={p.name: sha256(p) for p in EXTRACTS},
        tolerances=TOLERANCES,
        cases=[reference(case) for case in CASES],
    )
    Path(args.output).write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    for case in out["cases"]:
        print(case["id"], case["newton_iterations"], "iterations; T",
              [round(s["temperature_K"], 4) for s in case["profile"]],
              "Qc", round(case["condenser_duty_W"], 3), "Qr", round(case["reboiler_duty_W"], 3))


if __name__ == "__main__":
    main()
