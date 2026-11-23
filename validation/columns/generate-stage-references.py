"""Independent equilibrium-stage (MESH) references for the columns group; no Rust calls.

Two countercurrent cascades with prescribed stage flows, constant equilibrium
ratios K_ij and constant stage molar enthalpies:

- `distillation`: five stages, total condenser (liquid distillate drawn from
  stage 1, no overhead vapor), saturated liquid feed on stage 3, reboiler on
  stage 5, three components.
- `absorber`: four stages, lean liquid on stage 1, rich gas on stage 4, a vapor
  side draw on stage 2, two components.

Liquid rates follow from the prescribed vapor rates, feeds and draws by the
total stage balances. For each component the stage component balances with
y_ij = K_ij x_ij (Seader, Henley and Roper, "Separation Process Principles",
3rd ed., Eqs. 10-1 to 10-10) form one linear system. It is solved here by
exact rational Gauss-Jordan elimination on the full matrix, independently of
the Thomas algorithm in the Rust engine. The stage duties then close each
stage energy balance with duty positive into the stage. Values are written as
decimal strings rounded from exact fractions.

Python standard library only. Run manually:

    python validation/columns/generate-stage-references.py --output <candidate.json>

Review the candidate before replacing validation/columns/stages-v1.json.
"""
import argparse
import json
from fractions import Fraction as Fr


def F(v):
    return Fr(v)


CASES = {
    "distillation": {
        "description": "Total condenser, saturated liquid feed on stage 3, reboiler on stage 5.",
        "components": 3,
        "vapor": ["0", "120", "120", "120", "120"],
        "liquid_draw": ["40", "0", "0", "0", "0"],
        "vapor_draw": ["0", "0", "0", "0", "0"],
        "feeds": [
            {"stage": 3, "flow": "100", "z": ["3/10", "2/5", "3/10"], "h": "-41000"}
        ],
        "k": [
            ["21/10", "1", "2/5"],
            ["5/2", "6/5", "1/2"],
            ["3", "7/5", "3/5"],
            ["7/2", "8/5", "7/10"],
            ["4", "9/5", "4/5"],
        ],
        "h_liquid": ["-42000", "-41500", "-41000", "-40200", "-39500"],
        "h_vapor": ["-12000", "-11500", "-11000", "-10400", "-9800"],
    },
    "absorber": {
        "description": "Lean liquid on stage 1, rich gas on stage 4, vapor side draw on stage 2.",
        "components": 2,
        "vapor": ["90", "95", "98", "100"],
        "liquid_draw": ["0", "0", "0", "0"],
        "vapor_draw": ["0", "2", "0", "0"],
        "feeds": [
            {"stage": 1, "flow": "50", "z": ["1/50", "49/50"], "h": "-48000"},
            {"stage": 4, "flow": "100", "z": ["9/10", "1/10"], "h": "-9000"},
        ],
        "k": [["8", "1/20"], ["15/2", "1/25"], ["7", "3/100"], ["13/2", "1/40"]],
        "h_liquid": ["-47800", "-47500", "-47100", "-46600"],
        "h_vapor": ["-9500", "-9300", "-9200", "-9100"],
    },
}


def solve_exact(a, b):
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if m[r][col] != 0)
        m[col], m[pivot] = m[pivot], m[col]
        p = m[col][col]
        m[col] = [v / p for v in m[col]]
        for r in range(n):
            if r != col and m[r][col] != 0:
                f = m[r][col]
                m[r] = [vr - f * vc for vr, vc in zip(m[r], m[col])]
    return [m[i][n] for i in range(n)]


def build(case):
    n = len(case["vapor"])
    c = case["components"]
    V = [F(v) for v in case["vapor"]]
    U = [F(v) for v in case["liquid_draw"]]
    W = [F(v) for v in case["vapor_draw"]]
    K = [[F(v) for v in row] for row in case["k"]]
    hl = [F(v) for v in case["h_liquid"]]
    hv = [F(v) for v in case["h_vapor"]]
    feed = [F(0)] * n
    feed_comp = [[F(0)] * c for _ in range(n)]
    feed_h = [F(0)] * n
    for f in case["feeds"]:
        j = f["stage"] - 1
        flow = F(f["flow"])
        feed[j] += flow
        for i, z in enumerate(f["z"]):
            feed_comp[j][i] += flow * F(z)
        feed_h[j] += flow * F(f["h"])

    # Total balance of stage j: L_{j-1} + V_{j+1} + F_j = L_j + U_j + V_j + W_j.
    L = []
    for j in range(n):
        above = L[j - 1] if j > 0 else F(0)
        below = V[j + 1] if j + 1 < n else F(0)
        L.append(above + below + feed[j] - U[j] - V[j] - W[j])
        assert L[j] >= 0, f"negative liquid rate on stage {j + 1}"

    x = [[F(0)] * c for _ in range(n)]
    for i in range(c):
        a = [[F(0)] * n for _ in range(n)]
        b = [F(0)] * n
        for j in range(n):
            if j > 0:
                a[j][j - 1] = L[j - 1]
            a[j][j] = -(L[j] + U[j] + (V[j] + W[j]) * K[j][i])
            if j + 1 < n:
                a[j][j + 1] = V[j + 1] * K[j + 1][i]
            b[j] = -feed_comp[j][i]
        for j, v in enumerate(solve_exact(a, b)):
            x[j][i] = v
    y = [[K[j][i] * x[j][i] for i in range(c)] for j in range(n)]

    # Duty positive into stage j closes its energy balance.
    Q = []
    for j in range(n):
        inflow = feed_h[j]
        if j > 0:
            inflow += L[j - 1] * hl[j - 1]
        if j + 1 < n:
            inflow += V[j + 1] * hv[j + 1]
        outflow = (L[j] + U[j]) * hl[j] + (V[j] + W[j]) * hv[j]
        Q.append(outflow - inflow)

    feed_total = [sum(feed_comp[j][i] for j in range(n)) for i in range(c)]
    product = [
        V[0] * y[0][i]
        + L[n - 1] * x[n - 1][i]
        + sum(U[j] * x[j][i] + W[j] * y[j][i] for j in range(n))
        for i in range(c)
    ]
    assert feed_total == product, "overall component balance does not close"
    assert sum(Q) + sum(feed_h) == (
        V[0] * hv[0]
        + L[n - 1] * hl[n - 1]
        + sum(U[j] * hl[j] + W[j] * hv[j] for j in range(n))
    ), "overall energy balance does not close"

    def s(v):
        return repr(float(v))

    out = {k: v for k, v in case.items()}
    out["liquid"] = [s(v) for v in L]
    out["x"] = [[s(v) for v in row] for row in x]
    out["y"] = [[s(v) for v in row] for row in y]
    out["duty_w"] = [s(v) for v in Q]
    out["liquid_summation"] = [s(sum(row) - 1) for row in x]
    out["vapor_summation"] = [s(sum(row) - 1) for row in y]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    doc = {
        "schema": "columns-stages-v1",
        "generator": "validation/columns/generate-stage-references.py",
        "method": "exact rational Gauss-Jordan elimination of the stage component balances",
        "conventions": {
            "stages": "numbered 1..N from the top",
            "liquid": "L_j leaves stage j to stage j+1; L_N is the bottom liquid product",
            "vapor": "V_j leaves stage j to stage j-1; V_1 is the top vapor product",
            "draws": "U_j liquid and W_j vapor side draws leave at stage composition",
            "duty": "Q_j in W, positive into the stage",
        },
        "cases": {name: build(case) for name, case in CASES.items()},
    }
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        json.dump(doc, f, indent=1)
        f.write("\n")


if __name__ == "__main__":
    main()
