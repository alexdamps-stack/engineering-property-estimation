"""Independent references for product distillation curves and cut-point
specifications of the petroleum fractionator (#190).

Python standard library only; no Rust code or output is used.

TBP curve of a product (definition of `docs/a17-columns-petroleum.md`,
"Distillation curves"): the components are ordered by normal boiling point;
component i boils over [(Tb[i-1] + Tb[i])/2, (Tb[i] + Tb[i+1])/2], the end
intervals extending half the neighbouring spacing beyond the end components
(a single component boils at its Tb); within its interval the temperature
rises linearly with the cumulative standard liquid volume, and its share of
the curve is its share of the product's standard liquid volume. Intervals are
fixed by the whole component list, including absent components. The 0 % and
100 % points are the ends of the lightest and heaviest components present.
Database components (water) are left out; their mole fraction is reported.

ASTM D86 from TBP: the inverse of the API method TBP = a D86^b (degR) of
M. R. Riazi and T. E. Daubert (Hydrocarbon Processing 59(3), 1980; Ind. Eng.
Chem. Process Des. Dev. 25, 1986), constants as tabulated in Fahim,
Al-Sahhaf and Elkilani, Fundamentals of Petroleum Refining (2010), Table 3.2,
typed here independently of the Rust table. D86 between the tabulated 0, 10,
30, 50, 70, 90 and 95 vol% is linear in volume percent.

Inputs come from committed, Rust-independent references only:
`validation/columns/petroleum-v1.json` (teqp MESH solutions with their
compositions, boiling points and standard liquid volumes) and the printed
example E3.1 of `validation/pseudo-components/fahim-2010-worked-examples.json`.

The cut-point case takes `crude-two-draws-reboiled-PR` and replaces the molar
flows of its side products by cut points read from the reference products:
TBP 95 vol% of side draw 1 and D86 50 vol% of side draw 2. The column that
meets those cut points is the reference column itself, so the expected side
product flows, profile and products are the reference values. Each cut-point
and product-curve temperature carries a budget propagated from the
composition budget of petroleum-v1.json: sum_i |dT/dx_i| * 2e-7, by central
differences; the excluded (water) mole fraction carries 2e-7 per excluded
component. The curve functions themselves (synthetic cases, products evaluated
at the reference compositions) agree to 1e-9 K.

    python validation/columns/generate-petroleum-curve-references.py --output /tmp/petroleum-curves-candidate.json
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

PATH = Path(__file__).resolve()
ROOT = PATH.parents[2]
COLUMNS = ROOT / "validation/columns/petroleum-v1.json"
FAHIM = ROOT / "validation/pseudo-components/fahim-2010-worked-examples.json"

# Fahim et al. (2010) Table 3.2: vol%, a, b of TBP = a D86^b in degR.
API = [(0, 0.9167, 1.0019), (10, 0.5277, 1.0900), (30, 0.7429, 1.0425), (50, 0.8920, 1.0176),
       (70, 0.8705, 1.0226), (90, 0.9490, 1.0110), (95, 0.8008, 1.0355)]
REPORTED = [5.0, 10.0, 30.0, 50.0, 70.0, 90.0, 95.0]
COMPOSITION_BUDGET = 2e-7
# Budgets fixed before any Rust comparison ([absolute, relative]). Curve values
# are closed-form: 1e-9 K covers floating-point evaluation order.
TOLERANCES = {"curve_temperature_K": [1e-9, 0], "excluded_mole_fraction": [1e-14, 0]}


def tbp_knots(components):
    """Knots [(V0, T0, V1, T1)] of the present components; components = [(Tb, volume)]."""
    order = sorted(range(len(components)), key=lambda i: components[i][0])
    tb = [components[i][0] for i in order]
    vol = [components[i][1] for i in order]
    n = len(tb)
    if n == 1:
        edges = [tb[0], tb[0]]
    else:
        edges = ([tb[0] - (tb[1] - tb[0])/2] + [(tb[i] + tb[i + 1])/2 for i in range(n - 1)]
                 + [tb[-1] + (tb[-1] - tb[-2])/2])
    assert edges[0] > 0
    total = math.fsum(vol)
    assert total > 0
    knots, run = [], 0.0
    for i in range(n):
        if vol[i] > 0:
            start = 100*run/total
            run += vol[i]
            knots.append([start, edges[i], 100*run/total, edges[i + 1]])
    knots[-1][2] = 100.0
    return knots


def tbp_at(knots, x):
    assert 0 <= x <= 100
    for v0, t0, v1, t1 in knots:
        if x <= v1:
            share = min(max((x - v0)/(v1 - v0), 0.0), 1.0)
            return t0 + share*(t1 - t0)
    v0, t0, v1, t1 = knots[-1]
    return t1


def d86_points(knots):
    out = []
    for v, a, b in API:
        t_r = 1.8*tbp_at(knots, v)
        out.append((v, (t_r/a)**(1/b)/1.8))
    return out


def d86_at(knots, x):
    assert 0 <= x <= 95
    pts = d86_points(knots)
    for (v0, t0), (v1, t1) in zip(pts, pts[1:]):
        if x <= v1:
            return t0 + (t1 - t0)*(x - v0)/(v1 - v0)
    raise AssertionError


def product_components(case, z):
    """[(Tb, x V)] of the pseudo-components and the excluded mole fraction."""
    comps, excluded = [], 0.0
    for c, x in zip(case["components"], z):
        if "standard_liquid_volume_m3_per_mol" in c:
            comps.append((c["normal_boiling_point_K"], max(x, 0.0)*c["standard_liquid_volume_m3_per_mol"]))
        else:
            excluded += x
    return comps, excluded


def curve_record(components):
    knots = tbp_knots(components)
    return dict(tbp_K=[[v, tbp_at(knots, v)] for v in [0.0] + REPORTED + [100.0]],
                d86_K=[[v, d86_at(knots, v)] for v in [0.0] + REPORTED])


def cut_temperature(case, z, basis, x):
    comps, _ = product_components(case, z)
    knots = tbp_knots(comps)
    return tbp_at(knots, x) if basis == "tbp" else d86_at(knots, x)


def propagated_budget(case, z, basis, x):
    """sum_i |dT/dz_i| * composition budget, central differences."""
    total = 0.0
    for i in range(len(z)):
        h = 1e-6*max(z[i], 1e-6)
        up, down = list(z), list(z)
        up[i] += h
        down[i] = max(down[i] - h, 0.0)
        slope = (cut_temperature(case, up, basis, x) - cut_temperature(case, down, basis, x))/(up[i] - down[i])
        total += abs(slope)*COMPOSITION_BUDGET
    return total


def published_example():
    """E3.1 read backwards: the printed TBP points give the printed D86 points."""
    e = next(e for e in json.loads(FAHIM.read_text())["examples"] if e["id"] == "E3.1")
    tbp_c = e["expected"][0]["values"]
    d86_c = e["inputs"]["d86_temperature_degC"]
    computed = []
    for (v, a, b), t in zip(API, tbp_c):
        t_r = 1.8*(t + 273.15)
        computed.append((t_r/a)**(1/b)/1.8 - 273.15)
    # Half a printed digit (0.05 degC) of TBP input times dD86/dTBP <= 1.1,
    # plus half a printed digit of D86.
    atol = 0.05*1.1 + 0.05
    for c, d in zip(computed, d86_c):
        assert abs(c - d) <= atol, (c, d)
    return dict(id="E3.1", tbp_degC=tbp_c, printed_d86_degC=d86_c, computed_d86_degC=computed,
                atol_degC=atol, rationale="TBP printed to 0.1 degC (0.05 x max slope 1.1) plus D86 printed to 0.1 degC")


def synthetic():
    cases = [
        dict(id="three-even", components=[[400.0, 1.0], [420.0, 2.0], [440.0, 1.0]]),
        dict(id="absent-middle", components=[[400.0, 1.0], [420.0, 0.0], [440.0, 1.0]]),
        dict(id="unsorted-uneven", components=[[510.0, 0.3], [350.0, 0.2], [430.0, 0.4], [370.0, 0.1]]),
        dict(id="single", components=[[450.0, 2.0]]),
    ]
    for c in cases:
        c.update(curve_record([tuple(x) for x in c["components"]]))
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    columns = json.loads(COLUMNS.read_text())
    products = []
    for case in columns["cases"]:
        labels = ["overhead", "bottoms"] + [f"side draw {k + 1}" for k in range(len(case["draws"]))]
        for label, outlet in zip(labels, case["outlets"]):
            comps, excluded = product_components(case, outlet["z"])
            record = curve_record(comps)
            z = outlet["z"]
            products.append(dict(
                case=case["id"], product=label, excluded_mole_fraction=excluded,
                excluded_budget=COMPOSITION_BUDGET*(len(z) - len(comps)),
                tbp_K=[p for p in record["tbp_K"] if p[0] in REPORTED],
                d86_K=[p for p in record["d86_K"] if p[0] in REPORTED],
                tbp_budget_K=[propagated_budget(case, z, "tbp", v) for v in REPORTED],
                d86_budget_K=[propagated_budget(case, z, "astm_d86", v) for v in REPORTED]))
    base = next(c for c in columns["cases"] if c["id"] == "crude-two-draws-reboiled-PR")
    cut_points = []
    for k, (basis, x) in enumerate([("tbp", 95.0), ("astm_d86", 50.0)]):
        z = base["outlets"][2 + k]["z"]
        cut_points.append(dict(side_draw=k + 1, basis=basis, volume_percent=x,
                               temperature_K=cut_temperature(base, z, basis, x),
                               temperature_budget_K=propagated_budget(base, z, basis, x),
                               molar_flow_mol_s=base["outlets"][2 + k]["molar_flow_mol_s"]))
    result = dict(
        schema="a17-petroleum-curves-reference-v1",
        python="standard library only",
        sha256={str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (PATH, COLUMNS, FAHIM)},
        api_d86_tbp_constants=API,
        reported_volume_percent=REPORTED,
        tolerances=TOLERANCES,
        published_example=published_example(),
        synthetic_curves=synthetic(),
        products=products,
        cut_point_case=dict(base_case=base["id"], cut_points=cut_points,
                            note="the reference column meets these cut points; compare with "
                                 "petroleum-v1.json budgets and each cut point's temperature budget"),
    )
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
