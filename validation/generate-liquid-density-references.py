"""Independent references for the engineering liquid density (#228). Review-only.

Never run automatically to repair tests: a regenerated fixture is a reviewed
change. Writes `validation/fixtures/liquid-density-v1.json` with two separate
kinds of evidence:

* implementation verification: the COSTALD saturated-liquid correlation
  (Hankinson and Thomson 1979), its mixing rules and the COSTALD-Tait
  compressed-liquid correction with the generalized Riedel saturation pressure
  (Thomson, Brobst and Hankinson 1982; API TDB 6A3.4), evaluated by the
  independent `chemicals` 1.5.2 implementation with the engine's transcribed
  parameters and critical temperatures, plus domain-failure cases;
* experimental validation: saturated and compressed liquid densities of the
  reference multiparameter equations of state in CoolProp 7.2.0, which are
  correlations of measured densities with stated uncertainties of about
  0.01-0.5 % in the liquid. Pure fluids and hydrocarbon/nitrogen bubble-point
  liquids (CoolProp's GERG-type mixture departure functions).

Neither reference shares code with the Rust engine. Shared inputs: the
critical temperatures of the embedded component database (CoolProp 7.2 cubic
constants) and the COSTALD parameters. Budgets are declared below, before
comparison, and are written into the fixture.

    python -I validation/generate-liquid-density-references.py --output validation/fixtures/liquid-density-v1.json

Requires `chemicals==1.5.2` and `CoolProp==7.2.0` (validation/liquid-density-requirements.txt).
"""
import argparse
import hashlib
import json
import math
import pathlib
import platform

import chemicals
import CoolProp
from chemicals.dippr import EQ105, EQ106
from chemicals.volume import COSTALD, COSTALD_compressed, COSTALD_mixture_parameters
from CoolProp.CoolProp import AbstractState, PropsSI
import CoolProp.CoolProp as CP

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXTRACTS = [
    ROOT / "data/source-extracts/components.provisional.json",
    ROOT / "data/source-extracts/components-reacting.provisional.json",
    ROOT / "data/source-extracts/liquid-density-costald.provisional.json",
    ROOT / "data/source-extracts/liquid-density-chemsep.provisional.json",
]

# Declared engine domain (must match thermo-models/src/liquid_density.rs).
REDUCED_TEMPERATURE_RANGE = (0.25, 0.95)
MAXIMUM_PRESSURE_PA = 69e6
# Water/methanol route (restricted after the comparison: the generalized Tait
# correction under-predicts their compression by more than 2 % above ~10 MPa).
IDEAL_MIXING_MAXIMUM_PRESSURE_PA = 10e6

# Budgets, declared before any comparison.
IMPLEMENTATION_BUDGET = {
    "quantity": "molar volume and pseudo saturation pressure",
    "atol": 0.0,
    "rtol": 1e-11,
    "rationale": "Same closed-form expressions in double precision; R differs from fluids.constants by 2e-12 relative",
}
EXPERIMENTAL_BUDGET = [
    {
        "reduced_temperature_max": 0.90,
        "max_abs_relative_deviation": 0.02,
        "rationale": "Hankinson-Thomson report about 0.4 % average deviation for pure saturated liquids and Thomson-Brobst-Hankinson about 1 % for compressed liquids below Tr 0.95; 2 % bounds single points",
    },
    {
        "reduced_temperature_max": 0.95,
        "max_abs_relative_deviation": 0.05,
        "rationale": "Near-critical band 0.90 < Tr <= 0.95 where the corresponding-states fit and the reference EOS both lose accuracy",
    },
]

# CoolProp fluid names for the embedded components.
COOLPROP = {
    "methane": "Methane",
    "ethane": "Ethane",
    "propane": "n-Propane",
    "nitrogen": "Nitrogen",
    "carbon-dioxide": "CarbonDioxide",
    "oxygen": "Oxygen",
    "water": "Water",
    "methanol": "Methanol",
    "hydrogen": "Hydrogen",
    "carbon-monoxide": "CarbonMonoxide",
}
R = 8.31446261815324

# Engine qualification policy (must match thermo-models/src/liquid_density.rs).
# Set from the first comparison against CoolProp: water (about 24 %),
# methanol (2-4 %) and hydrogen (7-11 %) exceed the declared budgets as pure
# fluids, and CO2 + hydrocarbon bubble-point liquids exceed them (3-9 %).
# Budgets were not widened; the qualified scope was restricted instead.
# Second revision (#228 follow-up): water and methanol use their ChemSep v8.3
# saturated-liquid correlations (DIPPR 106 / 105) instead of COSTALD, with
# ideal (Amagat) mixing of pure saturated volumes whenever they exceed the
# trace fraction. The same 2 % / 5 % budgets apply; they were declared before
# this comparison. Water/methanol mixtures are not validated (unassessed).
QUALIFIED = ["methane", "water", "methanol", "ethane", "propane", "nitrogen", "carbon-dioxide", "oxygen", "carbon-monoxide"]
MIXTURE_QUALIFIED = ["methane", "ethane", "propane", "nitrogen", "oxygen", "carbon-monoxide"]
TRACE = 1e-3


def scope(ids, x):
    """'qualified' or 'unassessed' for a state the engine evaluates."""
    major = [i for i, xi in zip(ids, x) if xi > TRACE]
    if len(major) == 1 or all(i in MIXTURE_QUALIFIED for i in major):
        return "qualified"
    return "unassessed"


def classify(comps, ids, x, t, p):
    """Engine outcome at (T, p): 'refused', 'qualified' or 'unassessed'."""
    if expected_volume(comps, ids, x, t, p) is None:
        return "refused"
    return scope(ids, x)


def chemsep_density(record, t):
    """mol/m3 from the transcribed ChemSep record via chemicals.dippr, or None out of range."""
    lo, hi = record["temperature_range_K"]
    if not lo <= t <= hi:
        return None
    a, b, c, d, e = record["coefficients"]
    if record["equation"] == 105:
        return 1000.0 * EQ105(t, a, b, c, d)
    return 1000.0 * EQ106(t, record["chemsep_critical_temperature_K"], a, b, c, d, e)


def load():
    comps = {}
    for path in EXTRACTS[:2]:
        for c in json.loads(path.read_text(encoding="utf-8"))["components"]:
            comps[c["id"]] = {
                "tc": c["constants"]["critical_temperature_K"],
                "mw": c["constants"]["molar_mass_kg_per_mol"],
            }
    costald = json.loads(EXTRACTS[2].read_text(encoding="utf-8"))
    for c in costald["components"]:
        comps[c["id"]]["vstar"] = c["characteristic_volume_m3_per_mol"]
        comps[c["id"]]["omega"] = c["omega_srk"]
    for c in json.loads(EXTRACTS[3].read_text(encoding="utf-8"))["components"]:
        comps[c["id"]]["chemsep"] = c
    return comps


def riedel_psat(t, tcm, vm, omega):
    """Pseudo saturation pressure of API TDB 6A3.4 (independent transcription)."""
    trm = t / tcm
    log_tr = math.log10(trm)
    alpha = 35.0 - 36.0 / trm - 96.736 * log_tr + trm**6
    beta = log_tr + 0.03721754 * alpha
    log_pr = 5.8031817 * log_tr + 0.07608141 * alpha + omega * 4.86601 * beta
    pcm = (0.291 - 0.080 * omega) * R * tcm / vm
    return 10.0**log_pr * pcm


def expected_volume(comps, ids, x, t, p):
    """chemicals evaluation; returns (saturated, reported, psat, corrected, Tr)
    or None where the engine must return OutsideValidityDomain."""
    if sum(xi for i, xi in zip(ids, x) if i not in QUALIFIED) > TRACE:
        return None
    tcs = [comps[i]["tc"] for i in ids]
    vs = [comps[i]["vstar"] for i in ids]
    ws = [comps[i]["omega"] for i in ids]
    tcm, vm, omega = COSTALD_mixture_parameters(x, tcs, vs, ws)
    trm = t / tcm
    if not REDUCED_TEMPERATURE_RANGE[0] <= trm <= REDUCED_TEMPERATURE_RANGE[1] or p > MAXIMUM_PRESSURE_PA:
        return None
    v_sat = COSTALD(t, tcm, vm, omega)
    if sum(xi for i, xi in zip(ids, x) if "chemsep" in comps[i]) > TRACE:
        if p > IDEAL_MIXING_MAXIMUM_PRESSURE_PA:
            return None
        # Ideal mixing of pure saturated volumes.
        volume = included = excluded = 0.0
        for i, xi in zip(ids, x):
            if xi == 0.0:
                continue
            c = comps[i]
            if "chemsep" in c:
                rho = chemsep_density(c["chemsep"], t)
                if rho is None:
                    return None
                vi = 1.0 / rho
            elif REDUCED_TEMPERATURE_RANGE[0] <= t / c["tc"] <= REDUCED_TEMPERATURE_RANGE[1]:
                vi = COSTALD(t, c["tc"], c["vstar"], c["omega"])
            else:
                excluded += xi
                continue
            volume += xi * vi
            included += xi
        if excluded > TRACE:
            return None
        v_sat = volume / included
    psat = riedel_psat(t, tcm, vm, omega)
    if p > psat:
        pcm = (0.291 - 0.080 * omega) * R * tcm / vm
        v = COSTALD_compressed(t, p, psat, tcm, pcm, omega, v_sat)
        return v_sat, v, psat, True, trm
    return v_sat, v_sat, psat, False, trm


def implementation_cases(comps):
    cases = []
    pure_trs = [0.3, 0.5, 0.7, 0.85, 0.93]
    for cid in COOLPROP:
        tc = comps[cid]["tc"]
        for tr in pure_trs:
            t = round(tr * tc, 6)
            for p in [1e5, 5e6, 3e7]:
                cases.append((f"pure-{cid}-Tr{tr}-p{p:g}", [cid], [1.0], t, p))
    mixtures = [
        (["methane", "ethane", "propane"], [0.2, 0.3, 0.5]),
        (["nitrogen", "methane"], [0.1, 0.9]),
        (["propane", "carbon-dioxide"], [0.7, 0.3]),
        (["methanol", "water"], [0.4, 0.6]),
        (["water", "carbon-dioxide", "hydrogen"], [0.97, 0.02, 0.01]),
        (["water", "carbon-dioxide", "hydrogen"], [0.9992, 0.0004, 0.0004]),
        (["water", "methanol", "propane"], [0.5, 0.3, 0.2]),
        (["propane", "methane", "carbon-monoxide", "oxygen"], [0.85, 0.05, 0.05, 0.05]),
        (["propane", "water", "hydrogen"], [0.9992, 0.0004, 0.0004]),
        (["propane", "water", "hydrogen"], [0.998, 0.001, 0.001]),
    ]
    for ids, x in mixtures:
        tcs = [comps[i]["tc"] for i in ids]
        vs = [comps[i]["vstar"] for i in ids]
        ws = [comps[i]["omega"] for i in ids]
        tcm, _, _ = COSTALD_mixture_parameters(x, tcs, vs, ws)
        for tr in [0.4, 0.6, 0.8, 0.92]:
            t = round(tr * tcm, 6)
            for p in [2e5, 1e7, 6e7]:
                cases.append((f"mix-{'+'.join(ids)}-{x}-Tr{tr}-p{p:g}", ids, x, t, p))
    out = []
    for case_id, ids, x, t, p in cases:
        expected = expected_volume(comps, ids, x, t, p)
        base = {"id": case_id, "components": ids, "mole_fractions": x, "temperature_K": t, "pressure_Pa": p}
        if expected is None:
            out.append(base | {"engine_scope": "refused"})
            continue
        v_sat, v, psat, corrected, trm = expected
        mw = sum(xi * comps[i]["mw"] for xi, i in zip(x, ids))
        out.append(
            base
            | {
                "engine_scope": scope(ids, x),
                "reduced_temperature": trm,
                "saturated_molar_volume_m3_per_mol": v_sat,
                "pseudo_saturation_pressure_Pa": psat,
                "pressure_corrected": corrected,
                "molar_volume_m3_per_mol": v,
                "mass_density_kg_per_m3": mw / v,
            }
        )
    return out


def failure_cases(comps):
    tc = comps["propane"]["tc"]
    return [
        {"id": "propane-above-Tr-0.95", "components": ["propane"], "mole_fractions": [1.0],
         "temperature_K": 0.96 * tc, "pressure_Pa": 5e6, "error": "OutsideValidityDomain"},
        {"id": "propane-below-Tr-0.25", "components": ["propane"], "mole_fractions": [1.0],
         "temperature_K": 0.24 * tc, "pressure_Pa": 1e5, "error": "OutsideValidityDomain"},
        {"id": "propane-above-pressure-limit", "components": ["propane"], "mole_fractions": [1.0],
         "temperature_K": 0.6 * tc, "pressure_Pa": 7e7, "error": "OutsideValidityDomain"},
        {"id": "methane-supercritical", "components": ["methane"], "mole_fractions": [1.0],
         "temperature_K": 250.0, "pressure_Pa": 5e6, "error": "OutsideValidityDomain"},
    ]


def band(tr):
    for b in EXPERIMENTAL_BUDGET:
        if tr <= b["reduced_temperature_max"]:
            return b["max_abs_relative_deviation"]
    return None


def experimental_pure(comps):
    points = []
    for cid, fluid in COOLPROP.items():
        tc_db = comps[cid]["tc"]
        t_lo = max(PropsSI("Ttriple", fluid) + 0.5, REDUCED_TEMPERATURE_RANGE[0] * tc_db + 1e-3)
        t_hi = min(REDUCED_TEMPERATURE_RANGE[1] * tc_db, PropsSI("Tcrit", fluid) - 0.5) - 1e-3
        n = 12
        for k in range(n):
            t = round(t_lo + (t_hi - t_lo) * k / (n - 1), 4)
            rho_ref = PropsSI("Dmolar", "T", t, "Q", 0, fluid)
            psat_ref = PropsSI("P", "T", t, "Q", 0, fluid)
            points.append({"id": f"{cid}-sat-{k}", "kind": "saturated", "components": [cid],
                           "mole_fractions": [1.0], "temperature_K": t, "pressure_Pa": psat_ref,
                           "reference_molar_density_mol_per_m3": rho_ref,
                           "engine_scope": classify(comps, [cid], [1.0], t, psat_ref), "budget": band(t / tc_db)})
        for tr in [0.5, 0.7, 0.85]:
            t = round(tr * tc_db, 4)
            if not (t_lo <= t <= t_hi):
                continue
            psat_ref = PropsSI("P", "T", t, "Q", 0, fluid)
            for dp in [1e6, 5e6, 2e7, 5e7]:
                p = psat_ref + dp
                if p > MAXIMUM_PRESSURE_PA:
                    continue
                try:
                    rho_ref = PropsSI("Dmolar", "T", t, "P", p, fluid)
                except ValueError:
                    continue  # outside the reference EOS range (recorded by absence)
                points.append({"id": f"{cid}-compressed-Tr{tr}-dp{dp:g}", "kind": "compressed",
                               "components": [cid], "mole_fractions": [1.0], "temperature_K": t,
                               "pressure_Pa": p, "reference_molar_density_mol_per_m3": rho_ref,
                               "engine_scope": classify(comps, [cid], [1.0], t, p), "budget": band(tr)})
    return points


def experimental_mixtures(comps):
    """Bubble-point liquids of CoolProp HEOS mixtures (Kunz-Wagner departure functions)."""
    mixtures = [
        (["methane", "ethane"], [0.5, 0.5], [120.0, 150.0, 180.0]),
        (["methane", "propane"], [0.3, 0.7], [130.0, 170.0, 210.0]),
        (["ethane", "propane"], [0.4, 0.6], [200.0, 240.0, 280.0, 320.0]),
        (["nitrogen", "methane"], [0.1, 0.9], [100.0, 120.0, 140.0]),
        (["methane", "ethane", "propane"], [0.85, 0.10, 0.05], [115.0, 130.0, 145.0]),
        (["propane", "carbon-dioxide"], [0.8, 0.2], [230.0, 260.0, 290.0]),
        (["ethane", "carbon-dioxide"], [0.7, 0.3], [220.0, 250.0, 270.0]),
        (["nitrogen", "oxygen"], [0.79, 0.21], [80.0, 90.0, 100.0]),
        (["methane", "carbon-monoxide"], [0.9, 0.1], [110.0, 130.0, 150.0]),
        (["nitrogen", "methane", "ethane", "propane"], [0.02, 0.88, 0.07, 0.03], [115.0, 135.0, 155.0]),
    ]
    points = []
    for ids, x, temperatures in mixtures:
        tcs = [comps[i]["tc"] for i in ids]
        vs = [comps[i]["vstar"] for i in ids]
        ws = [comps[i]["omega"] for i in ids]
        tcm, _, _ = COSTALD_mixture_parameters(x, tcs, vs, ws)
        for t in temperatures:
            state = AbstractState("HEOS", "&".join(COOLPROP[i] for i in ids))
            state.set_mole_fractions(x)
            try:
                state.update(CP.QT_INPUTS, 0.0, t)
            except ValueError as e:
                raise SystemExit(f"CoolProp bubble point failed for {ids} at {t} K: {e}")
            points.append({"id": f"{'+'.join(ids)}-bubble-{t:g}", "kind": "bubble",
                           "components": ids, "mole_fractions": x, "temperature_K": t,
                           "pressure_Pa": state.p(),
                           "reference_molar_density_mol_per_m3": state.rhomolar(),
                           "engine_scope": classify(comps, ids, x, t, state.p()), "budget": band(t / tcm)})
    return points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    output = pathlib.Path(parser.parse_args().output)
    if chemicals.__version__ != "1.5.2" or CoolProp.__version__ != "7.2.0":
        raise SystemExit("requires chemicals 1.5.2 and CoolProp 7.2.0")
    comps = load()
    experimental = experimental_pure(comps) + experimental_mixtures(comps)
    if any(p["budget"] is None for p in experimental):
        raise SystemExit("experimental point outside the declared reduced-temperature bands: " + ", ".join(p["id"] for p in experimental if p["budget"] is None))
    fixture = {
        "schema": "liquid-density-reference-v1",
        "generator": "validation/generate-liquid-density-references.py",
        "references": {
            "implementation": f"chemicals {chemicals.__version__} (MIT) COSTALD, COSTALD_mixture_parameters, COSTALD_mixture_compressed; Riedel pseudo saturation pressure transcribed independently in the generator",
            "experimental": f"CoolProp {CoolProp.__version__} HEOS reference equations of state (saturated, compressed and bubble-point liquid molar density)",
            "python": platform.python_version(),
        },
        "inputs_sha256": {
            str(p.relative_to(ROOT)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in EXTRACTS
        },
        "gas_constant_J_per_mol_K": R,
        "declared_domain": {
            "reduced_temperature": list(REDUCED_TEMPERATURE_RANGE),
            "maximum_pressure_Pa": MAXIMUM_PRESSURE_PA,
            "ideal_mixing_maximum_pressure_Pa": IDEAL_MIXING_MAXIMUM_PRESSURE_PA,
            "mixing_rule": "Hankinson-Thomson COSTALD mixing; when water and methanol exceed the trace fraction, ideal (Amagat) mixing of pure saturated volumes (ChemSep DIPPR 105/106 for water and methanol, pure COSTALD for the others; components without a pure liquid volume at T may total at most the trace fraction and are left out); COSTALD-Tait with the Hankinson-Thomson mixture parameters in both cases",
            "pressure_rule": "Tait correction only above the generalized Riedel pseudo saturation pressure; otherwise the saturated value is reported with pressure_corrected=false",
        },
        "qualification_policy": {
            "qualified_components": QUALIFIED,
            "mixture_qualified_components": MIXTURE_QUALIFIED,
            "trace_mole_fraction": TRACE,
            "rule": "refused (OutsideValidityDomain) when non-qualified components exceed the trace fraction in total; qualified for a single component above trace or when every component above trace is mixture-qualified; otherwise unassessed. Budgets apply to qualified points only; other points are retained as evidence for the restriction.",
        },
        "tolerances": {"implementation": IMPLEMENTATION_BUDGET, "experimental": EXPERIMENTAL_BUDGET},
        "implementation": implementation_cases(comps),
        "failures": failure_cases(comps),
        "experimental": experimental,
    }
    output.write_text(json.dumps(fixture, indent=1) + "\n", encoding="utf-8")
    # Review summary (not part of the fixture): COSTALD deviation per data set.
    summary = {}
    for p in experimental:
        expected = expected_volume(comps, p["components"], p["mole_fractions"],
                                   p["temperature_K"], p["pressure_Pa"])
        if expected is None:
            continue
        v = expected[1]
        dev = 1.0 / v / p["reference_molar_density_mol_per_m3"] - 1.0
        key = ("+".join(p["components"]), p["kind"], p["engine_scope"])
        s = summary.setdefault(key, [0, 0.0, 0.0, 0])
        s[0] += 1
        s[1] += abs(dev)
        s[2] = max(s[2], abs(dev))
        s[3] += abs(dev) > p["budget"]
    for (name, kind, sc), (n, total, worst, failed) in sorted(summary.items()):
        print(f"{name:40s} {kind:10s} {sc:10s} n={n:2d} AAD={100 * total / n:6.3f}% max={100 * worst:6.3f}% over budget={failed}")


if __name__ == "__main__":
    main()
