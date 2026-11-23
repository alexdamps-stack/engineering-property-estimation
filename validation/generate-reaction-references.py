"""Independent A17 reactor references; never invokes the Rust engine.

Cantera evaluates the original NASA7 correlations, including the formation
integration constant a6, for an ideal-gas mixture defined from the pinned
source extract. Element counts are parsed here from the source formulas.
Reactor outlet flows are computed by plain arithmetic from the case
definitions below, and each case reports the ideal-gas part of the reactor
duty, H_out(T_out) - H_in(T_in), with formation enthalpies included.

The Rust test removes the PR/PR-BM residual enthalpy of the inlet and outlet
states before comparing, so this fixture checks stoichiometry, mass yields,
formation enthalpies and the reacting energy balance, not the cubic departure
(covered by the Clapeyron and A3 fixtures).

Generate to an explicit candidate path and review before adoption.
"""
import argparse
import hashlib
import importlib.metadata
import json
import platform
import re
from pathlib import Path

import cantera as ct

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/source-extracts/components.provisional.json"
DATA = {c["id"]: c for c in json.loads(SOURCE.read_text())["components"]}
COMPONENTS = ["methane", "ethane", "methanol", "oxygen", "nitrogen", "carbon-dioxide", "water"]
# Flows and duties are compared relative to the inlet flow and to sum |n h|.
TOLERANCES = {"component_flow_relative": 1e-12, "duty_relative": 1e-9,
              "formation_enthalpy_J_per_mol": 1e-6}

CASES = [
    {"id": "methane-combustion-conversion", "kind": "conversion",
     "feed": {"temperature_K": 450.0, "pressure_Pa": 1.5e5, "flow_mol_per_s": 2.0,
              "mole_fractions": {"methane": 0.08, "oxygen": 0.2, "nitrogen": 0.72}},
     "reactions": [{"coefficients": {"methane": -1, "oxygen": -2, "carbon-dioxide": 1, "water": 2},
                    "conversion": {"key": "methane", "fraction": 0.75}}],
     "outlet_temperature_K": 560.0, "pressure_drop_Pa": 2.0e4},
    {"id": "sequential-ethane-methane-combustion", "kind": "conversion",
     "feed": {"temperature_K": 420.0, "pressure_Pa": 2.0e5, "flow_mol_per_s": 1.0,
              "mole_fractions": {"methane": 0.05, "ethane": 0.03, "oxygen": 0.25, "nitrogen": 0.67}},
     "reactions": [{"coefficients": {"ethane": -1, "oxygen": -3.5, "carbon-dioxide": 2, "water": 3},
                    "conversion": {"key": "ethane", "fraction": 1.0}},
                   {"coefficients": {"methane": -1, "oxygen": -2, "carbon-dioxide": 1, "water": 2},
                    "extent_mol_per_s": 0.02}],
     "outlet_temperature_K": 590.0, "pressure_drop_Pa": 0.0},
    {"id": "methanol-mass-yield", "kind": "yield",
     "feed": {"temperature_K": 470.0, "pressure_Pa": 1.2e5, "flow_mol_per_s": 1.5,
              "mole_fractions": {"methanol": 0.2, "water": 0.1, "nitrogen": 0.7}},
     "mass_yields": {"methane": 1.0, "carbon-dioxide": 2.5, "water": 1.2},
     "inerts": ["nitrogen"],
     "outlet_temperature_K": 520.0, "pressure_drop_Pa": 1.0e4},
]


def elements(formula):
    counts = {}
    for symbol, digits in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        counts[symbol] = counts.get(symbol, 0) + (int(digits) if digits else 1)
    return counts


def gas():
    species = []
    for i in COMPONENTS:
        low, high = DATA[i]["ideal_gas"]["coefficients"]
        assert DATA[i]["ideal_gas"]["temperature_breaks_K"] == [200, 1000, 6000]
        s = ct.Species(i, elements(DATA[i]["formula"]))
        s.thermo = ct.NasaPoly2(200, 6000, 100000, [1000, *high, *low])
        species.append(s)
    return ct.Solution(thermo="ideal-gas", species=species)


def enthalpy_flow(g, flows, t):
    """Total ideal-gas enthalpy flow in W (formation basis), flows in mol/s."""
    total = sum(flows.values())
    g.TPX = t, 1e5, {k: v / total for k, v in flows.items() if v > 0}
    return total * g.enthalpy_mole / 1000.0


def scale(g, flows, t):
    total = sum(flows.values())
    return sum(abs(n * g.species(k).thermo.h(t) / 1000.0) for k, n in flows.items()), total


def run(g, case):
    feed = case["feed"]
    inlet = {i: feed["flow_mol_per_s"] * feed["mole_fractions"].get(i, 0.0) for i in COMPONENTS}
    flows = dict(inlet)
    extents = []
    if case["kind"] == "conversion":
        for r in case["reactions"]:
            nu = r["coefficients"]
            for e in {e for i in nu for e in elements(DATA[i]["formula"])}:
                assert sum(v * elements(DATA[i]["formula"]).get(e, 0) for i, v in nu.items()) == 0
            if "conversion" in r:
                key = r["conversion"]["key"]
                xi = r["conversion"]["fraction"] * flows[key] / -nu[key]
            else:
                xi = r["extent_mol_per_s"]
            for i, v in nu.items():
                flows[i] += v * xi
            assert all(n >= -1e-15 for n in flows.values()), case["id"]
            flows = {i: max(n, 0.0) for i, n in flows.items()}
            extents.append(xi)
    else:
        mass = {i: n * DATA[i]["constants"]["molar_mass_kg_per_mol"] for i, n in inlet.items()}
        reacting = sum(m for i, m in mass.items() if i not in case["inerts"])
        weights = case["mass_yields"]
        total = sum(weights.values())
        for i in COMPONENTS:
            if i not in case["inerts"]:
                flows[i] = reacting * weights.get(i, 0.0) / total / DATA[i]["constants"]["molar_mass_kg_per_mol"]
    h_in = enthalpy_flow(g, inlet, feed["temperature_K"])
    h_out = enthalpy_flow(g, flows, case["outlet_temperature_K"])
    generation = {i: flows[i] - inlet[i] for i in COMPONENTS}
    standard = sum(n * g.species(i).thermo.h(298.15) / 1000.0 for i, n in generation.items())
    atoms = {}
    for i, n in generation.items():
        for e, a in elements(DATA[i]["formula"]).items():
            atoms[e] = atoms.get(e, 0.0) + n * a
    return {**case, "outlet_flows_mol_per_s": flows, "extents_mol_per_s": extents,
            "standard_reaction_enthalpy_W": standard,
            "ideal_gas_duty_W": h_out - h_in,
            "duty_scale_W": scale(g, inlet, feed["temperature_K"])[0] + scale(g, flows, case["outlet_temperature_K"])[0],
            "element_generation_mol_per_s": atoms}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    g = gas()
    formation = {i: g.species(i).thermo.h(298.15) / 1000.0 for i in DATA if i in COMPONENTS}
    for i in DATA:
        if i not in COMPONENTS:
            low, high = DATA[i]["ideal_gas"]["coefficients"]
            formation[i] = ct.NasaPoly2(200, 6000, 100000, [1000, *high, *low]).h(298.15) / 1000.0
    output = {"schema": "a17-reactor-reference-v1", "generator": "validation/generate-reaction-references.py",
              "python": platform.python_version(), "packages": {"cantera": importlib.metadata.version("cantera")},
              "input_sha256": {"data/source-extracts/components.provisional.json": hashlib.sha256(SOURCE.read_bytes()).hexdigest()},
              "independence": "Cantera NASA7 ideal-gas mixture with formation constants; Python formula parsing and stoichiometric arithmetic; no Rust calls",
              "scope": "Ideal-gas reacting energy and material balances of fixed conversion/yield cases; the PR/PR-BM departure is removed by the Rust test and qualified separately",
              "components": COMPONENTS, "tolerances": TOLERANCES,
              "formation_enthalpy_J_per_mol": dict(sorted(formation.items())),
              "cases": [run(g, c) for c in CASES]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {len(output['cases'])} reactor cases and {len(formation)} formation enthalpies to {args.output}")


if __name__ == "__main__":
    main()
