"""Explicit review-only #224 source transcription. Never run automatically to repair tests.

Adds 1-butanol and the partially miscible water/1-butanol NRTL pair (plus
methanol/1-butanol) for the VLLE/LLE references. Everything comes from the DWSIM
commit already pinned by the A7 and #225 transcriptions:

- ChemSep v8.3 pure-component database (`chemsep1.xml`, Artistic License 2.0):
  critical constants, acentric factor, molar mass, the RPP ideal-gas heat
  capacity (ChemSep equation 100, a quartic polynomial), the ideal-gas
  formation enthalpy and absolute entropy, the vapour pressure (equation 101)
  and the saturated-liquid density (equation 105);
- ChemSep/DECHEMA NRTL interaction parameters at 1 atm (`nrtl.dat`), cross-checked
  against the thermo package's ChemSep table, as for methanol/water.

The RPP polynomial is represented exactly by NASA7 coefficients a1..a5 (cp/R);
a6 and a7 reproduce the ChemSep formation enthalpy and absolute entropy at
298.15 K. No value is fitted. Writes two extension extracts and evaluates no model.
Standard library only.

    python validation/prepare-chemsep-nrtl-source-data.py --output-directory data/source-extracts
"""
import argparse
import hashlib
import json
import math
import pathlib
import urllib.request
import xml.etree.ElementTree as ET

DWSIM = "1abf72d1b6b41d3e9a8cc770d3cc4e8fc76e5766"
THERMO = "2bb466e98439c2395a004e7095d07f9b97a5a0f0"
SPDX = "v3.24.0"
SOURCES = {
    "chemsep_pcd": f"https://raw.githubusercontent.com/DanWBR/dwsim/{DWSIM}/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml",
    "dwsim_chemsep_nrtl": f"https://raw.githubusercontent.com/DanWBR/dwsim/{DWSIM}/DWSIM.Thermodynamics/Assets/nrtl.dat",
    "thermo_chemsep_nrtl": f"https://raw.githubusercontent.com/CalebBell/thermo/{THERMO}/thermo/Interaction%20Parameters/ChemSep/nrtl.json",
    "artistic_license_2_0": f"https://raw.githubusercontent.com/spdx/license-list-data/{SPDX}/text/Artistic-2.0.txt",
}
COMPONENT_VERSION = "2026-11-27.1"
LIQUID_VERSION = "2026-11-27.1"
CAL = 4.184
R = 8.31446261815324
T0 = 298.15
# Database ID, ChemSep library index, CAS, formula.
BUTANOL = ("1-butanol", "1105", "71-36-3", "C4H10O")
# Pairs: (i, j, ChemSep row prefix, thermo CAS key i, thermo CAS key j).
PAIRS = [
    ("water", "1-butanol", "1921;1105;", "7732-18-5", "71-36-3"),
    ("methanol", "1-butanol", "1101;1105;", "67-56-1", "71-36-3"),
]
INFERRED = ("isobaric 1 atm data: interval between pure normal boiling points from the "
            "transcribed vapor-pressure correlations, rounded outward to 0.01 K")


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def compound(root, index):
    for c in root:
        node = c.find("LibraryIndex")
        if node is not None and node.get("value") == index:
            return c
    raise SystemExit(f"ChemSep compound {index} not found")


def value(c, tag, units):
    node = c.find(tag)
    if node is None or node.get("units") != units:
        raise SystemExit(f"unexpected ChemSep {tag}")
    return float(node.get("value"))


def correlation(c, tag, eqno, letters, units):
    node = c.find(tag)
    if node is None or node.get("units") != units or node.find("eqno").get("value") != str(eqno):
        raise SystemExit(f"unexpected ChemSep {tag} correlation")
    coefficients = [float(node.find(k).get("value")) if node.find(k) is not None else 0.0 for k in letters]
    return coefficients, [float(node.find("Tmin").get("value")), float(node.find("Tmax").get("value"))]


def nasa7_from_rpp(coefficients, hf, s0):
    """cp/R = a1 + a2 T + ... + a5 T^4 from cp = A + B T + ... + E T^4 (J/mol/K)."""
    a = [c / R for c in coefficients]
    h = sum(a[k] * T0 ** (k + 1) / (k + 1) for k in range(5))
    s = a[0] * math.log(T0) + sum(a[k] * T0 ** k / k for k in range(1, 5))
    return a + [hf / R - h, s0 / R - s]


def dippr101(c, t):
    a, b, cc, d, e = c
    return math.exp(a + b / t + cc * math.log(t) + d * t ** e)


def coolprop_pv(vp, t):
    tr, pr = vp["reducing_temperature_K"], vp["reducing_value"]
    theta = 1.0 - t / tr
    return pr * math.exp(tr / t * sum(n * theta ** e for n, e in zip(vp["n"], vp["t"])))


def normal_boiling(p, lo, hi, pressure=101325.0):
    """Bisection only to derive the documented isobaric regression interval."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if p(mid) < pressure else (lo, mid)
    return 0.5 * (lo + hi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-directory", required=True)
    out = pathlib.Path(parser.parse_args().output_directory)
    raw = {key: fetch(url) for key, url in SOURCES.items()}
    provenance = [{"id": k, "url": SOURCES[k], "sha256": hashlib.sha256(v).hexdigest()} for k, v in raw.items()]

    cid, index, cas, formula = BUTANOL
    c = compound(ET.fromstring(raw["chemsep_pcd"]), index)
    if c.find("CompoundID").get("value") != "1-butanol":
        raise SystemExit("unexpected ChemSep compound name")
    if c.find("CAS") is not None and c.find("CAS").get("value") != cas:
        raise SystemExit("ChemSep CAS mismatch")
    tc = value(c, "CriticalTemperature", "K")
    pc = value(c, "CriticalPressure", "Pa")
    omega = value(c, "AcentricityFactor", "_")
    mw = value(c, "MolecularWeight", "kg/kmol") / 1000.0
    hf = value(c, "HeatOfFormation", "J/kmol") / 1000.0
    s0 = value(c, "AbsEntropy", "J/kmol/K") / 1000.0
    rpp, rpp_range = correlation(c, "RPPHeatCapacityCp", 100, "ABCDE", "J/kmol/K")
    vp, vp_range = correlation(c, "VaporPressure", 101, "ABCDE", "Pa")
    rho, rho_range = correlation(c, "LiquidDensity", 105, "ABCD", "kmol/m3")
    if not (rpp_range[0] <= 200.0 and rpp_range[1] >= 1000.0):
        raise SystemExit("RPP heat capacity range does not cover 200..1000 K")
    nasa = nasa7_from_rpp([x / 1000.0 for x in rpp], hf, s0)
    component = {
        "id": cid, "cas": cas, "formula": formula, "source_name": "1-butanol",
        "constants": {"source": "chemsep_pcd", "chemsep_library_index": int(index),
                      "critical_temperature_K": tc, "critical_pressure_Pa": pc,
                      "acentric_factor": omega, "molar_mass_kg_per_mol": mw},
        "ideal_gas": {
            "source": "chemsep_pcd", "species": "1-butanol", "correlation": "NASA7",
            # The polynomial is valid to 1000 K only: the upper interval is empty.
            "temperature_breaks_K": [200.0, 1000.0, 1000.0],
            "coefficients": [nasa, nasa],
            "source_correlation": {"name": "ChemSep RPPHeatCapacityCp (Reid-Prausnitz-Poling)",
                                   "eqno": 100, "coefficients_J_per_kmol_K": rpp,
                                   "temperature_range_K": rpp_range},
            "source_note": "a1..a5 = RPP cp polynomial / R (exact); a6, a7 from ChemSep IG heat of formation "
                           f"({hf} J/mol) and IG absolute entropy ({s0} J/(mol K)) at 298.15 K",
            "coefficient_order": "a1,a2,a3,a4,a5,a6,a7",
            "reference_state_status": "engine anchoring h=s=0 at 298.15 K, 1 bar as for the base extract",
        },
    }
    components = {
        "schema": "a2-source-extract-provisional-v1", "version": COMPONENT_VERSION,
        "status": "source transcription for #224 VLLE/LLE implementation verification; no experimental accuracy qualification",
        "sources": [p for p in provenance if p["id"] in ("chemsep_pcd", "artistic_license_2_0")],
        "components": [component],
    }

    pure = [{
        "id": cid, "cas": cas, "source": "chemsep_pcd",
        "vapor_pressure": {"form": "chemsep-101",
                           "description": "ln(p/Pa) = A + B/T + C ln(T/K) + D T^E",
                           "coefficients": vp, "temperature_range_K": vp_range, "units": "Pa"},
        "saturated_liquid_density": {"form": "chemsep-105",
                                     "description": "rho' = A / B^(1 + (1 - T/C)^D)",
                                     "coefficients": rho, "temperature_range_K": rho_range,
                                     "units": "kmol/m^3"},
    }]
    base = json.loads((out / "liquid-nrtl.provisional.json").read_text(encoding="utf-8"))
    water_vp = next(p for p in base["pure_liquid"] if p["id"] == "water")["vapor_pressure"]
    methanol_vp = next(p for p in base["pure_liquid"] if p["id"] == "methanol")["vapor_pressure"]
    boiling = {
        "water": normal_boiling(lambda t: coolprop_pv(water_vp, t), 300.0, 640.0),
        "methanol": normal_boiling(lambda t: coolprop_pv(methanol_vp, t), 250.0, 500.0),
        cid: normal_boiling(lambda t: dippr101(vp, t), 300.0, 550.0),
    }

    table = json.loads(raw["thermo_chemsep_nrtl"])
    meta = table["metadata"]
    if meta["type"] != "NRTL original T" or meta["source"] != "ChemSep" or meta["symmetric"]:
        raise SystemExit("unexpected thermo NRTL metadata")
    lines = raw["dwsim_chemsep_nrtl"].decode("latin-1").splitlines()
    if "Units=cal/mol" not in lines[0] or "DECHEMA NRTL data @ 1atm" not in lines[0]:
        raise SystemExit("unexpected ChemSep IPD header")
    nrtl = []
    for i, j, prefix, cas_i, cas_j in PAIRS:
        rows = [line for line in lines if line.startswith(prefix)]
        reverse = ";".join(prefix.split(";")[1::-1]) + ";"
        if len(rows) != 1 or any(line.startswith(reverse) for line in lines):
            raise SystemExit(f"expected exactly one ChemSep row {prefix}")
        _, _, a12, a21, alpha, comment = rows[0].split(";")
        a12, a21, alpha = float(a12), float(a21), float(alpha)
        b12, b21 = a12 * CAL / R, a21 * CAL / R
        for b, key in ((b12, f"{cas_i} {cas_j}"), (b21, f"{cas_j} {cas_i}")):
            rec = table["data"].get(key)
            if rec is None or abs(b - rec["bij"]) > 1e-12 * abs(b) or rec["alphaij"] != alpha:
                raise SystemExit(f"thermo and ChemSep IPD values disagree for {key}")
        lo, hi = sorted([boiling[i], boiling[j]])
        nrtl.append({
            "i": i, "j": j, "source": "dwsim_chemsep_nrtl; cross-check thermo_chemsep_nrtl",
            "source_row": rows[0], "source_comment": comment,
            "source_units": "A_ij in cal/mol (thermochemical calorie 4.184 J); tau_ij = A_ij/(R T)",
            "A_ij_cal_per_mol": a12, "A_ji_cal_per_mol": a21, "alpha": alpha,
            "tau_form": "tau_ij = t0_ij + t1_ij/T + t2_ij ln(T/K) + t3_ij T; alpha_ij = a0_ij + a1_ij T",
            "tau_ij": {"t0": 0.0, "t1_K": b12, "t2": 0.0, "t3_per_K": 0.0},
            "tau_ji": {"t0": 0.0, "t1_K": b21, "t2": 0.0, "t3_per_K": 0.0},
            "alpha_ij": {"a0": alpha, "a1_per_K": 0.0},
            "regression_pressure_Pa": 101325.0,
            "regression_temperature_K": [math.floor(lo * 100) / 100, math.ceil(hi * 100) / 100],
            "regression_temperature_basis": INFERRED,
        })
    liquids = {
        "schema": "a7-liquid-nrtl-source-extract-v1", "version": LIQUID_VERSION,
        "status": "source transcription for #224 VLLE/LLE implementation verification; no experimental accuracy qualification",
        "gas_constant_J_per_mol_K": R, "calorie_J": CAL,
        "sources": provenance, "pure_liquid": pure, "nrtl_binary": nrtl,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "components-chemsep.provisional.json").write_text(json.dumps(components, indent=2) + "\n", encoding="utf-8")
    (out / "liquid-chemsep-nrtl.provisional.json").write_text(json.dumps(liquids, indent=2) + "\n", encoding="utf-8")
    print("Prepared 1-butanol and", len(nrtl), "NRTL pairs; no expected values generated.")


if __name__ == "__main__":
    main()
