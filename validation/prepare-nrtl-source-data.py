"""Explicit review-only A7 source transcription. Never run automatically to repair tests.

Downloads pinned CoolProp fluid files (vapor-pressure and saturated-liquid density
ancillaries), the MIT-licensed `thermo` redistribution of the ChemSep/DECHEMA NRTL
table and the DWSIM copy of the original ChemSep IPD file (cross-check of the
cal/mol source values). Writes a transcription; evaluates no thermodynamic model.
Standard library only.
"""
import argparse
import hashlib
import json
import math
import pathlib
import urllib.request

COOLPROP = "98b3523d5daa98454618d381d2ae53f7471d216b"
THERMO = "2bb466e98439c2395a004e7095d07f9b97a5a0f0"
DWSIM = "1abf72d1b6b41d3e9a8cc770d3cc4e8fc76e5766"
SOURCES = {
    "coolprop_water": f"https://raw.githubusercontent.com/CoolProp/CoolProp/{COOLPROP}/dev/fluids/Water.json",
    "coolprop_methanol": f"https://raw.githubusercontent.com/CoolProp/CoolProp/{COOLPROP}/dev/fluids/Methanol.json",
    "thermo_chemsep_nrtl": f"https://raw.githubusercontent.com/CalebBell/thermo/{THERMO}/thermo/Interaction%20Parameters/ChemSep/nrtl.json",
    "thermo_license": f"https://raw.githubusercontent.com/CalebBell/thermo/{THERMO}/LICENSE.txt",
    "dwsim_chemsep_nrtl": f"https://raw.githubusercontent.com/DanWBR/dwsim/{DWSIM}/DWSIM.Thermodynamics/Assets/nrtl.dat",
}
# Thermochemical calorie (exact) and the A1 gas constant.
CAL = 4.184
R = 8.31446261815324
PURE = [("water", "7732-18-5", "coolprop_water"), ("methanol", "67-56-1", "coolprop_methanol")]


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read()


def ancillary(record, kind):
    keys = {"T_r", "Tmax", "Tmin", "n", "t", "reducing_value", "type", "using_tau_r", "max_abserror_percentage", "description"}
    if not keys <= record.keys() or record["type"] != kind or len(record["n"]) != len(record["t"]):
        raise SystemExit(f"unexpected {kind} ancillary shape")
    return {
        "form": kind,
        "description": record["description"],
        "reducing_temperature_K": record["T_r"],
        "reducing_value": record["reducing_value"],
        "using_tau_r": record["using_tau_r"],
        "n": record["n"],
        "t": record["t"],
        "temperature_range_K": [record["Tmin"], record["Tmax"]],
        "source_max_abs_error_percent_vs_reference_eos": record["max_abserror_percentage"],
    }


def normal_boiling(vp, pressure=101325.0):
    """Bisection on the transcribed pV ancillary, only to derive the documented
    isobaric regression interval of the 1 atm NRTL data set."""
    tr, pr = vp["reducing_temperature_K"], vp["reducing_value"]

    def p(t):
        theta = 1.0 - t / tr
        return pr * math.exp(tr / t * sum(n * theta**e for n, e in zip(vp["n"], vp["t"])))

    lo, hi = vp["temperature_range_K"][0], vp["temperature_range_K"][1] - 1e-6
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

    pure = []
    for cid, cas, key in PURE:
        fluid = json.loads(raw[key])
        if fluid["INFO"]["CAS"] != cas:
            raise SystemExit(f"CAS mismatch for {cid}")
        anc = fluid["ANCILLARIES"]
        vp = ancillary(anc["pS"], "pV")
        rho = ancillary(anc["rhoL"], "rhoLnoexp")
        if not vp["using_tau_r"] or rho["using_tau_r"]:
            raise SystemExit("unexpected ancillary reduction flags")
        vp["units"] = {"reducing_value": "Pa"}
        rho["units"] = {"reducing_value": "mol/m^3"}
        pure.append({"id": cid, "cas": cas, "source": key, "reference_eos": fluid["EOS"][0]["BibTeX_EOS"],
                     "vapor_pressure": vp, "saturated_liquid_density": rho})

    table = json.loads(raw["thermo_chemsep_nrtl"])
    meta = table["metadata"]
    if meta["type"] != "NRTL original T" or meta["source"] != "ChemSep" or meta["symmetric"]:
        raise SystemExit("unexpected thermo NRTL metadata")
    forward = table["data"]["67-56-1 7732-18-5"]
    reverse = table["data"]["7732-18-5 67-56-1"]
    lines = raw["dwsim_chemsep_nrtl"].decode("latin-1").splitlines()
    if "Units=cal/mol" not in lines[0] or "DECHEMA NRTL data @ 1atm" not in lines[0]:
        raise SystemExit("unexpected ChemSep IPD header")
    row = [line for line in lines if line.startswith("1101;1921;")]
    if len(row) != 1 or any(line.startswith("1921;1101;") for line in lines):
        raise SystemExit("expected exactly one methanol/water ChemSep row")
    _, _, a12, a21, alpha, comment = row[0].split(";")
    a12, a21, alpha = float(a12), float(a21), float(alpha)
    b12, b21 = a12 * CAL / R, a21 * CAL / R
    # thermo stores A/R_cal; require agreement with the original cal/mol row.
    for b, rec in ((b12, forward), (b21, reverse)):
        if abs(b - rec["bij"]) > 1e-12 * abs(b) or rec["alphaij"] != alpha:
            raise SystemExit("thermo and ChemSep IPD values disagree")
    boiling = {p["id"]: normal_boiling(p["vapor_pressure"]) for p in pure}
    regression = [math.floor(boiling["methanol"] * 100) / 100, math.ceil(boiling["water"] * 100) / 100]
    nrtl = [{
        "i": "methanol", "j": "water", "source": "dwsim_chemsep_nrtl; cross-check thermo_chemsep_nrtl",
        "source_row": row[0], "source_comment": comment,
        "source_units": "A_ij in cal/mol (thermochemical calorie 4.184 J); tau_ij = A_ij/(R T)",
        "A_ij_cal_per_mol": a12, "A_ji_cal_per_mol": a21, "alpha": alpha,
        "tau_form": "tau_ij = t0_ij + t1_ij/T + t2_ij ln(T/K) + t3_ij T; alpha_ij = a0_ij + a1_ij T",
        "tau_ij": {"t0": 0.0, "t1_K": b12, "t2": 0.0, "t3_per_K": 0.0},
        "tau_ji": {"t0": 0.0, "t1_K": b21, "t2": 0.0, "t3_per_K": 0.0},
        "alpha_ij": {"a0": alpha, "a1_per_K": 0.0},
        "regression_pressure_Pa": 101325.0,
        "regression_temperature_K": regression,
        "regression_temperature_basis": "isobaric 1 atm data: interval between pure normal boiling points from the transcribed vapor-pressure ancillaries, rounded outward to 0.01 K",
    }]
    document = {
        "schema": "a7-liquid-nrtl-source-extract-v1", "version": "2026-11-08.1",
        "status": "source transcription for implementation verification; no experimental accuracy qualification",
        "gas_constant_J_per_mol_K": R, "calorie_J": CAL,
        "sources": provenance, "pure_liquid": pure, "nrtl_binary": nrtl,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "liquid-nrtl.provisional.json").write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    (out / "thermo-LICENSE.txt").write_bytes(raw["thermo_license"])
    print(f"Prepared {len(pure)} pure-liquid and {len(nrtl)} NRTL transcriptions; no expected values generated.")


if __name__ == "__main__":
    main()
