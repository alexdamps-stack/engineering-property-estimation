"""Explicit review-only #225 source transcription. Never run automatically to repair tests.

Downloads the ChemSep v8.3 pure-component database (Artistic License 2.0) from a
pinned DWSIM commit and the SPDX text of that licence. Transcribes, for the eight
qualified transport components, the critical volume and dipole moment used by
the Chung et al. (1988) gas method and the ChemSep temperature correlations of
saturated-liquid viscosity, liquid thermal conductivity and surface tension.
Adds the Chung association factors of water and methanol from the cited paper.
Writes a transcription; evaluates no transport model. Standard library only.
For #228 it also transcribes the ChemSep saturated-liquid density correlations
of water and methanol into a separate liquid-density-chemsep.provisional.json.

    python validation/prepare-transport-source-data.py --output-directory data/source-extracts
"""
import argparse
import hashlib
import json
import pathlib
import urllib.request
import xml.etree.ElementTree as ET

DWSIM = "1abf72d1b6b41d3e9a8cc770d3cc4e8fc76e5766"
SPDX = "v3.24.0"
SOURCES = {
    "chemsep_pcd": f"https://raw.githubusercontent.com/DanWBR/dwsim/{DWSIM}/DWSIM.Thermodynamics/Assets/Databases/chemsep1.xml",
    "artistic_license_2_0": f"https://raw.githubusercontent.com/spdx/license-list-data/{SPDX}/text/Artistic-2.0.txt",
}
VERSION = "2026-11-24.1"
SCHEMA = "transport-source-extract-v1"
# 1 debye = 1e-21 / c  C m (exact definition through the speed of light).
DEBYE_C_M = 1e-21 / 299792458.0
# Component ID, CAS, Chung association factor kappa. kappa is zero for
# non-associating species; the water and methanol values are those of
# Chung, Ajlan, Lee and Starling (1988), Ind. Eng. Chem. Res. 27, 671-679,
# Table I, as reproduced in Poling, Prausnitz and O'Connell (2001),
# The Properties of Gases and Liquids, 5th ed., Table 9-1.
COMPONENTS = [
    ("methane", "74-82-8", 0.0),
    ("ethane", "74-84-0", 0.0),
    ("propane", "74-98-6", 0.0),
    ("nitrogen", "7727-37-9", 0.0),
    ("carbon-dioxide", "124-38-9", 0.0),
    ("oxygen", "7782-44-7", 0.0),
    ("water", "7732-18-5", 0.075908),
    ("methanol", "67-56-1", 0.215175),
]
# Correlation records and the ChemSep equation numbers this transcription accepts:
# 101: Y = exp(A + B/T + C ln T + D T^E); 16: Y = A + exp(B/T + C + D T + E T^2).
CORRELATIONS = {
    "liquid_viscosity": ("LiquidViscosity", 101, "Pa.s", "Pa s"),
    "liquid_thermal_conductivity": ("LiquidThermalConductivity", 16, "W/m/K", "W/(m K)"),
    "surface_tension": ("SurfaceTension", 16, "N/m", "N/m"),
}


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def value(compound, tag, units):
    element = compound.find(tag)
    if element is None or element.get("units") != units:
        raise SystemExit(f"missing or unexpected {tag}")
    return float(element.get("value"))


def correlation(compound, tag, equation, units, label):
    element = compound.find(tag)
    if element is None or element.get("units") != units:
        raise SystemExit(f"missing or unexpected {tag}")
    if int(element.find("eqno").get("value")) != equation:
        raise SystemExit(f"{tag}: expected ChemSep equation {equation}")
    coefficients = [float(element.find(k).get("value")) for k in "ABCDE"]
    t_min = float(element.find("Tmin").get("value"))
    t_max = float(element.find("Tmax").get("value"))
    if element.find("Tmin").get("units") != "K" or element.find("Tmax").get("units") != "K" or not 0 < t_min < t_max:
        raise SystemExit(f"{tag}: invalid temperature range")
    return {"equation": equation, "coefficients": coefficients, "units": label, "temperature_range_K": [t_min, t_max]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-directory", required=True)
    out = pathlib.Path(parser.parse_args().output_directory)
    raw = {key: fetch(url) for key, url in SOURCES.items()}
    provenance = [{"id": k, "url": SOURCES[k], "sha256": hashlib.sha256(v).hexdigest()} for k, v in raw.items()]
    root = ET.fromstring(raw["chemsep_pcd"].decode("latin-1"))
    header = root.get("library", "")
    if "ChemSep v8.3" not in header or "artistic_license_2_0" not in header:
        raise SystemExit(f"unexpected ChemSep library header: {header!r}")
    by_cas = {}
    for compound in root.findall("compound"):
        cas = compound.find("CAS")
        if cas is not None:
            by_cas.setdefault(cas.get("value"), []).append(compound)

    records = []
    for cid, cas, kappa in COMPONENTS:
        matches = by_cas.get(cas, [])
        if len(matches) != 1:
            raise SystemExit(f"expected one ChemSep record for {cid} ({cas}), found {len(matches)}")
        c = matches[0]
        dipole = value(c, "DipoleMoment", "Coulomb.m")
        record = {
            "id": cid,
            "cas": cas,
            "chemsep_name": c.find("CompoundID").get("value"),
            "chemsep_index": int(c.find("LibraryIndex").get("value")),
            "critical_volume_m3_per_mol": value(c, "CriticalVolume", "m3/kmol") / 1000.0,
            "dipole_moment_C_m": dipole,
            "dipole_moment_debye": dipole / DEBYE_C_M,
            "association_factor": kappa,
        }
        for key, (tag, equation, units, label) in CORRELATIONS.items():
            record[key] = correlation(c, tag, equation, units, label)
        records.append(record)

    extract = {
        "schema": SCHEMA,
        "version": VERSION,
        "status": "source transcription for #225 transport properties; accuracy is qualified separately against CoolProp 7.2.0 in validation/transport",
        "sources": provenance,
        "chemsep_library": header,
        "equations": {
            "101": "Y = exp(A + B/T + C ln(T/K) + D T^E), T in K",
            "16": "Y = A + exp(B/T + C + D T + E T^2), T in K",
        },
        "association_factor_source": "Chung, Ajlan, Lee, Starling (1988), Ind. Eng. Chem. Res. 27, 671-679, Table I; reproduced in Poling, Prausnitz, O'Connell (2001), Table 9-1. Zero for non-associating species.",
        "unit_conversions": {"critical_volume": "ChemSep m3/kmol divided by 1000", "debye_C_m": DEBYE_C_M},
        "components": records,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "transport.provisional.json").write_text(json.dumps(extract, indent=2) + "\n", encoding="utf-8")
    (out / "liquid-density-chemsep.provisional.json").write_text(
        json.dumps(liquid_density_extract(root, header, provenance), indent=2) + "\n", encoding="utf-8"
    )
    (out / "ChemSep-LICENSE.txt").write_bytes(
        b"ChemSep v8.3 pure component data - Copyright (c) Harry Kooijman and Ross Taylor (2021).\n"
        + f"Distributed under the Artistic License 2.0; obtained from {SOURCES['chemsep_pcd']}\n".encode()
        + b"Only the transcribed values in transport.provisional.json and liquid-density-chemsep.provisional.json are redistributed here; they are unmodified apart from SI unit conversion.\n\n"
        + raw["artistic_license_2_0"]
    )



# #228: component-specific saturated-liquid densities for the components that
# COSTALD does not qualify. Hydrogen is excluded: its ChemSep range (13.8-33.2 K)
# lies below every package temperature interval.
LIQUID_DENSITY_VERSION = "2026-11-27.1"
LIQUID_DENSITY_COMPONENTS = [("water", "7732-18-5", 106), ("methanol", "67-56-1", 105)]


def liquid_density_extract(root, header, provenance):
    records = []
    for cid, cas, equation in LIQUID_DENSITY_COMPONENTS:
        matches = [c for c in root.findall("compound") if c.find("CAS") is not None and c.find("CAS").get("value") == cas]
        if len(matches) != 1:
            raise SystemExit(f"expected one ChemSep record for {cid} ({cas})")
        c = matches[0]
        element = c.find("LiquidDensity")
        if element is None or element.get("units") != "kmol/m3" or int(element.find("eqno").get("value")) != equation:
            raise SystemExit(f"{cid}: unexpected LiquidDensity record")
        coefficients = [float(element.find(k).get("value")) if element.find(k) is not None else 0.0 for k in "ABCDE"]
        t_min = float(element.find("Tmin").get("value"))
        t_max = float(element.find("Tmax").get("value"))
        if element.find("Tmin").get("units") != "K" or element.find("Tmax").get("units") != "K" or not 0 < t_min < t_max:
            raise SystemExit(f"{cid}: invalid temperature range")
        records.append({
            "id": cid,
            "cas": cas,
            "chemsep_name": c.find("CompoundID").get("value"),
            "chemsep_index": int(c.find("LibraryIndex").get("value")),
            "equation": equation,
            # SI: kmol/m3 is numerically mol/L; times 1000 gives mol/m3.
            "coefficients": coefficients,
            "coefficient_units": "kmol/m3 (A); multiply the result by 1000 for mol/m3",
            "temperature_range_K": [t_min, t_max],
            "chemsep_critical_temperature_K": value(c, "CriticalTemperature", "K"),
        })
    return {
        "schema": "liquid-density-chemsep-source-extract-v1",
        "version": LIQUID_DENSITY_VERSION,
        "status": "source transcription for the #228 engineering liquid density; accuracy is qualified against CoolProp 7.2.0 in validation/fixtures/liquid-density-v1.json",
        "sources": provenance,
        "chemsep_library": header,
        "equations": {
            "105": "rho = A / B^(1 + (1 - T/C)^D), kmol/m3, T in K",
            "106": "rho = A (1 - Tr)^(B + C Tr + D Tr^2 + E Tr^3), Tr = T/Tc with the ChemSep critical temperature, kmol/m3",
        },
        "components": records,
    }

if __name__ == "__main__":
    main()
