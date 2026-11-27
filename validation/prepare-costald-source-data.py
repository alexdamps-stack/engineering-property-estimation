"""Explicit review-only #228 source transcription. Never run automatically to repair tests.

Downloads the pinned MIT-licensed `chemicals` redistribution of the COSTALD
characteristic volumes and SRK acentric factors (Hankinson and Thomson 1979,
Table 1) for the embedded component set and writes a transcription with the
source URL, SHA-256 and licence notice. Evaluates no thermodynamic model.
Standard library only.

    python validation/prepare-costald-source-data.py --output-directory data/source-extracts
"""
import argparse
import csv
import hashlib
import io
import json
import pathlib
import urllib.request

CHEMICALS = "7fca9acf7708e2a646bc04d2b6f6996d0325e77d"  # tag v1.5.2
SOURCES = {
    "chemicals_costald": f"https://raw.githubusercontent.com/CalebBell/chemicals/{CHEMICALS}/chemicals/Density/COSTALD%20Parameters.tsv",
    "chemicals_license": f"https://raw.githubusercontent.com/CalebBell/chemicals/{CHEMICALS}/LICENSE.txt",
}
# Embedded component database identifiers and CAS numbers (base and reacting extracts).
COMPONENTS = [
    ("methane", "74-82-8"),
    ("ethane", "74-84-0"),
    ("propane", "74-98-6"),
    ("nitrogen", "7727-37-9"),
    ("carbon-dioxide", "124-38-9"),
    ("oxygen", "7782-44-7"),
    ("water", "7732-18-5"),
    ("methanol", "67-56-1"),
    ("hydrogen", "1333-74-0"),
    ("carbon-monoxide", "630-08-0"),
]
VERSION = "2026-11-24.1"


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-directory", required=True)
    out = pathlib.Path(parser.parse_args().output_directory)
    raw = {key: fetch(url) for key, url in SOURCES.items()}
    provenance = [
        {"id": k, "url": SOURCES[k], "sha256": hashlib.sha256(v).hexdigest()}
        for k, v in raw.items()
    ]
    rows = list(csv.DictReader(io.StringIO(raw["chemicals_costald"].decode("utf-8")), delimiter="\t"))
    if not rows or set(rows[0]) != {"CAS", "Chemical", "omega_SRK", "Vchar", "Z_RA"}:
        raise SystemExit("unexpected COSTALD table header")
    records = []
    for cid, cas in COMPONENTS:
        matches = [r for r in rows if r["CAS"] == cas]
        if len(matches) != 1:
            raise SystemExit(f"expected exactly one COSTALD row for {cid} ({cas})")
        row = matches[0]
        # Keep the decimal text so the transcription is exact.
        records.append(
            {
                "id": cid,
                "cas": cas,
                "source_name": row["Chemical"],
                "characteristic_volume_m3_per_mol": float(row["Vchar"]),
                "characteristic_volume_text": row["Vchar"],
                "omega_srk": float(row["omega_SRK"]),
                "omega_srk_text": row["omega_SRK"],
            }
        )
    extract = {
        "schema": "issue-228-costald-source-extract-v1",
        "version": VERSION,
        "status": "source transcription; experimental comparison in validation/fixtures/liquid-density-v1.json",
        "gas_constant_J_per_mol_K": 8.31446261815324,
        "method": {
            "saturated": "Hankinson, R. W.; Thomson, G. H. A new correlation for saturated densities of liquids and their mixtures. AIChE J. 1979, 25 (4), 653-663. doi:10.1002/aic.690250412",
            "compressed": "Thomson, G. H.; Brobst, K. R.; Hankinson, R. W. An improved correlation for densities of compressed liquids and liquid mixtures. AIChE J. 1982, 28 (4), 671-676. doi:10.1002/aic.690280420",
            "parameters": "V* and omega_SRK from Hankinson and Thomson (1979) Table 1 as redistributed by chemicals 1.5.2 (MIT, Caleb Bell); critical temperatures are taken from the embedded component database, not from this extract",
        },
        "sources": provenance,
        "components": records,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "liquid-density-costald.provisional.json").write_text(
        json.dumps(extract, indent=2) + "\n", encoding="utf-8"
    )
    (out / "chemicals-LICENSE.txt").write_bytes(raw["chemicals_license"])


if __name__ == "__main__":
    main()
