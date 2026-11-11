"""Compare a freshly generated reference candidate with its committed fixture.

Never writes or replaces fixtures. Structure, strings, input hashes and other
non-numeric values must match exactly. Numbers are compared as
|candidate - committed| <= atol + rtol*|committed| with:

* the fixture's own declared "tolerances" entry for the quantity, when the
  number's key (or an alias in DECLARED below) names one;
* a 1e-8 absolute floor for solver diagnostics that are near zero by
  construction (fugacity residuals, TPD minima, reduced Gibbs changes), the
  magnitude of the A4 fixtures' declared stable_tm_floor;
* otherwise a reproduction tolerance (default rtol 1e-12, atol 1e-15), which
  absorbs last-bit floating-point differences between machines and library
  builds (observed up to ~7e-18 absolute in cubic-v1.json) while staying three
  orders of magnitude tighter than the tightest declared budget (2e-9).

The recorded interpreter/package versions ("python", "packages") describe the
environment the committed fixture was generated in; differences there are
reported but do not by themselves fail the comparison.

usage: python validation/compare-reference-candidate.py CANDIDATE COMMITTED [--rtol R] [--atol A]
"""

import argparse
import json
from pathlib import Path
import sys

ENVIRONMENT_KEYS = {"python", "packages"}
# Number key -> candidate names in the fixture's "tolerances" block; the first
# one the fixture declares is used.
DECLARED = {
    "T": ["temperature_K", "temperature"],
    "target_h": ["enthalpy_J_per_mol"], "vapor_h": ["enthalpy_J_per_mol"], "liquid_h": ["enthalpy_J_per_mol"],
    "vapor_fraction": ["vapor_fraction", "phase_fraction"],
    "x": ["composition"], "z": ["composition"], "vapor": ["composition"], "liquid": ["composition"],
    "Z_vapor": ["Z"], "Z_liquid": ["Z"],
}
DIAGNOSTIC_FLOOR = 1e-8
DIAGNOSTICS = {"fugacity_residual", "feed_min_tm", "reduced_gibbs_change"}


def tolerance(key, declared, args):
    if key in DIAGNOSTICS:
        return DIAGNOSTIC_FLOOR, 0.0
    for name in [key] + DECLARED.get(key, []):
        value = declared.get(name)
        if isinstance(value, list) and len(value) == 2:
            return float(value[0]), float(value[1])
    return args.atol, args.rtol


def compare(candidate, committed, path, key, declared, args, failures, worst):
    if isinstance(committed, dict):
        if not isinstance(candidate, dict):
            failures.append(f"{path}: type differs")
            return
        if set(candidate) != set(committed):
            failures.append(f"{path}: keys differ ({sorted(set(candidate) ^ set(committed))})")
            return
        for k in committed:
            compare(candidate[k], committed[k], f"{path}.{k}", k, declared, args, failures, worst)
    elif isinstance(committed, list):
        if not isinstance(candidate, list) or len(candidate) != len(committed):
            failures.append(f"{path}: length differs")
            return
        for i, (a, b) in enumerate(zip(candidate, committed)):
            compare(a, b, f"{path}[{i}]", key, declared, args, failures, worst)
    elif isinstance(committed, (int, float)) and not isinstance(committed, bool):
        if isinstance(candidate, bool) or not isinstance(candidate, (int, float)):
            failures.append(f"{path}: type differs")
            return
        deviation = abs(candidate - committed)
        if deviation > worst[0]:
            worst[:] = [deviation, path]
        atol, rtol = tolerance(key, declared, args)
        if deviation > atol + rtol * abs(committed):
            failures.append(f"{path}: candidate {candidate!r} vs committed {committed!r} (atol {atol:g}, rtol {rtol:g})")
    elif candidate != committed:
        failures.append(f"{path}: candidate {candidate!r} vs committed {committed!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("committed", type=Path)
    parser.add_argument("--rtol", type=float, default=1e-12, help="reproduction tolerance for undeclared quantities")
    parser.add_argument("--atol", type=float, default=1e-15, help="reproduction tolerance for undeclared quantities")
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    committed = json.loads(args.committed.read_text(encoding="utf-8"))
    for key in sorted(ENVIRONMENT_KEYS & (set(candidate) | set(committed))):
        if candidate.get(key) != committed.get(key):
            print(f"note: {args.committed} recorded {key}={committed.get(key)!r}; "
                  f"this run used {candidate.get(key)!r}")
    declared = committed.get("tolerances", {}) if isinstance(committed.get("tolerances"), dict) else {}
    failures, worst = [], [0.0, None]
    compare({k: v for k, v in candidate.items() if k not in ENVIRONMENT_KEYS},
            {k: v for k, v in committed.items() if k not in ENVIRONMENT_KEYS},
            "$", "", declared, args, failures, worst)
    print(f"{args.committed}: maximum absolute numeric deviation {worst[0]:.3g} at {worst[1]}")
    if failures:
        for failure in failures[:50]:
            print("  " + failure, file=sys.stderr)
        sys.exit(f"{args.committed}: {len(failures)} difference(s) from the fresh candidate; "
                 "review explicitly and never replace the fixture automatically")
    print(f"{args.committed}: reproduced within declared tolerances (undeclared: rtol {args.rtol:g}, atol {args.atol:g})")


if __name__ == "__main__":
    main()
