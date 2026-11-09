"""Compare a freshly generated reference candidate with its committed fixture.

Never writes or replaces fixtures. Structure, strings, input hashes and integer
or exact inputs must match. Numbers must agree within a reproduction tolerance
(default: exact), applied as |candidate - committed| <= atol + rtol*|committed|.
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


def compare(candidate, committed, path, args, failures, worst):
    if isinstance(committed, dict):
        if not isinstance(candidate, dict):
            failures.append(f"{path}: type differs")
            return
        if set(candidate) != set(committed):
            failures.append(f"{path}: keys differ ({sorted(set(candidate) ^ set(committed))})")
            return
        for key in committed:
            compare(candidate[key], committed[key], f"{path}.{key}", args, failures, worst)
    elif isinstance(committed, list):
        if not isinstance(candidate, list) or len(candidate) != len(committed):
            failures.append(f"{path}: length differs")
            return
        for i, (a, b) in enumerate(zip(candidate, committed)):
            compare(a, b, f"{path}[{i}]", args, failures, worst)
    elif isinstance(committed, (int, float)) and not isinstance(committed, bool):
        if isinstance(candidate, bool) or not isinstance(candidate, (int, float)):
            failures.append(f"{path}: type differs")
            return
        deviation = abs(candidate - committed)
        if deviation > worst[0]:
            worst[:] = [deviation, path]
        if deviation > args.atol + args.rtol * abs(committed):
            failures.append(f"{path}: candidate {candidate!r} vs committed {committed!r}")
    elif candidate != committed:
        failures.append(f"{path}: candidate {candidate!r} vs committed {committed!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("committed", type=Path)
    parser.add_argument("--rtol", type=float, default=0.0)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args()
    candidate = json.loads(args.candidate.read_text(encoding="utf-8"))
    committed = json.loads(args.committed.read_text(encoding="utf-8"))
    for key in sorted(ENVIRONMENT_KEYS & (set(candidate) | set(committed))):
        if candidate.get(key) != committed.get(key):
            print(f"note: {args.committed} recorded {key}={committed.get(key)!r}; "
                  f"this run used {candidate.get(key)!r}")
    failures, worst = [], [0.0, None]
    compare({k: v for k, v in candidate.items() if k not in ENVIRONMENT_KEYS},
            {k: v for k, v in committed.items() if k not in ENVIRONMENT_KEYS}, "$", args, failures, worst)
    print(f"{args.committed}: maximum absolute numeric deviation {worst[0]:.3g} at {worst[1]}")
    if failures:
        for failure in failures[:50]:
            print("  " + failure, file=sys.stderr)
        sys.exit(f"{args.committed}: {len(failures)} difference(s) from the fresh candidate; "
                 "review explicitly and never replace the fixture automatically")
    print(f"{args.committed}: reproduced (rtol {args.rtol:g}, atol {args.atol:g})")


if __name__ == "__main__":
    main()
