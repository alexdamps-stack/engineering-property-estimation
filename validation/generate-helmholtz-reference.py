"""Explicit dilute Helmholtz reference generation using only Python's stdlib.

80-digit Decimal evaluation of the published PR residual Helmholtz expression;
no Rust calls or thermodynamic software/runtime dependencies. Never run by tests.
"""
import argparse
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = ROOT / "data/source-extracts/components.provisional.json"
    extract = json.loads(source.read_text(), parse_float=D, parse_int=D)
    constants = next(c for c in extract["components"] if c["id"] == "methane")["constants"]
    cases = []
    with localcontext() as ctx:
        ctx.prec = 80
        r = D("8.31446261815324")
        tc = constants["critical_temperature_K"]
        pc = constants["critical_pressure_Pa"]
        w = constants["acentric_factor"]
        t = D(300)
        tr = t / tc
        m = D("0.37464") + D("1.54226")*w - D("0.26992")*w*w
        b = D("0.07779607390388846")*r*tc/pc
        sqrt2 = D(2).sqrt()
        for variant in ["PR", "PR-BM"]:
            if variant == "PR":
                alpha = (1 + m*(1-tr.sqrt()))**2
            else:
                d = 1 + m/2
                alpha = (2*(1-1/d)*(1-(d*tr.ln()).exp())).exp()
            a = D("0.4572355289213822")*r*r*tc*tc/pc*alpha
            for volume, moles in [(D(1), D(1000)), (D("1e27"), D("1e15")), (D("1e30"), D("1e15"))]:
                v = volume/moles
                q = ((v+(1+sqrt2)*b)/(v+(1-sqrt2)*b)).ln()/(2*sqrt2*b)
                energy = moles*(-r*t*(1-b/v).ln()-a*q)
                cases.append({"variant": variant, "component": "methane", "temperature_K": str(t),
                              "volume_m3": str(volume), "moles": str(moles), "energy_J": str(energy)})
    output = {"schema": "helmholtz-decimal-reference-v1", "precision_decimal_digits": 80,
              "generator": "validation/generate-helmholtz-reference.py",
              "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "method": "Independent Decimal evaluation of untranslated PR Helmholtz expression, same pinned constants; no Rust calls",
              "cases": cases}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8")
    print(f"Wrote {len(cases)} 80-digit Helmholtz reference states")


if __name__ == "__main__":
    main()
