"""Independent GERG fugacity references from teqp 0.23.2 AD; no Rust calls.

Uses independently solved TP reference densities, adds retained-zero and trace
components. Defaults to read-only fixture verification; --output writes a candidate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import teqp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"validation/fixtures/gerg-density-v1.json"
ALIASES = {"carbon-dioxide":"carbondioxide", "carbon-monoxide":"carbonmonoxide", "hydrogen-sulfide":"hydrogensulfide"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    assert teqp.__version__ == "0.23.2"
    cases = []
    inputs = json.loads(SOURCE.read_text())["cases"]
    for case in inputs:
        cases.append(dict(id=case["id"],composition=case["composition"],temperature_k=case["temperature_k"],
                          molar_density=case["teqp_molar_density"]))
    for tiny in [0.,1e-12,1e-6]:
        cases.append(dict(id=f"infinite-dilution-{tiny}",temperature_k=300.,molar_density=1000.,
            composition={"methane":.9-tiny,"ethane":.1,"water":tiny,"hydrogen-sulfide":0.,"n-decane":0.}))
    for c in cases:
        c["component_ids"] = list(c["composition"])
        x = np.array(list(c["composition"].values()))
        model = teqp.make_model(dict(kind="GERG2008resid",model=dict(names=[ALIASES.get(n,n) for n in c["composition"]])))
        c["ln_fugacity_coefficients"] = np.log(model.get_fugacity_coefficients(c["temperature_k"],c["molar_density"]*x)).tolist()
        assert all(np.isfinite(c["ln_fugacity_coefficients"])), c
    output = dict(reference="teqp 0.23.2 GERG2008 residual Helmholtz AD chemical potentials",
        versions=dict(teqp=teqp.__version__,numpy=np.__version__),
        input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        tolerance=[2e-8,2e-9],cases=cases)
    if args.output:
        args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
        print(f"Wrote {len(cases)} GERG fugacity reference states to {args.output}")
    else:
        assert json.loads((ROOT/"validation/fixtures/gerg-fugacity-v1.json").read_text()) == output, "GERG fugacity reproduction differs; review a candidate explicitly"
        print(f"Verified {len(cases)} GERG fugacity reference states")


if __name__ == "__main__":
    main()
