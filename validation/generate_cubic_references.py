"""Generate fixtures with teqp, not the Rust engine. Run from the repository root.

Requires requirements-reference.txt. Output changes must be explicitly reviewed.
Uses teqp Helmholtz derivatives and numerical pressure bracketing; does not port
the Rust fugacity/departure expressions or polynomial root solver.
"""

import argparse
import hashlib
import json
import platform
from pathlib import Path

import CoolProp
import CoolProp.CoolProp as CP
import numpy as np
import scipy
from scipy.optimize import brentq
import teqp

ROOT = Path(__file__).resolve().parents[1]
R = 8.31446261815324
FLUIDS = {
    "methane": "Methane", "ethane": "Ethane", "propane": "Propane",
    "nitrogen": "Nitrogen", "carbon-dioxide": "CarbonDioxide", "oxygen": "Oxygen",
    "water": "Water", "methanol": "Methanol",
}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def generate_data():
    records = []
    for identifier, fluid in FLUIDS.items():
        raw = CP.get_fluid_param_string(fluid, "JSON")
        fluid_data = json.loads(raw)[0]
        eos = fluid_data["EOS"][0]
        reducing = eos["STATES"]["reducing"]
        records.append({
            "id": identifier, "cas": fluid_data["INFO"]["CAS"],
            "critical_temperature_k": reducing["T"],
            "critical_pressure_pa": reducing["p"],
            "molar_mass_kg_mol": eos["molar_mass"],
            "acentric_factor": eos["acentric"],
            "source_fluid": fluid,
            "source_json_sha256": hashlib.sha256(raw.encode()).hexdigest(),
        })
    data = {
        "version": "coolprop-8.0.0-initial-v1",
        "source": "CoolProp 8.0.0 bundled fluid JSON, first EOS reducing T/p, molar mass and acentric factor",
        "source_url": "https://github.com/CoolProp/CoolProp/tree/v8.0.0/dev/fluids",
        "usage_terms": "CoolProp MIT license; see data/NOTICE.md",
        "limitations": "Selected constants for implementation checks, not experimentally qualified mixtures. No interaction parameters or ideal caloric model supplied.",
        "components": records,
    }
    write_json(ROOT / "data/components.json", data)
    return data


def request(model, ids, fractions, temperature, pressure, root="largest", kij=0.0):
    n = len(ids)
    matrix = [[0.0 if i == j else kij for j in range(n)] for i in range(n)]
    return dict(model=model, components=ids, fractions=fractions, basis="mole",
                temperature_k=temperature, pressure_pa=pressure,
                binary_interactions=matrix, root=root)


def reference(req, records):
    cs = [records[c] for c in req["components"]]
    temperature, pressure = req["temperature_k"], req["pressure_pa"]
    x = np.array(req["fractions"], dtype=float)
    if req["basis"] == "mass":
        x /= np.array([c["molar_mass_kg_mol"] for c in cs])
    x /= x.sum()
    alphas = []
    for c in cs:
        w = c["acentric_factor"]
        k = 0.37464 + 1.54226*w - 0.26992*w*w
        if req["model"] == "PR-BM" and temperature > c["critical_temperature_k"]:
            d = 1 + k/2
            # Twu L=2*(1-1/d), M=1, N=d is exactly the BM continuation.
            alphas.append({"type": "Twu", "c": [2*(1-1/d), 1.0, d]})
        else:
            alphas.append({"type": "Mathias-Copeman", "c": [k, 0.0, 0.0]})
    model_spec = {"kind": "cubic", "model": {
        "type": "PR", "Tcrit / K": [c["critical_temperature_k"] for c in cs],
        "pcrit / Pa": [c["critical_pressure_pa"] for c in cs],
        "acentric": [c["acentric_factor"] for c in cs],
        "alpha": alphas, "kmat": req["binary_interactions"],
    }}
    model = teqp.make_model(model_spec)
    # A pressure function evaluated through teqp's independent Helmholtz AD.
    def residual(rho):
        return rho * model.get_R(x) * temperature * (1 + model.get_Ar01(temperature, rho, x)) - pressure
    # Search molar density directly, bounded by the EOS excluded volume.
    grid = np.geomspace(1e-8, (1 - 1e-10) / model.get_b(temperature, x), 12000)
    densities = []
    previous = grid[0]
    for rho in grid[1:]:
        if residual(previous) * residual(rho) < 0:
            densities.append(brentq(residual, previous, rho, xtol=1e-10, rtol=1e-14))
        previous = rho
    if not densities:
        raise RuntimeError("No independently bracketed density root")
    zs = sorted(pressure / (rho * R * temperature) for rho in densities)
    z = zs[0] if req["root"] == "smallest" else zs[-1]
    rho = pressure / (z * R * temperature)
    ar00 = model.get_Ar00(temperature, rho, x)
    ar10 = model.get_Ar10(temperature, rho, x)
    ar01 = model.get_Ar01(temperature, rho, x)
    return {"request": req, "teqp_model": model_spec, "expected": {
        "roots_z": zs, "z": z, "molar_density_mol_m3": rho,
        "mass_density_kg_m3": rho * sum(v*c["molar_mass_kg_mol"] for v, c in zip(x, cs)),
        "ln_fugacity_coefficients": np.log(model.get_fugacity_coefficients(temperature, rho*x)).tolist(),
        "residual_enthalpy_j_mol": R*temperature*(ar10+ar01),
        "residual_entropy_j_mol_k": R*(ar10-ar00+np.log(z)),
    }}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-data", action="store_true", help="also regenerate qualified constant selection")
    args = parser.parse_args()
    actual_versions = {"teqp": teqp.__version__, "CoolProp": CoolProp.__version__, "numpy": np.__version__, "scipy": scipy.__version__}
    expected_versions = {"teqp": "0.23.2", "CoolProp": "8.0.0", "numpy": "2.5.4", "scipy": "1.18.1"}
    if actual_versions != expected_versions:
        raise RuntimeError(f"Reference generation requires the pinned versions: {expected_versions}; got {actual_versions}")
    data = generate_data() if args.refresh_data else json.loads((ROOT / "data/components.json").read_text())
    records = {c["id"]: c for c in data["components"]}
    cases = []
    for model in ["PR", "PR-BM"]:
        cases.extend([
            request(model, ["methane"], [1.0], 300.0, 1e5),
            request(model, ["methane"], [1.0], 150.0, 1e6),
            request(model, ["methane"], [1.0], 150.0, 1e6, "smallest"),
            request(model, ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"], [0.9, 0.05, 0.02, 0.02, 0.01], 288.15, 5e6),
            request(model, ["methane", "ethane"], [0.7, 0.3], 250.0, 3e6, kij=0.02),
            request(model, ["methane", "ethane"], [0.7, 0.3], 500.0, 1e7, kij=0.02),
            request(model, ["propane"], [1.0], 300.0, 1e6, "smallest"),
            request(model, ["nitrogen", "methane"], [0.5, 0.5], 115.0, 1e6),
        ])
    fixtures = []
    for index, req in enumerate(cases):
        result = reference(req, records)
        result["id"] = f"cubic-{index:03d}"
        fixtures.append(result)
    write_json(ROOT / "validation/cubic-reference.json", {
        "purpose": "implementation verification, homogeneous roots only; no stability or experimental accuracy claim",
        "versions": {"teqp": teqp.__version__, "CoolProp": CoolProp.__version__, "numpy": np.__version__, "scipy": scipy.__version__, "python": platform.python_version()},
        "data_version": data["version"],
        "gas_constant_j_mol_k": R,
        "source": "https://pages.nist.gov/teqp-docs/en/latest/models/cubics.html",
        "generation": "python validation/generate_cubic_references.py" + (" --refresh-data" if args.refresh_data else ""),
        "tolerances": {"relative": 2e-8, "z_absolute": 1e-10, "density_absolute": 1e-7, "ln_phi_absolute": 1e-9, "enthalpy_absolute_j_mol": 1e-5, "entropy_absolute_j_mol_k": 1e-7},
        "cases": fixtures,
    })
    print(f"Generated {len(fixtures)} independent homogeneous reference cases")


if __name__ == "__main__":
    main()
