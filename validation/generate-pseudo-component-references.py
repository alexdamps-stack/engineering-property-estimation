"""Independent pseudo-component references (#229); never invokes the Rust engine.

Two kinds of reference values for the kerosene surrogate case of
docs/pseudo-components-status.md:

1. Surrogate species data from CoolProp's multiparameter (HEOS) fluids: normal
   boiling point, specific gravity at 60 degF/60 degF, molar mass, critical
   constants, acentric factor and ideal-gas cp0. These define the "surrogate"
   representation (defined compounds) and the model-comparison reference for
   the Kesler-Lee ideal-gas cp.
2. An independent plain-Python evaluation of the characterisation
   correlations (Kesler-Lee Tc/Pc/cp, Lee-Kesler omega, Riazi-Daubert M) at
   each species' Tb and SG, for implementation agreement with the Rust code.
   The coefficients are re-typed here from Fahim et al. (2010), chapter 3.

Generate to an explicit candidate path and review before adoption:

    python validation/generate-pseudo-component-references.py --output candidate.json
"""
import argparse
import importlib.metadata
import json
import math
import platform

import CoolProp
import CoolProp.CoolProp as CP

# 60 degF in K and the pressure of the specific-gravity definition.
T60F = (60.0 + 459.67) / 1.8
P_ATM = 101325.0
R = 8.31446261815324  # thermo_types::MOLAR_GAS_CONSTANT (CODATA 2018 exact value)

SPECIES = [
    # (pseudo-component id, CoolProp fluid, surrogate mole fraction)
    ("pseudo:surrogate-m-xylene", "m-Xylene", 0.2),
    ("pseudo:surrogate-n-decane", "n-Decane", 0.4),
    ("pseudo:surrogate-n-dodecane", "n-Dodecane", 0.4),
]
CP0_TEMPERATURES_K = [300.0, 400.0, 500.0, 600.0]

# Implementation agreement (same formulation, independent code), and the
# model-comparison budgets of the surrogate case, declared before the first
# comparison. Relative unless the name says otherwise.
TOLERANCES = {
    "correlation_relative": 1e-12,
    "kesler_lee_cp0_vs_coolprop_relative": 0.05,
    "critical_temperature_vs_surrogate_relative": 0.02,
    "critical_pressure_vs_surrogate_relative": 0.08,
    "acentric_factor_vs_surrogate_relative": 0.10,
    "molar_mass_vs_surrogate_relative": 0.05,
    "bubble_dew_temperature_vs_surrogate_K": 5.0,
    "liquid_mass_density_vs_surrogate_relative": 0.06,
    "specific_enthalpy_change_vs_surrogate_relative": 0.05,
    "recorded_exceedance_reproduction_absolute": 0.001,
}

# Exceedances of the declared budgets found by the first comparison. They are
# recorded, not hidden by a wider budget: the Rust test requires each to stay
# above its budget and reproduce the recorded deviation, and every other
# species/quantity must meet the budget.
RECORDED_EXCEEDANCES = [
    {"species": "pseudo:surrogate-m-xylene", "quantity": "molar_mass_vs_surrogate_relative",
     "representation": "pseudo", "observed": 0.0509,
     "reason": "Riazi-Daubert (1987) M for a pure aromatic with Watson K 10.4; the correlation is fitted to petroleum fractions"},
    {"species": "pseudo:surrogate-m-xylene", "quantity": "kesler_lee_cp0_vs_coolprop_relative",
     "representation": "pseudo", "temperatures_K": [300.0, 400.0, 500.0, 600.0],
     "observed": [0.0813, 0.0768, 0.0847, 0.0959],
     "reason": "Kesler-Lee cp overpredicts a pure aromatic (Watson K 10.4) when evaluated with the correlated omega and M; with the HEOS omega and M (surrogate representation) the deviation is 2.3-3.9 %"},
]


def kesler_lee_tc(tb, sg):
    return 189.8 + 450.6 * sg + (0.4244 + 0.1174 * sg) * tb + (14410.0 - 100688.0 * sg) / tb


def kesler_lee_pc(tb, sg):
    ln_bar = (5.689 - 0.0566 / sg
              - 1e-3 * tb * (0.436392 + 4.12164 / sg + 0.213426 / sg ** 2)
              + 1e-7 * tb ** 2 * (4.75794 + 11.819 / sg + 1.53015 / sg ** 2)
              - 1e-10 * tb ** 3 * (2.45055 + 9.901 / sg ** 2))
    return math.exp(ln_bar) * 1e5


def lee_kesler_omega(tb, tc, pc, k):
    tbr = tb / tc
    if tbr < 0.8:
        num = math.log(P_ATM / pc) - 5.92714 + 6.09648 / tbr + 1.28862 * math.log(tbr) - 0.169347 * tbr ** 6
        den = 15.2518 - 15.6875 / tbr - 13.4721 * math.log(tbr) + 0.43577 * tbr ** 6
        return num / den
    return -7.904 + 0.1352 * k - 0.007465 * k ** 2 + 8.359 * tbr + (1.408 - 0.01063 * k) / tbr


def riazi_daubert_m(tb, sg):
    return 42.965 * math.exp(2.097e-4 * tb - 7.78712 * sg + 2.08476e-3 * tb * sg) * tb ** 1.26007 * sg ** 4.98308 / 1000.0


def kesler_lee_cp_over_r(k, omega, molar_mass):
    cf = ((12.8 - k) * (10.0 - k) / (10.0 * omega)) ** 2
    per_mass = [
        -1.41779 + 0.11828 * k - cf * (1.09223 - 2.48245 * omega),
        -(6.99724 - 8.69326 * k + 0.27715 * k * k) * 1e-4 + cf * (3.434 - 7.14 * omega) * 1e-3,
        -2.2582e-6 + cf * (7.2661 - 9.2561 * omega) * 1e-7,
    ]
    return [a * molar_mass * 1000.0 / R for a in per_mass]


def species_record(pid, fluid, fraction):
    tb = CP.PropsSI("T", "P", P_ATM, "Q", 0, fluid)
    rho = CP.PropsSI("D", "T", T60F, "P", P_ATM, fluid)
    rho_w = CP.PropsSI("D", "T", T60F, "P", P_ATM, "Water")
    sg = rho / rho_w
    state = CP.AbstractState("HEOS", fluid)
    cp0 = []
    for t in CP0_TEMPERATURES_K:
        state.update(CP.DmolarT_INPUTS, 1e-6, t)
        cp0.append(state.cp0molar())
    k = (1.8 * tb) ** (1.0 / 3.0) / sg
    tc = kesler_lee_tc(tb, sg)
    pc = kesler_lee_pc(tb, sg)
    omega = lee_kesler_omega(tb, tc, pc, k)
    m = riazi_daubert_m(tb, sg)
    return {
        "id": pid,
        "coolprop_fluid": fluid,
        "cas": CP.get_fluid_param_string(fluid, "CAS"),
        "surrogate_mole_fraction": fraction,
        "normal_boiling_point_K": tb,
        "specific_gravity_60F": sg,
        "molar_mass_kg_per_mol": CP.PropsSI("molar_mass", fluid),
        "critical_temperature_K": CP.PropsSI("Tcrit", fluid),
        "critical_pressure_Pa": CP.PropsSI("pcrit", fluid),
        "acentric_factor": CP.PropsSI("acentric", fluid),
        "ideal_gas_cp_J_per_mol_K": cp0,
        "correlations": {
            "watson_k": k,
            "critical_temperature_K": tc,
            "critical_pressure_Pa": pc,
            "acentric_factor": omega,
            "molar_mass_kg_per_mol": m,
            "ideal_gas_cp_over_R": kesler_lee_cp_over_r(k, omega, m),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = {
        "schema": "pseudo-component-reference-v1",
        "generator": "validation/generate-pseudo-component-references.py",
        "python": platform.python_version(),
        "packages": {"CoolProp": importlib.metadata.version("CoolProp"), "CoolProp.__version__": CoolProp.__version__},
        "independence": "CoolProp HEOS multiparameter fluids (the reference equations of state cited by CoolProp 7.2.0 for m-xylene, n-decane and n-dodecane) for species data; correlations re-typed in plain Python; no Rust calls",
        "scope": "Kerosene surrogate case of docs/pseudo-components-status.md: species data of the defined-compound surrogate, an independent evaluation of the characterisation correlations, and declared model-comparison budgets",
        "specific_gravity_definition": "rho(T, 101325 Pa) / rho_water(T, 101325 Pa) at T = 60 degF = 288.7056 K, both from CoolProp HEOS",
        "ideal_gas_cp_temperatures_K": CP0_TEMPERATURES_K,
        "comparison_pressure_Pa": P_ATM,
        "tolerances": TOLERANCES,
        "recorded_exceedances": RECORDED_EXCEEDANCES,
        "species": [species_record(*s) for s in SPECIES],
    }
    with open(args.output, "w", encoding="utf-8") as f:
        f.write(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
