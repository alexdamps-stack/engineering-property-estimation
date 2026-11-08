"""Cross-check the independent NIST GERG fixture with teqp 0.23.2 AD.

No Rust calls, no output mutation. Uses the pinned reference-requirements.txt.
teqp's ideal model exposes its ideal Helmholtz energy through the Ar API.
"""
import json
import math
from pathlib import Path
import numpy as np
import teqp

assert teqp.__version__ == "0.23.2", teqp.__version__
fixture = json.loads((Path(__file__).resolve().parent/"fixtures/gerg-v1.json").read_text())
data = json.loads((Path(__file__).resolve().parents[1]/"data/gerg-2008.json").read_text())


def ideal_offsets(name):
    # teqp recomputes integration constants for exact h0=s0=0; NIST keeps
    # rounded constants. Map teqp to the NIST reference with A+B/T only.
    # This cannot change heat capacities or residual derivatives.
    i = data["components"].index(name)+1
    n, theta = data["tables"]["n0i"][i], data["tables"]["th0i"][i]
    t0 = 298.15
    d0 = 101.325/(8.314472*t0)
    a = math.log(d0)+n[1]+n[2]/t0-n[3]*math.log(t0)
    at = n[3]+n[2]/t0
    for j in range(4,8):
        if theta[j] == 0:
            continue
        v = theta[j]/t0
        if j in (4,6):
            a += n[j]*math.log(math.sinh(v))
            at += n[j]*v/math.tanh(v)
        else:
            a -= n[j]*math.log(math.cosh(v))
            at -= n[j]*v*math.tanh(v)
    return a-at, t0*(at+1)


offsets = {name: ideal_offsets(name) for name in data["components"]}
aliases = {"carbon-dioxide":"carbondioxide", "carbon-monoxide":"carbonmonoxide",
           "hydrogen-sulfide":"hydrogensulfide"}
keys = {"ar":"00", "ar_d":"01", "ar_dd":"02", "ar_t":"10", "ar_tt":"20", "ar_dt":"11"}
max_errors = {}
for case in fixture["cases"]:
    ids = list(case["composition"])
    x = np.array(list(case["composition"].values()))
    names = [aliases.get(name,name) for name in ids]
    t, rho = case["temperature_k"], case["molar_density"]
    residual = teqp.make_model(dict(kind="GERG2008resid",model=dict(names=names)))
    ideal = teqp.make_model(dict(kind="GERG2008idealgas",model=dict(names=names)))
    for key, derivative in keys.items():
        actual = getattr(residual,"get_Ar"+derivative)(t,rho,x)
        expected = case["expected"][key]
        error = abs(actual-expected)
        max_errors[key] = max(max_errors.get(key,0),error)
        assert error <= 2e-8+2e-10*abs(expected), (case["id"],key,actual,expected)
    for key, derivative in {"a0":"00","a0_t":"10","a0_tt":"20"}.items():
        actual = getattr(ideal,"get_Ar"+derivative)(t,rho,x)
        offset_a = sum(xx*offsets[name][0] for name,xx in zip(ids,x))
        offset_b = sum(xx*offsets[name][1] for name,xx in zip(ids,x))
        if key == "a0":
            actual += offset_a+offset_b/t
        elif key == "a0_t":
            actual += offset_b/t
        expected = case["expected"][key]
        error = abs(actual-expected)
        max_errors[key] = max(max_errors.get(key,0),error)
        assert error <= 2e-8+2e-10*abs(expected), (case["id"],key,actual,expected)
print(f"teqp {teqp.__version__}: {len(fixture['cases'])} ideal/residual derivative states passed")
print(json.dumps(max_errors,indent=2))
