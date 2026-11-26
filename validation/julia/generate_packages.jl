# Independent #226 package-operation references. No Rust calls and no Rust outputs.
#
# - GERG-2008: native Clapeyron GERG2008 bubble/dew pressures and temperatures
#   (Clapeyron's own solvers, no supplied solutions) inside the package domain
#   230..350 K, 0.1..10 MPa; isentropic-efficiency compressor and expander from
#   Clapeyron entropy and enthalpy. Enthalpies carry constant per-component
#   offsets to the NIST GERG ideal anchor (as check-gerg.jl); offsets cancel in
#   every entropy difference and in the isentropic outlet temperature.
# - NRTL-IG: the A7 Clapeyron gamma-phi adapter of generate_nrtl.jl (included,
#   unchanged) for pump, compressor and expander outlets; outlet temperatures
#   from bisection on the adapter's TP entropy and enthalpy.
# - EOS-CG-2021: native Clapeyron EOS_CG (Neumann et al. 2023) CO2/H2 gas-branch
#   isentropic outlets and PS states inside the package domain 330..350 K,
#   0.1..2 MPa. The package's caloric zero differs by constant per-component
#   offsets, so only temperatures and enthalpy/entropy differences are stored.
using Clapeyron, JSON3, SHA

include(joinpath(@__DIR__, "generate_nrtl.jl"))   # NRTL adapter; its main() does not run

const RGERG = 8.314472
const GERGDATA = joinpath(ROOT, "data/gerg-2008.json")
const NIST = JSON3.read(read(GERGDATA, String))
const GERG_ALIASES = Dict("carbon-dioxide" => "carbon dioxide")

# --- GERG-2008 ------------------------------------------------------------------
function nist_ideal_anchor(name)
    i = findfirst(==(name), String.(NIST.components)) + 1
    n, theta = NIST.tables.n0i[i], NIST.tables.th0i[i]
    t = 298.15
    a = log(101.325 / (RGERG * t)) + n[2] + n[3] / t - n[4] * log(t)
    at = n[4] + n[3] / t
    for j in 4:7
        theta[j+1] == 0 && continue
        v, c = theta[j+1] / t, n[j+1]
        if j == 4 || j == 6
            a += c * log(sinh(v)); at += c * v / tanh(v)
        else
            a -= c * log(cosh(v)); at -= c * v * tanh(v)
        end
    end
    return RGERG * t * (1 + at), RGERG * (at - a)
end

function gerg(ids)
    m = GERG2008([get(GERG_ALIASES, id, id) for id in ids])
    Clapeyron.Rgas(m) == RGERG || error("unmatched GERG gas constant")
    offsets = map(ids) do id
        ideal = Clapeyron.idealmodel(GERG2008([get(GERG_ALIASES, id, id)]))
        h, s = nist_ideal_anchor(id)
        v0 = RGERG * 298.15 / 101325.0
        (h - Clapeyron.VT_enthalpy(ideal, v0, 298.15, [1.0]), s - Clapeyron.VT_entropy(ideal, v0, 298.15, [1.0]))
    end
    return m, offsets
end

gerg_h(m, off, p, T, z) = enthalpy(m, p, T, z; phase=:vapour) + sum(z[i] * off[i][1] for i in eachindex(z))
gerg_s(m, off, p, T, z) = entropy(m, p, T, z; phase=:vapour) + sum(z[i] * off[i][2] for i in eachindex(z))

in_gerg(T, p) = 230.0 <= T <= 350.0 && 1e5 <= p <= 1e7

function chemical_residual(m, T, vl, vv, x, y)
    ml = Clapeyron.VT_chemical_potential(m, vl, T, x)
    mv = Clapeyron.VT_chemical_potential(m, vv, T, y)
    maximum(abs(ml[i] - mv[i]) for i in eachindex(x) if x[i] > 0 && y[i] > 0) / (RGERG * T)
end

function gerg_saturation(ids, z, kind, fixed)
    m, _ = gerg(ids)
    if kind == "bubble-pressure"
        p, vl, vv, y = bubble_pressure(m, fixed, z, Clapeyron.ChemPotBubblePressure(; atol=1e-12, rtol=1e-14)); T = fixed; x = z
    elseif kind == "dew-pressure"
        p, vl, vv, x = dew_pressure(m, fixed, z, Clapeyron.ChemPotDewPressure(; atol=1e-12, rtol=1e-14)); T = fixed; y = z
    elseif kind == "bubble-temperature"
        T, vl, vv, y = bubble_temperature(m, fixed, z, Clapeyron.ChemPotBubbleTemperature(; atol=1e-12, rtol=1e-14)); p = fixed; x = z
    else
        T, vl, vv, x = dew_temperature(m, fixed, z, Clapeyron.ChemPotDewTemperature(; atol=1e-12, rtol=1e-14)); p = fixed; y = z
    end
    x, y = collect(Float64, x), collect(Float64, y)
    in_gerg(T, p) || error("$ids $kind outside the GERG package domain: T=$T, p=$p")
    vv > 1.5 * vl || error("$ids $kind near-critical or trivial: vl=$vl, vv=$vv")
    r = chemical_residual(m, T, vl, vv, x, y)
    r <= 1e-8 || error("$ids $kind chemical-potential residual $r")
    for (v, w) in ((vl, x), (vv, y))
        q = Clapeyron.VT_pressure(m, v, T, w)
        abs(q - p) <= 1e-6 * p || error("$ids $kind pressure mismatch $q vs $p")
    end
    Dict("ids" => ids, "z" => z, "kind" => kind, "T" => T, "p" => p,
         "liquid" => x, "vapor" => y, "vl" => vl, "vv" => vv)
end

function bisect64(f, lo, hi)
    flo = f(lo)
    sign(flo) != sign(f(hi)) || error("reference outlet not bracketed in $lo..$hi K")
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        fm = f(mid)
        (sign(fm) == sign(flo)) ? (lo = mid; flo = fm) : (hi = mid)
    end
    return (lo + hi) / 2
end

# One isentropic-efficiency machine: h(T, p) and s(T, p) at the feed composition.
function machine(id, kind, h, s, T1, p1, p2, eta, n, bracket)
    s1, h1 = s(T1, p1), h(T1, p1)
    Ts = bisect64(T -> s(T, p2) - s1, bracket...)
    dhs = h(Ts, p2) - h1
    dh = kind == "expander" ? eta * dhs : dhs / eta
    T2 = bisect64(T -> h(T, p2) - (h1 + dh), bracket...)
    Dict("id" => id, "kind" => kind, "efficiency" => eta, "molar_flow_mol_s" => n,
         "inlet" => Dict("T" => T1, "p" => p1, "h" => h1, "s" => s1),
         "isentropic_outlet" => Dict("T" => Ts, "p" => p2, "h" => h(Ts, p2), "s" => s(Ts, p2)),
         "outlet" => Dict("T" => T2, "p" => p2, "h" => h(T2, p2), "s" => s(T2, p2)),
         "isentropic_power_W" => n * dhs, "power_W" => n * dh)
end

const NG = ["methane", "ethane", "propane", "nitrogen", "carbon-dioxide"]
const NG_Z = [0.88, 0.06, 0.03, 0.02, 0.01]

function gerg_units()
    m, off = gerg(NG)
    h(T, p) = gerg_h(m, off, p, T, NG_Z)
    s(T, p) = gerg_s(m, off, p, T, NG_Z)
    units = [machine("gerg-ng-compressor", "compressor", h, s, 260.0, 2e6, 4e6, 0.75, 10.0, (230.0, 350.0)),
             machine("gerg-ng-expander", "expander", h, s, 340.0, 6e6, 3e6, 0.8, 10.0, (230.0, 350.0))]
    for u in units
        u["ids"], u["z"] = NG, NG_Z
        for key in ("inlet", "isentropic_outlet", "outlet")
            T, p = u[key]["T"], u[key]["p"]
            in_gerg(T, p) || error("$(u["id"]) $key outside the GERG domain")
            Clapeyron.isstable(m, p, T, NG_Z) || error("$(u["id"]) $key is not a stable single phase")
            u[key]["phases"] = 1
        end
    end
    return units
end

# --- NRTL-IG (gamma-phi adapter of generate_nrtl.jl) -----------------------------
function nrtl_state(T, p, z)
    c = tp_case(T, p, z)
    s = sum(ph["fraction"] * ph["state"]["s"] for ph in c["phases"])
    beta = sum(ph["phase"] == "vapor" ? ph["fraction"] : 0.0 for ph in c["phases"])
    return Float64(c["h"]), Float64(s), length(c["phases"]), beta
end

function nrtl_units()
    z = [0.5, 0.5]
    cases = [("nrtl-pump", "pump", 320.0, 101325.0, 500000.0, 0.7, 5.0),
             ("nrtl-vapor-compressor", "compressor", 380.0, 100000.0, 200000.0, 0.8, 2.0),
             ("nrtl-vapor-expander", "expander", 420.0, 300000.0, 101325.0, 0.8, 2.0)]
    map(cases) do (id, kind, T1, p1, p2, eta, n)
        u = machine(id, kind, (T, p) -> nrtl_state(T, p, z)[1], (T, p) -> nrtl_state(T, p, z)[2],
                    T1, p1, p2, eta, n, (273.16, 450.0))
        u["ids"], u["z"] = IDS, z
        for key in ("inlet", "isentropic_outlet", "outlet")
            _, _, phases, beta = nrtl_state(u[key]["T"], u[key]["p"], z)
            u[key]["phases"] = phases
            u[key]["vapor_fraction"] = beta
        end
        u
    end
end

# --- EOS-CG-2021 CO2/H2 gas branch -------------------------------------------------
const EOSCG_IDS = ["carbon-dioxide", "hydrogen"]
const EOSCG_Z = [0.3, 0.7]

function eoscg_units()
    m = EOS_CG(["carbon dioxide", "hydrogen"])
    h(T, p) = enthalpy(m, p, T, EOSCG_Z; phase=:vapour)
    s(T, p) = entropy(m, p, T, EOSCG_Z; phase=:vapour)
    units = [machine("eoscg-compressor", "compressor", h, s, 331.0, 1.0e6, 1.15e6, 0.8, 1.0, (330.0, 350.0)),
             machine("eoscg-expander", "expander", h, s, 349.0, 1.5e6, 1.25e6, 0.8, 1.0, (330.0, 350.0))]
    for u in units
        for key in ("inlet", "isentropic_outlet", "outlet")
            T, p = u[key]["T"], u[key]["p"]
            (330.0 <= T <= 350.0 && 1e5 <= p <= 2e6) || error("$(u["id"]) $key outside the EOS-CG package domain")
            Clapeyron.isstable(m, p, T, EOSCG_Z) || error("$(u["id"]) $key is not a stable gas")
            # Absolute h and s use Clapeyron's reference: keep only differences.
            u[key]["h_minus_inlet"] = u[key]["h"] - u["inlet"]["h"]
            u[key]["s_minus_inlet"] = u[key]["s"] - u["inlet"]["s"]
        end
        for key in ("isentropic_outlet", "outlet", "inlet")
            delete!(u[key], "h"); delete!(u[key], "s")
        end
        u["ids"], u["z"] = EOSCG_IDS, EOSCG_Z
    end
    return units
end

function packages_main()
    length(ARGS) == 2 && ARGS[1] == "--output" || error("usage: julia --project=validation/julia validation/julia/generate_packages.jl --output candidate.json")
    @assert pkgversion(Clapeyron) == v"0.6.30"
    saturation = Any[
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "bubble-pressure", 240.0),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "bubble-pressure", 270.0),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "bubble-pressure", 300.0),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "bubble-temperature", 5e6),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "bubble-temperature", 7e6),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "dew-temperature", 1e6),
        gerg_saturation(["methane", "propane"], [0.3, 0.7], "dew-temperature", 3e6),
        gerg_saturation(["ethane", "propane"], [0.5, 0.5], "bubble-pressure", 260.0),
        gerg_saturation(["ethane", "propane"], [0.5, 0.5], "dew-pressure", 260.0),
        gerg_saturation(["ethane", "propane"], [0.5, 0.5], "bubble-temperature", 1.5e6),
        gerg_saturation(["ethane", "propane"], [0.5, 0.5], "dew-temperature", 1.5e6),
        gerg_saturation(["methane", "ethane", "propane"], [0.2, 0.45, 0.35], "bubble-temperature", 6e6),
        gerg_saturation(["methane", "ethane", "propane"], [0.2, 0.45, 0.35], "dew-temperature", 2.5e6),
    ]
    output = Dict("schema" => "clapeyron-package-operations-v1", "julia_version" => string(VERSION),
        "clapeyron_version" => string(pkgversion(Clapeyron)),
        "clapeyron_commit" => "3440a15ba224664cac5b24441ac8adacbaf326f9",
        "input_sha256" => Dict(replace(relpath(p, ROOT), '\\' => '/') => bytes2hex(sha256(read(p))) for p in [GAS, LIQ, GERGDATA]),
        "generator_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
        "formulation" => "Clapeyron GERG2008 (native, NIST ideal-anchor enthalpy offsets) and EOS_CG (EOS-CG-2021) with Clapeyron bubble/dew solvers and Float64 bisection on entropy/enthalpy; NRTL-IG via the generate_nrtl.jl gamma-phi adapter",
        "scope" => "Implementation agreement for matched formulations inside the package domains; no experimental accuracy claim.",
        # Chosen before any Rust comparison: the A13, A14 and A6 budgets.
        "tolerances" => Dict("temperature" => [2e-6, 0.0], "pressure" => [0.02, 2e-8],
            "composition" => [2e-8, 1e-7], "volume" => [1e-12, 2e-8], "h" => [2e-5, 2e-9],
            "s" => [2e-7, 2e-9], "power" => [2e-5, 2e-8], "fraction" => [2e-8, 0.0]),
        "gerg_saturation" => saturation, "gerg_units" => gerg_units(), "nrtl_units" => nrtl_units(),
        "eoscg_units" => eoscg_units())
    open(ARGS[2], "w") do io
        JSON3.pretty(io, output)
        println(io)
    end
    println("Wrote $(length(saturation)) GERG saturation points and GERG, NRTL and EOS-CG machines to ", ARGS[2])
end

if abspath(PROGRAM_FILE) == @__FILE__
    packages_main()
end
