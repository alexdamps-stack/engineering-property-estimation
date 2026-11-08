# Independent A7 NRTL / gamma-phi references. No Rust calls and no Rust outputs.
#
# Clapeyron.jl supplies the NRTL and aspenNRTL excess Gibbs energies and activity
# coefficients, the CoolProp-form PolExpSat vapor pressures and PolExpLiquid
# densities, and the ReidIdeal ideal-gas caloric functions. Temperature derivatives
# use ForwardDiff through those Clapeyron functions. Clapeyron 0.6.30's native
# gamma-phi VLE cannot use correlation-based pure fluids (no a_res for
# FluidCorrelation), so this script solves the modified-Raoult equations itself
# in 256-bit BigFloat: y_i p = x_i gamma_i psat_i with an ideal-gas vapor.
using Clapeyron, JSON3, SHA

const FD = Clapeyron.ForwardDiff
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const R = 8.31446261815324
const P0 = 100000.0
const T0 = 298.15
const GAS = joinpath(ROOT, "data/source-extracts/components.provisional.json")
const LIQ = joinpath(ROOT, "data/source-extracts/liquid-nrtl.provisional.json")
const GASDATA = Dict(String(c.id) => c for c in JSON3.read(read(GAS, String)).components)
const LIQDOC = JSON3.read(read(LIQ, String))
const LIQDATA = Dict(String(c.id) => c for c in LIQDOC.pure_liquid)
const BIN = only(LIQDOC.nrtl_binary)
const IDS = ["methanol", "water"]
setprecision(BigFloat, 256)

@assert Clapeyron.R̄ == R
@assert String(BIN.i) == "methanol" && String(BIN.j) == "water"

# --- Pure-component building blocks (Clapeyron models) -----------------------
anc(a, type, tau) = Clapeyron.GenericAncEvaluator(Float64.(a.n), Float64.(a.t),
    Float64(a.reducing_temperature_K), Float64(a.reducing_value), type, tau)
const SAT = Dict(id => PolExpSat(anc(LIQDATA[id].vapor_pressure, :exp, true)) for id in IDS)
const RHO = Dict(id => PolExpLiquid(anc(LIQDATA[id].saturated_liquid_density, :noexp, false)) for id in IDS)
psat(id, T) = first(saturation_pressure(SAT[id], T))
lnpsat(id, T) = log(psat(id, T))
dlnpsat(id, T) = FD.derivative(t -> lnpsat(id, t), T)
d2lnpsat(id, T) = FD.derivative(t -> dlnpsat(id, t), T)

function ideal_model(id)
    c = GASDATA[id]
    a = Float64.(c.ideal_gas.coefficients[1][1:5]) .* R
    m = ReidIdeal([id]; userlocations=(; a=[a[1]], b=[a[2]], c=[a[3]], d=[a[4]], e=[a[5]],
        Mw=[Float64(c.constants.molar_mass_kg_per_mol) * 1000]),
        reference_state=ReferenceState(:ideal_gas; T0=T0, P0=P0, H0=0.0, S0=0.0))
    @assert abs(enthalpy(m, P0, T0)) < 1e-8 && abs(entropy(m, P0, T0)) < 1e-8
    return m
end
const IDEAL = Dict(id => ideal_model(id) for id in IDS)
h_ig(id, T) = enthalpy(IDEAL[id], P0, T)
s_ig(id, T) = entropy(IDEAL[id], P0, T)           # pure, at p0
cp_ig(id, T) = isobaric_heat_capacity(IDEAL[id], P0, T)
mw(id) = Float64(GASDATA[id].constants.molar_mass_kg_per_mol)

# Clausius-Clapeyron pure liquid consistent with ideal-gas vapor, zero liquid volume.
h_l(id, T) = h_ig(id, T) - R * T^2 * dlnpsat(id, T)
s_l(id, T) = s_ig(id, T) - R * (lnpsat(id, T) - log(P0)) - R * T * dlnpsat(id, T)
cp_l(id, T) = cp_ig(id, T) - R * (2T * dlnpsat(id, T) + T^2 * d2lnpsat(id, T))

# --- NRTL models ---------------------------------------------------------------
const B12 = Float64(BIN.tau_ij.t1_K)
const B21 = Float64(BIN.tau_ji.t1_K)
const ALPHA = Float64(BIN.alpha_ij.a0)
@assert BIN.tau_ij.t0 == 0 && BIN.tau_ij.t2 == 0 && BIN.tau_ij.t3_per_K == 0 && BIN.alpha_ij.a1_per_K == 0
const MW_G = [mw(id) * 1000 for id in IDS]
const ACT = NRTL(IDS; puremodel=BasicIdeal, userlocations=(a=[0.0 0.0; 0.0 0.0],
    b=[0.0 B12; B21 0.0], c=[0.0 ALPHA; ALPHA 0.0], Mw=MW_G))

# Synthetic directed ternary exercising every temperature term (formulation check).
const TERNARY_IDS = ["comp-a", "comp-b", "comp-c"]
const TERNARY = (
    t0=[0.0 0.8 -0.4; 1.3 0.0 0.25; 0.6 -0.9 0.0],
    t1=[0.0 -150.0 220.0; 310.0 0.0 -75.0; -40.0 500.0 0.0],
    t2=[0.0 0.05 -0.02; -0.03 0.0 0.04; 0.01 -0.06 0.0],
    t3=[0.0 -0.001 0.0015; 0.002 0.0 -0.0005; 0.0007 0.0012 0.0],
    a0=[0.0 0.2 0.35; 0.2 0.0 0.3; 0.35 0.3 0.0],
    a1=[0.0 0.0004 -0.0002; 0.0004 0.0 0.0001; -0.0002 0.0001 0.0],
)
const ACT3 = aspenNRTL(TERNARY_IDS; puremodel=BasicIdeal, userlocations=TERNARY)

gE_over_RT(model, T, x) = Clapeyron.excess_gibbs_free_energy(model, 1e5, T, x) / (R * T)
lngamma(model, T, x) = log.(activity_coefficient(model, 1e5, T, x))

function excess_record(model, T, x)
    g(t) = gE_over_RT(model, t, x)
    d1 = FD.derivative(g, T)
    d2 = FD.derivative(t -> FD.derivative(g, t), T)
    hE = -R * T^2 * d1
    gE = R * T * g(T)
    # Cross-check Clapeyron's own excess-enthalpy path for the activity model.
    @assert isapprox(Clapeyron.excess(model, 1e5, T, x, enthalpy), hE; rtol=1e-9, atol=1e-8)
    Dict("T" => T, "x" => x, "ln_gamma" => lngamma(model, T, x),
        "d_ln_gamma_dT" => [FD.derivative(t -> lngamma(model, t, x)[i], T) for i in eachindex(x)],
        "g_excess" => gE, "h_excess" => hE, "s_excess" => (hE - gE) / T,
        "cp_excess" => -R * (2T * d1 + T^2 * d2))
end

# --- Gamma-phi phase properties -------------------------------------------------
function liquid_props(T, x)
    ex = excess_record(ACT, T, x)
    h = sum(x[i] * h_l(IDS[i], T) for i in 1:2) + ex["h_excess"]
    s = sum(x[i] * s_l(IDS[i], T) for i in 1:2) - R * sum(xi > 0 ? xi * log(xi) : 0.0 for xi in x) + ex["s_excess"]
    cp = sum(x[i] * cp_l(IDS[i], T) for i in 1:2) + ex["cp_excess"]
    v = sum(x[i] / (1 / volume(RHO[IDS[i]], P0, T, [1.0])) for i in 1:2)
    Dict("h" => h, "s" => s, "cp" => cp, "molar_volume" => v, "ln_gamma" => ex["ln_gamma"])
end

function vapor_props(T, p, y)
    h = sum(y[i] * h_ig(IDS[i], T) for i in 1:2)
    s = sum(y[i] * s_ig(IDS[i], T) for i in 1:2) - R * log(p / P0) - R * sum(yi > 0 ? yi * log(yi) : 0.0 for yi in y)
    cp = sum(y[i] * cp_ig(IDS[i], T) for i in 1:2)
    Dict("h" => h, "s" => s, "cp" => cp, "molar_volume" => R * T / p)
end

# --- Modified-Raoult solves in BigFloat ----------------------------------------
big_gamma(T, x) = lngamma(ACT, BigFloat(T), BigFloat.(x))
big_lnpsat(T) = [lnpsat(id, BigFloat(T)) for id in IDS]

function bubble_p(T, x)
    lg, lp = big_gamma(T, x), big_lnpsat(T)
    part = [x[i] > 0 ? x[i] * exp(lg[i] + lp[i]) : big(0.0) for i in 1:2]
    p = sum(part)
    return p, part ./ p
end

function dew_p(T, y)
    lp = big_lnpsat(T)
    lg = zeros(BigFloat, 2)
    for _ in 1:2000
        terms = [y[i] > 0 ? y[i] / exp(lg[i] + lp[i]) : big(0.0) for i in 1:2]
        p = 1 / sum(terms)
        x = terms .* p
        new = big_gamma(T, x)
        if maximum(abs.(new .- lg)) < big(1e-60)
            return p, x
        end
        lg = new
    end
    error("reference dew pressure did not converge at T=$T")
end

function bisect(f, lo, hi)
    lo, hi = BigFloat(lo), BigFloat(hi)
    flo = f(lo)
    @assert sign(flo) != sign(f(hi))
    for _ in 1:400
        mid = (lo + hi) / 2
        fm = f(mid)
        (sign(fm) == sign(flo)) ? (lo = mid; flo = fm) : (hi = mid)
    end
    return (lo + hi) / 2
end

bubble_T(p, x) = bisect(T -> log(first(bubble_p(T, x)) / p), 273.16, 450.0)
dew_T(p, y) = bisect(T -> log(first(dew_p(T, y)) / p), 273.16, 450.0)

function tp_split(T, p, z)
    lp = big_lnpsat(T)
    lg = big_gamma(T, z)
    beta = big(0.5)
    x = y = z
    for _ in 1:5000
        K = exp.(lg .+ lp .- log(BigFloat(p)))
        rr(b) = sum(z[i] * (K[i] - 1) / (1 + b * (K[i] - 1)) for i in 1:2)
        beta = bisect(rr, 1 / (1 - maximum(K)) + big(1e-70), 1 / (1 - minimum(K)) - big(1e-70))
        x = [z[i] / (1 + beta * (K[i] - 1)) for i in 1:2]
        x ./= sum(x)
        y = K .* x
        y ./= sum(y)
        new = big_gamma(T, x)
        if maximum(abs.(new .- lg)) < big(1e-60)
            break
        end
        lg = new
    end
    r = maximum(abs.(log.(y) .+ log(BigFloat(p)) .- log.(x) .- big_gamma(T, x) .- lp))
    @assert r < 1e-50 "reference TP split did not converge"
    return beta, x, y
end

f64(v) = Float64.(v)

function tp_case(T, p, z)
    pb, _ = bubble_p(T, z)
    pd, _ = dew_p(T, z)
    phases = Any[]
    if p >= pb
        push!(phases, Dict("phase" => "liquid", "fraction" => 1.0, "composition" => z, "state" => liquid_props(T, z)))
    elseif p <= pd
        push!(phases, Dict("phase" => "vapor", "fraction" => 1.0, "composition" => z, "state" => vapor_props(T, p, z)))
    else
        beta, x, y = tp_split(T, p, z)
        x, y = f64(x), f64(y)
        push!(phases, Dict("phase" => "vapor", "fraction" => Float64(beta), "composition" => y, "state" => vapor_props(T, p, y)))
        push!(phases, Dict("phase" => "liquid", "fraction" => Float64(1 - beta), "composition" => x, "state" => liquid_props(T, x)))
    end
    h = sum(ph["fraction"] * ph["state"]["h"] for ph in phases)
    Dict("T" => T, "p" => p, "z" => z, "bubble_pressure" => Float64(pb), "dew_pressure" => Float64(pd),
         "h" => h, "phases" => phases)
end

function pure_saturation_case(id, p, beta)
    z = id == "methanol" ? [1.0, 0.0] : [0.0, 1.0]
    T = Float64(bisect(t -> lnpsat(id, t) - log(BigFloat(p)), 273.16, 450.0))
    l, v = liquid_props(T, z), vapor_props(T, p, z)
    Dict("id" => "pure-$id-$(p)-$(beta)", "p" => p, "z" => z, "T" => T, "vapor_fraction" => beta,
         "target_h" => beta * v["h"] + (1 - beta) * l["h"], "h_liquid" => l["h"], "h_vapor" => v["h"])
end

function main()
    length(ARGS) == 2 && ARGS[1] == "--output" || error("usage: julia --project=validation/julia validation/julia/generate_nrtl.jl --output candidate.json")
    @assert pkgversion(Clapeyron) == v"0.6.30"
    xs = [[0.0, 1.0], [0.05, 0.95], [0.3, 0.7], [0.5, 0.5], [0.8, 0.2], [1.0, 0.0]]
    activity = [excess_record(ACT, T, x) for T in (300.0, 337.69, 350.0, 373.13, 420.0) for x in xs]
    ternary = [excess_record(ACT3, T, x) for T in (280.0, 330.0, 400.0)
               for x in ([0.2, 0.3, 0.5], [0.6, 0.4, 0.0], [0.1, 0.1, 0.8])]
    pure = Any[]
    for id in IDS, T in (273.16, 300.0, 337.69, 373.13, 450.0)
        push!(pure, Dict("id" => id, "T" => T, "psat" => psat(id, T), "dlnpsat_dT" => dlnpsat(id, T),
            "d2lnpsat_dT2" => d2lnpsat(id, T), "h_liquid" => h_l(id, T), "s_liquid" => s_l(id, T),
            "cp_liquid" => cp_l(id, T), "liquid_molar_density" => 1 / volume(RHO[id], P0, T, [1.0]),
            "h_ideal" => h_ig(id, T)))
    end
    saturation = Any[]
    for T in (330.0, 350.0, 370.0), x1 in (0.1, 0.5, 0.9)
        z = [x1, 1 - x1]
        pb, y = bubble_p(T, z)
        pd, xd = dew_p(T, z)
        push!(saturation, Dict("kind" => "bubble-pressure", "T" => T, "p" => Float64(pb), "liquid" => z, "vapor" => f64(y)))
        push!(saturation, Dict("kind" => "dew-pressure", "T" => T, "p" => Float64(pd), "liquid" => f64(xd), "vapor" => z))
    end
    for p in (50000.0, 101325.0, 300000.0), x1 in (0.1, 0.5, 0.9)
        z = [x1, 1 - x1]
        tb = bubble_T(p, z); _, y = bubble_p(tb, z)
        td = dew_T(p, z); _, xd = dew_p(td, z)
        push!(saturation, Dict("kind" => "bubble-temperature", "T" => Float64(tb), "p" => p, "liquid" => z, "vapor" => f64(y)))
        push!(saturation, Dict("kind" => "dew-temperature", "T" => Float64(td), "p" => p, "liquid" => f64(xd), "vapor" => z))
    end
    tp = Any[tp_case(T, p, z) for (T, p, z) in (
        (350.0, 70000.0, [0.5, 0.5]), (355.0, 101325.0, [0.3, 0.7]), (345.0, 101325.0, [0.7, 0.3]),
        (365.0, 101325.0, [0.1, 0.9]), (400.0, 300000.0, [0.4, 0.6]), (400.0, 450000.0, [0.4, 0.6]),
        (320.0, 101325.0, [0.5, 0.5]), (380.0, 101325.0, [0.5, 0.5]))]
    ph = Any[]
    for (i, c) in enumerate(tp)
        push!(ph, Dict("id" => "tp-$i", "p" => c["p"], "z" => c["z"], "T" => c["T"], "target_h" => c["h"],
            "vapor_fraction" => sum(ph["phase"] == "vapor" ? ph["fraction"] : 0.0 for ph in c["phases"])))
    end
    for id in IDS, beta in (0.0, 0.25, 0.75, 1.0)
        push!(ph, pure_saturation_case(id, 101325.0, beta))
    end
    output = Dict("schema" => "clapeyron-nrtl-reference-v1", "julia_version" => string(VERSION),
        "clapeyron_version" => string(pkgversion(Clapeyron)),
        "clapeyron_commit" => "3440a15ba224664cac5b24441ac8adacbaf326f9",
        "input_sha256" => Dict(relpath(p, ROOT) => bytes2hex(sha256(read(p)))  for p in [GAS, LIQ]),
        "generator_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
        "formulation" => "Clapeyron NRTL (tau=a+b/T, alpha=c) and aspenNRTL (tau=t0+t1/T+t2 lnT+t3 T, alpha=a0+a1 T); PolExpSat/PolExpLiquid CoolProp ancillaries; ReidIdeal NASA7 low-range Cp with pure ideal h=s=0 at 298.15 K, 1 bar; modified Raoult with ideal-gas vapor, no Poynting; Clausius-Clapeyron pure liquids; ForwardDiff temperature derivatives; BigFloat(256) VLE solves",
        "scope" => "Implementation agreement for matched formulation and parameters; no experimental accuracy claim.",
        # Chosen before Rust comparisons.
        "tolerances" => Dict("ln_gamma" => [1e-12, 1e-11], "d_ln_gamma_dT" => [1e-14, 1e-9],
            "g_excess" => [1e-9, 1e-11], "h_excess" => [1e-8, 1e-9], "s_excess" => [1e-10, 1e-9],
            "cp_excess" => [1e-7, 1e-7], "psat" => [0.0, 1e-12], "dlnpsat_dT" => [0.0, 1e-10],
            "d2lnpsat_dT2" => [0.0, 1e-8], "h" => [2e-5, 2e-9], "s" => [2e-7, 2e-9], "cp" => [2e-5, 1e-8],
            "liquid_molar_density" => [0.0, 1e-12], "molar_volume" => [0.0, 1e-9],
            "pressure" => [0.0, 1e-10], "temperature" => [1e-8, 0.0], "composition" => [1e-10, 1e-9],
            "fraction" => [1e-9, 0.0]),
        "activity" => activity, "ternary" => ternary, "ternary_parameters" => Dict(String(k) => [collect(v[i, :]) for i in 1:3] for (k, v) in pairs(TERNARY)),  # rows: [i][j] = parameter_ij
        "ternary_ids" => TERNARY_IDS, "pure" => pure, "saturation" => saturation, "tp" => tp, "ph" => ph)
    open(ARGS[2], "w") do io
        JSON3.pretty(io, output)
        println(io)
    end
    println("Wrote $(length(activity)) binary + $(length(ternary)) ternary activity, $(length(pure)) pure, $(length(saturation)) saturation, $(length(tp)) TP and $(length(ph)) PH NRTL states to ", ARGS[2])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
