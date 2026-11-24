# Independent #224 VLLE/LLE references. No Rust calls and no Rust outputs.
#
# PR/PR-BM: Clapeyron's native MultiPhaseTPFlash/MichelsenTPFlash supply initial
# phase sets; each candidate is solved to 256-bit precision by Newton's method on
# equal ln fugacities (Clapeyron fugacity coefficients) and component balances,
# and accepted only if an exhaustive composition-grid tangent-plane test (Clapeyron
# fugacities, minimum-Gibbs root) finds no negative TPD. Missing phases found by
# the grid are added and the system is re-solved. Clapeyron 0.6.30's native flash
# alone is not relied on for the phase count: it reports one phase for the
# 14 bar methane/propane/water case.
#
# NRTL: Clapeyron's NRTL supplies activity coefficients and PolExpSat/PolExpLiquid
# the CoolProp ancillaries (as in generate_nrtl.jl); the gamma-phi equations
# (ideal-gas vapor) are solved by the same 256-bit Newton method from grid
# tangent-plane minima. The NRTL pair is synthetic (see vlle-inputs.json).
#
# Invariant PH cases: the binary three-phase temperature at fixed p is solved by
# Newton in 256 bits; phase amounts follow from the stated vapor fraction and the
# component balance, and the enthalpy target from those amounts.
include("generate.jl")
include("generate_nrtl.jl")

const INPUT = joinpath(ROOT, "validation/vlle-inputs.json")
const VIN = JSON3.read(read(INPUT, String))
const NRTL_SYN = NRTL(String.(VIN.nrtl.ids); puremodel=BasicIdeal,
    userlocations=(a=[0.0 0.0; 0.0 0.0], b=[0.0 Float64(VIN.nrtl.b12_K); Float64(VIN.nrtl.b21_K) 0.0],
                   c=[0.0 Float64(VIN.nrtl.alpha); Float64(VIN.nrtl.alpha) 0.0], Mw=[mw(id)*1000 for id in IDS]))
@assert String.(VIN.nrtl.ids) == IDS

kij_matrix(ids) = [[i != j && ((ids[i] == "water") != (ids[j] == "water")) ? 0.5 : 0.0 for j in eachindex(ids)] for i in eachindex(ids)]

# ---------------------------------------------------------------------------
# Generic 256-bit Newton on u = ln n_ki: ln f_i(x_k) = ln f_i(x_1), Σ_k n_ki = z_i.
function polish(lnf, z, n0; maxit=200)
    P, C = size(n0)
    zb = BigFloat.(z)
    function F(u)
        n = exp.(reshape(u, P, C))
        x = [n[k, :] ./ sum(n[k, :]) for k in 1:P]
        L = [lnf(k, x[k]) for k in 1:P]
        r = BigFloat[]
        for k in 2:P, i in 1:C
            push!(r, L[k][i] - L[1][i])
        end
        for i in 1:C
            push!(r, sum(n[:, i]) - zb[i])
        end
        r
    end
    u = vec(log.(BigFloat.(n0)))
    for _ in 1:maxit
        f = F(u)
        nf = maximum(abs, f)
        nf < big"1e-60" && return reshape(exp.(u), P, C)
        h = big"1e-35"
        J = reduce(hcat, [(F(u .+ h .* ((1:length(u)) .== j)) .- F(u .- h .* ((1:length(u)) .== j))) ./ (2h) for j in eachindex(u)])
        du = -(J \ f)
        lambda = big(1.0)
        accepted = false
        while lambda > 1e-12
            cand = u .+ lambda .* du
            fc = try F(cand) catch; nothing end
            if fc !== nothing && all(isfinite, fc) && sum(abs2, fc) < (1 - 1e-4 * lambda) * sum(abs2, f)
                u = cand
                accepted = true
                break
            end
            lambda /= 2
        end
        accepted || error("Newton line search failed (residual $(Float64(nf)))")
    end
    error("Newton did not converge")
end

# Nonnegative least-squares phase amounts for compositions X (rows) and feed z.
function lever(X, z)
    keep = collect(1:size(X, 1))
    while true
        A = permutedims(X[keep, :])
        b = A \ z
        all(b .> 0) && return keep, b
        deleteat!(keep, argmin(b))
        isempty(keep) && error("no nonnegative lever solution")
    end
end

# Composition grid on the simplex, dense near the edges.
function grid_values()
    v = vcat(0.0, 10.0 .^ (-12:-3), collect(0.002:0.002:0.998), 1 .- 10.0 .^ (-3:-1:-12), 1.0)
    sort(unique(v))
end
function simplex_grid(C)
    v = grid_values()
    C == 2 && return [[a, 1 - a] for a in v]
    C == 3 || error("grid for C=$C not implemented")
    pts = Vector{Vector{Float64}}()
    for a in v, b in v
        c = 1 - a - b
        c >= -1e-15 && push!(pts, [a, b, max(c, 0.0)])
    end
    pts
end
tm_value(w, lnfw, d) = sum(w[i] > 0 ? w[i] * (lnfw[i] - d[i]) : 0.0 for i in eachindex(w))

# Deepest grid points with negative TPD, at least 0.05 apart (by kind).
function negative_minima(tms; limit=3)
    sorted = sort(tms; by=t -> t[1])
    out = Any[]
    for t in sorted
        t[1] < -1e-9 || break
        any(o -> o[3] == t[3] && maximum(abs.(o[2] .- t[2])) < 0.05, out) && continue
        push!(out, t)
        length(out) == limit && break
    end
    out
end

# ---------------------------------------------------------------------------
# PR / PR-BM.
function pr_lnf(model, p, T, x, phase)
    log.(x) .+ log.(fugacity_coefficient(model, p, T, x; phase))
end
classify(model, p, T, x, v) = abs(v - volume(model, p, T, x; phase=:liquid)) <= abs(v - volume(model, p, T, x; phase=:vapour)) ? :liquid : :vapour

function pr_grid(model, p, T, d)
    [(tm_value(w, log.(max.(w, 1e-300)) .+ log.(fugacity_coefficient(model, p, T, w)), d), w, :fluid)
     for w in simplex_grid(length(d)) if all(w .>= 0)]
end

function pr_solve(model, p, T, z, phases, n0)
    P = length(phases)
    P == 1 && return [1.0], [Float64.(z)]
    pb, Tb = BigFloat(p), BigFloat(T)
    n = polish((k, x) -> pr_lnf(model, pb, Tb, x, phases[k]), z, n0)
    beta = [Float64(sum(n[k, :])) for k in 1:P]
    x = [Float64.(n[k, :] ./ sum(n[k, :])) for k in 1:P]
    return beta, x
end

function pr_case(case)
    ids, z, variant = String.(case.ids), Float64.(case.z), String(case.variant)
    model = matched_model(ids, variant, kij_matrix(ids))
    T, p = Float64(case.T), Float64(case.p)
    candidates = Any[]
    for method in (MultiPhaseTPFlash(), MultiPhaseTPFlash(full_tpd=true), MichelsenTPFlash())
        r = Clapeyron.tp_flash2(model, p, T, z, method)
        push!(candidates, (Float64.(r.fractions), [Float64.(x) for x in r.compositions], Float64.(r.volumes)))
    end
    for (b0, x0, v0) in candidates
        accepted = try
            pr_candidate(model, p, T, z, b0, x0, v0)
        catch e
            println("  candidate rejected: ", sprint(showerror, e)[1:min(200, end)])
            nothing
        end
        accepted === nothing || return pr_record(case, model, ids, z, T, p, accepted...)
    end
    error("no candidate passed the grid stability test for $(case.id)")
end

# Solves a Clapeyron candidate; while the grid finds a negative TPD, splits one of
# the phases with Clapeyron's two-phase MichelsenTPFlash seeded by the grid minimum.
function pr_candidate(model, p, T, z, b0, x0, v0)
    phases = [classify(model, p, T, x0[k], v0[k]) for k in eachindex(x0)]
    beta, x = pr_solve(model, p, T, z, phases, permutedims(reduce(hcat, [b0[k] .* x0[k] for k in eachindex(x0)])))
    for round in 1:3
        ref = argmax(beta)
        d = pr_lnf(model, p, T, x[ref], phases[ref])
        minima = negative_minima(pr_grid(model, p, T, d); limit=1)
        isempty(minima) && return (beta, x, phases)
        length(beta) >= 3 && error("grid indicates a fourth phase")
        w = minima[1][2]
        split = nothing
        for k in sortperm(beta; rev=true)
            try
                r = Clapeyron.tp_flash2(model, p, T, x[k], MichelsenTPFlash(K0=max.(w, 1e-30) ./ x[k]))
                length(r.fractions) == 2 || continue
                others = [j for j in eachindex(x) if j != k]
                nb = vcat(beta[others], beta[k] .* Float64.(r.fractions))
                nx = vcat(x[others], [Float64.(c) for c in r.compositions])
                np = vcat(phases[others], [classify(model, p, T, Float64.(r.compositions[j]), Float64(r.volumes[j])) for j in 1:2])
                split = (nb, nx, np)
                break
            catch
            end
        end
        split === nothing && error("no phase split reproduces the grid instability")
        nb, nx, np = split
        beta, x = pr_solve(model, p, T, z, np, permutedims(reduce(hcat, [nb[j] .* nx[j] for j in eachindex(nx)])))
        phases = np
    end
    error("phase set did not settle")
end

function pr_record(case, model, ids, z, T, p, beta, x, phases)
    for k in eachindex(x)
        # The solved root must be the minimum-Gibbs root at the phase composition.
        v = volume(model, p, T, x[k]; phase=phases[k])
        @assert isapprox(volume(model, p, T, x[k]), v; rtol=1e-9) "$(case.id): phase $k is not on the minimum-Gibbs root"
    end
    out = [Dict("fraction" => beta[k], "composition" => x[k],
                "clapeyron_phase" => String(Clapeyron.VT_identify_phase(model, volume(model, p, T, x[k]; phase=phases[k]), T, x[k])),
                "state" => phase_properties(model, p, T, x[k], phases[k])) for k in eachindex(x)]
    sort!(out; by=ph -> ph["state"]["mass_density"])
    residual = length(x) == 1 ? 0.0 : maximum(maximum(pr_lnf(model, p, T, x[k], phases[k]) .- pr_lnf(model, p, T, x[1], phases[1]); init=-Inf) for k in eachindex(x))
    h = sum(ph["fraction"] * ph["state"]["h"] for ph in out)
    println("PR ", case.id, ": ", length(out), " phases, fractions ", [ph["fraction"] for ph in out])
    Dict("id" => String(case.id), "variant" => String(case.variant), "ids" => ids, "z" => z,
         "kij" => kij_matrix(ids), "T" => T, "p" => p, "h" => h, "phases" => out)
end

# ---------------------------------------------------------------------------
# NRTL gamma-phi (synthetic pair).
syn_lngamma(T, x) = log.(activity_coefficient(NRTL_SYN, 1e5, T, x))
function g_lnf(T, p, x, kind)
    kind == :vapour && return log.(x) .+ log(p)
    log.(x) .+ syn_lngamma(T, x) .+ [lnpsat(id, T) for id in IDS]
end
function syn_liquid_state(T, p, x)
    ex = excess_record(NRTL_SYN, T, x)
    h = sum(x[i] * h_l(IDS[i], T) for i in 1:2) + ex["h_excess"]
    v = sum(x[i] / (1 / volume(RHO[IDS[i]], P0, T, [1.0])) for i in 1:2)
    Dict("h" => h, "molar_volume" => v, "mass_density" => sum(x[i] * mw(IDS[i]) for i in 1:2) / v,
         "ln_phi" => ex["ln_gamma"] .+ [lnpsat(id, T) for id in IDS] .- log(p))
end
function syn_vapor_state(T, p, y)
    s = vapor_props(T, p, y)
    Dict("h" => s["h"], "molar_volume" => s["molar_volume"],
         "mass_density" => sum(y[i] * mw(IDS[i]) for i in 1:2) / s["molar_volume"], "ln_phi" => [0.0, 0.0])
end
syn_state(T, p, x, kind) = kind == :vapour ? syn_vapor_state(T, p, x) : syn_liquid_state(T, p, x)

# Negative local minima of the TPD along the binary composition grid, per kind.
function g_grid(T, p, d)
    out = Any[]
    for kind in (:liquid, :vapour)
        ws = [w for w in simplex_grid(2) if w[1] > 0 && w[2] > 0]
        tms = [tm_value(w, g_lnf(T, p, w, kind), d) for w in ws]
        for i in eachindex(ws)
            left = i == 1 ? Inf : tms[i-1]
            right = i == length(ws) ? Inf : tms[i+1]
            tms[i] <= left && tms[i] <= right && push!(out, (tms[i], ws[i], kind))
        end
    end
    out
end

function g_min_kind(T, p, x)
    gl = sum(x .* g_lnf(T, p, x, :liquid)); gv = sum(x .* g_lnf(T, p, x, :vapour))
    gl <= gv ? :liquid : :vapour
end

function nrtl_case(case)
    z, T, p = Float64.(case.z), Float64.(case.T), Float64.(case.p)
    kinds = [g_min_kind(T, p, z)]
    x = [z]; beta = [1.0]
    for round in 1:4
        ref = argmax(beta)
        d = g_lnf(T, p, x[ref], kinds[ref])
        minima = negative_minima(g_grid(T, p, d))
        isempty(minima) && return nrtl_record(case, z, T, p, beta, x, kinds)
        println("  ", case.id, " round ", round, ": TPD minima ", [(m[3], m[2][1], m[1]) for m in minima])
        round == 4 && break
        # Re-initialise from the current phases and the grid minima.
        # Two or more minima bracket the feed on their own; one minimum joins the current phases.
        allx = vcat(length(minima) >= 2 ? [] : x, [m[2] for m in minima])
        allk = vcat(length(minima) >= 2 ? Symbol[] : kinds, [m[3] for m in minima])
        if length(allx) > 3
            allx = allx[end-2:end]; allk = allk[end-2:end]
        end
        if length(minima) >= 2
            keep, b = lever(permutedims(reduce(hcat, allx)), z)
            allx, allk = allx[keep], allk[keep]
            n0 = permutedims(reduce(hcat, [b[j] .* allx[j] for j in eachindex(keep)]))
        else
            # One new phase: start it with 10 % of the feed; Newton closes the balances.
            amounts = vcat(0.9 .* beta, [0.1])
            n0 = permutedims(reduce(hcat, [amounts[j] .* allx[j] for j in eachindex(allx)]))
        end
        if length(allx) == 1
            x, kinds, beta = allx, allk, [1.0]
            continue
        end
        Tb, pb = BigFloat(T), BigFloat(p)
        n = polish((k, xx) -> g_lnf(Tb, pb, xx, allk[k]), z, n0)
        beta = [Float64(sum(n[k, :])) for k in eachindex(allk)]
        x = [Float64.(n[k, :] ./ sum(n[k, :])) for k in eachindex(allk)]
        kinds = allk
    end
    error("NRTL case $(case.id) did not settle")
end

function nrtl_record(case, z, T, p, beta, x, kinds)
    out = [Dict("fraction" => beta[k], "composition" => x[k], "kind" => String(kinds[k]),
                "state" => syn_state(T, p, x[k], kinds[k])) for k in eachindex(x)]
    sort!(out; by=ph -> (ph["kind"] == "liquid", ph["state"]["mass_density"]))
    h = sum(ph["fraction"] * ph["state"]["h"] for ph in out)
    println("NRTL ", case.id, ": ", [(ph["kind"], ph["fraction"]) for ph in out])
    Dict("id" => String(case.id), "z" => z, "T" => T, "p" => p, "h" => h, "phases" => out)
end

# ---------------------------------------------------------------------------
# Binary three-phase temperature at fixed p: unknowns ln x of the minor component
# in each phase and T; ln f equalities of phases 2, 3 with phase 1.
function three_phase(lnf_all, p, T0, x0)
    # x0: initial mole fraction of component 1 in each of 3 phases.
    u = vcat(BigFloat.(x0), BigFloat(T0))
    F(u) = begin
        xs = [[u[k], 1 - u[k]] for k in 1:3]
        L = [lnf_all(k, u[4], xs[k]) for k in 1:3]
        vcat(L[2] .- L[1], L[3] .- L[1])
    end
    for _ in 1:200
        f = F(u)
        maximum(abs, f) < big"1e-60" && return u
        h = big"1e-35"
        J = reduce(hcat, [(F(u .+ h .* ((1:4) .== j)) .- F(u .- h .* ((1:4) .== j))) ./ (2h) for j in 1:4])
        du = -(J \ f)
        lambda = big(1.0)
        while lambda > 1e-12
            c = u .+ lambda .* du
            if all(0 .< c[1:3] .< 1) && sum(abs2, F(c)) < sum(abs2, f)
                u = c
                break
            end
            lambda /= 2
        end
        lambda > 1e-12 || error("three-phase Newton line search failed")
    end
    error("three-phase Newton did not converge")
end

function invariant_amounts(x1, z, vapor, beta_v)
    # Phases 1, 2 liquids, `vapor` index of the vapor; solve the other two amounts.
    liquids = [k for k in 1:3 if k != vapor]
    a, b = liquids
    rest = z[1] - beta_v * x1[vapor]
    bb = (rest - (1 - beta_v) * x1[a]) / (x1[b] - x1[a])
    beta = zeros(3); beta[vapor] = beta_v; beta[b] = bb; beta[a] = 1 - beta_v - bb
    @assert all(beta .> 0) "invariant amounts not positive: $beta"
    beta
end

function pr_invariant(case)
    ids, z, variant = String.(case.ids), Float64.(case.z), String(case.variant)
    model = matched_model(ids, variant, kij_matrix(ids))
    p = Float64(case.p)
    phases = [:liquid, :liquid, :vapour]   # propane-rich liquid, aqueous liquid, vapor
    lnf_all(k, T, x) = pr_lnf(model, BigFloat(p), T, x, phases[k])
    u = three_phase(lnf_all, p, Float64(case.T_guess), [0.999, 1e-9, 0.99])
    T = Float64(u[4]); x1 = Float64.(u[1:3])
    beta = invariant_amounts(x1, z, 3, Float64(case.vapor_fraction))
    xs = [[x1[k], 1 - x1[k]] for k in 1:3]
    out = [Dict("fraction" => beta[k], "composition" => xs[k],
                "state" => phase_properties(model, p, T, xs[k], phases[k])) for k in 1:3]
    for k in 1:3
        @assert isapprox(volume(model, p, T, xs[k]), out[k]["state"]["molar_volume"]; rtol=1e-9) "invariant phase $k not on the minimum-Gibbs root"
    end
    sort!(out; by=ph -> ph["state"]["mass_density"])
    h = sum(ph["fraction"] * ph["state"]["h"] for ph in out)
    println("PR invariant ", case.id, ": T=", T)
    Dict("id" => String(case.id), "variant" => variant, "ids" => ids, "z" => z, "kij" => kij_matrix(ids),
         "p" => p, "vapor_fraction" => Float64(case.vapor_fraction), "T" => T, "target_h" => h, "phases" => out)
end

function nrtl_invariant(case)
    z, p = Float64.(case.z), Float64.(case.p)
    kinds = [:liquid, :liquid, :vapour]
    lnf_all(k, T, x) = g_lnf(T, BigFloat(p), x, kinds[k])
    u = three_phase(lnf_all, p, Float64(case.T_guess), [0.85, 0.12, 0.6])
    T = Float64(u[4]); x1 = Float64.(u[1:3])
    beta = invariant_amounts(x1, z, 3, Float64(case.vapor_fraction))
    xs = [[x1[k], 1 - x1[k]] for k in 1:3]
    out = [Dict("fraction" => beta[k], "composition" => xs[k], "kind" => String(kinds[k]),
                "state" => syn_state(T, p, xs[k], kinds[k])) for k in 1:3]
    sort!(out; by=ph -> (ph["kind"] == "liquid", ph["state"]["mass_density"]))
    h = sum(ph["fraction"] * ph["state"]["h"] for ph in out)
    println("NRTL invariant ", case.id, ": T=", T)
    Dict("id" => String(case.id), "z" => z, "p" => p, "vapor_fraction" => Float64(case.vapor_fraction),
         "T" => T, "target_h" => h, "phases" => out)
end

function vlle_main()
    VERSION == v"1.11.7" || error("VLLE references pin Julia 1.11.7")
    pkgversion(Clapeyron) == v"0.6.30" || error("VLLE references pin Clapeyron 0.6.30")
    length(ARGS) == 2 && ARGS[1] == "--output" || error("usage: generate_vlle.jl --output candidate.json")
    cubic = [pr_case(c) for c in VIN.cubic_tp]
    nrtl = [nrtl_case(c) for c in VIN.nrtl_tp]
    cubic_inv = [pr_invariant(c) for c in VIN.cubic_ph_invariant]
    nrtl_inv = [nrtl_invariant(c) for c in VIN.nrtl_ph_invariant]
    paths = [SOURCE, LIQ, INPUT, joinpath(@__DIR__, "generate.jl"), joinpath(@__DIR__, "generate_nrtl.jl")]
    output = Dict("schema" => "vlle-clapeyron-v1", "julia_version" => string(VERSION),
        "clapeyron_version" => string(pkgversion(Clapeyron)),
        "input_sha256" => Dict(replace(relpath(p, ROOT), '\\' => '/') => bytes2hex(sha256(read(p))) for p in paths),
        "generator_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
        "method" => "PR/PR-BM: Clapeyron MultiPhaseTPFlash/MichelsenTPFlash initial phase sets, 256-bit Newton on Clapeyron ln fugacities and balances, composition-grid TPD acceptance (minimum-Gibbs root). NRTL: Clapeyron NRTL activity coefficients and PolExpSat ancillaries, ideal-gas vapor, 256-bit Newton from grid TPD minima with grid acceptance. Invariant PH: 256-bit Newton for the binary three-phase temperature; amounts from the stated vapor fraction. No Rust calls or Rust outputs.",
        "scope" => "Implementation agreement for matched formulations: explicit water-hydrocarbon kij = 0.5 (approximation) and a synthetic NRTL pair on methanol/water pure data. Not experimental accuracy; the grid test is evidence of stability at its resolution, not a proof.",
        # Chosen before comparison; acceptance requires compositions and fugacity residuals <= 1e-7.
        "tolerances" => Dict("temperature" => [2e-6, 0.0], "fraction" => [2e-8, 0.0], "composition" => [1e-8, 1e-7],
            "molar_volume" => [1e-12, 2e-9], "mass_density" => [1e-7, 2e-9], "h" => [2e-5, 2e-9],
            "ln_phi" => [2e-9, 2e-9], "fugacity_residual" => [1e-7, 0.0],
            # Remaining PR phase properties (live-reproduction budgets of clapeyron-v1.json).
            "Z" => [2e-10, 2e-9], "s" => [2e-7, 2e-9], "h_res" => [2e-5, 2e-9], "s_res" => [2e-7, 2e-9],
            "cp" => [2e-5, 2e-8], "cv" => [2e-5, 2e-8], "sound" => [2e-4, 2e-8]),
        "cubic_tp" => cubic, "nrtl_tp" => nrtl, "cubic_ph_invariant" => cubic_inv, "nrtl_ph_invariant" => nrtl_inv)
    @assert finite_tree(output)
    open(ARGS[2], "w") do io
        JSON3.pretty(io, output)
        println(io)
    end
    println("Wrote $(length(cubic)) cubic, $(length(nrtl)) NRTL TP and $(length(cubic_inv) + length(nrtl_inv)) invariant PH references to ", ARGS[2])
end

if abspath(PROGRAM_FILE) == @__FILE__
    vlle_main()
end
