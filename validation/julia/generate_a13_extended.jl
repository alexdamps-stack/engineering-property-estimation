# Independent #227 extended bubble/dew references. No Rust calls and no Rust outputs.
#
# Model construction reuses the qualified adapter of generate.jl (exact PR factors,
# PRAlpha/BMAlpha, vdW1fRule, NoTranslation, explicit kij and l = 0). The
# component records of the reacting and ChemSep extracts are added to its lookup;
# pseudo-components use the critical constants of the independent Python
# characterisation fixture (validation/fixtures/pseudo-components-v1.json) with
# Clapeyron's BasicIdeal (the ideal part does not enter bubble/dew equations).
#
# Two-phase boundaries: Clapeyron's native bubble_pressure, dew_pressure,
# bubble_temperature and dew_temperature. For a water-bearing dew point the native
# solve is repeated from aqueous and hydrocarbon incipient-liquid guesses (no
# temperature or pressure guess). Every candidate must satisfy fugacity equality,
# a vapor/liquid volume ratio > 1.5, the stated search bounds, and a tangent-plane
# test of both phases (a composition grid for up to three components, Clapeyron's
# own `tpd` for more); exactly one distinct accepted candidate is required.
#
# Three-phase (heterogeneous) bubble points: Clapeyron's MultiPhaseTPFlash at the
# stated liquid-side guess supplies the two liquids; the 3C+1 equations (two
# liquids with Σ n = z, equal ln fugacities of all three phases, Σ y = 1, and the
# unknown ln T or ln p) are solved by Newton's method at 256 bits on Clapeyron
# fugacity coefficients, and accepted only if the composition grid finds no
# negative TPD from the liquid tangent plane.
include("generate.jl")

const EXTRA_SOURCES = [joinpath(ROOT, "data/source-extracts/components-reacting.provisional.json"),
                       joinpath(ROOT, "data/source-extracts/components-chemsep.provisional.json")]
for path in EXTRA_SOURCES, c in JSON3.read(read(path, String)).components
    haskey(DATA, String(c.id)) && error("duplicate component $(c.id)")
    DATA[String(c.id)] = c
end
const PSEUDO_SOURCE = joinpath(ROOT, "validation/fixtures/pseudo-components-v1.json")
const PSEUDO = Dict(String(s.id) => s.correlations for s in JSON3.read(read(PSEUDO_SOURCE, String)).species)
const INPUT = joinpath(ROOT, "validation/a13-extended-inputs.json")

# PR/PR-BM with pseudo-component constants; same formulation checks as matched_model.
function pseudo_model(ids, variant, kij)
    tc = Float64[]; pc = Float64[]; w = Float64[]; mw = Float64[]
    for id in ids
        c = startswith(id, "pseudo:") ? PSEUDO[id] : DATA[id].constants
        push!(tc, Float64(c.critical_temperature_K)); push!(pc, Float64(c.critical_pressure_Pa))
        push!(w, Float64(c.acentric_factor)); push!(mw, Float64(c.molar_mass_kg_per_mol) * 1000)
    end
    params = (; Tc=tc, Pc=pc, Mw=mw, acentricfactor=w,
                k=hcat([Float64.(row) for row in kij]...), l=zeros(length(ids), length(ids)))
    alpha = variant == "PR" ? PRAlpha : variant == "PR-BM" ? BMAlpha : error("unknown variant $variant")
    model = PR(String.(ids); idealmodel=BasicIdeal, alpha, translation=NoTranslation,
        mixing=vdW1fRule, userlocations=params)
    @assert Clapeyron.Rgas(model) == R
    for i in eachindex(ids)
        @assert isapprox(model.params.a.values[i,i], OA*(R*tc[i])^2/pc[i]; rtol=5e-14)
        @assert isapprox(model.params.b.values[i,i], OB*R*tc[i]/pc[i]; rtol=5e-14)
    end
    model
end

model_for(c) = any(id -> startswith(String(id), "pseudo:"), c.ids) ?
    pseudo_model(String.(c.ids), String(c.variant), c.kij) :
    matched_model(String.(c.ids), String(c.variant), c.kij)

# ---------------------------------------------------------------------------
# Tangent-plane tests.
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
lnf(model, p, T, x, phase) = log.(max.(x, 1e-300)) .+ log.(fugacity_coefficient(model, p, T, x; phase))
# Minimum modified TPD over the grid (minimum-Gibbs root at each trial) from the
# tangent plane d of a reference phase.
function grid_min_tpd(model, p, T, d)
    best = Inf
    for w in simplex_grid(length(d))
        lw = log.(max.(w, 1e-300)) .+ log.(fugacity_coefficient(model, p, T, w))
        best = min(best, sum(w[i] > 0 ? w[i] * (lw[i] - d[i]) : 0.0 for i in eachindex(w)))
    end
    best
end
# True when no trial phase lies below the tangent plane of composition x on `phase`.
function tangent_plane_stable(model, p, T, x, phase)
    count(>(0), x) < 2 && return true
    # The selected root must be the minimum-Gibbs root of its composition.
    g(ph) = sum(x[i] > 0 ? x[i] * log(fugacity_coefficient(model, p, T, x; phase=ph)[i]) : 0.0 for i in eachindex(x))
    g(phase) <= min(g(:liquid), g(:vapour)) + 1e-10 || return false
    if length(x) <= 3
        return grid_min_tpd(model, p, T, lnf(model, p, T, x, phase)) >= -1e-9
    end
    _, tpds = Clapeyron.tpd(model, p, T, x)
    return isempty(tpds) || minimum(tpds) >= -1e-9
end

# ---------------------------------------------------------------------------
# Two-phase boundaries.
function native(m, c, fixed, z, x0)
    bubble = c.kind == "bubble"
    if c.unknown == "pressure"
        bubble ? bubble_pressure(m, fixed, z) :
            x0 === nothing ? dew_pressure(m, fixed, z) : dew_pressure(m, fixed, z, ChemPotDewPressure(x0=x0))
    else
        bubble ? bubble_temperature(m, fixed, z) :
            x0 === nothing ? dew_temperature(m, fixed, z) : dew_temperature(m, fixed, z, ChemPotDewTemperature(x0=x0))
    end
end

function incipient_guesses(c, z)
    ids = String.(c.ids)
    w = findfirst(==("water"), ids)
    (c.kind == "dew" && w !== nothing && count(>(0), z) > 1) || return Any[nothing]
    dry = sum(z[i] for i in eachindex(z) if i != w)
    guesses = Any[nothing]
    for xw in (0.999, 0.001)
        push!(guesses, [i == w ? xw : (1 - xw) * z[i] / dry for i in eachindex(z)])
    end
    guesses
end

# Refines a native boundary at 256 bits: Newton on ln x_i + ln φ_i^L − ln y_i − ln φ_i^V = 0
# and Σ w − 1 = 0 for the unknown ln T or ln p and the incipient ln w. Native
# Clapeyron solves stop near 1e-8 in ln f, which is too loose for the A13 budgets.
function polish_two_phase(m, c, fixed, z, u0, w0)
    pressure_unknown = c.unknown == "pressure"; bubble = c.kind == "bubble"
    act = [i for i in eachindex(z) if z[i] > 0]
    setprecision(BigFloat, 256) do
        zb = BigFloat.(z); fb = BigFloat(fixed)
        function F(q)
            T, p = pressure_unknown ? (fb, exp(q[1])) : (exp(q[1]), fb)
            w = zeros(BigFloat, length(z)); w[act] .= exp.(q[2:end])
            x, y = bubble ? (zb, w) : (w, zb)
            l = log.(fugacity_coefficient(m, p, T, x ./ sum(x); phase=:liquid))
            v = log.(fugacity_coefficient(m, p, T, y ./ sum(y); phase=:vapour))
            vcat([log(x[i]) + l[i] - log(y[i]) - v[i] for i in act], [sum(w) - 1])
        end
        q = vcat(BigFloat(u0), log.(BigFloat.(w0[act])))
        for _ in 1:100
            f = F(q)
            if maximum(abs, f) < big"1e-60"
                T, p = pressure_unknown ? (fixed, Float64(exp(q[1]))) : (Float64(exp(q[1])), fixed)
                w = zeros(length(z)); w[act] .= Float64.(exp.(q[2:end]))
                return T, p, w ./ sum(w)
            end
            h = big"1e-35"
            J = reduce(hcat, [(F(q .+ h .* ((1:length(q)) .== j)) .- F(q .- h .* ((1:length(q)) .== j))) ./ (2h) for j in eachindex(q)])
            q = q .- J \ f
        end
        error("$(c.id): 256-bit refinement did not converge")
    end
end

function two_phase_case(c)
    m = model_for(c)
    z = Float64.(c.z); fixed = Float64(c.fixed)
    tb = Float64.(c.bounds.T); pb = Float64.(c.bounds.p)
    accepted = Any[]
    for x0 in incipient_guesses(c, z)
        r = try native(m, c, fixed, z, x0) catch e
            println("  native solve failed: ", sprint(showerror, e)[1:min(160, end)]); continue
        end
        all(isfinite, r[1]) && all(isfinite, r[4]) || (println("  native solve not finite"); continue)
        t, p, w = try
            count(>(0), z) == 1 ? ((c.unknown == "pressure" ? (fixed, Float64(r[1])) : (Float64(r[1]), fixed))..., z) :
                polish_two_phase(m, c, fixed, z, log(Float64(r[1])), Float64.(r[4]))
        catch e
            println("  refinement failed: ", sprint(showerror, e)[1:min(160, end)]); continue
        end
        x, y = c.kind == "bubble" ? (z, w) : (w, z)
        vl, vv = volume(m, p, t, x; phase=:liquid), volume(m, p, t, y; phase=:vapour)
        ok = all(isfinite, [t, p, vl, vv]) && vl > 0 && vv / vl > 1.5 &&
             tb[1] <= t <= tb[2] && pb[1] <= p <= pb[2]
        ok || (println("  candidate outside bounds or near-critical: T=$t p=$p"); continue)
        l = log.(fugacity_coefficient(m, p, t, x; phase=:liquid))
        v = log.(fugacity_coefficient(m, p, t, y; phase=:vapour))
        act = [i for i in eachindex(z) if x[i] > 0 && y[i] > 0]
        residual = maximum(abs.(log.(x[act]) .+ l[act] .- log.(y[act]) .- v[act]))
        residual < 1e-8 || (println("  candidate fugacity residual $residual"); continue)
        stable = tangent_plane_stable(m, p, t, x, :liquid) && tangent_plane_stable(m, p, t, y, :vapour)
        stable || (println("  candidate rejected by the tangent-plane test: T=$t p=$p"); continue)
        any(a -> abs(a["T"] - t) < 1e-5 && abs(a["p"] / p - 1) < 1e-7, accepted) && continue
        push!(accepted, Dict("id"=>String(c.id), "T"=>t, "p"=>p, "x"=>x, "y"=>y, "vl"=>vl, "vv"=>vv,
                             "fugacity_residual"=>residual))
    end
    length(accepted) == 1 || error("$(c.id): $(length(accepted)) accepted boundaries")
    accepted[1]
end

# ---------------------------------------------------------------------------
# Three-phase bubble points.
function three_phase_case(c)
    m = model_for(c)
    z = Float64.(c.z); C = length(z); fixed = Float64(c.fixed)
    pressure_unknown = c.unknown == "pressure"
    g = Float64(c.guess)
    T0, p0 = pressure_unknown ? (fixed, g) : (g, fixed)
    flash = Clapeyron.tp_flash2(m, p0, T0, z, MultiPhaseTPFlash())
    length(flash.fractions) == 2 || error("$(c.id): guess state is not two-phase: $(flash.fractions)")
    xs = [Float64.(x) for x in flash.compositions]
    n0 = [flash.fractions[k] .* xs[k] for k in 1:2]
    phi1 = fugacity_coefficient(m, p0, T0, xs[1]; phase=:liquid)
    y0 = xs[1] .* phi1 ./ sum(xs[1] .* phi1)
    setprecision(BigFloat, 256) do
        zb = BigFloat.(z); fb = BigFloat(fixed)
        coords(u) = pressure_unknown ? (fb, exp(u)) : (exp(u), fb)
        function F(q)
            T, p = coords(q[1])
            n1 = exp.(q[2:C+1]); n2 = exp.(q[C+2:2C+1]); y = exp.(q[2C+2:3C+1])
            x1 = n1 ./ sum(n1); x2 = n2 ./ sum(n2)
            l1 = lnf(m, p, T, x1, :liquid); l2 = lnf(m, p, T, x2, :liquid)
            v = log.(y) .+ log.(fugacity_coefficient(m, p, T, y ./ sum(y); phase=:vapour))
            vcat(n1 .+ n2 .- zb, l2 .- l1, v .- l1, [sum(y) - 1])
        end
        q = vcat(log(BigFloat(g)), log.(BigFloat.(n0[1])), log.(BigFloat.(n0[2])), log.(BigFloat.(y0)))
        converged = false
        for _ in 1:200
            f = F(q)
            if maximum(abs, f) < big"1e-60"
                converged = true; break
            end
            h = big"1e-35"
            J = reduce(hcat, [(F(q .+ h .* ((1:length(q)) .== j)) .- F(q .- h .* ((1:length(q)) .== j))) ./ (2h) for j in eachindex(q)])
            dq = -(J \ f)
            lambda = big(1.0); accepted = false
            while lambda > 1e-12
                cand = q .+ lambda .* dq
                fc = try F(cand) catch; nothing end
                if fc !== nothing && all(isfinite, fc) && sum(abs2, fc) < (1 - 1e-4 * lambda) * sum(abs2, f)
                    q = cand; accepted = true; break
                end
                lambda /= 2
            end
            accepted || error("$(c.id): three-phase Newton line search failed")
        end
        converged || error("$(c.id): three-phase Newton did not converge")
        T, p = Float64.(coords(q[1]))
        n1 = Float64.(exp.(q[2:C+1])); n2 = Float64.(exp.(q[C+2:2C+1])); y = Float64.(exp.(q[2C+2:3C+1]))
        b1, b2 = sum(n1), sum(n2); x1 = n1 ./ b1; x2 = n2 ./ b2
        tb = Float64.(c.bounds.T); pb = Float64.(c.bounds.p)
        @assert tb[1] <= T <= tb[2] && pb[1] <= p <= pb[2] "$(c.id): outside bounds T=$T p=$p"
        @assert maximum(abs.(x1 .- x2)) > 0.03 "$(c.id): liquids coincide"
        @assert 1e-9 < b1 < 1 - 1e-9 "$(c.id): liquid amount vanishes"
        vl1 = volume(m, p, T, x1; phase=:liquid); vl2 = volume(m, p, T, x2; phase=:liquid)
        vv = volume(m, p, T, y; phase=:vapour)
        @assert vv / max(vl1, vl2) > 1.5 "$(c.id): near-critical vapor"
        for (x, ph) in ((x1, :liquid), (x2, :liquid), (y, :vapour))
            @assert tangent_plane_stable(m, p, T, x, ph) "$(c.id): phase on a metastable root"
        end
        # No phase below the common tangent plane (a fourth phase).
        @assert grid_min_tpd(m, p, T, lnf(m, p, T, x1, :liquid)) >= -1e-9 "$(c.id): grid finds a fourth phase"
        r = maximum(maximum(abs.(lnf(m, p, T, a, pa) .- lnf(m, p, T, x1, :liquid))) for (a, pa) in ((x2, :liquid), (y, :vapour)))
        # Lighter liquid (lower mass density) first, as reported by the engine.
        ρ1 = mass_density(m, p, T, x1; phase=:liquid); ρ2 = mass_density(m, p, T, x2; phase=:liquid)
        light, heavy = ρ1 <= ρ2 ? ((x1, vl1, b1), (x2, vl2, b2)) : ((x2, vl2, b2), (x1, vl1, b1))
        Dict("id"=>String(c.id), "T"=>T, "p"=>p, "x"=>light[1], "vl"=>light[2], "x2"=>heavy[1],
             "vl2"=>heavy[2], "beta"=>heavy[3], "y"=>y, "vv"=>vv, "fugacity_residual"=>r)
    end
end

function a13_extended_main()
    VERSION == v"1.11.7" || error("#227 references pin Julia 1.11.7")
    pkgversion(Clapeyron) == v"0.6.30" || error("#227 references pin Clapeyron 0.6.30")
    length(ARGS) == 2 && ARGS[1] == "--output" || error("usage: generate_a13_extended.jl --output candidate.json")
    cases = Any[]
    for c in JSON3.read(read(INPUT, String)).cases
        println("#227 ", c.id); flush(stdout)
        r = get(c, :three_phase, false) ? three_phase_case(c) : two_phase_case(c)
        println("  T=", r["T"], " p=", r["p"]); flush(stdout)
        push!(cases, r)
    end
    paths = vcat([SOURCE], EXTRA_SOURCES, [PSEUDO_SOURCE, INPUT, joinpath(@__DIR__, "generate.jl")])
    output = Dict("schema"=>"a13-extended-clapeyron-v1", "julia_version"=>string(VERSION),
        "clapeyron_version"=>string(pkgversion(Clapeyron)),
        "input_sha256"=>Dict(replace(relpath(p, ROOT), '\\'=>'/')=>bytes2hex(sha256(read(p))) for p in paths),
        "generator_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml")))),
        "method"=>"Two-phase: native Clapeyron bubble/dew solvers (water-bearing dew points also from aqueous and hydrocarbon incipient-liquid guesses), accepted after fugacity, bound, volume-ratio and tangent-plane tests, exactly one accepted root. Three-phase bubble: Clapeyron MultiPhaseTPFlash liquids at a liquid-side guess, 256-bit Newton on Clapeyron ln fugacities, grid TPD acceptance. No Rust calls or Rust outputs.",
        "scope"=>"Matched exact-factor PR/PR-BM, NoTranslation, explicit kij of each case (synthetic zero, or the explicit water/nonpolar 0.5 approximation). Independent implementation verification, not experimental accuracy or a global-stability proof.",
        # Declared before any comparison with the engine: the unchanged A13 budgets,
        # extended to the second liquid and its amount.
        "tolerances"=>Dict("temperature"=>[2e-6, 0.0], "pressure"=>[0.02, 2e-8], "composition"=>[2e-8, 1e-7],
            "volume"=>[1e-12, 2e-8], "x2"=>[2e-8, 1e-7], "vl2"=>[1e-12, 2e-8], "beta"=>[2e-8, 0.0]),
        "cases"=>cases)
    @assert finite_tree(output)
    open(ARGS[2], "w") do io; JSON3.pretty(io, output); println(io); end
    println("Wrote ", length(cases), " #227 saturation references to ", ARGS[2])
end
if abspath(PROGRAM_FILE) == @__FILE__; a13_extended_main(); end
