# Independent Clapeyron.jl references. No Rust calls and no Rust output inputs.
# The existing fixtures supply case inputs/enthalpy targets, not Julia solutions.
# NASA7 low-range Cp/R coefficients are represented exactly by ReidIdeal over
# the tested 200..600 K range; Clapeyron integrates and differentiates them.
using Clapeyron, JSON3, SHA

const ROOT = normpath(joinpath(@__DIR__, "../.."))
const R = 8.31446261815324
const SOURCE = joinpath(ROOT, "data/source-extracts/components.provisional.json")
const DATA = Dict(String(c.id) => c for c in JSON3.read(read(SOURCE, String)).components)
const OA = 0.4572355289213822
const OB = 0.07779607390388846

function matched_model(ids, variant, kij)
    records = [DATA[String(id)] for id in ids]
    coeffs = [Float64.(c.ideal_gas.coefficients[1][1:5]) .* R for c in records]
    constants = [c.constants for c in records]
    mw = [Float64(c.molar_mass_kg_per_mol)*1000 for c in constants]
    ideal = ReidIdeal(String.(ids); userlocations=(;
        a=getindex.(coeffs, 1), b=getindex.(coeffs, 2), c=getindex.(coeffs, 3),
        d=getindex.(coeffs, 4), e=getindex.(coeffs, 5), Mw=mw))
    params = (; Tc=[Float64(c.critical_temperature_K) for c in constants],
                Pc=[Float64(c.critical_pressure_Pa) for c in constants],
                Mw=mw, acentricfactor=[Float64(c.acentric_factor) for c in constants],
                k=hcat([Float64.(row) for row in kij]...), l=zeros(length(ids), length(ids)))
    alpha = variant == "PR" ? PRAlpha : variant == "PR-BM" ? BMAlpha : error("unknown variant $variant")
    model = PR(String.(ids); idealmodel=ideal, alpha, translation=NoTranslation,
        mixing=vdW1fRule, userlocations=params,
        reference_state=ReferenceState(:ideal_gas; T0=298.15, P0=100000.0, H0=0.0, S0=0.0))
    @assert Clapeyron.Rgas(model) == R
    # Refuse hidden formulation changes rather than widening comparison budgets.
    for i in eachindex(ids)
        tc, pc = params.Tc[i], params.Pc[i]
        @assert isapprox(model.params.a.values[i,i], OA*(R*tc)^2/pc; rtol=5e-14)
        @assert isapprox(model.params.b.values[i,i], OB*R*tc/pc; rtol=5e-14)
    end
    for pure in split_model(Clapeyron.idealmodel(model))
        @assert abs(enthalpy(pure, 100000.0, 298.15)) < 1e-8
        @assert abs(entropy(pure, 100000.0, 298.15)) < 1e-8
    end
    return model
end

function phase_properties(model, p, t, x, phase)
    v = volume(model, p, t, x; phase)
    s = entropy(model, p, t, x; phase)
    # Clapeyron's entropy_res is at fixed T,V; our departure is at fixed T,p.
    # Subtract the independently evaluated ideal entropy at the same T,p,x.
    s_res = s - entropy(Clapeyron.idealmodel(model), p, t, x)
    properties = Dict("Z"=>p*v/(R*t), "molar_volume"=>v,
        "mass_density"=>mass_density(model, p, t, x; phase),
        "h"=>enthalpy(model, p, t, x; phase), "s"=>s,
        "h_res"=>enthalpy_res(model, p, t, x; phase), "s_res"=>s_res,
        "cp"=>isobaric_heat_capacity(model, p, t, x; phase),
        "cv"=>isochoric_heat_capacity(model, p, t, x; phase),
        "sound"=>speed_of_sound(model, p, t, x; phase),
        "ln_phi"=>log.(fugacity_coefficient(model, p, t, x; phase)))
    return Dict(key => (value isa AbstractVector ? Float64.(value) : Float64(value)) for (key,value) in properties)
end

function ph_case(case)
    ids, z, variant = String.(case.ids), Float64.(case.z), String(case.variant)
    model = matched_model(ids, variant, case.kij)
    p, h = Float64(case.p), Float64(case.target_h)
    # No temperature or phase-composition solution is supplied as an initial guess.
    flash = ph_flash(model, p, h, z; atol=1e-12, rtol=1e-14, max_iters=200)
    t = flash.data.T
    phases = Any[]
    for i in eachindex(flash.fractions)
        x = flash.compositions[i]
        vl = volume(model, p, t, x; phase=:liquid)
        vv = volume(model, p, t, x; phase=:vapour)
        # At a pure saturation endpoint a fresh minimum-Gibbs root request is
        # ambiguous. Honor the actual volume returned by Clapeyron's PH solve.
        phase = abs(flash.volumes[i]-vl) <= abs(flash.volumes[i]-vv) ? :liquid : :vapour
        push!(phases, Dict("fraction"=>flash.fractions[i], "composition"=>x,
                          "state"=>phase_properties(model, p, t, x, phase)))
    end
    # Normalize ordering by actual molar volume, vapor first for two phases.
    sort!(phases; by=ph->ph["state"]["molar_volume"], rev=true)
    actual_h = sum(ph["fraction"]*ph["state"]["h"] for ph in phases)
    @assert isfinite(t) && abs(actual_h-h)/max(abs(h),1.0) <= 1e-8 "$(case.variant) $(case.id): T=$t h=$actual_h target=$h"
    Dict("id"=>String(case.id), "ids"=>ids, "z"=>z, "variant"=>variant,
         "kij"=>case.kij, "p"=>p, "target_h"=>h, "T"=>t, "h"=>actual_h, "phases"=>phases)
end

function homogeneous_case(case)
    ids, x, variant = String.(case.ids), Float64.(case.x), String(case.variant)
    model = matched_model(ids, variant, case.kij)
    t, p = Float64(case.T), Float64(case.p)
    roots = Any[]
    for (i, root) in enumerate(case.roots)
        root.mechanically_stable || continue
        selection = i == 1 ? "smallest-volume" : "largest-volume"
        phase = i == 1 ? :liquid : :vapour
        # Clapeyron 0.6.30's Float64 cubic volume loses precision in the
        # 0.01 Pa diagnostic limit. Evaluate the *same package* at 256 bits;
        # do not relax the regular-state budgets or replace its EOS/roots.
        state = if p < 1.0
            setprecision(BigFloat,256) do
                phase_properties(model, BigFloat(p), BigFloat(t), BigFloat.(x), phase)
            end
        else
            phase_properties(model, p, t, x, phase)
        end
        push!(roots, Dict("selection"=>selection, "state"=>state))
    end
    Dict("ids"=>ids, "x"=>x, "variant"=>variant, "T"=>t, "p"=>p, "kij"=>case.kij,
         "evaluation_precision_bits"=>(p < 1.0 ? 256 : 53), "roots"=>roots)
end

function finite_tree(x)
    x isa AbstractFloat && return isfinite(x)
    x isa AbstractDict && return all(finite_tree, values(x))
    x isa AbstractVector && return all(finite_tree, x)
    return true
end

function main()
    length(ARGS) == 2 && ARGS[1] == "--output" || error("usage: julia --project=validation/julia validation/julia/generate.jl --output candidate.json")
    @assert pkgversion(Clapeyron) == v"0.6.30"
    ph_input = joinpath(ROOT,"validation/fixtures/ph-v1.json")
    homogeneous_input = joinpath(ROOT,"validation/fixtures/cubic-v1.json")
    ph = Any[]
    for case in JSON3.read(read(ph_input,String)).cases
        println("Clapeyron PH: ", case.variant, " ", case.id)
        push!(ph, ph_case(case))
    end
    homogeneous = homogeneous_case.(JSON3.read(read(homogeneous_input,String)).cases)
    output = Dict("schema"=>"clapeyron-reference-v1", "julia_version"=>string(VERSION),
        "clapeyron_version"=>string(pkgversion(Clapeyron)),
        "clapeyron_commit"=>"3440a15ba224664cac5b24441ac8adacbaf326f9",
        "input_sha256"=>Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in [SOURCE,ph_input,homogeneous_input]),
        "generator_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))),
        "formulation"=>"Exact PR factors; PRAlpha or BMAlpha; vdW1fRule; explicit symmetric kij and l=0; NoTranslation; R=8.31446261815324; NASA7 low Cp coefficients via ReidIdeal; pure ideal h=s=0 at 298.15 K, 1 bar",
        "scope"=>"Implementation agreement for matched parameters. Native Clapeyron PH solve without supplied T/compositions. Shared input parameters, independent derivatives and solvers. No experimental accuracy or global stability proof.",
        # Chosen before comparisons; reuse A3/A5 budgets for the same regular cases.
        "tolerances"=>Dict("temperature"=>[2e-6,0.0], "fraction"=>[2e-8,0.0], "composition"=>[2e-8,1e-7],
            "Z"=>[2e-10,2e-9], "mass_density"=>[1e-7,2e-9], "molar_volume"=>[1e-12,2e-9],
            "h"=>[2e-5,2e-9], "s"=>[2e-7,2e-9], "h_res"=>[2e-5,2e-9], "s_res"=>[2e-7,2e-9],
            "cp"=>[2e-5,2e-8], "cv"=>[2e-5,2e-8], "sound"=>[2e-4,2e-8], "ln_phi"=>[2e-9,2e-9]),
        "ph"=>ph, "homogeneous"=>homogeneous)
    @assert length(ph) == 28 && length(homogeneous) == 60
    @assert finite_tree(output)
    open(ARGS[2],"w") do io
        JSON3.pretty(io, output)
        println(io)
    end
    println("Wrote 28 PH and 60 homogeneous Clapeyron states to ", ARGS[2])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
