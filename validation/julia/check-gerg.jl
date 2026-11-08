# Native Clapeyron GERG-2008 EOS/AD equilibrium-residual comparisons on the
# independently solved physical exterior branches. No Rust calls or values.
using Clapeyron, JSON3, SHA

VERSION == v"1.11.7" || error("GERG check pins Julia 1.11.7")
pkgversion(Clapeyron) == v"0.6.30" || error("GERG check pins Clapeyron 0.6.30")
const ROOT = normpath(joinpath(@__DIR__,"../.."))
const INPUT = joinpath(ROOT,"validation/fixtures/gerg-equilibrium-v1.json")
const FIXTURE = JSON3.read(read(INPUT,String))
const ALIASES = Dict("carbon-dioxide"=>"carbon dioxide")
const RGERG = 8.314472
const NIST = JSON3.read(read(joinpath(ROOT,"data/gerg-2008.json"),String))
const OFFSETS = Dict{String,Tuple{Float64,Float64}}()

function nist_ideal_anchor(name)
    i = findfirst(==(name),String.(NIST.components))+1
    n,theta = NIST.tables.n0i[i],NIST.tables.th0i[i]
    t = 298.15
    a = log(101.325/(RGERG*t))+n[2]+n[3]/t-n[4]*log(t)
    at = n[4]+n[3]/t
    for j in 4:7
        theta[j+1] == 0 && continue
        v,c = theta[j+1]/t,n[j+1]
        if j == 4 || j == 6
            a += c*log(sinh(v));at += c*v/tanh(v)
        else
            a -= c*log(cosh(v));at -= c*v*tanh(v)
        end
    end
    return RGERG*t*(1+at),RGERG*(at-a)
end

offset(ids,x,index) = sum(xi*OFFSETS[String(id)][index] for (id,xi) in zip(ids,x))

function check(label, actual, expected, key)
    atol, rtol = Float64.(FIXTURE.tolerances[Symbol(key)])
    abs(actual-expected) <= atol+rtol*abs(expected) || error("$label $key: $actual != $expected")
end

function model(case)
    ids = [get(ALIASES,String(id),String(id)) for id in case.ids]
    m = GERG2008(ids)
    Clapeyron.Rgas(m) == RGERG || error("unmatched GERG gas constant")
    for (original,id) in zip(String.(case.ids),ids)
        haskey(OFFSETS,original) && continue
        ideal = Clapeyron.idealmodel(GERG2008([id]))
        h,s = nist_ideal_anchor(original)
        # Reference-only constant shifts. No EOS, cp, fugacity or density
        # parameter changes: Clapeyron ideal h0 and s0 differ from rounded NIST.
        # Evaluate the ideal anchor at the identical native GERG molar density,
        # avoiding any separate ideal-model PT volume convention.
        v0=RGERG*298.15/101325.
        OFFSETS[original] = (h-Clapeyron.VT_enthalpy(ideal,v0,298.15,[1.]),s-Clapeyron.VT_entropy(ideal,v0,298.15,[1.]))
    end
    return m
end

function state(m,p,t,x,v,ids)
    actual_p = Clapeyron.VT_pressure(m,v,t,x)
    abs(actual_p-p) <= 1e-4+2e-10*p || error("Clapeyron GERG pressure mismatch: $actual_p != $p")
    Dict("rho"=>1/v,"Z"=>actual_p*v/(RGERG*t),
        "h"=>Clapeyron.VT_enthalpy(m,v,t,x)+offset(ids,x,1),"s"=>Clapeyron.VT_entropy(m,v,t,x)+offset(ids,x,2),
        "cp"=>Clapeyron.VT_isobaric_heat_capacity(m,v,t,x))
end

function verify_equilibrium(c,t,p,compositions,densities,fractions,target_h)
    active = findall(>(0),Float64.(c.z))
    ids = String.(c.ids)[active]
    m = model((;ids))
    states = []
    potentials = []
    for (composition,rho) in zip(compositions,densities)
        x = Float64.(composition)[active]
        v = 1/rho
        push!(states,state(m,p,t,x,v,ids))
        push!(potentials,Clapeyron.VT_chemical_potential(m,v,t,x))
    end
    if length(potentials) == 2
        residual = maximum(abs.(potentials[1].-potentials[2]))/(RGERG*t)
        isfinite(residual) && residual <= 1e-8 || error("$(c.id): Clapeyron chemical-potential residual $residual")
        balance = sum(fractions[i].*Float64.(compositions[i]) for i in eachindex(fractions))
        maximum(abs.(balance.-Float64.(c.z))) <= 1e-8 || error("$(c.id): component balance mismatch")
    end
    h = sum(fractions[i]*states[i]["h"] for i in eachindex(fractions))
    check(String(c.id),h,target_h,"h")
    return states
end

function main()
    FIXTURE.source_sha256 == bytes2hex(sha256(read(joinpath(ROOT,"data/source-extracts/gerg/GERG2008.cpp")))) || error("GERG source changed")
    FIXTURE.generator_sha256 == bytes2hex(sha256(read(joinpath(ROOT,"validation/generate-gerg-equilibrium-references.py")))) || error("GERG generator changed; review references explicitly")
    FIXTURE.a4_adapter_sha256 == bytes2hex(sha256(read(joinpath(ROOT,"validation/generate-flash-references.py")))) || error("A4 adapter changed")
    for c in FIXTURE.tp
        println("Native Clapeyron GERG TP state/residual: ",c.id);flush(stdout)
        p,t,z = Float64(c.p),Float64(c.T),Float64.(c.z)
        compositions = c.phases == 1 ? [c.z] : [c.vapor,c.liquid]
        densities = c.phases == 1 ? [Float64(c.rho)] : [Float64(c.vapor_properties.rho),Float64(c.liquid_properties.rho)]
        fractions = c.phases == 1 ? [1.] : [Float64(c.vapor_fraction),1-Float64(c.vapor_fraction)]
        states = verify_equilibrium(c,t,p,compositions,densities,fractions,Float64(c.enthalpy))
        for (i,ph) in enumerate(states)
            expected = c.phases == 1 ? c : i == 1 ? c.vapor_properties : c.liquid_properties
            for key in ["Z","h","s","cp"]
                check(String(c.id),ph[key],Float64(expected[Symbol(key)]),key)
            end
        end
        if String(c.id) in ("binary-split","binary-boundary-high")
            # Additional fresh native flash benchmarks, with Wilson guesses only.
            # Other cases are covered by NIST/teqp/SciPy solved references above;
            # nonfinite results from the package's general GERG flash initialization
            # are not reclassified as valid single phases or accepted as evidence.
            m = model(c)
            flash = Clapeyron.tp_flash2(m,p,t,z,MichelsenTPFlash(K_tol=1e-12,ss_iters=2000))
            order = sortperm(flash.volumes;rev=true)
            length(order) == 2 || error("$(c.id): native binary flash failed")
            check(String(c.id),flash.fractions[order[1]],Float64(c.vapor_fraction),"vapor_fraction")
            for (i,j) in enumerate(order), (actual,want) in zip(flash.compositions[j],compositions[i])
                check(String(c.id),actual,Float64(want),"composition")
            end
        end
    end
    for c in FIXTURE.ph
        println("Native Clapeyron GERG PH state/residual: ",c.id);flush(stdout)
        p,t,h = Float64(c.p),Float64(c.T),Float64(c.target_h)
        if haskey(c,:saturation)
            sat = c.saturation
            compositions = [c.z,c.z]
            densities = [Float64(sat.vapor_properties.rho),Float64(sat.liquid_properties.rho)]
            fractions = [Float64(c.vapor_fraction),1-Float64(c.vapor_fraction)]
        else
            expected = c.expected_tp
            compositions = expected.phases == 1 ? [expected.z] : [expected.vapor,expected.liquid]
            densities = expected.phases == 1 ? [Float64(expected.rho)] : [Float64(expected.vapor_properties.rho),Float64(expected.liquid_properties.rho)]
            fractions = expected.phases == 1 ? [1.] : [Float64(expected.vapor_fraction),1-Float64(expected.vapor_fraction)]
        end
        verify_equilibrium(c,t,p,compositions,densities,fractions,h)
    end
    println("Native Clapeyron GERG EOS/AD check passed: $(length(FIXTURE.tp)) TP and $(length(FIXTURE.ph)) PH equilibrium states plus two fresh binary TP solves; no reference files modified.")
end

main()
