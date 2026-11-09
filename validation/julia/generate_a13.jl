# Native Clapeyron saturation solves; no Rust calls, outputs or fitted values.
# Reuse only the already qualified independent model-construction adapter.
include("generate.jl")

function a13_main()
    VERSION == v"1.11.7" || error("A13 pins Julia 1.11.7")
    pkgversion(Clapeyron) == v"0.6.30" || error("A13 pins Clapeyron 0.6.30")
    length(ARGS)==2 && ARGS[1]=="--output" || error("usage: generate_a13.jl --output candidate.json")
    input=joinpath(ROOT,"validation/a13-inputs.json")
    cases=Any[]
    for c in JSON3.read(read(input,String)).cases
        m=matched_model(String.(c.ids),String(c.variant),c.kij)
        z=Float64.(c.z); fixed=Float64(c.fixed)
        pressure_unknown=c.unknown=="pressure"; bubble=c.kind=="bubble"
        f=pressure_unknown ? (bubble ? bubble_pressure : dew_pressure) : (bubble ? bubble_temperature : dew_temperature)
        r=f(m,fixed,z)
        t,p=pressure_unknown ? (fixed,r[1]) : (r[1],fixed)
        x,y=bubble ? (z,Float64.(r[4])) : (Float64.(r[4]),z)
        vl,vv=Float64(r[2]),Float64(r[3])
        @assert all(isfinite,[t,p,vl,vv]) && vl>0 && vv/vl>1.5 "nonphysical/near-critical $(c.id): $r"
        l=log.(fugacity_coefficient(m,p,t,x;phase=:liquid))
        v=log.(fugacity_coefficient(m,p,t,y;phase=:vapour))
        residual=maximum(abs.(log.(x).+l.-log.(y).-v))
        @assert residual<1e-8 "reference fugacity failure $(c.id): $residual"
        push!(cases,Dict("id"=>String(c.id),"T"=>t,"p"=>p,"x"=>x,"y"=>y,"vl"=>vl,"vv"=>vv,"fugacity_residual"=>residual))
        println("A13 ",c.id," T=",t," p=",p," residual=",residual)
        flush(stdout)
    end
    paths=[SOURCE,input,joinpath(@__DIR__,"generate.jl")]
    output=Dict("schema"=>"a13-clapeyron-v1","julia_version"=>string(VERSION),"clapeyron_version"=>string(pkgversion(Clapeyron)),
        "input_sha256"=>Dict(replace(relpath(p,ROOT),'\\'=>'/')=>bytes2hex(sha256(read(p))) for p in paths),
        "generator_sha256"=>bytes2hex(sha256(read(@__FILE__))),"manifest_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))),
        "method"=>"Native Clapeyron bubble_pressure/dew_pressure/bubble_temperature/dew_temperature; no supplied solution, Rust calls or Rust outputs",
        "scope"=>"Matched exact-factor PR/PR-BM, NoTranslation, explicit zero kij. Independent implementation verification, not experimental accuracy or global-stability proof.",
        "tolerances"=>Dict("temperature"=>[2e-6,0.],"pressure"=>[0.02,2e-8],"composition"=>[2e-8,1e-7],"volume"=>[1e-12,2e-8]),"cases"=>cases)
    @assert finite_tree(output)
    open(ARGS[2],"w") do io;JSON3.pretty(io,output);println(io);end
end
if abspath(PROGRAM_FILE)==@__FILE__;a13_main();end
