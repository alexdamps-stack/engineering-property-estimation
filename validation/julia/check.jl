# Full live standard check. Instantiates the locked Julia environment, validates
# the committed fixtures' provenance, generates temporary candidates and tests
# the current Rust engine against them. Never replaces committed expected values.
# Covers PR/PR-BM (clapeyron-v1.json) and the A7 NRTL/gamma-phi path (nrtl-v1.json).
using Pkg, SHA

VERSION == v"1.11.7" || error("Clapeyron standard check pins Julia 1.11.7; found $VERSION")
Pkg.instantiate()
using JSON3

const ROOT = normpath(joinpath(@__DIR__, "../.."))

function verify_provenance(fixture_name, generator_name)
    fixture = JSON3.read(read(joinpath(ROOT, "validation/fixtures", fixture_name), String))
    for (path, expected) in pairs(fixture.input_sha256)
        actual = bytes2hex(sha256(read(joinpath(ROOT, String(path)))))
        actual == expected || error("$fixture_name input changed: $path; regenerate and review references explicitly")
    end
    fixture.generator_sha256 == bytes2hex(sha256(read(joinpath(@__DIR__, generator_name)))) || error("$generator_name changed; regenerate and review $fixture_name explicitly")
    fixture.manifest_sha256 == bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))) || error("Julia environment changed; regenerate and review $fixture_name explicitly")
end

# Number leaves use the committed fixture's own declared comparison tolerances,
# |fresh - committed| <= atol + rtol*|committed|. Keys that are not named in the
# fixture's "tolerances" block resolve through these aliases; any other number
# (input parameters such as kij) must reproduce exactly.
const TOLERANCE_ALIASES = Dict(
    "T"=>"temperature", "p"=>"pressure", "dew_pressure"=>"pressure", "bubble_pressure"=>"pressure",
    "x"=>"composition", "y"=>"composition", "z"=>"composition", "liquid"=>"composition", "vapor"=>"composition",
    "vl"=>"volume", "vv"=>"volume", "vapor_fraction"=>"fraction",
    "target_h"=>"h", "h_liquid"=>"h", "h_vapor"=>"h", "h_ideal"=>"h", "s_liquid"=>"s", "cp_liquid"=>"cp")
# Solver diagnostics recorded as bounds rather than reference values; the
# generators assert these bounds themselves (generate_a13.jl: residual < 1e-8).
const DIAGNOSTIC_BOUNDS = Dict("fugacity_residual"=>1e-8)

function leaf_tolerance(tolerances, key)
    haskey(DIAGNOSTIC_BOUNDS, key) && return (DIAGNOSTIC_BOUNDS[key], 0.0)
    name = haskey(tolerances, Symbol(key)) ? key : get(TOLERANCE_ALIASES, key, nothing)
    (name === nothing || !haskey(tolerances, Symbol(name))) && return (0.0, 0.0)
    atol, rtol = tolerances[Symbol(name)]
    return (Float64(atol), Float64(rtol))
end

# Recursively compares a fresh reference with the committed fixture. Structure,
# strings, provenance hashes and exact inputs must be identical; numbers must
# agree within the tolerance declared for their quantity.
function compare_reference!(failures, fresh, committed, tolerances, path, key)
    if committed isa AbstractDict
        fresh isa AbstractDict || return push!(failures, "$path: type differs")
        Set(keys(fresh)) == Set(keys(committed)) || return push!(failures, "$path: keys differ")
        for k in keys(committed)
            compare_reference!(failures, fresh[k], committed[k], tolerances, "$path.$k", String(k))
        end
    elseif committed isa AbstractVector
        (fresh isa AbstractVector && length(fresh) == length(committed)) || return push!(failures, "$path: length differs")
        for i in eachindex(committed)
            compare_reference!(failures, fresh[i], committed[i], tolerances, "$path[$i]", key)
        end
    elseif committed isa Real && !(committed isa Bool)
        (fresh isa Real && !(fresh isa Bool)) || return push!(failures, "$path: type differs")
        atol, rtol = leaf_tolerance(tolerances, key)
        abs(Float64(fresh) - Float64(committed)) <= atol + rtol * abs(Float64(committed)) ||
            push!(failures, "$path: fresh $fresh vs committed $committed (atol $atol, rtol $rtol)")
    else
        isequal(fresh, committed) || push!(failures, "$path: fresh $(repr(fresh)) vs committed $(repr(committed))")
    end
    return failures
end

function assert_matches_committed(candidate, fixture_name)
    committed = JSON3.read(read(joinpath(ROOT, "validation/fixtures", fixture_name), String))
    fresh = JSON3.read(read(candidate, String))
    failures = compare_reference!(String[], fresh, committed, committed.tolerances, "\$", "")
    if !isempty(failures)
        foreach(f -> println(stderr, "  ", f), first(failures, 50))
        error("Fresh Clapeyron reference differs from committed $fixture_name in $(length(failures)) value(s); review explicitly and never replace the fixture automatically")
    end
    println("Fresh Clapeyron reference reproduces committed $fixture_name within its declared tolerances.")
end

function live(generator_name, fixture_name, candidate, test, variable; minimal_compilation)
    generator = joinpath(@__DIR__, generator_name)
    project = @__DIR__
    # The PR/PR-BM generator avoids minutes of specializing reference-only flash
    # calls. The NRTL generator's BigFloat solves are far slower interpreted
    # (over 15 min) than compiled (about 1 min), so it runs compiled. Algorithms
    # and AD remain Clapeyron's own implementation either way.
    flags = minimal_compilation ? `--compile=min` : ``
    run(`$(Base.julia_cmd()) $flags --project=$project $generator --output $candidate`)
    assert_matches_committed(candidate, fixture_name)
    cd(ROOT) do
        # Environment inheritance preserves caller CARGO_TARGET_DIR/toolchain settings.
        run(addenv(`cargo test -p thermo-solvers --test $test --locked -- --nocapture`, variable=>candidate))
    end
end

function main()
    verify_provenance("clapeyron-v1.json", "generate.jl")
    verify_provenance("nrtl-v1.json", "generate_nrtl.jl")
    verify_provenance("a13-v1.json", "generate_a13.jl")
    mktempdir() do dir
        live("generate.jl", "clapeyron-v1.json", joinpath(dir,"clapeyron-current.json"), "clapeyron", "THERMO_CLAPEYRON_REFERENCE"; minimal_compilation=true)
        live("generate_nrtl.jl", "nrtl-v1.json", joinpath(dir,"nrtl-current.json"), "nrtl_clapeyron", "THERMO_NRTL_REFERENCE"; minimal_compilation=false)
        live("generate_a13.jl", "a13-v1.json", joinpath(dir,"a13-current.json"), "a13_clapeyron", "THERMO_A13_REFERENCE"; minimal_compilation=false)
    end
    println("Live Clapeyron standard check passed (PR/PR-BM, NRTL and A13 saturation); fresh references reproduce the committed fixtures, which were not modified.")
end

main()
