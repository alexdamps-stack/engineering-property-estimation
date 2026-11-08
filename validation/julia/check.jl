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

function live(generator_name, candidate, test, variable; minimal_compilation)
    generator = joinpath(@__DIR__, generator_name)
    project = @__DIR__
    # The PR/PR-BM generator avoids minutes of specializing reference-only flash
    # calls. The NRTL generator's BigFloat solves are far slower interpreted
    # (over 15 min) than compiled (about 1 min), so it runs compiled. Algorithms
    # and AD remain Clapeyron's own implementation either way.
    flags = minimal_compilation ? `--compile=min` : ``
    run(`$(Base.julia_cmd()) $flags --project=$project $generator --output $candidate`)
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
        live("generate.jl", joinpath(dir,"clapeyron-current.json"), "clapeyron", "THERMO_CLAPEYRON_REFERENCE"; minimal_compilation=true)
        live("generate_nrtl.jl", joinpath(dir,"nrtl-current.json"), "nrtl_clapeyron", "THERMO_NRTL_REFERENCE"; minimal_compilation=false)
        live("generate_a13.jl", joinpath(dir,"a13-current.json"), "a13_clapeyron", "THERMO_A13_REFERENCE"; minimal_compilation=false)
    end
    println("Live Clapeyron standard check passed (PR/PR-BM, NRTL and A13 saturation); committed reference fixtures unchanged.")
end

main()
