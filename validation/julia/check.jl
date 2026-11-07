# Full live standard check. Instantiates the locked Julia environment, validates
# the committed fixture's provenance, generates a temporary candidate and tests
# the current Rust engine against it. Never replaces committed expected values.
using Pkg, SHA

VERSION == v"1.11.7" || error("Clapeyron standard check pins Julia 1.11.7; found $VERSION")
Pkg.instantiate()
using JSON3

const ROOT = normpath(joinpath(@__DIR__, "../.."))

function main()
    fixture = JSON3.read(read(joinpath(ROOT, "validation/fixtures/clapeyron-v1.json"), String))
    for (path, expected) in pairs(fixture.input_sha256)
        actual = bytes2hex(sha256(read(joinpath(ROOT, String(path)))))
        actual == expected || error("Clapeyron fixture input changed: $path; regenerate and review references explicitly")
    end
    fixture.generator_sha256 == bytes2hex(sha256(read(joinpath(@__DIR__,"generate.jl")))) || error("Clapeyron generator changed; regenerate and review references explicitly")
    fixture.manifest_sha256 == bytes2hex(sha256(read(joinpath(@__DIR__,"Manifest.toml")))) || error("Julia environment changed; regenerate and review references explicitly")
    mktempdir() do dir
        candidate = joinpath(dir,"clapeyron-current.json")
        generator = joinpath(@__DIR__,"generate.jl")
        project = @__DIR__
        # Avoid spending minutes specializing the reference-only flash call;
        # numerical algorithms and AD remain Clapeyron's own implementation.
        run(`$(Base.julia_cmd()) --compile=min --project=$project $generator --output $candidate`)
        cd(ROOT) do
            # Environment inheritance preserves caller CARGO_TARGET_DIR/toolchain settings.
            run(addenv(`cargo test -p thermo-solvers --test clapeyron --locked -- --nocapture`, "THERMO_CLAPEYRON_REFERENCE"=>candidate))
        end
    end
    println("Live Clapeyron standard check passed; committed reference fixture unchanged.")
end

main()
