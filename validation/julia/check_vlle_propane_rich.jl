# Independent reference for the propane-rich tray regression in issue #240.
# Run with: julia --project=validation/julia validation/julia/check_vlle_propane_rich.jl
include("generate_vlle.jl")

VERSION == v"1.11.7" || error("VLLE references pin Julia 1.11.7")
pkgversion(Clapeyron) == v"0.6.30" || error("VLLE references pin Clapeyron 0.6.30")

case = (id = "pr-propane-rich-three-phase-20bar", ids = ["propane", "water"],
        variant = "PR", z = [0.9957, 0.0043], p = 2e6, T_guess = 330.0,
        vapor_fraction = 0.43)
result = redirect_stdout(stderr) do
    pr_invariant(case)
end
JSON3.pretty(stdout, result)
println()
