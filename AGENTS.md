# Guidance for AI contributors

For PR/PR-BM model, caloric-data, or equilibrium-solver changes, run the
Clapeyron standard check in addition to ordinary Rust checks:
`julia --project=validation/julia validation/julia/check.jl` with Julia 1.11.7.
The locked Julia environment and offline reference tests are described in
`validation/clapeyron-checks.md`. Keep reference generation independent of the
Rust engine. Do not automatically replace committed fixtures or widen tolerances
to resolve a failure; review formulation, data, and reference-state compatibility.

For flowsheet, PFD/P&ID, equipment-symbol, or process-interface design, consult
`docs/design-reference/iso-10628/README.md` and `manifest.json` before choosing
equipment graphics. The library contains all 296 symbols from the linked Commons
reference as editable SVG assets, plus searchable HTML and PNG contact sheets.
Reuse the traced source geometry and retain its attribution/license. Define
connection ports separately; the reference does not imply implemented simulator
capabilities. Follow `docs/flowsheet-and-ui.md` for the current interface proposal.
