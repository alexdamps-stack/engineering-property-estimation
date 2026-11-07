# Independent numerical validation

This directory contains A2 preparation tooling and [pending tolerance proposals](tolerances.provisional.json). No independently generated thermodynamic expected outputs or numerical model validation results are supplied yet.

The [testing strategy](../docs/testing-strategy.md) defines the reference-case requirements, independence rules, numerical scenarios, tolerances, and reporting conventions. Keep small qualified fixtures in Git for offline crate integration tests. Reference-generation tools may use external software; ordinary Rust tests must not require it. Record and review every regeneration of expected values.

Record the reference implementation and version, model formulation, composition, parameters, units, reference states, domain, and per-property absolute/relative tolerances. Keep comparison code separate from the production implementation to avoid comparing a model with itself.

Report implementation verification, model-to-model comparisons, and experimental validation separately. Historical prototype test counts are not current repository results. The initial scaffold CI only checks buildability, formatting, linting, and documentation.

Run `powershell -NoProfile -ExecutionPolicy Bypass -File validation/check-preparation.ps1` from the repository root for offline artifact-integrity checks. The explicitly invoked `prepare-source-data.ps1 -OutputDirectory <directory>` downloads pinned sources and transcribes selected records; it does not implement or evaluate any model. Review regeneration before adopting changes. Exact commands, hashes, checks and limitations are in [preparation checks](preparation-checks.md) and [reference qualification](../docs/reference-qualification.md). These files are provisional and do not finalize A1 integration or A2 acceptance.
