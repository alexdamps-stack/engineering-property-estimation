# Contributing

Use English for new documentation, source comments, and commit messages. Keep the German original plan as a historical baseline; update the English development plan when scope or acceptance criteria change.

## Local workflow

1. Install the toolchain selected by `rust-toolchain.toml` and your platform's native linker prerequisites.
2. Make a focused change on a feature branch.
3. Run the formatting, Clippy, test, and documentation commands in the root README.
4. Describe behavior, evidence, and remaining limits in the pull request.

Commit the root `Cargo.lock`. Do not commit build output, credentials, local environment files, or unreviewed datasets. Dependency and toolchain updates should be explicit changes with the workspace checks rerun.

## Numerical development

The initial homogeneous API is documented in the [calculation contract](docs/calculation-contract.md). Extend that contract with new capabilities rather than assuming homogeneous roots are phase-equilibrium results. Keep GUI dependencies out of the calculation crates.

Follow the [testing strategy](docs/testing-strategy.md) for test layers, placement, independent reference cases, tolerances, CI, and completion criteria. Add tests with each implemented capability and regression cases with reproducible bug fixes. Review regeneration of expected values explicitly; do not replace them with our own engine's output to resolve failures. Generated tests must be reproducible using recorded seeds and inputs.

Every future reference case must record its composition basis, units, model variant, parameters, provenance, reference version, and tolerances. Distinguish implementation verification, model comparisons, and experimental validation. Do not copy historical pass counts into current test results.

Qualify data provenance and usage terms before adding third-party parameters or code. Do not silently fill missing interaction parameters with zero. Document any deliberately selected zero-interaction approximation in the case that uses it.
