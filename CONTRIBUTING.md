# Contributing

Use English for new documentation, source comments, and commit messages. Keep the German original plan as a historical baseline; update the English development plan when scope or acceptance criteria change.

## Local workflow

1. Install the toolchain selected by `rust-toolchain.toml` and your platform's native linker prerequisites.
2. Make a focused change on a feature branch.
3. Run the formatting, Clippy, test, and documentation commands in the root README.
4. Refresh the installed local desktop prototype and check its source version.
5. Describe behavior, evidence, and remaining limits in the pull request.

Commit the root `Cargo.lock`. Do not commit build output, credentials, local environment files, or unreviewed datasets. Dependency and toolchain updates should be explicit changes with the workspace checks rerun.

### Local desktop update

After codebase changes on the local Windows/Ubuntu machine, run
`powershell -File scripts/update-local-app.ps1`, followed by the same command
with `-Check`. The check verifies the source fingerprint, installed binary and
Windows shortcuts; it returns exit code 1 for a missing or stale installation.
Hosts without this environment should report the update as unavailable.

The updater builds with the pinned Rust toolchain in WSL Ubuntu, self-tests a
candidate and atomically switches the installed `current` symlink after success.
Failed builds and self-tests preserve the previous installation. An update lock
prevents competing installers; edits during a build prevent publication. The
source fingerprint includes tracked and nonignored untracked files, including
local edits. `build.json` records the source checkout, commit and binary hash.
Previous releases remain in
`~/.local/share/property-estimation-prototype/releases/` for recovery.

Save open projects and reopen **Property Estimation Prototype** from the desktop
to use the new build. The updater leaves running windows open. This is a required
workflow step; there is no file watcher or automatic restart.

Use `-Distribution` to select another WSL distribution. From Ubuntu, the updater
can also run as `python3 scripts/update-local-app.py`; use PowerShell to set up or
repair the Windows shortcuts. Run its regression checks in Ubuntu with
`python3 -m unittest discover -s scripts/tests -v`.

## Numerical development

The PR/PR-BM numerical APIs and TP/PH flash contracts are documented in the
milestone status records. Keep GUI dependencies out of the calculation crates.

The Clapeyron comparison is a standard numerical check for supported PR/PR-BM
changes. Run `julia --project=validation/julia validation/julia/check.jl` with
Julia 1.11.7; the locked environment instantiates automatically. Workspace
tests also run the committed Clapeyron fixtures offline. Both checks must pass;
generation or comparison failures must remain visible. Review reference,
tolerance and Julia dependency updates explicitly. See
[Clapeyron checks](validation/clapeyron-checks.md).

Follow the [testing strategy](docs/testing-strategy.md) for test layers, placement, independent reference cases, tolerances, CI, and completion criteria. Add tests with each implemented capability and regression cases with reproducible bug fixes. Review regeneration of expected values explicitly; do not replace them with our own engine's output to resolve failures. Generated tests must be reproducible using recorded seeds and inputs.

Every future reference case must record its composition basis, units, model variant, parameters, provenance, reference version, and tolerances. Distinguish implementation verification, model comparisons, and experimental validation. Do not copy historical pass counts into current test results.

Qualify data provenance and usage terms before adding third-party parameters or code. Do not silently fill missing interaction parameters with zero. Document any deliberately selected zero-interaction approximation in the case that uses it.
