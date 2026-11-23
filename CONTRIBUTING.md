# Contributing

Use English for new documentation, source comments, and commit messages. Keep the German original plan as a historical baseline; update the English development plan when scope or acceptance criteria change.

## Local workflow

1. Install the toolchain selected by `rust-toolchain.toml` and your platform's native linker prerequisites.
2. Make a focused change on a feature branch.
3. Run the formatting, Clippy, test, and documentation commands in the root README.
4. Refresh the local desktop prototype and verify it using the checklist below.
5. Describe behavior, evidence, and remaining limits in the pull request.

### Local desktop update checklist

On the local Windows/Ubuntu development machine, complete these steps after
codebase changes (including engine, GUI, data, dependencies and examples):

- Run the appropriate source checks from the README and testing strategy.
- Run `powershell -File scripts/update-local-app.ps1` from this checkout.
- Run `powershell -File scripts/update-local-app.ps1 -Check`; a stale or missing
  installation returns exit code 1. Report update errors instead of claiming
  the installed application was refreshed.
- Save work in any open app window, close it, and reopen **Property Estimation
  Prototype** on the desktop to use the new build. Agents must not force-close
  app windows or discard unsaved projects.

The updater builds with the pinned toolchain in Ubuntu, runs the desktop
acceptance self-test on the candidate, and atomically switches the installed
`current` symlink only after success. Build or self-test failure leaves the
previous installation selected. An update lock prevents competing installers.
Running windows continue using their existing executable. No file watcher or
automatic restart is installed.

`build.json` beside the installed executable records the source checkout,
commit, a SHA-256 of Git-tracked and nonignored untracked files (including local
edits), and the installed binary hash. The check verifies both source and binary
hashes. If files change during the build, the updater refuses to publish it;
run it again after editing finishes. Previous releases remain under
`~/.local/share/property-estimation-prototype/releases/` for recovery.
The PowerShell check also verifies that both Windows shortcuts target this
installation. Run updater regression checks in Ubuntu with
`python3 -m unittest discover -s scripts/tests -v`.

The entry point is [update-local-app.ps1](scripts/update-local-app.ps1);
the WSL implementation is [update-local-app.py](scripts/update-local-app.py).
From Ubuntu, run `python3 scripts/update-local-app.py` directly. Use the
PowerShell entry point for initial installation or to repair the Windows
shortcuts. `-Distribution` selects a WSL distribution other than Ubuntu.
This local installation step does not replace CI or the numerical reference
checks; hosts without this desktop environment should report it as unavailable.

### AI-assisted workflow with pstack

The repository uses a local [pstack skill](.agents/skills/pstack/SKILL.md),
adapted from [cursor/plugins pstack](https://github.com/cursor/plugins/tree/main/pstack).
`AGENTS.md` makes it part of the normal AI contributor workflow. In clients
that discover project skills, invoke `$pstack`; otherwise read the linked
`SKILL.md` directly. No global plugin installation is required for this adapter.

- **Investigate or sweep:** trace behavior, reproduce suspected failures, and
  report actionable issues with file/line evidence, impact, and verification.
  An issue-finding request produces a report; apply fixes when requested.
- **Fix a bug:** reproduce the failure, identify its cause, make a focused fix,
  and rerun the original reproduction plus relevant regression checks.
- **Build or refactor:** trace callers and domain types, define observable
  acceptance criteria, implement in verifiable units, and check the actual result.

The local adapter retains the project's Rust, numerical-reference, and UI checks.
It uses available session tools and models, without requiring Cursor commands,
external review panels, or automatic PR creation. Read its provenance section
before updating it; the adapted workflow is maintained locally rather than
downloaded from upstream on every task.

Commit the root `Cargo.lock`. Do not commit build output, credentials, local environment files, or unreviewed datasets. Dependency and toolchain updates should be explicit changes with the workspace checks rerun.

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
