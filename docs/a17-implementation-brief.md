# A17 implementation brief (shared by all contributors)

Work package A17 (#170) implements the 48 missing palette items listed in
[unit-operation-gap-analysis.md](unit-operation-gap-analysis.md), one GitHub
sub-issue each (#171–#218). Design inspiration: DWSIM (open source, readable
reference for formulations and specification sets) and the Aspen suite
(HYSYS/Aspen Plus specification sets and user expectations). Never copy
proprietary text or code; cite open literature for every correlation.

Solids handling (#198–#211) is deferred from this palette effort. The
foundation's project records for existing solids units remain supported.

## Rules

1. Read `AGENTS.md`, `CONTRIBUTING.md`, `docs/flowsheet-and-ui.md`,
   `docs/testing-strategy.md` and, once it exists, `docs/a17-architecture.md`
   (the extension conventions written by the foundation step). Follow them.
2. **Build and test in WSL2 Ubuntu**, not Windows: Smart App Control blocks
   cargo on Windows (os error 4551). Use
   `wsl -e bash -lc 'cd "$(wslpath "<worktree>")" && CARGO_TARGET_DIR=$HOME/targets/<your-branch> cargo test -p <crate>'`.
   Use your own `CARGO_TARGET_DIR` so parallel agents do not block each other.
   Before finishing run `cargo fmt --all`, `cargo clippy --workspace --all-targets -- -D warnings`
   and `cargo test --workspace` (or the crates you touched plus their dependents
   if the full run exceeds your time; say which).
3. For changes to PR/PR-BM, NRTL, caloric data or equilibrium solvers also run
   `julia --project=validation/julia validation/julia/check.jl` (Julia 1.11.7,
   in WSL). Keep reference generation independent of the Rust engine
   (Python/Julia scripts under `validation/`). Never replace committed fixtures
   or widen tolerances to make a check pass.
4. **Ownership:** edit only your group's files plus the small registration
   hook points named in `docs/a17-architecture.md`. If you must change a shared
   file (engine core, `PortSpec`, project format version, package API), keep
   the change minimal and list it in your report.
5. Every unit: ports, mutually exclusive specifications as enums, engine
   calculation with component and energy balance, explicit failures (no silent
   fallbacks), tests of behaviour and failure paths, independent reference
   values where a reference exists, project-format persistence (shared format
   version bump is done once by the integrator; add your records behind the
   current version and note it), palette entry with an ISO 10628 symbol from
   `docs/design-reference/iso-10628` (keep attribution), and a spec form in
   the desktop app.
6. A unit that cannot be implemented honestly in steady state with the
   available models is not faked: keep it drawable, `calculation_planned`
   false, document why in `docs/a17-status.md` under your group, and say so.
7. Commit on your branch with messages that reference the sub-issue numbers
   (`Refs #186`), ending with
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Do not push,
   open PRs or close issues; the orchestrator integrates.
8. Report back in at most 250 words: what was implemented per issue number,
   what is drawn-only and why, shared files touched, check results (exact
   commands and pass/fail), open risks.

## Columns and reactors as packed beds

Columns and reactors must support packed beds, with configurable modelling:

- **Columns** (distillation, absorber/stripper, extractor, 3-phase): each
  section is either *trayed* or *packed*. Calculation mode per section:
  - *Equilibrium (ideal theoretical stages)*: number of stages, optional
    Murphree/overall efficiency; for packing, theoretical stages from packed
    height / HETP (user HETP or a correlation).
  - *Rate-based (mass transfer)*: packed height, column diameter, packing
    type and data (specific area, void fraction, nominal size, packing factor;
    a small library of random and structured packings with cited data),
    mass-transfer correlation (Onda 1968, Billet–Schultes 1999, Bravo–Rocha–Fair
    for structured), film model, segment count; hydraulics (pressure drop,
    flooding fraction) reported. DWSIM's rate-based absorber and Aspen
    RateSep are the user-facing references.
- **Reactors** (PFR, CSTR, Gibbs/equilibrium where relevant): homogeneous or
  *packed (catalytic) bed*: bed length/diameter or catalyst mass, bed void
  fraction, particle diameter and density, Ergun pressure drop, rate basis
  per catalyst mass or per reactor volume, optional effectiveness factor,
  adiabatic / isothermal / specified-duty / heat-transfer-coefficient
  thermal modes, multiple tubes. Kinetics: power law and LHHW, Arrhenius
  temperature dependence, basis (concentration or partial pressure).
