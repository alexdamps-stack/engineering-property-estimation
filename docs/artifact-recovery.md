# A0 artifact recovery record

Search date: November 7, 2026. Scope: one bounded local recovery session against
the dispatched scaffold baseline `bf86da9`. Writable scope was limited to this
file and `artifact-inventory.md`; other workers' concurrent changes were not
treated as recovered historical artifacts.

## Result and decision

`ThermoCore_v0_3_0_Source.zip`, an identifiable ThermoCore Python source tree,
original formulation/parameter bundles, original test logs and inputs, and the
full earlier Rust solver/research reports were **not recovered**. This is a
statement about the inspected locations, not proof that the artifacts no longer
exist elsewhere or that historical download links are dead.

Proceed with **fresh Rust implementation, not a Python prototype port**. A port
cannot be reviewed without source, exact model conventions, parameters, usage
terms, and executable tests. Continue through A1 contract review and A2 source
qualification before model acceptance. An additional local Rust implementation
was found and recorded below as a review candidate; no adoption is authorized by
this recovery record. Stop searching now unless a concrete new artifact or
location is supplied; the plan's approximately three-day maximum is not a reason
to extend an unproductive search.

## Search boundaries and method

All paths below are Windows paths. Directory depths are relative to each listed
root (root = depth zero). Enumeration used PowerShell `Get-ChildItem -Force` and
an explicit breadth-first queue, a 10,000-file cap per root, and no traversal of
reparse points. `.git`, `target`, `node_modules`, `.venv`, `venv`, `.orca`, and
`.codex` directories were excluded; the first three searches also excluded
`.agents`, `ChatGPT`, `Codex`, and `deepseek-harness` directory names. Counts are
observed snapshots, not assertions about a shared checkout's later contents.

Filename matching was case-insensitive:
`thermo|gerg|peng.?robinson|nrtl|solver|formulation|research.*report|AGENTS\.md`.
The exact archive name also falls under `thermo`. Files without informative
names, nested archives, unrelated binaries, and unrelated document bodies were
not exhaustively searched. No archive was extracted and no discovered program
or generator was executed.

| Location inspected | Bounds and observed evidence | Recovery outcome |
|---|---|---|
| `C:\Users\PCUser\orca\workspaces` | Immediate directory names only, to identify nearby project storage | The named Property Estimation project was selected; unrelated projects were not traversed |
| `C:\Users\PCUser\orca\workspaces\Engineering Tools - Property Estimation` | Depth 8; 49 directories and 77 files listed; no depth truncation; includes `octopus` and `UI-Flowsheet-Design` | Thermo-named matches were scaffold crate directories; no archive/source bundle/report match |
| `C:\Users\PCUser\Downloads` | Depth 3; 3 directories and 109 files listed; no depth truncation | No ThermoCore/solver/formulation name match; one generic research-report candidate inspected as described below |
| `C:\Users\PCUser\Documents` | Depth 3; 13 directories and 14 files listed, with the exclusions above | No matching names; depth stopped at the 11 game-storage directories listed below |
| `C:\Users\PCUser\Documents\ChatGPT\Engineering Tools - Property Estimation` | Exact project directory discovered through `git rev-parse --git-common-dir`; separately enumerated to depth 8, 27 directories and 45 files; no depth truncation; `.git` contents excluded | Separate Rust implementation and reference files found; no ThermoCore archive or original reports |
| `C:\Users\PCUser\OneDrive\Documents` | `Test-Path` returned false | Ordinary alternative Documents location absent |
| `C:\workspace\scratch\2a408b01b5b5` | `Test-Path` returned false for this local analogue only | Does not test or access the historical remote Linux filesystem |

The Documents depth limit stopped below
`C:\Users\PCUser\Documents\My Games\Tabletop Simulator` at
`Mods\Assetbundles`, `Mods\Audio`, `Mods\Images`, `Mods\Images Raw`,
`Mods\Models`, `Mods\Models Raw`, `Mods\PDF`, `Mods\Text`,
`Mods\Translations`, `Mods\Workshop`, and `Saves\Saved Objects`.
These directories were not searched further.

The parent Documents `ChatGPT` directory was excluded from the ordinary search.
Only the exact Git project directory subsequently identified by Git was opened;
no sibling private chat exports, chat stores, account databases, secrets,
application profiles, or entire drives were searched. Downloads included an
unrelated thermodynamic-database-named installer archive (`DDB2025_Update1.zip`)
and a DWSIM installer; their names alone establish no recoverable ThermoCore
content or redistribution rights, and their contents were not opened. Generic
archives were not recursively unpacked. No network/download/private-conversation
access was attempted.

## Content evidence inspected

- Read the assigned `docs/parallel-work-plan.md`, `docs/development-plan.md`,
  `CONTRIBUTING.md`, and original `docs/artifact-inventory.md`. The English plan
  supplies the archive name, remote scratch path, historical counts, and private
  conversation links; these are secondary planning records, not recovered source.
- `git ls-tree -r --name-only bf86da9` listed scaffold files and documentation.
  Read this worktree's `crates/thermo-models/src/lib.rs` and
  `crates/thermo-solvers/src/lib.rs`; both explicitly say scaffold only.
- Searched immediate documentation files in `octopus/docs` and
  `UI-Flowsheet-Design/docs` with `Select-String` for `ThermoCore`,
  `2a408b01b5b5`, `53 tests`, `111 GERG`, and `solver report`. Hits repeated the
  planning/inventory claims and parallel assignment, not original evidence.
- In the exact original project directory, searched text files with extensions
  `.md`, `.txt`, `.py`, `.rs`, `.toml`, `.json`, `.csv`, each smaller than 2 MiB,
  for `ThermoCore|2a408b01b5b5|53 tests|111|solver report|Solver Research`.
  Hits were the German/English plans, old inventory, and incidental numeric
  `111` substrings in reference JSON; none recovered the historical artifacts.
- Read only the first 35 lines of
  `C:\Users\PCUser\Downloads\deep-research-report.md`; its heading and summary
  concern Hamburg consulting-office locations, so it was rejected as unrelated.
  No unrelated content is reproduced here.
- Checked `AGENTS.md` at the workspace and each ancestor through `C:\`, and
  searched its filename in the enumerated project trees; none was found.
- Loaded the installed orchestration skill and `orca skills get orchestration`
  for worker obligations. No helper workers were started.

## Separate Rust candidate: provenance and limits

The original project directory is
`C:\Users\PCUser\Documents\ChatGPT\Engineering Tools - Property Estimation`.
Its read-only Git inspection reported a clean working tree and HEAD
`0985d1ac85cdbddc386467286044f20489964fc6`, with subject
"Implement homogeneous PR and PR-BM engine with independent reference tests"
and recorded author date `2026-11-01T20:29:03+01:00`. Git dates and descriptions
are local provenance metadata, not independent proof of correctness or authorship.
The last four commit summaries were inspected; no history checkout was made.

Inspected `data/NOTICE.md`, `validation/README.md`, and
`validation/requirements-reference.txt` in full; the first 45 lines of
`crates/thermo-models/src/lib.rs`; first 55 lines of
`validation/generate_cubic_references.py`; first 38 lines and parsed case count
of `validation/cubic-reference.json`. Observed:

- The Rust module identifies homogeneous PR76 and Boston–Mathias supercritical
  continuation, model version `pr76-exact-constants-bm-v1`, and explicitly omits
  phase stability, flash, and ideal caloric contribution.
- The JSON contains 16 cases, labels its purpose homogeneous implementation
  verification only, and declares teqp 0.23.2, CoolProp 8.0.0, NumPy 2.5.4,
  SciPy 1.18.1, and Python 3.14.4. These are file declarations, not verified
  installed versions or regenerated reference results.
- The notice declares selected CoolProp constants and an upstream MIT notice;
  it explicitly does not select a license for original repository code. Usage
  terms and source matching require review before any import.
- The validation README reports earlier ordinary tests and a 200-case grid.
  No such tests were rerun here, and these claims do not qualify A2/A3, establish
  experimental accuracy, or replace the unavailable ThermoCore results.

SHA-256 digests of the inspected candidate files anchor this observation:

| Relative path in original checkout | SHA-256 |
|---|---|
| `crates/thermo-models/src/lib.rs` | `2CC2AAAE16B0EDBBDFAECB955F7361F114146B86C2713F566CF9A78720B891E7` |
| `validation/generate_cubic_references.py` | `CE32F14F810D51D5A9488F777BC112700FD05974617DF03A9EDA31AB10591385` |
| `validation/cubic-reference.json` | `5FAC3485E03B65CDE261B3ECBEEC6E3C8944E4BD7781AAB9972EBB05E8641244` |
| `data/NOTICE.md` | `8C6FFD7FAD9199610E748476146441D78EFCF83D17F04036FF36DCDF5B64D10F` |

No source, constants, fixtures, or original license text were copied into this
worktree. Coordinator guidance is required before importing third-party code or
making any change outside the two assigned files. The separate candidate may
be useful for later review without changing the decision that an original
ThermoCore port is presently unsupported.

## Later recovery compatibility and open evidence

If original material is supplied later, preserve its original archive bytes,
hash, filename, acquisition date, source/owner, version, and permission evidence
before extraction. Inspect archive paths and file types safely without running
embedded code; arrange coordinator-approved storage before copying it here.
Keep original tests/reports alongside their inputs, dependency versions, and
unmodified expected outputs. Missing licenses or parameter usage terms remain
blockers to reuse, not reasons to infer permission.

Map recovered formulations to the reviewed Rust contract: PR versus PR-BM alpha
variant and derivatives, cubic constants, mixing and directed interaction
conventions, units and composition basis, caloric reference states, ideal/residual
contributions, liquid-volume treatment, validity domain, and solver termination
criteria. Preserve model/data version identifiers and use explicit adapters or
conversion records; do not silently change established conventions to match a
prototype. A teqp GERG adapter remains a reference path, not a native Rust GERG
implementation. Recovered expected values require independence and formulation
review before becoming accepted fixtures; historical counts stay historical.

The coordinator can review A0 on the explicit-unavailability branch. Missing
ThermoCore artifacts do not block A1 contract work under this decision. Missing
qualified formulas, datasets, independent references, and executable current
checks still gate later numerical acceptance. Candidate reuse is a separate
coordinator decision; no model implementation or acceptance is claimed here.

## Checks and command outcomes

The initial discovery command attempted `rg`, which is not installed: the
PowerShell tool result exited 1 with `CommandNotFoundException`. This was a
**failed discovery check**, not evidence of no matches. Subsequent bounded
PowerShell enumeration and `Select-String` replaced that check successfully.
The same initial command's `Get-Item -ErrorAction SilentlyContinue` probes were
not used as absence evidence; the successful explicit `Test-Path` checks above
replaced them. All subsequent search/inspection command tool results exited 0.

Git tree/status/log inspection, JSON parsing/counting, SHA-256 hashing, and
`git diff --check` completed successfully. Documentation-only changes require
no model test run; no numerical validation, source execution, package install,
archive extraction, staging, commit, reset, or checkout was performed.
