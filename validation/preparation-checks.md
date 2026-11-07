# A2 preparation evidence

Date 2026-11-07; baseline `bf86da9`; Windows PowerShell 5.1.26100.9711.
No stage, commit, reset, checkout, engine/runtime dependency or source-crate
mutation was performed. Other workers' concurrent changes are outside this report.

| Check / command | Outcome |
|---|---|
| `orca skills get orchestration` | Passed; worker obligations loaded; no helpers spawned |
| Initial `rg --files ...` documentation discovery | Failed: rg CommandNotFoundException; replaced with bounded PowerShell Get-ChildItem and Get-Content successfully |
| `python --version` and attempted `python -c` package import | Deliberate reference-availability probes: Windows Store alias reports Python unavailable; no Python reference execution |
| `Get-Command julia,py,node -ErrorAction SilentlyContinue` | Deliberate availability probe: none found, no Julia evaluation |
| Pinned GitHub source and license HTTP retrieval | Passed for final four generator URLs; tag resolution pinned CoolProp and Cantera commits |
| PR/BM candidate source reading | Passed; exact formulas inspected; no execution or copying of local candidate code/data |
| `https://pypi.org/pypi/teqp/0.23.2/json` and `https://pypi.org/pypi/CoolProp/8.0.0/json` | Deliberate availability probes: both HTTP 200 with matching version fields; not installed/generation evidence |
| Clapeyron license discovery at `.../v0.6.30/LICENSE` | HTTP 404; corrected to pinned `LICENSE.md`, successfully inspected; no inference of permission from the failed URL |
| `powershell -NoProfile -ExecutionPolicy Bypass -File validation/prepare-source-data.ps1 -OutputDirectory data/source-extracts` | Passed, exit 0; eight source records and two notices |
| Same generator targeting a fresh temporary directory; `Get-FileHash -Algorithm SHA256` comparison of all three output files | Passed, exit 0; byte-identical outputs; temporary files and empty directory removed |
| `powershell -NoProfile -ExecutionPolicy Bypass -File validation/check-preparation.ps1` | Passed, exit 0; source records, NASA7 shape/ranges, positive/finite constants, IDs, source pins, two 9-state grids/model, dry feeds and unit-bearing proposals; three malformed compositions rejected |
| `git diff --check -- data validation examples docs/reference-qualification.md` | Passed, exit 0; rerun after final documentation changes |

The temporary-directory reproduction command used
`$reproDirectory=Join-Path $env:TEMP ('a2-source-reproduction-'+[guid]::NewGuid().ToString())`,
ran the generator with `-OutputDirectory $reproDirectory`, then compared each
named file's SHA-256 with `data/source-extracts`. No random thermodynamic inputs
were generated; the GUID chooses an isolated temporary directory only.

Reproduced file hashes:

| File | SHA-256 |
|---|---|
| `data/source-extracts/components.provisional.json` | `CF685B510B3C5A319B2D2BC866A3E0FC67F222888617CFFCE9CEEE5169A7108C` |
| `data/source-extracts/CoolProp-LICENSE.txt` | `9BF835333EF602AF4CB19338B9F9D43671E174FA029B00280E9BDBA6EA4719B2` |
| `data/source-extracts/Cantera-License.txt` | `E92980B9712CE20E73898A97B0116889E84E07F548D6BE8591E87DCAD79C41BB` |

Changed files owned by A2:

- `data/README.md`
- `data/source-extracts/components.provisional.json`
- `data/source-extracts/CoolProp-LICENSE.txt`
- `data/source-extracts/Cantera-License.txt`
- `examples/README.md`
- `examples/cases.provisional.json`
- `examples/process.provisional.json`
- `validation/README.md`
- `validation/prepare-source-data.ps1`
- `validation/check-preparation.ps1`
- `validation/tolerances.provisional.json`
- `validation/preparation-checks.md`
- `docs/reference-qualification.md`

Exact source URLs/commits, formulas, permission evidence, independence limitations,
candidate review, domains and next actions are in
`docs/reference-qualification.md` and the source JSON. No source is qualified as
experimental validation. No expected thermodynamic values were generated.

Remaining gates: A1/coordinator review, exact versus rounded cubic factors,
shared R/caloric offsets, approved external reference runtime and source audit,
independent expected values for both PR and PR-BM, phase/stability/PH qualification,
saturation state selection and evidence-based tolerance approval. The original
local candidate's redistribution terms are unresolved, so it was not imported.
These are reported open gaps, not passing numerical checks or A2 acceptance.
