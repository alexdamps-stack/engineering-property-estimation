# Artifact inventory

Updated: November 1, 2026.

| Artifact | Status | Next action |
|---|---|---|
| German development plan | Preserved at repository root | Retain as historical baseline |
| English development plan | Included in `docs/development-plan.md` | Maintain as current planning document |
| Rust workspace | First homogeneous PR/PR-BM engine and JSON CLI | Extend the documented calculation contract toward total caloric properties and phase equilibrium |
| ThermoCore v0.3.0 Python source | Referenced in historical conversation; not recovered locally | Recover `ThermoCore_v0_3_0_Source.zip` and inspect contents under A0 |
| Earlier formulations and parameter sources | Described in earlier reports; not qualified in this repository | Recover, identify versions, and review provenance |
| Earlier test results | Historical report: 53 tests and 111 GERG comparisons over 37 gas states | Reproduce only after recovering source and reference inputs |
| Rust solver/research reports | Some assistant content was not readable beyond references | Recover full reports if available; do not assume implementation exists |
| Component constants | Eight selected records from CoolProp 8.0.0 with source hashes and notice | Qualify physical applicability and add required interaction/ideal caloric data under A2 |
| Independent numerical references | Sixteen teqp 0.23.2 homogeneous cases and reproducible generator included | Broaden model/state coverage; preserve independence and reviewed regeneration |

The historical `/workspace/scratch/2a408b01b5b5` location is not a path available in this local repository. This inventory does not assert that linked archives or conversation content remain downloadable.

Fresh-implementation decision for A0: no earlier source archive is locally available. The current Rust implementation was developed from documented formulations and checked against independently generated teqp fixtures. It does not claim to port or reproduce the unavailable prototype's entire feature set or historical test suite.
