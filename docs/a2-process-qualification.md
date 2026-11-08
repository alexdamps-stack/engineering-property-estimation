# A2 demonstration phase and state qualification

Reviewed 2026-11-08 for issues #15 and #9. The historical design input is
`examples/process.provisional.json`; its null expectations are superseded by
`validation/fixtures/process-v1.json`, not used as the current qualification.
The original feed rates, compositions and setpoints are retained unchanged.

The new independent adapter uses pinned teqp 0.23.2, Cantera 3.2.0, NumPy 2.3.4
and SciPy 1.16.3. It reuses the reviewed A3 caloric and A4 equilibrium adapters:
48 seeded random TPD starts plus Wilson starts, Powell phase-split solves and
Brent energy roots. No Rust calculation supplies expected values. The fixture
records input and adapter hashes, phase compositions, molar enthalpies, flows,
stability evidence and cooler duty for both exact PR and PR-BM separately.
Reference h uses pure ideal h=0 at 298.15 K; entropy uses 100000 Pa and the
same source/reference convention as A3. All interactions are explicit synthetic
zeros and volumes are untranslated. This qualifies implementation agreement,
not measured natural-gas accuracy or global stability.

| Stream | PR T (K) | PR-BM T (K) | Absolute p (MPa) | Expected fluid phase | Flow (mol/s) |
|---|---:|---:|---:|---|---:|
| A | 330 | 330 | 5 | single gas-like fluid | 100 |
| B | 350 | 350 | 5 | single gas-like fluid | 20 |
| Mixer outlet | 333.698837612 | 333.697581743 | 5 | single gas-like fluid | 120 |
| Cooler outlet | 240 | 240 | 5 | vapor + liquid | 120 |
| Valve outlet | 214.110342594 | 213.981036124 | 1.5 | vapor + liquid | 120 |
| Separator vapor | same as valve | same as valve | 1.5 | saturated vapor | 115.799201135 / 115.754245346 |
| Separator liquid | same as valve | same as valve | 1.5 | saturated liquid | 4.200798865 / 4.245754654 |

Paired flows are PR / PR-BM. Cooler vapor fractions are 0.993206358932 /
0.994149910191; valve vapor fractions are 0.964993342793 / 0.964618711217.
The mixer temperature solves the flow-weighted inlet enthalpy balance. The
valve independently solves constant molar enthalpy at 1.5 MPa. Separator
compositions and flows use the valve equilibrium split with no further heat,
work or pressure drop. The separated streams are phase endpoints, not interior
two-phase feeds; finite-start TPD confirms no lower fluid split at those inputs.

For every stream, the allowed computational/search envelope is 200..400 K and
1..5 MPa, intersected with the cubic caloric domain and explicit solver
diagnostics. The independently qualified set consists of the exact T/p/x
points above, with the fixture's 2e-6 K comparison budget for energy roots;
the envelope is **not** a validated continuous phase map. Feed and cooler inputs
remain fixed; changing them requires fresh qualification. Separator endpoint
compositions are specific to each variant. Do not extrapolate phase labels
over the search envelope. Water, methanol, solids, additional liquid phases,
reactions and critical-neighborhood claims are excluded. A13 owns the future
bubble/dew and phase-envelope range qualification; A8 owns execution of the
connected flowsheet through unit-operation APIs.

Reproduce without replacing the committed fixture:

```sh
python validation/generate-process-references.py --output /tmp/process-candidate.json
cargo test -p thermo-solvers --test process_reference --locked
```

The Rust regression compares all 14 states, including phase counts, phase
compositions/fractions and enthalpies; PH recovery of mixer/valve temperatures;
mixing and separator component/energy balances; valve enthalpy and cooler duty.
Budgets are the existing A5 budgets, fixed before these comparisons. No earlier
fixture is replaced and no tolerance is widened. See the final acceptance and
check record in `docs/a1-a2-acceptance.md`.
