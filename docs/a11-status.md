# A11 — Tables and splines

Implementation and measurement record, November 8, 2026. Issues #59 (A11),
#60 (profile), #61 (bounded single-phase tabulation experiment), #62 (fallback).

A11 adds the experimental crate `property-tables`. For **one fixed mixture** of
the exact PR model it builds a thermodynamically consistent Gibbs-energy table
over an explicitly qualified single-phase domain, measures value and derivative
errors separately against the direct model, measures the actual speedup, and
uses the direct model or a clear error outside the approved domain. It does not
change the models, solvers, flowsheet engine, package or project format; no
flowsheet path uses the table.

## Reproduce

```sh
cargo run --release -p property-tables --example a11_benchmark --locked
cargo run --release -p property-tables --example a11_benchmark --locked -- --quick   # smoke run
cargo test -p property-tables --locked
```

The benchmark prints the profile, table build/qualification evidence, timing
(min/median/max of 11 repetitions), accuracy over fixed-seed samples, a grid
refinement study and the fallback outcomes. All inputs are fixed seeds; timings
are wall-clock and machine-dependent, and are **not** a pull-request gate
(testing strategy, section 6). Accuracy targets are gated in the tests.

## Targets (fixed before the grid was chosen)

| class | quantity | target |
|---|---|---|
| value | molar enthalpy h | ≤ 0.1 J/mol absolute |
| value | molar entropy s | ≤ 1e-3 J/(mol K) absolute |
| value | molar volume / mass density | ≤ 1e-6 relative |
| derivative | cp, cv, (∂v/∂T)_p, (∂v/∂p)_T, speed of sound | ≤ 1e-3 relative |
| speed | table vs `CubicModel::evaluate`, same property set | ≥ 10× median |

The first pilot grid (29×31) missed the density target (1.4e-5). The grid was
refined until the targets held, not the targets relaxed (refinement study below).

## Experiment definition

- **Mixture:** dry natural gas of `examples/homogeneous-natural-gas.json`:
  methane/ethane/propane/nitrogen/CO2 = 0.88/0.06/0.03/0.02/0.01, **exact PR**
  (`PR-exact-1976`), explicit synthetic zero kij (not fitted VLE data).
  PR-BM is rejected by `GibbsTable::build`: its alpha second derivative is
  one-sided at each component Tc (ethane 305 K, CO2 304 K, propane 370 K lie in
  the domain), which would put cp kinks inside cells. Qualifying PR-BM would need
  node alignment at those temperatures and is not done.
- **Approved domain:** 280–420 K, 1–80 bar (`experiment::DOMAIN`), 57 temperature
  nodes (2.5 K) × 91 pressure nodes uniform in ln p (Δln p = ln 80/90 ≈ 0.0487).
- **Formulation:** the table stores molar Gibbs energy g(T, x = ln p) and its exact
  derivatives ∂g/∂T = −s, ∂g/∂x = p v, ∂²g/∂T∂x = p (∂v/∂T)_p from
  `CubicModel::evaluate` at every node, and interpolates with a tensor-product
  cubic Hermite surface (C¹ across cells). All properties derive from that one
  surface: s = −g_T, v = g_x/p, h = g + Ts, u = h − pv, cp = −T g_TT,
  (∂v/∂T)_p = g_Tx/p, (∂v/∂p)_T = (g_xx − g_x)/p², cv = cp + T(∂v/∂T)²/(∂v/∂p),
  w² = −v² cp /(M cv (∂v/∂p)_T). Hence h = g + Ts, u = h − pv, dh/dT = cp and the
  Maxwell relation (∂s/∂p)_T = −(∂v/∂T)_p hold on the table itself (tested),
  which independent property splines would not guarantee. cp, cv and (∂v/∂p)_T
  are second derivatives of the surface and are therefore only C⁰ across cell
  edges; their errors are reported separately.
- **Qualification (excludes phase boundaries and the critical region):** every
  point of a 4× refined lattice (81,225 points, including all nodes and edges)
  must have exactly one real cubic root, be stable in the default TPD trial set,
  and have reduced bulk modulus −v²(∂p/∂v)_T/(RT) ≥ 0.25 (1 for an ideal gas,
  0 at a critical point). Any failure aborts the build with
  `TableError::QualificationFailed` and the location/reason. Measured: minimum
  reduced bulk modulus 0.617; no TPD trial converged to a nontrivial stationary
  point. Finite TPD trials are evidence, not proof, of stability (as in A4);
  qualification between lattice points relies on the lattice density.

## Measured evidence

Machine: AMD Ryzen 7 PRO 4750U, 29.7 GB RAM, Windows 11 Pro; rustc/cargo 1.99.1;
`--release`, single thread. 5000 fixed-seed points for cheap calls, 500 for
stability/flash/PH; 11 repetitions; ns per call.

### Direct-calculation profile (#60)

| operation | min ns | median ns | max ns |
|---|---|---|---|
| `CubicModel::compressibility_roots` | 1578 | 1788 | 2257 |
| `CubicModel::ln_fugacity` | 2376 | 2416 | 2577 |
| `CubicModel::evaluate` (all properties, no phase check) | 4596 | 4876 | 5378 |
| `assess_stability` (default TPD trials) | 152116 | 163959 | 199431 |
| `tp_flash` (stability + single-phase properties) | 152167 | 174632 | 279038 |
| `ph_flash` (single-phase target h) | 6682738 | 6946639 | 8079159 |

Findings: in single-phase states the TP flash cost is ~94 % stability analysis;
`evaluate` itself is ~5 µs. PH flash averages 35.0 outer iterations, each a full
TP flash with stability, ≈ 7 ms per call. **Tabulation pays off most where a
caller repeatedly needs phase-qualified single-phase properties (TP flash) or PH
inversions for one composition in a known single-phase region**, e.g. repeated
heater/exchanger evaluations in recycle passes. Plain `evaluate` calls gain far
less in absolute time.

### Table cost and speed (#61)

Build including qualification: 13.5 s (dominated by 81,225 stability
assessments); node data 165,984 bytes (5187 nodes × 32 bytes). Break-even
against direct TP flashes is ≈ 13.5 s / 175 µs ≈ 77,000 calls; against PH flash
≈ 1,900 calls.

| operation | min ns | median ns | max ns |
|---|---|---|---|
| `GibbsTable::evaluate` (same property set as `evaluate`) | 127 | 155 | 183 |
| `TabulatedMixture::evaluate` (domain check + table) | 110 | 128 | 167 |
| `TabulatedMixture::ph` (table PH inversion) | 601 | 696 | 925 |

Median speedups: **31.5×** vs `CubicModel::evaluate` (target ≥ 10× met), 1129×
vs `tp_flash`, ≈ 10,000× for table PH vs `ph_flash`. A quick run on the same
machine gave 42.8× vs `evaluate`; timings vary between runs by tens of percent.

### Accuracy against the direct model (#61)

Values: g, h absolute J/mol; s absolute J/(mol K); v, ρ relative. Derivatives: relative.

| grid / points | n | g | h | s | v | ρ | cp | cv | (∂v/∂T)_p | (∂v/∂p)_T | w |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **57×91 approved, random** | 200,000 | 1.73e-5 | 1.78e-3 | 6.29e-6 | 6.39e-7 | 6.39e-7 | 1.11e-4 | 2.25e-4 | 3.15e-6 | 1.00e-4 | 1.04e-4 |
| 57×91 approved, cell centres | 5,040 | 1.30e-5 | 1.43e-4 | 4.64e-7 | 2.46e-9 | 2.46e-9 | 6.28e-5 | 6.99e-5 | 2.50e-8 | 5.38e-5 | 2.33e-5 |
| 29×31 (pilot) | 20,000 | 1.19e-3 | 1.88e-2 | 6.40e-5 | 1.41e-5 | 1.41e-5 | 4.08e-4 | 8.99e-4 | 5.62e-5 | 6.04e-4 | 5.50e-4 |
| 29×61 | 20,000 | 7.36e-5 | 1.26e-2 | 4.46e-5 | 1.93e-6 | 1.93e-6 | 4.11e-4 | 7.51e-4 | 1.25e-5 | 1.31e-4 | 2.21e-4 |
| 57×61 | 20,000 | 7.94e-5 | 2.18e-3 | 7.52e-6 | 1.92e-6 | 1.92e-6 | 9.70e-5 | 2.18e-4 | 6.87e-6 | 1.31e-4 | 1.31e-4 |
| 113×181 | 20,000 | 1.03e-6 | 1.83e-4 | 6.50e-7 | 7.57e-8 | 7.57e-8 | 2.57e-5 | 5.27e-5 | 3.49e-7 | 2.37e-5 | 2.57e-5 |

All approved-grid maxima meet the targets; density has the smallest margin
(0.64 of target). Volume error is governed by the ln p spacing, cp by the T
spacing (≈ second order), as expected for cubic Hermite second derivatives.
Maxima occur away from cell centres for first derivatives (v, s), so random
sampling, not only cell centres, is used for acceptance.

PH (500 targets from direct `evaluate`): max |T_table − T_direct| = 2.44e-5 K;
max |h_direct(T_table) − h_target| = 1.25e-3 J/mol.

## Outside the approved domain (#62)

`OutsidePolicy` is explicit per `TabulatedMixture`; results carry `Source::Table`
or `Source::DirectModel`, never a mix.

- `Error`: `TableError::OutsideApprovedDomain { T, p, domain }` ("… outside the
  approved table domain (280..=420 K, 100000..=8000000 Pa); use the direct model"),
  `EnthalpyOutsideApprovedDomain { p, h, min, max }` for PH, `CompositionMismatch`
  for any other composition (per-entry tolerance 1e-12).
- `DirectModel`: the stability-checked direct `tp_flash` / `ph_flash` (default
  options) at the requested state and composition. A single-phase result is
  returned as `Source::DirectModel`; a split returns
  `DirectResultNotSinglePhase { phases }` (e.g. 200 K, 30 bar → 2 phases).
  Direct solver errors propagate unchanged.
- Non-finite inputs are an input error under both policies. Domain bounds are
  inclusive; the table never extrapolates.

Measured outcomes (benchmark): 350 K/40 bar → table (both policies);
450 K/20 bar and 300 K/120 bar → error vs direct; 200 K/30 bar → error vs
`DirectResultNotSinglePhase { phases: 2 }`.

## Tests (`crates/property-tables/tests/tables.rs`)

Accuracy targets on 1500 fixed-seed points and all cell centres; node
reproduction; table consistency (h = g + Ts, u = h − pv, finite-difference
dh/dT = cp, Maxwell relation, finite-difference (∂v/∂p)_T); inclusive edges and
clear outside errors; non-finite inputs; direct fallback equality with
`evaluate` and source labelling; two-phase fallback error; composition mismatch
under both policies; table PH against `ph_flash`, PH range errors and PH
fallback; qualification rejection across the dew curve (200–300 K, 10–50 bar)
and by the critical-region bound; invalid domains and PR-BM rejection. The suite
builds the fully qualified table once (≈ 60 s in a debug test build).

No independent reference changes: the table is compared only with the existing
Rust direct model, whose Clapeyron qualification is unchanged; no fixture or
tolerance was modified.

## Limits and open work

- One fixed composition, exact PR, zero kij, one rectangle. Variable composition,
  PR-BM, GERG, NRTL and liquid or near-critical domains are not tabulated.
- Qualification checks a lattice; it is evidence, not a proof, of single-phase
  stability between lattice points, and shares A4's finite-trial caveat.
- Second-derivative properties are C⁰ across cell edges; a quintic Hermite or
  node-stored second derivatives would be needed for smooth cp.
- Not integrated into `property-packages`, the flowsheet or the project format;
  integration (and invalidation of tables on package changes) is separate work.
- Timings come from one laptop; repeat the benchmark on target hardware before
  relying on a speedup figure.
