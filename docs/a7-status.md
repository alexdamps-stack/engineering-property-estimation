# A7 implementation and verification: NRTL and liquid properties

Date: 2026-11-08. Starting revision: `cfdcb01` (A5 merged). Scope: an NRTL
liquid with an explicitly selected ideal-gas vapor model (modified Raoult),
pure-liquid correlations, caloric coupling to the existing NASA7 reference, and
gamma–phi bubble/dew points and TP/PH flash. The first qualified system is
methanol/water. Executed checks are recorded in
[the A7 verification report](../validation/a7-checks.md).

## Data and provenance

`data/source-extracts/liquid-nrtl.provisional.json` is produced only by the
explicitly invoked, standard-library
[`prepare-nrtl-source-data.py`](../validation/prepare-nrtl-source-data.py). It
downloads pinned files, records their SHA-256 hashes and transcribes:

| Data | Source (pinned) | Form |
|---|---|---|
| Vapor pressure, water and methanol | CoolProp `98b3523` `dev/fluids/{Water,Methanol}.json`, `pS` ancillary (MIT) | `p = p_r exp((T_r/T) Σ n_k θ^t_k)`, θ = 1 − T/T_r; declared 273.16–647.096 K and 175.61–512.49 K |
| Saturated-liquid density | same files, `rhoL` ancillary | `ρ' = ρ_r (1 + Σ n_k θ^t_k)` mol/m³ |
| Methanol/water NRTL | ChemSep IPD "DECHEMA NRTL data @ 1atm", row `1101;1921;-189.0469;792.802;0.2999;Methanol/Water p61 1/1a`, via DWSIM `1abf72d` `nrtl.dat`; cross-checked with `thermo` `2bb466e` `ChemSep/nrtl.json` (MIT) | `τ_ij = A_ij/(R T)`, A in cal/mol |

The cal/mol values are converted with the thermochemical calorie (4.184 J) and
the A1 gas constant: `t1_ij = A_ij · 4.184 / R`, giving `t1(methanol→water) =
-95.132 K` and `t1(water→methanol) = 398.953 K`. The script refuses the extract
unless `thermo`'s stored `bij` reproduces this conversion to 1e-12 relative and
the ChemSep header states 1 atm and cal/mol. The CoolProp ancillaries report
0.014 % maximum deviation from their reference EOS for vapor pressure; that is a
fit statement from the source, not a qualification by this project.

The regression interval of the 1 atm data is recorded as 337.69–373.13 K: the
pure normal boiling points from the transcribed vapor-pressure ancillaries,
rounded outward. The source gives no explicit interval.

`thermo-data::LiquidDatabase` loads the extract, validates its header, shapes,
ranges and duplicate pairs, and orients a stored (i, j) record for a request
(j, i). A pair or pure-liquid record that is absent returns `MissingParameter`.
No zero interaction is ever substituted.

## Formulation

**NRTL (`thermo-models::NrtlModel`).**

```
g^E/(RT) = Σ_i x_i (Σ_j x_j τ_ji G_ji) / (Σ_k x_k G_ki),   G_ji = exp(−α_ji τ_ji)
τ_ij = t0_ij + t1_ij/T + t2_ij ln(T/K) + t3_ij T             (directed, τ_ij ≠ τ_ji)
α_ij = α_ji = a0_ij + a1_ij T
```

This matches Clapeyron's `aspenNRTL` (and `NRTL` when only t0, t1 and a0 are
used). Diagonal terms must be exactly zero and α must be symmetric. A
nonpositive α(T) returns `OutsideValidityDomain`. `ln γ_i`, `∂ln γ_i/∂T`,
`g^E`, `h^E = −RT² ∂(g^E/RT)/∂T`, `s^E` and `c_p^E` are propagated analytically
with second-order temperature jets. There is no finite differencing.
Explicit parameter sets need a nonblank source.

Validity reports whether T lies inside every interacting pair's regression
interval (`WithinDeclaredDomain`) or outside one (`Extrapolated`). A pair
interacts when both components are present. Explicit parameters without an
interval report `Unassessed`. So does a liquid with no interacting pair (a
single active component, where g^E = 0), because no interval applies.

**Gamma–phi package (`thermo-models::ActivityModel`).** The vapor model is a
required `VaporModel` argument with no default. `IdealGas` is the only variant.
Equilibrium is `y_i p = x_i γ_i p_i^sat`, with no Poynting factor and no φ^sat.
To match this selection, liquid volume is neglected in the Gibbs energy, so the
pure-liquid standard state is `g_i^L(T) = g_i^ig(T, p0) + RT ln(p_i^sat/p0)`,
which gives:

```
h_i^L  = h_i^ig − R T² dln p^sat/dT                   (Clausius–Clapeyron)
s_i^L  = s_i^ig°(T) − R ln(p^sat/p0) − R T dln p^sat/dT
c_p,i^L = c_p,i^ig − R (2T L' + T² L'')
h^L = Σ x_i h_i^L + h^E,  s^L = Σ x_i s_i^L − R Σ x_i ln x_i + s^E
```

Liquid h, s and g are pressure independent. The caloric reference is the A3/A5
one: pure ideal gas h = s = 0 at 298.15 K and 1 bar, NASA7 integrals. Flowsheet
streams can therefore move between the cubic and gamma–phi packages without a
reference shift. The latent heat implied by the vapor-pressure slope reflects
the selected ideal-gas/zero-volume model. For water at 373 K it is about
41.3 kJ/mol, against the measured 40.7 kJ/mol. This is model error, not an
implementation discrepancy.

Liquid molar volume comes from the saturated-liquid density ancillaries with
ideal (Amagat) mixing. It is reported as a separately labelled property and is
**not** coupled into fugacities or caloric properties, as the A1 contract
requires.

Computational envelope: 273.16–450 K and 100 Pa–500 kPa, and every pure
correlation must cover T. Outside it, the package returns `OutsideValidityDomain`.

## Solvers (`thermo-solvers::activity`)

- **Bubble pressure** is explicit: `p = Σ x_i γ_i p_i^sat`.
- **Dew pressure** solves the incipient liquid by successive substitution on γ.
- **Bubble/dew temperature** use safeguarded Newton on `ln p_sat,mix(T) − ln p`
  over the option bracket. The dew slope `Σ x_i (∂ln γ_i/∂T + dln p_i^sat/dT)`
  is exact at fixed y by Gibbs–Duhem. If no root is bracketed, the result is
  `OutsideValidityDomain`.
- **TP flash.** With an ideal-gas vapor and pressure-independent liquid
  fugacities, the feed is liquid for `p ≥ p_bubble(T, z)`, vapor for
  `p ≤ p_dew(T, z)`, and two-phase between them. The split iterates γ(x) with
  the existing negative-flash Rachford–Rice, converging on
  `max |ln(y_i p) − ln(x_i γ_i p_i^sat)| ≤ 1e-13`. That is near machine
  precision, so derived phase properties such as the liquid volume inherit
  only ~1e-14 composition error. Every liquid result passes
  finite near-pure NRTL tangent-plane trials. A negative TPD (liquid–liquid
  split) returns `UnsupportedCapability`; VLLE is not supported.
- **PH flash** locates the feed's bubble and dew temperatures at p relative to
  the option bracket. They split it into a liquid branch (`h ≤ h_L(T_b)`), a
  two-phase range and a vapor branch (`h ≥ h_V(T_d)`). A saturation temperature
  outside the bracket truncates or removes the corresponding branch and is
  reported as `A7-PH-SATURATION-OUTSIDE-BRACKET`. For example, at 1 kPa the
  bubble temperature of most mixtures lies below 273.16 K, so there is no
  liquid branch. A target that cannot be reached inside the bracket returns
  `OutsideValidityDomain`. The single-phase branches use Newton with exact c_p
  and the A5 criterion `|h − h_t|/max(|h_t|, 1 J/mol) ≤ 1e-10`. In the
  two-phase range, a mixture uses Illinois regula falsi on the split enthalpy
  (stability trials are skipped during the search, and the accepted
  temperature is flashed again with them). Its energy scale is
  `max(|h_t|, h(T_hi) − h(T_lo), 1 J/mol)` over the two-phase range. The split's
  composition tolerance and the temperature resolution limit the attainable
  enthalpy precision to a small fraction of that span. With a 1 J/mol floor, two-phase targets near
  h = 0 (for example 95 % vaporized at 1 atm with this reference) were
  unreachable. A single active component uses the lever rule at the
  correlation saturation temperature, retaining both endpoint phases. The
  scale used is reported (`energy_scale_J_per_mol` in the CLI).
- A split that converges to a vapor fraction within 1e-10 outside (0, 1),
  which is roundoff at p ≈ p_dew or p_bubble, is reported as the saturated
  single phase with `A7-SPLIT-AT-SATURATION-BOUNDARY`. Larger excursions
  remain `NumericalFailure`.

Every loop has a documented limit (500 composition, 200 temperature
iterations). Tests deliberately exhaust each limit. Invalid options, a
nonfinite target, an out-of-envelope state and an unbracketed point are
explicit errors.

The CLI accepts `"model": "NRTL-IG"` for `tp-flash`, `ph-flash`,
`bubble-pressure`, `dew-pressure`, `bubble-temperature` and `dew-temperature`;
see `examples/nrtl-*.json`. It rejects `kij`, `interaction_source` and `root`
for this model. Its provenance records the component and liquid/NRTL data
identifiers and versions. Bubble and dew outputs report the NRTL validity of
the saturated liquid.

## Independent evidence

`validation/fixtures/nrtl-v1.json` comes from `validation/julia/generate_nrtl.jl`
in the locked Julia 1.11.7 / Clapeyron 0.6.30 environment, which never calls
Rust:

| Group | States | Reference path |
|---|---|---|
| Methanol/water activity and excess | 30 (5 T incl. extrapolated 300 and 420 K × 6 x incl. both infinite-dilution limits) | Clapeyron `NRTL`, ForwardDiff for T derivatives; `excess(…, enthalpy)` cross-check |
| Directed ternary, every τ and α temperature term nonzero | 9 | Clapeyron `aspenNRTL` |
| Pure vapor pressure, derivatives, liquid h/s/c_p, density | 10 | `PolExpSat`, `PolExpLiquid`, `ReidIdeal` (NASA7 low-range c_p, anchored at 298.15 K/1 bar) |
| Bubble/dew pressure and temperature | 36 | modified-Raoult equations solved in 256-bit BigFloat with Clapeyron γ and p^sat |
| TP flash (3 liquid/vapor single phase, 5 two-phase) with phase h, s, c_p, v | 8 | BigFloat split; caloric assembled from the Clapeyron pieces above |
| PH recovery: 8 TP states and pure water/methanol at β = 0, 0.25, 0.75, 1 | 16 | targets from the reference states |

Clapeyron 0.6.30's own gamma–phi VLE cannot be used with correlation-based pure
fluids: `FluidCorrelation` defines no `a_res`. The generator therefore solves
the same equations itself, with Clapeyron supplying every property function.
The A7 Clapeyron check now runs alongside the PR/PR-BM check in `check.jl`
(push/PR and monthly CI). Ordinary workspace tests compare the committed
fixture offline.

Rust-only tests cover the ideal-solution limit; `g^E/RT = Σ x ln γ`;
Gibbs–Duhem; `ln γ_i = ∂(n g^E/RT)/∂n_i`; finite-difference T derivatives;
permutation of a directed ternary; directed parameter swaps; dataset
orientation; and every missing-data and invalid-parameter failure. They also
check caloric identities (`dh/dT = c_p`, `ds/dT = c_p/T`, Gibbs–Helmholtz),
pure-liquid/vapor Gibbs equality at p^sat, pressure independence of liquid
properties, and equilibrium identities with Gibbs decrease of the split.
Further tests cover TP→PH round trips, component reordering, inactive
components, pure lever-rule PH, a narrow-boiling mixture, LLE detection with a
hypothetical immiscible parameter set, and iteration exhaustion.

## Supported domain and remaining gates

- Implementation agreement is established for the methanol/water parameter set
  and the formulation above within 273.16–450 K, 100 Pa–500 kPa. It is not an
  experimental accuracy claim. Results outside 337.69–373.13 K are flagged as
  extrapolated NRTL parameters.
- Only water and methanol have pure-liquid correlations. Other components fail
  with `MissingParameter`, as does any pair other than methanol/water.
- The ideal-gas vapor and the neglected liquid volume restrict the package to
  low pressure. Vapor association, a nonideal vapor model, Poynting and φ^sat
  corrections are not implemented.
- At most one liquid phase. Finite liquid TPD trials are evidence, not proof,
  of liquid stability. Supercritical or noncondensable components (Henry's law)
  are not supported.
- Mixing of cubic and gamma–phi streams in a flowsheet needs A8 to select one
  property package per stream. The shared caloric reference only makes that
  possible.
- Hosted CI and a code review were completed on PR #75. The review findings and
  their resolution are recorded in [the A7 verification report](../validation/a7-checks.md).
