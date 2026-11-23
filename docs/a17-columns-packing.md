# A17 columns: packed sections (ideal theoretical-stage path)

Part of the columns group (#186, #187, #190); see
[a17-implementation-brief.md](a17-implementation-brief.md) § "Columns and
reactors as packed beds". Code: `crates/flowsheet-core/src/units/columns/packing.rs`
(`flowsheet_core::units::columns::packing`). Status: foundation only. No column
kind is calculated through it yet, and no record, form or palette change exists.

## Inputs

| Field | Unit | Validation |
|---|---|---|
| HETP | m | Required, finite, > 0 |
| Basis: packed height | m | Required, finite, > 0; with a rounding rule |
| Basis: theoretical stages | – | 1 to `MAX_PACKED_STAGES` (500) |

The two bases are mutually exclusive (`PackedSectionBasis`). Missing or
invalid values are reported through `PackedSectionSpec::issues` as ordinary
`SpecIssue`s (`Missing`, `NotFinite`, `NotPositive`, `CountOutOfRange`).

## Stage derivation and rounding

Height basis: \(N_\text{raw} = Z/\text{HETP}\) (Kister, *Distillation Design*,
1992, ch. 11; Seader, Henley and Roper, *Separation Process Principles*, 3rd
ed., §6.8). The equilibrium-stage solver needs a whole number of stages:

1. If \(N_\text{raw}\) lies within a relative 1e-9 of a whole number, that
   number is used. Floating-point quotients such as 1.2/0.4 =
   2.9999999999999996 therefore give 3 under every rule.
2. Otherwise the selected `StageRounding` applies: **Nearest** (default, halves
   up, deviation at most half an HETP), **Down** (conservative rating; never
   credits more separation than the bed height provides) or **Up**.
3. A result of 0 stages, or more than 500, fails with an explicit message.
   Nothing is clamped.

The resolved `PackedStages` reports the specified height, the raw quotient,
the modelled height \(N\cdot\text{HETP}\), its deviation from \(Z\), and the
effective HETP \(Z/N\). The stage basis derives \(Z = N\cdot\text{HETP}\)
exactly (as in RadFrac or DWSIM, where the stage count is primary).

Fractional stages are not represented by a partial-efficiency stage. That
would be a different approximation and would mix HETP with a Murphree
efficiency.

## Pairing with the equilibrium-stage solver

- `PackedSectionSpec::resolve()` → `PackedStages { stages, hetp_m, … }` is the
  interface. A packed section contributes `stages` consecutive ideal stages
  to the stage stack with **efficiency 1** (`stage_efficiency()`). HETP already
  includes non-ideal contacting, so a Murphree or overall efficiency must not
  also be applied to packed stages. Trayed sections keep their own
  efficiencies.
- `stage_at_depth(depth)` maps an elevation below the top of the bed to a
  zero-based stage of the section. It uses the effective HETP so that the full
  bed maps onto the modelled stages. Feeds and draws given by elevation can be
  located with it.
- Condenser and reboiler are separate equilibrium stages and are not part of
  a packed section's count.

## Pairing with a later rate-based model

The rate-based batch adds a calculation-mode choice alongside this path.
It reuses the same packed height, adds column diameter, packing data (specific
area, void fraction, nominal size, packing factor), the mass-transfer
correlation (Onda 1968, Billet–Schultes 1999, Bravo–Rocha–Fair), film model and
segment count, and reports hydraulics (pressure drop, flooding fraction). Rate-based
segments are not theoretical stages, and no HETP applies there. An HETP
correlation would be a further `hetp_m` source in this module. It belongs to
that batch because it needs the same packing data. None of these correlations
are implemented here.

## Integration needs

- A section choice (trayed/packed) and calculation mode in the column specs
  (`columns.rs`), records (`records/columns.rs`, fields behind format
  version 4), and the column spec form.
- The stage solver consumes `PackedStages`, with its stage efficiency fixed at
  1 for packed stages.
- A `docs/a17-status.md` Columns row once a column kind is calculated with
  packed sections.

## Placement in the stage cascade

`stages.rs::PlacedPackedSection` converts the section-local zero-based
`stage_at_depth` to one-based column stages. It also checks condenser and
reboiler exclusion and efficiency 1 against a `StageCascade`
(see [a17-columns-stages.md](a17-columns-stages.md#packed-sections-joint-contract-with-packingrs)).
