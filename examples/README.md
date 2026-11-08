# Reproducible examples

The [case inputs](cases.provisional.json) define a dry natural-gas mixture and a methane/propane binary. The [process proposal](process.provisional.json) defines two feeds, mixing, cooling, pressure reduction and separation, with explicit mole/flow bases, SI units and assumptions.

These provisional files are engineering design inputs, with deliberately explicit zero-interaction approximations. Their null process-phase expectations remain intentional gaps. See [reference qualification](../docs/reference-qualification.md).

The executable homogeneous CLI cases are [natural gas](homogeneous-natural-gas.json)
and [methane/propane](homogeneous-binary.json):

```sh
cargo run -p thermo-cli --locked -- examples/homogeneous-natural-gas.json
cargo run -p thermo-cli --locked -- examples/homogeneous-binary.json
```

They produce homogeneous properties with explicit SI quantities and provenance.
Root selection compares fixed-composition roots; phase stability and flashes
remain unassessed in that mode. Two-phase TP flashes with stability evidence use
`"calculation": "tp-flash"`, see [binary](tp-flash-binary.json) and
[mixed natural-gas separator state](tp-flash-natural-gas.json); see [A4 status](../docs/a4-status.md).
The provisional process file is not a CLI input format.

PH calculations use `"calculation": "ph-flash"` with `enthalpy_j_per_mol`
instead of `temperature_k`. The target is molar enthalpy on the reported
caloric reference; temperature is an output. Root requests are rejected.
Run the [natural-gas example](ph-flash-natural-gas.json):

```sh
cargo run -p thermo-cli --locked -- examples/ph-flash-natural-gas.json
```

TP output now includes phase-weighted `enthalpy_J_per_mol`. To calculate an
ideal adiabatic valve, use this inlet enthalpy as the PH target at the outlet
pressure, keeping composition, model, interactions and caloric reference fixed.
PH output includes the energy residual, phase properties and equilibrium
evidence. See [A5 status](../docs/a5-status.md) for supported cases and limits.

The A7 gamma–phi package uses `"model": "NRTL-IG"` (NRTL liquid with an
explicitly selected ideal-gas vapor; embedded methanol/water data only, no
`kij`). It supports `tp-flash`, `ph-flash`, `bubble-pressure`/`dew-pressure`
(input `temperature_k`) and `bubble-temperature`/`dew-temperature` (input
`pressure_pa`). See the [TP example](nrtl-tp-flash-methanol-water.json), the
[bubble-temperature example](nrtl-bubble-temperature-methanol-water.json) and
[A7 status](../docs/a7-status.md).
