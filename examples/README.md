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
remain unassessed. The provisional process file is not a CLI input format.
