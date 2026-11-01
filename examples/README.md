# Reproducible examples

Two JSON requests are runnable with the Rust CLI:

```sh
cargo run -p thermo-cli --locked -- calculate --input examples/natural-gas.json
cargo run -p thermo-cli --locked -- calculate --input examples/binary-hydrocarbon.json
```

`natural-gas.json` uses 90/5/2/2/1 mol% methane/ethane/propane/nitrogen/CO₂ at 288.15 K and 5 MPa absolute, PR-BM, and explicitly zero kij. `binary-hydrocarbon.json` uses 70/30 mol% methane/ethane at 250 K and 3 MPa absolute with illustrative kij=0.02. Both select the largest homogeneous root. Their reference outputs are in the independent [fixture](../validation/cubic-reference.json).

Neither request establishes equilibrium phase behavior. The [calculation contract](../docs/calculation-contract.md) defines input/output and numerical limitations. Add a mixed/cooled/expanded/separated gas process only after TP/PH flash is qualified; no flowsheet project format is defined yet.
