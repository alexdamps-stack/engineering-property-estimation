# Component and interaction data

`components.json` supplies eight selected constant records from CoolProp 8.0.0. Records include IDs, CAS numbers, EOS reducing temperature/pressure used as PR critical inputs, molar mass, acentric factor, source-fluid names, and source JSON hashes. See [the notice](NOTICE.md) for attribution and [the calculation contract](../docs/calculation-contract.md) for the version and usage limits.

Before adding data, record component identifiers, units, source and version, usage terms, applicable model, validity limits, and missing parameters. Initial planned components are methane, ethane, propane, nitrogen, carbon dioxide, oxygen, water, and methanol.

The initial constant schema is consumed by `thermo-data`; it supplies no binary interaction parameters or ideal caloric coefficients. Matrices must be supplied explicitly in each request. The sample matrices are illustrative choices, not fitted data. These constants have not been qualified as accurate property packages across all fluid/phase domains, especially polar liquids.
