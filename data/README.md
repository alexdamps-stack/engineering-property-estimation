# Component and interaction data

The provisional [source extract](source-extracts/components.provisional.json) contains eight CAS-matched component constant records from pinned CoolProp 7.2.0 data and gas-phase NASA7 records from pinned Cantera 3.2.0 data. Upstream license notices are retained alongside it.

Before adding data, record component identifiers, units, source and version, usage terms, applicable model, validity limits, and missing parameters. Initial planned components are methane, ethane, propane, nitrogen, carbon dioxide, oxygen, water, and methanol.

`thermo-data` now checks and loads this exact pinned extract and evaluates NASA7 using the reconciled A1 caloric convention. The historical filename/schema is retained for source reproduction. Source transcription and permission evidence do not establish experimental accuracy. See [A1–A3 status](../docs/a1-a3-status.md) and [reference qualification](../docs/reference-qualification.md) for formulations, executed comparisons and remaining gates. No fitted binary interactions or implicit numerical defaults are supplied.

A7 adds [`liquid-nrtl.provisional.json`](source-extracts/liquid-nrtl.provisional.json): CoolProp 7.2 vapor-pressure and saturated-liquid density ancillaries for water and methanol, and the ChemSep/DECHEMA 1 atm methanol/water NRTL row (via pinned DWSIM and the MIT-licensed `thermo` redistribution, whose [notice](source-extracts/thermo-LICENSE.txt) is retained). It is produced only by the explicit [`prepare-nrtl-source-data.py`](../validation/prepare-nrtl-source-data.py). Other NRTL pairs and pure-liquid records are deliberately absent and fail as missing data. See [A7 status](../docs/a7-status.md).
