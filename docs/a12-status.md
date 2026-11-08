# A12 status: external integration (DWSIM)

Issues: #63 (work package), #64 (define the use case and interface), #65 (demonstrate it).
Excel stays deferred (development plan, section 3).

## 1. Scope decision (#64)

**Use case.** *Independent ThermoCore equilibrium check of a DWSIM separator.* A
process engineer has a DWSIM flowsheet with a feed stream and a separator vessel.
The feed state, the compound identities and the Peng–Robinson interaction
parameters DWSIM actually uses are taken from the DWSIM objects and evaluated by
the ThermoCore cubic TP flash. ThermoCore returns the phase split, phase
properties, convergence evidence and full provenance, and the host reports it
next to DWSIM's own split.

**Interface.** Two layers, each versioned and independent of DWSIM:

| Layer | Artifact | Contract |
|---|---|---|
| Rust boundary | `crates/thermo-bridge` → `thermocore_bridge.dll` (cdylib) | C ABI v1, four functions (`include/thermocore_bridge.h`); UTF-8 JSON schemas `thermocore.bridge.{request,response,info}.v1` |
| .NET client | `integrations/dwsim/src/ThermoCoreBridge.cs` → `ThermoCore.Bridge.dll` | .NET Framework 4.x, C# 5 (DWSIM on Windows needs .NET Framework 4.6.2+); P/Invoke, ABI check, typed `ThermoCoreException` |
| DWSIM host | `integrations/dwsim/src/DwsimThermoCoreCheck.cs` | DWSIM `Automation3` (headless), DWSIM 9.0.x |

Why this boundary: DWSIM's own extension points (CAPE-OPEN Thermo 1.1 COM
servers, `PropertyPackage` subclasses loaded by DWSIM, `IExternalUnitOperation`
plugins with editor and drawing code) each need a large, GUI-coupled surface and
would make ThermoCore results replace DWSIM's internal ones without the
provenance survival this project requires. A JSON-over-C ABI is the smallest
stable boundary that any of those later hosts (and a later Excel add-in) can
reuse unchanged; the DWSIM Automation API is DWSIM's supported way to drive
real flowsheets from .NET without the GUI.

**Bridge v1 scope.** `PR` and `PR-BM`, `tp-flash` and `ph-flash`, at most two
fluid phases, through `property_packages::CubicPackage` (200–600 K search
interval). GERG-2008, NRTL, saturation points and flowsheets are not exposed in v1.
No solver or model code was changed.

## 2. Contract (bridge v1)

C functions (all reentrant, no global mutable state):

```c
uint32_t tc_bridge_abi_version(void);           /* == 1 */
char *tc_bridge_info(void);                       /* info.v1 document */
char *tc_bridge_calculate(const char *request);   /* request.v1 -> response.v1 */
void tc_bridge_free(char *document);              /* frees any returned document */
```

Request (`integrations/dwsim/requests/natural-gas-separator-tp.json`):

```json
{"schema": "thermocore.bridge.request.v1", "request_id": "…", "calculation": "tp-flash",
 "package": {"model": "PR-BM", "components": ["methane", "…"], "kij": [[…]], "interaction_source": "…"},
 "state": {"mole_fractions": […], "temperature_K": 220.0, "pressure_Pa": 1500000.0}}
```

- **SI only**, with the unit in every field name; pressure is absolute; enthalpy
  (`enthalpy_j_per_mol`, PH input) is on the reported caloric reference.
- **No physics defaults.** Unknown fields are rejected (`pressure_bar` fails).
  Component order, a complete symmetric `kij` and a non-empty interaction source
  are mandatory; mole fractions must sum to one within 1e-10.
- **Response envelope**: `schema`, `abi_version`, `bridge_version`, echoed
  `request_id`, `status` `ok` with `result`, or `error` with
  `{kind, message}` and no result. Kinds: `invalid-request`, `invalid-input`,
  `missing-parameter`, `unsupported-capability`, `outside-validity-domain`,
  `non-convergence`, `numerical-failure`, `internal-panic` (panics are caught at
  the ABI and never unwind into the host). A null or non-UTF-8 request yields an
  `invalid-request` document, never a null pointer.
- **Result**: `temperature_K`, `pressure_Pa`, mixture `enthalpy_J_per_mol` and
  `entropy_J_per_mol_K`, `molar_vapor_fraction` (null for a supercritical single
  phase), and per phase: label, molar phase fraction, mole fractions, molar mass,
  Z, molar volume, mass density, h, s, cp, cv, speed of sound, ln φ and validity
  status. `convergence` holds iterations and the fugacity, component-balance and
  energy residuals. `provenance` holds model id/implementation/formulation (the
  formulation contains the interaction source and the kij matrix), dataset
  id/version/source and the caloric reference id/definition. Solver
  diagnostics are passed through with code and severity.
- **Component mapping** is explicit: `info` lists the catalogue with CAS numbers;
  the DWSIM host maps DWSIM compounds by CAS and stops on any unmapped compound.
- **Numbers.** The Rust side writes shortest round-trip decimals and parses with
  serde_json `float_roundtrip`, so values cross the boundary exactly
  (`tests/c_abi.rs` compares bit-for-bit with a direct `tp_flash`). The .NET
  client writes doubles with `G17`. .NET Framework number parsing is not
  correctly rounded, so host-side values may differ by about one ulp.

## 3. Acceptance evidence (#65)

Reproduce (Windows; Rust toolchain; `csc.exe` of .NET Framework 4.x, present on
Windows 10/11):

```powershell
powershell -ExecutionPolicy Bypass -File integrations/dwsim/run-checks.ps1 -DwsimPath "<dir with DWSIM.Automation.dll>"
```

Without `-DwsimPath` only the DWSIM-independent CLR checks run.

**DWSIM used:** official release
[v9.0.4](https://github.com/DanWBR/dwsim/releases/tag/v9.0.4),
`DWSIM_v904_win64_portable.7z`, SHA-256
`3a2afccf887018882c4ca3e6779d39e24c1040ef3bf6606c48ea179964258f87`; reported
version "DWSIM version 9.0.4.0 (28.07.2025 13:19:16)". The Automation3 and
interface members used were read by reflection from the shipped assemblies.

Run log: `integrations/dwsim/evidence/run-checks-dwsim-9.0.4.txt`; report with
the full request, ThermoCore response and DWSIM values:
`integrations/dwsim/evidence/dwsim-9.0.4-separator-check.json`. The script
also saves the solved DWSIM flowsheet as `target/a12-dwsim/dwsim-separator-check.dwxmz`.

| Check | Result |
|---|---|
| CLR (.NET Framework 4.x) bridge client, 19 checks: ABI, info/CAS, SI echo, provenance, residual, typed errors (150 K, unknown component, non-SI field, malformed JSON, NaN), 16 concurrent calls | all pass |
| Bridge result vs `thermo-cli` (separate process, same state): vapor fraction | 0.98749222958091643 vs 0.98749222958091631 (1 ulp, .NET parsing) |
| DWSIM flowsheet FEED → SEPARATOR → VAPOR/LIQUID solved by DWSIM 9.0.4 "Peng-Robinson (PR)", 220 K, 1.5 MPa, 1 mol/s | solved, mole balance closed |
| kij taken from the DWSIM package (e.g. CH4–CO2 0.0793, C2H6–CO2 0.1326) and carried into ThermoCore provenance | pass |
| Vapor fraction DWSIM / ThermoCore PR | 0.989676 / 0.989579 (Δ = 9.7e-5) |
| max \|Δy\| / max \|Δx\| | 6.9e-5 / 1.05e-4 |
| Vapor Z DWSIM / ThermoCore | 0.89334 / 0.88984 (Δ = 3.5e-3) |
| ThermoCore ln-fugacity residual | ≤ 1e-10 |
| 150 K DWSIM state | typed `unsupported-capability` |
| DWSIM n-butane stream | rejected before the call (no ThermoCore component for CAS 106-97-8) |

The gates (|Δβ| ≤ 0.02, |Δx|,|Δy| ≤ 0.02, |ΔZ| ≤ 0.01) were fixed before the
first DWSIM run. They detect integration faults (units, component order, kij
mapping); DWSIM is a cross-check, not a qualified reference. The remaining
differences are expected from different critical constants and model details.
The vapor-Z difference is larger than the split difference and was not
investigated further; it does not affect the integration verdict.

Rust checks: `cargo test -p thermo-bridge` (4 unit tests, 3 C-ABI tests),
`cargo clippy -p thermo-bridge --all-targets` clean, `cargo fmt --check`, and
`cargo test --workspace`.

## 4. Limits and follow-ups

- Not a DWSIM property package: DWSIM still calculates its own streams with its
  own package. Registering ThermoCore inside DWSIM (CAPE-OPEN Thermo 1.1 or a
  DWSIM `PropertyPackage`) and use from the DWSIM GUI or an IronPython script are
  follow-ups that can reuse bridge v1 unchanged.
- DWSIM is not available in CI; the DWSIM step runs only with `-DwsimPath`. The
  Rust and CLR steps need no DWSIM.
- The kij read uses a non-public DWSIM field (`m_pr.InteractionParameters`) and
  fails explicitly if a later DWSIM version changes it.
- `float_roundtrip` is a workspace-unified serde_json feature: builds that
  include `thermo-bridge` parse JSON floats exactly in every crate. Inputs then
  differ from builds without it by at most one ulp; no workspace test changed.
- Licensing: DWSIM is GPL-3.0. Neither `thermocore_bridge.dll` nor
  `ThermoCore.Bridge.dll` references DWSIM; only the check program compiles
  against a DWSIM installation provided by the user, and it is not part of the
  Windows package.
