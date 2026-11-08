# ThermoCore ↔ DWSIM (A12)

Status, contract and acceptance evidence: [docs/a12-status.md](../../docs/a12-status.md).

| Path | Content |
|---|---|
| `src/ThermoCoreBridge.cs` | .NET Framework 4.x client of the bridge C ABI v1 (no DWSIM dependency) |
| `src/DwsimThermoCoreCheck.cs` | Headless DWSIM 9 separator check through `DWSIM.Automation.Automation3` |
| `tests/BridgeClientTests.cs` | CLR foreign-interface tests of the C ABI |
| `requests/` | v1 request documents |
| `evidence/` | Run log and report from DWSIM 9.0.4 |
| `run-checks.ps1` | Builds everything and runs the checks |

```powershell
# CLR checks only
powershell -ExecutionPolicy Bypass -File integrations/dwsim/run-checks.ps1
# plus the real DWSIM flowsheet (installed or portable DWSIM 9.0.x for Windows)
powershell -ExecutionPolicy Bypass -File integrations/dwsim/run-checks.ps1 -DwsimPath "C:\Program Files\DWSIM"
```

Outputs go to `target/a12-dwsim/`. The C# sources target the C# 5 compiler in
`%WINDIR%\Microsoft.NET\Framework64\v4.0.30319\csc.exe`, so no .NET SDK is needed.

Using the client from another .NET host:

```csharp
ThermoCoreBridge.Load(@"path\to\thermocore_bridge.dll");   // checks ABI == 1
var result = ThermoCoreBridge.Calculate(ThermoCoreBridge.TpFlashRequest(
    "id", "PR", components, kij, "source of kij", moleFractions, 220.0 /* K */, 1.5e6 /* Pa */));
// throws ThermoCoreException { Kind = "unsupported-capability", ... } on failure
```
