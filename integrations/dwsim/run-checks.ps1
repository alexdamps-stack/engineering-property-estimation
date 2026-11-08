<#
.SYNOPSIS
  A12 foreign-interface checks: builds thermocore_bridge.dll, compiles the
  .NET Framework client with the csc.exe shipped in Windows, runs the CLR
  bridge tests and, with -DwsimPath, the real DWSIM separator check.

.PARAMETER DwsimPath
  Directory containing DWSIM.Automation.dll (installed or portable DWSIM
  9.0.x for Windows). Omit to run only the DWSIM-independent CLR checks.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File integrations/dwsim/run-checks.ps1 -DwsimPath "C:\Program Files\DWSIM"
#>
param(
    [string]$DwsimPath = "",
    [string]$OutDir = ""
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = Resolve-Path (Join-Path $here "..\..")
if ($OutDir -eq "") { $OutDir = Join-Path $root "target\a12-dwsim" }
New-Item -ItemType Directory -Force $OutDir | Out-Null

$csc = Join-Path $env:WINDIR "Microsoft.NET\Framework64\v4.0.30319\csc.exe"
if (-not (Test-Path $csc)) { throw "csc.exe for .NET Framework 4.x not found at $csc" }

Push-Location $root
try {
    Write-Host "== cargo build (release): thermo-bridge, thermo-cli"
    & cargo build --release -p thermo-bridge -p thermo-cli
    if ($LASTEXITCODE -ne 0) { throw "cargo build failed" }
} finally { Pop-Location }
$dll = Join-Path $root "target\release\thermocore_bridge.dll"
Copy-Item $dll $OutDir -Force
$bridgeDll = Join-Path $OutDir "thermocore_bridge.dll"

Write-Host "== csc: ThermoCore.Bridge.dll (.NET Framework 4.x, C# 5)"
$client = Join-Path $OutDir "ThermoCore.Bridge.dll"
& $csc /nologo /target:library /platform:x64 /optimize+ /warnaserror+ "/out:$client" `
    /reference:System.Web.Extensions.dll (Join-Path $here "src\ThermoCoreBridge.cs")
if ($LASTEXITCODE -ne 0) { throw "csc ThermoCore.Bridge failed" }

$tests = Join-Path $OutDir "BridgeClientTests.exe"
& $csc /nologo /target:exe /platform:x64 /optimize+ /warnaserror+ "/out:$tests" `
    "/reference:$client" /reference:System.Web.Extensions.dll (Join-Path $here "tests\BridgeClientTests.cs")
if ($LASTEXITCODE -ne 0) { throw "csc BridgeClientTests failed" }

# Independent process path through the existing CLI for the same state.
$cli = Join-Path $root "target\release\thermo-cli.exe"
$cliOut = & $cli (Join-Path $root "examples\tp-flash-natural-gas.json") | Out-String
if ($LASTEXITCODE -ne 0) { throw "thermo-cli failed" }
$cliBeta = [regex]::Match($cliOut, '"molar_vapor_fraction":\s*([0-9eE.+-]+)').Groups[1].Value
Write-Host "thermo-cli molar_vapor_fraction = $cliBeta"

Write-Host "== CLR bridge tests"
& $tests $bridgeDll (Join-Path $here "requests\natural-gas-separator-tp.json") $cliBeta
if ($LASTEXITCODE -ne 0) { throw "CLR bridge tests failed ($LASTEXITCODE)" }

if ($DwsimPath -eq "") {
    Write-Host "DWSIM check skipped: pass -DwsimPath <dir containing DWSIM.Automation.dll>"
    exit 0
}
$DwsimPath = (Resolve-Path $DwsimPath).Path
if (-not (Test-Path (Join-Path $DwsimPath "DWSIM.Automation.dll"))) { throw "DWSIM.Automation.dll not in $DwsimPath" }

Write-Host "== csc: DWSIM separator check against $DwsimPath"
$demo = Join-Path $OutDir "DwsimThermoCoreCheck.exe"
$refs = @("DWSIM.Automation.dll", "DWSIM.Interfaces.dll", "DWSIM.Thermodynamics.dll", "DWSIM.SharedClasses.dll",
          "DWSIM.GlobalSettings.dll", "DWSIM.UnitOperations.dll") | ForEach-Object { "/reference:" + (Join-Path $DwsimPath $_) }
& $csc /nologo /target:exe /platform:x64 /optimize+ /warnaserror+ "/out:$demo" "/reference:$client" `
    /reference:System.Web.Extensions.dll @refs (Join-Path $here "src\DwsimThermoCoreCheck.cs")
if ($LASTEXITCODE -ne 0) { throw "csc DwsimThermoCoreCheck failed" }
# DWSIM's own binding redirects apply to the host process.
Copy-Item (Join-Path $DwsimPath "DWSIM.exe.config") "$demo.config" -Force

$report = Join-Path $OutDir "dwsim-separator-check.json"
& $demo $DwsimPath $bridgeDll $report
if ($LASTEXITCODE -ne 0) { throw "DWSIM separator check failed ($LASTEXITCODE)" }
Write-Host "Report: $report"
