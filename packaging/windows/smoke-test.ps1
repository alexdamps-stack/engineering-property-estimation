<#
.SYNOPSIS
  Installs the built setup silently for the current user, runs the installed
  executable's headless acceptance self-test, checks the Start-menu entry,
  and uninstalls again. Exits non-zero on any failure.
#>
param(
    [Parameter(Mandatory = $true)][string]$Setup,
    [string]$ReportDir = ""
)

$ErrorActionPreference = 'Stop'
if (-not $ReportDir) { $ReportDir = Split-Path -Parent $Setup }
$dir = Join-Path $env:TEMP "pe-flowsheet-smoke-$PID"

function Invoke-Checked([string]$File, [string[]]$Arguments) {
    $p = Start-Process -FilePath $File -ArgumentList $Arguments -Wait -PassThru
    if ($p.ExitCode -ne 0) { throw "$File $($Arguments -join ' ') exited with $($p.ExitCode)" }
}

Invoke-Checked $Setup @('/VERYSILENT', '/SUPPRESSMSGBOXES', '/NORESTART', '/CURRENTUSER', "/DIR=`"$dir`"", "/LOG=`"$ReportDir\install.log`"")
$exe = Join-Path $dir 'thermo-desktop.exe'
if (-not (Test-Path $exe)) { throw "installed executable missing: $exe" }
foreach ($f in 'README.txt', 'notices\NOTICE.txt', 'notices\THIRD-PARTY-CRATES.txt', 'notices\CoolProp-LICENSE.txt') {
    if (-not (Test-Path (Join-Path $dir $f))) { throw "installed file missing: $f" }
}
$menu = Join-Path ([Environment]::GetFolderPath('Programs')) 'Property Estimation Flowsheet'
if (-not (Get-ChildItem -Recurse -Filter '*.lnk' $menu -ErrorAction SilentlyContinue)) {
    throw "Start-menu shortcuts missing in $menu"
}

$version = Join-Path $ReportDir 'installed-version.txt'
Invoke-Checked $exe @('--version', "`"$version`"")
Get-Content $version

$report = Join-Path $ReportDir 'installed-self-test.json'
Invoke-Checked $exe @('--self-test', "`"$report`"")
$r = Get-Content $report -Raw | ConvertFrom-Json
if ($r.status -ne 'passed') { throw "self-test status $($r.status)" }
$r.report.steps | ForEach-Object { "  $_" }

$uninstaller = Get-ChildItem $dir -Filter 'unins*.exe' | Select-Object -First 1
Invoke-Checked $uninstaller.FullName @('/VERYSILENT', '/SUPPRESSMSGBOXES', '/NORESTART')
# The uninstaller relaunches itself from a temporary copy; wait for it.
for ($i = 0; $i -lt 60 -and (Test-Path $exe); $i++) { Start-Sleep -Milliseconds 500 }
if (Test-Path $exe) { throw "uninstall left $exe behind" }
'Installer smoke test passed.'
