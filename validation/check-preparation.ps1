$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$root = Split-Path $PSScriptRoot -Parent
function Read-Json([string]$path) { Get-Content -Raw -LiteralPath (Join-Path $root $path) | ConvertFrom-Json }
function Require([bool]$condition, [string]$message) { if (-not $condition) { throw $message } }
function Check-Composition($ids, $values, $known) {
    Require ($ids.Count -eq $values.Count) 'Composition length mismatch'
    Require (@($ids | Select-Object -Unique).Count -eq $ids.Count) 'Duplicate component'
    $sum = 0.0
    foreach ($i in 0..($values.Count - 1)) {
        Require ($ids[$i] -in $known) 'Unknown component'
        Require ([double]::IsNaN($values[$i]) -eq $false -and [double]::IsInfinity($values[$i]) -eq $false -and $values[$i] -ge 0) 'Invalid composition'
        $sum += $values[$i]
    }
    Require ([Math]::Abs($sum - 1) -le 1e-12) 'Unnormalized composition'
}
$data = Read-Json 'data/source-extracts/components.provisional.json'
Require ($data.components.Count -eq 8) 'Expected eight initial components'
$ids = @($data.components.id)
Require (@($ids | Select-Object -Unique).Count -eq 8) 'Duplicate identities'
Require (@($data.components.cas | Select-Object -Unique).Count -eq 8) 'Duplicate CAS'
foreach ($source in $data.sources) {
    Require ($source.url -match '/[a-f0-9]{40}/') 'Source URL not pinned to commit'
    Require ($source.sha256_utf8_response_text -match '^[a-f0-9]{64}$') 'Invalid source digest'
}
foreach ($component in $data.components) {
    $c = $component.constants
    foreach ($v in @($c.critical_temperature_K,$c.critical_pressure_Pa,$c.molar_mass_kg_per_mol)) {
        Require ($v -gt 0 -and -not [double]::IsInfinity($v)) 'Invalid positive constant'
    }
    Require (-not [double]::IsNaN($c.acentric_factor) -and -not [double]::IsInfinity($c.acentric_factor)) 'Invalid acentric factor'
    $ig = $component.ideal_gas
    Require ($ig.correlation -eq 'NASA7') 'Unexpected ideal model'
    Require ($ig.temperature_breaks_K.Count -eq 3) 'Invalid temperature intervals'
    Require ($ig.temperature_breaks_K[0] -eq 200 -and $ig.temperature_breaks_K[1] -eq 1000 -and $ig.temperature_breaks_K[2] -eq 6000) 'Unexpected source intervals'
    Require ($ig.coefficients.Count -eq 2) 'Invalid polynomial count'
    foreach ($coefficients in $ig.coefficients) {
        Require ($coefficients.Count -eq 7) 'Invalid NASA7 coefficient count'
        foreach ($v in $coefficients) { Require (-not [double]::IsNaN($v) -and -not [double]::IsInfinity($v)) 'Nonfinite coefficient' }
    }
}
$cases = Read-Json 'examples/cases.provisional.json'
$count = 0
foreach ($case in $cases.cases) {
    Check-Composition $case.components $case.mole_fractions $ids
    Require ($case.temperature_K -ge 200 -and $case.temperature_K -le 6000 -and $case.pressure_Pa -gt 0 -and $case.flow_mol_per_s -gt 0) 'Invalid case state'
    Require ($null -eq $case.expected_phase) 'Unqualified phase claim'
    $count += $case.sweep_temperature_K.Count * $case.sweep_pressure_Pa.Count
}
Require ($count -eq 18) 'Unexpected comparison grid count'
Require ($null -eq $cases.independent_expected_outputs) 'Unexpected unqualified numerical fixture'
$process = Read-Json 'examples/process.provisional.json'
Require ('water' -notin $process.components -and 'methanol' -notin $process.components) 'Wet demonstration feed'
foreach ($feed in $process.feeds) { Check-Composition $process.components $feed.mole_fractions $ids }
foreach ($phase in $process.expected_phase_by_stream.PSObject.Properties) { Require ($null -eq $phase.Value) 'Unqualified process phase claim' }
Require ($null -eq $process.independent_expected_outputs) 'Unexpected process expectations'
$tol = Read-Json 'validation/tolerances.preparation-history.json'
foreach ($q in $tol.quantities) { Require ($q.atol -gt 0 -and $q.rtol -gt 0 -and $q.unit.Length -gt 0 -and $q.minimum_states_per_model -gt 0) 'Invalid tolerance proposal' }
# Exercise the checker on actual malformed inputs, with no model implementation.
$rejections = 0
foreach ($bad in @(@(0.7,0.7),@(1.1,-0.1),@([double]::NaN,0.0))) {
    try { Check-Composition @('methane','propane') $bad $ids } catch { $rejections++ }
}
Require ($rejections -eq 3) 'Composition checks did not reject malformed input'
Write-Output "PASS: eight source records, NASA7 shapes/ranges, source pins, two cases/$count grid states per model, dry process, pending expectations, quantity tolerances, three invalid-composition rejections."
Write-Output 'These are artifact integrity checks, not thermodynamic verification or A2 acceptance.'
