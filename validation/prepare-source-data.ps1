param([Parameter(Mandatory=$true)][string]$OutputDirectory)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
# Explicit review-only generator: never run this automatically to repair tests.
$cp = '98b3523d5daa98454618d381d2ae53f7471d216b'
$ct = '4a8358eb80cfeb50474386b5f9ec0b3a83519889'
$sources = [ordered]@{
    constants = "https://raw.githubusercontent.com/CoolProp/CoolProp/$cp/dev/cubics/all_cubic_fluids.json"
    nasa7 = "https://raw.githubusercontent.com/Cantera/cantera/$ct/data/nasa_gas.yaml"
    coolprop_license = "https://raw.githubusercontent.com/CoolProp/CoolProp/$cp/LICENSE"
    cantera_license = "https://raw.githubusercontent.com/Cantera/cantera/$ct/License.txt"
}
$downloaded = @{}
$provenance = @()
foreach ($key in $sources.Keys) {
    $response = Invoke-WebRequest -UseBasicParsing -Uri $sources[$key]
    $body = [string]$response.Content
    $downloaded[$key] = $body
    $hash = [System.Security.Cryptography.SHA256]::Create()
    try { $digest = [BitConverter]::ToString($hash.ComputeHash([Text.Encoding]::UTF8.GetBytes($body))).Replace('-','').ToLowerInvariant() } finally { $hash.Dispose() }
    $provenance += [ordered]@{id=$key;url=$sources[$key];sha256_utf8_response_text=$digest}
}
$rows = @(
    @('methane','74-82-8','CH4'), @('ethane','74-84-0','C2H6'),
    @('propane','74-98-6','C3H8'), @('nitrogen','7727-37-9','N2'),
    @('carbon-dioxide','124-38-9','CO2'), @('oxygen','7782-44-7','O2'),
    @('water','7732-18-5','H2O'), @('methanol','67-56-1','CH3OH')
)
$all = ConvertFrom-Json $downloaded.constants
$components = @()
foreach ($row in $rows) {
    $matches = @($all | Where-Object CAS -eq $row[1])
    if ($matches.Count -ne 1) { throw "Expected unique CAS $($row[1])" }
    $fluid = $matches[0]
    if ($fluid.Tc_units -ne 'K' -or $fluid.pc_units -ne 'Pa' -or $fluid.molemass_units -ne 'kg/mol') { throw 'Unexpected source units' }
    $pattern = '(?ms)^- name: ' + [regex]::Escape($row[2]) + '\r?\n.*?(?=^- name: |\z)'
    $blocks = [regex]::Matches($downloaded.nasa7, $pattern)
    if ($blocks.Count -ne 1) { throw "Expected unique gas species $($row[2])" }
    $block = $blocks[0].Value
    if ($block -notmatch 'model: NASA7') { throw 'Unexpected correlation' }
    $ranges = [regex]::Match($block, 'temperature-ranges: \[([^\]]+)\]').Groups[1].Value
    $range = @($ranges.Split(',') | ForEach-Object {[double]::Parse($_.Trim(), [Globalization.CultureInfo]::InvariantCulture)})
    $coefficients = @()
    foreach ($match in [regex]::Matches($block, '(?s)    - \[([^\]]+)\]')) {
        $coefficients += ,@($match.Groups[1].Value.Split(',') | ForEach-Object {[double]::Parse($_.Trim(), [Globalization.CultureInfo]::InvariantCulture)})
    }
    if ($range.Count -ne 3 -or $coefficients.Count -ne 2 -or $coefficients[0].Count -ne 7 -or $coefficients[1].Count -ne 7) { throw 'Unexpected NASA7 shape' }
    $note = [regex]::Match($block, '(?m)^    note: (.+)').Groups[1].Value.Trim()
    $components += [ordered]@{
        id=$row[0];cas=$row[1];formula=$row[2];source_name=$fluid.name
        constants=[ordered]@{source='constants';critical_temperature_K=$fluid.Tc;critical_pressure_Pa=$fluid.pc;acentric_factor=$fluid.acentric;molar_mass_kg_per_mol=$fluid.molemass}
        ideal_gas=[ordered]@{source='nasa7';species=$row[2];correlation='NASA7';temperature_breaks_K=$range;coefficients=$coefficients;source_note=$note;coefficient_order='a1,a2,a3,a4,a5,a6,a7';reference_state_status='source integration constants retained; engine anchoring pending'}
    }
}
$document = [ordered]@{
    schema='a2-source-extract-provisional-v1';version='2026-11-07.1';status='source-transcription-only; A1 reconciliation and model qualification pending'
    sources=$provenance;components=$components
}
[IO.Directory]::CreateDirectory([IO.Path]::GetFullPath($OutputDirectory)) | Out-Null
$utf8 = New-Object Text.UTF8Encoding($false)
[IO.File]::WriteAllText((Join-Path $OutputDirectory 'components.provisional.json'), ($document | ConvertTo-Json -Depth 14) + "`n", $utf8)
foreach ($pair in @(@('coolprop_license','CoolProp-LICENSE.txt'),@('cantera_license','Cantera-License.txt'))) {
    [IO.File]::WriteAllText((Join-Path $OutputDirectory $pair[1]), $downloaded[$pair[0]], $utf8)
}
Write-Output "Prepared $($components.Count) component source transcriptions; no expected thermodynamic values generated."
