<#
.SYNOPSIS
  Updates the installed Ubuntu/WSL desktop prototype, or checks its source version.
#>
param(
    [switch]$Check,
    [string]$Distribution = 'Ubuntu'
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$arguments = @('-d', $Distribution, '--cd', $root, '--', 'python3', 'scripts/update-local-app.py')
if ($Check) { $arguments += '--check' }
& wsl.exe @arguments
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Update shortcuts only after the candidate has passed its self-test and been published.
$gui = (Get-Command wslg.exe -ErrorAction Stop).Source
$linuxHome = (& wsl.exe -d $Distribution -- printenv HOME).Trim()
if ($LASTEXITCODE -ne 0 -or -not $linuxHome.StartsWith('/')) {
    throw 'Cannot determine the installed app path in WSL.'
}
$documents = [Environment]::GetFolderPath('MyDocuments')
if ($documents -notmatch '^([A-Za-z]):[\\/]') { throw 'Expected a local Windows Documents folder.' }
$linuxDocuments = '/mnt/' + $documents.Substring(0, 1).ToLower() + $documents.Substring(2).Replace('\', '/')
$app = $linuxHome + '/.local/share/property-estimation-prototype/current/thermo-desktop'
$launchArguments = '-d "' + $Distribution + '" --cd "' + $linuxDocuments + '" -- "' + $app + '" --demo'
$shell = New-Object -ComObject WScript.Shell
foreach ($folder in @([Environment]::GetFolderPath('Desktop'), [Environment]::GetFolderPath('Programs'))) {
    $shortcutPath = Join-Path $folder 'Property Estimation Prototype.lnk'
    if ($Check -and -not (Test-Path -LiteralPath $shortcutPath)) {
        Write-Output 'An app shortcut is missing; run the updater.'
        exit 1
    }
    $shortcut = $shell.CreateShortcut($shortcutPath)
    if ($Check) {
        if ($shortcut.TargetPath -ne $gui -or $shortcut.Arguments -ne $launchArguments) {
            Write-Output 'An app shortcut is stale; run the updater.'
            exit 1
        }
        continue
    }
    $shortcut.TargetPath = $gui
    $shortcut.Arguments = $launchArguments
    $shortcut.WorkingDirectory = $documents
    $shortcut.Description = 'Property Estimation desktop prototype (Ubuntu / WSL GUI)'
    $shortcut.IconLocation = 'C:\Windows\System32\shell32.dll,21'
    $shortcut.Save()
}
if ($Check) {
    Write-Output 'Desktop and Start menu shortcuts are current.'
    exit 0
}
Write-Output 'Desktop and Start menu shortcuts point to the verified app. Reopen the app to use the update.'
