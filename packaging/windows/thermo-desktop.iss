; Inno Setup 6 script for the Property Estimation Flowsheet desktop editor.
; Build through packaging/windows/build.ps1, which stages the files and passes
; AppVersion, StageDir and OutputDir:
;   ISCC.exe /DAppVersion=0.1.0 /DStageDir=... /DOutputDir=... thermo-desktop.iss
;
; Per-user installation by default (no administrator rights needed); the
; user may choose an all-users installation in the dialog or with /ALLUSERS.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#ifndef StageDir
  #error StageDir must be defined
#endif
#ifndef OutputDir
  #define OutputDir "."
#endif

#define AppName "Property Estimation Flowsheet"
#define AppExe "thermo-desktop.exe"

[Setup]
AppId={{8C2F1E4A-6B1D-4E7A-9A51-3F0C2D7B9E10}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=Property Estimation project
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog commandline
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
OutputDir={#OutputDir}
OutputBaseFilename=PropertyEstimationFlowsheet-{#AppVersion}-x64-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
UninstallDisplayIcon={app}\{#AppExe}
UninstallDisplayName={#AppName}

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#StageDir}\{#AppExe}"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#StageDir}\README.txt"; DestDir: "{app}"; Flags: ignoreversion
Source: "{#StageDir}\notices\*"; DestDir: "{app}\notices"; Flags: ignoreversion recursesubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{group}\{#AppName} (demonstration)"; Filename: "{app}\{#AppExe}"; Parameters: "--demo"
Name: "{group}\Third-party notices"; Filename: "{app}\notices\NOTICE.txt"
Name: "{group}\{cm:UninstallProgram,{#AppName}}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Parameters: "--demo"; Description: "Open the demonstration project"; Flags: nowait postinstall skipifsilent
