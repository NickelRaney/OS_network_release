param([switch]$VerboseBuild)
$ErrorActionPreference = 'Stop'
$standaloneRoot = $PSScriptRoot
$vsWhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio/Installer/vswhere.exe'
if (-not (Test-Path -LiteralPath $vsWhere)) { throw 'Install Visual Studio C++ build tools first.' }
$vsInstall = & $vsWhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
if (-not $vsInstall) { throw 'Visual Studio x64 C++ compiler not found.' }
$vcVars = Join-Path $vsInstall 'VC/Auxiliary/Build/vcvars64.bat'
$buildDirectory = Join-Path $standaloneRoot 'build'
New-Item -ItemType Directory -Path $buildDirectory -Force | Out-Null
$buildBatch = Join-Path $buildDirectory 'compile.cmd'
# /MT links the C++ runtime statically; the executable needs no MATLAB DLLs.
$batchText = @"
@echo off
call "$vcVars" >nul
if errorlevel 1 exit /b %errorlevel%
cd /d "$buildDirectory"
cl /nologo /O2 /std:c++17 /EHsc /fp:strict /arch:AVX2 /MT /W4 /DNDEBUG /DMINIZ_NO_ARCHIVE_APIS /DMINIZ_NO_DEFLATE_APIS /I"$standaloneRoot" "$standaloneRoot/main.cpp" "$standaloneRoot/v1_simulator.cpp" "$standaloneRoot/v1_l6.cpp" "$standaloneRoot/mat5.cpp" "$standaloneRoot/third_party/miniz/miniz_tinfl.c" /Fe:"$standaloneRoot/v1_4_standalone.exe" /link /INCREMENTAL:NO
exit /b %errorlevel%
"@
Set-Content -LiteralPath $buildBatch -Value $batchText -Encoding Ascii
if ($VerboseBuild) { Write-Output $batchText }
& $env:ComSpec /d /c $buildBatch
if ($LASTEXITCODE -ne 0) { throw "C++ compilation failed ($LASTEXITCODE)." }
Write-Output "Built $(Join-Path $standaloneRoot 'v1_4_standalone.exe')"
