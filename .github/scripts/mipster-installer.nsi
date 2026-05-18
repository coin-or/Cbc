; mipster-installer.nsi — NSIS installer script for MIPster on Windows x86_64.
;
; Produces: mipster-@VERSION@-setup-windows-x86_64.exe
;
; Features:
;   - Installs to C:\Program Files\MIPster (user-changeable)
;   - Adds bin\ to the system PATH (Machine scope, requires admin)
;   - Registers in Add/Remove Programs with a working uninstaller
;   - Uninstall removes bin\ from PATH cleanly
;
; Build with makensis (called from build-windows.sh after substituting placeholders):
;   @VERSION@      — version string (e.g. devel, 3.0.0)
;   @DIST_DIR@     — Windows path to dist/ directory
;   @DIST_SUBDIR@  — name of the dist subdirectory (mipster-windows-x86_64)

Unicode True

!define APPNAME     "MIPster"
!define APPVER      "@VERSION@"
!define PUBLISHER   "h-g-s"
!define INSTDIR_KEY "Software\MIPster"
!define UNINST_KEY  "Software\Microsoft\Windows\CurrentVersion\Uninstall\MIPster"
!define DIST_DIR    "@DIST_DIR@"
!define DIST_SUBDIR "@DIST_SUBDIR@"

Name            "${APPNAME} ${APPVER}"
OutFile         "${DIST_DIR}\mipster-${APPVER}-setup-windows-x86_64.exe"
InstallDir      "$PROGRAMFILES64\MIPster"
InstallDirRegKey HKLM "${INSTDIR_KEY}" "InstallDir"
RequestExecutionLevel admin
SetCompressor   /SOLID lzma

; ── Pages ─────────────────────────────────────────────────────────────────────
!include "MUI2.nsh"
!define MUI_ABORTWARNING
!define MUI_ICON    "${NSISDIR}\Contrib\Graphics\Icons\modern-install.ico"
!define MUI_UNICON  "${NSISDIR}\Contrib\Graphics\Icons\modern-uninstall.ico"

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"

; ── PATH helpers (via PowerShell — handles broadcast notification) ─────────────

!macro AddToPath
  DetailPrint "Adding $INSTDIR\bin to system PATH..."
  nsExec::ExecToLog 'powershell.exe -NoProfile -NonInteractive -Command \
    "$p = [Environment]::GetEnvironmentVariable(\"Path\", \"Machine\"); \
     $b = \"$INSTDIR\bin\"; \
     if (($p -split \";\") -notcontains $b) { \
       [Environment]::SetEnvironmentVariable(\"Path\", $p + \";\" + $b, \"Machine\") \
     }"'
!macroend

!macro RemoveFromPath
  DetailPrint "Removing $INSTDIR\bin from system PATH..."
  nsExec::ExecToLog 'powershell.exe -NoProfile -NonInteractive -Command \
    "$b = \"$INSTDIR\bin\"; \
     $p = [Environment]::GetEnvironmentVariable(\"Path\", \"Machine\"); \
     $p2 = ($p -split \";\" | Where-Object { $_ -ne $b }) -join \";\"; \
     [Environment]::SetEnvironmentVariable(\"Path\", $p2, \"Machine\")"'
!macroend

; ── Install section ───────────────────────────────────────────────────────────
Section "MIPster" SecMain
  SectionIn RO  ; mandatory

  SetOutPath "$INSTDIR\bin"
  File "${DIST_DIR}\${DIST_SUBDIR}\bin\mipster.exe"
  File "${DIST_DIR}\${DIST_SUBDIR}\bin\mipster-generic.exe"
  File "${DIST_DIR}\${DIST_SUBDIR}\bin\mipster-avx2.exe"

  SetOutPath "$INSTDIR\lib\generic"
  File /nonfatal "${DIST_DIR}\${DIST_SUBDIR}\lib\generic\*.dll"
  File /nonfatal "${DIST_DIR}\${DIST_SUBDIR}\lib\generic\*.dll.a"

  SetOutPath "$INSTDIR\lib\avx2"
  File /nonfatal "${DIST_DIR}\${DIST_SUBDIR}\lib\avx2\*.dll"
  File /nonfatal "${DIST_DIR}\${DIST_SUBDIR}\lib\avx2\*.dll.a"

  SetOutPath "$INSTDIR\include\coin-or"
  File /r "${DIST_DIR}\${DIST_SUBDIR}\include\coin-or\*.*"

  ; Registry: save install dir and version
  WriteRegStr HKLM "${INSTDIR_KEY}" "InstallDir" "$INSTDIR"
  WriteRegStr HKLM "${INSTDIR_KEY}" "Version"    "${APPVER}"

  ; Add/Remove Programs entry
  WriteRegStr   HKLM "${UNINST_KEY}" "DisplayName"          "${APPNAME} ${APPVER}"
  WriteRegStr   HKLM "${UNINST_KEY}" "DisplayVersion"        "${APPVER}"
  WriteRegStr   HKLM "${UNINST_KEY}" "Publisher"             "${PUBLISHER}"
  WriteRegStr   HKLM "${UNINST_KEY}" "UninstallString"       "$INSTDIR\uninstall.exe"
  WriteRegStr   HKLM "${UNINST_KEY}" "QuietUninstallString"  "$INSTDIR\uninstall.exe /S"
  WriteRegStr   HKLM "${UNINST_KEY}" "InstallLocation"       "$INSTDIR"
  WriteRegDWORD HKLM "${UNINST_KEY}" "NoModify"              1
  WriteRegDWORD HKLM "${UNINST_KEY}" "NoRepair"              1

  ; Add to system PATH
  !insertmacro AddToPath

  ; Write uninstaller last
  WriteUninstaller "$INSTDIR\uninstall.exe"
SectionEnd

; ── Uninstall section ─────────────────────────────────────────────────────────
Section "Uninstall"
  ; Remove from PATH first (InstallDir is still known from registry)
  !insertmacro RemoveFromPath

  Delete "$INSTDIR\uninstall.exe"
  Delete "$INSTDIR\bin\mipster.exe"
  Delete "$INSTDIR\bin\mipster-generic.exe"
  Delete "$INSTDIR\bin\mipster-avx2.exe"
  RMDir  "$INSTDIR\bin"
  RMDir /r "$INSTDIR\lib"
  RMDir /r "$INSTDIR\include"
  RMDir  "$INSTDIR"

  DeleteRegKey HKLM "${UNINST_KEY}"
  DeleteRegKey HKLM "${INSTDIR_KEY}"
SectionEnd
