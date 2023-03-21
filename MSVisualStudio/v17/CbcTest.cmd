@REM This file will run command line tests, comparable to test\makefile.am

@echo off
SET "BINDIR=%~1"
SET "SAMPLEDIR=%~2"
SET "NETLIBDIR=%~3"
SET "MIPLIBDIR=%~4"

echo INFO: Running %0 %*
echo INFO: Using bindir '%BINDIR%' and sampledir '%SAMPLEDIR%'
echo INFO: Netlibdir '%NETLIBDIR%' and miplibdir '%MIPLIBDIR%' are ignored for these tests.

if "%BINDIR%"=="" echo ERROR: No bindir given. && goto :usage
if not exist "%BINDIR%" echo ERROR: Folder bindir %BINDIR% does not exist. && goto :usage

if "%SAMPLEDIR%"=="" echo ERROR: No sampledir given. && goto :usage
if not exist %SAMPLEDIR% echo ERROR: Folder sampledir %SAMPLEDIR% does not exist. && goto :usage
if not %errorlevel%==0 echo ERROR: %SAMPLEDIR% cannot contain spaces. && goto :usage

goto :test

:usage
echo INFO: Usage %0 ^<bindir^> ^<sampledir^>
echo INFO: where ^<bindir^> contains the executables, sampledir the sample files.
echo INFO: This script runs automated test.
echo INFO: The ^<sampledir^> must not contain spaces!
echo INFO: For example: %0 "D:\Some Directory\" ..\..\samples\
goto :error

:error
echo ERROR: An error occurred while running the tests!
echo INFO: Finished Tests but failed
exit /b 1

:test
echo INFO: Starting Tests
@echo on

"%BINDIR%\gamsTest.exe"
@if not %errorlevel%==0 goto :error

"%BINDIR%\osiUnitTest.exe" -mpsDir=%SAMPLEDIR%
@if not %errorlevel%==0 goto :error

"%BINDIR%\cbc.exe" -dirSample %SAMPLEDIR% -unitTest
@if not %errorlevel%==0 goto :error

@echo off
echo INFO: Finished Tests successfully (%ERRORLEVEL%)