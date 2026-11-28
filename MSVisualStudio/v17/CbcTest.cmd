@REM This file will run command line tests, comparable to test\makefile.am

@echo off
SET "BINDIR=%~1"
SET "SAMPLEDIR=%~2"
SET "NETLIBDIR=%~3"
SET "MIPLIBDIR=%~4"

echo INFO: Running %0 %*
echo INFO: Using bindir '%BINDIR%' and sampledir '%SAMPLEDIR%'.
echo INFO: Netlibdir '%NETLIBDIR%' is ignored for these tests.
echo INFO: Miplibdir '%MIPLIBDIR%' is ignored for these tests.

if "%BINDIR%"=="" echo ERROR: No bindir given. && goto :usage
if not exist "%BINDIR%" echo ERROR: Folder bindir %BINDIR% does not exist. && goto :usage

if "%SAMPLEDIR%"=="" echo ERROR: No sampledir given. && goto :usage
if not exist %SAMPLEDIR% echo ERROR: Folder sampledir %SAMPLEDIR% does not exist. && goto :usage
if not %errorlevel%==0 echo ERROR: %SAMPLEDIR% cannot contain spaces. && goto :usage

@REM if "%NETLIBDIR%"=="" echo ERROR: No netlibdir given. && goto :usage
@REM if not exist %NETLIBDIR% echo ERROR: Folder netlibdir %NETLIBDIR% does not exist. && goto :usage
@REM if not %errorlevel%==0 echo ERROR: %NETLIBDIR% cannot contain spaces. && goto :usage

@REM if "%MIPLIBDIR%"=="" echo ERROR: No miplibdir given. && goto :usage
@REM if not exist %MIPLIBDIR% echo ERROR: Folder miplibdir %MIPLIBDIR% does not exist. && goto :usage
@REM if not %errorlevel%==0 echo ERROR: %MIPLIBDIR% cannot contain spaces. && goto :usage

goto :test

:usage
echo INFO: Usage %0 ^<bindir^> ^<sampledir^>
echo INFO: where ^<bindir^> contains the executables and sampledir the sample files.
echo INFO: This script runs automated test.
echo INFO: The sampledir must not contain spaces!
echo INFO: For example: %0 "D:\Some Directory\" ..\..\samples\
goto :error

:error
echo ERROR: An error occurred while running the tests!
echo INFO: Finished Tests but failed
exit /b 1

:test
echo INFO: Starting Tests
@echo on
@REM Runs make test tests
"%BINDIR%\cbcmaketest.exe"
@if not %errorlevel%==0 echo ERROR: cbcmaketest.exe exited with code %errorlevel% && goto :error

@REM Solve p0033 from the sample data; the full regression suite is make sanity-test.
"%BINDIR%\cbc.exe" %SAMPLEDIR%\p0033.mps -solve -quit
@if not %errorlevel%==0 echo ERROR: cbc.exe exited with code %errorlevel% && goto :error

@echo off
echo INFO: Finished Tests successfully (%ERRORLEVEL%)