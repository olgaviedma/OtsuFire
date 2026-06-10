@echo off
REM ===========================================================================
REM LONG BALANCED — RESUME ORCHESTRATOR (skip the ~52-min SHARED rebuild)
REM CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
REM
REM Use when SHARED upstream is already done (SHARED_SUCCESS.marker present) and
REM you only need to (re)run the cheap tail: FULL train+score+map+EFFIS, then
REM NO-HOTSPOT. SHARED is reused untouched.
REM
REM CANONICAL route architecture (G2.3): the SCORE step finds the SHARED features
REM purely from the EXPLICIT cfg route (labelled_features = SHARED_FEATURES_GPKG
REM in 20_full.R / 30_no_hotspot.R). NO directory junctions are created or
REM required (the old junction-based RESUME workaround is GONE). The runners'
REM validate_shared_inputs() fails fast if any SHARED product is unresolvable.
REM
REM Same failure policy as MASTER_RUN.bat: on any failure STOP, write
REM FAILURE_REPORT.txt, do NOT launch downstream. Writes RESUME_LOG.txt and
REM per-stage *_resume logs (does NOT clobber the original run logs).
REM ===========================================================================
REM EDIT the four <PLACEHOLDER> paths below for your run, then launch.
setlocal enableextensions

set "RSCRIPT=<PATH_TO_Rscript.exe>"
set "ORCH=<PATH_TO_ORCHESTRATION>"
set "ROOT=<PATH_TO_LONG_RUN_ROOT>"
set "SNAPID=<SNAPSHOT_ID>"
set "RLOG=%ORCH%\RESUME_LOG.txt"
set "SHARED_MARK=%ROOT%\SHARED\SHARED_SUCCESS.marker"
set "FULL_MARK=%ROOT%\LONG_2017_BALANCED_FULL\FULL_SUCCESS.marker"
set "NOHS_MARK=%ROOT%\LONG_2017_BALANCED_NO_HOTSPOT\NO_HOTSPOT_SUCCESS.marker"
set "SHARED_FEAT=%ROOT%\SHARED\ENGINE_ROUTES\2017\Min_Min\SUPERVISED\balanced\03_FEATURES\features_geometry.gpkg"

echo [%date% %time%] RESUME START >> "%RLOG%"

REM ---- PRECHECK A: SHARED upstream must already be done (skip the rebuild) -----
if not exist "%SHARED_MARK%" goto :fail_no_shared
echo [%date% %time%] PRECHECK SHARED_SUCCESS.marker present (rebuild skipped) >> "%RLOG%"

REM ---- PRECHECK B: the SHARED features GPKG (the explicit labelled_features) ---
REM     must resolve. NO junction is needed; the runners pass this path directly.
if not exist "%SHARED_FEAT%" goto :fail_no_features
echo [%date% %time%] PRECHECK SHARED features_geometry.gpkg resolves (explicit route, no junction) >> "%RLOG%"

REM ---- Clear any stale FULL marker so the gate is honest -----------------------
if exist "%FULL_MARK%" del /q "%FULL_MARK%"

REM ---- STAGE 1: FULL TAIL (re-run 20_full.R; SCORE finds SHARED features) ------
echo [%date% %time%] STAGE FULL (resume; deterministic A/B/C retrain + SCORE+MAP+EFFIS) >> "%RLOG%"
"%RSCRIPT%" "%ORCH%\20_full.R" >> "%ROOT%\LONG_2017_BALANCED_FULL\full_resume_run.log" 2>&1
if errorlevel 1 goto :fail_full
if not exist "%FULL_MARK%" goto :fail_full
echo [%date% %time%] FULL OK >> "%RLOG%"

REM ---- STAGE 2: NO-HOTSPOT (ONLY if FULL closed cleanly) -----------------------
echo [%date% %time%] STAGE NO_HOTSPOT (train 38-feature model + SCORE+MAP+EFFIS) >> "%RLOG%"
"%RSCRIPT%" "%ORCH%\30_no_hotspot.R" >> "%ROOT%\LONG_2017_BALANCED_NO_HOTSPOT\no_hotspot_resume_run.log" 2>&1
if errorlevel 1 goto :fail_nohs
if not exist "%NOHS_MARK%" goto :fail_nohs
echo [%date% %time%] NO_HOTSPOT OK >> "%RLOG%"

echo [%date% %time%] RESUME COMPLETE (FULL + NO_HOTSPOT) >> "%RLOG%"
echo RESUME_COMPLETE %date% %time% > "%ORCH%\RESUME_COMPLETE.marker"
goto :eof

:fail_no_shared
echo [%date% %time%] FAILURE: SHARED_SUCCESS.marker MISSING; refusing to resume >> "%RLOG%"
echo FAILURE stage=precheck reason=SHARED_SUCCESS.marker missing (cannot resume without SHARED upstream) > "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: nothing launched. Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:fail_no_features
echo [%date% %time%] FAILURE: SHARED features_geometry.gpkg does not resolve >> "%RLOG%"
echo FAILURE stage=precheck reason=SHARED 03_FEATURES features_geometry.gpkg missing > "%ORCH%\FAILURE_REPORT.txt"
echo expected=%SHARED_FEAT% >> "%ORCH%\FAILURE_REPORT.txt"
echo fix=re-run 10_shared_upstream.R; do NOT patch package code >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:fail_full
echo [%date% %time%] FAILURE: FULL >> "%RLOG%"
echo FAILURE stage=FULL see LONG_2017_BALANCED_FULL\full_resume_run.log + full.log > "%ORCH%\FAILURE_REPORT.txt"
echo last_valid_output=SHARED upstream (intact) >> "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: NO_HOTSPOT NOT launched (FULL did not close). Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:fail_nohs
echo [%date% %time%] FAILURE: NO_HOTSPOT >> "%RLOG%"
echo FAILURE stage=NO_HOTSPOT see LONG_2017_BALANCED_NO_HOTSPOT\no_hotspot_resume_run.log + no_hotspot.log > "%ORCH%\FAILURE_REPORT.txt"
echo last_valid_output=FULL profile (INTACT, kept). >> "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: FULL kept intact. Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:stop
echo [%date% %time%] RESUME STOPPED ON FAILURE >> "%RLOG%"
endlocal
exit /b 1
