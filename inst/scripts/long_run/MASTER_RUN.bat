@echo off
REM ===========================================================================
REM LONG BALANCED — MASTER ORCHESTRATOR (unattended, sequential+conditional)
REM CANONICAL VERSIONED TEMPLATE (inst/scripts/long_run/).
REM   env snapshot -> SHARED upstream -> FULL -> (gate) -> NO-HOTSPOT
REM Failure policy: if a stage fails, STOP, write a FAILURE report, do NOT
REM launch downstream. NO-HOTSPOT only runs if FULL_SUCCESS.marker exists.
REM Each Rscript writes its own log; this master writes MASTER_LOG.txt.
REM
REM Route architecture (G2.3): SHARED upstream products are exposed to each
REM profile as EXPLICIT cfg routes (00_common.R SHARED_* + labelled_features =
REM SHARED_FEATURES_GPKG in the SCORE call). NO directory junctions are used or
REM required. validate_shared_inputs() fails fast before any heavy compute.
REM ===========================================================================
REM EDIT the four <PLACEHOLDER> paths below for your run, then launch.
setlocal enableextensions

set "RSCRIPT=<PATH_TO_Rscript.exe>"
set "ORCH=<PATH_TO_ORCHESTRATION>"
set "ROOT=<PATH_TO_LONG_RUN_ROOT>"
set "SNAPID=<SNAPSHOT_ID>"
set "MLOG=%ORCH%\MASTER_LOG.txt"
set "SHARED_MARK=%ROOT%\SHARED\SHARED_SUCCESS.marker"
set "FULL_MARK=%ROOT%\LONG_2017_BALANCED_FULL\FULL_SUCCESS.marker"
set "NOHS_MARK=%ROOT%\LONG_2017_BALANCED_NO_HOTSPOT\NO_HOTSPOT_SUCCESS.marker"

echo [%date% %time%] MASTER START >> "%MLOG%"

REM ---- STAGE 0: ENV SNAPSHOT -------------------------------------------------
echo [%date% %time%] STAGE env_snapshot >> "%MLOG%"
"%RSCRIPT%" "%ORCH%\05_env_snapshot.R" >> "%ORCH%\env_snapshot_run.log" 2>&1
if errorlevel 1 (
  echo [%date% %time%] FAILURE: env_snapshot ^(exit %errorlevel%^) >> "%MLOG%"
  echo FAILURE stage=env_snapshot exit=%errorlevel% see env_snapshot_run.log > "%ORCH%\FAILURE_REPORT.txt"
  echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
  echo fix=resolve env/path error before relaunch; do NOT patch package code >> "%ORCH%\FAILURE_REPORT.txt"
  goto :stop
)

REM ---- STAGE 1: SHARED UPSTREAM ---------------------------------------------
echo [%date% %time%] STAGE shared_upstream >> "%MLOG%"
"%RSCRIPT%" "%ORCH%\10_shared_upstream.R" >> "%ROOT%\SHARED\shared_upstream_run.log" 2>&1
if errorlevel 1 goto :fail_shared
if not exist "%SHARED_MARK%" goto :fail_shared
echo [%date% %time%] SHARED OK >> "%MLOG%"

REM ---- STAGE 2: FULL PROFILE ------------------------------------------------
echo [%date% %time%] STAGE FULL >> "%MLOG%"
"%RSCRIPT%" "%ORCH%\20_full.R" >> "%ROOT%\LONG_2017_BALANCED_FULL\full_run.log" 2>&1
if errorlevel 1 goto :fail_full
if not exist "%FULL_MARK%" goto :fail_full
echo [%date% %time%] FULL OK >> "%MLOG%"

REM ---- STAGE 3: NO-HOTSPOT (ONLY if FULL closed cleanly) --------------------
echo [%date% %time%] STAGE NO_HOTSPOT >> "%MLOG%"
"%RSCRIPT%" "%ORCH%\30_no_hotspot.R" >> "%ROOT%\LONG_2017_BALANCED_NO_HOTSPOT\no_hotspot_run.log" 2>&1
if errorlevel 1 goto :fail_nohs
if not exist "%NOHS_MARK%" goto :fail_nohs
echo [%date% %time%] NO_HOTSPOT OK >> "%MLOG%"

echo [%date% %time%] MASTER COMPLETE (FULL + NO_HOTSPOT) >> "%MLOG%"
echo MASTER_COMPLETE %date% %time% > "%ORCH%\MASTER_COMPLETE.marker"
goto :eof

:fail_shared
echo [%date% %time%] FAILURE: shared_upstream >> "%MLOG%"
echo FAILURE stage=shared_upstream see SHARED\shared_upstream_run.log + SHARED\shared_upstream.log > "%ORCH%\FAILURE_REPORT.txt"
echo last_valid_output=none ^(upstream is first heavy stage^) >> "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: FULL and NO_HOTSPOT NOT launched. Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:fail_full
echo [%date% %time%] FAILURE: FULL >> "%MLOG%"
echo FAILURE stage=FULL see LONG_2017_BALANCED_FULL\full_run.log + full.log > "%ORCH%\FAILURE_REPORT.txt"
echo last_valid_output=SHARED upstream ^(intact^) >> "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: NO_HOTSPOT NOT launched (FULL did not close). Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:fail_nohs
echo [%date% %time%] FAILURE: NO_HOTSPOT >> "%MLOG%"
echo FAILURE stage=NO_HOTSPOT see LONG_2017_BALANCED_NO_HOTSPOT\no_hotspot_run.log + no_hotspot.log > "%ORCH%\FAILURE_REPORT.txt"
echo last_valid_output=FULL profile (INTACT, kept). >> "%ORCH%\FAILURE_REPORT.txt"
echo snapshot=%SNAPID% >> "%ORCH%\FAILURE_REPORT.txt"
echo NOTE: FULL kept intact. Do NOT patch/relaunch with changes. >> "%ORCH%\FAILURE_REPORT.txt"
goto :stop

:stop
echo [%date% %time%] MASTER STOPPED ON FAILURE >> "%MLOG%"
endlocal
exit /b 1
