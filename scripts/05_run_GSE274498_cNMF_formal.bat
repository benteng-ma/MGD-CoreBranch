@echo off
setlocal

REM ============================================================
REM Formal cNMF analysis for GSE274498 meibomian-gland cells
REM
REM Validated analysis parameters:
REM   Input cells: 10,307
REM   Genes in input: 27,267
REM   K: 5-20
REM   NMF replicates: 100 per K
REM   HVGs: 3,000
REM   Master seed: 14
REM   beta loss: frobenius
REM   initialization: random
REM   max NMF iterations: 1000
REM   selected solution: K=8
REM   consensus local-density threshold: 0.1
REM
REM IMPORTANT:
REM This workflow can take many hours because it performs
REM 16 K values x 100 replicates = 1,600 NMF runs.
REM ============================================================


REM ------------------------------------------------------------
REM Resolve project directory relative to this script
REM ------------------------------------------------------------

set "SCRIPT_DIR=%~dp0"

for %%I in ("%SCRIPT_DIR%..") do (
    set "PROJECT_ROOT=%%~fI"
)


REM ------------------------------------------------------------
REM Input and output locations
REM ------------------------------------------------------------

set "COUNTS=%PROJECT_ROOT%\source_data\GSE274498\GSE274498_MG_all_rawcounts.h5ad"

set "OUTPUT_ROOT=%PROJECT_ROOT%\results\cnmf"

set "RUN_NAME=GSE274498_allMG_formal"


echo.
echo ============================================================
echo GSE274498 FORMAL cNMF
echo ============================================================
echo Project root:
echo %PROJECT_ROOT%
echo.
echo Counts:
echo %COUNTS%
echo.
echo Output root:
echo %OUTPUT_ROOT%
echo.
echo Run name:
echo %RUN_NAME%
echo ============================================================
echo.


REM ------------------------------------------------------------
REM Check input
REM ------------------------------------------------------------

if not exist "%COUNTS%" (
    echo ERROR: Input file not found:
    echo %COUNTS%
    exit /b 1
)

if not exist "%OUTPUT_ROOT%" (
    mkdir "%OUTPUT_ROOT%"
)


REM ------------------------------------------------------------
REM STEP 1
REM Prepare cNMF inputs and allocate NMF runs
REM ------------------------------------------------------------

echo.
echo [1/5] PREPARE
echo.

cnmf prepare ^
    --output-dir "%OUTPUT_ROOT%" ^
    --name "%RUN_NAME%" ^
    -c "%COUNTS%" ^
    -k 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ^
    --n-iter 100 ^
    --numgenes 3000 ^
    --seed 14 ^
    --max-nmf-iter 1000 ^
    --beta-loss frobenius ^
    --init random ^
    --total-workers 1

if errorlevel 1 (
    echo.
    echo ERROR during cNMF prepare.
    exit /b 1
)


REM ------------------------------------------------------------
REM STEP 2
REM Run all 1,600 NMF factorizations on one worker
REM ------------------------------------------------------------

echo.
echo [2/5] FACTORIZE
echo This is the long-running step.
echo.

cnmf factorize ^
    --output-dir "%OUTPUT_ROOT%" ^
    --name "%RUN_NAME%" ^
    --worker-index 0 ^
    --total-workers 1

if errorlevel 1 (
    echo.
    echo ERROR during cNMF factorize.
    exit /b 1
)


REM ------------------------------------------------------------
REM STEP 3
REM Combine replicate spectra for every K
REM ------------------------------------------------------------

echo.
echo [3/5] COMBINE
echo.

cnmf combine ^
    --output-dir "%OUTPUT_ROOT%" ^
    --name "%RUN_NAME%"

if errorlevel 1 (
    echo.
    echo ERROR during cNMF combine.
    exit /b 1
)


REM ------------------------------------------------------------
REM STEP 4
REM Calculate K-selection stability and reconstruction error
REM ------------------------------------------------------------

echo.
echo [4/5] K SELECTION
echo.

cnmf k_selection_plot ^
    --output-dir "%OUTPUT_ROOT%" ^
    --name "%RUN_NAME%"

if errorlevel 1 (
    echo.
    echo ERROR during cNMF K selection.
    exit /b 1
)


REM ------------------------------------------------------------
REM STEP 5
REM Final consensus solution
REM
REM K=8 was selected after inspection of K=5-20 because it
REM showed the highest consensus stability while retaining a
REM comparatively parsimonious solution.
REM
REM Formal density threshold = 0.1
REM ------------------------------------------------------------

echo.
echo [5/5] K=8 CONSENSUS
echo.

cnmf consensus ^
    --output-dir "%OUTPUT_ROOT%" ^
    --name "%RUN_NAME%" ^
    --components 8 ^
    --local-density-threshold 0.1 ^
    --show-clustering

if errorlevel 1 (
    echo.
    echo ERROR during cNMF consensus.
    exit /b 1
)


echo.
echo ============================================================
echo FORMAL cNMF COMPLETE
echo ============================================================
echo.
echo Results:
echo %OUTPUT_ROOT%\%RUN_NAME%
echo.
echo Final solution:
echo K = 8
echo local-density threshold = 0.1
echo.
echo ============================================================

endlocal