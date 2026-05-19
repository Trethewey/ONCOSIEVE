@echo off
REM ============================================================================
REM ONCOSIEVE v1.0 -- Windows One-Click Installer & Runner
REM
REM Double-click this file to:
REM   1. Check Python 3.10+ is installed
REM   2. Install Python dependencies
REM   3. Check for bcftools/tabix (via conda or manual)
REM   4. Run the full ONCOSIEVE pipeline
REM
REM Requirements:
REM   - Python 3.10+ on PATH
REM   - Git Bash or WSL for bcftools/tabix (or conda with bioconda channel)
REM   - Data files populated in data/ directory
REM ============================================================================

setlocal enabledelayedexpansion
title ONCOSIEVE v1.0 -- Pan-Cancer Whitelist Curation

echo.
echo   ================================================================
echo                         O N C O S I E V E  v1.0
echo        Pan-Cancer Somatic Variant Whitelist Curation Tool
echo   ================================================================
echo.

REM Navigate to script directory
cd /d "%~dp0\.."
set "ONCOSIEVE_DIR=%cd%"
echo   Working directory: %ONCOSIEVE_DIR%
echo.

REM ── Step 1: Check Python ──────────────────────────────────────────
echo   [1/5] Checking Python installation...
python --version >nul 2>&1
if %errorlevel% neq 0 (
    python3 --version >nul 2>&1
    if %errorlevel% neq 0 (
        echo.
        echo   ERROR: Python not found on PATH.
        echo   Please install Python 3.10+ from https://www.python.org/downloads/
        echo   Make sure to check "Add Python to PATH" during installation.
        echo.
        pause
        exit /b 1
    )
    set "PYTHON=python3"
) else (
    set "PYTHON=python"
)

for /f "tokens=2 delims= " %%v in ('%PYTHON% --version 2^>^&1') do set PYVER=%%v
echo         Python %PYVER% found.

REM Check Python version >= 3.10
for /f "tokens=1,2 delims=." %%a in ("%PYVER%") do (
    if %%a LSS 3 (
        echo   ERROR: Python 3.10+ required, found %PYVER%
        pause
        exit /b 1
    )
    if %%a EQU 3 if %%b LSS 10 (
        echo   ERROR: Python 3.10+ required, found %PYVER%
        pause
        exit /b 1
    )
)

REM ── Step 2: Create virtual environment ────────────────────────────
echo   [2/5] Setting up virtual environment...
if not exist "%ONCOSIEVE_DIR%\.venv\Scripts\python.exe" (
    echo         Creating virtual environment...
    %PYTHON% -m venv "%ONCOSIEVE_DIR%\.venv"
    if %errorlevel% neq 0 (
        echo   ERROR: Failed to create virtual environment.
        pause
        exit /b 1
    )
)
set "PYTHON=%ONCOSIEVE_DIR%\.venv\Scripts\python.exe"
set "PIP=%ONCOSIEVE_DIR%\.venv\Scripts\pip.exe"
echo         Virtual environment ready.

REM ── Step 3: Install Python dependencies ───────────────────────────
echo   [3/5] Installing Python dependencies...
"%PIP%" install --quiet --upgrade pip >nul 2>&1
"%PIP%" install --quiet -r "%ONCOSIEVE_DIR%\requirements.txt"
if %errorlevel% neq 0 (
    echo   ERROR: pip install failed. Check requirements.txt.
    pause
    exit /b 1
)
echo         Dependencies installed.

REM ── Step 4: Check system tools ────────────────────────────────────
echo   [4/5] Checking system tools...

set "HAS_BCFTOOLS=0"
set "HAS_TABIX=0"

where bcftools >nul 2>&1 && set "HAS_BCFTOOLS=1"
where tabix >nul 2>&1 && set "HAS_TABIX=1"

if "%HAS_BCFTOOLS%"=="0" (
    echo.
    echo   WARNING: bcftools not found on PATH.
    echo   Options:
    echo     a) Install via conda: conda install -c bioconda bcftools htslib
    echo     b) Use Git Bash / WSL where bcftools is available
    echo     c) Download from https://www.htslib.org/download/
    echo.
)
if "%HAS_TABIX%"=="0" (
    echo   WARNING: tabix not found on PATH.
    echo   tabix is part of htslib. Install alongside bcftools.
    echo.
)

if "%HAS_BCFTOOLS%"=="0" (
    echo   Cannot proceed without bcftools and tabix.
    echo   Install them and re-run this script.
    pause
    exit /b 1
)

REM ── Step 5: Check data files ──────────────────────────────────────
echo   [5/5] Checking data files...
"%PYTHON%" pre_check.py --config config.yaml
if %errorlevel% neq 0 (
    echo.
    echo   Pre-check failed. Fix the errors above and re-run.
    pause
    exit /b 1
)
echo.

REM ── Run Pipeline ──────────────────────────────────────────────────
echo   ================================================================
echo                     Starting ONCOSIEVE Pipeline
echo   ================================================================
echo.

set "TIMESTAMP=%date:~-4%%date:~4,2%%date:~7,2%_%time:~0,2%%time:~3,2%%time:~6,2%"
set "TIMESTAMP=%TIMESTAMP: =0%"
set "LOG_FILE=logs\run_%TIMESTAMP%.log"
if not exist logs mkdir logs

echo   Log file: %LOG_FILE%
echo.

REM Run build_whitelist.py
echo   Building whitelist...
"%PYTHON%" build_whitelist.py --config config.yaml 2>&1 | tee "%LOG_FILE%"
if %errorlevel% neq 0 (
    echo   ERROR: build_whitelist.py failed. Check %LOG_FILE%
    pause
    exit /b 1
)

REM Index VCF
echo   Indexing VCF...
tabix -f -p vcf "output\pan_cancer_whitelist_GRCh38.vcf.gz"

REM Post-process
echo   Post-processing...
"%PYTHON%" post_process_whitelist.py ^
    --whitelist "output\pan_cancer_whitelist_GRCh38.tsv.gz" ^
    --out "output\pan_cancer_whitelist_GRCh38_annotated.tsv.gz"

REM Post-pipeline (REVEL, Excel, HTML)
if exist "data\REVEL\revel_with_transcript_ids" (
    echo   Running post-pipeline (REVEL, Excel, HTML report)...
    "%PYTHON%" tools\post_pipeline.py ^
        --whitelist "output\pan_cancer_whitelist_GRCh38_annotated.tsv.gz" ^
        --vcf "output\pan_cancer_whitelist_GRCh38.vcf.gz" ^
        --revel "data\REVEL\revel_with_transcript_ids" ^
        --out-dir output\
) else (
    echo   REVEL file not found -- skipping post-pipeline.
)

echo.
echo   ================================================================
echo                     Pipeline Complete!
echo   ================================================================
echo.
echo   Output files in: %ONCOSIEVE_DIR%\output\
echo.
echo   Key files:
echo     - pan_cancer_whitelist_GRCh38.vcf.gz      (VCF)
echo     - pan_cancer_whitelist_GRCh38_full.tsv.gz  (Full TSV)
echo     - pan_cancer_whitelist_GRCh38_full.xlsx    (Excel)
echo     - oncosieve_report_*.html                  (HTML report)
echo.
pause
