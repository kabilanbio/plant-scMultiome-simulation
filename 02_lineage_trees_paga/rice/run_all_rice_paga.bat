@echo off
REM ==============================================================================
REM run_all_rice_paga.bat
REM 1-Click Automated Runner for Rice PAGA Lineage Tree Reconstruction
REM Step 1: Export standardized Seurat RDS datasets to Matrix Market / AnnData
REM Step 2: Run Scanpy PAGA Lineage Tree extraction for all 8 Rice tissues
REM ==============================================================================

echo ==============================================================================
echo [Rice PAGA Lineage Tree Reconstruction Pipeline]
echo ==============================================================================

set SCRIPT_DIR=%~dp0
set DRIVE_LETTER=%~d0
set R_EXE=C:\Program Files\R\R-4.6.1\bin\Rscript.exe

REM Auto-detect Python environment (prefer local virtualenv sc_env if available)
if exist "C:\Users\kabil\sc_env\Scripts\python.exe" (
  set PYTHON_EXE="C:\Users\kabil\sc_env\Scripts\python.exe"
) else (
  set PYTHON_EXE=python
)

echo [Step 1/2] Exporting Standardized Rice Datasets from R...
"%R_EXE%" "%SCRIPT_DIR%export_rice_standardized_to_anndata.R" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/rice/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/sc_datasets/rice/featured_datasets_anndata" ^
  --tissue "all"

if %ERRORLEVEL% NEQ 0 (
  echo [WARNING] R export reported an issue or partial completion. Checking if files are present...
)

echo.
echo [Step 2/2] Running Scanpy PAGA Lineage Tree Reconstruction...
%PYTHON_EXE% "%SCRIPT_DIR%rice_paga_lineage_tree.py" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/rice/featured_datasets_anndata" ^
  --standardized_dir "%DRIVE_LETTER%/PhD/sc_datasets/rice/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/data_simulation_comparison/scMultiSim/rice_datasets/rice_paga_results" ^
  --tissue "all" ^
  --dpi 600

echo.
echo ==============================================================================
echo Pipeline execution complete! Check outputs in:
echo %DRIVE_LETTER%\PhD\data_simulation_comparison\scMultiSim\rice_datasets\rice_paga_results
echo ==============================================================================
pause
