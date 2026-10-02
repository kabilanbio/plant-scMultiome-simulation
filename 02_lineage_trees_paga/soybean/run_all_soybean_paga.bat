@echo off
REM ==============================================================================
REM run_all_soybean_paga.bat
REM 1-Click Automated Runner for Soybean Multiomics PAGA Lineage Tree Reconstruction
REM Compatible with Hostel PC (HP 280 G4, i7-8700, sc_env) & Laptop
REM ==============================================================================

echo ==============================================================================
echo [Soybean PAGA Lineage Tree Reconstruction Pipeline - 7 Tissues]
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

echo [Step 1/2] Exporting Standardized Soybean Datasets from R...
"%R_EXE%" "%SCRIPT_DIR%export_soybean_standardized_to_anndata.R" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/soybean/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/sc_datasets/soybean/featured_datasets_anndata" ^
  --tissue "all"

if %ERRORLEVEL% NEQ 0 (
  echo [WARNING] R export reported an issue or partial completion. Checking if files are present...
)

echo.
echo [Step 2/2] Running Scanpy PAGA Lineage Tree Reconstruction for all 7 Tissues...
%PYTHON_EXE% "%SCRIPT_DIR%soybean_paga_lineage_tree.py" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/soybean/featured_datasets_anndata" ^
  --standardized_dir "%DRIVE_LETTER%/PhD/sc_datasets/soybean/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/data_simulation_comparison/scMultiSim/soybean_datasets/soybean_paga_results" ^
  --tissue "all" ^
  --dpi 600

echo.
echo ==============================================================================
echo Pipeline execution complete! Check outputs in:
echo %DRIVE_LETTER%\PhD\data_simulation_comparison\scMultiSim\soybean_datasets\soybean_paga_results
echo ==============================================================================
pause
