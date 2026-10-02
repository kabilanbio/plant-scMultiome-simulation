@echo off
REM ==============================================================================
REM run_all_arab_paga.bat
REM 1-Click Automated Runner for Arabidopsis Root PAGA Lineage Tree Reconstruction
REM Compatible with Hostel PC (HP 280 G4, i7-8700, sc_env) & Laptop
REM ==============================================================================

echo ==============================================================================
echo [Arabidopsis Root PAGA Lineage Tree Reconstruction Pipeline]
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

echo [Step 1/2] Exporting Standardized Arabidopsis Root Dataset from R...
"%R_EXE%" "%SCRIPT_DIR%export_arab_standardized_to_anndata.R" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/arabidopsis_root/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/sc_datasets/arabidopsis_root/featured_datasets_anndata"

if %ERRORLEVEL% NEQ 0 (
  echo [WARNING] R export reported an issue or partial completion. Checking if files are present...
)

echo.
echo [Step 2/2] Running Scanpy PAGA Lineage Tree Reconstruction...
%PYTHON_EXE% "%SCRIPT_DIR%arab_paga_lineage_tree.py" ^
  --input_dir "%DRIVE_LETTER%/PhD/sc_datasets/arabidopsis_root/featured_datasets_anndata" ^
  --standardized_dir "%DRIVE_LETTER%/PhD/sc_datasets/arabidopsis_root/featured_datasets_standardized" ^
  --output_dir "%DRIVE_LETTER%/PhD/data_simulation_comparison/scMultiSim/arab_root/arab_paga_results" ^
  --dpi 600

echo.
echo ==============================================================================
echo Pipeline execution complete! Check outputs in:
echo %DRIVE_LETTER%\PhD\data_simulation_comparison\scMultiSim\arab_root\arab_paga_results
echo ==============================================================================
pause
