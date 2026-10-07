# Multiaxial fatigue methods

Python workflow for calibrating and comparing multiaxial fatigue methods against experimental material datasets. The implementation already separates critical-plane methods, invariant methods, and energy-based methods through a central registry.

## Setup

Use Python 3.10 or newer from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

## Run

1. Edit `MATERIAL_NAME` and `MODEL_NAME` near the top of `1_run_material.py`.
2. Run `python 1_run_material.py`. It reads a workbook in `Experimental_Data/` and writes a CSV in `Results/`.
3. Run `python 2_compare_results.py` to compare CSV results and generate figures and tables in `Results/Comparison/`.
4. Optionally run `python 0_calibrate_material_params.py` to print suggested material parameters. This does not modify the workbooks.

The existing CSV files in `Results/` are example outputs from the supplied archive. `Results/Comparison/` contains generated plots and is excluded from Git.

## Add a fatigue method

Implement a method in `fatigue/models/cp_methods/`, `fatigue/models/invariant_methods/`, or `fatigue/models/energy_based_methods/`, then register it in `fatigue/models/registry.py`. See [adding a method](docs/adding-methods.md) for the expected interfaces and naming conventions. The [architecture notes](docs/architecture.md) describe the data flow and CSV fields.

## Layout

| Path | Purpose |
| --- | --- |
| `fatigue/imp/` | Excel data import |
| `fatigue/planes/` | Critical-plane sampling |
| `fatigue/models/` | Method implementations, adapters, and registry |
| `fatigue/eval/` | Calibration and evaluation |
| `fatigue/report/` | CSV export |
| `Experimental_Data/` | Supplied material workbooks |
| `Results/` | Supplied result CSVs and new output |

This is the supplied v7.2 implementation with repository metadata added. Numerical methods and experimental data have not been changed.
