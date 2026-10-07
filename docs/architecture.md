# Architecture and data flow

`1_run_material.py` reads one material workbook from `Experimental_Data/`, selects uniaxial design points, fits a power-law relation between a method's damage parameter and cycles to failure, evaluates the remaining cases, and writes `Results/<METHOD>_<MATERIAL>.csv`. `2_compare_results.py` reads those CSVs and generates plots and tables under `Results/Comparison/`. `0_calibrate_material_params.py` prints suggested material parameters without editing the workbooks.

| Path | Role |
| --- | --- |
| `fatigue/imp/experimental_data.py` | Workbook import, material properties, stress/strain tensors and harmonics |
| `fatigue/planes/sampling.py` | Candidate plane rotations; each plane normal is the third column of its rotation matrix |
| `fatigue/models/base.py` | Model interfaces and result types |
| `fatigue/models/cp_adapter.py` | Case-level wrapper around critical-plane methods |
| `fatigue/models/registry.py` | Names exposed to the runner |
| `fatigue/eval/` | Critical-plane evaluation, calibration and fitting |
| `fatigue/report/export_csv.py` | Per-method CSV output |

The importer reads material properties from the `Summary` sheet and cases from the `Dataset` sheet. It filters invalid/notched/run-out cases. Optional radial and torsional phase columns produce stress and strain harmonics (`mean`, `sin`, `cos`) for single-frequency sinusoidal loading. The critical-plane evaluator uses these for non-proportional loading while retaining the legacy two-loadstep interface (`S0`/`E0` and `S1`/`E1`). Direct methods receive harmonics but must use them explicitly.

The main CSV fields include `Material`, `Method`, `Method family`, `Loading type`, `Flag`, `Metric_mode`, `Nf_exp`, `Nf_expected`, `CP_value`, `DP_fit`, `Error_ln`, and `Error_ln_dp`. `Flag=0` marks calibration points. `Metric_mode=Nf` compares predicted and experimental cycles; `DP_DIFF` compares damage parameters when the calibration slope is nearly flat (the runner switches when its slope parameter exceeds 15). Preserve these field names when extending reporting.
