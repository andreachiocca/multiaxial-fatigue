# Adding a fatigue method

The user-facing method name is registered in `fatigue/models/registry.py`. The runner selects `MODEL_NAME` and calibrates the names in `MODEL_NAMES` near the top of `1_run_material.py`. The `_ext` suffix selects the extreme critical-plane variant; calibrate its base method name.

## Critical-plane method

1. Add a class in `fatigue/models/cp_methods/`. Follow an existing class such as `FatemiSocie` in `fs.py` and the `CriticalPlaneMethod` protocol in `fatigue/models/base.py`.
2. Provide `name`, `required_params()`, and `evaluate_on_plane(*, S0r, E0r, S1r, E1r, params) -> CPPlaneResult`. Return finite `damage` and `metric`; the evaluator selects the plane with the highest metric. The local plane normal is the `z` axis, and the relevant in-plane shear components are `xz` and `yz`.
3. Import and register the class in `_available_cp_methods()` in `fatigue/models/registry.py`. The `CPModelAdapter` exposes the standard and `_ext` values.
4. Add the name to `MODEL_NAMES` in `1_run_material.py`, select it with `MODEL_NAME`, and run it on a material. Check finite positive values on the uniaxial design subset and compare proportional and non-proportional cases. The evaluator handles phase shifts for critical-plane methods.

## Direct method

1. Add a class in `fatigue/models/invariant_methods/` or `fatigue/models/energy_based_methods/`. Follow an existing model and the `FatigueModel` protocol in `fatigue/models/base.py`.
2. Provide `name`, `required_params()`, and `evaluate_case(*, S0, E0, S1, E1, params, R_list=None, harmonics=None) -> CaseResult`. The result uses `values={"": damage}`. A direct model must explicitly use `harmonics` if its calculation needs phase information.
3. Import and register the class in `_available_direct_models()` in `fatigue/models/registry.py`, then add and select it in `1_run_material.py`.

If a method needs a new material parameter, add it to the Excel `Summary` sheet, update the mapping in `fatigue/imp/experimental_data.py`, pass it in the runner's `params` dictionary, and declare it in `required_params()`.

For a new method family, update `_family_cmap_name()` in `2_compare_results.py` if the plots need a distinct color palette. If calibration produces NaN, check the method's values on the design subset and its required parameters first.
