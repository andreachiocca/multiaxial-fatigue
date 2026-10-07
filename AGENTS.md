# Working on this repository

- Preserve the existing numerical behavior unless the task explicitly asks to change it. Check the experimental workbook schema and CSV output columns before editing shared data flow.
- Add a critical-plane method in `fatigue/models/cp_methods/`; add a direct method in `fatigue/models/invariant_methods/` or `fatigue/models/energy_based_methods/`. Register new methods in `fatigue/models/registry.py`.
- Follow `docs/adding-methods.md` for model interfaces and naming conventions, including `_ext` variants. See `docs/architecture.md` for the data flow.
- The three numbered scripts are current entry points. Run them from the repository root because data and output paths are relative to it.
- Do not commit Python caches or `Results/Comparison/`. Treat the existing `Results/*.csv` files as reference outputs; review their diff if a run overwrites them.
- For a new method or numerical change, add focused tests or a reproducible comparison against existing results. Do not claim validation solely from import or syntax checks.
