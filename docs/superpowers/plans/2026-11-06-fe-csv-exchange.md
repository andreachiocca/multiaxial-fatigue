# FE CSV exchange format Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a CSV, exported by the user or a solver script, the recommended way to bring FE results in: a specified format with a header, a streaming importer into the HDF5 history file (extended to version 1.1), a standard-library writer helper, an example, and a node assessment driver joined to the load-history bridge.

**Architecture:** Extend `multifat.fe` additively (model, store, tensors, select), add `fe/csvio.py` (header, columns, importer), `fe/csvwriter.py` (stdlib only) and a CSV flavour of the synthetic tube generator; add `history/nodes.py` (node histories, unit-case superposition, `assess_nodes`). Nothing existing changes behaviour; version 1.0 files keep reading.

**Tech Stack:** Python 3.10+, numpy, pandas (CSV), optional `h5py` (extra `fe`), pytest, black (100), ruff, `uv`.

**Spec:** `docs/superpowers/specs/2026-11-06-fe-csv-exchange-design.md` (read first; its "Decisions from the builder check" section fixes the open rules).

## Global Constraints

- Internal units: MPa, mm, tensor strain (shear = engineering / 2), global Cartesian, six components `xx, yy, zz, xy, yz, xz`; the CSV may declare `Pa, kPa, MPa, GPa, psi, ksi, N/mm^2`, `m, cm, mm, in, ft`, a `strain_unit` (`1, mm/mm, m/m, percent, microstrain`, required when strain exists) and `tensor` or `engineering` shear, and is converted on import. Units are explicit tables, not `pint` or `xarray`.
- Nothing is inferred: required header keys, units, `coordinate_system: global_cartesian`, `strain_convention` when strain exists, `cycle: stabilised` for `elastic_plastic`; an unknown column, a canonical column together with its alias, a missing required column, a truncated or ragged grid, duplicate keys, non-finite or asymmetric tensors, `time` not strictly increasing along ascending `step` are refusals with a diagnostic.
- A **location** is the triple `(instance, node_id, material)`; keys `(instance, node_id, material, step)` are unique; a missing `instance` or `material` column is the empty string; `n_nodes` counts locations.
- Strain rules: `Eel` and `Epl` (total = sum); total `E_xx..E_xz` (Abaqus `E11..`) is total mechanical strain; `Eel` alone is elastic data; `Epl` alone needs `elastic_constants` (E and nu per material, MPa) and derives the elastic part by Hooke (recorded `derived_elastic_part`); `linear_elastic` with no strain column may derive (recorded `derived_linear_elastic`); elastic-plastic with no strain column is refused.
- `nodal_average_global` with more than one material or instance is refused unless `allow_global_average=True` (recorded).
- Step order is ascending `step` value (non-contiguous allowed, duplicates refused); step-major files stream (peak memory one step); node-major files are sorted in memory up to a row limit.
- `step_role: unit_case` needs `analysis_type: linear_elastic` and a load-function CSV with `time` (or `index`) and `f_<step>` columns matched by name.
- `assess_nodes` records a per-node `status` and continues by default (`on_refusal="record"`); unexpected exceptions propagate.
- Notation is data: rows in `docs/notation.md` before new public names; no retired synonyms (`ratio`, `mat`, `spec`, `dataset`, `tolerance`, ...); run `tests/test_notation.py` after every task. Numpydoc, module docstring on line 1, black 100, ruff correctness only. `h5py` lazy and optional; the header parser and `csvwriter.py` use only the standard library.
- Do not change any existing model, calibration record, workbook or existing CSV; version 1.0 HDF5 files and the VTK importer keep working; existing `fe`/`history` tests stay green. Commits: Conventional Commits, local, no push, no attribution lines.
- Software regression and scientific validation stay separate in every doc: nothing here validates a fatigue method.

## Review Focus

1. Silent convention errors: Abaqus `S13 -> xz`, `S23 -> yz`, Ansys order, engineering shear halving only for strain (Eel, Epl, E), never for stress; tests must use distinct values per component. Task 3, Task 4, Task 6.
2. Truncated, ragged or duplicate grids, unsorted steps, `time` not increasing, `n_nodes`/`n_steps` mismatch. Task 4.
3. Interface nodes (one node, two materials) and node ids repeated across instances: import, selection, assessment. Task 1, Task 4, Task 7.
4. Strain composition and derivation: `Epl` alone, elastic-plastic without a stabilised cycle, derived strain flagged, no derivation in plasticity. Task 4.
5. `assess_nodes`: per-node statuses, `unit_case` on elastic-plastic data, missing or extra load-function columns, a material without a setup. Task 7.

## File Structure

| File | Responsibility |
|---|---|
| `src/multifat/fe/model.py`, `store.py`, `select.py`, `tensors.py` | version 1.1 additions (instance, material, split strain, temperature, metadata, steps), Hooke helper |
| `src/multifat/fe/csvio.py` | header, column resolution with aliases, `import_fe_csv` |
| `src/multifat/fe/csvwriter.py` | stdlib-only `FeCsvWriter` |
| `src/multifat/fe/synthetic.py` | `write_tube_csv`, `--csv` in `main` |
| `src/multifat/history/nodes.py` | `MaterialSetup`, `node_histories`, `read_load_functions`, `assess_nodes` |
| `src/multifat/history/signal.py` | column-name aliases for `from_csv` |
| `tests/test_fecsv_*.py`, `tests/test_hist_nodes.py` | tests per module and an end-to-end file |
| `docs/fe-csv-format.md`, `docs/examples/fe-csv/*.csv`, `docs/fe-input-format.md`, `docs/load-histories.md`, `docs/notation.md`, `docs/decision-log.md`, `AGENTS.md`, `README.md`, `docs/audit-checkpoint.md` | documentation and examples |

---

### Task 1: File format 1.1 (model, store, selection)

**Files:** Modify `src/multifat/fe/model.py`, `store.py`, `select.py`, `tensors.py`, `docs/notation.md`. Create `tests/test_fecsv_store.py`.

**Interfaces:**
- Produces (`model.py`): `AVERAGING` gains `"nodal_average_within_material"`, `"nodal_average_global"`; `HistoryMeta` gains defaults `analysis_type: str = "unspecified"` (`unspecified`, `linear_elastic`, `elastic_plastic`), `cycle: str = "unspecified"`, `step_role: str = "time_sample"` (`time_sample`, `unit_case`), `load_description: str = ""`, validated; `strain_kind` may also be `"derived_linear_elastic"` or `"derived_elastic_part"`; `LoadCases` gains `steps: np.ndarray | None` (integer step values); `HistoryBlock` gains optional `instance`, `material` (`(M,)` string arrays), `E_elastic_hist`, `E_plastic_hist` (`(M, N, 3, 3)`), `T_hist` (`(M, N)`); location uniqueness is over `(location_id, instance, material)` when `instance` or `material` is given, else over `location_id`.
- Produces (`store.py`): `FORMAT_VERSION = "1.1"` (reader accepts major version 1); `HistoryWriter(..., instance=None, material=None, with_strain_split=False, with_temperature=False)` and `write_step(step, S, E, *, E_elastic=None, E_plastic=None, T=None)`; datasets `/locations/instance`, `/locations/material` (UTF-8 strings), `/fields/strain_elastic`, `/fields/strain_plastic` (`(P, N, 6)`), `/fields/temperature` (`(P, N)`), `/load_cases/step`; `HistoryFile` attributes `instance`, `material`, `has_strain_split`, `has_temperature`, and `read` returns the extra fields in the block.
- Produces: `select.select_by_ids(f, ids, *, instance=None, material=None)` returning every matching position when `instance` or `material` is omitted; `tensors.elastic_strain(S, E, nu)` (Hooke, tensor shear, `S` in MPa, `E` in MPa); the unit tables of `to_mpa` and `to_mm` gain `psi`, `ksi`, `N/mm^2` and `in`, `ft`, and a new `to_strain(a, unit)` converts `1, mm/mm, m/m, percent, microstrain` to dimensionless. Two existing tests use `psi` as an example of an unknown unit (`tests/test_fe_tensors.py::test_unit_scaling`, `tests/test_hist_signal.py::test_from_arrays_requires_stress_unit_and_knows_units`): change them to an unambiguous unknown unit such as `furlong`.

- [ ] **Step 1: Write failing tests** (`tests/test_fecsv_store.py`, `importorskip("h5py")`): `test_version_1_1_roundtrip_of_instance_material_split_strain_temperature_and_steps`; `test_same_node_label_in_two_instances_and_an_interface_node_in_two_materials_are_distinct_locations`; `test_duplicate_location_triple_rejected`; `test_node_id_only_duplicates_still_rejected_without_instance_or_material`; `test_version_1_0_file_still_opens` (craft a 1.0 file as in `tests/test_fe_store.py::_craft`); `test_select_by_ids_with_instance_and_material_and_ambiguous_without`; `test_meta_new_fields_roundtrip_and_defaults_for_older_dicts`; `test_averaging_values_validated`; `test_elastic_strain_matches_hooke_by_hand` (uniaxial and shear values); `test_us_customary_and_strain_units_convert` (`psi`, `ksi`, `N/mm^2`, `in`, `ft`, `percent`, `microstrain`, an unknown unit raises).
- [ ] **Step 2: Run** `uv run pytest tests/test_fecsv_store.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the additions; add notation rows (`E_elastic_hist`, `E_plastic_hist`, `T_hist`, `instance`, `material`, `step`, `analysis_type`, `step_role`); run `tests/test_notation.py` and the existing `tests/test_fe_*.py`.
- [ ] **Step 4: Run** all `tests/test_fe_*.py tests/test_fecsv_store.py`. Expected: PASS. black, ruff.
- [ ] **Step 5: Commit** `feat(fe): extend the history file to version 1.1`.

---

### Task 2: Header and column resolution

**Files:** Create `src/multifat/fe/csvio.py` (first part), `tests/test_fecsv_header.py`.

**Interfaces:**
- Produces: `CSV_FORMAT = "multifat-fe-csv 1.0"`; `CsvHeader` frozen dataclass (fields `format, solver, solver_version, stress_unit, length_unit, coordinate_system, averaging, analysis_type, step_role`, optional `strain_convention, strain_unit, cycle, time_unit, temperature_unit, n_nodes, n_steps`, free text `load_description, model, source, notes`, and `extra: dict` for unknown keys); `read_header(path) -> tuple[CsvHeader, int]` returning the header and the number of header lines; value sets and cross rules enforced with `ValueError` naming the key. `COLUMN_ALIASES: dict[str, str]` (lower-case alias to canonical) and `resolve_columns(columns, *, ignore=()) -> ColumnMap` (`ColumnMap.names`: canonical to actual column; `ColumnMap.has_strain`, `has_split`, `has_total`, `has_temperature`, `has_time`, `has_instance`, `has_material`). Canonical: `node_id, instance, material, step, time, x, y, z, S_xx..S_xz, Eel_xx..Eel_xz, Epl_xx..Epl_xz, E_xx..E_xz, T`. Aliases: Abaqus `S11 S22 S33 S12 S23 S13` to `S_xx S_yy S_zz S_xy S_yz S_xz`, `EE..` to `Eel_*`, `PE..` to `Epl_*`, `E11 E22 E33 E12 E23 E13` to `E_*`; Ansys `SX SY SZ SXY SYZ SXZ`, `EPELX EPELY EPELZ EPELXY EPELYZ EPELXZ` to `Eel_*`, `EPPLX..` to `Epl_*`, `EPTOX..` to `E_*`, `TEMP` to `T`; short `Sxx..`, `node` to `node_id`.

- [ ] **Step 1: Write failing tests:** `test_valid_header_parses_all_fields_and_keeps_unknown_keys`; one parametrised `test_missing_required_key_names_the_key` over every required key; `test_value_sets_enforced` (units, `coordinate_system` other than `global_cartesian`, `averaging`, `analysis_type`, `step_role`); `test_strain_convention_and_unit_required_when_strain_columns_exist_is_checked_by_the_importer_contract` (header parse alone does not need them); `test_elastic_plastic_requires_stabilised_cycle`; `test_time_unit_required_with_time_column` (checked in `resolve`/import task); `test_abaqus_names_resolve_with_s13_to_xz_and_s23_to_yz` (distinct sentinel mapping asserted for all six components, strains and plastic strains); `test_ansys_names_resolve`; `test_canonical_and_alias_together_rejected`; `test_unknown_column_rejected_unless_ignored`; `test_missing_required_columns_listed`; `test_partial_strain_group_rejected`; `test_case_insensitive_names`.
- [ ] **Step 2-5:** run (FAIL), implement, run (PASS plus notation), commit `feat(fe): add the CSV header and column resolution`.

---

### Task 3: Tube CSV generator and writer helper

**Files:** Create `src/multifat/fe/csvwriter.py`, `tests/test_fecsv_writer.py`. Modify `src/multifat/fe/synthetic.py`.

**Interfaces:**
- Produces (`csvwriter.py`, standard library only): `FeCsvWriter(path, header: Mapping[str, object], *, columns: Sequence[str])` context manager with `write_row(row: Mapping[str, object])` and `write_rows(rows)`; writes the header lines (`# key: value`, required keys checked), the column line, and floats with 17 significant digits; `REQUIRED_HEADER_KEYS`.
- Produces (`synthetic.py`): `write_tube_csv(path, geom, material, loads, *, n_radii, n_theta, flavour="canonical" | "abaqus" | "ansys", shear="tensor" | "engineering", stress_unit="MPa", length_unit="mm", order="step_major" | "node_major", with_instance=False, with_material=False, with_time=False, strain="split" | "total" | "none", analysis_type="linear_elastic", header=None) -> Path` using `FeCsvWriter`; `Eel` is the Hooke strain and `Epl` a known fraction of it for `elastic_plastic` tests; `main` accepts `--csv DIR` writing `tube_two_steps.csv` (2 radii, 2 angles, two load cases).

- [ ] **Step 1: Write failing tests:** `test_writer_uses_only_the_standard_library` (AST of `csvwriter.py`); `test_writer_header_and_rows_exact_text` (small fixed file compared by text, 17 significant digits); `test_writer_rejects_missing_required_header_key`; `test_write_tube_csv_flavours_use_the_solver_names` (column lines of the three flavours); `test_engineering_shear_doubles_only_strain_columns` (stress unchanged, strain off-diagonals doubled, checked against the tensor-shear oracle); `test_units_scaled` (Pa and m); `test_node_major_order_groups_by_location`; `test_example_flag_writes_the_two_step_example`.
- [ ] **Step 2-5:** run (FAIL), implement, run (PASS), commit `feat(fe): add the CSV writer helper and the tube CSV generator`.

---

### Task 4: Importer

**Files:** Modify `src/multifat/fe/csvio.py`. Create `tests/test_fecsv_import.py`.

**Interfaces:**
- Consumes: Tasks 1-3; `fe.store.HistoryWriter`, `fe.tensors.to_mpa`, `to_mm`, `components_to_sym3`, `elastic_strain`; `fe.model.sha256_of`.
- Produces: `import_fe_csv(csv_path, out_path, *, elastic_constants: Mapping[str, tuple[float, float]] | None = None, allow_global_average: bool = False, ignore_columns: Sequence[str] = (), chunk_rows: int = 200_000, max_node_major_rows: int = 2_000_000) -> HistoryMeta`. Behaviour per the Global Constraints: validate header and columns before reading rows; stream step-major files one step at a time, sort node-major files in memory under the row limit (a larger one is refused with advice to export step-major); check the grid against `n_nodes` and `n_steps` and for completeness; refuse duplicate keys, non-finite or asymmetric tensors, `time` not strictly increasing, steps not ascending in a step-major file; convert units and shear; compose strain by the rules; apply the averaging rule; remove the partial output on failure; record in the meta `solver`, `solver_version`, units, `source_frame="global_cartesian"`, `averaging`, `analysis_type`, `cycle`, `step_role`, `load_description`, the file's SHA-256 in `source_files`, and `source_conventions` (strain convention, alias mapping used, component order); `LoadCases` carry `time` (if present) and `steps`.

- [ ] **Step 1: Write failing tests** (`importorskip("h5py")`; the tube oracle `tube_cartesian` and its strain are the expectation): `test_roundtrip_matches_oracle` parametrised over flavour, shear, units, order, with and without instance/material/time, strain `split`/`total`/`none`; `test_engineering_shear_not_applied_to_stress`; `test_instance_repeated_labels_and_interface_node_import_as_distinct_locations`; `test_truncated_file_refused_against_n_nodes_and_n_steps`; `test_ragged_grid_refused`; `test_duplicate_keys_refused`; `test_unsorted_steps_in_step_major_file_refused`; `test_time_not_increasing_refused_and_never_resorted`; `test_non_contiguous_steps_accepted_and_recorded`; `test_asymmetric_and_nan_tensors_refused`; `test_missing_strain_convention_or_strain_unit_refused_when_strain_present`; `test_strain_in_percent_and_microstrain_converted`; `test_psi_and_inch_file_roundtrips`; `test_non_global_frame_refused`; `test_elastic_plastic_requires_stabilised_cycle`; `test_total_strain_is_eel_plus_epl`; `test_epl_alone_needs_elastic_constants_and_then_derives_the_elastic_part` (flagged `derived_elastic_part`, equals Hooke by hand); `test_linear_elastic_without_strain_derives_only_with_constants_and_is_flagged`; `test_elastic_plastic_without_any_strain_refused`; `test_global_average_across_materials_refused_unless_allowed_and_recorded`; `test_partial_output_removed_on_failure`; `test_node_major_row_limit_gives_advice`; `test_streaming_peak_memory_independent_of_node_count` (tracemalloc, fixed step size, `chunk_rows` small); `test_provenance_records_hash_header_and_alias_mapping`.
- [ ] **Step 2: Run** `uv run pytest tests/test_fecsv_import.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `import_fe_csv` (a first pass over the key columns to validate the grid and collect the location table, a second pass to stream the values; builder's choice of algorithm).
- [ ] **Step 4: Run** the same command plus earlier fe tests. Expected: PASS.
- [ ] **Step 5: Commit** `feat(fe): add the streaming CSV importer`.

---

### Task 5: Aligned single-location CSV and the example

**Files:** Modify `src/multifat/history/signal.py`. Create `docs/examples/fe-csv/tube_two_steps.csv` (generated), tests in `tests/test_hist_signal.py` and `tests/test_fecsv_import.py`.

- [ ] **Step 1: Write failing tests:** `test_from_csv_accepts_the_node_csv_names` (`time`, `S_xx..S_xz`, `Eel_*`/`E_*` as total strain) alongside the old `t`, `sxx`; `test_committed_example_equals_the_generated_one` (`python -m multifat.fe.synthetic --csv` output equals the committed file byte for byte); `test_committed_example_imports_and_matches_the_oracle`.
- [ ] **Step 2-5:** run (FAIL), implement the alias lists in `TimeHistory.from_csv` (the first present name wins, both present is an error), generate the example with the command, run (PASS), commit `feat(fe): align time-history CSV names and add the example`.

---

### Task 6: Node histories, superposition and node assessment

**Files:** Create `src/multifat/history/nodes.py`, `tests/test_hist_nodes.py`.

**Interfaces:**
- Consumes: `fe.store.HistoryFile`, `history.signal.TimeHistory`, `history.assess.assess_history`, `history.accumulate.SNCurve`, `history.eligibility.IneligibleHistoryError`.
- Produces: `MaterialSetup(model, params, curve, E=None, nu=None)`; `NodeRef(position, location_id, instance, material)`; `read_load_functions(path, steps) -> tuple[np.ndarray, np.ndarray, str]` returning the axis, `F` of shape `(Nt, K)` and the axis kind, columns `time` or `index` and `f_<step>` matched by name (missing or extra function columns refused); `node_histories(path, *, positions=None, load_functions=None, chunk_locations=4096) -> Iterator[tuple[NodeRef, TimeHistory]]` (time axis from `time` else the step index with `axis_kind="cycle"`; `unit_case` files need `analysis_type == "linear_elastic"` and load functions, superposed as `S(t) = sum_k f_k(t) S_k`, strain alike); `assess_nodes(path, materials, *, R_list, positions=None, load_functions=None, variant="", on_refusal="record" | "raise", **assess_options) -> list[dict]` with keys `position, location_id, instance, material, status, damage_per_block, blocks_to_failure, cycles_to_failure, n_events, approximation, notes` (`materials` is a `MaterialSetup` for all nodes or a mapping by material label).

- [ ] **Step 1: Write failing tests** (`importorskip("h5py")`; files made with `write_tube_csv` and `import_fe_csv`): `test_two_steps_give_one_cycle_and_the_constant_amplitude_value` (the node result equals `assess_history` on the two-sample history and equals the direct two-state evaluation); `test_time_sample_history_equals_direct_assess_history_for_every_node`; `test_axis_is_time_when_present_else_step_index_as_cycle`; `test_unit_case_superposition_equals_the_directly_built_history`; `test_unit_case_on_elastic_plastic_refused`; `test_load_function_columns_matched_by_name_missing_and_extra_refused`; `test_per_node_status_records_ineligible_and_continues`; `test_on_refusal_raise_aborts_at_the_first`; `test_material_without_setup_gets_status_and_run_continues`; `test_interface_node_assessed_once_per_material_and_repeated_ids_across_instances_distinct`; `test_positions_subset_and_chunking_do_not_change_results`; `test_nonproportional_node_refused_then_allowed_with_the_pass_through_flag`.
- [ ] **Step 2: Run** `uv run pytest tests/test_hist_nodes.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `nodes.py`; add notation rows (`MaterialSetup` vocabulary).
- [ ] **Step 4: Run** the same command and `tests/test_hist_*.py tests/test_notation.py`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(history): assess the nodes of an FE import`.

---

### Task 7: Documentation, isolation and final gate

**Files:** Create `docs/fe-csv-format.md`. Modify `docs/fe-input-format.md`, `docs/load-histories.md`, `docs/notation.md`, `docs/decision-log.md` (`d086` onward), `AGENTS.md`, `README.md`, `docs/audit-checkpoint.md`. Extend `tests/test_fe_isolation.py` and `tests/test_hist_isolation.py`.

- [ ] **Step 1: Write failing tests:** `test_fe_csv_modules_import_no_workbook_calibration_or_history_code` (`csvio.py` and `csvwriter.py`: nothing from `multifat.imp`, `eval`, `models`, `history`, `cli`, `report`); `test_header_parser_and_writer_need_no_h5py_or_numpy_at_import_of_the_writer`; `test_doc_lists_every_header_key_alias_and_the_example_command` (`docs/fe-csv-format.md` mentions each required header key, every alias group, `python -m multifat.fe.synthetic --csv`).
- [ ] **Step 2: Write** `docs/fe-csv-format.md`: purpose and the layered picture, the header (required, conditional, free keys), the columns and aliases, the rules (named columns, shear convention, global frame, averaging per material, node identity, units, strains and derivation, steps and `step_role`), size guidance (node sets, step-major order, row limits, Parquet later), **solver checklists** for Ansys and Abaqus (what to set so the export is global, averaged within materials and in the declared convention; stated as checklists, not tested code, with a request for tested scripts), the writer helper, the example, limits and the later list, questions for the maintainer. Update the other docs; decision-log entries `d086`-`d092` (CSV as user format; location triple; strain composition and derivation; step order and `step_role`; `assess_nodes` statuses and setups; header discipline and integrity keys; writer helper and checklists); checkpoint section.
- [ ] **Step 3: Run the gate:** `uv run pytest -q`, `uv run black --check . && uv run ruff check .`, `uv run python tools/check_docs.py`, `uv run python tools/check_scientific_reviews.py`, `uv run python tools/generate_method_docs.py --check`, `uv run python tools/generate_workbook_docs.py --check`. Expected: all pass. `git diff origin/main --stat -- src/multifat ':!src/multifat/fe' ':!src/multifat/history'` shows nothing.
- [ ] **Step 4: Regression snapshot** on UniTN as before (bundle, `scp -O`, Task Scheduler, `--skip-calibrate`), `origin/main` against the branch; record in the checkpoint.
- [ ] **Step 5: Commit** `docs(fe): publish the CSV exchange format, decisions and example`.

---

## Self-review

- **Spec coverage:** goals 1-9 map to Tasks 7 (specification, docs), 3 and 5 (example, writer, alignment), 2 and 4 (importer, aliases, rules), 1 (format 1.1, derivation helper), 6 (step role, superposition, node assessment), 4 and 6 (tests of refusals, streaming, end to end); decisions from the builder check are in the Global Constraints and tested in Tasks 4 and 6.
- **Types:** `CsvHeader`, `ColumnMap`, `MaterialSetup`, `NodeRef`, the `HistoryBlock` and `HistoryMeta` additions are defined once and used with the same names.
- **Known limits:** Ansys and Abaqus exporters are checklists; the maintainer's tested scripts and one real file are requested. The first-harmonic reduction of a sampled non-proportional cycle is not in this plan.
