# FE local-history interface Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a solver-independent local-history object, an HDF5 "universal" file, a VTK XML (`.vtu`/`.pvd`) importer, a synthetic thick-tube generator that serves as test oracle, bounded-memory selection, range strategies and a bridge to the existing two-state models, without changing any existing number.

**Architecture:** One new package `src/multifat/fe/`, imported by nothing that exists today. Pure-numpy modules (tensor conventions, data model, ranges, bridge, VTK reader, tube generator) do not need `h5py`; only `store.py` does, imported lazily. Tests use the analytical tube as oracle.

**Tech Stack:** Python 3.10+, numpy, standard-library `xml.etree`, optional `h5py` (extra `fe`), pytest, black (100), ruff, `uv`.

**Spec:** `docs/superpowers/specs/2026-11-05-fe-local-history-interface-design.md` (read it first; this plan only records decisions the spec leaves to the builder).

## Global Constraints

- Units: MPa and mm. Strain is dimensionless **tensor** strain (shear = engineering shear / 2). Stress/strain stored as six components in the fixed order `xx, yy, zz, xy, yz, xz`, shape `(P, N, 6)`, float64 default.
- Strain is the **total mechanical strain** (elastic plus plastic); attribute `strain_kind = "total_mechanical"`.
- History is global Cartesian, ordered load cases, straight segments between consecutive ones; `local_only = True`.
- `h5py` is an optional extra `fe` (`h5py>=3.10`), imported lazily with a clear error; the package must install and run without it.
- Notation is data: new public parameter names get a row in `docs/notation.md` before use; `E` is Young's modulus, strain tensors carry the `E_` prefix; do not use retired synonyms.
- Numpydoc docstrings, module docstring on line 1 (before `from __future__`, no leading newline); black line length 100; ruff correctness-only.
- No GPL code (SeAn-CP is GPL-3.0; the repo is MIT, `d071`). The tube is written from the textbook Lame/Saint-Venant equations.
- Do not edit any existing file under `src/multifat/` or any workbook, calibration record, archive or CSV. New code is not imported by existing modules.
- Commits: Conventional Commits, local only on branch `feat/fe-input-interface`, never push, no `Co-authored-by` or AI attribution lines. Confirm with the user before the first commit.
- Keep software regression and scientific validation separate in every doc and docstring: nothing here validates a fatigue method.

## Review Focus

Failure modes the spec implies but a straight reading of the tasks would not test; each has its test in the owning task.

1. Duplicate or mismatched location IDs between load steps (silently pairs wrong points). Task 3, Task 5.
2. `N = 1` and constant histories (zero range, tie-breaking, no NaN). Task 6.
3. Missing strain flowing through the bridge into a stress-only and a strain-needing model. Task 7.
4. Component-order permutation and engineering-vs-tensor shear going unnoticed because the test tensor is too symmetric. Task 1, Task 5.
5. Unsupported VTK encoding (zlib-compressed, appended) misread instead of rejected. Task 5.

## File Structure

| File | Responsibility |
|---|---|
| `pyproject.toml`, `uv.lock`, `.github/workflows/quality.yml` | extra `fe`; CI installs it for the test job |
| `src/multifat/fe/__init__.py` | exports only; no `h5py` import |
| `src/multifat/fe/tensors.py` | component orders, shear conversion, unit scaling, frame rotation |
| `src/multifat/fe/model.py` | `HistoryMeta`, `SourceFile`, `LoadCases`, `HistoryBlock` with validation |
| `src/multifat/fe/store.py` | HDF5 streaming writer, reader, `iter_chunks` |
| `src/multifat/fe/select.py` | selection helpers and screening |
| `src/multifat/fe/synthetic.py` | tube solution, universal-file writer, `.vtu`/`.pvd` writer, `main` |
| `src/multifat/fe/vtk.py` | `.vtu`/`.pvd` reader and `import_pvd` |
| `src/multifat/fe/ranges.py` | `RangeStrategy`, `ExplicitPair`, `LongestChord` |
| `src/multifat/fe/bridge.py` | `to_case_inputs`, eligibility |
| `tests/test_fe_*.py` | one file per module above, plus `test_fe_isolation.py` |
| `docs/fe-input-format.md`, `docs/notation.md`, `docs/decision-log.md`, `AGENTS.md`, `README.md`, `docs/audit-checkpoint.md` | documentation |

---

### Task 1: Branch, extra, tensor conventions

**Files:** Create `src/multifat/fe/__init__.py`, `src/multifat/fe/tensors.py`, `tests/test_fe_tensors.py`. Modify `pyproject.toml`, `.github/workflows/quality.yml`, `uv.lock`.

**Interfaces:**
- Produces (`tensors.py`):
  - `CANONICAL = ("xx","yy","zz","xy","yz","xz")`, `VTK9 = ("xx","xy","xz","yx","yy","yz","zx","zy","zz")`
  - `sym3_to_six(T: ndarray[..., 3, 3]) -> ndarray[..., 6]`, `six_to_sym3(v: ndarray[..., 6]) -> ndarray[..., 3, 3]`
  - `components_to_sym3(a: ndarray[M, k], order: Sequence[str], *, kind: Literal["stress","strain"], shear: Literal["tensor","engineering"] = "tensor") -> ndarray[M, 3, 3]` (engineering shear divided by 2 for strain only; for a 9-component order the two off-diagonal copies must agree within tolerance, else `ValueError`)
  - `to_mpa(a, unit: str)`, `to_mm(a, unit: str)` for `Pa, kPa, MPa, GPa` and `m, mm, cm`; unknown unit raises `ValueError`
  - `rotate_to_global(T: ndarray[M,3,3], R: ndarray[M,3,3]) -> ndarray[M,3,3]` computing `R T R^T` (columns of `R` are the local axes in global coordinates)
  - `cylindrical_rotations(points: ndarray[M,3], axis=(0,0,1), origin=(0,0,0)) -> ndarray[M,3,3]` with columns `e_r, e_theta, e_z`; a point on the axis raises `ValueError`

- [ ] **Step 1:** Create branch `feat/fe-input-interface` from `main`. Add `[project.optional-dependencies] fe = ["h5py>=3.10"]` to `pyproject.toml`; run `uv lock`; change the test job's sync in `quality.yml` to `uv sync --locked --extra fe`. Expected: `uv lock` succeeds; if it cannot reach the network, stop and report.
- [ ] **Step 2: Write failing tests** in `tests/test_fe_tensors.py`: `test_six_roundtrip` (random symmetric tensors); `test_component_order_detects_permutation` (a tensor with six distinct values, e.g. `xx=1, yy=2, zz=3, xy=4, yz=5, xz=6`, read with the wrong order must differ); `test_engineering_shear_halved_for_strain_only`; `test_vtk9_asymmetric_pair_rejected` (`xy != yx` raises); `test_unit_scaling` (`Pa -> MPa` is `1e-6`); `test_cylindrical_to_cartesian_analytic` (point `(r, 0, 0)` and `(0, r, 0)`: a stress `diag(s_r, s_t, s_z)` in `(r,theta,z)` maps to the known Cartesian diagonal; the invariants trace and von Mises are preserved at 30 random points); `test_point_on_axis_rejected`.
- [ ] **Step 3: Run** `uv run pytest tests/test_fe_tensors.py -q`. Expected: FAIL (module missing).
- [ ] **Step 4: Implement** `tensors.py` as specified. `__init__.py` holds the module docstring only.
- [ ] **Step 5: Run** the same command. Expected: PASS. Then `uv run black --check . && uv run ruff check .`.
- [ ] **Step 6: Commit** `feat(fe): add tensor conventions for FE results`.

---

### Task 2: Data model and notation rows

**Files:** Create `src/multifat/fe/model.py`, `tests/test_fe_model.py`. Modify `docs/notation.md`.

**Interfaces:**
- Consumes: `tensors.sym3_to_six`.
- Produces (`model.py`):
  - `SourceFile(name: str, sha256: str)` frozen dataclass; `sha256_of(path) -> str`
  - `LoadCases(time: ndarray[N] | None, labels: tuple[str, ...] | None)` frozen dataclass
  - `HistoryMeta(material_id, solver, solver_version, analysis_step, adapter, adapter_version, source_files: tuple[SourceFile, ...], source_units: dict[str, str], source_frame: str, averaging: Literal["none","nodal_average","extrapolated"], location_kind: Literal["node","element_centroid","integration_point","element_nodal"], periodic: bool = False, local_only: bool = True, strain_kind: str | None = "total_mechanical", kind: str = "fe_prediction", format_version: str = "1.0", created: str = "", notes: str = "")` frozen dataclass; `to_dict()` / `from_dict()` JSON-safe
  - `HistoryBlock(location_id: ndarray[M], coordinates: ndarray[M,3], S_hist: ndarray[M,N,3,3], E_hist: ndarray[M,N,3,3] | None)` frozen dataclass whose `__post_init__` raises `ValueError` for wrong shapes, non-finite values, `N < 1`, duplicate IDs, or asymmetry beyond `1e-9 * max|T|`

- [ ] **Step 1: Write failing tests** `test_fe_model.py`: `test_block_accepts_valid`; `test_block_rejects_nan`; `test_block_rejects_asymmetric`; `test_block_rejects_duplicate_ids`; `test_block_rejects_shape_mismatch_between_S_and_E`; `test_meta_dict_roundtrip`; `test_sha256_of_known_bytes` (known digest of `b"abc"` written to a temp file).
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_model.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `model.py`. Add rows to `docs/notation.md` for `S_hist`, `E_hist`, `S_range`, `E_range`, `location_id`, `pair`; run `uv run pytest tests/test_notation.py -q`.
- [ ] **Step 4: Run** both test files. Expected: PASS.
- [ ] **Step 5: Commit** `feat(fe): add history data model and notation rows`.

---

### Task 3: HDF5 store, selection, bounded iteration

**Files:** Create `src/multifat/fe/store.py`, `src/multifat/fe/select.py`, `tests/test_fe_store.py`, `tests/test_fe_select.py`.

**Interfaces:**
- Consumes: `model.HistoryMeta`, `LoadCases`, `HistoryBlock`; `tensors.sym3_to_six`, `six_to_sym3`.
- Produces (`store.py`):
  - `HistoryWriter(path, meta: HistoryMeta, location_id: ndarray[P], coordinates: ndarray[P,3], load_cases: LoadCases, n_load_cases: int, *, with_strain: bool, dtype: str = "float64", chunk_locations: int = 4096)`, a context manager; `write_step(step: int, S: ndarray[P,3,3], E: ndarray[P,3,3] | None) -> None` writes the slab `[:, step, :]`; leaving the context without all steps written raises `ValueError`
  - `HistoryFile(path)`, a context manager with attributes `meta`, `load_cases`, `location_id`, `coordinates`, `n_locations`, `n_load_cases`, `has_strain`; `read(index: ndarray | slice) -> HistoryBlock`; `iter_chunks(*, max_locations: int | None = None, max_bytes: int | None = None, index: ndarray | None = None) -> Iterator[HistoryBlock]`
  - file layout and names: the spec's, with root attributes `format = "multifat-fe-history"` and `format_version`; reserved empty group `/mesh`; datasets `/fields/stress`, `/fields/strain`, chunked `(chunk_locations, N, 6)`
  - `require_h5py()` raising `ImportError("... install multifat[fe]")` when absent
- Produces (`select.py`): `select_by_ids(f: HistoryFile, ids) -> ndarray[int]`, `select_by_mask(f, mask) -> ndarray[int]`, `select_in_box(f, lo, hi) -> ndarray[int]`, `von_mises_max(block: HistoryBlock) -> ndarray[M]`, `screen_by(f, scalar: Callable[[HistoryBlock], ndarray], *, top: int | None = None, threshold: float | None = None, max_locations: int = 4096) -> ndarray[int]`; all return sorted unique positions, and unknown IDs raise `KeyError`

- [ ] **Step 1: Write failing tests** (start each with `pytest.importorskip("h5py")`): `test_roundtrip_equals_input` (write 5 locations x 4 steps stepwise, read back, `array_equal` in float64); `test_float32_storage_within_tolerance`; `test_missing_step_raises_on_close`; `test_provenance_and_sha_roundtrip`; `test_duplicate_location_ids_rejected_by_writer`; `test_iter_chunks_covers_all_locations_once`; `test_iter_chunks_empty_selection_yields_nothing`; `test_peak_memory_independent_of_P` (build files with `P = 2000` and `P = 8000`, same `N`, iterate with `max_locations=100` under `tracemalloc`; peaks differ by less than 25 %); `test_without_h5py_message` (monkeypatch `require_h5py` import to fail; message names the extra). In `test_fe_select.py`: `test_select_by_ids_unknown_raises`, `test_select_in_box`, `test_screen_by_top_k_matches_bruteforce`.
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_store.py tests/test_fe_select.py -q`. Expected: FAIL. If `h5py` is missing locally, `uv sync --extra fe` first.
- [ ] **Step 3: Implement** both modules. `iter_chunks` derives `max_locations` from `max_bytes` using the per-location size `N * 6 * 8 * (1 + has_strain)` times 3 for the expanded blocks.
- [ ] **Step 4: Run** the same command. Expected: PASS; `uv run black --check . && uv run ruff check .`.
- [ ] **Step 5: Commit** `feat(fe): add HDF5 history store with chunked reading and selection`.

---

### Task 4: Thick-tube solution and universal-file writer

**Files:** Create `src/multifat/fe/synthetic.py` (first part), `tests/test_fe_synthetic.py`.

**Interfaces:**
- Consumes: `store.HistoryWriter`, `model.*`, `tensors.cylindrical_rotations`, `rotate_to_global`.
- Produces (`synthetic.py`):
  - `TubeGeometry(a: float, b: float)`, `TubeMaterial(E: float, nu: float)`, `TubeLoad(p_i: float, F: float, T: float)` frozen dataclasses (MPa, N, N mm; inner radius `a`, outer `b`, mm); outer pressure is zero
  - `tube_stress_cyl(r: ndarray, load: TubeLoad, geom: TubeGeometry, *, closed_ends: bool = True) -> ndarray[..., 3, 3]` in `(r, theta, z)`: Lame `A = p_i a^2 / (b^2 - a^2)`, `B = a^2 b^2 p_i / (b^2 - a^2)`, `s_rr = A - B/r^2`, `s_tt = A + B/r^2`, `s_zz = F/area + (A if closed_ends else 0)`, `t_tz = T r / J`
  - `tube_strain_cyl(stress: ndarray, mat: TubeMaterial) -> ndarray` (isotropic Hooke, tensor shear)
  - `tube_points(geom, n_radii: int, n_theta: int, z: float = 0.0) -> tuple[ndarray, ndarray, ndarray, ndarray]` returning `location_id, coordinates[P,3], r, theta`
  - `tube_cartesian(points..., load, geom, mat, *, closed_ends=True) -> tuple[ndarray, ndarray]` global Cartesian `S, E`
  - `write_tube_history(path, geom, mat, loads: Sequence[TubeLoad], *, n_radii: int, n_theta: int, closed_ends: bool = True, dtype: str = "float64") -> None`, streaming one load case at a time through `HistoryWriter`

- [ ] **Step 1: Write failing tests** `test_fe_synthetic.py`: independent oracles, not reproductions of the formulas: `test_boundary_conditions` (`s_rr(a) == -p_i`, `s_rr(b) == 0` to `1e-9`); `test_radial_equilibrium` (finite-difference `ds_rr/dr + (s_rr - s_tt)/r` is small, tolerance `1e-5 * p_i`); `test_torque_resultant` (numerical `int t_tz * r dA` over the section equals `T` within `1e-3`); `test_axial_resultant` (`int s_zz dA == F + A * area` for closed ends); `test_hooke_plane_stress_invariants` (trace of strain equals `(1 - 2 nu)/E` times trace of stress); `test_cartesian_rotation_preserves_invariants` (von Mises identical in both frames); `test_write_tube_history_roundtrip` (importorskip `h5py`; file contents equal `tube_cartesian` output).
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_synthetic.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the functions above from the stated equations.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: Commit** `feat(fe): add analytical thick-tube generator`.

---

### Task 5: VTK XML writer for fixtures and importer

**Files:** Modify `src/multifat/fe/synthetic.py` (add the writer). Create `src/multifat/fe/vtk.py`, `tests/test_fe_vtk.py`.

**Interfaces:**
- Consumes: Tasks 1-4.
- Produces (`synthetic.py`): `write_tube_pvd(directory, geom, mat, loads, *, n_radii, n_theta, encoding: Literal["ascii","base64"] = "ascii", frame: Literal["cartesian","cylindrical"] = "cartesian", order: Sequence[str] = CANONICAL, shear: Literal["tensor","engineering"] = "tensor", stress_unit: str = "MPa", length_unit: str = "mm", closed_ends: bool = True) -> Path` returning the `.pvd` path; one `.vtu` per load case, a point cloud of `VTK_VERTEX` cells, point arrays `Stress`, `Strain`, `NodeID` (and `Rotation` is not written: the cylindrical frame is declared to the importer)
- Produces (`vtk.py`):
  - `CylindricalFrame(axis: tuple[float,float,float] = (0,0,1), origin: tuple[float,float,float] = (0,0,0))`, `RotationArray(name: str)`
  - `TensorSpec(array: str, order: Sequence[str], unit: str = "MPa", shear: Literal["tensor","engineering"] = "tensor", location: Literal["point","cell"] = "point")`
  - `read_vtu(path) -> VtuData` (`points[P,3]`, `point_data: dict`, `cell_data: dict`, `cell_centroids[C,3]`); supports `ascii` and inline `binary` (base64, `UInt32` header, uncompressed); any `compressor`, `appended` data or other encoding raises `ValueError` naming it
  - `read_pvd(path) -> list[tuple[float, Path]]` in timestep order
  - `import_pvd(pvd_path, out_path, *, stress: TensorSpec, strain: TensorSpec | None, id_array: str | None, length_unit: str = "mm", frame: Literal["cartesian"] | CylindricalFrame | RotationArray = "cartesian", material_id: str, solver: str, solver_version: str, analysis_step: str, averaging: str = "none", location_kind: str = "node", periodic: bool = False) -> HistoryMeta` streaming one `.vtu` at a time into `HistoryWriter`; records the SHA-256 of every source file, the declared frame and units; checks that IDs and coordinates agree across steps and raises otherwise

- [ ] **Step 1: Write failing tests** `test_fe_vtk.py` (reader tests need no `h5py`): parametrised `test_reader_recovers_oracle` over `encoding in {ascii, base64}`, `frame in {cartesian, cylindrical}`, `order in {CANONICAL, VTK9}`, `shear in {tensor, engineering}`, `unit/length in {(MPa, mm), (Pa, m)}`: read the steps, convert with `components_to_sym3`, `to_mpa`, rotate, and compare to `tube_cartesian` within `1e-9` relative; `test_import_pvd_to_universal_file_matches_oracle` (importorskip `h5py`; whole path through `import_pvd` then `HistoryFile.read`); `test_cell_data_uses_centroids`; `test_rejects_compressed_vtu` (hand-edit a fixture header to `vtkZLibDataCompressor`; expect `ValueError` mentioning it); `test_rejects_appended_data`; `test_mismatched_ids_between_steps_raises`; `test_mismatched_point_count_between_steps_raises`; `test_missing_array_name_lists_available`; `test_provenance_records_sha256_of_each_file`.
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_vtk.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the writer and reader. Use `xml.etree.ElementTree` and `numpy.frombuffer`; the reader must not guess: unknown `type`, `Name` or `NumberOfComponents` mismatches raise.
- [ ] **Step 4: Run** the same command. Expected: PASS; run black and ruff.
- [ ] **Step 5: Commit** `feat(fe): add VTK XML importer with explicit tensor conventions`.

---

### Task 6: Range strategies

**Files:** Create `src/multifat/fe/ranges.py`, `tests/test_fe_ranges.py`.

**Interfaces:**
- Produces: `RangeResult(S_range: ndarray[3,3], E_range: ndarray[3,3] | None, pair: tuple[int, int])`; `RangeStrategy` protocol `__call__(S_hist: ndarray[N,3,3], E_hist: ndarray[N,3,3] | None) -> RangeResult`; `ExplicitPair(a: int, b: int)` (range `X[a] - X[b]`, indices validated); `LongestChord(metric: Literal["deviatoric_norm","principal_spread"] = "deviatoric_norm", on: Literal["stress","strain"] = "stress")` returning the pair `(a, b)` with `a < b` maximising the metric of `X[a] - X[b]`, ties to the lexicographically smallest pair; `N = 1` returns pair `(0, 0)` and a zero range; a constant history returns pair `(0, 1)` and a zero range. Memory is `O(N)` (loop over `a`).

- [ ] **Step 1: Write failing tests** `test_fe_ranges.py`: `test_explicit_pair_range_and_bounds`; `test_longest_chord_two_states_returns_that_pair`; `test_longest_chord_matches_bruteforce` (random `N = 12`, both metrics, both `on` values); `test_single_load_case_zero_range`; `test_constant_history_zero_range_no_nan`; `test_tie_break_is_smallest_pair`; `test_longest_chord_matches_d069_closed_form` (pure-shear `xz`, `yz` ellipse sampled at 720 points; the deviatoric Frobenius norm divided by `sqrt(2)` equals `2 sqrt(max(0, 0.5 (aa + bb + hypot(aa - bb, 2 ab))))` of `multifat.eval.evaluator` within `1e-4` relative, where `aa, bb, ab` are the squared norms and dot of the sine and cosine shear vectors); `test_strain_missing_with_on_strain_raises`.
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_ranges.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the module; document the metrics as provisional, with scientific review and the remaining rules in issue #20.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: Commit** `feat(fe): add range strategies for load-case histories`.

---

### Task 7: Bridge to existing models and eligibility

**Files:** Create `src/multifat/fe/bridge.py`, `tests/test_fe_bridge.py`.

**Interfaces:**
- Consumes: `model.HistoryBlock`, `ranges.RangeStrategy`; existing `multifat.models.registry.get_models`, `stress_only_models`, `multifat.imp.experimental_data.extract_tensor`, `multifat.planes.sampling.build_R_list`.
- Produces: `IneligibleMethodError(ValueError)`; `CaseInputs(location_id, pair, S0, E0, S1, E1, harmonics: dict[str, ndarray])` frozen dataclass; `check_eligibility(meta: HistoryMeta, method) -> None` raising when `meta.local_only` and `getattr(method, "requires_nonlocal", False)`; `to_case_inputs(block: HistoryBlock, strategy: RangeStrategy, *, meta: HistoryMeta | None = None, method=None) -> list[CaseInputs]`. Convention copied from `extract_tensor`: `S0 = S_hist[a]`, `S1 = S_hist[b]`, `S_mean = (S0 + S1)/2`, `S_sin = (S0 - S1)/2`, `S_cos = 0`, likewise for `E`. Without strain, `E0`, `E1` and all `E_*` are `NaN` `(3,3)` arrays, the convention `stress_only_models()` probes with.

- [ ] **Step 1: Write failing tests** `test_fe_bridge.py`: `test_bridge_equals_extract_tensor_convention` (build the nested workbook dict `data[mat]["cases"][0][1][1|2]["stress"|"strain"]` by hand from the same two tensors; `extract_tensor(..., include_harmonics=True)` equals the bridge outputs exactly); `test_two_state_history_runs_existing_model_identically` (`get_models(["FS"])[0].evaluate_case` on bridge output and on `extract_tensor` output, with `params = {"Sy": 700.0, "k_FS": 1.0}` and `R_list = build_R_list(n1=3, n2=6)`, gives equal `values`); `test_stress_only_model_runs_without_strain` (pick one name from `stress_only_models()`, finite value) and `test_strain_model_without_strain_fails_clearly` (a model outside that set; assert the observed error type, not silence); `test_stub_requires_nonlocal_refused_with_reason`; `test_all_registered_methods_pass_eligibility` (`get_models` over every name accepted by the registry, skipping disabled ones via `ensure_enabled`); `test_nonproportional_reduction_is_documented_lossy` (a quarter-period out-of-phase history: the two-state pair differs from the harmonic `S_sin`/`S_cos` description; the assertion encodes the stated limit).
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_bridge.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the module.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: Commit** `feat(fe): add bridge from local histories to case-level models`.

---

### Task 8: Isolation and type separation

**Files:** Create `tests/test_fe_isolation.py`.

- [ ] **Step 1: Write tests:** `test_no_existing_module_imports_fe` (walk `src/multifat` with `ast`; no module outside `fe/` imports `multifat.fe`); `test_fe_imports_no_calibration_or_workbook_modules` (modules in `fe/` import nothing from `multifat.imp`, `multifat.eval` or `multifat.report`; only `bridge.py` may import from `multifat.models`, and only for type hints); `test_workbook_and_calibration_entry_points_reject_fe_objects` (pass a `HistoryBlock` and an `.h5` path to `load_material_workbook`, `import_experimental_data_grouped_by_material` and the comparison reader; assert the exact error each raises today, not silence; if any entry silently accepts, add the smallest type guard in that function and report it to the user).
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_isolation.py -q`. Fix to PASS.
- [ ] **Step 3: Verify no existing source changed:** `git diff main --stat -- src/multifat ':!src/multifat/fe'` prints nothing (unless a guard was added in step 2, in which case list it for the user).
- [ ] **Step 4: Commit** `test(fe): guard isolation from workbook and calibration paths`.

---

### Task 9: Large-input test

**Files:** Create `tests/test_fe_large.py`.

- [ ] **Step 1: Write tests:** `test_stream_many_locations_bounded` (default run: `P = 40_000`, `N = 8`, `write_tube_history` then `iter_chunks(max_bytes=4_000_000)` under `tracemalloc`; peak below 3x the budget; every location seen once); `test_large_scale_opt_in` (skipped unless `MULTIFAT_SLOW=1`: `P = 1_000_000`, `N = 16`, `chunk` budget 64 MB, peak below 4x budget, reads every location).
- [ ] **Step 2: Run** `uv run pytest tests/test_fe_large.py -q` and, once, `MULTIFAT_SLOW=1 uv run pytest tests/test_fe_large.py -q`. Expected: PASS; report timings and disk use of the opt-in run.
- [ ] **Step 3: Commit** `test(fe): add bounded-memory large-input tests`.

---

### Task 10: Documentation, decision log and final gate

**Files:** Create `docs/fe-input-format.md`. Modify `docs/decision-log.md`, `docs/notation.md` (if rows are missing), `AGENTS.md` (layout line for `fe/`), `README.md` (one paragraph and link), `docs/audit-checkpoint.md`. Add `python -m multifat.fe.synthetic --out DIR` via `main(argv=None) -> int` in `synthetic.py`, writing the universal file and a `.pvd` series.

- [ ] **Step 1: Write failing test** `test_documented_example_runs` in `tests/test_fe_synthetic.py`: runs `main(["--out", tmp])` and checks that the universal file and the `.pvd` exist and import back.
- [ ] **Step 2: Implement** `main`; write `docs/fe-input-format.md` by hand: purpose, contract table, HDF5 layout, conventions (units, tensor shear, component order, frame), VTK importer arguments and supported encodings, limits (one load step must fit in memory; straight-segment assumption; local-only; two-state bridge is lossy for non-proportional paths), the example command, and the follow-up list. Add decision-log entries (`d072` onward, one per spec decision: scope boundary, HDF5 extra, VTK first, strain meaning, frames, eligibility mechanism, tube from textbook) in the existing table and detailed-record format.
- [ ] **Step 3: Run the gate:** `uv run pytest -q`, `uv run black --check . && uv run ruff check .`, `uv run python tools/check_docs.py`, `uv run python tools/check_scientific_reviews.py`, `uv run python tools/generate_method_docs.py --check`, `uv run python tools/generate_workbook_docs.py --check`. Expected: all pass.
- [ ] **Step 4: Regression snapshot (no existing source changed, so this confirms rather than discovers):** `uv run python tools/regression_snapshot.py run <dir-main> --skip-calibrate` at `main` and `... run <dir-branch> --skip-calibrate` on the branch, then `compare`. This takes a while; ask the user whether to run it locally or on a workstation (see global instructions) before starting. Expected: no differences.
- [ ] **Step 5: Update `docs/audit-checkpoint.md`** (top section: result, files, remaining work, exact next action) and **commit** `docs(fe): publish FE input specification and decisions`.

---

## Self-review

- **Spec coverage:** goals 1-7 map to Tasks 2-3 (object, file), 5 (VTK), 4 (generator), 3 (selection, chunks), 6-7 (ranges, bridge, eligibility), 4-9 (tests); non-goals are in Global Constraints and Task 8; acceptance criteria 1-7 are Tasks 10, 5, 4-5, 9, 8-10, 8, 7.
- **Types:** `HistoryBlock`, `HistoryMeta`, `LoadCases`, `HistoryWriter`, `HistoryFile`, `RangeResult`, `CaseInputs` are defined once and used with the same names later.
- **Known limit:** the opt-in large test is not part of CI; the default-run bounded-memory tests are.
