# Load-history bridge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bridge a general (proportional) multiaxial load history to damage parameters and an accumulated life: a time-history contract, scalar rainflow and longest-chord events, damage per event through the existing models, linear Miner on a caller-supplied Basquin curve, eligibility diagnostics, and reproducible results.

**Architecture:** A new package `src/multifat/history/`, imported by nothing existing. It reuses `multifat.fe` (`HistoryBlock`, `ExplicitPair`, `LongestChord`, `to_case_inputs`, `IneligibleMethodError`) and the existing models unchanged. Exact for proportional histories; a non-proportional history is refused unless an explicit, recorded approximation flag is set.

**Tech Stack:** Python 3.10+, numpy, pandas (CSV), optional `h5py` (extra `fe`), pytest, black (100), ruff, `uv`.

**Spec:** `docs/superpowers/specs/2026-11-05-load-history-bridge-design.md` (read it first; it holds the reasoning and the issue's acceptance criteria). Upstream `main` is merged (5 November): criteria may raise `LoadingScopeError`, models may declare `design_loading`.

## Global Constraints

- MPa and mm; tensor shear; symmetric 3x3 tensors; one location per history in this round.
- The history is **one block that repeats until failure**: damage per block `D = sum n_i / N_f,i`, life `1/D` blocks, life in cycles `(sum of counts per block) / D`.
- Default scalar channel: projection of the mean-removed, shear-weighted (`sqrt(2)`) six-component stress path on its dominant singular direction (sign fixed by the largest component).
- Counting: ASTM E1049 three-point rainflow; `residue="half"` (default) or `"rearrange"`; the policy is recorded.
- The S-N curve is **caller-supplied**: `SNCurve(A, b, method, variant, design_loading, source)`, `DP = A N^b`, built from `A_surv` and `b` of the pipeline's in-memory fit; `b >= 0` or non-finite is refused. Life from a damage parameter is `multifat.eval.fat_cycles.fatigue_cycles(DP, A, b)`.
- Native-life methods (`getattr(model, "native_life", None) is not None`) are refused for accumulation; `LoadingScopeError` on an event is reported as `IneligibleHistoryError` naming method and event; strain-needing methods need strain.
- A non-proportional history raises `IneligibleHistoryError` unless `allow_nonproportional_approximation=True`, and the flag is written into the result.
- Notation is data: new public names get rows in `docs/notation.md` before use; `E` is Young's modulus, strain tensors carry `E_`; no retired synonyms (`mat`, `spec`, `dataset`, ... see the guard `tests/test_notation.py`; run it after every task).
- Numpydoc docstrings, module docstring on line 1; black line 100; ruff correctness only; `h5py` stays optional and lazily imported; no GPL code (MIT, `d071`); no equation implemented from memory beyond ASTM E1049 rainflow and Miner.
- Do not edit any existing source file, record, workbook or CSV; `fe` and the rest must not import `history`. Commits: Conventional Commits, local, never push, no `Co-authored-by` or AI attribution.
- Software regression and scientific validation stay separate in every docstring and doc: nothing here validates a fatigue method.

## Review Focus

1. Flat or tiny histories: `N = 2`, constant history, plateaus (equal consecutive values), zero-range reversals: no events, no NaN, no division by zero. Task 3, Task 7.
2. Bad time axes and inputs: non-monotonic or duplicate times, NaN, asymmetric tensors, missing units in a CSV, a generator returning the wrong shape. Task 1.
3. Proportionality verdict: noise at 1e-9 and a static hydrostatic offset must not flip a proportional history; a visibly elliptical path must be refused. Task 2, Task 7.
4. Events with `DP <= 0`, events below the endurance limit, a history with zero total damage (blocks to failure infinite): no `log(0)`, no `inf * 0`. Task 5, Task 7.
5. Curve and model mismatch: other method name, `_ext` requested on a model without it, `LoadingScopeError` mid-history, curve with `b >= 0`. Task 5, Task 6, Task 7.

## File Structure

| File | Responsibility |
|---|---|
| `src/multifat/history/__init__.py` | docstring only |
| `src/multifat/history/signal.py` | `TimeHistory`, constructors, validation, `record` |
| `src/multifat/history/channels.py` | scalar channels |
| `src/multifat/history/proportionality.py` | non-proportionality measure |
| `src/multifat/history/rainflow.py` | `turning_points`, `Cycle`, `rainflow` |
| `src/multifat/history/events.py` | `Event`, events from cycles / longest chord, merging |
| `src/multifat/history/accumulate.py` | `SNCurve`, `Endurance`, `miner` |
| `src/multifat/history/eligibility.py` | `IneligibleHistoryError`, checks, method matrix |
| `src/multifat/history/assess.py` | `assess_history`, `HistoryAssessment` |
| `src/multifat/history/example.py` | `python -m multifat.history.example` |
| `tests/test_hist_*.py` | one file per module, plus `test_hist_isolation.py` |
| `docs/load-histories.md`, `docs/notation.md`, `docs/decision-log.md`, `AGENTS.md`, `README.md`, `docs/audit-checkpoint.md` | documentation |

---

### Task 1: Time-history contract

**Files:** Create `src/multifat/history/__init__.py`, `src/multifat/history/signal.py`, `tests/test_hist_signal.py`. Modify `docs/notation.md`.

**Interfaces:**
- Consumes: `multifat.fe.model.HistoryBlock`, `HistoryMeta`, `SourceFile`; `multifat.fe.tensors.to_mpa`, `components_to_sym3`; `multifat.fe.store.HistoryFile` (lazy).
- Produces (`signal.py`): `TimeHistory(t, S_hist, E_hist, axis_kind, units, provenance)` frozen dataclass (`t` `(N,)` strictly increasing, `N >= 2`; `S_hist` `(N,3,3)` MPa; `E_hist` `(N,3,3)` or `None`; `axis_kind` `"time"` or `"cycle"`; `units` dict with `stress`; `provenance` JSON-safe dict) with
  - `sha256() -> str` over `t`, `S_hist`, `E_hist` as little-endian float64;
  - `to_block(location_id: int = 0) -> HistoryBlock`; `meta() -> HistoryMeta` (`adapter="multifat.history"`, `local_only=True`, material from provenance or `"unspecified"`);
  - `record(path) -> Path` writing `<path>.npz` (arrays) and `<path>.json` (units, provenance, `sha256`); classmethod `from_record(path)`;
  - classmethods `from_arrays(t, S_hist, E_hist=None, *, stress_unit: str, axis_kind="time", provenance=None)`, `from_csv(path, *, stress_unit: str, time_unit: str, delimiter=",", columns=None, strain_shear="tensor")` (default columns `t, sxx, syy, szz, sxy, syz, sxz` and optional `exx, eyy, ezz, exy, eyz, exz`; units are required, never assumed), `from_file(path, location_id=None)` (HDF5 history with one location, needs `h5py`), `from_function(fn, t, *, name=None, **params)` where `fn(t, **params)` returns `(N,3,3)` or a dict `{"S": ..., "E": ...}`; provenance records `{"source": "function", "function": "<module>.<qualname>", "params": params}` and `params` must be JSON-serialisable.

- [ ] **Step 1: Write failing tests** (`tests/test_hist_signal.py`): `test_rejects_non_increasing_or_duplicate_times`; `test_rejects_nan_and_asymmetric_tensors`; `test_rejects_single_sample`; `test_from_arrays_requires_stress_unit` (omitting it is a `TypeError`; `"psi"` raises `ValueError`); `test_stress_unit_pa_is_converted_to_mpa`; `test_from_csv_roundtrip_and_missing_unit_rejected` (write with pandas, read back equal; engineering-shear strain halved); `test_from_function_records_name_and_params`; `test_from_function_wrong_shape_rejected`; `test_from_function_non_json_params_rejected`; `test_record_roundtrip_reproduces_arrays_and_hash` (`from_record(record(p))` has equal `sha256()`); `test_to_block_matches_fe_contract` (`HistoryBlock` shapes `(1, N, 3, 3)`); `test_from_file_reads_one_location_history` (`importorskip("h5py")`, written with `fe.synthetic.write_tube_history` using `n_radii=1, n_theta=1`).
- [ ] **Step 2: Run** `uv run pytest tests/test_hist_signal.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** `signal.py`; add notation rows for `t`, `axis_kind`, `TimeHistory` vocabulary (`S_hist`, `E_hist` exist), then `uv run pytest tests/test_notation.py -q`.
- [ ] **Step 4: Run** both test files. Expected: PASS. `uv run black --check . && uv run ruff check .`.
- [ ] **Step 5: Commit** `feat(history): add the time-history contract`.

---

### Task 2: Channels and proportionality

**Files:** Create `src/multifat/history/channels.py`, `src/multifat/history/proportionality.py`, `tests/test_hist_channels.py`.

**Interfaces:**
- Consumes: `signal.TimeHistory`; `multifat.fe.tensors.sym3_to_six`.
- Produces (`channels.py`): `dominant(S_hist) -> ndarray[N]`, `component(S_hist, i: int, j: int)`, `normal_stress(S_hist, n)`, `signed_von_mises(S_hist)` (sign of the trace, sign of the largest principal value where the trace is zero), `resolve_channel(channel, history) -> ndarray[N]` accepting `"dominant"`, `"von_mises"`, `("component", i, j)`, `("normal", n)` or a callable `f(S_hist) -> (N,)`; unknown names raise `ValueError`.
- Produces (`proportionality.py`): `Proportionality(ratio, proportional, tol)` and `proportionality(S_hist, *, tol: float = 1e-6)` where `ratio` is the second over the first singular value of the mean-removed shear-weighted six-component path (0 for a rank-one or constant path) and `proportional = ratio <= tol`.

- [ ] **Step 1: Write failing tests:** `test_dominant_channel_is_affine_in_the_load_factor` (`S(t) = lam(t) A + S_static`, random `A`, nonzero `S_static`, `lam` a sum of sines: a linear fit of the channel on `lam` has residual below 1e-10 and a nonzero slope); `test_dominant_channel_sign_is_fixed` (negating `A` leaves the channel's sign convention unchanged: the largest component of the singular direction is positive); `test_von_mises_is_signed_for_uniaxial_tension_and_compression`; `test_von_mises_zero_trace_uses_principal_sign`; `test_component_and_normal_channels_match_direct_projection`; `test_resolve_channel_forms_and_unknown_rejected`; `test_proportional_with_static_hydrostatic_offset_is_proportional`; `test_noise_1e9_does_not_flip_proportionality`; `test_elliptical_shear_path_is_not_proportional` (`ratio > 0.1`); `test_constant_history_is_proportional_ratio_zero`.
- [ ] **Step 2: Run** `uv run pytest tests/test_hist_channels.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** the two modules; document the sign rules and the tolerance.
- [ ] **Step 4: Run** the same command plus `tests/test_notation.py`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(history): add scalar channels and a proportionality measure`.

---

### Task 3: Rainflow

**Files:** Create `src/multifat/history/rainflow.py`, `tests/test_hist_rainflow.py`.

**Interfaces:**
- Produces: `turning_points(x: ndarray, gate: float = 0.0) -> ndarray[int]` (indices of the endpoints and every local peak or valley, plateaus collapsed to their first sample, reversals smaller than `gate` removed by the usual hysteresis filter); `Cycle(span: float, mean: float, count: float, i_start: int, i_end: int)` (indices into `x`, in time order); `rainflow(x, *, gate=0.0, residue: Literal["half","rearrange"]="half") -> list[Cycle]`: ASTM E1049 three-point algorithm, full cycles `count = 1.0`, residue and start-point half cycles `count = 0.5`; `"rearrange"` rotates the cyclic reversal sequence to start and end at the global extremum of largest magnitude so that no residue remains (every cycle `count = 1.0`).

- [ ] **Step 1: Write failing tests:** `test_astm_e1049_example_sequence` (reversals `-2, 1, -3, 5, -1, 3, -4, 4, -2`: summed counts by span `{3: 0.5, 4: 1.5, 6: 0.5, 8: 1.0, 9: 0.5}`; this is the commonly reproduced example of the standard; cite it in the test and verify the numbers against ASTM E1049 or Lee and Tjhung (Zotero key `A8XFUU86`) before merging, ledger the source); `test_three_point_equals_independent_four_point_on_random_sequences` (a four-point implementation written in the test file, McInnes and Meehan 2008; 200 random sequences of lengths 3 to 40; equal multisets of `(span, mean, count)`); `test_conservation_of_variation` (sum over cycles of `span * count` equals the total variation of the reversal sequence for the half-cycle policy); `test_rearranged_closed_history_has_only_full_cycles`; `test_cycle_multiset_is_invariant_to_the_starting_reversal_for_rearrange`; `test_brute_force_range_pairs_on_short_sequences`; `test_indices_refer_to_original_samples` (`x[i_start]` and `x[i_end]` reproduce `mean +- span/2`); `test_constant_series_gives_no_cycles`; `test_two_samples_give_one_half_cycle`; `test_plateaus_and_equal_neighbours_are_collapsed`; `test_gate_removes_small_reversals`; `test_ca_block_valley_to_valley_gives_m_full_cycles` (`m` repetitions of valley-peak ending in a valley: `m` cycles of count 1 with `"half"`).
- [ ] **Step 2: Run** `uv run pytest tests/test_hist_rainflow.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** with a stack; no equation beyond the standard algorithm.
- [ ] **Step 4: Run** the same command. Expected: PASS.
- [ ] **Step 5: Commit** `feat(history): add ASTM E1049 rainflow counting`.

---

### Task 4: Events

**Files:** Create `src/multifat/history/events.py`, `tests/test_hist_events.py`.

**Interfaces:**
- Consumes: `rainflow.Cycle`, `signal.TimeHistory`, `multifat.fe.ranges.LongestChord`, `RangeResult`.
- Produces: `Event(pair: tuple[int,int], count: float, S_range: ndarray, E_range: ndarray | None, span: float | None)` frozen dataclass; `events_from_cycles(history, cycles) -> list[Event]` (`pair = (i_start, i_end)`, `S_range = S[i_start] - S[i_end]`); `events_from_longest_chord(history, *, metric="deviatoric_norm", on="stress") -> list[Event]` (one event, `count = 1.0`); `merge_identical(history, events) -> list[Event]` merging events whose two end tensors (stress and strain) are equal, summing counts and keeping the first pair; `check_event_count(events, max_events: int = 2000)` raising `ValueError` naming the limit.

- [ ] **Step 1: Write failing tests:** `test_events_carry_full_tensors_at_the_reversals`; `test_ranges_match_difference_of_end_tensors`; `test_longest_chord_event_matches_fe_longest_chord`; `test_merge_identical_sums_counts_and_keeps_distinct_means` (same span, different mean stay separate); `test_merge_is_idempotent_and_conserves_total_count`; `test_event_count_guard_names_the_limit`.
- [ ] **Step 2-5:** run (FAIL), implement, run (PASS, plus notation), commit `feat(history): add events from counted cycles and longest chord`.

---

### Task 5: Accumulation

**Files:** Create `src/multifat/history/accumulate.py`, `tests/test_hist_accumulate.py`.

**Interfaces:**
- Consumes: `multifat.eval.fat_cycles.fatigue_cycles`.
- Produces: `SNCurve(A, b, method, variant="", design_loading="tension", source="")` (validation: finite, `A > 0`, `b < 0`, else `ValueError` mentioning "no usable life"); classmethod `from_fit(fit: Mapping, *, method, variant="", design_loading="tension", source="")` reading `fit["A_surv"]` and `fit["b"]`; `cycles_to_failure(dp: ndarray) -> ndarray` (`inf` where `dp <= 0`), `dp_at(N)`; `Endurance(mode: Literal["extrapolate","ignore_below"]="extrapolate", dp_limit: float | None=None)` with classmethods `ignore_below(dp_limit)` and `from_knee(curve, N_knee)` (`dp_limit = A N_knee^b`); `miner(dp, counts, curve, endurance=Endurance()) -> MinerResult(damage: ndarray, n_f: ndarray, total: float)`; damage `count / N_f`, zero where `dp <= 0` or (`ignore_below` and `dp < dp_limit`).

- [ ] **Step 1: Write failing tests:** `test_basquin_roundtrip_dp_and_cycles` (`cycles_to_failure(dp_at(N)) == N` to 1e-9); `test_curve_rejects_flat_positive_or_nonfinite_slope`; `test_from_fit_uses_a_surv`; `test_miner_equals_hand_sum_on_analytic_curve` (`A = 1000, b = -0.1`, dp values and counts chosen by hand); `test_zero_and_negative_dp_give_zero_damage_without_warnings`; `test_ignore_below_differs_from_extrapolate_exactly_by_the_events_below_the_limit`; `test_from_knee_limit_matches_curve`; `test_total_damage_zero_is_finite_and_reported_as_zero`.
- [ ] **Step 2-5:** run (FAIL), implement, run (PASS), commit `feat(history): add Basquin curve and Miner accumulation`.

---

### Task 6: Eligibility

**Files:** Create `src/multifat/history/eligibility.py`, `tests/test_hist_eligibility.py`.

**Interfaces:**
- Consumes: `multifat.models.registry` (`get_models`, `stress_only_models`, `_available_cp_methods`, `_available_direct_models`), `multifat.models.references.ensure_enabled`, `multifat.fe.bridge.check_eligibility`, `accumulate.SNCurve`, `proportionality.Proportionality`.
- Produces: `IneligibleHistoryError(ValueError)`; `is_native_life(model) -> bool`; `check_accumulation_eligibility(model, curve: SNCurve, *, variant: str, has_strain: bool)` raising for: curve method or variant different from the model's name and the requested variant, native-life model, strain-needing model without strain (`stress_only_models()`), model declaring `requires_nonlocal`; `check_history_kind(p: Proportionality, *, allow_approximation: bool) -> str | None` raising unless proportional or allowed, returning the approximation note when allowed; `method_matrix() -> list[dict]` one row per enabled registered method with keys `method`, `strain_needed`, `native_life`, `requires_nonlocal`, `accumulation`, `reason`; `python -m multifat.history.eligibility` prints it as a markdown table.

- [ ] **Step 1: Write failing tests:** `test_matrix_has_a_row_for_every_enabled_registered_method`; `test_native_life_methods_refused_for_accumulation` (every row with `native_life` true raises); `test_strain_method_refused_without_strain_and_accepted_with_it`; `test_curve_method_and_variant_mismatch_rejected`; `test_stub_requires_nonlocal_refused`; `test_nonproportional_refused_unless_allowed_and_note_returned`; `test_markdown_table_lists_every_method`.
- [ ] **Step 2-5:** run (FAIL), implement, run (PASS, plus notation), commit `feat(history): add eligibility checks and the method matrix`.

---

### Task 7: Assessment

**Files:** Create `src/multifat/history/assess.py`, `tests/test_hist_assess.py`.

**Interfaces:**
- Consumes: all earlier tasks; `multifat.fe.bridge.to_case_inputs`, `fe.ranges.ExplicitPair`; `multifat.models.errors.LoadingScopeError`; `multifat.planes.sampling.build_R_list`.
- Produces: `assess_history(history: TimeHistory, model, params: Mapping[str, float], curve: SNCurve, *, R_list, variant: str = "", counting: Literal["rainflow","longest_chord"] = "rainflow", channel="dominant", gate: float = 0.0, residue="half", endurance: Endurance = Endurance(), allow_nonproportional_approximation: bool = False, max_events: int = 2000) -> HistoryAssessment`; `HistoryAssessment` frozen dataclass with `method`, `variant`, `curve`, `settings` (counting, channel name, gate, residue, endurance, plane-grid description `{"n_planes": len(R_list)}`, proportionality ratio, `approximation` note or `None`), `history_sha256`, `history_provenance`, `events` (list of dicts: `pair`, `count`, `S_range`, `dp`, `n_f`, `damage`), `damage_per_block`, `blocks_to_failure` (`inf` when zero damage), `cycles_per_block`, `cycles_to_failure`, `to_dict()`, `to_json(path)`. Flow: check eligibility, resolve channel, count, build and merge events, guard the count, evaluate each unique event with `evaluate_case` through `to_case_inputs(history.to_block(), ExplicitPair(a, b), meta=history.meta(), method=model)` reading `values[variant]` (a missing variant raises `IneligibleHistoryError`; `LoadingScopeError` is re-raised as `IneligibleHistoryError` naming method and event), then `miner`.

- [ ] **Step 1: Write failing tests:** `test_constant_amplitude_history_reproduces_the_ca_value_and_miner_sum` (history of `m = 5` repetitions of one valley-peak-valley cycle of a tension-torsion proportional load, `FIN` or `FS` with nominal constants and `R_list = build_R_list(n1=3, n2=6)`; the single merged event's `dp` equals `model.evaluate_case(...)` called directly on the two end tensors, `cycles_per_block == 5`, `damage_per_block == 5 / curve.cycles_to_failure(dp)` to 1e-12, `blocks_to_failure == 1 / damage_per_block`); `test_variable_amplitude_history_sums_event_damages` (a three-level block; total equals the hand sum of per-event damage); `test_residue_half_and_rearrange_both_recorded_and_differ_as_documented`; `test_longest_chord_counting_gives_one_event`; `test_nonproportional_history_refused_then_allowed_with_recorded_note`; `test_native_life_method_refused`; `test_strain_method_refused_without_strain`; `test_loading_scope_error_reported_with_method_and_event` (GOLOS on a non-fully-reversed history); `test_missing_variant_rejected` (`variant="_ext"` on a model without it); `test_curve_for_another_method_rejected`; `test_zero_damage_history_has_infinite_life_without_warnings`; `test_result_json_is_reproducible_and_carries_the_history_hash` (two runs give equal `to_dict()`; the hash equals `history.sha256()`); `test_event_count_guard`.
- [ ] **Step 2: Run** `uv run pytest tests/test_hist_assess.py -q`. Expected: FAIL.
- [ ] **Step 3: Implement** `assess.py`.
- [ ] **Step 4: Run** the same command, the earlier `tests/test_hist_*.py`, `tests/test_notation.py`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(history): add the history assessment and result record`.

---

### Task 8: Isolation, example, documentation, gate

**Files:** Create `tests/test_hist_isolation.py`, `src/multifat/history/example.py`, `docs/load-histories.md`. Modify `docs/decision-log.md` (`d079` onward), `docs/notation.md` (if rows missing), `AGENTS.md` (layout line for `history/`), `README.md` (software-map row), `docs/audit-checkpoint.md`.

- [ ] **Step 1: Write failing tests:** `test_no_existing_module_imports_history` (AST walk like `tests/test_fe_isolation.py`); `test_history_imports_no_cli_or_report` (`history/` imports nothing from `multifat.cli` or `multifat.report`); `test_fe_does_not_import_history`; `test_documented_example_runs` (`example.main(["--out", tmp])` returns 0, writes the recorded history and the result JSON; the result's `blocks_to_failure` is finite and positive); `test_doc_lists_every_matrix_method` (every `method_matrix()` name appears in backticks in `docs/load-histories.md`).
- [ ] **Step 2: Implement** `example.py` (a documented generator function building a proportional tension-torsion block, an analytic curve, one method, prints the result); write `docs/load-histories.md` by hand: purpose, the contract (inputs, units, constructors), the proportional-domain rule and the approximation flag, counting and residue policies, the curve and accumulation, the eligibility matrix (from `python -m multifat.history.eligibility`), reproducibility record, limits, benchmarks B1-B4 and follow-ups F1-F8 with their Zotero keys (copy from the spec), the example command. Add decision-log entries `d079`-`d085` (scope, proportional domain, caller-supplied curve, channel, residue and block semantics, native-life refusal, package name) in the existing index and record format; AGENTS.md and README rows; checkpoint section.
- [ ] **Step 3: Run the gate:** `uv run pytest -q`, `uv run black --check . && uv run ruff check .`, `uv run python tools/check_docs.py`, `uv run python tools/check_scientific_reviews.py`, `uv run python tools/generate_method_docs.py --check`, `uv run python tools/generate_workbook_docs.py --check`. Expected: all pass. Verify `git diff origin/main --stat -- src/multifat ':!src/multifat/fe' ':!src/multifat/history'` shows nothing new beyond the already merged `fe` work (no existing file changed by this plan).
- [ ] **Step 4: Regression snapshot** (no existing source changed, so this confirms): run `tools/regression_snapshot.py` on the workstation as before (bundle, `scp -O`, Task Scheduler) comparing `origin/main` with the branch, `--skip-calibrate`; record the result in the checkpoint.
- [ ] **Step 5: Commit** `docs(history): publish the load-history specification, decisions and example`.

---

## Self-review

- **Spec coverage:** goals 1-8 map to Tasks 1 (contract, generator, record), 3-4 (counting, events, longest chord), 7 (damage per event, result), 5 (accumulation), 6 (eligibility), 8 (docs, benchmarks B1-B4 and follow-ups F1-F8, isolation, preservation); acceptance criteria to Tasks 1, 6, 8 (benchmarks and tasks), 7 and 8 (CA preservation, snapshot).
- **Types:** `TimeHistory`, `Cycle`, `Event`, `SNCurve`, `Endurance`, `MinerResult`, `Proportionality`, `HistoryAssessment`, `IneligibleHistoryError` are defined once and used with the same names.
- **Known limits:** the ASTM E1049 example counts come from memory of the standard's worked example and are verified against a source in Task 3 before merge; multiaxial counting and equivalent-cycle extractors are the follow-ups, not part of this plan.
