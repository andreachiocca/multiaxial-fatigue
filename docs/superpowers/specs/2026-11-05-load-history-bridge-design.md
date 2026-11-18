# General load history to damage parameters and life: design

Issue: [andreachiocca/multiaxial-fatigue#20](https://github.com/andreachiocca/multiaxial-fatigue/issues/20),
"Define time-history inputs and a variable-amplitude assessment pathway" (assignee: Tommaso Grossi).
Date: 2026-11-05. Status: draft for approval, nothing built. Builds on the FE local-history interface
([design](2026-11-05-fe-local-history-interface-design.md), `src/multifat/fe/`).

## Intent and why

The assessment path represents constant-amplitude, single-harmonic loading. A user who has a
**general** multiaxial history (synchronised tensor components over time) cannot feed it to the
repository. This work is the **bridge** from such a history to the damage parameters of the
existing criteria and on to a fatigue life:

    general history -> reduction (counting / equivalent cycle) -> events -> damage parameter per event
                    -> life per event from the method's calibrated curve -> accumulated damage and life

The formalism of the exact critical-plane paper (Grossi, Chiocca and Sgamma, Sec. 2.2) fixes the
vocabulary: a history is N load cases; a **range strategy** maps them to ranges; peak quantities are
maxima over load cases. Counting algorithms, equivalent-cycle extractors and accumulation rules are
the missing layers around it.

## Decided with the user (2026-11-05)

- The bridge covers the whole chain above, **including accumulation into a life**.
- Reduction algorithms are selectable and recorded, never silently substituted for each other
  (rainflow, Wang-Brown, PDMR, longest chord, minimum circumscribed circle, and so on).
- A literature corpus was built (below). **This round implements only the basic parts**; the
  multiaxial counting and equivalent-cycle rules are specified as follow-up tasks and implemented
  after the user (or a cheaper model) has fetched and the sources have been read.
- Work continues on the branch `feat/fe-and-time-history-inputs`.

## Scientific frame for the basic implementation (my call, open to change)

For a **proportional** history, `sigma(t) = lambda(t) A + sigma_static`, a scalar channel carries all
the information: rainflow counting of `lambda(t)` is exact, and each counted cycle is a straight
two-state path between two samples of the history, which the existing case-level models and the
two-state bridge (`fe.bridge`) evaluate exactly (the bridge is exact for two-state paths). So the
basic implementation is rigorous in a clearly stated domain: **proportional (single-channel)
histories of any length and sequence**. For a non-proportional history the same steps are only an
approximation: it is detected, refused by default with a diagnostic, and allowed only with an explicit
`allow_nonproportional_approximation=True`, which is written into the result.

## Goals

1. A **time-history input contract**: one location, `N` synchronised stress tensors (and optional
   strain), a time or cycle axis, units, provenance. Built from arrays, from a CSV file, from an
   HDF5 history file (`fe`), or **generated from a documented Python function**. The sampled history
   actually used is recorded (arrays, SHA-256, generator name and parameters) so a result is
   reproducible.
2. **Reduction**, basic set: turning-point extraction with an optional gate; **scalar-channel rainflow**
   (ASTM E1049 three-point form with half-cycle residue; four-point equivalence tested); selectable
   scalar channels; `LongestChord` (exists) as the single-equivalent-cycle option. Events keep the
   sample indices of their two reversals so the full tensors are available.
3. **Damage parameter per event** through the existing models and the two-state bridge (no model
   change), with explicit eligibility diagnostics (below).
4. **Accumulation into a life**: damage per event from the method's calibrated Basquin master curve
   `DP = A N^b`, linear Miner sum, and life in blocks and in cycles (below).
   - **Where the curve comes from.** The repository does not persist the damage-parameter curve:
     `A`, `b` and `A_surv` are fitted in memory by the pipeline (`Evaluation.fit` in `cli/run.py`,
     `fit_power_law_with_survival_std`) and no record or CSV stores them, and the spec forbids changing
     any record or CSV. So `SNCurve(A, b, method, variant, source)` is **supplied by the caller**, with a
     constructor from the in-memory curve mapping the pipeline produces (`cal["fit"]` or
     `cal["fit_ext"]`, keys `A`, `b`, `A_surv`). It uses `A_surv`, as the existing life prediction does
     (`run.py`). No dependency on `cli`; a loader that re-runs a workbook calibration is a follow-up.
   - **Variant and plane set.** `assess_history` takes the damage-parameter variant (`""` or `"_ext"`),
     the plane set (`R_list`, or the grid sizes it is built from), and refuses a curve whose method or
     variant differs. It records the plane-grid description; that the curve was calibrated on the same
     grid is the caller's responsibility, stated in the docs.
   - **Block semantics.** The history is **one block that repeats until failure**: damage per block
     `D = sum n_i / N_f,i`, life `1 / D` blocks, life in cycles `(sum of counts per block) / D`.
   - **Residue.** `residue="rearrange"` (default after review: the block repeats, so every cycle must
     close; the block is rotated to start at its global extremum so that no residue remains, exact for
     strictly periodic repetition) or `"half"` (ASTM E1049: residue as half cycles; it under-counts a
     repeating block and the result says so). The policy is recorded in the result.
   - **Endurance.** `extrapolate` (the curve is used at every damage-parameter level) or
     `ignore_below(dp_limit)`: an event whose damage parameter is below `dp_limit` contributes zero,
     per event, whatever its count. A helper converts a knee life to a limit, `dp_limit = A N_knee^b`.
5. A **result object** that records history, counting algorithm and settings, curve, accumulation
   rule and every event, serialisable to JSON.
6. **Eligibility and diagnostics**: per registered method and per history kind; `IneligibleHistoryError`
   naming the reason: non-proportional history, missing strain for a strain-needing method, method
   declaring `requires_nonlocal`, curve missing or mismatched, a curve with `b >= 0` (the pipeline
   treats it as "no usable life"), and **methods that define their own life law** (`native_life`,
   for example MDC, MMCCM, Morel): accumulation through a single Basquin curve of the damage parameter
   is not their formulation, so they are refused for accumulation (they remain usable for damage
   parameters per event). A method is a native-life method when `getattr(model, "native_life", None)`
   is not `None` (the test `cli/run.py` uses). After the merge with upstream `main` (5 November) two
   more cases apply: a criterion whose loading scope excludes an event (`LoadingScopeError`, for
   example GOLOS, defined for fully reversed in-phase plane stress only) is reported as ineligible for
   that history, naming the method and the event, never skipped silently; and a curve records the
   model's `design_loading` (`tension` or `torsion`) as provenance. Life from a damage parameter is
   `fatigue_cycles(DP, A_surv, b)` of `eval/fat_cycles.py`, the function the pipeline uses.
7. **Documentation**: `docs/load-histories.md` with the contract, the eligibility matrix, the defined
   benchmarks and the follow-up task list; decision-log entries; notation rows.
8. Tests that preserve the constant-amplitude results (a history made of repeated constant-amplitude
   blocks reproduces the constant-amplitude damage parameter and life) and leave every existing
   number unchanged (regression snapshot).

## Non-goals and anti-goals

- No multiaxial counting (Wang-Brown, enclosing-surface, minimum-circumscribed-circle rainflow), no
  PDMR, MCC, MCE, prismatic hull, maximum rectangular hull or moment-of-inertia extractors, no
  maximum-normal-stress-over-cycle extension of the criteria (exact-cp Eqs. 8-15), no nonlinear or
  Haibach accumulation, no sequence-effect models: all are **specified as follow-ups**, because they
  need the papers to be read. No equation is implemented from memory.
- No application of the harmonic criteria to a non-proportional history without the explicit,
  recorded approximation flag. No silent substitution of one algorithm for another.
- No change to any model, calibration record, fingerprint or existing CSV; the new package is imported
  by nothing existing; no GPL code (repository is MIT, `d071`).
- No claim of validation: software regression tests and analytical oracles are not experimental
  validation of a fatigue method; experimental benchmarks are defined, not run.
- Not a notch, size, surface or mean-stress-correction layer: the damage parameters are whatever the
  existing criteria compute, and the curve is whatever the calibration stored.

## Constraints

- Same repository conventions as the FE work: notation rows before new public names (`E` is Young's
  modulus, strain tensors `E_`), numpydoc, module docstring on line 1, black 100, ruff correctness,
  MPa and mm, tensor shear, decisions in `docs/decision-log.md` (after `d078`), checkpoint updated,
  regression snapshot unchanged, commits local, no push, no attribution lines.
- `h5py` stays optional (extra `fe`); CSV and array inputs need only numpy and pandas.
- The two-state bridge convention (`S0 = S_hist[a]`, `S1 = S_hist[b]`, `S_sin = (S0 - S1)/2`, strain
  NaN when absent) is reused, not redefined. `fe` must not import the new package.

## Design (names are the builder's within the notation rule; behaviour is fixed)

New package `src/multifat/history/`:

| Module | Content |
|---|---|
| `signal.py` | `TimeHistory`: `t` (N,), `axis_kind` (`time` or `cycle`, strictly increasing), `S_hist` (N,3,3) MPa, optional `E_hist`, units and provenance; constructors from arrays, `from_csv`, `from_file` (HDF5), `from_function(fn, t, **params)` (records the function's qualified name and parameters); `to_block()` giving the `fe` `HistoryBlock` of one location; `record(path)` writing the sampled arrays and metadata |
| `channels.py` | scalar channels of a tensor history, each returning `(N,)`. The **default** is the projection of the mean-removed six-component stress path on its dominant singular direction (sign fixed by the largest component): for a proportional history it is exactly linear in the load factor `lambda(t)`, so rainflow on it counts the true load reversals even with a static stress (the signed von Mises channel is not linear in `lambda` when `sigma_static != 0`). Also offered, explicitly: a component, the normal stress on a given direction, the signed von Mises stress (sign from the hydrostatic part, stated), a user callable |
| `rainflow.py` | `turning_points(x, gate)`; `rainflow(x, residue="half")` returning cycles with `range`, `mean`, `count` (1.0 or 0.5) and the two sample indices; residue handling as in ASTM E1049, with the alternative of closing the residue by repeating the block recorded in the result |
| `events.py` | `Event(pair, count, S_range, E_range)`; events from rainflow cycles, from `LongestChord`, or from an explicit pair; caching of identical pairs; a guard on the number of events |
| `proportionality.py` | a measure of non-proportionality of the history (rank of the centred deviatoric stress path by singular values; threshold stated) and the verdict used by the eligibility check |
| `accumulate.py` | `SNCurve(A, b, source)` and `cycles_to_failure(dp)`; `miner(events, curve, endurance=...)` returning damage per event and total; loading a curve from a calibration record for a method and variant |
| `assess.py` | `assess_history(history, model, params, curve, *, R_list, variant="", counting="rainflow", channel="dominant", gate=0.0, residue="half", endurance=..., allow_nonproportional_approximation=False)` and the result object |
| `eligibility.py` | `IneligibleHistoryError`; the eligibility check and the matrix of methods |

The damage parameter of an event is `model.evaluate_case(S0, E0, S1, E1, params, R_list, harmonics)`
through `fe.bridge` with `ExplicitPair`, reading the `""` (or `"_ext"`) value; the plane set is the
caller's, as in the existing pipeline.

## Issue acceptance criteria (verbatim from #20) and where each is met

- *Document and validate the time-history input contract* -> goal 1, verification 1, `docs/load-histories.md`.
- *Identify supported methods and give explicit diagnostics for unsupported ones* -> goal 6,
  verification 4, the eligibility matrix.
- *Define source-supported benchmarks and separate implementation tasks* -> benchmarks B1-B4 and
  follow-ups F1-F8 below.
- *Preserve current constant-amplitude results and calibration compatibility* -> goal 8,
  verification 3 and 7 (no record or curve format changes; the curve is an input).
- Issue note: *do not apply existing harmonic criteria to arbitrary histories without a justified
  formulation* -> the proportional-domain rule and the recorded approximation flag.

## Verification

1. **Contract documented and validated**: tests reject non-monotonic time, non-finite or asymmetric
   tensors, mismatched shapes, missing units, a generator that returns the wrong shape; generator and
   CSV round trips equal the arrays; `record()` reproduces the history exactly (SHA-256).
2. **Counting against independent oracles**: the ASTM E1049 example sequence with its published
   counts (to be taken from the standard or a secondary source in the corpus and cited in the test),
   equality of three-point and four-point results on random sequences (McInnes and Meehan 2008), the
   conservation property of rainflow (every reversal-to-reversal variation is accounted once),
   sequence-order independence of the cycle multiset for a closed history, a brute-force range-pair
   check on short sequences.
3. **Constant-amplitude preservation**: a history that is `m` repetitions of one constant-amplitude
   cycle, starting and ending at the same valley (so the residue is empty), gives `m` identical full
   events whose damage parameter equals the constant-amplitude value of the existing path, and
   accumulated damage `m / N_f`. With `residue="rearrange"` they are `m` full events; with `"half"` the
   standard half-cycle rule at the start makes them half events whose counts sum to `m` (observed
   during the build; the test asserts the sum).
4. **Eligibility**: a proportional history is accepted for every registered method (strain-needing
   ones only with strain); a non-proportional history raises `IneligibleHistoryError` unless the
   approximation flag is set, and the flag is recorded; a stub `requires_nonlocal` method is refused.
5. **Accumulation**: Miner sum on an analytic curve equals the hand sum; `ignore_below` and
   `extrapolate` differ exactly by the events below the limit; blocks-to-failure is `1 / D`.
6. **Reproducibility**: the same history and settings give the same result JSON; the JSON contains the
   sampled-history hash.
7. **Isolation**: no existing module imports `history`; regression snapshot unchanged; full suite,
   black, ruff, `check_docs` pass.

## Benchmarks defined for later (not run here)

- **B1** ASTM E1049 example counts; **B2** three-point/four-point equivalence (McInnes and Meehan);
  **B3** constant-amplitude reproduction (this round).
- **B4** experimental variable-amplitude data sets that must be read first: Gates and Fatemi
  (un-notched and notched, Parts I and II, 2017-2018), Susmel and Tovo (2011), Wang and Susmel
  (2015, 2016), Morettini et al. (random loading on ASTM A105, 2019 and 2021), Sonsino and Kueppers
  (welded joints, 2001), Kim and Park (1999). Each needs its load histories or enough to regenerate them,
  the calibration constants and the lives.

## Follow-up tasks (separate issues, each blocked on reading the cited sources)

Sources are in the Zotero collection "multifat #20 load histories" (key `4J6FG4BS`, under
"Multiaxial Fatigue"); PDFs are to be fetched later.

- **F1** Enclosing-surface methods and the Modified Wang-Brown multiaxial rainflow: Meggiolaro and de
  Castro 2012 Parts I and II; Wang and Brown 1996 Parts 1 and 2.
- **F2** Rainflow with minimum-circumscribed-circle shear ranges: Bibbo et al. 2022; Anes et al. 2014.
- **F3** Equivalent-range extractors (PDMR, MCC, MCE, prismatic hull, maximum rectangular hull, moment
  of inertia): Papadopoulos 1998; Carpinteri et al. 2014; Araujo et al. 2014; Meggiolaro and de Castro 2014.
- **F4** Maxima over the cycle's load cases and the criteria that need them (exact-cp Eqs. 8-15).
- **F5** Accumulation variants: Haibach modified rule (source still to be located), Lv et al. 2014,
  Benkabouche et al. 2015, Okamura et al. 1979.
- **F6** Critical-plane counting and damage under variable amplitude: Susmel and Tovo 2011, Carpinteri,
  Spagnoli and Vantadori 2017 (review), Langlais 2003, Kim and Park 1999.
- **F7** Residue handling and sequence effects: Marsh et al. 2016, Anthes 1997, Glinka and Kam 1987,
  Amzallag et al. 1994, Rychlik 1987, Socie 1992, Lee and Tjhung 2012.
- **F8** Benchmark reproductions B4 and comparison with the reviews of Chen et al. 2012 and
  Bibbo et al. 2022.

## Assistant's calls (open to change)

The proportional-domain scope and the approximation flag; reuse of `fe` types and the two-state bridge;
the curve being caller-supplied (with a constructor from the pipeline's in-memory curve) using
`A_surv`; refusing native-life methods for accumulation; the dominant-direction default channel;
block-repeats-until-failure semantics with `half` and `rearrange` residue policies; per-event
`dp_limit` endurance; variant and plane-grid handling; package name `history`; one location per
history in this round.

## Left to the builder

Module and function names within the notation rule; the numerical threshold of the proportionality
measure and its justification; the default gate; CSV column naming; the JSON schema of the result;
the exact ASTM E1049 residue treatment (taken from the standard or a corpus paper and cited);
the event-count guard value; where tests keep their random seeds.
