# Load histories: from a general history to damage parameters and life

A bridge from a **general** multiaxial history (synchronised stress tensors over time) to the damage
parameters of the existing criteria and on to an accumulated fatigue life. This page is the
specification of the contract and of the basic implementation (`multifat.history`), the eligibility
matrix, the benchmarks that are defined and the follow-up tasks. It builds on the FE local-history
interface ([fe-input-format.md](fe-input-format.md)).

Software regression and scientific validation are kept apart. The tests show that the counting agrees
with independent algorithms and with a worked example of the standard as reproduced in secondary
sources, that constant-amplitude results are reproduced and that no existing number moved. None of this validates a fatigue criterion against
experiments, and the example uses an illustrative curve, not a calibration.

## The chain

    history -> scalar channel -> rainflow (or longest chord) -> events -> damage parameter per event
            -> life per event from the method's Basquin curve -> Miner sum -> blocks to failure

The vocabulary is that of the exact critical-plane paper (Grossi, Chiocca and Sgamma, Sec. 2.2): a
history is `N` load cases; a counting or extraction rule turns them into **events**, each a pair of
samples `(a, b)` with a count; the criteria see the full tensors at the two samples.

## What is exact, and what is not

For a **proportional** history, `S(t) = lambda(t) A + S_static`, a scalar channel carries all the
information. Rainflow on a channel that is affine in `lambda(t)` counts the true load reversals, and
each counted cycle is a straight two-state path between two samples, which the existing criteria and
the two-state bridge evaluate exactly. So the **reduction** (counting and event construction) is exact for
**proportional histories of any length and sequence**. The damage parameter is the criterion's own, and
the Basquin curve and Miner's sum remain modelling choices; nothing here validates them.

A non-proportional history is detected (`proportionality`: the second over the first singular value of
the mean-removed shear-weighted six-component stress path, default tolerance `1e-6`, a numerical
floor that `proportionality_tol` changes, for example for noisy measured data) and **refused** with
`IneligibleHistoryError`. For criteria that read strain the **strain path** is checked too, and
stress and strain must share one affine scalar load function: separate straight paths are not
enough if their reversals differ. The joint rank test normalizes the two paths independently, so
their units cannot hide one behind the other. Constant stress with varying strain is refused
because the stress channel cannot count its cycles. `allow_nonproportional_approximation=True` accepts it as an
approximation, and the flag and a note are written into the result. Multiaxial counting for
non-proportional paths is a follow-up (below).

## The time-history contract

`TimeHistory` holds one location: `N >= 2` stress tensors (MPa, tensor components), optional total
strain tensors (dimensionless, tensor shear) and a strictly increasing time or cycle axis, with units
and provenance. Units are always stated, never assumed.

| Constructor | Input |
|---|---|
| `TimeHistory.from_arrays(time, S_hist, E_hist, stress_unit=...)` | arrays; `stress_unit` required (Pa, kPa, MPa, GPa) |
| `TimeHistory.from_csv(path, stress_unit=..., time_unit=...)` | columns `t, sxx, syy, szz, sxy, syz, sxz` and optionally `exx ... exz`; `columns=` renames them; engineering shear via `strain_shear="engineering"` |
| `TimeHistory.from_file(path, location_id=None)` | a time-sample `multifat.fe` file with one uniquely identified location (needs `h5py`, extra `fe`); repeated IDs across instances/materials and unit cases are refused; use `node_histories` with complete location selection or load functions |
| `TimeHistory.from_function(fn, time, **params)` | a documented Python function returning `(N, 3, 3)` stress in MPa or a dict with `"S"` and `"E"`; the function name and `params` (JSON-serialisable) are recorded |

Validation rejects non-increasing or duplicate times, non-finite values, asymmetric tensors, shape
mismatches, a single sample and unknown units. `history.record(path)` writes the sampled arrays
(`.npz`) and the metadata with the SHA-256 (`.json`); `TimeHistory.from_record` reads them back and
verifies the hash, so a result can be reproduced from the exact history that was assessed.

## Counting and events

- **Channel.** The default is the projection of the mean-removed, shear-weighted six-component stress
  path on its dominant singular direction (sign fixed by its largest component). For a proportional
  history it is exactly affine in `lambda(t)`, even with a static stress. Also available: a
  component, the normal stress on a direction, the von Mises stress signed by the hydrostatic part
  (not affine in `lambda` when `S_static != 0`, so never the default) and a callable.
- **Rainflow.** ASTM E1049 three-point form (Downing and Socie 1982): a range that contains the
  starting point is a half cycle, other closed ranges are full cycles, the residue is counted as half
  cycles. `residue="rearrange"` (the default of `assess_history`, because a block repeats until failure)
  rotates the reversals to start at the highest peak so that no residue remains and every cycle is
  full (exact for strictly periodic repetition; the end of the block is taken to connect to its
  start). `"half"` is the standard rule for a single history; for a repeating block it counts the
  cycle that spans consecutive blocks as halves, which under-counts (about a factor 2 on the shipped
  example), and the result says so in `notes`. `gate` is a hysteresis
  filter. Cycles keep the sample indices of their two reversals.
- **Longest chord.** `counting="longest_chord"` gives one equivalent event, the pair of load cases
  farthest apart (`multifat.fe.ranges.LongestChord`).
- Events with identical end tensors (in either order) are merged and their counts summed. A guard
  (`max_events`, default 2000) stops runaway cases.

## Damage parameter and accumulation

The damage parameter of an event is the criterion's value on the two-state path (`""` or `"_ext"`
variant), through `multifat.fe.bridge`; no model changes. The plane set `R_list` is the caller's.

The life uses `DP = A N^b`, an **`SNCurve` supplied by the caller**: the repository does not persist
the damage-parameter curve (the pipeline fits `A`, `b`, `A_surv` in memory), and nothing here changes
a calibration record. `SNCurve.from_fit(fit, method=..., variant=...)` reads `A_surv` and `b` of the
pipeline's in-memory fit, as the existing life prediction does; the curve must come from the same
method, variant and plane set as the assessment. A curve with `b >= 0` is refused.

- The history is **one block that repeats until failure**: `D = sum n_i / N_f,i` per block, life
  `1 / D` blocks, life in cycles `(sum of counts per block) / D`. Zero damage gives infinite lives.
- Events with a non-positive damage parameter contribute nothing. `Endurance.ignore_below(dp_limit)`
  ignores events below a limit (`Endurance.from_knee(curve, N_knee)` converts a knee life); the default
  extrapolates the curve at every level.
- With `residue="half"` a constant-amplitude block of `m` cycles is counted as halves that sum to `m`
  (the half-cycle rule at the start); with `"rearrange"` it is `m` full cycles. A life below one
  cycle (a damage parameter above the curve's `A`) is noted in the result.

The result (`HistoryAssessment`, `to_dict()` / `to_json()`) records the settings, the curve, the
history hash and provenance, every event with its damage parameter, life and damage, and the lives.
Infinite values are written as the string `"inf"`.

## Eligibility

`IneligibleHistoryError` names the reason. Refused: a curve of another method or variant; a method
with a **native life law** (`native_life`), because accumulation through a single Basquin curve of its
damage parameter is not its formulation; a strain-needing method on a history without strain; a method
declaring `requires_nonlocal`; an event outside a criterion's loading scope (`LoadingScopeError`, for
example GOLOS, defined for fully reversed in-phase plane stress only), reported with the method and
the event; a missing variant; a non-proportional history (unless allowed). Native-life methods stay
usable for damage parameters per event through `event_damage_parameters(history, model, params,
R_list=...)`, which needs no curve.

A method marked `scope restricted` rejects a generic load with a mean stress through
`LoadingScopeError` (found by evaluation, like the stress-only list); expect `IneligibleHistoryError` on
general events for it.

Matrix of the enabled registered methods (`python -m multifat.history.eligibility`):

| method | strain needed | native life | scope restricted | accumulation | reason |
|---|---|---|---|---|---|
| `FS` | yes | no | no | yes |  |
| `FIN` | no | no | no | yes |  |
| `SWT` | yes | no | no | yes |  |
| `MATAKE` | no | no | no | yes |  |
| `MCD` | no | no | no | yes |  |
| `MDC` | no | yes | no | refused | native life law |
| `MWCM` | no | no | no | yes |  |
| `MMCCM` | yes | yes | no | refused | native life law |
| `VSE` | yes | no | no | yes |  |
| `MGSE_YU` | yes | no | no | yes |  |
| `MGSE_ZHU` | yes | no | no | yes |  |
| `MSWT` | yes | no | no | yes |  |
| `MKBM` | yes | no | no | yes |  |
| `LIU2021` | yes | no | no | yes |  |
| `LI` | yes | no | no | yes |  |
| `GSE` | yes | no | no | yes |  |
| `GSA` | yes | no | no | yes |  |
| `NINIC` | no | no | no | yes |  |
| `MWHS` | yes | no | no | yes |  |
| `LIU2025` | no | no | yes | yes |  |
| `ZHANGYAO` | no | no | no | yes |  |
| `KUJAWSKI2021` | yes | no | no | yes |  |
| `BP` | no | no | no | yes |  |
| `CAIM` | no | no | no | yes |  |
| `SWTD` | yes | no | no | yes |  |
| `OTT` | no | no | no | yes |  |
| `TVEIT2025` | no | yes | no | refused | native life law |
| `LIUMAHADEVAN` | no | yes | no | refused | native life law |
| `DANGVAN` | no | no | no | yes |  |
| `CROSSLAND` | no | no | no | yes |  |
| `WANG2022` | no | no | yes | yes |  |
| `GAM` | no | no | no | yes |  |
| `GOLOS` | no | no | yes | yes |  |
| `SINES` | no | no | no | yes |  |
| `PAPADOPOULOS` | no | no | no | yes |  |
| `PCR` | no | no | no | yes |  |
| `PI` | no | no | no | yes |  |
| `RFB` | no | no | no | yes |  |
| `VHN` | no | no | no | yes |  |
| `PBP` | no | no | no | yes |  |
| `COIN` | no | no | no | yes |  |
| `MMMH` | no | no | no | yes |  |
| `MMP` | yes | no | no | yes |  |
| `SSF` | no | no | yes | yes |  |
| `MOREL` | no | yes | no | refused | native life law |
| `EK` | yes | no | no | yes |  |

`MMP` is listed as needing strain only because the registry's probe constants (used to find the
stress-only methods) do not satisfy its validity condition `Sigm1 < Sig0 < 2 Sigm1`, so its probe
evaluation fails; with valid constants it is a stress criterion. The probe belongs to the registry,
which this work does not change; the row is conservative (a stress-only history is refused for it).

## FE nodes: from an import to a table per node

`multifat.history.nodes` joins the FE import to this chain. `node_histories(path)` turns each
location of a `multifat.fe` history file (for example from the CSV importer,
[fe-csv-format.md](fe-csv-format.md)) into a `TimeHistory`: for a `time_sample` file the load steps
are the samples (axis: the file's time, else the solver step numbers as a cycle axis); for a
`unit_case` file (linear elastic only) the unit cases are superposed with a table of load functions,
`S(t) = sum_k f_k(t) S_k`, matched to the steps by name. `assess_nodes(path, materials, R_list=...)`
then runs `assess_history` on every selected node with the `MaterialSetup` (criterion, constants,
Basquin curve) of the node's material and returns one row per node. A refusal (non-proportional
history, native-life method, missing strain, a material without a setup) is recorded as the node's
`status` and the run continues; `on_refusal="raise"` stops at the first. The checks that make this
exact (proportional node histories) and the approximations are those above; nothing here validates a
criterion.

## Limits

- Proportional histories only for an exact result; non-proportional ones are an explicit approximation.
- One location per history; no notch, size, surface or mean-stress correction beyond what a criterion
  itself contains.
- The criteria receive two load steps, so the maximum normal stress over all the load cases of a cycle
  is not carried (follow-up F4).
- Linear damage accumulation, no sequence effects beyond what counting keeps.
- `max_events` bounds the number of distinct events; evaluating each costs one plane search.

## Worked example

```bash
python -m multifat.history.example --out out
```

builds a proportional tension-torsion block from the documented function `tension_torsion_block`,
records it (`out/history.npz`, `out/history.json`), counts it, assesses it with the Findley criterion
on an illustrative curve and writes `out/result.json`. The curve and the constant are placeholders,
not a calibration.

## Benchmarks (defined)

- **B1** The ASTM E1049 worked example: reversals `-2, 1, -3, 5, -1, 3, -4, 4, -2` give ranges 3
  (0.5), 4 (1.5), 6 (0.5), 8 (1.0), 9 (0.5). Tested. The series is the one in the MathWorks rainflow
  documentation (citing E1049) and the counts those reported in `iamlikeme/rainflow` issue 4; both
  secondary sources, to be checked against the standard or Lee and Tjhung (Zotero `A8XFUU86`).
- **B2** Three-point and four-point equivalence (McInnes and Meehan 2008): equal counts on random
  sequences against an independent four-point implementation. Tested. Also tested: conservation of the
  total variation, and agreement of `residue="rearrange"` with Rychlik's definition on periodic
  histories.
- **B3** Constant-amplitude reproduction: a history of `m` identical cycles gives the constant-amplitude
  damage parameter and `m / N_f`. Tested.
- **B4** Experimental variable-amplitude data sets, defined but not run (they need the papers read):
  Gates and Fatemi (un-notched and notched, 2017-2018), Susmel and Tovo (2011), Wang and Susmel
  (2015, 2016), Morettini et al. (random loading, 2019 and 2021), Sonsino and Kueppers (welded
  joints, 2001), Kim and Park (1999). Each needs the load histories (or enough to regenerate them),
  the calibration constants and the lives.

## Follow-up tasks (each blocked on reading the cited sources)

The sources are in the Zotero collection "multifat #20 load histories" (key `4J6FG4BS`, under
"Multiaxial Fatigue"); PDFs are to be fetched.

- **F1** Enclosing-surface methods and the Modified Wang-Brown multiaxial rainflow: Meggiolaro and de
  Castro 2012 Parts I and II; Wang and Brown 1996 Parts 1 and 2.
- **F2** Rainflow with minimum-circumscribed-circle shear ranges: Bibbo et al. 2022; Anes et al. 2014.
- **F3** Equivalent-range extractors (PDMR, MCC, MCE, prismatic hull, maximum rectangular hull, moment
  of inertia): Papadopoulos 1998; Carpinteri et al. 2014; Araujo et al. 2014; Meggiolaro and de Castro 2014.
- **F4** Maxima over the cycle's load cases and the criteria that need them (exact-cp Eqs. 8-15).
- **F5** Accumulation variants: Haibach modified rule (source still to be located), Lv et al. 2014,
  Benkabouche et al. 2015, Okamura et al. 1979.
- **F6** Critical-plane counting and damage under variable amplitude: Susmel and Tovo 2011,
  Carpinteri, Spagnoli and Vantadori 2017 (review), Langlais 2003, Kim and Park 1999.
- **F7** Residue handling and sequence effects: Marsh et al. 2016, Anthes 1997, Glinka and Kam 1987,
  Amzallag et al. 1994, Rychlik 1987, Socie 1992, Lee and Tjhung 2012.
- **F8** Benchmark reproductions B4 and comparison with the reviews of Chen et al. 2012 and Bibbo et
  al. 2022.

Decisions: [decision-log.md](decision-log.md), `d079` to `d085`.
