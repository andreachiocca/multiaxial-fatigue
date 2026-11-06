# Development continuation checkpoint

## Milestone, solver importers and review fixes (6 November 2026)

Approved scope: `docs/superpowers/plans/2026-11-06-solver-importers-handoff.md`;
execution plan: `docs/superpowers/plans/2026-11-06-solver-importers.md`.
Started on clean branch `feat/fe-and-time-history-inputs`, HEAD `7ea8a7d`, no open PR.
Ruling: use this checkpoint as the only ledger and retain the existing branch, as the repository
and handoff require; no separate skill ledger/worktree or approval pause. The approved handoff
is the design authority; write the small execution plan and proceed inline.
Pre-flight: all three adapters share canonical row validation/atomic output; no numerical
module outside fe/history is changed. Result-state ordinals map solver step/increment identities
into ascending CSV steps; original identities and times remain in notes. Step-time resets are
refused instead of guessed into total time. Flat explicit INP nodes only; transformed part/assembly
coordinates and shell section points are refused.

Task 1: Abaqus adapter and shared output validation added test-first (initial collection failed
on missing adapter). Real MIT axis-fem RPT/INP preserved with license/source, eight nodes checked
against file values and the axial expectation. Component sentinels, separate/concatenated frames,
unit and strain declarations, conflicting element-node values and atomic refusal tests included.
Focused tests plus notation: 627 passed, 22 existing skips; black and ruff pass.
Native probe succeeds: ansys-mapdl-reader 0.56.0 imports
under Python 3.12 with numpy 2.5.4 in uv's isolated environment. No project numpy downgrade.
UniTN reachable for final full suite. Public MITC_python Ansys shell output inspected but unsuitable
as a PRNSOL stress fixture (section-point strains only); text parser will be explicitly sample-tested.
Task 1 committed as `3bed7ec`.
Task 2: PRNSOL/NLIST and complete Workbench tensor-table conversion implemented test-first;
targeted tests plus notation: 646 passed, 22 existing skips. Text Ansys support uses hand-written
samples only; the inspected real shell output is intentionally unsupported. Native source/data
inspection continues (first RST has no stresses; second is modal, to be refused by the adapter).
Next exact action: native RST compatibility/frame verification, tests, documentation and final gate.

Task 3: native adapter added after successful isolated Python 3.12/numpy 2.5.4 probe.
Real MIT pymapdl-reader `tests/testfiles/vol_test.rst` (uncommitted, 2.6 MB) converted and
re-imported; 11 native tests passed including independent global element-stress averaging.
Default suite of adapters plus notation: 661 passed, 23 skips (22 existing; one opt-in real RST).
Dependency emits one upstream PyVista offset deprecation warning. Native scope: default-option
SOLID185/186/187, one material, real static/transient, small displacement, stress only. Refuse
modal/complex/cyclic, shell/beam, layered/nondefault options, non-finite unselected midside data
and multi-material global averaging. Source verifies global stress rotation before averaging;
strain is not rotated the same way and therefore not exported. Extra `ansys` bounds the reader
to the inspected 0.56 API; project numpy remains unconstrained upward. No binary fixture committed.
Next exact action: documentation, mutation test, final quality gate on UniTN and fresh review.

Task 4 in progress: guides, README, AGENTS and decisions d093-d097 updated, with explicit real
versus sample-only evidence and unsupported DAT/ODB/Workbench variants. S13/S23 process-local
mutation fails the sentinel test (expected); clean process restored automatically. Two additional
refusals were RED then GREEN: mixed Step/Total Time and flat INP reuse across multiple instances.
User approved a local merge of origin/main PR #62; merge `1695f5b` is conflict-free. Its new native
LIUMAHADEVAN needed one method-matrix documentation row (isolation test RED then GREEN).
Protected-source diff against origin/main `d0e11c1` is empty. No regression snapshot is required:
all branch source changes remain inside fe/history, including after the upstream merge.
Focused adapter/notation/history integration checks: 682 passed, 23 skips. Documentation links,
scientific dossiers and both generators pass. A downloaded upstream Python source inside ignored
scratch was included by black's scan; retain it with a text suffix and repeat the final check.
Full suite and whole-branch reviewer remain outstanding; do not report completion yet.

Resumed after model quota: the first UniTN suite reached 100% without failure markers, but pytest
raised WinError 5 while cleaning the shared `pytest-current` temporary link; exit 1, no valid green
full-suite claim. Peak capped job memory 1,948,860,416 bytes. Repeat using a fresh private basetemp,
without deleting or altering the workstation's shared temporary directory.
The first fresh reviewer (gpt-6-astra) hit quota before completing. It found two Important history
bridge bugs: repeated node labels selected the first location silently, and timed unit cases were
read as time samples. Both are reproduced RED then GREEN with identity-provenance coverage;
`history/signal.py` now refuses ambiguous IDs and unit cases, directing users to `node_histories`.
Focused signal/nodes/notation: 661 passed, 22 existing skips; black/ruff/docs/diff checks pass.
These guard-only changes are uncommitted while a replacement fresh whole-branch reviewer
(gpt-6.1-sol, the available model after quota) completes the review. No existing-main numerical
source changed; protected-source diff remains empty. Next action: remaining reviewer findings,
one test-first fix pass, isolated full-suite rerun, final checkpoint/local commit; no push.

Whole-branch replacement review completed on `d0e11c1..d7ad4b3`: no Critical findings;
five further Important findings, seven total including the first review's two history guards.
All entered one test-first fix pass: strict CSV row width before pandas projection; atomic CSV/PVD
output and refusal of all source aliases; report-wide Abaqus local/section/material checks before
the first frame; a normalized joint stress/strain scalar-load test (also refuses constant-stress
cycling strain); and descending numeric solver-step refusal. Each new failing test was observed
before its fix. The new atomic helper preserves sources and preexisting destinations.
Ruling: the review declined scientific criterion validity, unsupported solver layouts and native
reader compatibility beyond recorded evidence; those remain explicitly outside claimed scope.
No new deferred Minor findings were reported. Focused combined checks running; full rerun pending.
The static-history regression exposed mean-subtraction roundoff in the new joint gate; subtracting
the first sample before centering preserves exact constant histories. No valid calculation changes
outside fe/history. The fix pass remains uncommitted until the isolated full-suite gate succeeds.
Focused fix-pass verification: 883 passed, 23 skips; black, ruff, documentation links,
scientific dossiers, both generated-document checks and protected-source diff pass.
Next exact action: archive the staged tree, run the full suite with private basetemp on UniTN,
then record its actual result and commit this verified fix pass locally. No further review
findings remain to implement; no second reviewer pass is needed after these covering tests.
The private-basetemp full rerun completed: 2 failed, 3066 passed, 26 skipped, 120 xfailed,
61 subtests passed (291.93 s). Failures were the two node-history life comparisons, reproducible
on UniTN alone. Investigation traced lost small-strain decimal digits to pandas' default float
parser. A new bit-exact strain round-trip test failed first; `float_precision="round_trip"` fixes
it without relaxing assertions. Local CSV/node tests: 81 passed; all 14 node tests on UniTN now
pass under the same cap. NumPy/pandas versions match across machines (2.5.4/3.0.6).
CLAUDE.md now carries the solver/history guidance corresponding to AGENTS.md (cross-agent parity).
The final frozen-tree confirmation run on WorkstationUniTN (fresh basetemp `pytest-confirm-6be012`, peak capped memory 1.91 GB) passed completely:
3069 passed, 26 skipped, 120 xfailed, 7 warnings, 61 subtests passed in 289.27 s (0:04:49).
The transient calibration callback failure in unchanged core code did not reproduce on the isolated rerun or across the full confirmation suite.
All quality gates pass: `black --check .`, `ruff check .`, `tools/check_docs.py` (1832 links OK), `tools/check_scientific_reviews.py`, `tools/generate_method_docs.py --check`, `tools/generate_workbook_docs.py --check`, and `git diff origin/main --stat -- src/multifat ':!src/multifat/fe' ':!src/multifat/history'` is empty.
Repository purged of external repository fixtures (`axis-fem` deleted). Ansys MAPDL 2026 R1 batch mode on WorkstationUniTN generated an in-house 2D plane stress holed-plate fatigue benchmark (`PLANE183`, `KEYOPT(3)=3`) under a 2-step alternating fatigue cycle ($Y$-tension 100 MPa then $X$-tension 100 MPa, rotating principal stress axes by 90 degrees). Refined FE mesh (64 hole divisions, 1.5 mm ligament sizing, 2,864 corner nodes) verified against Howland (1930) finite-width analytical solution: peak stress matches within +0.41% (304.48 MPa vs 303.22 MPa), sub-1% error across entire notch influence zone, global force equilibrium verified to +0.009% (10000.87 N vs 10000.00 N). `rst.py` updated to support PLANE182/183; both `.rst` and `.out` convert into exchange CSV and HDF5. All 46 criteria tested on FE history (38 ok, 8 safely refused). Gallery script `examples/plot_fe_benchmark_maps.py` and 6-panel 300 DPI damage map with inset zooms added (`docs/images/hole_plate_fatigue_damage_maps.png`). Real Abaqus output tests marked with xfail TODO pending Abaqus license/installation. Decision d099 recorded. Local branch `feat/fe-and-time-history-inputs` ready for user review/push instruction.

## Milestone, FE CSV exchange format (6 November 2026)

Branch `feat/fe-and-time-history-inputs`, local and not pushed. Design spec and plan:
`docs/superpowers/specs/2026-11-06-fe-csv-exchange-design.md` and
`docs/superpowers/plans/2026-11-06-fe-csv-exchange.md`; documentation `docs/fe-csv-format.md`;
decisions `d086`-`d092`. This answers the maintainer's CSV proposal on issue 19: the CSV is now the
recommended user-facing format, the HDF5 file (version 1.1) the internal store.

New or extended, all inside `multifat.fe` and `multifat.history` (no existing module changed):
`fe/model.py`, `store.py`, `select.py`, `tensors.py` (file version 1.1: instance and material labels,
split strain, temperature, step numbers, analysis metadata, US-customary and strain unit tables, Hooke
helper); `fe/csvio.py` (header, column aliases for Abaqus and Ansys, streaming `import_fe_csv`);
`fe/csvwriter.py` (standard library only); `fe/synthetic.py` (`write_tube_csv`, `--csv`);
`history/signal.py` (`from_csv` accepts the node-CSV names); `history/nodes.py` (`node_histories`,
unit-case superposition, `assess_nodes` with per-node status). Committed example
`docs/examples/fe-csv/tube_two_steps.csv`.

Verified: full pytest 2640 passed, 25 skipped, 76 expected failures (10 min 42 s); black, ruff, `check_docs`, `check_scientific_reviews` and both generators `--check` pass.
`git diff 84ea401 --stat -- src/multifat ':!src/multifat/fe' ':!src/multifat/history'` is empty (no
existing source file changed since the merge-base), so no regression snapshot was run for this
milestone; the last one (UniTN, 0 differences) covers the existing code. No whole-branch review of the
CSV work was run (model quota); a fresh review is the first task of the next session. `origin/main` has
moved on (`96de6ad`, upstream PRs #55-#59 adding criteria): merged into the branch on 6 November (`b1df4d8`, no conflicts); full pytest on the merged tree 2942 passed, 25 skipped, 111 expected failures. Two history tests needed adapting to the new criteria: the method table in `docs/load-histories.md` now lists ZHANGYAO, KUJAWSKI2021, WANG2022, MMP and SSF, and the all-methods end-to-end test gives MMP valid limits. Finding for upstream: the registry probe constants (`_PROBE_PARAMS`) break MMP's condition Sigm1 < Sig0 < 2 Sigm1, so `stress_only_models()` lists MMP as needing strain (conservative; not changed here).
Software checks only: the importer reproduces the analytical tube for
every flavour, shear convention, unit, row order and strain mode; two deliberate mutations (swapped
`S13`/`S23`, halved stress shear) fail the tests. Nothing here validates a fatigue criterion.

Deviations from the plan, recorded: x, y, z are required columns (the store needs coordinates); the
asymmetric-tensor test became a non-finite test (six components cannot be asymmetric); `TimeHistory.
from_csv` accepts `E_*` or `Eel_*` together with `Epl_*`, and refuses `Eel_*` alone (it is not the total
strain under plasticity); `Eel_*` alone is accepted by the importer for linear-elastic files only.

Open: reply to the maintainer (questions in `docs/fe-csv-format.md`); tested Abaqus and Ansys export
scripts and one real file each; automatic solver importers (Abaqus reports, Ansys result files) tested
on public data; harmonic reduction of a sampled non-proportional cycle; the pull request decision.

Next exact action: run the handoff prompt `docs/superpowers/plans/2026-11-06-solver-importers-handoff.md`
in a new session (Abaqus reports, Ansys listings, optional native `.rst`, tested on public data), a fresh
whole-branch review, then merge `origin/main` and decide on the pull request.

## Milestone, load-history bridge, issue 20 (5 November 2026)

Branch `feat/fe-and-time-history-inputs`, local and not pushed, merged with `origin/main` (`84ea401`). Design
spec and plan: `docs/superpowers/specs/2026-11-05-load-history-bridge-design.md` and
`docs/superpowers/plans/2026-11-05-load-history-bridge.md`; documentation `docs/load-histories.md`;
decisions `d079`-`d085`.

New package `src/multifat/history/`, imported by no existing module: `signal.py` (time-history contract
from arrays, CSV, `fe` file or a documented function, with a hashed record), `channels.py` and
`proportionality.py`, `rainflow.py` (ASTM E1049 three-point, `half` or `rearrange` residue),
`events.py`, `accumulate.py` (caller-supplied Basquin `SNCurve`, Miner, endurance), `eligibility.py`
(method matrix, `IneligibleHistoryError`), `assess.py` (result record) and `example.py`
(`python -m multifat.history.example`). Exact for proportional histories; a non-proportional history is
refused unless an explicit, recorded approximation flag is set. Multiaxial counting and the
equivalent-range extractors are follow-up tasks F1-F8 (sources: Zotero collection `4J6FG4BS`, 42 papers,
PDFs not fetched). The ASTM E1049 example counts are confirmed by two secondary sources only.

Verified: full pytest 2476 passed, 25 skipped, 76 expected failures, plus one isolation test of the FE
work that forbade any outside import of `multifat.fe` and was corrected to name `history` as its
intended consumer (the two isolation files rerun green); black, ruff, `check_docs`,
`check_scientific_reviews` and both generators pass. Regression snapshot (`--skip-calibrate`,
WorkstationUniTN): `origin/main` (`84ea401`) against the branch at `1be7d8b`: 2,062 output files and 481
log files, 0 differences; later commits change only a test; no existing source file changed. A
whole-branch review (0 critical, 3 important, 12 minor) led to fixes with failing tests first: settable
proportionality tolerance, `rearrange` as the default residue for a repeating block (`half` under-counts
it by about a factor 2 and is noted), `event_damage_parameters` for native-life methods, a strain-path
proportionality check, a scope-restricted column in the method matrix, and hardening of the record,
counting guards and overflow handling. Deferred minors: strain units other than dimensionless are not
converted; `record()` with a `.json` suffix writes `x.json.json`. The ASTM E1049 example counts rest on
secondary sources only.

Next exact action: the pull request decision, and the FE CSV exchange format proposed by the maintainer
(issue 19 comment), see the discussion with the maintainer.

## Milestone, FE local-history interface, issue 19 (5 November 2026)

Branch `feat/fe-and-time-history-inputs` (renamed from `feat/fe-input-interface`) from `main` (`538ebc3`), merged with `origin/main` (`84ea401`) on 5 November, local and not pushed. Design spec and
plan: `docs/superpowers/specs/2026-11-05-fe-local-history-interface-design.md` and
`docs/superpowers/plans/2026-11-05-fe-local-history-interface.md`; format specification
`docs/fe-input-format.md`; decisions `d072`-`d078`.

New package `src/multifat/fe/`, imported by no existing module: `tensors.py` (component orders,
shear convention, units, frames), `model.py`, `store.py` (HDF5 file, chunked reading), `select.py`,
`vtk.py` (`.vtu`/`.pvd` importer), `ranges.py` (explicit pair, longest chord), `bridge.py` (two-state
bridge, `requires_nonlocal` eligibility), `synthetic.py` (analytical thick-tube oracle and
`python -m multifat.fe.synthetic`). `h5py` is the optional extra `fe` (`pyproject.toml`, `uv.lock`;
the CI test job installs it). No existing source, workbook, calibration record or result changed:
`git diff main --stat -- src/multifat ':!src/multifat/fe'` is empty.

Verified: the new tests (tensors, model, store, select, synthetic, vtk, ranges, bridge, isolation,
large input) pass; the opt-in test of 1,000,000 locations by 16 load cases (float32, 64 MB chunk
budget) passed in 81 s with peak memory below four times the budget. The tube is checked against
its boundary conditions, radial equilibrium and torque and axial resultants, and every importer
convention is recovered to 1e-10 relative. These are software checks, not validation of a fatigue
method. Full pytest after the review fix pass and the minor fixes: 1835 passed, 28 skipped, 41 expected failures (the documented
misprints); black, ruff, `check_docs`, `check_scientific_reviews` and both generators pass. A first full
run caught three retired-synonym findings in the new modules (`mat`, `spec`, `dataset`), renamed
before the first docs commit. A whole-branch review then led to five fixes with failing tests first:
the VTK base64 layout VTK actually writes (one stream; the separate-header layout is still read),
a `complete` attribute so a crashed export is refused instead of read as zeros, refusal of
strain-needing methods on a stress-only history, required stress and length units, and an HDF5 chunk
size that follows the number of load cases (the former default failed above about 21,800 load cases).
The review's minor findings were fixed too (rotation arrays validated, multi-part `.pvd` rejected,
source conventions recorded, newer 1.x metadata fields tolerated, unchunked foreign files read,
duplicate IDs refused on open, NaN screening values ignored, finite component check); the byte
budget wording, the hydrostatic blind spot of the longest chord and the dropped per-location kind and
surface flag are documented instead.

Regression snapshot (`tools/regression_snapshot.py run --skip-calibrate --jobs 4`, run on
WorkstationUniTN): `main` (`538ebc3`) against the branch at `8d02c5c` (the later commit changes only
`fe/`, tests and docs): 1,708 output files and 409 log files, 0 differences; no existing source file
changed. The calibration stage was skipped, so calibration records are not re-derived by that run.
Remaining: the pull request, decided by the maintainer. Follow-ups, not started: PyAnsys/DPF, Abaqus `.odb`, CalculiX and Exodus adapters;
issue 20 (time-history contract, remaining range rules, maxima over all load cases).

Next exact action: review the branch, then decide on the snapshot run and the pull request.

## Issue #13 methods merged; integration audit (5 November 2026)

PRs #45 to #49 (GAM, NINIC, GOLOS, MWHS, LIU2025) and #50 (homepage figures) are merged into
`main`.

An integration audit, `--all` with the new models, found four problems, all fixed before merging:
1. GOLOS was misclassified as needing strains.
2. Load-scaling searches failed on trial loads beyond a criterion's domain. Such loads now count
   as +inf damage and the bracketing accepts an infinite end.
3. LIU2025 could divide by zero in the Serensen branch. It now raises `LoadingScopeError`.
4. MWHS was fitted on tension tests although its source calibrates on torsion. It now uses the
   new `design_loading = "torsion"` option.

Checks: full suite green and CI lint green. The regression snapshot shows no change to existing
model outputs.

Next exact action: implement further #13 candidates as their papers arrive, refreshing the
homepage figures each time (rule in `docs/adding-methods.md`). Meanwhile, continue the #32
guides.

## Issue #13 new methods from maintainer-supplied papers (5 November 2026)

Five criteria were implemented from the original PDFs, each with a guide, a paper test, a catalogue
record and an audit row. They are stacked PRs; merge them in the order listed.

| Order | Method | PR | Branch | Paper check |
| --- | --- | --- | --- | --- |
| 1 | GAM (N08) | #45 | `feat/gam-criterion` | Tables 1-4: in-phase rows exact; out-of-phase rows used a phase-free peak, contrary to Eq. (11) (xfail) |
| 2 | NINIC (N11) | #46 | `feat/ninic-criterion` | Tables 2-4: 29/32 rows within 0.08; Eq. (4) valid only for S_e/T_e >= sqrt(2) |
| 3 | GOLOS (N13) | stacked | `feat/golos-criterion` | Eq. (30) identities |
| 4 | MWHS (N28) | stacked | `feat/mwhs-criterion` | Fig. 14 reproduced, except 90 deg of (b) |
| 5 | LIU2025 (new record N58) | stacked | `feat/liu2025-criterion` | Fig. 3 reproduced, except 2 points |

Notes on the GOLOS and LIU2025 branches:
- The GOLOS branch adds `LoadingScopeError`. The pipeline reports a model that raises it as
  skipped, with eligibility `out_of_scope`.
- GOLOS is out of scope on every supplied workbook.
- LIU2025 is out of scope on 42CrMo4_QT and Ti6Al4V_Berto because of the estimated Tu.

Checks: full suite green, apart from one shear-rule registry test, which was fixed and re-run.
Regression snapshots show no change to existing model outputs.

Next exact action: implement the remaining #13 candidates as their PDFs arrive. Requested:
- Zhang & Yao 2010;
- Wang et al. 2022;
- Kujawski 2021;
- Papuga & Fojtík 2017;
- Anes et al. 2014;
- Lu & Martínez-Esnaola 2021;
- Stefanov 1993;
- Liu & Mahadevan 2005.

Meanwhile, continue the #32 guides (LI, LIU2021, MGSE_YU, MGSE_ZHU, MSWT, MMCCM).

## Issue #32 guides MKBM, VSE, GSE, GSA; issue #13 candidates (5 November 2026)

- PR #44 (`docs/guides-strain-energy`) adds four guides: MKBM, VSE, GSE and GSA.
  - It completes the `mkbm`, `liu_vse` and `gsa` docstrings; 19 modules remain in
    `PENDING_DOCSTRINGS`.
  - Checks: guide, notation and method tests pass (536); black, ruff and `check_docs` are clean.
  - Documentation only; the regression snapshot was not run.
- Issue #13: the prioritised candidate list is posted as a comment on the issue.
  - Batch 1: N08, N06, N11, N05, N13, N28, N26.
  - Batch 2: G05, N24, N32, G32, N16.
  - Blocker: coding waits for the original PDFs from the maintainer.

Next exact action: #32 guides for LI, LIU2021, MGSE_YU, MGSE_ZHU, MSWT and MMCCM, then the
invariant and integral families. Implement the #13 candidates as their PDFs arrive. TVEIT2025
stays untouched.

## Open work after PR #33 (5 November 2026)

Open PRs, each from `main` at `2efee9f` and independent:

- #34 (`claude/multiaxial-fatigue-repo-ml3b4w`): calibration integrity. Design flags of the tests
  behind calibrated plane constants; `k_MGSE` from Zhu Eq. (12) regardless of the model selection;
  L_LI searched in [0, 20]; constants stopped by an arbitrary search bound are `unidentified` and
  not used; run-report notes for flat master curves. Changes calibration outputs (see CHANGELOG,
  *Changed numerically*); `regression_snapshot.py` outputs identical, `calibrate.log` differs as
  documented.
- #35 (`docs/method-guides-batch2`): guides for BP, CAIM, PCR, PI, MOREL (15/35 methods).
- #36 (`feat/survival-probability`): issue #16, `--survival-probability`; closes #16.

Each branch was checked with the full suite, black, ruff, `check_docs` and the regression snapshot
against `5cd5124`. Issues #14 and #32 carry progress comments; #16 is closed by #36.

Open scientific points, not changed:
- `zhu_k_mgse` averages Eq. (12) over N = 1e4..2e6, whereas `docs/calibration-source-audit.md`
  records 5,000-50,000 cycles for metals in Zhu et al. (2018), Section 4.1. Needs the paper.
- Flat scatter minima (e.g. Ti6Al4V k_MGSE on [2, 5], a_MSWT on [0.5, 1]) still take the smallest
  point of the plateau; they could be treated as `unidentified` like bound hits.
- TVEIT2025 within-cycle failure interpolation (see below): owned by another contributor.

Next exact action: after the three PRs are merged, continue #32 with GSE, GSA and VSE, then the
remaining critical-plane energy methods; then issue #18 (fatigue-limit vs fatigue-life workflows).

## Issue #32 method guides, first ten methods (5 November 2026)

Branch `claude/multiaxial-fatigue-repo-ml3b4w`, from `main` at `5cd5124`. Documentation only;
TVEIT2025 is deliberately untouched (another contributor is working on it).

Done: `docs/method-guides/` with a template, an index of the shared conventions and guides for
FIN, FS, SWT, MATAKE, MCD, MWCM, CROSSLAND, SINES, DANGVAN and PAPADOPOULOS. Each separates the
published formulation from the repository implementation (`Repository choices`), states the
original and the repository calibration, and has a worked example whose numbers are asserted
against hand calculations and whose printed output is compared by `tests/test_method_guides.py`.
Module docstrings were added or completed for the same ten modules. The new test also ratchets the
`Repository choices:` paragraph of the method modules: 25 modules are listed as pending.
README, user guide (FAQ and test counts), `adding-methods.md` and CHANGELOG link the guides.

Verified: the guide CLI commands run on `Al7075_T6` (all ten models `ok`; `--constants` reuse of
the sidecar works); full test suite, black, ruff, `check_docs`, `generate_method_docs --check`;
`regression_snapshot.py` against `5cd5124`: 0 differing outputs and logs.

Design flags of calibrated plane constants (affects issue #14), fixed after PR #33 was merged:
the tests of the scatter objective behind `k_FI`, `k_FS`, `L_LI`, `k_MGSE`, `a_MSWT`, and the
Basquin tests behind a Zhu `k_MGSE`, were exported as `Flag = 1`. They are now design rows; 410
rows of the five life datasets change flag, nothing else (6.7e-14 relative rounding). The guides
no longer carry the warning.

Remaining (issue #32): guides for the other registered methods (next: BP, CAIM, PCR, PI, GSE,
GSA, VSE, MOREL, then the rest; TVEIT2025 only after its owner's work lands), and the 25 pending
docstrings. Next exact action: write the BP and CAIM guides from the template.

## Main after Tveit integration (29 October 2026)

PR #28 was merged into PR #27, and PR #27 was merged into `main` at
`3555826`. PR #21 was closed without merge; no overlapping PR remains open.
The documentation follow-up starts from this merge. It records the current
CLI, registry and comparison-capability flow, distinguishes the historical
`output_quantity` proposal from the implemented API, and corrects generated
catalogue wording for active, partially implemented TVEIT2025. No model,
calibration record, experimental dataset, benchmark or numerical result is
changed. The within-cycle Tveit failure-interpolation question recorded below
still requires scientific approval before any numerical correction.

## PR #27 native-life reporting review (29 October 2026)

The focused follow-up branch `fix/pr27-native-reporting-review` starts from
PR #27 head `05fa153`. It evaluates native life even when no reversed-axial
design set is available, while retaining Morel's explicit insufficient-design
diagnostic. The native-life summary now distinguishes calibration/design rows,
assessment failures and run-outs. Focused regression tests cover both cases.
The Tveit continuum equations, numerical constants, scientific benchmarks,
material records and dataset files are unchanged.

The Tveit within-cycle failure interpolation in
`src/multifat/models/invariant_methods/tveit2025.py` appears to add one sampled
increment to the elapsed fraction. Correcting it would change native numerical
results and requires Scientific Reviewer or numerical-owner approval and
updated reference checks. PR #27 and the overlapping original PR #21 must not
both be merged.

## Issue #17 comparison-error reporting (29 October 2026)

The focused follow-up starts from merged `main` at `3555826`. It adds
`summary_error_stats_by_role.csv` to distinguish design residuals from other
assessed failures while leaving the existing aggregate statistics and CSV
fields unchanged. Analytical tests cover the DP error's base and sign,
nonpositive/nonfinite masking, role means and counts, and run-out exclusion.
`Flag=1` is an assessment row, not proof of source- or batch-independent
validation. Issue #14 retains the held-out validation protocol. No fatigue
equation, calibration record, experimental data or numerical benchmark changes.

## Janssens final workbook (28 October 2026)

The final author-supplied AISI316L workbook replaces the previous file. The R = -1 tensile fatigue limit is 278.5 MPa, the first uniaxial run-out is 17,225,000 cycles, and the non-proportional torsional phases are 90 degrees. The overlay has been reduced to 15 independent operations. The updated dataset fingerprint is `a63f9325d41155bb7611b2c6272e0bbee9a46e33fefe4bc33faa47f8470f4750`. Older checkpoint notes describe superseded source revisions.


## Issue #22, approved fatigue-limit selection policy (28 October 2026)

The maintainer confirmed that fatigue-limit values measured at the same stress ratio
should be averaged, separately for tensile and torsional loading. The implementation
on PR #23 already follows this policy. For each available R = -1 parameter,
`Sigm1` or `Taum1`, the revised Summary importer calculates an equal-weight
arithmetic mean of eligible positive, finite source values. It does not average
across different stress ratios or loading modes. Entries representing bounds rather
than measured values are not used as exact limits. The selected value retains its
contributing Excel rows, references, authors, stress ratios and reference cycle
counts. Missing usable R = -1 data cause a warning and leave the property
unavailable. The original workbook and calibration JSON records are not edited.

On the final Janssens-supplied AISI316L workbook, the selected values are
278.5 MPa in tension and 144 MPa in torsion at R = -1. Braun's R = 0 stress
range (387.8 MPa at 10^6 cycles) is excluded from these averages. The author-
corrected workbook and its 15-operation overlay yield 39 failures and 13
run-outs. Both push and pull-request CI runs passed on commit `2178636`,
including lint and the complete test matrix (Ubuntu 3.10-3.13 and Windows
3.12). The previous AISI316L overlay remains available as a historical record.

The source identity and batch-compatibility limitations remain separate data
review questions. In particular, the CX18 Summary names 30NCD16, and a Ti6Al4V
source uses a different material label. The importer retains their source
provenance and reports identity conflicts but does not establish metallurgical
equivalence. Fatigue strengths reported at incompatible reference lives or as
stress ranges must not be interpreted as independently validated endurance
limits. No unverified Goodman conversion or other mean-stress correction is
automatically applied. Issue #18 remains separate and concerns distinct
fatigue-limit and fatigue-life assessment workflows.

## Milestone, multifat architecture integration (22 October 2026)

Branch `dev/tg-developments`, local and not pushed. Upstream `main` at `ba84152` (PRs #10 and #11)
is merged as `aca31b0`. Further commits on this branch: the overlay layer for the five workbooks
edited after supply (`5b1a852`), the explicit Findley assessment for FIN and FIN_ext (`994981e`),
and this combined `AGENTS.md` with the checkpoint. The default paths are unchanged.

Carried over: the N04 scientific dossier (five files, check tool, CI step). Susmel N04 stays
implemented here with the MWCM record (enabled, tested), as decided by the user; upstream records
it as `not_implemented`, so the catalogue differs from upstream on that record. The explicit
Findley option `--shear-amplitude-algorithm` for FIN and FIN_ext, with its five tests rewritten
against this architecture. Not carried: the upstream calibration store
(`fatigue/eval/calibration_store.py`), because this code calibrates Findley inside each run.

Verified on the final code (base `206748c`, the pre-round HEAD; new `6243442`, code-identical to the
final commit `84edf55`, which changes only this file):

- Full pytest at `84edf55`: 1633 passed, 27 skipped, 41 expected failures (documented misprints).
  Black, ruff, `check_docs`, `check_scientific_reviews` and both generators pass at `84edf55`.
- Calibration through the multifat-calibrate CLI, including single-workbook copies: identical to the
  base for all 15 workbooks. Overlays follow the workbook by content, so a copied `8_AISI316L` gives the
  base constants.
- Pipeline over all 15 workbooks (folder import): every calibrated constant and run report identical.
  The differences are the expected categories only: timestamps, git revision, the five overlaid
  workbooks' hashes, the new `assumptions_applied` field, timing, output paths and the new
  `[--no-overlay]` usage text.
- Regression snapshot: 1,021 run logs identical after traceback-frame normalisation; exit codes
  identical; the five overlaid workbooks identical cell by cell.
- Plots: 78 run outcomes identical after path normalisation.
- Differential harness: 34,070 cases, no value, warning or message differences.

An earlier parity run on `882d662` (before the overlay relocation) gave the same results.

Done in this round: MIT licence (step F, `f590efa`), the Calibrations loading test (step E, `bb546e1`),
the explicit Findley option with rewritten tests (step B, `994981e`), the [project decision log](decision-log.md)
recording `d001`–`d071`.

Still open: pushing and the draft PR (decided by the maintainer, not done here);
`docs/workbook-migration-evidence.json` is not regenerated (no generator found in the tree).

Next exact action: review the draft PR text and decide whether to push
dev/tg-developments and open the PR.

Path map: `fatigue/` is `src/multifat/` (see `AGENTS.md`).

## Current status, Findley Phase 3B and Susmel N04 documentation

Updated 20 October 2026. Findley Phase 3B was merged into `main` through
[PR #11](https://github.com/andreachiocca/multiaxial-fatigue/pull/11) at
`356357bf162991fb24583b3376490049c9655ef0`. Its numerical implementation,
regression tests and checkpoint remain intact. The Susmel N04 scientific dossier
is in [PR #10](https://github.com/andreachiocca/multiaxial-fatigue/pull/10),
which remains open for human review. This merge reconciles the two checkpoint
entries without changing the five dossier files, Findley numerical source,
calibration records, experimental data or reference results. Susmel remains
`not_implemented` and `not_tested` in the method catalogue. No Phase 4 or Susmel
implementation is started here.

The next action is human review of PR #10 after its documentation, regression
and GitHub quality checks pass. Do not merge automatically.

## Completed phase, Findley production amplitude assessment (Phase 3B)

Updated 20 October 2026. Branch `feat/findley-amplitude-assessment` starts from
main `ab0e828a216622bc8a26ed39c04b1377ebf8aeae`, the merge of PR #9. PRs
#7, #8 and #9 are merged. The baseline 101-test suite passed locally, with
three pre-existing skips. PR #10 is open for the independent Susmel N04
documentation and has no overlapping Findley source ownership. Direct Git
transport is unavailable in this container; source files were reconstructed
from verified GitHub blobs, with exact-byte checks on the numerical files and
four archived workbooks needed by the tests. No Susmel file is modified.

Phase 3B integrates Script 1 only for FIN/FIN_ext. Its omitted
`--shear-amplitude-algorithm` option uses the unchanged production FIN model
and verified legacy MCE calibration path. Explicit `mcc_harmonic`,
`mce_harmonic` and `rss` use the Phase 2 evaluator through the established
Findley calibration callable. The original stress sine/cosine tensors are
required before endpoint extraction; unsupported cases are excluded with a
diagnostic. RSS resolves to harmonic MCE. Other methods reject an explicit
override. No automatic calibration or general registry override is introduced.

Assessment lookup now optionally requires the recorded evaluation path as well
as the resolved algorithm. Default requires `legacy`; explicit choices require
`harmonic`. This is an additional strict assessment constraint, not a
fingerprint exception or a change to the calibration API's existing default
lookup. Material, failure-input, grid, fitting, procedure, version and source
checks remain in force. The production FIN, harmonic evaluator, shared
evaluator, adapter, registry and amplitude library remain byte-identical to
main, preserving their stored formulation identities.

Default CSV names, columns and numerical values are unchanged. Explicit CSVs
use `FIN.mcc_harmonic`/`FIN.mce_harmonic` in the filename and `Method` field;
the `_ext` suffix remains available. Each Findley CSV has a `.meta.json`
sidecar with selected option, resolved algorithm/version, evaluation path,
calibration ID, input/formulation provenance, settings and CSV checksum.
The existing comparison reader accepts these files without schema changes.

Five new focused assessment tests pass. They cover path-specific record
selection and ID rejection, legacy-record isolation, MCC/MCE/RSS output
metadata, FIN_ext labels, unsupported methods and unknown identifiers, raw
history rejection before fallback, comparison-reader compatibility and an
independent per-plane SVD reference for non-proportional production assessment.
The complete suite passes 106 tests, with the same three expected skips.
Default AISI316L FIN CSV and unrelated FS CSV each compare byte-for-byte with
the original main runner on the archived workbooks. Documentation checks pass
1,396 local links; the generated method catalogue check passes 161 records,
117 counted candidates and 51/51 review entries. `git diff --check` passes;
a protected-source diff confirms unchanged model, evaluator, adapter, registry,
experimental, calibration, result and scientific-review paths.

Temporary Script 0 and Script 1 CLI runs on AISI316L and Al7075_T6 produced
separate compatible MCC/MCE records and result files at the production
522-plane grid. AISI316L exported 44 rows and Al7075_T6 101 rows. Both
algorithms selected `k_FI=0.75` for each material under the coarse
five-point, one-refinement calibration search. Their production-grid selected
damage values coincided on all exported rows, despite distinct algorithm and
record identities. A matched `n1=3, n2=6` Al7075_T6 integration check found
19 differing rows; `Plain_multi_32` gave MCC 170.000000 MPa and MCE
195.830052 MPa. The latter case is independently checked by per-plane SVD.
These are grid-dependent numerical observations, not experimental predictive
validation or evidence that the amplitude definitions are interchangeable.

Ruff 0.14.2 is not installed and package download was blocked by the container's
HTTP proxy. Local Ruff was not run; the GitHub quality workflow passed Ruff,
tests and catalogue/documentation checks. No source workbook, historical
calibration JSON or saved result CSV was modified. Changed tracked files are
Script 1, the calibration store, eligibility, CSV exporter, README, the
shear-amplitude, calibration-workflow and architecture guides, this checkpoint,
and the new `tests/test_findley_assessment.py`.

The first verified commit is `3cbcb3722b6930d6bdbf175e910efbc18c035097`
on [PR #11](https://github.com/andreachiocca/multiaxial-fatigue/pull/11).
The published tree matches all ten tested local file blobs and differs from
main only in those files. Both GitHub quality workflows passed for this
commit, push run 39387667205 and PR run 39387737462. This following
checkpoint-only status update changes no numerical source, test or calibration
fingerprint. PR #11 was subsequently merged into `main` at
`356357bf162991fb24583b3376490049c9655ef0`. Phase 4 has not started.

## Documentation phase, Susmel N04 scientific review

Updated 20 October 2026. This documentation branch starts from verified `main`
`ab0e828a216622bc8a26ed39c04b1377ebf8aeae`, the merge of Findley Phase 3A
PR #9. No open repository PR was found at the baseline. The earlier Phase 3A
section below is retained as a historical record; its statement that PR #9 is
open is superseded by the merge. Findley Phase 3B and other numerical development
are outside this task.

The supplied N04 dossier is integrated as five byte-identical files under
`docs/scientific-reviews/N04-susmel/`. The new reviews index explains the common
LaTeX, Markdown, BibTeX, JSON and PDF convention. The catalogue source links N04
to the dossier, records the original-paper review and independent B1-B8 analytical
benchmarks, and retains `not_implemented` and `not_tested`. Only its derived
catalogue documentation is regenerated. A dependency-free dossier check and CI
step validate required files, method identity and distinct benchmark IDs. Temporary
LaTeX outputs are ignored. No fatigue criterion, calibration record, benchmark
value, experimental data or numerical reference result is modified.

Publication screening found no private PragTic code, third-party article PDF,
private experimental data, embedded PDF attachments, internal paths or credentials
in the five public files. The four citations and DOI metadata were checked against
the supplied source register and publisher records. The source dossier describes
the repository at an earlier commit; its dated implementation assessment is
preserved, not presented as a current code audit.

Local checks completed before publication: `python tools/generate_method_docs.py`
and `--check` passed for 161 records; `python tools/check_scientific_reviews.py`
passed for one dossier; a structured JSON check confirmed B1-B8, both B7 outcomes
and undefined B8 values. `latexmk -pdf -interaction=nonstopmode -halt-on-error
-outdir=/tmp/wg43-n04-build N04-susmel.tex` completed with BibTeX and no
unresolved citations or references. The supplied and rebuilt five-page PDFs have
identical extracted text and rendered pixels, so the supplied PDF is retained.
The first published documentation commit is
`bdd72e0657c78b142bafbe3b225c3702c1dabd65`. Its branch
[Quality checks run](https://github.com/andreachiocca/multiaxial-fatigue/actions/runs/39385413635)
passed both jobs. GitHub Actions ran `python tools/generate_method_docs.py
--check`, `python tools/check_docs.py`, `python tools/check_scientific_reviews.py`,
`python -m pip check`, `python -m ruff check .` and `python tools/run_tests.py`.
These are documentation and software regression checks, not independent Susmel
predictive validation. A protected-path review of the 16-file commit confirmed
that numerical source, calibration records, experimental inputs and reference
results were untouched.

PR #10 was subsequently opened for review. Its documentation branch is now
updated against the Findley Phase 3B merge. Do not merge automatically or begin
Susmel development.

## Previous phase, configurable Findley calibration (Phase 3A)

Updated 19 October 2026. Branch `feat/findley-amplitude-calibration` starts at
verified main `1ffb8896d1447e9f9f40f37f0c8c7dd08a563251`, the merge of PR #8.
The initial working tree was clean and no Phase 3A PR was open. Phase 1 library,
Phase 2 interface and experimental-data changes are present. PR #8 is merged;
its historical open-PR status below is superseded. Baseline suite passed all
83 tests before editing (12.8 s, no skips).

The objective is selectable Findley calibration in Script 0, method-scoped JSON
provenance and compatibility, while production assessment remains MCE. Reviewed
Script 0, calibration/fitting, store, eligibility, all loader call sites and the
verified Findley evaluator. FatPy contributing, testing and code-style pages and
current repository were reviewed. Adopt explicit typed configuration, focused
components, unit/integration checks and concise API/theory documentation. No
FatPy source, dependency or package/test-framework migration is introduced.

Implemented at this milestone:

- `fatigue/eval/findley_calibration.py` resolves default/RSS to canonical MCE,
  validates method scope and original stored histories, and adapts the Phase 2
  evaluator for both coefficient optimisation and master-curve fitting.
- Script 0 accepts `--shear-amplitude-algorithm`. Default numerical branches and
  observation selection remain intact. Explicit requests check raw case histories
  before `extract_tensor` can synthesise proportional coefficients. Batch preflight
  skips unsupported methods or missing axial histories; selected coefficient cases
  are also checked before optimisation.
- `calibration.py` uses the configured model and persists its effective identity.
  `calibration_store.py` adds optional FIN-scoped amplitude metadata within schema
  1, includes it in record signatures, validates optimisation/curve consistency
  and checks algorithm compatibility before choosing the newest record.
- The original formulation hashes are preserved. Explicit harmonic records add
  dependencies on the calibration adapter, Phase 2 evaluator and shared amplitude
  library. Metadata-free MCE is accepted only for the verified historical FIN hash
  `f65f1330...953acb7b`, matching current source plus all previous provenance/input/
  settings checks. Removing MCC metadata cannot confer legacy MCE compatibility.
- Ordinary loader calls require MCE. Script 1, production FIN, evaluator, adapter,
  registry, other criteria, historical JSON and Excel inputs are unchanged.

The first verified milestone is published as
`36ff0a12225ca9bfeb6ebb4d146f3dfdf4112132`. Actual final checks: all 18 new
`tests/test_findley_calibration.py` tests pass. The full suite passes all 101 tests
(83 existing plus 18 new, 15.7 s, no skips). Repository Ruff, documentation links
(1,394), generated method catalogue, CLI help and diff checks pass. New checks include actual historical JSON acceptance,
unknown-provenance rejection, canonical metadata/record reuse, coefficient and
curve consistency, independent non-proportional objective SVD, dependency/version
changes, batch skips, and exact Script 1 CSV comparison with a newer MCC present.
An explicitly selected MCC ID and an MCC-only store cannot supply MCE assessment.
Missing coefficient-subset histories also produce a safe batch skip. Explicit
coefficient reuse checks amplitude compatibility and equivalent reruns reuse IDs.
The first attempt at the new missing-torsion test used an obsolete case number
from a historical numerical fixture; the test now selects an actual failed torsion
observation. No scientific selection rule or saved reference value was changed.

Complete Script 0 runs used current `Al7075_T6`, n1=3, n2=6, n_grid=5, n_refine=1,
plane_stride=1, parameter_stdnum=2 and master_stdnum=0 in a temporary directory.
All three selections (default, explicit MCE, MCC) gave k_FI=0.75; final master
A=A_surv=724.2113871168674, b=-0.10692522368903569 and STD=0.5061625852437134.
Coefficient fitting used 23 failures (11 axial, 12 torsional); the master curve
used 11 axial failures. Fourteen run-outs were excluded. Selected observations
and master-curve damage series agree exactly between algorithms on these
proportional calibration cases, while their metadata/record identities differ.
The coefficient-stage default diagnostic fit has A=634.6482916208016,
b=-0.10977780510589384, A_surv=445.01621466900093 and STD=1.6167214068576232.
Original-main Script 0 and calibration collector source were separately executed
on these same data/settings. Coefficient result, selected cases, master fit and
damage series match the new default exactly. Production lookup selects MCE even
with a newer MCC in that temporary archive; an explicit MCC ID is rejected.

Modified files are Script 0, calibration.py, calibration_store.py, eligibility.py,
the new Findley calibration adapter, its tests, docs/shear-amplitude.md,
docs/calibration-workflow.md, docs/architecture.md and this existing checkpoint.
A protected-path diff confirms unchanged Script 1, other production/review scripts,
all models, shared evaluator, importer, plane sampling, historical calibrations,
Excel inputs, results and numerical fixtures. All method formulation fingerprints
on the original path remain unchanged. The documentation records exact metadata,
compatibility rules, coarse example and the FatPy development/interoperability note.

Phase 3A is complete in [PR #9](https://github.com/andreachiocca/multiaxial-fatigue/pull/9),
open and unmerged on `feat/findley-amplitude-calibration`. The verified implementation
and documentation commit is `0105b695ca42b6f73995b089cced0bd8f3797230`. The published
tree matches the tested local tree exactly. Both GitHub workflows passed for this
commit: push run 39268373284 and PR run 39268498442. This following checkpoint-only
status commit does not change numerical source, tests or metadata identities.

Next action is human review of PR #9. No Phase 3A tests remain outstanding. Keep
the PR open without merge. Only complete single-harmonic histories and finite
plane grids are supported; existing float64 normal-extrema limitations remain.
The example uses a coarse verification grid and does not establish predictive
validation. No stored calibration or reference result is replaced. Phase 3B,
including MCC production life prediction, requires separate explicit approval.

## Previous phase, controlled Findley harmonic amplitudes

Updated 19 October 2026. Branch `feat/findley-shear-amplitude` starts from verified
main `71eaac38ef619144c7ea055edf3b5ef3835ed7fe`, the merge of Phase 1 PR #7.
PR #7 is merged; the historical open-PR status below is superseded. There was
no existing Phase 2 implementation or open PR. The initial working tree was clean.
The dataset work from PR #6 and independent amplitude library remain present.

The objective is Findley-only numerical damage evaluation with explicit MCC/MCE,
without changing production life prediction or calibration compatibility. Review
of `fin.py`, `cp_adapter.py`, `evaluator.py` and calibration fingerprints found
that editing any shared production source would invalidate stored curves, even
for unchanged algorithms. The smallest safe interface is therefore the separate
`fatigue/models/cp_methods/fin_harmonic.py::evaluate_findley` function. It uses the
existing case arguments and returns `CaseResult`. No registry entry is added.

Omitted selection delegates exactly to the unchanged FIN adapter. Explicit MCC,
MCE or RSS selections require both harmonic stress tensors. They rotate the actual
coefficients on every candidate plane, call the shared amplitude library, reuse
the existing signed harmonic normal-stress extrema, and maximise Findley damage.
FIN and FIN_ext consequently return the same maximum. No algorithm is hidden in
material parameters. No other method can receive this option through its adapter.

Verification at this milestone: all 70 baseline tests passed before editing.
The 13 new tests in `tests/test_findley_harmonic.py` first failed on the missing
module, then passed; focused Ruff checks also passed. Tests include analytical
circles/ellipses, a two-plane ranking reversal, independent SVD/projection checks
on 30 random planes, proportional equality, phase/basis invariance, mean shear
and signed mean normal stress, endpoint fallback and explicit-request rejection,
extreme shear scales, identical FIN variants, captured main outputs for FIN/FS/
Matake/McDiarmid and unchanged formulation fingerprints for these four methods.

The first verified milestone is published as
`176ce53600a427865b9822e6294a641af087381b`. The changed files are the new Findley
diagnostic module, its tests, `docs/shear-amplitude.md` and this checkpoint.
The complete post-change suite passed all 83 tests (70 existing plus 13 new),
without skips, including the final float64 conversion rejection tests. Focused
tests, repository-wide Ruff, documentation links, the generated method catalogue
and diff checks passed. A protected-path diff confirms that these four files are
the entire change: JSON records, experimental inputs, numbered scripts, registry,
other criteria, production FIN, evaluator, adapter and calibration store are
unchanged. Production FIN/FS/Matake/McDiarmid fingerprint snapshots pass exactly.

The numerical note documents inputs, selection, example damage values and the
calibration boundary. Circular/elliptical prescribed-plane damage is 12.4/16.542136
and 12.4/13.580340 (MCC/MCE). In the ranking example MCC selects normal x with
13.806248 and MCE selects normal z with 15.542136. Default proportional and
non-proportional snapshots match main exactly; explicitly selected MCE agrees
within float64 rounding. FIN and FIN_ext agree for both algorithms.
Known limitations are complete single-harmonic histories and the supplied finite
plane grid; normal-extrema float64 arithmetic retains existing limitations. No
experimental validation or compatible MCC damage-life calibration is claimed.
Phase 2 is complete in [PR #8](https://github.com/andreachiocca/multiaxial-fatigue/pull/8),
open and unmerged on `feat/findley-shear-amplitude`. Verified implementation and
documentation commit is `50dabaf1c7f7c38b0af34ec8074d1f60758d0093`. The published
Git tree matches the tested local tree exactly. Both GitHub quality workflows
passed for that commit: push run 39248633505 and PR run 39248694751. This later
checkpoint-only status commit does not alter numerical code or test definitions.

Next action is human review of PR #8. Keep it open without merge. No tests remain
outstanding for Phase 2. Phase 3 requires separate approval and resolved algorithm/
definition-version metadata with dependency-scoped calibration compatibility
before production life prediction can use alternative amplitudes. Do not edit
JSON records, enable a global option or apply MCE curves to MCC damage.

## Previous phase, independent harmonic shear amplitudes

Updated 19 October 2026. Branch `feat/shear-amplitude-harmonic` starts at
`1154f86`, the stable main merge of dataset PR #6. The working tree was clean.
Experimental inputs, JSON records and previous implementation are preserved.
No existing independent two-component amplitude API was found. Existing private
RSS range and BP major-semiaxis expressions remain untouched.

Current objective is a standalone MCC/MCE library, analytical and reference-code
tests, and one short technical note. No production integration or configuration
change is authorised in this phase.

Reference inspection found active MCC and commented MCE expressions in the two
algebraic MATLAB Ca.m files. MOI/Ca.m uses path-length quadrature and is a different
definition. PragTic MCCM_IM/MCEM_IM process sampled polygons iteratively and are
not the harmonic solver to port. The full PDF is absent; supplied equation
transcriptions and MATLAB sources are available. Do not claim a full PDF review
or native MATLAB execution.

Baseline `1154f86` passed all 58 existing tests before implementation. The new
`fatigue/models/utils/shear_amplitude.py` supplies vectorised `shear_amplitude`
with explicit mcc_harmonic/mce_harmonic identifiers and the rss alias. It uses a
scaled closed-form Gram eigenvalue for MCC and the scaled Frobenius norm for MCE.
It rejects incompatible shapes, non-real/non-finite inputs and unrepresentable
outputs; exact zero is the only zero threshold. No production import was added.

The 12 tests in `tests/test_shear_amplitude.py` passed. They cover analytical
values, degeneracy, rotations/reflections and phase shifts, mean translation,
time sampling and LCM, batched inputs, SVD cross-checks, float64 extreme scales,
input rejection, existing RSS/engineering-strain conventions, BP equivalence,
and two MATLAB-equation load reconstructions. Ruff passed for both new files.
Tests first failed on the absent module, then passed after implementation.

The tested numerical milestone is published as `88229e93e2a2fb3527f23164537e2ab0c513351f`.
The technical note is [Independent harmonic shear amplitudes](shear-amplitude.md).
The complete suite passed all 70 tests (58 existing plus 12 new), without skips.
Ruff, documentation links, generated method catalogue and diff checks passed.
A final input-validation check also rejects nonzero extended-precision values
that would silently become zero when converted to float64.

A protected-path diff against `1154f86` confirmed unchanged production helpers,
numbered scripts, methods, criteria, requirements, calibrations, Excel inputs and
reference results. Only the new module, tests, note and this checkpoint belong to
this phase. The new module is not imported by any production file. Thus the
existing numerical implementation and its calibration fingerprints are preserved.
No new independent experimental fatigue validation is claimed.

The algebraic MATLAB files have identical SHA-256 `e96f70827fc0ce5c158ab56ade21af3a9eb8e4b72e47f205c5f3c93ab08092e5`.
Their active expression and commented alternative were reconstructed in two test
cases. Neither native MATLAB/Octave nor the PragTic application was executed.
The appendix MAT comparison array was inspected for scope, not used as a pure
amplitude benchmark. A 100,000-path timing comparison measured about 0.020 s for
scaled MCC versus 0.062 s for NumPy SVD; it is a local observation only.

Phase 1 is complete in [PR #7](https://github.com/andreachiocca/multiaxial-fatigue/pull/7),
open without merge. Verified implementation commit is
`2a17791d06e951a557c3d8707c2d156dbbb49f09`; a later commit only records this status
and clarifies a test comment. Both GitHub workflows passed, push run 39242032913
and PR run 39242046649. The published tree matches the tested local tree exactly.
The PR contains only four files: the amplitude module, its tests, the technical
note and this existing checkpoint. The previous experimental-data work remains
preserved on main and is not part of the PR diff.

Next action is human review of PR #7. Keep it open without merge. Do not integrate
these algorithms into a criterion until that separate phase is authorised. Phase 2 must separately review method-specific
integration, effective-algorithm metadata and dependency-scoped compatibility.
The following dataset checkpoint is a retained record of completed earlier work.

## Previous dataset workflow checkpoint

Updated 19 October 2026. Branch `feat/dataset-eligibility`, based on `022b6c0`.
PR #5 was merged into main before this session. The recovered JSON/run-out
implementation is preserved. No uncommitted changes existed at session start.

## Current objective and actual status

The revised dataset importer, method preflight, safe batch skips and observation
subset reporting are implemented. Draft PR [#6](https://github.com/andreachiocca/multiaxial-fatigue/pull/6) is open without merge. JSON history, formulations and family labels
are preserved. The remote branch exists; its first importer milestone
is `dfb514da82bd2792dd1d4da86fc0bb9dc9891e79`. Direct Git push has no credentials;
publication uses the GitHub connector. Local milestone history is retained.

- The tested local migration used 12 revised workbooks plus unchanged Al7075_T6, absent from
  the revised collection. The user explicitly approved publication on 19 October 2026; the revised collection is now installed locally for publication on the same PR. Four previous active files remain byte-identical in
  `Experimental_Data/archive/pre-dataset-migration` for reproducible regressions.
- Header aliases are explicit. Stress columns remain authoritative; force and
  moment columns are preserved metadata, never silently substituted. Phases are
  relative to the axial phase. Missing elastic properties leave strain undefined.
- Per-method input, conditional approximation, calibration and case checks are in
  `fatigue/eval/eligibility.py`. Script 0 skips each failed pair independently;
  Script 1 issues a controlled warning; Script 3 exports the matrix and can run
  compatible combinations. No automatic recalibration occurs during assessment.
- `fatigue/report/subsets.py` and Script 2 intersect identified failure
  observations separately for life and damage metrics. Individual results remain
  intact. Dataset fingerprints prevent mixing changed inputs; ambiguous IDs and
  duplicate results cannot enter a common comparison.
- New tests cover missing inputs/results, run-out-only data, conditional method
  requirements, batch continuation, and pairwise/group matched observations.

## Verification executed

`MPLBACKEND=Agg python tools/run_tests.py` passed 58 tests. Ruff, documentation
links, method catalogue and `git diff --check` passed. The original numerical
fixtures run against the archived original active workbooks, not substituted
new data. The attached previous collection was compared with the new collection
using `tools/audit_dataset_migration.py`; the detailed structured audit was
generated locally and can be regenerated from the supplied source files. All 56 shared 42CrMo4_QT and 44 AISI316L
observations retain the supplied stresses. AISI316L properties changed.

Script 0 with n1=3, n2=6, grid=7, refine=2 explicitly generated/reused all 15
method calibrations on revised 42CrMo4_QT under `/tmp/wg43-new-calibrations`.
Script 3 assessed all 15 base methods plus FS_ext and skipped 48 unsupported
pairs in the selected four-material batch. BP/SWT common life comparison used
45 failure observations per method. Eleven run-outs per method remain in
individual outputs with undefined conventional errors. These are software and
regression checks, not independent scientific validation or final calibration.
Temporary example records are not installed as default production calibrations.

The exact published implementation `904c61994e1bb10e5808ac3dd73f99ac794b53b7`
was checked in `/tmp/wg43-public-check`. The suite passed with 58 discovered
checks, 55 executed successfully and 3 full-collection tests skipped because the
withheld workbooks were absent. Ruff, documentation links and catalogue passed.
GitHub Actions succeeded for both push run 39235075215 and PR run 39235140578.
The archived-input end-to-end runner completed 60 combinations and 4,290 rows
(3,646 failures, 644 run-outs); four GSA/LIU2021 combinations on Al/Ti were
explicitly skipped because the steel empirical estimate is unsupported there.
This runner checks execution, counts and recomputed errors. The frozen numerical
fixtures separately check unchanged coefficients and method results; a new exact
full-table comparison against main was not performed this session.
The common parity plot was inspected visually. The installer dry run, byte
preservation, legacy snapshot retirement and conflicting archive checks passed.

## Scientific limitations and precise next action

Do not silently correct measured source data. CX18 declares 30NCD16 in Summary;
Ti6Al4V includes a Summary name Ti64Al4V. Both identity conflicts are flagged.
AISI316L averages fatigue strengths at different R, so methods requiring its
reversed axial fatigue strength are blocked until that input is resolved.
En3B/GJS600 have only notched observations in the revised files and remain in
storage outside the current plain-specimen scope. Missing nu blocks strain
methods; run-out-only observations do not define an S–N curve. Steel empirical
estimates are blocked for non-steels unless the required fatigue properties are
supplied. Published historical results are unchanged.

The previous raw-workbook publication block was resolved by the user's explicit
"si autorizzo" on 19 October 2026, answering the request to publish all twelve
workbooks, raw observations and embedded content in the public repository.
The same PR #6 is used without merge. The supplied ZIP was installed with
`tools/install_revised_dataset.py --apply`; the detailed migration audit and
eligibility snapshots are regenerated. Original workbook bytes are preserved.
The complete pre-publication recovery history remains at local branch
`recovery/revised-dataset-local-20261019` (`4111686`).

Automatic review accepted workbook publication under this approval but rejected
the additional detailed audit of source-property data. That optional JSON export
is excluded from the publication; its generator remains available. Do not publish
that extra report without separate approval. This does not block publishing the
explicitly authorised twelve workbooks.

The fresh 58-test run passed without skips after installation of the authorised
collection. Ruff, documentation links, the method catalogue and diff checks
passed. The eligibility snapshot was regenerated. Commit
`52b061b87d351be310d4e10127c085a5c985b24b` publishes the twelve workbooks on PR #6.
Every published workbook was compared byte for byte with `New_Exp_data/` in the
supplied ZIP. The public migration excludes only the additional detailed
source-property audit. Historical workbooks, JSON records and reference results
remain preserved. GitHub CI passed for the data commit on both the push
(run 39237587054) and PR (run 39237595131) workflows, with the complete revised
collection available.

Current publication is complete and verified. PR #6 remains open without merge.
Next action is human review of PR #6; no further dataset publication step is pending. Future scientific
work should address the source identity and mixed-fatigue-strength restrictions
above before interpreting affected method/material comparisons. No repeated
implementation or further permission for the twelve workbooks is required.
Resume from here, not the historical audit.

## Historical scientific audit, 7 October 2026

## Completed

- Resumed the existing working tree and preserved the previously published Böhme–Papuga and CAIM implementations and their tests. BP and CAIM are distinct formulations and remain separately selectable.
- Reviewed all 20 registered base methods. Corrected engineering shear factors, Findley plane selection, Li mean shear, Liu's removable zero-shear singularity, harmonic strain decomposition, zero-shear cancellation, and critical-plane ties. See [method-audit.md](method-audit.md) for equation-by-equation scope and limitations.
- Disabled CARSPA, DANGVAN, OTT, SWTD and ZHU_EDP, including direct evaluation guards, with reasons and restoration requirements. Removed obsolete ZHU_EDP calibration.
- Added source metadata and a [bibliography](bibliography.md) for every method. CSVs identify the method reference and implementation revision; comparisons exclude disabled and pre-audit results by default.
- Added `Error_log10_dp = log10(DP_e/DP)` to export, statistics and plots for all slopes, including exactly flat curves. Fitting no longer inverts nearly flat curves. Invalid fits/calibration objectives fail explicitly.
- Preserved experimental workbooks and historical result CSVs. Corrected results must be regenerated; fitted method constants need recalibration before quantitative model rankings.

## Verification completed

- `MPLBACKEND=Agg python -m unittest discover -s tests -q`: **20 tests passed**, including the existing BP/CAIM reference cases and new numerical/reporting regressions.
- `MPLBACKEND=Agg python tests/run_material_validation.py /tmp/fatigue-audit`: **64 runs / 3,936 rows passed**, covering 15 retained base models plus FS_ext on all four supplied materials. Checks include finite damage/log ratios, source metadata, ratio recomputation and suppressed life predictions in DP_DIFF mode.
- `MPLBACKEND=Agg python 2_compare_results.py --results_dir /tmp/fatigue-audit --out_dir /tmp/fatigue-audit-plots`: passed, **21 PNG plots**, plus TikZ and statistics. All four materials have log10 damage statistics for 16 variants; the three with usable negative slopes also have life statistics. AISI316L has 496 damage-space rows and no life predictions. Inspected its damage boxplot visually.
- An exactly-flat synthetic curve passed CSV export, plotting and statistics; legacy CSVs missing DP_fit passed the explicit opt-in path without crashing.
- Five retained parameter calibrations passed a coarse-grid smoke check on 42CrMo4_QT: FS, FIN, LI, MGSE_ZHU, MSWT. This is an execution check, not final recalibration of workbook constants.
- Python parsing/compilation and `git diff --check` passed. Final diff reviewed; no raw data or archived Results changes.
- Non-failing warnings: empirical steel-property fallbacks where material coefficients are missing; existing Matplotlib deprecation/marker warnings.

## Unresolved evidence and exact next action

Full original-equation/plane-rule verification is still incomplete for **LIU2021, MGSE_ZHU, MKBM and LI**. MSWT is associated with the reported scalar modified-SWT form and is not claimed to implement Jiang's complete incremental model. These limitations are explicit in the audit; passing implementation tests is not independent scientific validation.

**Next action:** obtain an accessible full original Liu et al. article (DOI `10.1016/j.ijpvp.2021.104532`), inspect Eq. (27) and the definitions of its normal strain and shear ranges, and compare them with `fatigue/models/cp_methods/liu2021.py`. Add an independently calculated published load case, then correct or disable the model if necessary. Continue in order with Zhu 2018 (MGSE_ZHU/MSWT), Li 2011 (MKBM), and Li 2021 (LI), using the DOI mapping already saved in the bibliography. Do not repeat completed BP/CAIM integration or the general audit. Then recalibrate affected parameters and regenerate comparisons before interpreting rankings.

No pending test failure is known. The unresolved item is source completeness, not a failed test. The shared elastic, harmonic and finite-plane-grid limitations are documented in the audit.

## Files changed in this audit

- `0_calibrate_material_params.py`
- `1_run_material.py`
- `2_compare_results.py`
- `README.md`
- `docs/adding-methods.md`
- `docs/architecture.md`
- `docs/audit-checkpoint.md`
- `docs/bibliography.md`
- `docs/method-audit.md`
- `fatigue/eval/evaluator.py`
- `fatigue/eval/fitting.py`
- `fatigue/eval/metrics.py`
- `fatigue/models/cp_methods/car_spa.py`
- `fatigue/models/cp_methods/dangvan.py`
- `fatigue/models/cp_methods/fin.py`
- `fatigue/models/cp_methods/fs.py`
- `fatigue/models/cp_methods/gsa.py`
- `fatigue/models/cp_methods/gse.py`
- `fatigue/models/cp_methods/li.py`
- `fatigue/models/cp_methods/liu2021.py`
- `fatigue/models/cp_methods/mgse_yu.py`
- `fatigue/models/cp_methods/mgse_zhu.py`
- `fatigue/models/cp_methods/mkbm.py`
- `fatigue/models/cp_methods/mswt.py`
- `fatigue/models/energy_based_methods/swtd.py`
- `fatigue/models/energy_based_methods/zhu_edp.py`
- `fatigue/models/invariant_methods/ottosen.py`
- `fatigue/models/references.py`
- `fatigue/models/registry.py`
- `fatigue/models/utils/plane_history.py`
- `fatigue/report/export_csv.py`
- `tests/run_material_validation.py`
- `tests/test_method_audit.py`


