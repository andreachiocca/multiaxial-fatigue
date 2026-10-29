# Development continuation checkpoint

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

