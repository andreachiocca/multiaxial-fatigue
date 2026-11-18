# Contributing to multifat

`multifat` is an open-source platform for multiaxial fatigue assessment. It began as the comparison tool behind the **FABER WG4.3 scientific review article** and stays public afterwards for evaluating and comparing published multiaxial fatigue criteria. The present phase covers **constant-amplitude loading of plain/smooth geometries**, P/NP where the selected law supports it, without explicit notches, cracks or significant gradients; local finite-element and time-history input bridges extend this scope under active development.

## Check external code before starting from scratch

Consult the [external-code inventory](docs/EXTERNAL_IMPLEMENTATIONS.md) and the method’s discovery metadata. An external implementation stays external until independently reviewed, integrated and tested. It does not change WG4.3 completeness or validation.

For any reuse or translation, record the exact file/revision, method/version correspondence, original paper, loading/material/geometry assumptions, dependencies, tests and license grant. Public code without an explicit license is not permission to copy, translate or redistribute. This repository uses the MIT License ([LICENSE](LICENSE)), and contributions are accepted on the same terms. Permissively licensed code (MIT, BSD, ISC, Apache-2.0) is compatible if its notices are kept. GPL-licensed and unlicensed code is not. Without permission, link the source and ask for it, or derive an independent implementation from the published formulation and review it. Do not relicense code of unknown origin.

Compare intermediate quantities and source-specific plane/amplitude conventions, not only final life. External tests and published claims are evidence leads, not validation of a port. Attach provenance and notices to the PR. For geometry/gradient/spectral software, first establish whether a faithful smooth-point CA subcase exists. Do not import a broader physical model silently.

Start with [METHODS.md](METHODS.md), the [suggested roadmap](ROADMAP.md) and [prepared Issues](docs/contribution-issues/README.md). Contributions can implement a method, resolve a source ambiguity, build a reproducible benchmark, validate existing code or improve documentation. Priority is advisory, not an assignment.

## Ten steps from a method to a PR

1. **Choose an ID.** Browse [all missing methods by family](docs/METHOD_INDEX.md#missing-methods-by-family), repairs or source-review tasks. The roadmap is a suggested sequence, though another catalogued method is also acceptable. A genuinely new candidate needs its relationship to existing entries established before it gets a permanent ID.
2. **Find or propose the task.** Use the entry’s Search Issues link, check aliases and closed Issues, then claim a task or open a scoped proposal. The Propose task link opens the method-specific template with its ID and title. Prepared briefs are not live Issues.
3. **Confirm the primary formulation.** Obtain the exact original article, relevant errata or supplements, and a traceable citation. DOI metadata confirms identity, not equations. An unresolved variant needs reference work before implementation. Do not upload copyrighted papers without redistribution rights.
4. **Write a bounded specification.** State CA P/NP/mean-stress scope, material domain, equations, amplitude, range and shear conventions, required measurements and the output (limit, damage parameter, native life or probability). Record repository adaptations separately from the original law.
5. **Inspect reusable components.** Read the [model interface guide](docs/adding-methods.md), one nearby implementation and the relevant projection/search helpers. Reuse must preserve the source plane rule and the physical simultaneous path. List any input, native-life or spatial/statistical prerequisite.
6. **Implement and integrate.** Place the law in its model family, declare parameters, register it and add runtime citation metadata. Follow the paths and interfaces below, and propose major prerequisites separately.
7. **Document inputs and limits.** Add material fields only when needed. State units, provenance, measured versus estimated constants, supported loading and explicitly unsupported cases. Never omit work terms silently or replace a native life law with a generic fit.
8. **Numerically verify.** Add focused analytical, degenerate, phase, mean-stress and frame checks, unsupported-input checks and search or integration convergence. Compute at least one reference result independently of the implementation.
9. **Validate and report evidence.** Reproduce a published case with its material and loading data and, where available, its intermediate values. Keep calibration separate from validation, and state inaccessible data or remaining predictive gaps. Run the regression and documentation checks below.
10. **Update and submit.** Update the JSON catalogue, source, tests, evidence and actual Issue links, then regenerate the documentation. Submit a PR from the template, citing the permanent ID and acceptance checklist. A maintainer reviews completeness and validation independently; passing unit tests alone does not promote either.

## Develop locally

Use Python 3.10+ and `uv sync` as in the README (or `pip install -e .`). Branch from current `main` and run entry scripts from the repository root. Preserve the experimental workbooks and `Results/archive/`. Fresh outputs go to `Results/current/` (ignored by Git), test output to a temporary directory. Do not commit generated caches, private data or article PDFs.

Follow [the architecture](docs/architecture.md), [adding methods](docs/adding-methods.md) and the
[notation table](docs/notation.md):

- Plane-level criteria belong in `src/multifat/models/cp_methods/`, implement `CriticalPlaneMethod.evaluate_on_plane` and return `CPPlaneResult(damage, metric)`. A scan may be reused only when it matches the source's plane rule.
- Direct criteria belong in `invariant_methods/`, `energy_based_methods/` or `integral_methods/`. Implement `FatigueModel.evaluate_case` and return `CaseResult`. Consume the physical `harmonics` when timing or phase matters.
- Register in `src/multifat/models/registry.py` and add runtime citation metadata in `src/multifat/models/references.py`. Declare every constant the model reads: required ones in `required_params()`, those read with `params.get` in `optional_params`. Do not add an optional model constant to `src/multifat/imp/schema.py` just to read it from the Summary (changing the schema alters every dataset fingerprint): a Summary column named like a declared `optional_params` key is already forwarded to the model. Follow existing names and `_ext` conventions.
- A new physical model or data structure needs its own proposal. The plane/scalar interface does not justify dropping source requirements. Keep BP/CAIM distinct and preserve all-slope `log10(DP_e/DP)` reporting and the life-suppression rules.

Declare a CP method's out-of-phase shear rule with `np_shear_rule` (`"rss"` legacy default or `"longest_chord"`). BP/CAIM use the major-axis/MCC convention. Strain shear is tensorial in tensors and engineering shear in the relevant scalar ranges. Record amplitude versus range, cycles versus reversals, units, stress sign, plane normal, coordinate basis and every empirical material fallback.

Ruff checks syntax errors, unused imports, redefinitions and undefined names. Black owns formatting (`uv run black .`). The one-off formatting commit is in `.git-blame-ignore-revs` (use `git config blame.ignoreRevsFile .git-blame-ignore-revs`). The documentation checker covers maintained local links and fragments; external availability and frozen research prose need separate review. Runtime and `dev` dependencies live in `pyproject.toml`, pinned by `uv.lock`.

## Evidence required in a PR

Provide the source and version, supported scope, equations implemented, parameter provenance, numerical choices, limitations and a reproducible validation command. Cite the original equation numbers without copying excessive text. Include focused checks for pure axial/torsion, degeneracy, mean stress, phase, and coordinate or time-origin invariance, and show search or quadrature convergence where it matters.

Reference cases should include the material data, the actual repeated path, intermediate quantities and the expected output, with at least one result calculated independently. Keep calibration and validation data separate, and record measurement uncertainty and inaccessible data honestly. A passing test of the coded equation is numerical verification, not predictive accuracy.

Run:

```bash
uv sync
uv run ruff check .
uv run black --check .
MPLBACKEND=Agg uv run pytest
uv run python tools/generate_method_docs.py
uv run python tools/generate_method_docs.py --check
uv run python tools/generate_workbook_docs.py --check
uv run python tools/check_docs.py
git diff --check
```

Include commands, results and known warnings in the PR. Routine documentation changes do not require regenerating experimental results.

### Changes that can alter numerical outputs

Shared code (plane search, calibration, import, export, a model's equations) must be checked against the previous outputs before merge:

```bash
git stash   # or check out the base commit in a worktree
uv run python tools/regression_snapshot.py run /tmp/snap-base --jobs 4
git stash pop
uv run python tools/regression_snapshot.py run /tmp/snap-new --jobs 4
uv run python tools/regression_snapshot.py compare /tmp/snap-base /tmp/snap-new
```

The snapshot runs every registered model and `_ext` variant on every fatigue-life workbook, the calibration and the comparison, then `multifat-pipeline --models all` on every workbook (calibrated plane constants, fatigue-limit workbooks) and one cross-batch transfer, and hashes every output, leaving out the volatile fields of the pipeline sidecars and run reports (about 40 minutes on a laptop; `--skip-pipeline` halves it). A refactor must be output-identical. An intended numerical change goes in its own commit, with the changed outputs and the reason in the PR (which models, which materials, by how much). Bump `IMPLEMENTATION_REVISION` in `src/multifat/models/references.py` (using a new dated string such as `"YYYY-MM-DD"`) whenever a change can alter numerical outputs (model equations, plane search, calibration, data import or overlays, master-curve fitting, or load-error/endurance evaluation); doing so stamps new result CSVs, flags older results as legacy in `multifat-compare`, rejects stale pipeline calibration sidecars, and triggers staleness checks in `examples/` manifests. The end-to-end sweep `uv run python tests/run_material_validation.py /tmp/wg43-validation` checks finite outputs and row counts on all materials.

### Workbook contributions

Workbooks follow the [workbook format](docs/workbook-format.md). Record the source of every value (paper, table or figure, page, tabulated or digitised) in [Experimental_Data/PROVENANCE.md](Experimental_Data/PROVENANCE.md). Supplied workbooks are never edited in place: keep their original bytes intact (recorded in [Experimental_Data/MANIFEST.json](Experimental_Data/MANIFEST.json)) and declare corrections as cell-level operations in `src/multifat/imp/overlay_files/<stem>.overlay.json`, matched by the workbook's SHA-256 and applied in memory on load (`--no-overlay` returns the supplied data). When adding or replacing a supplied workbook, update [Experimental_Data/MANIFEST.json](Experimental_Data/MANIFEST.json), and treat any workbook or overlay change as numerical (snapshot above).

## Update the catalogue

Experimental inputs and fitted coefficients use separate sources. Preserve the Excel observation-status convention and keep run-outs out of coefficient optimisation, master curves and ordinary error statistics (they are exported as `Flag = 2` rows). Calibrated constants belong in the calibration sidecar, never the workbook. Changing the calibration data flow needs a failure-only versus failure-plus-run-out regression; see the [calibration workflow](docs/calibration-workflow.md) for the records in `Calibrations/`.

`docs/catalogue/methods.json` is the source of truth. Keep permanent IDs, add canonical relationships and optional secondary classifications, and count each identifiable scoped method once by its primary family. Exclude unresolved bundles and overlaps from progress totals until reconciled. The archived research CSV is provenance, not a live status source.

Update implementation status, declared scope, module and symbol, tests, remaining work, validation evidence and the actual Issue URL, then generate `METHODS.md`, `ROADMAP.md`, the CSV, the README counts and the prepared briefs with `tools/generate_method_docs.py`. CI checks output freshness and evidence paths against the registry. Do not edit generated files independently. Roadmap selection lives in each JSON record as `roadmap_rank`, `development_track` and `development_reason`, and every record needs an explicit disposition. Unscheduled candidates stay in the catalogue. The generator’s Issue-ID list controls which prepared briefs are emitted, not which methods contributors may choose.

### Ready for “fully implemented”

- [ ] Exact supported formulation/version identified, with explicit adaptations and loading/material restrictions.
- [ ] Every essential component of that declared formulation implemented; no known essential omission.
- [ ] Required inputs/parameters supported and checked; measured and estimated properties distinguished.
- [ ] Original plane, path, integration and native life rules faithfully handled within declared scope.
- [ ] Registered, documented and covered by meaningful automated tests, including unsupported-input behavior.
- [ ] Source → catalogue → code → tests → evidence → Issue/PR links complete.
- [ ] Existing numerical regression checks pass and generated docs are current.
- [ ] Completeness reviewed by a maintainer; uncertain completeness remains `requires_review`.

Validation is a separate promotion: `not_tested` → `basic_numerical_tests` → `reference_examples_reproduced` → `independently_validated`. Reference reproduction needs independently calculated published examples. Independent validation needs traceable published or experimental data beyond self-consistency tests. A complete scoped implementation may still rest on basic numerical evidence.

## Submit

Open a PR linking the permanent method ID and the relevant Issue. Summarize the scientific capability, scope, departures and evidence. Split large workstreams into separate reviewable prerequisites. Do not re-enable disabled proxies until their recorded audit gaps are resolved. By submitting, you agree to license your contribution under the MIT License.
