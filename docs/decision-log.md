# Decision Log

This log records architectural, numerical, and data decisions made during the development of
`multifat` (the carry-forward of `andreachiocca/multiaxial-fatigue`). Each decision is assigned a stable
reference identifier (`d0NN`) cited across commit messages, code docstrings, data overlays, and audit
ledgers.

---

## Decision Index

| Ref | Date | Status | Topic | Summary |
|---|---|---|---|---|
| [`d001`](#d001) | 2026-10-15 | **Active** | Fork baseline | Established `dev/tg-carry-forward` based on upstream `e70a8af`. |
| [`d011`](#d011) | 2026-10-15 | **Active** | IP boundary | Third-party GPL code (CPOpt, FS-Sol, etc.) may not be ported directly; implementations must be independently derived from papers. |
| [`d036`](#d036) | 2026-10-18 | **Active** | Provenance N1 | Corrected 42CrMo4 `Plain_multiaxial_9` torsional phase to $0^\circ$ (proportional PP loading). |
| [`d037`](#d037) | 2026-10-18 | **Active** | MSWT weighting | Converted parameter $a_{\mathrm{MSWT}}$ via $a' = a / (2 - a)$ to match Jiang (2000) / Zhu (2018) printed weighting without altering predicted lives. |
| [`d038`](#d038) | 2026-10-18 | **Active** | Provenance N6/N8 | Corrected AISI 316L non-proportional phase to $90^\circ$ and appended 8 omitted Janssens Table 4 test specimens. |
| [`d040`](#d040) | 2026-10-18 | **Active** | Reference outputs | Archived historical reference CSVs into `Results/archive/pre-audit/` rather than committing 170+ regenerated files. |
| [`d041`](#d041) | 2026-10-18 | *Reversed* | Maintenance logs | Originally deleted `docs/audit-checkpoint.md`; reversed when upstream maintainer requested tracking it. |
| [`d042`](#d042) | 2026-10-18 | *Superseded* | MMMH mean stress | Kept MMMH unregistered pending Mertens references; superseded when endurance limit formulation was verified and registered. |
| [`d043`](#d043) | 2026-10-18 | **Active** | Banvillet G26 | Implemented closed-form tension-torsion index; kept unregistered pending rotating-bending limits and specimen geometry in workbooks. |
| [`d044`](#d044) | 2026-10-18 | *Superseded* | Project licence | Held licence choice open between MIT, BSD, Apache, and GPL-3.0; superseded by `d071`. |
| [`d046`](#d046) | 2026-10-18 | **Active** | Output paths | Configured `multifat-run` default results directory to `Results/current/`. |
| [`d047`](#d047) | 2026-10-18 | **Active** | Point kernels | Retained uncalibrated models (G24, G26, G29, G30, N34) as modular functions with missing inputs recorded in the catalogue. |
| [`d050`](#d050) | 2026-10-18 | **Active** | Release readiness | Closed Phase 4 public-release preparation milestone. |
| [`d052`](#d052) | 2026-10-19 | **Active** | Pipeline design | Retained unified pipeline architecture (evaluating missing constants in one pass) rather than upstream stateful multi-script workflow. |
| [`d054`](#d054) | 2026-10-19 | **Active** | Dataset merge | Integrated upstream revised collection (`71eaac3`, 12 numbered workbooks) while retaining unnumbered `Al7075_T6.xlsx`. |
| [`d055`](#d055) | 2026-10-19 | **Active** | Batch workbooks | Created `14_42CrMo4_Reis.xlsx` and `15_Ti6Al4V_Berto.xlsx` to isolate distinct metallurgical batches and printed constants. |
| [`d056`](#d056) | 2026-10-19 | **Active** | AISI 316L data | Restored bulk $E = 193\,\mathrm{GPa}$, $\nu = 0.30$, cleared XRD lattice constants, fixed $90^\circ$ phase, and appended Table 4 tests. |
| [`d057`](#d057) | 2026-10-19 | *Refined* | Fatigue index error | Defined fatigue index error for run-out-only datasets ($FI = DP/DP_{\mathrm{ref}} - 1$, Böhme & Papuga 2024); refined by `d063`. |
| [`d058`](#d058) | 2026-10-19 | **Active** | Calibration records | Supported upstream `Calibrations/` JSON records via `--constants FILE`. |
| [`d059`](#d059) | 2026-10-19 | **Active** | Steel Poisson ratio | Set $\nu = 0.30$ (Dubbel 2004) for steels missing Poisson's ratio in Summary tables (`2_CX18`, `3_25CrMo4`, `4_30NCD16`, `7_982FA`). |
| [`d061`](#d061) | 2026-10-19 | **Active** | Braun mean stress | Restored $\sigma_m = \sigma_a$ on Braun $R = 0$ plain rows in `8_AISI316L.xlsx` to prevent exclusion under revised-format loader. |
| [`d062`](#d062) | 2026-10-19 | **Active** | Compute platform | Shifted numerical snapshot benchmarks to WorkstationUniPi due to i9-14900K hardware instability on WorkstationUniTN. |
| [`d063`](#d063) | 2026-10-19 | **Active** | Load-scaled index | Implemented load-scaled index $FI_{\mathrm{load}} = 1/s - 1$ where $DP(s \cdot \text{case}) = DP_{\mathrm{ref}}$ for quadratic/energy criteria (SWT). |
| [`d065`](#d065) | 2026-10-19 | **Active** | Branch management | Fast-forwarded main development line to `0e87e95` and pruned transient synchronization branches. |
| [`d066`](#d066) | 2026-10-19 | **Active** | Pipeline edge cases | Resolved deferred edge cases (Excel lock files `~$*`, duplicate material IDs, non-finite values). |
| [`d067`](#d067) | 2026-10-20 | **Active** | Upstream PRs #8-#9 | Integrated harmonic Findley evaluation (`fin_harmonic.py`) and configurable calibration options. |
| [`d068`](#d068) | 2026-10-20 | **Active** | Findley amplitude | Added `--shear-amplitude-algorithm` to `multifat-calibrate`; default remains legacy production MCE. |
| [`d069`](#d069) | 2026-10-20 | **Active** | MCC plane rule | Routed Findley MCC calibration through the refined plane search (`FindleyMCC`, `longest_chord` rule). |
| [`d070`](#d070) | 2026-10-20 | **Active** | Calibration guard | Enforced rejection of non-MCE calibration records when loading through `--constants` in production pipeline. |
| [`d071`](#d071) | 2026-10-22 | **Active** | MIT licensing | Formally relicensed `multifat` under the MIT License, superseding `d044`. |
| [`d072`](#d072) | 2026-11-05 | **Active** | FE interface scope | Issue 19 delivers the local-history contract, a file format, one importer and the range-strategy interface; further range rules and criteria for general histories belong to issue 20. |
| [`d073`](#d073) | 2026-11-05 | **Active** | FE file format | HDF5 is the universal FE history file; `h5py` is the optional extra `fe`. |
| [`d074`](#d074) | 2026-11-05 | **Active** | First FE adapter | VTK XML (`.vtu`, `.pvd`) is the first third-party importer, built on the standard library and numpy; ANSYS and Abaqus adapters are follow-up issues. |
| [`d075`](#d075) | 2026-11-05 | **Active** | FE conventions | MPa, mm, tensor shear, total mechanical strain, global Cartesian axes, six components in a fixed order; nothing is inferred by an adapter. |
| [`d076`](#d076) | 2026-11-05 | **Active** | FE eligibility | A method declares `requires_nonlocal`; the bridge refuses it on a local-only history; no registered method needs more than local data. |
| [`d077`](#d077) | 2026-11-05 | **Active** | FE test oracle | The synthetic thick tube is written from the textbook Lame and Saint-Venant equations, not from GPL code, and serves as analytical oracle. |
| [`d078`](#d078) | 2026-11-05 | **Active** | FE two-state bridge | Selecting a pair of load cases reduces a history to the two-step model inputs; this is exact only for two-state paths and is documented as lossy otherwise. |
| [`d079`](#d079) | 2026-11-05 | **Active** | History bridge scope | Issue 20 delivers the time-history contract, scalar rainflow, longest-chord events, damage per event through the existing models and Miner accumulation into a life; multiaxial counting and equivalent-cycle extractors are specified follow-ups. |
| [`d080`](#d080) | 2026-11-05 | **Active** | Proportional domain | The basic implementation is exact for proportional histories; a non-proportional history is refused unless an explicit, recorded approximation flag is set. |
| [`d081`](#d081) | 2026-11-05 | **Active** | Caller-supplied curve | The Basquin curve of the damage parameter is supplied by the caller (`SNCurve`, from the pipeline's in-memory `A_surv` and `b`); no record or CSV changes. |
| [`d082`](#d082) | 2026-11-05 | **Active** | Counting channel | The default scalar channel is the dominant singular direction of the mean-removed stress path, affine in the load factor; signed von Mises is explicit only. |
| [`d083`](#d083) | 2026-11-05 | **Active** | Block semantics | The history is one block repeating until failure; residue policy `rearrange` is the default, `half` (ASTM E1049) is selectable with a note; the choice is recorded. |
| [`d084`](#d084) | 2026-11-05 | **Active** | Native-life refusal | Methods with their own life law (`native_life`) are refused for accumulation; loading-scope errors are reported as ineligibility, never skipped. |
| [`d085`](#d085) | 2026-11-05 | **Active** | History package | New package `multifat.history`, imported by nothing existing, reusing `multifat.fe`; literature corpus in the Zotero collection `4J6FG4BS`, PDFs to be fetched later. |
| [`d086`](#d086) | 2026-11-06 | **Active** | FE CSV exchange format | A CSV file with a commented header is the recommended way to bring FE results in (the maintainer's proposal on issue 19, tightened): one row per location and step, named columns matched case-insensitively with Abaqus... |
| [`d087`](#d087) | 2026-11-06 | **Active** | Location identity and file version 1.1 | A location is the triple `(instance, node_id, material)`; keys `(instance, node_id, material, step)` are unique and the grid is complete. |
| [`d088`](#d088) | 2026-11-06 | **Active** | Strain composition and derivation | Total `E_*` is the total mechanical strain; `Eel_*` with `Epl_*` are stored and summed; `Eel_*` alone is the total strain for a linear-elastic file and refused for an elastic-plastic one; `Epl_*` alone needs... |
| [`d089`](#d089) | 2026-11-06 | **Active** | Step order and step role | Steps ascend; non-contiguous numbers are allowed and recorded; `time` must be strictly increasing along ascending steps and is never re-sorted. |
| [`d090`](#d090) | 2026-11-06 | **Active** | Node assessment | `assess_nodes` takes one `MaterialSetup(model, params, curve, E=None, nu=None)` or a mapping by material label, runs `assess_history` on each selected node and returns one row per node with a status: `ok`,... |
| [`d091`](#d091) | 2026-11-06 | **Active** | Header discipline and integrity | Required header keys, the value sets, `n_nodes` and `n_steps` checked against the file, undeformed coordinates required and checked across steps, `nodal_average_global` across several materials or instances refused... |
| [`d092`](#d092) | 2026-11-06 | **Active** | Writer helper and solver checklists | `FeCsvWriter` uses only the standard library so that a solver script can copy it; `write_tube_csv` writes the analytical tube in the canonical, Abaqus and Ansys column names as the test oracle. |
| [`d093`](#d093) | 2026-11-06 | **Active** | Solver adapters target CSV gate | `fe/solvers` writes atomic FE CSV through FeCsvWriter after canonical validation; no direct HDF5 output. |
| [`d094`](#d094) | 2026-11-06 | **Active** | Abaqus nodal reports | Read Field Output Reports with explicit nodal positions; collapse element-nodal duplicates only on exact agreement. |
| [`d095`](#d095) | 2026-11-06 | **Active** | Ansys text listings | Read complete PRNSOL tensor tables with global NLIST coordinates or complete named Workbench tables. |
| [`d096`](#d096) | 2026-11-06 | **Active** | Optional native RST adapter | Support solid and plane elements (SOLID185/186/187, PLANE182/183) with `ansys-mapdl-reader`. |
| [`d097`](#d097) | 2026-11-06 | **Active** | Solver evidence and untrusted files | Document real versus hand-written coverage; no third-party code execution. |
| [`d098`](#d098) | 2026-11-06 | **Active** | Whole-branch safeguards | Validate CSV row width before pandas projection; atomic output; require shared scalar load functions. |
| [`d099`](#d099) | 2026-11-06 | **Active** | In-house Ansys fatigue benchmark | Repository purged of external repository fixtures; in-house PLANE183 2-step fatigue cycle generated on WorkstationUniTN. |

---

## Detailed Records

### d001
* **Date:** 2026-10-15
* **Topic:** Carry-forward fork baseline
* **Status:** Active
* **Decision:** Establish `dev/tg-carry-forward` as an independent development branch starting from upstream commit `e70a8af`.
* **Rationale:** Preserves clean upstream provenance while restructuring into a modern `src/` layout with comprehensive method verification.

### d011
* **Date:** 2026-10-15
* **Topic:** Intellectual property and clean-room implementation
* **Status:** Active
* **Decision:** Code from GPL-3.0 repositories authored by collaborators (including CPOpt, FS-Sol, FI-Sol, SeAn-CP, and ECP) may not be ported directly into `multifat`.
* **Rationale:** Ensures clean IP provenance. All multiaxial fatigue criteria must be independently derived and verified directly from original scientific publications.

### d036
* **Date:** 2026-10-18
* **Topic:** 42CrMo4 proportional loading correction (Finding N1)
* **Status:** Active
* **Decision:** Set torsional phase to $0^\circ$ for row `Plain_multiaxial_9` in `42CrMo4_QT.xlsx`.
* **Rationale:** Anes et al. (IJF 58 Table 3, PP case) and Anes et al. (FFEMS 37 Table 2, case 1) confirm this test was proportional (in-phase). It had been erroneously tagged as $90^\circ$ out-of-phase.

### d037
* **Date:** 2026-10-18
* **Topic:** MSWT parameter weighting definition
* **Status:** Active
* **Decision:** Convert workbook parameter $a_{\mathrm{MSWT}}$ to $a' = a / (2 - a)$ to match the printed Jiang (2000) Eq. (6) and Zhu (2018) Eq. (13) weighting.
* **Rationale:** The original code implementation used twice the shear weight relative to the published formulation. The algebraic mapping $a' = a / (2 - a)$ harmonizes the definition with the published literature while leaving predicted lives exactly unchanged.

### d038
* **Date:** 2026-10-18
* **Topic:** AISI 316L audit corrections (Findings N6, N7, N8)
* **Status:** Active
* **Decision:** In `8_AISI316L.xlsx`, set torsional phase to $90^\circ$ on non-proportional rows (N6); set bulk elastic properties to $E = 193\,\mathrm{GPa}$ and $\nu = 0.30$ (N7); append the 8 missing test specimens from Janssens (2018) Table 4 (N8).
* **Rationale:** Janssens (2018) explicitly specifies $90^\circ$ phase shift for the non-proportional series. Literature values from Zhu et al. (2019) provide verified bulk elastic properties. The 8 omitted tests in Table 4 are valid fatigue specimens.

### d040
* **Date:** 2026-10-18
* **Topic:** Historical reference output preservation
* **Status:** Active
* **Decision:** Archive pre-audit CSV results into `Results/archive/pre-audit/` rather than committing 170+ regenerated CSV files to the root results folder.
* **Rationale:** Keeps git history lean, avoids massive diff noise across numerical commits, and preserves historical pre-audit outputs for audit traceability.

### d041
* **Date:** 2026-10-18
* **Topic:** Maintenance log retirement
* **Status:** Reversed
* **Decision:** Originally retired `docs/audit-checkpoint.md` as internal maintenance scratch. Reversed during upstream sync when Andrea Chiocca requested maintaining `docs/audit-checkpoint.md` as the unified progress record.

### d042
* **Date:** 2026-10-18
* **Topic:** MMMH criterion registration fallback
* **Status:** Superseded by Batch P4-2
* **Decision:** Originally deferred registering MMMH because Mertens & Liebich (2017) references [16, 23] for the Haigh diagram $M_N$ curve were not yet retrieved.
* **Outcome:** Dubbel C10 Section 10 verified the endurance limit relation $M = (\sigma_A(R=-1) - \sigma_A(R=0)) / \sigma_m(R=0)$, enabling full registration of `MMMH` in `src/multifat/models/registry.py`.

### d043
* **Date:** 2026-10-18
* **Topic:** Banvillet criterion (G26) scope
* **Status:** Active
* **Decision:** Implement Banvillet (2003) closed-form tension-torsion index with 15 verified tests, but keep the criterion unregistered in the general pipeline.
* **Rationale:** Full general evaluation requires rotating-bending fatigue limits and cross-section geometry (solid vs. thin-walled tube) not present in standard workbooks.

### d044
* **Date:** 2026-10-18
* **Topic:** Licence selection options
* **Status:** Superseded by `d071`
* **Decision:** Enumerated licence choices (MIT, BSD-3-Clause, Apache-2.0, GPL-3.0) and held decision open for author alignment. Formally superseded by `d071` selecting MIT.

### d046
* **Date:** 2026-10-18
* **Topic:** Default results output path
* **Status:** Active
* **Decision:** Set default output path for `multifat-run` and `multifat-pipeline` to `Results/current/`.

### d047
* **Date:** 2026-10-18
* **Topic:** Unregistered point kernel disposition
* **Status:** Active
* **Decision:** Retain unregistered point kernels (G24, G26, G29, G30, N34) as tested standalone modules, documenting their unresolved workbook inputs in `docs/catalogue/methods.json` rather than exposing partial criteria.

### d050
* **Date:** 2026-10-18
* **Topic:** Phase 4 milestone completion
* **Status:** Active
* **Decision:** Concluded public-release documentation restructuring and verification sweep on WorkstationUniTN.

### d052
* **Date:** 2026-10-19
* **Topic:** Pipeline architecture vs. upstream calibration store
* **Status:** Active
* **Decision:** Maintain `multifat`'s single-command pipeline architecture with inline missing-constant calibration and sidecar JSON records, declining upstream's stateful multi-script workflow (`fatigue/eval/calibration_store.py`).

### d054
* **Date:** 2026-10-19
* **Topic:** October 2026 dataset migration
* **Status:** Active
* **Decision:** Adopt upstream revised dataset collection (`71eaac3`, 12 numbered workbooks `1_`–`8_`, `10_`–`13_`), while retaining unnumbered `Al7075_T6.xlsx` which had no upstream replacement.

### d055
* **Date:** 2026-10-19
* **Topic:** Dedicated batch workbooks
* **Status:** Active
* **Decision:** Split secondary test batches into dedicated workbooks:
  - `14_42CrMo4_Reis.xlsx`: 29 Anes/Reis tests ($S_u = 1100\,\mathrm{MPa}$ batch, printed strain-life constants).
  - `15_Ti6Al4V_Berto.xlsx`: 48 Berto tests ($S_u = 978\,\mathrm{MPa}$, printed fatigue limits).
* **Rationale:** Prevents averaging static and fatigue properties across distinct heats and metallurgical conditions.

### d056
* **Date:** 2026-10-19
* **Topic:** AISI 316L data overlay
* **Status:** Active
* **Decision:** Declare corrections for `8_AISI316L.xlsx` in [`src/multifat/imp/overlay_files/8_AISI316L.overlay.json`](../src/multifat/imp/overlay_files/8_AISI316L.overlay.json):
  - Summary: Janssens $E = 193\,\mathrm{GPa}$, $\nu = 0.30$; Braun $E$ and $\nu$ set to `Null` (published values were XRD lattice constants for residual stress evaluation).
  - Dataset: 11 non-proportional cases set to $90^\circ$ phase shift; 8 omitted Table 4 cases appended.

### d057
* **Date:** 2026-10-19
* **Topic:** Fatigue index error for fatigue-limit datasets
* **Status:** Refined by `d063`
* **Decision:** Evaluate run-out-only datasets via the fatigue index error (Böhme & Papuga 2024, Eqs. 8–9):
  $$FI_{\mathrm{axial}} = \frac{DP(\text{case})}{DP(\text{axial } R=-1 \text{ at } \sigma_{-1})} - 1$$
* **Rationale:** Provides quantitative model error assessment when cycle counts are censored run-outs at the endurance limit.

### d058
* **Date:** 2026-10-19
* **Topic:** Loading external calibration records
* **Status:** Active
* **Decision:** Support upstream `Calibrations/` JSON records as an explicit pipeline input via `--constants <file>`, validating input compatibility before assessment.

### d059
* **Date:** 2026-10-19
* **Topic:** Steel Poisson's ratio fallback
* **Status:** Active
* **Decision:** Set Poisson's ratio $\nu = 0.30$ (Dubbel 2004, Part C) in Summary rows for steel workbooks lacking reported $\nu$ (`2_CX18`, `3_25CrMo4`, `4_30NCD16`, `7_982FA`).
* **Rationale:** Enables elastic strain tensor reconstruction for strain- and energy-based criteria without inventing material-specific constants.

### d061
* **Date:** 2026-10-19
* **Topic:** Braun mean stress preservation
* **Status:** Active
* **Decision:** Set mean axial stress $\sigma_m = \sigma_a$ for the 17 Braun $R = 0$ plain rows in `8_AISI316L.xlsx`.
* **Rationale:** The revised dataset loader interprets empty stress cells as unavailable loading components. Setting $\sigma_m = \sigma_a$ correctly reflects the $R = 0$ loading condition and prevents dropping all 17 tests.

### d062
* **Date:** 2026-10-19
* **Topic:** Benchmark execution host
* **Status:** Active
* **Decision:** Execute full regression snapshot sweeps on WorkstationUniPi.
* **Rationale:** Recurring Windows kernel / OpenBLAS memory faults on WorkstationUniTN (Intel Core i9-14900K Raptor Lake instability).

### d063
* **Date:** 2026-10-19
* **Topic:** Load-scaled fatigue index ($FI_{\mathrm{load}}$)
* **Status:** Active
* **Decision:** Implement load-scaled fatigue index $FI_{\mathrm{load}} = 1/s - 1$ where $DP(s \cdot \text{case}) = DP_{\mathrm{ref}}$ (via root-finding in `src/multifat/eval/endurance.py`).
* **Rationale:** Direct damage parameter ratios ($DP / DP_{\mathrm{ref}} - 1$) distort relative stress capacity for criteria where damage scales quadratically with stress or strain (e.g., SWT, Smith–Watson–Topper). $FI_{\mathrm{load}}$ accurately recovers the equivalent load-scaling margin.

### d065
* **Date:** 2026-10-19
* **Topic:** Branch consolidation
* **Status:** Active
* **Decision:** Fast-forward `dev/tg-carry-forward` to `0e87e95` and prune temporary integration branches.

### d066
* **Date:** 2026-10-19
* **Topic:** Pipeline edge case resolution
* **Status:** Active
* **Decision:** Resolve nine deferred review items: ignore Excel lock files (`~$*.xlsx`), reject duplicate material IDs, guard against non-finite values in summaries, and preserve case identity metadata.

### d067
* **Date:** 2026-10-20
* **Topic:** Integration of upstream PRs #8 and #9
* **Status:** Active
* **Decision:** Port isolated harmonic Findley damage evaluation (`fin_harmonic.py`) and configurable calibration amplitudes from upstream `ab0e828`.

### d068
* **Date:** 2026-10-20
* **Topic:** Configurable Findley calibration amplitude
* **Status:** Active
* **Decision:** Provide `--shear-amplitude-algorithm` option in `multifat-calibrate` (`mce_harmonic`, `mcc_harmonic`, `rss`). Production default remains legacy MCE.

### d069
* **Date:** 2026-10-20
* **Topic:** Findley MCC calibration plane search
* **Status:** Active
* **Decision:** Execute MCC calibration via the `longest_chord` plane rule (`FindleyMCC`) using `multifat`'s refined plane search.

### d070
* **Date:** 2026-10-20
* **Topic:** Calibration record algorithm guard
* **Status:** Active
* **Decision:** Reject non-MCE calibration records when loading via `--constants` into the production pipeline.

### d071
* **Date:** 2026-10-22
* **Topic:** MIT licence adoption
* **Status:** Active (supersedes `d044`)
* **Decision:** Formally relicense `multifat` under the **MIT License** (Copyright 2026 Andrea Chiocca, Michele Sgamma, and Tommaso Grossi) in commit `f590efa`.
* **Rationale:** Permissive licensing provides maximum research utility, industry adoption, and post-processor integration. All incorporated third-party code must be permissively licensed (MIT, BSD, ISC, Apache-2.0).

### d072
* **Date:** 2026-11-05
* **Topic:** FE interface scope (issue 19)
* **Status:** Active
* **Decision:** Issue 19 delivers a solver-independent local-history contract, a universal file format, one third-party importer and a range-strategy interface with two rules (explicit pair, longest chord). The remaining range rules, cycle counting and the extension of the criteria to maxima over all load cases (exact-cp Eqs. 8-15) belong to issue 20.
* **Rationale:** Both issues need the same object, an ordered sequence of load cases at a material point; the damage formalism does not depend on how load cases are obtained or how ranges are extracted.

### d073
* **Date:** 2026-11-05
* **Topic:** FE file format
* **Status:** Active
* **Decision:** The universal file is HDF5 (format 1.0, `docs/fe-input-format.md`). `h5py>=3.10` is the optional extra `fe`, imported lazily; the package installs and runs without it.
* **Rationale:** HDF5 underlies Exodus and XDMF, supports chunked partial reads and is one file. The alternative, a folder of `.npy` arrays with a JSON manifest, avoids the dependency but is not a recognisable FE container.

### d074
* **Date:** 2026-11-05
* **Topic:** First FE adapter
* **Status:** Active
* **Decision:** VTK XML (`.vtu`, `.pvd` series) is imported with the standard library and numpy: ASCII and uncompressed inline base64, point or cell data. Appended and compressed data are rejected. Adapters for ANSYS (PyAnsys/DPF), Abaqus (`.odb`), CalculiX (`.frd`) and Exodus are follow-up issues that write the universal file.
* **Rationale:** ParaView, CalculiX, Code_Aster, FEniCS and Elmer write VTK XML, and ANSYS and Abaqus reach it through ParaView; the proprietary adapters need their own installs to test.

### d075
* **Date:** 2026-11-05
* **Topic:** FE conventions
* **Status:** Active
* **Decision:** Stress in MPa, lengths in mm, strain dimensionless with tensor shear and meaning the total mechanical strain, global Cartesian axes; stored as six components in the order `xx, yy, zz, xy, yz, xz`. Adapters take units, component order, shear convention and source frame (global, cylindrical or per-point rotation) as arguments and never infer them. IDs are integers.
* **Rationale:** The existing models take total strain `E0, E1` and derive the elastic part from the stress; a wrongly guessed convention would change numbers silently.

### d076
* **Date:** 2026-11-05
* **Topic:** FE eligibility
* **Status:** Active
* **Decision:** A history is `local_only`. A method may declare `requires_nonlocal = True`; `fe.bridge` refuses it with a diagnostic. No registry, base or model file changed, and no registered method declares it: every one takes the tensors of one material point.
* **Rationale:** A local stress result is not sufficient for nonlocal or structural methods (issue text). Integral methods integrate over planes at a point, not over volume, so they stay eligible.

### d077
* **Date:** 2026-11-05
* **Topic:** FE test oracle
* **Status:** Active
* **Decision:** The synthetic thick tube (Lame solution, axial membrane stress, Saint-Venant torsion, isotropic Hooke) is written from the textbook equations. It is checked by boundary conditions, radial equilibrium and torque and axial resultants, then used as oracle for every importer convention.
* **Rationale:** The upstream tube code is GPL-3.0 and the repository is MIT (`d071`, `d011`); the equations are textbook.

### d078
* **Date:** 2026-11-05
* **Topic:** FE two-state bridge
* **Status:** Active
* **Decision:** `fe.bridge.to_case_inputs` takes the raw tensors at the pair of load cases a range strategy selects as `S0, S1` (and `E0, E1`), with `S_sin = (S0 - S1)/2` and `S_cos = 0` as in `extract_tensor`. It is exact for two-state paths; for curved or phase-shifted histories it discards information and the maximum normal stress over all load cases is not carried. The longest-chord rule is cross-checked against the `d069` plane rule, not reused.
* **Rationale:** Keeps every existing model and fingerprint unchanged; handling general histories is issue 20.

### d079
* **Date:** 2026-11-05
* **Topic:** History bridge scope (issue 20)
* **Status:** Active
* **Decision:** Issue 20 delivers a time-history contract, ASTM E1049 scalar rainflow, longest-chord events, damage parameters per event through the existing models and the two-state bridge, and Miner accumulation into a life. Multiaxial counting (enclosing surface, Modified Wang-Brown, minimum-circumscribed-circle rainflow), the equivalent-range extractors (PDMR, MCC, MCE, prismatic hull, maximum rectangular hull, moment of inertia), the criteria that need maxima over a cycle's load cases and the other accumulation rules are specified follow-up tasks F1-F8.
* **Rationale:** The multiaxial rules need their papers read; no equation is implemented from memory.

### d080
* **Date:** 2026-11-05
* **Topic:** Proportional domain
* **Status:** Active
* **Decision:** Counting a scalar channel and evaluating two-state events is exact for proportional histories. A non-proportional history (second over first singular value of the stress path above 1e-6) raises `IneligibleHistoryError` unless `allow_nonproportional_approximation=True`; the flag and a note are recorded in the result.
* **Rationale:** Issue 20: do not apply the existing criteria to arbitrary histories without a justified formulation.

### d081
* **Date:** 2026-11-05
* **Topic:** Caller-supplied S-N curve
* **Status:** Active
* **Decision:** `SNCurve(A, b, method, variant, design_loading, source)` is supplied by the caller, with a constructor from the pipeline's in-memory fit (`A_surv`, `b`); `b >= 0` is refused. The curve must come from the same method, variant and plane set.
* **Rationale:** The repository does not persist the damage-parameter curve, and the work may not change calibration records or CSVs.

### d082
* **Date:** 2026-11-05
* **Topic:** Counting channel
* **Status:** Active
* **Decision:** The default channel is the projection of the mean-removed shear-weighted stress path on its dominant singular direction. The von Mises stress signed by the hydrostatic part, a component, a normal stress and a callable are explicit options.
* **Rationale:** Only the dominant direction is exactly affine in the load factor when a static stress is present.

### d083
* **Date:** 2026-11-05
* **Topic:** Block semantics and residue
* **Status:** Active
* **Decision:** The history is one block repeating until failure: damage per block `D`, life `1/D` blocks, `(sum of counts)/D` cycles. Residue policy `rearrange` (default: start at the highest peak, all cycles full) or `half` (ASTM E1049, with a note) is recorded. A constant-amplitude block counted with `half` gives half cycles that sum to `m`.
* **Rationale:** Rearranging is exact for strictly periodic repetition; the standard's half-cycle rule under-counts the cycle that spans consecutive blocks (about a factor 2 on the shipped example, found by review), so it is not the default for a repeating block. This reverses the first draft of the spec, which made `half` the default.

### d084
* **Date:** 2026-11-05
* **Topic:** Eligibility
* **Status:** Active
* **Decision:** Native-life methods (MDC, MMCCM, MOREL, TVEIT2025 at this date) are refused for accumulation; a criterion whose loading scope excludes an event (`LoadingScopeError`) is reported as ineligible with the method and the event; strain-needing methods need strain; unknown or missing variants are refused.
* **Rationale:** A single Basquin curve of the damage parameter is not a native-life method's formulation; nothing is skipped silently.

### d085
* **Date:** 2026-11-05
* **Topic:** History package and literature corpus
* **Status:** Active
* **Decision:** New package `multifat.history`, imported by no existing module and reusing `multifat.fe`. The literature for the follow-ups is the Zotero collection "multifat #20 load histories" (key `4J6FG4BS`, 42 papers); PDFs are not fetched yet and the ASTM E1049 example counts are confirmed only by secondary sources.
* **Rationale:** Keeps every existing number unchanged; the corpus is the source list for the follow-up tasks.

### d086
* **Date:** 2026-11-06
* **Topic:** FE CSV exchange format
* **Status:** Active
* **Decision:** A CSV file with a commented header is the recommended way to bring FE results in (the maintainer's proposal on issue 19, tightened): one row per location and step, named columns matched case-insensitively with Abaqus and Ansys aliases, a header that declares the solver, the units, the frame (`global_cartesian` only), the averaging, the analysis type, the step role, and the strain convention and unit when strain exists. The HDF5 file stays as the internal store that the importer fills step by step. Nothing is inferred: a missing key or unit, an unknown column, a column with its alias, an incomplete strain group, blank, `nan` or `inf` cells are refused with a diagnostic.
* **Rationale:** Layers instead of a choice: a text format is easy to write from any solver and to read by eye, the chunked file is what a large model needs. The maintainer's format is better for users; the earlier design is better for storage.

### d087
* **Date:** 2026-11-06
* **Topic:** Location identity and file version 1.1
* **Status:** Active
* **Decision:** A location is the triple `(instance, node_id, material)`; keys `(instance, node_id, material, step)` are unique and the grid is complete. History file version 1.1 adds optional `instance` and `material` labels, the split `strain_elastic` and `strain_plastic`, `temperature`, solver `step` numbers, `analysis_type`, `cycle`, `step_role` and `load_description`, and two averaging values (`nodal_average_within_material`, `nodal_average_global`). Version 1.0 files still open; `select_by_ids` takes optional `instance` and `material` and returns every match when they are omitted.
* **Rationale:** With nodal averaging within materials, an interface node has one value per material; the same label in two instances is two nodes. Collapsing them would silently mix materials.

### d088
* **Date:** 2026-11-06
* **Topic:** Strain composition and derivation
* **Status:** Active
* **Decision:** Total `E_*` is the total mechanical strain; `Eel_*` with `Epl_*` are stored and summed; `Eel_*` alone is the total strain for a linear-elastic file and refused for an elastic-plastic one; `Epl_*` alone needs `elastic_constants` per material and the elastic part is derived by Hooke's law and flagged `derived_elastic_part`; a linear-elastic file without strain may derive it from `elastic_constants` (`derived_linear_elastic`); an elastic-plastic file without any strain is refused and must declare `cycle: stabilised`. `TimeHistory.from_csv` also accepts the node-CSV names and applies the same rules (`Eel_*` alone refused). Strain units `1`, `mm/mm`, `m/m`, `percent`, `microstrain` and the US-customary `psi`, `ksi`, `in`, `ft` are explicit tables.
* **Rationale:** Elastic strain is not the total strain under plasticity, and a derived quantity must not be mistaken for solver output. Units are tables rather than `pint` or `xarray` because those do not check a declared unit against the data.

### d089
* **Date:** 2026-11-06
* **Topic:** Step order and step role
* **Status:** Active
* **Decision:** Steps ascend; non-contiguous numbers are allowed and recorded; `time` must be strictly increasing along ascending steps and is never re-sorted. A step-major file streams with the memory of one step; a node-major file is sorted in memory under a row limit and larger ones are refused with advice. `step_role: unit_case` (linear elastic only) is superposed with a load-function CSV of an axis column (`time` or `index`) and one `f_<step>` column per unit case, matched by name. The axis of a node history is the time if present, else the solver step numbers as a cycle axis.
* **Rationale:** Sorting hides export errors; streaming is the only way to import large models; name matching removes the position errors of unit-case tables.

### d090
* **Date:** 2026-11-06
* **Topic:** Node assessment
* **Status:** Active
* **Decision:** `assess_nodes` takes one `MaterialSetup(model, params, curve, E=None, nu=None)` or a mapping by material label, runs `assess_history` on each selected node and returns one row per node with a status: `ok`, `ineligible: <reason>` or `no_material_setup`. `on_refusal="record"` (default) continues, `"raise"` stops at the first refusal; unexpected exceptions are never caught. With `E` and `nu` a stress-only linear-elastic file gets its Hooke strain at assessment, flagged in the provenance. No setup is resolved from workbooks in this round.
* **Rationale:** A node-level table must survive a few nodes that a criterion cannot treat; the reason has to be visible, not a gap.

### d091
* **Date:** 2026-11-06
* **Topic:** Header discipline and integrity
* **Status:** Active
* **Decision:** Required header keys, the value sets, `n_nodes` and `n_steps` checked against the file, undeformed coordinates required and checked across steps, `nodal_average_global` across several materials or instances refused unless `allow_global_average=True` (recorded), stress shear never scaled, engineering shear halved only for strain. The importer removes the partial output on any failure and records the CSV's SHA-256, the units, the conventions and the column aliases used in the file.
* **Rationale:** A convention error is silent and large; every convention is therefore a declared header key or a refusal.

### d092
* **Date:** 2026-11-06
* **Topic:** Writer helper and solver checklists
* **Status:** Active
* **Decision:** `FeCsvWriter` uses only the standard library so that a solver script can copy it; `write_tube_csv` writes the analytical tube in the canonical, Abaqus and Ansys column names as the test oracle. The Abaqus and Ansys instructions in `docs/fe-csv-format.md` are checklists, not tested code, with a pure-shear check for the shear convention; tested exporter scripts and one real file per solver are requested from the maintainer.
* **Rationale:** No solver is available to test against here; the documentation says what has and has not been checked.

### d093
* **Date:** 2026-11-06
* **Topic:** Solver adapters target the CSV gate
* **Status:** Active
* **Decision:** `fe/solvers` writes atomic FE CSV through FeCsvWriter after canonical alias/header/grid validation; no direct HDF5 output. Explicit units, analysis, step role, global frame and relevant strain/cycle/time declarations are mandatory. Original solver identities and provenance remain in the header; CSV states are ascending ordinals.
* **Rationale:** One inspectable exchange format and one internal importer prevent competing conventions.

### d094
* **Date:** 2026-11-06
* **Topic:** Abaqus nodal reports
* **Status:** Active
* **Decision:** Read complete Field Output Reports with explicit nodal position, tensor columns, instance, state/time and averaging metadata; undeformed coordinates come from flat global explicit INP nodes or COORD columns. Collapse element-node duplicates only on component agreement (1e-10 relative, 1e-14 absolute in source units); never average. Refuse ambiguous coordinates, local frames, shell section points and material conflicts. DAT and ODB are outside scope.
* **Rationale:** Flattening interfaces or section points silently would destroy the physical meaning of a nodal tensor.

### d095
* **Date:** 2026-11-06
* **Topic:** Ansys text subset
* **Status:** Active
* **Decision:** Read complete PRNSOL tensor tables with global NLIST coordinates or complete named Workbench tab/comma tables. Require averaging explicitly; check stated units and local-frame markers. Disjoint strain tables join by node. Refuse scalar output, shell/element section points, duplicate component rows and packed fixed-width numbers. No suitable real PRNSOL stress file was located; tests are hand-written and documentation states that limitation.
* **Rationale:** A strict supported subset is reviewable without implying that every Mechanical export layout was tested.

### d096
* **Date:** 2026-11-06
* **Topic:** Optional native RST stress adapter
* **Status:** Active
* **Decision:** Optional `ansys-mapdl-reader>=0.56,<0.57` installed/imported under Python 3.12/numpy 2.5.4. Accept single-material default-option SOLID185/186/187, real static/transient small-displacement linear elasticity. Global stress rotation is verified in source and against public RST element tensors. Export stress only because EEL rotation differs. Refuse unsupported states, materials, elements and non-finite stresses; node filtering must be explicit. No project numpy downgrade.
* **Rationale:** Native reader averaging and private API assumptions must remain bounded and visible.

### d097
* **Date:** 2026-11-06
* **Topic:** Solver evidence and untrusted public files
* **Status:** Active
* **Decision:** Each download uses a fresh empty directory after license inspection; no downloaded code/input is executed. Commit only small MIT Abaqus RPT/INP with license/source. Native RST is opt-in external test data; ordinary tests never fetch. Document real versus hand-written coverage, provenance hashes and component mutation results. No fatigue validation is claimed.
* **Rationale:** Import fidelity, regression stability and experimental validation are different evidence.

### d098
* **Date:** 2026-11-06
* **Topic:** Whole-branch import and history safeguards
* **Status:** Active
* **Decision:** Validate CSV row width before pandas projection and use round-trip float parsing; CSV/PVD imports publish atomically and refuse source aliases. Validate report-wide local/section/material declarations even before a frame. Numeric solver steps ascend. `TimeHistory.from_file` refuses ambiguous node labels and unit cases. Strain-reading criteria require one shared scalar stress/strain load function, tested after independent path normalization, and refuse constant-stress cycling strain unless the explicit approximation flag is set.
* **Rationale:** The whole-branch review reproduced source deletion and silent changes of physical meaning that were absent from ordinary round trips. Each guard has a failing-then-passing regression test.

### d099
* **Date:** 2026-11-06
* **Topic:** In-house Ansys fatigue benchmark and external data purge
* **Status:** Active
* **Decision:** Purge repository of third-party repository fixtures (`axis-fem` removed). Generate an in-house 2D plane stress holed-plate benchmark using Ansys MAPDL 2026 R1 batch mode on WorkstationUniTN (`PLANE183`, `KEYOPT(3)=3`) under a 2-step alternating fatigue cycle ($Y$-tension 100 MPa then $X$-tension 100 MPa, rotating principal stress directions by $90^\circ$). Support PLANE182/183 elements in `rst.py`. Ingest both `.out` and `.rst` into exchange CSV and HDF5, verified against Howland's finite-width analytical solution ($K_{tg} \approx 3.032$). Mark real Abaqus output tests as xfail TODO pending Abaqus license/installation on workstation.
* **Rationale:** Eliminates external repository dependencies, guarantees clean provenance, and directly exercises multi-step biaxial fatigue loading.
