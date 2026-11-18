# FE local-history interface: design

Issue: [andreachiocca/multiaxial-fatigue#19](https://github.com/andreachiocca/multiaxial-fatigue/issues/19),
"Define an interface for importing finite-element results" (assignee: Tommaso Grossi).
Date: 2026-11-05. Status: draft for approval, nothing built.

## Intent and why

The framework evaluates criteria on experimental workbooks only. Finite-element (FE) analyses give
local stress and strain tensors at many material points and load steps, which is the input a
critical-plane assessment of a real component needs. FE solvers differ in file format, tensor
component order, shear convention, units, coordinate frame and data volume, so each criterion cannot
be taught to read each solver.

The goal is to be **prepared** for FE data, not to run a component assessment yet: define one
solver-independent contract for local tensor histories, a recommended "universal" file format that
carries it, one working third-party importer, and a bounded-memory way to select and stream
locations. Solver-specific adapters (ANSYS, PyAnsys, Abaqus, others) later only have to produce
that file.

This work is deliberately coupled to issue #20 (time-history inputs and a variable-amplitude
pathway). Both need the same object: an ordered sequence of load cases at a material point.

## Scientific frame (decided)

The history at a material point is a finite sequence of N load cases `(sigma_i, epsilon_i)`, as in
Grossi, Chiocca and Sgamma, *Exact Critical Plane Search for a Class of Damage Criteria under General
Load Histories* (`exact-cp`, Sec. 2.2). Two parts are kept separate, as in that paper:

- **Load cases**: where they come from (FE steps, strain gauges, a cycle counter) is not part of the
  damage formalism.
- **Range strategy**: a rule that maps the load cases to the ranges `Delta sigma = sigma_a - sigma_b`,
  `Delta epsilon = epsilon_a - epsilon_b`. Peak quantities (maximum normal stress and the like) are
  maxima over the load cases (exact-cp Eqs. 8, 9, 15).

Between consecutive load cases the tensor is taken to move along a straight segment (exact-cp
Eq. 6). The contract states this as an assumption the data producer must satisfy (dense enough
load cases, e.g. after plastic shakedown or ratcheting); the importer cannot verify it.

The existing case-level models take two load steps `S0, S1, E0, E1`. That is the current model
interface, not a limit of the formalism. #19 does not change any model.

## Goals

1. A solver-independent in-memory object for local tensor histories: P locations x N load cases,
   stress and optional strain, with units, load-case axis, location identity, material identity,
   analysis step and provenance (below).
2. A documented universal file format (HDF5) that stores it, with a written specification and a
   worked example.
3. One representative third-party importer: VTK XML (`.vtu`, time series as `.pvd`), with explicit
   tensor-convention arguments. Standard library plus numpy only.
4. A synthetic FE-export generator based on a thick-walled tube under pressure, axial force and
   torque (textbook Lame and Saint-Venant solution), used as test oracle and as the documented
   example.
5. Location selection and bounded-memory chunked iteration, so peak memory does not grow with the
   number of locations.
6. A range-strategy interface with two implementations: explicit pair and longest chord (see
   below), plus a thin bridge that turns the selected pair into the existing `S0, S1, E0, E1` model
   inputs.
7. Tests for tensor transformations, units, location identities, provenance and large-input
   handling, and proof that existing workbook results do not change.

## Non-goals and anti-goals

- No solver-specific adapters for ANSYS (PyAnsys/DPF), Abaqus (`.odb`) or CalculiX in this issue.
  They are specified as follow-up issues; they need proprietary or separate installs to test. They
  must target the universal file.
- No component-level fatigue assessment. A local stress result is not sufficient for nonlocal,
  gradient, volume-integral or structural methods (issue text). The object carries no neighbour or
  gradient information; methods that need it must be reported as ineligible with a diagnostic, never
  evaluated silently on local data. The file format reserves an optional mesh-connectivity group
  (not read in v1) so adding it later does not break files.
- FE predictions never enter the experimental dataset, calibration or the comparison tables. They
  have their own type, which the workbook importer, calibration and `multifat-compare` do not accept.
- No inference of units, tensor convention or frame. Every adapter takes them as arguments or reads
  them from the source; a missing one is an error.
- No claim of validation. Regression tests show numbers did not change; the synthetic oracle shows
  the importer is faithful. Neither validates a fatigue method against experiments.
- No GPL code. SeAn-CP (GPL-3.0) and the repository's other upstream tools are not copied; the
  repository is MIT (d071). The tube solution is re-derived from textbook equations.
- No change to any model, calibration record, fingerprint or existing CSV. The new code is not
  imported by existing pipelines.
- No optional-dependency leakage: `h5py` is an optional extra (`fe`), imported lazily with a clear
  error; the rest of the package installs and runs without it.

## Constraints from the repository

- `AGENTS.md` conventions apply: numpydoc docstrings with the module docstring on line 1; black
  (line length 100) and ruff; **notation is data**: every new public parameter name needs a row in
  `docs/notation.md` (or a reasoned entry in `tests/data/notation_allow.txt`) before it is used;
  `E` is Young's modulus, strain tensors carry the `E_` prefix; MPa and mm.
- Workbooks and archived results are never edited. Decisions are logged in `docs/decision-log.md`
  (next free ID after `d071`). `docs/audit-checkpoint.md` is updated at each milestone.
- The regression snapshot (`tools/regression_snapshot.py`) must run unchanged over the new branch as
  the evidence that no existing number moved.
- A longest-chord plane rule already exists (decision `d069`, `FindleyMCC`). The longest-chord range
  strategy must reuse or reference it, not duplicate it.

## Design

### 1. The history object

One immutable description of FE-derived local histories, read lazily from a file:

| Part | Content |
|---|---|
| Locations | `P` locations; stable source ID (node or element/integration-point ID), coordinates `(P,3)` in mm in the global frame; the `kind` (`node`, `element_centroid`, `integration_point`, `element_nodal`) is one value per file; a per-location kind and a surface flag are left to a later format version (dropped from 1.0 after the final review) |
| Load cases | ordered `N`; index, optional time or cycle coordinate, solver step/increment label; flag `periodic` (last load case connects to first) |
| Stress | `(P, N, 3, 3)` symmetric, MPa, global Cartesian frame |
| Strain | optional, same shape, dimensionless **tensor** strain (shear = engineering shear / 2). It is the **total mechanical strain** (elastic plus plastic), the quantity the existing models take as `E0, E1` (they derive the elastic part from the stress). A `strain_kind` attribute records this; an optional separate `strain_plastic` field is reserved for later, not read in v1 |
| Material | `material_id` (free text; matches a workbook material ID when one exists) |
| Averaging | `none`, `nodal_average` or `extrapolated`: how the solver produced values at the locations (nodal averaging smooths discontinuities and must be visible) |
| Provenance | solver name and version, source file names with SHA-256, adapter name and version, source units and frame, creation time, free-text notes |

Elastic-only FE data may omit strain; strain-based criteria are then ineligible unless the material
supplies elastic constants, using the repository's existing elastic strain conversion.

### 2. Universal file format (HDF5)

Recommended and chosen: HDF5, because it is the de facto FE-results container (Exodus, XDMF), supports
chunked partial reads, and is one file. Layout (version `1.0`; exact names are the builder's, the
spec document is the authority):

- root attributes: `format = "multifat-fe-history"`, `format_version`;
- `/meta`: provenance and material attributes above;
- `/locations`: IDs, kinds, coordinates, flags;
- `/load_cases`: index, time or cycle, labels;
- `/fields/stress`, `/fields/strain`: the six unique components in a **fixed documented order**
  (`xx, yy, zz, xy, yz, xz`, tensor shear), shape `(P, N, 6)`, float64 (float32 allowed as a
  documented storage option), chunked along locations so a location's full history is contiguous;
- reserved empty group `/mesh` for future connectivity.

Writers validate symmetry and finiteness; readers expand to `(.,.,3,3)`. The format specification
(`docs/fe-input-format.md`, hand-written) lists every field, unit, convention and the validation
rules, and includes the worked example.

### 3. VTK XML importer (first third-party adapter)

Reads `.vtu` and a `.pvd` time series (one `.vtu` per load step). Explicit arguments:
point or cell data, array names, component order of 6- or 9-component tensors, strain-shear
convention (`tensor` or `engineering`), source stress units, and the source frame. Supported
source frames: global Cartesian (default); a cylindrical frame about a declared axis and origin
(needed because the tube oracle is natural in `(r, theta, z)`); and a per-point rotation array
(the general case; the cylindrical option just builds it). Everything is rotated to global
Cartesian before it enters the history, and the source frame is recorded in the provenance. Handles the
encodings the tests can actually produce and verify (ASCII, inline base64 binary; ParaView's default
appended/zlib output only if a test fixture for it is built and passes). It must reject, with a
diagnostic, anything it does not support, never guess.

### 4. Synthetic generator and oracle

Thick-walled tube: inner radius `a`, outer radius `b`, harmonic or stepped pressure, axial force
and torque, evaluated at a set of radii and angles. Stress from the textbook Lame solution plus
axial membrane stress plus Saint-Venant torsion, in `(r, theta, z)`, rotated to global Cartesian
(this exercises the frame transformation). Strain from isotropic Hooke. The generator writes
(a) the universal file directly, (b) a `.pvd`/`.vtu` series with chosen tensor conventions and
units (MPa and Pa, 6- and 9-component orders), and the analytical tensors are the oracle: a round
trip through either route must reproduce them within a stated tolerance.

### 5. Selection and bounded memory

- Selection by location IDs, boolean mask, coordinate box, and an optional scalar screening
  function (e.g. rank by maximum von Mises over the history, to shortlist locations; a screening
  aid, not a damage claim).
- `iter_chunks(...)`: yields contiguous blocks of locations bounded by a number of locations or a
  byte budget, each with all N load cases. Peak memory depends on the chunk, not on `P`.
- The whole path is bounded, not only reading. The `.pvd` import is step-major while the file is
  location-major, so conversion streams one load step at a time: read one `.vtu`, rotate, write the
  slab into `[:, step, :]`. Peak memory is one load step of the fields (order `P x 6` per field),
  never `P x N`. Screening scans the file chunk by chunk. A single load step is assumed to fit in
  memory; the documentation states this limit.

### 6. Range strategy and the bridge to existing models

Interface: given the `(N, 3, 3)` stress (and strain) of one location, return `Delta sigma`,
`Delta epsilon` and the selected pair of load-case indices `(a, b)`. Implementations in #19:

- **explicit pair** (user or cycle counter supplies `a, b`);
- **longest chord** over pairs of load cases under a selectable distance (deviatoric norm of the
  difference; spread of principal values). Documented as a provisional definition; its scientific
  review and the other rules (rainflow-based) belong to issue #20.

Bridge: from a selected pair build `S0, S1, E0, E1` for the existing case-level models, where
`S0, S1, E0, E1` are the raw tensors at load cases `a` and `b` (the existing interface derives
`S_mean` and `S_sin` from them, see `docs/notation.md`; the builder confirms the convention in the
code). It is exact only when the history is a two-state path; for curved or phase-shifted histories
the two-state reduction discards information, and the maximum normal stress over all load cases is
not carried by the two-step interface. The documentation says so explicitly. Handling of those
cases is #20.

Longest chord considers all pairs of load cases, so the `periodic` flag does not change it; the
flag matters only for explicit consecutive-pair selection and is recorded for #20. The existing
`longest_chord` plane rule (`d069`, in `eval/evaluator.py`) is a private vectorised function on
rotated harmonic arrays and cannot be called on tensors as it stands, so the new rule is
cross-checked against its closed form in a test. Refactoring the existing code is out of scope.

Eligibility diagnostic (a deliverable, with the narrowest mechanism): every registered method
takes the tensors at one material point (`evaluate_case` receives only `S0, E0, S1, E1` and
constants), so none can use gradients or volume data and all are eligible on local FE data
(integral methods integrate over planes at a point, not over volume). The mechanism is therefore a
capability check: a method may declare `requires_nonlocal = True` (read with `getattr`, default
`False`), the history declares itself `local_only`, and the bridge refuses such a method with a
diagnostic naming the reason. No registry, base or model file is edited and no registered method
sets the attribute; a stub method tests the refusal. (Corrected after reading the code: an earlier
draft named the integral family as ineligible.)

## Verification (maps to the acceptance criteria)

1. **Specification and example published**: `docs/fe-input-format.md` plus the runnable synthetic
   example; `tools/check_docs.py` passes.
2. **Importer**: `.vtu`/`.pvd` series imported and written to the universal file; read back.
3. **Transformations, units, identities, provenance**: round trip against the analytical tube for
   each convention (component order, shear convention, Pa vs MPa, cylindrical vs Cartesian frame);
   location IDs and coordinates preserved; provenance and SHA-256 recorded and checked; deliberately
   inconsistent inputs (missing units, asymmetric tensor, NaN, mismatched step counts) raise.
4. **Large input**: a generated history with many locations is processed in chunks with a peak
   memory bound that is independent of `P` (opt-in slow test).
5. **No change to existing results**: `uv run python tools/regression_snapshot.py` shows identical
   output to the base commit; full `pytest`, `black --check`, `ruff`, `check_docs`,
   `check_scientific_reviews` pass. For a two-state synthetic history the bridge reproduces the
   inputs the harmonic importer would give for the same loading.
6. Type separation: a test shows that an FE history object is rejected by the workbook, calibration
   and comparison entry points.
7. Eligibility: a test shows that the bridge refuses a stub method declaring
   `requires_nonlocal = True` with a diagnostic, and accepts every registered method.

## Follow-ups (separate issues, not built here)

- PyAnsys/DPF adapter; Abaqus `.odb` adapter (runs in Abaqus Python, writes the universal file);
  CalculiX `.frd`; Exodus; long-format CSV.
- Issue #20: time-history contract details, the remaining range rules, cycle counting, and extending
  the criteria to use maxima over all load cases (exact-cp Eqs. 8-15).
- Mesh connectivity and gradient support for nonlocal methods.

## Decisions and who made them

Decided with the user (2026-11-05): scope boundary above (contract plus file format plus range-rule
interface in #19; remaining range rules and criteria extension in #20); HDF5 with an optional
`h5py` extra; synthetic fixtures first, one real export later; `.vtu`/`.pvd` as the first adapter;
ANSYS and Abaqus adapters as follow-ups.

Assistant's calls (open to change), including those added after the builder check: strain is
total mechanical strain; the VTK importer supports global, cylindrical and per-point-rotation
source frames; conversion is bounded by one load step; the bridge uses raw tensors at the selected
pair and cross-checks (does not reuse) `d069`; the eligibility diagnostic lives in the bridge and
reads a `requires_nonlocal` method attribute only; tensor-shear and six-component storage; float64 default;
`periodic` flag; the `averaging` field; fixed component order; selection and screening helpers;
longest chord with two distance options; the bridge to `S0, S1, E0, E1` and its stated limits;
reserved `/mesh` group; strain optional; lazy import of `h5py`; slow-test marker for the large-input
test; textbook re-derivation of the tube instead of reuse of SeAn-CP code.

## Left to the builder

Module and function names (subject to the notation rule), exact HDF5 dataset names and chunk sizes,
the byte budget and large-input sizes, the tolerance values, VTU encodings beyond ASCII and base64,
where the example script lives, and the split into commits.
