# FE input format

Local stress and strain histories from finite-element (FE) analyses, in one solver-independent
form. This page is the specification of the contract (`multifat.fe`), of the universal file, of the
importers (the CSV exchange format, [fe-csv-format.md](fe-csv-format.md), and VTK XML) and of
the bridge to the existing models. It prepares the framework for FE data; it does not run a component
assessment. **Users bring results in through the CSV file**; the HDF5 file below is the internal
store that the CSV importer fills.

Software regression and scientific validation are kept apart: the tests below show that the
importer is faithful to an analytical solution and that no existing number moved. Neither validates
a fatigue method against experiments.

## What a history is

An FE history is the sequence of `N` ordered **load cases** `(sigma_i, epsilon_i)` of `P` material
**locations**, following Sec. 2.2 of the exact critical-plane paper (Grossi, Chiocca and Sgamma):
where the load cases come from (FE steps, strain gauges, a cycle counter) is separate from the
**range strategy** that turns them into ranges. Between consecutive load cases the tensor is
assumed to move along a straight segment; the producer must supply load cases dense enough for
that (for example after plastic shakedown or ratcheting). The importer cannot verify it.

| Part | Content |
|---|---|
| Locations | `P` of them; stable integer `location_id` (node, element or integration-point label), coordinates `(P, 3)` in mm, global frame; one `location_kind` for the whole file (`node`, `element_centroid`, `integration_point`, `element_nodal`); no per-location kind or surface flag. Since 1.1 a location is the triple `(instance, location_id, material)` when the optional `instance` and `material` labels exist (the same label in two instances, or an interface node once per material, are distinct locations), else the `location_id` alone |
| Load cases | ordered, `N >= 1`; optional time or cycle coordinate and solver labels; flag `periodic` (last connects to first) |
| Stress | `(P, N, 3, 3)` symmetric, MPa, global Cartesian axes |
| Strain | optional, same shape, dimensionless **tensor** strain (shear = engineering shear / 2); the **total mechanical strain** (elastic plus plastic), as the models take `E0, E1`; `strain_kind = "total_mechanical"`, or `derived_linear_elastic` / `derived_elastic_part` when the importer derived the elastic strain by Hooke's law. Since 1.1 the elastic and plastic parts may be stored separately (`E_elastic_hist`, `E_plastic_hist`) and a temperature `T_hist` `(M, N)` |
| Material | `material_id` (free text; equal to a workbook material ID when one exists) |
| Averaging | `none`, `nodal_average`, `extrapolated`, `nodal_average_within_material` or `nodal_average_global`: how the solver produced the values at the locations |
| Analysis (1.1) | `analysis_type` (`linear_elastic`, `elastic_plastic`, `unspecified`), `cycle`, `step_role` (`time_sample`, `unit_case`), `load_description` |
| Provenance | solver and version, analysis step, source files with SHA-256, adapter and version, source units and frame, source conventions (component order, strain shear, point or cell data, column aliases used), creation time, notes |

Symbols follow [notation.md](notation.md): `S_hist`, `E_hist`, `location_id`, `S_range`, `E_range`,
`pair`.

An FE history is a **prediction**. It has its own types (`HistoryBlock`, `HistoryMeta`) that the
workbook importer, the calibration and the comparison do not accept, and nothing here changes how
experimental observations are imported or calibrated.

A local stress result is not sufficient for nonlocal, gradient, volume-integral or structural
methods. A history declares `local_only = True`, and a method that needs more declares
`requires_nonlocal = True`; the bridge then refuses it with a diagnostic. All registered models take
the tensors of one material point, so none is refused today.

## The universal file (HDF5, format 1.1)

HDF5 is the de facto container of FE results (Exodus, XDMF), allows chunked partial reads and is one
file. `h5py` is an optional dependency: `uv sync --extra fe` (or `pip install multifat[fe]`). The
rest of the package, including the VTK reader, runs without it.

| Path | Content |
|---|---|
| root attributes | `format = "multifat-fe-history"`, `format_version = "1.1"`, `complete` (true only after every load case is written; a file without it is refused on open) |
| `/meta` | attribute `json`: the metadata above |
| `/locations/id`, `/locations/coordinates` | integer IDs `(P,)`, coordinates `(P, 3)` in mm; optional `/locations/instance`, `/locations/material` (UTF-8 strings, 1.1) |
| `/load_cases/index` | `0 .. N-1`; optional `time` and `step` (solver step numbers, 1.1) datasets and `labels` attribute |
| `/fields/stress`, `/fields/strain` | `(P, N, 6)`, float64 (float32 is a storage option), chunked along locations (about 4 MiB per chunk, whatever `N`) so one location's history is contiguous |
| `/fields/strain_elastic`, `/fields/strain_plastic`, `/fields/temperature` | optional (1.1): `(P, N, 6)` and `(P, N)` |
| `/mesh` | empty group, reserved for connectivity (not read in 1.0) |

The six components are stored in the fixed order `xx, yy, zz, xy, yz, xz`, with tensor shear. The
writer checks that every tensor is symmetric and finite, that IDs are unique integers (unique as a
triple when instance or material labels exist) and that every load case was written. A version 1.0
file still opens: the new fields are simply absent.

```python
from multifat.fe.store import HistoryFile
from multifat.fe.select import select_in_box, screen_by, von_mises_max

with HistoryFile("tube_history.h5") as f:
    print(f.meta.solver, f.n_locations, f.n_load_cases)
    hot = screen_by(f, von_mises_max, top=10)           # shortlist, not a damage claim
    for block in f.iter_chunks(max_bytes=64_000_000, index=hot):
        ...                                              # block.S_hist: (M, N, 3, 3)
```

`iter_chunks` yields blocks of locations with all load cases, bounded by a number of locations or a
byte budget; peak memory depends on the chunk, not on `P`. The byte budget is an estimate, not a hard cap: the numpy temporaries of one chunk (copies, the expanded tensors and the symmetry check) stay below three times it, as the tests check; HDF5's own buffers are not counted. Files written by other tools with unchunked (contiguous) datasets are read too. Selection helpers (`select_by_ids`,
`select_by_mask`, `select_in_box`, `screen_by`) return sorted file positions; `select_by_ids` takes
optional `instance` and `material` filters and returns every match when they are omitted (an
interface node gives one position per material).

## CSV importer

`multifat.fe.csvio.import_fe_csv` reads the FE CSV exchange format
([fe-csv-format.md](fe-csv-format.md)) one load step at a time, validates the header, the columns and
the grid, converts units and shear, composes the strain by stated rules and writes this file. It is
the recommended entry point for solver results.

## VTK XML importer

`multifat.fe.vtk` reads unstructured grids (`.vtu`) and `.pvd` time series (one file per load
step). Supported: ASCII and inline base64 `binary` arrays (uncompressed, `UInt32` or `UInt64`
header, little endian), point data or cell data (a cell value sits at the cell centroid).
**Appended data and compressed data are rejected** with a diagnostic naming the feature; other
variants are not claimed. In ParaView, save with "Data Mode: Binary" or "ASCII" and compression
off.

Nothing is inferred. The caller states:

| Argument | Meaning |
|---|---|
| `TensorSpec(array, order, unit, shear, location)` | array name, component order of the columns (`CANONICAL`, `VTK9` or any permutation of `xx, yy, zz, xy, yz, xz`), stress unit (`Pa`, `kPa`, `MPa`, `GPa`; required for stress, ignored for strain), shear convention of the strain (`tensor` or `engineering`), `point` or `cell` |
| `length_unit` | `m`, `cm` or `mm`; required |
| `frame` | `"cartesian"`; `CylindricalFrame(axis, origin)` (tensors in `(r, theta, z)`; origin in the source length unit); or `RotationArray(name)` (a 9-component per-point rotation, columns are the local axes) |
| `id_array` | integer array of location IDs; if omitted, the VTK point or cell index is used and the notes say so |

```python
from multifat.fe.tensors import CANONICAL
from multifat.fe.vtk import TensorSpec, import_pvd

meta = import_pvd(
    "vtk/tube.pvd", "from_vtk.h5",
    stress=TensorSpec("Stress", CANONICAL, unit="MPa"),
    strain=TensorSpec("Strain", CANONICAL, shear="tensor"),
    id_array="NodeID", material_id="42CrMo4_QT", length_unit="mm",
    solver="my-solver", solver_version="1.0", analysis_step="step 1",
)
```

Base64 `binary` arrays may carry the byte-count header in the same stream as the data (the layout VTK writes) or as a block of its own.

The conversion streams one load step at a time, so peak memory is one load step of the fields
(order `P x 6` values per field), never `P x N`; a single load step is assumed to fit in memory.
Location IDs and coordinates must agree across steps and the partial output is removed on failure. A `.pvd` with several parts or repeated timesteps is rejected, and a `RotationArray` must hold proper rotations (orthonormal, determinant +1). The component order and strain shear convention are recorded in the file.
The SHA-256 of every source file is recorded in the file.

Solver-specific adapters for ANSYS (PyAnsys/DPF), Abaqus (`.odb`), CalculiX (`.frd`) and Exodus are
follow-ups; they need their own installs to test and must write this file or the CSV.

## Worked example

```bash
python -m multifat.fe.synthetic --out out
```

writes `out/tube_history.h5` and `out/vtk/tube.pvd` for the same three load cases of an elastic
thick-walled tube under internal pressure, axial force and torque (Lame and Saint-Venant solution,
closed ends, inner radius 6 mm, outer radius 7.5 mm), and the two import back to the same tensors
(`tests/test_fe_synthetic.py::test_documented_example_runs`). The analytical field is the test
oracle: boundary conditions, radial equilibrium, torque and axial resultants are checked
independently of the formulas, and every importer convention (component order, shear convention,
units, cylindrical to Cartesian frame, ASCII and binary encodings, point and cell data) is
recovered to `1e-10` relative. The solution is textbook and written from scratch; it excludes end
effects, plasticity and residual stress.

## Ranges and the bridge to the existing models

A **range strategy** maps the `N` load cases of one location to `S_range = S[a] - S[b]` (and
`E_range`) for a selected pair `(a, b)`; peak quantities such as the maximum normal stress are
maxima over all load cases and do not come from it.

| Strategy | Pair |
|---|---|
| `ExplicitPair(a, b)` | supplied by the user or a cycle counter |
| `LongestChord(metric, on)` | the pair farthest apart in the space of symmetric tensors (both distances ignore hydrostatic change, so a pure pressure cycle is a tie and gives pair `(0, 1)`); `metric` is the Frobenius norm of the deviatoric part (`deviatoric_norm`) or the spread of the principal values (`principal_spread`), on stress or strain; ties go to the smallest pair; one load case gives pair `(0, 0)` and a zero range |

The longest-chord distances are provisional definitions; their scientific review and the other rules
(rainflow-based ones) belong to the variable-amplitude work (issue 20). For a pure-shear ellipse the
chord equals the existing `longest_chord` plane rule of `eval/evaluator.py` (tested).

`multifat.fe.bridge.to_case_inputs(block, strategy, meta=..., method=...)` turns the selected pair
into `S0, E0, S1, E1` and the harmonic dictionary of the existing models, with the convention of
`imp.experimental_data.extract_tensor` (`S_sin = (S0 - S1) / 2`, `S_cos = 0`). **The reduction is
exact only for a two-state path.** For curved or phase-shifted histories it discards information
(a circular shear path becomes a straight line through it, with no quadrature component), and the
maximum normal stress over all load cases does not reach the two-step interface. Extending the
criteria to those histories is issue 20. Without strain, `E0`, `E1` and every `E_*` are NaN; when a
`method` is given, the bridge refuses every method that is not stress-only (`IneligibleMethodError`
naming the missing strain), so no model is left to return NaN or fail on its own.

## Limits

- Straight segments between consecutive load cases; dense enough load cases are the producer's job.
- Single-point (local) data only; no gradients, volume or connectivity (reserved `/mesh` group).
- One load step must fit in memory during VTK conversion.
- Elastic-only data can omit strain; strain-based criteria then need it from the material.
- Integer location IDs only.

Decisions: [decision-log.md](decision-log.md), `d072` to `d078` (interface) and `d086` to `d092`
(CSV format and file version 1.1).
