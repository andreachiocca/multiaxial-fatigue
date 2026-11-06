# FE CSV exchange format

The recommended way to bring finite-element (FE) results into `multifat`: one CSV file with a
commented header, written by a solver script or by hand, imported by `multifat.fe.csvio.import_fe_csv`
into the universal HDF5 history file ([fe-input-format.md](fe-input-format.md)) and then assessed node
by node ([load-histories.md](load-histories.md)). The format follows the proposal of the maintainer
(issue 19: one row per node and step, named columns, a header with the units and conventions) and
tightens it so that nothing is guessed.

Software checks and scientific validation are separate here too: the round trips below show that the
importer is faithful to an analytical solution under every convention it accepts. They do not validate
any fatigue criterion against experiments.

## The layers

```
solver script / hand  ->  FE CSV (this page)  ->  import_fe_csv  ->  HDF5 history file
                                                                      -> node_histories -> TimeHistory
                                                                      -> assess_history -> table per node
```

The CSV is what a user produces; the HDF5 file is the internal, chunked, completeness-flagged store
(version 1.1) that the importer fills one load step at a time. Both stay in the repository: the CSV is
easy to write from any solver and to read by eye, the HDF5 file is what a million-node model needs.

## Quick start

```bash
python -m multifat.fe.synthetic --csv docs/examples/fe-csv
```

writes [`docs/examples/fe-csv/tube_two_steps.csv`](examples/fe-csv/tube_two_steps.csv): four nodes of
an elastic thick-walled tube at two load steps (the unloaded state and the maximum), with the total
strain. The committed file is exactly the generated one (`tests/test_fecsv_import.py`). Its values are
the output of the analytical oracle, so figures such as `7.3e-16` are rounding noise, not data.

```python
import numpy as np
from multifat.fe.csvio import import_fe_csv
from multifat.history.accumulate import SNCurve
from multifat.history.nodes import MaterialSetup, assess_nodes
from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

import_fe_csv("docs/examples/fe-csv/tube_two_steps.csv", "tube.h5")     # needs h5py (extra fe)
(model,) = get_models(["FS"])
setup = MaterialSetup(model, params, SNCurve(A=2000.0, b=-0.12, method="FS"))   # params: your constants
rows = assess_nodes("tube.h5", setup, R_list=build_R_list(n1=3, n2=6))
```

Each row has `position`, `location_id`, `instance`, `material`, `status`, `damage_per_block`,
`blocks_to_failure`, `cycles_to_failure`, `n_events`, `approximation` and `notes`. Two steps (a maximum
and a minimum) are the two-state case and give exactly one cycle. The constants, the plane set and the
Basquin curve are the caller's, exactly as for `assess_history`; the curve must come from the same
method, variant and plane set.

## The file

Header lines first, then one column line, then the rows. UTF-8, `\n` or `\r\n`, decimal point, no
quoting (a text value may not contain a comma). A header line is `# key: value`; only the first colon
separates key and value, so free text may contain colons. A repeated key, a `#` line without a colon
and a file without a column line are refused.

### Header

Required keys (a missing or invalid one is refused, naming the key):

| Key | Value |
|---|---|
| `format` | `multifat-fe-csv 1.0` (a later 1.x minor is read) |
| `solver`, `solver_version` | free text, recorded in the provenance |
| `stress_unit` | `Pa`, `kPa`, `MPa`, `GPa`, `N/mm^2`, `psi`, `ksi` |
| `length_unit` | `m`, `cm`, `mm`, `in`, `ft` |
| `coordinate_system` | `global_cartesian` only; the exporter must write the tensors in the global frame |
| `averaging` | `none`, `nodal_average_within_material`, `nodal_average_global`, `extrapolated` |
| `analysis_type` | `linear_elastic` or `elastic_plastic` |
| `step_role` | `time_sample` (each step is a state of one load history) or `unit_case` (each step is a unit load case) |

Required when relevant:

| Key | When |
|---|---|
| `strain_convention` (`tensor` or `engineering`) and `strain_unit` (`1`, `mm/mm`, `m/m`, `percent`, `microstrain`) | the file has any strain column |
| `cycle` | `elastic_plastic` files must say `stabilised` (the shakedown or stabilised cycle); otherwise a positive integer or `not_applicable` |
| `time_unit` | the file has a `time` column |

Optional: `n_nodes` (distinct locations) and `n_steps` (distinct steps), checked against the file, so a
truncated export is refused; `temperature_unit`; free text `load_description`, `model`, `source`,
`notes`. Unknown keys are kept in the header record and otherwise ignored.

`unit_case` needs `analysis_type: linear_elastic`: superposing unit cases is meaningless under
plasticity.

### Columns

Columns are matched **by name**, case-insensitively, never by position; the order in the file is free.

| Canonical column | Meaning | Required |
|---|---|---|
| `node_id` | integer node label | yes |
| `instance`, `material` | text; absent means the empty string | no |
| `step` | integer solver step or increment number; need not start at 1 or be contiguous | yes |
| `time` | one value per step, strictly increasing along ascending steps | no |
| `x`, `y`, `z` | **undeformed** coordinates of the node, in `length_unit` | yes |
| `S_xx S_yy S_zz S_xy S_yz S_xz` | stress components, global frame, `stress_unit` | yes |
| `Eel_xx ... Eel_xz`, `Epl_xx ... Epl_xz` | elastic and plastic strain | no |
| `E_xx ... E_xz` | total mechanical strain | no |
| `T` | temperature, passed through in `temperature_unit` | no |

Solver names are accepted as aliases (matching is case-insensitive; the spelling below is the
solver's own):

| Canonical | Abaqus | Ansys |
|---|---|---|
| `node_id` | `NODE` | `NODE` |
| `S_xx` `S_yy` `S_zz` `S_xy` `S_yz` `S_xz` | `S11` `S22` `S33` `S12` `S23` `S13` | `SX` `SY` `SZ` `SXY` `SYZ` `SXZ` |
| `Eel_xx` ... `Eel_xz` | `EE11` `EE22` `EE33` `EE12` `EE23` `EE13` | `EPELX` `EPELY` `EPELZ` `EPELXY` `EPELYZ` `EPELXZ` |
| `Epl_xx` ... `Epl_xz` | `PE11` `PE22` `PE33` `PE12` `PE23` `PE13` | `EPPLX` `EPPLY` `EPPLZ` `EPPLXY` `EPPLYZ` `EPPLXZ` |
| `E_xx` ... `E_xz` | `E11` `E22` `E33` `E12` `E23` `E13` | `EPTOX` `EPTOY` `EPTOZ` `EPTOXY` `EPTOYZ` `EPTOXZ` |
| `T` | | `TEMP` |

Abaqus numbers the components 1, 2, 3 for the global axes `x`, `y`, `z`, so `S12` is `xy`, `S23` is
`yz` and `S13` is `xz`; the mapping is tested with a distinct value per component. The short names
`Sxx` ... `Sxz` are accepted too.

A file that has both a canonical column and an alias of it, an unknown column (unless listed in
`ignore_columns`), a missing required column, an incomplete strain group (all six components or none)
or a total strain `E_*` together with `Eel_*`/`Epl_*` is refused with a message that names the columns.

### Rows

One row per **location** and step. A location is the triple `(instance, node_id, material)`: the same
node label in two instances is two locations, and with `nodal_average_within_material` a node on a
material interface appears once per material, each with its own values. The keys
`(instance, node_id, material, step)` must be unique and the grid complete (every location at every
step); `n_nodes` counts locations.

* **Step-major** order (all locations of one step, then the next step) is streamed: peak memory is one
  step. Steps must ascend and each step must form one block; the order of the nodes within a step is
  free (rows are matched by key). Prefer this order for large models.
* **Node-major** order (all steps of one node, then the next) is read into memory and sorted by step;
  more than `max_node_major_rows` rows (default 2 000 000) is refused with advice to export
  step-major. A file whose first two rows are the same location is read as node-major.
* `time` must be one finite value per step and strictly increasing along ascending steps; a file
  where it is not is refused, never re-sorted. The coordinates of a node must not change between steps.
* Blank numeric cells, `nan` and `inf` are refused.
* Every row must have exactly the declared number of fields, including ignored columns. Width is
  checked before column projection so an extra field cannot be interpreted as an implicit index.
* Floating-point text uses round-trip parsing. With unchanged units and tensor shear, the writer's
  17-digit values survive import exactly; the default decimal parser can lose small strain digits.

CSV and VTK imports write to a temporary sibling file and replace the destination only after
success. Refusal preserves an existing destination. Using any source (including a referenced VTU)
as the destination, through the same path or a filesystem alias, is refused.

## Conventions the importer applies

| Item | Rule |
|---|---|
| Stress | converted to MPa; components are never scaled by the shear convention |
| Length | converted to mm |
| Strain | `strain_unit` to dimensionless; `engineering` shear strain is halved to the tensor convention (only for strain, never for stress) |
| Frame | global Cartesian; a non-global tensor is the exporter's job to rotate (refused by the header otherwise) |
| Averaging | recorded; `nodal_average_global` with more than one material or instance is refused unless `allow_global_average=True` (averaged values at an interface mix materials) |

### Strain composition and derivation

| File has | Result | Recorded as |
|---|---|---|
| `E_*` | total mechanical strain | `total_mechanical` |
| `Eel_*` and `Epl_*` | both stored; total is the sum | `total_mechanical` |
| `Eel_*` only, `linear_elastic` | total strain (no plasticity) | `total_mechanical` |
| `Eel_*` only, `elastic_plastic` | **refused**: the elastic strain is not the total under plasticity | |
| `Epl_*` only | needs `elastic_constants={material: (E, nu)}`; the elastic part is derived by Hooke's law from the stress | `derived_elastic_part` |
| nothing, `linear_elastic` | stress-only file; with `elastic_constants` the Hooke strain is derived | none / `derived_linear_elastic` |
| nothing, `elastic_plastic` | **refused**: elastic strain cannot be derived under plasticity | |

`E` is in MPa; the dictionary key is the `material` column value, or the empty string for a file
without one. A derived strain is flagged in the file's metadata so that nobody mistakes it for solver
output. Criteria that need strain use the total mechanical strain.

### Steps and what they mean

* `step_role: time_sample`: each node's series over the steps is a `TimeHistory` (axis: `time` if the
  file has it, else the solver step numbers as a cycle axis). Two steps, maximum and minimum, give
  exactly one cycle under the default residue policy; a sampled cycle goes through rainflow counting.
  Proportional node histories are exact; a non-proportional one is refused unless
  `allow_nonproportional_approximation=True`, which is recorded.
* `step_role: unit_case`: each step is a unit load case `S_k`; a **load-function** CSV with one axis
  column (`time` or `index`) and one column `f_<step>` per unit case, matched to the steps **by name**
  (a missing or extra function column is refused), gives `S(t) = sum_k f_k(t) S_k` per node and the
  same for the strain. Linear elastic only.

```python
rows = assess_nodes("unit_cases.h5", setup, R_list=R_list, load_functions="functions.csv")
```

### Refusals in a node run

`assess_nodes` records a `status` per node and continues (`on_refusal="record"`, the default):
`ok`, `ineligible: <reason>` (for example a non-proportional history, a native-life method, missing
strain) or `no_material_setup` (the node's material has no entry in the `materials` mapping).
`on_refusal="raise"` stops at the first. Unexpected exceptions are never caught. `MaterialSetup(model,
params, curve, E=None, nu=None)` carries the criterion, its constants and the Basquin curve of a
material; with `E` and `nu` a stress-only linear-elastic file gets its Hooke strain at assessment
(flagged in the history provenance).

## Units

Units are explicit tables at the import edge (stress, length and strain, above), not `pint` or
`xarray`: those do not check a declared unit against the data, and unit-aware arrays would add heavy
dependencies for a closed problem. An optional labelled `xarray` view of a history (dimensions
location, step, component) is a possible later extra; it would not replace the HDF5 store.

## Size

Export the nodes you need (a node set at the critical region, the surface) rather than the whole
model, and write step-major. A row is about 20 numbers, a few hundred bytes of text; ten million rows
is a few gigabytes, which the importer streams. For models beyond that, write the HDF5 file directly (the
`HistoryWriter` is documented in [fe-input-format.md](fe-input-format.md)) or a binary exchange format
later; Parquet is a candidate and is not implemented.

## Solver checklists

These are **checklists, not tested code**. They say what must hold for the exported numbers to be
what the header declares; the importer cannot see inside the solver. Tested export scripts from users
of each solver are requested (see the questions below), and a pure-shear check (below) is the
solver-independent way to verify a shear convention.

**Abaqus**

1. Request nodal output of the stress `S` (and, for strain-based criteria, `E`, `EE` and `PE`)
   with the node label and the instance name; export the **undeformed** nodal coordinates.
2. Check that the components are in the **global** frame: a local orientation (`*ORIENTATION`) or a
   shell or beam section changes the frame of the output. If it does, transform before writing.
3. Choose the averaging on purpose and write it in `averaging`: unaveraged element-nodal values,
   nodal values averaged within each material (preferred for components with several materials), or
   averaged across everything. Averaging across a material interface mixes materials.
4. Write one row per node and step; for an elastic-plastic analysis export the **stabilised** cycle
   and say `cycle: stabilised`. Write the time of each step if the load history has a time scale.
5. Declare `strain_convention` after checking your version's documentation and the pure-shear check.

**Ansys (Mechanical APDL or Mechanical)**

1. Export the nodal component stress `SX ... SXZ` (and `EPEL`, `EPPL`, `EPTO` for strain) with the
   node number; use the **global Cartesian** result coordinate system (`RSYS,0`), not a local or
   cylindrical one.
2. Check the averaging: nodal averaging across bodies of different materials mixes them; select the
   nodes by material or export per body and set `material` accordingly.
3. Component strains in Ansys are documented with engineering shear; check your version and
   declare `strain_convention` accordingly.
4. Export undeformed coordinates, one row per node and result set; for plasticity use the stabilised
   result set and say so in `cycle`.

**Pure-shear check of a convention.** Run a one-element model under simple shear `gamma = 0.002` and
read the shear strain component the solver reports: `0.002` means `engineering`, `0.001` means
`tensor`. The same check applies to the writer you use.

## Writer helper

`multifat.fe.csvwriter.FeCsvWriter` uses only the standard library, so a solver script can copy the
single file. It writes the header (`format` first), the column line and the rows with 17 significant
digits, and checks that the required header keys are present, that every row has exactly the declared
columns and that numbers are finite. It does not validate units or conventions; the importer does.

```python
from multifat.fe.csvwriter import FeCsvWriter

header = {"format": "multifat-fe-csv 1.0", "solver": "MyFE", "solver_version": "2.1",
          "stress_unit": "MPa", "length_unit": "mm", "coordinate_system": "global_cartesian",
          "averaging": "nodal_average_within_material", "analysis_type": "linear_elastic",
          "step_role": "time_sample"}
columns = ["node_id", "material", "step", "x", "y", "z", "S_xx", "S_yy", "S_zz", "S_xy", "S_yz", "S_xz"]
with FeCsvWriter("out.csv", header, columns=columns) as writer:
    writer.write_row({"node_id": 1, "material": "steel", "step": 1, "x": 0.0, "y": 0.0, "z": 0.0,
                      "S_xx": 100.0, "S_yy": 0.0, "S_zz": 0.0, "S_xy": 0.0, "S_yz": 0.0, "S_xz": 0.0})
```

The synthetic generator `multifat.fe.synthetic.write_tube_csv` produces the analytical tube in the
canonical, Abaqus and Ansys column names, with engineering or tensor shear, other units, either row
order and optional instance, material and time: it is the test oracle of every convention above.

## Solver importers

The Python adapters in `multifat.fe.solvers` write this CSV, never HDF5 directly. Each returns
the validated `CsvHeader` and replaces the destination atomically only after every frame passes.
An existing output survives a refusal. No solver input is executed.

```python
from multifat.fe.solvers.abaqus import abaqus_report_to_csv
from multifat.fe.solvers.ansys import ansys_listing_to_csv

declarations = dict(stress_unit="MPa", length_unit="mm", analysis_type="linear_elastic",
                    step_role="time_sample", coordinate_system="global_cartesian", time_unit="s")
abaqus_report_to_csv(["frame1.rpt", "frame2.rpt"], "abaqus.csv", inp="flat_nodes.inp",
                     material="steel", strain_convention="tensor", strain_unit="1", **declarations)
ansys_listing_to_csv("stress.out", "ansys.csv", nodes="global_nodes.lis", material="steel",
                     averaging="nodal_average_within_material", **declarations)
```

Then call `import_fe_csv` and `assess_nodes` as above. Units, analysis type, step role and global
coordinate system are required declarations; any strain requires `strain_unit` and
`strain_convention`, any time requires `time_unit`, and elastic-plastic text requires
`cycle="stabilised"`. Choose these from the solver setup; the example values are not defaults.
`solver_version` defaults to `unknown` for text, or can be supplied explicitly.

**Abaqus RPT.** Read complete Field Output Report tables with `Node Label`, optional `Element Label`,
six `S.S11 ... S.S23` components and optional complete `E.`, `EE.` or `PE.` groups. The original
step, increment, Step Time/Total Time, ODB name, instance, averaging criteria and region-boundary
text are retained. Components reuse `COLUMN_ALIASES`; `S.S13` is xz and `S.S23` is yz.
Coordinates come from flat explicit global `*NODE` INP data or report `COORD.COOR1/2/3` columns
asserted to be undeformed. Part/assembly transforms, generated nodes and `*INCLUDE` are refused.
Supply a material-specific selection and its `material` label; no INP section/material mapping
is inferred. Conflicting material declarations, local output and shell section points are refused.

The report must state Unique Nodal or Element Nodal. Element-node duplicates collapse only if
all values agree (`rel_tol=1e-10`, `abs_tol=1e-14` in source units); no averaging is performed.
Otherwise export Unique Nodal within material or separate labelled element/material selections.
Region-boundary averaging records `nodal_average_within_material`; unqualified nodal averaging
records `nodal_average_global`; unaveraged output needs an explicit `averaging` argument.
Conflicting arguments are refused. `.dat` and native `.odb` are unsupported; ODB needs Abaqus.

**Ansys text.** Read `PRNSOL,S,COMP` tables with `NODE SX SY SZ SXY SYZ SXZ`, optionally accompanied
by disjoint complete `EPEL*`, `EPPL*`, `EPTO*` strain tables for the same state. Coordinates come
from global Cartesian NLIST (`DSYS=0`) or table columns. Load step, substep, time and coordinate/
averaging text are recorded; `averaging` must always be supplied. Use `RSYS,0` and export separately
per material. Shell section points, element tables, repeated node/component values and fixed-width
numbers that run into each other are refused. Complete tab/comma-delimited Workbench tensor tables
are supported with `Node Number`, X/Y/Z and all six named components. Bracketed units must match
the arguments; provide `step` and optionally `time` when absent. Scalar equivalent-stress or
single-component exports cannot define a tensor. This is a strict subset of Mechanical layouts,
tested only on hand-written text samples, not a claim to support all Workbench exports.

For both text adapters the complete column set, location grid and undeformed coordinates must
remain consistent; finite values, ascending states and ascending time are required. CSV `step`
is a result-state ordinal starting at 1; solver identities are recorded in `notes`. Step-time
resets are refused: export total time instead. No time offset is guessed. Global averaging across
several materials or instances is refused. Reports stream by frame with disk spooling before
publication. No `.dat` parser or command-line entry point is provided.

**Native Ansys RST.** Install `uv sync --extra ansys --extra fe` and call
`multifat.fe.solvers.rst.ansys_rst_to_csv("solid.rst", "solid.csv", material="steel", **declarations)`.
The MIT optional reader is bounded to the inspected 0.56 API; it installed/imported successfully
with Python 3.12 and numpy 2.5.4, without downgrading project numpy. `nodal_stress` returns
X,Y,Z,XY,YZ,XZ after global rotation and nodal averaging. Source inspection confirms that
`read_element_result(as_global=1)` rotates ENS using EUL before accumulation. EEL is not rotated
the same way, so native **strain is not exported**. The caller must declare linear elasticity.
Only default-option SOLID185/186/187, one material, real static/transient states and small
displacements are supported. Modal, complex, cyclic, shell/beam, nondefault/layered elements and
multiple materials are refused. Averaging is `nodal_average_global`; no interface splitting is
available. SI units are checked against Pa/m; unsupported known systems are refused. User Defined
units require caller declarations. Missing/non-finite stress is refused, including midside nodes
without stored stress: pass `nodes=[...]` explicitly to select corner nodes; none are silently dropped.

**Evidence, 6 November 2026.** Fetch-free tests exercise refusals, units/conventions, separate and
concatenated frames, atomic output and distinct component sentinels. A process-local alias swap
of S13/S23 fails the sentinel test; production files remain intact.

| Adapter | Real output checked | Evidence and limits |
|---|---|---|
| Abaqus | Synthetic in-memory reports and parser unit tests | Hand-written sentinels, multi-frame composition, element-nodal agreement and all refusals verified. Real Abaqus output generator marked as an expected failure TODO pending Abaqus license/installation on workstation. No external repository data is included. |
| Ansys text | In-house 2D plane stress holed-plate biaxial fatigue benchmark (`hole_plate_fatigue.out`) | Self-generated with Ansys MAPDL 2026 R1 batch mode on WorkstationUniTN using `PLANE183` elements under a 2-step alternating fatigue cycle ($Y$-tension 100 MPa, then $X$-tension 100 MPa, rotating principal axes by $90^\circ$). Verified against Howland finite-width analytical baseline ($K_{tg} \approx 3.032$). See [benchmark record](../tests/data/solvers/ansys/README.md). |
| Native RST | In-house 2D plane stress holed-plate biaxial fatigue benchmark (`hole_plate_fatigue.rst`) | Self-generated with Ansys MAPDL 2026 R1 batch mode on WorkstationUniTN. Converted into exchange CSV and HDF5 with 2,864 corner nodes across 2 load steps. Matches text `PRNSOL` values to within text printing precision. Optional `ansys` extra verified on Python 3.12 / NumPy 2.5.4. |

Licenses were checked before use; each downloaded file went to a new empty scratch directory,
was inspected as text where appropriate, and inspection/conversion Python ran with `-I`.
No downloaded script or input was executed. Native repository commit
`90dd3a7ba54c8a0efb408ff59956dbef6f8b1f44`:
[MIT license](https://github.com/ansys/pymapdl-reader/blob/90dd3a7ba54c8a0efb408ff59956dbef6f8b1f44/LICENSE),
[RST source](https://github.com/ansys/pymapdl-reader/blob/90dd3a7ba54c8a0efb408ff59956dbef6f8b1f44/tests/testfiles/vol_test.rst),
[rotation code](https://github.com/ansys/pymapdl-reader/blob/v0.56.0/ansys/mapdl/reader/cython/_binary_reader.pyx).
RST SHA-256: `5e261e7dd46a8d5f28a8cde8e90acb2eaa8e89e1c91f60053229dc7bfa735727`.
Obtain it separately after license inspection, set `MULTIFAT_TEST_RST` to its absolute path and
run `uv run --extra ansys --extra fe python -I -m pytest tests/test_solver_rst.py -q`.
Without that path the ordinary suite skips the single real-RST test; it never downloads test data.
Format references: [Abaqus report positions](https://docs.software.vt.edu/abaqusv2025/English/SIMACAECAERefMap/simacae-t-reportfieldsetuptabbtn.htm),
[Ansys PRNSOL components](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_cmd/Hlp_C_PRNSOL.html).
These checks establish software fidelity; nothing here validates a fatigue criterion.

## Verification

Round trips against the analytical tube for every flavour, shear convention, unit (`MPa`/`mm`, `Pa`/`m`,
`psi`/`in`, `ksi`/`ft`, `N/mm^2`/`cm`), row order and strain mode, including a distinct value for each of
the six components; sentinel tests of the Abaqus and Ansys component mapping; refusal tests (truncated
and ragged grids, duplicate keys, unsorted steps, time not increasing, non-finite and blank cells,
header omissions, a non-global frame, a missing stabilised cycle, global averaging across materials,
unit cases on elastic-plastic data); interface nodes and repeated node labels; derived strain against
Hooke's law by hand; a tracemalloc test that the memory of the import does not grow with the number of
steps; and the node assessment against `assess_history` on the same histories. Two deliberate mutations
(swapping `S13` and `S23`, and halving the stress shear) were tried and fail the tests.

## Limits

* Local (single-point) data only; the step-to-cycle choice is the user's: two states, a sampled cycle or
  unit cases with load functions. A non-proportional sampled cycle is not reduced to harmonics yet.
* Global Cartesian components only; the exporter rotates.
* Undeformed coordinates; geometric nonlinearity is not followed.
* Native Abaqus `.odb` needs Abaqus and is outside scope; optional Ansys RST conversion writes CSV.
* Nothing here validates a fatigue criterion.

## Questions for the maintainer

1. Is `material` a workbook material ID (so a later resolver can fill the criterion constants)?
2. What should `instance` be for Ansys (a body, a part, an element-type group)?
3. Is `step` always in time order in your exports, and does it restart per load case?
4. Can you share one real Abaqus and one real Ansys export (a small model), and the script that wrote
   them? The checklists above should be replaced by tested scripts.

Decisions: [decision-log.md](decision-log.md), `d086` to `d098`.
