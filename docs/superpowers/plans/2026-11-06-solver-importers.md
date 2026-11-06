# Solver importer implementation plan

Authority: [approved handoff](2026-11-06-solver-importers-handoff.md). Execute inline on
`feat/fe-and-time-history-inputs`; local commits only. The handoff is the approved specification.
Progress and rulings live only in `docs/audit-checkpoint.md` (repository instruction overrides
the skill's separate ledger). No existing numerical source outside `fe`/`history` changes.

## Design decisions

Use `fe/solvers/{common,abaqus,ansys,rst}.py`. Shared code validates the existing CsvHeader,
ColumnMap and COLUMN_ALIASES, checks frame grids, and writes atomically using FeCsvWriter.
Abaqus parses frame-by-frame; Ansys text accepts complete six-component tables and NLIST
coordinates. Neither executes solver input. Reject unsupported output instead of guessing.
Require caller declaration of global Cartesian and undeformed coordinates when not established
by the file. Preserve source metadata in notes, including original step/frame identifiers.
CSV step is the ascending result-state ordinal; solver identifiers are recorded separately.
Do not synthesize total time from step time resets. Multiple source steps with resetting time
are refused unless the caller explicitly supplies a consistent time axis.
Native RST is optional and conditional on successful compatibility and public-data checks.
No CLI is needed: the Python functions are the thin public interface.

## Task 1: Abaqus reports

Files: shared adapter support, Abaqus parser, licensed small fixtures and provenance,
`tests/test_solver_abaqus.py`. Read public files in isolated scratch directories after license
inspection. Add failing round-trip tests (real axial report/INP and component sentinels), then
implement. Exercise concatenated frames, separate files, units/conventions, duplicate element
nodes, missing stress/coordinates/units, nonascending frames, changed columns and metadata.
Check targeted tests and notation; record results and commit.

## Task 2: Ansys text

Files: `ansys.py`, `tests/test_solver_ansys.py`, hand-written listings with documented source
format; inspect public APDL output for an applicable real table. Add failing PRNSOL/NLIST and
Workbench tests then implement, using the same aliases and CSV validation. Refuse missing
components, coordinate frame ambiguity, conflicting rows, changed columns and nonascending time.
Check targeted tests and notation; record results and commit.

## Task 3: Optional native Ansys

Probe Python 3.12/numpy 2 compatibility with `ansys-mapdl-reader` in an isolated environment;
never constrain project numpy downward. Inspect official reader source for global rotation and
component order, verify against a public licensed RST. If viable, add `rst.py`, optional `ansys`
extra and tests (real data as opt-in external path; fetch-free contract tests). Explicitly refuse
unsupported analyses/material averaging. Otherwise document the failed probe and defer the adapter.
Check targeted tests and notation; record results and commit.

## Task 4: Documentation and final gate

Update FE CSV guide, README, AGENTS, decision log d093 onward and checkpoint. Document tested
formats, licenses, actual real-file evidence, refusal limits, ODB exclusion and use examples.
Run all focused checks, full pytest once on a capped workstation job without concurrent edits,
black/ruff, documentation and scientific checks, generated-doc checks, and protected-source diff.
Run a distinct reviewer subagent on the whole branch diff. Fix important findings test-first,
rerun affected checks (and full suite after code fixes), checkpoint and commit locally. Do not push.

## Pre-flight interfaces and review focus

Tasks 1/2/3 consume shared canonical row validation and atomic CSV output; test component mapping,
time ordering, complete location grids, metadata consistency and no partial output on refusal.
Review pagination, multi-instance input, material interfaces, local frames, shell section points,
incomplete tables, unit declarations and native reader material averaging. Distinguish unsupported
formats from successful conversion. Nothing here validates fatigue criteria.
