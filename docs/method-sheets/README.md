# Method sheets (LaTeX)

One LaTeX document that states every registered criterion as its paper defines it, for reuse in
the FABER WG4.3 review article.

- `methods.tex`: the master file (preamble, shared macros, `\input` of each sheet).
- `notation.tex`: the notation shared by all sheets (frame, harmonic loading, invariants, plane
  quantities, material constants, master curve, load-equivalent error, plane search).
- `sheets/<NAME>.tex`: one sheet per method. The first line is `% registry: NAME` (several names
  separated by commas when one sheet covers variants). Every sheet has the blocks Source, Idea and
  scope, Inputs, Formulation, Constants and calibration, Life, Repository choices, Worked example
  and Verification, in that order.
- `references.bib`: the original papers.

Build:

```bash
cd docs/method-sheets
latexmk -pdf -interaction=nonstopmode -halt-on-error methods.tex
```

`tests/test_method_sheets.py` checks that every registered model has a sheet or is listed as
pending, that every sheet is input by `methods.tex`, that the blocks are present and in order,
that every `\cite` key is in `references.bib`, and that braces and environments balance. Every
registered model has a sheet; the pending list is empty, so a new method arrives with its sheet.

The worked example of a sheet repeats numbers from the method guide in `docs/methods/`, which the
tests execute.
