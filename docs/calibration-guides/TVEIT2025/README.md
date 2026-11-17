# TVEIT2025 calibration guide

A step-by-step guide to calibrating TVEIT2025 (Tveit et al. 2025) for a new material: the data it
needs, the six steps of the authors' procedure with the rules this repository adds, how to read the
calibration confidence, and worked examples. The examples are the calibrations of this repository:

- Al7075-T6;
- Ti6Al4V, two batches (Berto, Meneghetti);
- 42CrMo4, two batches (Santus, Reis);
- Al7075-T651 (Wang).

The method guide [`TVEIT2025.md`](../../method-guides/TVEIT2025.md) is the reference for the
procedure; this guide explains it with examples.

| File | Content |
|---|---|
| [`tveit2025_calibration_guide.pdf`](tveit2025_calibration_guide.pdf) | The guide, compiled: a reading copy of the LaTeX source. |
| [`tveit2025_calibration_guide.tex`](tveit2025_calibration_guide.tex) | The guide (LaTeX). |
| [`figures/`](figures) | The figures, in PDF. |
| [`make_figures.py`](make_figures.py) | Draws the figures from the calibration records. |

The LaTeX source is authoritative; rebuild the PDF when the source changes. Build the guide
(pdflatex, twice, for the references):

```bash
cd docs/calibration-guides/TVEIT2025
pdflatex tveit2025_calibration_guide.tex
pdflatex tveit2025_calibration_guide.tex
```

Redraw the figures from the repository root:

```bash
python docs/calibration-guides/TVEIT2025/make_figures.py
```

The script uses the pipeline's calibration: it reuses the record committed in
`Calibrations/<material>/` and calibrates a workbook without one (Al7075_T651_Wang, about 20 s).
The numbers quoted in the guide come from the same records; recalibrated records require the guide
to be updated.
