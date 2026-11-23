"""Build the fatigue-limit workbook 20_34CrNiMo6_Pallares from Pallares-Santasmartas et al. (2018).

Sources (one steel, one supplier and heat treatment, tested by the same research group):

  [A] Pallares-Santasmartas L, Albizuri J, Aviles A, Aviles R. Mean stress effect on the axial
      fatigue strength of DIN 34CrNiMo6 quenched and tempered steel. Metals 2018;8(4):213,
      doi:10.3390/met8040213. Table 2 (monotonic properties) and Table 3 (axial fatigue strengths
      at 2e6 cycles for R = -2, -1, -0.5 and 0.05).
  [T] Pallares-Santasmartas L, Albizuri J, Aviles A, Saintier N, Merzeau J. Influence of mean
      shear stress on the torsional fatigue behaviour of 34CrNiMo6 steel. Int J Fatigue
      2018;113:54-68, doi:10.1016/j.ijfatigue.2018.04.008. Table 2 (the same monotonic
      properties) and Table 3 (torsional fatigue strengths at 2e6 cycles for a mean shear stress
      of 0, 150, 250, 350 and 500 MPa, staircase method of ISO 12107).

Both tables are transcribed in Experimental_Data/sources/Pallares2018/34crnimo6_limits.csv; this
script turns that transcription into the multifat workbook.

Scope (maintainer request of 2026-11-23, decision d120): smooth hourglass specimens, fatigue
strengths at 2e6 cycles, uniaxial loading only (axial with a mean normal stress, torsion with a
mean shear stress):

  - Summary: Su = 1210 and Sy = 1084 MPa (Table 2 of [A] and [T]); Sigm1 = 615 MPa, the push-pull
    limit at R = -1 of [A]; Taum1 = 432.5 MPa, the torsion limit at zero mean stress of [T]. [T]
    itself combines these two values (sigma_-1 = 615 MPa, Section 3 and Table 5). Neither paper
    gives E or nu: they are Null, so the models that need strains are skipped.
  - Dataset: one limit state per row (Run out = Yes at N = 2e6, Observation type = fatigue
    limit, as in the other fatigue-limit workbooks): the 4 axial and 5 torsional strengths.
    The pure R = -1 axial and torsional rows are the reference values Sigm1 and Taum1.

The repeated torsion limit tau_0 = 765 MPa of [T] Table 5 is not a tested load case (Table 3 has
no tau_m = tau_a test), and the R = 0.05 axial strength is not an R = 0 limit, so Sig0 and Tau0
are left empty. The stresses are the nominal (elastic) values of the sources; for tau_m = 250,
350 and 500 MPa the maximum shear stress exceeds the torsional yield strength (Sy/sqrt(3) =
626 MPa), as [T] states in Section 5.3.

Usage: uv run python tools/build_pallares_34crnimo6.py [transcription.csv] [output.xlsx]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Experimental_Data" / "sources" / "Pallares2018" / "34crnimo6_limits.csv"
OUTPUT = ROOT / "Experimental_Data" / "20_34CrNiMo6_Pallares.xlsx"
REFS = {
    "A": ("doi:10.3390/met8040213", "Pallares-Santasmartas L, Albizuri J, Aviles A, Aviles R.", 2018),
    "T": (
        "doi:10.1016/j.ijfatigue.2018.04.008",
        "Pallares-Santasmartas L, Albizuri J, Aviles A, Saintier N, Merzeau J.",
        2018,
    ),
}  # fmt: skip
MATERIAL = "34CrNiMo6_Pallares"
N_REF = 2e6  # cycle basis of every strength in [A] Table 3 and [T] Table 3

SUMMARY_HEADERS = [
    "Article DOI / Reference", "Authors", "Article year", "Material name",
    "Young's modulus (MPa)", "Poisson's ratio (-)", "Tensile strength (MPa)",
    "Yield strength (MPa)", "Tensile fatigue limit (MPa)", "Torsional fatigue limit (MPa)",
    "Number of Cycles (cycles)", "Stress ratio R", "Axial fatigue limit loading",
]  # fmt: skip
DATASET_HEADERS = [
    "Authors", "Article year", "Article DOI / Reference", "Specimen name",
    "Number of cycles to failure", "Run out (Yes/No)", "Failure type (complete rupture/stiffness decrease)",
    "Geometry type (flat/cylindrical/hollow cylinder)", "Notched", "Diameter (mm)",
    "Phase RADIAL load (°)", "Phase AXIAL load (°)", "Phase TORSIONAL load (°)",
    "Mean RADIAL stress  (MPa)", "Amplitute RADIAL stress  (MPa)",
    "Mean AXIAL stress  (MPa)", "Amplitute AXIAL stress  (MPa)", "derived from (AXIAL/BENDING) load",
    "Mean SHEAR stress  (MPa)", "Amplitute SHEAR stress  (MPa)",
    "Observation type (failure/run-out/fatigue limit)",
]  # fmt: skip
DIAMETER = {
    "A": "Null",
    "T": 10.0,
}  # [T] Fig. 8: 10 mm minimum diameter; [A] Fig. 4 is not transcribed


def _num(text: str) -> float | None:
    return float(text) if text.strip() else None


def read_transcription(path: Path = SOURCE) -> tuple[dict[str, float], list[dict]]:
    """Monotonic properties and fatigue strengths of the transcription."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    props = {r["quantity"]: float(r["value"]) for r in rows if r["kind"] == "property"}
    limits = [
        dict(
            source=r["source"],
            load=r["quantity"],
            R=_num(r["R"]),
            sm=_num(r["sigma_m_MPa"]) or 0.0,
            sa=_num(r["sigma_a_MPa"]) or 0.0,
            tm=_num(r["tau_m_MPa"]) or 0.0,
            ta=_num(r["tau_a_MPa"]) or 0.0,
            s=_num(r["s_MPa"]),
            n=_num(r["n_specimens"]),
        )  # fmt: skip
        for r in rows
        if r["kind"] == "limit"
    ]
    return props, limits


def stress_ratio(lim: dict) -> float:
    """R = (mean - amplitude) / (mean + amplitude) of the active component."""
    m, a = (lim["sm"], lim["sa"]) if lim["load"] == "axial" else (lim["tm"], lim["ta"])
    return (m - a) / (m + a)


def build(source: Path = SOURCE, output: Path = OUTPUT) -> dict:
    props, limits = read_transcription(source)
    for lim in limits:  # the printed R of [A] Table 3 must follow from the printed stresses
        if lim["R"] is not None and abs(stress_ratio(lim) - lim["R"]) > 0.01:
            raise ValueError(f"R = {stress_ratio(lim):.3f} from the stresses, {lim['R']} printed")
    sigm1 = next(lim["sa"] for lim in limits if lim["load"] == "axial" and lim["sm"] == 0.0)
    taum1 = next(lim["ta"] for lim in limits if lim["load"] == "torsion" and lim["tm"] == 0.0)

    out = openpyxl.Workbook()
    summary = out.active
    summary.title = "Summary"
    summary.append(SUMMARY_HEADERS)
    values = ["Null", "Null", props["Su_MPa"], props["Sy_MPa"], sigm1, taum1]
    for key in ("A", "T"):
        doi, authors, year = REFS[key]
        row_s, row_t = (sigm1, "Null") if key == "A" else ("Null", taum1)
        summary.append([doi, authors, year, MATERIAL, "Null", "Null", props["Su_MPa"],
                        props["Sy_MPa"], row_s, row_t, N_REF, -1, "push-pull" if key == "A" else None])  # fmt: skip
    summary.append([None] * 3 + ["Average", *values, N_REF, -1, "push-pull"])
    ds = out.create_sheet("Dataset")
    ds.append(DATASET_HEADERS)
    for k, lim in enumerate(limits):
        doi, authors, year = REFS[lim["source"]]
        axial = lim["load"] == "axial"
        mean = lim["sm"] if axial else lim["tm"]
        name = f"{'AXIAL' if axial else 'TORSION'}_MEAN_{mean:+.0f}MPA_LIMIT_{k + 1}"
        ds.append([
            authors, year, doi, name, N_REF, "Yes", "stiffness decrease", "cylindrical", "No",
            DIAMETER[lim["source"]],
            "Null", 0 if axial else "Null", "Null" if axial else 0,
            0, 0, lim["sm"], lim["sa"], "AXIAL", lim["tm"], lim["ta"],
            "fatigue limit",  # strength at 2e6 cycles (staircase in [T]) (d113)
        ])  # fmt: skip
    out.save(output)
    kappa = sigm1 / taum1
    return dict(Sigm1=sigm1, Taum1=taum1, kappa=kappa, rows=ds.max_row - 1,
                c_T=sigm1 / (3.0 - math.sqrt(3.0) * kappa))  # fmt: skip


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", nargs="?", type=Path, default=SOURCE)
    parser.add_argument("output", nargs="?", type=Path, default=OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for k, v in build(args.source, args.output).items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
