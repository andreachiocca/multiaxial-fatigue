"""Build the fatigue-limit workbook 17_Ti6Al4V_Delahay from Delahay & Palin-Luc (2006).

Source: Delahay T, Palin-Luc T. Estimation of the fatigue strength distribution in high-cycle
multiaxial fatigue taking into account the stress-strain gradient effect. Int J Fatigue
2006;28(5-6):474-484, doi:10.1016/j.ijfatigue.2005.06.048. Table 1 (p. 478) and the Ti-6Al-4V
rows of Table 2 (p. 479) are transcribed in Experimental_Data/sources/Delahay2006/
ti6al4v_tables.csv; this script turns that transcription into the multifat workbook.

Scope (maintainer decision d109, 2026-11-11): annealed Ti-6Al-4V, smooth specimens, median
staircase endurance limits (P_f = 0.5, >= 15 specimens, 1e6 or 1e7 cycles), all R = -1:

  - Summary: E, nu, Su, Sy of Table 1; Sigm1 = the plane-bending limit 652 MPa (loading recorded
    as "plane bending", FatPy fat_lim_ben_m1) and Taum1 = the torsion limit 411 MPa. The assessed
    points are plane bending + torsion, so calibrating on the same load types lets the stress
    gradient largely cancel in a pointwise criterion. The tension (583 MPa) and rotating-bending
    (602 MPa) limits stay in the transcription as information; they are not in the workbook.
  - Dataset: one limit state per row (Run out = Yes at N = 1e6, the convention of the other
    fatigue-limit workbooks; Observation type = fatigue limit, decision d113): plane bending 652, torsion 411, plane bending + torsion 442/255
    in phase and 567/328 at 90 deg. The bending stress is the nominal surface amplitude, entered
    in the axial column with "derived from" = BENDING.

The paper does not state the specimen diameter of each test type (the size-effect figures name
7.98 mm for the reference specimens); the diameter is left Null.

Usage: uv run python tools/build_delahay_ti6al4v.py [transcription.csv] [output.xlsx]
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Experimental_Data" / "sources" / "Delahay2006" / "ti6al4v_tables.csv"
OUTPUT = ROOT / "Experimental_Data" / "17_Ti6Al4V_Delahay.xlsx"
DOI = "10.1016/j.ijfatigue.2005.06.048"
AUTHORS = "Delahay T, Palin-Luc T."
YEAR = 2006
MATERIAL = "Ti6Al4V_Delahay"
N_REF = 1e6  # cycle basis of the limits (Table 2: "1e6 or 1e7 cycles")
DESIGN = {"Sigm1": "plane bending", "Taum1": "torsion"}  # d109
ASSESSED = ("plane bending", "torsion", "plane bending + torsion")

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


def read_transcription(path: Path = SOURCE) -> tuple[dict[str, float], list[dict]]:
    """Material properties (Table 1) and endurance limits (Table 2) of the transcription."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    props = {r["quantity"]: float(r["value"]) for r in rows if r["kind"] == "property"}
    limits = []
    for r in rows:
        if r["kind"] != "limit":
            continue
        limits.append(
            dict(
                load=r["load"],
                sa=float(r["sigma_a_MPa"] or 0.0),
                ta=float(r["tau_a_MPa"] or 0.0),
                s=float(r["s_MPa"]),
                phase=float(r["phase_deg"] or 0.0),
            )
        )
    return props, limits


def threshold_stress(limits: list[dict]) -> float:
    """sigma* of Eq. (6) from the tension and rotating-bending limits (a transcription check)."""
    by = {lim["load"]: lim["sa"] for lim in limits}
    return math.sqrt(2.0 * by["tension"] ** 2 - by["rotating bending"] ** 2)


def build(source: Path = SOURCE, output: Path = OUTPUT) -> dict:
    props, limits = read_transcription(source)
    sigma_star = threshold_stress(limits)
    if abs(sigma_star - props["sigma_star_MPa"]) > 1.0:  # Table 1 gives 563 MPa
        raise ValueError(
            f"Eq. (6) gives sigma* = {sigma_star:.1f}, Table 1 {props['sigma_star_MPa']}"
        )
    single = {lim["load"]: lim for lim in limits if lim["load"] in DESIGN.values()}
    sigm1, taum1 = single[DESIGN["Sigm1"]]["sa"], single[DESIGN["Taum1"]]["ta"]

    out = openpyxl.Workbook()
    summary = out.active
    summary.title = "Summary"
    summary.append(SUMMARY_HEADERS)
    values = [props["E_MPa"], props["nu"], props["Su_MPa"], props["Sy_MPa"], sigm1, taum1]
    summary.append([f"doi:{DOI}", AUTHORS, YEAR, MATERIAL, *values, N_REF, -1, DESIGN["Sigm1"]])
    summary.append([None] * 3 + ["Average", *values, N_REF, -1, DESIGN["Sigm1"]])
    ds = out.create_sheet("Dataset")
    ds.append(DATASET_HEADERS)
    for k, lim in enumerate(lim for lim in limits if lim["load"] in ASSESSED):
        bending, shear = lim["sa"] != 0.0, lim["ta"] != 0.0
        name = lim["load"].replace(" + ", "_").replace(" ", "_").upper()
        if bending and shear:
            name += f"_{lim['phase']:.0f}DEG"
        ds.append([
            AUTHORS, YEAR, f"doi:{DOI}", f"{name}_LIMIT_{k + 1}", N_REF, "Yes",
            "Null", "cylindrical", "No", "Null",
            "Null", 0 if bending else "Null", lim["phase"] if shear else "Null",
            0, 0, 0, lim["sa"], "BENDING" if bending else "AXIAL", 0, lim["ta"],
            "fatigue limit",  # median staircase limit (d113)
        ])  # fmt: skip
    out.save(output)
    return dict(Sigm1=sigm1, Taum1=taum1, sigma_star=sigma_star, rows=ds.max_row - 1)


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
