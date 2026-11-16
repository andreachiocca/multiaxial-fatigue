"""Build the fatigue-life workbook 18_Al7075_T651_Wang from Wang, Shang, Sun & Chen (2018).

Source: Wang X-W, Shang D-G, Sun Y-J, Chen H. Multiaxial high-cycle fatigue life prediction
model considering mean shear stress effect under constant and variable amplitude loading.
Theor Appl Fract Mech (2018), doi:10.1016/j.tafmec.2017.10.007 (accepted manuscript). Tables 3-5
(52 constant-amplitude tests on 7075-T651 thin-walled tubes, load controlled, sine 10 Hz) are
transcribed in tests/data/papers/wang2018_tables3_5.csv; this script turns that transcription,
with Tables 1 and 2, into the multifat workbook (maintainer decision, 16 November 2026).

  - Summary: E = 71.7 GPa and G = 27.5 GPa (Table 1), so nu = E/(2G) - 1 = 0.3036; Sy = 501 and
    Su = 561 MPa. Sigm1 and Taum1 are the fully reversed axial and torsional S-N curves of Table 2
    at N = 1e6 cycles (the first of the paper's N_0 = 1e6, 2e6, 1e7): 888.90 N^-0.10042 = 222.0
    and 688.85 N^-0.12233 = 127.1 MPa (axial loading push-pull). No R = 0 limits are given.
  - Dataset: one row per test, all failures (Run out = No, Observation type = failure). sigma is
    the axial stress, tau the shear stress; tau lags sigma by delta (Eq. 8), entered as the
    torsional phase. The radial stress is 0. Thin-walled tube, outer diameter 12.5 mm (Fig. 1).

Usage: uv run python tools/build_wang_al7075_t651.py [transcription.csv] [output.xlsx]
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests" / "data" / "papers" / "wang2018_tables3_5.csv"
OUTPUT = ROOT / "Experimental_Data" / "18_Al7075_T651_Wang.xlsx"
DOI = "10.1016/j.tafmec.2017.10.007"
AUTHORS = "Wang X-W, Shang D-G, Sun Y-J, Chen H."
YEAR = 2018
MATERIAL = "Al7075_T651_Wang"
N_REF = 1e6
E_MPA, G_MPA, SY, SU = 71700.0, 27500.0, 501.0, 561.0  # Table 1
SIGM1 = 888.90 * N_REF**-0.10042  # Table 2, axial
TAUM1 = 688.85 * N_REF**-0.12233  # Table 2, torsional

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


def read_transcription(path: Path = SOURCE) -> list[dict]:
    """The 52 tests of Tables 3-5."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(line for line in fh if not line.startswith("#")))
    return [{k: (v if k == "table" else float(v)) for k, v in r.items()} for r in rows]


def build(source: Path = SOURCE, output: Path = OUTPUT) -> dict:
    tests = read_transcription(source)
    nu = E_MPA / (2.0 * G_MPA) - 1.0
    out = openpyxl.Workbook()
    summary = out.active
    summary.title = "Summary"
    summary.append(SUMMARY_HEADERS)
    values = [E_MPA, round(nu, 4), SU, SY, round(SIGM1, 2), round(TAUM1, 2)]
    summary.append([f"doi:{DOI}", AUTHORS, YEAR, MATERIAL, *values, N_REF, -1, "push-pull"])
    summary.append([None] * 3 + ["Average", *values, N_REF, -1, "push-pull"])
    ds = out.create_sheet("Dataset")
    ds.append(DATASET_HEADERS)
    for t in tests:
        axial, shear = t["sigma_a"] != 0 or t["sigma_m"] != 0, t["tau_a"] != 0 or t["tau_m"] != 0
        kind = "AXIAL_TORSION" if axial and shear else ("AXIAL" if axial else "TORSION")
        ds.append([
            AUTHORS, YEAR, f"doi:{DOI}", f"{kind}_{int(t['specimen'])}", int(t["N_f"]), "No",
            "complete rupture", "hollow cylinder", "No", 12.5,
            "Null", 0 if axial else "Null", t["delta"] if shear else "Null",
            0, 0, t["sigma_m"], t["sigma_a"], "AXIAL", t["tau_m"], t["tau_a"],
            "failure",
        ])  # fmt: skip
    out.save(output)
    return dict(Sigm1=round(SIGM1, 2), Taum1=round(TAUM1, 2), nu=round(nu, 4), rows=ds.max_row - 1)


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
