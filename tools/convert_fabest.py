"""Convert the FABEST database (42CrMo4 +QT) into the multifat workbook 42CrMo4_FABEST.

Source: The FABEST database: A benchmark fatigue dataset on 42CrMo4 +QT steel, version
FAB_v20261014, Zenodo, doi:10.5281/zenodo.23359038, CC BY 4.0. The original file is kept unchanged
in Experimental_Data/sources/FABEST/ and checked against SHA256 below before conversion.

Scope (maintainer decision, 2026-11-09): only the PBS series of plain hollow hourglass specimens
(outer/inner diameter 8.1/7 mm, turned), which holds every torsion and multiaxial curve of the
database, so that the design set and the assessed tests share lab, geometry and finish:
axial FAB060, FAB070; torsion FAB062, FAB072; tension-torsion FAB064, FAB074 (0 deg) and FAB065,
FAB075 (90 deg). All stresses are R = -1, as reported (nominal, net section).

Reading rules:
  - one row per test point of a curve sheet; the "<curve>_Rm" row (static strength, N < 1) is
    skipped;
  - "finished" is a failure with N = "till break"; "interrupted" is a run-out with N = "till END"
    (Observation type "run-out": an interrupted test of an S-N curve, decision d113);
  - channel 1 = axial (tension-compression) or torsion, channel 2 = torsion; the channel-2 phase
    shift is the torsional phase relative to the axial load.

Summary (maintainer decision): E, nu, Su and the lower yield stress ReL are the averages of the
four PBS tensile tests; Sigm1 and Taum1 are the authors' Basquin strengths at N = 1e6 (SN_FL,
SN_NFL) of the axial and torsion curves, averaged per load type (a finite-life strength, not a
staircase limit); sigma_f, b_fat, epsilon_f, c_fat are the strain-life constants of sheet FAB501.

Usage: uv run python tools/convert_fabest.py [source.xlsm] [output.xlsx]
"""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Experimental_Data" / "sources" / "FABEST" / "FAB_v20261014.xlsm"
OUTPUT = ROOT / "Experimental_Data" / "16_42CrMo4_FABEST.xlsx"
SHA256 = "720f637f5a1c74f579c62197137db54d0d2a495d6efdc607ddc4548a7d4cefff"
DOI = "10.5281/zenodo.23359038"
REFERENCE = f"FABEST database FAB_v20261014, doi:{DOI}"
AUTHORS = "FABEST consortium (COST CA23109 FABER)"
YEAR = 2026
MATERIAL = "42CrMo4_FABEST"
AXIAL = ("FAB060", "FAB070")
TORSION = ("FAB062", "FAB072")
MULTIAXIAL = ("FAB064", "FAB065", "FAB074", "FAB075")
CURVES = AXIAL + TORSION + MULTIAXIAL
OD, ID = 8.1, 7.0  # mm, PBS hollow hourglass specimens (Master Index)

SUMMARY_HEADERS = [
    "Article DOI / Reference", "Authors", "Article year", "Material name",
    "Young's modulus (MPa)", "Poisson's ratio (-)", "Tensile strength (MPa)",
    "Yield strength (MPa)", "Tensile fatigue limit (MPa)", "Torsional fatigue limit (MPa)",
    "Number of Cycles (cycles)", "Stress ratio R", "Axial fatigue limit loading",
    "sigma_f (MPa)", "b_fat (-)", "epsilon_f (-)", "c_fat (-)",
]  # fmt: skip
DATASET_HEADERS = [
    "Authors", "Article year", "Article DOI / Reference", "Specimen name",
    "Number of cycles to failure", "Run out (Yes/No)", "Failure type (complete rupture/stiffness decrease)",
    "Geometry type (flat/cylindrical/hollow cylinder)", "Notched", "Diameter (mm)",
    "Inner diameter (mm)", "Curve ID",
    "Phase RADIAL load (°)", "Phase AXIAL load (°)", "Phase TORSIONAL load (°)",
    "Freq RADIAL load [Hz]", "Freq AXIAL load[Hz]", "Freq TORSIONAL load [Hz]",
    "Mean RADIAL stress  (MPa)", "Amplitute RADIAL stress  (MPa)",
    "Mean AXIAL stress  (MPa)", "Amplitute AXIAL stress  (MPa)", "derived from (AXIAL/BENDING) load",
    "Mean SHEAR stress  (MPa)", "Amplitute SHEAR stress  (MPa)",
    "Observation type (failure/run-out/fatigue limit)",
]  # fmt: skip


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _value(rows, key):
    for r in rows:
        if r and r[0] == key:
            return r[1]
    raise KeyError(key)


def curve_points(ws) -> list[dict]:
    """Test points of one curve sheet (failures and run-outs)."""
    rows = list(ws.iter_rows(values_only=True))
    two = rows[0][11] == "load channel No. 2"
    ch1 = rows[1][7]
    n_end, n_break, status = (15, 16, 17) if two else (11, 12, 13)
    out = []
    for r in rows[3:]:
        pid = r[5] if len(r) > 5 else None
        if not isinstance(pid, str) or not pid.startswith(ws.title + "_") or pid.endswith("_Rm"):
            if pid is None and out:
                break
            continue
        if not isinstance(r[7], (int, float)):
            continue
        finished = r[status] == "finished"
        if r[status] not in ("finished", "interrupted"):
            raise ValueError(f"{pid}: unknown status {r[status]!r}")
        N = r[n_break] if finished else r[n_end]
        if not (isinstance(N, (int, float)) and N >= 1):
            raise ValueError(f"{pid}: no cycle count")
        p = dict(id=pid, N=float(N), runout=not finished, freq=r[9])
        amp1, mean1 = float(r[7]), float(r[8] or 0.0)
        if ch1 == "tension-compression":
            p.update(sa=amp1, sm=mean1, ta=0.0, tm=0.0, phase=0.0)
        elif ch1 == "torsion":
            p.update(sa=0.0, sm=0.0, ta=amp1, tm=mean1, phase=0.0)
        else:
            raise ValueError(f"{ws.title}: channel 1 {ch1!r}")
        if two:
            if rows[1][11] != "torsion":
                raise ValueError(f"{ws.title}: channel 2 {rows[1][11]!r}")
            p.update(ta=float(r[11]), tm=float(r[12] or 0.0), phase=float(r[14] or 0.0))
            if r[13] != r[9]:
                raise ValueError(f"{pid}: channels at different frequencies")
        if p["sm"] or p["tm"]:
            raise ValueError(f"{pid}: mean stress outside the R = -1 scope")
        out.append(p)
    return out


def convert(source: Path = SOURCE, output: Path = OUTPUT) -> dict:
    if sha256(source) != SHA256:
        raise ValueError(f"{source}: SHA-256 differs from FAB_v20261014 ({SHA256})")
    wb = openpyxl.load_workbook(source, read_only=True, data_only=True, keep_vba=False)
    tens = list(wb["Tensile test"].iter_rows(values_only=True))
    avg = next(r for r in tens if r and r[0] == "averaged")
    ReL, Rm, E, nu = avg[2], avg[3], avg[5], avg[6]
    fl = {c: _value(list(wb[c].iter_rows(values_only=True)), "SN_FL") for c in AXIAL + TORSION}
    nfl = {c: _value(list(wb[c].iter_rows(values_only=True)), "SN_NFL") for c in AXIAL + TORSION}
    if set(nfl.values()) != {1000000}:
        raise ValueError(f"SN_NFL is not 1e6 for every curve: {nfl}")
    sigm1 = sum(fl[c] for c in AXIAL) / len(AXIAL)
    taum1 = sum(fl[c] for c in TORSION) / len(TORSION)
    sl = list(wb["FAB501"].iter_rows(values_only=True))[3]
    sigma_f, b_fat, eps_f, c_fat = sl[4], sl[5], sl[7], sl[8]

    out = openpyxl.Workbook()
    summary = out.active
    summary.title = "Summary"
    summary.append(SUMMARY_HEADERS)
    summary.append([REFERENCE, AUTHORS, YEAR, MATERIAL, E, nu, Rm, ReL, sigm1, taum1, 1e6, -1,
                    "push-pull", sigma_f, b_fat, eps_f, c_fat])  # fmt: skip
    summary.append([None] * 3 + ["Average", E, nu, Rm, ReL, sigm1, taum1, 1e6, -1, "push-pull",
                    sigma_f, b_fat, eps_f, c_fat])  # fmt: skip
    ds = out.create_sheet("Dataset")
    ds.append(DATASET_HEADERS)
    counts = {}
    for curve in CURVES:
        points = curve_points(wb[curve])
        counts[curve] = (sum(not p["runout"] for p in points), sum(p["runout"] for p in points))
        for p in points:
            ax = p["sa"] != 0.0
            sh = p["ta"] != 0.0
            ds.append([
                AUTHORS, YEAR, REFERENCE, p["id"], p["N"], "Yes" if p["runout"] else "No",
                "complete rupture", "hollow cylinder", "No", OD, ID, curve,
                "Null", 0 if ax else "Null", p["phase"] if sh else "Null",
                "Null", p["freq"] if ax else "Null", p["freq"] if sh else "Null",
                0, 0, p["sm"], p["sa"], "AXIAL", p["tm"], p["ta"],
                "run-out" if p["runout"] else "failure",  # interrupted test of an S-N curve (d113)
            ])  # fmt: skip
    out.save(output)
    return dict(
        Sigm1=sigm1,
        Taum1=taum1,
        fl=fl,
        E=E,
        nu=nu,
        Su=Rm,
        Sy=ReL,
        counts=counts,
        strain_life=(sigma_f, b_fat, eps_f, c_fat),
        rows=ds.max_row - 1,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", nargs="?", type=Path, default=SOURCE)
    parser.add_argument("output", nargs="?", type=Path, default=OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    info = convert(args.source, args.output)
    for k, v in info.items():
        print(f"{k}: {v}")
    assert math.isfinite(info["Sigm1"]) and math.isfinite(info["Taum1"])


if __name__ == "__main__":
    main()
