"""Abaqus report conversion: real axial output, mapping sentinels and fail-closed parsing."""

from pathlib import Path
import csv

import numpy as np
import pytest

from multifat.fe.solvers.abaqus import abaqus_report_to_csv
from multifat.fe.csvio import import_fe_csv
from multifat.fe.store import HistoryFile

DATA = Path(__file__).parent / "data/solvers/axis-fem"
REPORT = DATA / "problema_1_strain_stress.rpt"
INP = DATA / "problema_1.inp"
ABAQUS_TODO = pytest.mark.xfail(
    not REPORT.exists(),
    reason="Real Abaqus output generator pending Abaqus installation/license on workstation (TODO)",
)
OPTIONS = dict(
    stress_unit="Pa",
    length_unit="m",
    analysis_type="linear_elastic",
    step_role="time_sample",
    coordinate_system="global_cartesian",
    strain_convention="tensor",
    strain_unit="1",
    time_unit="s",
)


def rows(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def report(increment=1, time=1, values="11 22 33 12 13 23", extra="", position="Unique Nodal"):
    return f"""Field Output Report
ODB: test.odb
Step: Step-1
Frame: Increment {increment}: Step Time = {time}
{position} values (Not averaged across region boundaries)
Field Output reported at nodes for part: PART-1-1
{extra}
Node Label S.S11 S.S22 S.S33 S.S12 S.S13 S.S23 COORD.COOR1 COORD.COOR2 COORD.COOR3
Label @Loc 1
-------------------
1 {values} 0 0 0
"""


def convert(tmp_path, text, **options):
    path = tmp_path / "input.rpt"
    path.write_text(text, encoding="utf-8")
    return abaqus_report_to_csv(path, tmp_path / "out.csv", **(OPTIONS | options))


@ABAQUS_TODO
def test_real_report_roundtrip(tmp_path):
    if not REPORT.exists():
        pytest.fail("Abaqus solver not available on workstation; real output generator is a TODO")
    header = abaqus_report_to_csv(REPORT, tmp_path / "out.csv", inp=INP, **OPTIONS)
    assert "75%" in header.notes and "Step-1" in header.notes
    assert "problema_1.odb" in header.source
    assert header.averaging == "nodal_average_within_material"
    got = rows(tmp_path / "out.csv")
    assert len(got) == 8
    assert float(got[0]["S_yz"]) == 8.06554e-12
    assert float(got[0]["S_xz"]) == 2.97323e-12
    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        np.testing.assert_allclose(block.S_hist[:, 0, 2, 2], 0.05)
        np.testing.assert_allclose(block.S_hist[:, 0, 0, 0], 0, atol=2e-17)
        np.testing.assert_allclose(block.E_hist[:, 0, 2, 2], 250e-9)


def test_component_sentinels_and_frames(tmp_path):
    header = convert(tmp_path, report() + report(2, 2))
    assert header.n_steps == 2
    got = rows(tmp_path / "out.csv")
    assert [float(got[0][f"S_{c}"]) for c in ("xx", "yy", "zz", "xy", "yz", "xz")] == [
        11,
        22,
        33,
        12,
        23,
        13,
    ]
    assert [int(row["step"]) for row in got] == [1, 2]


@ABAQUS_TODO
def test_separate_real_frames(tmp_path):
    if not REPORT.exists():
        pytest.fail("Abaqus solver not available on workstation; real output generator is a TODO")
    other = tmp_path / "second.rpt"
    other.write_text(
        REPORT.read_text().replace("Increment      1", "Increment      2").replace("1.000", "2.000")
    )
    header = abaqus_report_to_csv([REPORT, other], tmp_path / "out.csv", inp=INP, **OPTIONS)
    assert header.n_steps == 2 and len(rows(tmp_path / "out.csv")) == 16


@pytest.mark.parametrize(
    "text, options, match",
    [
        (report(), {"stress_unit": None}, "stress_unit"),
        (report(), {"length_unit": None}, "length_unit"),
        (report(), {"coordinate_system": None}, "coordinate_system"),
        (report(), {"time_unit": None}, "time_unit"),
        (report() + report(), {}, "ascend"),
        (report() + report(2, 0), {}, "time"),
        (report() + report(2, 2).replace("S.S23", "S.Mises"), {}, "column|stress"),
        (report().replace("S.S", "U.U"), {}, "column|stress"),
        (report().replace("COORD.COOR3", "BAD"), {}, "column|coordinate"),
        (report(extra="Local coordinate system: 3"), {}, "global|local"),
        (report(extra="Section Point 1"), {}, "section"),
        (report(), {"averaging": "nodal_average_global"}, "averaging"),
        (report().replace("1 11 22", "1 nan 22"), {}, "finite"),
    ],
)
def test_refusals_preserve_existing_output(tmp_path, text, options, match):
    out = tmp_path / "out.csv"
    out.write_text("existing")
    with pytest.raises(ValueError, match=match):
        convert(tmp_path, text, **options)
    assert out.read_text() == "existing"


def test_element_nodal_agreement_and_disagreement(tmp_path):
    text = (
        report(position="Element Nodal")
        .replace("Node Label S.S11", "Element Node S.S11")
        .replace("1 11 22", "1 1 11 22")
    )
    convert(tmp_path, text + "2 1 11 22 33 12 13 23 0 0 0\n")
    assert len(rows(tmp_path / "out.csv")) == 1
    with pytest.raises(ValueError, match="Unique Nodal"):
        convert(tmp_path, text + "2 1 99 22 33 12 13 23 0 0 0\n")


@ABAQUS_TODO
def test_real_strain_unit_and_convention(tmp_path):
    if not REPORT.exists():
        pytest.fail("Abaqus solver not available on workstation; real output generator is a TODO")
    abaqus_report_to_csv(
        REPORT,
        tmp_path / "out.csv",
        inp=INP,
        **(OPTIONS | {"strain_unit": "microstrain", "strain_convention": "engineering"}),
    )
    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        np.testing.assert_allclose(
            block.E_hist[0, 0, 0, 2], 29.7323e-24 * 1e-6 / 2, rtol=1e-14, atol=0
        )


def test_refuse_mixed_time_basis(tmp_path):
    with pytest.raises(ValueError, match="time basis"):
        convert(tmp_path, report() + report(2, 2).replace("Step Time", "Total Time"))


@pytest.mark.parametrize(
    "declaration,match",
    [
        ("Local coordinate system: 3", "local|global"),
        ("Section Point 1", "section"),
        ("Material: other", "material"),
    ],
)
def test_report_wide_declarations_checked_before_first_frame(tmp_path, declaration, match):
    with pytest.raises(ValueError, match=match):
        convert(tmp_path, report().replace("Frame:", declaration + "\nFrame:"))


@ABAQUS_TODO
def test_refuse_flat_inp_reused_for_multiple_instances(tmp_path):
    if not REPORT.exists():
        pytest.fail("Abaqus solver not available on workstation; real output generator is a TODO")
    text = REPORT.read_text()
    path = tmp_path / "multiple.rpt"
    # One frame, two instances with overlapping labels: a flat INP cannot identify their positions.
    path.write_text(text + text[text.index("Loc 1") :].replace("PART-1-1", "PART-2-1"))
    with pytest.raises(ValueError, match="instance"):
        abaqus_report_to_csv(path, tmp_path / "out.csv", inp=INP, **OPTIONS)


@ABAQUS_TODO
def test_missing_strain_declaration_and_changed_grid(tmp_path):
    if not REPORT.exists():
        pytest.fail("Abaqus solver not available on workstation; real output generator is a TODO")
    with pytest.raises(ValueError, match="strain_convention"):
        abaqus_report_to_csv(
            REPORT, tmp_path / "out.csv", inp=INP, **(OPTIONS | {"strain_convention": None})
        )
    with pytest.raises(ValueError, match="grid|coordinates"):
        convert(tmp_path, report() + report(2, 2).replace("1 11 22", "2 11 22"))
