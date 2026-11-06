"""Fetch-free documented PRNSOL/NLIST and complete Workbench tensor tables."""

import pytest
import numpy as np

from multifat.fe.solvers.ansys import ansys_listing_to_csv
from multifat.fe.csvio import import_fe_csv
from multifat.fe.store import HistoryFile

OPTIONS = dict(
    stress_unit="MPa",
    length_unit="mm",
    analysis_type="linear_elastic",
    step_role="time_sample",
    coordinate_system="global_cartesian",
    averaging="nodal_average_within_material",
    time_unit="s",
)
NLIST = """LIST ALL SELECTED NODES. DSYS= 0
NODE X Y Z THXY THYZ THZX
1 10 20 30 0 0 0
2 40 50 60 0 0 0
"""


def listing(substep=1, time=1):
    return f"""***** POST1 NODAL STRESS LISTING *****
LOAD STEP= 1 SUBSTEP= {substep}
TIME= {time} LOAD CASE= 0
THE FOLLOWING X,Y,Z VALUES ARE IN GLOBAL COORDINATES
NODE SX SY SZ SXY SYZ SXZ
1 11 22 33 12 23 13
2 41 52 63 42 53 43
"""


def convert(tmp_path, text, coordinates=NLIST, **options):
    path = tmp_path / "result.out"
    path.write_text(text, encoding="utf-8")
    nodes = tmp_path / "nodes.lis"
    nodes.write_text(coordinates, encoding="utf-8")
    return ansys_listing_to_csv(path, tmp_path / "out.csv", nodes=nodes, **(OPTIONS | options))


def test_prnsol_nlist_roundtrip(tmp_path):
    header = convert(tmp_path, listing() + listing(2, 2))
    assert header.n_steps == 2
    assert "SUBSTEP" in header.notes
    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        np.testing.assert_array_equal(
            block.S_hist[0, 0], [[11, 12, 13], [12, 22, 23], [13, 23, 33]]
        )
        np.testing.assert_array_equal(block.coordinates, [[10, 20, 30], [40, 50, 60]])


def test_workbench_complete_tensor_table(tmp_path):
    text = "Node Number\tX [mm]\tY [mm]\tZ [mm]\tSX [MPa]\tSY [MPa]\tSZ [MPa]\tSXY [MPa]\tSYZ [MPa]\tSXZ [MPa]\n1\t0\t0\t0\t11\t22\t33\t12\t23\t13\n"
    header = convert(tmp_path, text, step=3, time=5)
    assert header.n_nodes == 1
    assert "3" in header.notes


@pytest.mark.parametrize(
    "text, coordinates, options, match",
    [
        (listing().replace("GLOBAL COORDINATES", "LOCAL COORDINATES"), NLIST, {}, "global|local"),
        (listing(), NLIST.replace("DSYS= 0", "DSYS= 1"), {}, "DSYS|global"),
        (listing(), NLIST.replace("1 10 20 30", "1 99 20 30") + "1 0 0 0 0 0 0\n", {}, "duplicate"),
        (listing() + listing(), NLIST, {}, "ascend"),
        (listing() + listing(2, 0), NLIST, {}, "time"),
        (listing().replace(" SXZ", " BAD"), NLIST, {}, "column|stress"),
        (listing().replace("1 11 22", "1 11"), NLIST, {}, "ragged"),
        (listing(), NLIST, {"stress_unit": None}, "stress_unit"),
        (listing(), NLIST, {"averaging": None}, "averaging"),
        (listing(), "", {}, "coordinate"),
        (
            listing().replace("NODAL STRESS", "ELEMENT NODAL STRESS")
            + "SHELL RESULTS FOR TOP/BOTTOM\n",
            NLIST,
            {},
            "shell|element",
        ),
    ],
)
def test_ansys_refusals(tmp_path, text, coordinates, options, match):
    with pytest.raises(ValueError, match=match):
        convert(tmp_path, text, coordinates, **options)
    assert not (tmp_path / "out.csv").exists()


def test_workbench_units_checked(tmp_path):
    text = "Node Number\tX [m]\tY [m]\tZ [m]\tSX [Pa]\tSY [Pa]\tSZ [Pa]\tSXY [Pa]\tSYZ [Pa]\tSXZ [Pa]\n1\t0\t0\t0\t1\t2\t3\t4\t5\t6\n"
    with pytest.raises(ValueError, match="unit"):
        convert(tmp_path, text, step=1, time=1)


def test_numeric_load_steps_must_ascend_even_without_time(tmp_path):
    text = listing().replace("LOAD STEP= 1", "LOAD STEP= 2") + listing()
    text = "\n".join(line for line in text.splitlines() if not line.startswith("TIME="))
    with pytest.raises(ValueError, match="ascend"):
        convert(tmp_path, text)


def test_inhouse_listing_fatigue_cycle(tmp_path):
    from pathlib import Path

    data_dir = Path(__file__).parent / "data/solvers/ansys"
    out_file = data_dir / "hole_plate_fatigue.out"
    if not out_file.exists():
        pytest.skip("in-house holed plate fatigue fixture not found")

    header = ansys_listing_to_csv(
        out_file,
        tmp_path / "out.csv",
        nodes=out_file,
        stress_unit="MPa",
        length_unit="mm",
        time_unit="s",
        coordinate_system="global_cartesian",
        averaging="nodal_average_global",
        analysis_type="linear_elastic",
        step_role="time_sample",
    )
    assert header.n_nodes == 2864
    assert header.n_steps == 2

    pytest.importorskip("h5py")
    import_fe_csv(tmp_path / "out.csv", tmp_path / "out.h5")
    with HistoryFile(tmp_path / "out.h5") as history:
        block = next(history.iter_chunks())
        assert block.n_locations == 2864
        assert block.S_hist.shape[1] == 2
        s1 = block.S_hist[:, 0]
        s2 = block.S_hist[:, 1]
        # Tensile concentration peaks match by 90-deg symmetry
        np.testing.assert_allclose(s1[:, 1, 1].max(), s2[:, 0, 0].max(), rtol=0.03)
