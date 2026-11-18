"""CSV exchange format: header parsing and column-name resolution (no h5py needed)."""

from __future__ import annotations

import pytest

from multifat.fe.csvio import (
    COLUMN_ALIASES,
    CSV_FORMAT,
    REQUIRED_HEADER_KEYS,
    CsvHeader,
    read_header,
    resolve_columns,
)

GOOD = {
    "format": CSV_FORMAT,
    "solver": "Abaqus",
    "solver_version": "2024",
    "stress_unit": "MPa",
    "length_unit": "mm",
    "coordinate_system": "global_cartesian",
    "averaging": "none",
    "analysis_type": "linear_elastic",
    "step_role": "time_sample",
}
STRESS = ["S_xx", "S_yy", "S_zz", "S_xy", "S_yz", "S_xz"]
BASE = ["node_id", "step", "x", "y", "z"] + STRESS


def _write(tmp_path, header=None, columns="node_id,step,x,y,z", extra_lines=()):
    header = dict(GOOD if header is None else header)
    lines = [f"# {k}: {v}" for k, v in header.items()] + list(extra_lines)
    path = tmp_path / "f.csv"
    path.write_text("\n".join(lines + [columns, "1,1,0,0,0"]) + "\n", encoding="utf-8")
    return path


def test_valid_header_parses_all_fields_and_keeps_unknown_keys(tmp_path):
    header = dict(
        GOOD,
        strain_convention="engineering",
        strain_unit="microstrain",
        cycle="stabilised",
        time_unit="s",
        temperature_unit="C",
        n_nodes="12",
        n_steps="3",
        load_description="max then min",
        model="tube",
        source="script v1",
        notes="a: b",
        vendor_key="kept",
    )
    path = _write(tmp_path, header)
    h, n_lines = read_header(path)
    assert isinstance(h, CsvHeader)
    assert n_lines == len(header)
    assert h.solver == "Abaqus" and h.solver_version == "2024"
    assert h.strain_convention == "engineering" and h.strain_unit == "microstrain"
    assert (h.n_nodes, h.n_steps) == (12, 3)
    assert h.cycle == "stabilised" and h.time_unit == "s" and h.temperature_unit == "C"
    assert h.notes == "a: b"  # only the first colon splits key and value
    assert h.extra == {"vendor_key": "kept"}


def test_header_tolerates_bom_crlf_and_empty_comment_lines(tmp_path):
    text = "﻿# format: multifat-fe-csv 1.0\r\n#\r\n" + "\r\n".join(
        f"# {k}: {v}" for k, v in GOOD.items() if k != "format"
    )
    path = tmp_path / "f.csv"
    path.write_bytes((text + "\r\nnode_id,step\r\n").encode("utf-8"))
    h, n_lines = read_header(path)
    assert h.solver == "Abaqus" and n_lines == len(GOOD) + 1


@pytest.mark.parametrize("key", sorted(REQUIRED_HEADER_KEYS))
def test_missing_required_key_names_the_key(tmp_path, key):
    header = {k: v for k, v in GOOD.items() if k != key}
    with pytest.raises(ValueError, match=key):
        read_header(_write(tmp_path, header))


def test_required_keys_are_exactly_the_documented_set():
    assert set(REQUIRED_HEADER_KEYS) == set(GOOD)


@pytest.mark.parametrize(
    "key, value",
    [
        ("format", "other-csv 1.0"),
        ("format", "multifat-fe-csv 2.0"),
        ("stress_unit", "furlong"),
        ("length_unit", "furlong"),
        ("coordinate_system", "cylindrical"),
        ("averaging", "smooth"),
        ("analysis_type", "plastic"),
        ("step_role", "sample"),
        ("strain_convention", "voigt"),
        ("strain_unit", "furlong"),
        ("n_nodes", "twelve"),
        ("n_steps", "0"),
        ("cycle", "later"),
    ],
)
def test_value_sets_enforced(tmp_path, key, value):
    with pytest.raises(ValueError, match=key):
        read_header(_write(tmp_path, dict(GOOD, **{key: value})))


def test_format_minor_version_accepted(tmp_path):
    h, _ = read_header(_write(tmp_path, dict(GOOD, format="multifat-fe-csv 1.3")))
    assert h.format == "multifat-fe-csv 1.3"


def test_duplicate_key_and_malformed_line_refused(tmp_path):
    with pytest.raises(ValueError, match="solver"):
        read_header(_write(tmp_path, extra_lines=["# solver: Ansys"]))
    with pytest.raises(ValueError, match="line"):
        read_header(_write(tmp_path, extra_lines=["# no colon here"]))


def test_file_without_column_line_refused(tmp_path):
    path = tmp_path / "f.csv"
    path.write_text("".join(f"# {k}: {v}\n" for k, v in GOOD.items()), encoding="utf-8")
    with pytest.raises(ValueError, match="column"):
        read_header(path)


def test_elastic_plastic_requires_stabilised_cycle(tmp_path):
    ep = dict(GOOD, analysis_type="elastic_plastic")
    with pytest.raises(ValueError, match="cycle"):
        read_header(_write(tmp_path, ep))
    with pytest.raises(ValueError, match="cycle"):
        read_header(_write(tmp_path, dict(ep, cycle="3")))
    with pytest.raises(ValueError, match="cycle"):
        read_header(_write(tmp_path, dict(ep, cycle="not_applicable")))
    h, _ = read_header(_write(tmp_path, dict(ep, cycle="stabilised")))
    assert h.cycle == "stabilised"
    read_header(_write(tmp_path, dict(GOOD, cycle="3")))  # integers and not_applicable are valid
    read_header(_write(tmp_path, dict(GOOD, cycle="not_applicable")))


def test_unit_case_requires_linear_elastic(tmp_path):
    bad = dict(GOOD, analysis_type="elastic_plastic", cycle="stabilised", step_role="unit_case")
    with pytest.raises(ValueError, match="unit_case"):
        read_header(_write(tmp_path, bad))
    read_header(_write(tmp_path, dict(GOOD, step_role="unit_case")))


def test_strain_convention_and_unit_required_when_strain_columns_exist_is_checked_by_the_importer_contract(
    tmp_path,
):
    h, _ = read_header(_write(tmp_path))  # the header alone does not need them
    assert h.strain_convention is None and h.strain_unit is None
    cmap = resolve_columns(BASE + ["Eel_xx", "Eel_yy", "Eel_zz", "Eel_xy", "Eel_yz", "Eel_xz"])
    with pytest.raises(ValueError, match="strain_convention"):
        cmap.check_header(h)
    with pytest.raises(ValueError, match="strain_unit"):
        cmap.check_header(
            CsvHeader.from_mapping(dict(GOOD, strain_convention="tensor"))  # no strain_unit
        )
    cmap.check_header(
        CsvHeader.from_mapping(dict(GOOD, strain_convention="tensor", strain_unit="1"))
    )


def test_time_unit_required_with_time_column():
    h = CsvHeader.from_mapping(GOOD)
    cmap = resolve_columns(BASE + ["time"])
    assert cmap.has_time
    with pytest.raises(ValueError, match="time_unit"):
        cmap.check_header(h)
    cmap.check_header(CsvHeader.from_mapping(dict(GOOD, time_unit="s")))


def test_elastic_plastic_without_strain_columns_refused():
    h = CsvHeader.from_mapping(dict(GOOD, analysis_type="elastic_plastic", cycle="stabilised"))
    with pytest.raises(ValueError, match="strain"):
        resolve_columns(BASE).check_header(h)


# Distinct sentinel per component: any swapped mapping fails.
ABAQUS_S = {
    "S11": "S_xx",
    "S22": "S_yy",
    "S33": "S_zz",
    "S12": "S_xy",
    "S23": "S_yz",
    "S13": "S_xz",
}
COMPONENT = {"xx": "11", "yy": "22", "zz": "33", "xy": "12", "yz": "23", "xz": "13"}
ANSYS_COMPONENT = {"xx": "X", "yy": "Y", "zz": "Z", "xy": "XY", "yz": "YZ", "xz": "XZ"}


def test_abaqus_names_resolve_with_s13_to_xz_and_s23_to_yz():
    for alias, canonical in ABAQUS_S.items():
        assert COLUMN_ALIASES[alias.lower()] == canonical
    for comp, num in COMPONENT.items():
        assert COLUMN_ALIASES[f"s{num}"] == f"S_{comp}"
        assert COLUMN_ALIASES[f"ee{num}"] == f"Eel_{comp}"
        assert COLUMN_ALIASES[f"pe{num}"] == f"Epl_{comp}"
        assert COLUMN_ALIASES[f"e{num}"] == f"E_{comp}"
    cmap = resolve_columns(
        ["NODE", "Step", "X", "Y", "Z", "S11", "S22", "S33", "S12", "S23", "S13"]
        + [f"EE{n}" for n in ("11", "22", "33", "12", "23", "13")]
        + [f"PE{n}" for n in ("11", "22", "33", "12", "23", "13")]
    )
    assert cmap.names["S_xz"] == "S13" and cmap.names["S_yz"] == "S23"
    assert cmap.names["Eel_xz"] == "EE13" and cmap.names["Epl_yz"] == "PE23"
    assert cmap.names["node_id"] == "NODE"
    assert cmap.has_split and not cmap.has_total


def test_ansys_names_resolve():
    for comp, suffix in ANSYS_COMPONENT.items():
        assert COLUMN_ALIASES[f"s{suffix.lower()}"] == f"S_{comp}"
        assert COLUMN_ALIASES[f"epel{suffix.lower()}"] == f"Eel_{comp}"
        assert COLUMN_ALIASES[f"eppl{suffix.lower()}"] == f"Epl_{comp}"
        assert COLUMN_ALIASES[f"epto{suffix.lower()}"] == f"E_{comp}"
    assert COLUMN_ALIASES["temp"] == "T"
    cmap = resolve_columns(
        ["NODE", "STEP", "X", "Y", "Z", "SX", "SY", "SZ", "SXY", "SYZ", "SXZ", "TEMP"]
    )
    assert cmap.names["S_xy"] == "SXY" and cmap.names["T"] == "TEMP"
    assert cmap.has_temperature and not cmap.has_strain


def test_short_stress_names_and_case_insensitive_names():
    cmap = resolve_columns(
        ["Node_ID", "STEP", "x", "Y", "z", "sxx", "SYY", "Szz", "sxy", "syz", "sxz", "Instance",
         "MATERIAL", "eel_xx", "EEL_YY", "Eel_zz", "eel_xy", "eel_yz", "eel_xz"]
    )  # fmt: skip
    assert cmap.has_instance and cmap.has_material and cmap.has_strain and not cmap.has_split
    assert cmap.names["instance"] == "Instance" and cmap.names["Eel_yy"] == "EEL_YY"
    assert cmap.has_elastic_only


def test_canonical_and_alias_together_rejected():
    with pytest.raises(ValueError, match="S_xx"):
        resolve_columns(BASE + ["S11"])
    with pytest.raises(ValueError, match="S_xx"):
        resolve_columns(BASE + ["SX"])
    with pytest.raises(ValueError, match="node_id"):
        resolve_columns(BASE + ["NODE"])


def test_unknown_column_rejected_unless_ignored():
    with pytest.raises(ValueError, match="Whatever"):
        resolve_columns(BASE + ["Whatever"])
    cmap = resolve_columns(BASE + ["Whatever"], ignore=("whatever",))
    assert cmap.ignored == ("Whatever",)
    assert "Whatever" not in cmap.names.values()


def test_missing_required_columns_listed():
    with pytest.raises(ValueError, match=r"S_yz.*node_id|node_id.*S_yz"):
        resolve_columns(["step", "x", "y", "z", "S_xx", "S_yy", "S_zz", "S_xy", "S_xz"] + ["time"])
    with pytest.raises(ValueError, match="step"):
        resolve_columns([c for c in BASE if c != "step"])


def test_partial_strain_group_rejected():
    with pytest.raises(ValueError, match="Eel"):
        resolve_columns(BASE + ["Eel_xx", "Eel_yy"])
    with pytest.raises(ValueError, match="E_"):
        resolve_columns(BASE + ["E11"])


def test_total_strain_with_a_split_part_rejected():
    total = [f"E_{c}" for c in ("xx", "yy", "zz", "xy", "yz", "xz")]
    eel = [f"Eel_{c}" for c in ("xx", "yy", "zz", "xy", "yz", "xz")]
    with pytest.raises(ValueError, match="total"):
        resolve_columns(BASE + total + eel)
    cmap = resolve_columns(BASE + total)
    assert cmap.has_total and not cmap.has_split and cmap.has_strain
    pl = [f"Epl_{c}" for c in ("xx", "yy", "zz", "xy", "yz", "xz")]
    assert resolve_columns(BASE + pl).has_plastic_only
    assert resolve_columns(BASE + eel + pl).has_split
