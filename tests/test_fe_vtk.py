"""VTK XML importer: conventions recovered from synthetic exports, rejections, provenance."""

from __future__ import annotations

import itertools
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from multifat.fe.model import sha256_of
from multifat.fe.synthetic import (
    TubeGeometry,
    TubeLoad,
    TubeMaterial,
    tube_cartesian,
    tube_points,
    write_tube_pvd,
)
from multifat.fe.tensors import CANONICAL, VTK9
from multifat.fe.vtk import (
    CylindricalFrame,
    TensorSpec,
    import_pvd,
    read_pvd,
    read_step,
    read_vtu,
)

GEOM = TubeGeometry(a=6.0, b=7.5)
MAT = TubeMaterial(E=210000.0, nu=0.3)
LOADS = [TubeLoad(0.0, 0.0, 0.0), TubeLoad(50.0, 12000.0, 0.0), TubeLoad(30.0, -5000.0, 8e4)]
ORDER_ODD = ("xx", "yy", "zz", "xz", "xy", "yz")  # a non-default order with distinct shear slots


def _series(tmp_path, **kw):
    kw.setdefault("n_radii", 3)
    kw.setdefault("n_theta", 5)
    return write_tube_pvd(tmp_path, GEOM, MAT, LOADS, **kw)


@pytest.mark.parametrize("source_kind", ["pvd", "vtu"])
def test_import_never_overwrites_a_source(tmp_path, source_kind):
    pvd = _series(tmp_path, n_radii=1, n_theta=1)
    source = pvd if source_kind == "pvd" else read_pvd(pvd)[1][1]
    original = source.read_bytes()
    with pytest.raises(ValueError, match="source|differ"):
        import_pvd(
            pvd,
            source,
            stress=TensorSpec("Stress", CANONICAL, "MPa"),
            strain=None,
            id_array="NodeID",
            material_id="steel",
            solver="test",
            solver_version="1",
            analysis_step="test",
            length_unit="mm",
        )
    assert source.read_bytes() == original


def _specs(order, shear, unit, location="point"):
    return (
        TensorSpec("Stress", order, unit=unit, location=location),
        TensorSpec("Strain", order, shear=shear, location=location),
    )


@pytest.mark.parametrize(
    "encoding,header,frame,order,shear,unit,length",
    [
        (e, h, f, o, s, u, ln)
        for e, h, f, o, s, (u, ln) in itertools.product(
            ["ascii", "base64"],
            ["UInt32"],
            ["cartesian", "cylindrical"],
            [CANONICAL, VTK9, ORDER_ODD],
            ["tensor", "engineering"],
            [("MPa", "mm"), ("Pa", "m")],
        )
        for h in ["UInt64" if e == "base64" and f == "cylindrical" and o == VTK9 else "UInt32"]
    ],
)
def test_reader_recovers_oracle(tmp_path, encoding, header, frame, order, shear, unit, length):
    pvd = _series(
        tmp_path, encoding=encoding, header_type=header, frame=frame, order=order, shear=shear,
        stress_unit=unit, length_unit=length,
    )  # fmt: skip
    stress, strain = _specs(order, shear, unit)
    frm = CylindricalFrame() if frame == "cylindrical" else "cartesian"
    ids, _, r, theta = tube_points(GEOM, 3, 5)
    steps = read_pvd(pvd)
    assert len(steps) == len(LOADS)
    for (_, path), load in zip(steps, LOADS):
        got = read_step(
            read_vtu(path), stress, strain, frame=frm, length_unit=length, id_array="NodeID"
        )
        S, E = tube_cartesian(r, theta, load, GEOM, MAT)
        assert np.array_equal(got.location_id, ids)
        assert np.allclose(got.S, S, rtol=1e-10, atol=1e-9)
        assert np.allclose(got.E, E, rtol=1e-10, atol=1e-15)


def test_cell_data_uses_centroids(tmp_path):
    pvd = _series(tmp_path, data_location="cell")
    stress, strain = _specs(CANONICAL, "tensor", "MPa", location="cell")
    _, xyz, r, theta = tube_points(GEOM, 3, 5)
    got = read_step(
        read_vtu(read_pvd(pvd)[1][1]), stress, strain, length_unit="mm", id_array="NodeID"
    )
    S, _ = tube_cartesian(r, theta, LOADS[1], GEOM, MAT)
    assert np.allclose(got.coordinates, xyz) and np.allclose(got.S, S, rtol=1e-10, atol=1e-9)


def _first_vtu(pvd):
    return read_pvd(pvd)[0][1]


def test_rejects_compressed_vtu(tmp_path):
    pvd = _series(tmp_path, encoding="base64")
    path = _first_vtu(pvd)
    path.write_text(
        path.read_text().replace("<VTKFile ", '<VTKFile compressor="vtkZLibDataCompressor" ', 1)
    )
    with pytest.raises(ValueError, match="vtkZLibDataCompressor"):
        read_vtu(path)


def test_rejects_appended_data(tmp_path):
    pvd = _series(tmp_path)
    path = _first_vtu(pvd)
    path.write_text(path.read_text().replace('format="ascii"', 'format="appended"', 1))
    with pytest.raises(ValueError, match="appended"):
        read_vtu(path)


def test_missing_array_name_lists_available(tmp_path):
    vtu = read_vtu(_first_vtu(_series(tmp_path)))
    with pytest.raises(KeyError, match="Stress"):
        read_step(vtu, TensorSpec("Nope", CANONICAL, unit="MPa"), None, length_unit="mm")


def test_wrong_column_count_rejected(tmp_path):
    vtu = read_vtu(_first_vtu(_series(tmp_path)))
    with pytest.raises(ValueError, match="columns"):
        read_step(vtu, TensorSpec("Stress", VTK9, unit="MPa"), None, length_unit="mm")


def _import(pvd, out, **kw):
    stress, strain = _specs(CANONICAL, "tensor", "MPa")
    args = dict(
        stress=stress, strain=strain, id_array="NodeID", material_id="tube", solver="synthetic",
        solver_version="1", analysis_step="s", length_unit="mm",
    )  # fmt: skip
    args.update(kw)
    return import_pvd(pvd, out, **args)


def test_import_pvd_to_universal_file_matches_oracle(tmp_path):
    pytest.importorskip("h5py")
    from multifat.fe.store import HistoryFile

    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in", encoding="base64")
    meta = _import(pvd, tmp_path / "out.h5")
    ids, xyz, r, theta = tube_points(GEOM, 3, 5)
    with HistoryFile(tmp_path / "out.h5") as f:
        block = f.read(slice(None))
        assert f.meta == meta and f.load_cases.labels[0] == "step_0000.vtu"
        assert list(f.load_cases.time) == [0.0, 1.0, 2.0]
    assert np.array_equal(block.location_id, ids) and np.allclose(block.coordinates, xyz)
    for k, load in enumerate(LOADS):
        S, E = tube_cartesian(r, theta, load, GEOM, MAT)
        assert np.allclose(block.S_hist[:, k], S, rtol=1e-10, atol=1e-9)
        assert np.allclose(block.E_hist[:, k], E, rtol=1e-10, atol=1e-15)


def test_cell_import_derives_element_centroid_kind(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in", data_location="cell")
    stress, strain = _specs(CANONICAL, "tensor", "MPa", location="cell")
    meta = _import(pvd, tmp_path / "out.h5", stress=stress, strain=strain)
    assert meta.location_kind == "element_centroid"


def test_provenance_records_sha256_of_each_file(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in")
    meta = _import(pvd, tmp_path / "out.h5", frame="cartesian")
    files = {s.name: s.sha256 for s in meta.source_files}
    assert files[pvd.name] == sha256_of(pvd)
    for _, path in read_pvd(pvd):
        assert files[path.name] == sha256_of(path)
    assert (
        meta.source_units == {"stress": "MPa", "length": "mm"} and meta.source_frame == "cartesian"
    )
    assert meta.created and meta.adapter == "multifat.fe.vtk"


def _edit(path, fn):
    tree = ET.parse(path)
    fn(tree.getroot())
    tree.write(path)


def test_mismatched_ids_between_steps_raises(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in")

    def swap(root):
        arr = next(a for a in root.iter("DataArray") if a.get("Name") == "NodeID")
        vals = arr.text.split()
        vals[0], vals[1] = vals[1], vals[0]
        arr.text = " ".join(vals)

    _edit(read_pvd(pvd)[1][1], swap)
    with pytest.raises(ValueError, match="IDs"):
        _import(pvd, tmp_path / "out.h5")
    assert not (tmp_path / "out.h5").exists()


def test_mismatched_point_count_between_steps_raises(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    pvd = _series(tmp_path / "a")
    other = _series(tmp_path / "b", n_theta=4)
    shutil.copy(read_pvd(other)[1][1], read_pvd(pvd)[1][1])
    with pytest.raises(ValueError, match="point count"):
        _import(pvd, tmp_path / "out.h5")


def test_duplicate_ids_rejected(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in")

    def dup(root):
        arr = next(a for a in root.iter("DataArray") if a.get("Name") == "NodeID")
        vals = arr.text.split()
        vals[1] = vals[0]
        arr.text = " ".join(vals)

    for _, path in read_pvd(pvd):
        _edit(path, dup)
    with pytest.raises(ValueError, match="duplicate"):
        _import(pvd, tmp_path / "out.h5")


@pytest.mark.parametrize("layout", ["joint", "separate"])
@pytest.mark.parametrize("header", ["UInt32", "UInt64"])
def test_both_base64_layouts_are_read(tmp_path, layout, header):
    pvd = _series(tmp_path, encoding="base64", header_type=header, base64_layout=layout)
    stress, strain = _specs(CANONICAL, "tensor", "MPa")
    _, _, r, theta = tube_points(GEOM, 3, 5)
    got = read_step(
        read_vtu(read_pvd(pvd)[1][1]), stress, strain, length_unit="mm", id_array="NodeID"
    )
    S, _ = tube_cartesian(r, theta, LOADS[1], GEOM, MAT)
    assert np.allclose(got.S, S, rtol=1e-10, atol=1e-9)


def test_jointly_encoded_stream_built_by_hand_is_read(tmp_path):
    """The layout VTK writes for uncompressed inline data: header and data in one base64 stream."""
    import base64

    values = np.arange(12, dtype="<f8") * 1.5
    raw = values.tobytes()
    body = base64.b64encode(len(raw).to_bytes(4, "little") + raw).decode()
    xml = (
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" '
        'header_type="UInt32"><UnstructuredGrid><Piece NumberOfPoints="4" NumberOfCells="0">'
        '<PointData/><Points><DataArray type="Float64" NumberOfComponents="3" format="binary">'
        f"{body}</DataArray></Points></Piece></UnstructuredGrid></VTKFile>"
    )
    path = tmp_path / "joint.vtu"
    path.write_text(xml)
    assert np.array_equal(read_vtu(path).points, values.reshape(4, 3))


def test_garbled_binary_reports_layout_problem(tmp_path):
    xml = (
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" '
        'header_type="UInt32"><UnstructuredGrid><Piece NumberOfPoints="1" NumberOfCells="0">'
        '<PointData/><Points><DataArray type="Float64" NumberOfComponents="3" format="binary">'
        "AAAAAAAAAAAAAAAA</DataArray></Points></Piece></UnstructuredGrid></VTKFile>"
    )
    path = tmp_path / "bad.vtu"
    path.write_text(xml)
    with pytest.raises(ValueError, match="header"):
        read_vtu(path)


def test_stress_unit_is_required(tmp_path):
    vtu = read_vtu(_first_vtu(_series(tmp_path)))
    with pytest.raises(ValueError, match="unit"):
        read_step(vtu, TensorSpec("Stress", CANONICAL), None, length_unit="mm")


def test_length_unit_is_required(tmp_path):
    pvd = _series(tmp_path)
    vtu = read_vtu(_first_vtu(pvd))
    with pytest.raises(TypeError, match="length_unit"):
        read_step(vtu, TensorSpec("Stress", CANONICAL, unit="MPa"), None)
    pytest.importorskip("h5py")
    stress, strain = _specs(CANONICAL, "tensor", "MPa")
    with pytest.raises(TypeError, match="length_unit"):
        import_pvd(
            pvd, tmp_path / "o.h5", stress=stress, strain=strain, id_array="NodeID",
            material_id="m", solver="s", solver_version="1", analysis_step="a",
        )  # fmt: skip


def _add_array(path, name, values, n_comp):
    def add(root):
        node = root.find("UnstructuredGrid/Piece/PointData")
        arr = ET.SubElement(
            node,
            "DataArray",
            type="Float64",
            Name=name,
            NumberOfComponents=str(n_comp),
            format="ascii",
        )
        arr.text = " ".join(f"{v:.17g}" for v in np.asarray(values).ravel())

    _edit(path, add)


_RZ90 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


def test_rotation_array_rotates_tensors_and_identity_is_a_noop(tmp_path):
    from multifat.fe.vtk import RotationArray

    path = _first_vtu(_series(tmp_path))
    n = 15
    _add_array(path, "Rot", np.tile(_RZ90, (n, 1, 1)), 9)
    stress, strain = _specs(CANONICAL, "tensor", "MPa")
    plain = read_step(read_vtu(path), stress, strain, length_unit="mm")
    got = read_step(read_vtu(path), stress, strain, frame=RotationArray("Rot"), length_unit="mm")
    assert np.allclose(got.S, _RZ90 @ plain.S @ _RZ90.T) and np.allclose(
        got.E, _RZ90 @ plain.E @ _RZ90.T
    )


@pytest.mark.parametrize(
    "matrix,message",
    [(2.0 * _RZ90, "orthonormal"), (np.diag([1.0, 1.0, -1.0]), "determinant")],
)
def test_rotation_array_must_be_a_proper_rotation(tmp_path, matrix, message):
    from multifat.fe.vtk import RotationArray

    path = _first_vtu(_series(tmp_path))
    _add_array(path, "Rot", np.tile(matrix, (15, 1, 1)), 9)
    stress, strain = _specs(CANONICAL, "tensor", "MPa")
    with pytest.raises(ValueError, match=message):
        read_step(read_vtu(path), stress, strain, frame=RotationArray("Rot"), length_unit="mm")


def test_multi_part_and_repeated_timesteps_rejected(tmp_path):
    pvd = _series(tmp_path)
    _edit(pvd, lambda root: next(root.iter("DataSet")).set("part", "1"))
    with pytest.raises(ValueError, match="part"):
        read_pvd(pvd)
    pvd2 = _series(tmp_path / "again")
    _edit(pvd2, lambda root: list(root.iter("DataSet"))[1].set("timestep", "0"))
    with pytest.raises(ValueError, match="timestep"):
        read_pvd(pvd2)


def test_conventions_recorded_and_strain_kind_none_without_strain(tmp_path):
    pytest.importorskip("h5py")
    (tmp_path / "in").mkdir()
    pvd = _series(tmp_path / "in", order=ORDER_ODD, shear="engineering")
    stress, strain = _specs(ORDER_ODD, "engineering", "MPa")
    meta = _import(pvd, tmp_path / "a.h5", stress=stress, strain=strain)
    assert meta.source_conventions == {
        "stress_order": "xx,yy,zz,xz,xy,yz",
        "strain_order": "xx,yy,zz,xz,xy,yz",
        "strain_shear": "engineering",
        "location": "point",
    }
    assert meta.strain_kind == "total_mechanical"
    meta2 = _import(pvd, tmp_path / "b.h5", stress=stress, strain=None)
    assert meta2.strain_kind is None and "strain_order" not in meta2.source_conventions
