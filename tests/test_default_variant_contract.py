"""Regression contract for the unqualified FS and FIN plane-selection rules.

The _ext suffix is an explicit alternative requested by the caller.
These checks intentionally do not claim that the pipeline already supports
a single-variant selection or immutable JSON calibration history.
"""

import numpy as np
import pytest

from multifat.cli.run import _split_model_variant
from multifat.models.cp_methods.fin import Findley
from multifat.models.cp_methods.fs import FatemiSocie


@pytest.mark.parametrize(
    ("method_tag", "expected"),
    [
        ("FS", ("FS", "")),
        ("FS_ext", ("FS", "_ext")),
        ("FIN", ("FIN", "")),
        ("FIN_ext", ("FIN", "_ext")),
    ],
)
def test_plain_method_and_explicit_ext_tags(method_tag, expected):
    assert _split_model_variant(method_tag) == expected


def test_fs_default_plane_metric_is_shear_strain_amplitude():
    S0r = np.zeros((3, 3))
    S1r = np.zeros((3, 3))
    E0r = np.zeros((3, 3))
    E1r = np.zeros((3, 3))
    S0r[2, 2] = 100.0
    S1r[2, 2] = -100.0
    E0r[0, 2] = 0.01
    E1r[0, 2] = -0.01

    result = FatemiSocie().evaluate_on_plane(
        S0r=S0r,
        S1r=S1r,
        E0r=E0r,
        E1r=E1r,
        params={"Sy": 500.0, "k_FS": 0.5},
    )
    # FS default selects by shear strain, not by the full damage parameter.
    assert result.metric == pytest.approx(0.02)
    assert result.damage == pytest.approx(0.022)


def test_findley_default_metric_is_the_complete_damage_parameter():
    S0r = np.zeros((3, 3))
    S1r = np.zeros((3, 3))
    E0r = np.zeros((3, 3))
    E1r = np.zeros((3, 3))
    S0r[0, 2] = 50.0
    S1r[0, 2] = -50.0
    S0r[2, 2] = 100.0
    S1r[2, 2] = -100.0

    result = Findley().evaluate_on_plane(
        S0r=S0r,
        S1r=S1r,
        E0r=E0r,
        E1r=E1r,
        params={"k_FI": 0.2},
    )
    # tau_a = 50 MPa, k_FI * sigma_n,max = 20 MPa.
    assert result.damage == pytest.approx(70.0)
    assert result.metric == pytest.approx(result.damage)
