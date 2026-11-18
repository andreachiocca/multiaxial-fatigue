"""Material constants outside a criterion's validity range skip the model, not fail it."""

from pathlib import Path

import pytest

from multifat.cli.pipeline import _eligibility_status, run_workbook
from multifat.models.errors import LoadingScopeError, MaterialScopeError
from multifat.models.integral_methods.robert_fogue_bahuaud import rfb_parameters

DATA = Path(__file__).resolve().parents[1] / "Experimental_Data"


@pytest.mark.parametrize("ratio", [0.517, 0.57, 0.87, 1.0])
def test_rfb_outside_its_validity_raises_material_scope(ratio):
    """Robert et al., Appendix: 0.577 < t_-1/s_-1 < 0.866 (4/3 < k < 3)."""
    with pytest.raises(MaterialScopeError, match="0.577 < t_-1/s_-1 < 0.866"):
        rfb_parameters(400.0, 400.0 * ratio, 600.0)


def test_material_scope_is_not_a_loading_scope():
    # load-scaling searches treat LoadingScopeError as an out-of-scope trial load
    assert not issubclass(MaterialScopeError, LoadingScopeError)
    assert issubclass(MaterialScopeError, ValueError)


def test_rfb_inside_its_validity_runs():
    a, b, d = rfb_parameters(400.0, 400.0 * 0.7, 600.0)
    assert a > 0 and b > 0


def test_pipeline_skips_rfb_on_aisi316l(tmp_path):
    """AISI316L: t_-1/s_-1 = 0.517, below the RFB range; RFB was reported as failed."""
    (outcome,) = run_workbook(
        DATA / "8_AISI316L.xlsx", models=["RFB"], results_dir=tmp_path, calibrate=False
    )
    assert outcome.status == "skipped"
    assert outcome.mode == "out_of_scope"
    assert "received 0.517" in outcome.reason
    assert _eligibility_status(outcome, calibrate=False) == "out_of_scope"
    assert not (tmp_path / "RFB_AISI316L.csv").exists()
