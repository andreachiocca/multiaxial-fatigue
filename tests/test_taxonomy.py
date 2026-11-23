"""The runtime taxonomy (issue #61) agrees with the catalogue and covers every registered model."""

from __future__ import annotations

import json
from pathlib import Path

from multifat.models import registry
from multifat.models.taxonomy import APPROACHES, QUANTITIES, TAXONOMY, approach_label, taxonomy

CATALOGUE = Path(__file__).resolve().parents[1] / "docs/catalogue/methods.json"


def _records():
    methods = json.loads(CATALOGUE.read_text(encoding="utf-8"))["methods"]
    return {m["registry_name"]: m for m in methods if m["registry_name"] and m["enabled"]}


def test_every_registered_model_is_classified():
    assert set(TAXONOMY) == set(registry.base_models())


def test_runtime_table_matches_the_catalogue():
    for name, m in _records().items():
        t = TAXONOMY[name]
        assert (t.approach, t.quantity, t.np_sensitivity, list(t.tags)) == (
            m["approach"],
            m["quantity"],
            m["np_sensitivity"],
            m["taxonomy_tags"],
        ), name


def test_vocabulary():
    for t in TAXONOMY.values():
        assert t.approach in APPROACHES and t.quantity in QUANTITIES
        assert t.np_sensitivity in ("emergent", "explicit")


def test_reduction_rule_cases():
    """The cases the reduction rule decides against the former directory or catalogue family."""
    assert taxonomy("DANGVAN").approach == "critical_plane"
    assert taxonomy("MOREL").approach == "critical_plane"  # max over planes of T_s
    assert taxonomy("PAPADOPOULOS").approach == "critical_plane"  # 2001 form: max_n T_a
    assert taxonomy("PCR").approach == "critical_plane"
    assert taxonomy("PI").approach == "integral"
    assert taxonomy("TVEIT2025").approach == "incremental"
    assert taxonomy("BOLCHOUN").np_sensitivity == "explicit"
    assert approach_label("FIN_ext") == "Critical plane"
    assert approach_label("NOT_A_MODEL") == "Unknown"
