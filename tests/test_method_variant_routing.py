"""Default and explicit method variants use a common public selector."""

import pytest

from multifat.cli.method_selection import resolve_method_tag, select_method_variants
from multifat.cli.run import _split_model_variant


@pytest.mark.parametrize(
    ("tag", "expected"),
    [
        ("FS", ("FS", "")),
        ("FS_ext", ("FS", "_ext")),
        ("FIN", ("FIN", "")),
        ("FIN_ext", ("FIN", "_ext")),
        ("CROSSLAND", ("CROSSLAND", "")),
    ],
)
def test_resolver_matches_single_model_runner(tag, expected):
    result = resolve_method_tag(tag)
    assert (result.base_name, result.suffix) == expected
    assert result.tag == tag
    assert _split_model_variant(tag) == expected


def test_request_expansion_preserves_order_and_deduplicates():
    selected = select_method_variants(["FS", "FIN", "CROSSLAND", "FS_ext"])
    assert [request.tag for request in selected] == ["FS", "FIN", "CROSSLAND", "FS_ext"]
    expanded = select_method_variants(["FS", "FIN", "CROSSLAND", "FIN"], all_variants=True)
    assert [request.tag for request in expanded] == ["FS", "FS_ext", "FIN", "FIN_ext", "CROSSLAND"]


@pytest.mark.parametrize("tag", ["CROSSLAND_ext", "SINES_ext", "MISSING", "", "FS_ext_ext"])
def test_unknown_or_unsupported_variants_fail_early(tag):
    with pytest.raises((ValueError, KeyError)):
        resolve_method_tag(tag)


def test_only_explicit_ext_does_not_add_plain_variant():
    assert [entry.tag for entry in select_method_variants(["FS_ext"])] == ["FS_ext"]
