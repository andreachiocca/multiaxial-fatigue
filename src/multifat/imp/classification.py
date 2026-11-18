"""Classification of an imported material block into an assessment workflow (issue #18).

One importer, one workbook per material; the classification decides which assessment the data
support and gives every plain observation an explicit role::

    fatigue_life     >= 2 uniaxial tension R = -1 failures with distinct lives (design set):
                     master curve on the design set, other failures assessed, run-outs censored
    fatigue_limit    no plain failure, run-outs present and both R = -1 fatigue limits usable:
                     the run-outs are taken as limit states and judged by the fatigue index
    ambiguous        data that support neither workflow (failures without a design set; run-outs
                     without usable limits); nothing is assessed, the reasons are reported
    out_of_scope     no plain (unnotched) observation

Rules that keep the scientific distinction explicit:

- a run-out is a censored life observation (a lower bound on life), never a failure;
- the Dataset column ``Observation type`` (decision d113) declares each run-out row either an
  interrupted test (``run-out``, censored) or an established limit state (``fatigue limit``,
  from a staircase or derived from an S-N curve); a declared type is never overridden;
- an undeclared run-out (blank cell, or a workbook without the column) follows the rule before
  d113: a limit state in a workbook without failures, a censored row otherwise, because
  finite-cycle survival and an established limit need different evidence;
- a ``fatigue_life`` workbook with declared fatigue limits is mixed: its failures go through the
  life assessment and its limit states through the fatigue index (``roles["limit_state"]``);
- a limit state stopped below the reference life of the Summary limits is reported in
  ``diagnostics``.

Repository choices:
  - the classification reproduces the routing of the pipeline before issue #18 for every
    unambiguous workbook (failures -> life path, run-outs only with limits -> limit path), and
    the d113 types of the supplied workbooks reproduce it too;
  - the reference life of a limit is the smallest ``reference_cycles`` among its R = -1 sources.
"""

from __future__ import annotations

from typing import Any

WORKFLOWS = ("fatigue_life", "fatigue_limit", "ambiguous", "out_of_scope")


def classify_dataset(block: dict, material: str) -> dict[str, Any]:
    """Workflow, per-observation roles, counts and diagnostics of one imported material block."""
    from multifat.imp import experimental_data as exp
    from multifat.report.subsets import observation_index

    cases = block.get("cases", {})
    runouts = block.get("runouts", {})
    props = block.get("material", {})
    flagged = block.get("property_issues", {})
    limits = {
        key: props.get(field)
        for key, field in (
            ("Sigm1", "Tensile_fatigue_limit_MPa"),
            ("Taum1", "Torsional_fatigue_limit_MPa"),
        )
    }
    usable_limits = all(v is not None and v > 0 for v in limits.values()) and not (
        {"Sigm1", "Taum1"} & set(flagged)
    )
    points = exp.extract_uniaxial_tension({material: block}).get(material, {}).get("points", [])
    design_ids = sorted(int(pt["case_id"]) for pt in points)
    design_lives = {float(pt["Nf_cycles"]) for pt in points}
    reasons: list[str] = []
    diagnostics: list[str] = []

    declared = {int(k): (r.get("meta") or {}).get("observation_type") for k, r in runouts.items()}
    limit_ids = sorted(
        k for k, kind in declared.items() if kind == "fatigue limit" or (kind is None and not cases)
    )
    censored_ids = sorted(set(declared) - set(limit_ids))

    if not cases and not runouts:
        workflow = "out_of_scope"
        reasons.append("no plain (unnotched) observation")
    elif cases and len(design_lives) >= 2:
        workflow = "fatigue_life"
        reasons.append(
            f"{len(design_ids)} uniaxial tension R = -1 failures with {len(design_lives)} "
            "distinct lives form the design set"
        )
        if censored_ids:
            diagnostics.append(
                f"{len(censored_ids)} run-outs are censored life observations, not fatigue limits"
            )
        if limit_ids and usable_limits:
            reasons.append(
                f"mixed: {len(limit_ids)} declared fatigue limits are also judged by the fatigue "
                f"index at Sigm1 = {limits['Sigm1']:g} and Taum1 = {limits['Taum1']:g} MPa"
            )
        elif limit_ids:
            diagnostics.append(
                f"{len(limit_ids)} declared fatigue limits are not assessed: the Summary lacks "
                "usable R = -1 fatigue limits"
            )
            limit_ids = []
    elif cases:
        workflow = "ambiguous"
        reasons.append(
            "failures present, but fewer than 2 uniaxial tension R = -1 failures with distinct "
            f"lives (found {len(design_ids)}): no master curve can be fitted"
        )
    elif not limit_ids:
        workflow = "ambiguous"
        reasons.append(
            f"{len(runouts)} run-outs declared as interrupted tests and no failure: neither a "
            "master curve nor a limit state"
        )
    elif usable_limits:
        workflow = "fatigue_limit"
        reasons.append(
            f"{len(limit_ids)} run-outs and no failure: limit states judged at the R = -1 limits "
            f"Sigm1 = {limits['Sigm1']:g} and Taum1 = {limits['Taum1']:g} MPa"
        )
        if censored_ids:
            diagnostics.append(
                f"{len(censored_ids)} run-outs declared as interrupted tests are not assessed "
                f"(ids {censored_ids})"
            )
    else:
        workflow = "ambiguous"
        missing = [k for k, v in limits.items() if v is None or v <= 0 or k in flagged]
        reasons.append(
            "run-outs only, and the Summary lacks usable R = -1 fatigue limits "
            f"({', '.join(missing)}): neither a master curve nor a fatigue index can be built"
        )
    if workflow not in ("fatigue_life", "fatigue_limit"):
        limit_ids = []
    reference = _reference_life(block) if limit_ids else None
    if reference is not None:
        short = [k for k in limit_ids if float(runouts[k]["meta"]["Nf_cycles"]) < reference]
        if short:
            diagnostics.append(
                f"{len(short)} run-outs stopped below the {reference:g}-cycle reference of the "
                f"Summary limits (ids {short}): finite-cycle survivals, not established limits"
            )

    ids = observation_index(block)

    def entries(kind: str, keys) -> list[dict]:
        return [{"id": int(k), **ids[(kind, int(k))]} for k in keys]

    design_set = set(design_ids) if workflow == "fatigue_life" else set()
    roles = {
        "design": entries("case", sorted(design_set)),
        "assessment": (
            entries("case", sorted(int(k) for k in cases if int(k) not in design_set))
            if workflow == "fatigue_life"
            else []
        ),
        "censored": entries("runout", censored_ids) if workflow == "fatigue_life" else [],
        "limit_state": entries("runout", limit_ids),
    }
    return {
        "material": material,
        "workflow": workflow,
        "reasons": reasons,
        "diagnostics": diagnostics,
        "counts": {
            "failures": len(cases),
            "runouts": len(runouts),
            "design": len(roles["design"]),
            "limit_states": len(roles["limit_state"]),
            "excluded": len(block.get("excluded_observations", [])),
        },
        "fatigue_limits_usable": usable_limits,
        "source_workbook": block.get("source", {}).get("workbook"),
        "roles": roles,
    }


def limit_state_runouts(block: dict, material: str) -> dict:
    """The run-outs of `block` that the classification takes as limit states, by run-out id."""
    ids = {entry["id"] for entry in classify_dataset(block, material)["roles"]["limit_state"]}
    return {k: v for k, v in block.get("runouts", {}).items() if int(k) in ids}


def _reference_life(block: dict) -> float | None:
    """Smallest reference life of the R = -1 limit sources, or None when none is recorded."""
    lives = [
        float(source["reference_cycles"])
        for sources in (block.get("fatigue_limit_sources") or {}).values()
        for source in sources
        if source.get("reference_cycles")
    ]
    return min(lives) if lives else None
