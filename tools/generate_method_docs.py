#!/usr/bin/env python3
"""Generate contributor documentation from docs/catalogue/methods.json (stdlib only)."""
import argparse
import ast
import collections
import csv
import io
import json
import re
from pathlib import Path
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/catalogue/methods.json"
STATUSES = ["fully_implemented", "partially_implemented", "not_implemented", "requires_review"]
LABEL = dict(
    zip(
        STATUSES,
        ["Fully implemented", "Partially implemented", "Not implemented", "Requires review"],
    )
)
ISSUE_IDS = [
    "I_LIU2021",
    "I_MGSE_ZHU",
    "I_MKBM",
    "I_LI",
    "I_MSWT",
    "I_CARSPA",
    "I_DANGVAN",
    "I_OTT",
    "I_SWTD",
    "I_ZHU_EDP",
    "N01",
    "N02",
    "G01",
    "G02",
    "N03",
    "N04",
]


def cell(s):
    return str("—" if s is None or s == "" else s).replace("|", "\\|").replace("\n", " ")


def table(headers, rows):
    return (
        "\n".join(
            ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
            + ["| " + " | ".join(cell(v) for v in row) + " |" for row in rows]
        )
        + "\n"
    )


def pathlink(p, prefix=""):
    return f"[{p}]({prefix}{p})" if p else "None"


def literal_assignment(file, name):
    tree = ast.parse((ROOT / file).read_text(encoding="utf-8"))
    for n in tree.body:
        if isinstance(n, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == name for t in n.targets
        ):
            return ast.literal_eval(n.value)
    raise ValueError(name)


def validate(data):
    entries = data["methods"]
    ids = [m["id"] for m in entries]
    assert len(ids) == len(set(ids)), "Duplicate permanent ID"
    active = set()
    tree = ast.parse((ROOT / "src/multifat/models/registry.py").read_text(encoding="utf-8"))
    for fn in tree.body:
        if isinstance(fn, ast.FunctionDef) and fn.name in {
            "_available_cp_methods",
            "_available_direct_models",
        }:
            for n in ast.walk(fn):
                if isinstance(n, ast.Return) and isinstance(n.value, ast.Dict):
                    active.update(ast.literal_eval(k) for k in n.value.keys)
    disabled = set(literal_assignment("src/multifat/models/references.py", "DISABLED_METHODS"))
    impl = {m["registry_name"]: m for m in entries if m["record_type"] == "implementation"}
    assert set(impl) == active | disabled, "Catalogue/registry drift"
    for m in entries:
        assert m["implementation_status"] in STATUSES
        assert m["feasibility"] in "ABCD"
        assert m["validation_status"] in {
            "not_tested",
            "basic_numerical_tests",
            "reference_examples_reproduced",
            "independently_validated",
        }
        assert m["canonical_id"] is None or m["canonical_id"] in ids
        assert not m["canonical_id"] or not m["count_in_progress"]
        assert not m["doi"] or m["doi"].startswith("https://doi.org/")
        for p in ([m["module"]] if m["module"] else []) + m["test_files"]:
            assert (ROOT / p).is_file(), f"Missing evidence path: {p}"
        if m["registry_name"]:
            assert m["enabled"] == (m["registry_name"] in active)
        if m["implementation_status"] == "fully_implemented":
            assert (
                m["enabled"]
                and m["module"]
                and m["test_files"]
                and m["loading_scope"]
                and m["limitations"]
            )
            assert m["completeness_evidence"] == "complete_for_declared_repository_scope"
    coverage = list(
        csv.DictReader(
            (ROOT / "docs/research/Review_51_Coverage.csv").open(encoding="utf-8", newline="")
        )
    )
    assert sorted(int(r["number"]) for r in coverage) == list(range(1, 52))
    metadata = json.loads((ROOT / "docs/research/doi-records.json").read_text(encoding="utf-8"))
    dois = {r["DOI"].lower() for r in metadata}
    for m in entries:
        if "saved metadata in" in m["reference_status"]:
            assert m["doi"].removeprefix("https://doi.org/").lower() in dois
    inherited = list(
        csv.DictReader(
            (ROOT / "docs/research/Multiaxial_Fatigue_Master_Catalogue.csv").open(
                encoding="utf-8", newline=""
            )
        )
    )
    assert set(ids) == {
        r["ID"] for r in inherited
    }, "Inherited research record lost or added without reconciliation"
    ranks = [m["roadmap_rank"] for m in entries if m["roadmap_rank"] is not None]
    assert sorted(ranks) == list(range(1, len(ranks) + 1)), "Duplicate or missing roadmap rank"
    for m in entries:
        assert m["development_reason"] and m["reference_action"]
        assert (m["roadmap_rank"] is not None) == (m["development_track"] == "prioritized_addition")
        if m["roadmap_rank"] is not None:
            assert m["count_in_progress"] and m["implementation_status"] == "not_implemented"
    sources = {s["id"]: s for s in data["external_sources"]}
    assert len(sources) == len(data["external_sources"]), "Duplicate external source ID"
    searches = json.loads(
        (ROOT / "docs/research/external-code-search.json").read_text(encoding="utf-8")
    )
    assert {q["id"] for q in searches["queries"]} == {
        m["id"] for m in entries if m["count_in_progress"]
    }
    for m in entries:
        assert set(m["external_source_ids"]) <= set(sources)
        assert m["current_project_scope"]
        for source_id in m["external_source_ids"]:
            assert m["id"] in sources[source_id]["method_ids"]
        linked = [sources[sid] for sid in m["external_source_ids"]]
        expected = (
            "confirmed_source_available"
            if any(s["discovery_status"] == "confirmed" for s in linked)
            else (
                "possible_or_related_source"
                if linked
                else (
                    "no_accessible_source_identified"
                    if m["count_in_progress"]
                    else "not_individually_searched"
                )
            )
        )
        assert m["external_code_search_status"] == expected, "External discovery status drift"
        if m["count_in_progress"]:
            assert m["external_code_search_status"] != "not_individually_searched"
    return entries


def overview(entries):
    counted = [m for m in entries if m["count_in_progress"]]
    rows = []
    for family in sorted({m["primary_family"] for m in counted}):
        c = collections.Counter(
            m["implementation_status"] for m in counted if m["primary_family"] == family
        )
        rows.append([family, *[c[s] for s in STATUSES], sum(c.values())])
    c = collections.Counter(m["implementation_status"] for m in counted)
    rows.append(["**TOTAL**", *[c[s] for s in STATUSES], len(counted)])
    return table(
        ["Primary family", "Fully implemented¹", "Partial", "Missing", "Requires review", "Total"],
        rows,
    )


def catalogue(entries, detailed=False, family_filter=None):
    txt = "# WG4.3 method catalogue\n\nGenerated by `python tools/generate_method_docs.py`; edit [methods.json](docs/catalogue/methods.json), then regenerate. The JSON is the current source of truth. The CSV is a generated export.\n\n"
    txt += "## Status and scope\n\nThis catalogue supports preparation of the FABER WG4.3 scientific review article and collaborative computational comparison. The present phase is CA multiaxial fatigue of plain/smooth geometries, without explicit notches, cracks, stress concentrations or significant gradients. Wider geometry/loading entries remain future/reference tracks, not current capabilities or commitments. See the independent [external-code inventory](docs/EXTERNAL_IMPLEMENTATIONS.md); external availability never changes WG4.3 completeness or validation.\n\n**Fully implemented** means essential components are present, integrated and covered by numerical tests **for the named repository scope**. Some names explicitly identify repository adaptations; this does not certify reproduction of the unrestricted original criterion. **Partially implemented** means code exists with known essential gaps and may be enabled for a documented stage. **Not implemented** means no usable corresponding implementation exists. **Requires review** means active code exists but source completeness cannot yet be established; it is not silently counted as complete or as a proven coding omission.\n\n"
    txt += "Validation is independent: `not_tested`, `basic_numerical_tests`, `reference_examples_reproduced`, `independently_validated`. No current entry claims published-dataset reproduction or independent experimental validation. CAIM’s independent calculation of supplied equations is numerical verification. An execution sweep is not scientific validation.\n\n"
    txt += "CA = constant amplitude; P/NP = proportional/non-proportional. Repeating arbitrary paths are within the research scope; current workbook import supports a single-frequency harmonic cycle. Life from a fitted damage curve is distinct from a published native life law. [Audit](docs/method-audit.md) records shear conventions, clipping, material-property approximations and `_ext` variants. `_ext` is a repository maximum-damage selection, not another counted theory.\n\n"
    txt += (
        "## Progress\n\n"
        + overview(entries)
        + "\n¹Complete only for the explicitly declared scope. Primary families are counted once; secondary tags do not change totals. Counts refer to identifiable method/scoped-candidate records, not a proof of independent theoretical novelty.\n\n"
    )
    excluded = [m for m in entries if not m["count_in_progress"]]
    txt += f"All {len(entries)} inherited records remain traceable. {len(excluded)} bundles, overlap checks, related publications, unresolved leads and broad workstreams are excluded from progress totals. Absence from those totals is not a scientific rejection. The [51-publication coverage register](docs/research/Review_51_Coverage.csv) preserves every supplied review entry.\n\n"
    txt += "Browse [contribution opportunities and status/loading/validation views](docs/METHOD_INDEX.md), or inspect the [complete count reconciliation and record dispositions](docs/COUNT_RECONCILIATION.md). Roadmap selection never removes an entry here.\n\n"
    txt += "## Browse by family\n\n"
    families = sorted({m["primary_family"] for m in entries})
    if family_filter is not None:
        families = [family_filter]
    for f in families:
        txt += f'- [{f}](#family-{f.lower().replace(" / ", "-").replace(" ", "-")})\n'
    for f in families:
        anchor = f.lower().replace(" / ", "-").replace(" ", "-")
        txt += f'\n<a id="family-{anchor}"></a>\n\n## {f}\n\n'
        group = [m for m in entries if m["primary_family"] == f]
        txt += table(
            ["ID", "Method / scope", "Implementation", "Validation", "Priority", "Counted"],
            [
                [
                    f'[{m["id"]}](#{m["id"].lower()})',
                    m["name"],
                    LABEL[m["implementation_status"]],
                    m["validation_status"],
                    m["priority"],
                    "Yes" if m["count_in_progress"] else m["record_type"],
                ]
                for m in group
            ],
        )
        for m in group:
            txt += f'\n### {m["id"]}\n\n**{m["name"]}** — {m["authors_year"]}\n\n'
            details = [
                (
                    "Primary / secondary family",
                    m["primary_family"] + "; " + ", ".join(m["secondary_classifications"]),
                ),
                ("Record / counting", m["record_type"] + "; " + m["count_reason"]),
                ("Reference", m["reference"]),
                ("DOI", m["doi"] or "No verified original DOI available"),
                ("Reference evidence", m["reference_status"]),
                ("Formulation evidence", m["formulation_evidence"]),
                ("Principle", m["principle"]),
                ("Loading scope", m["loading_scope"]),
                ("Inputs / material parameters", m["required_inputs"]),
                ("Output", m["prediction_capability"]),
                (
                    "Implementation / evidence",
                    LABEL[m["implementation_status"]] + "; " + m["completeness_evidence"],
                ),
                ("Validation / evidence", m["validation_status"] + "; " + m["validation_evidence"]),
                ("Python implementation", pathlink(m["module"]) + "; " + m["symbol"]),
                (
                    "Automated tests",
                    ", ".join(pathlink(p) for p in m["test_files"]) or "None for this candidate",
                ),
                ("Limitations", m["limitations"]),
                ("Current FABER phase", m["current_project_scope"]),
                ("External source discovery", external_links(m)),
                ("Missing tasks", m["missing_tasks"]),
                ("Reference retrieval action", m["reference_action"]),
                (
                    "Development disposition",
                    m["development_track"] + "; " + m["development_reason"],
                ),
                ("Roadmap rank", m["roadmap_rank"]),
                ("Reuse", m["reuse"]),
                (
                    "Feasibility / effort / priority",
                    m["feasibility"] + " / " + m["complexity"] + " / " + m["priority"],
                ),
                ("Provenance / prior status", m["source_record_id"] + "; " + m["review_status"]),
                (
                    "Related catalogue ID",
                    (
                        f'[{m["canonical_id"]}](#{m["canonical_id"].lower()})'
                        if m["canonical_id"]
                        else "None recorded"
                    ),
                ),
                ("GitHub Issue / prepared brief", issue_links(m)),
            ]
            if not detailed:
                page = (
                    "docs/methods/"
                    + m["primary_family"].lower().replace(" / ", "-").replace(" ", "-")
                    + ".md"
                )
                details = [
                    (
                        "Reference / DOI",
                        m["reference"] + "; " + (m["doi"] or "Original DOI unverified"),
                    ),
                    (
                        "WG4.3 completeness / validation",
                        LABEL[m["implementation_status"]] + "; " + m["validation_status"],
                    ),
                    ("Current FABER phase", m["current_project_scope"]),
                    ("External sources", external_links(m)),
                    (
                        "Complete scientific and implementation record",
                        f'[Inputs, formulation, scope, code, tests, validation, missing work and Issues]({page}#{m["id"].lower()})',
                    ),
                ]
            txt += table(["Field", "Evidence / task"], details)
    return txt


def issue(m):
    source = "../../"
    return f"""# {m['id']}: {m['name']}

**Prepared Issue body; not an open GitHub Issue.** Before posting, search open and closed Issues for the permanent ID and aliases. Do not create a duplicate. Suggested labels: `method:{m['primary_family'].lower().replace(' ', '-')}`, `status:{'missing' if m['implementation_status']=='not_implemented' else 'review' if m['implementation_status']=='requires_review' else 'partial'}`, `help-wanted`, `validation-needed`. No `good-first-issue` claim is made.

## Scope and scientific reference

- Primary family: {m['primary_family']}; secondary: {', '.join(m['secondary_classifications']) or 'none recorded'}.
- Status: {LABEL[m['implementation_status']]}; enabled: {m['enabled']}.
- Authors/year: {m['authors_year']}.
- Reference: {m['reference']}.
- DOI: {m['doi'] or 'Original DOI not verified; retrieve the traceable primary source before coding.'}
- Reference evidence: {m['reference_status']}.
- Mathematical evidence: {m['formulation_evidence']}.
- Principle: {m['principle']}.
- Loading scope: {m['loading_scope']}.

## Available components and remaining work

FABER phase: {m['current_project_scope']}

External source discovery: {external_links(m, '../')}. Consult the linked inventory for exact variant, language, license, evidence and reuse restrictions. External code does not count as integrated WG4.3 implementation.

{m['missing_tasks']}

Existing implementation: {pathlink(m['module'], source)}; symbol `{m['symbol'] or 'not implemented'}`. Existing tests: {', '.join(pathlink(p, source) for p in m['test_files']) or 'none for this candidate'}.

Reuse: {m['reuse']}

Expected inputs: {m['required_inputs']}. Expected output: {m['prediction_capability']}.

Known restrictions: {m['limitations']}

## Acceptance criteria

- [ ] Identify exact original formulation, version, supplements and errata; record equations and departures. If source identity/equivalence is unresolved, resolve it before implementation.
- [ ] Agree bounded CA scope: proportional, phase-shifted NP, mean stresses and applicable material domain. Reject unsupported inputs explicitly.
- [ ] Implement missing components faithfully using physical simultaneous tensors and the original amplitude/plane rule; declare material fields, units, calibration and missing-data behavior.
- [ ] Integrate using [the existing interfaces](../../docs/adding-methods.md), registry and runtime source metadata. Do not enable an archived proxy before repair.
- [ ] Check axial/torsion limits, zero/degenerate cases, means, phase/time-origin and frame invariance as applicable; demonstrate grid/quadrature/optimizer convergence where needed.
- [ ] Obtain at least one independently calculated original-paper numerical case with material/loading data and intermediate quantities. Keep calibration and validation separate; record inaccessible data rather than inventing it.
- [ ] Preserve BP/CAIM, disabled guards, reporting and all-slope `log10(DP_e/DP)` behavior; run `MPLBACKEND=Agg python -m unittest discover -s tests -v`.
- [ ] Update `docs/catalogue/methods.json`, attach Issue/PR and validation evidence, regenerate documentation and pass `python tools/generate_method_docs.py --check`.
- [ ] Review implementation completeness separately from validation; never promote validation merely because numerical tests pass.

Validation lead: {m['validation_evidence']}

Dependencies: confirm source/version and measured parameters first; consult [ROADMAP.md](../../ROADMAP.md) for shared cycle, quadrature, native-life and data prerequisites. G01/G02 share a paper but retain distinct laws; other candidate dependencies require specification, not automatic coupling.

Use separate reviewable PRs for source/benchmark specification, prerequisites and numerical implementation when the scope warrants it. Repository entry: [METHODS.md#{m['id'].lower()}](../../METHODS.md#{m['id'].lower()}).
"""


def roadmap(entries):
    byid = {m["id"]: m for m in entries}
    roadmap_ids = [
        m["id"]
        for m in sorted(
            (e for e in entries if e["roadmap_rank"] is not None), key=lambda e: e["roadmap_rank"]
        )
    ]
    txt = "# Suggested implementation roadmap\n\nThese are research recommendations, not formally agreed WG4.3 commitments. No new fatigue law is implemented by this documentation change. Preserve the audit baseline and select one original formulation per PR.\n\n## Existing-code evidence and restoration\n\nFirst resolve LIU2021 Eq. (27), then MGSE_ZHU, MKBM, LI and scalar MSWT scope. Their code is active but completeness remains unestablished. Restore CARSPA plane orientation, Dang Van simultaneous shakedown/shear–hydrostatic history, the bounded OTT reduction/full-model distinction, physically paired SWTD directions and the actual ZHU_EDP work/path model in separate tasks. Prepared briefs specify evidence and acceptance criteria.\n\n"
    txt += roadmap_explanation(entries)
    txt += "## FABER WG4.3 phase and external-code prerequisites\n\nThis programme supports a review article in preparation. Near-term implementation targets plain/smooth CA fatigue, P/NP only where the chosen law supports it. Gradient, volume, notch, structural and VA/random targets below are retained as potential future/reference tracks, not current capabilities or formal commitments. Each row now states that boundary.\n\nBefore writing a missing law from scratch, inspect [external sources](docs/EXTERNAL_IMPLEMENTATIONS.md): openfat may help Crossland/Sines; FLife provides a COIN-LiWI amplitude kernel; ASU shares G35 model/case scripts, but no reuse grant was identified. Source availability changes preparation tasks, not scientific validation or priority automatically. FS-Sol/FI-Sol can help benchmark existing proportional branches. ECP remains future geometry work. No licenses or third-party code are adopted here.\n\n"
    txt += "## Ordered programme across families\n\n"
    txt += table(
        [
            "Order",
            "ID / method",
            "Added capability / principle",
            "Exact paper to retrieve",
            "DOI",
            "Prerequisite / missing work",
        ],
        [
            [
                i,
                f'[{id}](METHODS.md#{id.lower()}) — {byid[id]["name"]}',
                byid[id]["development_reason"],
                byid[id]["reference"],
                byid[id]["doi"] or "Original DOI unverified; use traceable citation",
                byid[id]["current_project_scope"]
                + " "
                + byid[id]["missing_tasks"]
                + " "
                + external_links(byid[id]),
            ]
            for i, id in enumerate(roadmap_ids, 1)
        ],
    )
    txt += "\n**First new-discovery paper:** Papuga & Růžička (2008), *Two new multiaxial criteria for high cycle fatigue computation*, [10.1016/j.ijfatigue.2007.02.015](https://doi.org/10.1016/j.ijfatigue.2007.02.015), for G01/G02. Keep the 2008 laws distinct from later PCRN/PIN/QCP versions. Consult their saved relationship notes; the 2019 mean-stress development and erratum are separate version requirements, not substitutions for the original paper.\n\n"
    txt += "## Reusable prerequisites\n\n"
    txt += table(
        ["Component", "Existing reuse", "Beneficiaries / boundary"],
        [
            [
                "Physical sampled repeating cycle",
                "Mean/sine/cosine tensors; direct model interface",
                "Periodic nonsinusoidal paths, energy integration and simultaneous invariant evaluation. No arbitrary-cycle importer currently exists.",
            ],
            [
                "Explicit source-specific amplitude rules",
                "BP major-axis helpers and legacy CP RSS convention",
                "Papuga, Papadopoulos, invariant/path methods; RSS and MCC must not be silently interchanged.",
            ],
            [
                "Orientation quadrature and convergence",
                "BP sphere grid and projections",
                "PI/PIN, Robert–Fogue–Bahuaud, Liu–Zenner; source-specific weights required.",
            ],
            [
                "Native fatigue resistance/life result",
                "Current empirical fitting only where mathematically appropriate",
                "MDC, MWCM, MMCCM, Morel and energy-life laws; do not replace native equations with one generic fit.",
            ],
            [
                "Prescribed stabilized strain response",
                "Tensor projection and plane search",
                "Energy/strain criteria; elastic workbook strains do not reproduce plastic work or NP hardening.",
            ],
            [
                "Gradient/spatial fields with units and weights",
                "Direct evaluator can be extended; no mesh/volume processor",
                "G29 bounded gradient input first; volumetric/nonlocal methods subsequently.",
            ],
            [
                "Probability calibration with run-outs",
                "Raw workbook provenance; processed importer filters run-outs",
                "Probability/weakest-link criteria; no claim of probabilistic validity from current point fits.",
            ],
        ],
    )
    txt += "\nThe full catalogue retains every other A/B candidate and wider C workstream. Bibliographic title/DOI identity does not establish reproducible equations. Conditional candidates and unresolved overlaps need source review before assignment. Future benchmarks should provide units, cyclic conventions, material coefficients, actual paths and intermediate results.\n"
    return txt


def external_links(m, doc_prefix="docs/"):
    state = m["external_code_search_status"]
    links = ", ".join(
        f"[{sid}]({doc_prefix}EXTERNAL_IMPLEMENTATIONS.md#{sid.lower()})"
        for sid in m["external_source_ids"]
    )
    return state + ("; " + links if links else "") + "; not integrated or independently validated"


def external_inventory(entries):
    counted = [m for m in entries if m["count_in_progress"]]
    stats = collections.Counter(m["external_code_search_status"] for m in counted)
    txt = "# Public external implementations and code leads\n\nGenerated from the external-source records in [catalogue/methods.json](catalogue/methods.json). Search/check date: 15 October 2026. Availability and limited static inspection are not equation verification, successful execution, predictive validation or permission to reuse. No third-party code was copied, translated, executed or integrated.\n\n"
    txt += "## Search coverage and meaning\n\n"
    txt += table(
        [
            "Counted catalogue method/scoped records searched",
            "Confirmed corresponding/bounded-variant source",
            "Possible or related source",
            "No accessible source identified",
        ],
        [
            [
                len(counted),
                stats["confirmed_source_available"],
                stats["possible_or_related_source"],
                stats["no_accessible_source_identified"],
            ]
        ],
    )
    txt += "\nThe 117 records are not a definitive count of independent theories. Two additional future workstreams, C04/ECP and C09/random-spectrum HCF, have confirmed author-associated source links and are reported separately, not added to the 117-method denominator. Other noncounted leads were not individually searched.\n\n"
    txt += "Confirmed means an accessible file and a traceable, explicitly bounded method correspondence were inspected. Possible includes related variants, supporting helpers, a command interface without numerical-kernel inspection, documented software without public source, and an inaccessible archive lead. One record can have several sources; confirmed takes precedence for the per-method count. No source found is a limited search outcome, not proof of absence. Source-code availability never changes WG4.3 implementation or validation status.\n\n"
    txt += "## License and scope gate\n\nSix source projects offer potentially useful licensed current-scope implementations or supporting helpers: FS-Sol, FI-Sol, openfat, FLife, Fatlab and fatpack. ECP is a seventh licensed project for future geometry/gradient work. These are project counts, not numbers of independently validated laws. GPLv3/GPLv2-or-later, MIT and ISC have different obligations. This repository is licensed MIT: MIT, BSD, ISC and Apache-2.0 code can be combined with it, provided notices are kept and the exact file scope is recorded; GPL code cannot be incorporated and may serve only as a reference for an independent implementation. No integration is pre-approved. Unknown-license ASU repositories are reference-only pending permission. Public visibility or a copyright/citation notice is not a reuse grant.\n\n"
    txt += "The current FABER phase is CA fatigue of plain/smooth geometries without explicit notches, cracks or significant gradients. Existing FE wrappers can sometimes supply point kernels, but that does not import geometry processing into this scope. FLife’s COIN-LiWI kernel uses complex amplitude spectra; its spectral Carpinteri–Spagnoli variant is not a replacement for the original CA plane construction. ASU’s 2006 strain-based and 2020 energy sources must not be substituted for differently cited catalogue laws.\n\n"
    txt += "## Source inventory\n\n"
    for s in EXTERNAL_SOURCES:
        txt += f"### {s['id']}\n\n**{s['project']}** — {s['discovery_status']}\n\n"
        fields = [
            (
                "Catalogue correspondence",
                ", ".join(f"[{mid}](../METHODS.md#{mid.lower()})" for mid in s["method_ids"])
                or "Unmapped related paper; no new criterion promoted",
            ),
            ("Project", s["url"]),
            ("Direct source / lead", s["source_url"]),
            (
                "Checked revision",
                s["revision"]
                or "Page/branch checked on stated date; not a pinned numerical release",
            ),
            ("Language / dependencies", s["language"] + "; " + s["dependencies"]),
            ("License / evidence", s["license"] + "; " + s["license_evidence"]),
            ("Formulation / loading / geometry", s["correspondence"]),
            ("Tests / scientific evidence", s["quality"]),
            ("Potential use and limitations", s["potential_use"]),
            ("Verification boundary", s["verification"]),
        ]
        txt += table(["Field", "Evidence / restriction"], fields) + "\n"
    txt += "## Method-by-method discovery register\n\nThis covers every counted record plus the two source-linked future workstreams. Follow the catalogue for original reference, parameters, WG4.3 implementation and validation.\n\n"
    txt += table(
        ["Method ID / name", "Search disposition", "External links", "Current phase"],
        [
            [
                f"[{m['id']}](../METHODS.md#{m['id'].lower()}) — {m['name']}",
                m["external_code_search_status"],
                external_links(m, ""),
                m["current_project_scope"],
            ]
            for m in entries
            if m["count_in_progress"] or m["id"] in {"C04", "C09"}
        ],
    )
    txt += "\n## Search method and unresolved leads\n\nThe [query register](research/external-code-search.json) records one named source-code query per counted record, with targeted stronger searches and citation/repository chaining for useful leads. GitHub, GitLab, MathWorks File Exchange, SourceForge, author/research sites, Zenodo and publication-linked software were screened. Useful links were checked against primary repository metadata, actual files, documentation and license notices. Generic matches and homonyms were rejected. This is a targeted first pass, not exhaustive code indexing or a saturation claim.\n\nPragTic provides valuable author documentation but public solver source/reuse rights were not established. A historical source-on-request statement is not an open-source license. The inspected FatPy repository was empty: an implementation Issue is not source code. MathWorks LoadPathMATLAB visualizes FE load flow, not the fatigue Moment of Load Path law; sports/physiological fatigue repositories were rejected. Zenodo modified-Findley archive files could not be opened and remain unverified. The current official Code_Aster GitLab page returned a loading shell; a historical GitHub mirror exposed the command interface, not an audited standalone numerical kernel.\n\nRecommendations: first compare proportional FS/Findley against the author MATLAB algorithms; specify source-faithful Crossland/Sines with openfat as a checked-but-unvalidated comparator; investigate the MIT COIN-LiWI kernel and its wrapper independently; ask for permission before porting ASU G35 or its data. Preserve ECP for a future geometry extension. Do not change method validation/completeness solely from this inventory. No live method Issues exist; relevant prepared briefs now point here, and actual Issues should inherit these references when posted.\n"
    return txt


def selected(entries):
    return [m for m in entries if m["roadmap_rank"] is not None]


def issue_links(m, prefix=""):
    if m["issue_url"]:
        return f'[Open Issue]({m["issue_url"]})'
    search = "https://github.com/andreachiocca/multiaxial-fatigue/issues?" + urlencode(
        {"q": f'{m["id"]} in:title,body'}
    )
    proposal = "https://github.com/andreachiocca/multiaxial-fatigue/issues/new?" + urlencode(
        {"template": "method.md", "title": f'[{m["id"]}] {m["name"]}'}
    )
    bits = [f"[Search Issues]({search})", f"[Propose task]({proposal})"]
    if m["id"] in ISSUE_IDS:
        bits.insert(0, f'[Prepared brief]({prefix}docs/contribution-issues/{m["id"]}.md)')
    return " · ".join(bits) + " (no live Issue linked)"


def roadmap_explanation(entries):
    missing = [
        m
        for m in entries
        if m["count_in_progress"] and m["implementation_status"] == "not_implemented"
    ]
    chosen = selected(entries)
    deferred = [m for m in missing if m["roadmap_rank"] is None]
    return f"""## Why the roadmap is shorter than the catalogue

This roadmap selects **{len(chosen)} of {len(missing)} missing method/scoped candidates** for an initial programme. **{len(deferred)} remain outside the sequence**: {sum(m['feasibility'] in 'AB' for m in deferred)} with preliminary A/B feasibility and {sum(m['feasibility']=='C' for m in deferred)} requiring major extensions. They remain contribution opportunities, including medium/low-priority candidates. The five existing repairs and five existing-source reviews above form separate tracks.

Selection balances recognized CA fatigue relevance, additional scientific capability, methodological diversity, reuse of current utilities, reference/benchmark leads, material data and shared prerequisites. Numerical simplicity is not the sole criterion: volumetric, mesoscopic and probabilistic extensions are deliberately retained. Order is advisory and may change when original papers/data establish better feasibility.

The earlier research shortlist had 23 entries. G14 (Glinka–Shen–Plumtree) was subsequently held for overlap review against GSE/GSA, leaving 22 in this sequence; it remains in the full catalogue and is not asserted equivalent.

No candidate is rejected just because it is absent here. Deferred rows have explicit requirements in the [complete per-record reconciliation](docs/COUNT_RECONCILIATION.md) and [contribution index](docs/METHOD_INDEX.md). For work outside the sequence, choose its permanent ID, check source identity and related criteria, then open a scoped proposal using the method template. A different sequence is welcome when justified by a source-faithful formulation and available validation data.

Full VA/random-only workflows are deferred from CA priorities. Bundles and uncertain identities need specification first; application/comparison papers are not automatically new laws. Those dispositions do not remove their source records from the catalogue.

"""


def reconciliation(entries):
    counted = [m for m in entries if m["count_in_progress"]]
    missing = [m for m in counted if m["implementation_status"] == "not_implemented"]
    chosen = selected(entries)
    outside = [m for m in missing if m["roadmap_rank"] is None]
    txt = "# Catalogue and roadmap count reconciliation\n\nGenerated from `catalogue/methods.json`; every inherited permanent ID is checked against the frozen research CSV.\n\n"
    txt += f"## Answer to “about 100 reviewed, but about 20 in the roadmap?”\n\nThe **{len(missing)} missing method/scoped candidates** describe the implementation backlog. The **{len(chosen)} roadmap entries** are an intentionally selected initial development sequence. The other **{len(outside)}** are retained, individually accounted for and open to proposals. The prior 23-step research shortlist was itself selective; the published initial roadmap had 22 additions. Neither was the literature census. The only difference between the earlier 23-entry shortlist and the 22-entry published sequence is **G14, Glinka–Shen–Plumtree**: it remains in the catalogue as an overlap-review record pending an equation-level comparison with GSE/GSA. This is not proof that the laws are equivalent. The missing explanation is corrected here.\n\n"
    txt += "There is no evidence of an inherited record being lost: all **161/161 IDs** survive, and all **51/51 numbered FABER publications** have dispositions. The weakness was that the selection policy and the nonselected backlog were not prominent enough. This update adds explicit development dispositions and browse links; it does not fabricate additional fatigue laws.\n\n"
    txt += "## Exact inventory arithmetic\n\n"
    prefix_rows = []
    for prefix, desc in [
        ("I_", "Current and archived implementation records"),
        ("N", "Inherited review-derived candidate packages"),
        ("G", "Extended-search publication/formulation packages"),
        ("R", "Restoration/variant notes"),
        ("C", "Wider workflows/model/application packages"),
        ("L", "Historical bibliographic/identity leads"),
    ]:
        group = [m for m in entries if m["id"].startswith(prefix)]
        prefix_rows.append(
            [
                prefix,
                desc,
                len(group),
                sum(m["count_in_progress"] for m in group),
                sum(not m["count_in_progress"] for m in group),
            ]
        )
    txt += table(
        [
            "ID prefix",
            "What it represents",
            "All records",
            "In method progress total",
            "Outside denominator",
        ],
        prefix_rows + [["TOTAL", "", len(entries), len(counted), len(entries) - len(counted)]],
    )
    txt += f'\n**161 = 117 + 44; 117 = 20 implementation records + 97 missing candidates; 97 = {len(chosen)} prioritized additions + {len(outside)} nonselected candidates.** Of the 97 missing candidates, {sum(m["feasibility"] in "AB" for m in missing)} have preliminary A/B feasibility and {sum(m["feasibility"]=="C" for m in missing)} have C feasibility. Feasibility does not certify available measured data, equations or a reproducible algorithm.\n\n'
    txt += "## What these numbers do and do not count\n\nA publication can introduce several laws (G01/PCR and G02/PI share one 2008 article), compare existing laws, extend one family or supply an algorithm. Different papers can describe the same formulation. N14 contains Gough–Pollard and Nishihara–Kawamoto and is a bundle, not one reconciled independent law. R01/R02/R03 refer to existing archived implementations and are linked to their canonical IDs. `_ext` variants are documented under their parent implementation. Meaningful later versions remain distinguished in scope/relationship notes; unresolved equivalence is not asserted.\n\n"
    txt += "**117 is a reconciled count of identifiable method/scoped-candidate records, not a scientifically final count of independent theories.** Exact independence remains unresolved for some bundles, variants and broad C workstreams; expanding those could increase the denominator, while finding equivalences could reduce it. Resolving that would require original-paper work beyond this documentation task. All candidates remain visible irrespective of difficulty or family.\n\n"
    c = collections.Counter(m["record_type"] for m in entries if not m["count_in_progress"])
    txt += table(["Outside progress denominator", "Records"], sorted(c.items()))
    txt += "\nBroader workstreams are excluded from the narrow method-progress denominator because their boundaries are not reconciled, not because major extensions are unwanted. Twelve **counted** missing candidates already have C feasibility, and three are on the initial roadmap. Full VA/random workflows are outside the CA priority scope. Nguyen’s manufacturing comparison does not establish an independent law. The per-record reasons below preserve these distinctions.\n\n"
    tracks = collections.Counter(m["development_track"] for m in entries)
    txt += "## Development dispositions\n\n" + table(["Track", "Records"], sorted(tracks.items()))
    txt += "\n## Every record accounted for\n\nUse the catalogue link for full references, loading/data scope and scientific limitations. “Outside denominator” is not a scientific exclusion. Development dispositions and reasons are fields in the JSON and generated CSV.\n\n"
    txt += table(
        [
            "ID / name",
            "Record type",
            "Counted",
            "Implementation",
            "Development disposition",
            "Selection / deferral / relationship reason",
        ],
        [
            [
                f'[{m["id"]}](../METHODS.md#{m["id"].lower()}) — {m["name"]}',
                m["record_type"],
                "Yes" if m["count_in_progress"] else "No",
                LABEL[m["implementation_status"]],
                f'#{m["roadmap_rank"]}' if m["roadmap_rank"] else m["development_track"],
                m["development_reason"],
            ]
            for m in entries
        ],
    )
    txt += "\nThe frozen [research checkpoint](research/Multiaxial_Fatigue_Research_Checkpoint.md), [original 51-publication review](research/Review-FABER-WG43.md), [coverage register](research/Review_51_Coverage.csv) and [prior assessment](research/feasibility-assessment.md) remain source evidence. Some inherited workstream numbering is inconsistent and remains flagged; it is not promoted to a verified method identity. No new broad search or original-equation verification was performed.\n"
    return txt


def readme_reconciliation(entries):
    counted = [m for m in entries if m["count_in_progress"]]
    missing = [m for m in counted if m["implementation_status"] == "not_implemented"]
    n = len(selected(entries))
    return f"The catalogue contains **{len(missing)} missing method/scoped candidates**; the roadmap selects **{n}** for an initial, dependency-aware programme. The other **{len(missing)-n} remain available for proposals**. The complete index retains **{len(entries)} research records**, including **{len(counted)} identifiable method/scoped records** and **{len(entries)-len(counted)} bundles, related sources, unresolved overlaps/leads and wider workstreams** outside that progress denominator. These are not {len(entries)} independent fatigue theories.\n"


def readme_navigation(entries):
    missing = [
        m
        for m in entries
        if m["count_in_progress"] and m["implementation_status"] == "not_implemented"
    ]
    txt = f'**{len(missing)} missing method candidates** are available to explore: {sum(m["feasibility"] in "AB" for m in missing)} have preliminary A/B feasibility (a new criterion or manageable extensions). Pick by scientific interest and available references/data; an off-roadmap proposal is welcome.\n\n'
    families = sorted({m["primary_family"] for m in missing})
    txt += (
        "Browse contribution opportunities by family: "
        + " · ".join(
            f'[{f}](docs/METHOD_INDEX.md#family-{f.lower().replace(" / ","-").replace(" ","-")})'
            for f in families
        )
        + ".\n\n"
    )
    active = [m for m in entries if m["enabled"]]
    txt += (
        "**Available in the current API:** "
        + ", ".join(f'[{m["registry_name"]}](METHODS.md#{m["id"].lower()})' for m in active)
        + ". Each link states the implemented scope and evidence; active does not automatically mean source completeness is established.\n"
    )
    return txt


def browse(entries):
    counted = [m for m in entries if m["count_in_progress"]]
    txt = "# Find a method to contribute\n\nGenerated from [catalogue/methods.json](catalogue/methods.json). The [full catalogue](../METHODS.md) includes every research record; [count reconciliation](COUNT_RECONCILIATION.md) explains the denominator. A/B/C/D and priorities are preliminary scientific/engineering judgments, not promises about available equations or measured material data.\n\n"
    txt += "[Available API methods](#available-api-methods) · [Repairs](#repairs) · [Source reviews](#source-reviews) · [Missing methods by family](#missing-methods-by-family) · [Validation and benchmarks](#validation-and-benchmarks) · [Loading scope](#loading-scope) · [Unresolved/related work](#unresolved-and-related-work)\n\n"
    for title, group in [
        ("Available API methods", [m for m in entries if m["enabled"]]),
        ("Repairs", [m for m in counted if m["implementation_status"] == "partially_implemented"]),
        ("Source reviews", [m for m in counted if m["implementation_status"] == "requires_review"]),
    ]:
        txt += f"## {title}\n\n"
        txt += (
            table(
                ["ID / method", "Completeness", "Validation evidence", "Task / Issue"],
                [
                    [
                        f'[{m["id"]}](../METHODS.md#{m["id"].lower()}) — {m["name"]}',
                        LABEL[m["implementation_status"]],
                        m["validation_status"] + "; " + m["validation_evidence"],
                        issue_links(m, "../"),
                    ]
                    for m in group
                ],
            )
            + "\n"
        )
    txt += "## Missing methods by family\n\nThese are all missing counted candidates, not only the prioritized roadmap. A numbered row is selected in the initial sequence; other rows remain open to proposals. Follow the catalogue link for the principle, required material parameters and exact outstanding tasks. DOI identity does not establish a fully reproducible algorithm.\n\n"
    missing = [m for m in counted if m["implementation_status"] == "not_implemented"]
    for family in sorted({m["primary_family"] for m in missing}):
        anchor = family.lower().replace(" / ", "-").replace(" ", "-")
        txt += f'<a id="family-{anchor}"></a>\n\n### {family}\n\n'
        txt += (
            table(
                [
                    "ID / method",
                    "Authors/year; reference",
                    "Loading scope / expected output",
                    "Feasibility; plan",
                    "Find or propose work",
                ],
                [
                    [
                        f'[{m["id"]}](../METHODS.md#{m["id"].lower()}) — {m["name"]}',
                        m["authors_year"]
                        + "; "
                        + (
                            f'[DOI]({m["doi"]})'
                            if m["doi"]
                            else "No verified original DOI; see catalogue citation"
                        ),
                        m["loading_scope"] + " / " + m["prediction_capability"],
                        m["current_project_scope"]
                        + " "
                        + external_links(m, "")
                        + " "
                        + m["feasibility"]
                        + "; "
                        + (
                            f'roadmap #{m["roadmap_rank"]}'
                            if m["roadmap_rank"]
                            else m["development_track"]
                        ),
                        issue_links(m, "../"),
                    ]
                    for m in missing
                    if m["primary_family"] == family
                ],
            )
            + "\n"
        )
    txt += "## Validation and benchmarks\n\n"
    txt += table(
        ["Evidence level", "Counted records", "What contributors can do"],
        [
            [level, sum(m["validation_status"] == level for m in counted), task]
            for level, task in [
                (
                    "not_tested",
                    "Implement or repair the law, then supply meaningful numerical checks. Disabled guards are not criterion tests.",
                ),
                (
                    "basic_numerical_tests",
                    "Supply independently calculated published cases, intermediate results and independent datasets.",
                ),
                (
                    "reference_examples_reproduced",
                    "Extend independent validation across material/loading domains.",
                ),
                (
                    "independently_validated",
                    "Extend evidence to new documented domains without changing the source formulation.",
                ),
            ]
        ],
    )
    txt += "\nExisting evidence is linked per active method above and in the catalogue. BP has analytical loading limits and phase checks; CAIM also has a separate dense calculation of supplied MATLAB equations. Those are numerical verification, not reproduced experimental datasets. Other active methods have the documented analytical/regression tests. The material workbooks are available for future comparisons; old result CSVs are historical and must not be treated as current validation.\n\n"
    txt += "## Loading scope\n\nCurrent executable input: physical single-frequency CA tensors `S_mean + S_sin sin(t) + S_cos cos(t)`, with strains inferred elastically by the workbook importer. API availability and each law’s validity domain are separate. A repeating NP path can be nonsinusoidal, but the current importer cannot represent arbitrary periodic waveforms.\n\nUse the loading/output columns above to screen candidate P/NP, mean-stress, strength-limit and life capabilities. These are source-domain assessments, not newly implemented capabilities. For proposed strain/plastic-work, gradient or probability use, inspect the required data and restrictions in the linked catalogue entry before coding.\n\n"
    txt += "## Unresolved and related work\n\nAll noncounted bundles, possible equivalences, wider physical models, applications and bibliographic leads have an individual [reconciliation row](COUNT_RECONCILIATION.md#every-record-accounted-for). They are valid reference/specification contribution opportunities. They need source/scope clarification before a claim that a new independent criterion has been implemented. They are not hidden rejections or restricted to critical-plane work.\n"
    return txt


def main():
    args = argparse.ArgumentParser()
    args.add_argument("--check", action="store_true")
    opt = args.parse_args()
    data = json.loads(DATA.read_text(encoding="utf-8"))
    entries = validate(data)
    global EXTERNAL_SOURCES
    EXTERNAL_SOURCES = data["external_sources"]
    outputs = {
        "METHODS.md": catalogue(entries),
        "ROADMAP.md": roadmap(entries),
        "docs/COUNT_RECONCILIATION.md": reconciliation(entries),
        "docs/METHOD_INDEX.md": browse(entries),
    }
    for family in sorted({m["primary_family"] for m in entries}):
        name = family.lower().replace(" / ", "-").replace(" ", "-")
        content = catalogue(entries, detailed=True, family_filter=family)

        def rebase(match):
            target = match.group(1)
            if target.startswith(("https:", "http:")):
                return match.group(0)
            return (
                "](" + ("../../METHODS.md" if target.startswith("#") else "../../") + target + ")"
            )

        outputs["docs/methods/" + name + ".md"] = re.sub(r"\]\(([^)]+)\)", rebase, content)
    outputs["docs/EXTERNAL_IMPLEMENTATIONS.md"] = external_inventory(entries)
    extbuf = io.StringIO(newline="")
    ew = csv.DictWriter(extbuf, fieldnames=list(EXTERNAL_SOURCES[0]), lineterminator="\n")
    ew.writeheader()
    ew.writerows(
        {
            k: json.dumps(v, ensure_ascii=False) if isinstance(v, list) else v
            for k, v in source.items()
        }
        for source in EXTERNAL_SOURCES
    )
    outputs["docs/catalogue/external_sources.csv"] = extbuf.getvalue()
    rows = []
    for m in entries:
        rows.append(
            {
                k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                for k, v in m.items()
            }
        )
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    outputs["docs/catalogue/methods.csv"] = buf.getvalue()
    index = "# Prepared contribution Issues\n\nThese versioned briefs are ready for maintainer review and posting. **No live Issues or labels have been created.** The connector available during preparation had no Issue/label creation operation. Search existing Issues by permanent ID and aliases before posting; add the returned Issue URL to `docs/catalogue/methods.json` and regenerate.\n\n"
    index += table(
        ["ID", "Scope", "Status", "Brief"],
        [
            [m["id"], m["name"], LABEL[m["implementation_status"]], f'[{m["id"]}]({m["id"]}.md)']
            for m in entries
            if m["id"] in ISSUE_IDS
        ],
    )
    index += "\nUse the repository Issue template for additional methods after checking reference identity and scope. Broad unresolved bundles should be specified before a numerical implementation Issue is opened. Suggested families/labels are advisory; no simple numerical task is automatically a `good-first-issue`.\n"
    outputs["docs/contribution-issues/README.md"] = index
    for m in entries:
        if m["id"] in ISSUE_IDS:
            outputs[f'docs/contribution-issues/{m["id"]}.md'] = issue(m)
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    before, rest = readme.split("<!-- METHOD_PROGRESS_START -->")
    _, after = rest.split("<!-- METHOD_PROGRESS_END -->")
    outputs["README.md"] = (
        before
        + "<!-- METHOD_PROGRESS_START -->\n\n"
        + overview(entries)
        + "\n<!-- METHOD_PROGRESS_END -->"
        + after
    )
    readme = outputs["README.md"]
    before, rest = readme.split("<!-- METHOD_NAVIGATION_START -->")
    _, after = rest.split("<!-- METHOD_NAVIGATION_END -->")
    outputs["README.md"] = (
        before
        + "<!-- METHOD_NAVIGATION_START -->\n\n"
        + readme_navigation(entries)
        + "\n<!-- METHOD_NAVIGATION_END -->"
        + after
    )
    readme = outputs["README.md"]
    before, rest = readme.split("<!-- METHOD_RECONCILIATION_START -->")
    _, after = rest.split("<!-- METHOD_RECONCILIATION_END -->")
    outputs["README.md"] = (
        before
        + "<!-- METHOD_RECONCILIATION_START -->\n\n"
        + readme_reconciliation(entries)
        + "\n<!-- METHOD_RECONCILIATION_END -->"
        + after
    )
    stale = []
    for path, content in outputs.items():
        p = ROOT / path
        if opt.check:
            if not p.exists() or p.read_text(encoding="utf-8") != content:
                stale.append(path)
        else:
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(content, encoding="utf-8", newline="\n")
    if stale:
        raise SystemExit("Stale generated files: " + ", ".join(stale))
    print(
        f'{len(entries)} records; {sum(m["count_in_progress"] for m in entries)} counted method/scoped candidates; 51/51 review entries; registry, paths, DOI evidence and generated files OK'
    )


if __name__ == "__main__":
    main()
