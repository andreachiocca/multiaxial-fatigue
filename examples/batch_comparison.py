"""Cross-batch statistics of every method on the batches of one material (issue #71).

A family of batches of one material is declared in ``examples/data/material_families.json``.
Every method is run twice on each batch:

- **own batch**: calibrated and assessed on the same batch (the ordinary pipeline run);
- **transfer**: the batch predicted from each other batch of the family
  (``multifat-pipeline --design-material``), with the Summary, constants, design set and master
  curve of the design batch.

The statistic is the load-equivalent error e = Error_load = ln(1/s*) of the assessed failures
(Flag 1), grouped by loading path (A, T, AT-IP, AT-NP, with "+mean" for a non-zero mean stress).
The analysis has a method-free part and a method part:

- method-free (``examples/batch_reference.py``): Basquin strength S(1e5) of each batch, path
  and test series, strength ratios between paths, and the scatter between the series of one
  steel in the FABEST database, the yardstick for a batch difference;
- per method, batch and path: n, bias = mean(e) with a 95 % percentile-bootstrap interval, SD,
  RMS and the share within +-20 % of the load;
- per path: additive decomposition of the cell biases,
  bias_mb = mu + alpha_m + beta_b + (alpha beta)_mb. The batch offsets beta_b get intervals from
  a joint bootstrap (the same specimens resampled for all methods), the interaction is compared
  with what sampling alone produces, and each offset is set against its method-free prediction;
- run-outs: the offsets recomputed with the run-outs as censored observations (sensitivity);
- transfer: bias and RMS per design -> target pair, the RMS change on the tests the own-batch run
  assesses (paired by Observation_id), the method-free prediction of the transfer bias, and the
  load factor that covers 95 % of the tests, own batch against transfer.

Effect sizes with intervals are reported, not p-values: the specimens of one S-N curve share
the lot, the machine and the laboratory, so a test at specimen level would only say that two
test series differ, which the intervals show already; the number of batches (two or three) is
too small for any inference on a population of batches.

The results are written to ``examples/data/batch_comparison/<family>/`` (CSV, LaTeX tables and a
manifest) and ``docs/images/batch_comparison_<family>_*.png``
(figures: ``examples/batch_plots.py``). The manifest records the registered models, the
implementation revision and the dataset fingerprints; ``check`` (run by
``tests/test_batch_comparison.py``) fails when a method or a workbook has changed since, so the
results are regenerated with every new method or data change:

    uv run python examples/batch_comparison.py all <tmp_dir> --jobs 4   # run + analyse
    uv run python examples/batch_comparison.py run <tmp_dir> --jobs 4   # pipeline runs only
    uv run python examples/batch_comparison.py analyse <tmp_dir>        # statistics and figures
    uv run python examples/batch_comparison.py check                    # committed results current?
"""

from __future__ import annotations

import argparse
import datetime
import json
import subprocess
import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.style import (  # noqa: E402
    BATCH_COLOURS,
    BATCH_MARKERS,
    DATA_DIR,
    DOCS_IMG_DIR,
    apply_academic_style,
)

FAMILIES_FILE = DATA_DIR / "material_families.json"
OUT_ROOT = DATA_DIR / "batch_comparison"
SEED = 20261110  # bootstrap seed, fixed for reproducible intervals
N_BOOT = 4000
MIN_N = 3  # smallest group entering a test or a decomposition
WITHIN = 0.2  # share of failures within +-20 % of the load
PATH_NAMES = {
    "Tensile": "A",
    "Torsion": "T",
    "Multiaxial proportional": "AT-IP",
    "Multiaxial non-proportional": "AT-NP",
}
PATH_ORDER = ["A", "A+mean", "T", "T+mean", "AT-IP", "AT-IP+mean", "AT-NP", "AT-NP+mean"]


# ---- families and inputs -------------------------------------------------------------------
def families() -> dict:
    data = json.loads(FAMILIES_FILE.read_text(encoding="utf-8"))
    return {k: v for k, v in data.items() if not k.startswith("_")}


def registered_models() -> list[str]:
    from multifat.models.registry import model_names

    return sorted(model_names())


def dataset_fingerprints(batches) -> dict:
    from multifat.cli import run as runmod
    from multifat.report.subsets import dataset_fingerprint

    out = {}
    for b in batches:
        data, mat = runmod.load_material(b)
        out[b] = dataset_fingerprint(data[mat])
    return out


def loading_paths(batches) -> pd.DataFrame:
    """Loading path of every failure and run-out, keyed by (Batch, Observation_id)."""
    from multifat.cli import run as runmod
    from multifat.report.subsets import observation_index

    rows = []
    for b in batches:
        data, mat = runmod.load_material(b)
        block = data[mat]
        index = observation_index(block)
        for kind, group in (("case", block["cases"]), ("runout", block.get("runouts", {}))):
            for cid, case in group.items():
                H = case["harmonics"]
                label = PATH_NAMES.get(
                    runmod.classify_loading_type(H=H, meta=case.get("meta", {})), "?"
                )
                if np.max(np.abs(np.asarray(H["S_mean"], float))) > 1e-6:
                    label += "+mean"
                rows.append(
                    dict(Batch=b, Observation_id=index[(kind, int(cid))]["Observation_id"],
                         Path=label)
                )  # fmt: skip
    return pd.DataFrame(rows)


# ---- pipeline runs -------------------------------------------------------------------------
def _pipeline(args: list[str]) -> int:
    cmd = [sys.executable, "-m", "multifat.cli.pipeline", *args, "--no-compare"]
    print("$", " ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=REPO_ROOT).returncode


def run(out: Path, family: str, jobs: int, models: str) -> None:
    """Own-batch runs, then every ordered transfer, into out/<family>/{own,transfer}.

    A batch listed under the family's ``constants`` (a published-constants record, e.g. the SSF
    surface of 42CrMo4_Reis) is run with ``--constants`` on its own and as a design batch; those
    constants then reach the other batches only through the transfer.
    """
    spec = families()[family]
    batches, constants = spec["batches"], spec.get("constants", {})
    base = out / family
    common = ["--models", models, "--jobs", str(jobs)]

    def materials(names):
        return [arg for b in names for arg in ("--material", b)]

    def supplied(batch):
        return ["--constants", str(REPO_ROOT / constants[batch])] if batch in constants else []

    plain = [b for b in batches if b not in constants]
    codes = []
    if plain:
        codes.append(_pipeline([*materials(plain), *common, "--results-dir", str(base / "own")]))
    for b in (b for b in batches if b in constants):
        codes.append(
            _pipeline([*materials([b]), *supplied(b), *common, "--results-dir", str(base / "own")])
        )
    for design in batches:
        others = [b for b in batches if b != design]
        codes.append(
            _pipeline(
                [*materials(others), "--design-material", design, *supplied(design), *common,
                 "--results-dir", str(base / "transfer")]
            )
        )  # fmt: skip
    if any(codes):  # a failed model is reported in the run reports; the others are usable
        print(f"warning: pipeline exit codes {codes}; see the *_run_report.tsv files")


def read_results(folder: Path, materials: dict[str, tuple[str, str | None]]) -> pd.DataFrame:
    """Rows of the base-variant CSVs of `materials` {material id: (target, design or None)}."""
    frames = []
    for path in sorted(folder.glob("*.csv")):
        head = pd.read_csv(path, nrows=1)
        if not {"Method", "Material", "Flag", "Error_load"} <= set(head.columns):
            continue
        material, method = str(head["Material"].iloc[0]), str(head["Method"].iloc[0])
        if material not in materials or method.endswith("_ext"):
            continue
        frame = pd.read_csv(path)
        target, design = materials[material]
        frame["Batch"], frame["Design"] = target, design or target
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    d = pd.concat(frames, ignore_index=True)
    d["Error_load"] = pd.to_numeric(d["Error_load"], errors="coerce")
    return d


# ---- statistics ----------------------------------------------------------------------------
def bootstrap_mean_ci(x: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95 % percentile-bootstrap interval of the mean (NaN below MIN_N values)."""
    if x.size < MIN_N:
        return (np.nan, np.nan)
    means = x[rng.integers(0, x.size, size=(N_BOOT, x.size))].mean(axis=1)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo), float(hi)


def cell_statistics(d: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """n, bias (with bootstrap interval), SD, RMS and share within +-20 % per group of `keys`."""
    rng = np.random.default_rng(SEED)
    rows = []
    for key, g in d.groupby(keys, sort=True):
        e = g["Error_load"].to_numpy(float)
        e = e[np.isfinite(e)]
        if e.size == 0:
            continue
        lo, hi = bootstrap_mean_ci(e, rng)
        rows.append(
            {**dict(zip(keys, key)), "n": int(e.size), "bias": e.mean(), "bias_lo": lo,
             "bias_hi": hi, "sd": e.std(ddof=1) if e.size > 1 else np.nan,
             "rms": float(np.sqrt(np.mean(e**2))),
             "within_20pct": float(np.mean(np.abs(np.expm1(e)) <= WITHIN))}
        )  # fmt: skip
    return pd.DataFrame(rows)


def _complete_matrix(g: pd.DataFrame, present: list[str]):
    """Errors as {batch: array (specimens, methods)} on the methods finite for every specimen."""
    wide = (
        g[g["Batch"].isin(present)]
        .set_index(["Batch", "Observation_id", "Method"])["Error_load"]
        .unstack("Method")
    )
    wide = wide.loc[:, np.isfinite(wide.to_numpy(float)).all(axis=0)]
    blocks = {b: wide.loc[b].to_numpy(float) for b in present if b in wide.index.levels[0]}
    return list(wide.columns), blocks


def _additive(means: np.ndarray):
    """Additive fit of means[..., methods, batches]: (mu, alpha, beta, residual)."""
    mu = means.mean(axis=(-2, -1), keepdims=True)
    alpha = means.mean(axis=-1, keepdims=True) - mu
    beta = means.mean(axis=-2, keepdims=True) - mu
    return mu, alpha, beta, means - mu - alpha - beta


def two_way(assessed: pd.DataFrame, batches: list[str]):
    """Additive decomposition of the cell biases per path, with joint bootstrap intervals.

    bias_mb = mu + alpha_m + beta_b + (alpha beta)_mb on the methods that are finite for every
    assessed failure of the batches compared, for every subset of >= 2 batches of a path. The
    same specimens are resampled for all methods (within each batch), which keeps the
    correlation between methods: every method sees the same tests, so the batch offsets are far
    better determined than independent cell intervals would suggest.

    Returns (summary, offsets, residuals):

    - summary, one row per path and subset: n_methods, mu, method_sd, interaction_rms,
      noise_rms (interaction expected from sampling alone, from the bootstrap), and
      interaction_excess = sqrt(max(interaction_rms^2 - noise_rms^2, 0));
    - offsets, long: Path, batches, Batch, n, offset with its 95 % interval;
    - residuals of the largest subset: Method, Batch, Path, residual.
    """
    rng = np.random.default_rng(SEED + 3)
    summary, offsets, residuals = [], [], []
    for path, g in assessed.groupby("Path", sort=False):
        sizes = g[np.isfinite(g["Error_load"])].groupby("Batch")["Observation_id"].nunique()
        available = [b for b in batches if sizes.get(b, 0) >= MIN_N]
        subsets = [c for k in range(len(available), 1, -1) for c in combinations(available, k)]
        for present in map(list, subsets):
            methods, blocks = _complete_matrix(g, present)
            if len(methods) < MIN_N or len(blocks) < len(present):
                continue
            means = np.stack([blocks[b].mean(axis=0) for b in present], axis=-1)
            mu, alpha, beta, res = _additive(means)
            boot = np.stack(
                [
                    blocks[b][rng.integers(0, len(blocks[b]), (N_BOOT, len(blocks[b])))].mean(1)
                    for b in present
                ],
                axis=-1,
            )
            _, _, beta_b, res_b = _additive(boot)
            rms2 = float(np.mean(res**2))
            noise2 = max(float(np.mean(res_b**2)) - rms2, 0.0)
            tag = "+".join(present)
            row = dict(Path=path, batches=tag, n_methods=len(methods), mu=float(mu.squeeze()),
                       method_sd=float(alpha.std(ddof=1)), interaction_rms=np.sqrt(rms2),
                       noise_rms=np.sqrt(noise2),
                       interaction_excess=np.sqrt(max(rms2 - noise2, 0.0)))  # fmt: skip
            for k, b in enumerate(present):
                lo, hi = np.percentile(beta_b[:, 0, k], [2.5, 97.5])
                row[f"offset_{b}"] = float(beta[0, k])
                offsets.append(
                    dict(Path=path, batches=tag, Batch=b, n=len(blocks[b]),
                         offset=float(beta[0, k]), offset_lo=lo, offset_hi=hi)
                )  # fmt: skip
            summary.append(row)
            if len(present) == len(available):
                long = pd.DataFrame(res, index=methods, columns=present).rename_axis("Method")
                long = long.reset_index().melt(
                    id_vars="Method", var_name="Batch", value_name="residual"
                )
                residuals.append(long.assign(Path=path))
    return (
        pd.DataFrame(summary),
        pd.DataFrame(offsets),
        pd.concat(residuals) if residuals else pd.DataFrame(),
    )


def predicted_offsets(offsets: pd.DataFrame, relative: pd.DataFrame) -> pd.DataFrame:
    """Add the method-free prediction of each batch offset (column ``predicted``).

    For a method that is homogeneous of degree one in the loads and calibrated on the axial
    curve of its own batch, the error of a test is e_L = ln c_m + ln[amp / S_A(N)], c_m a
    constant of the method on the path. Its batch offset is then the mean of ln[amp / S_A(N)]
    over the tests of the batch, centred over the batches compared: a property of the S-N data
    alone (``batch_reference.relative_strengths`` with design = target).
    """
    if offsets.empty:
        return offsets
    own = relative[relative["Design"] == relative["Target"]] if len(relative) else relative
    q = {(r.Target, r.Path): r.q for r in own.itertuples()}
    out = offsets.copy()
    out["predicted"] = np.nan
    for _, idx in out.groupby(["Path", "batches"]).groups.items():
        values = [q.get((out.at[i, "Batch"], out.at[i, "Path"])) for i in idx]
        if all(v is not None for v in values):
            out.loc[idx, "predicted"] = np.asarray(values) - np.mean(values)
    return out


def censored_mean(
    failures: np.ndarray, bounds: np.ndarray, *, tol: float = 1e-10, max_iter: int = 10_000
) -> float:
    """Mean of a normal sample with right-censored observations (true value > bound), by MLE.

    Expectation-maximisation for the normal distribution with right censoring: each censored
    value is replaced by its conditional moments given that it exceeds its bound,
    E[X | X > c] = mu + sd * lambda(z) and E[X^2 | X > c] = mu^2 + sd^2 + sd (mu + c) lambda(z),
    with z = (c - mu) / sd and lambda the inverse Mills ratio phi(z) / (1 - Phi(z)); the mean and
    variance are then re-estimated. The iteration converges to the maximum-likelihood estimate.
    Standard library and numpy only (scipy is not a project dependency).
    """
    from math import erfc, exp, pi, sqrt

    failures = np.asarray(failures, float)
    bounds = np.asarray(bounds, float)
    if bounds.size == 0 or failures.size < MIN_N:
        return float(failures.mean()) if failures.size else float("nan")

    def mills(z: float) -> float:
        tail = 0.5 * erfc(z / sqrt(2.0))
        if tail < 1e-300:  # far tail: lambda(z) -> z
            return z
        return exp(-0.5 * z * z) / sqrt(2.0 * pi) / tail

    n = failures.size + bounds.size
    mu = float(np.concatenate([failures, bounds]).mean())
    sd = max(float(failures.std(ddof=1)), 1e-6)
    for _ in range(max_iter):
        lam = np.array([mills((c - mu) / sd) for c in bounds])
        first = mu + sd * lam
        second = mu**2 + sd**2 + sd * (mu + bounds) * lam
        mu_new = (failures.sum() + first.sum()) / n
        var = ((failures**2).sum() + second.sum()) / n - mu_new**2
        sd_new = sqrt(max(var, 1e-24))
        converged = abs(mu_new - mu) < tol and abs(sd_new - sd) < tol
        mu, sd = mu_new, sd_new
        if converged:
            break
    return float(mu)


def runout_sensitivity(own: pd.DataFrame, batches: list[str]) -> pd.DataFrame:
    """Batch offsets with the run-outs entered as censored observations, against failures only.

    A run-out at N_ro bounds the error from below: the test would have failed later, so its
    e_L exceeds the value computed at N_ro. The bound is available where the damage parameter
    is homogeneous in the loads (e_L = ln(DP / DP_fit) / d, with d read from the failures of
    the same file); the other methods are left out of this check. Both offsets are computed on
    the same methods.
    """
    rows = []
    for (method, batch), g in own.groupby(["Method", "Batch"]):
        fail = g[(g["Flag"] == 1) & np.isfinite(g["Error_load"]) & (g["Error_load"] != 0)]
        with np.errstate(divide="ignore", invalid="ignore"):
            degree = (fail["Error_ln_dp"] / fail["Error_load"]).to_numpy(float)
        degree = degree[np.isfinite(degree)]
        if degree.size < MIN_N or np.ptp(degree) > 1e-6 * abs(np.median(degree)):
            continue  # not homogeneous: no closed-form bound
        d = float(np.median(degree))
        ro = g[g["Flag"] == 2]
        bound = np.log(ro["CP_value"].astype(float) / ro["DP_fit"].astype(float)) / d
        ro = ro.assign(bound=bound)[np.isfinite(bound)]
        for path, f in g[g["Flag"] == 1].groupby("Path"):
            e = f["Error_load"].to_numpy(float)
            e = e[np.isfinite(e)]
            c = ro[ro["Path"] == path]["bound"].to_numpy(float)
            if e.size < MIN_N:
                continue
            rows.append(dict(Method=method, Batch=batch, Path=path, n=e.size, n_runouts=c.size,
                             bias=e.mean(), bias_censored=censored_mean(e, c)))  # fmt: skip
    cells = pd.DataFrame(rows)
    out = []
    if cells.empty:
        return cells
    for path, g in cells.groupby("Path", sort=False):
        table = g.pivot(index="Method", columns="Batch", values=["bias", "bias_censored"])
        present = [b for b in batches if b in table["bias"].columns]
        table = table.dropna()
        if len(present) < 2 or len(table) < MIN_N:
            continue
        counts = g.groupby("Batch")["n_runouts"].max()
        for column in ("bias", "bias_censored"):
            means = table[column][present].to_numpy(float)
            beta = means.mean(axis=0) - means.mean()
            for b, value in zip(present, beta):
                out.append(dict(Path=path, Batch=b, n_methods=len(table), kind=column,
                                n_runouts=int(counts.get(b, 0)), offset=float(value)))  # fmt: skip
    wide = pd.DataFrame(out)
    if wide.empty:
        return wide
    wide = wide.pivot_table(
        index=["Path", "Batch", "n_methods", "n_runouts"], columns="kind", values="offset"
    ).reset_index()
    return wide.rename(columns={"bias": "offset", "bias_censored": "offset_censored"})


def transfer_penalty(own: pd.DataFrame, moved: pd.DataFrame) -> pd.DataFrame:
    """RMS of the own-batch and transfer errors on the same tests (own-batch Flag 1, paired)."""
    a = own[own["Flag"] == 1][["Method", "Batch", "Observation_id", "Error_load"]]
    b = moved[moved["Flag"] == 1][["Method", "Batch", "Design", "Observation_id", "Error_load"]]
    j = b.merge(a, on=["Method", "Batch", "Observation_id"], suffixes=("_transfer", "_own"))
    j = j[np.isfinite(j["Error_load_transfer"]) & np.isfinite(j["Error_load_own"])]
    rows = []
    for (method, design, target), g in j.groupby(["Method", "Design", "Batch"], sort=True):
        if len(g) < MIN_N:
            continue
        rms_own = float(np.sqrt(np.mean(g["Error_load_own"] ** 2)))
        rms_tr = float(np.sqrt(np.mean(g["Error_load_transfer"] ** 2)))
        rows.append(
            dict(Method=method, Design=design, Target=target, n=len(g),
                 bias_own=float(g["Error_load_own"].mean()),
                 bias_transfer=float(g["Error_load_transfer"].mean()), rms_own=rms_own,
                 rms_transfer=rms_tr, penalty=rms_tr - rms_own)
        )  # fmt: skip
    return pd.DataFrame(rows)


def transfer_prediction(own_cells, transfer_cells, relative: pd.DataFrame) -> pd.DataFrame:
    """Transfer bias against its method-free prediction, per design -> target pair and path.

    For a method calibrated on batch A, the bias on the tests of batch B minus its bias on the
    same path of A is, for a damage parameter homogeneous of degree one, the difference of the
    relative strengths q(B | A) - q(A | A) of ``batch_reference.relative_strengths``: how the
    tests of B and of A lie against the axial curve of A. The observed shift is the median over
    the methods of bias(A -> B, p) - bias(A, p); the bias of A on its own design path is zero by
    construction.
    """
    q = {(r.Design, r.Target, r.Path): r.q for r in relative.itertuples()}
    own = {(r.Method, r.Batch, r.Path): r.bias for r in own_cells.itertuples()}
    rows = []
    for (design, target, path), g in transfer_cells.groupby(["Design", "Target", "Path"]):
        g = g[g["n"] >= MIN_N]
        shifts = []
        for r in g.itertuples():
            base = own.get((r.Method, design, path))
            if base is None and path in ("A", "T"):
                base = 0.0  # design path of the method: not assessed on its own batch
            if base is not None:
                shifts.append(r.bias - base)
        if len(shifts) < MIN_N:
            continue
        a, b = q.get((design, design, path)), q.get((design, target, path))
        q1, med, q3 = np.percentile(shifts, [25, 50, 75])
        rows.append(
            dict(Design=design, Target=target, Path=path, n_methods=len(shifts), shift=med,
                 shift_q1=q1, shift_q3=q3,
                 predicted=b - a if a is not None and b is not None else np.nan,
                 bias_median=float(g["bias"].median()))
        )  # fmt: skip
    return pd.DataFrame(rows)


COVERAGE_Z = 1.645  # one-sided 95 % under a normal model of e_L


def designer_summary(own: pd.DataFrame, moved: pd.DataFrame) -> pd.DataFrame:
    """Load factor that covers 95 % of the tests, own batch against transfer, same tests.

    gamma = exp(-(mean - 1.645 SD)) of e_L: the factor by which the load allowed by the method
    must be divided so that 95 % of the tests are on the safe side, under a normal model of
    e_L. gamma < 1 means the method is conservative even at that level.
    """
    a = own[own["Flag"] == 1][["Method", "Batch", "Observation_id", "Error_load"]]
    b = moved[moved["Flag"] == 1][["Method", "Batch", "Design", "Observation_id", "Error_load"]]
    j = b.merge(a, on=["Method", "Batch", "Observation_id"], suffixes=("_transfer", "_own"))
    j = j[np.isfinite(j["Error_load_transfer"]) & np.isfinite(j["Error_load_own"])]
    rows = []
    for (method, design, target), g in j.groupby(["Method", "Design", "Batch"], sort=True):
        if len(g) < 2 * MIN_N:
            continue

        def gamma(e):
            return float(np.exp(-(e.mean() - COVERAGE_Z * e.std(ddof=1))))

        rows.append(dict(Method=method, Design=design, Target=target, n=len(g),
                         gamma_own=gamma(g["Error_load_own"]),
                         gamma_transfer=gamma(g["Error_load_transfer"])))  # fmt: skip
    return pd.DataFrame(rows)


# ---- analysis, figures, tables ---------------------------------------------------------------
def analyse(out: Path, family: str) -> Path:
    spec = families()[family]
    batches = spec["batches"]
    from multifat.cli.run import transfer_material_id

    base = out / family
    own = read_results(base / "own", {b: (b, None) for b in batches})
    moved = read_results(
        base / "transfer",
        {transfer_material_id(t, dsg): (t, dsg) for t, dsg in permutations(batches, 2)},
    )
    if own.empty:
        raise SystemExit(f"no own-batch results in {base / 'own'}; run first")
    paths = loading_paths(batches)
    own = own.merge(paths, on=["Batch", "Observation_id"], how="left", validate="m:1")
    if not moved.empty:
        moved = moved.merge(paths, on=["Batch", "Observation_id"], how="left", validate="m:1")
    dest = OUT_ROOT / family
    dest.mkdir(parents=True, exist_ok=True)

    from examples import batch_plots, batch_reference

    assessed = own[own["Flag"] == 1]
    write = dict(index=False, float_format="%.5g")
    reference = batch_reference.write_reference(
        dest, batches, fabest=bool(spec.get("fabest_benchmark"))
    )
    cells = cell_statistics(assessed, ["Method", "Batch", "Path"])
    decomposition, offsets, residuals = two_way(assessed, batches)
    offsets = predicted_offsets(offsets, reference["relative"])
    if len(residuals):
        cells = cells.merge(residuals, on=["Method", "Batch", "Path"], how="left")
    cells.to_csv(dest / "own_cells.csv", **write)
    decomposition.to_csv(dest / "own_two_way.csv", **write)
    offsets.to_csv(dest / "own_offsets.csv", **write)
    runout_sensitivity(own, batches).to_csv(dest / "own_runout_sensitivity.csv", **write)
    penalty = prediction = designer = pd.DataFrame()
    if not moved.empty:
        tcells = cell_statistics(moved[moved["Flag"] == 1], ["Method", "Design", "Batch", "Path"])
        tcells = tcells.rename(columns={"Batch": "Target"})
        tcells.to_csv(dest / "transfer_cells.csv", **write)
        penalty = transfer_penalty(own, moved)
        penalty.to_csv(dest / "transfer_penalty.csv", **write)
        prediction = transfer_prediction(cells, tcells, reference["relative"])
        prediction.to_csv(dest / "transfer_prediction.csv", **write)
        designer = designer_summary(own, moved)
        designer.to_csv(dest / "transfer_designer.csv", **write)
    for stale in ("own_batch_effect.csv", "batch_reference.csv"):  # replaced outputs
        (dest / stale).unlink(missing_ok=True)
    # a figure that cannot be drawn any more (e.g. no own-batch result) must not survive
    for old in DOCS_IMG_DIR.glob(f"batch_comparison_{family}_*.png"):
        old.unlink()
    _latex_tables(dest, decomposition, penalty, batches)
    _figure_bias(cells, batches, family)
    if len(penalty):
        _figure_transfer(penalty, batches, family)
    batch_plots.all_figures(
        family, batches, reference=reference, offsets=offsets, cells=cells,
        penalty=penalty, designer=designer,
    )  # fmt: skip
    manifest = dict(
        family=family, batches=batches, notes=spec.get("notes", ""),
        constants=spec.get("constants", {}),
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        git_commit=_git("rev-parse", "HEAD"), registered_models=registered_models(),
        implementation_revision=_revision(), dataset_fingerprints=dataset_fingerprints(batches),
        methods_with_results=sorted(own["Method"].unique()),
        statistic="Error_load = ln(1/s*), assessed failures (Flag 1)", bootstrap=[N_BOOT, SEED],
    )  # fmt: skip
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"results -> {dest}")
    return dest


def _git(*flags) -> str | None:
    proc = subprocess.run(["git", *flags], cwd=REPO_ROOT, capture_output=True, text=True)
    return proc.stdout.strip() or None


def _revision() -> str:
    from multifat.models.references import IMPLEMENTATION_REVISION

    return IMPLEMENTATION_REVISION


def _short(batch: str, family: str) -> str:
    return batch.removeprefix(family + "_") or batch


def _latex_tables(dest: Path, decomposition, penalty, batches) -> None:
    family = dest.name
    lines = [
        r"\begin{tabular}{llrr" + "r" * len(batches) + "rr}",
        r"\hline",
        "Path & Batches & Methods & "
        + " & ".join(rf"$\beta_{{\mathrm{{{_short(b, family)}}}}}$" for b in batches)
        + r" & Interaction RMS & Sampling noise & Method SD \\",
        r"\hline",
    ]
    for _, r in decomposition.iterrows():
        offsets = " & ".join(
            f"{r[f'offset_{b}']:+.3f}" if pd.notna(r.get(f"offset_{b}")) else "--" for b in batches
        )
        lines.append(
            f"{r['Path']} & {r['batches'].count('+') + 1} & {r['n_methods']} & {offsets} & "
            f"{r['interaction_rms']:.3f} & {r['noise_rms']:.3f} & {r['method_sd']:.3f} \\\\"
        )
    lines += [r"\hline", r"\end{tabular}"]
    (dest / "own_two_way.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if penalty.empty:
        return
    s = (
        penalty.groupby("Method")
        .apply(
            lambda g: pd.Series(
                dict(pairs=len(g), rms_own=np.sqrt(np.average(g.rms_own**2, weights=g.n)),
                     rms_transfer=np.sqrt(np.average(g.rms_transfer**2, weights=g.n)))
            ),
            include_groups=False,
        )
        .sort_values("rms_transfer")
    )  # fmt: skip
    s["penalty"] = s["rms_transfer"] - s["rms_own"]
    s.to_csv(dest / "transfer_summary.csv", float_format="%.5g")
    rows = [r"\begin{tabular}{lrrrr}", r"\hline",
            r"Method & Pairs & RMS own & RMS transfer & $\Delta$ \\", r"\hline"]  # fmt: skip
    for m, r in s.iterrows():
        rows.append(
            f"{m} & {int(r.pairs)} & {r.rms_own:.3f} & {r.rms_transfer:.3f} & {r.penalty:+.3f} \\\\"
        )
    rows += [r"\hline", r"\end{tabular}"]
    (dest / "transfer_summary.tex").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _figure_bias(cells: pd.DataFrame, batches, family: str) -> None:
    """Own-batch bias with 95 % interval per method, one panel per path shared by >= 2 batches."""
    import matplotlib.pyplot as plt

    apply_academic_style()
    shared = [
        p
        for p in PATH_ORDER
        if cells[(cells["Path"] == p) & (cells["n"] >= MIN_N)]["Batch"].nunique() >= 2
    ]
    if not shared:
        return
    methods = sorted(cells["Method"].unique())
    fig, axes = plt.subplots(
        1, len(shared), figsize=(2.6 * len(shared) + 1.2, 0.2 * len(methods) + 1.4), sharey=True
    )
    axes = np.atleast_1d(axes)
    y = {m: i for i, m in enumerate(methods)}
    handles = {}
    for ax, path in zip(axes, shared):
        sub = cells[(cells["Path"] == path) & (cells["n"] >= MIN_N)]
        counts = []
        for k, b in enumerate(batches):
            g = sub[sub["Batch"] == b]
            if g.empty:
                continue
            counts.append(f"{_short(b, family)} {int(g['n'].max())}")
            yy = np.array([y[m] for m in g["Method"]]) + (k - (len(batches) - 1) / 2) * 0.22
            handles[b] = ax.errorbar(
                g["bias"], yy,
                xerr=[g["bias"] - g["bias_lo"], g["bias_hi"] - g["bias"]],
                fmt=BATCH_MARKERS[k], ms=3.2, lw=0.8, color=BATCH_COLOURS[k],
            )  # fmt: skip
        ax.axvline(0.0, color="#52514e", lw=0.8)
        for x in (np.log(1 / 1.2), np.log(1 / 0.8)):  # +-20 % of the load
            ax.axvline(x, color="#888888", lw=0.5, ls=":")
        ax.set_title(f"{path}\n" + r"$n$: " + ", ".join(counts), fontsize=9)
        ax.set_xlabel(r"bias $\bar e_L = \overline{\ln(1/s^*)}$")
        ax.grid(axis="x", color="#dddddd", lw=0.5)
    fig.legend(
        [handles[b] for b in batches if b in handles],
        [_short(b, family) for b in batches if b in handles],
        loc="upper center", ncol=len(handles), bbox_to_anchor=(0.5, 0.965), frameon=False,
    )  # fmt: skip
    axes[0].set_yticks(range(len(methods)), methods, fontsize=7)
    axes[0].invert_yaxis()
    fig.suptitle(
        f"{family}: own-batch load error by loading path (mean, 95 % bootstrap interval; "
        "dotted: $\\pm$20 % of the load)",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(DOCS_IMG_DIR / f"batch_comparison_{family}_bias.png", dpi=200)
    plt.close(fig)


def _figure_transfer(penalty: pd.DataFrame, batches, family: str) -> None:
    """RMS of the load error, own batch and predicted from each other batch, on the same tests."""
    import matplotlib.pyplot as plt

    apply_academic_style()
    cols, labels = [], []
    for t in batches:
        cols.append((t, None))
        labels.append(f"{_short(t, family)}\nown")
        for dsg in batches:
            if dsg != t:
                cols.append((t, dsg))
                labels.append(f"{_short(t, family)}\nfrom {_short(dsg, family)}")
    methods = sorted(penalty["Method"].unique())
    grid = np.full((len(methods), len(cols)), np.nan)
    for i, m in enumerate(methods):
        for j, (t, dsg) in enumerate(cols):
            g = penalty[(penalty["Method"] == m) & (penalty["Target"] == t)]
            if dsg is None:
                g = g.drop_duplicates("Target")
                value = g["rms_own"]
            else:
                value = g[g["Design"] == dsg]["rms_transfer"]
            if len(value):
                grid[i, j] = float(value.iloc[0])
    fig, ax = plt.subplots(figsize=(1.0 * len(cols) + 2.0, 0.2 * len(methods) + 1.6))
    vmax = float(np.nanpercentile(grid, 98)) if np.isfinite(grid).any() else 1.0
    image = ax.imshow(grid, cmap="Blues", vmin=0.0, vmax=vmax, aspect="auto")
    for i in range(len(methods)):
        for j in range(len(cols)):
            if np.isfinite(grid[i, j]):
                ax.text(
                    j, i, f"{grid[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if grid[i, j] > 0.6 * vmax else "#0b0b0b",
                )  # fmt: skip
    ax.set_xticks(range(len(cols)), labels, fontsize=6.5)
    ax.set_yticks(range(len(methods)), methods, fontsize=7)
    for j in range(1, len(cols)):
        if cols[j][1] is None:
            ax.axvline(j - 0.5, color="white", lw=2)
    fig.colorbar(image, ax=ax, shrink=0.6, label=r"RMS of $e_L$ (same tests)")
    ax.set_title(f"{family}: own-batch versus cross-batch prediction")
    fig.tight_layout()
    fig.savefig(DOCS_IMG_DIR / f"batch_comparison_{family}_transfer.png", dpi=200)
    plt.close(fig)


# ---- staleness check -----------------------------------------------------------------------
def check(family: str | None = None) -> list[str]:
    """Reasons the committed results are out of date (empty list: current)."""
    problems = []
    for name, spec in families().items():
        if family and name != family:
            continue
        path = OUT_ROOT / name / "manifest.json"
        if not path.exists():
            problems.append(f"{name}: no results ({path.relative_to(REPO_ROOT)})")
            continue
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if manifest["batches"] != spec["batches"]:
            problems.append(f"{name}: batches changed")
        if manifest.get("constants", {}) != spec.get("constants", {}):
            problems.append(f"{name}: supplied constants changed")
        added = sorted(set(registered_models()) - set(manifest["registered_models"]))
        removed = sorted(set(manifest["registered_models"]) - set(registered_models()))
        if added or removed:
            problems.append(f"{name}: models added {added}, removed {removed}")
        if manifest["implementation_revision"] != _revision():
            problems.append(f"{name}: implementation revision changed")
        now = dataset_fingerprints(spec["batches"])
        changed = [b for b in spec["batches"] if manifest["dataset_fingerprints"].get(b) != now[b]]
        if changed:
            problems.append(f"{name}: datasets changed {changed}")
    return problems


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=["run", "analyse", "all", "check"])
    parser.add_argument("out", type=Path, nargs="?", help="folder of the pipeline runs")
    parser.add_argument("--family", default=None, help="one family (default: all)")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--models", default="all")
    args = parser.parse_args(argv)
    if args.action == "check":
        problems = check(args.family)
        for p in problems:
            print(p)
        if problems:
            raise SystemExit("rerun: uv run python examples/batch_comparison.py all <tmp_dir>")
        print("batch comparison results are current")
        return
    if args.out is None:
        parser.error(f"{args.action} needs the run folder")
    for family in [args.family] if args.family else list(families()):
        if args.action in ("run", "all"):
            run(args.out, family, args.jobs, args.models)
        if args.action in ("analyse", "all"):
            analyse(args.out, family)


if __name__ == "__main__":
    main()
