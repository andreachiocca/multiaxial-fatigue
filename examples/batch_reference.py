"""Method-free comparison of the batches of one material (issue #71).

Before any fatigue criterion is applied, the batches of a family are compared through their own
S-N data, per loading path and per test series:

- ``series_strengths``: Basquin fit of the R = -1 failures of each batch, path and series, and
  the strength S(N_REF) with a 95 % bootstrap interval. A series is one S-N curve of the source
  (the workbook's ``Curve ID`` when present, otherwise authors and year). The pooled fit of a
  batch and path is the series ``all``.
- ``strength_ratios``: strength of each path relative to the axial strength of the same batch,
  S_path(N_REF) / S_A(N_REF): tau_a/sigma_a for torsion, and the axial amplitude under combined
  loading over the pure axial one.
- ``relative_strengths``: the tests of each batch and path measured against the axial Basquin
  curve of a design batch, at their own lives: the method-free prediction of the batch offsets
  and of the transfer shifts.
- ``variance_components``: the spread between the datasets of a family split into a material
  batch part and a test-series part, the latter from the FABEST series (one material batch,
  maintainer, d105).
- ``fabest_benchmark``: every unnotched, R = -1 axial and torsion series of the FABEST database
  (one steel, five laboratories), read from the committed source file. The scatter of S(N_REF)
  between these series is the yardstick for a difference between two batches: a difference that
  small can come from the specimen, the machining or the laboratory alone.

Repository choices:
  - Basquin regression of log amplitude on log life, as in ``multifat.eval.material_fits``;
  - N_REF = 1e5 cycles, inside the tested range of every 42CrMo4 batch; a strength outside the
    tested lives of its series is flagged ``extrapolated``;
  - intervals: percentile bootstrap of the specimens of one series (pairs), fixed seed. The
    specimens of one series share the lot, the machine and the laboratory, so an interval
    describes that series, not the population of series.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

N_REF = 1e5  # reference life of the strength comparison, cycles
N_BOOT = 4000
SEED = 20261110
MIN_N = 3
PATH_NAMES = {
    "Tensile": "A",
    "Torsion": "T",
    "Multiaxial proportional": "AT-IP",
    "Multiaxial non-proportional": "AT-NP",
}
VON_MISES = 1.0 / np.sqrt(3.0)  # tau/sigma of the von Mises criterion


def basquin_strength(amp, N, N_ref: float = N_REF) -> tuple[float, float]:
    """(S(N_ref), exponent b) of amp = A N^b fitted on log amplitude against log life."""
    b, a = np.polyfit(np.log10(np.asarray(N, float)), np.log10(np.asarray(amp, float)), 1)
    return float(10.0 ** (a + b * np.log10(N_ref))), float(b)


def bootstrap_strength(amp, N, rng, N_ref: float = N_REF) -> np.ndarray:
    """Bootstrap sample of ln S(N_ref), resampling the specimens of one series."""
    amp, N = np.log10(np.asarray(amp, float)), np.log10(np.asarray(N, float))
    idx = rng.integers(0, amp.size, size=(N_BOOT, amp.size))
    x, y = N[idx], amp[idx]
    xm, ym = x.mean(axis=1, keepdims=True), y.mean(axis=1, keepdims=True)
    sxx = ((x - xm) ** 2).sum(axis=1)
    slope = np.where(sxx > 0, ((x - xm) * (y - ym)).sum(axis=1) / np.where(sxx > 0, sxx, 1), 0.0)
    return np.log(10.0) * (ym[:, 0] + slope * (np.log10(N_ref) - xm[:, 0]))


def reversed_tests(batch: str) -> pd.DataFrame:
    """R = -1 failures of a batch: path, series, reference amplitude (MPa), life, tau/sigma."""
    from multifat.cli import run as runmod

    data, mat = runmod.load_material(batch)
    block = data[mat]
    sheet = pd.read_excel(runmod.workbook_path(batch), "Dataset")
    if "Curve ID" in sheet:
        label = sheet["Curve ID"].astype(str)
    else:
        label = sheet["Authors"].astype(str) + " " + sheet["Article year"].astype(str)
    series = dict(zip(sheet["Specimen name"].astype(str), label))
    rows = []
    for case in block["cases"].values():
        H = case["harmonics"]
        if np.max(np.abs(np.asarray(H["S_mean"], float))) > 1e-6:
            continue
        amp = np.hypot(H["S_sin"], H["S_cos"])
        path = PATH_NAMES.get(runmod.classify_loading_type(H=H, meta=case["meta"]), "?")
        sigma, tau = float(amp[1, 1]), float(amp[0, 1])
        if path == "?" or float(amp[0, 0]) > 1e-6:
            continue
        rows.append(
            dict(Batch=batch, Path=path, Series=series.get(str(case["meta"]["specimen_id"]), "?"),
                 amp=tau if path == "T" else sigma, N=float(case["meta"]["Nf_cycles"]),
                 tau_over_sigma=tau / sigma if sigma > 0 else np.nan)
        )  # fmt: skip
    return pd.DataFrame(rows)


def _strength_row(g: pd.DataFrame, rng) -> dict | None:
    if len(g) < MIN_N or g["N"].nunique() < 2:
        return None
    S, b = basquin_strength(g["amp"], g["N"])
    lo, hi = np.exp(np.percentile(bootstrap_strength(g["amp"], g["N"], rng), [2.5, 97.5]))
    return dict(n=len(g), N_min=g["N"].min(), N_max=g["N"].max(), b=b, S_ref=S, S_lo=lo, S_hi=hi,
                decades=float(np.log10(g["N"].max() / g["N"].min())),
                extrapolated=not (g["N"].min() <= N_REF <= g["N"].max()))  # fmt: skip


def series_strengths(batches, Su: dict | None = None) -> pd.DataFrame:
    """S(N_REF) of every batch, path and series (and of the pooled series ``all``)."""
    rng = np.random.default_rng(SEED)
    rows = []
    for batch in batches:
        tests = reversed_tests(batch)
        for path, g in tests.groupby("Path", sort=False):
            groups = [("all", g)] + (list(g.groupby("Series")) if g["Series"].nunique() > 1 else [])
            for name, sub in groups:
                row = _strength_row(sub, rng)
                if row is None:
                    continue
                ratio = sub["tau_over_sigma"].dropna()
                row.update(Batch=batch, Path=path, Series=name,
                           tau_over_sigma=float(ratio.median()) if len(ratio) else np.nan)  # fmt: skip
                if Su and Su.get(batch):
                    row["S_ref_over_Su"] = row["S_ref"] / Su[batch]
                rows.append(row)
    lead = ["Batch", "Path", "Series"]
    out = pd.DataFrame(rows)
    return out[lead + [c for c in out.columns if c not in lead]] if len(out) else out


def strength_ratios(batches) -> pd.DataFrame:
    """S_path(N_REF) / S_A(N_REF) per batch and path, with a 95 % bootstrap interval."""
    rng = np.random.default_rng(SEED + 1)
    rows = []
    for batch in batches:
        tests = reversed_tests(batch)
        axial = tests[tests["Path"] == "A"]
        if len(axial) < MIN_N:
            continue
        ln_axial = bootstrap_strength(axial["amp"], axial["N"], rng)
        S_axial, _ = basquin_strength(axial["amp"], axial["N"])
        for path, g in tests[tests["Path"] != "A"].groupby("Path", sort=False):
            if len(g) < MIN_N or g["N"].nunique() < 2:
                continue
            S, _ = basquin_strength(g["amp"], g["N"])
            sample = bootstrap_strength(g["amp"], g["N"], rng) - ln_axial
            lo, hi = np.exp(np.percentile(sample, [2.5, 97.5]))
            rows.append(
                dict(Batch=batch, Path=path, n=len(g), n_axial=len(axial), ratio=S / S_axial,
                     ratio_lo=lo, ratio_hi=hi, ln_ratio=float(np.log(S / S_axial)),
                     ln_ratio_se=float(sample.std(ddof=1)),
                     extrapolated=not (
                         max(g["N"].min(), axial["N"].min())
                         <= N_REF
                         <= min(g["N"].max(), axial["N"].max())
                     ))
            )  # fmt: skip
    return pd.DataFrame(rows)


def relative_strengths(batches) -> pd.DataFrame:
    """Mean ln[amp_i / S_A(N_i)] of the tests of each batch and path, per design batch.

    S_A is the axial Basquin curve of the design batch, evaluated at the life of each test. For
    a damage parameter homogeneous of degree one, calibrated on that axial curve, the load error
    of a test is ln c_m plus this quantity, c_m a constant of the method on the path. Its mean
    over the tests of a batch therefore predicts, without any method, the batch offsets of the
    own-batch errors (design = target) and the shift of a transfer (design != target). Unlike
    the ratio at N_REF, it accounts for the lives at which each batch was actually tested.
    """
    tests = {b: reversed_tests(b) for b in batches}
    rows = []
    for design in batches:
        axial = tests[design][tests[design]["Path"] == "A"]
        if len(axial) < MIN_N or axial["N"].nunique() < 2:
            continue
        S_ref, b = basquin_strength(axial["amp"], axial["N"])
        for target in batches:
            for path, g in tests[target].groupby("Path", sort=False):
                q = np.log(g["amp"] / (S_ref * (g["N"] / N_REF) ** b))
                rows.append(
                    dict(Design=design, Target=target, Path=path, n=len(g), q=float(q.mean()),
                         q_se=float(q.std(ddof=1) / np.sqrt(len(g))) if len(g) > 1 else np.nan)
                )  # fmt: skip
    return pd.DataFrame(rows)


def batch_Su(batches) -> dict:
    """Tensile strength of each batch from its Summary (MPa), None when absent."""
    from multifat.cli import run as runmod

    out = {}
    for b in batches:
        data, mat = runmod.load_material(b)
        out[b] = data[mat]["material"].get("Tensile_strength_MPa")
    return out


# ---- FABEST reproducibility benchmark ----------------------------------------------------------
def fabest_benchmark() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Unnotched R = -1 axial and torsion series of the FABEST database, and their scatter.

    Returns the per-series table and a summary of the scatter of ln S(N_REF) between series:
    all unnotched series of one load case, and the series of one drawing (v07, hollow 8.1/7 mm
    hourglass). Series that cannot be read or fitted are listed with the reason.
    """
    import openpyxl

    from tools.convert_fabest import SHA256, SOURCE, curve_points, sha256

    if sha256(SOURCE) != SHA256:
        raise ValueError(f"{SOURCE}: SHA-256 differs from the recorded FABEST file")
    wb = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    index = [list(r) for r in wb["Master Index"].iter_rows(values_only=True)]
    header = next(r for r in index if "Curve ID" in r)
    col = {name: header.index(name) for name in header if name}
    rng = np.random.default_rng(SEED + 2)
    rows = []
    for r in index[index.index(header) + 1 :]:
        # rows of the index are shifted by one cell for most series: locate the curve id
        curve = next((c for c in r if isinstance(c, str) and c.startswith("FAB")), None)
        if curve is None or curve not in wb.sheetnames:
            continue
        shift = r.index(curve) - col["Curve ID"]
        cell = lambda name: r[col[name] + shift]  # noqa: E731
        kind, load, R = str(cell("Type")), cell("Load case"), cell("R")
        if not kind.startswith("Unnotched") or load not in ("Axial", "Torsion") or R != -1:
            continue
        if int(curve[3:]) >= 500:  # FAB5xx: strain-controlled low-cycle tests
            continue
        base = dict(Curve=curve, Load=load, Lab=cell("Tested at"), Type=kind,
                    Drawing=cell("Drawing"), D_mm=cell("Minimum outer diameter"),
                    d_mm=cell("Inner diameter"),
                    Finish=cell("Final surface operation on outer surface"),
                    Machine=cell("Test machine"))  # fmt: skip
        try:
            points = [p for p in curve_points(wb[curve]) if not p["runout"]]
        except (ValueError, KeyError, TypeError, IndexError) as exc:
            rows.append({**base, "skipped": f"{type(exc).__name__}: {exc}"})
            continue
        frame = pd.DataFrame(
            dict(amp=[p["sa"] or p["ta"] for p in points], N=[p["N"] for p in points])
        )
        row = _strength_row(frame, rng)
        if row is None:
            rows.append({**base, "skipped": f"{len(frame)} failures"})
            continue
        rows.append({**base, **row, "freq_Hz": points[0]["freq"], "skipped": ""})
    table = pd.DataFrame(rows)
    ok = table[(table["skipped"] == "") & ~table["extrapolated"].astype(bool)]
    summary = []
    for load, g in ok.groupby("Load"):
        groups = [("all unnotched series", g)]
        same = g[g["Drawing"].astype(str) == "v07"]
        if len(same) >= 3:
            groups.append(("drawing v07 (hollow 8.1/7 mm hourglass)", same))
        for name, sub in groups:
            ln = np.log(sub["S_ref"].to_numpy(float))
            summary.append(
                dict(Load=load, group=name, series=len(sub), labs=sub["Lab"].nunique(),
                     S_median=float(np.exp(np.median(ln))), S_min=float(np.exp(ln.min())),
                     S_max=float(np.exp(ln.max())), sd_ln=float(ln.std(ddof=1)),
                     range_ln=float(ln.max() - ln.min()))
            )  # fmt: skip
    return table, pd.DataFrame(summary)


def _chi2_cdf(x: float, df: int) -> float:
    """Chi-square distribution function: regularized lower incomplete gamma P(df/2, x/2)."""
    from math import exp, lgamma, log

    a, x = 0.5 * df, 0.5 * x
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:  # series
        term = total = 1.0 / a
        for n in range(1, 10_000):
            term *= x / (a + n)
            total += term
            if term < total * 1e-16:
                break
        return total * exp(-x + a * log(x) - lgamma(a))
    # continued fraction for the upper tail (Lentz)
    tiny = 1e-300
    b, c, d = x + 1.0 - a, 1.0 / tiny, 1.0 / (x + 1.0 - a)
    h = d
    for i in range(1, 10_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return 1.0 - exp(-x + a * log(x) - lgamma(a)) * h


def chi2_ppf(p: float, df: int) -> float:
    """Quantile of the chi-square distribution, by bisection on ``_chi2_cdf`` (no scipy)."""
    lo, hi = 0.0, max(1.0, float(df))
    while _chi2_cdf(hi, df) < p:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if _chi2_cdf(mid, df) < p else (lo, mid)
        if hi - lo < 1e-14 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


def variance_components(series: pd.DataFrame, fabest_summary: pd.DataFrame) -> pd.DataFrame:
    """Split the spread of ln S(N_REF) between the datasets of a family into batch and series.

    A dataset (workbook) confounds the material batch with the test series (specimen,
    machining, laboratory, machine). The FABEST database is one material batch tested in
    several series and laboratories, so the variance of ln S between its unnotched series is the
    series component alone, s_series^2. The batch component is then estimated as
    s_batch^2 = max(s_datasets^2 - s_series^2, 0), with s_datasets^2 the variance of ln S between
    the pooled datasets of the family (raw, and divided by Su). With two or three datasets the
    estimate has one or two degrees of freedom: ``datasets_sd_hi`` is the upper end of the 95 %
    chi-square interval of s_datasets, which bounds what the data can exclude.
    """
    if fabest_summary is None or fabest_summary.empty:
        return pd.DataFrame()
    pooled = series[series["Series"] == "all"]
    rows = []
    for path, load in (("A", "Axial"), ("T", "Torsion")):
        ref = fabest_summary[
            (fabest_summary["Load"] == load) & (fabest_summary["group"] == "all unnotched series")
        ]
        g = pooled[pooled["Path"] == path]
        if ref.empty or len(g) < 2:
            continue
        s_series = float(ref["sd_ln"].iloc[0])
        for scale in ("raw", "over_Su"):
            column = "S_ref" if scale == "raw" else "S_ref_over_Su"
            if column not in g or g[column].isna().any():
                continue
            ln = np.log(g[column].to_numpy(float))
            k = ln.size
            var = float(ln.var(ddof=1))
            hi = np.sqrt(var * (k - 1) / chi2_ppf(0.025, k - 1))
            rows.append(
                dict(Path=path, scale=scale, datasets=k, datasets_sd=np.sqrt(var),
                     datasets_sd_hi=float(hi), series_sd=s_series,
                     series_count=int(ref["series"].iloc[0]),
                     batch_sd=float(np.sqrt(max(var - s_series**2, 0.0))))
            )  # fmt: skip
    return pd.DataFrame(rows)


def write_reference(dest: Path, batches, *, fabest: bool = False) -> dict:
    """Write the method-free tables of a family into `dest`; returns them."""
    write = dict(index=False, float_format="%.5g")
    out = {
        "series": series_strengths(batches, batch_Su(batches)),
        "ratios": strength_ratios(batches),
        "relative": relative_strengths(batches),
    }
    out["series"].to_csv(dest / "reference_series.csv", **write)
    out["relative"].to_csv(dest / "reference_relative.csv", **write)
    out["ratios"].to_csv(dest / "reference_ratios.csv", **write)
    if fabest:
        out["fabest"], out["fabest_summary"] = fabest_benchmark()
        out["fabest"].to_csv(dest / "fabest_series.csv", **write)
        out["fabest_summary"].to_csv(dest / "fabest_summary.csv", **write)
        out["components"] = variance_components(out["series"], out["fabest_summary"])
        out["components"].to_csv(dest / "reference_components.csv", **write)
    return out
