"""Explicit end-to-end check (not part of the quick unit-test suite).

Run from repository root: python tests/run_material_validation.py /tmp/fatigue-audit
This exercises all registered base models plus FS_ext with the standard plane grid, on
workbook constants only: a model whose required constants the workbook lacks (MDC, MMCCM need the
Basquin and strain-life constants that multifat-pipeline supplies) is reported as SKIP.
Outputs are new files in the requested directory; archived Results are untouched.
"""

from pathlib import Path
import contextlib
import sys
import os

os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import pandas as pd
from multifat.models.registry import model_names
from multifat.eval.metrics import log10_dp_ratio
from multifat.cli import run as runner
from multifat.imp.experimental_data import load_material_workbook
from multifat.models.errors import MissingParamsError


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/fatigue-audit").resolve()
    out.mkdir(parents=True, exist_ok=True)
    runner.RESULTS_DIR = str(out)
    runner.MAKE_PLOTS = False
    # Plain failures of the materials with life data (the fatigue-limit sets need the pipeline)
    counts = {
        "AISI316L": 39,
        "42CrMo4_Santus": 45,
        "42CrMo4_Reis": 25,
        "Al7075_T6": 87,
        "Ti6Al4V_Meneghetti": 10,
        "Ti6Al4V_Berto": 48,
    }
    runout_counts = {
        m: len(load_material_workbook(runner.resolve_workbook(m))[m]["runouts"]) for m in counts
    }
    names = model_names() + ["FS_ext"]
    rows = runs_ok = 0
    for name in names:
        for material, count in counts.items():
            runner.MODEL_NAME, runner.MATERIAL_NAME = name, material
            try:
                with (out / "validation.log").open("a") as log, contextlib.redirect_stdout(log):
                    runner.main()
            except MissingParamsError as e:
                print(f"SKIP {name} / {material}: {e}", flush=True)
                continue
            df = pd.read_csv(out / f"{name}_{material}.csv")
            runs = df[df.Flag == 2]
            assert len(runs) == runout_counts[material], (name, material, "run-out rows")
            assert runs.Error_log10_dp.isna().all() and runs.Error_ln.isna().all()
            assert np.isfinite(runs.CP_value).all(), (name, material, "nonfinite run-out DP")
            df = df[df.Flag != 2]
            assert len(df) == count, (name, material, len(df), count)
            assert np.isfinite(df.CP_value).all(), (name, material, "nonfinite DP")
            assert np.isfinite(df.Error_log10_dp).all(), (name, material, "missing log ratio")
            np.testing.assert_allclose(
                df.Error_log10_dp, log10_dp_ratio(df.DP_fit, df.CP_value), atol=1e-12
            )
            assert df.Method_reference.notna().all(), (name, material, "missing reference")
            if df.Metric_mode.eq("DP_DIFF").all():
                assert df.Nf_expected.isna().all(), (name, material, "flat curve life prediction")
            rows += len(df)
            runs_ok += 1
            print(f"PASS {name} / {material}: {len(df)} rows", flush=True)
    print(f"PASS {runs_ok} of {len(names) * len(counts)} runs, {rows} rows", flush=True)


if __name__ == "__main__":
    main()
