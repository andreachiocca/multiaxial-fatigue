"""Explicit end-to-end check (not part of the quick unit-test suite).

Run from repository root: python tests/run_material_validation.py /tmp/fatigue-audit
This exercises all retained base models plus FS_ext with the standard plane grid.
Outputs are new files in the requested directory; archived Results are untouched.
"""
from pathlib import Path
import contextlib
import runpy
import sys
import os

os.environ.setdefault('MPLBACKEND','Agg')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from fatigue.models.registry import _available_cp_methods, _available_direct_models
from fatigue.eval.metrics import log10_dp_ratio


def main():
    out = Path(sys.argv[1] if len(sys.argv)>1 else '/tmp/fatigue-audit').resolve()
    out.mkdir(parents=True,exist_ok=True)
    runner=runpy.run_path('1_run_material.py')
    g=runner['main'].__globals__
    g['RESULTS_DIR']=str(out)
    g['MAKE_PLOTS']=False
    counts={'AISI316L':31,'42CrMo4_QT':70,'Al7075_T6':87,'Ti6Al4V_andrea_(AXIAL)':58}
    names=list(_available_cp_methods())+list(_available_direct_models())+['FS_ext']
    rows=0
    for name in names:
        for material,count in counts.items():
            g['MODEL_NAME'],g['MATERIAL_NAME']=name,material
            with (out/'validation.log').open('a') as log, contextlib.redirect_stdout(log):
                runner['main']()
            df=pd.read_csv(out/f'{name}_{material}.csv')
            assert len(df)==count,(name,material,len(df),count)
            assert np.isfinite(df.CP_value).all(),(name,material,'nonfinite DP')
            assert np.isfinite(df.Error_log10_dp).all(),(name,material,'missing log ratio')
            np.testing.assert_allclose(df.Error_log10_dp,log10_dp_ratio(df.DP_fit,df.CP_value),atol=1e-12)
            assert df.Method_reference.notna().all(),(name,material,'missing reference')
            if df.Metric_mode.eq('DP_DIFF').all():
                assert df.Nf_expected.isna().all(),(name,material,'flat curve life prediction')
            rows+=len(df)
            print(f'PASS {name} / {material}: {len(df)} rows',flush=True)
    print(f'PASS {len(names)*len(counts)} runs, {rows} rows',flush=True)


if __name__=='__main__': main()
