"""Damage comparison independent of master-curve slope."""
import numpy as np


def log10_dp_ratio(expected, observed):
    """log10(DP_e/DP); nonpositive or nonfinite inputs yield NaN."""
    expected, observed = np.broadcast_arrays(np.asarray(expected, dtype=float),
                                             np.asarray(observed, dtype=float))
    out = np.full(expected.shape, np.nan)
    ok = np.isfinite(expected) & np.isfinite(observed) & (expected > 0) & (observed > 0)
    out[ok] = np.log10(expected[ok])-np.log10(observed[ok])
    return out
