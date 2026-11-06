"""Worked example: a generated multiaxial load history assessed with the Findley criterion.

``python -m multifat.history.example --out DIR`` builds a proportional tension-torsion history of
variable amplitude from the documented function ``tension_torsion_block``, records the sampled
history, counts it by rainflow on the dominant channel, evaluates the Findley damage parameter of each
event and accumulates the damage with Miner's rule on an illustrative Basquin curve. The curve and the
material constant are placeholders chosen for the example: they are not a calibration, and the output
is not a validation of the criterion.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

from .accumulate import SNCurve
from .assess import assess_history
from .signal import TimeHistory

LEVELS = (0.0, 1.0, -0.4, 0.7, -1.0, 0.3, -0.6, 0.0)


def tension_torsion_block(
    time: np.ndarray,
    *,
    sigma_amplitude: float = 180.0,
    tau_ratio: float = 0.5,
    sigma_mean: float = 40.0,
    levels: Sequence[float] = LEVELS,
) -> np.ndarray:
    """Proportional history that visits the load-factor ``levels`` in order, one level per segment.

    The axial stress is ``sigma_mean + sigma_amplitude * lambda(t)`` and the shear stress
    ``tau_ratio * sigma_amplitude * lambda(t)``, both in MPa; ``lambda`` is piecewise linear in time.
    Returns the stress tensors, shape ``(len(time), 3, 3)``.
    """
    knots = np.linspace(time[0], time[-1], len(levels))
    lam = np.interp(time, knots, levels)
    S = np.zeros((len(time), 3, 3))
    S[:, 2, 2] = sigma_mean + sigma_amplitude * lam
    S[:, 0, 2] = S[:, 2, 0] = tau_ratio * sigma_amplitude * lam
    return S


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example; writes ``history.npz``, ``history.json`` and ``result.json`` to ``--out``."""
    parser = argparse.ArgumentParser(
        prog="python -m multifat.history.example", description=main.__doc__
    )
    parser.add_argument("--out", required=True, help="output directory (created if missing)")
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    time = np.linspace(0.0, 70.0, 281)
    history = TimeHistory.from_function(
        tension_torsion_block, time, sigma_amplitude=180.0, tau_ratio=0.5
    )
    history.record(out / "history")

    (model,) = get_models(["FIN"])
    curve = SNCurve(A=900.0, b=-0.12, method="FIN", source="illustrative curve, not a calibration")
    result = assess_history(
        history,
        model,
        {"k_FI": 0.3},
        curve,
        R_list=build_R_list(n1=3, n2=6),
        residue="rearrange",
    )
    result.to_json(out / "result.json")
    print(
        f"{len(result.events)} distinct events, {result.cycles_per_block:g} cycles per block, "
        f"damage per block {result.damage_per_block:.3e}, blocks to failure {result.blocks_to_failure:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
