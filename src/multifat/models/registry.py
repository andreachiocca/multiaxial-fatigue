"""Central registry for all *case-level* fatigue models.

Design goal: keep a **single** place where new methods are wired into the
user-facing API (``get_models``), while keeping the implementation files
themselves isolated in their respective subfolders:

  - ``multifat/models/cp_methods``
  - ``multifat/models/invariant_methods``
  - ``multifat/models/energy_based_methods``
  - ``multifat/models/integral_methods``
"""

from __future__ import annotations

import functools
from typing import Iterable

from .cp_adapter import CPModelAdapter
from .base import FatigueModel
from .comparison import DEFAULT_CAPABILITIES
from .references import ensure_enabled, method_reference

# CP (critical-plane) plane-level methods
from .cp_methods.fs import FatemiSocie
from .cp_methods.fin import Findley
from .cp_methods.swt import SWT
from .cp_methods.matake import Matake
from .cp_methods.mcdiarmid import McDiarmid
from .cp_methods.mdc import MDC
from .cp_methods.mwcm import MWCM
from .cp_methods.mmccm import MMCCM
from .cp_methods.liu_vse import LiuVSE

# Additional critical-plane formulations
from .cp_methods.mgse_yu import MGSE_Yu
from .cp_methods.mgse_zhu import MGSE_Zhu
from .cp_methods.mswt import MSWT
from .cp_methods.mkbm import MKBM
from .cp_methods.liu2021 import Liu2021

# Li (2021), Ince & Glinka (2014)
from .cp_methods.li import LI
from .cp_methods.gse import GSE
from .cp_methods.gsa import GSA

# Case-level integral and hybrid methods
from .integral_methods.bohme_papuga import BohmePapuga
from .integral_methods.caim import CAIM
from .energy_based_methods.swtd import SWTD
from .invariant_methods.ottosen import OttosenReduced
from .invariant_methods.tveit2025 import Tveit2025
from .cp_methods.dangvan import DangVan
from .invariant_methods.crossland import Crossland
from .invariant_methods.sines import Sines
from .integral_methods.papadopoulos import Papadopoulos
from .integral_methods.papuga_ruzicka import PapugaCriticalPlane, PapugaIntegral
from .integral_methods.robert_fogue_bahuaud import RobertFogueBahuaud
from .invariant_methods.vu_halm_nadot import VuHalmNadot
from .invariant_methods.pbp import ProjectionByProjection
from .invariant_methods.coin_liwi import CoinLiwi
from .invariant_methods.mmmh import MMMH
from .integral_methods.morel import Morel
from .energy_based_methods.ellyin_kujawski import EllyinKujawski


def _available_cp_methods():
    """Instantiate all available critical-plane methods."""

    return {
        "FS": FatemiSocie(),
        "FIN": Findley(),
        "SWT": SWT(),
        "MATAKE": Matake(),
        "MCD": McDiarmid(),
        "MDC": MDC(),
        "MWCM": MWCM(),
        "MMCCM": MMCCM(),
        "VSE": LiuVSE(),
        "MGSE_YU": MGSE_Yu(),
        "MGSE_ZHU": MGSE_Zhu(),
        "MSWT": MSWT(),
        "MKBM": MKBM(),
        "LIU2021": Liu2021(),
        "LI": LI(),
        "GSE": GSE(),
        "GSA": GSA(),
    }


def _available_direct_models():
    """Instantiate all available non-CP (direct / invariant) models."""

    return {
        "BP": BohmePapuga(),
        "CAIM": CAIM(),
        "SWTD": SWTD(),
        "OTT": OttosenReduced(),
        "TVEIT2025": Tveit2025(),
        "DANGVAN": DangVan(),
        "CROSSLAND": Crossland(),
        "SINES": Sines(),
        "PAPADOPOULOS": Papadopoulos(),
        "PCR": PapugaCriticalPlane(),
        "PI": PapugaIntegral(),
        "RFB": RobertFogueBahuaud(),
        "VHN": VuHalmNadot(),
        "PBP": ProjectionByProjection(),
        "COIN": CoinLiwi(),
        "MMMH": MMMH(),
        "MOREL": Morel(),
        "EK": EllyinKujawski(),
    }


def declared_param_names() -> set[str]:
    """Every material constant a registered model requires or optionally reads."""
    names: set[str] = set()
    for m in list(_available_cp_methods().values()) + list(_available_direct_models().values()):
        names |= set(m.required_params()) | set(getattr(m, "optional_params", {}))
    return names


def param_table() -> list[tuple[str, list[str], dict[str, str]]]:
    """(model, required params, {optional param: default or fallback}) for each registered model."""
    rows = []
    for name, m in {**_available_cp_methods(), **_available_direct_models()}.items():
        rows.append((name, sorted(m.required_params()), dict(getattr(m, "optional_params", {}))))
    return rows


def get_models(names: Iterable[str]):
    """Build a list of models by name.

    Parameters
    ----------
    names:
        Iterable of model names. Can contain:
          - CP method names: "FS", "FIN", "SWT", "MATAKE"
          - direct model names: "BP", "CAIM"

    Returns
    -------
    list
        A list of case-level model instances.

    Notes
    -----
    CP names are wrapped into a case-level :class:`~multifat.models.base.FatigueModel`
    via :class:`~multifat.models.cp_adapter.CPModelAdapter`.
    """

    direct = _available_direct_models()
    cp = _available_cp_methods()

    out: list[FatigueModel] = []
    for n in names:
        ensure_enabled(n)
        if n in direct:
            out.append(direct[n])
        elif n in cp:
            out.append(CPModelAdapter(cp[n]))
        else:
            raise KeyError(
                f"Unknown model '{n}'. " f"Direct={list(direct.keys())}, CP={list(cp.keys())}"
            )

    for model in out:
        model.reference = method_reference(model.name)
        if not hasattr(model, "comparison_capabilities"):
            model.comparison_capabilities = DEFAULT_CAPABILITIES
    return out


# Nominal constants for probing a model; only finiteness of its value is used.
_PROBE_PARAMS = {
    "Sy": 700.0, "Su": 900.0, "Sigm1": 400.0, "Taum1": 250.0, "Sig0": 300.0, "Tau0": 200.0,
    "k_FS": 1.0, "k_FI": 0.3, "L_LI": 3.0, "k_MGSE": 1.0, "a_MSWT": 0.3, "alpha_w": 0.0,
    "xi_ZHU": 1.0, "sigma_f": 1200.0, "b_fat": -0.09, "epsilon_f": 0.5, "c_fat": -0.6,
    "tau_f": 700.0, "gamma_f": 0.9, "b0": -0.09, "c0": -0.6,
    "Sigm1_TVEIT": 400.0, "Taum1_TVEIT": 250.0,
    "npts_TVEIT": 1001.0, "min_blocks_TVEIT": 8.0, "max_blocks_TVEIT": 10000.0,
    "tol_alpha_TVEIT": 1e-6, "tol_damage_TVEIT": 1e-5,
}  # fmt: skip


@functools.lru_cache(maxsize=1)
def stress_only_models() -> frozenset[str]:
    """Base models whose value needs no strain, so they run on a workbook without E or nu.

    Found by evaluation rather than listed: each model sees one 90-degree tension-torsion case
    with NaN strains and E = nu = None, and every other declared constant at a nominal value.
    """
    import warnings

    import numpy as np

    from multifat.planes.sampling import build_R_list

    nan = np.full((3, 3), np.nan)
    z = np.zeros((3, 3))
    s_sin = np.array([[0.0, 0.0, 0.0], [0.0, 200.0, 0.0], [0.0, 0.0, 0.0]])
    s_cos = np.array([[0.0, 100.0, 0.0], [100.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    harm = {"S_mean": z, "S_sin": s_sin, "S_cos": s_cos, "E_mean": nan, "E_sin": nan, "E_cos": nan}
    R_list = build_R_list(n1=6, n2=12)
    out = set()
    for name, model in {**_available_cp_methods(), **_available_direct_models()}.items():
        m = CPModelAdapter(model) if name in _available_cp_methods() else model
        params = {k: _PROBE_PARAMS.get(k, 1.0) for k in declared_param_names()}
        params.update(E=None, nu=None)
        try:
            with warnings.catch_warnings(), np.errstate(all="ignore"):
                warnings.simplefilter("ignore")
                res = m.evaluate_case(
                    S0=s_sin, E0=nan, S1=-s_sin, E1=nan, params=params, R_list=R_list,
                    harmonics=harm,
                )  # fmt: skip
            v = res.values.get("")
            if v is not None and np.isfinite(v):
                out.add(name)
        except Exception:
            pass
    return frozenset(out)
