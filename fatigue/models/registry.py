# fatigue/models/registry.py
from __future__ import annotations

"""Central registry for all *case-level* fatigue models.

Design goal: keep a **single** place where new methods are wired into the
user-facing API (``get_models``), while keeping the implementation files
themselves isolated in their respective subfolders:

  - ``fatigue/models/cp_methods``
  - ``fatigue/models/invariant_methods``
  - ``fatigue/models/energy_based_methods`` (empty for now)
"""

from typing import Iterable

from .cp_adapter import CPModelAdapter
from .base import FatigueModel

# CP (critical-plane) plane-level methods
from .cp_methods.fs import FatemiSocie
from .cp_methods.fin import Findley
from .cp_methods.swt import SWT
from .cp_methods.dangvan import DangVan
from .cp_methods.car_spa import CarSpa
from .cp_methods.matake import Matake
from .cp_methods.mcdiarmid import McDiarmid

# New CP methods from attached papers
from .cp_methods.mgse_yu import MGSE_Yu
from .cp_methods.mgse_zhu import MGSE_Zhu
from .cp_methods.mswt import MSWT
from .cp_methods.mkbm import MKBM
from .cp_methods.liu2021 import Liu2021

# New CP methods from attached papers (Li 2021, Ince & Glinka 2014)
from .cp_methods.li import LI
from .cp_methods.gse import GSE
from .cp_methods.gsa import GSA

# Invariant methods
from .invariant_methods.ottosen import OttosenReduced

# Energy based methods
from .energy_based_methods.swtd import SWTD
from .energy_based_methods.zhu_edp import ZhuEDP

def _available_cp_methods():
    """Instantiate all available critical-plane methods."""

    return {
        "FS": FatemiSocie(),
        "FIN": Findley(),
        "SWT": SWT(),
        "DANGVAN": DangVan(),
        "CARSPA": CarSpa(),
        "MATAKE": Matake(),
        "MCD": McDiarmid(),

        # Added methods
        "MGSE_YU": MGSE_Yu(),
        "MGSE_ZHU": MGSE_Zhu(),
        "MSWT": MSWT(),
        "MKBM": MKBM(),
        "LIU2021": Liu2021(),

        # Added in v6+ patch (critical-plane)
        "LI": LI(),
        "GSE": GSE(),
        "GSA": GSA(),
    }


def _available_direct_models():
    """Instantiate all available non-CP (direct / invariant) models."""

    return {
        "OTT": OttosenReduced(),
        "SWTD": SWTD(),
        "ZHU_EDP": ZhuEDP(),
        # "TRESCA": TrescaEqvStressAmp(),
        # "W_MAX": ElasticStrainEnergyMax(),
    }


def get_models(names: Iterable[str]):
    """Build a list of models by name.

    Parameters
    ----------
    names:
        Iterable of model names. Can contain:
          - CP method names: "FS", "FIN", "SWT", "DANGVAN", "CARSPA"
          - direct model names: "VM", "SWTD" (and future direct models)

    Returns
    -------
    list
        A list of case-level model instances.

    Notes
    -----
    CP names are wrapped into a case-level :class:`~fatigue.models.base.FatigueModel`
    via :class:`~fatigue.models.cp_adapter.CPModelAdapter`.
    """

    direct = _available_direct_models()
    cp = _available_cp_methods()

    out: list[FatigueModel] = []
    for n in names:
        if n in direct:
            out.append(direct[n])
        elif n in cp:
            out.append(CPModelAdapter(cp[n]))
        else:
            raise KeyError(
                f"Unknown model '{n}'. "
                f"Direct={list(direct.keys())}, CP={list(cp.keys())}"
            )

    return out
