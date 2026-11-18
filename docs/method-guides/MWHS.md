# Modified Wu–Hu–Song (`MWHS`)

> **Registry name:** `MWHS` (and `MWHS_ext`) · **Catalogue entry:** [N28](../../METHODS.md#n28) ·
> **Source:** Wu, Li, Fang & Song (2017), <https://doi.org/10.1177/0954406217708214> ·
> **Implementation:** [`cp_methods/mwhs.py`](../../src/multifat/models/cp_methods/mwhs.py)

## Overview

- **Type and family:** strain-based critical-plane parameter that combines the Fatemi–Socie and
  Smith–Watson–Topper ideas.
- **Idea:** on the plane of maximum shear strain, the shear strain amplitude, which drives shear
  cracking, is added to an SWT-like normal term, $\sqrt{\sigma_{n,\max}\Delta\varepsilon_n/E}$, which
  drives tensile opening. The normal term is weighted by $\tfrac12(1 + \sigma_{n,\max}/\sigma_y)$.
  This replaces the life-dependent weight $k$ of the earlier WHS parameter, so that only torsional
  fatigue properties are needed.
- **Proposed for:** low-cycle multiaxial fatigue of the Ni-based superalloy GH4169 under
  proportional and non-proportional tension–torsion. It is also applied to pure titanium,
  titanium BT9, S45C steel and 1045 steel.
- **Status:** implemented for the declared repository scope. Eq. (19) was checked on the paper,
  and Fig. 14 is reproduced (see [Verification](#verification-and-limitations)). The life
  comparisons of the paper are figures only.

## Original formulation

Wu et al. (2017, Eq. 19):

```math
\frac{\Delta\gamma_{\max}}{2} + \frac12\left(1 + \frac{\sigma_{n,\max}}{\sigma_y}\right)
\left(\frac{\sigma_{n,\max}\,\Delta\varepsilon_n}{E}\right)^{0.5}
= \frac{\tau_f'}{G}(2N_f)^{b_0} + \gamma_f'(2N_f)^{c_0},
```

evaluated on the maximum shear strain plane, where:

- $\Delta\gamma_{\max}/2$ is the maximum engineering shear strain amplitude;
- $\sigma_{n,\max}$ and $\Delta\varepsilon_n$ are the maximum normal stress and the normal strain
  range on that plane;
- $\sigma_y$ is the yield strength.

The earlier WHS parameter, Eq. (13), weighted the normal term by a constant $k$. Eq. (18) shows
that $k$ must vary with life (Fig. 12); the paper suggests 0.6 for GH4169. The modification removes
$k$.

## Repository implementation

On every candidate plane (normal = local $z$ axis):

```math
\mathrm{DP} = \gamma_a + \frac12\left(1 + \frac{\sigma_{n,\max}}{S_y}\right)\sqrt{\frac{\sigma_{n,\max}\,\Delta\varepsilon_n}{E}},
\qquad \sigma_{n,\max} = \max(\sigma_{n,m} + \sigma_{n,a}, 0).
```

DP is dimensionless.

- `MWHS`: DP on the plane of largest shear strain range; among tied planes, the larger DP.
- `MWHS_ext`: the largest DP over all planes. This is a repository variant.

**Repository choices**

- *Out-of-phase shear range.* The longest chord of the shear-strain path
  (`np_shear_rule = "longest_chord"`). With it, the maximum shear strain plotted in Fig. 14(a) is
  reproduced at every phase angle; the `rss` rule is not.
- *Compression clipped.* $\sigma_{n,\max} < 0$ is set to 0. The square root needs a non-negative
  product, and a compressed plane then carries no normal term.
- *Monotonic yield strength and elastic strains.* The workbook `Sy` is used for $\sigma_y$, and
  strains are elastic. The source applies the parameter to measured elastic–plastic strains in the
  low-cycle regime.
- *Life.* A power-law master curve $\mathrm{DP} = A N_f^b$ replaces the torsional Coffin–Manson
  curve of Eq. (19). As in the source, it is fitted on the **fully reversed pure torsion failures**
  (`design_loading = "torsion"`), not on the tension tests used by the other methods
  ([README](README.md#from-damage-parameter-to-life-the-master-curve)). The single power law
  has no separate plastic term.

## Material constants and inputs

| Symbol | Python / workbook name | Unit | Source | Required |
| --- | --- | --- | --- | --- |
| $\sigma_y$ | `Sy` / Summary yield strength | MPa | workbook | yes, $> 0$ |
| $E$, $\nu$ | `E`, `nu` / Summary | MPa, – | workbook | yes |
| $\boldsymbol\sigma(t)$ | `S_mean`, `S_sin`, `S_cos` | MPa | Dataset sheet | yes |

## Calibration

MWHS has no fitted criterion constant; this is its stated advantage over WHS. The life follows
from torsion tests alone.

1. **Original procedure.** The torsional Coffin–Manson constants $\tau_f'$, $b_0$, $\gamma_f'$, $c_0$,
   fitted on pure torsion tests, give the life (Eq. 19). The uniaxial tests then validate the
   parameter (Fig. 13).
2. **Repository.** The master curve $A$, $b$ is fitted on the fully reversed pure torsion failures
   of the workbook, which are reported as design rows (`Flag = 0`). The uniaxial tension tests are
   predictions (`Flag = 1`), as in the paper. MWHS is skipped, with the reason in the run report,
   when the workbook has no such torsion failures or when their master curve is not identifiable
   under the shared rule of `eval/fitting.py` (at least 5 tests, at least one decade of life, and
   a negative slope beyond two standard errors).

MWHS therefore uses a design set different from the other methods. Comparisons restricted to
common observations (`--common-methods`) report this, and should be run with `--flag_filter other`.

With elastic strains, the uniaxial parameter on the 45° plane is

```math
\mathrm{DP} = (1+\nu)\frac{\sigma_a}{E} + \frac12\left(1 + \frac{\sigma_a}{2S_y}\right)\frac{\sigma_a}{E}\sqrt{\frac{1-\nu}{2}},
```

and the torsional parameter is $\mathrm{DP} = 2(1+\nu)\tau_a/E$. The implied torsion-to-tension
limit ratio therefore depends on $\sigma_a/S_y$.

## Worked example

GH4169 elastic constants and yield strength (paper's Table 3) with $\nu = 0.3$; synthetic loads.

```python
import numpy as np

from multifat.models.registry import get_models
from multifat.planes.sampling import build_R_list

E, nu, Sy = 198500.0, 0.3, 1083.1  # MPa, -, MPa: GH4169 E and sigma_y (paper's Table 3)
mwhs = get_models(["MWHS"])[0]
params = {"Sy": Sy, "E": E, "nu": nu}
planes = build_R_list(n1=20, n2=40)
zero = np.zeros((3, 3))


def hooke(S):
    """Isotropic linear-elastic strain tensor (tensorial shear components)."""
    return ((1 + nu) * S - nu * np.trace(S) * np.eye(3)) / E


def evaluate(S_sin, S_cos=zero):
    """MWHS and MWHS_ext of S(t) = S_sin sin(wt) + S_cos cos(wt), elastic strains."""
    result = mwhs.evaluate_case(
        S0=S_sin, S1=-S_sin, E0=hooke(S_sin), E1=-hooke(S_sin), params=params, R_list=planes,
        harmonics={"S_mean": zero, "S_sin": S_sin, "S_cos": S_cos,
                   "E_sin": hooke(S_sin), "E_cos": hooke(S_cos)},
    )  # fmt: skip
    return result.values[""], result.values["_ext"]


s_a, t_a = 600.0, 400.0
tension, tension_ext = evaluate(np.diag([s_a, 0.0, 0.0]))
torsion_tensor = np.zeros((3, 3))
torsion_tensor[0, 1] = torsion_tensor[1, 0] = t_a
torsion, _ = evaluate(torsion_tensor)

# 45 deg plane in tension: gamma_a = (1+nu) s/E, sigma_n,max = s/2, d_eps_n = (1-nu) s/E
hand = (1 + nu) * s_a / E + 0.5 * (1 + s_a / 2 / Sy) * (s_a / E) * np.sqrt((1 - nu) / 2)
assert abs(tension - hand) < 1e-9 * hand
# torsion: gamma_a = tau/G on the plane of tau_xy, no normal stress there
assert abs(torsion - 2 * (1 + nu) * t_a / E) < 1e-9 * torsion

# 90 deg out-of-phase tension-torsion, sigma_a = 2 tau_a
oop, oop_ext = evaluate(np.diag([400.0, 0.0, 0.0]), 200.0 * (torsion_tensor / t_a))

print(f"MWHS tension 600 MPa: {tension:.6f} (MWHS_ext {tension_ext:.6f})")
print(f"MWHS torsion 400 MPa: {torsion:.6f}")
print(f"MWHS 90 deg out of phase: {oop:.6f} (MWHS_ext {oop_ext:.6f})")
```

Output:

```text
MWHS tension 600 MPa: 0.005071 (MWHS_ext 0.005556)
MWHS torsion 400 MPa: 0.005239
MWHS 90 deg out of phase: 0.004571 (MWHS_ext 0.004571)
```

In tension, `MWHS_ext` finds a larger DP on a plane closer to the load axis than the 45° plane of
maximum shear strain, as for [FS](FS.md) and [MKBM](MKBM.md).

**On a supplied workbook:**

```bash
uv run multifat-pipeline --material Al7075_T6 --models MWHS --all-variants --results-dir Results/current/guides
```

## Comparison capabilities and outputs

- **Life:** master curve (`Nf_expected`) fitted on the torsion failures, not the torsional
  strain–life curve with its plastic term.
- **DP comparison:** supported. DP is not homogeneous in the load, so use `Error_load` to compare
  with other methods.
- **`_ext`:** available; maximum DP over planes.

## Verification and limitations

[`tests/methods/test_mwhs.py`](../../tests/methods/test_mwhs.py) checks Fig. 14 for the GH4169 test
with $\varepsilon_a = 0.612$ %, $\gamma_a = 0.88$ %, $\sigma_a = 788.6$ MPa and $\tau_a = 404$ MPa.
The evaluation uses the measured strain amplitudes, $\nu = 0.48$ (Table 3) and the stresses of
the test:

- (a), the maximum shear strain amplitude, is reproduced within 2 % at all seven phase angles,
  0°–90°;
- (b), the normal term of Eq. (19), is reproduced within 5 % from 0° to 75°. At 90° the
  evaluation gives 0.0039 against 0.0049 read from the figure; this is unresolved and kept as a
  strict expected failure.

The module also checks the uniaxial and torsional closed forms with elastic strains, and the
clipping of compression.

Not reproduced: the life correlations of Figs. 15–19, which are figures only and use measured
elastic–plastic strains and the torsional Coffin–Manson curve.

Limitations:

- elastic strains on the workbooks, while the parameter was proposed for low-cycle fatigue;
- monotonic instead of cyclic yield strength.

## References

- Z. R. Wu, X. Li, L. Fang, Y. D. Song (2017), Evaluation of multiaxial fatigue life prediction
  criteria for Ni-based superalloy GH4169, *Proceedings of the IMechE, Part C: Journal of
  Mechanical Engineering Science*. <https://doi.org/10.1177/0954406217708214>
- Z. R. Wu, X. T. Hu, Y. D. Song (2014), Multiaxial fatigue life prediction for titanium alloy TC4
  under proportional and nonproportional loading, *International Journal of Fatigue* 59, 170–175.
