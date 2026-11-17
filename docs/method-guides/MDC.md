# MDC — Margetin, Ďurka & Chmelko (2016)

**Registry name:** `MDC`  
**Family:** Critical plane, stress, native life  
**Source file:** [`src/multifat/models/cp_methods/mdc.py`](../../src/multifat/models/cp_methods/mdc.py)

## Overview

`MDC` is a critical-plane stress criterion that defines an equivalent shear stress $\tau_r$ on the maximum-shear-stress plane and maps it directly onto the torsional Basquin $S$–$N$ curve for native life prediction. A stress-level-dependent normal-stress weighting factor $k$ interpolates smoothly between the axial and torsional Basquin curves over the entire life regime without a single fixed material ratio.

## Original formulation

In Margetin, Ďurka & Chmelko (2016), Eqs. (4)–(8), the equivalent shear stress and native life relation on the maximum-shear-stress-amplitude plane $(\sigma_{a,\mathrm{cr}}, \tau_{a,\mathrm{cr}})$ are

```math
\tau_r = \sqrt{k\,\sigma_{a,\mathrm{cr}}^2 + \tau_{a,\mathrm{cr}}^2} = \tau'_f (2 N_f)^{b_0}
```

with uniaxial and torsional Basquin curves $\sigma_a = \sigma'_f (2 N_f)^b$ and $\tau_a = \tau'_f (2 N_f)^{b_0}$ and stress-level factor

```math
k = \left(\frac{2\tau'_f}{\sigma'_f}\right)^{\!2}
    \left(\frac{2\,\sigma_{a,\mathrm{cr}}}{\sigma'_f}\right)^{\!\frac{2(b_0 - b)}{b}} - 1.
```

Because pure uniaxial stress $\sigma_a$ yields $\sigma_{a,\mathrm{cr}} = \tau_{a,\mathrm{cr}} = \sigma_a / 2$ on $45^\circ$ planes, the factor $k$ forces the criterion to match $\sigma_a = \sigma'_f (2 N_f)^b$ exactly; pure torsion ($\sigma_{a,\mathrm{cr}} = 0$) recovers $\tau_a = \tau'_f (2 N_f)^{b_0}$ independently of $k$.

## Repository implementation

**Critical plane selection.** The plane maximizing $\tau_{a,\mathrm{cr}}$ is chosen over the half-sphere grid; equal shear-amplitude ties are broken in favour of the larger $\tau_r$. Out-of-phase cycles inherit the default `rss` shear-amplitude projection (`np_shear_rule = "rss"`). For extended evaluation (`MDC_ext`), maximum-damage plane selection is also registered (`ext_variant = True`).

**Repository choices and departures.** Eqs. (7)–(8) in the paper print $\sigma_{\max}$ in place of the $\sigma_{a,\mathrm{cr}}$ of Eq. (4). Because all test series in Margetin et al. (2016) are fully reversed and using the applied axial stress misses Table 4 by factors of 2–70, the repository evaluates $\sigma_{\max} = \sigma_{a,\mathrm{cr}}$ (zero mean-stress sensitivity), reproducing every $N_{f,\mathrm{com}}$ entry of Table 4 within $0.3\%$ in $\tau_r$. When $k < 0$ at low normal stress yields a negative radicand, $k\,\sigma_{a,\mathrm{cr}}^2 + \tau_{a,\mathrm{cr}}^2$ is clamped to $0$.

## Material constants and inputs

| Symbol | Workbook key | Python name | Unit | Notes |
|---|---|---|---|---|
| $\sigma'_f$ | `sigma_f` | `sigma_f` | MPa | Axial fatigue strength coefficient |
| $b$ | `b_fat` | `b_fat` | — | Axial Basquin exponent ($< 0$) |
| $\tau'_f$ | `tau_f` | `tau_f` | MPa | Torsional fatigue strength coefficient |
| $b_0$ | `b0` | `b0` | — | Torsional Basquin exponent ($< 0$) |

All four stress-life Basquin parameters are required (`required_params`). Strains `E0`, `E1` are not used.

## Calibration

No additional optional parameters are exposed or calibrated (`optional_params = {}`). Material Basquin parameters are read directly from the workbook, or fitted by `multifat-pipeline` on the workbook's $R=-1$ tension and torsion tests when absent; all four are required, otherwise the run stops with an explicit error. The Uniform Material Law fallback would set $b_0 = b$ and make $k$ constant, so it is not used.

## Worked example

Evaluate the combined bending–torsion state $\sigma_a = 220\,\text{MPa}$, $\tau_a = 120\,\text{MPa}$ with the steel constants of Margetin et al. (2016) Table 2:

```python
import numpy as np
from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.mdc import MDC, mdc_life
from multifat.planes.sampling import build_R_list

p = {"sigma_f": 636.0, "b_fat": -0.0531, "tau_f": 550.0, "b0": -0.0736}
S = np.array([[220.0, 120.0, 0.0], [120.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
Z = np.zeros((3, 3))
res = best_cp_for_case(
    S0=S, S1=-S, E0=Z, E1=Z, R_list=build_R_list(), method=MDC(), params=p
)
tau_r = res.by_metric
nf = mdc_life(tau_r, p["tau_f"], p["b0"])
assert np.isclose(tau_r, 174.2036, atol=1e-3)
print(f"tau_r = {tau_r:.3f} MPa, Nf = {nf:.0f}")
```

Output:

```text
tau_r = 174.204 MPa, Nf = 3040907
```

Batch CLI evaluation:

```bash
multifat-run --model MDC
```

## Comparison capabilities and outputs

`MDC` sets `ComparisonCapabilities(life_prediction="native", dp_comparison=True, dp_reference="model_law")`. The reported damage parameter is $\tau_r$ in MPa, and predicted life $N_f$ comes directly from inverting $\tau_r = \tau'_f (2 N_f)^{b_0}$ via `mdc_life` without fitting a separate master curve.

## Verification and limitations

- **Software verification (paper reproduction).** Pure axial and torsional boundary conditions match exact Basquin identity to `1e-12`, and every combined bending–torsion row in Table 4 of Margetin et al. (2016) matches within $0.3\%$ in $\tau_r$ (life ratios $0.96$–$1.005$ at slope $1/b_0 = -13.6$) in [`tests/methods/test_mdc.py`](../../tests/methods/test_mdc.py).
- **Mean-stress limitation.** Taking $\sigma_{\max} = \sigma_{a,\mathrm{cr}}$ suppresses mean normal stress sensitivity; asymmetric ($R \neq -1$) tests are not differentiated by normal mean stress.

## References

- Margetin, M., Ďurka, R., & Chmelko, V. (2016). Multiaxial fatigue criterion based on parameters from torsion and axial S-N curve. *Frattura ed Integrità Strutturale*, 37, 146–152. [doi:10.3221/IGF-ESIS.37.20](https://doi.org/10.3221/IGF-ESIS.37.20)
- Shared conventions and symbol definitions: [`docs/method-guides/README.md`](README.md), [`docs/notation.md`](../notation.md), [`docs/method-audit.md`](../method-audit.md).
