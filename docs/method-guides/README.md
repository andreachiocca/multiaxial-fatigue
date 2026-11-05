# Method guides

One guide per implemented fatigue criterion. Each guide states the published formulation, what
the repository computes and where it departs from the source, which material constants are needed
and how they are obtained, a worked example that can be checked by hand, and the comparison
outputs the method supports.

The guides complement, and do not replace:

- the generated [method catalogue](../../METHODS.md) (status, evidence and DOI of every catalogued
  method, generated from `docs/catalogue/methods.json`);
- the [numerical method audit](../method-audit.md) (one line per implemented method and its audit
  decisions);
- [comparison capabilities](../comparison-capabilities.md) and the
  [calibration workflow](../calibration-workflow.md).

When a guide and the source code disagree, the code is what runs: please open an issue. Every
Python example in these guides is executed by `tests/test_method_guides.py`, so the printed
numbers are checked on every CI run.

## Available guides

The most widely used criteria are documented first.

| Guide | Registry name | Family | Damage parameter (DP) | Material constants |
| --- | --- | --- | --- | --- |
| [Findley](FIN.md) | `FIN`, `FIN_ext` | Critical plane, stress | $\tau_a + k\,\sigma_{n,\max}$ (MPa) | `k_FI` |
| [Fatemi–Socie](FS.md) | `FS`, `FS_ext` | Critical plane, strain | $\gamma_a\,(1 + k\,\sigma_{n,\max}/\sigma_y)$ (–) | `k_FS`, `Sy`, `E`, `nu` |
| [Smith–Watson–Topper](SWT.md) | `SWT`, `SWT_ext` | Critical plane, energy-like | $\sigma_{n,\max}\,\varepsilon_{n,a}$ (MPa) | `E`, `nu` |
| [Matake](MATAKE.md) | `MATAKE`, `MATAKE_ext` | Critical plane, stress | $\tau_a + (2t_{-1}/f_{-1} - 1)\,\sigma_{n,\max}$ (MPa) | `Sigm1`, `Taum1` |
| [McDiarmid](MCD.md) | `MCD`, `MCD_ext` | Critical plane, stress | $\tau_a + t_{-1}\sigma_{n,\max}/(2S_u)$ (MPa) | `Taum1`, `Su` |
| [Modified Wöhler Curve Method](MWCM.md) | `MWCM` | Critical plane, stress | $\tau_a/\tau_{A,\mathrm{Ref}}(\rho)$ (index) | `Sigm1`, `Taum1` |
| [Crossland](CROSSLAND.md) | `CROSSLAND` | Stress invariants | $(\sqrt{J_{2,a}} + \kappa\,\sigma_{H,\max})/t_{-1}$ (index) | `Sigm1`, `Taum1` |
| [Sines](SINES.md) | `SINES` | Stress invariants | $(\sqrt{J_{2,a}} + \kappa\,\sigma_{H,m})/t_{-1}$ (index) | `Sigm1`, `Taum1`, `Su` |
| [Dang Van](DANGVAN.md) | `DANGVAN` | Mesoscopic | $\sigma_{eq,DV}/f_{-1}$ (index) | `Sigm1`, `Taum1` |
| [Papadopoulos](PAPADOPOULOS.md) | `PAPADOPOULOS` | Mesoscopic, integral | $\gamma_N/t_{-1}$, Eq. (41) (index) | `Sigm1`, `Taum1` |
| [Böhme–Papuga](BP.md) | `BP` | Integral | $\sqrt{\langle a\tau_a^2 + b\sigma_{n,a}^2 + c\tau_a\tau_m + d\sigma_{n,a}\sigma_{n,m}\rangle}/s_{-1}$ (index) | `Sigm1`, `Taum1`, `Sig0`/`Tau0` or `Su` |
| [CAIM](CAIM.md) | `CAIM` | Critical plane + integral | $\sqrt{\max(a\tau_a^2 + b\sigma_{n,a}^2) + \langle c\tau_a\tau_m + d\sigma_{n,a}\sigma_{n,m}\rangle}/s_{-1}$ (index) | as BP |
| [Papuga–Růžička PCr](PCR.md) | `PCR` | Critical plane, stress | $\max\sqrt{a_c C_a^2 + b_c(N_a + t_{-1}N_m/f_0)}/f_{-1}$ (index) | `Sigm1`, `Taum1`, `Sig0` or `Su` |
| [Papuga–Růžička PI](PI.md) | `PI` | Integral | $\sqrt{\langle a_I C_a^2 + b_I(N_a + \kappa N_m)\rangle}/f_{-1}$ (index) | `Sigm1`, `Taum1` |
| [Morel](MOREL.md) | `MOREL` | Mesoscopic, native life | $T_S/T_{S,\lim}$ (index); life from Eq. (38) | `Sigm1`, `Taum1`; $p, q, r$ fitted |

Guides still to be written follow the same [template](TEMPLATE.md). Until a guide exists, use the
method's entry in the [catalogue](../../METHODS.md), its line in the
[method audit](../method-audit.md) and its module docstring under `src/multifat/models/`.

## What every method shares

### Load history and sign conventions

Every observation is a constant-amplitude, single-harmonic stress history

```math
\boldsymbol\sigma(t) = \boldsymbol\sigma_m + \boldsymbol\sigma_s \sin\omega t + \boldsymbol\sigma_c \cos\omega t
```

stored as the harmonic tensors `S_mean`, `S_sin` and `S_cos` (MPa). Strain tensors `E_mean`,
`E_sin`, `E_cos` follow from the stresses by isotropic linear elasticity with the workbook `E` and
`nu`; the workflow does not simulate cyclic plasticity. Tensor shear strains are tensorial
($\varepsilon_{ij}$); engineering shear strain is $\gamma_{ij} = 2\varepsilon_{ij}$. Symbols follow
[notation.md](../notation.md).

### Critical-plane search

Critical-plane criteria (`FIN`, `FS`, `SWT`, `MATAKE`, `MCD`, `MWCM`, …) are evaluated on a
522-plane grid (`build_R_list(n1=20, n2=40)`), completed with analytical seed planes and a local
refinement of the plane normal. On each plane the local $z$ axis is the normal, and the components
$xz$, $yz$ are the two in-plane shear components. For an out-of-phase history the shear vector on a
plane describes an ellipse with semi-axes $a \ge b$; the amplitude used is

| `np_shear_rule` | Shear amplitude | Used by |
| --- | --- | --- |
| `rss` (default) | $\sqrt{a^2 + b^2}$ | methods whose source gives no rule (e.g. FIN, MATAKE, MCD) |
| `longest_chord` | $a$ (half the longest chord, radius of the minimum circumscribed circle) | FS, MWCM, MKBM, LI, LIU2021, VSE, MMCCM |

For proportional loading $b = 0$ and the two rules coincide. See
[shear-amplitude.md](../shear-amplitude.md).

Each critical-plane method returns two values:

- `NAME`: the DP on the plane that maximises the method's **selection metric** (for example the
  shear strain range for FS); ties are resolved towards the larger DP;
- `NAME_ext`: the largest DP over all planes. This is a repository variant, not necessarily the
  author's plane rule. For methods whose metric is the DP itself (FIN) the two coincide.

### From damage parameter to life: the master curve

Unless a method declares its own life law (MDC, MMCCM, MOREL, TVEIT2025), the life is obtained
from a power-law master curve fitted to the **design set**: the uniaxial tension failures at
$R = -1$ (`Flag = 0` in the CSV, together with any other test that fitted a constant of the method),

```math
\mathrm{DP} = A\,N_f^{\,b}, \qquad \hat N_f = \left(\mathrm{DP}/A\right)^{1/b},
```

by least squares in $\ln \mathrm{DP}$ versus $\ln N_f$ (`multifat.eval.fitting`). All other failures
(`Flag = 1`) are predicted with the inverted curve. Run-outs (`Flag = 2`) are evaluated but never
enter a fit or an error statistic. This is an **assessment-layer** life law shared by all methods,
not the life law of the original author. The life error is
$\mathrm{Error\_ln} = \ln(N_{f,\exp}/\hat N_f)$ and the damage-space error is
$\mathrm{Error\_log10\_dp} = \log_{10}(\mathrm{DP\_fit}/\mathrm{CP\_value})$ (see
[comparison capabilities](../comparison-capabilities.md#error-statistics-and-observation-roles)).

### Calibration of method constants

Method constants come, in this order of precedence, from

1. the workbook `Summary` sheet (a column named after the constant, e.g. `k FI`);
2. a calibration record passed with `multifat-pipeline --constants FILE` (a file under
   `Calibrations/` or an earlier `<MATERIAL>_calibration.json`);
3. the calibration that `multifat-pipeline` runs for the constants that are still missing. It is
   written to `<MATERIAL>_calibration.json` in the results folder, never into the workbook.

The repository calibrates the one-parameter plane constants (`k_FS`, `k_FI`, `L_LI`, `k_MGSE`,
`a_MSWT`) by **minimising the scatter** of the uniaxial tension and torsion failures at $R=-1$
around one master curve: for each trial value, DP is evaluated on every selected test, the power
law above is fitted, and the standard deviation of the residuals in $\ln N_f$ is minimised by a
grid search with refinement. This makes tension and torsion lives collapse on a single curve at
finite life. It is generally **not** the calibration printed by the original authors, which
usually matches two fatigue limits; each guide states the difference. On a fatigue-limit dataset
(run-outs only), a missing plane constant is chosen so that the R = -1 tension and torsion limits
give the same DP (`multifat.eval.endurance`).

The tests of the scatter objective (tension and torsion) are reported as design rows
(`Flag = 0`) of every model that adopts the calibrated constant, as are the Basquin tests behind
constants derived from Basquin fits. With `--constants`, a sidecar of this pipeline carries its
design tests only when it was made on the same failures; for constants from the workbook or an
upstream record the tests are unknown and stay `Flag = 1`. The held-out validation protocol is
tracked in [issue #14](https://github.com/andreachiocca/multiaxial-fatigue/issues/14).

The default master curve is Basquin; `--sn-model kohout_vechet` or `stromeyer` adds a fatigue-limit
asymptote, with the caveats listed in the
[method audit](../method-audit.md#master-curve-s-n-models).

To rank methods against each other use the load-equivalent error `Error_load`, which does not
depend on the slope of the master curve nor on the degree of the DP
([comparison capabilities](../comparison-capabilities.md#load-equivalent-error-the-cross-method-ranking-metric)).

Calibration is not validation: an error statistic computed on tests that fed a constant or the
master curve measures the fit, not the predictive accuracy.

### Running one method

```bash
# Calibrate the missing constants, run FIN and FIN_ext on one material, compare.
uv run multifat-pipeline --material Al7075_T6 --models FIN --all-variants --results-dir Results/current/guides

# Reuse constants from a calibration record instead of calibrating.
uv run multifat-pipeline --material Al7075_T6 --models FIN --all-variants \
    --constants Results/current/guides/Al7075_T6_calibration.json --results-dir Results/current/guides_reuse
```

The results folder then holds:

| File | Content |
| --- | --- |
| `<MODEL>_<MATERIAL>.csv` | one row per test: `CP_value` (DP), `Nf_exp`, `Nf_expected`, `DP_fit`, `Flag`, errors, identities |
| `<MATERIAL>_calibration.json` | calibrated constants, objective, selected test IDs, provenance (workbook hash, git revision) |
| `<MATERIAL>_run_report.tsv` | status of each model (`ok`, `skipped` with the reason) |
| `<MATERIAL>_validation.txt` | workbook checks |
| `Comparison/` | error statistics, plots and LaTeX tables |
