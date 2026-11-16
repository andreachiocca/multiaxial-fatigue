# Harmonic shear amplitudes and controlled Findley evaluation

Provenance: upstream pull requests #7 through #9 (commits 88229e9, 176ce53, 36ff0a1), carried into `multifat` together with their tests: the utility (88229e9, 2a17791), the standalone Findley evaluation `fin_harmonic.py` (176ce53, 50dabaf), and the Findley calibration amplitude (36ff0a1, 0105b69). Import statements were altered, and both the calibration integration and the snapshot values in the tests were adapted, as detailed below.

Phase 1 introduces a numerical utility; no fatigue criterion, plane-selection rule, calibration record, experimental input or reporting format is changed. Neither the current evaluator nor the plane-history utilities or the Böhme–Papuga implementation import the new module. Phase 2 supplies a distinct entry point for the Findley damage parameter, intended for numerical comparisons; production predictions remain unaffected. Phase 3A makes the Findley calibration amplitude configurable (`multifat-calibrate`, discussed below). The production Findley assessment keeps evaluating MCE and keeps accepting only MCE-compatible constants.

## Definition and API

For a complete single-frequency trajectory expressed in a common two-dimensional tangent basis,

\[
\boldsymbol\tau(t)=\boldsymbol\tau_m+\mathbf a\sin t+\mathbf b\cos t,
\qquad A=[\mathbf a\;\mathbf b],
\]

the singular values \(s_1\geq s_2\) of \(A\) are the ellipse semiaxes. An amplitude expressed in the units supplied by the caller is what the API returns.

| Identifier | Returned amplitude | Harmonic equivalence |
|---|---|---|
| `mcc_harmonic` | \(s_1\) | Radius of the minimum circumscribed circle, also half the longest chord |
| `mce_harmonic` | \(\sqrt{s_1^2+s_2^2}=\sqrt{\|a\|^2+\|b\|^2}\) | RSS definition already used for harmonic amplitudes |
| `rss` | Alias for `mce_harmonic` | Identical numerical implementation |

```python
from multifat.models.utils.shear_amplitude import shear_amplitude

mcc = shear_amplitude([10, 0], [0, 5], algorithm="mcc_harmonic")  # 10
mce = shear_amplitude([10, 0], [0, 5], algorithm="mce_harmonic")  # sqrt(125)
```

Accepted inputs take shape `(2,)` or identical batch shapes `(..., 2)`. One coefficient pair yields a `float`; a batch yields an array whose last axis has been removed. Coefficients must be finite, real and representable in float64. Broadcasting between shapes, complex values and unrecognised identifiers are refused. `method_default` is intentionally excluded as an algorithm identifier, and no global user option is added.

Because the mean vector locates the ellipse centre, it is not an input to the amplitude. Its norm is a separate scalar quantity. In general, neither amplitude equals the largest instantaneous magnitude of the translated trajectory. Take a circle of radius 10 whose mean is `(30, -40)`: MCC is 10, MCE is \(10\sqrt2\), the mean magnitude is 50 and the maximum total magnitude is 60. Twice the selected amplitude gives a shear range. Strain conventions are not converted by the function; a caller needing an engineering shear measure must convert tensorial shear strain beforehand.

## Numerical formulation

MCC takes the larger eigenvalue of the 2 by 2 Gram matrix and never subtracts nearly equal quantities. Each input pair is first divided by its largest absolute coefficient; the square root of

\[
\lambda_+=\tfrac12\left(a\cdot a+b\cdot b+
\operatorname{hypot}(a\cdot a-b\cdot b,2a\cdot b)\right)
\]

is then multiplied by that scale. MCE instead uses the scaled Frobenius norm. Squaring the original coefficients is thereby avoided, and paths at very small or very large scales remain safe. A path counts as zero only when it is exactly zero. No absolute physical cutoff, eigenvalue clipping tolerance or iterative optimisation is involved.

Calculations run in float64. After normalisation, contributions lost to underflow cannot affect the leading amplitude at float64 precision. Subnormal outcomes keep float64 rounding, while a result outside the finite range raises `OverflowError`. Below the subnormal spacing, no relative accuracy can be guaranteed. Complete harmonic ellipses, together with their line and point limits, are covered by these algorithms. They are not enclosing-circle or enclosing-ellipse solvers for arbitrary sampled paths, partial cycles or multiple frequencies, and no estimate of the minor semiaxis is returned.

In a local comparison over 100,000 random matrices, with the median taken across five batches of three evaluations, this MCC implementation measured approximately 0.020 s, against 0.062 s for NumPy SVD. That observation is what motivates the simple closed form; it should not be read as a portable timing guarantee. Independently of the analytical examples, the tests compare both amplitudes with SVD on 1,000 matrices.

## Source interpretation and verification

The scientific source is the 2021 thesis by Eva Cízová, *Parameters of Stress Tensor History and Their Impact on Fatigue Life* (Czech Technical University in Prague), supervised by Jan Papuga. Equation numbers follow the searchable reference extracted from that thesis and supplied for this work. The complete PDF could not be obtained in this session; consequently no independent page-level or image-level verification of the equations is claimed.

- The longest-chord amplitude is defined in Eq. (2.1), while Eqs. (2.5)-(2.6) define MCC.
- The MCE semiaxis combination comes from Eq. (2.9). The supplied summary cites Eqs. (3.26), (3.27) and (3.29) for the harmonic ellipse identities.
- Section 2.5 presents a moment-of-inertia definition weighted by path length. That is a different amplitude definition and is not implemented here.

The MATLAB appendix was inspected directly. The algebraic `Ca.m` files it contains, both under `Programy_DP/Nová kritéria` and in the circle/ellipse directory, are byte-identical (SHA-256 `e96f70827fc0ce5c158ab56ade21af3a9eb8e4b72e47f205c5f3c93ab08092e5`). They define `Cl=f*sin(t)+g*cos(t)` and `Cr=p*sin(t)+q*cos(t)`; consequently the Python inputs map as `a=(f,p)` and `b=(g,q)`. The expression actually used returns the major semiaxis (MCC), and the commented `sqrt(a1^2+a2^2)` alternative corresponds to MCE. `MOI/Ca.m` works differently, using temporal quadrature weighted by path speed.

Two reconstructions of the source equations take `phi=0`, `theta=pi/4`, torsional phase `pi/2`, zero radial and hoop amplitudes, an axial amplitude of 20, and a torsional amplitude of \(10\sqrt2\) or \(5\sqrt2\). The corresponding coefficients are `a=(0,-10), b=(-10,0)` and `a=(0,-10), b=(-5,0)`. The expected MCC/MCE amplitudes are 10 / \(\sqrt{200}\) and 10 / \(\sqrt{125}\). Automated tests re-implement the MATLAB algebra in Python and match it to 13 decimal places. No MATLAB or Octave runtime was available, so the comparison is made directly at equation level instead of executing MATLAB. The appendix result file `FatLimResultsAll.mat` holds a 286 by 60 comparison array, whose criterion-level results were not treated as isolated amplitude benchmarks. None of the thesis experimental data are redistributed.

The PragTic `solve.cpp` file was inspected with respect to `MCCM_IM` and `MCEM_IM`. Those routines take sampled polygons and determine centres and amplitudes iteratively, with `MCEM_IM` combining the fitted axes through their quadratic norm. This phase neither ports nor validates their tolerances, convergence behaviour or general-path algorithms. Jan Papuga supplied the files for WG4.3 study and adaptation. The new implementation follows the harmonic identities stated above, and no external solver implementation or source archive has been copied into the package.

Focused checks are run with `uv run pytest tests/test_shear_amplitude.py`. They cover exact zero, proportional and degenerate paths, circles and ellipses, rotations and reflections, phase shifts, translated means, direct time sampling, longest chords, extreme scales, invalid inputs, batch evaluation, reconstructions of the source equations, and agreement with the existing RSS and BP helpers. The full suite is invoked by `uv run pytest`.

## Relation to the critical-plane shear rules

On every plane, the critical-plane evaluator already applies an ellipse-aware shear range, chosen per method through `np_shear_rule` (`multifat.eval.evaluator.SHEAR_RULES`). The `longest_chord` rule is twice the `mcc_harmonic` amplitude, whereas `rss` (the default) is twice `mce_harmonic`. `tests/test_shear_amplitude.py::test_longest_chord_rule_is_twice_mcc` verifies both on a single plane, and `tests/test_np_shear_rule.py` records which method uses which rule, and why. As in the upstream design, this utility remains independent of that code path.

## Controlled Findley evaluation

The method-specific function `evaluate_findley`, located in `src/multifat/models/cp_methods/fin_harmonic.py`, takes the case arguments already used by the CP adapter (`S0`, `E0`, `S1`, `E1`, `params`, `R_list`, `harmonics`) plus the optional keyword `shear_amplitude_algorithm`. It returns the existing `CaseResult`, whose damage parameters are stress values, not fatigue lives.

```python
from multifat.models.cp_methods.fin_harmonic import evaluate_findley

# case contains the existing CP adapter arguments, including params['k_FI'].
legacy = evaluate_findley(**case)
mce = evaluate_findley(**case, shear_amplitude_algorithm="mce_harmonic")
mcc = evaluate_findley(**case, shear_amplitude_algorithm="mcc_harmonic")
print(mcc.values[""], mcc.values["_ext"])  # FIN and FIN_ext damage parameters
```

When the selection is omitted (`None`), the call delegates exactly to the production FIN adapter, which is unchanged. The existing harmonic MCE/RSS convention is thus preserved, together with the existing two-endpoint calculation whenever harmonic coefficients are unavailable. `None` acts as an interface sentinel, not as a further numerical amplitude algorithm. Any explicit selection, including MCE and its `rss` alias, needs both `harmonics['S_sin']` and `harmonics['S_cos']` as finite, real 3 by 3 stress tensors. Coefficients that are partial or absent raise `ValueError`, since two endpoints cannot carry phase information. The identifier `method_default` is not accepted.

A single global basis is shared by the stress tensors and the coefficient tensors. Every orthonormal matrix in `R_list` stores the candidate-plane axes as its columns, the third column being the normal. Once transformed by \(R^T S R\), the local xz and yz entries of each coefficient tensor become the two input vectors passed to the shared amplitude function. Strain tensors keep the adapter interface but play no role in Findley. As in the production evaluator, the mean stress is defined by the midpoint `(S0 + S1)/2`; a distinct `S_mean` entry is never read. Explicit harmonic amplitudes are determined by the coefficient tensors, not by the synthetic endpoint shear range.

For **each** candidate plane, the function first evaluates \(P=\tau_a+k_{\mathrm{FI}}\sigma_{n,\max}\) and then maximises that complete expression. The normal-extrema helper already in place is reused unmodified, signed compression included. It returns \(\sigma_{n,\max}=\sigma_{n,m}+\sqrt{\sigma_{n,s}^2+\sigma_{n,c}^2}\) regardless of which shear algorithm was chosen. Mean shear plays no part in the amplitude, and no factor of two multiplies the amplitude inside the Findley expression. Since Findley itself selects planes according to damage, both output keys hold the same maximum; no artificial difference between FIN and FIN_ext is created.

With a prescribed plane whose mean normal stress is 7, whose normal sine and cosine coefficients are 3 and 4, and with \(k_{\mathrm{FI}}=0.2\), the normal contribution amounts to 2.4.

| Shear semiaxes | MCC amplitude | MCE amplitude | MCC Findley | MCE Findley |
|---|---:|---:|---:|---:|
| 10, 10 | 10 | 14.142136 | 12.4 | 16.542136 |
| 10, 5 | 10 | 11.180340 | 12.4 | 13.580340 |
| 10, 0 | 10 | 10 | 12.4 | 12.4 |

The automated two-plane example goes further and reverses the ranking. Along normal z, taking a circle of radius 10 with mean normal 7 yields MCC/MCE damage of 11.4/15.542136. Along normal x, a line of amplitude \(\sqrt{164}\) with mean normal 5 gives 13.806248 under either algorithm. MCC therefore selects normal x, while MCE selects normal z. Substituting the amplitude only after an MCE search would produce an incorrect MCC result. The tests also compare 30 random planes against independent SVD and normal-traction projections, and they verify default snapshots, phase and basis invariance, mean stresses, endpoint rejection, small and large shear scales, and isolation from other methods. Tests are run with `uv run pytest tests/test_findley_harmonic.py`. Between the supplied planes, the production adapter additionally refines the critical plane (F2); its default is therefore compared against the supplied-plane values (`refine=False`) and never falls below them. The non-proportional FS snapshot follows the longest-chord rule (F3) instead of the RSS value used upstream.

The search is confined to the finite plane set supplied. The input assumptions and arbitrary-path restrictions of the shared library continue to hold. At ordinary stress scales, explicit MCE and legacy harmonic MCE agree to within float64 rounding; the explicit route bypasses the legacy synthetic-shear endpoint cutoff and the mean subtraction. The omitted-selection path is left exactly as it was. Arithmetic in the normal-extrema helper keeps the existing float64 limits, overflow included, when normal coefficients are squared at extreme scales. These checks are numerical verifications; they do not constitute independent experimental validation of Findley with alternative amplitudes.

The separate entry point is deliberate. Upstream calibration fingerprints hash whole production source files; keeping the entry point separate leaves production results unchanged here. Neither the general adapter, nor the production registry, nor the material dictionary gains an algorithm parameter (the module is unregistered, and production code does not import it). The calibration and assessment options described below now cover every criterion that exposes a plane shear amplitude (see [Selection for every method](#selection-for-every-method)). The diagnostic function neither loads nor applies a damage-life curve. In particular, an MCE-calibrated curve must not be used with MCC damage parameters.

## Findley calibration selection

`multifat-calibrate --shear-amplitude-algorithm method_default|mcc_harmonic|mce_harmonic` (with `rss` also accepted as a canonicalised MCE alias) sets the shear amplitude used by the Findley criterion during k_FI calibration, through `calibrate_material_params(..., shear_amplitude_algorithm=...)`. The default follows the established MCE path, whose results are unchanged. An explicit selection requires the stored stress sine and cosine tensors for every uniaxial tension and torsion failure employed; inspecting the original case in this way avoids the proportional endpoint fallback of the general tensor extractor. Since November 2026 the selection applies to every calibrated critical-plane criterion (FS, FIN, LI, MGSE_ZHU, MSWT). With `method_default` the result is unchanged and records `amplitude_settings` for `FIN` only; an explicit selection records, for each calibrated criterion, the resolved algorithm, `DEFINITION_VERSION` and the evaluation path.

Upstream commits 36ff0a1 and 0105b69 perform calibration with a separate Findley model and a versioned JSON store. In the present implementation the choice is instead the plane shear rule of the established search (`np_shear_rule` in `multifat.eval.evaluator.SHEAR_RULES`). As in production Findley, MCE is `rss`, the root-sum-square of the shear semiaxes. MCC corresponds to `longest_chord`, half the longest chord of the shear path, that is, the major semiaxis (`FindleyMCC` in `multifat/eval/findley_calibration.py`). Both choices therefore share the refined plane search, the plane pools and the bounds of the calibration. `tests/test_findley_calibration.py` verifies the rules against `evaluate_findley` (MCC and MCE) on prescribed planes, including the reversed plane ranking reported above.

For axial and torsion calibration failures the shear path degenerates to a line; MCC and MCE are then equal and the fitted k_FI is the same (as in the upstream example). The choice becomes relevant for non-proportional calibration data. Production Findley (`multifat-run`, `multifat-pipeline`) evaluates MCE, so a k_FI calibrated with MCC must not be used there. The option is not exposed by `multifat-pipeline`, and `--constants` rejects any record whose `amplitude_settings.FIN.algorithm` is neither absent nor `mce_harmonic`. Predicting life with MCC in production would require an integration approved separately. Upstream's record format and numerical example are given in the [calibration workflow](calibration-workflow.md#findley-amplitude-selection). Neither the calibration agreement nor the software checks presented here establish independent experimental validation of Findley with alternative amplitudes.
## Explicit amplitude in production assessment

`multifat-run` takes `--shear-amplitude-algorithm {method_default,longest_chord,mcc_harmonic,mce_harmonic,rss}`
for every model that exposes a plane shear amplitude (see below). `method_default` (the default)
is the production path and its output is unchanged. An explicit selection calibrates a model's
constants with the same criterion it assesses with, so the two cannot disagree; `rss` resolves to
`mce_harmonic`, `longest_chord` to `mcc_harmonic`. Results are written as, e.g.,
`FIN.mcc_harmonic` / `FS.mce_harmonic` (`_ext` kept), so
`multifat-compare` reads them by name, and a `.meta.json` next to each CSV records the selected option,
the resolved algorithm, the evaluation path, the definition version and the SHA-256 of the CSV.

Cases without stored stress sine and cosine tensors cannot be assessed on the explicit path; they are
skipped with a `RuntimeWarning`. The option is not exposed by `multifat-pipeline`, which evaluates
production Findley (MCE) only.

## Selection for every method

Maintainer decision (16 November 2026): the user selects the plane shear amplitude of every method
that uses one; each method's default is the amplitude of its source, otherwise the longest chord.
`multifat.eval.shear_amplitude_option` implements the selection:

| Selection | Plane rule (`np_shear_rule`) | Single harmonic |
| --- | --- | --- |
| `method_default` | the method's own | as documented in each guide |
| `longest_chord`, `mcc_harmonic` | `longest_chord` | MCC radius = major semiaxis $a$ = half the longest chord |
| `mce_harmonic`, `rss` | `rss` | MCE amplitude $\sqrt{a^2 + b^2}$ |

For one harmonic the longest chord and the MCC coincide, so a method has at most two results. The
selection reaches every critical-plane method (through the plane preparation) and every case-level
model declaring `np_shear_rule` (SIH). Models whose amplitude is fixed by their definition (for
example BP and RFB use the MCC, PAPADOPOULOS the quadratic mean of its resolved shear, LIUMAHADEVAN
the rss by maintainer decision) or that use none (invariant criteria, WANG2018 with its scalar shear)
run unchanged, with a `RuntimeWarning`, and their outputs keep the plain name.

`multifat-calibrate --shear-amplitude-algorithm` fits FS, FIN, LI, MGSE_ZHU and MSWT with the
selection; `multifat-run --shear-amplitude-algorithm` assesses with it, taking calibrated constants
(e.g. k_FS) from the workbook, which must have been calibrated with the same option.
`multifat-pipeline` evaluates `method_default` only, and its `--constants` refuses a record whose
`amplitude_settings` hold another amplitude for any method. `tests/test_shear_amplitude_option.py`
checks that every critical-plane method and SIH follow the selection, that a selection does not
leak into other instances, and that in-phase values are unchanged while out-of-phase values move.

