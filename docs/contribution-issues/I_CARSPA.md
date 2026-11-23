# I_CARSPA: Carpinteri–Spagnoli

**Prepared Issue body; not an open GitHub Issue.** Before posting, search open and closed Issues for the permanent ID and aliases. Do not create a duplicate. Suggested labels: `method:critical-plane`, `status:partial`, `help-wanted`, `validation-needed`. No `good-first-issue` claim is made.

## Scope and scientific reference

- Primary family: Critical plane; secondary: Critical plane.
- Status: Fully implemented; enabled: True.
- Authors/year: A CARPINTERI, A SPAGNOLI (2001).
- Reference: Multiaxial high-cycle fatigue criterion for hard metals, International Journal of Fatigue.
- DOI: https://doi.org/10.1016/s0142-1123(00)00075-x
- Reference evidence: Exact-title Crossref work record; implementation attribution subject to audit caveats; saved metadata in docs/research/doi-records.json.
- Mathematical evidence: Checked against the PDFs of Carpinteri & Spagnoli (2001) Eqs. (22)-(30) and Tables 1-4, Carpinteri et al. (1999) Eqs. (1)-(9), Fig. 3 and Table 1 (Euler-angle reduction), and Carpinteri, Brighenti & Spagnoli (2000) Fig. 3 and cases (a1)-(c2) (notation) (2026-11-23).
- Principle: Weighted mean principal directions from the reduced Euler angles (Eqs. 22-23, W = (sigma_1/sigma_af)^m_sigma for sigma_1 >= c sigma_af); critical-plane normal w at delta = 67.5 [1 - (tau_af/sigma_af)^2] deg from 1^ in the plane 1^3^ (Eq. 24); sqrt(N_max^2 + (sigma_af/tau_af)^2 C_a^2) <= sigma_af, N_max = N_m + N_a, C_a the major semi-axis of the shear ellipse (Eqs. 25-27)..
- Loading scope: Single-harmonic CA, proportional and out of phase, mean stresses; N_m < -N_a on the critical plane out of scope (source); hard metals 1/sqrt(3) <= tau_af/sigma_af <= 1 (MaterialScopeError outside)..

## Available components and remaining work

FABER phase: Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply.

External source discovery: possible_or_related_source; [EXT06](../EXTERNAL_IMPLEMENTATIONS.md#ext06); not integrated or independently validated. Consult the linked inventory for exact variant, language, license, evidence and reuse restrictions. External code does not count as integrated WG4.3 implementation.

Unexplained eta_cal = 0 of three beta = 90 deg rows of Tables 2-4; the later ductile-metal extensions (Carpinteri et al. 2003 onwards) are separate methods.

Existing implementation: [src/multifat/models/cp_methods/car_spa.py](../../src/multifat/models/cp_methods/car_spa.py); symbol `CarSpa`. Existing tests: [tests/methods/test_carspa.py](../../tests/methods/test_carspa.py).

Reuse: Physical harmonic tensors; direct FatigueModel interface; original-law-specific postprocessing. Tensor rotations and candidate plane search where the source plane rule permits.

Expected inputs: Sigm1, Taum1; m_sigma (or b_fat) when the principal directions rotate within the weighted window; c_CARSPA = 0.5 by default. Expected output: Equivalent stress amplitude sigma_eq (Eq. 27) and error index; fitted life through the shared master curve.

Known restrictions: Hard metals only (source scope); mirror-and-return orientation of the fixed frame (repository choice); m_sigma from -1/b_fat when not supplied; right-limit directions at equal principal stresses; peak-instant directions when sigma_1 < c sigma_af everywhere.

## Acceptance criteria

- [ ] Identify exact original formulation, version, supplements and errata; record equations and departures. If source identity/equivalence is unresolved, resolve it before implementation.
- [ ] Agree bounded CA scope: proportional, phase-shifted NP, mean stresses and applicable material domain. Reject unsupported inputs explicitly.
- [ ] Implement missing components faithfully using physical simultaneous tensors and the original amplitude/plane rule; declare material fields, units, calibration and missing-data behavior.
- [ ] Integrate using [the existing interfaces](../../docs/adding-methods.md), registry and runtime source metadata. Do not enable an archived proxy before repair.
- [ ] Check axial/torsion limits, zero/degenerate cases, means, phase/time-origin and frame invariance as applicable; demonstrate grid/quadrature/optimizer convergence where needed.
- [ ] Obtain at least one independently calculated original-paper numerical case with material/loading data and intermediate quantities. Keep calibration and validation separate; record inaccessible data rather than inventing it.
- [ ] Preserve BP/CAIM, disabled guards, reporting and all-slope `log10(DP_e/DP)` behavior; run `MPLBACKEND=Agg python -m unittest discover -s tests -v`.
- [ ] Update `docs/catalogue/methods.json`, attach Issue/PR and validation evidence, regenerate documentation and pass `python tools/generate_method_docs.py --check`.
- [ ] Review implementation completeness separately from validation; never promote validation merely because numerical tests pass.

Validation lead: tests/methods/test_carspa.py: Tables 2-4 (30 Nishihara-Kawamoto limits): error index I from the printed eta_cal for all 30 rows; eta_cal within 1.5 deg for 27 rows, 3 beta = 90 deg rows with tau_a/sigma_a = 0.2 strict xfail (table prints 0, the procedure gives 2.3-3.5 deg); end-to-end I for the 27 rows; cases (a1)-(c2) of the 2000 paper; mirror invariance; closed forms at tau_af/sigma_af = 1/sqrt(3) and 1.

Dependencies: confirm source/version and measured parameters first; consult [ROADMAP.md](../../ROADMAP.md) for shared cycle, quadrature, native-life and data prerequisites. G01/G02 share a paper but retain distinct laws; other candidate dependencies require specification, not automatic coupling.

Use separate reviewable PRs for source/benchmark specification, prerequisites and numerical implementation when the scope warrants it. Repository entry: [METHODS.md#i_carspa](../../METHODS.md#i_carspa).
