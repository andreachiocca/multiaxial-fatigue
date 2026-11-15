# Public external implementations and code leads

Generated from the external-source records in [catalogue/methods.json](catalogue/methods.json). Search/check date: 15 October 2026. Availability and limited static inspection are not equation verification, successful execution, predictive validation or permission to reuse. No third-party code was copied, translated, executed or integrated.

## Search coverage and meaning

| Counted catalogue method/scoped records searched | Confirmed corresponding/bounded-variant source | Possible or related source | No accessible source identified |
| --- | --- | --- | --- |
| 119 | 6 | 10 | 103 |

The 117 records are not a definitive count of independent theories. Two additional future workstreams, C04/ECP and C09/random-spectrum HCF, have confirmed author-associated source links and are reported separately, not added to the 117-method denominator. Other noncounted leads were not individually searched.

Confirmed means an accessible file and a traceable, explicitly bounded method correspondence were inspected. Possible includes related variants, supporting helpers, a command interface without numerical-kernel inspection, documented software without public source, and an inaccessible archive lead. One record can have several sources; confirmed takes precedence for the per-method count. No source found is a limited search outcome, not proof of absence. Source-code availability never changes WG4.3 implementation or validation status.

## License and scope gate

Six source projects offer potentially useful licensed current-scope implementations or supporting helpers: FS-Sol, FI-Sol, openfat, FLife, Fatlab and fatpack. ECP is a seventh licensed project for future geometry/gradient work. These are project counts, not numbers of independently validated laws. GPLv3/GPLv2-or-later, MIT and ISC have different obligations. This repository is licensed MIT: MIT, BSD, ISC and Apache-2.0 code can be combined with it, provided notices are kept and the exact file scope is recorded; GPL code cannot be incorporated and may serve only as a reference for an independent implementation. No integration is pre-approved. Unknown-license ASU repositories are reference-only pending permission. Public visibility or a copyright/citation notice is not a reuse grant.

The current FABER phase is CA fatigue of plain/smooth geometries without explicit notches, cracks or significant gradients. Existing FE wrappers can sometimes supply point kernels, but that does not import geometry processing into this scope. FLife’s COIN-LiWI kernel uses complex amplitude spectra; its spectral Carpinteri–Spagnoli variant is not a replacement for the original CA plane construction. ASU’s 2006 strain-based and 2020 energy sources must not be substituted for differently cited catalogue laws.

## Source inventory

### EXT01

**FS-Sol** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_FS](../METHODS.md#i_fs) |
| Project | https://github.com/achiocca1/FS-Sol |
| Direct source / lead | https://github.com/achiocca1/FS-Sol/blob/bbab33d9d6d9979c8ea21e5543c4cf023818a6f8/FSSol_v1/FSSol_v1.m |
| Checked revision | bbab33d9d6d9979c8ea21e5543c4cf023818a6f8 |
| Language / dependencies | MATLAB; MATLAB; scripts read ANSYS-exported files and use user-set material constants. |
| License / evidence | GPL-3.0; license file inspected |
| Formulation / loading / geometry | Bounded extended FS, isotropic linear elasticity and proportional loading; not the original maximum-shear plane search or arbitrary NP implementation. |
| Tests / scientific evidence | Author-associated repository links 10.1111/ffe.14153; example/validation data advertised; script inspected, not executed. |
| Potential use and limitations | Adapt or benchmark the proportional elastic _ext branch after GPL-compatible licensing and independent verification; not a drop-in NP solver. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT02

**FI-Sol** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_FIN](../METHODS.md#i_fin) |
| Project | https://github.com/achiocca1/FI-Sol |
| Direct source / lead | https://github.com/achiocca1/FI-Sol/blob/fd4406dff6092dcb54590296043193301a3aa583/FISol_v1.0/MAIN_FI.m |
| Checked revision | fd4406dff6092dcb54590296043193301a3aa583 |
| Language / dependencies | MATLAB; MATLAB; reads ANSYS-exported stress/coordinate CSV files. |
| License / evidence | GPL-3.0; license file inspected |
| Formulation / loading / geometry | Closed-form Findley evaluation for proportional linear-elastic loading; stress-point calculation applies to smooth geometries although supplied examples include components. |
| Tests / scientific evidence | Author-associated repository links 10.1016/j.euromechsol.2024.105274; supplied examples/validation described; no external run performed. |
| Potential use and limitations | Analytical reference or future acceleration of proportional Findley; respect domain and GPL obligations. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT03

**openfat** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [N01](../METHODS.md#n01), [N02](../METHODS.md#n02) |
| Project | https://github.com/cloestro/openfat |
| Direct source / lead | https://github.com/cloestro/openfat/blob/782bbdf93080358aebf9decea4eb127fd3fe9265/crossland.c |
| Checked revision | 782bbdf93080358aebf9decea4eb127fd3fe9265 |
| Language / dependencies | C; C/math library, matrix.c and utilities.c; archived Google Code mirror. |
| License / evidence | GPL-2.0-or-later; explicit source-file notices inspected |
| Formulation / loading / geometry | Crossland and Sines numerical routines accepting stress histories; exact NP amplitude, normalization and interpolation choices still need comparison with selected originals. |
| Tests / scientific evidence | crossland.c, sines.c and progtest.c read; contains driver/example, not a verified predictive validation suite. TODO comments and historical numerical choices remain. |
| Potential use and limitations | Candidate translation/adaptation or independent comparison after GPL-compatible license decision; not automatically correct. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT04

**ASU time-derivative energy model** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [G35](../METHODS.md#g35) |
| Project | https://github.com/ymlasu/Energy-based-time-derivative-damage-accumulation-model-under-uniaxial-and-multiaxial-random-loadings |
| Direct source / lead | https://github.com/ymlasu/Energy-based-time-derivative-damage-accumulation-model-under-uniaxial-and-multiaxial-random-loadings/blob/6cf6d82e607f15e6cd0bb8045105608cdaf810fe/Fatigue_Model_D.m |
| Checked revision | 6cf6d82e607f15e6cd0bb8045105608cdaf810fe |
| Language / dependencies | MATLAB; MATLAB; supplied material/path data and case-specific calibration scripts. |
| License / evidence | Not established; repository root and README contain copyright notice but no explicit license grant identified |
| Formulation / loading / geometry | Author-associated implementation for constant and random uniaxial/multiaxial loading; constant_*.m cases are present. A faithful repeating-cycle integration still needs specification. |
| Tests / scientific evidence | README cites original model and arXiv 2106.02128; MATLAB model and constant/random case listings inspected. Hard-coded material ratio in Fatigue_Model_D.m needs review; no execution or equation audit. |
| Potential use and limitations | Reference only pending permission/license clarification; obtain original paper and reproduce a CA example before considering translation. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT05

**FLife COIN-LiWI** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [G08](../METHODS.md#g08) |
| Project | https://github.com/ladisk/FLife |
| Direct source / lead | https://github.com/ladisk/FLife/blob/dfe660987296abd5670d87f16e709065da3315db/FLife/multiaxial/criteria.py |
| Checked revision | dfe660987296abd5670d87f16e709065da3315db |
| Language / dependencies | Python; NumPy/SciPy plus FLife spectral interfaces; wrapper expects 3D amplitude spectrum. |
| License / evidence | MIT; LICENSE inspected |
| Formulation / loading / geometry | _coin_liwi uses complex tensor amplitude spectra; unlike most package methods it is not PSD-only. A harmonic smooth-point adapter may fit current scope, but means, spectrum normalization and wrapper behavior require tests. |
| Tests / scientific evidence | eq_stress.coin_liwi cites Schmidt and Kraft original paper. Package contains tests/multiaxial_test.py and comparison notebooks; not run here and coverage of this exact kernel not established. |
| Potential use and limitations | High-value kernel/reference candidate, not wholesale spectral-framework adoption. Independently verify branch handling and wrapper output before reuse. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT06

**FLife spectral Carpinteri–Spagnoli** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_CARSPA](../METHODS.md#i_carspa) |
| Project | https://github.com/ladisk/FLife |
| Direct source / lead | https://github.com/ladisk/FLife/blob/dfe660987296abd5670d87f16e709065da3315db/FLife/multiaxial/criteria.py |
| Checked revision | dfe660987296abd5670d87f16e709065da3315db |
| Language / dependencies | Python; NumPy/SciPy and multiaxial PSD/spectral data. |
| License / evidence | MIT; LICENSE inspected |
| Formulation / loading / geometry | Accessible _cs is a frequency-domain Carpinteri–Spagnoli variant, not established equivalent to the original CA plane-orientation procedure required by WG4.3. |
| Tests / scientific evidence | Source names and package references inspected; package tests/notebooks exist, but exact CA original correspondence not verified. |
| Potential use and limitations | Related algorithm reference only; do not use it to re-enable the disabled WG4.3 original formulation. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT07

**Code_Aster fatigue commands** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_FS](../METHODS.md#i_fs), [I_MATAKE](../METHODS.md#i_matake), [I_DANGVAN](../METHODS.md#i_dangvan), [N01](../METHODS.md#n01), [N03](../METHODS.md#n03) |
| Project | https://github.com/ralic/Code_Aster |
| Direct source / lead | https://github.com/ralic/Code_Aster/blob/master/code_aster/Cata/Commands/post_fatigue.py |
| Checked revision | Page/branch checked on stated date; not a pinned numerical release |
| Language / dependencies | Python / Fortran ecosystem; Full Code_Aster FE/runtime stack; numerical kernels not reviewed in this pass. |
| License / evidence | GPL-3.0 (repository metadata); GitHub license metadata; inspected command file, not underlying numerical routines |
| Formulation / loading / geometry | Command schema names modified periodic Matake/Dang Van, Crossland/Papadopoulos and nonperiodic FS variants. Numerical source correspondence and exact versions still require inspection. |
| Tests / scientific evidence | Command source accessible in a historical mirror; current official GitLab page exposed only a loading shell in this environment. R7.04.04 documents periodic/VA formulations and amplitude conventions. |
| Potential use and limitations | Strong reference/documentation lead; not a verified standalone criterion source and not ready for copying/translation. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT08

**Fatlab equivalent-stress helpers** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [N38](../METHODS.md#n38), [N39](../METHODS.md#n39), [N40](../METHODS.md#n40) |
| Project | https://sourceforge.net/p/fatlab/code/HEAD/tree/trunk/calc_fatigue_stress.m |
| Direct source / lead | https://sourceforge.net/p/fatlab/code/HEAD/tree/trunk/calc_fatigue_stress.m |
| Checked revision | SVN r3 |
| Language / dependencies | MATLAB; MATLAB; full GUI/FE environment is not needed conceptually for point helpers but dependencies require isolation. |
| License / evidence | GPL-3.0 (project declaration); SourceForge project licensing declaration |
| Formulation / loading / geometry | calc_fatigue_stress.m includes von Mises and principal/Tresca stress-history helpers; these are not automatically correct fatigue-amplitude reductions for reversed or NP cycles. |
| Tests / scientific evidence | SVN r3 source inspected; benchmarks directory exists. In particular unsigned equivalent-history half-ranges must not replace tensor-amplitude definitions. |
| Potential use and limitations | Potential licensed helper/reference adaptation only after defining signed/amplitude policy and testing degeneracy; not a complete verified CA baseline evaluator. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT09

**fatpack SWT stress correction** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_SWT](../METHODS.md#i_swt) |
| Project | https://github.com/Gunnstein/fatpack |
| Direct source / lead | https://github.com/Gunnstein/fatpack/blob/master/fatpack/stresscorrection.py |
| Checked revision | Page/branch checked on stated date; not a pinned numerical release |
| Language / dependencies | Python; Python/NumPy; stress-range/mean conventions must match. |
| License / evidence | ISC; repository license metadata; stress-correction source inspected |
| Formulation / loading / geometry | find_swt_equivalent_stress is a uniaxial stress-life mean-stress correction; it is not the repository multiaxial sigma_max*epsilon_a plane search. |
| Tests / scientific evidence | Source cites mean-stress literature; package helper not executed or mathematically validated here. |
| Potential use and limitations | Supporting licensed reference for mean-stress checks, not a complete multiaxial SWT implementation. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT10

**PragTic Wang–Brown** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [N46](../METHODS.md#n46) |
| Project | https://www.pragtic.com/PragTicHelp/files/%7BB1530853-7C95-4C02-B465-FEB43C9456D0%7D.htm |
| Direct source / lead | https://www.pragtic.com/PragTicHelp/files/%7BB1530853-7C95-4C02-B465-FEB43C9456D0%7D.htm |
| Checked revision | Page/branch checked on stated date; not a pinned numerical release |
| Language / dependencies | Not established for inspected distribution; PragTic application; no standalone source/dependencies verified. |
| License / evidence | Not established / access restrictions; official help and project status; historical source access by request is not a public reuse grant |
| Formulation / loading / geometry | Official help describes Wang–Brown 1993/1996 and alternative normal-excursion/plane selections; original version must be selected. |
| Tests / scientific evidence | Public method documentation inspected; source code not publicly retrieved. Freeware executable availability does not imply open-source reuse. |
| Potential use and limitations | Documented software/reference lead only. Ask rights holders if source access or reuse is sought; no automatic request sent. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT11

**ASU strain-based characteristic-plane model** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [N16](../METHODS.md#n16) |
| Project | https://github.com/ymlasu/strain-based-model |
| Direct source / lead | https://github.com/ymlasu/strain-based-model/blob/0cb1b31d1a30484f4b7490b40a38f5634f94be3f/mafa304sin.m |
| Checked revision | 0cb1b31d1a30484f4b7490b40a38f5634f94be3f |
| Language / dependencies | MATLAB; MATLAB; cyclic material and strain data; sinusoidal and stepped examples. |
| License / evidence | Not established; README copyright notice, no explicit root license identified |
| Formulation / loading / geometry | Code is for the separate 2006 strain-based Liu–Mahadevan publication (10.1111/j.1460-2695.2005.00957.x), not proven equivalent to N16's 2007 stress-based formulation. |
| Tests / scientific evidence | Author-associated README and sinusoidal CA strain/phase script inspected; one material supplied; related source only. |
| Potential use and limitations | Reference only pending permission and exact variant mapping. Do not label N16 implemented from this differently cited source. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT12

**Author ECP implementation** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [C04](../METHODS.md#c04) |
| Project | https://github.com/achiocca1/ECP |
| Direct source / lead | https://github.com/achiocca1/ECP/blob/9d5115258d5a13b926788a0c636a417d94e85efd/MAIN.m |
| Checked revision | 9d5115258d5a13b926788a0c636a417d94e85efd |
| Language / dependencies | MATLAB; MATLAB and spatial stress/strain/material-radius data; examples involve notches. |
| License / evidence | GPL-3.0; LICENSE inspected |
| Formulation / loading / geometry | Nonlocal volume-averaged effective critical plane for geometry/stress concentrations; outside present plain/no-gradient phase. |
| Tests / scientific evidence | Author repository links 10.1016/j.ijfatigue.2024.108565 and 108716; source tree, examples and result files accessible. README advertises some newer paths absent from inspected tree; inspect revision before relying on those paths. |
| Potential use and limitations | Future geometry/gradient reference or GPL-compliant adaptation, not a current smooth-point drop-in. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT13

**ASU random-spectrum HCF model** — confirmed

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [C09](../METHODS.md#c09) |
| Project | https://github.com/ymlasu/Multiaxial-high-cycle-fatigue-life-prediction-model |
| Direct source / lead | https://github.com/ymlasu/Multiaxial-high-cycle-fatigue-life-prediction-model/blob/96e5768cbb8d071c307b27f2d1b78666da086ab3/stress_based_prediction_constant.m |
| Checked revision | 96e5768cbb8d071c307b27f2d1b78666da086ab3 |
| Language / dependencies | MATLAB; MATLAB and constant_loading.mat/material calibration; random workflow beyond current scope. |
| License / evidence | Not established; README copyright notice; no explicit license grant identified |
| Formulation / loading / geometry | Author code for the specific Wei et al. random-spectrum publication; includes an explicit constant-loading driver. Distinct from automatically claiming the base Liu–Mahadevan original is implemented. |
| Tests / scientific evidence | README/source agree on 10.1016/j.ijfatigue.2019.105462; source material constants and phase/path call inspected; tests not run. |
| Potential use and limitations | Reference CA driver and data pending permission; fuller random-spectrum implementation is a future extension. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT14

**ASU biaxial energy model** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | Unmapped related paper; no new criterion promoted |
| Project | https://github.com/ymlasu/An-energy-based-model-under-multiaxial-biaxial-loading |
| Direct source / lead | https://github.com/ymlasu/An-energy-based-model-under-multiaxial-biaxial-loading/blob/f8b9f47125df2406e41f075632abde86cd837d0f/Fatigue_Model.m |
| Checked revision | f8b9f47125df2406e41f075632abde86cd837d0f |
| Language / dependencies | MATLAB; MATLAB and supplied loading/material workbooks. |
| License / evidence | Not established; README copyright notice; no explicit license grant identified |
| Formulation / loading / geometry | 2020 energy paper 10.1016/j.ijfatigue.2020.105858; not matched to an existing exact catalogue criterion. In particular this is not G52's 2017 critical-plane energy paper. |
| Tests / scientific evidence | Author README/tree and material example listings inspected; primary-method correspondence remains unassigned. |
| Potential use and limitations | Retain as unmapped reference lead, not a new independent criterion or a substitute for G52. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

### EXT15

**Zenodo modified Findley software** — possible

| Field | Evidence / restriction |
| --- | --- |
| Catalogue correspondence | [I_FIN](../METHODS.md#i_fin) |
| Project | https://zenodo.org/records/10370862 |
| Direct source / lead | https://zenodo.org/records/10370862 |
| Checked revision | Page/branch checked on stated date; not a pinned numerical release |
| Language / dependencies | Unverified; Unverified. |
| License / evidence | Unverified; search record only; record/file listing could not be opened |
| Formulation / loading / geometry | FKM-based modified Findley, not the unrestricted original criterion. |
| Tests / scientific evidence | Record 10370862 located by search; direct fetch inaccessible. Files, license and numerical implementation not verified. |
| Potential use and limitations | Unverified lead only; no reuse until source and license can be inspected. |
| Verification boundary | Availability/identity/limited static inspection only; no third-party code executed, ported or integrated. |

## Method-by-method discovery register

This covers every counted record plus the two source-linked future workstreams. Follow the catalogue for original reference, parameters, WG4.3 implementation and validation.

| Method ID / name | Search disposition | External links | Current phase |
| --- | --- | --- | --- |
| [I_FS](../METHODS.md#i_fs) — Fatemi–Socie (repository elastic/harmonic scope) | confirmed_source_available | confirmed_source_available; [EXT01](EXTERNAL_IMPLEMENTATIONS.md#ext01), [EXT07](EXTERNAL_IMPLEMENTATIONS.md#ext07); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_FIN](../METHODS.md#i_fin) — Findley | confirmed_source_available | confirmed_source_available; [EXT02](EXTERNAL_IMPLEMENTATIONS.md#ext02), [EXT15](EXTERNAL_IMPLEMENTATIONS.md#ext15); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_SWT](../METHODS.md#i_swt) — Smith–Watson–Topper (repository plane extension) | possible_or_related_source | possible_or_related_source; [EXT09](EXTERNAL_IMPLEMENTATIONS.md#ext09); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MATAKE](../METHODS.md#i_matake) — Matake (tension-only repository adaptation) | possible_or_related_source | possible_or_related_source; [EXT07](EXTERNAL_IMPLEMENTATIONS.md#ext07); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MCD](../METHODS.md#i_mcd) — McDiarmid (surface-crack repository scope) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MGSE_YU](../METHODS.md#i_mgse_yu) — Yu generalized strain energy (repository damage-plane adaptation) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MGSE_ZHU](../METHODS.md#i_mgse_zhu) — Zhu weighted generalized strain energy | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MSWT](../METHODS.md#i_mswt) — Scalar modified SWT reported by Zhu (2018) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_MKBM](../METHODS.md#i_mkbm) — Li modified Kandil–Brown–Miller | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_LIU2021](../METHODS.md#i_liu2021) — Liu et al. (2021) modified shear amplitude | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_LI](../METHODS.md#i_li) — Li et al. (2021) improved strain energy density | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_GSE](../METHODS.md#i_gse) — Ince–Glinka generalized strain energy | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_GSA](../METHODS.md#i_gsa) — Ince–Glinka generalized strain amplitude | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_BP](../METHODS.md#i_bp) — Böhme–Papuga | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_CAIM](../METHODS.md#i_caim) — Critical–integral approach (CAIM) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_CARSPA](../METHODS.md#i_carspa) — Carpinteri–Spagnoli | possible_or_related_source | possible_or_related_source; [EXT06](EXTERNAL_IMPLEMENTATIONS.md#ext06); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_DANGVAN](../METHODS.md#i_dangvan) — Dang Van | possible_or_related_source | possible_or_related_source; [EXT07](EXTERNAL_IMPLEMENTATIONS.md#ext07); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_OTT](../METHODS.md#i_ott) — Ottosen/Tveit proportional reduction (archived scope) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_SWTD](../METHODS.md#i_swtd) — Deviatoric SWT | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [I_ZHU_EDP](../METHODS.md#i_zhu_edp) — Zhu equivalent damage parameter | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N01](../METHODS.md#n01) — Crossland | confirmed_source_available | confirmed_source_available; [EXT03](EXTERNAL_IMPLEMENTATIONS.md#ext03), [EXT07](EXTERNAL_IMPLEMENTATIONS.md#ext07); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N02](../METHODS.md#n02) — Sines | confirmed_source_available | confirmed_source_available; [EXT03](EXTERNAL_IMPLEMENTATIONS.md#ext03); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N03](../METHODS.md#n03) — Papadopoulos, specified 2001 formulation | possible_or_related_source | possible_or_related_source; [EXT07](EXTERNAL_IMPLEMENTATIONS.md#ext07); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N04](../METHODS.md#n04) — Modified Wöhler Curve Method / Susmel–Lazzarin | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N05](../METHODS.md#n05) — Wang–Sun–Liu–Liu–Shang maximum-shear-plane HCF criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N06](../METHODS.md#n06) — Zhang–Yao nonlinear critical-plane HCF criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N07](../METHODS.md#n07) — Sectional critical-plane criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N08](../METHODS.md#n08) — Gonçalves–Araújo–Mamiya prismatic-hull criterion for hard metals | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N09](../METHODS.md#n09) — Unified equivalent shear stress | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N10](../METHODS.md#n10) — Modified Marin-type mean-shear criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N11](../METHODS.md#n11) — Ninić Damage Function (quadratic damage-plane criterion) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N12](../METHODS.md#n12) — Revised multiaxial fatigue damage function | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N13](../METHODS.md#n13) — Golos–Dębski–Dębski modified Gough limit curve | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N15](../METHODS.md#n15) — Original Liu–Zenner integral criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N16](../METHODS.md#n16) — Liu–Mahadevan critical-plane criterion and life model | possible_or_related_source | possible_or_related_source; [EXT11](EXTERNAL_IMPLEMENTATIONS.md#ext11); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N17](../METHODS.md#n17) — Mean-shear-sensitive weighted-plane criterion, CA branch | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N18](../METHODS.md#n18) — Piecewise ruled S–N surface | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N19](../METHODS.md#n19) — Improved Gough / Nishihara–Kawamoto criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N20](../METHODS.md#n20) — Stress-invariant endurance criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N21](../METHODS.md#n21) — Equivalent Second Moment of Load Path (ESMLP) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N22](../METHODS.md#n22) — Moment of Load Path (MLP), closed CA stress-path branch | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N23](../METHODS.md#n23) — Equal-life surface with path correction | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N24](../METHODS.md#n24) — Stress Scale Factor (SSF) equivalent shear stress | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N25](../METHODS.md#n25) — A90(T) deviatoric tensile critical-plane model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N26](../METHODS.md#n26) — Kujawski negative-mean-stress modified SWT (mSWT) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N27](../METHODS.md#n27) — Wu–Hu–Song (WHS) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N28](../METHODS.md#n28) — Modified Wu–Hu–Song (MWHS) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N29](../METHODS.md#n29) — Original Kandil–Brown–Miller (KBM) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N30](../METHODS.md#n30) — Normalized energy critical-plane parameter | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N31](../METHODS.md#n31) — Load-path-sensitive modified fatigue parameter | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N32](../METHODS.md#n32) — Multiaxial Fatigue Space (MFS) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N34](../METHODS.md#n34) — Probabilistic slip-site model, homogeneous prescribed-volume scope | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N35](../METHODS.md#n35) — (V)HCF-DAM, prescribed proportional blocks | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Full VA/random/block workflow is beyond present CA scope; any faithful repeating-cycle subcase needs separate specification. |
| [N36](../METHODS.md#n36) — Tanaka–Mura–Wu HCF crack-nucleation model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N37](../METHODS.md#n37) — Projection-by-Projection (PbP), harmonic CA scope | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N38](../METHODS.md#n38) — von Mises equivalent-amplitude baseline | possible_or_related_source | possible_or_related_source; [EXT08](EXTERNAL_IMPLEMENTATIONS.md#ext08); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N39](../METHODS.md#n39) — Tresca equivalent-amplitude baseline | possible_or_related_source | possible_or_related_source; [EXT08](EXTERNAL_IMPLEMENTATIONS.md#ext08); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N40](../METHODS.md#n40) — Maximum principal stress / Rankine baseline | possible_or_related_source | possible_or_related_source; [EXT08](EXTERNAL_IMPLEMENTATIONS.md#ext08); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N41](../METHODS.md#n41) — Elastic strain-energy equivalent-stress baseline | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N42](../METHODS.md#n42) — IIW normal/shear interaction assessment | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [N43](../METHODS.md#n43) — Eurocode 3 normal/shear fatigue interaction | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [N44](../METHODS.md#n44) — Effective Equivalent Stress Hypothesis (EESH) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [N46](../METHODS.md#n46) — Wang–Brown strain-based assessment, selected original version | possible_or_related_source | possible_or_related_source; [EXT10](EXTERNAL_IMPLEMENTATIONS.md#ext10); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N47](../METHODS.md#n47) — Original Marin mean-stress criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N48](../METHODS.md#n48) — Froustey–Lasserre energetic criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N50](../METHODS.md#n50) — Grubisic–Simbürger formulation | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N52](../METHODS.md#n52) — Itoh NP strain correction / criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N53](../METHODS.md#n53) — Substitute-amplitude/NP-damage method, CA branch | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N55](../METHODS.md#n55) — Gradient-corrected local evaluator with externally supplied gradients | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [N57](../METHODS.md#n57) — Combined deviatoric strain amplitudes | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [C04](../METHODS.md#c04) — Chiocca–Frendo ECP workstream | confirmed_source_available | confirmed_source_available; [EXT12](EXTERNAL_IMPLEMENTATIONS.md#ext12); not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [C09](../METHODS.md#c09) — Wei random-spectrum extension | confirmed_source_available | confirmed_source_available; [EXT13](EXTERNAL_IMPLEMENTATIONS.md#ext13); not integrated or independently validated | Full VA/random/block workflow is beyond present CA scope; any faithful repeating-cycle subcase needs separate specification. |
| [G01](../METHODS.md#g01) — Papuga–Růžička critical-plane criterion (PCR; later PCRN) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G02](../METHODS.md#g02) — Papuga–Růžička integral criterion (PI; later PIN) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G03](../METHODS.md#g03) — Papuga QCP mean-stress development | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G04](../METHODS.md#g04) — Kakuno–Kawada | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G05](../METHODS.md#g05) — Manson–McKnight–Papuga (MMP) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G06](../METHODS.md#g06) — Margetin–Ďurka–Chmelko (MDC) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G07](../METHODS.md#g07) — Modified Mohr–Mises hypothesis (MMMH) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G08](../METHODS.md#g08) — Complex-invariant linear wave interference (COIN LiWI) | confirmed_source_available | confirmed_source_available; [EXT05](EXTERNAL_IMPLEMENTATIONS.md#ext05); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G09](../METHODS.md#g09) — Moment-adapted COIN LiWI (MA-COIN LiWI) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G10](../METHODS.md#g10) — Abasolo–Pallares-Santasmartas–Eizmendi exponent criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G11](../METHODS.md#g11) — Yankin–Mugatarov–Wildemann loading-path criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G12](../METHODS.md#g12) — Twin-shear unified-yield fatigue model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G13](../METHODS.md#g13) — Braccesi–Cianetti–Lori–Pioli equivalent uniaxial process | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G15](../METHODS.md#g15) — Liu virtual strain energy (VSE modes I/II) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G16](../METHODS.md#g16) — Ellyin–Kujawski multiaxial energy criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G17](../METHODS.md#g17) — Garud plastic-work approach | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G18](../METHODS.md#g18) — Łagoda–Macha–Będkowski critical-plane energy | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G20](../METHODS.md#g20) — Lu–Wu–Zhong constant-amplitude energy model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G21](../METHODS.md#g21) — Lu–Melendez–Martínez-Esnaola energy/plane family | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G22](../METHODS.md#g22) — Modified Manson–Coffin Curve Method (MMCCM) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G23](../METHODS.md#g23) — Morel mesoscopic constant-amplitude life model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G24](../METHODS.md#g24) — Flacelière–Morel probabilistic fatigue family | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Full spatial surface/volume framework is future work; a homogeneous plain-volume restriction requires original-source review. |
| [G25](../METHODS.md#g25) — Morel–Huyen heterogeneous shakedown model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G26](../METHODS.md#g26) — Palin-Luc–Lasserre / Banvillet volumetric energy | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G27](../METHODS.md#g27) — Morel–Palin-Luc nonlocal HCF theory | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G28](../METHODS.md#g28) — Delahay–Palin-Luc probabilistic gradient-energy model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G29](../METHODS.md#g29) — Papadopoulos–Panoskaltsis invariant-gradient criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G30](../METHODS.md#g30) — Berto–Lazzarin averaged strain energy density (SED) | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G31](../METHODS.md#g31) — Wormsen–Härkegård weakest-link / extremes model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Future geometry/gradient/structural extension; not part of the present plain/smooth, no-gradient phase. |
| [G32](../METHODS.md#g32) — Stefanov curvilinear integral / IDD family | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G33](../METHODS.md#g33) — Zhai–Zou–Hu–Chen–Jin shakedown-constrained Walker model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G35](../METHODS.md#g35) — Tien–Wei–Chen–Liu energy time-derivative damage | confirmed_source_available | confirmed_source_available; [EXT04](EXTERNAL_IMPLEMENTATIONS.md#ext04); not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G36](../METHODS.md#g36) — Chaudonneret multiaxial continuum damage | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G37](../METHODS.md#g37) — Lemaitre–Sermage–Desmorat two-scale damage | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G38](../METHODS.md#g38) — Monchiet–Charkaluk–Kondo PSB damage criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G39](../METHODS.md#g39) — Pessard–Morel–Morel–Bellett competing-mechanism model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G40](../METHODS.md#g40) — Shang–Wang critical-plane damage model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G41](../METHODS.md#g41) — Robert–Fogue–Bahuaud integral life criterion | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G42](../METHODS.md#g42) — Malnati mesoplastic cycle-indicator method | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G43](../METHODS.md#g43) — Xu et al. mean-stress/biaxiality modification | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G44](../METHODS.md#g44) — Feng et al. elevated-temperature octahedral-stress model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G45](../METHODS.md#g45) — Shao et al. compression–torsion bearing-steel model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G46](../METHODS.md#g46) — Probabilistic FS/Chaboche material-variability framework | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G47](../METHODS.md#g47) — Gates–Fatemi revised FS normalization | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G48](../METHODS.md#g48) — Sharifimehr–Fatemi normal/shear interaction extension | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G49](../METHODS.md#g49) — Leis total octahedral-energy parameter | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [G52](../METHODS.md#g52) — Wei–Liu critical-plane energy model | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N58](../METHODS.md#n58) — Liu–Xie–Hao equivalent stress on the maximum shear stress amplitude plane | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |
| [N59](../METHODS.md#n59) — Albinmousa–Al Hussain polar damage sum, Smith–Watson–Topper parameter | no_accessible_source_identified | no_accessible_source_identified; not integrated or independently validated | Plain/smooth CA candidate or scoped point formulation; original loading, material and response requirements still apply. |

## Search method and unresolved leads

The [query register](research/external-code-search.json) records one named source-code query per counted record, with targeted stronger searches and citation/repository chaining for useful leads. GitHub, GitLab, MathWorks File Exchange, SourceForge, author/research sites, Zenodo and publication-linked software were screened. Useful links were checked against primary repository metadata, actual files, documentation and license notices. Generic matches and homonyms were rejected. This is a targeted first pass, not exhaustive code indexing or a saturation claim.

PragTic provides valuable author documentation but public solver source/reuse rights were not established. A historical source-on-request statement is not an open-source license. The inspected FatPy repository was empty: an implementation Issue is not source code. MathWorks LoadPathMATLAB visualizes FE load flow, not the fatigue Moment of Load Path law; sports/physiological fatigue repositories were rejected. Zenodo modified-Findley archive files could not be opened and remain unverified. The current official Code_Aster GitLab page returned a loading shell; a historical GitHub mirror exposed the command interface, not an audited standalone numerical kernel.

Recommendations: first compare proportional FS/Findley against the author MATLAB algorithms; specify source-faithful Crossland/Sines with openfat as a checked-but-unvalidated comparator; investigate the MIT COIN-LiWI kernel and its wrapper independently; ask for permission before porting ASU G35 or its data. Preserve ECP for a future geometry extension. Do not change method validation/completeness solely from this inventory. No live method Issues exist; relevant prepared briefs now point here, and actual Issues should inherit these references when posted.
