# Provenance of the experimental workbooks (draft)

Sources, row counts and known issues of the workbooks. The active folder holds the revised
October 2026 collection (upstream 71eaac3, numbered files) plus two batch workbooks and the edits listed in
[October 2026 collection](#october-2026-collection-upstream-71eaac3); the sections after it describe the
previous four workbooks (now in `archive/pre-dataset-migration/` as upstream archived them, and with our
2026-10-18 edits in commit 23e1dfe). Counts are for the `Dataset`
sheet. "Plain" means `Notched = No`; "failures" means `Run out = No`. Page numbers are physical PDF pages of
the versions of record. Rows are **digitised** when they were read off a figure: their cycle counts are
non-integer, and this holds for every row in such a group.

Open decisions are listed at the end of each section; every changed cell carries an Excel comment.

**Overlay layer.** Five workbooks carry corrections that are not in their supplied files: `2_CX18`,
`3_25CrMo4`, `4_30NCD16`, `7_982FA` and `8_AISI316L`. Those files are kept byte-identical to the
version their source supplied (`upstream` 71eaac3 and later). Each correction is a machine-readable
operation in `src/multifat/imp/overlay_files/<stem>.overlay.json`, with its decision and source, the value it replaces
(checked on every application) and the SHA-256 of the file it applies to. The importer applies it in
memory when the file is read, so the data are those of the revised workbooks; `--no-overlay` gives
the supplied data. The cell-level changes below are the human-readable record of the same operations.
The equivalence (as supplied plus overlay equals the revised workbook, same dataset fingerprint) is
tested in `tests/test_dataset_overlays.py`.

## October 2026 collection (upstream 71eaac3)

Merged on 2026-10-19 (decision d054). Twelve revised workbooks (`1_` to `13_`, upstream 71eaac3) plus the retained
`Al7075_T6.xlsx`; the material ID is the file name without its number. The revised Summary has no fitted
constants (they are calibrated by the pipeline, or read from `Calibrations/` with `--constants`).

What the revised files contain compared with the previous workbooks (`tools/audit_dataset_migration.py`):

| Material | Revised file | Change from the previous workbook |
|---|---|---|
| 42CrMo4_Santus | `6_42CrMo4_Santus.xlsx` | Santus and Benedetti rows only (45 plain failures, 11 run-outs); the Anes/Reis rows (Su 1100 batch) are not in it |
| Ti6Al4V | `1_Ti6Al4V_Meneghetti.xlsx` | Meneghetti 10 failures and 4 Delahay run-outs (removed by overlay, d109); the 48 Berto failures are not in it |
| AISI316L | `8_AISI316L.xlsx` | Janssens and Braun rows as originally supplied (without N6, N8 and our E, nu) |
| Al7075_T6 | `Al7075_T6.xlsx` | measured data unchanged; Summary coefficient cells removed |
| eight new materials | `2_`-`5_`, `7_`, `10_`-`13_` | En3B and GJS600 notched only; CX18, 30NCD16, 25CrMo4, 982FA, 34Cr4, GGG60 and S65A have only fatigue-limit rows (`Run out = Yes` at the endurance basis) |

### Edits of 2026-10-19 (Excel COM, comment on every changed cell)

| Workbook | Change | Source | Decision |
|---|---|---|---|
| `14_42CrMo4_Reis.xlsx` (new) | the 29 Anes/Reis rows of the previous 42CrMo4_Santus.xlsx (25 failures, 4 run-outs, including N1), in the revised layout; Summary rows Anes IJF 58 and FFEMS 37 with E 206000, nu 0.3, Su 1100, Sy 980 and the strain-life constants sigma_f 1154 MPa, b -0.061, epsilon_f 0.180, c -0.53 (new Summary columns) | Anes et al., IJF 58 (2014) Table 1 p. 2 | d055 |
| `15_Ti6Al4V_Berto.xlsx` (new) | the 48 Berto rows of the previous Ti6Al4V_andrea_(AXIAL).xlsx; Summary E 110000, nu 0.3, Su 978, Sy 894, fatigue limits 475.74 / 388.25 MPa at 2e6 cycles, R = -1 | Berto et al., FFEMS 38 (2015): Table 1 p. 3, text p. 9, Table 3 p. 5 | d055 |
| `8_AISI316L.xlsx` | N6: torsional phase 90 on Non_prop_1-11; N8: the 8 Janssens Table 4 rows (Uniaxial_10-14, Torsional_8-9, Non_prop_12) | Janssens 2018 Table 4 p. 4 | d056 |
| `8_AISI316L.xlsx` | Summary: Janssens E 193000, nu 0.3; Braun E, nu set to Null (were 219000 / 0.29) | Zhu et al. 2019 p. 8; Braun 2024 p. 8 gives 219 GPa / 0.29 as the X-ray elastic constants of the (220) planes for residual stresses, not bulk values | d056 |
| `8_AISI316L.xlsx` | Braun plain rows: Mean AXIAL stress Null -> sigma_a (17 rows) | R = 0 tests, Braun 2024 Fig. 13(i); values of the previous workbook (load ratio 0) | ruling, see below |
| `2_CX18`, `4_30NCD16`, `3_25CrMo4`, `7_982FA` | nu Null -> 0.3 on every Summary source row | Dubbel, 21st ed. (2004), Part C, Hooke's law, PDF p. 123 ("nu_Stahl = 0,30"); not given in the FatLim tables (Boehme & Papuga 2024 Table 1) | d059 |

Rows are copied by header, with the previous workbook's cached values; a legacy `Null` stress or phase is
written as 0 (the legacy convention for an inactive component or an in-phase load). Every copied row was
checked to give the same stress tensors and phases as before. With the Braun mean stresses restored,
AISI316L has 39 plain failures and 13 run-outs, as before the merge.

The Braun mean stresses are not covered by d056: in the revised format an empty or `Null` stress is
unavailable loading and excludes the row (upstream rule, ported), so the 17 Braun plain rows were dropped.
They are R = 0 tests whose mean stress the previous workbook computed from the load ratio; restoring it was
confirmed by T. Grossi on 2026-10-19 (decision d061).

### Open questions on revised workbooks

- 42CrMo4 and Ti6Al4V: confirm whether the Anes/Reis and Berto rows should remain separated as distinct batch workbooks
  (`14_42CrMo4_Reis.xlsx` and `15_Ti6Al4V_Berto.xlsx`) or consolidated into the primary files.
  For 42CrMo4 they stay separate: the three workbooks are the batches compared in issue #71.
- AISI316L: the Braun Summary row gives 387.8 MPa as tensile fatigue limit at R = 0; Braun Fig. 13(i)
  prints it as the stress range Delta sigma_R,50% at 1e6 cycles (amplitude 193.9 MPa). Averaged with
  Janssens' 241 MPa (R = -1) it gives Sigm1 = 314.4 MPa, which the importer does not use (stress-ratio
  check), so models needing Sigm1 skip AISI316L. Its reference cell also holds a pasted Delahay title.
- AISI316L: the Janssens non-proportional rows have torsional phase 0, Table 4 says 90 deg (N6), and the 8
  Table 4 points of N8 are missing.
- CX18: 3 of its 5 rows are identical to 30NCD16 rows (Lasserre & Froustey 1992), and its Summary names
  30NCD16. Same steel?
- CX18: Summary cell M2 (K) refers to another workbook (`[1]!Tabella1`); its cached 1.578 is stale and a
  full recalculation gives `#REF!` (the importer does not read K; the cell was left as supplied).
- 25CrMo4: E = 127000 MPa (Troost 1987 row) is low for a steel; FatLim lists no E.
- Ti6Al4V: the Delahay row (Ti64Al4V, Su 1090, Sy 980) is averaged with Meneghetti (Su 980, Sy 915).
  Resolved 2026-11-11 (decision d108): `src/multifat/imp/overlay_files/1_Ti6Al4V_Meneghetti.overlay.json` sets the
  Delahay row to `Null` and the Average and
  standard-deviation rows to the Meneghetti values: E 110000, nu 0.3, Su 980, Sy 915, Sigm1 500,
  Taum1 422.5 MPa (2e6 cycles, R = -1). `--no-overlay` returns the supplied averages.
- Ti6Al4V: the 4 Delahay rows (Dataset rows 58-61, `Run out = Yes`) are median staircase endurance
  limits of annealed Ti-6Al-4V (Delahay & Palin-Luc 2006, Table 2), not run-outs; 652 and the two
  combined rows are plane bending, not axial, and the torsion limit is 411, not 441 MPa. Resolved
  2026-11-11 (decision d109): the overlay clears the 4 rows (`drop_case`).
- Fatigue-limit datasets: are the `Run out = Yes` rows the fatigue-limit states at the stated cycle basis?
  They are assessed as such (fatigue index error).
- Fatigue-limit datasets: the pure R = -1 rows do not match the Summary limits used as reference, so the
  fatigue index errors mostly measure that mismatch until it is resolved:
  - CX18: rows 2-3 give 310 / 179 MPa, as FatLim's XC18 (Boehme & Papuga 2024 Table 1: 332 / 186 MPa),
    but the Summary holds the 30NCD16 values 710 / 450 MPa.
  - 30NCD16: rows 2-3 give 710 / 450 MPa (Summary average 685 / 430 over the 1992 and 1989 papers); rows
    12-13 give 310 / 179 MPa, the XC18 values again.
  - S65A: rows 2 and 7 give 866 / 541 MPa against a Summary of 583.8 / 370.7 MPa.
  - 34Cr4: rows 2-3 give 410 / 256 MPa; the Summary averages three sources (439.3 / 290 MPa).
  - GGG60: the axial amplitudes (e.g. 549.8 MPa with a static shear) are about twice the Summary 275 MPa;
    are they ranges?

## 16_42CrMo4_FABEST.xlsx (FABEST database, added 2026-11-09)

**Source.** The FABEST database: A benchmark fatigue dataset on 42CrMo4 +QT steel, version
FAB_v20261014 (Zenodo record v3, published 2026-10-14), doi:10.5281/zenodo.23359038, licence
CC BY 4.0. The original file is kept unchanged in `sources/FABEST/FAB_v20261014.xlsm`. The importer
does not read subfolders.
- SHA-256: `720f637f5a1c74f579c62197137db54d0d2a495d6efdc607ddc4548a7d4cefff`
- MD5: `b5ec39bc3df17148015305655e9ebc3c`, as published by Zenodo.

**Conversion.** The workbook is generated by `tools/convert_fabest.py`, which checks the SHA-256
first. `tests/test_fabest_conversion.py` checks that the committed workbook equals the conversion.
Reading rules:
- each test point of a curve sheet is one row;
- `finished` is a failure at "till break"; `interrupted` is a run-out at "till END";
- the static `<curve>_Rm` row is skipped;
- the torsional phase is the channel-2 phase shift relative to the axial channel;
- stresses are as reported, nominal and related to the net section.

**Scope (maintainer decision, 2026-11-09).** One batch only: PBS, plain hollow hourglass specimens
(outer/inner diameter 8.1/7 mm), turned. This series holds every torsion and multiaxial curve of
the database, so the design set and the assessed tests share lab, geometry and finish.

| Curve | Loading | Failures | Run-outs |
|---|---|---|---|
| FAB060, FAB070 | axial, R = -1 (design set) | 20 | 2 |
| FAB062, FAB072 | torsion, R = -1 | 20 | 4 |
| FAB064, FAB074 | axial + torsion, in phase, R = -1 | 18 | 2 |
| FAB065, FAB075 | axial + torsion, 90 deg out of phase, R = -1 | 21 | 2 |
| **Total** | | **79** | **10** |

**Not converted yet.** The other 54 curves of the database (the PBS bar series FAB110-114 are now
`19_42CrMo4_FABEST_bar`): axial curves on other sizes, geometries,
finishes and labs, the mean-stress series (R = 0.1, 0, 0.3, -2, -3, -5), notched curves (V-notches,
transverse hole) and fretting curves. Their axial R = -1 strengths at 1e6 cycles range from 302 to
566 MPa, so pooling them would mix size and surface effects. They are left for later workbooks
(the notched ones for #19).

**Summary (maintainer decision).**
- E = 206362 MPa, nu = 0.296, Su = Rm = 1110.9 MPa, Sy = ReL = 1034.8 MPa: averages of the four PBS
  tensile tests (sheet `Tensile test`).
- Sigm1 = 517.3 MPa and Taum1 = 366.5 MPa: the authors' Basquin strengths at N = 1e6 (`SN_FL`,
  `SN_NFL`) of FAB060/FAB070 and FAB062/FAB072, averaged.
- Strain-life constants sigma_f = 1312.6 MPa, b = -0.0636, epsilon_f = 3.0235, c = -0.888: sheet
  `FAB501`, strain-controlled R = -1, LabControl rig.

**Open points.**
- **The 1e6-cycle strengths lie below the run-out levels.** They are extrapolations of a Basquin fit,
  not fatigue limits, so Sigm1 and Taum1 are probably conservative for the limit-based criteria.
  - Axial: FAB060 has run-outs at 600 MPa (2.6e6 and 1.3e7 cycles), its lowest failure is at 610 MPa,
    and FAB070 has failures down to 530 MPa; Sigm1 is 517 MPa.
  - Torsion: run-outs at 398–419 MPa against Taum1 = 366.5 MPa.
- The two axial curves of the same series differ: FAB060 runs out at 600 MPa, while FAB070 fails at
  530 MPa.
- The strain-life constants come from LCF tests (lives of 358–11,867 cycles). Models that read
  `sigma_f` and `b_fat` from the Summary (e.g. MDC) use them instead of a fit to the HCF data.
- Whether the shear stress of the hollow specimens refers to the outer surface or to the mean wall
  radius ("related to net section") is not stated. Values are used as given.

## 19_42CrMo4_FABEST_bar.xlsx (FABEST database, PBS solid bars, added 2026-11-23)

**Source.** As `16_42CrMo4_FABEST.xlsx`: FABEST database FAB_v20261014, doi:10.5281/zenodo.23359038,
CC BY 4.0, same SHA-256. Converted by `tools/convert_fabest.py --bar` (`convert_bar`), with the
same reading rules; `tests/test_fabest_conversion.py` checks the committed workbook.

**Scope (maintainer decision d119, 2026-11-23).** The PBS mean-stress series of the database are on
solid hourglass bars (drawing v05, d = 6 mm, turned), not on the tube of `16_42CrMo4_FABEST`. They
form a dataset of their own, with the fully reversed curve of the same bar as design set, so that the
assessed tests differ from the design tests only in the mean stress:

| Curve | Loading | Failures | Run-outs |
|---|---|---|---|
| FAB110 | axial, R = -1 (design set) | 13 | 3 |
| FAB111 | axial, R = 0 | 10 | 3 |
| FAB112 | axial, R = 0.3 | 8 | 1 |
| FAB113 | axial, R = -2 | 10 | 0 |
| FAB114 | axial, R = -5 | 10 | 1 |
| **Total** | | **51** | **8** |

**Summary.**
- E, nu, Su, Sy and the strain-life constants: as `16_42CrMo4_FABEST` (same material batch).
- Sigm1 = 412.3 MPa: the authors' Basquin strength of FAB110 at 1e6 cycles (`SN_FL`).
- Taum1 = 366.5 MPa: the torsion strength of the PBS tube curves FAB062/FAB072, as in
  `16_42CrMo4_FABEST`. FABEST has no torsion curve on the bar; the value is from another specimen of
  the same batch and laboratory (maintainer decision d119).

**Open points.**
- **Taum1/Sigm1 = 0.89 is high for a steel** (the tube gives 0.71). The bar is weaker than the tube in
  axial loading (412 against 517 MPa at 1e6 cycles, both Basquin extrapolations below the run-out
  levels), and the torsion value comes from the tube. Criteria that use the ratio are affected: RFB
  and SIH leave their validity range (0.577-0.866) and are skipped; criteria whose mean-stress
  sensitivity follows from the ratio (e.g. GAM, NINIC) show large errors on this dataset.
- FAB111_2 (R = 0, amplitude and mean 600 MPa) has a maximum stress of 1200 MPa, above Su = 1111
  MPa, and failed after 19 cycles: a quasi-static failure, kept as supplied. It weighs heavily in
  the RMS error of the high-cycle criteria on this dataset.
- No torsion and no multiaxial tests: the dataset assesses the uniaxial mean-stress response only.
- The other unnotched mean-stress series stay out: CTU R = -3 (FAB010, solid bar v01; a separate
  dataset with FAB001/014/081 as design set is possible) and UNIOVI R = 0.1 (FAB090-093, dog-bone
  bars, no fully reversed curve to design on).

## 18_Al7075_T651_Wang.xlsx (Wang, Shang, Sun & Chen 2018, added 2026-11-16)

**Source.** Wang X-W, Shang D-G, Sun Y-J, Chen H. Multiaxial high-cycle fatigue life prediction
model considering mean shear stress effect under constant and variable amplitude loading. Theor
Appl Fract Mech (2018), doi:10.1016/j.tafmec.2017.10.007 (accepted manuscript). Tables 3-5 (pp.
34-36) are transcribed in `tests/data/papers/wang2018_tables3_5.csv`, also used by the WANG2018
paper-reproduction tests.

**Conversion.** `tools/build_wang_al7075_t651.py` builds the workbook from the transcription and
Tables 1-2. `tests/test_wang_al7075_dataset.py` checks that:
- the committed workbook equals the build;
- the table sizes are 13, 16 and 23 tests;
- the transcribed Table 4 amplitudes and means add up to its printed tau_max.

**Content (maintainer decision, 16 November 2026).** 7075-T651 aluminium alloy, thin-walled tubes
(outer diameter 12.5 mm, inner 10.5 mm), load controlled, sine 10 Hz, constant amplitude.
52 failures:
- uniaxial: R = -1 (5 tests) and sigma_m about 1.13 sigma_a (8);
- pure torsion: R = -1 (8) and R = 0.06 (8);
- tension-torsion with both means, at 0, 30, 45, 90 and 180 deg (23).

The shear lags the axial stress by delta (Eq. 8).

| Summary value | Source |
| --- | --- |
| E = 71 700 MPa, Sy = 501, Su = 561 MPa | Table 1 |
| nu = 0.3036 | Table 1: E/(2G) - 1 with G = 27.5 GPa |
| Sigm1 = 221.99 MPa (push-pull) | Table 2 axial curve 888.90 N^-0.10042 at N = 1e6 |
| Taum1 = 127.10 MPa | Table 2 torsional curve 688.85 N^-0.12233 at N = 1e6 |

**Open points.**
- The fatigue limits are S-N curve values at 1e6 cycles, the first of the paper's three choices
  (1e6, 2e6, 1e7). Aluminium has no true fatigue limit.
- No R = 0 limits are given. Criteria needing Sig0/Tau0 use their own estimates or are skipped.
- Whether the shear stresses refer to the outer surface or to the mean wall radius is not stated.
  The values are used as given.

## 17_Ti6Al4V_Delahay.xlsx (Delahay & Palin-Luc 2006, added 2026-11-11)

**Source.** Delahay T, Palin-Luc T. Estimation of the fatigue strength distribution in high-cycle
multiaxial fatigue taking into account the stress-strain gradient effect. Int J Fatigue 28 (2006)
474-484, doi:10.1016/j.ijfatigue.2005.06.048. Table 1 (p. 478) and the Ti-6Al-4V rows of Table 2
(p. 479) are transcribed in `sources/Delahay2006/ti6al4v_tables.csv`.

**Conversion.** `tools/build_delahay_ti6al4v.py` builds the workbook from the transcription;
`tests/test_delahay_dataset.py` checks that the committed workbook equals the build and that the
transcription agrees with quantities the paper derives from it: sigma* = 563 MPa from Eq. (6), and
tau*/tau_ref = 384/411 = 0.934 (Section 4.2).

**Scope (maintainer decision d109).** Annealed Ti-6Al-4V (E 110 GPa, nu 0.3, Sy 980, Su 1090 MPa),
smooth specimens, median staircase endurance limits (P_f = 0.5, at least 15 specimens, 1e6 or
1e7 cycles; entered at N = 1e6), all R = -1:

| Load | sigma_a (MPa) | tau_a (MPa) | s (MPa) | Phase | Role |
|---|---|---|---|---|---|
| tension | 583 | - | 17 | - | information (not in the workbook) |
| rotating bending | 602 | - | 34 | - | information (not in the workbook) |
| plane bending | 652 | - | 20 | - | Sigm1 (loading "plane bending") |
| torsion | - | 411 | 27 | - | Taum1 |
| plane bending + torsion | 442 | 255 | 47 | 0 deg | assessed |
| plane bending + torsion | 567 | 328 | 17 | 90 deg | assessed |

Sigm1 is the plane-bending limit, not the tension one: the assessed points are plane bending +
torsion, so calibrating on the same load types lets the stress gradient largely cancel in a
pointwise criterion. Every Dataset row is a limit state (`Run out = Yes`), as in the other
fatigue-limit workbooks; the bending stress is the nominal surface amplitude.

**Open points.**
- The paper does not give the specimen diameter per test type (the size-effect figures name 7.98
  mm for the reference specimens); the diameter is Null.
- The cycle basis is "1e6 or 1e7 cycles" without saying which applies to each limit.
- The authors note that the 90 deg limit was measured on specimens of a different roughness.
- These limits were formerly in `1_Ti6Al4V_Meneghetti.xlsx` as 4 run-outs, mislabelled (removed there, d109).

## 14_42CrMo4_Reis.xlsx: R = -1 fatigue limits (overlay, 2026-11-10)

**Maintainer decision d100 (issue #71).** The Summary had no fatigue limits, so the 22 limit-based
methods and WANG2022 skipped this batch. As for FABEST, the limits are the authors' strengths at
N = 1e6 from their own S-N trend lines: Anes, Reis, Li, Fonte, de Freitas, Int. J. Fatigue 62 (2014)
21-33, doi:10.1016/j.ijfatigue.2013.05.004, Table 5 p. 28, case 1 at N = 1e6:

| Limit | Value | Trend line (Table 4, p. 27) |
| --- | --- | --- |
| Sigm1 (push-pull, R = -1) | 452 MPa | sigma_a = 1204.3 N^-0.071 (case 1, pure axial) |
| Taum1 (R = -1) | 372 MPa | tau_a = 864.78 N^-0.061 (case 2, pure shear) |

The IJF 62 paper analyses the same Reis tests as the two Summary sources (the SSF surface is fitted
on them, `Calibrations/42CrMo4_Reis/SSF_published.json`).

- The supplied workbook is unchanged. The values are a new Summary source row (row 4, with the
  column `Axial fatigue limit loading` = push-pull), added by
  `src/multifat/imp/overlay_files/14_42CrMo4_Reis.overlay.json`; `--no-overlay` returns the
  workbook without limits.
- The longest R = -1 failures reach 2.7e5 (axial) and 3.2e5 cycles (torsion), so the values at
  1e6 extrapolate the authors' trend lines; they are finite-life strengths, not staircase limits.
- Ratio Taum1/Sigm1 = 0.823, against 0.794 (QT) and 0.709 (FABEST) (issue #71).

## 42CrMo4_Santus.xlsx

| Source (as cited in workbook) | Plain failures | Plain run-outs | Notched | Where the data are printed | How read |
|---|---|---|---|---|---|
| Santus et al. 2018, IJF 113, 113-125, [10.1016/j.ijfatigue.2018.03.037](https://doi.org/10.1016/j.ijfatigue.2018.03.037) | 17 (axial, R = -1) | 0 | 14 (sharp V, axial) | Fig. 5(a), p. 5 | digitised |
| Benedetti et al. 2022, IJF 159, 106787, [10.1016/j.ijfatigue.2022.106787](https://doi.org/10.1016/j.ijfatigue.2022.106787) | 28 (axial R = -1: 10; torsion R = -1: 9; torsion R = 0.1: 9) | 11 | 23 | per-specimen values not printed; Fig. 11(a) and Table 12 (p. 12) | author or supplementary data |
| Anes et al. "2013", [10.1016/j.ijfatigue.2013.03.014](https://doi.org/10.1016/j.ijfatigue.2013.03.014) (IJF 58 (2014) 102-113) | 11 (axial 6, torsion 5) | 2 | 0 | Table 3, p. 6 (PT, PS) | tabulated |
| Anes et al. 2014, FFEMS 37, 1338-1354, [10.1111/ffe.12192](https://doi.org/10.1111/ffe.12192) | 14 (in-phase 7, 90° out-of-phase 7) | 2 | 0 | Table 2 cases 1-2, p. 11; the same numbers are IJF 58 Table 3 (PP, OP) | tabulated |

Both Anes papers republish data from L. Reis, PhD thesis, IST Lisbon (2004) (ref. 27 of the FFEMS paper).
The Anes stresses are amplitudes at R = -1: the FFEMS paper normalises the waveforms by σa and √3·τa (Fig. 6).

Material constants printed in the sources:

- Santus 2018 Table 1 (p. 4): E = 206 GPa, Sy = 727 MPa, Su = 875 MPa (42CrMo4+QT, Trento bars).
  Table 2 (p. 6), plain axial, σa = k·N_f^b: R = -1: k = 1163 MPa, b = -0.080139; R = 0.1: k = 578.9 MPa,
  b = -0.039655.
- Anes, IJF 58 Table 1 (p. 2), measured by the authors on the tested 42CrMo4: E = 206 GPa, ν = 0.3,
  Su = 1100 MPa, Sy = 980 MPa, σ'_f = 1154 MPa, b = -0.061, ε'_f = 0.180, c = -0.53.
- Benedetti 2022 Table 12 (p. 12), plain R = -1: torsion τa = 681.6·N_f^-0.048 MPa; axial
  σa = 103 + 859.9·N_f^-0.066 MPa; strengths at 5×10^6 cycles: 412 MPa (axial), 328 MPa (torsion).
  The paper says its supplies are distinct from those of its ref. [17]; it does not say they are the Santus
  2018 bars. Checked 2026-11-11: ref. [17] is Benedetti, Santus, Berto, IJF 126 (2019) 306-318, not
  Santus 2018 (which is its ref. [15]); the statement concerns the 2019 material only.
- **Maintainer decision d101 (2026-11-10, issue #71):** the Santus 2018 and Benedetti 2022 tests of this
  workbook are one batch, and `42CrMo4_Santus` is treated as one batch in the cross-batch comparison
  (`examples/data/material_families.json`).

Known issues:

1. The Santus row of the Summary sheet carries k_FS = 2.07, k_FI = 0.83, L_LI = 5, k_MGSE = 1.15,
   a_MSWT = 0.53, alpha_w = 0 and xi_ZHU = 1. Santus 2018 reports axial tests only, and none of these values
   appears in it. ν = 0.33 is not printed there either.
2. The Summary average mixes two heat treatments: Su 987.5 = (875 + 1100)/2 and Sy 853.5 = (727 + 980)/2
   (Santus QT vs Anes/Reis). The E standard deviation (5900 MPa) is Santus' printed st. dev.
3. The Benedetti 2022 fatigue limits in the Summary (433 / 344 MPa) differ from Table 12 (412 / 328 MPa).
4. The 9 plain torsion R = 0.1 failures and 2 run-outs credited to Benedetti 2022 are not in that paper. The
   points at 600 MPa / 3620 cycles (axial) and 480 MPa / 2113 cycles (torsion) lie outside Fig. 11 (its axis
   starts at 10^4 cycles). The paper's supplementary files are MATLAB scripts for the control-radius
   search and contain no test data, so these rows and the 433 / 344 MPa limits remain unsourced.
5. Row `Plain_multiaxial_9` (σa 610, τa 352.2 MPa, N 4114) was tagged 90° but is a proportional test (PP /
   case 1). Fixed on 2026-10-18 (phase 0).
6. The Anes 2013 Summary row leaves Su and Sy empty, although Table 1 prints 1100 and 980 MPa.

## AISI316L.xlsx

### PSI stress amplitudes corrected from original records (29 October and 3 November 2026)

Following Arda Akay's correspondence, the eight Janssens PSI tests introduced
by the in-memory overlay have been reconciled with the original records
rather than using the rounded stresses printed by Janssens et al. (2018).
Five uniaxial axial amplitudes require the **1.1558 load-cell multiplier**.
Their new stress amplitudes in MPa are:

| Specimen | Published stress amplitude (MPa) | Corrected axial amplitude (MPa) |
| --- | ---: | ---: |
| `Uniaxial_10` | 240 | 277.4 |
| `Uniaxial_11` | 260 | 300.5 |
| `Uniaxial_12` | 290 | 335.2 |
| `Uniaxial_13` | 300 | 346.7 |
| `Uniaxial_14` | 330 | 381.4 |

For `Non_prop_12`, the **unrounded original axial amplitude is 233.35 MPa**
(233 MPa as printed). Applying 1.1558 yields **269.70593 MPa**, stored
at 0.1 MPa precision as **269.7 MPa**. Its original shear amplitude is
**134.7 MPa** (135 MPa as printed). The original-record torsional amplitudes
of `Torsional_8` and `Torsional_9` are **150.1 and 155.9 MPa**
(respectively 150 and 156 MPa as printed). These three shear adjustments
reflect **publication rounding only**; no axial load-cell correction is
applied to the shear stress amplitudes. Phase angles and fatigue lives are
unchanged.

Source: Arda Akay, 29 October 2026 email and follow-up clarification provided
3 November 2026, based on original PSI test records. The in-memory overlay
records both published amplitudes and source-precision values. The original
`8_AISI316L.xlsx` file remains unchanged and is identified through its
existing `base_sha256`. The corrected post-overlay dataset fingerprint is
`e245074653e53c978ee11a94812ba9d9fd6e48c0721b98a6f80c5f7a6ee2d9f5` (regenerated from the importer); it supersedes the intermediate
`d5a42243370f4944595891f69a87d1e47d9532aeca45594749ab15edf7836fef`
fingerprint. Historical calibration records are not automatically migrated
to the new dataset fingerprint.

This correction changes experimental stress histories for all fatigue
criteria, not only TVEIT2025. The distinct Janssens tube and Braun
wrought-sheet material/batch sources remain subject to Issue #24.

**Janssens final revision, supplied 28 October 2026.** The final author workbook (`8_AISI316L(1).xlsx`, original SHA-256 `a70a09f76fb8fa22f0d6ef9843512eefd20db8fae980639067ef6d207187a35a`) replaces the previous revision. Its reported R = -1 tensile fatigue limit is 278.5 MPa (previously 241); the torsional limit remains 144 MPa. Its first uniaxial run-out is 17,225,000 cycles (previously 11,000,000), and non-proportional torsional phases in Dataset rows 18-28 are now 90 degrees. The Summary averages and derived values are also updated. All underlying ZIP members of the workbook are preserved verbatim, although the repository copy is recompressed at ZIP-container level. The repository XLSX SHA-256 is `26a97aec4603aa19df29bbd6576c3b72259cc3792eaa6e60df1d5e977f52017c`. Because the 11 phase corrections are now part of the source, the active overlay has 15 remaining operations. The new post-overlay fingerprint is `a63f9325d41155bb7611b2c6272e0bbee9a46e33fefe4bc33faa47f8470f4750` (39 failure observations and 13 run-outs). The previous overlay mapping remains in the historical overlay JSON. No other experimental workbooks are changed.


**Author-corrected workbook received on 28 October 2026.** The maintainer confirmed
that the new `8_AISI316L.xlsx` is the authoritative Janssens-corrected experimental input.
The workbook is preserved byte-for-byte (SHA-256
`b637deefb415243c5a071a38f188c02202501d2c6e85b87c30af41523a211949`). Relative to the previous input, it
contains revised geometry and stress values; these observations are not reverted or
recalculated. Its `Summary` still reports Janssens' R = -1 tensile and torsional
limits as 241 and 144 MPa, respectively, alongside Braun's R = 0 tensile stress
range of 387.8 MPa. The importer selects the R = -1 source values and does not use
the mixed-ratio `Average` value as an endurance limit.

The in-memory overlay retains the documented independent corrections (bulk elastic
properties, non-proportional phases and eight additional published Janssens cases).
On the author-corrected workbook it yields 39 failures and 13 run-outs. The overlay
now targets the corrected workbook hash and records the regenerated dataset fingerprint
`b139cef997c825a12b8f7c0e9f4e72402addbf193222f4506db8b7a4caba666b`. The prior overlay is preserved without
modification as `8_AISI316L.before-author-correction.overlay.json` to retain the
original input-to-correction mapping. Its stored fingerprint refers to the earlier
importer behaviour and should not be interpreted as a validation of current results.
Neither source workbook values nor existing calibration records were modified.
Material-batch comparability and the interpretation of the original 241 MPa fatigue
limit after the corrected stress amplitudes remain matters for scientific review.


| Source | Plain failures | Plain run-outs | Notched | Where printed | How read |
|---|---|---|---|---|---|
| Janssens 2018, IJF 110, 42-48, [10.1016/j.ijfatigue.2017.12.013](https://doi.org/10.1016/j.ijfatigue.2017.12.013) | 25 (axial 11, torsion 7, 90° non-proportional 7) | 10 (Nf > 11×10^6, stored as 11e6) | 0 | Table 4, p. 4 | tabulated |
| Braun et al. 2024, Mater. Des. 244, 113122, [10.1016/j.matdes.2024.113122](https://doi.org/10.1016/j.matdes.2024.113122) | 14 (axial, R = 0, wrought sheet) | 3 | 42 (V and U notches) | Fig. 13(i), p. 10 | digitised, σa = σm = Δσ_nom/2 |

Material constants printed in the sources:

- Janssens Table 3 (p. 3), SA-312 316L tube: Rp0.2 = 265 ± 5 MPa, Rm = 585 ± 5 MPa, σ-1 = 241 MPa,
  τ-1 = 144 MPa. No S-N or strain-life constants.
- Braun Fig. 13(i) (p. 10), wrought plain, R = 0: Δσ_R,50% = 387.8 MPa at 10^6 cycles, k = 8.33, for
  N = 10^6·(Δσ/Δσ_R)^-k (Eq. (2)).

Changes made on 2026-10-18 (each changed cell carries an Excel comment):

- **8 Janssens points added** (`Uniaxial_10` to `Uniaxial_14`, `Torsional_8`, `Torsional_9`, `Non_prop_12`):
  uniaxial 240 MPa / 37k, 260 / 14k, 290 / 16k, 300 / 21k, 330 / 4k cycles; torsional 150 / 318k,
  156 / 29k; non-proportional 233 / 135 MPa, 74k. They are printed in Table 4 but were absent from the
  original workbook, and the paper states no exclusion criterion (tests with σa ≥ 300 MPa ran at 1 Hz).
  Each row copies its group's geometry and formulas; all other rows are unchanged.
- **E and ν taken from the literature**, because Janssens prints neither: E = 193,000 MPa and ν = 0.3
  (Zhu et al. 2019, IJF, [10.1016/j.ijfatigue.2018.11.025](https://doi.org/10.1016/j.ijfatigue.2018.11.025),
  p. 8, 316L). Consistent with Liu et al. 2021 (IJPVP, Table 1: 197 GPa, 0.3) and Hwang et al. 2020
  (IJF: 195 GPa measured on SA-312 TP316). The previous values, 174,000 MPa and 0.27, had no source.

Known issues:

1. k_FS, k_FI, L_LI, k_MGSE and a_MSWT in the Summary are not in Janssens 2018.
2. The non-proportional Janssens rows had an empty phase, but the paper states a 90° phase shift. Fixed on
   2026-10-18 (phase 90).
3. Braun's text gives the run-out limit as 2×10^6 cycles, but Fig. 13(i) draws the run-outs near 10^7, where
   they were digitised.
4. Janssens (tube) and Braun (hot-rolled sheet) are different batches, and the Summary has Janssens values
   only.

## Al7075_T6.xlsx

| Source | Plain failures | Plain run-outs | Notched | Where printed | How read |
|---|---|---|---|---|---|
| Pedranz et al. 2025, Int. J. Mech. Sci. 288, 109961, [10.1016/j.ijmecsci.2025.109961](https://doi.org/10.1016/j.ijmecsci.2025.109961) | 87 (axial R = -1/0.1: 11/11; torsion R = -1/0.1: 12/11; multiaxial 42) | 14 | 88 | supplementary Appendix B, Tables 6-22 | tabulated |
| Benedetti et al. 2021, FFEMS 44, 2975-2989, [10.1111/ffe.13530](https://doi.org/10.1111/ffe.13530) | 0 (15 axial R = 0.1 tests on the phi60 bar, plotted only in Fig. 9B; not imported, d111) | 0 | 0 | Summary sheet only | - |

All 189 rows match Appendix B point by point (cycles and amplitude, in the right series), except one notched
blunt-torsion run-out stored as 5,000,000 cycles where Table 13 prints 4,800,000.

Material constants printed in the sources (same 60 mm extruded bar; Pedranz Table 2 takes its properties from
Benedetti 2021):

- Benedetti 2021 Table 3 (p. 9), φ60: E = 70.5 ± 0.2 GPa, Sy = 531 ± 7.3 MPa, Su = 595 ± 5.9 MPa;
  ν = 0.33 (p. 4).
- Benedetti 2021 Table 4 (p. 10), φ60: σ'_f = 1067 MPa, b = -0.0893, ε'_f = 0.131, c = -0.674 (from
  strain-controlled axial LCF). Cyclic Ramberg-Osgood in the form ε = σ/E + (σ/H)^(1/n): H = 617 MPa,
  n = 0.0233. This n is not a conventional n' and needs a convention check before use.
- Pedranz Appendix B Table 5, plain series, with Eq. (9) σa = kσ1/N_f^kσ2 + kσ3 and Eq. (10)
  τa = kτ1/N_f^kτ2: axial R = -1: 1344.8 MPa, 0.206, 88.2 MPa; axial R = 0.1: 4761 MPa, 0.406, 94.5 MPa;
  torsion R = -1: 355.9 MPa, 0.0817; torsion R = 0.1: 359.3 MPa, 0.0913. Strengths at 5×10^6 cycles:
  143.8 MPa (axial, R = -1) and 100.9 MPa (torsion, R = -1).
- Benedetti 2022 Table 12 (p. 12) prints the same plain axial fit (1345, 0.206, 88.2), so that series is
  shared with Pedranz. Its torsion fit differs: τa = 403.6·N_f^-0.090 MPa, strength 99.2 MPa.

Known issues:

1. The Benedetti 2021 Summary row carries k_FS = 2.09, k_FI = 0.39, L_LI = 3.42, k_MGSE = 0.64 and
   a_MSWT = 0.39. That paper reports axial tests only, and none of these values appears in it.
2. The Pedranz fatigue limits in the Summary (150 / 100 MPa) differ from Appendix B Table 5 (143.8 /
   100.9 MPa). Resolved 2026-11-11 (decision d111): `src/multifat/imp/overlay_files/Al7075_T6.overlay.json`
   sets the Table 5 strengths at N = 5e6, R = -1, push-pull (the paper's basis, Tables 3-4); 150 / 100
   is printed in neither source and fits neither 5e6 nor 2e6 (155.9 / 108.8 MPa from the same fits).
   The same overlay adds the measured plastic strain-life terms of Benedetti 2021 Table 4 (eps'_f =
   0.131, c = -0.674; gamma'_f = sqrt(3) eps'_f, c0 = c by the rule of `strain_life_fill`), which
   replace the steel uniform-material-law estimates. sigma'_f and b stay the Basquin fit of the HCF
   tests: the LCF values (1067 MPa, -0.0893) give 311 MPa at 2N = 1e6, against about 145 MPa
   measured. The Summary name becomes Al 7075-T6, as both papers write it (T6511 once in Pedranz).
3. Workbook group R = -1, φ = 0° (19 rows) merges two biaxiality ratios (λ = 1/√3 and λ = 1).

## Ti6Al4V_andrea_(AXIAL).xlsx

| Source | Plain failures | Plain run-outs | Notched | Where printed | How read |
|---|---|---|---|---|---|
| Meneghetti et al. 2018, TAFM 96, 509-533, [10.1016/j.tafmec.2018.06.010](https://doi.org/10.1016/j.tafmec.2018.06.010) | 10 (axial 6, torsion 4) | 0 (figure shows 3) | 16 (sharp V, axial 6, torsion 10) | Fig. 10(a), p. 17; counts in Table 3 | digitised, verified |
| Berto et al. 2014, FFEMS 38 (2015) 503-517, [10.1111/ffe.12272](https://doi.org/10.1111/ffe.12272) | 48 (axial R = -1: 9; torsion R = -1, 0, 0.25, 0.5, -2, -3: 39) | 0 | 0 | Figs 4-6; series in Table 3, p. 5 | digitised |

Material constants printed in the sources:

- Meneghetti Table 1 (p. 6): σy = 915 MPa, σu = 980 MPa, E = 110,000 MPa, ν = 0.3, A = 13 %, Z = 29 %.
  Table 3: plain axial Δσ_A = 1014 MPa, k = 50; plain torsion
  Δτ_A = 894 MPa, k = 22 (at 2×10^6 cycles, Ps = 50%). The text (p. 19) gives 1016 and 847 MPa instead.
- Berto Table 1 (p. 3): Su = 978 MPa, Sy = 894 MPa; text (p. 9): E = 110 GPa, ν = 0.3. Table 3 (p. 5)
  prints the inverse slope k and the amplitude at 10^6 / 2×10^6 cycles for each series, e.g. tension R = -1:
  k = 9.25, 512.76 / 475.74 MPa; torsion R = -1: k = 22.13, 400.60 / 388.25 MPa.
- No strain-life constants are printed in either paper.

Known issues:

1. The Summary torsional limit 422 MPa is the level of the plain torsion run-out in Meneghetti Fig. 10(a)
   (τa = 422.3 MPa at 2×10^6 cycles, consistent with the text's Δτ_A = 847 MPa). The axial 500 MPa matches
   no printed value: Table 3 gives 507, the text 508, and the two axial run-outs sit at 503 and 507 MPa.
2. Berto Table 3 lists 72 plain specimens; 48 are in the workbook, with no run-outs.
3. The Meneghetti plain rows agree with the Fig. 10(a) markers (read from the vector graphics) within about
   2 % in cycles and 1 % in stress. The figure shows 4 torsion failures, as in the workbook, although Table 3
   implies 5. The 3 plain run-outs (2 axial, 1 torsion) are not in the workbook.
4. k_FS, k_FI, L_LI, k_MGSE, a_MSWT and xi_ZHU are not in either paper.
5. The Berto Summary row leaves E and ν empty, although the paper prints them.
6. τa/σa ≈ 1 at 10^5 cycles in the Meneghetti plain data (checked 2026-11-11, issue #24) is a property
   of the data, not an error. The torsion values are nominal net-section shear amplitudes (Eq. 12b,
   p. 11). The axial curve is almost flat (k = 50) and the torsion curve steeper (k = 22), so the two
   curves converge at short lives; at 2×10^6 cycles Table 3 gives τA/σA = 894/1014 = 0.88 (0.83 with
   the text's values).

## Further edits applied on 2026-10-18 (approved by T. Grossi; Excel comment on each cell)

| Workbook | Cell | Old | New | Reason |
|---|---|---|---|---|
| 42CrMo4_Santus | Dataset O103 (`Plain_multiaxial_9`, phase torsional) | 90 | 0 | Proportional test: PP row of Anes IJF 58 Table 3 = case 1 of Anes FFEMS 37 Table 2 |
| AISI316L | Dataset O25-O36 (`Non_prop_1`-`12`, phase torsional) | Null | 90 | Janssens Table 4, "Non-Prop 90° Phase Shift" |
| 42CrMo4_Santus | Summary P2 (a_MSWT, Santus row) | 0.53 | 0.3605 | a' = a/(2 - a): printed MSWT weighting (Jiang 2000 Eq. 6 / Zhu 2018 Eq. 13); lives unchanged |
| AISI316L | Summary P2 (a_MSWT, Janssens row) | 0.36 | 0.2195 | same |
| Al7075_T6 | Summary P2 (a_MSWT, Benedetti 2021 row) | 0.39 | 0.2422 | same |
| Ti6Al4V_andrea_(AXIAL) | Summary P2 (a_MSWT, Meneghetti row) | 0.63 | 0.4599 | same |

The Average rows are formulas and now read the converted a_MSWT values.

## TODO

- [ ] **Benedetti 2022 data in `42CrMo4_Santus.xlsx` (source verification needed):** The paper and its
  supplementary files (two MATLAB scripts, no test data) do not contain: the 9 plain torsion R = 0.1
  failures and 2 run-outs; the points at 600 MPa / 3620 cycles (axial) and 480 MPa / 2113 cycles (torsion);
  the Summary fatigue limits 433 / 344 MPa (Table 12 prints 412 / 328 MPa). All are kept unchanged until the
  source is known.
- [ ] Open numerical decisions on the other workbooks (A7 Su/Sy mix, misattributed calibration constants,
  Pedranz and Meneghetti fatigue limits,
  untranscribed run-outs): listed in the issues above, with evidence in the local audit notes
  (documented in the literature audit records).

## Constants for decision D (source first)

Strain-life constants (σ'_f, b, ε'_f, c) are printed for **42CrMo4 (Anes/Reis batch)** and **7075-T6 (φ60
bar)** only. No source prints the torsional constants τ'_f, b0, γ'_f or c0 for any material, so those must
be fitted or skipped.


## Observation type column, decision d113 (12 November 2026)

At the maintainer's request every Dataset sheet gained the column
`Observation type (failure/run-out/fatigue limit)` (issue #18). Run out `No` rows are `failure`;
Run out `Yes` rows are `run-out` in workbooks that hold S-N curves (at least one failure) and
`fatigue limit` in workbooks that hold fatigue limits only. `tools/add_observation_type.py` added
the column to fifteen workbooks by editing only the XML of the Dataset sheet and of its Excel table
(formulas, cached values, styles and all other ZIP members unchanged; the table is extended to the
new column); `16_42CrMo4_FABEST.xlsx` and `17_Ti6Al4V_Delahay.xlsx` were rebuilt by their
conversion scripts. Every other cell value is unchanged, as are the dataset fingerprints. The active
overlays are re-pinned to the new SHA-256 with unchanged operations; the eight AISI316L `add_case`
rows also declare `failure`. `MANIFEST.json` holds the new hashes.

| Workbook | Previous SHA-256 | SHA-256 |
| --- | --- | --- |
| `10_En3B.xlsx` | `16a56e868fec91f5…` | `0f21e6c41a65a83afbe25dc359d1e0b9fe04c48ee04640568e106275f4e7bf94` |
| `11_GJS600.xlsx` | `d4c82e5c80ca604f…` | `39c07d6f7e9d15d5c3b4416f1c3a581bb2e7ee619fe777408880aaa9794b89b7` |
| `12_GGG60.xlsx` | `e53345f03ea4dffa…` | `367bc0e8f28bc134f694356812fcd3f0f639e6ff7160c0c9ba71cc4e01d1c45a` |
| `13_S65A.xlsx` | `dd03c3b0a89906fa…` | `efa0778484ad51afc4e1186953b4e480316c007043de27cf703b042b4c91aab0` |
| `14_42CrMo4_Reis.xlsx` | `8a2f87a325e2925e…` | `3fa6c5d8a20e54a39dd8272d17285fcb6f4708c39a85dd917c0d46c530b053e3` |
| `15_Ti6Al4V_Berto.xlsx` | `4cb8262d51574e38…` | `500fb0dc62d63afafb99b9a632824b063b6ee9000759527840ace10cb0caf765` |
| `16_42CrMo4_FABEST.xlsx` | `04afa35fc8b3c953…` | `242b6098c4e8fce2065c1116c1c0c885d84d7deeb2dbcb7731b075c399e87a25` |
| `17_Ti6Al4V_Delahay.xlsx` | `a67fbf23ff5ecd1e…` | `fd0e79f7845f76540d57319fa07b94be8b7c5ba47c6f4f5f38da3a2dd703df30` |
| `1_Ti6Al4V_Meneghetti.xlsx` | `86f0f71c50dffa3b…` | `78808b28fc8546bf3361b121c2f768dc5b4291be53b68b79d90e8cab1ece04dd` |
| `2_CX18.xlsx` | `0f8040cc6c00c363…` | `f9fd79fdfb657e41050be1c867eb7125f05795c45638d00d40f7226612ce9fc1` |
| `3_25CrMo4.xlsx` | `ba08322b49592323…` | `3fcc36daf35c12806b467b1005d4d5450d4284b618b758597705cc34c3afe17b` |
| `4_30NCD16.xlsx` | `debf822173d4914a…` | `0deaa15c8d667ec409483c78e9e28cd85a4af767f977ae9b2412e40ff049595e` |
| `5_34Cr4.xlsx` | `9352b72a54ddb872…` | `7f3e60c008e6449beb3f13ce0306f004c0e413ae2c8211f14ca895b648987edf` |
| `6_42CrMo4_Santus.xlsx` | `87f40c72a2e4d766…` | `1327f64f68305d4e4453e784803a5bccce7596dc17db0ab032125598a1fd93d6` |
| `7_982FA.xlsx` | `f786dc8e99f598dd…` | `cb91fc0b0c2a2f7ca88c44b7686cf68a389cc6b2a9d5c382fb754af731f8ac34` |
| `8_AISI316L.xlsx` | `26a97aec4603aa19…` | `f822e4c39deb6bacb26c2019bda4744ebaf4c7a327c808ccfc19c9b7358821ee` |
| `Al7075_T6.xlsx` | `2c250c383cfe6476…` | `2a378823426639986ba534f6fda4adefce621d1757def2531f157818306f94fc` |

## AISI316L tensile fatigue limit adjustment, decision d114 (18 November 2026)

Following Issue #102, the tensile fatigue limit $\sigma_{-1}$ in `8_AISI316L.xlsx` (Janssens row) was adjusted from the published $241 \times 1.1558 = 278.5\text{ MPa}$ to $271.6\text{ MPa}$ (the staircase load level $235\text{ MPa} \times 1.1558$, as used in the PSI calibration of Tveit) via a `set_summary` overlay operation in `8_AISI316L.overlay.json`.

The published value ($278.5\text{ MPa}$) exceeded 5 axial $R = -1$ failures (at $271.6$ and $277.4\text{ MPa}$). Adjusting to $271.6\text{ MPa}$ places 3 failures above the threshold and 2 directly on it.

The base workbook SHA-256 (`f822e4c39deb6bacb26c2019bda4744ebaf4c7a327c808ccfc19c9b7358821ee`) is unchanged. The updated post-overlay dataset fingerprint is:
- `5e2e99b9e3456f7d96a5dde8e66e4213e241cf88592a068d679506490824445d`

