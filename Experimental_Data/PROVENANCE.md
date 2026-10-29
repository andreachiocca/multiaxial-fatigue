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
| 42CrMo4_QT | `6_42CrMo4_QT.xlsx` | Santus and Benedetti rows only (45 plain failures, 11 run-outs); the Anes/Reis rows (Su 1100 batch) are not in it |
| Ti6Al4V | `1_Ti6Al4V.xlsx` | Meneghetti 10 failures and 4 Delahay run-outs; the 48 Berto failures are not in it |
| AISI316L | `8_AISI316L.xlsx` | Janssens and Braun rows as originally supplied (without N6, N8 and our E, nu) |
| Al7075_T6 | `Al7075_T6.xlsx` | measured data unchanged; Summary coefficient cells removed |
| eight new materials | `2_`-`5_`, `7_`, `10_`-`13_` | En3B and GJS600 notched only; CX18, 30NCD16, 25CrMo4, 982FA, 34Cr4, GGG60 and S65A have only fatigue-limit rows (`Run out = Yes` at the endurance basis) |

### Edits of 2026-10-19 (Excel COM, comment on every changed cell)

| Workbook | Change | Source | Decision |
|---|---|---|---|
| `14_42CrMo4_Reis.xlsx` (new) | the 29 Anes/Reis rows of the previous 42CrMo4_QT.xlsx (25 failures, 4 run-outs, including N1), in the revised layout; Summary rows Anes IJF 58 and FFEMS 37 with E 206000, nu 0.3, Su 1100, Sy 980 and the strain-life constants sigma_f 1154 MPa, b -0.061, epsilon_f 0.180, c -0.53 (new Summary columns) | Anes et al., IJF 58 (2014) Table 1 p. 2 | d055 |
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

## 42CrMo4_QT.xlsx

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
  2018 bars.

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

### PSI axial amplitude correction (29 October 2026, pending confirmation)

Arda Akay reported that six of the eight Janssens PSI observations added by the
in-memory overlay still carried the uncorrected, published axial stress
amplitudes. The overlay now applies the reported **1.1558 load-cell factor**
to `Uniaxial_10` through `Uniaxial_14` and `Non_prop_12`. The applied
amplitudes (MPa) are 277.4, 300.5, 335.2, 346.7, 381.4 and **269.3
(provisional)**, respectively, replacing 240, 260, 290, 300, 330 and 233.
The two torsional-only observations and the 135 MPa shear amplitude in
`Non_prop_12` are unchanged. The original `8_AISI316L.xlsx` file is
preserved; the correction is applied at import time. The overlay records
the original values, multiplier, correspondent's clarification and the
updated dataset fingerprint.

**Open check before merging:** Arda's email of 29 October gives 269.7 MPa
for `Non_prop_12`, whereas the stated multiplier applied to the published
233 MPa yields 269.3014 MPa. This PR provisionally uses **269.3 MPa** after
rounding to one decimal. Verify the original unrounded PSI stress or obtain
Arda's confirmation before treating the non-proportional correction as final.
The historical fingerprints and fatigue-life calibrations based on the
uncorrected dataset must not be silently reused.

This correction changes experimental stress histories for all criteria,
not Tveit alone, and does not resolve the separate Janssens/Braun material
and batch distinction tracked in Issue #24.

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
| Benedetti et al. 2021, FFEMS 44, 2975-2989, [10.1111/ffe.13530](https://doi.org/10.1111/ffe.13530) | 0 | 0 | 0 | Summary sheet only | - |

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
   100.9 MPa).
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

## Further edits applied on 2026-10-18 (approved by T. Grossi; Excel comment on each cell)

| Workbook | Cell | Old | New | Reason |
|---|---|---|---|---|
| 42CrMo4_QT | Dataset O103 (`Plain_multiaxial_9`, phase torsional) | 90 | 0 | Proportional test: PP row of Anes IJF 58 Table 3 = case 1 of Anes FFEMS 37 Table 2 |
| AISI316L | Dataset O25-O36 (`Non_prop_1`-`12`, phase torsional) | Null | 90 | Janssens Table 4, "Non-Prop 90° Phase Shift" |
| 42CrMo4_QT | Summary P2 (a_MSWT, Santus row) | 0.53 | 0.3605 | a' = a/(2 - a): printed MSWT weighting (Jiang 2000 Eq. 6 / Zhu 2018 Eq. 13); lives unchanged |
| AISI316L | Summary P2 (a_MSWT, Janssens row) | 0.36 | 0.2195 | same |
| Al7075_T6 | Summary P2 (a_MSWT, Benedetti 2021 row) | 0.39 | 0.2422 | same |
| Ti6Al4V_andrea_(AXIAL) | Summary P2 (a_MSWT, Meneghetti row) | 0.63 | 0.4599 | same |

The Average rows are formulas and now read the converted a_MSWT values.

## TODO

- [ ] **Benedetti 2022 data in `42CrMo4_QT.xlsx` (source verification needed):** The paper and its
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
