"""Regression cases with analytical values, not just import/syntax checks."""

from multifat.models.cp_methods.mmccm import _validated_constants as _mmccm_constants, mmccm_life
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from multifat.models.registry import get_models, available_cp_methods
from multifat.models.references import DISABLED_METHODS, REFERENCES
from multifat.eval.evaluator import best_cp_for_case
from multifat.models.base import CPPlaneResult
from multifat.models.utils.plane_history import split_strain_ranges, shear_maximum
from multifat.report.export_csv import save_fatigue_csv
from multifat.eval.fitting import fit_power_law_with_survival_std
from multifat.eval.metrics import log10_dp_ratio
from multifat.cli import calibrate, compare


class MethodAuditTests(unittest.TestCase):
    def test_analytical_plane_values_all_retained_cp_models(self):
        # Fully reversed, engineering shear amplitude .02; tensor amplitude .01.
        s = np.array([[0.0, 0.0, 40.0], [0.0, 0.0, 0.0], [40.0, 0.0, 30.0]])
        e = np.array([[0.0, 0.0, 0.01], [0.0, 0.0, 0.0], [0.01, 0.0, 0.003]])
        p = dict(
            E=200000.0,
            nu=0.3,
            Su=600.0,
            Sy=250.0,
            Sigm1=240.0,
            Taum1=145.0,
            k_FS=0.5,
            k_FI=0.2,
            k_MGSE=0.4,
            a_MSWT=0.3,
            L_LI=2.0,
            sigma_f=200.0,
            tau_f=100.0,
            b_fat=-0.1,
            b0=-0.08,
            epsilon_f=0.3,
            c_fat=-0.6,
            gamma_f=0.5,
            c0=-0.55,
        )
        ga, ea, ta, sn = 0.02, 0.003, 40.0, 30.0
        ge = ta / (p["E"] / (2 * (1 + p["nu"])))
        ne = sn / p["E"]
        expected = dict(
            FS=ga * (1 + 0.5 * sn / 250),
            FIN=ta + 0.2 * sn,
            SWT=sn * ea,
            KUJAWSKI2021=sn * ea,  # no mean: identical to SWT
            MATAKE=ta + (2 * 145 / 240 - 1) * sn,
            MCD=ta + 145 / (2 * 600) * sn,
            MGSE_YU=ta * ga + sn * ea,
            MGSE_ZHU=ta * ga + 0.4 * sn * ea,
            MSWT=4 * 0.3 * sn * ea
            + 2 * 0.7 * ta * ga,  # Zhu 2018 Eq. (13): 2a de s + (1-a)/2 dg dt
            MKBM=ga + (1 + sn / 250) * ea,
            LI=ta * ga + 2 * sn * ea,
            GSE=ta * ga + sn * ea,
            GSA=ta / 100 * ge + (ga - ge) + sn / 200 * ne + (ea - ne),
            LIU2021=ga + sn * np.sqrt(sn * p["E"] * 2 * ea) / (4 * (p["E"] / 2.6) * 100),
            MDC=np.sqrt(
                ((2 * 100 / 200) ** 2 * (2 * sn / 200) ** (2 * 0.02 / -0.1) - 1) * sn**2 + ta**2
            ),
            MWCM=ta / (145 - sn / ta * (145 - 240 / 2)),
            MMCCM=1.0 / float(mmccm_life(ga, sn / ta, _mmccm_constants(p))),
            VSE=(2 * sn) * (2 * ea) + (2 * ta) * ga,
            LIU2025=np.sqrt(sn**2 + 3 * ta**2),
            ZHANGYAO=np.sqrt(ta**2 - ta * sn + (145 / 240) ** 2 * (240 + 2 * ta) * sn) / 145,
            MWHS=0.02 + 0.5 * (1 + sn / 250) * np.sqrt(sn * 2 * ea / p["E"]),
            NINIC=np.sqrt(
                (ta / 145) ** 2 + (240 / 145) ** 2 * (1 - (120 / 145) ** 2) * (sn / 240) ** 2
            ),
        )
        self.assertEqual(set(expected), set(available_cp_methods()))
        for name, value in expected.items():
            with self.subTest(name=name):
                model = available_cp_methods()[name]
                r = model.evaluate_on_plane(S0r=s, S1r=-s, E0r=e, E1r=-e, params=p)
                self.assertAlmostEqual(r.damage, value, places=11)

    def test_li_mean_shear_cancels_in_reversed_loading(self):
        from multifat.models.cp_methods.li import LI

        e = np.zeros((3, 3))
        e[0, 2] = e[2, 0] = 0.01
        for mean in (0.0, 10.0, -10.0):
            s0, s1 = np.zeros((3, 3)), np.zeros((3, 3))
            s0[0, 2] = s0[2, 0] = mean + 40
            s1[0, 2] = s1[2, 0] = mean - 40
            r = LI().evaluate_on_plane(S0r=s0, S1r=s1, E0r=e, E1r=-e, params={"L_LI": 1})
            self.assertAlmostEqual(r.damage, (40 + abs(mean)) * 0.02)

    def test_plane_ties_are_independent_of_order(self):
        class Tied:
            name = "TEST"

            def evaluate_on_plane(self, **kw):
                return CPPlaneResult(float(kw["S0r"][2, 2]), 1.0)

        r = np.array([np.eye(3), [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]])
        s = np.diag([2.0, 0.0, 1.0])
        z = np.zeros((3, 3))
        for planes in (r, r[::-1]):
            result = best_cp_for_case(
                S0=s, S1=-s, E0=z, E1=z, R_list=planes, method=Tied(), params={}
            )
            self.assertEqual(result.by_metric, 2.0)

    def test_cancelled_harmonic_shear_overrides_legacy_endpoints(self):
        from multifat.eval.evaluator import _replace_shear_components

        s = np.zeros((3, 3))
        s[0, 2] = s[2, 0] = 40.0
        mean = np.zeros((3, 3))
        mean[0, 2] = mean[2, 0] = 7.0
        zero = np.zeros((3, 3))
        hi, lo = _replace_shear_components(
            mean + s, mean - s, delta_eq=0.0, T_sin_r=zero, T_cos_r=zero
        )
        np.testing.assert_allclose(hi, mean)
        np.testing.assert_allclose(lo, mean)

    def test_liu_zero_shear_limit_is_continuous(self):
        from multifat.models.cp_methods.liu2021 import Liu2021

        model = Liu2021()
        p = dict(E=200000.0, nu=0.3, tau_f=100.0)
        s = np.diag([0.0, 0.0, 30.0])
        e = np.diag([0.0, 0.0, 0.003])
        vals = []
        for shear in (0.0, 1e-12):
            strain = e.copy()
            strain[0, 2] = strain[2, 0] = shear
            vals.append(
                model.evaluate_on_plane(S0r=s, S1r=-s, E0r=strain, E1r=-strain, params=p).damage
            )
        self.assertGreater(vals[0], 0.0)
        self.assertAlmostEqual(vals[0], vals[1], places=10)

    def test_findley_selects_damage_plane(self):
        from multifat.models.cp_methods.fin import Findley

        s = np.zeros((3, 3))
        s[0, 2] = s[2, 0] = 10
        s[2, 2] = 20
        r = Findley().evaluate_on_plane(S0r=s, S1r=-s, E0r=s, E1r=-s, params={"k_FI": 0.3})
        self.assertEqual(r.metric, r.damage)

    def test_actual_harmonic_split_has_no_fictitious_plastic_strain(self):
        s = np.array([[100.0, 50.0, 20.0], [50.0, -40.0, 10.0], [20.0, 10.0, 80.0]])
        c = np.array([[0.0, 25.0, -20.0], [25.0, 80.0, 5.0], [-20.0, 5.0, -10.0]])

        def elastic(t):
            return (1.3 * t - 0.3 * np.trace(t) * np.eye(3)) / 200000

        ge, gp, ne, np_ = split_strain_ranges(s, c, elastic(s), elastic(c), 200000, 0.3)
        self.assertEqual(gp, 0.0)
        self.assertEqual(np_, 0.0)
        from multifat.models.cp_methods.gse import GSE

        r = best_cp_for_case(
            S0=s,
            S1=-s,
            E0=elastic(s),
            E1=-elastic(s),
            R_list=np.array([np.eye(3)]),
            refine=False,
            method=GSE(),
            params={"E": 200000.0, "nu": 0.3},
            S_sin=s,
            S_cos=c,
            E_sin=elastic(s),
            E_cos=elastic(c),
        )
        expected = shear_maximum(0 * s, s, c) * ge / 2 + np.hypot(s[2, 2], c[2, 2]) * ne / 2
        self.assertAlmostEqual(r.ext, expected, places=12)

    def test_shifted_ellipse_maximum(self):
        m, s, c = (np.zeros((3, 3)) for _ in range(3))
        m[0, 2] = 3
        s[0, 2] = 2
        c[1, 2] = 2
        self.assertAlmostEqual(shear_maximum(m, s, c), 5.0)

    def test_retired_models_and_references(self):
        for name in DISABLED_METHODS:
            for variant in (name, name + "_ext"):
                with self.assertRaisesRegex(NotImplementedError, "disabled by numerical audit"):
                    get_models([variant])
        for model in get_models(list(available_cp_methods()) + ["BP", "CAIM"]):
            self.assertIn(model.name, REFERENCES)
            self.assertTrue(model.reference["citation"])


class DamageMetricTests(unittest.TestCase):
    def test_flat_and_shallow_fits(self):
        n = np.array([1e3, 1e4, 1e5, 1e6])
        for b in (0.0, -1e-10, -0.1, 0.01):
            with self.subTest(b=b):
                dp = 20 * n**b
                with np.errstate(all="raise"):
                    fit = fit_power_law_with_survival_std(n, dp, stdnum=0)
                self.assertTrue(np.isfinite(fit["A_surv"]))
                np.testing.assert_allclose(fit["A_surv"] * n ** fit["b"], dp, rtol=1e-10)
        # Life inversion used to overflow for a noisy almost-flat curve.
        dp = np.exp(np.array([0.1, -0.1, -0.1, 0.1]))
        fit = fit_power_law_with_survival_std(n, dp, stdnum=2)
        self.assertEqual(fit["b"], 0.0)
        self.assertTrue(np.isfinite(fit["A_surv"]))
        self.assertAlmostEqual(
            fit["A_surv"], np.exp(-2 * np.std(np.log(dp), ddof=2))
        )  # n - 2 dof (#74)
        self.assertTrue(np.isnan(fit["STD"]))

    def test_unidentifiable_calibration_is_not_reported_as_an_optimum(self):
        with self.assertRaisesRegex(ValueError, "No finite calibration objective"):
            calibrate.minimize_1d_grid_refine(lambda p: np.inf, (0.0, 1.0), n_grid=3, n_refine=1)

    def test_exactly_flat_plot_pipeline(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            results = base / "results"
            plots = base / "plots"
            save_fatigue_csv(
                out_dir=str(results),
                method_name="BP",
                material_name="flat",
                N_exp=[1e3, 1e4, 1e5],
                CP=[8.0, 10.0, 12.0],
                N_pred=[np.nan] * 3,
                basquin_meta={"A_surv": 10.0, "b": 0.0},
                metric_mode="DP_DIFF",
            )
            compare.main(["--results_dir", str(results), "--out_dir", str(plots)])
            self.assertTrue((plots / "flat" / "flat_dp_boxplot_error.png").is_file())
            tex = (plots / "flat" / "flat_dp_error_pdf.tex").read_text()
            self.assertIn("DP_e < DP", tex)
            self.assertIn("DP_e > DP", tex)
            self.assertNotIn("Not-safe", tex)
            stats = pd.read_csv(plots / "summary_error_stats.csv")
            self.assertEqual(stats.iloc[0]["Metric"], "Error_log10_dp")
            self.assertEqual(stats.iloc[0]["n"], 3)

    def test_damage_log_base_sign_and_invalid_inputs(self):
        np.testing.assert_allclose(
            log10_dp_ratio([100, 10, 1, 0, -1, np.nan], [10, 10, 10, 1, 1, 1]),
            [1, 0, -1, np.nan, np.nan, np.nan],
            equal_nan=True,
        )
        self.assertEqual(float(log10_dp_ratio(1e300, 1e-300)), 600.0)

    def test_export_and_statistics_for_all_slopes(self):
        for mode in ("Nf", "DP_DIFF"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as d:
                path = save_fatigue_csv(
                    out_dir=d,
                    method_name="FS",
                    material_name="flat",
                    N_exp=[1e5],
                    CP=[10],
                    N_pred=[np.nan],
                    DP_fit=[100],
                    N_design=[1e4],
                    F_design=[100],
                    N_pred_design=[np.nan],
                    DP_fit_design=[10],
                    metric_mode=mode,
                )
                df = pd.read_csv(path)
                np.testing.assert_allclose(df.Error_log10_dp, [-1, 1])
                self.assertTrue(df.Method_reference.str.contains("doi.org").all())
                stats = compare.compute_stats_table(df)
                row = stats[stats.Metric == "Error_log10_dp"].iloc[0]
                self.assertEqual(row["n"], 2)
                self.assertAlmostEqual(row.Error_mean, 0.0)


if __name__ == "__main__":
    unittest.main()
