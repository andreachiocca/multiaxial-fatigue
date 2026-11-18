"""Findley damage-parameter checks, without alternative-amplitude life prediction.

Ported from upstream tests/test_findley_harmonic.py (A. Chiocca, 176ce53, 50dabaf); imports only
changed. Not ported: the output_quantity assertion and the formulation-fingerprint test (the
versioned calibration store is not used here); the latter is replaced by an isolation test.
"""

import unittest

import numpy as np

from multifat.eval.evaluator import best_cp_for_case
from multifat.models.cp_methods.fin_harmonic import evaluate_findley
from multifat.models.registry import get_models


IDENTITY = np.eye(3)
NORMAL_X = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
ALGORITHMS = ("mcc_harmonic", "mce_harmonic", "rss")


def prescribed_case(minor=10.0):
    """On normal z, shear semiaxes are 10/minor and sigma_max = 7 + 5."""
    sine = np.array([[0.0, 0.0, 10.0], [0.0, 0.0, 0.0], [10.0, 0.0, 3.0]])
    cosine = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, minor], [0.0, minor, 4.0]])
    mean = np.diag([0.0, 0.0, 7.0])
    return dict(
        S0=mean + sine,
        S1=mean - sine,
        E0=np.zeros((3, 3)),
        E1=np.zeros((3, 3)),
        params={"k_FI": 0.2},
        R_list=np.array([IDENTITY]),
        harmonics=dict(S_sin=sine, S_cos=cosine),
    )


def ranking_case(proportional=False):
    """Two planes with analytically known opposite MCC/MCE rankings."""
    sine = np.array([[0.0, 8.0, 10.0], [8.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    cosine = (
        0.5 * sine
        if proportional
        else np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 10.0], [0.0, 10.0, 0.0]])
    )
    mean = np.diag([5.0, 0.0, 7.0])
    return dict(
        S0=mean + sine,
        S1=mean - sine,
        E0=sine / 200000,
        E1=-sine / 200000,
        R_list=np.array([IDENTITY, NORMAL_X]),
        params=dict(k_FI=0.2, k_FS=0.5, Sy=250.0, Sigm1=240.0, Taum1=145.0, Su=600.0),
        harmonics=dict(S_sin=sine, S_cos=cosine, E_sin=sine / 200000, E_cos=cosine / 200000),
    )


def supplied_planes_only(name, case):
    """Production criterion on exactly the supplied planes (grid-only, no refinement).

    The production adapter also refines the maximum between planes (F2), so it can exceed
    the supplied-plane value; the upstream snapshots and the harmonic variant are
    supplied-plane values.
    """
    method = get_models([name])[0]._m
    harmonics = case.get("harmonics") or {}
    best = best_cp_for_case(
        S0=case["S0"],
        E0=case["E0"],
        S1=case["S1"],
        E1=case["E1"],
        R_list=case["R_list"],
        method=method,
        params=case["params"],
        S_sin=harmonics.get("S_sin"),
        S_cos=harmonics.get("S_cos"),
        E_sin=harmonics.get("E_sin"),
        E_cos=harmonics.get("E_cos"),
        refine=False,
    )
    return {"": float(best.by_metric), "_ext": float(best.ext)}


class FindleyHarmonicTests(unittest.TestCase):
    def test_default_matches_unmodified_main_snapshots_and_adapter(self):
        # Captured on upstream main 71eaac3 (supplied planes only) before the diagnostic entry point.
        for proportional, expected in ((True, 15.317821063276353), (False, 15.542135623730951)):
            case = ranking_case(proportional)
            self.assertEqual(supplied_planes_only("FIN", case), {"": expected, "_ext": expected})
            actual = evaluate_findley(**case)
            self.assertEqual(actual, get_models(["FIN"])[0].evaluate_case(**case))
            # The production adapter refines between planes: never below the supplied planes.
            for variant in ("", "_ext"):
                self.assertGreaterEqual(actual.values[variant], expected)
            explicit = evaluate_findley(**case, shear_amplitude_algorithm="mce_harmonic")
            for variant in ("", "_ext"):
                np.testing.assert_allclose(explicit.values[variant], expected, rtol=2e-15)

    def test_endpoint_default_and_explicit_selection_without_harmonics(self):
        case = prescribed_case()
        case.pop("harmonics")
        self.assertEqual(evaluate_findley(**case), get_models(["FIN"])[0].evaluate_case(**case))
        self.assertAlmostEqual(supplied_planes_only("FIN", case)[""], 12.0)
        self.assertGreaterEqual(evaluate_findley(**case).values[""], 12.0)
        for algorithm in ALGORITHMS:
            with self.assertRaisesRegex(ValueError, "S_sin and S_cos"):
                evaluate_findley(**case, shear_amplitude_algorithm=algorithm)

    def test_circular_and_elliptical_complete_findley_expression(self):
        for minor in (10.0, 5.0, 0.0):
            case = prescribed_case(minor)
            for algorithm, amplitude in (
                ("mcc_harmonic", 10.0),
                ("mce_harmonic", np.sqrt(100 + minor**2)),
            ):
                result = evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
                for variant in ("", "_ext"):
                    self.assertAlmostEqual(
                        result.values[variant], amplitude + 0.2 * 12.0, places=13
                    )

    def test_proportional_same_plane_and_full_search(self):
        case = ranking_case(True)
        for planes in (case["R_list"], case["R_list"][:1], case["R_list"][1:]):
            one = dict(case, R_list=planes)
            mcc = evaluate_findley(**one, shear_amplitude_algorithm="mcc_harmonic")
            mce = evaluate_findley(**one, shear_amplitude_algorithm="mce_harmonic")
            for variant in ("", "_ext"):
                self.assertAlmostEqual(mcc.values[variant], mce.values[variant], places=13)

    def test_algorithm_changes_plane_ranking_and_search(self):
        case = ranking_case()
        # Plane z is circular, plane x is a line of length amplitude sqrt(164).
        expected = {
            "mcc_harmonic": [10.0 + 1.4, np.sqrt(164.0) + 1.0],
            "mce_harmonic": [np.sqrt(200.0) + 1.4, np.sqrt(164.0) + 1.0],
        }
        winners = []
        for algorithm, values in expected.items():
            actual = [
                evaluate_findley(
                    **dict(case, R_list=np.array([r])), shear_amplitude_algorithm=algorithm
                ).values[""]
                for r in case["R_list"]
            ]
            np.testing.assert_allclose(actual, values, rtol=2e-15)
            winners.append(int(np.argmax(actual)))
            for planes in (case["R_list"], case["R_list"][::-1]):
                result = evaluate_findley(
                    **dict(case, R_list=planes), shear_amplitude_algorithm=algorithm
                )
                for variant in ("", "_ext"):
                    self.assertAlmostEqual(result.values[variant], max(values), places=13)
        self.assertEqual(winners, [1, 0])

    def test_mean_normal_and_signed_compression_are_not_algorithm_dependent(self):
        for mean_normal in (7.0, -20.0):
            case = prescribed_case()
            for key in ("S0", "S1"):
                case[key][2, 2] += mean_normal - 7.0
            for algorithm, amplitude in (("mcc_harmonic", 10.0), ("mce_harmonic", np.sqrt(200.0))):
                result = evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
                self.assertAlmostEqual(
                    (result.values[""] - amplitude) / 0.2, mean_normal + 5.0, places=12
                )

    def test_mean_shear_translation_does_not_change_amplitude(self):
        case = prescribed_case()
        for algorithm in ALGORITHMS:
            expected = evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
            shifted = dict(case, S0=case["S0"].copy(), S1=case["S1"].copy())
            for key in ("S0", "S1"):
                for i, value in ((0, 30.0), (1, -40.0)):
                    shifted[key][i, 2] += value
                    shifted[key][2, i] += value
            self.assertEqual(
                evaluate_findley(**shifted, shear_amplitude_algorithm=algorithm), expected
            )

    def test_rss_alias_and_tiny_harmonics_without_endpoint_reconstruction(self):
        case = prescribed_case()
        self.assertEqual(
            evaluate_findley(**case, shear_amplitude_algorithm="rss"),
            evaluate_findley(**case, shear_amplitude_algorithm="mce_harmonic"),
        )
        case["params"] = {"k_FI": 0.0}
        for scale in (0.0, 1e-200, 1e-16, 1e150):
            scaled = dict(
                case, harmonics={key: value * scale for key, value in case["harmonics"].items()}
            )
            for algorithm, amplitude in (("mcc_harmonic", 10.0), ("mce_harmonic", np.sqrt(200.0))):
                value = evaluate_findley(**scaled, shear_amplitude_algorithm=algorithm).values[""]
                np.testing.assert_allclose(value, amplitude * scale, rtol=2e-15, atol=0.0)

    def test_phase_and_tangent_basis_invariance(self):
        case = prescribed_case(5.0)
        angle = 0.63
        a, b = (case["harmonics"][key] for key in ("S_sin", "S_cos"))
        shifted = dict(
            case,
            harmonics=dict(
                S_sin=a * np.cos(angle) - b * np.sin(angle),
                S_cos=a * np.sin(angle) + b * np.cos(angle),
            ),
        )
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        rotated = dict(case, R_list=np.array([rotation]))
        for algorithm in ALGORITHMS:
            expected = evaluate_findley(**case, shear_amplitude_algorithm=algorithm).values[""]
            for changed in (shifted, rotated):
                value = evaluate_findley(**changed, shear_amplitude_algorithm=algorithm).values[""]
                self.assertAlmostEqual(value, expected, places=13)

    def test_random_planes_against_independent_svd_and_normal_projection(self):
        rng = np.random.default_rng(403)
        mean, sine, cosine = (rng.normal(size=(3, 3)) for _ in range(3))
        mean, sine, cosine = (t + t.T for t in (mean, sine, cosine))
        planes = np.array([np.linalg.qr(rng.normal(size=(3, 3)))[0] for _ in range(30)])
        case = dict(
            S0=mean + sine,
            S1=mean - sine,
            E0=np.zeros((3, 3)),
            E1=np.zeros((3, 3)),
            params={"k_FI": 0.3},
            R_list=planes,
            harmonics=dict(S_sin=sine, S_cos=cosine),
        )
        for algorithm in ("mcc_harmonic", "mce_harmonic"):
            scores = []
            for r in planes:
                n = r[:, 2]
                matrix = r[:, :2].T @ np.column_stack((sine @ n, cosine @ n))
                sv = np.linalg.svd(matrix, compute_uv=False)
                amplitude = sv[0] if algorithm == "mcc_harmonic" else np.hypot(*sv)
                normal = n @ mean @ n + np.hypot(n @ sine @ n, n @ cosine @ n)
                score = amplitude + 0.3 * normal
                scores.append(score)
                one = evaluate_findley(
                    **dict(case, R_list=np.array([r])), shear_amplitude_algorithm=algorithm
                )
                np.testing.assert_allclose(one.values[""], score, rtol=3e-15)
            result = evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
            np.testing.assert_allclose(result.values[""], max(scores), rtol=3e-15)
            self.assertEqual(result.values[""], result.values["_ext"])

    def test_invalid_algorithm_missing_coefficients_and_bad_inputs(self):
        case = prescribed_case()
        for algorithm in ("method_default", "mcc", "", 1, []):
            with self.assertRaisesRegex(ValueError, "Unknown shear amplitude algorithm"):
                evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
        for h in ({}, {"S_sin": np.eye(3)}, {"S_cos": np.eye(3)}):
            with self.assertRaisesRegex(ValueError, "S_sin and S_cos"):
                evaluate_findley(
                    **dict(case, harmonics=h), shear_amplitude_algorithm="mcc_harmonic"
                )
        for key in ("S_sin", "S_cos"):
            for bad in (np.zeros((2, 2)), np.full((3, 3), np.nan), np.full((3, 3), np.inf)):
                h = dict(case["harmonics"], **{key: bad})
                with self.assertRaises(ValueError):
                    evaluate_findley(
                        **dict(case, harmonics=h), shear_amplitude_algorithm="mcc_harmonic"
                    )
            h = dict(case["harmonics"], **{key: np.eye(3, dtype=complex)})
            with self.assertRaises(TypeError):
                evaluate_findley(
                    **dict(case, harmonics=h), shear_amplitude_algorithm="mcc_harmonic"
                )
        for key in ("S0", "S1"):
            with self.assertRaises(ValueError):
                evaluate_findley(
                    **dict(case, **{key: np.full((3, 3), np.nan)}),
                    shear_amplitude_algorithm="mcc_harmonic",
                )
        if np.finfo(np.longdouble).tiny < np.finfo(float).tiny:
            tiny = np.eye(3, dtype=np.longdouble) * np.longdouble("1e-400")
            h = dict(case["harmonics"], S_sin=tiny)
            with self.assertRaisesRegex(ValueError, "float64 conversion"):
                evaluate_findley(
                    **dict(case, harmonics=h), shear_amplitude_algorithm="mcc_harmonic"
                )
        for planes in ([], None, np.zeros((1, 3, 3))):
            with self.assertRaises(ValueError):
                evaluate_findley(
                    **dict(case, R_list=planes), shear_amplitude_algorithm="mcc_harmonic"
                )
        for k in (-0.1, np.inf, np.nan):
            with self.assertRaises(ValueError):
                evaluate_findley(
                    **dict(case, params={"k_FI": k}), shear_amplitude_algorithm="mcc_harmonic"
                )

    def test_other_methods_unchanged_and_do_not_accept_algorithm(self):
        # Upstream 71eaac3 gave FS 0.00014340125522463185 for the non-proportional case; here FS
        # uses the longest-chord shear range (F3), hence the different value. MATAKE and MCD use the
        # circumscribed circle (longest chord) since 17 Nov 2026 (maintainer decision; upstream and
        # RSS gave 15.600468957064283 and 14.987968957064284 out of phase).
        snapshots = {
            True: {
                "FS": 0.00014460999273909115,
                "MATAKE": 15.35948772994302,
                "MCD": 14.92198772994302,
            },
            False: {
                "FS": 0.00012934310959614355,
                "MATAKE": 13.847915141532363,
                "MCD": 13.410415141532363,
            },
        }
        for proportional, values in snapshots.items():
            case = ranking_case(proportional)
            for algorithm in ALGORITHMS:
                evaluate_findley(**case, shear_amplitude_algorithm=algorithm)
                for name, value in values.items():
                    model = get_models([name])[0]
                    self.assertEqual(supplied_planes_only(name, case), {"": value, "_ext": value})
                    model.evaluate_case(**case)
                    with self.assertRaises(TypeError):
                        model.evaluate_case(**case, shear_amplitude_algorithm=algorithm)
        with self.assertRaises(TypeError):
            get_models(["FIN"])[0].evaluate_case(
                **ranking_case(), shear_amplitude_algorithm="mcc_harmonic"
            )

    def test_production_code_does_not_import_the_harmonic_variant(self):
        # The registry, adapter and FIN module must stay independent of the opt-in evaluator.
        import importlib
        import inspect

        for name in (
            "multifat.models.registry",
            "multifat.models.cp_adapter",
            "multifat.models.cp_methods.fin",
        ):
            source = inspect.getsource(importlib.import_module(name))
            self.assertNotIn("fin_harmonic", source)
            self.assertNotIn("findley_calibration", source)


if __name__ == "__main__":
    unittest.main()
