"""Unit tests for gaia-lensing-dr4 (standard library only): python -m unittest -v test_lensfit

MV2, Munim, Inc. Apache-2.0. Made with AI agents working in MV2's name.
"""
import math
import random
import unittest

import events
import lensfit as L

ALPHA, DELTA = 217.4, -62.7
GEOM = dict(alpha=ALPHA, delta=DELTA, lens_offset=(1500.0, 300.0), lens_pm=(-600.0, -100.0),
            lens_plx=200.0)
SRC = (0.0, 0.0, 0.5, -2.0, 3.0)


def synth(theta2, n=90, noise=0.3, seed=1, geom=GEOM, src=SRC):
    rng = random.Random(seed)
    pts = []
    for i in range(n):
        t = -2.9 + 5.4 * (i + rng.random()) / n
        psi = rng.uniform(0.0, 2.0 * math.pi)
        p = dict(t=t, s=math.sin(psi), c=math.cos(psi), fr=0.0, sig=noise)
        p["f"] = L.model_parallax_factor(geom["alpha"], geom["delta"], p)
        D = L.separation_track([p], geom, src)[0]
        r2 = D[0] ** 2 + D[1] ** 2
        w = (src[0] + src[3] * t) * p["s"] + (src[1] + src[4] * t) * p["c"] + src[2] * p["f"]
        w += theta2 * (D[0] * p["s"] + D[1] * p["c"]) / r2
        p["w"] = w + rng.gauss(0.0, noise)
        pts.append(p)
    return pts


class TestGeometry(unittest.TestCase):
    def test_equinox(self):
        # 2017 Mar 20 12:00 TT (JD 2457833.0): Sun near lambda 0, Gaia near -x.
        t = (2457833.0 - L.JD_J2000) / 365.25 + 2000.0 - 2017.5
        x, y, z = L.observer_xyz(t)
        self.assertAlmostEqual(x, -1.01 * 0.996, delta=0.02)
        self.assertLess(abs(y), 0.03)
        self.assertLess(abs(z), 0.02)

    def test_solstice(self):
        # 2017 Jun 21 (JD 2457926.0): Sun at lambda 90, so Gaia at (0, -cos eps, -sin eps) * 1.026
        t = (2457926.0 - L.JD_J2000) / 365.25 + 2000.0 - 2017.5
        x, y, z = L.observer_xyz(t)
        self.assertLess(abs(x), 0.03)
        self.assertAlmostEqual(y, -1.026 * math.cos(math.radians(23.44)), delta=0.03)
        self.assertAlmostEqual(z, -1.026 * math.sin(math.radians(23.44)), delta=0.03)

    def test_parallax_direction(self):
        # observer on +x, star at alpha=90 deg on the equator: shift toward -x = +alpha
        pa, pd = L.parallax_offsets(90.0, 0.0, (1.0, 0.0, 0.0))
        self.assertAlmostEqual(pa, 1.0)
        self.assertAlmostEqual(pd, 0.0)
        # observer on +z, star on the equator: shift toward -z = south
        pa, pd = L.parallax_offsets(0.0, 0.0, (0.0, 0.0, 1.0))
        self.assertAlmostEqual(pd, -1.0)

    def test_locate(self):
        pts = synth(0.0, n=40)
        a, d, rms = L.locate_from_parallax_factors(pts, coarse_deg=3.0, fine_deg=0.2, span_deg=4.0)
        self.assertLess(rms, 0.01)
        self.assertLess(abs(d - DELTA), 1.0)
        self.assertLess(abs(((a - ALPHA + 180.0) % 360.0) - 180.0) * math.cos(math.radians(DELTA)), 1.0)

    def test_geometry_residuals(self):
        pts = synth(0.0, n=20)
        self.assertLess(L.median(L.geometry_residuals(ALPHA, DELTA, pts)), 1e-9)
        self.assertGreater(L.median(L.geometry_residuals(ALPHA + 90.0, -DELTA, pts)), L.GEOMETRY_TOL)

    def test_tangent_offset(self):
        da, dd = L.tangent_offset_mas(10.0 + 1.0 / 3600.0, 60.0, 10.0, 60.0)
        self.assertAlmostEqual(da, 500.0, delta=0.5)
        self.assertAlmostEqual(dd, 0.0)
        da, _ = L.tangent_offset_mas(0.0001, 0.0, 359.9999, 0.0)
        self.assertAlmostEqual(da, 720.0, delta=0.01)

    def test_closest_approach(self):
        t, r = L.closest_approach(GEOM, SRC)
        # straight-line part: lens offset (1500,300) + (-600,-100) t - source pm; min near t=2.5
        # (parallax of 199.5 mas against a 607 mas/yr relative motion shifts it by < 0.35 yr)
        self.assertAlmostEqual(t, 2.52, delta=0.45)
        self.assertLess(r, 260.0)


class TestFit(unittest.TestCase):
    def test_recovers_injected_deflection(self):
        theta2 = L.KAPPA * 0.3 * (200.0 - 0.5)
        fr = L.fit_event(synth(theta2, noise=0.1), GEOM, (0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertEqual(fr["status"], "OK")
        self.assertLess(abs(fr["thetaE2"] - theta2), 4.0 * fr["thetaE2_err"])
        self.assertGreater(fr["snr"], 3.0)
        self.assertEqual(L.event_verdict(fr), "DEFLECTION_DETECTED")
        M, sM = L.mass_from_fit(fr, 200.0, 0.05)
        self.assertLess(abs(M - 0.3), 4.0 * sM)
        self.assertGreaterEqual(fr["sep_min_used_mas"], L.R_CUT_MAS)

    def test_null(self):
        fr = L.fit_event(synth(0.0, seed=7), GEOM, (0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertLess(abs(fr["snr"]), 4.0)
        self.assertIn(L.event_verdict(fr), ("NOT_DETECTED", "DEFLECTION_DETECTED", "NEGATIVE_DEFLECTION"))
        ul = L.upper_limit_mass(fr, 200.0)
        self.assertGreater(ul, 0.0)

    def test_negative(self):
        theta2 = -L.KAPPA * 0.5 * 200.0
        fr = L.fit_event(synth(theta2, noise=0.1, seed=3), GEOM, (0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertEqual(L.event_verdict(fr), "NEGATIVE_DEFLECTION")

    def test_source_parameters(self):
        fr = L.fit_event(synth(0.0, noise=0.1, seed=11), GEOM, (0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertAlmostEqual(fr["params"]["parallax_mas"][0], 0.5, delta=5 * fr["params"]["parallax_mas"][1])
        self.assertAlmostEqual(fr["params"]["pmra_masyr"][0], -2.0, delta=5 * fr["params"]["pmra_masyr"][1])

    def test_too_few(self):
        pts = synth(0.0, n=10)
        self.assertEqual(L.fit_event(pts, GEOM, SRC)["status"], "TOO_FEW")
        self.assertEqual(L.event_verdict(dict(status="TOO_FEW")), "TOO_FEW")

    def test_cut_removes_close_points(self):
        theta2 = L.KAPPA * 0.3 * 199.5
        pts = synth(theta2)
        D = L.separation_track(pts, GEOM, SRC)
        n_close = sum(1 for d in D if math.hypot(*d) < L.R_CUT_MAS)
        self.assertGreater(n_close, 0)
        fr = L.fit_event(pts, GEOM, (0.0, 0.0, 0.0, 0.0, 0.0))
        self.assertLessEqual(fr["n_points"], len(pts) - n_close)
        fr2 = L.fit_event(pts, GEOM, (0.0, 0.0, 0.0, 0.0, 0.0), r_cut=L.R_CUT_ALT_MAS)
        self.assertLess(fr2["n_points"], fr["n_points"])

    def test_wls_exact(self):
        cols = [[1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 2.0, 3.0]]
        y = [1.0, 3.0, 5.0, 7.0]
        f = L.solve_wls(cols, y, [1.0] * 4)
        self.assertAlmostEqual(f["beta"][0], 1.0)
        self.assertAlmostEqual(f["beta"][1], 2.0)
        self.assertAlmostEqual(f["chi2"], 0.0)

    def test_invert(self):
        A = [[4.0, 1.0], [1.0, 3.0]]
        B = L.invert(A)
        self.assertAlmostEqual(B[0][0], 3.0 / 11.0)
        self.assertAlmostEqual(B[0][1], -1.0 / 11.0)


class TestGrading(unittest.TestCase):
    def test_mass_propagation(self):
        fr = dict(thetaE2=100.0, thetaE2_err=10.0, plx_S=0.0, plx_S_err=0.0, cov_thetaE2_plx=0.0)
        M, sM = L.mass_from_fit(fr, 100.0, 0.0)
        self.assertAlmostEqual(M, 100.0 / (L.KAPPA * 100.0))
        self.assertAlmostEqual(sM, 10.0 / (L.KAPPA * 100.0))
        self.assertEqual(L.mass_from_fit(fr, -1.0, 0.0), (None, None))

    def test_predicted(self):
        self.assertAlmostEqual(L.predicted_thetaE2(1.0, 100.0, 0.0), 814.4)

    def test_validation(self):
        self.assertEqual(L.validation_check(0.6, 0.1, 0.56, 0.08), "PASS")
        self.assertEqual(L.validation_check(1.5, 0.1, 0.56, 0.08), "FAIL")

    def test_overall(self):
        self.assertEqual(L.overall(False, [], []), "GEOMETRY_FAIL")
        self.assertEqual(L.overall(True, ["PASS"], ["DEFLECTION_DETECTED", "NOT_DETECTED"]),
                         "VALIDATED_1_HEADLINE_DETECTED")
        self.assertEqual(L.overall(True, ["PASS", "FAIL"], []), "SYSTEMATICS_0_HEADLINE_DETECTED")
        self.assertEqual(L.overall(True, [], ["TOO_FEW"]), "UNVALIDATED_0_HEADLINE_DETECTED")


class TestJobHelpers(unittest.TestCase):
    def test_propagate(self):
        import lens_job
        ra, dec = lens_job.propagate(10.0, 60.0, 3600.0, 1800.0, 2015.5, 2017.5)
        self.assertAlmostEqual(dec - 60.0, 3600.0 / 3.6e6)
        self.assertAlmostEqual((ra - 10.0) * 0.5, 7200.0 / 3.6e6, places=9)
        self.assertEqual(lens_job.propagate(1.0, 2.0, None, None, 2016.0, 2017.5), (1.0, 2.0))

    def test_num(self):
        import lens_job
        self.assertIsNone(lens_job.num(""))
        self.assertIsNone(lens_job.num("nan"))
        self.assertEqual(lens_job.num("1.5"), 1.5)

    def test_rehearsal_grade(self):
        import rehearsal
        g = rehearsal.grade([0.01] * 12, [0.5] * 9 + [5.0] * 3)
        self.assertEqual(g["outcome"], "REHEARSAL_PASS")
        g = rehearsal.grade([0.01] * 11 + [0.2], [0.5] * 12)
        self.assertEqual(g["R1_geometry"], "FAIL")
        g = rehearsal.grade([0.01] * 12, [0.5] * 8 + [5.0] * 4)
        self.assertEqual(g["R2_null"], "FAIL")


class TestEvents(unittest.TestCase):
    def test_table(self):
        ev = events.as_dicts()
        self.assertEqual(len(ev), 43)
        self.assertEqual([e["num"] for e in ev], list(range(1, 44)))
        self.assertTrue(all(e["t_ca"] < 2020.5 for e in ev))
        self.assertTrue(all(isinstance(e["lens_dr2"], int) and isinstance(e["source_dr2"], int) for e in ev))
        nums = {e["num"] for e in ev}
        self.assertTrue(set(events.HEADLINE) <= nums)
        self.assertTrue(set(events.VALIDATION) <= nums)
        # headline rule: sM5 / M_in <= 0.20 (plus LAWD 37)
        by = {e["num"]: e for e in ev}
        for n in events.HEADLINE:
            if n != 37:
                self.assertLessEqual(by[n]["sm5"] / by[n]["m_in"], 0.20 + 1e-9)
        rule = {e["num"] for e in ev if e["sm5"] and e["sm5"] / e["m_in"] <= 0.20 + 1e-9}
        self.assertEqual(rule | {37}, set(events.HEADLINE))


if __name__ == "__main__":
    unittest.main()
