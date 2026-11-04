import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from numgeo import floats, geom, linalg, poly


class TestFloats(unittest.TestCase):
    def test_eps(self):
        self.assertEqual(floats.machine_eps(), 2.0 ** -52)

    def test_kahan_beats_naive(self):
        xs = [0.1] * 1_000_000
        exact = float(floats.exact_sum(xs))
        self.assertLess(abs(floats.kahan_sum(xs) - exact),
                        abs(floats.naive_sum(xs) - exact))

    def test_stable_quadratic(self):
        # racines 1e-8 et 1e8 : la formule naïve perd la petite racine
        a, b, c = 1.0, -(1e8 + 1e-8), 1.0
        small = floats.quad_roots_stable(a, b, c)[0]
        self.assertAlmostEqual(small / 1e-8, 1.0, places=9)
        naive = floats.quad_roots_naive(a, b, c)[0]
        self.assertGreater(abs(naive / 1e-8 - 1.0), abs(small / 1e-8 - 1.0))


class TestPoly(unittest.TestCase):
    def test_horner_derivative(self):
        p = [1.0, -3.0, 2.0]  # x^2 - 3x + 2
        self.assertEqual(poly.horner(p, 2.0), 0.0)
        self.assertEqual(poly.derivative(p), [2.0, -3.0])

    def test_newton_sqrt2(self):
        r, it = poly.newton([1.0, 0.0, -2.0], 1.0)
        self.assertAlmostEqual(r, math.sqrt(2), places=14)
        self.assertLess(it, 8)

    def test_real_roots(self):
        p = poly.from_roots([0.1, 0.5, 0.9])
        r = poly.real_roots(p, 0.0, 1.0)
        for got, want in zip(r, [0.1, 0.5, 0.9]):
            self.assertAlmostEqual(got, want, places=10)

    def test_real_roots_none(self):
        self.assertEqual(poly.real_roots([1.0, 0.0, 1.0], -5, 5), [])


class TestLinalg(unittest.TestCase):
    def test_solve(self):
        A = [[2, 1, 1], [1, 3, 2], [1, 0, 0]]
        x = linalg.solve(A, [4, 5, 6])
        for got, want in zip(linalg.matvec(A, x), [4, 5, 6]):
            self.assertAlmostEqual(got, want, places=12)

    def test_pivot_needed(self):
        A = [[1e-20, 1.0], [1.0, 1.0]]
        good = linalg.solve(A, [1.0, 2.0])
        bad = linalg.solve(A, [1.0, 2.0], pivot=False)
        self.assertAlmostEqual(good[0], 1.0, places=12)
        self.assertGreater(abs(bad[0] - 1.0), 0.1)

    def test_det_and_singular(self):
        self.assertAlmostEqual(linalg.det([[1, 2], [3, 4]]), -2.0)
        with self.assertRaises(ZeroDivisionError):
            linalg.lu([[1, 2], [2, 4]])

    def test_cond_hilbert_grows(self):
        self.assertGreater(linalg.cond1(linalg.hilbert(8)), 1e9)

    def test_lstsq_exact_line(self):
        xs = [0, 1, 2, 3]
        A = linalg.vandermonde(xs, 1)
        y = [1.0, 3.0, 5.0, 7.0]  # y = 2x + 1
        for fn in (linalg.lstsq_qr, linalg.lstsq_normal):
            a, b = fn(A, y)
            self.assertAlmostEqual(a, 2.0, places=10)
            self.assertAlmostEqual(b, 1.0, places=10)


class TestGeom(unittest.TestCase):
    def test_hull_square(self):
        pts = [(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0.5), (0.5, 0)]
        h = geom.convex_hull(pts)
        self.assertEqual(sorted(h), [(0, 0), (0, 1), (1, 0), (1, 1)])
        self.assertAlmostEqual(geom.polygon_area(h), 1.0)

    def test_orient_exact(self):
        # trois points quasi alignés : le filtre flottant doit basculer en exact
        self.assertEqual(geom.orient((0, 0), (1, 1), (2, 2)), 0)
        self.assertEqual(geom.orient((0, 0), (1, 0), (0, 1)), 1)
        self.assertEqual(geom.orient((0.1, 0.1), (0.3, 0.3), (0.5, 0.5)),
                         geom.orient((0.1, 0.1), (0.3, 0.3), (0.5, 0.5)))

    def test_point_in_polygon(self):
        sq = [(0, 0), (2, 0), (2, 2), (0, 2)]
        self.assertTrue(geom.point_in_polygon((1, 1), sq))
        self.assertTrue(geom.point_in_polygon((2, 1), sq))
        self.assertFalse(geom.point_in_polygon((3, 1), sq))

    def test_bezier_line(self):
        ctrl = [(0, 0), (1, 2), (2, -2), (3, 0)]
        ts = geom.bezier_line_intersections(ctrl, (0, 0), (3, 0))
        self.assertEqual(len(ts), 3)  # t = 0, 0.5, 1
        for t, want in zip(ts, [0.0, 0.5, 1.0]):
            self.assertAlmostEqual(t, want, places=8)
            self.assertAlmostEqual(geom.bezier_point(ctrl, t)[1], 0.0, places=8)


if __name__ == "__main__":
    unittest.main()
