"""Géométrie plane : orientation robuste, enveloppe convexe, aire, Bézier."""
from fractions import Fraction
from math import comb

from . import poly


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def orient(o, a, b):
    """Signe exact de cross(o,a,b). Le calcul flottant sert de filtre ; on
    bascule en Fraction (exact) quand le résultat est dans le bruit d'arrondi."""
    v = cross(o, a, b)
    scale = abs(a[0] - o[0]) * abs(b[1] - o[1]) + abs(a[1] - o[1]) * abs(b[0] - o[0])
    if abs(v) > 1e-14 * scale:
        return (v > 0) - (v < 0)
    e = cross(*[tuple(map(Fraction, p)) for p in (o, a, b)])
    return (e > 0) - (e < 0)


def convex_hull(points):
    """Monotone chain (Andrew), O(n log n), sans points colinéaires."""
    pts = sorted(set(map(tuple, points)))
    if len(pts) <= 2:
        return pts

    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and orient(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h

    return half(pts)[:-1] + half(reversed(pts))[:-1]


def polygon_area(pts):
    """Aire signée (formule du lacet), positive en sens trigonométrique."""
    n = len(pts)
    return sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1]
               for i in range(n)) / 2


def point_in_polygon(pt, pts):
    """Lancer de rayon ; True si à l'intérieur ou sur le bord."""
    x, y = pt
    n, inside = len(pts), False
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        if (orient(a, b, pt) == 0 and min(a[0], b[0]) <= x <= max(a[0], b[0])
                and min(a[1], b[1]) <= y <= max(a[1], b[1])):
            return True
        if (a[1] > y) != (b[1] > y):
            if x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]):
                inside = not inside
    return inside


def bezier_point(ctrl, t):
    """De Casteljau (stable numériquement)."""
    pts = [tuple(p) for p in ctrl]
    while len(pts) > 1:
        pts = [((1 - t) * p[0] + t * q[0], (1 - t) * p[1] + t * q[1])
               for p, q in zip(pts, pts[1:])]
    return pts[0]


def _bernstein_to_power(vals):
    """Valeurs de contrôle (base de Bernstein) -> coefficients de puissance,
    degré décroissant."""
    n = len(vals) - 1
    coef = [comb(n, k) * sum((-1) ** (k - i) * comb(k, i) * vals[i]
                             for i in range(k + 1)) for k in range(n + 1)]
    return coef[::-1]


def bezier_line_intersections(ctrl, p, q):
    """Paramètres t dans [0,1] où la courbe coupe la droite (pq) : on injecte
    la courbe dans l'équation de la droite, ce qui donne un polynôme en t."""
    nx, ny = -(q[1] - p[1]), q[0] - p[0]
    vals = [nx * (c[0] - p[0]) + ny * (c[1] - p[1]) for c in ctrl]
    return poly.real_roots(_bernstein_to_power(vals), 0.0, 1.0)
