"""TP6 - Géométrie du plan."""
import math
import random

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def vecteur(A, B):
    return (B[0] - A[0], B[1] - A[1])


def produit_scalaire(u, v):
    return u[0] * v[0] + u[1] * v[1]


def norme(u):
    return math.hypot(*u)


def angle(u, v):
    c = produit_scalaire(u, v) / (norme(u) * norme(v))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))  # clamp : évite acos(1.0000000002)


def determinant(u, v):
    return u[0] * v[1] - u[1] * v[0]


def distance_point_droite(P, A, B):
    return abs(determinant(vecteur(A, B), vecteur(A, P))) / norme(vecteur(A, B))


assert norme((3, 4)) == 5.0 and produit_scalaire((1, 2), (3, 4)) == 11
assert math.isclose(angle((1, 0), (1, 1)), 45) and math.isclose(angle((1, 0), (0, 1)), 90)
assert math.isclose(distance_point_droite((0, 0), (0, 1), (1, 0)), math.sqrt(2) / 2)
assert distance_point_droite((3, 4), (0, 0), (1, 0)) == 4.0


def orientation(A, B, C):
    d = determinant(vecteur(A, B), vecteur(A, C))
    return (d > 0) - (d < 0)


assert orientation((0, 0), (1, 0), (1, 1)) == 1
assert orientation((0, 0), (1, 0), (1, -1)) == -1
assert orientation((0, 0), (1, 1), (2, 2)) == 0


def aire(polygone):
    n = len(polygone)
    # (i + 1) % n relie le dernier sommet au premier
    s = sum(determinant(polygone[i], polygone[(i + 1) % n]) for i in range(n))
    return abs(s) / 2


MAISON = [(0, 0), (2, 0), (2, 1.5), (1, 2.5), (0, 1.5)]
assert aire([(0, 0), (1, 0), (1, 1), (0, 1)]) == 1.0
assert aire([(0, 0), (4, 0), (0, 3)]) == 6.0
assert aire(MAISON) == 4.0


def se_coupent(A, B, C, D):
    """Segments [AB] et [CD] qui se croisent strictement. Les cas alignés/tangents renvoient False."""
    return (orientation(A, B, C) * orientation(A, B, D) < 0
            and orientation(C, D, A) * orientation(C, D, B) < 0)


assert se_coupent((0, 0), (1, 1), (0, 1), (1, 0))
assert not se_coupent((0, 0), (1, 0), (0, 1), (1, 1))
assert not se_coupent((0, 0), (1, 1), (2, 2), (3, 3))


def intersection(A, B, C, D):
    """Point d'intersection par Cramer, ou None si les segments ne se croisent pas."""
    if not se_coupent(A, B, C, D):
        return None
    u, v, w = vecteur(A, B), vecteur(C, D), vecteur(A, C)
    t = determinant(w, v) / determinant(u, v)
    return (A[0] + t * u[0], A[1] + t * u[1])


assert intersection((0, 0), (2, 2), (0, 2), (2, 0)) == (1.0, 1.0)


def dedans(P, polygone):
    """Lancer de rayon horizontal vers la droite. Inégalité stricte d'un côté (a[1] > y)
    et large de l'autre : un sommet exactement à la hauteur de P n'est compté qu'une fois."""
    x, y = P
    n, inside = len(polygone), False
    for i in range(n):
        a, b = polygone[i], polygone[(i + 1) % n]
        if (a[1] > y) != (b[1] > y):
            if x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]):
                inside = not inside
    return inside


assert dedans((1, 1), MAISON) and dedans((1, 2.4), MAISON)
assert not dedans((3, 1), MAISON) and not dedans((0.2, 2.2), MAISON)

random.seed(0)
N = 200000
hits = sum(dedans((random.uniform(0, 2), random.uniform(0, 2.5)), MAISON) for _ in range(N))
est = hits / N * 2 * 2.5
print("Monte-Carlo aire maison :", est, "(exact 4.0)")  # erreur ~ 1/sqrt(N) : lent

ETOILE = [((2 if i % 2 == 0 else 0.8) * math.cos(math.pi / 2 + i * math.pi / 5),
           (2 if i % 2 == 0 else 0.8) * math.sin(math.pi / 2 + i * math.pi / 5)) for i in range(10)]
assert dedans((0, 0), ETOILE) and dedans((0, 1.5), ETOILE) and not dedans((1.2, 1.2), ETOILE)
print("aire étoile :", aire(ETOILE))

pts = [(random.uniform(-2.5, 2.5), random.uniform(-2.5, 2.5)) for _ in range(60)]
plt.figure()
plt.fill(*zip(*ETOILE), alpha=0.2)
for p in pts:
    plt.plot(*p, "go" if dedans(p, ETOILE) else "ro", ms=4)
plt.axis("equal")
plt.savefig("tp6_etoile.png", dpi=100)
plt.close()


def zone_de(P, zones):
    for i, z in enumerate(zones):
        if dedans(P, z):
            return i
    return None


ZONES = [MAISON, [(3, 0), (4, 0), (4, 1), (3, 1)]]
assert zone_de((1, 1), ZONES) == 0 and zone_de((3.5, 0.5), ZONES) == 1 and zone_de((2.5, 0.5), ZONES) is None


def enveloppe(points):
    """Andrew : tri par abscisse puis chaînes inférieure et supérieure (virages à gauche)."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def chaine(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and orientation(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h

    return chaine(pts)[:-1] + chaine(reversed(pts))[:-1]


assert sorted(enveloppe([(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0.5)])) == [(0, 0), (0, 1), (1, 0), (1, 1)]
gauss = [(random.gauss(0, 1), random.gauss(0, 1)) for _ in range(200)]
env = enveloppe(gauss)
print("sommets de l'enveloppe de 200 points gaussiens :", len(env))
plt.figure()
plt.plot(*zip(*gauss), "k.", ms=3)
plt.plot(*zip(*(env + env[:1])), "r-")
plt.savefig("tp6_enveloppe.png", dpi=100)
plt.close()

print("\n== Réflexion IA : un point-dans-polygone naïf ==")


def dedans_naif(P, poly):
    """Version typique : inégalités larges des deux côtés."""
    x, y = P
    c = False
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        if min(a[1], b[1]) <= y <= max(a[1], b[1]) and a[1] != b[1]:
            if x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]):
                c = not c
    return c


print("(1, 1.5) naïf :", dedans_naif((1, 1.5), MAISON), "| correct :", dedans((1, 1.5), MAISON))
print("(0, 1.5) naïf :", dedans_naif((0, 1.5), MAISON), "| correct :", dedans((0, 1.5), MAISON))
# À la hauteur d'un sommet, le rayon touche deux côtés : avec <= partout, le sommet est
# compté deux fois et la parité est faussée. Les points sur un côté restent indéterminés
# (dedans ou dehors selon le côté) : il faut un test d'appartenance au bord à part.
