"""Expériences du projet. Produit le rapport `resultats.md`.

Lancer : python experiences.py
Fil rouge : à chaque étage (flottants -> polynômes -> matrices -> géométrie),
on mesure où les flottants trahissent l'algorithme « naïf » et ce qui le répare.
"""
import math
import os
import random
import sys

sys.stdout.reconfigure(encoding="utf-8")

from numgeo import floats, geom, linalg, poly

OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def sci(x):
    return f"{x:.3e}"


def exp_sommes():
    say("## 1. Sommation : naïve vs Kahan vs exacte")
    say("| n | erreur naïve | erreur Kahan |")
    say("|---|---|---|")
    for n in (10**3, 10**4, 10**5, 10**6):
        xs = [0.1] * n
        ex = floats.exact_sum(xs)
        say(f"| {n} | {sci(abs(float(floats.naive_sum(xs) - ex)))} "
            f"| {sci(abs(float(floats.kahan_sum(xs) - ex)))} |")
    say()


def exp_quadratique():
    say("## 2. Équation du second degré x² - (1e k + 1e -k)x + 1 : petite racine 1e-k")
    say("| k | erreur relative naïve | erreur relative stable |")
    say("|---|---|---|")
    for k in (2, 4, 6, 8, 10, 12):
        b = -(10.0**k + 10.0**-k)
        want = 10.0**-k
        n = floats.quad_roots_naive(1.0, b, 1.0)[0]
        s = floats.quad_roots_stable(1.0, b, 1.0)[0]
        say(f"| {k} | {sci(abs(n - want) / want)} | {sci(abs(s - want) / want)} |")
    say()


def exp_wilkinson():
    say("## 3. Polynôme de Wilkinson : racines 1..n, coefficients arrondis")
    say("Les coefficients de ∏(x-i) sont calculés en flottants, puis on cherche "
        "les racines (bissection sur [0, n+1]). Écart max à la vraie racine :")
    say("| n | écart max |")
    say("|---|---|")
    for n in (5, 10, 13, 15, 17, 20):
        p = poly.from_roots([float(i) for i in range(1, n + 1)])
        roots = poly.real_roots(p, 0.0, n + 1.0)
        if len(roots) == n:
            err = max(abs(r - i) for r, i in zip(roots, range(1, n + 1)))
            say(f"| {n} | {sci(err)} |")
        else:
            say(f"| {n} | seulement {len(roots)}/{n} racines réelles retrouvées |")
    say()


def exp_newton():
    say("## 4. Newton : convergence quadratique vs bissection")
    p = [1.0, 0.0, -2.0]
    r, it = poly.newton(p, 1.0)
    say(f"- Newton sur x²-2 depuis 1.0 : {it} itérations, erreur "
        f"{sci(abs(r - math.sqrt(2)))}")
    # racine double : Newton devient linéaire
    p2 = poly.from_roots([1.0, 1.0, 3.0])
    r2, it2 = poly.newton(p2, 0.0, tol=1e-12)
    say(f"- Newton sur (x-1)²(x-3) depuis 0.0 (racine double) : {it2} itérations "
        f"pour 1e-12 (convergence linéaire, facteur 1/2)")
    say()


def exp_lu():
    say("## 5. LU avec / sans pivot ; conditionnement de Hilbert")
    A = [[1e-20, 1.0], [1.0, 1.0]]
    say(f"- Système [[1e-20,1],[1,1]] x = [1,2] (solution ≈ [1,1]) : "
        f"avec pivot {linalg.solve(A, [1.0, 2.0])}, "
        f"sans pivot {linalg.solve(A, [1.0, 2.0], pivot=False)}")
    say()
    say("| n | cond₁(Hilbert) | erreur relative sur x (b = H·1) |")
    say("|---|---|---|")
    for n in (3, 6, 9, 12):
        H = linalg.hilbert(n)
        b = linalg.matvec(H, [1.0] * n)
        x = linalg.solve(H, b)
        err = max(abs(v - 1.0) for v in x)
        say(f"| {n} | {sci(linalg.cond1(H))} | {sci(err)} |")
    say()


def exp_moindres_carres():
    say("## 6. Ajustement polynomial : équations normales vs QR")
    say("Données y = polynôme de degré d (coefficients 1), 30 points sur [0,1], "
        "bruit nul. Erreur max sur les coefficients retrouvés :")
    say("| degré | équations normales | QR (Gram-Schmidt modifié) |")
    say("|---|---|---|")
    xs = [i / 29 for i in range(30)]
    for d in (3, 5, 7, 9):
        A = linalg.vandermonde(xs, d)
        y = [poly.horner([1.0] * (d + 1), x) for x in xs]
        res = []
        for fn in (linalg.lstsq_normal, linalg.lstsq_qr):
            try:
                c = fn(A, y)
                res.append(sci(max(abs(v - 1.0) for v in c)))
            except ZeroDivisionError:
                res.append("échec (matrice singulière en flottants)")
        say(f"| {d} | {res[0]} | {res[1]} |")
    say()


def exp_geometrie():
    say("## 7. Géométrie : orientation exacte et enveloppe convexe")
    random.seed(0)
    # test classique (Kettner et al.) : p décrit une grille 256x256 de pas 2^-53
    # autour de (0.5, 0.5), q=(12,12) et r=(24,24) : tous quasi colinéaires.
    u = 2.0**-53
    a, b = (12.0, 12.0), (24.0, 24.0)
    wrong = total = 0
    for i in range(256):
        for j in range(256):
            p = (0.5 + i * u, 0.5 + j * u)
            raw = geom.cross(p, a, b)
            wrong += geom.orient(p, a, b) != (raw > 0) - (raw < 0)
            total += 1
    say(f"- Test d'orientation flottant pur vs exact sur {total} triplets "
        f"quasi colinéaires : {wrong} signes faux ({100 * wrong / total:.1f} %)")
    pts = [(random.random(), random.random()) for _ in range(5000)]
    h = geom.convex_hull(pts)
    say(f"- Enveloppe convexe de 5000 points aléatoires dans [0,1]² : "
        f"{len(h)} sommets, aire {geom.polygon_area(h):.6f}")
    say()


def exp_bezier():
    say("## 8. Application : intersection courbe de Bézier cubique / droite")
    ctrl = [(0, 0), (1, 3), (2, -3), (3, 0)]
    ts = geom.bezier_line_intersections(ctrl, (0, 0), (3, 0))
    say("Courbe de points de contrôle (0,0),(1,3),(2,-3),(3,0) coupée par l'axe Ox :")
    for t in ts:
        x, y = geom.bezier_point(ctrl, t)
        say(f"- t = {t:.12f} -> point ({x:.9f}, {y:.2e})")
    say()


if __name__ == "__main__":
    say("# Résultats des expériences")
    say()
    for f in (exp_sommes, exp_quadratique, exp_wilkinson, exp_newton, exp_lu,
              exp_moindres_carres, exp_geometrie, exp_bezier):
        f()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultats.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
