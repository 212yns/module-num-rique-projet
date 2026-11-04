"""TP4 - Matrices et systèmes linéaires."""
import math
import random
import time
from fractions import Fraction

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def chrono(f, *args):
    t = time.perf_counter()
    r = f(*args)
    return r, time.perf_counter() - t


def transposee(A):
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


print("== 1.1 Transposée ==")
assert transposee([[1, 2, 3], [4, 5, 6]]) == [[1, 4], [2, 5], [3, 6]]


def produit(A, B):
    if len(A[0]) != len(B):
        raise ValueError("dimensions incompatibles")
    # [[0]*q for _ in range(n)] : une liste par ligne. [[0]*q]*n répéterait la MÊME ligne n fois.
    C = [[0] * len(B[0]) for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] += A[i][k] * B[k][j]
    return C


print("\n== 2.1 Produit ==")
I2 = [[1, 0], [0, 1]]
assert produit([[1, 2], [3, 4]], I2) == [[1, 2], [3, 4]]
print(produit([[1, 2], [3, 4]], [[0, 1], [1, 0]]))  # [[2, 1], [4, 3]]
try:
    produit([[1, 2]], [[1, 2]])
except ValueError as e:
    print("refusé :", e)

print("\n== 2.2 Coût du produit ==")
random.seed(0)
prev = None
for n in (40, 80, 120):
    A = [[random.random() for _ in range(n)] for _ in range(n)]
    _, t = chrono(produit, A, A)
    print(n, f"{t:.3f}s", f"x{t / prev:.1f}" if prev else "")
    prev = t
# O(n^3) : doubler n multiplie le temps par ~8.

print("\n== 3.1 numpy ==")
n = 120
A = [[random.random() for _ in range(n)] for _ in range(n)]
_, t1 = chrono(produit, A, A)
An = np.array(A)
_, t2 = chrono(lambda: An @ An)
print(f"python {t1:.4f}s, numpy {t2:.6f}s, rapport {t1 / t2:.0f}")


def resoudre(A, b):
    """Gauss avec pivot partiel ; copies locales, arguments inchangés."""
    n = len(A)
    M = [list(map(float, A[i])) + [float(b[i])] for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))  # plus grand pivot en valeur absolue
        if M[p][k] == 0:
            raise ZeroDivisionError("système singulier")
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


print("\n== 4.1 Gauss ==")
A = [[1, 2, 2], [1, 3, -2], [3, 5, 8]]
b = [2, -1, 8]
x = resoudre(A, b)
print(x)
assert all(math.isclose(a, c, abs_tol=1e-12) for a, c in zip(x, [3, -1, 0.5]))
assert all(math.isclose(a, c) for a, c in zip(np.array(A) @ np.array(x), b))

print("\n== 4.2 Comparaison à numpy.linalg.solve (n = 50, graine 7) ==")
rng = np.random.default_rng(7)
An, bn = rng.random((50, 50)), rng.random(50)
mine = resoudre(An.tolist(), bn.tolist())
print(np.max(np.abs(np.array(mine) - np.linalg.solve(An, bn))))
assert np.allclose(mine, np.linalg.solve(An, bn), atol=1e-9)

print("\n== 5.1 Hilbert : solve vs inv ==")
for n in (10, 12):
    H = np.array([[1 / (i + j + 1) for j in range(n)] for i in range(n)])
    bh = H @ np.ones(n)
    e1 = np.max(np.abs(np.linalg.solve(H, bh) - 1))
    e2 = np.max(np.abs(np.linalg.inv(H) @ bh - 1))
    print(f"n={n} cond={np.linalg.cond(H):.2e} erreur solve={e1:.2e} erreur inv@b={e2:.2e}")
# inv calcule n systèmes (toutes les colonnes de H^-1), chacun mal conditionné, puis on re-multiplie :
# plus de calcul, plus d'erreurs accumulées. solve résout directement un seul système.

print("\n== 6.1 Rotations ==")


def R(deg):
    t = math.radians(deg)
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]])


assert np.allclose(R(60) @ R(30), R(90))
maison = np.array([(0, 0), (2, 0), (2, 1.5), (1, 2.5), (0, 1.5), (0, 0)]).T
plt.figure()
for deg in (0, 30, 60, 90):
    p = R(deg) @ maison
    plt.plot(p[0], p[1], label=f"{deg}°")
plt.axis("equal")
plt.legend()
plt.savefig("tp4_maison.png", dpi=100)
plt.close()
# Homothétie 2 puis symétrie / rotation : l'ordre compte (produit matriciel non commutatif).
S = np.array([[1, 0], [0, -1]])
print("symétrie∘rotation == rotation∘symétrie ?", np.allclose(S @ R(30), R(30) @ S))


def determinant(A):
    n = len(A)
    M = [list(map(float, r)) for r in A]
    d = 1.0
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        if M[p][k] == 0:
            return 0.0
        if p != k:
            M[k], M[p] = M[p], M[k]
            d = -d  # chaque échange de lignes change le signe
        d *= M[k][k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n):
                M[i][j] -= f * M[k][j]
    return d


print("\n== 7.1 Déterminant ==")
print(determinant([[1, 2, 2], [1, 3, -2], [3, 5, 8]]), np.linalg.det(np.array([[1, 2, 2], [1, 3, -2], [3, 5, 8]])))
assert math.isclose(determinant([[1, 2, 2], [1, 3, -2], [3, 5, 8]]), -2)
assert determinant([[1, 2], [2, 4]]) == 0

print("\n== 8.1 Hilbert exact ==")
n = 12
H = [[Fraction(1, i + j + 1) for j in range(n)] for i in range(n)]
bf = [sum(row) for row in H]


def resoudre_exact(A, b):
    """Même algorithme sur des Fractions. La version flottante convertissait en float()
    (d'où la perte) ; et `M[p][k] == 0` testait un flottant arrondi. Ici : exact."""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for k in range(n):
        p = next(i for i in range(k, n) if M[i][k] != 0)  # en exact, tout pivot non nul convient
        M[k], M[p] = M[p], M[k]
        for i in range(k + 1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= f * M[k][j]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


sol = resoudre_exact(H, bf)
print(sol[:3], "...")
assert all(v == Fraction(1, 1) for v in sol)
# `resoudre` (flottants) convertit les Fractions avec float() : sur Hilbert 12 elle
# ne plante pas ici mais renvoie une solution fausse (la perte vient du conditionnement).
# Un test de singularité du type `abs(pivot) < tol` (ou un pivot qui tombe à 0 selon
# l'implémentation) est ce qui déclenche l'erreur « système singulier » de l'énoncé.
print("erreur de resoudre en flottants sur Hilbert 12 :",
      max(abs(v - 1) for v in resoudre(H, bf)))


def inverse(A):
    """Gauss-Jordan sur [A | I]."""
    n = len(A)
    M = [list(map(float, A[i])) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k]))
        M[k], M[p] = M[p], M[k]
        piv = M[k][k]
        M[k] = [v / piv for v in M[k]]
        for i in range(n):
            if i != k:
                f = M[i][k]
                M[i] = [a - f * c for a, c in zip(M[i], M[k])]
    return [r[n:] for r in M]


Ainv = inverse([[1, 2, 2], [1, 3, -2], [3, 5, 8]])
print("inverse, 1re ligne :", [round(v, 6) for v in Ainv[0]])
assert np.allclose(Ainv[0], [-17, 3, 5])
assert np.allclose(np.array(Ainv) @ np.array([[1, 2, 2], [1, 3, -2], [3, 5, 8]]), np.eye(3))


def translation(tx, ty):
    return np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]], dtype=float)


def rotation_h(deg):
    t = math.radians(deg)
    return np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1]])


print("\n== 9.1 Coordonnées homogènes ==")
P = np.array([1, 0, 1])
a = rotation_h(90) @ translation(3, 1) @ P  # d'abord la translation, puis la rotation
b = translation(3, 1) @ rotation_h(90) @ P  # d'abord la rotation, puis la translation
print(np.round(a[:2], 6), np.round(b[:2], 6))
assert np.allclose(a[:2], [-1, 4]) and np.allclose(b[:2], [3, 2])
