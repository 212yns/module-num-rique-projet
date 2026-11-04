"""TP3 - Polynômes et fonctions. Convention : coefficients du terme constant au plus haut degré."""
import math
import time

import matplotlib

matplotlib.use("Agg")  # on écrit des PNG au lieu d'ouvrir une fenêtre
import matplotlib.pyplot as plt


def chrono(f, *args):
    t = time.perf_counter()
    r = f(*args)
    return r, time.perf_counter() - t


def evaluer_naif(P, x):
    return sum(c * x ** i for i, c in enumerate(P))


def evaluer(P, x):
    """Horner : a0 + x(a1 + x(a2 + ...)), une multiplication et une addition par coefficient."""
    r = 0
    for c in reversed(P):
        r = r * x + c
    return r


print("== 1.1 / 2.1 ==")
assert evaluer_naif([2, 3, 0, 5], 2) == 48 == evaluer([2, 3, 0, 5], 2)
assert evaluer_naif([2, 3, 0, 5], 0) == 2
for x in (-2, -0.5, 0, 1.5, 3):
    assert math.isclose(evaluer_naif([2, 3, 0, 5], x), evaluer([2, 3, 0, 5], x))

print("\n== 2.2 Temps sur flottants (x = 1.0001) ==")
for d in (1000, 5000, 20000):
    P = [1] * (d + 1)
    a, t1 = chrono(evaluer_naif, P, 1.0001)
    b, t2 = chrono(evaluer, P, 1.0001)
    print(d, f"naïf {t1:.5f}s Horner {t2:.5f}s rapport {t1 / t2:.1f}  écart {abs(a - b):.2e}")

print("\n== 2.3 Temps sur entiers (x = 3) ==")
for d in (200, 400, 800):
    P = [1] * (d + 1)
    a, t1 = chrono(evaluer_naif, P, 3)
    b, t2 = chrono(evaluer, P, 3)
    assert a == b
    print(d, f"rapport {t1 / t2:.1f}")


def deriver(P):
    return [i * c for i, c in enumerate(P)][1:]


print("\n== 3.1 deriver ==")
assert deriver([2, 3, 0, 5]) == [3, 0, 15] and deriver([7]) == []
assert evaluer(deriver([7]), 5) == 0


def dichotomie(f, a, b, tolerance=1e-12):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("f(a) et f(b) doivent être de signes opposés")
    while b - a > tolerance:
        m = (a + b) / 2
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return (a + b) / 2


print("\n== 4.1 Dichotomie ==")
r = dichotomie(lambda x: x * x - 2, 0, 2)
print(r)
assert abs(r - math.sqrt(2)) < 1e-11
try:
    dichotomie(lambda x: x * x + 1, 0, 2)
except ValueError as e:
    print("refusé :", e)


def newton(P, depart, tolerance=1e-12, max_iter=50, trace=False):
    """Retourne (x, itérations) ; garde-fous : tangente horizontale et plafond d'itérations."""
    dP, x = deriver(P), depart
    for k in range(1, max_iter + 1):
        d = evaluer(dP, x)
        if d == 0:
            raise ZeroDivisionError("tangente horizontale")
        nx = x - evaluer(P, x) / d
        if trace:
            print("  x =", nx)
        if abs(nx - x) <= tolerance * max(1, abs(nx)):  # critère relatif
            return nx, k
        x = nx
    return x, max_iter


print("\n== 4.2 Newton sur x^3 - 2x - 5 ==")
P = [-5, -2, 0, 1]
r, k = newton(P, 2)
print(r, k)
assert abs(r - 2.0945514815) < 1e-9 and k <= 5

print("\n== 4.3 Newton qui échoue : x^3 - 2x + 2 depuis 0 ==")
r, k = newton([2, -2, 0, 1], 0, max_iter=8, trace=True)
print("pas de convergence après", k, "itérations (cycle 0 -> 1 -> 0)")

print("\n== 4.4 Coût du point de départ ==")
for x0 in (2, 1, 0):
    print(x0, newton(P, x0))

print("\n== 5.1 Tracé ==")
xs = [-3 + 6 * i / 600 for i in range(601)]
plt.figure()
plt.plot(xs, [evaluer(P, x) for x in xs])
plt.axhline(0, color="k", lw=0.5)
plt.ylim(-20, 20)
r, _ = newton(P, 2)
plt.plot([r], [0], "ro")
# itérés de Newton depuis 2
x = 2.0
its = [x]
for _ in range(4):
    x = x - evaluer(P, x) / evaluer(deriver(P), x)
    its.append(x)
plt.plot(its, [evaluer(P, v) for v in its], "g.")
plt.title("x^3 - 2x - 5 : une seule racine réelle (~2.0945)")
plt.savefig("tp3_trace.png", dpi=100)
plt.close()
print("tp3_trace.png écrit : un seul passage par l'axe ; les deux autres racines sont complexes conjuguées")


def somme(P, Q):
    n = max(len(P), len(Q))
    P, Q = P + [0] * (n - len(P)), Q + [0] * (n - len(Q))
    return [a + b for a, b in zip(P, Q)]


def produit_poly(P, Q):
    R = [0] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q):
            R[i + j] += a * b
    return R


print("\n== 6.1 / 6.2 Algèbre ==")
assert produit_poly([1, 1], [1, 1]) == [1, 2, 1]
assert produit_poly([-1, 1], [1, 1]) == [-1, 0, 1]
assert somme([1, 2], [0, 0, 3]) == [1, 2, 3]
A, B = [1, 2, 3], [4, 0, -1, 2]
assert evaluer(produit_poly(A, B), 1.7) == evaluer(A, 1.7) * evaluer(B, 1.7) or math.isclose(
    evaluer(produit_poly(A, B), 1.7), evaluer(A, 1.7) * evaluer(B, 1.7))
pascal = [1]
for _ in range(10):
    pascal = produit_poly(pascal, [1, 1])
print(pascal, sum(pascal))
assert pascal == [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1] and sum(pascal) == 1024


def derivee_avant(f, x, h):
    return (f(x + h) - f(x)) / h


def derivee_centree(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)


print("\n== 7.1 Meilleur pas ==")
print(f"{'h':>8} {'avant':>12} {'centrée':>12}")
for e in (1, 2, 4, 6, 8, 10, 12):
    h = 10.0 ** -e
    print(f"1e-{e:<5} {abs(derivee_avant(math.sin, 1, h) - math.cos(1)):>12.2e} "
          f"{abs(derivee_centree(math.sin, 1, h) - math.cos(1)):>12.2e}")
# Erreur = troncature (~h ou ~h^2) + arrondi (~eps/h) : la somme est minimale pour un h intermédiaire.


def rectangles(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + i * h) for i in range(n))


def trapezes(f, a, b, n):
    h = (b - a) / n
    return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))


def simpson(f, a, b, n):
    n += n % 2  # n pair
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4 * sum(f(a + i * h) for i in range(1, n, 2))
    s += 2 * sum(f(a + i * h) for i in range(2, n, 2))
    return s * h / 3


print("\n== 8.1 Intégration de x^2 sur [0,1] (exact 1/3) ==")
for n in (10, 100, 1000):
    print(n, *(f"{abs(m(lambda x: x * x, 0, 1, n) - 1 / 3):.2e}" for m in (rectangles, trapezes, simpson)))
print("sin sur [0, pi] avec Simpson :")
for n in (10, 100, 1000):
    print(n, f"{abs(simpson(math.sin, 0, math.pi, n) - 2):.2e}")

print("\n== 8.2 Calcul de pi ==")
g = lambda x: 4 / (1 + x * x)
print(abs(trapezes(g, 0, 1, 1000) - math.pi), abs(simpson(g, 0, 1, 1000) - math.pi))


def newton_general(f, depart, tolerance=1e-12, max_iter=50):
    x = depart
    for k in range(1, max_iter + 1):
        d = derivee_centree(f, x, 1e-6)
        if d == 0:
            raise ZeroDivisionError("tangente horizontale")
        nx = x - f(x) / d
        if abs(nx - x) < tolerance * max(1, abs(nx)):
            return nx, k
        x = nx
    return x, max_iter


print("\n== Newton général sur cos x - x ==")
print(newton_general(lambda x: math.cos(x) - x, 1.0))
# Départ où il diverge : f'(x) = -sin x - 1 s'annule en x = -pi/2 (tangente horizontale).
try:
    print(newton_general(lambda x: math.cos(x) - x, -math.pi / 2))
except ZeroDivisionError as e:
    print("départ -pi/2 :", e)
