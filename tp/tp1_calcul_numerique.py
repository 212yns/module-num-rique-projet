"""TP1 - Calcul numérique. Chaque exercice s'exécute et vérifie son résultat attendu."""
import math
import sys
from fractions import Fraction
from decimal import Decimal, getcontext

print("== Démonstrations ==")
print(2 ** 10, 2 ** 200, len(str(2 ** 200)))
print(0.1 + 0.2, 0.1 + 0.2 == 0.3)
a, b = 0.1 + 0.2, 0.3
print(abs(a - b) < 1e-9, math.isclose(a, b))

print("\n== 1.1 Division euclidienne ==")
q, r = 2026 // 7, 2026 % 7
print(q, r)
assert (q, r) == (289, 3) and 289 * 7 + 3 == 2026

print("\n== 1.2 Dernier chiffre ==")
n = 123456789
print(n % 10, n // 10)
assert (n % 10, n // 10) == (9, 12345678)

print("\n== 2.1 Erreur qui s'accumule ==")
s = 0
for _ in range(1000):
    s += 0.1
print(s, s - 100)
assert s != 100 and abs(s - 100) < 1e-10


def presque_egaux(a, b, tol=1e-9):
    """Égalité à une tolérance absolue près (inadaptée aux très grands nombres)."""
    return abs(a - b) < tol


print("\n== 2.2 presque_egaux ==")
print(presque_egaux(0.1 + 0.2, 0.3), 0.1 + 0.2 == 0.3)
# Pour 1e20, l'écart entre deux flottants voisins vaut ~16384 : la tolérance fixe
# 1e-9 ne reconnaît jamais deux nombres de cet ordre comme égaux.
print(presque_egaux(1e20, 1e20 + 16384), math.isclose(1e20, 1e20 + 16384))

print("\n== 3.1 Table de sept ==")
for i in range(1, 11):
    print(f"7 x {i} = {7 * i}")

print("\n== 3.2 Somme 1..100 ==")
total = 0
for i in range(1, 101):
    total += i
print(total, 100 * 101 // 2)
assert total == 5050 == 100 * 101 // 2


def nombre_de_chiffres(n):
    """Nombre de chiffres de n >= 0, sans conversion en chaîne."""
    c = 1
    while n >= 10:
        n //= 10
        c += 1
    return c


print("\n== 4.1 nombre_de_chiffres ==")
assert nombre_de_chiffres(2026) == 4 and nombre_de_chiffres(7) == 1
print(nombre_de_chiffres(2026), nombre_de_chiffres(7))


def heron(a, tol=1e-15):
    """Racine carrée par Héron ; renvoie (valeur, nombre de tours)."""
    x, tours = a, 0
    while True:
        y = (x + a / x) / 2
        tours += 1
        if abs(y - x) <= tol * y:
            return y, tours
        x = y


print("\n== 5.1 Héron ==")
v, t = heron(2)
print(v, math.sqrt(2), abs(v - math.sqrt(2)), t)
assert abs(v - math.sqrt(2)) < 1e-12 and t < 10
for a in (2, 1e6, 1e30):
    print(a, heron(a)[1], "tours")  # plus a est grand, plus il faut de tours


def vers_binaire(n):
    """Écriture binaire par divisions successives par 2."""
    if n == 0:
        return "0"
    bits = ""
    while n > 0:
        bits = str(n % 2) + bits
        n //= 2
    return bits


print("\n== 6.1 vers_binaire ==")
assert vers_binaire(37) == "100101" and vers_binaire(0) == "0"
print(vers_binaire(37), vers_binaire(0))

print("\n== 6.2 epsilon machine ==")
eps = 1.0
while 1.0 + eps / 2 != 1.0:
    eps /= 2
print(eps, sys.float_info.epsilon)
assert eps == sys.float_info.epsilon == 2.0 ** -52


def racines_naif(a, b, c):
    d = math.sqrt(b * b - 4 * a * c)
    return (-b + d) / (2 * a), (-b - d) / (2 * a)


print("\n== 7.1 Formule naïve sur x^2 - 1e8 x + 1 ==")
petite, grande = racines_naif(1, -1e8, 1)
print(petite, grande, petite * grande)  # petite ~ 7.45e-9 au lieu de 1e-8


def racines(a, b, c):
    """Évite la soustraction de nombres proches : x1 sans annulation, x2 = c/(a x1)."""
    d = math.sqrt(b * b - 4 * a * c)
    x1 = (-b - d) / (2 * a) if b > 0 else (-b + d) / (2 * a)
    return x1, c / (a * x1)


print("\n== 7.2 Formule stable ==")
x1, x2 = racines(1, -1e8, 1)
print(x1, x2, x1 * x2)
for x in (x1, x2):
    res = x * x - 1e8 * x + 1
    print("résidu relatif :", abs(res) / (x * x + 1e8 * x + 1))
assert math.isclose(x2, 1e-8) and math.isclose(x1, 1e8)

print("\n== 7.3 Calcul exact ==")
s = Fraction(0)
for _ in range(1000):
    s += Fraction(1, 10)
print(s, Fraction(1, 10) + Fraction(2, 10) == Fraction(3, 10))
assert s == 100

getcontext().prec = 50
print(Decimal(1) / Decimal(7))
