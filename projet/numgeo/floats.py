"""Flottants : erreurs d'arrondi, sommation, équation du second degré."""
import math
from fractions import Fraction


def machine_eps():
    """Plus petit eps (puissance de 2) tel que 1 + eps > 1, calculé."""
    eps = 1.0
    while 1.0 + eps / 2 > 1.0:
        eps /= 2
    return eps


def naive_sum(xs):
    s = 0.0
    for x in xs:
        s += x
    return s


def kahan_sum(xs):
    """Sommation compensée de Kahan : erreur O(eps) au lieu de O(n*eps)."""
    s = c = 0.0
    for x in xs:
        y = x - c
        t = s + y
        c = (t - s) - y
        s = t
    return s


def exact_sum(xs):
    """Référence exacte : chaque float est un rationnel dyadique."""
    return sum((Fraction(x) for x in xs), Fraction(0))


def quad_roots_naive(a, b, c):
    """Racines réelles de ax^2+bx+c par la formule du lycée (annulation possible)."""
    d = b * b - 4 * a * c
    if d < 0:
        return ()
    r = math.sqrt(d)
    return tuple(sorted(((-b - r) / (2 * a), (-b + r) / (2 * a))))


def quad_roots_stable(a, b, c):
    """Évite l'annulation : q = -(b + sign(b) sqrt(d))/2, racines q/a et c/q."""
    d = b * b - 4 * a * c
    if d < 0:
        return ()
    q = -(b + math.copysign(math.sqrt(d), b)) / 2
    if q == 0:
        return (0.0, 0.0)
    return tuple(sorted((q / a, c / q)))
