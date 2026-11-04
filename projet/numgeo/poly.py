"""Polynômes (coefficients par degré décroissant) et recherche de racines."""


def horner(p, x):
    r = 0.0
    for c in p:
        r = r * x + c
    return r


def derivative(p):
    n = len(p) - 1
    return [c * (n - i) for i, c in enumerate(p[:-1])] or [0.0]


def strip(p):
    """Retire les zéros de tête exacts."""
    i = 0
    while i < len(p) - 1 and p[i] == 0:
        i += 1
    return list(p[i:])


def mul(p, q):
    r = [0.0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return r


def from_roots(roots):
    p = [1.0]
    for r in roots:
        p = mul(p, [1.0, -r])
    return p


def bisection(p, a, b, tol=1e-14, maxit=200):
    fa, fb = horner(p, a), horner(p, b)
    if fa == 0:
        return a
    if fb == 0:
        return b
    if fa * fb > 0:
        raise ValueError("pas de changement de signe")
    for _ in range(maxit):
        m = (a + b) / 2
        fm = horner(p, m)
        if fm == 0 or (b - a) / 2 < tol:
            return m
        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return (a + b) / 2


def newton(p, x0, tol=1e-14, maxit=100):
    """Retourne (racine, nb d'itérations)."""
    dp = derivative(p)
    x = x0
    for k in range(maxit):
        d = horner(dp, x)
        if d == 0:
            raise ZeroDivisionError("dérivée nulle")
        dx = horner(p, x) / d
        x -= dx
        if abs(dx) < tol * max(1.0, abs(x)):
            return x, k + 1
    return x, maxit


def real_roots(p, a, b, tol=1e-14):
    """Racines réelles de p dans [a,b] : on borne les racines de p' puis on
    bissecte sur chaque intervalle où p est monotone (donc <= 1 racine)."""
    p = strip(p)
    if len(p) == 1:
        return []
    if len(p) == 2:
        r = -p[1] / p[0]
        return [r] if a <= r <= b else []
    crit = real_roots(derivative(p), a, b, tol)
    pts = [a] + crit + [b]
    out = []
    for lo, hi in zip(pts, pts[1:]):
        flo, fhi = horner(p, lo), horner(p, hi)
        if flo == 0:
            out.append(lo)
        elif flo * fhi < 0:
            out.append(bisection(p, lo, hi, tol))
    if horner(p, pts[-1]) == 0:
        out.append(pts[-1])
    dedup = []
    for r in sorted(out):
        if not dedup or abs(r - dedup[-1]) > 1e-12:
            dedup.append(r)
    return dedup
