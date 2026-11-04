"""TP2 - Arithmétique et nombres premiers."""
import time


def chrono(f, *args):
    t = time.perf_counter()
    r = f(*args)
    return r, time.perf_counter() - t


def diviseurs(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def est_parfait(n):
    return n > 1 and sum(diviseurs(n)[:-1]) == n


print("== 1.1 Diviseurs ==")
assert diviseurs(28) == [1, 2, 4, 7, 14, 28] and diviseurs(97) == [1, 97]
print(diviseurs(28), diviseurs(97))
print("parfaits < 10000 :", [n for n in range(1, 10000) if est_parfait(n)])


def pgcd_naif(a, b):
    for d in range(min(a, b), 0, -1):
        if a % d == 0 and b % d == 0:
            return d


def pgcd(a, b):
    """Euclide : pgcd(a, b) = pgcd(b, a % b) ; le reste décroît strictement, donc terminaison."""
    while b:
        a, b = b, a % b
    return a


print("\n== 2.1 PGCD ==")
assert pgcd(1071, 462) == pgcd_naif(1071, 462) == 21
print(pgcd(1071, 462), pgcd_naif(1071, 462))


def est_premier_naif(n):
    if n < 2:
        return False
    for d in range(2, n):
        if n % d == 0:
            return False
    return True


def est_premier(n):
    """Teste les diviseurs jusqu'à racine de n : un diviseur > racine a un cofacteur < racine."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


print("\n== 3.1 Tableau des temps ==")
print(f"{'n':>10} {'naïf (s)':>12} {'racine (s)':>12}")
for n in (10007, 100003, 1000003):
    assert est_premier_naif(n) and est_premier(n)
    _, t1 = chrono(est_premier_naif, n)
    _, t2 = chrono(est_premier, n)
    print(f"{n:>10} {t1:>12.6f} {t2:>12.6f}")
# naïf : O(n) -> x10 à chaque ligne ; racine : O(sqrt n) -> x sqrt(10) ~ 3.

print("\n== Encadré IA : cas limites d'un test de primalité ==")
print("n=1 :", est_premier(1), "| n=2 :", est_premier(2), "| n=-7 :", est_premier(-7))
assert not est_premier(1) and est_premier(2) and not est_premier(-7)


def crible(n):
    """Ératosthène : on barre les multiples de p à partir de p*p (les plus petits
    multiples 2p, 3p... ont déjà été barrés par un facteur premier plus petit)."""
    est = [True] * (n + 1)
    est[0] = est[1] = False
    p = 2
    while p * p <= n:
        if est[p]:
            for m in range(p * p, n + 1, p):
                est[m] = False
        p += 1
    return [i for i, v in enumerate(est) if v]


print("\n== 4.1 Crible ==")
pr, t1 = chrono(crible, 100000)
cnt, t2 = chrono(lambda: sum(1 for i in range(100000) if est_premier(i)))
print(len(pr), f"{t1:.4f}s", "| un par un :", cnt, f"{t2:.4f}s")
assert len(pr) == cnt == 9592


def facteurs(n):
    f, d = [], 2
    while d * d <= n:
        while n % d == 0:
            f.append(d)
            n //= d
        d += 1
    if n > 1:  # reste final : un facteur premier > racine (sans lui, facteurs(14) = [2])
        f.append(n)
    return f


print("\n== 5.1 Facteurs ==")
assert facteurs(360) == [2, 2, 2, 3, 3, 5] and facteurs(14) == [2, 7]
print(facteurs(360), facteurs(14))
p, q = 999983, 999979  # deux premiers de six chiffres
f, t = chrono(facteurs, p * q)
print(f, f"{t:.3f}s")  # ~ sqrt(10^12) = 10^6 essais ; pour 300 chiffres : 10^150 essais


def fibonacci_recursif(n, compteur=None):
    if compteur is not None:
        compteur[n] = compteur.get(n, 0) + 1
    return n if n < 2 else fibonacci_recursif(n - 1, compteur) + fibonacci_recursif(n - 2, compteur)


print("\n== 6.1 Fibonacci récursif ==")
for n in (20, 25, 30):
    v, t = chrono(fibonacci_recursif, n)
    print(n, v, f"{t:.4f}s")

print("\n== 6.2 Compteur d'appels ==")
cpt = {}
fibonacci_recursif(30, cpt)
print("F(1) calculé", cpt[1], "fois ; appels totaux", sum(cpt.values()))
assert cpt[1] == 832040 and sum(cpt.values()) == 2692537


def fibonacci(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


print("\n== 6.3 Fibonacci itératif ==")
assert all(fibonacci(i) == fibonacci_recursif(i) for i in range(26))
v, t = chrono(fibonacci, 10000)
print(len(str(v)), "chiffres", f"{t:.5f}s")
assert len(str(v)) == 2090


def temps_de_vol(n, plafond=10 ** 6):
    c = 0
    while n != 1 and c < plafond:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        c += 1
    return c


print("\n== Syracuse ==")
best = max(range(1, 10000), key=temps_de_vol)
print(temps_de_vol(27), best, temps_de_vol(best))
assert temps_de_vol(27) == 111 and best == 6171 and temps_de_vol(6171) == 261


def euclide_etendu(a, b):
    """Renvoie (g, u, v) avec a*u + b*v = g. Récursion adaptée : une seule branche
    d'appel (profondeur ~ log b), donc pas d'explosion comme pour Fibonacci."""
    if b == 0:
        return a, 1, 0
    g, u, v = euclide_etendu(b, a % b)
    return g, v, u - (a // b) * v


print("\n== 7.1 Euclide étendu ==")
assert euclide_etendu(1071, 462) == (21, -3, 7)
assert euclide_etendu(5, 26) == (1, -5, 1) and (-5) % 26 == 21
print(euclide_etendu(1071, 462), euclide_etendu(5, 26))


def puissance_mod_naive(a, e, n):
    r = 1
    for _ in range(e):
        r = r * a % n
    return r


def puissance_mod(a, e, n):
    """Exponentiation rapide : e = bits ; on élève au carré à chaque bit. O(log e)."""
    r, a = 1, a % n
    while e:
        if e & 1:
            r = r * a % n
        a = a * a % n
        e >>= 1
    return r


print("\n== 7.2 Puissance modulaire ==")
r1, t1 = chrono(puissance_mod_naive, 7, 10 ** 7, 1000003)
r2, t2 = chrono(puissance_mod, 7, 10 ** 7, 1000003)
assert r1 == r2 == pow(7, 10 ** 7, 1000003)
print(r1, f"naïf {t1:.3f}s, rapide {t2 * 1e6:.1f} µs")

print("\n== Fermat : menteurs en base 2 sous 10000 ==")
print(pow(2, 340, 341))
menteurs = [n for n in range(3, 10000) if pow(2, n - 1, n) == 1 and not est_premier(n)]
print(len(menteurs), menteurs)
