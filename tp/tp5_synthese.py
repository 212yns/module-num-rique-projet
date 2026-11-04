"""TP5 - Synthèse : chiffrement affine, brute force, analyse fréquentielle, Vigenère."""
import random
import unicodedata
from collections import Counter
from math import gcd


def lettre_vers_nombre(c):
    return ord(c) - ord("A")


def nombre_vers_lettre(n):
    return chr(n % 26 + ord("A"))


def nettoyer(texte):
    """Majuscules, accents retirés (décomposition NFD), seulement A-Z."""
    sans = unicodedata.normalize("NFD", texte.upper())
    return "".join(c for c in sans if "A" <= c <= "Z")


assert lettre_vers_nombre("C") == 2 and nombre_vers_lettre(2) == "C"
assert nettoyer("Été, à Paris !") == "ETEAPARIS"


def chiffrer(texte, a, b):
    return "".join(nombre_vers_lettre(a * lettre_vers_nombre(c) + b) for c in nettoyer(texte))


assert chiffrer("ATTAQUEALAUBE", 5, 8) == "IZZIKECILIENC"
print("Chiffré (5, 8) :", chiffrer("ATTAQUEALAUBE", 5, 8))
print("César (1, 3)   :", chiffrer("ATTAQUEALAUBE", 1, 3))

# 3.1 Clefs valides : a premier avec 26.
valides = [a for a in range(1, 26) if gcd(a, 26) == 1]
print("\na valides :", valides, len(valides), "| clefs :", len(valides) * 26)
assert len(valides) == 12 and len(valides) * 26 == 312


def inverse_modulaire(a, n):
    """Inverse de a mod n par Euclide étendu ; erreur si pgcd(a, n) != 1."""
    def eucl(a, b):
        if b == 0:
            return a, 1, 0
        g, u, v = eucl(b, a % b)
        return g, v, u - (a // b) * v
    g, u, _ = eucl(a % n, n)
    if g != 1:
        raise ValueError(f"{a} n'est pas inversible modulo {n} (pgcd = {g})")
    return u % n


assert inverse_modulaire(5, 26) == 21
try:
    inverse_modulaire(4, 26)
except ValueError as e:
    print("Erreur attendue :", e)


def dechiffrer(texte, a, b):
    inv = inverse_modulaire(a, 26)  # refuse les clefs invalides
    return "".join(nombre_vers_lettre(inv * (lettre_vers_nombre(c) - b)) for c in texte)


msg = nettoyer("Attaque a l aube")
assert dechiffrer(chiffrer(msg, 5, 8), 5, 8) == msg
try:
    dechiffrer(chiffrer(msg, 4, 7), 4, 7)
except ValueError as e:
    print("Clef (4, 7) refusée :", e)
# Avec a = 4 : x et x + 13 donnent la même lettre chiffrée (4*13 = 52 = 0 mod 26),
# donc deux lettres claires distinctes ont le même chiffré : l'information est perdue.
print("4*0+7 ->", nombre_vers_lettre(7), "| 4*13+7 ->", nombre_vers_lettre(4 * 13 + 7))

PARAGRAPHE = nettoyer(
    "Le chiffrement affine remplace chaque lettre par une autre selon une formule simple. "
    "Cette methode est tres ancienne et elle ne protege plus rien aujourd hui car le nombre "
    "de clefs possibles est tres petit. Une attaque par force brute essaie toutes les clefs "
    "et l analyse des frequences des lettres permet de retrouver le message en quelques "
    "instants, meme sans connaitre la clef utilisee pour chiffrer le texte."
)
SECRET = chiffrer(PARAGRAPHE, 7, 12)


def toutes_les_clefs(secret):
    return [((a, b), dechiffrer(secret, a, b)) for a in valides for b in range(26)]


cands = toutes_les_clefs(SECRET)
print("\nForce brute :", len(cands), "candidats")
lisibles = [(k, t) for k, t in cands if "CHIFFREMENT" in t]
print("clef lisible :", [k for k, _ in lisibles], "|", lisibles[0][1][:40], "...")
assert [k for k, _ in lisibles] == [(7, 12)]

# 6.1 Analyse fréquentielle : la lettre la plus fréquente est l'image de E.
top = Counter(SECRET).most_common(1)[0][0]
y = lettre_vers_nombre(top)
essais = [(a, b) for a in valides for b in range(26) if (a * 4 + b) % 26 == y]
print("\nLettre la plus fréquente :", top, "->", len(essais), "clefs à tester au lieu de 312")
trouve = [(a, b) for a, b in essais if "CHIFFREMENT" in dechiffrer(SECRET, a, b)]
print("clef retrouvée :", trouve)
assert trouve == [(7, 12)]
# Sur « ATTAQUEALAUBE » (13 lettres) : A et T dominent par hasard, rien ne garantit que
# la lettre la plus fréquente soit l'image de E : l'hypothèse statistique n'a pas assez
# de données pour s'appliquer (loi des grands nombres).

# Extension : score = somme(fréquence observée * fréquence française attendue).
FR = dict(zip("ESAINTRULODCPMVQFBGHJXYZKW",
              [14.7, 7.9, 7.6, 7.5, 7.2, 6.9, 6.5, 6.3, 5.5, 5.4, 3.7, 3.4, 3.0, 2.9, 1.6, 1.4,
               1.1, 0.9, 0.9, 0.7, 0.6, 0.4, 0.3, 0.1, 0.05, 0.05]))


def score(t):
    n = len(t)
    return sum(c / n * FR[l] for l, c in Counter(t).items())


meilleur = max(cands, key=lambda kt: score(kt[1]))
print("meilleur score fréquentiel :", meilleur[0])


def indice_coincidence(t):
    N = len(t)
    return sum(n * (n - 1) for n in Counter(t).values()) / (N * (N - 1))


print("\n== Indice de coincidence ==")
random.seed(1)
alea = "".join(random.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for _ in range(len(PARAGRAPHE)))
print("clair", round(indice_coincidence(PARAGRAPHE), 4),
      "| affine", round(indice_coincidence(SECRET), 4),
      "| aléatoire", round(indice_coincidence(alea), 4))
assert abs(indice_coincidence(PARAGRAPHE) - indice_coincidence(SECRET)) < 1e-12
# La substitution ne change que les noms des lettres, pas les effectifs n_i : l'IC est invariant.


def vigenere(texte, cle):
    return "".join(nombre_vers_lettre(lettre_vers_nombre(c) + lettre_vers_nombre(cle[i % len(cle)]))
                   for i, c in enumerate(texte))


V = vigenere(PARAGRAPHE, "CLE")
print("IC Vigenère (clef CLE) :", round(indice_coincidence(V), 4))
print("L : IC moyen des colonnes")
for L in range(1, 8):
    cols = [V[i::L] for i in range(L)]
    print(L, round(sum(indice_coincidence(c) for c in cols) / L, 4))
