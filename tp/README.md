# TP

Un script par TP, avec les exercices dans l'ordre de l'énoncé (étapes d'entraînement
comprises). Chaque script vérifie par `assert` les résultats attendus de l'énoncé et
affiche ses mesures.

| Fichier | Sujet |
|---|---|
| `tp1_calcul_numerique.py` | entiers, flottants, Héron, binaire, formule stable, `Fraction` |
| `tp2_arithmetique.py` | PGCD, primalité, crible, Fibonacci, Euclide étendu, puissance modulaire |
| `tp3_polynomes_fonctions.py` | Horner, dichotomie, Newton, dérivation et intégration numériques |
| `tp4_matrices_gauss.py` | produit, Gauss, Hilbert, déterminant, calcul exact, coordonnées homogènes |
| `tp5_synthese.py` | chiffrement affine, force brute, fréquences, indice de coïncidence, Vigenère |
| `tp6_geometrie_plan.py` | produit scalaire, déterminant, segments, point dans polygone, enveloppe convexe |
| `ia.md` | notes sur les encadrés « L'IA vous le donne » |

Ce que chaque script affiche est enregistré dans `sorties/` (un fichier `.txt` par TP).

Lancer un script : `python tp1_calcul_numerique.py`. Les TP3, TP4 et TP6 utilisent
`numpy` et `matplotlib` (`pip install numpy matplotlib`) et écrivent leurs figures en PNG
dans ce dossier.

Écarts avec l'énoncé :

- TP4 5.1 : sur Hilbert 12 j'obtiens des erreurs de 0,28 (`solve`) et 5,7 (`inv @ b`), contre
  0,7 et 21 dans l'énoncé. Le conditionnement (1,6e16) est le même ; les erreurs dépendent
  de la version de LAPACK.
- TP4 8.1 : mon `resoudre` en flottants ne lève pas d'erreur de matrice singulière sur
  Hilbert 12, il renvoie une solution fausse (erreur 0,006).
- TP5 : le paragraphe chiffré avec la clef (7, 12) est un texte que j'ai écrit, pas celui
  de l'énoncé. L'indice de coïncidence du clair vaut 0,093 au lieu de 0,087.
- TP3 2.2 : le rapport naïf/Horner mesuré est ~3,5 au lieu de ~2.
