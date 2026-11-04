# numgeo

Petite bibliothèque Python, sans dépendance, qui reprend les thèmes du module : flottants,
arithmétique, polynômes, matrices et géométrie plane. Pour chaque thème, j'ai écrit
l'algorithme naïf, mesuré où il échoue à cause des flottants, puis écrit la version corrigée.

## Contenu

```
projet/
├── numgeo/
│   ├── floats.py   epsilon machine, somme de Kahan, somme exacte, second degré stable
│   ├── poly.py     Horner, dérivée, bissection, Newton, racines réelles sur un intervalle
│   ├── linalg.py   LU avec pivot, déterminant, inverse, conditionnement, moindres carrés
│   └── geom.py     orientation exacte, enveloppe convexe, aire, Bézier, intersection avec une droite
├── tests/
│   └── test_numgeo.py   16 tests
├── experiences.py       lance les mesures et écrit resultats.md
└── resultats.md         résultats chiffrés
```

## Lancer

Python 3 suffit, il n'y a rien à installer.

```bash
cd projet
python -m unittest discover -s tests
python experiences.py
```

`experiences.py` régénère `resultats.md`.

## Ce que montrent les expériences

- La somme naïve de 10⁶ fois 0,1 est fausse de 1,3e-6. Avec la somme de Kahan, l'erreur est nulle.
- La formule du lycée pour le second degré perd toute précision sur la petite racine dès
  que les deux racines sont très éloignées. La formule stable reste exacte.
- Les racines du polynôme de Wilkinson de degré 20 sont fausses à 3,5e-3 près.
- Sans pivot, un système 2×2 donne [0, 1] au lieu de [1, 1]. Sur la matrice de Hilbert de
  taille 12, dont le conditionnement vaut environ 4e16, l'erreur monte à 0,08.
- Le test d'orientation en flottants seuls se trompe sur environ 18 % des triplets quasi
  alignés de mon jeu d'essai. Le filtre flottant suivi d'un calcul exact avec `Fraction`
  donne toujours le bon signe.

L'intersection d'une courbe de Bézier et d'une droite rassemble ces thèmes : on obtient un
polynôme en t, dont on isole les racines dans [0, 1] avant de les affiner.

## Limites

- Les moindres carrés utilisent Gram-Schmidt modifié, moins stable que Householder. Aux
  degrés élevés, QR ne fait pas beaucoup mieux que les équations normales dans mes essais.
- L'isolation des racines par la dérivée est fragile avec les racines multiples : je fusionne
  les racines à moins de 1e-12 l'une de l'autre.
- Le seuil du filtre de `orient` (1e-14 fois l'échelle des termes) est heuristique, pas prouvé.
