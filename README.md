# Outils fondamentaux : module-num-rique-projet

**Nom :** Yanis El Moussafer

## Contenu

- `tp/` : mes solutions des six TP, un script par TP (voir `tp/README.md`).
- `projet/` : le projet personnel, `numgeo`.

## Le projet : numgeo

`numgeo` est une petite bibliothèque Python, sans dépendance, qui reprend les thèmes du
module (flottants, arithmétique, polynômes, matrices, géométrie plane). Pour chacun, j'ai
écrit l'algorithme naïf, mesuré où les flottants le font échouer, puis écrit la version
corrigée.

- `floats.py` : epsilon machine, somme de Kahan, somme exacte avec `Fraction`, équation du
  second degré stable. La formule du lycée perd toute précision sur la petite racine de
  x² − (10¹⁰ + 10⁻¹⁰)x + 1 ; la somme de 10⁶ fois 0,1 est fausse de 1,3e-6 sans Kahan, exacte avec.
- `poly.py` : Horner, dérivée, bissection, Newton, racines réelles sur un intervalle en
  isolant les racines avec la dérivée. Le polynôme de Wilkinson de degré 20 donne des racines
  fausses à 3,5e-3 près.
- `linalg.py` : décomposition LU avec pivot partiel, déterminant, inverse, conditionnement,
  Hilbert, moindres carrés par équations normales et par QR. Sans pivot, un système 2×2 donne
  [0, 1] au lieu de [1, 1] ; sur Hilbert 12 (conditionnement 4e16), l'erreur monte à 0,08.
- `geom.py` : test d'orientation exact (calcul flottant, puis `Fraction` quand le résultat
  est trop proche de zéro), enveloppe convexe, aire, point dans un polygone, courbes de
  Bézier, intersection Bézier/droite. Le test d'orientation en flottants seuls se trompe
  sur environ 18 % des triplets quasi alignés de mon jeu d'essai.

L'intersection d'une courbe de Bézier et d'une droite regroupe ces thèmes : on obtient un
polynôme en t, dont on isole les racines dans [0, 1] avant de les affiner.

Lancer les tests et les expériences :
cd projet
python -m unittest discover -s tests
python experiences.py
Les chiffres cités ici viennent de `projet/resultats.md`.

### Limites

- Pour les moindres carrés, j'ai utilisé Gram-Schmidt modifié, moins stable que Householder.
  Aux degrés élevés, QR ne fait pas beaucoup mieux que les équations normales dans mes essais.
- L'isolation des racines par la dérivée est fragile avec les racines multiples : je
  fusionne les racines à moins de 1e-12 l'une de l'autre.
- Le seuil du filtre de `orient` (1e-14 fois l'échelle des termes) est heuristique, il n'est
  pas prouvé comme dans les prédicats de Shewchuk.
- Je n'ai pas tracé de figures pour le projet.

## Usage de l'IA

J'ai utilisé Claude (Anthropic) pour écrire le code du projet. Les tests et les scripts ont été exécutés et passent. Ce que j'ai relu et
vérifié moi-même : Tous le Projet. Les notes sur les encadrés « L'IA vous le donne » sont
dans `tp/ia.md`.
