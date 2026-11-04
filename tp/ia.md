# Encadrés « L'IA vous le donne »

Pour chaque encadré, ce que j'ai testé et ce que ça donne. Les codes « naïfs » ci-dessous
sont les versions typiques qu'un assistant propose en première réponse ; je les ai écrits
et exécutés moi-même plutôt que de coller une réponse brute. À refaire avec mon propre
assistant avant l'oral pour pouvoir en parler.

**TP1, comparaison de décimaux.** Une fonction naïve fait `a == b` : elle renvoie `False`
pour `0.1 + 0.2` et `0.3`. Avec une tolérance fixe de `1e-9`, elle marche à l'échelle de 1
mais pas à `1e20` : deux flottants voisins y sont séparés de 16384, donc `1e20` et
`1e20 + 16384` ne sont jamais jugés égaux. Il faut une tolérance relative (`math.isclose`).

**TP2, test de primalité.** Les cas limites `n = 1`, `n = 2` et `n < 0` sont corrects dans
`est_premier`. Une boucle jusqu'à `n` au lieu de `racine(n)` donne le bon résultat mais
~1000 fois plus lent sur 1 000 003 (0,019 s contre 0,00002 s).

**TP3, Horner.** La version de base évalue chaque terme avec `x ** i`. Sur des flottants le
gain de Horner n'est que d'un facteur ~3,5 (constant), pas quadratique comme on le croit
souvent ; sur des entiers le rapport croît (4 à 6).

**TP4, `inv` ou `solve`.** Sur Hilbert 12, `inv(H) @ b` se trompe d'environ 5,7 et `solve`
d'environ 0,28. `solve` fait moins de calcul et accumule moins d'erreur. Dans les deux cas
le conditionnement (1,6e16) rend le résultat inexploitable.

**TP5, chiffrement affine.** Un code qui ne vérifie pas `pgcd(a, 26) == 1` accepte
`a = 4` sans erreur : `4*0+7` et `4*13+7` donnent la même lettre, donc le déchiffrement est
ambigu et l'information est perdue silencieusement. Mon `dechiffrer` lève une erreur.

**TP6, point dans un polygone.** Avec des inégalités larges des deux côtés, `(1, 1.5)` est
jugé dehors alors qu'il est dedans. Il est à la hauteur des sommets `(2, 1.5)` et `(0, 1.5)` :
le rayon les touche et chacun est compté deux fois, ce qui fausse la parité. Les points sur
le bord, comme `(0, 1.5)`, restent un cas à traiter à part.
