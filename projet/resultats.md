# Résultats des expériences

## 1. Sommation : naïve vs Kahan vs exacte
| n | erreur naïve | erreur Kahan |
|---|---|---|
| 1000 | 1.407e-12 | 0.000e+00 |
| 10000 | 1.588e-10 | 0.000e+00 |
| 100000 | 1.885e-08 | 0.000e+00 |
| 1000000 | 1.333e-06 | 0.000e+00 |

## 2. Équation du second degré x² - (1e k + 1e -k)x + 1 : petite racine 1e-k
| k | erreur relative naïve | erreur relative stable |
|---|---|---|
| 2 | 1.990e-13 | 0.000e+00 |
| 4 | 2.023e-09 | 0.000e+00 |
| 6 | 7.614e-06 | 0.000e+00 |
| 8 | 4.901e-01 | 0.000e+00 |
| 10 | 1.000e+00 | 0.000e+00 |
| 12 | 1.000e+00 | 0.000e+00 |

## 3. Polynôme de Wilkinson : racines 1..n, coefficients arrondis
Les coefficients de ∏(x-i) sont calculés en flottants, puis on cherche les racines (bissection sur [0, n+1]). Écart max à la vraie racine :
| n | écart max |
|---|---|
| 5 | 1.599e-14 |
| 10 | 2.860e-10 |
| 13 | 2.629e-08 |
| 15 | 6.708e-07 |
| 17 | 2.387e-05 |
| 20 | 3.538e-03 |

## 4. Newton : convergence quadratique vs bissection
- Newton sur x²-2 depuis 1.0 : 6 itérations, erreur 2.220e-16
- Newton sur (x-1)²(x-3) depuis 0.0 (racine double) : 32 itérations pour 1e-12 (convergence linéaire, facteur 1/2)

## 5. LU avec / sans pivot ; conditionnement de Hilbert
- Système [[1e-20,1],[1,1]] x = [1,2] (solution ≈ [1,1]) : avec pivot [1.0, 1.0], sans pivot [0.0, 1.0]

| n | cond₁(Hilbert) | erreur relative sur x (b = H·1) |
|---|---|---|
| 3 | 7.480e+02 | 2.220e-16 |
| 6 | 2.907e+07 | 6.821e-11 |
| 9 | 1.100e+12 | 1.207e-05 |
| 12 | 3.832e+16 | 7.619e-02 |

## 6. Ajustement polynomial : équations normales vs QR
Données y = polynôme de degré d (coefficients 1), 30 points sur [0,1], bruit nul. Erreur max sur les coefficients retrouvés :
| degré | équations normales | QR (Gram-Schmidt modifié) |
|---|---|---|
| 3 | 1.510e-14 | 6.417e-14 |
| 5 | 6.279e-10 | 2.276e-11 |
| 7 | 1.292e-07 | 8.558e-08 |
| 9 | 1.824e-04 | 1.614e-04 |

## 7. Géométrie : orientation exacte et enveloppe convexe
- Test d'orientation flottant pur vs exact sur 65536 triplets quasi colinéaires : 11972 signes faux (18.3 %)
- Enveloppe convexe de 5000 points aléatoires dans [0,1]² : 18 sommets, aire 0.996079

## 8. Application : intersection courbe de Bézier cubique / droite
Courbe de points de contrôle (0,0),(1,3),(2,-3),(3,0) coupée par l'axe Ox :
- t = 0.000000000000 -> point (0.000000000, 0.00e+00)
- t = 0.500000000000 -> point (1.500000000, 0.00e+00)
- t = 1.000000000000 -> point (3.000000000, 0.00e+00)

