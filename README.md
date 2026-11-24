# Solveur de matrices bandes – C++

Mini-projet de L3 Ingénierie Mathématique (Université de Bordeaux), repris et corrigé en 2026.

- LU avec pivotage partiel (bande élargie à ku + kl), Cholesky, SOR.
- Normes 1, 2, infinie et Frobenius ; déterminant et log|det|.
- Lecture du format bande et de MatrixMarket (`.mtx`).

`mat_vec_io` (lecture et écriture des matrices) a été fourni par l'enseignant.

## Lancer
```sh
make
./mini_prog sherman1.mtx resultat.txt    # matrices : https://math.nist.gov/MatrixMarket
make test                                # comparaison à un calcul dense, matrices aléatoires
```

## Résultats (LU, second membre A·1)

| Matrice | n | kl = ku | Temps LU (s) | Erreur max | Résidu max |
|---|---:|---:|---:|---:|---:|
| sherman1 | 1000 | 100 | 0,0015 | 7,8·10⁻¹³ | 3,8·10⁻¹⁵ |
| sherman2 | 1080 | 221 | 0,033 | 3,3·10⁻⁷ | 6,0·10⁻⁷ |
| sherman3 | 5005 | 385 | 0,086 | 6,4·10⁻¹² | 8,7·10⁻¹⁰ |
| sherman4 | 1104 | 368 | 0,014 | 3,3·10⁻¹⁴ | 3,8·10⁻¹⁴ |
| sherman5 | 3312 | 1106 | 0,21 | 9,3·10⁻¹⁴ | 1,9·10⁻¹² |

Temps : une seule exécution, Apple M1, `-O2`. SOR (ω = 1) ne converge en 10 000 itérations que pour sherman4.
