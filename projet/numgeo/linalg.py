"""Algèbre linéaire dense : LU pivoté, conditionnement, moindres carrés."""
import math


def matvec(A, x):
    return [sum(a * b for a, b in zip(row, x)) for row in A]


def transpose(A):
    return [list(r) for r in zip(*A)]


def matmul(A, B):
    Bt = transpose(B)
    return [[sum(a * b for a, b in zip(r, c)) for c in Bt] for r in A]


def lu(A, pivot=True):
    """Retourne (LU compact, permutation, signe). ZeroDivisionError si singulière."""
    n = len(A)
    M = [list(map(float, r)) for r in A]
    perm = list(range(n))
    sign = 1
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k])) if pivot else k
        if M[p][k] == 0.0:
            raise ZeroDivisionError("matrice singulière")
        if p != k:
            M[k], M[p] = M[p], M[k]
            perm[k], perm[p] = perm[p], perm[k]
            sign = -sign
        for i in range(k + 1, n):
            M[i][k] /= M[k][k]
            f = M[i][k]
            for j in range(k + 1, n):
                M[i][j] -= f * M[k][j]
    return M, perm, sign


def lu_solve(fact, b):
    M, perm, _ = fact
    n = len(M)
    y = [b[perm[i]] for i in range(n)]
    for i in range(n):
        y[i] -= sum(M[i][j] * y[j] for j in range(i))
    for i in reversed(range(n)):
        y[i] = (y[i] - sum(M[i][j] * y[j] for j in range(i + 1, n))) / M[i][i]
    return y


def solve(A, b, pivot=True):
    return lu_solve(lu(A, pivot), b)


def det(A):
    M, _, s = lu(A)
    return s * math.prod(M[i][i] for i in range(len(M)))


def inverse(A):
    f = lu(A)
    n = len(A)
    cols = [lu_solve(f, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
    return transpose(cols)


def norm1(A):
    return max(sum(abs(A[i][j]) for i in range(len(A))) for j in range(len(A[0])))


def cond1(A):
    """Conditionnement en norme 1 : ||A|| * ||A^-1||."""
    return norm1(A) * norm1(inverse(A))


def hilbert(n):
    return [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)]


def vandermonde(xs, deg):
    """Colonnes x^deg, ..., x, 1 (même ordre que poly.horner)."""
    return [[x ** k for k in range(deg, -1, -1)] for x in xs]


def lstsq_normal(A, y):
    """Équations normales A^T A x = A^T y : cond(A)^2, à éviter."""
    At = transpose(A)
    return solve(matmul(At, A), matvec(At, y))


def lstsq_qr(A, y):
    """Moindres carrés par Gram-Schmidt modifié (A = QR, puis R x = Q^T y)."""
    n = len(A[0])
    cols = transpose(A)
    Q, R = [], [[0.0] * n for _ in range(n)]
    for j in range(n):
        v = cols[j][:]
        for i, q in enumerate(Q):
            R[i][j] = sum(a * b for a, b in zip(q, v))
            v = [a - R[i][j] * b for a, b in zip(v, q)]
        R[j][j] = math.sqrt(sum(a * a for a in v))
        Q.append([a / R[j][j] for a in v])
    qty = [sum(a * b for a, b in zip(q, y)) for q in Q]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (qty[i] - sum(R[i][j] * x[j] for j in range(i + 1, n))) / R[i][i]
    return x
