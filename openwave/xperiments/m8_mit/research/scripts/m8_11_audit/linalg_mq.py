"""Exact linear algebra over MQ."""
from mq import MQ, ZERO, ONE


def rref(A):
    rows = [list(r) for r in A]
    nr = len(rows)
    nc = len(rows[0]) if nr else 0
    piv = []
    r = 0
    for c in range(nc):
        p = None
        for i in range(r, nr):
            if not rows[i][c].is_zero():
                p = i
                break
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        inv = rows[r][c].inv()
        rows[r] = [x * inv for x in rows[r]]
        for i in range(nr):
            if i != r and not rows[i][c].is_zero():
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == nr:
            break
    return rows, piv


def nullspace(A):
    nc = len(A[0])
    rows, piv = rref(A)
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fc in free:
        v = [ZERO] * nc
        v[fc] = ONE
        for i, pc in enumerate(piv):
            v[pc] = -rows[i][fc]
        basis.append(v)
    return basis


def rank(A):
    return len(rref(A)[1])


def solve(A, b):
    """solve square nonsingular A x = b."""
    n = len(A)
    aug = [list(A[i]) + [b[i]] for i in range(n)]
    rows, piv = rref(aug)
    assert piv == list(range(n)), "singular"
    return [rows[i][n] for i in range(n)]


def matvec(A, x):
    return [sum((a * b for a, b in zip(row, x) if not a.is_zero() and not b.is_zero()), ZERO) for row in A]


def transpose(A):
    return [list(c) for c in zip(*A)]
