#!/usr/bin/env python3
"""
Route (a) generalized to a finite effective action in SU(2)_L x SU(2)_R.

NON-EVIDENTIARY pilot code. Tuning-set groups only. No adjudication-case parameters.

WHAT CHANGES AND WHAT DOES NOT
  Unchanged: the RBF-FD weights, the left-invariant frame, and the Hodge assembly.
  Those are properties of S^3, not of Gamma. Only the ORBIT REDUCTION changes.

  Deck elements are ACTION PAIRS [u, v] acting as F(x) = u x v, applied as 4x4
  rotations. The scalar reduction is plain orbit summation. The one-form reduction
  is Ad-EQUIVARIANT, from the pre-registration section 2 law

      f(gamma . x) = Ad(v_gamma^-1) f(x)

  which was fixed by a manufactured pointwise pullback test, not by inference. So for
  node j in orbit l, reached from the representative by gamma_j,

      R[aM+m, cM+l] = sum_{j in orbit l} sum_b A_ab[i_m, j] [Ad(v_{gamma_j}^-1)]_bc

  The index order in that accumulation is itself a hazard: Ad is orthogonal, so
  Ad(v)^T and Ad(v^-1) are the SAME matrix and no pointwise test can separate them.
  What separates them is where the rotation sits in the sum, which is why
  assembly_selftest below compares the reduced operator against the full one.
"""

import numpy as np

from route_a_tuned import build_L
from route_a_oneform import left_frame, build_dirderiv
from route_a_hodge import hodge_matrix

ONE = np.array([1.0, 0.0, 0.0, 0.0])
IMAG = [np.array([0., 1, 0, 0]), np.array([0., 0, 1, 0]), np.array([0., 0, 0, 1])]


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2])


def qconj(a):
    return np.array([a[0], -a[1], -a[2], -a[3]])


def act_matrix(u, v):
    """4x4 rotation of x -> u x v."""
    E = np.eye(4)
    return np.array([qmul(qmul(u, E[i]), v) for i in range(4)]).T


def Ad(v):
    """3x3 rotation with columns v u_j v^-1, in the (i,j,k) basis."""
    vi = qconj(v)
    return np.array([qmul(qmul(v, u), vi)[1:] for u in IMAG]).T


def rot(t):
    return np.array([np.cos(t), np.sin(t), 0.0, 0.0])


def pairs_lens(q, s1, s2):
    """Deck pairs of L(q; s1, s2). Requires q odd so 2 is invertible mod q."""
    i2 = pow(2, -1, q)
    al, be = ((s1 + s2) * i2) % q, ((s1 - s2) * i2) % q
    return [(rot(2*np.pi*al*m/q), rot(2*np.pi*be*m/q)) for m in range(q)]


def close_pairs(gens, tol=1e-9, cap=512):
    """Close a deck group from GENERATOR pairs alone.

    The routes derive the effective action and its order from the generators supplied
    in Packet I. The declared case parameters are NOT used as an expected-order
    fixture: an order that is told to you is not an order you checked. Composition is
    (u1,v1) o (u2,v2) = (u1 u2, v2 v1), and pairs are identified modulo the central
    kernel [u, v] ~ [-u, -v], so what is counted is the EFFECTIVE element.
    """
    out = [(ONE.copy(), ONE.copy())]
    frontier = list(out)

    def same(a, b):
        return ((np.abs(a[0] - b[0]).max() < tol and np.abs(a[1] - b[1]).max() < tol) or
                (np.abs(a[0] + b[0]).max() < tol and np.abs(a[1] + b[1]).max() < tol))

    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = (qmul(g[0], a[0]), qmul(a[1], g[1]))
                if not any(same(b, e) for e in out):
                    out.append(b); nxt.append(b)
                    if len(out) > cap:
                        raise ValueError("closure cap exceeded")
        frontier = nxt
    return out


def pairs_left(elems):
    """A left action written in pair form, identity first."""
    out = [(ONE.copy(), ONE.copy())]
    for g in elems:
        if np.abs(g - ONE).max() > 1e-12:
            out.append((g.copy(), ONE.copy()))
    return out


def pairs_swapped(pairs):
    """[u, v] -> [v, u]. Used only to exercise the left/right swap mutation."""
    return [(v.copy(), u.copy()) for u, v in pairs]


def cloud(seeds, pairs):
    """Orbit cloud. Node order per seed follows `pairs`, so gid 0 is the rep."""
    mats = [act_matrix(u, v) for u, v in pairs]
    X, oid, gid = [], [], []
    for l, sd in enumerate(seeds):
        for g, Mg in enumerate(mats):
            X.append(Mg @ sd); oid.append(l); gid.append(g)
    return np.array(X), np.array(oid), np.array(gid), len(seeds)


def reduce_scalar(A, oid, M):
    reps = [np.where(oid == m)[0][0] for m in range(M)]
    R = np.zeros((M, M))
    for mm, i in enumerate(reps):
        row = A[i]
        for l in range(M):
            R[mm, l] = row[oid == l].sum()
    return R


def reduce_oneform(A, oid, gid, pairs, M, N, order="correct"):
    """Ad-equivariant reduction. `order` exists so the accumulation can be mutated."""
    reps = [np.where(oid == m)[0][0] for m in range(M)]
    AdI = [Ad(qconj(v)) for _, v in pairs]      # Ad(v^-1), one per group element
    R = np.zeros((3 * M, 3 * M))
    for l in range(M):
        idx = np.where(oid == l)[0]
        for mm, i in enumerate(reps):
            for a in range(3):
                acc = np.zeros(3)
                for j in idx:
                    W = AdI[gid[j]]
                    row = A[a*N + i, :]
                    blk = np.array([row[b*N + j] for b in range(3)])
                    if order == "correct":
                        acc += blk @ W                 # sum_b A_ab W_bc
                    elif order == "transposed":
                        acc += W @ blk                 # wrong index order
                    elif order == "omitted":
                        acc += blk
                R[a*M + mm, np.arange(3)*M + l] = acc
    return R


def assembly_selftest(A, oid, gid, pairs, M, N, order="correct", trials=3, seed=11):
    """Does the reduced operator agree with the full one on an invariant field?

    Build f on representatives, extend it by the frozen equivariance law, apply the
    FULL operator, restrict to representatives, and compare against the reduced
    operator applied to the representative values. This is the assembly-level check
    that the pointwise pullback test cannot do, because it is sensitive to the index
    order in the accumulation.
    """
    rng = np.random.default_rng(seed)
    reps = [np.where(oid == m)[0][0] for m in range(M)]
    AdI = [Ad(qconj(v)) for _, v in pairs]
    R = reduce_oneform(A, oid, gid, pairs, M, N, order=order)
    worst = 0.0
    for _ in range(trials):
        fr = rng.normal(size=(M, 3))
        full = np.zeros((N, 3))
        for j in range(N):
            full[j] = AdI[gid[j]] @ fr[oid[j]]
        vec = np.concatenate([full[:, a] for a in range(3)])
        got_full = A @ vec
        lhs = np.concatenate([got_full[a*N + np.array(reps)] for a in range(3)])
        rhs = R @ np.concatenate([fr[:, a] for a in range(3)])
        worst = max(worst, np.abs(lhs - rhs).max() / max(1e-30, np.abs(lhs).max()))
    return worst


def operators(X, k=110, m=7, p=4):
    """Scalar positive Laplacian and the Hodge operator on the cloud."""
    Lpos = -build_L(X, k, m, p)
    D = build_dirderiv(X, left_frame(X), k, m, p)
    return Lpos, hodge_matrix(Lpos, D, len(X))


def relax(seeds, pairs, iters=40, step=0.12):
    """Riesz relaxation of seeds against the full orbit cloud."""
    mats = [act_matrix(u, v) for u, v in pairs]
    S = seeds.copy()
    for _ in range(iters):
        P = np.vstack([np.array([Mg @ s for Mg in mats]) for s in S])
        F = np.zeros_like(S)
        for i, s in enumerate(S):
            d = P - s
            r2 = (d ** 2).sum(1)
            keep = r2 > 1e-12
            F[i] = -(d[keep] / r2[keep][:, None] ** 1.5).sum(0)
        S = S + step * F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-12)
        S /= np.linalg.norm(S, axis=1, keepdims=True)
    return S
