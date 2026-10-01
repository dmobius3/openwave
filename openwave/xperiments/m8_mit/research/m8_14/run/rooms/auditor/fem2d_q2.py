"""Method A: biquadratic (Q2) tensor-product finite elements directly in (y, w)
on the rectangle [0, pi] x [-W, W], with the anti-periodic reflected seam
phi(pi, -w) = -phi(0, w) imposed as a nodal identification, area weight
|cos y|, zero data on w = +-W and on the collapsed fiber y = pi/2.

y-mesh: nodes at distance s_k = (pi/2)(k/N)^G from the fiber on each side
(G = 3, cubic grading), Q2 midpoints at element centres; refinement doubles N
(nested).  w-mesh: M uniform Q2 elements across [-W, W] (symmetric, so the
reflection w -> -w maps nodes to nodes); M = N/2.

Bilinear forms (R = 1):  a(u,v) = int |cos y| u_y v_y + u_w v_w / |cos y|,
m(u,v) = int |cos y| u v.  With the tensor structure,
K = A_y (x) M_w + B_y (x) K_w,   M = M_y (x) M_w.
"""
import math
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from common import WIDTHS, NEIG, element_matrices_s, graded_s, richardson, save_json

G = 3
NQUAD = 8
LEVELS = [8, 16, 32, 64, 128, 256]


def y_matrices(N, g=G, nquad=NQUAD):
    """Global 1D Q2 matrices in y on [0, pi]; node index 0 is y=0, 2N is the
    fiber y=pi/2, 4N is y=pi.  Returns (A, B, Mm, y) as dense-free sparse."""
    s = graded_s(N, g)
    n = 4 * N + 1
    rows, cols, va, vb, vm = [], [], [], [], []
    for k in range(N):
        A, B, Mm = element_matrices_s(s[k], s[k + 1], 2, nquad, None)
        hs = [2 * k, 2 * k + 1, 2 * k + 2]  # half-indices (s-order)
        for side in (-1, +1):
            idx = [2 * N + side * h for h in hs]
            for a in range(3):
                for b in range(3):
                    rows.append(idx[a]); cols.append(idx[b])
                    va.append(A[a, b]); vb.append(B[a, b]); vm.append(Mm[a, b])
    mk = lambda v: sp.csr_matrix((v, (rows, cols)), shape=(n, n))
    return mk(va), mk(vb), mk(vm)


def w_matrices(W, M):
    """Global 1D Q2 matrices on [-W, W] with M uniform elements.
    Returns (K_w, M_w); node j at w = -W + j W / M, j = 0..2M."""
    xi, wq = np.polynomial.legendre.leggauss(6)
    V = np.array([xi * (xi - 1) / 2, 1 - xi ** 2, xi * (xi + 1) / 2])
    D = np.array([xi - 0.5, -2 * xi, xi + 0.5])
    h = 2 * W / M
    jac = h / 2
    Ke = (D * (wq / jac)) @ D.T
    Me = (V * (wq * jac)) @ V.T
    n = 2 * M + 1
    rows, cols, vk, vm = [], [], [], []
    for e in range(M):
        idx = [2 * e, 2 * e + 1, 2 * e + 2]
        for a in range(3):
            for b in range(3):
                rows.append(idx[a]); cols.append(idx[b])
                vk.append(Ke[a, b]); vm.append(Me[a, b])
    mk = lambda v: sp.csr_matrix((v, (rows, cols)), shape=(n, n))
    return mk(vk), mk(vm)


def seam_prolongation(N, M, seam_sign=-1.0, reflect=True):
    """P maps reduced dofs to full nodal values.  Full index = i*(2M+1)+j with
    y-node i in 0..4N and w-node j in 0..2M.  Reduced dofs: i in 0..4N-1 except
    the fiber i = 2N, j in 1..2M-1.  Row i = 4N (y = pi) is
    phi(pi, w_j) = seam_sign * phi(0, w_{2M-j})  (reflect=True), i.e.
    phi(pi, -w) = seam_sign * phi(0, w)."""
    ny, nw = 4 * N + 1, 2 * M + 1
    red = {}
    for i in range(4 * N):
        if i == 2 * N:
            continue
        for j in range(1, 2 * M):
            red[(i, j)] = len(red)
    rows, cols, vals = [], [], []
    for (i, j), r in red.items():
        rows.append(i * nw + j); cols.append(r); vals.append(1.0)
    for j in range(1, 2 * M):
        jj = (2 * M - j) if reflect else j
        rows.append(4 * N * nw + j); cols.append(red[(0, jj)]); vals.append(seam_sign)
    P = sp.csr_matrix((vals, (rows, cols)), shape=(ny * nw, len(red)))
    return P


def assemble(W, N, M, seam_sign=-1.0, reflect=True, g=G, nquad=NQUAD):
    Ay, By, My = y_matrices(N, g, nquad)
    Kw, Mw = w_matrices(W, M)
    K = sp.kron(Ay, Mw) + sp.kron(By, Kw)
    Mm = sp.kron(My, Mw)
    P = seam_prolongation(N, M, seam_sign, reflect)
    Kr = (P.T @ K @ P).tocsc()
    Mr = (P.T @ Mm @ P).tocsc()
    return Kr, Mr, P


def solve(W, N, M=None, k=NEIG, seam_sign=-1.0, reflect=True, g=G, nquad=NQUAD,
          return_vectors=False):
    if M is None:
        M = N // 2
    Kr, Mr, P = assemble(W, N, M, seam_sign, reflect, g, nquad)
    vals, vecs = spla.eigsh(Kr, k=k, M=Mr, sigma=0.0, which="LM", tol=1e-13)
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    # relative residual of every returned pair at this level
    res = [float(np.linalg.norm(Kr @ vecs[:, i] - vals[i] * (Mr @ vecs[:, i]))
                 / np.linalg.norm(Kr @ vecs[:, i])) for i in range(k)]
    out = {"vals": vals, "res": res, "ndof": Kr.shape[0]}
    if return_vectors:
        out.update(vecs=vecs, P=P, Kr=Kr, Mr=Mr, N=N, M=M)
    return out


def main():
    results = {"method": "A: Q2 tensor FEM in (y,w), grading exponent %d, M=N/2" % G,
               "levels_N": LEVELS, "widths": {}}
    for label, W in WIDTHS:
        levels = []
        for N in LEVELS:
            t0 = time.time()
            r = solve(W, N)
            levels.append({"N": N, "M": N // 2, "ndof": r["ndof"],
                           "lowest": [float(v) for v in r["vals"]],
                           "max_rel_residual": max(r["res"]),
                           "seconds": time.time() - t0})
            print("A W=%-5s N=%4d ndof=%7d  %s  res=%.1e" % (
                label, N, r["ndof"], " ".join("%.12f" % v for v in r["vals"]),
                max(r["res"])), flush=True)
        rich = []
        for i in range(NEIG):
            seq = [lv["lowest"][i] for lv in levels]
            rr = richardson(seq)
            # second estimate from the preceding triple, to test stability of p
            rr_prev = richardson(seq[:-1])
            rr["extrap_prev_triple"] = rr_prev["extrap"]
            rr["p_prev_triple"] = rr_prev["p"]
            rr["err_combined"] = max(rr["err"], abs(rr["extrap"] - rr_prev["extrap"]))
            rich.append(rr)
        results["widths"][label] = {"W": W, "levels": levels, "richardson": rich}
        print("A W=%-5s bottom extrap=%.12f p=%s err=%.2e (combined %.2e)" % (
            label, rich[0]["extrap"], rich[0]["p"], rich[0]["err"],
            rich[0]["err_combined"]), flush=True)
    save_json("results_A.json", results)


if __name__ == "__main__":
    main()
