"""Method B: exact Fourier-sine expansion in w, high-order 1D finite elements
in y for each transverse sector.

Sector k uses S_k(w) = sin(k pi (w+W)/(2W)), the Dirichlet eigenfunctions of
-d^2/dw^2 on [-W, W] with eigenvalue mu_k = (k pi/(2W))^2, and
S_k(-w) = (-1)^(k+1) S_k(w).  For phi = u(y) S_k(w) the seam condition
phi(pi, -w) = -phi(0, w) becomes u(pi) = (-1)^k u(0), imposed as a nodal
identification of y = pi with y = 0 (the derivative matching is the natural
condition of the form).  The fiber y = pi/2 carries zero data.  Each sector:
  a_k(u,u) = int_0^pi |cos y| u'^2 + mu_k u^2/|cos y| dy,  m(u,u) = int |cos y| u^2.
Degree-Q Lagrange elements on Gauss-Lobatto nodes, nodes graded toward the
fiber as s_k = (pi/2)(k/N)^G on each side; refinement doubles N (nested).

The sectors are invariant subspaces of the 2D problem and the sines are
complete on [-W, W], so the 2D spectrum is the union over k of the sector
spectra.  Sector k has every eigenvalue >= mu_k (1/|cos y| >= |cos y|), so
sectors are added until mu_k exceeds the current sixth-lowest value.
"""
import math

import numpy as np
import scipy.linalg as sla

from common import WIDTHS, NEIG, element_matrices_s, graded_s, richardson, save_json

G = 4
Q = 8
LEVELS = [2, 4, 8, 16, 32]
Q_ALT = 6  # second polynomial degree, for a p-direction cross-check


def sector_eigs(W, k, N, q=Q, g=G, nev=NEIG, seam_sign=-1.0, mu_scale=1.0, nquad=None):
    """Lowest nev eigenvalues of sector k.  seam_sign is the sign in
    phi(pi,-w) = seam_sign*phi(0,w); mu_scale is for planted-defect tests."""
    if nquad is None:
        nquad = 2 * q + 6
    mu = mu_scale * (k * math.pi / (2 * W)) ** 2
    s = graded_s(N, g)
    nh = q * N + 1               # nodes per half, half-index 0 at the fiber
    n = 2 * (nh - 1) + 1         # global y-nodes; index nh-1 is the fiber
    fib = nh - 1
    A = np.zeros((n, n)); Mm = np.zeros((n, n))
    for e in range(N):
        Ae, Be, Me = element_matrices_s(s[e], s[e + 1], q, nquad, None)
        hs = np.arange(q * e, q * e + q + 1)
        for side in (-1, +1):
            idx = fib + side * hs
            A[np.ix_(idx, idx)] += Ae + mu * Be
            Mm[np.ix_(idx, idx)] += Me
    # y = 0 is index 0, y = pi is index n-1.  Seam: u(pi) = sigma u(0),
    # sigma = seam_sign * (-1)^(k+1)   [from S_k(-w) = (-1)^(k+1) S_k(w)].
    sigma = seam_sign * (-1) ** (k + 1)
    keep = [i for i in range(n - 1) if i != fib]
    P = np.zeros((n, len(keep)))
    for c, i in enumerate(keep):
        P[i, c] = 1.0
    P[n - 1, 0] = sigma
    Kr = P.T @ A @ P
    Mr = P.T @ Mm @ P
    # On the graded mesh Mr is nearly singular (tip elements of size ~N^-G),
    # while Kr is well scaled; so solve Mr x = theta Kr x for the largest
    # theta = 1/lambda (Cholesky of Kr), which is the same pencil.
    nn = len(keep)
    m = min(nev, nn)
    theta = sla.eigh(Mr, Kr, eigvals_only=True, subset_by_index=[nn - m, nn - 1])
    vals = np.sort(1.0 / theta)
    return mu, vals


def spectrum(W, N, q=Q, g=G, seam_sign=-1.0, mu_scale=1.0, nquad=None):
    """Lowest NEIG eigenvalues of the 2D problem, with sector labels."""
    pairs = []
    k = 1
    sectors_used = []
    while True:
        mu = mu_scale * (k * math.pi / (2 * W)) ** 2
        if len(pairs) >= NEIG and mu >= sorted(pairs)[NEIG - 1][0]:
            break
        _, vals = sector_eigs(W, k, N, q, g, NEIG, seam_sign, mu_scale, nquad)
        pairs += [(float(v), k, m) for m, v in enumerate(vals)]
        sectors_used.append(k)
        k += 1
    pairs.sort()
    return pairs[:NEIG], sectors_used


def study(W, q=Q, levels=LEVELS, g=G):
    lv = []
    for N in levels:
        low, used = spectrum(W, N, q, g)
        lv.append({"N": N, "lowest": [p[0] for p in low],
                   "labels_k_m": [[p[1], p[2]] for p in low], "sectors": used})
    rich = []
    for i in range(NEIG):
        seq = [x["lowest"][i] for x in lv]
        rr = richardson(seq)
        rr_prev = richardson(seq[:-1])
        rr["extrap_prev_triple"] = rr_prev["extrap"]
        rr["p_prev_triple"] = rr_prev["p"]
        rich.append(rr)
    return lv, rich


def main():
    results = {"method": "B: sine modes in w (exact), degree-%d GLL Lagrange FEM in y, "
                         "grading exponent %d, N per side doubled" % (Q, G),
               "levels_N": LEVELS, "q": Q, "q_alt": Q_ALT, "widths": {}}
    for label, W in WIDTHS:
        lv, rich = study(W, Q)
        lv6, rich6 = study(W, Q_ALT)
        best, err = [], []
        for i in range(NEIG):
            r = rich[i]
            fin = lv[-1]["lowest"][i]
            # error estimate: the larger of |extrap - finest| (h-direction, q=Q),
            # |extrap - extrap from the preceding triple|, |finest(Q) - finest(Q_ALT)|
            # and a roundoff floor of 1e-12 relative.
            e = max(r["err"], abs(r["extrap"] - r["extrap_prev_triple"]),
                    abs(fin - lv6[-1]["lowest"][i]), 1e-12 * abs(fin))
            best.append(r["extrap"]); err.append(e)
        results["widths"][label] = {"W": W, "levels": lv, "richardson": rich,
                                    "levels_q_alt": lv6, "richardson_q_alt": rich6,
                                    "best": best, "err": err}
        for x in lv:
            print("B W=%-5s N=%3d  %s" % (label, x["N"], " ".join("%.14f" % v for v in x["lowest"])))
        for x in lv6:
            print("B6 W=%-5s N=%3d  %s" % (label, x["N"], " ".join("%.14f" % v for v in x["lowest"])))
        print("B W=%-5s labels(k,m)=%s" % (label, lv[-1]["labels_k_m"]))
        print("B W=%-5s best=%s" % (label, " ".join("%.14f" % v for v in best)))
        print("B W=%-5s err =%s  p(bottom)=%s" % (label, " ".join("%.1e" % v for v in err), rich[0]["p"]),
              flush=True)
    save_json("results_B.json", results)


if __name__ == "__main__":
    main()
