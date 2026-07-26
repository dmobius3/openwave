#!/usr/bin/env python3
"""m8_1_eigensolve_xcheck.py

T3 cross-check of the T1 (realization A) spectrum at W = 1.0 and W = 2.2 with a
methodologically different discretization: a DIRECT 2D sparse FEM on the (y, w)
rectangle (bilinear quads), no separation of variables.

Operator and realization implemented (same spec as m8_1_eigensolve.py):
  energy form  q(psi) = int |f| psi_y^2 + (1/|f|) psi_w^2  dy dw,
  mass         m(psi) = int |f| psi^2 dy dw,   f = cos y, R = 1.
  - Twisted seam: node (y=pi, w_j) identified with -1 * node (y=0, -w_j).
  - Neumann at w = +-W: natural in the weak form (nothing imposed).
  - Cone line y = pi/2: the mesh is split into two blocks (left/right of the
    cone) with NO continuity across it; on each block the cone-edge DOFs are
    constrained to a single constant c_L (resp. c_R):
      * the transverse-CONSTANT trace stays free with the natural condition
        => u_N = 0 on each side (log branch has infinite energy) = realization A
      * every mu>0 transverse component gets zero trace = regular branch.
    The constrained cone basis function is constant along the edge, so its
    w-derivative vanishes identically; this is imposed EXACTLY by zeroing the
    d/dw of cone-edge local shape functions during assembly of the (1/|f|)
    psi_w psi_w term (the constrained space has no other w-variation there).

Also re-lists, side by side:
  - method 1: 1D FEM + Richardson (imported from m8_1_eigensolve.py)
  - method 3 (mu=0 channel only): exact Frobenius-series zeros (W-independent).

Output: ../data/m8_1_xcheck.json
"""

import json
import os
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

import m8_1_eigensolve as m1

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
PI = np.pi

G1, GW1 = np.polynomial.legendre.leggauss(3)


def solve_2d(W, ny_half, nw, k=10, p=2):
    """Lowest k eigenvalues of the 2D FEM (realization A). Returns sorted array."""
    yl = m1.graded_nodes(ny_half, p)          # left block nodes: 0 .. pi/2
    yr = (PI - yl)[::-1]                      # right block nodes: pi/2 .. pi
    wv = np.linspace(-W, W, nw + 1)           # symmetric uniform w grid

    nwn = nw + 1
    # global DOF numbering
    #   left block rows i=0..ny_half-1 (i=ny_half is the cone edge -> c_L)
    #   c_L, c_R
    #   right block rows i=1..ny_half-1 (i=0 -> c_R ; i=ny_half -> seam map)
    c_L = ny_half * nwn
    c_R = c_L + 1
    ndof = (2 * ny_half - 1) * nwn + 2

    def node_map(block, i, j):
        """(dof, sign) for node (i, j) of block 'L'/'R'."""
        if block == "L":
            if i == ny_half:
                return c_L, 1.0
            return i * nwn + j, 1.0
        # right block
        if i == 0:
            return c_R, 1.0
        if i == ny_half:
            # seam: psi(pi, w_j) = -psi(0, -w_j) = -psi(0, w_{nw-j})
            return (nw - j), -1.0
        return c_R + 1 + (i - 1) * nwn + j, 1.0

    rows, cols, kdat, mdat = [], [], [], []
    for block, ynod in (("L", yl), ("R", yr)):
        cone_i = ny_half if block == "L" else 0   # element-row index touching cone
        for i in range(ny_half):
            y0, y1 = ynod[i], ynod[i + 1]
            hy = y1 - y0
            touches_cone = (block == "L" and i == ny_half - 1) or \
                           (block == "R" and i == 0)
            # local nodes: 0=(i,j) 1=(i+1,j) 2=(i,j+1) 3=(i+1,j+1)
            cone_locals = ((1, 3) if block == "L" else (0, 2)) if touches_cone else ()
            for j in range(nw):
                w0, w1 = wv[j], wv[j + 1]
                hw = w1 - w0
                Kl = np.zeros((4, 4))
                Ml = np.zeros((4, 4))
                for gx, wx in zip(G1, GW1):
                    xi = 0.5 * (1 + gx)
                    y = y0 + hy * xi
                    a = abs(np.cos(y))
                    Ny = np.array([1 - xi, xi, 1 - xi, xi])
                    dNy = np.array([-1, 1, -1, 1]) / hy
                    for gy, wy in zip(G1, GW1):
                        eta = 0.5 * (1 + gy)
                        Nw = np.array([1 - eta, 1 - eta, eta, eta])
                        dNw = np.array([-1, -1, 1, 1]) / hw
                        N = Ny * Nw
                        By = dNy * Nw                    # d/dy of shapes
                        Bw = Ny * dNw                    # d/dw of shapes
                        # exact constrained space: cone-edge basis is constant
                        # along the edge -> its w-derivative is identically 0
                        for cl in cone_locals:
                            Bw = Bw.copy()
                            Bw[cl] = 0.0
                        wq = wx * wy * 0.25 * hy * hw
                        Kl += wq * (a * np.outer(By, By) +
                                    (1.0 / a) * np.outer(Bw, Bw))
                        Ml += wq * a * np.outer(N, N)
                nd = [node_map(block, i, j), node_map(block, i + 1, j),
                      node_map(block, i, j + 1), node_map(block, i + 1, j + 1)]
                for aa in range(4):
                    da, sa = nd[aa]
                    for bb in range(4):
                        db, sb = nd[bb]
                        rows.append(da)
                        cols.append(db)
                        kdat.append(sa * sb * Kl[aa, bb])
                        mdat.append(sa * sb * Ml[aa, bb])
    K = sp.coo_matrix((kdat, (rows, cols)), shape=(ndof, ndof)).tocsc()
    M = sp.coo_matrix((mdat, (rows, cols)), shape=(ndof, ndof)).tocsc()
    vals = eigsh(K, k=k, M=M, sigma=-0.5, which="LM", return_eigenvectors=False)
    return np.sort(vals)


def main():
    os.makedirs(DATA, exist_ok=True)
    out = {"task": "T3 cross-check of realization (A)",
           "method2": "direct 2D sparse bilinear FEM (no separation of variables)",
           "per_W": {}}
    # method 3: exact mu=0-channel eigenvalues (series zeros; W-independent)
    grid = np.linspace(1e-6, 45.0, 4501)
    vr = m1.mu0_seam_series(grid)[0]
    dvr = m1.mu0_seam_series(grid)[1]
    zr = m1.refine_zeros(lambda x: float(m1._vreg(np.array([x]))[0]), grid, vr)
    zdr = m1.refine_zeros(lambda x: float(m1._dvreg(np.array([x]))[0]), grid, dvr)
    mu0_exact = sorted([0.0] + zr + zdr)
    out["mu0_exact_series_zeros"] = mu0_exact[:8]

    for W in (1.0, 2.2):
        # method 1: 1D FEM + Richardson (recomputed; identical to T1 pipeline)
        sec = m1.run_fem_all(W)
        merged = []
        for n, d in sec.items():
            for e in d["eigs"]:
                merged.append({"lambda": e["lambda"], "err": e["err"],
                               "sector_n": n})
        merged.sort(key=lambda e: e["lambda"])
        m1_low = merged[:8]

        # method 2: 2D FEM at two resolutions + 2-point Richardson (q=2 assumed)
        res = [(96, 40), (144, 60)]
        raw = {f"{ny}x{nw}": [float(v) for v in solve_2d(W, ny, nw, k=10)]
               for ny, nw in res}
        r1 = np.array(raw[f"{res[0][0]}x{res[0][1]}"])
        r2 = np.array(raw[f"{res[1][0]}x{res[1][1]}"])
        rr = (res[1][0] / res[0][0]) ** 2               # h ratio^2 (q=2)
        extr = (rr * r2 - r1) / (rr - 1.0)
        out["per_W"][f"{W:g}"] = {
            "method1_1dfem_richardson_lowest8": m1_low,
            "method2_2dfem_raw": raw,
            "method2_resolutions_ny_half_x_nw": res,
            "method2_2dfem_extrapolated_q2": [float(v) for v in extr],
        }
        print(f"W={W:g}")
        for i, e in enumerate(m1_low):
            v2 = extr[i] if i < len(extr) else float("nan")
            print(f"  #{i}: 1D {e['lambda']:.8f}+-{e['err']:.1e} (n={e['sector_n']})"
                  f"   2D {v2:.8f}   diff {v2 - e['lambda']:+.2e}")

    with open(os.path.join(DATA, "m8_1_xcheck.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print("written", os.path.join(DATA, "m8_1_xcheck.json"))


if __name__ == "__main__":
    main()
