"""T4: P1 finite elements for -Delta on the admissible class, R = 1.

Domain: rectangle [0, pi] x [-W, W] in (y, w), metric dy^2 + cos^2 y dw^2.
  stiffness  int ( |cos y| phi_y^2 + phi_w^2 / |cos y| ) dy dw
  mass       int |cos y| phi^2 dy dw
Seam: node (pi, -w_j) is the dof of node (0, w_j) with sign -1.
Zero data on w = +-W and on the fiber y = pi/2.
Mesh: y-nodes pi/2 -+ (pi/2)(k/N)^2, k = 0..N, N = 16; 8 uniform cells in w.
Each rectangular cell is cut into two triangles by its (y0,w0)-(y1,w1) diagonal.
Uniform refinement = bisect every cell in y and in w (midpoints in (y, w)),
applied four times; this is the reading of "uniform refinement" taken here.
Quadrature: on each triangle, 12-point Gauss-Legendre in y and 3-point
Gauss-Legendre in w over the w-section; the gradient is constant per
triangle, so the stiffness only needs the weight integrals.

Pass --plant to plant defects: (a) periodic instead of anti-periodic seam,
(b) weight cos y instead of |cos y|. The PASS checks must then print FAIL.
"""
import json
import os
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

PLANT = "--plant" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
NEIG = 6
WIDTHS = {"1/4": 0.25, "1/2": 0.5, "1": 1.0, "7/5": 1.4, "3/2": 1.5,
          "pi/2": np.pi/2, "17/10": 1.7}

GY, GYW = np.polynomial.legendre.leggauss(12)
GW, GWW = np.polynomial.legendre.leggauss(3)


def weight(y):
    c = np.cos(y)
    return c if PLANT else np.abs(c)


def base_mesh(N=16, nw=8):
    k = np.arange(N + 1)
    left = np.pi/2 - (np.pi/2)*(k[::-1]/N)**2      # 0 ... pi/2
    right = np.pi/2 + (np.pi/2)*(k[1:]/N)**2       # ... pi
    return np.concatenate([left, right]), nw


def refine(ys):
    mid = 0.5*(ys[:-1] + ys[1:])
    out = np.empty(2*len(ys) - 1)
    out[0::2], out[1::2] = ys, mid
    return out


def dof_map(ys, ws):
    """dof index (or -1) and sign for every node (i, j)."""
    ny, nw = len(ys), len(ws)
    ifib = int(np.argmin(np.abs(ys - np.pi/2)))
    assert abs(ys[ifib] - np.pi/2) < 1e-14
    idx = -np.ones((ny, nw), dtype=int)
    sgn = np.ones((ny, nw))
    n = 0
    for i in range(ny - 1):              # last row y = pi is slaved to row 0
        if i == ifib:
            continue
        for j in range(1, nw - 1):
            idx[i, j] = n
            n += 1
    for j in range(nw):                  # (pi, -w_j) = -(0, w_j); -w_j is node nw-1-j
        idx[ny - 1, nw - 1 - j] = idx[0, j]
        sgn[ny - 1, nw - 1 - j] = 1.0 if PLANT else -1.0
    return idx, sgn, n


def tri_integrals(P, wt):
    """For triangle P (3x2 array of (y,w)) return area-weight integrals:
    I_c = int wt(y), I_ic = int 1/wt(y), Mloc = int wt(y) l_i l_j."""
    y0, y1 = P[:, 0].min(), P[:, 0].max()
    A = np.array([[1, P[0, 0], P[0, 1]], [1, P[1, 0], P[1, 1]], [1, P[2, 0], P[2, 1]]])
    C = np.linalg.inv(A)                 # barycentric l_i = C[0,i] + C[1,i] y + C[2,i] w
    yq = 0.5*(y1 - y0)*GY + 0.5*(y1 + y0)
    wyq = 0.5*(y1 - y0)*GYW
    Ic = Iic = 0.0
    M = np.zeros((3, 3))
    for yy, wy in zip(yq, wyq):
        # w-section of the right triangle at height yy
        ws = []
        for a in range(3):
            b = (a + 1) % 3
            ya, yb = P[a, 0], P[b, 0]
            if (ya - yy)*(yb - yy) <= 0 and ya != yb:
                s = (yy - ya)/(yb - ya)
                ws.append(P[a, 1] + s*(P[b, 1] - P[a, 1]))
            elif ya == yb == yy:
                ws += [P[a, 1], P[b, 1]]
        wl, wr = min(ws), max(ws)
        L = wr - wl
        c = weight(yy)
        Ic += wy*c*L
        Iic += wy*L/c
        wq = 0.5*L*GW + 0.5*(wl + wr)
        lam = C[0][None, :] + C[1][None, :]*yy + C[2][None, :]*wq[:, None]
        M += wy*c*0.5*L*(lam.T*GWW) @ lam
    return Ic, Iic, M, C


def assemble(ys, ws):
    idx, sgn, n = dof_map(ys, ws)
    rows, cols, kv, mv = [], [], [], []
    for i in range(len(ys) - 1):
        for j in range(len(ws) - 1):
            nodes = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            for tri in ([0, 1, 2], [0, 2, 3]):
                nd = [nodes[t] for t in tri]
                P = np.array([[ys[a], ws[b]] for a, b in nd])
                Ic, Iic, M, C = tri_integrals(P, weight)
                gy, gw = C[1], C[2]
                # gradients are exact zeros where coordinates coincide; on a triangle
                # with an edge on the fiber the free node has gw = 0 exactly, so the
                # (large) 1/|cos| integral is multiplied by 0.
                K = Ic*np.outer(gy, gy) + Iic*np.outer(gw, gw)
                for a in range(3):
                    da = idx[nd[a]]
                    if da < 0:
                        continue
                    for b in range(3):
                        db = idx[nd[b]]
                        if db < 0:
                            continue
                        s = sgn[nd[a]]*sgn[nd[b]]
                        rows.append(da); cols.append(db)
                        kv.append(s*K[a, b]); mv.append(s*M[a, b])
    Kmat = sps.csr_matrix((kv, (rows, cols)), shape=(n, n))
    Mmat = sps.csr_matrix((mv, (rows, cols)), shape=(n, n))
    return Kmat, Mmat, idx, sgn


def lowest(K, M, k=NEIG):
    vals = spla.eigsh(K, k=k, M=M, sigma=0.0, which="LM", return_eigenvectors=False)
    return np.sort(vals)


def checks(W):
    """Checks on the coarse mesh that can fail; see RETURN.md."""
    ys, nw = base_mesh()
    ws = np.linspace(-W, W, nw + 1)
    K, M, idx, sgn = assemble(ys, ws)
    lines = []
    # (1) mass of the constant-free test: total weighted area of the rectangle = 4W
    #     computed from the element weight integrals over all triangles.
    tot = 0.0
    for i in range(len(ys) - 1):
        for j in range(nw):
            P1 = np.array([[ys[i], ws[j]], [ys[i+1], ws[j]], [ys[i+1], ws[j+1]]])
            P2 = np.array([[ys[i], ws[j]], [ys[i+1], ws[j+1]], [ys[i], ws[j+1]]])
            tot += tri_integrals(P1, weight)[0] + tri_integrals(P2, weight)[0]
    ok = abs(tot - 4*W) < 1e-12
    lines.append(f"{'PASS' if ok else 'FAIL'}  W={W:.6f} weighted area {tot:.15f} equals 4W={4*W:.15f}")
    # (2) seam consistency: interpolate an admissible phi (anti-periodic by construction)
    #     through the dof map and read it back at the slaved row y = pi.
    phi = lambda y, w: np.where(y <= np.pi/2, 1.0, -1.0)*np.abs(np.cos(y))*np.cos(np.pi*w/(2*W))*(1 + 0.3*np.sin(w + y))
    # check that this phi really is anti-periodic: phi(pi,-w) = -phi(0,w)
    assert np.allclose(phi(np.pi, -ws), -phi(0.0, ws))
    u = np.zeros(int(idx.max()) + 1)
    for i in range(len(ys) - 1):
        for j in range(nw + 1):
            if idx[i, j] >= 0:
                u[idx[i, j]] = phi(ys[i], ws[j])
    back = np.array([sgn[-1, j]*u[idx[-1, j]] if idx[-1, j] >= 0 else 0.0 for j in range(nw + 1)])
    err = np.max(np.abs(back - phi(np.pi, ws)))
    lines.append(f"{'PASS' if err < 1e-14 else 'FAIL'}  W={W:.6f} seam dof map reproduces phi(pi,w) (max err {err:.2e})")
    # (3) symmetry and positive definiteness of K, M
    sym = max(abs(K - K.T).max(), abs(M - M.T).max())
    mineig = np.linalg.eigvalsh(M.toarray()).min()
    lines.append(f"{'PASS' if sym < 1e-13 and mineig > 0 else 'FAIL'}  W={W:.6f} K, M symmetric (asym {sym:.1e}) and M positive definite (min eig {mineig:.3e})")
    return lines


def run():
    results, extra = {}, {}
    for lab, W in WIDTHS.items():
        ys, nw = base_mesh()
        ws = np.linspace(-W, W, nw + 1)
        levels = []
        for lev in range(5):
            K, M, _, _ = assemble(ys, ws)
            ev = lowest(K, M)
            levels.append([float(v) for v in ev])
            print(f"W={lab:6s} level {lev} ({len(ys)-1}x{len(ws)-1} cells, {K.shape[0]} dofs): "
                  + " ".join(f"{v:.10f}" for v in ev), flush=True)
            ys, ws = refine(ys), refine(ws)
        L = np.array(levels)
        def est(col):
            a, b, c = L[2, col], L[3, col], L[4, col]
            p = float(np.log2(abs(a - b)/abs(b - c)))
            ex = float(c + (c - b)/(2**p - 1))
            return p, ex, float(abs(ex - c))
        p, ex, er = est(0)
        results[lab] = {"levels": levels, "p": p, "extrap": ex, "err": er}
        extra[lab] = [dict(zip(("p", "extrap", "err"), est(cidx))) for cidx in range(NEIG)]
        print(f"W={lab:6s} p={p:.6f} extrap={ex:.12f} err={er:.3e}", flush=True)
    return results, extra


if __name__ == "__main__":
    if "--checks" in sys.argv:
        for W in (0.25, np.pi/2, 1.7):
            for line in checks(W):
                print(line)
    else:
        res, extra = run()
        with open(os.path.join(HERE, "t4.json"), "w") as fh:
            json.dump({"t4": res, "t4_extra_all_six": extra}, fh, indent=1)
