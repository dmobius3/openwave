"""T4: P1 finite elements in (y, w) for -Delta on the admissible class (R = 1).

Protocol (spec sheet T4):
  rectangle [0, pi] x [-W, W], anti-periodic seam (0, w) ~ (pi, -w) with phi(pi,-w) = -phi(0,w),
  area weight |cos y|, zero data on w = +-W and on the fiber y = pi/2.
  Level 0: y nodes at pi/2 -+ (pi/2)(k/N)^2, k = 0..N, N = 16 per side; 8 uniform cells across in w.
  Four uniform (midpoint / red) refinements -> five levels.
  p from the last three levels, Richardson extrapolation of the bottom from the last two using that p,
  error estimate |lambda_extrap - lambda_finest|.
Readings (see RETURN.md): each (y,w) cell is cut by the diagonal (y0,w0)-(y1,w1); integrals of the
weights are done by 16-point Gauss-Legendre in y and exact (Simpson) integration in w.
Checks (each can fail; see --plant):
  C4.1  sum of all entries of the unconstrained mass matrix = int |cos y| dy dw = 4W   (rel 1e-12)
  C4.2  y-stiffness of phi = y (unconstrained) = 4W                                    (rel 1e-12)
  C4.3  K and M (constrained) symmetric, M positive definite (Cholesky succeeds)
  C4.4  on every dropped divergent element term the free-node w-gradients are exactly zero
  C4.5  finest-level Rayleigh quotient of the interpolated closed-form bottom eigenfunction
        lies within 5% of the closed-form value (a gross-error check of seam and weights;
        not used to tune anything)
Plants: abscos (cos instead of |cos|), seam (periodic seam sign +1), ygrad (drop 1/dy in d/dy),
        asym (perturb one local stiffness entry), fibidx (check the wrong node in C4.4).
Use --only LABEL to run a single width (used for the planted runs).
Writes t4_out.json.
"""
import sys, os, json, math
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
here = os.path.dirname(os.path.abspath(__file__))

WIDTHS = {"1/4": 0.25, "1/2": 0.5, "1": 1.0, "7/5": 1.4, "3/2": 1.5, "pi/2": math.pi / 2, "17/10": 1.7}
N0, M0, NLEV, NEIG = 16, 8, 5, 6
GX, GW = np.polynomial.legendre.leggauss(16)

ok_all = True
fails = []


def report(name, ok):
    global ok_all
    ok_all &= bool(ok)
    if not ok:
        fails.append(name)
    print(f"  {name}: {'PASS' if ok else 'FAIL'}")


def weight(yv):
    return np.cos(yv) if plant == "abscos" else np.abs(np.cos(yv))


def ynodes(level):
    k = np.arange(N0 + 1)
    side = (math.pi / 2) * (k / N0) ** 2          # distance from fiber
    yy = np.concatenate([math.pi / 2 - side[::-1], math.pi / 2 + side[1:]])
    for _ in range(level):                         # uniform refinement = midpoint insertion
        mid = 0.5 * (yy[:-1] + yy[1:])
        new = np.empty(2 * len(yy) - 1)
        new[0::2], new[1::2] = yy, mid
        yy = new
    return yy


def assemble(W, level, constrained=True):
    yy = ynodes(level)
    Mw = M0 * 2**level
    ww = np.linspace(-W, W, Mw + 1)
    ny = len(yy)                                   # indices 0..2N, fiber index N
    ifib = (ny - 1) // 2
    h = ww[1] - ww[0]
    # dof map: node (i, j) -> (dof, sign)
    nfree_y = [i for i in range(ny - 1) if i != ifib]
    pos = {i: r for r, i in enumerate(nfree_y)}
    def dof(i, j):
        if j == 0 or j == Mw or i == ifib:
            return -1, 0.0
        if i == ny - 1:                            # (pi, w_j) ~ (0, -w_j) = (0, w_{Mw-j}), sign -1
            return pos[0] * (Mw - 1) + (Mw - j - 1), (1.0 if plant == "seam" else -1.0)
        return pos[i] * (Mw - 1) + (j - 1), 1.0
    if not constrained:
        def dof(i, j):
            return i * (Mw + 1) + j, 1.0
        ndof = ny * (Mw + 1)
    else:
        ndof = len(nfree_y) * (Mw - 1)
    rows, cols, kv, mv = [], [], [], []
    dropped_ok = True
    for i in range(ny - 1):
        y0, y1 = yy[i], yy[i + 1]
        dy = y1 - y0
        yg = 0.5 * (y0 + y1) + 0.5 * dy * GX
        wg = 0.5 * dy * GW
        t = (yg - y0) / dy                         # in (0,1)
        cy = weight(yg)
        sy = np.abs(np.cos(yg))
        fib_lo, fib_hi = (i == ifib), (i + 1 == ifib)
        ddy = 1.0 if plant == "ygrad" else 1.0 / dy
        for typ in (1, 2):
            # T1: a=(y0,w0) b=(y1,w0) c=(y1,w1); width(y) = h t, w-segment [w0, w0 + h t]
            # T2: a=(y0,w0) d=(y0,w1) c=(y1,w1); width(y) = h (1-t), w-segment [w0 + h t, w1]
            if typ == 1:
                loc = [(0, 0), (1, 0), (1, 1)]
                grads = np.array([[-ddy, 0.0], [ddy, -1 / h], [0.0, 1 / h]])
                width = h * t
                # local coords s = (w-w0)/h on the segment s in [0, t]
                s_lo, s_hi = np.zeros_like(t), t
                def basis(s):
                    return [1 - t, t - s, s]
                wid_at_fiber_nonzero = fib_hi   # width h at y1
            else:
                loc = [(0, 0), (0, 1), (1, 1)]
                grads = np.array([[0.0, -1 / h], [-ddy, 1 / h], [ddy, 0.0]])
                width = h * (1 - t)
                s_lo, s_hi = t, np.ones_like(t)
                def basis(s):
                    return [1 - s, s - t, t]
                wid_at_fiber_nonzero = fib_lo   # width h at y0
            Ic = np.sum(wg * cy * width)
            if wid_at_fiber_nonzero:
                Is = 0.0                       # divergent integral; must multiply zero free gradients
                free_on_fiber_edge = [0] if typ == 1 else [2]
                if plant == "fibidx":
                    free_on_fiber_edge = [1]
                dropped_ok &= all(grads[q, 1] == 0.0 for q in free_on_fiber_edge)
            else:
                Is = np.sum(wg * width / sy)
            # mass: Simpson in w (exact for quadratics), Gauss in y
            sm = 0.5 * (s_lo + s_hi)
            B = [basis(s_lo), basis(sm), basis(s_hi)]
            Mloc = np.zeros((3, 3))
            for p_ in range(3):
                for q_ in range(3):
                    inner = width * (B[0][p_] * B[0][q_] + 4 * B[1][p_] * B[1][q_] + B[2][p_] * B[2][q_]) / 6
                    Mloc[p_, q_] = np.sum(wg * cy * inner)
            Kloc = Ic * np.outer(grads[:, 0], grads[:, 0]) + Is * np.outer(grads[:, 1], grads[:, 1])
            if plant == "asym":
                Kloc[0, 1] *= 1.01
            for j in range(Mw):
                ds = [dof(i + di, j + dj) for (di, dj) in loc]
                for p_ in range(3):
                    dp, sp_ = ds[p_]
                    if dp < 0:
                        continue
                    for q_ in range(3):
                        dq, sq_ = ds[q_]
                        if dq < 0:
                            continue
                        rows.append(dp); cols.append(dq)
                        kv.append(sp_ * sq_ * Kloc[p_, q_]); mv.append(sp_ * sq_ * Mloc[p_, q_])
    K = sps.csc_matrix((kv, (rows, cols)), shape=(ndof, ndof))
    M = sps.csc_matrix((mv, (rows, cols)), shape=(ndof, ndof))
    return K, M, yy, ww, dof, dropped_ok


def exact_bottom_fn(W):
    nu = math.pi / (2 * W)
    def phi(yv, wv):                           # y-picture of the lune function cos(psi)^nu sin(nu(w+W))
        if yv <= math.pi / 2:
            psi, ww_ = yv, wv
            sgn = 1.0
        else:
            psi, ww_ = yv - math.pi, -wv
            sgn = -1.0
        return sgn * abs(math.cos(psi)) ** nu * math.sin(nu * (ww_ + W))
    return nu * (nu + 1), phi


results = {}
for lab, W in WIDTHS.items():
    if only and lab != only:
        continue
    print(f"W = {lab}  ({W!r})")
    # unconstrained consistency checks at level 0
    Ku, Mu, yy, ww, _, _ = assemble(W, 0, constrained=False)
    tot = Mu.sum()
    report(f"C4.1 [W={lab}] mass total {tot:.15f} vs 4W={4*W:.15f}", abs(tot - 4 * W) < 1e-12 * 4 * W)
    Mw = len(ww) - 1
    vy = np.repeat(yy, Mw + 1)
    # y-stiffness of phi=y: phi has zero w-gradient, so the full K applies (Is terms vanish)
    Ky = vy @ (Ku @ vy)
    report(f"C4.2 [W={lab}] y-stiffness of phi=y {Ky:.15f} vs 4W", abs(Ky - 4 * W) < 1e-12 * 4 * W)
    levels, sizes = [], []
    for lev in range(NLEV):
        K, M, yy, ww, dof, dropped_ok = assemble(W, lev)
        if lev == 0:
            symm = (abs(K - K.T).max() == 0) and (abs(M - M.T).max() == 0)
            try:
                np.linalg.cholesky(M.toarray()); spd = True
            except np.linalg.LinAlgError:
                spd = False
            report(f"C4.3 [W={lab}] K, M symmetric and M SPD", symm and spd)
        report(f"C4.4 [W={lab}, level {lev}] dropped divergent terms meet zero gradients", dropped_ok)
        vals = spla.eigsh(K, k=NEIG, M=M, sigma=0.0, which="LM", return_eigenvectors=False)
        vals = np.sort(vals)
        levels.append([float(v) for v in vals])
        sizes.append(K.shape[0])
        print(f"   level {lev} ndof={K.shape[0]:6d} :", " ".join(f"{v:.10f}" for v in vals))
        if lev == NLEV - 1:
            lam_ex, phi = exact_bottom_fn(W)
            Mw = len(ww) - 1
            v = np.zeros(K.shape[0])
            for i in range(len(yy)):
                for j in range(Mw + 1):
                    d, s = dof(i, j)
                    if d >= 0 and i < len(yy) - 1:
                        v[d] = s * phi(yy[i], ww[j])
            rq = (v @ (K @ v)) / (v @ (M @ v))
            report(f"C4.5 [W={lab}] interpolant Rayleigh quotient {rq:.10f} vs closed form {lam_ex:.10f}",
                   abs(rq - lam_ex) < 0.05 * lam_ex)
    l3, l4, l5 = levels[2][0], levels[3][0], levels[4][0]
    d1, d2 = l3 - l4, l4 - l5
    ratio = d1 / d2 if d2 != 0 else float("nan")
    p = math.log2(ratio) if ratio > 0 else float("nan")
    extrap = l5 + (l5 - l4) / (2**p - 1) if p == p else float("nan")
    err = abs(extrap - l5) if extrap == extrap else float("nan")
    print(f"   p = {p:.6f}  extrap = {extrap:.12f}  err = {err:.3e}  (ratio {ratio:.6f})")
    results[lab] = {"levels": levels, "p": p, "extrap": extrap, "err": err, "ndof": sizes, "ratio": ratio}

print("T4_ALL", "PASS" if ok_all else "FAIL", ("failed: " + "; ".join(fails)) if fails else "")
if plant is None and only is None:
    with open(os.path.join(here, "t4_out.json"), "w") as fh:
        json.dump(results, fh, indent=1)
