"""M5.32 R29-1: the vacuum gap of the projected kinetic Hessian about the relaxed hedgehog (the D8 row), with the
two known-answer controls of R29-0 (d).

EQUATIONS FIRST
---------------
Static energy E[M] of the R23-1 instrument at c = 0 (the R22-1 lattice: certified quartic + V4 with the targets
(-g, 1, delta, delta), g = 8, delta = 0.3, w = 25 W1, h = 1.5, the shell pinned at the seed), M_0i = 0 and M_00
slaved on the relaxed field; the full symmetric gradient G = dE/dM over all ten entries per cell.
A small oscillation M = M* + eps(t) D about the relaxed field M* has, under the certified action,
    E = E* + 1/2 S2(D) eps^2 + K(D) eps_dot^2,
    S2(D) = D . Hess E . D   (central difference of the gradient, eps scaled to |D|_max),
    K(D)  = 4 h^3 sum_x sum_br wt sum_i <[D, A_i]_eta, [D, A_i]_eta>_eta,   A_i = d_i M*   (R4's kin at lambda = 0),
so omega^2 = S2 / (2 K).  K is LOCAL in D (no derivative of D): per free cell it is a 10 x 10 symmetric form k_x,
whose eigenvectors U_x split the ten directions into K-positive (K+), K-negative (K-) and K-null (static-stiff).
In the K-eigenbasis coordinates D = U y the weight is diagonal, B = diag(kappa_x) > 0 on K+, so the lowest
omega^2 on K+ is the smallest eigenvalue of the pencil (U+^T H U+, B+), found by LOBPCG with the generator span
projected out (the constraint Y). On K- the GHOST word asks whether S2 can be negative: the smallest eigenvalue of
U-^T H U- (Euclidean) is computed the same way. The trace direction eta is in the K-null space of every cell (R29-0).
Pre-registered words (the task doc): GAPPED / GAPLESS / STATIC-STIFF / GHOST, read on the 1/L fit across rungs.

CONTROLS (R29-0 (d); the same solver, a known functional):
  (i) the free massive field: E_c[D] = 1/2 sum_x [ sum_i q_h(d_i D) + m^2 q_h(D) ] with q_h = K_h's form at the
      vacuum (positive on all ten directions) and the same q_h as the omega^2 weight, so every mode has
      omega^2 = k^2 + m^2 exactly in the continuum; in the pinned box the lowest mode sits at k ~ pi/L_free:
      m^2 = 1 makes that 0.4 percent at n 36 (the 2 percent bar);
  (ii) the Goldstone control: E_g[M] = V4(M) + 1/2 sum_x sum_i q_h(d_i M) about the uniform R23 vacuum, omega^2
      weight q_h: the internal rotation modes of V4 are flat, so their lowest omega^2 = (pi/L_free)^2 -> 0 in 1/L.

Modes (python3 m5_32_r29_1_gap.py MODE ...):
    census  TAG                 the K census of a stored field (counts per cell, the trace mode)
    gap     TAG [--nev 6]       the pencil on a stored relaxed field (K+ lowest, K- lowest, generators, static-stiff)
    control-mass n L            control (i) on an n^3 box of side L
    control-gold n L            control (ii)
    relax   n L [--from TAG]    relax a c = 0 hedgehog on a new rung through the R23-1 instrument (the pinned rad
                                seed, or the interpolated stored field TAG as the warm start), under the R23-1 gate
    collect                     the 1/L fit over the rungs and the word
Set OMP_NUM_THREADS=1 (one thread per job; the pool shape is the machine's, never written here).
Out: ../data/m5_32_r29_1/*.json (per mode), ../data/m5_32_r29_1_gap.json (collect); fields in ../data/m5_32_r23_1/.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.ndimage import zoom
from scipy.sparse.linalg import LinearOperator, eigsh, lobpcg

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_DIR = os.path.join(DATA, "m5_32_r29_1")
OUT_JSON = os.path.join(DATA, "m5_32_r29_1_gap.json")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


C23 = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, R0, W1, G = C23.R21, C23.R20, C23.R0, C23.W1, C23.G
B3 = R21.B3
DELTA = 0.3
W1S = 25.0
L_OVERRIDE = None                                    # set from --L for a rung whose spacing is not 1.5


def sym_basis():
    out = []
    for a in range(4):
        X = np.zeros((4, 4)); X[a, a] = 1.0; out.append(X)
    for j in range(1, 4):
        X = np.zeros((4, 4)); X[0, j] = X[j, 0] = 1.0; out.append(X)
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        X = np.zeros((4, 4)); X[i, j] = X[j, i] = 1.0; out.append(X)
    return np.array(out)


BX = sym_basis()                                    # (10, 4, 4); Frobenius norms 1 (diag) and sqrt 2 (off)
BW = np.array([np.sum(b * b) for b in BX])          # the metric of the basis coordinates


# ================= the lattice objects =================
class Rung:
    """one box: cfg, the free mask, the energy / gradient, the K weight, the K eigenbasis."""

    def __init__(self, n, L, M=None, functional="cert", m2=1.0):
        self.n, self.L = n, float(L)
        self.cfg = R21.cfg_of(n, L, G, DELTA)
        self.h = self.cfg["h"]
        self.p = R21.params_of(G, DELTA)
        self.pot = ("v4", R0.roots_of(self.cfg, degenerate=True), W1 * W1S)
        self.free = ~B3.pin_shell(n, self.h)
        self.nfree = int(self.free.sum())
        self.functional, self.m2 = functional, m2
        self.M = M
        self.X, self.Y, self.Z = B3.coords(n, self.h)
        self.r = np.sqrt(self.X ** 2 + self.Y ** 2 + self.Z ** 2)
        self.ngrad = 0

    # ---- energy and gradient of the three functionals ----
    def grad(self, M):
        self.ngrad += 1
        if self.functional == "cert":
            E, Gm, info = C23.energy_grad(M, self.cfg, self.p, self.pot, 0.0)
            return E, Gm
        # the control functionals: quadratic q_h gradient energy (+ V4 for the Goldstone control, + m^2 q_h for the mass control)
        h3, h = self.h ** 3, self.h
        E, Gm = 0.0, np.zeros_like(M)
        for br, wt in B3.branches(self.cfg["stencil"]):
            for ax in range(3):
                A = B3.d1(M, ax, h, br)
                E += wt * 0.5 * h3 * float(np.sum(A * A))             # q_h(X) = tr(X X) at the vacuum (h = 1)
                Gm += wt * h3 * B3.d1_adj(A, ax, h, br)
        if self.functional == "mass":
            E += 0.5 * h3 * self.m2 * float(np.sum(M * M))
            Gm += h3 * self.m2 * M
        else:                                                           # "gold": V4 about the vacuum
            ev, gv = R0.v4_energy_grad(M, self.cfg, self.pot[1], self.pot[2], True)
            E += float(ev); Gm += gv
        return E, Gm

    # ---- the omega^2 weight (certified: R4's kin at lambda 0; controls: q_h) ----
    def K_cell_matrices(self):
        """k_x[a, b] on the free cells, (nfree, 10, 10)."""
        n = self.n
        if self.functional != "cert":
            k = np.zeros((self.nfree, 10, 10))
            k[:] = np.diag(0.5 * self.h ** 3 * BW)                      # 1/2 h^3 q_h(X) = 1/2 h^3 |X|_F^2
            return k
        Ab = {br: [B3.d1(self.M, ax, self.h, br) for ax in range(3)] for br, _ in B3.branches(self.cfg["stencil"])}

        def kc(Xc):
            tot = np.zeros((n, n, n))
            for br, wt in B3.branches(self.cfg["stencil"]):
                for i in range(3):
                    F = B3.comm_eta(Xc, Ab[br][i])
                    tot += wt * B3.inner_eta(F, F)
            return 4.0 * self.h ** 3 * tot
        kx = np.zeros((n, n, n, 10, 10))
        diag = [kc(np.broadcast_to(BX[a], self.M.shape)) for a in range(10)]
        for a in range(10):
            kx[..., a, a] = diag[a]
            for b in range(a + 1, 10):
                kx[..., a, b] = kx[..., b, a] = 0.5 * (kc(np.broadcast_to(BX[a] + BX[b], self.M.shape)) - diag[a] - diag[b])
        return kx[self.free]

    def build_basis(self, null_tol=1e-8):
        """the K eigenbasis per free cell in ORTHONORMAL coordinates (the basis metric BW folded in)."""
        k = self.K_cell_matrices()                                      # the form in basis coordinates x
        sw = np.sqrt(BW)
        ko = k / sw[None, :, None] / sw[None, None, :]                  # the form in orthonormal coordinates y = sw x
        ev, U = np.linalg.eigh(ko)                                      # (nfree, 10), (nfree, 10, 10)
        kscale = float(np.max(np.abs(ev)))
        tol = null_tol * kscale
        self.kev, self.kU, self.kscale = ev, U, kscale
        self.pos = ev > tol
        self.neg = ev < -tol
        self.nul = ~(self.pos | self.neg)
        self.sw = sw
        return {"K_scale": kscale, "null_tol": tol, "n_pos": int(self.pos.sum()), "n_neg": int(self.neg.sum()), "n_null": int(self.nul.sum()),
                "per_cell_signature": sorted({f"{int(a)}+/{int(b)}-/{int(c)}0" for a, b, c in zip(self.pos.sum(1), self.neg.sum(1), self.nul.sum(1))})}

    # ---- the block-Jacobi preconditioner: the exact per-cell Hessian of V4 (polarisation of the local density) plus
    #      the gradient term's diagonal scale, in the K-eigen coordinates, inverted per cell on the selected subspace ----
    def v4_cell_hessian(self):
        n = self.n
        if self.functional == "mass":
            return np.broadcast_to(np.diag(self.h ** 3 * self.m2 * BW), (self.nfree, 10, 10)).copy()
        roots, w = self.pot[1], self.pot[2]
        cp = [sum(q ** k for q in roots) for k in range(1, 5)]

        def dens(M):
            N = M @ ETA; tot = np.zeros(M.shape[:3]); P = np.broadcast_to(np.eye(4), M.shape).copy()
            for k in range(1, 5):
                P = P @ N; tot += (np.einsum("...kk->...", P) - cp[k - 1]) ** 2
            return w * self.h ** 3 * tot
        M0 = self.M; e = 1e-3
        Hc = np.zeros((n, n, n, 10, 10))
        f0 = dens(M0)
        fa = [(dens(M0 + e * BX[a]) + dens(M0 - e * BX[a]) - 2 * f0) / (e * e) for a in range(10)]
        for a in range(10):
            Hc[..., a, a] = fa[a]
            for b in range(a + 1, 10):
                fab = (dens(M0 + e * (BX[a] + BX[b])) + dens(M0 - e * (BX[a] + BX[b])) - 2 * f0) / (e * e)
                Hc[..., a, b] = Hc[..., b, a] = 0.5 * (fab - fa[a] - fa[b])
        return Hc[self.free]

    def build_preconditioner(self, sel):
        Hc = self.v4_cell_hessian()                                     # basis coordinates
        lap = 6.0 * self.h / 1.0                                        # the gradient term's diagonal scale (h^3 x 6/h^2 per unit |X|_F^2) times the term's weight
        wgt = 1.0 if self.functional != "cert" else 4.0 * 2.0 * float(np.mean(np.abs(self.kev)))   # a scale guess for the certified curvature term
        Hc = Hc + (lap * wgt) * np.diag(BW)[None]
        Ho = Hc / self.sw[None, :, None] / self.sw[None, None, :]
        Hk = np.einsum("fij,fjk,fkl->fil", self.kU.transpose(0, 2, 1), Ho, self.kU)   # K-eigen coordinates
        # the selected block per cell, regularised and inverted
        selc = sel                                                      # (nfree, 10) boolean
        inv = np.zeros_like(Hk)
        for f in range(Hk.shape[0]):
            idx = np.where(selc[f])[0]
            if len(idx) == 0:
                continue
            B = Hk[f][np.ix_(idx, idx)]
            B = 0.5 * (B + B.T) + 1e-6 * np.eye(len(idx)) * max(1.0, np.abs(B).max())
            wB, VB = np.linalg.eigh(B)
            wB = np.maximum(np.abs(wB), 1e-6 * np.abs(wB).max())
            inv[f][np.ix_(idx, idx)] = (VB / wB) @ VB.T
        self.precond_cells = inv

        def apply(z):
            zz = np.zeros(self.kU.shape[:2]); zz[sel] = z
            out = np.einsum("fij,fj->fi", self.precond_cells, zz)
            return out[sel]
        return apply

    # ---- maps between a field direction D (n,n,n,4,4) and the coordinates on a K-subspace ----
    def field_to_coords(self, D, sel):
        """D -> the K-eigenbasis coordinates on the selected (per-cell boolean) directions; the complement is dropped."""
        x = np.einsum("fab,cab->fc", D[self.free], BX) / BW             # basis coordinates per free cell (f, 10)
        y = x * self.sw                                                 # orthonormal
        z = np.einsum("fij,fi->fj", self.kU, y)                         # K-eigen coordinates
        return z[sel]

    def coords_to_field(self, zsel, sel):
        z = np.zeros(self.kU.shape[:2]); z[sel] = zsel
        y = np.einsum("fij,fj->fi", self.kU, z)
        x = y / self.sw
        D = np.zeros((self.n, self.n, self.n, 4, 4))
        D[self.free] = np.einsum("fc,cab->fab", x, BX)
        return D

    # ---- the second variation as a matrix-free operator on a subspace ----
    def S2_apply(self, D, eps=1e-3):
        s = eps / max(float(np.max(np.abs(D))), 1e-300)
        _, Gp = self.grad(self.M + s * D); _, Gm = self.grad(self.M - s * D)
        return (Gp - Gm) / (2 * s)

    def S2_apply_richardson(self, D, eps=1e-3):
        a = self.S2_apply(D, eps); b = self.S2_apply(D, eps / 2)
        return (4 * b - a) / 3.0, float(np.max(np.abs(a - b)) / max(1e-300, float(np.max(np.abs(b)))))

    def quad(self, D):
        """S2(D) = D . H D (two gradients), K(D), and the ratio."""
        HD = self.S2_apply(D)
        s2 = float(np.sum(D * HD))
        kd = self.K_of(D)
        return s2, kd

    def K_of(self, D):
        x = np.einsum("fab,cab->fc", D[self.free], BX) / BW
        y = x * self.sw
        z = np.einsum("fij,fi->fj", self.kU, y)
        return float(np.sum(self.kev * z * z))

    def smooth_starts(self, sel, nev, seed=29):
        """nev smooth random field directions (low-k Fourier sums times a Gaussian envelope on the free cells), in the
        coordinates of sel: the lowest modes are long-wavelength exterior textures, which white noise reaches slowly."""
        rng = np.random.default_rng(seed)
        kmax = 3 * 2 * np.pi / self.L
        env = np.exp(-(self.r / (0.35 * self.L)) ** 2)
        free = self.free.astype(float)[..., None, None]
        cols = []
        for _ in range(nev):
            D = np.zeros((self.n, self.n, self.n, 4, 4))
            for _m in range(6):
                k = rng.uniform(-kmax, kmax, 3); ph = rng.uniform(0, 2 * np.pi)
                S = rng.standard_normal((4, 4)); S = 0.5 * (S + S.T)
                D += np.cos(k[0] * self.X + k[1] * self.Y + k[2] * self.Z + ph)[..., None, None] * S
            D *= env[..., None, None] * free
            cols.append(self.field_to_coords(D, sel))
        X0 = np.array(cols).T
        X0 += 1e-3 * rng.standard_normal(X0.shape) * np.abs(X0).max()
        return X0

    def null_block_inverse(self):
        """the per-cell inverse of the Hessian on the K-null directions (the exact per-cell V4 Hessian plus the gradient
        term's diagonal scale, in the K-eigen coordinates): the auxiliary-variable elimination of the author's note of
        2026-10-28 ((K_null)^-1 H_null-null, done one step block-Jacobi; the null block is 1e5 times stiffer than the
        soft modes on this field, so the neglected inter-cell couplings are a second-order correction)."""
        if getattr(self, "_null_inv", None) is not None:
            return self._null_inv
        Hc = self.v4_cell_hessian()
        lap = 6.0 * self.h; wgt = 1.0 if self.functional != "cert" else 4.0 * 2.0 * float(np.mean(np.abs(self.kev)))
        Hc = Hc + (lap * wgt) * np.diag(BW)[None]
        Ho = Hc / self.sw[None, :, None] / self.sw[None, None, :]
        Hk = np.einsum("fij,fjk,fkl->fil", self.kU.transpose(0, 2, 1), Ho, self.kU)
        inv = np.zeros_like(Hk)
        for f in range(Hk.shape[0]):
            idx = np.where(self.nul[f])[0]
            if len(idx) == 0:
                continue
            B = Hk[f][np.ix_(idx, idx)]; B = 0.5 * (B + B.T)
            inv[f][np.ix_(idx, idx)] = np.linalg.pinv(B, rcond=1e-10)
        self._null_inv = inv
        return inv

    def S2_apply_schur(self, D, eps=1e-3):
        """H_eff D = (H (D + U0 z))  with z = -(H00)^-1 (H D)_0 : the K-null components respond to D through their own
        stationarity equation instead of being held at zero. Four gradients instead of two."""
        HD = self.S2_apply(D, eps)
        inv = self.null_block_inverse()
        x = np.einsum("fab,cab->fc", HD[self.free], BX) / BW
        yv = x * self.sw
        zc = np.einsum("fij,fi->fj", self.kU, yv)                       # (HD) in K-eigen coordinates
        z0 = -np.einsum("fij,fj->fi", inv, np.where(self.nul, zc, 0.0))   # the null response
        Dz = self.coords_to_field(z0[self.nul], self.nul)
        self.last_schur_norm = float(np.sqrt(np.sum(Dz * Dz)) / max(1e-300, np.sqrt(np.sum(D * D))))
        return self.S2_apply(D + Dz, eps)

    def S2_schur_cg(self, D, eps=1e-3, iters=30, tol=1e-6):
        """the K-null auxiliary elimination with the TRUE null block: solve H00 z = -(H D)_0 by preconditioned conjugate
        gradients (matrix-free, two gradients per step, the per-cell block inverse as the preconditioner), then
        S2_eff = D . H (D + U0 z). Returns (S2_eff, cg_residual_ratio, steps, |z|/|D|)."""
        HD = self.S2_apply(D, eps)
        inv = self.null_block_inverse()
        def to_null(F):
            x = np.einsum("fab,cab->fc", F[self.free], BX) / BW
            return np.einsum("fij,fi->fj", self.kU, x * self.sw)[self.nul]
        def from_null(zn):
            return self.coords_to_field(zn, self.nul)
        def H00(zn):
            return to_null(self.S2_apply(from_null(zn), eps))
        def Pinv(rn):
            zz = np.zeros(self.kU.shape[:2]); zz[self.nul] = rn
            return np.einsum("fij,fj->fi", inv, zz)[self.nul]
        b = -to_null(HD)
        z = Pinv(b); r = b - H00(z); y = Pinv(r); pdir = y.copy(); ry = float(r @ y); b2 = float(np.sqrt(b @ b)); steps = 0
        for steps in range(1, iters + 1):
            Ap = H00(pdir); alpha = ry / float(pdir @ Ap)
            z = z + alpha * pdir; r = r - alpha * Ap
            if np.sqrt(float(r @ r)) < tol * b2:
                break
            y = Pinv(r); ry_new = float(r @ y); pdir = y + (ry_new / ry) * pdir; ry = ry_new
        Dz = from_null(z)
        s2 = float(np.sum(D * self.S2_apply(D + Dz, eps)))
        return s2, float(np.sqrt(float(r @ r)) / b2), steps, float(np.sqrt(np.sum(Dz * Dz)) / np.sqrt(np.sum(D * D)))

    def pencil(self, sel, nev, Y_fields=(), maxiter=400, tol=1e-5, seed=29, label="", X0=None, warm=None, schur=False):
        """the smallest eigenvalues of (U_sel^T H U_sel, B_sel) by LOBPCG; B = the K eigenvalues on sel (sign-corrected)."""
        dim = int(sel.sum())
        bvals = self.kev[sel]
        sign = 1.0 if bvals.min() > 0 else -1.0                          # K- : use -K as the (positive) metric for the GHOST read
        Bdiag = sign * bvals
        S2fun = self.S2_apply_schur if schur else self.S2_apply
        A = LinearOperator((dim, dim), matvec=lambda z: self.field_to_coords(S2fun(self.coords_to_field(z, sel)), sel), dtype=float)
        if label == "ghost":                                             # GHOST: the Euclidean minimum of S2 on K- (B = I)
            Bop = None
        else:
            Bop = LinearOperator((dim, dim), matvec=lambda z: Bdiag * z, dtype=float)
        if warm is not None and os.path.exists(warm):
            Z = np.load(warm)["modes"]
            if Z.shape[1] != self.n:                                   # a smaller rung's modes: interpolate (order 1) to this grid
                f = self.n / Z.shape[1]
                Z = np.array([np.stack([np.stack([zoom(Zi[..., a, b], f, order=1) for b in range(4)], -1) for a in range(4)], -2) for Zi in Z])
                Z *= self.free.astype(float)[None, ..., None, None]
            X0 = np.array([self.field_to_coords(Z[i], sel) for i in range(min(nev, Z.shape[0]))]).T
            if X0.shape[1] < nev:
                X0 = np.hstack([X0, self.smooth_starts(sel, nev - X0.shape[1], seed)])
        elif X0 is None:
            X0 = self.smooth_starts(sel, nev, seed)
        Mpre = LinearOperator((dim, dim), matvec=self.build_preconditioner(sel), dtype=float)
        Yc = None
        if Y_fields:
            Yc = np.array([self.field_to_coords(D, sel) for D in Y_fields]).T
            # keep only the independent constraints
            u, s, _ = np.linalg.svd(Yc, full_matrices=False)
            Yc = u[:, s > 1e-8 * s[0]]
        t0 = time.time(); self.ngrad = 0
        if label == "gap-arpack":
            # Lanczos on the generalised pencil (M = B diagonal, Minv explicit), the constraint span projected out of A
            Pc = (lambda v: v) if Yc is None else (lambda v, Yq=np.linalg.qr(Yc)[0]: v - Yq @ (Yq.T @ v))
            Aop = LinearOperator((dim, dim), matvec=lambda v: Pc(A.matvec(Pc(v))), dtype=float)
            Minv = LinearOperator((dim, dim), matvec=lambda v: v / Bdiag, dtype=float)
            v0 = X0[:, 0] / np.linalg.norm(X0[:, 0])
            vals, vecs = eigsh(Aop, k=nev, M=Bop, Minv=Minv, which="SA", v0=v0, ncv=max(4 * nev, 24), maxiter=maxiter, tol=tol)
            hist = [vals.copy(), vals.copy()]; n_iter = -1
        elif label == "ghost":
            # ARPACK Lanczos for the smallest algebraic eigenvalues of the indefinite Euclidean problem on K- (LOBPCG ended
            # worse than its start on this spectrum, which spans -0.5 to 1e5); the constraint span is projected out of
            # the operator (P A P with P = 1 - Y Y^T), the start vector is the best saved witness or a smooth direction.
            if Yc is not None:
                Yq, _ = np.linalg.qr(Yc)
                Pc = lambda v: v - Yq @ (Yq.T @ v)
                Aop = LinearOperator((dim, dim), matvec=lambda v: Pc(A.matvec(Pc(v))), dtype=float)
            else:
                Aop = A
            v0 = X0[:, 0] / np.linalg.norm(X0[:, 0])
            vals, vecs = eigsh(Aop, k=nev, which="SA", v0=v0, ncv=max(4 * nev, 20), maxiter=maxiter, tol=tol)
            hist = [vals.copy(), vals.copy()]; n_iter = -1
        else:
            vals, vecs, hist = lobpcg(A, X0, B=Bop, M=Mpre, Y=Yc, largest=False, maxiter=maxiter, tol=tol, retLambdaHistory=True)
            n_iter = len(hist)
        order = np.argsort(vals); vals, vecs = vals[order], vecs[:, order]
        out = []
        for i in range(nev):
            D = self.coords_to_field(vecs[:, i], sel)
            D /= np.sqrt(np.sum(D * D))
            HD, rich = self.S2_apply_richardson(D)
            s2 = float(np.sum(D * HD)); kd = self.K_of(D)
            s2_schur = float(np.sum(D * self.S2_apply_schur(D)))        # the auxiliary-eliminated value beside the bare one, always
            # the mode's profile: the radial distribution of |D|^2 and its dominant wavelength
            w2 = np.sum(D * D, axis=(-1, -2))
            rbins = np.linspace(0, 0.5 * self.L, 12)
            prof = [float(w2[(self.r >= rbins[j]) & (self.r < rbins[j + 1])].sum()) for j in range(len(rbins) - 1)]
            out.append({"eigenvalue": float(vals[i]), "S2": s2, "S2_schur_(K-null eliminated)": s2_schur, "null_response_norm_ratio": getattr(self, "last_schur_norm", None),
                        "K": kd, "omega2_S2_over_2K": s2 / (2 * kd) if abs(kd) > 1e-300 else None,
                        "omega2_schur": s2_schur / (2 * kd) if abs(kd) > 1e-300 else None,
                        "richardson_rel_err_HD": rich, "radial_weight_profile_(bins to L/2)": prof,
                        "r_weighted_mean": float(np.sum(w2 * self.r) / w2.sum())})
        conv = [float(abs(hist[-1][i] - hist[-2][i]) / max(1e-300, abs(hist[-1][i]))) for i in range(nev)] if len(hist) > 1 else None
        self.last_mode_fields = np.array([self.coords_to_field(vecs[:, i], sel) for i in range(nev)])
        h10 = [float(hist[j][0]) for j in range(0, len(hist), max(1, len(hist) // 20))]
        return {"dim": dim, "nev": nev, "schur_in_the_iteration": bool(schur), "wall_s": round(time.time() - t0, 1), "gradients": self.ngrad, "iterations": n_iter,
                "lowest_eigenvalue_history_(every_1/20th_iteration)": h10,
                "last_step_rel_change": conv, "converged_(all_last_steps_below_1e-4)": bool(conv is not None and max(conv) < 1e-4),
                "modes": out, "constraints": 0 if Yc is None else int(Yc.shape[1])}


# ================= generators on a field =================
def generators(R: Rung):
    M, h, free = R.M, R.h, R.free.astype(float)[..., None, None]
    dc = lambda ax: 0.5 * (B3.d1(M, ax, h, "fwd") + B3.d1(M, ax, h, "bwd"))
    A = [dc(ax) for ax in range(3)]
    co = [R.X, R.Y, R.Z]
    g = {}
    for k in range(3):
        g[f"translation_{k}"] = A[k] * free
        J = np.zeros((4, 4)); i, j = [(2, 3), (3, 1), (1, 2)][k]; J[i, j], J[j, i] = -1.0, 1.0
        i0, j0 = [(1, 2), (2, 0), (0, 1)][k]
        g[f"orbital_rotation_{k}"] = (co[i0][..., None, None] * A[j0] - co[j0][..., None, None] * A[i0]) * free
        g[f"internal_rotation_{k}"] = (J @ M - M @ J) * free
        Kb = np.zeros((4, 4)); Kb[0, k + 1] = Kb[k + 1, 0] = 1.0
        g[f"internal_boost_{k}"] = (Kb @ M + M @ Kb) * free
    return g


# ================= modes =================
def mode_census(tag):
    M = np.load(os.path.join(C23.OUT_NPZ, tag + ".npz"))["M"]
    n = M.shape[0]; L = L_OVERRIDE if L_OVERRIDE else n * 1.5
    R = Rung(n, L, M)
    cen = R.build_basis()
    # the trace mode in the null space: the null eigenvectors' span contains eta in orthonormal coordinates
    eta_y = (np.einsum("ab,cab->c", ETA, BX) / BW) * R.sw
    zz = np.einsum("fij,i->fj", R.kU, eta_y)
    outside_null = float(np.max(np.abs(zz[~R.nul]))) if (~R.nul).any() else 0.0
    cen["eta_component_outside_null_space_max"] = outside_null
    cen["field"] = tag
    return cen


def mode_gap(tag, nev, constrain=True, warm_dir=None, maxiter=400, schur=False):
    M = np.load(os.path.join(C23.OUT_NPZ, tag + ".npz"))["M"]
    n = M.shape[0]; L = L_OVERRIDE if L_OVERRIDE else n * 1.5
    R = Rung(n, L, M)
    E0, G0 = R.grad(M)
    rec = {"field": tag, "n": n, "L": L, "E": E0, "max|G| spatial free": float(np.abs(G0[..., 1:, 1:] * R.free[..., None, None]).max()),
           "census": R.build_basis()}
    gens = generators(R)
    grows = {}
    for name, D in gens.items():
        nn = float(np.sqrt(np.sum(D * D)))
        Dn = D / nn
        HD, rich = R.S2_apply_richardson(Dn)
        grows[name] = {"S2": float(np.sum(Dn * HD)), "richardson_rel_err": rich, "K": R.K_of(Dn)}
    rec["generators"] = grows
    log(f"{tag}: generators done ({R.ngrad} gradients)")
    # K+ : the gap pool with the generator span projected out
    Yf = list(gens.values()) if constrain else []
    wp = (warm_dir if warm_dir.endswith(".npz") else os.path.join(warm_dir, f"{tag}_Kpos_modes.npz")) if warm_dir else None
    rec["generator_span_projected_out"] = bool(constrain)
    rec["K_positive"] = R.pencil(R.pos, nev, Y_fields=Yf, label=("gap-arpack" if os.environ.get("R29_ARPACK") else "gap"), warm=wp, maxiter=maxiter, schur=schur)
    np.savez_compressed(os.path.join(OUT_DIR, f"{tag}_Kpos_modes.npz"), modes=R.last_mode_fields)
    log(f"{tag}: K+ pencil done, lowest omega2 {rec['K_positive']['modes'][0]['omega2_S2_over_2K']}")
    # K- : the GHOST read (the Euclidean minimum of S2 on K-)
    wn = (warm_dir.replace("_Kpos_", "_Kneg_") if warm_dir.endswith(".npz") else os.path.join(warm_dir, f"{tag}_Kneg_modes.npz")) if warm_dir else None
    rec["K_negative_ghost_read"] = R.pencil(R.neg, min(nev, 4), Y_fields=Yf, label="ghost", warm=wn, maxiter=300, tol=1e-4)
    np.savez_compressed(os.path.join(OUT_DIR, f"{tag}_Kneg_modes.npz"), modes=R.last_mode_fields)
    log(f"{tag}: K- read done, lowest S2 {rec['K_negative_ghost_read']['modes'][0]['S2']}")
    rec["static_stiff_directions"] = rec["census"]["n_null"]
    rec["S2_scale"] = None
    rec["wall_s"] = round(time.time() - T0, 1)
    return rec


def mode_control(kind, n, L):
    if kind == "mass":
        R = Rung(n, L, np.zeros((n, n, n, 4, 4)), functional="mass", m2=1.0)
    else:
        Mv = np.diag([G, 1.0, DELTA, DELTA])
        R = Rung(n, L, np.broadcast_to(Mv, (n, n, n, 4, 4)).copy(), functional="gold")
    cen = R.build_basis()
    rec = {"control": kind, "n": n, "L": L, "census": cen}
    if kind == "mass":
        res = R.pencil(R.pos, 4, label="gap")
        rec["pencil"] = res
        wc = int(np.ceil(1.6 / R.h)); Lf = (n - 2 * wc + 1) * R.h
        rec["expected_lowest_omega2"] = 1.0 + 3.0 * (2.0 / R.h * np.sin(np.pi / Lf * R.h / 2)) ** 2
        rec["L_free"] = Lf
        rec["lowest_omega2"] = res["modes"][0]["omega2_S2_over_2K"]
        rec["pass_2_percent"] = bool(abs(rec["lowest_omega2"] - 1.0) < 0.02)
    else:
        # restrict to the internal-rotation directions of the vacuum: the Goldstone sector (rot12, rot13 about (g,1,d,d); rot23 is a symmetry)
        res = R.pencil(R.pos, 6, label="gap")
        rec["pencil"] = res
        wc = int(np.ceil(1.6 / R.h)); Lf = (n - 2 * wc + 1) * R.h
        rec["expected_lowest_omega2"] = 3.0 * (2.0 / R.h * np.sin(np.pi / Lf * R.h / 2)) ** 2
        rec["L_free"] = Lf
        rec["pass_within_fit_error_proxy_10_percent"] = None
        rec["lowest_omega2"] = res["modes"][0]["omega2_S2_over_2K"]
        rec["pass_within_fit_error_proxy_10_percent"] = bool(abs(rec["lowest_omega2"] / rec["expected_lowest_omega2"] - 1.0) < 0.1)
    rec["wall_s"] = round(time.time() - T0, 1)
    return rec


def mode_relax(n, L, src, maxit=6000):
    """a c = 0 hedgehog on a new rung through the R23-1 instrument: the job dict of R23-1 with n, L; the warm start
    (if src) is the stored field interpolated to the new grid (zoom), then re-pinned by the seed on the shell."""
    j = {"delta": DELTA, "w1s": W1S, "seed": "rad", "src": None, "c": 0.0, "n": n, "L": float(L), "part": "R29-1 rung"}
    C23.MAX_ITER[n] = maxit                                             # the stage file resumes, so a higher cap continues a stopped rung
    if src:
        Ms = np.load(os.path.join(C23.OUT_NPZ, src + ".npz"))["M"]
        f = n / Ms.shape[0]
        Mi = np.stack([np.stack([zoom(Ms[..., a, b], f, order=1) for b in range(4)], -1) for a in range(4)], -2)
        cfg = R21.cfg_of(n, L, G, DELTA)
        seed = C23.CORES.seed_core(cfg, "rad", DELTA)
        pin = B3.pin_shell(n, cfg["h"])
        Mi[pin] = seed[pin]
        tag = C23.job_tag(j)
        os.makedirs(C23.OUT_NPZ, exist_ok=True)
        stage = os.path.join(C23.OUT_NPZ, tag + "_stage.npz")
        if not os.path.exists(stage):
            np.savez_compressed(stage, M=Mi, done=0, chunks=json.dumps([{"note": f"warm start: {src} interpolated (order 1) and re-pinned", "E": float("nan")}]))
        j["warm_start"] = src
    row = C23.run_job(j)
    row.pop("chunks", None)
    return row


def mode_collect():
    rows = {}
    for f in sorted(os.listdir(OUT_DIR)):
        if f.startswith("gap_") and f.endswith(".json"):
            rows[f] = json.load(open(os.path.join(OUT_DIR, f)))
    pts = []
    for f, r in rows.items():
        pts.append((r["L"], r["n"], r["K_positive"]["modes"][0]["omega2_S2_over_2K"], r["K_negative_ghost_read"]["modes"][0]["S2"], r["static_stiff_directions"]))
    pts.sort()
    out = {"rungs": [{"L": p[0], "n": p[1], "lowest_Kpos_omega2": p[2], "lowest_Kneg_S2": p[3], "static_stiff": p[4]} for p in pts]}
    if len(pts) >= 3:
        x = np.array([1.0 / p[0] for p in pts]); y = np.array([p[2] for p in pts])
        A = np.vstack([np.ones_like(x), x]).T
        coef, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
        resid = float(np.sqrt(res[0] / max(1, len(x) - 2))) if len(res) else 0.0
        out["fit_omega2_vs_1_over_L"] = {"intercept": float(coef[0]), "slope": float(coef[1]), "residual": resid}
    json.dump(out, open(OUT_JSON, "w"), indent=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("args", nargs="*")
    ap.add_argument("--nev", type=int, default=6); ap.add_argument("--from", dest="src", default=None)
    ap.add_argument("--no-constraint", action="store_true"); ap.add_argument("--warm", default=None); ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--suffix", default=""); ap.add_argument("--schur", action="store_true"); ap.add_argument("--maxit", type=int, default=6000)
    ap.add_argument("--L", type=float, default=None)
    a = ap.parse_args()
    global L_OVERRIDE
    L_OVERRIDE = a.L
    os.makedirs(OUT_DIR, exist_ok=True)
    if a.mode == "census":
        rec = mode_census(a.args[0]); name = f"census_{a.args[0]}.json"
    elif a.mode == "gap":
        rec = mode_gap(a.args[0], a.nev, constrain=not a.no_constraint, warm_dir=a.warm, maxiter=a.maxiter, schur=a.schur); name = f"gap_{a.args[0]}{a.suffix}.json"
    elif a.mode == "ghost":
        tag = a.args[0]
        M = np.load(os.path.join(C23.OUT_NPZ, tag + ".npz"))["M"]; n = M.shape[0]
        R = Rung(n, L_OVERRIDE if L_OVERRIDE else n * 1.5, M); R.build_basis()
        gens = generators(R)
        wn = os.path.join(a.warm, f"{tag}_Kneg_modes.npz") if a.warm else None
        rec = {"field": tag, "K_negative_ghost_read": R.pencil(R.neg, a.nev, Y_fields=list(gens.values()) if not a.no_constraint else [], label="ghost", warm=wn, maxiter=a.maxiter, tol=1e-4)}
        np.savez_compressed(os.path.join(OUT_DIR, f"{tag}_Kneg_modes_arpack.npz"), modes=R.last_mode_fields)
        name = f"ghost_{tag}{a.suffix}.json"
    elif a.mode == "schur-post":
        # the K-null elimination applied AFTER convergence to the saved bare K+ modes of a rung: one-step model vs CG truth
        tag = a.args[0]
        M = np.load(os.path.join(C23.OUT_NPZ, tag + ".npz"))["M"]; n = M.shape[0]
        R = Rung(n, L_OVERRIDE if L_OVERRIDE else n * 1.5, M); R.build_basis()
        Z = np.load(os.path.join(OUT_DIR, f"{tag}_Kpos_modes.npz"))["modes"]
        rows = []
        for i in range(Z.shape[0]):
            D = Z[i] / np.sqrt(np.sum(Z[i] ** 2))
            s2 = float(np.sum(D * R.S2_apply(D))); k = R.K_of(D)
            s2_one = float(np.sum(D * R.S2_apply_schur(D)))
            s2_cg, res, steps, zr = R.S2_schur_cg(D)
            rows.append({"mode": i, "S2_bare": s2, "S2_one_step": s2_one, "S2_cg": s2_cg, "cg_residual_ratio": res, "cg_steps": steps, "null_response_norm_ratio": zr,
                         "K": k, "omega2_bare": s2 / (2 * k), "omega2_one_step": s2_one / (2 * k), "omega2_cg": s2_cg / (2 * k)})
            log(f"{tag} mode {i}: bare {s2:.4g} one-step {s2_one:.4g} cg {s2_cg:.4g} (res {res:.1e}, {steps} steps)")
        rec = {"field": tag, "modes": rows}; name = f"schurpost_{tag}{a.suffix}.json"
    elif a.mode in ("control-mass", "control-gold"):
        n, L = int(a.args[0]), float(a.args[1])
        rec = mode_control("mass" if a.mode == "control-mass" else "gold", n, L); name = f"{a.mode}_n{n}_L{L:g}.json"
    elif a.mode == "relax":
        n, L = int(a.args[0]), float(a.args[1])
        rec = mode_relax(n, L, a.src, a.maxit); name = f"relax_n{n}_L{L:g}.json"
    elif a.mode == "collect":
        rec = mode_collect(); name = None
    else:
        raise SystemExit(f"unknown mode {a.mode}")
    if name:
        json.dump(rec, open(os.path.join(OUT_DIR, name), "w"), indent=1, default=float)
        log(f"wrote {name}")
    print(json.dumps({k: v for k, v in rec.items() if k not in ("K_positive", "K_negative_ghost_read", "pencil", "generators")}, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
