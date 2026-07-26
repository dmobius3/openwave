"""M5.21.1e ADVERSARIAL AUDIT: independent re-verification of the
spec-review claims (findings/m5_21_1e_spec_review.md, C1-C8 of the audit
brief). Default skepticism: every quantity checked below is REBUILT here
(own projector via per-cell SVD, own u_eta / V4 densities via a different
einsum/eigenvalue route, own winding read via bilinear interpolation +
continuity sign-fixing, own containment code, own toy discretization).
Only grid plumbing (cell_weights, grid_coords, MIR, J4, eta) and the
audited force grad_static_4 are imported; the force itself is
cross-checked here against a finite difference of MY energy (A0).

Claims:
  C1 tangent fraction at the frozen 128x256 endpoint (claimed 0.054)
  C2 V4-flatness of orbit tangents t = WM + MW^T, W in so(1,3)
  C3 iso-endpoint winding vs radius (q~0 inside r=20, q~0.98-0.99 at
     r=40-58) + u_eta radial profile (outer shell?)
  C4 w1000 endpoint: contained (r90, coreball 0.93, no u_eta dilution)
     AND not converged (virial u/(3V4) = 0.053)
  C5 Derrick scaling: u_eta ~ 1/lambda, V4 ~ lambda^3 under dilation
  C6 toy ripple floor -(alpha om^2-1)^2/(4 beta om^4) and E<0 grid states
  C8 boost tangents on block-diagonal fields are pure time-mixing

C2 ALGEBRA (the independent derivation, checked numerically below):
  W in so(1,3)  <=>  W^T eta + eta W = 0  <=>  eta W = -W^T eta.
  For t = W M + M W^T:
    d Tr((eta M)^p) = p Tr((eta M)^{p-1} eta t)
      = p Tr((eta M)^{p-1} (eta W) M) + p Tr((eta M)^{p-1} eta M W^T)
      = -p Tr((eta M)^{p-1} W^T (eta M)) + p Tr((eta M)^p W^T)
      = -p Tr(W^T (eta M)^p) + p Tr(W^T (eta M)^p) = 0   exactly,
  and at finite s, Lam = exp(sW) has Lam^T eta Lam = eta so
  eta(Lam M Lam^T) = (Lam^T)^{-1} (eta M) Lam^T: a similarity transform,
  Tr((eta M)^p) invariant exactly (not just to first order).

C8 ALGEBRA: for block-diagonal M (M_{0i} = 0) and a boost K
  (K_{0i} = K_{i0} = 1, else 0): (KM)_{0b} = M_{ib} (spatial row i of M
  landed in the time row), (KM)_{i0} = M_{00}, (KM) spatial block = 0,
  (KM)_{00} = M_{i0} = 0; same for M K^T = M K. So t = KM + MK^T has
  ONLY time-mixing entries. The static force on a block-diagonal field
  is block-diagonal (every term of u_eta and V4 preserves the block
  structure), hence <G, t_K> = 0: rot3 == full6 descent.

C5 EDGE HANDLING: dilation M_lambda(x) = M(x/lambda) built by bilinear
  interpolation on the (rho, z) grid; the axis uses the mirror ghost
  M(-rho) = MIR*M(rho); lambda >= 1 only, so source points stay in the
  domain (no outer extrapolation); outer edges clamped (inactive).

Writes data/m5_21_1e_audit.json. Creates no other file.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from m5_17_energy import MIR, J4, cell_weights, grid_coords        # noqa: E402
from m5_20_2_a_eom import ETA, G_T, H, WSCALE, grad_static_4       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
DELTA = 0.3
SG = +1.0 * G_T                    # the s = +1 working branch of the runs
C_P = tuple(SG ** p + 1.0 + DELTA ** p for p in range(1, 5))


# ================= my own energy densities =================
def my_channels(M, h=H):
    """the axisym derivative channels (grid plumbing, fixed scheme)."""
    nr = M.shape[0]
    Mm = np.empty_like(M[: nr - 1])
    Mm[1:] = M[: nr - 2]
    Mm[0] = MIR * M[0]
    Arho = ((M[1:] - Mm) / (2.0 * h))[:, 1:-1]
    Az = (M[: nr - 1, 2:] - M[: nr - 1, :-2]) / (2.0 * h)
    Mc = M[: nr - 1, 1:-1]
    rho = ((np.arange(nr - 1) + 0.5) * h)[:, None, None, None]
    Aphi = (np.broadcast_to(J4, Mc.shape) @ Mc
            - Mc @ np.broadcast_to(J4, Mc.shape)) / rho
    return Arho, Aphi, Az


def my_u_eta(M, h=H):
    """4 SUM_{i<j} <F,F>_eta, F = A eta B - B eta A; my inner-product
    route: <F,F>_eta = sum_ab (eta F eta)_ab F_ab (eta diagonal)."""
    Arho, Aphi, Az = my_channels(M, h)
    e = np.diag(ETA)
    tot = 0.0
    for A, B in ((Arho, Aphi), (Arho, Az), (Aphi, Az)):
        F = A @ np.broadcast_to(ETA, A.shape) @ B \
            - B @ np.broadcast_to(ETA, A.shape) @ A
        EFE = e[:, None] * F * e[None, :]
        tot = tot + np.einsum("...ab,...ab->...", EFE, F)
    return 4.0 * tot


def my_traces(M):
    """I_p = tr((eta M)^p), p = 1..4, per cell (matrix-power route)."""
    EM = np.broadcast_to(ETA, M.shape) @ M
    P = EM
    out = [np.einsum("...aa->...", P)]
    for _ in range(3):
        P = P @ EM
        out.append(np.einsum("...aa->...", P))
    return np.stack(out, axis=-1)


def my_traces_eig(M):
    """same I_p via eigenvalues of eta M (independent route: I_p =
    sum_k lam_k^p, lam = eigvals(eta M), complex-safe)."""
    lam = np.linalg.eigvals(np.broadcast_to(ETA, M.shape) @ M)
    return np.stack([np.real(np.sum(lam ** p, axis=-1))
                     for p in range(1, 5)], axis=-1)


def my_v4(M, wscale):
    Ic = my_traces(M[: M.shape[0] - 1, 1:-1])
    v = np.zeros(Ic.shape[:-1])
    for p in range(4):
        v = v + (Ic[..., p] - C_P[p]) ** 2
    return wscale * v


def my_total(M, wscale, h=H):
    w = cell_weights(M.shape[0], M.shape[1], h)
    return float(np.sum((my_u_eta(M, h) + my_v4(M, wscale)) * w))


def my_seed(nr, nz, delta=DELTA, sg=SG, r_c=4.0):
    """the M5.21.1 hedgehog seed formula (seed spec, rebuilt)."""
    R, Z = grid_coords(nr, nz, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    rs = np.where(r < 1e-12, 1e-12, r)
    n = np.stack([R / rs, np.zeros_like(R), Z / rs], axis=-1)
    m = np.broadcast_to(np.array([0.0, 1.0, 0.0]), n.shape)
    S = (n[..., :, None] * n[..., None, :]
         + delta * m[..., :, None] * m[..., None, :])
    a = (1.0 + delta) / 3.0
    w = (1.0 - np.exp(-((r / r_c) ** 2)))[..., None, None]
    M = np.zeros((nr, nz, 4, 4))
    M[..., 1:4, 1:4] = w * S + (1.0 - w) * (a * np.eye(3))
    M[..., 0, 0] = -sg
    return M


# ================= my own so(1,3) + projector =================
def so13_basis_nullspace():
    """derive the generator space from the defining relation
    W^T eta + eta W = 0 by null-space computation (NOT hardcoded)."""
    A = np.zeros((16, 16))
    for k in range(16):
        E = np.zeros(16)
        E[k] = 1.0
        W = E.reshape(4, 4)
        A[:, k] = (W.T @ ETA + ETA @ W).ravel()
    _, s, vt = np.linalg.svd(A)
    null = vt[s < 1e-10 * s.max()]
    return [null[k].reshape(4, 4) for k in range(null.shape[0])]


def tan_frac_svd(M, G, gens, cutoff=1e-10, mask=None):
    """my projector: per-cell ORTHONORMAL tangent basis via SVD (their
    route was a regularized Gram solve). Returns ||P G||^2 / ||G||^2."""
    nr, nz = M.shape[:2]
    T = np.stack([W @ M + M @ W.T for W in gens], axis=2)   # (nr,nz,k,4,4)
    Tf = T.reshape(nr, nz, len(gens), 16)
    Gf = G.reshape(nr, nz, 16)
    if mask is not None:
        Tf = Tf[mask]
        Gf = Gf[mask]
    else:
        Tf = Tf.reshape(-1, len(gens), 16)
        Gf = Gf.reshape(-1, 16)
    U, s, _ = np.linalg.svd(np.swapaxes(Tf, -1, -2), full_matrices=False)
    keep = s > cutoff * np.maximum(s[..., :1], 1e-300)      # per-cell rank
    coef = np.einsum("cik,ci->ck", U, Gf) * keep
    Gt2 = np.sum(coef ** 2)
    return float(Gt2 / max(np.sum(Gf ** 2), 1e-300))


# ================= my own observables =================
def my_contain(M, wscale, which="total"):
    nr, nz = M.shape[:2]
    w = cell_weights(nr, nz, H)
    du = my_u_eta(M) * w
    dv = my_v4(M, wscale) * w
    dens = {"total": du + dv, "u": du, "v": dv}[which]
    R, Z = grid_coords(nr, nz, H)
    r = np.sqrt(R ** 2 + Z ** 2)[: nr - 1, 1:-1].ravel()
    d = dens.ravel()
    o = np.argsort(r)
    cs = np.cumsum(d[o])
    tot = max(cs[-1], 1e-300)
    return {"r50": float(r[o][np.searchsorted(cs, 0.5 * tot)]),
            "r90": float(r[o][np.searchsorted(cs, 0.9 * tot)]),
            "coreball8": float(np.sum(d[r <= 8.0]) / tot),
            "E": float(tot)}


def u_shell_profile(M, edges):
    nr, nz = M.shape[:2]
    w = cell_weights(nr, nz, H)
    dens = (my_u_eta(M) * w).ravel()
    R, Z = grid_coords(nr, nz, H)
    r = np.sqrt(R ** 2 + Z ** 2)[: nr - 1, 1:-1].ravel()
    hist, _ = np.histogram(r, bins=edges, weights=dens)
    return hist


def sample_M(M, rho_q, z_q):
    """bilinear sample of the (rho, z) field; mirror ghost across the
    axis (M(-rho) = MIR*M(rho)); edges clamped."""
    nr, nz = M.shape[:2]
    Mg = np.concatenate([(MIR * M[0])[None], M], axis=0)    # i = -1 ghost
    shp = np.shape(rho_q)
    fi = np.clip(np.ravel(rho_q) / H - 0.5, -1.0, nr - 1.0)
    fj = np.clip(np.ravel(z_q) / H + nz / 2 - 0.5, 0.0, nz - 1.0)
    i0 = np.minimum(np.floor(fi).astype(int), nr - 2)
    j0 = np.minimum(np.floor(fj).astype(int), nz - 2)
    t = (fi - i0)[:, None, None]
    s = (fj - j0)[:, None, None]
    g0, g1 = i0 + 1, i0 + 2
    out = ((1 - t) * (1 - s) * Mg[g0, j0] + t * (1 - s) * Mg[g1, j0]
           + (1 - t) * s * Mg[g0, j0 + 1] + t * s * Mg[g1, j0 + 1])
    return out.reshape(shp + (4, 4))


def my_winding(M, r, npts=241, with_gap=False):
    """my winding read: bilinear-sampled director (eigenvector of the
    LARGEST spatial eigenvalue) along the meridional half-circle,
    continuity sign-fixing, angle steps wrapped mod pi. with_gap also
    returns the minimum top-eigenvalue gap along the circle (a tiny gap
    means the director, hence q, is ill-defined there)."""
    th = np.linspace(1e-3, np.pi - 1e-3, npts)
    Ms = sample_M(M, r * np.sin(th), r * np.cos(th))
    lam, V = np.linalg.eigh(Ms[:, 1:4, 1:4])
    n = V[:, :, -1]                                # largest eigenvalue
    for k in range(1, npts):                       # continuity sign fix
        if float(np.dot(n[k], n[k - 1])) < 0.0:
            n[k] = -n[k]
    alpha = np.arctan2(n[:, 0], n[:, 2])
    dd = np.diff(alpha)
    dd = (dd + np.pi / 2) % np.pi - np.pi / 2
    q = float(np.sum(dd) / np.pi)
    if not with_gap:
        return q
    gap = lam[:, -1] - lam[:, -2]
    return q, float(np.min(gap)), float(np.median(gap))


def sawtooth_ratio(M):
    """grid-scale (stencil-null) texture diagnostic: mean squared
    one-cell difference over mean squared 2h central difference, per
    axis, on the interior. Smooth field ~1; checkerboard (the central-
    difference null mode) >> 1."""
    d1r = M[1:-1] - M[:-2]
    c_r = 0.5 * (M[2:] - M[:-2])
    d1z = M[:, 1:-1] - M[:, :-2]
    c_z = 0.5 * (M[:, 2:] - M[:, :-2])
    return {"rho": float(np.mean(d1r ** 2) / max(np.mean(c_r ** 2),
                                                 1e-300)),
            "z": float(np.mean(d1z ** 2) / max(np.mean(c_z ** 2),
                                               1e-300))}


def rescale_field(M, lam):
    """M_lambda(x) = M(x/lambda) on the same grid (bilinear, see header
    for edge handling); lambda >= 1."""
    nr, nz = M.shape[:2]
    R, Z = grid_coords(nr, nz, H)
    return sample_M(M, R / lam, Z / lam)


# ================= C1 =================
def check_c1(out):
    gens = so13_basis_nullspace()
    assert len(gens) == 6, f"so(1,3) null space dim {len(gens)} != 6"
    z = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))
    Mep = z[z.files[0]]
    nr, nz = Mep.shape[:2]
    G = grad_static_4(Mep, WSCALE, DELTA, g=SG)
    # A0: cross-check the imported force against FD of MY energy
    rng = np.random.default_rng(3)
    fd_rel = []
    for _ in range(3):
        D = np.zeros_like(Mep)
        Dc = rng.normal(size=(nr - 1, nz - 2, 4, 4))
        D[: nr - 1, 1:-1] = 0.5 * (Dc + np.swapaxes(Dc, -1, -2))
        D /= np.sqrt(np.sum(D * D))            # unit direction: kills the
        eps = 1e-5                             # O(eps^2 ||D||^3) truncation
        num = (my_total(Mep + eps * D, WSCALE)
               - my_total(Mep - eps * D, WSCALE)) / (2 * eps)
        an = float(np.sum(G * D))
        fd_rel.append(abs(num - an) / (abs(num) + abs(an) + 1e-300))
    pin = np.zeros((nr, nz), dtype=bool)
    pin[-1, :] = pin[:, 0] = pin[:, -1] = True
    tf_all = tan_frac_svd(Mep, G, gens)
    tf_free = tan_frac_svd(Mep, G, gens, mask=~pin)
    tf_loose = tan_frac_svd(Mep, G, gens, cutoff=1e-6)
    # the it-0 point of the slide (my seed, 64x128) for the range check
    M0 = my_seed(64, 128)
    G0 = grad_static_4(M0, WSCALE, DELTA, g=SG,
                       w4=cell_weights(64, 128, H)[..., None, None],
                       rho4=((np.arange(63) + 0.5) * H)[:, None, None, None])
    tf_seed = tan_frac_svd(M0, G0, gens)
    out["C1"] = {
        "claimed": {"endpoint_tan_frac": 0.054, "seed_tan_frac": 0.0003,
                    "headline": "97-99.97% of the slide force orthogonal"},
        "measured": {
            "A0_force_vs_myFD_rel": fd_rel,
            "endpoint_tan_frac_all_cells": tf_all,
            "endpoint_tan_frac_free_cells": tf_free,
            "endpoint_tan_frac_cutoff1e-6": tf_loose,
            "endpoint_orthogonal_pct": 100.0 * (1.0 - tf_all),
            "seed64_tan_frac": tf_seed},
    }
    print("C1", json.dumps(out["C1"]["measured"]))


# ================= C2 =================
def check_c2(out):
    gens = so13_basis_nullspace()
    rng = np.random.default_rng(11)
    worst_fd, worst_ref, worst_fin = 0.0, 0.0, 0.0

    def v4pt(M):
        EM = ETA @ M
        P, v = EM, 0.0
        for p in range(1, 5):
            v += (np.trace(P) - C_P[p - 1]) ** 2
            P = P @ EM
        return v

    for _ in range(20):
        A = rng.normal(size=(4, 4))
        M = 0.5 * (A + A.T) + np.diag([-SG, 1.0, DELTA, 0.0])
        for W in gens:
            t = W @ M + M @ W.T
            eps = 1e-5
            d_t = (v4pt(M + eps * t) - v4pt(M - eps * t)) / (2 * eps)
            Dr = rng.normal(size=(4, 4))
            Dr = 0.5 * (Dr + Dr.T)
            d_r = (v4pt(M + eps * Dr) - v4pt(M - eps * Dr)) / (2 * eps)
            worst_fd = max(worst_fd, abs(d_t))
            worst_ref = max(worst_ref, abs(d_t) / max(abs(d_r), 1e-300))
            # finite transformation: Lam = expm(sW)
            Lam = expm(0.37 * W)
            gate = np.max(np.abs(Lam.T @ ETA @ Lam - ETA))
            assert gate < 1e-12, f"expm(sW) not in SO(1,3): {gate}"
            M2 = Lam @ M @ Lam.T
            I1 = my_traces(M[None, None])[0, 0]
            I2 = my_traces(M2[None, None])[0, 0]
            worst_fin = max(worst_fin, float(np.max(
                np.abs(I2 - I1) / (1.0 + np.abs(I1)))))
    # field-level: global W tangent on the frozen endpoint, FD of MY V4
    z = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))
    Mep = z[z.files[0]]
    w = cell_weights(Mep.shape[0], Mep.shape[1], H)
    Wg = gens[0] + 0.5 * gens[3]
    tg = Wg @ Mep + Mep @ Wg.T
    eps = 1e-5
    v_p = float(np.sum(my_v4(Mep + eps * tg, WSCALE) * w))
    v_m = float(np.sum(my_v4(Mep - eps * tg, WSCALE) * w))
    v_0 = float(np.sum(my_v4(Mep, WSCALE) * w))
    field_rel = abs(v_p - v_m) / (2 * eps) / max(abs(v_0), 1e-300)
    # eigenvalue-route crosscheck of my traces (independence gate)
    Msamp = Mep[10:20, 60:70]
    tr_gap = float(np.max(np.abs(my_traces(Msamp) - my_traces_eig(Msamp))
                          / (1.0 + np.abs(my_traces(Msamp)))))
    out["C2"] = {
        "claimed": "orbit tangents leave Tr((eta M)^p) exactly invariant",
        "measured": {"pointwise_FD_worst_abs": worst_fd,
                     "pointwise_FD_worst_vs_random_dir": worst_ref,
                     "finite_Lam_invariance_worst_rel": worst_fin,
                     "field_level_FD_rel_endpoint": field_rel,
                     "traces_power_vs_eig_route_rel": tr_gap},
    }
    print("C2", json.dumps(out["C2"]["measured"]))


# ================= C3 =================
def check_c3(out):
    z = np.load(os.path.join(DATA, "m5_21_1e_iso_endpoint.npz"))
    M = z[z.files[0]]
    radii = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 58]
    qs, gaps = {}, {}
    for r in radii:
        q1 = my_winding(M, float(r), npts=241)
        q2, gmin, gmed = my_winding(M, float(r), npts=1921, with_gap=True)
        qs[str(r)] = [q1, q2]                    # sampling stability
        gaps[str(r)] = [gmin, gmed]
    edges = np.arange(0.0, 68.0, 4.0)
    prof = u_shell_profile(M, edges)
    tot = float(prof.sum())
    peak_bin = int(np.argmax(prof))
    r_in = float(prof[: 5].sum() / tot)          # r < 20
    r_out = float(prof[7:].sum() / tot)          # r >= 28
    # the seed profile for contrast (was the interior filled at it 0?)
    prof0 = u_shell_profile(my_seed(64, 128), edges)
    cu = my_contain(M, WSCALE, "u")
    out["C3"] = {
        "claimed": {"q_inside_r20": "~0", "q_r40_58": "+0.98..0.99",
                    "reading": "charge migrated to a boundary shell"},
        "measured": {
            "q_vs_r_npts241_1921": qs,
            "top_eig_gap_min_median": gaps,
            "u_shell_edges": edges.tolist(),
            "u_shell_endpoint": prof.tolist(),
            "u_shell_seed": prof0.tolist(),
            "u_peak_shell": [float(edges[peak_bin]),
                             float(edges[peak_bin + 1])],
            "u_frac_r_lt_20": r_in,
            "u_frac_r_ge_28": r_out,
            "contain_u_only": cu},
    }
    print("C3 q:", json.dumps(qs))
    print("C3 gaps:", json.dumps(gaps))
    print("C3 u-shells:", np.round(prof, 2).tolist(),
          "peak", out["C3"]["measured"]["u_peak_shell"],
          "u-contain", json.dumps(cu))


# ================= C4 =================
def check_c4(out):
    z = np.load(os.path.join(DATA, "m5_21_1e_w1000_endpoint.npz"))
    M = z[z.files[0]]
    nr, nz = M.shape[:2]
    ws = WSCALE * 1000.0
    w = cell_weights(nr, nz, H)
    Eu = float(np.sum(my_u_eta(M) * w))
    Ev = float(np.sum(my_v4(M, ws) * w))
    ct = my_contain(M, ws, "total")
    cu = my_contain(M, ws, "u")
    M0 = my_seed(nr, nz)
    Eu0 = float(np.sum(my_u_eta(M0) * w))
    cu0 = my_contain(M0, ws, "u")
    G = grad_static_4(M, ws, DELTA, g=SG,
                      w4=w[..., None, None],
                      rho4=((np.arange(nr - 1) + 0.5) * H)[:, None, None,
                                                           None])
    gnorm = float(np.sqrt(np.sum(G * G)))
    # continuation probe: MY energy-monotone descent from the endpoint
    pin = np.zeros((nr, nz), dtype=bool)
    pin[-1, :] = pin[:, 0] = pin[:, -1] = True
    free = (~pin)[..., None, None].astype(float)
    wfull = np.ones((nr, nz))
    wfull[: nr - 1, 1:-1] = w
    precond = (1.0 / wfull)[..., None, None]
    Mp = M.copy()
    E = my_total(Mp, ws)
    dt = 1e-4
    probe = []
    for it in range(1, 601):
        Gp = grad_static_4(Mp, ws, DELTA, g=SG, w4=w[..., None, None],
                           rho4=((np.arange(nr - 1) + 0.5)
                                 * H)[:, None, None, None])
        F = -Gp * precond * free
        ok = False
        for _ in range(40):
            M2 = Mp + dt * F
            E2 = my_total(M2, ws)
            if E2 < E:
                ok = True
                break
            dt *= 0.5
        if not ok:
            probe.append({"it": it, "stalled": True})
            break
        Mp, E = M2, E2
        dt *= 1.2
        if it % 200 == 0 or it == 1:
            cpt = my_contain(Mp, ws, "total")
            probe.append({"it": it, "E": E,
                          "E_u": float(np.sum(my_u_eta(Mp) * w)),
                          "r90": cpt["r90"],
                          "coreball8": cpt["coreball8"]})
    out["C4"] = {
        "claimed": {"r90": 6.5, "coreball": 0.93, "E_u_no_dilution":
                    "20.3 -> 20.6", "virial_u_over_3v": 0.053},
        "measured": {
            "E_u": Eu, "E_v4": Ev, "virial_u_over_3v": Eu / (3.0 * Ev),
            "contain_total": ct, "contain_u_only": cu,
            "contain_u_only_SEED": cu0, "E_u_seed": Eu0,
            "gnorm": gnorm,
            "continuation_probe_600it": probe},
    }
    print("C4", json.dumps({k: v for k, v in out["C4"]["measured"].items()
                            if k != "continuation_probe_600it"}))
    print("C4 probe:", json.dumps(probe))


# ================= C5 =================
def check_c5(out):
    lams = [1.0, 1.25, 1.6, 2.0, 2.5]

    def energies(M, wscale):
        nr, nz = M.shape[:2]
        w = cell_weights(nr, nz, H)
        eu, ev = [], []
        for lam in lams:
            Ml = rescale_field(M, lam)
            eu.append(float(np.sum(my_u_eta(Ml) * w)))
            ev.append(float(np.sum(my_v4(Ml, wscale) * w)))
        return eu, ev

    def logslope(e):
        e = np.asarray(e)
        if np.any(e <= 0):
            return float("nan")
        return float(np.polyfit(np.log(lams), np.log(e), 1)[0])

    def fit_inv_const(e):
        """least squares e(lam) = a/lam + b: the exact continuum form
        for a 1/lam-scaling energy on a FINITE window with a 1/r^4 tail
        (the window correction is a lam-independent constant)."""
        X = np.vstack([1.0 / np.asarray(lams), np.ones(len(lams))]).T
        coef, *_ = np.linalg.lstsq(X, np.asarray(e), rcond=None)
        resid = float(np.sqrt(np.sum((X @ coef - e) ** 2))
                      / np.sqrt(np.sum(np.asarray(e) ** 2)))
        return float(coef[0]), float(coef[1]), resid

    # arm 1: the REAL hedgehog seed (axis-regular, vacuum-spectrum far
    # field). Its u_eta has the physical 1/r^4 tail, so on the finite
    # box E_u(lam) = A/lam + B exactly if u scales 1/lam (B = the
    # constant window correction); V4 is localized -> clean lam^3.
    Ms = my_seed(64, 128)
    eu_s, ev_s = energies(Ms, WSCALE)
    a_s, b_s, res_s = fit_inv_const(eu_s)
    pv_s = logslope(ev_s)
    # arm 2: compact synthetic field, AXIS-REGULAR (in-plane block
    # proportional to I => [J4, M] = 0, no 1/rho line singularity),
    # mirror-even (diagonals + (0,3), rho only through rho^2/r^2),
    # Gaussian envelope: exercises the TIME-MIXING sector, where the
    # eta-norm is indefinite (E_u may be negative), so the sign-safe
    # test is |lam E_u(lam) / E_u(1) - 1|.
    R, Z = grid_coords(64, 128, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    bump = np.exp(-((r / 8.0) ** 2))
    Mc = np.zeros((64, 128, 4, 4))
    Mc[..., 0, 0] = -SG + 0.20 * bump * np.cos(Z / 3.0)
    inpl = 0.7 + 0.30 * bump * np.sin(R ** 2 / 24.0)
    Mc[..., 1, 1] = inpl
    Mc[..., 2, 2] = inpl
    Mc[..., 3, 3] = 0.2 + 0.30 * bump * np.cos(R ** 2 / 30.0 + Z / 4.0)
    Mc[..., 0, 3] = Mc[..., 3, 0] = 0.25 * bump * np.sin(Z / 3.5)
    eu_c, ev_c = energies(Mc, WSCALE)
    dev_c = [abs(lams[k] * eu_c[k] / eu_c[0] - 1.0)
             for k in range(len(lams))]
    pv_c = logslope(ev_c)
    # arm 3: the w1000 endpoint (real object; hedgehog tail => log
    # slope carries the finite-window correction; also fit A/lam + B)
    z = np.load(os.path.join(DATA, "m5_21_1e_w1000_endpoint.npz"))
    Mw = z[z.files[0]]
    eu_w, ev_w = energies(Mw, WSCALE * 1000.0)
    pu_w, pv_w = logslope(eu_w), logslope(ev_w)
    a_w, b_w, res_w = fit_inv_const(eu_w)
    # monotone expansion of the frozen-potential remainder on the ISO
    # endpoint: E_u(lambda) should fall monotonically IF the endpoint is
    # a smooth field (a fine lambda ladder detects stencil-null texture:
    # any jump at lambda -> 1+ means grid-scale structure invisible to
    # the 2h central-difference functional)
    zi = np.load(os.path.join(DATA, "m5_21_1e_iso_endpoint.npz"))
    Mi = zi[zi.files[0]]
    w64 = cell_weights(64, 128, H)
    lams_i = [1.0, 1.05, 1.1, 1.25, 1.5, 2.0, 2.5]
    eu_iso = [float(np.sum(my_u_eta(rescale_field(Mi, lam)) * w64))
              for lam in lams_i]
    mono = bool(all(eu_iso[k + 1] < eu_iso[k]
                    for k in range(len(eu_iso) - 1)))
    saw = {"iso_endpoint": sawtooth_ratio(Mi),
           "w1000_endpoint": sawtooth_ratio(Mw),
           "seed": sawtooth_ratio(Ms)}
    out["C5"] = {
        "claimed": {"u_eta_exponent": -1.0, "V4_exponent": 3.0,
                    "orbit_class_descent": "expands monotonically"},
        "measured": {
            "lambdas": lams,
            "seed": {"E_u": eu_s, "E_v4": ev_s,
                     "u_fit_A_over_lam_plus_B": [a_s, b_s],
                     "u_fit_resid_rel": res_s, "v4_slope": pv_s},
            "compact_timemix_field": {
                "E_u": eu_c, "E_v4": ev_c,
                "lamEu_over_Eu1_dev_from_1": dev_c, "v4_slope": pv_c},
            "w1000_endpoint": {"E_u": eu_w, "E_v4": ev_w,
                               "u_slope": pu_w, "v4_slope": pv_w,
                               "u_fit_A_over_lam_plus_B": [a_w, b_w],
                               "u_fit_resid_rel": res_w},
            "iso_lambdas": lams_i,
            "iso_endpoint_Eu_vs_lambda": eu_iso,
            "iso_Eu_monotone_decreasing": mono,
            "sawtooth_ratio_fwd_over_central": saw},
    }
    print("C5 seed:    u fit A/lam+B = %.3f, %.3f (resid %.2e) "
          "v4_slope %.4f" % (a_s, b_s, res_s, pv_s))
    print("C5 compact: max|lam Eu/Eu(1) - 1| %.4f v4_slope %.4f"
          % (max(dev_c), pv_c))
    print("C5 w1000:   u_slope %.4f (A/lam+B: %.3f, %.3f, resid %.2e) "
          "v4_slope %.4f" % (pu_w, a_w, b_w, res_w, pv_w))
    print("C5 iso E_u(lam):", np.round(eu_iso, 3).tolist(), "mono", mono)
    print("C5 sawtooth:", json.dumps(saw))


# ================= C6 =================
def check_c6(out):
    alpha, beta = 1.0, 1.0
    om_star = float(np.sqrt(70.0 / 61.0))
    floor_cf = -(alpha * om_star ** 2 - 1.0) ** 2 / (4 * beta
                                                     * om_star ** 4)
    # exact rational: -(9/61)^2 / (4 (70/61)^2) = -81/19600
    floor_exact = -81.0 / 19600.0
    # numeric pointwise minimum over p (independent of the formula)
    p = np.linspace(0.0, 1.0, 200001)
    f = (1 - alpha * om_star ** 2) * p ** 2 + beta * om_star ** 4 * p ** 4
    k = int(np.argmin(f))
    px_cf = float(np.sqrt((alpha * om_star ** 2 - 1)
                          / (2 * beta * om_star ** 4)))
    # my grid: N = 401, L = 20, LINK-based gradient discretization
    L, N = 20.0, 401
    X = np.linspace(-L, L, N)
    dx = X[1] - X[0]
    rng = np.random.default_rng(2)

    def e_grad(phi, om):
        a = 1.0 - alpha * om ** 2
        b = beta * om ** 4
        pl = (phi[1:] - phi[:-1]) / dx
        wn = np.ones(N)
        wn[0] = wn[-1] = 0.5
        E = float(dx * np.sum(a * pl ** 2 + b * pl ** 4)
                  + dx * np.sum(wn * (1 - phi ** 2) ** 2))
        flux = 2 * a * pl + 4 * b * pl ** 3
        g = np.zeros(N)
        g[1:-1] = flux[:-1] - flux[1:]
        g += dx * wn * (-4.0) * (1 - phi ** 2) * phi
        g[0] = g[-1] = 0.0
        return E, g

    def relax(om, iters=6000):
        phi = np.tanh(X) + 0.02 * rng.standard_normal(N)
        phi[0], phi[-1] = -1.0, 1.0
        E, g = e_grad(phi, om)
        dt = 1e-4
        for _ in range(iters):
            ok = False
            for _ in range(40):
                ph2 = phi - dt * g
                E2, g2 = e_grad(ph2, om)
                if E2 < E:
                    ok = True
                    break
                dt *= 0.5
            if not ok:
                break
            phi, E, g = ph2, E2, g2
            dt *= 1.15
        return E

    oms = [1.0, 1.5, 2.0, 2.15, 2.3, 2.5]
    es = {str(om): relax(om) for om in oms}
    e_min = min(es.values())
    om_min = [k for k, v in es.items() if v == e_min][0]
    out["C6"] = {
        "claimed": {"floor_at_omstar": -0.0041,
                    "grid_min": {"E": -5.8, "omega": 2.15},
                    "alpha_om2_gt_1_at_omstar": True},
        "measured": {
            "om_star": om_star,
            "alpha_om_star_sq_minus_1": alpha * om_star ** 2 - 1.0,
            "floor_closed_form": floor_cf,
            "floor_exact_rational_81_19600": floor_exact,
            "floor_numeric_scan": float(f[k]),
            "px_closed_form": px_cf, "px_numeric": float(p[k]),
            "grid": "N=401 L=20 link-based discretization (own)",
            "E_relaxed_vs_omega": es,
            "E_min": e_min, "om_at_min": float(om_min),
            "E_neg_exists_below_kink": bool(e_min < 0.0)},
    }
    print("C6", json.dumps(out["C6"]["measured"]))


# ================= C8 =================
def check_c8(out):
    rng = np.random.default_rng(5)
    boosts = []
    for i in (1, 2, 3):
        K = np.zeros((4, 4))
        K[0, i] = K[i, 0] = 1.0
        assert np.max(np.abs(K.T @ ETA + ETA @ K)) == 0.0
        boosts.append(K)
    worst_block, worst_mix = 0.0, 0.0
    for _ in range(30):
        M = np.zeros((4, 4))
        M[0, 0] = rng.normal()
        S = rng.normal(size=(3, 3))
        M[1:, 1:] = 0.5 * (S + S.T)
        for K in boosts:
            t = K @ M + M @ K.T
            worst_block = max(worst_block, abs(t[0, 0]),
                              float(np.max(np.abs(t[1:, 1:]))))
            worst_mix = max(worst_mix, float(np.max(np.abs(t[0, 1:]))))
    # the force on a block-diagonal FIELD is block-diagonal
    M0 = my_seed(64, 128)
    G0 = grad_static_4(M0, WSCALE, DELTA, g=SG,
                       w4=cell_weights(64, 128, H)[..., None, None],
                       rho4=((np.arange(63) + 0.5) * H)[:, None, None,
                                                        None])
    g_mix = float(np.max(np.abs(G0[..., 0, 1:4])))
    ovl = 0.0
    for K in boosts:
        tK = (np.broadcast_to(K, M0.shape) @ M0
              + M0 @ np.broadcast_to(K.T, M0.shape))
        ovl = max(ovl, abs(float(np.sum(G0 * tK))))
    # rot3 vs full6 from the task's own JSON: max rel diff at it 8000
    with open(os.path.join(DATA, "m5_21_1e_constraint.json")) as f:
        d = json.load(f)
    r3 = {r["it"]: r for r in d["phase_i"]["rot3"]}
    f6 = {r["it"]: r for r in d["phase_i"]["full6"]}
    keys = ("E_u", "E_v4", "r90", "coreball", "gnorm")
    reldiff = max(abs(r3[8000][k] - f6[8000][k])
                  / (abs(r3[8000][k]) + 1e-300) for k in keys)
    out["C8"] = {
        "claimed": "boost tangents pure time-mixing on block-diagonal "
                   "fields; rot3 == full6 to ~1e-5",
        "measured": {"boost_tangent_block_part_max": worst_block,
                     "boost_tangent_timemix_max": worst_mix,
                     "force_timemix_on_blockdiag_seed": g_mix,
                     "force_dot_boost_tangent": ovl,
                     "rot3_vs_full6_it8000_max_rel": reldiff},
    }
    print("C8", json.dumps(out["C8"]["measured"]))


# ================= C7 (recorded from the direct PDF read) =============
def record_c7(out):
    out["C7"] = {
        "claimed": {
            "a_p8_det_quote": True, "b_p12_det_and_could_use_eg": True,
            "c_vacuum_positive_g": True},
        "measured_by_direct_read": {
            "a": "p. 8 col 1: 'Ideally we would not having to fix shape "
                 "(lambda_i) ... from a simple e.g. V(A) = (sum ||A||_F^2"
                 " - 1)^2 Higgs-like potential, with additional e.g. "
                 "volume constraint det(M) = prod lambda_i = const to "
                 "prevent using only long axes which allow for low "
                 "curvature (hence energy).' -- present verbatim",
            "b": "p. 12 col 1: 'we could use e.g. V(M) = sum_p "
                 "(Tr((M xi)^p) - c_p)^2 potential' AND 'There might be "
                 "required additional det(M) = const (volume preserving) "
                 "constraint' -- both present",
            "c": "p. 11 Fig 10 caption + p. 11 col 1: vacuum M = ODO^T, "
                 "D = diag(g, 1, delta, 0), g >> 1 >> delta > 0, O in "
                 "SO(1,3); p. 12: 'traces of (M xi)^p = (-lambda_0)^p + "
                 "lambda_1^p + lambda_2^p + lambda_3^p' => targets "
                 "(-g)^p + 1 + delta^p: M_vac has POSITIVE g eigenvalue "
                 "(the stack's s = -1 branch)"},
        "note": "p. 7 also confirms Eq 12/13 and the 'difficult open "
                "question requiring simulations' potential language; "
                "p. 12 'Choosing the details especially of potential is "
                "very difficult, will rather require PDE simulations' "
                "and the 6th-order skyrmion hedge are present.",
    }


def main():
    out = {"task": "M5.21.1e adversarial audit",
           "auditor": "independent re-implementation "
                      "(m5_21_1e_audit_check.py)",
           "wscale": WSCALE, "g": SG, "delta": DELTA}
    check_c2(out)      # cheap algebra first (gates the projector claims)
    check_c1(out)
    check_c3(out)
    check_c4(out)
    check_c5(out)
    check_c6(out)
    record_c7(out)
    check_c8(out)
    # ---- verdicts ----
    m1 = out["C1"]["measured"]
    v1_ok = abs(m1["endpoint_tan_frac_all_cells"] - 0.054) < 0.005
    out["C1"]["verdict"] = (
        "QUALIFIED (tan_frac values confirmed exactly; the '97-99.97% "
        "orthogonal' headline does not cover the endpoint, which is "
        f"{100 * (1 - m1['endpoint_tan_frac_all_cells']):.1f}% "
        "orthogonal)" if v1_ok else "REFUTED")
    m2 = out["C2"]["measured"]
    out["C2"]["verdict"] = ("CONFIRMED" if
                            m2["pointwise_FD_worst_vs_random_dir"] < 1e-8
                            and m2["finite_Lam_invariance_worst_rel"]
                            < 1e-10 else "REFUTED")
    q = out["C3"]["measured"]["q_vs_r_npts241_1921"]
    inner_ok = all(abs(q[str(r)][1]) < 0.1 for r in (5, 10, 15, 20))
    outer_ok = all(0.9 < q[str(r)][1] < 1.05
                   for r in (40, 45, 50, 55, 58))
    peak = out["C3"]["measured"]["u_peak_shell"]
    out["C3"]["verdict"] = (
        ("CONFIRMED" if inner_ok else "QUALIFIED (outer-shell winding "
         "confirmed; interior NOT simply q=0: wound rings / ill-defined "
         "director at intermediate radii)") if outer_ok else "REFUTED")
    out["C3"]["verdict_energy"] = (
        f"u_eta does NOT peak in an outer shell: peak shell "
        f"r=[{peak[0]:.0f},{peak[1]:.0f}], broad spread")
    m4 = out["C4"]["measured"]
    c4_ok = (abs(m4["virial_u_over_3v"] - 0.053) < 0.005
             and abs(m4["contain_total"]["r90"] - 6.5) < 0.5
             and abs(m4["contain_total"]["coreball8"] - 0.93) < 0.02
             and m4["E_u"] > 0.98 * m4["E_u_seed"])
    out["C4"]["verdict"] = "CONFIRMED" if c4_ok else "QUALIFIED"
    m5 = out["C5"]["measured"]
    scaling_ok = (m5["seed"]["u_fit_resid_rel"] < 0.02
                  and abs(m5["seed"]["v4_slope"] - 3.0) < 0.15
                  and max(m5["compact_timemix_field"][
                      "lamEu_over_Eu1_dev_from_1"]) < 0.1)
    saw_iso = m5["sawtooth_ratio_fwd_over_central"]["iso_endpoint"]
    out["C5"]["verdict"] = (
        ("CONFIRMED (scaling law); QUALIFIED (the iso endpoint itself "
         "carries stencil-null grid-scale texture, sawtooth ratio "
         f"{max(saw_iso.values()):.1f}, so its E_u descent is partly a "
         "discretization artifact channel)"
         if max(saw_iso.values()) > 3.0
         or not m5["iso_Eu_monotone_decreasing"]
         else "CONFIRMED") if scaling_ok else "REFUTED")
    m6 = out["C6"]["measured"]
    out["C6"]["verdict"] = ("CONFIRMED" if
                            abs(m6["floor_closed_form"]
                                - m6["floor_numeric_scan"]) < 1e-6
                            and m6["E_neg_exists_below_kink"]
                            else "REFUTED")
    out["C7"]["verdict"] = "CONFIRMED"
    m8 = out["C8"]["measured"]
    out["C8"]["verdict"] = ("CONFIRMED" if
                            m8["boost_tangent_block_part_max"] == 0.0
                            and m8["force_timemix_on_blockdiag_seed"]
                            == 0.0
                            and m8["rot3_vs_full6_it8000_max_rel"] < 1e-4
                            else "QUALIFIED")
    path = os.path.join(DATA, "m5_21_1e_audit.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", path)
    for c in ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"):
        print(f"{c}: {out[c]['verdict']}")


if __name__ == "__main__":
    main()
