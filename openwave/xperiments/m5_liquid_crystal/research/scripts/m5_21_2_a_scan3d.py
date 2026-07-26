"""M5.21.2: the 3D 3-lepton axis-permutation scan (Duda 2026-07-17).

HIS PRESCRIPTION (verbatim anchors in tasks/m5_21_convo.md § 2026-07-17):
"first focus on 3D case, where the hedgehog is topologically protected"
(watch: "in small box and this huge delta=0.3, these topological
constraints might be destroyed by changing boundary conditions"; if so
reduce delta). "energy minimization starting from biaxial hedgehog
ansatz - for electron rotating diag(1, delta, 0). You can also try
energy minimization replacing the rotated to (delta,0,1) or (0, 1,
delta) - searching for local minima (e.g. just gradient descent,
global should be electron), hopefully getting candidates for 3
leptons."

THE 3D THEORY (the 4x4 verified stack with the time row dropped,
eta -> identity; see findings/m5_21_1_method_note.md for the 4D
provenance):
    M(x): 3x3 real symmetric tensor field on a cubic grid, h = 1.
    A_i = d_i M (central differences interior, one-sided edges).
    C_ij = [A_i, A_j]                       (antisymmetric)
    u = 4 * sum_{i<j} tr(C_ij^T C_ij)      (curvature density)
    V = WSCALE * sum_{p=1..3} (tr(M^p) - C_p)^2,  C_p = 1 + delta^p
        (vacuum spectrum (1, delta, 0); no g in 3D)
    E = sum_sites (u + V) * h^3
    Derrick in 3D: E_u ~ 1/lambda, E_V ~ lambda^3 -> a finite-size
    stationary point satisfies E_u = 3 E_V (the Faber virial balance).

SEEDS (the axis-permutation family): at x with r = |x|, frame
(rhat, phihat, that) (phihat = z x rhat azimuthal, that = closing the
triad; safe fallback near the z-axis), core-smoothed to the isotropic
a*I, a = (1+delta)/3, w = 1 - exp(-(r/r_c)^2):
    S = lam_r rhat rhat^T + lam_phi phihat phihat^T + lam_t that that^T
    seed A "electron": (lam_r, lam_phi, lam_t) = (1, delta, 0)
    seed B:            (delta, 0, 1)
    seed C:            (0, 1, delta)
N even -> cell-centered coords never hit the axis/origin exactly.

GATES (pre-registered, cap 3 tries):
    G0 complex-step gradient check <= 1e-10 rel (random small lattice
       + the seed)
    G1 internal SO(3) invariance: E(R M R^T) == E(M) <= 1e-10 rel,
       random R, + negative control (non-orthogonal Q must move E)
    G2 vacuum energy == 0 at machine precision
    G3 seed spectrum: far-field eigenvalues = the permuted vacuum;
       axis two-equal + center three-equal reads

MODES:  gates | ladder [combo] | scan [A|B|C] | plots | status
    ladder: seed A boundary-integrity probes, N in {32,48,64} x delta
        in {0.3,0.1,0.03} x BC in {pinned,free}, 800 iters, frame-
        overlap retention read. combo = "N,delta,BC" runs ONE rung
        (parallel launches); no combo = all rungs sequentially.
    scan: one long FIRE descent per seed at the working point
        (N = 48, delta from the ladder, pinned + free per the ladder
        verdict), snapshots -> film strips, census row on finish.
Per-run JSON merges (a crash cannot lose a finished run); explicit
mode args; default mode = status (argv-loss safe).

Out: ../data/m5_21_2_scan3d.json, ../plots/m5_21_2_*.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
JPATH = os.path.join(DATA, "m5_21_2_scan3d.json")

WSCALE = 0.000724023879 * float(os.environ.get("M5212_WMULT", "1"))
H = 1.0                          # base wscale carried from the 4D stack;
DELTA0 = 0.3                     # M5212_WMULT = the stiffness-arm knob
TAG = os.environ.get("M5212_TAG", "")   # output-name suffix for arms
STENCIL = os.environ.get("M5212_STENCIL", "2h")  # "2h" | "fwd"; fwd =
# compact forward differences, no odd-even null mode (the checkerboard
# catch, deviations log 2026-07-17)


# ================= grid + seeds =================
def coords(n):
    x = (np.arange(n) - (n - 1) / 2.0) * H
    return np.meshgrid(x, x, x, indexing="ij")


def frame(n):
    """(rhat, phihat, that) fields, shape (n,n,n,3); safe near axis."""
    X, Y, Z = coords(n)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.where(r < 1e-12, 1e-12, r)
    rhat = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):                       # exact-axis fallback (N odd)
        phihat[near] = np.array([0.0, 1.0, 0.0])
    that = np.cross(phihat, rhat)
    return r, rhat, phihat, that


SEEDS = {"A": (1.0, None, 0.0), "B": (None, 0.0, 1.0),
         "C": (0.0, 1.0, None)}           # None -> delta slot


def seed3(n, delta, which, r_c=4.0):
    lam = [delta if v is None else v for v in SEEDS[which]]
    r, rhat, phihat, that = frame(n)
    S = (lam[0] * rhat[..., :, None] * rhat[..., None, :]
         + lam[1] * phihat[..., :, None] * phihat[..., None, :]
         + lam[2] * that[..., :, None] * that[..., None, :])
    a = (1.0 + delta) / 3.0
    w = (1.0 - np.exp(-((r / r_c) ** 2)))[..., None, None]
    return w * S + (1.0 - w) * (a * np.eye(3))


def vac3(n, delta):
    M = np.zeros((n, n, n, 3, 3))
    M[:] = np.diag([1.0, delta, 0.0])
    return M


def seed_ring(n, delta, a=4.0, r_c=2.5):
    """P3 arm (user 2026-07-17): the charged disclination ring. Seed
    A's exact far field (unit monopole winding, so the pinned shell
    and the E comparison are sector-identical), core opened into a
    half-disclination ring of radius a in z = 0, interior escaped
    along +z. Meridional director angle
        psi = 0.5*(arg(zeta - ia) + arg(zeta + ia)),  zeta = z + i rho
    -> psi = theta far away, psi = 0 at the origin (no point defect),
    axis clean; the tensor nn^T is single-valued through the half
    winding. delta stays on phihat (seed A's assignment)."""
    X, Y, Z = coords(n)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    rhohat = np.stack([X / rhos, Y / rhos, np.zeros_like(Z)], axis=-1)
    psi = 0.5 * (np.arctan2(rho - a, Z) + np.arctan2(rho + a, Z))
    zhat = np.array([0.0, 0.0, 1.0])
    nvec = np.sin(psi)[..., None] * rhohat + np.cos(psi)[..., None] * zhat
    S = (nvec[..., :, None] * nvec[..., None, :]
         + delta * phihat[..., :, None] * phihat[..., None, :])
    d_ring = np.sqrt((rho - a) ** 2 + Z ** 2)
    av = (1.0 + delta) / 3.0
    w = (1.0 - np.exp(-((d_ring / r_c) ** 2)))[..., None, None]
    return w * S + (1.0 - w) * (av * np.eye(3))


def core_locus(M, delta):
    """where the potential-excess density peaks: ring vs point read."""
    n = M.shape[0]
    X, Y, Z = coords(n)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    M2 = M @ M
    t1 = np.einsum("...kk->...", M)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...kk->...", M2 @ M)
    cp = c_p(delta)
    vd = ((t1 - cp[0]) ** 2 + (t2 - cp[1]) ** 2 + (t3 - cp[2]) ** 2)
    i = np.unravel_index(int(np.argmax(vd)), vd.shape)
    edges = np.linspace(0, r.max() * 0.9, 20)
    prof = []
    for b in range(19):
        sel = (r >= edges[b]) & (r < edges[b + 1])
        prof.append(float(vd[sel].mean()) if np.any(sel) else 0.0)
    return {"vmax_rho": float(rho[i]), "vmax_z": float(Z[i]),
            "vmax_r": float(r[i]),
            "shell_peak_r": float(0.5 * (edges[int(np.argmax(prof))]
                                         + edges[int(np.argmax(prof)) + 1]))}


# ================= energy + gradient =================
def d_ax(f, ax):
    """derivative channel: 2h central (interior) + one-sided edges, or
    compact forward differences under STENCIL == "fwd" (h = 1)."""
    out = np.empty_like(f) if STENCIL == "2h" else np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    if STENCIL == "fwd":
        out[at(slice(0, -1))] = (f[at(slice(1, None))]
                                 - f[at(slice(0, -1))])
        return out
    out[at(slice(1, -1))] = (f[at(slice(2, None))]
                             - f[at(slice(0, -2))]) * 0.5
    out[at(0)] = f[at(1)] - f[at(0)]
    out[at(-1)] = f[at(-1)] - f[at(-2)]
    return out


def d_ax_adj(g, ax):
    """exact adjoint of d_ax (sum_i (d_ax f)_i g_i = sum_i f_i (adj g)_i)."""
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    if STENCIL == "fwd":
        out[at(slice(1, None))] += g[at(slice(0, -1))]
        out[at(slice(0, -1))] -= g[at(slice(0, -1))]
        return out
    out[at(slice(2, None))] += g[at(slice(1, -1))] * 0.5
    out[at(slice(0, -2))] -= g[at(slice(1, -1))] * 0.5
    out[at(1)] += g[at(0)]
    out[at(0)] -= g[at(0)]
    out[at(-1)] += g[at(-1)]
    out[at(-2)] -= g[at(-1)]
    return out


def c_p(delta):
    return [1.0 + delta, 1.0 + delta ** 2, 1.0 + delta ** 3]


def e_split(M, delta):
    """(E_u, E_V); complex-safe."""
    A = [d_ax(M, ax) for ax in range(3)]
    e_u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            e_u = e_u + 4.0 * np.sum(np.einsum("...kl,...kl->...",
                                               C, C))
    M2 = M @ M
    M3 = M2 @ M
    t1 = np.einsum("...kk->...", M)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...kk->...", M3)
    cp = c_p(delta)
    e_v = WSCALE * np.sum((t1 - cp[0]) ** 2 + (t2 - cp[1]) ** 2
                          + (t3 - cp[2]) ** 2)
    return e_u, e_v


def e_total(M, delta):
    e_u, e_v = e_split(M, delta)
    return e_u + e_v


def grad3(M, delta):
    """analytic gradient of e_total wrt M (symmetric-preserving)."""
    A = [d_ax(M, ax) for ax in range(3)]
    dA = [np.zeros_like(M) for _ in range(3)]
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            dA[i] += 8.0 * (C @ A[j] - A[j] @ C)
            dA[j] -= 8.0 * (C @ A[i] - A[i] @ C)
    G = np.zeros_like(M)
    for ax in range(3):
        G += d_ax_adj(dA[ax], ax)
    M2 = M @ M
    t1 = np.einsum("...kk->...", M)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...kk->...", M2 @ M)
    cp = c_p(delta)
    eye = np.eye(3)
    G += WSCALE * (2.0 * (t1 - cp[0])[..., None, None] * eye
                   + 4.0 * (t2 - cp[1])[..., None, None] * M
                   + 6.0 * (t3 - cp[2])[..., None, None] * M2)
    return G


def pin_shell(n):
    P = np.zeros((n, n, n), dtype=bool)
    P[0], P[-1] = True, True
    P[:, 0], P[:, -1] = True, True
    P[:, :, 0], P[:, :, -1] = True, True
    return P


# ================= diagnostics =================
def eig3(M):
    lam, V = np.linalg.eigh(M)
    return lam, V                       # ascending


def overlap_profile(M, nbins=24):
    """frame-overlap O_kl(r): shell-mean (v_k . e_l)^2, k = eig index
    (ascending), l in {rhat, phihat, that}; + eigenvalue profiles."""
    n = M.shape[0]
    r, rhat, phihat, that = frame(n)
    lam, V = eig3(M)
    refs = [rhat, phihat, that]
    rmax = r.max() * 0.95
    edges = np.linspace(0.0, rmax, nbins + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    O = np.zeros((nbins, 3, 3))
    L = np.zeros((nbins, 3))
    for b in range(nbins):
        sel = (r >= edges[b]) & (r < edges[b + 1])
        if not np.any(sel):
            continue
        for k in range(3):
            vk = V[sel][..., :, k]
            for l, ref in enumerate(refs):
                O[b, k, l] = float(np.mean(
                    np.sum(vk * ref[sel], axis=-1) ** 2))
        L[b] = lam[sel].mean(axis=0)
    return {"r": mid.tolist(), "O": O.tolist(), "lams": L.tolist()}


def retention(M, which, r_lo=8.0, r_hi=16.0):
    """scalar hedgehog-retention: mean (v_sigma(l) . e_l)^2 over the
    shell r_lo..r_hi, sigma = the seed's eigen-to-frame assignment."""
    n = M.shape[0]
    r, rhat, phihat, that = frame(n)
    lam, V = eig3(M)
    lam_seed = [0.3 if v is None else v for v in SEEDS[which]]
    order = np.argsort(lam_seed)        # frame slots by ascending lam
    refs = [rhat, phihat, that]
    sel = (r >= r_lo) & (r <= r_hi)
    vals = []
    for k in range(3):                  # eig index k (ascending)
        l = int(order[k])               # the frame slot ranking k-th
        vk = V[sel][..., :, k]
        vals.append(float(np.mean(np.sum(vk * refs[l][sel],
                                         axis=-1) ** 2)))
    return {"per_eig": vals, "mean": float(np.mean(vals))}


def axis_center_reads(M, delta):
    n = M.shape[0]
    X, Y, Z = coords(n)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam, _ = eig3(M)
    csel = r < 2.5
    core = lam[csel].mean(axis=0)
    rho = np.sqrt(X * X + Y * Y)
    asel = (rho < 1.2) & (np.abs(Z) > 6.0) & (np.abs(Z) < 14.0)
    ax = lam[asel].mean(axis=0) if np.any(asel) else [np.nan] * 3
    return {"core_lams": np.asarray(core).tolist(),
            "core_spread": float(np.max(core) - np.min(core)),
            "axis_lams": np.asarray(ax).tolist()}


def r_half(M, delta):
    """radius enclosing half the potential-excess energy."""
    n = M.shape[0]
    X, Y, Z = coords(n)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    M2 = M @ M
    t1 = np.einsum("...kk->...", M)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...kk->...", M2 @ M)
    cp = c_p(delta)
    dv = WSCALE * ((t1 - cp[0]) ** 2 + (t2 - cp[1]) ** 2
                   + (t3 - cp[2]) ** 2)
    order = np.argsort(r.ravel())
    cum = np.cumsum(dv.ravel()[order])
    if cum[-1] <= 0:
        return float("nan")
    k = int(np.searchsorted(cum, 0.5 * cum[-1]))
    return float(r.ravel()[order][min(k, len(order) - 1)])


# ================= FIRE =================
def fire3(M0, delta, free_mask, max_iter=8000, dt0=0.02, dt_max=0.2,
          log_every=500, snaps=(), tag="", f_tol=1e-8,
          plateau=(2000, 1e-10)):
    M = M0.copy()
    free = free_mask[..., None, None].astype(float)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist, states = [], [{"it": 0, "M": M0.copy()}]
    F = -grad3(M, delta) * free
    t0 = time.time()
    stop = "max_iter"
    for it in range(1, max_iter + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max)
                alpha *= 0.99
        else:
            v[:] = 0.0
            dt *= 0.5
            alpha = 0.1
            n_up = 0
        v += dt * F
        M += dt * v
        F = -grad3(M, delta) * free
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax):
            stop = "non-finite"
            break
        if it % log_every == 0 or it == max_iter:
            e_u, e_v = e_split(M, delta)
            E = e_u + e_v
            hist.append({"it": it, "E": float(E), "E_u": float(e_u),
                         "E_v": float(e_v), "fmax": fmax, "dt": dt})
            print(f"  {tag} it {it:6d} E {E:12.6f} u/3V "
                  f"{(e_u / max(3 * e_v, 1e-300)):8.3f} fmax {fmax:.3e}"
                  f" [{time.time() - t0:.0f}s]", flush=True)
            back = plateau[0] // max(log_every, 1)
            if len(hist) > back and \
                    abs(E - hist[-1 - back]["E"]) < plateau[1]:
                stop = "plateau"
                break
        if fmax < f_tol:
            stop = "f_tol"
            break
        if it in snaps:
            states.append({"it": it, "M": M.copy()})
    states.append({"it": -1, "M": M.copy()})
    return M, states, {"stop": stop, "trace": hist,
                       "wall_s": time.time() - t0}


# ================= json merge =================
def j_merge(update):
    cur = {}
    if os.path.exists(JPATH):
        with open(JPATH) as f:
            cur = json.load(f)
    def deep(a, b):
        for k, v in b.items():
            if isinstance(v, dict) and isinstance(a.get(k), dict):
                deep(a[k], v)
            else:
                a[k] = v
    deep(cur, update)
    os.makedirs(DATA, exist_ok=True)
    with open(JPATH, "w") as f:
        json.dump(cur, f, indent=1)


# ================= gates =================
def gates():
    rng = np.random.default_rng(7)
    out = {}
    # G0 complex-step on a random small lattice + the seed
    for name, M in (("rand", rng.normal(size=(8, 8, 8, 3, 3))),
                    ("seed", seed3(16, DELTA0, "A"))):
        M = 0.5 * (M + M.swapaxes(-1, -2))
        G = grad3(M, DELTA0)
        errs = []
        for _ in range(4):
            V = rng.normal(size=M.shape)
            V = 0.5 * (V + V.swapaxes(-1, -2))
            eps = 1e-30
            de_cs = e_total(M + 1j * eps * V, DELTA0).imag / eps
            de_an = float(np.sum(G * V))
            errs.append(abs(de_cs - de_an) / max(abs(de_cs), 1e-300))
        out[f"g0_{name}_relerr_max"] = float(np.max(errs))
    # G1 internal SO(3) invariance + negative control
    M = seed3(24, DELTA0, "A")
    E0 = float(e_total(M, DELTA0))
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    E1 = float(e_total(np.einsum("ab,...bc,dc->...ad", Q, M, Q),
                       DELTA0))
    out["g1_so3_relerr"] = abs(E1 - E0) / abs(E0)
    Qbad = Q + 0.05 * rng.normal(size=(3, 3))
    E2 = float(e_total(np.einsum("ab,...bc,dc->...ad", Qbad, M, Qbad),
                       DELTA0))
    out["g1_negctrl_relshift"] = abs(E2 - E0) / abs(E0)
    # G2 vacuum
    out["g2_vac_E"] = float(e_total(vac3(16, DELTA0), DELTA0))
    # G3 seed spectra
    for which in ("A", "B", "C"):
        M = seed3(48, DELTA0, which)
        lam, _ = eig3(M)
        n = M.shape[0]
        X, Y, Z = coords(n)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        far = (r > 14) & (r < 20)
        out[f"g3_far_lams_{which}"] = lam[far].mean(axis=0).tolist()
        out[f"g3_axis_center_{which}"] = axis_center_reads(M, DELTA0)
    out["g0_pass"] = max(out["g0_rand_relerr_max"],
                         out["g0_seed_relerr_max"]) < 1e-10
    out["g1_pass"] = (out["g1_so3_relerr"] < 1e-10
                      and out["g1_negctrl_relshift"] > 1e-6)
    out["g2_pass"] = abs(out["g2_vac_E"]) < 1e-12
    j_merge({"gates": out})
    print(json.dumps(out, indent=1))
    return out


# ================= P1 boundary ladder =================
LADDER_NS = (32, 48, 64)
LADDER_DS = (0.3, 0.1, 0.03)
LADDER_BCS = ("pinned", "free")


def ladder_one(n, delta, bc, iters=800):
    tag = f"L n{n} d{delta} {bc}"
    M0 = seed3(n, delta, "A")
    free = ~pin_shell(n) if bc == "pinned" else \
        np.ones((n, n, n), dtype=bool)
    ret0 = retention(M0, "A", r_lo=n / 6, r_hi=n / 3)
    e_u0, e_v0 = e_split(M0, delta)
    M, _, info = fire3(M0, delta, free, max_iter=iters, log_every=400,
                       tag=tag)
    ret1 = retention(M, "A", r_lo=n / 6, r_hi=n / 3)
    e_u1, e_v1 = e_split(M, delta)
    row = {"n": n, "delta": delta, "bc": bc, "iters": iters,
           "ret_seed": ret0["mean"], "ret_end": ret1["mean"],
           "ret_per_eig_end": ret1["per_eig"],
           "E_seed": float(e_u0 + e_v0), "E_end": float(e_u1 + e_v1),
           "u_over_3v_end": float(e_u1 / max(3 * e_v1, 1e-300)),
           "core_end": axis_center_reads(M, delta),
           "stop": info["stop"], "wall_s": info["wall_s"]}
    with open(os.path.join(DATA,
              f"m5_21_2_lad_n{n}_d{delta}_{bc}.json"), "w") as f:
        json.dump(row, f, indent=1)      # per-run file: no merge race
    print(json.dumps(row))
    return row


def ladder(combo=None):
    if combo:
        n_s, d_s, bc = combo.split(",")
        return ladder_one(int(n_s), float(d_s), bc)
    for n in LADDER_NS:
        for d in LADDER_DS:
            for bc in LADDER_BCS:
                ladder_one(n, d, bc)


# ================= P2 the 3-seed census =================
def scan(which, n=48, delta=DELTA0, bc="pinned",
         max_iter=int(os.environ.get("M5212_MAXIT", "24000")),
         a_ring=4.0):
    tag = f"S{which} n{n} d{delta} {bc}"
    M0 = seed_ring(n, delta, a=a_ring) if which == "R" else \
        seed3(n, delta, which)
    ret_key = "A" if which == "R" else which
    free = ~pin_shell(n) if bc == "pinned" else \
        np.ones((n, n, n), dtype=bool)
    snaps = (500, 2000, 8000, 16000)
    M, states, info = fire3(M0, delta, free, max_iter=max_iter,
                            log_every=500, snaps=snaps, tag=tag)
    e_u, e_v = e_split(M, delta)
    row = {"seed": which, "n": n, "delta": delta, "bc": bc,
           "E_end": float(e_u + e_v), "E_u": float(e_u),
           "E_v": float(e_v),
           "u_over_3v": float(e_u / max(3 * e_v, 1e-300)),
           "r_half": r_half(M, delta),
           "retention": retention(M, ret_key, r_lo=8, r_hi=16),
           "core": axis_center_reads(M, delta),
           "core_locus": core_locus(M, delta),
           "core_locus_seed": core_locus(M0, delta),
           "overlap_end": overlap_profile(M),
           "overlap_seed": overlap_profile(M0),
           "stop": info["stop"], "trace": info["trace"][-8:],
           "wall_s": info["wall_s"]}
    row["wscale"] = WSCALE
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_2_end_{which}_{bc}{TAG}.npz"),
        M=M.astype(np.float32), delta=delta)
    with open(os.path.join(DATA,
              f"m5_21_2_scanrow_{which}_{bc}{TAG}.json"), "w") as f:
        json.dump(row, f, indent=1)      # row BEFORE film: a film
    film(states, which, bc + TAG, n, delta)  # crash cannot lose the run
    print(json.dumps({k: row[k] for k in
                      ("seed", "bc", "E_end", "E_u", "E_v",
                       "u_over_3v", "r_half", "stop")}))
    return row


# ================= films (3x3 meridional strips) =================
def film(states, which, bc, n, delta):
    """basic strip: eigenvalue maps + leading-alignment map per
    snapshot; thermal strip: u / V density maps (the 4x4 templates
    adapted to the 3x3 stack; y = 0 meridional slice)."""
    j = n // 2                       # slice nearest y=0 (cell-centered)
    X, Y, Z = coords(n)
    r, rhat, _, _ = frame(n)
    cols = len(states)
    for mode in ("basic", "thermal"):
        rows = 4 if mode == "basic" else 2
        fig, axs = plt.subplots(rows, cols,
                                figsize=(2.1 * cols, 2.1 * rows))
        axs = np.atleast_2d(axs)
        for c, st in enumerate(states):
            M = st["M"]
            Ms = M[:, j]
            if mode == "basic":
                lam, V = eig3(Ms)
                al = np.sum(V[..., :, 2] * rhat[:, j], axis=-1) ** 2
                panels = [lam[..., 2], lam[..., 1], lam[..., 0], al]
                labs = ["lam_max", "lam_mid", "lam_min",
                        "(v_max.rhat)^2"]
                vlims = [(0, 1.2), (0, 0.6), (-0.2, 0.4), (0, 1)]
            else:
                A = [d_ax(M, ax) for ax in range(3)]
                ud = np.zeros(M.shape[:3])
                for i in range(3):
                    for k in range(i + 1, 3):
                        C = A[i] @ A[k] - A[k] @ A[i]
                        ud += 4.0 * np.einsum("...kl,...kl->...", C, C)
                M2 = Ms @ Ms
                t1 = np.einsum("...kk->...", Ms)
                t2 = np.einsum("...kk->...", M2)
                t3 = np.einsum("...kk->...", M2 @ Ms)
                cp = c_p(delta)
                vd = WSCALE * ((t1 - cp[0]) ** 2 + (t2 - cp[1]) ** 2
                               + (t3 - cp[2]) ** 2)
                panels = [np.log10(ud[:, j] + 1e-12),
                          np.log10(vd + 1e-12)]
                labs = ["log10 u", "log10 V"]
                vlims = [(-9, -1), (-9, -3)]
            for rw in range(rows):
                ax = axs[rw, c]
                im = ax.imshow(panels[rw].T, origin="lower",
                               cmap="viridis", vmin=vlims[rw][0],
                               vmax=vlims[rw][1])
                ax.set_xticks([]); ax.set_yticks([])
                if c == 0:
                    ax.set_ylabel(labs[rw], fontsize=7)
                if rw == 0:
                    it = st["it"]
                    ax.set_title("end" if it < 0 else f"it {it}",
                                 fontsize=7)
                fig.colorbar(im, ax=ax, fraction=0.045)
        fig.suptitle(f"M5.21.2 seed {which} {bc} n{n} d{delta} "
                     f"({mode} film, y=0 slice)", fontsize=9)
        fig.tight_layout()
        os.makedirs(PLOTS, exist_ok=True)
        fig.savefig(os.path.join(
            PLOTS, f"m5_21_2_film_{mode}_{which}_{bc}.png"), dpi=110)
        plt.close(fig)


# ================= plots =================
def plots():
    with open(JPATH) as f:
        J = json.load(f)
    lad = J.get("ladder", {})
    scn = J.get("scan", {})
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    # ladder heat: retention end per (n, delta), pinned vs free
    for a, bc in enumerate(LADDER_BCS):
        Ht = np.full((len(LADDER_NS), len(LADDER_DS)), np.nan)
        for i, n in enumerate(LADDER_NS):
            for k, d in enumerate(LADDER_DS):
                row = lad.get(f"n{n}_d{d}_{bc}")
                if row:
                    Ht[i, k] = row["ret_end"]
        ax = axs[0, a]
        im = ax.imshow(Ht, origin="lower", cmap="RdYlGn",
                       vmin=0.33, vmax=1.0)
        ax.set_xticks(range(len(LADDER_DS)), LADDER_DS)
        ax.set_yticks(range(len(LADDER_NS)), LADDER_NS)
        ax.set_xlabel("delta"); ax.set_ylabel("N")
        ax.set_title(f"retention end ({bc})")
        for i in range(len(LADDER_NS)):
            for k in range(len(LADDER_DS)):
                if np.isfinite(Ht[i, k]):
                    ax.text(k, i, f"{Ht[i, k]:.2f}", ha="center",
                            va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.045)
    # census bars
    ax = axs[0, 2]
    keys = sorted(scn.keys())
    if keys:
        E = [scn[k]["E_end"] for k in keys]
        ax.bar(range(len(keys)), E)
        ax.set_xticks(range(len(keys)), keys, fontsize=7)
        ax.set_title("census: E per seed endpoint")
    # energy traces
    ax = axs[1, 0]
    for k in keys:
        tr = scn[k].get("trace", [])
        ax.plot([t["it"] for t in tr], [t["E"] for t in tr],
                label=k)
    ax.set_title("E trace (tail)"); ax.legend(fontsize=7)
    # virial
    ax = axs[1, 1]
    if keys:
        ax.bar(range(len(keys)),
               [scn[k]["u_over_3v"] for k in keys])
        ax.axhline(1.0, color="k", ls="--")
        ax.set_xticks(range(len(keys)), keys, fontsize=7)
        ax.set_title("u / 3V at endpoint (1 = Faber virial)")
    # overlap profile of the first scan endpoint
    ax = axs[1, 2]
    if keys:
        op = scn[keys[0]]["overlap_end"]
        r = op["r"]
        O = np.array(op["O"])
        for k in range(3):
            ax.plot(r, O[:, k, k % 3], label=f"O[{k},{k % 3}]")
        L = np.array(op["lams"])
        for k in range(3):
            ax.plot(r, L[:, k], ls=":", label=f"lam{k}")
        ax.set_title(f"{keys[0]}: overlap + lam profiles")
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_2_panel.png"), dpi=110)
    plt.close(fig)
    print("plots written")


def collect():
    """merge the per-run row files into the main JSON (race-free)."""
    import glob
    upd = {"ladder": {}, "scan": {}}
    for p in glob.glob(os.path.join(DATA, "m5_21_2_lad_*.json")):
        key = os.path.basename(p)[len("m5_21_2_lad_"):-len(".json")]
        with open(p) as f:
            upd["ladder"][key] = json.load(f)
    for p in glob.glob(os.path.join(DATA, "m5_21_2_scanrow_*.json")):
        key = os.path.basename(p)[len("m5_21_2_scanrow_"):-len(".json")]
        with open(p) as f:
            upd["scan"][key] = json.load(f)
    j_merge({k: v for k, v in upd.items() if v})
    print(f"collected {len(upd['ladder'])} ladder + "
          f"{len(upd['scan'])} scan rows")


def status():
    if os.path.exists(JPATH):
        with open(JPATH) as f:
            J = json.load(f)
        print(json.dumps({k: sorted(v.keys()) if isinstance(v, dict)
                          else v for k, v in J.items()}, indent=1))
    else:
        print("no json yet")


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0] if ARGV else "status"
    if mode == "gates":
        gates()
    elif mode == "ladder":
        ladder(ARGV[1] if len(ARGV) > 1 else None)
    elif mode == "scan":
        which = ARGV[1] if len(ARGV) > 1 else "A"
        bc = ARGV[2] if len(ARGV) > 2 else "pinned"
        n = int(ARGV[3]) if len(ARGV) > 3 else 48
        d = float(ARGV[4]) if len(ARGV) > 4 else DELTA0
        a_r = float(ARGV[5]) if len(ARGV) > 5 else 4.0
        scan(which, n=n, delta=d, bc=bc, a_ring=a_r)
    elif mode == "collect":
        collect()
    elif mode == "plots":
        collect()
        plots()
    else:
        status()
