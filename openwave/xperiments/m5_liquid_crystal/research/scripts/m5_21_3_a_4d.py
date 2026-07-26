"""M5.21.3: the 4D lift: does energy minimization select nonzero
field time derivatives? (his 2026-07-18 prescription).

THE STACK (the audited conventions, joined):
  - grid/stencil layer = M5.21.2b (sym stencil = 1/2(fwd+bwd) with
    exact adjoints, h-aware, physical box L = N*h, pinned shell)
  - eta layer = the M5.20.2 production convention (audited M5.21.1):
    ETA = diag(-1,1,1,1); comm_eta(A,B) = A eta B - B eta A;
    <F,G>_eta = tr(eta F eta G^T);
    u = 4 sum_{mu<nu} <F_munu, F_munu>_eta
  - vacuum (branch s): M_vac = diag(-sg, 1, delta, 0), sg = s*g;
    eta.M_vac spectrum = (sg, 1, delta, 0)
  - potential V4 = W1 * sum_{p=1..4} (tr((M eta)^p) - C_p)^2 with
    C_p = sg^p + 1 + delta^p (the audited trace-target; complex-step
    safe end-to-end; the 2b T2-eigenvalue lesson is a documented
    caveat, its Lorentz-invariant lift = future arm)

THE TIME SECTOR (this task's new piece): the velocity-field reading
of "energy minimization including field time derivatives": for a
chosen normalized clock direction field a0(x) (generator G, envelope
w(r)), the config M(x, t) with dM/dt = omega * a0 contributes
  E_kin(omega) = omega^2 * kin(M; G),
  kin = 4 sum_i <comm_eta(a0, A_i), comm_eta(a0, A_i)>_eta
(pure quadratic in omega: no linear term exists, F_0i ~ omega).
The eta-metric makes kin SIGN-INDEFINITE: kin < 0 is the signature
channel that would drive nonzero time derivatives; kin > 0 for all
generators means static wins at toy parameters. The z-twist channel
(theta = k z) enters A_z -> A_z + k a0: E(k) = E + b k + d k^2 with
a LINEAR term b (spontaneous twist selection possible).

MODES: gates | p1 s=+1|-1 [maxit] | hess tagbase | p2 tagbase |
       p3 tagbase gen=<name> [omegas] | collect | status
Out: ../data/m5_21_3_*.json / npz (race-free rows).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

W1 = 0.000724023879                  # the 2b/M5.21.2 WSCALE, carried
G_T = 8.0                            # the engine LC_G (m5_20_2_a_eom)
DELTA0 = 0.3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
SEED3_NPZ = os.path.join(DATA, "m5_21_2b_end_i2_A_T2.npz")


# ================= stencils (copied from m5_21_2b, provenance) ======
def d1(f, ax, h, st):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    if st == "fwd":
        out[at(slice(0, -1))] = (f[at(slice(1, None))]
                                 - f[at(slice(0, -1))]) / h
    elif st == "bwd":
        out[at(slice(1, None))] = (f[at(slice(1, None))]
                                   - f[at(slice(0, -1))]) / h
    else:
        out[at(slice(1, -1))] = (f[at(slice(2, None))]
                                 - f[at(slice(0, -2))]) / (2 * h)
        out[at(0)] = (f[at(1)] - f[at(0)]) / h
        out[at(-1)] = (f[at(-1)] - f[at(-2)]) / h
    return out


def d1_adj(g, ax, h, st):
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    if st == "fwd":
        out[at(slice(1, None))] += g[at(slice(0, -1))] / h
        out[at(slice(0, -1))] -= g[at(slice(0, -1))] / h
    elif st == "bwd":
        out[at(slice(1, None))] += g[at(slice(1, None))] / h
        out[at(slice(0, -1))] -= g[at(slice(1, None))] / h
    else:
        out[at(slice(2, None))] += g[at(slice(1, -1))] / (2 * h)
        out[at(slice(0, -2))] -= g[at(slice(1, -1))] / (2 * h)
        out[at(1)] += g[at(0)] / h
        out[at(0)] -= g[at(0)] / h
        out[at(-1)] += g[at(-1)] / h
        out[at(-2)] -= g[at(-1)] / h
    return out


def branches(st):
    return [("fwd", 0.5), ("bwd", 0.5)] if st == "sym" else [(st, 1.0)]


def coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def pin_shell(n, h, depth=1.6):
    wc = max(1, int(np.ceil(depth / h)))
    P = np.zeros((n, n, n), dtype=bool)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = slice(0, wc); P[tuple(sl)] = True
        sl[ax] = slice(n - wc, n); P[tuple(sl)] = True
    return P


# ================= eta algebra =================
def comm_eta(A, B):
    return A @ ETA @ B - B @ ETA @ A


def inner_eta(F, G):
    """<F,G>_eta per cell = tr(eta F eta G^T)."""
    return np.einsum("...ab,...cd,ac,bd->...", F, G, ETA, ETA,
                     optimize=True)


def sym4(X):
    return 0.5 * (X + X.swapaxes(-1, -2))


# ================= config / seeds =================
def base_cfg(**kw):
    cfg = {"s": 1.0, "g": G_T, "n": 32, "L": 48.0, "delta": DELTA0,
           "stencil": "sym", "maxit": 6000, "tag": "", "renv": 10.0}
    cfg.update(kw)
    cfg["h"] = cfg["L"] / cfg["n"]
    cfg["sg"] = cfg["s"] * cfg["g"]
    return cfg


def vac4(cfg):
    return np.diag([-cfg["sg"], 1.0, cfg["delta"], 0.0])


def c4_of(cfg):
    sg, d = cfg["sg"], cfg["delta"]
    return tuple(sg ** p + 1.0 + d ** p for p in range(1, 5))


def load_seed3():
    Z = np.load(SEED3_NPZ)
    return Z["M"].astype(np.float64)


def embed34(M3, cfg):
    n = M3.shape[0]
    M4 = np.zeros((n, n, n, 4, 4))
    M4[..., 1:4, 1:4] = M3
    M4[..., 0, 0] = -cfg["sg"]
    return M4


def envelope(cfg):
    X, Y, Z = coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return np.exp(-((r / cfg["renv"]) ** 4))


# ================= energy: spatial u + V4 =================
def a_fields(M, cfg, st=None):
    st = st or cfg["stencil"]
    return {br: ([d1(M, ax, cfg["h"], br) for ax in range(3)], wt)
            for br, wt in branches(st)}


def e_parts(M, cfg, st=None):
    """(E_u, E_V) h^3-weighted (no eps arm in 4D: eps = 0 era)."""
    h3 = cfg["h"] ** 3
    e_u = 0.0
    for br, (A, wt) in a_fields(M, cfg, st).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                e_u = e_u + wt * 4.0 * np.sum(inner_eta(F, F))
    Me = M @ ETA
    P = Me
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    cp = c4_of(cfg)
    vd = sum((t[p] - cp[p]) ** 2 for p in range(4))
    return h3 * e_u, h3 * W1 * np.sum(vd)


def e_total(M, cfg, st=None):
    a, b = e_parts(M, cfg, st)
    return a + b


def grad(M, cfg):
    """exact gradient of e_total wrt symmetric M (complex-step gated)."""
    h3 = cfg["h"] ** 3
    G = np.zeros_like(M)
    for br, wt in branches(cfg["stencil"]):
        A = [d1(M, ax, cfg["h"], br) for ax in range(3)]
        dA = [np.zeros_like(M) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                WF = 8.0 * (ETA @ F @ ETA)          # dU/dF (F antisym)
                # dF = dAi eta Aj - Aj eta dAi  (+ Ai eta dAj - dAj eta Ai)
                dA[i] += WF @ (ETA @ A[j]).swapaxes(-1, -2) \
                    - (A[j] @ ETA).swapaxes(-1, -2) @ WF
                dA[j] += (A[i] @ ETA).swapaxes(-1, -2) @ WF \
                    - WF @ (ETA @ A[i]).swapaxes(-1, -2)
        for ax in range(3):
            G += wt * d1_adj(dA[ax], ax, cfg["h"], br)
    # V4 part
    Me = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for p in range(1, 4):
        pows.append(pows[-1] @ Me)
    t = [np.einsum("...kk->...", P @ Me) for P in pows]
    cp = c4_of(cfg)
    GV = np.zeros_like(M)
    for p in range(1, 5):
        coef = 2.0 * W1 * (t[p - 1] - cp[p - 1]) * p
        X = ETA @ pows[p - 1]                     # d tr((M eta)^p)/dM^T
        GV += coef[..., None, None] * X.swapaxes(-1, -2)
    return h3 * G + h3 * sym4(GV)


# ================= kinetic sector =================
def gen_catalog(cfg, M):
    """named generator/velocity fields a0(x) (normalized below)."""
    n = cfg["n"]
    w = envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    out = {}

    def local_rot(vhat):
        W = np.zeros(vhat.shape[:-1] + (4, 4))
        n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
        W[..., 1, 2], W[..., 1, 3] = -n3, n2
        W[..., 2, 1], W[..., 2, 3] = n3, -n1
        W[..., 3, 1], W[..., 3, 2] = -n2, n1
        return W

    # clock_local: rotation about the local LEADING spatial eigenvector
    out["clock_local"] = local_rot(V[..., :, 2])
    # plane_1d: rotation about the local LOWEST eigenvector (rotates
    # the (1, delta)-eigenplane into itself)
    out["plane_1d"] = local_rot(V[..., :, 0])
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4)); Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0
    for nm, Gm in (("rot_z", Jz), ("rot_x", Jx),
                   ("boost_z", Kz), ("boost_x", Kx)):
        out[nm] = np.broadcast_to(Gm, M.shape)
    a0s = {}
    for nm, Gm in out.items():
        a0 = w * (Gm @ M - M @ Gm.swapaxes(-1, -2))
        nrm = np.sqrt(np.sum(a0 * a0))
        a0s[nm] = a0 / max(nrm, 1e-300)          # unit Frobenius norm
    return a0s


def kin_of(M, a0, cfg, st=None):
    """kin(M; a0) = 4 sum_i <comm_eta(a0, A_i)>^2_eta, h^3-weighted.
    E_kin(omega) = omega^2 * kin."""
    h3 = cfg["h"] ** 3
    k = 0.0
    for br, (A, wt) in a_fields(M, cfg, st).items():
        for i in range(3):
            F = comm_eta(a0, A[i])
            k = k + wt * 4.0 * np.sum(inner_eta(F, F))
    return h3 * k


def kin_grad(M, a0, cfg, da0_of=None):
    """gradient of kin wrt M with a0 FROZEN (the velocity-field read;
    the a0(M) dependence is deliberately not chained: the clock
    direction is the ansatz, complex-step gated in this frozen form)."""
    h3 = cfg["h"] ** 3
    G = np.zeros_like(M)
    for br, wt in branches(cfg["stencil"]):
        A = [d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            F = comm_eta(a0, A[i])
            WF = 8.0 * (ETA @ F @ ETA)
            dAi = (a0 @ ETA).swapaxes(-1, -2) @ WF \
                - WF @ (ETA @ a0).swapaxes(-1, -2)
            G += wt * d1_adj(dAi, i, cfg["h"], br)
    return h3 * sym4(G)


def twist_read(M, a0, cfg, ks=(-0.02, -0.01, 0.0, 0.01, 0.02)):
    """E(k): A_z -> A_z + k*a0 (the axial-twist channel); fit
    E = e0 + b k + d k^2."""
    h3 = cfg["h"] ** 3
    es = []
    for k in ks:
        e_u = 0.0
        for br, (A, wt) in a_fields(M, cfg).items():
            Az = A[2] + k * a0
            AA = [A[0], A[1], Az]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = comm_eta(AA[i], AA[j])
                    e_u = e_u + wt * 4.0 * np.sum(inner_eta(F, F))
        es.append(h3 * e_u)
    co = np.polyfit(np.array(ks), np.array(es), 2)
    d, b, e0 = float(co[0]), float(co[1]), float(co[2])
    kstar = -b / (2.0 * d) if abs(d) > 1e-300 else float("nan")
    return {"b": b, "d": d, "k_star": kstar, "E_at_k": es,
            "dE_at_kstar": float(b * kstar + d * kstar ** 2)
            if np.isfinite(kstar) else float("nan")}


# ================= FIRE (4x4; omega-augmented optional) ============
def fire(M0, cfg, free_mask, max_iter, a0=None, omega=0.0,
         log_every=500, tag="", f_tol=1e-8, plateau=(2000, 1e-10),
         dt0=0.02, dt_max=0.2, dive_floor=None):
    M = M0.copy()
    free = free_mask[..., None, None].astype(float)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist = []

    def tot_grad(Mx):
        Gt = grad(Mx, cfg)
        if a0 is not None and omega != 0.0:
            Gt = Gt + (omega ** 2) * kin_grad(Mx, a0, cfg)
        return Gt

    def tot_e(Mx):
        e = e_total(Mx, cfg)
        if a0 is not None and omega != 0.0:
            e = e + (omega ** 2) * kin_of(Mx, a0, cfg)
        return e
    F = -tot_grad(M) * free
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
        F = -tot_grad(M) * free
        fmax = float(np.max(np.abs(F)))
        if not np.isfinite(fmax):
            stop = "non-finite"
            break
        if it % log_every == 0 or it == max_iter:
            E = float(tot_e(M))
            hist.append({"it": it, "E": E, "fmax": fmax, "dt": dt})
            print(f"  {tag} it {it:6d} E {E:14.6f} fmax {fmax:.3e} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
            if dive_floor is not None and E < dive_floor:
                stop = "dive"
                break
            back = plateau[0] // max(log_every, 1)
            if len(hist) > back and \
                    abs(E - hist[-1 - back]["E"]) < plateau[1]:
                stop = "plateau"
                break
        if fmax < f_tol:
            stop = "f_tol"
            break
    return M, {"stop": stop, "trace": hist,
               "wall_s": time.time() - t0}


# ================= gates =================
def gates():
    rng = np.random.default_rng(21)
    out = {}
    cfg = base_cfg(n=10, L=15.0)
    # G0: complex-step gradient of the full static functional
    M = rng.normal(size=(10, 10, 10, 4, 4))
    M = sym4(M)
    G = grad(M, cfg)
    errs = []
    for _ in range(4):
        V = sym4(rng.normal(size=M.shape))
        de_an = float(np.sum(G * V))
        t = 1e-30
        de = e_total(M + 1j * t * V, cfg).imag / t
        errs.append(abs(de - de_an) / max(abs(de), 1e-300))
    out["g0_static"] = float(np.max(errs))
    # G0k: complex-step on the kinetic gradient (a0 frozen)
    a0 = rng.normal(size=M.shape)
    a0 = a0 / np.sqrt(np.sum(a0 * a0))
    Gk = kin_grad(M, a0, cfg)
    errs = []
    for _ in range(4):
        V = sym4(rng.normal(size=M.shape))
        de_an = float(np.sum(Gk * V))
        t = 1e-30
        de = kin_of(M + 1j * t * V, a0, cfg).imag / t
        errs.append(abs(de - de_an) / max(abs(de), 1e-300))
    out["g0_kin"] = float(np.max(errs))
    # G1: SO(1,3) invariance of E_static + negative control
    from scipy.linalg import expm
    cfgB = base_cfg(n=14, L=21.0)
    M3 = load_seed3()
    sub = M3[::2, ::2, ::2][:14, :14, :14]
    MB = np.zeros((14, 14, 14, 4, 4))
    MB[..., 1:4, 1:4] = sub
    MB[..., 0, 0] = -cfgB["sg"]
    E0 = e_total(MB, cfgB)
    Gm = np.zeros((4, 4))
    Gm[0, 1] = Gm[1, 0] = 0.11          # boost x
    Gm[2, 3], Gm[3, 2] = -0.23, 0.23    # rotation
    L = expm(Gm)
    assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-12)
    ML = np.einsum("ab,...bc,dc->...ad", L, MB, L)
    out["g1_so13"] = abs(float(e_total(ML, cfgB)) - E0) / abs(E0)
    Lb = L + 0.05 * rng.normal(size=(4, 4))
    MLb = np.einsum("ab,...bc,dc->...ad", Lb, MB, Lb)
    out["g1_negctrl"] = abs(float(e_total(MLb, cfgB)) - E0) / abs(E0)
    # G2: vacuum == 0 (both branches)
    for s in (1.0, -1.0):
        c2 = base_cfg(s=s, n=8, L=12.0)
        Mv = np.zeros((8, 8, 8, 4, 4))
        Mv[:] = vac4(c2)
        eu, ev = e_parts(Mv, c2)
        out[f"g2_vac_s{int(s):+d}"] = [float(eu), float(ev)]
    # G3: 3D regression: block-diag embed => E_u(4D) == E_u(3D read)
    cfg3 = base_cfg(n=32, L=48.0)
    M3f = load_seed3()
    M4 = embed34(M3f, cfg3)
    eu4, _ = e_parts(M4, cfg3)
    # 3D read with the 2b machinery: plain commutators on the 3x3 block
    h = cfg3["h"]
    e3 = 0.0
    for br, wt in branches("sym"):
        A = [d1(M3f, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                e3 += wt * 4.0 * np.sum(
                    np.einsum("...kl,...kl->...", C, C))
    e3 *= h ** 3
    out["g3_regression"] = abs(float(eu4) - float(e3)) \
        / max(abs(float(e3)), 1e-300)
    out["pass"] = bool(out["g0_static"] < 5e-9 and out["g0_kin"] < 5e-9
                       and out["g1_so13"] < 1e-9
                       and out["g1_negctrl"] > 1e-6
                       and max(abs(x) for v in
                               (out["g2_vac_s+1"], out["g2_vac_s-1"])
                               for x in v) < 1e-16
                       and out["g3_regression"] < 1e-12)
    with open(os.path.join(DATA, "m5_21_3_gates.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


# ================= phases =================
def p1(s, maxit=6000):
    cfg = base_cfg(s=float(s), maxit=int(maxit),
                   tag=f"p1_s{int(float(s)):+d}")
    M0 = embed34(load_seed3(), cfg)
    free = ~pin_shell(cfg["n"], cfg["h"])
    e0u, e0v = e_parts(M0, cfg)
    M, info = fire(M0, cfg, free, cfg["maxit"], tag=cfg["tag"])
    eu, ev = e_parts(M, cfg)
    offb = float(np.max(np.abs(M[..., 0, 1:4])))
    lam_c = np.linalg.eigvalsh(
        M[cfg["n"] // 2, cfg["n"] // 2, cfg["n"] // 2, 1:4, 1:4])
    row = {"tag": cfg["tag"], "s": cfg["s"], "E_end": float(eu + ev),
           "E_u": float(eu), "E_v": float(ev),
           "E_seed": float(e0u + e0v), "offblock_max": offb,
           "core_lam3": lam_c.tolist(), "stop": info["stop"],
           "trace": info["trace"][-5:], "wall_s": info["wall_s"]}
    np.savez_compressed(os.path.join(DATA, f"m5_21_3_{cfg['tag']}.npz"),
                        M=M.astype(np.float32), s=cfg["s"],
                        delta=cfg["delta"], h=cfg["h"])
    with open(os.path.join(DATA, f"m5_21_3_row_{cfg['tag']}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in ("tag", "E_end", "E_u", "E_v",
                                          "offblock_max", "stop")}))
    return row


def load_p1(tagbase):
    Z = np.load(os.path.join(DATA, f"m5_21_3_{tagbase}.npz"))
    cfg = base_cfg(s=float(Z["s"]), tag=str(tagbase))
    return Z["M"].astype(np.float64), cfg


def hess(tagbase, ndir=24, t=1e-3):
    """time-mixing directional curvatures at the P1 endpoint (central
    FD, unit-norm localized directions with ONLY 0i entries)."""
    M, cfg = load_p1(tagbase)
    rng = np.random.default_rng(31)
    w = envelope(cfg)
    E0 = e_total(M, cfg)
    curv = []
    for _ in range(ndir):
        P = np.zeros_like(M)
        for i in range(1, 4):
            f = rng.normal(size=M.shape[:3]) * w
            P[..., 0, i] = f
            P[..., i, 0] = f
        P /= np.sqrt(np.sum(P * P))
        c = (e_total(M + t * P, cfg) - 2 * E0
             + e_total(M - t * P, cfg)) / t ** 2
        curv.append(float(c))
    # spatial-sector control directions (block-diagonal)
    ctrl = []
    for _ in range(8):
        P = np.zeros_like(M)
        B = sym4(rng.normal(size=M.shape[:3] + (3, 3))
                 if False else rng.normal(size=M.shape[:3] + (3, 3)))
        B = 0.5 * (B + B.swapaxes(-1, -2))
        P[..., 1:4, 1:4] = B * w[..., None, None]
        P /= np.sqrt(np.sum(P * P))
        c = (e_total(M + t * P, cfg) - 2 * E0
             + e_total(M - t * P, cfg)) / t ** 2
        ctrl.append(float(c))
    out = {"tag": f"hess_{tagbase}", "E0": float(E0),
           "curv_timemix": curv, "n_neg": int(sum(c < 0 for c in curv)),
           "min_curv": float(np.min(curv)),
           "curv_spatial_ctrl": ctrl,
           "n_neg_ctrl": int(sum(c < 0 for c in ctrl))}
    with open(os.path.join(DATA, f"m5_21_3_row_hess_{tagbase}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ("tag", "n_neg", "min_curv",
                                          "n_neg_ctrl")}))
    return out


def p2(tagbase):
    """the generator catalog: kin sign + the twist channel per a0."""
    M, cfg = load_p1(tagbase)
    a0s = gen_catalog(cfg, M)
    rows = {}
    for nm, a0 in a0s.items():
        k = kin_of(M, a0, cfg)
        tw = twist_read(M, a0, cfg)
        rows[nm] = {"kin": float(k), "twist": tw}
        print(f"  {tagbase} {nm:12s} kin {k:+.6e} "
              f"k* {tw['k_star']:+.4f} dE(k*) {tw['dE_at_kstar']:+.3e}",
              flush=True)
    out = {"tag": f"p2_{tagbase}", "rows": rows,
           "min_kin": float(min(r["kin"] for r in rows.values())),
           "argmin": min(rows, key=lambda nm: rows[nm]["kin"])}
    with open(os.path.join(DATA, f"m5_21_3_row_p2_{tagbase}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"tag": out["tag"], "min_kin": out["min_kin"],
                      "argmin": out["argmin"]}))
    return out


def p3(tagbase, gen_name, omegas=(0.05, 0.1, 0.2, 0.4, 0.8),
       maxit=3000):
    """omega-ladder with profile re-relaxation: E*(omega) =
    min_M [E_static(M) + omega^2 kin(M; a0)] (a0 re-derived from the
    warm-start endpoint per rung, then FROZEN within the rung)."""
    M, cfg = load_p1(tagbase)
    free = ~pin_shell(cfg["n"], cfg["h"])
    E_stat = float(e_total(M, cfg))
    ladder = [{"omega": 0.0, "E": E_stat, "kin": None, "stop": "ref"}]
    dive_floor = E_stat - 10.0 * abs(E_stat) - 100.0
    for om in omegas:
        a0 = gen_catalog(cfg, M)[gen_name]
        k0 = float(kin_of(M, a0, cfg))
        tg = f"p3_{tagbase}_{gen_name}_w{om:g}"
        M, info = fire(M, cfg, free, maxit, a0=a0, omega=om, tag=tg,
                       dive_floor=dive_floor)
        eu, ev = e_parts(M, cfg)
        kin_end = float(kin_of(M, a0, cfg))
        E = float(eu + ev + om ** 2 * kin_end)
        ladder.append({"omega": om, "E": E, "E_u": float(eu),
                       "E_v": float(ev), "kin": kin_end,
                       "kin_pre": k0, "stop": info["stop"]})
        print(json.dumps(ladder[-1]), flush=True)
        if info["stop"] in ("dive", "non-finite"):
            break
    out = {"tag": f"p3_{tagbase}_{gen_name}", "E_static": E_stat,
           "ladder": ladder}
    with open(os.path.join(DATA,
                           f"m5_21_3_row_p3_{tagbase}_{gen_name}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_3_p3_{tagbase}_{gen_name}.npz"),
        M=M.astype(np.float32), s=cfg["s"], delta=cfg["delta"],
        h=cfg["h"])
    return out


def collect():
    import glob
    rows = {}
    for p in sorted(glob.glob(os.path.join(DATA, "m5_21_3_row_*.json"))):
        key = os.path.basename(p)[len("m5_21_3_row_"):-len(".json")]
        with open(p) as f:
            rows[key] = json.load(f)
    with open(os.path.join(DATA, "m5_21_3_all.json"), "w") as f:
        json.dump(rows, f, indent=1)
    print(f"collected {len(rows)} rows")


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0] if ARGV else "status"
    if mode == "gates":
        gates()
    elif mode == "p1":
        p1(ARGV[1], *(int(a) for a in ARGV[2:3]))
    elif mode == "hess":
        hess(ARGV[1])
    elif mode == "p2":
        p2(ARGV[1])
    elif mode == "p3":
        om = tuple(float(x) for x in ARGV[3].split(",")) if \
            len(ARGV) > 3 else (0.05, 0.1, 0.2, 0.4, 0.8)
        p3(ARGV[1], ARGV[2], om)
    elif mode == "collect":
        collect()
    else:
        collect()
