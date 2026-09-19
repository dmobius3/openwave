"""M5.32 R22-2: the Coulomb certificate pair (the pair gate of the author's
2026-09-15 post on the coordination thread, section II).

EQUATIONS FIRST
---------------
Field: a unit director n(x) on the degenerate-pair vacuum,
    M_sp(n) = delta I + (1 - delta) n n^T,  M_00 constant,  M_0i = 0,
energy = the certified static quartic of the stack on M(n) (the potential
vanishes identically, the eigenvalues never leave the vacuum):
    E[n] = 4 h^3 sum_br wt sum_cells sum_{i<j} <F_ij, F_ij>,  F_ij = [A_i, A_j],
    A_i = d_i M_sp with the stack's forward / backward branches (wt 1/2 each).
By R22-0 (a) the continuum density is k |E_top|^2 with k = 8 (1 - delta)^4 in
the stack normalization (the author's 4 (1 - delta)^4, c4 = 2), so
    single core outside r_c:   4 pi k / r_c      = 32 pi (1 - delta)^4 / r_c
    pair cross term, free space: -/+ 8 pi k / d  = -/+ 64 pi (1 - delta)^4 / d
(48.27 / d at delta 0.3; the plan post's 24.14 / d is the same number in the
author's normalization).

The box. Every cell outside the pinned cores relaxes, the faces included. The
natural boundary condition of int |E_top|^2 is E_top normal to the wall
(the boundary term of the variation is a . (E x nu) with a_i = n . (dn x d_i n)),
the grounded-conductor condition. The box-aware prediction is therefore the
Dirichlet Green function of the cube, not 1 / d:
    U_box(d) = -/+ 8 pi k G_D(x_1, x_2),   G_D = 1 / d + psi_1(x_2),
    psi_1 harmonic in the cube, psi_1 = -1 / |x - x_1| on the walls
(solved on the same cell-centered lattice with ghost-cell walls by a DST-II).
Both predictions are reported; the free-space slope is the claim under test,
the box form says how much of any departure is the wall. Audit note (C11.1):
the wall condition holds (E x nu = 0, no multiplier from the fixed flux), but
the BULK Euler-Lagrange equation is only (curl E) x E = 0, weaker than
curl E = 0. So 8 pi k G_D is the lower bound over divergence-free fields,
attained only if the box Coulomb field is realizable as a director texture:
a bound-type prediction, which is exactly the author's claim under test.

Objects (centers on the x axis at +/- d / 2, the stereographic pole on z, so
perpendicular to the pole; no quotient-map seed):
    like    w = w_1 w_2            total degree 2 (the author's product ansatz)
    unlike  w = w_1 conj(w_2)      total degree 0 (hedgehog times a reflected one)
    ref     w = w_1                one core at + d / 2 in the same box
    w_k = (x_k + i y) / (r_k + z), n = (2 Re w, 2 Im w, 1 - |w|^2) / (1 + |w|^2).
Cores: cells with r_k < R_C are PINNED to the exact single-hedgehog texture,
ball 1 to r-hat_1, ball 2 to diag(-1, -1, 1) r-hat_2 (like) or
diag(-1, 1, 1) r-hat_2 (unlike), the constant O(3) images the seeds approach
at the centers. The energy is invariant under a constant O(3) on n, so the
pinned balls carry exactly the reference's self-energy, and
    U(d) = E_pair(d) - 2 E_ref(d)
(the core at - d / 2 is the INVERSION image of the reference, x -> -x, an exact
symmetry of the stack's functional; a single-axis mirror is not, see R22-1).
Pin artifact, stated in advance: a ball of fixed uniform flux excludes the
other core's field E_0 = 1 / d^2. Outside, the induced dipole p = -E_0 r_c^3 / 2
adds k (2 pi / 3) r_c^3 E_0^2; inside, the pinned texture never carries E_0 and
k (4 pi / 3) r_c^3 E_0^2 is missing; the cross terms vanish by the angular
average. Net, per pair, a CHARGE-EVEN  -(4 pi / 3) k r_c^3 / d^4
(-0.097 at d 6, -0.031 at d 8, -0.006 at d 12: 1.2 to 0.15 percent of U).
Audit note (C12.2): that is the l = 1 term; the full multipole sum
net_l = -(l / (l + 1)) (4 pi / (2 l + 1)) k r_c^(2l+1) / d^(2l+2) gives
-0.1127 / -0.0333 / -0.0140 / -0.00628 at d 6 / 8 / 10 / 12 (collect uses it).

Descent: scipy L-BFGS on the unconstrained m, n = m / |m| (the radial
direction is flat and harmless), free cells only; the residual is the max
over free cells of the tangential gradient per cell volume; stop at
RES_GATE or MAX_ITER, the attained residual is recorded per row.

Reads: E, U, the flux of the central-difference E_top through a cube around
each core and through the outer box (degree = flux / 4 pi), the energy
of E_top - E_pred in the free region (the transverse energy, sign of E_pred
chosen by the overlap), the central-stencil energy of the end field (the
second discretization, R22-0 (a) B).

Labels (pre-registered in the plan post):
    COULOMB_CERTIFIED  U(+-) and U(++) follow the prediction, U_grav inside
                       the stated residuals
    LIKE_EXCESS        the like pair keeps transverse energy
    UNLIKE_TUBE        a flip tube sits between the unlike cores
    PAIR_BOX_LIMITED   the two boxes disagree

Modes: smoke | run [workers] | collect
Regenerate: python m5_32_r22_2_pair.py run 10   (detached; see the task doc)
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_2_pair.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r22_2")
G, DELTA = 8.0, 0.3
R_C = 2.5
RES_GATE = 1e-7
MAX_ITER = 4000
CKPT_EVERY = 200
MEM_GB = {32: 0.5, 48: 1.5, 64: 3.5}
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def ups_guard():
    try:
        sys.path.insert(0, os.path.expanduser("~/Library/Application Support/UPSSentinel"))
        import ups_guard as ug

        ug.register()
        return ug
    except Exception:  # noqa: BLE001
        return None


# ================= the director energy =================
def m3_of(n, delta):
    return delta * np.eye(3) + (1.0 - delta) * n[..., :, None] * n[..., None, :]


def energy_grad_n(m, h, delta, need_grad=True, stencil="sym"):
    """E and dE/dm for n = m / |m| on the stack's stencil. m: (n, n, n, 3)."""
    nm = np.linalg.norm(m, axis=-1, keepdims=True)
    n = m / nm
    M3 = m3_of(n, delta)
    e = 0.0
    GM = np.zeros_like(M3) if need_grad else None
    for br, wt in B3.branches(stencil):
        A = [B3.d1(M3, ax, h, br) for ax in range(3)]
        GA = [np.zeros_like(M3) for _ in range(3)] if need_grad else None
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ A[j] - A[j] @ A[i]
                e += wt * 4.0 * float(np.sum(F * F))
                if need_grad:
                    GA[i] += wt * 8.0 * (F @ A[j] - A[j] @ F)
                    GA[j] += wt * 8.0 * (A[i] @ F - F @ A[i])
        if need_grad:
            for ax in range(3):
                GM += B3.d1_adj(GA[ax], ax, h, br)
    E = h**3 * e
    if not need_grad:
        return E, None
    GM = 0.5 * (GM + GM.swapaxes(-1, -2)) * h**3
    gn = 2.0 * (1.0 - delta) * np.einsum("...ab,...b->...a", GM, n)
    gt = gn - n * np.sum(gn * n, axis=-1, keepdims=True)
    return E, gt / nm


def density_n(n, h, delta, stencil="sym"):
    M3 = m3_of(n, delta)
    e = np.zeros(n.shape[:3])
    for br, wt in B3.branches(stencil):
        A = [B3.d1(M3, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ A[j] - A[j] @ A[i]
                e += wt * 4.0 * np.einsum("...ab,...ab->...", F, F)
    return e


# ================= seeds and pins =================
def n_of_w(w):
    a = np.abs(w) ** 2
    return np.stack([2 * w.real, 2 * w.imag, 1 - a], axis=-1) / (1 + a)[..., None]


def w_single(X, Y, Z, cx):
    x = X - cx
    r = np.sqrt(x * x + Y * Y + Z * Z)
    return (x + 1j * Y) / np.maximum(r + Z, 1e-12)


def build(kind, d, n_, L):
    h = L / n_
    X, Y, Z = B3.coords(n_, h)
    c1, c2 = +d / 2.0, -d / 2.0
    r1 = np.sqrt((X - c1) ** 2 + Y * Y + Z * Z)
    r2 = np.sqrt((X - c2) ** 2 + Y * Y + Z * Z)
    rh1 = np.stack([X - c1, Y, Z], axis=-1) / np.maximum(r1, 1e-300)[..., None]
    rh2 = np.stack([X - c2, Y, Z], axis=-1) / np.maximum(r2, 1e-300)[..., None]
    w1, w2 = w_single(X, Y, Z, c1), w_single(X, Y, Z, c2)
    pin = r1 < R_C
    if kind == "ref":
        n = rh1.copy()
    else:
        n = n_of_w(w1 * w2) if kind == "like" else n_of_w(w1 * np.conj(w2))
        img = np.array([-1.0, -1.0, 1.0]) if kind == "like" else np.array([-1.0, 1.0, 1.0])
        b2 = r2 < R_C
        n[b2] = rh2[b2] * img
        pin = pin | b2
    n[r1 < R_C] = rh1[r1 < R_C]
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return n, pin, (X, Y, Z, r1, r2, h)


# ================= the box-aware prediction =================
def psi_dirichlet(n_, L, x1):
    """psi harmonic in the cube [-L/2, L/2]^3 with psi = -1/|x - x1| on the walls; cell-centered, ghost-cell walls."""
    from scipy.fft import dstn, idstn

    h = L / n_
    c = (np.arange(n_) - (n_ - 1) / 2.0) * h
    rhs = np.zeros((n_, n_, n_))
    for ax in range(3):
        for side, idx in ((-L / 2.0, 0), (L / 2.0, n_ - 1)):
            P = np.meshgrid(c, c, indexing="ij")
            pts = [None, None, None]
            o = [a for a in range(3) if a != ax]
            pts[o[0]], pts[o[1]] = P[0], P[1]
            pts[ax] = np.full_like(P[0], side)
            g = -1.0 / np.sqrt(
                (pts[0] - x1[0]) ** 2 + (pts[1] - x1[1]) ** 2 + (pts[2] - x1[2]) ** 2
            )
            sl = [slice(None)] * 3
            sl[ax] = idx
            rhs[tuple(sl)] += -2.0 * g / h**2
    k = np.arange(1, n_ + 1)
    lam = -(2.0 - 2.0 * np.cos(np.pi * k / n_)) / h**2
    den = lam[:, None, None] + lam[None, :, None] + lam[None, None, :]
    return idstn(dstn(rhs, type=2) / den, type=2), c


def interp3(F, c, x):
    h = c[1] - c[0]
    out = F
    for ax in range(3):
        t = (x[ax] - c[0]) / h
        i0 = int(np.clip(np.floor(t), 0, len(c) - 2))
        f = t - i0
        out = np.take(out, i0, axis=0) * (1 - f) + np.take(out, i0 + 1, axis=0) * f
    return float(out)


def green_box(d, n_, L):
    x1 = np.array([d / 2.0, 0.0, 0.0])
    x2 = np.array([-d / 2.0, 0.0, 0.0])
    psi, c = psi_dirichlet(n_, L, x1)
    return 1.0 / d + interp3(psi, c, x2), psi, c


# ================= reads =================
def etop_central(n, h):
    dn = [np.gradient(n, h, axis=ax) for ax in range(3)]
    cr = lambda a, b: np.einsum("...a,...a->...", n, np.cross(dn[a], dn[b]))  # noqa: E731
    return np.stack([cr(1, 2), cr(2, 0), cr(0, 1)], axis=-1)


def cube_flux(E, c, h, center, half):
    """flux of E through the cube |x - center|_inf < half (cell-centered faces, midpoint rule)."""
    idx = [np.where(np.abs(c - center[a]) < half)[0] for a in range(3)]
    tot = 0.0
    for ax in range(3):
        o = [a for a in range(3) if a != ax]
        for sgn, i_face in ((+1, idx[ax][-1]), (-1, idx[ax][0])):
            sl = [None] * 3
            sl[ax] = i_face
            sl[o[0]], sl[o[1]] = idx[o[0]][:, None], idx[o[1]][None, :]
            # the face sits between cell i_face and its outer neighbor: average the two
            sl2 = list(sl)
            sl2[ax] = i_face + sgn
            val = 0.5 * (E[tuple(sl) + (ax,)] + E[tuple(sl2) + (ax,)])
            tot += sgn * float(np.sum(val)) * h * h
    return tot


def reads(n, pin, geo, kind, d, n_, L):
    X, Y, Z, r1, r2, h = geo
    k = 8.0 * (1.0 - DELTA) ** 4
    out = {
        "E": energy_grad_n(n, h, DELTA, need_grad=False)[0],
        "E_central": energy_grad_n(n, h, DELTA, need_grad=False, stencil="cen")[0],
        "E_pinned_cells": float(density_n(n, h, DELTA)[pin].sum() * h**3),
    }
    c = (np.arange(n_) - (n_ - 1) / 2.0) * h
    E = etop_central(n, h)
    half = R_C + 2.0 * h
    out["degree_core1"] = cube_flux(E, c, h, (d / 2.0, 0, 0), half) / (4 * np.pi)
    if kind != "ref":
        out["degree_core2"] = cube_flux(E, c, h, (-d / 2.0, 0, 0), half) / (4 * np.pi)
    out["degree_box"] = cube_flux(E, c, h, (0, 0, 0), L / 2.0 - 2.5 * h) / (4 * np.pi)
    # the transverse energy against the box Coulomb field
    q2 = {"ref": 0.0, "like": 1.0, "unlike": -1.0}[kind]
    Ep = np.zeros_like(E)
    for q, cx in ((1.0, d / 2.0), (q2, -d / 2.0)):
        if q == 0.0:
            continue
        psi, _ = psi_dirichlet(n_, L, np.array([cx, 0.0, 0.0]))
        rr = np.sqrt((X - cx) ** 2 + Y * Y + Z * Z)
        free = np.stack([X - cx, Y, Z], axis=-1) / np.maximum(rr, 1e-300)[..., None] ** 3
        Ep += q * (free - np.stack(np.gradient(psi, h), axis=-1))
    m = (~pin) & (r1 > R_C + 2 * h) & ((r2 > R_C + 2 * h) | (kind == "ref"))
    inner = np.zeros_like(m)
    inner[2:-2, 2:-2, 2:-2] = True
    m &= inner
    best = None
    for sg in (+1.0, -1.0):
        dE = sg * E - Ep
        val = k * float(np.sum(dE[m] ** 2)) * h**3
        best = val if best is None else min(best, val)
    out["E_transverse_vs_box_coulomb"] = best
    out["E_coulomb_pred_same_mask"] = k * float(np.sum(Ep[m] ** 2)) * h**3
    out["E_top_energy_same_mask"] = k * float(np.sum(E[m] ** 2)) * h**3
    return out


# ================= the job =================
def job_tag(j):
    return f"{j['kind']}_d{j['d']:g}_n{j['n']}_L{j['L']:g}"


def run_job(j):
    from scipy.optimize import minimize

    t0 = time.time()
    tag = job_tag(j)
    ug = ups_guard()
    n0, pin, geo = build(j["kind"], j["d"], j["n"], j["L"])
    h = geo[-1]
    free = ~pin
    os.makedirs(OUT_NPZ, exist_ok=True)
    ck = os.path.join(OUT_NPZ, tag + "_ckpt.npz")
    m = n0.copy()
    it0 = 0
    if os.path.exists(ck):
        Zc = np.load(ck)
        m, it0 = Zc["m"], int(Zc["it"])
    row = dict(j)
    row.update(
        {
            "tag": tag,
            "h": h,
            "r_c": R_C,
            "pinned_cells": int(pin.sum()),
            "seed_reads": reads(n0, pin, geo, j["kind"], j["d"], j["n"], j["L"]),
            "resumed_at_iter": it0,
        }
    )
    st = {"it": it0, "trace": [], "stop": None}
    cache = {}

    def fun(x):
        mm = m.copy()
        mm[free] = x.reshape(-1, 3)
        E, g = energy_grad_n(mm, h, DELTA)
        res = float(np.max(np.abs(g[free]))) / h**3
        cache["last"] = (E, res)
        return E, g[free].ravel()

    class _Stop(Exception):
        pass

    def cb(xk):
        st["it"] += 1
        E, res = cache["last"]
        if st["it"] % 50 == 0 or st["it"] == it0 + 1:
            st["trace"].append([st["it"], E, res])
            log(f"{tag} it {st['it']:5d} E {E:.8f} res {res:.3e}")
        if st["it"] % CKPT_EVERY == 0:
            mm = m.copy()
            mm[free] = xk.reshape(-1, 3)
            np.savez_compressed(ck, m=mm, it=st["it"])
        if ug is not None and st["it"] % 25 == 0 and ug.wrap_up():
            st["stop"] = "UPS wrap-up"
            mm = m.copy()
            mm[free] = xk.reshape(-1, 3)
            np.savez_compressed(ck, m=mm, it=st["it"])
            raise _Stop()

    x = m[free].ravel().copy()
    try:
        left = max(1, j.get("max_iter", MAX_ITER) - it0)
        r = minimize(
            fun,
            x,
            jac=True,
            method="L-BFGS-B",
            callback=cb,
            options={
                "maxiter": left,
                "maxfun": 3 * left,
                "maxcor": 20,
                "ftol": 0.0,
                "gtol": RES_GATE * h**3,
            },
        )
        x = r.x
        st["stop"] = st["stop"] or str(r.message)
    except _Stop:
        pass
    m[free] = x.reshape(-1, 3)
    n_end = m / np.linalg.norm(m, axis=-1, keepdims=True)
    E_end, g_end = energy_grad_n(n_end, h, DELTA)
    row.update(
        {
            "iters": st["it"],
            "stop": st["stop"],
            "trace": st["trace"],
            "residual_end": float(np.max(np.abs(g_end[free]))) / h**3,
            "end_reads": reads(n_end, pin, geo, j["kind"], j["d"], j["n"], j["L"]),
        }
    )
    row["E"] = row["end_reads"]["E"]
    row["label"] = (
        "AT_GATE" if row["residual_end"] < RES_GATE else f"FALLING (res {row['residual_end']:.2e})"
    )
    np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), n=n_end)
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} E {row['E']:.6f} res {row['residual_end']:.2e} it {row['iters']} wall {row['wall_s']}"
    )
    return row


def jobs_all():
    J = []
    for n_, L in ((48, 48.0), (32, 48.0), (64, 64.0)):
        for d in (6.0, 8.0, 10.0, 12.0):
            if n_ == 64 and d not in (8.0, 12.0):
                continue
            for kind in ("ref", "unlike", "like"):
                J.append({"kind": kind, "d": d, "n": n_, "L": L})
    J.sort(key=lambda j: -j["n"])
    return J


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R22-2", "rows": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers, mem_budget=14.0):
    rows = load_json()["rows"]
    pending = [j for j in jobs if job_tag(j) not in rows]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    running = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        while pending or running:
            used = sum(MEM_GB[j["n"]] for j in running.values())
            for j in list(pending):
                if len(running) >= workers:
                    break
                if used + MEM_GB[j["n"]] <= mem_budget:
                    running[ex.submit(run_job, j)] = j
                    pending.remove(j)
                    used += MEM_GB[j["n"]]
            done = next(as_completed(list(running)))
            row = done.result()
            running.pop(done)
            J = load_json()
            J["rows"][row["tag"]] = row
            save_json(J)
    log("pool done")


# ================= smoke =================
def smoke():
    R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
    out = {}
    n_, L = 24, 36.0
    n, pin, geo = build("like", 8.0, n_, L)
    h = geo[-1]
    rng = np.random.default_rng(3)
    m = n * (1.0 + 0.1 * rng.normal(size=n.shape[:3] + (1,))) + 0.05 * rng.normal(size=n.shape)
    # 1 the energy against the production stack
    nn = m / np.linalg.norm(m, axis=-1, keepdims=True)
    cfg = B3.base_cfg(s=-1.0, g=G, n=n_, L=L, delta=DELTA)
    Mfull = B3.embed34(m3_of(nn, DELTA), cfg)
    pot = ("v4", R20.R0.roots_of(cfg, degenerate=True), R20.W1)
    Ep, Gp, _ = R20.energy_grad(Mfull, cfg, R20.params_of(G), pot)
    E, g = energy_grad_n(m, h, DELTA)
    out["E_rel_vs_production"] = abs(E - Ep) / abs(Ep)
    # 2 the gradient against the production gradient through the chain rule
    nm = np.linalg.norm(m, axis=-1, keepdims=True)
    gn = 2.0 * (1.0 - DELTA) * np.einsum("...ab,...b->...a", B3.sym4(Gp)[..., 1:, 1:], nn)
    gp = (gn - nn * np.sum(gn * nn, axis=-1, keepdims=True)) / nm
    out["grad_rel_vs_production"] = float(np.max(np.abs(g - gp)) / np.max(np.abs(gp)))
    # 3 a directional finite difference
    v = rng.normal(size=m.shape)
    t = 1e-5
    fd = (
        energy_grad_n(m + t * v, h, DELTA, False)[0] - energy_grad_n(m - t * v, h, DELTA, False)[0]
    ) / (2 * t)
    out["grad_fd_rel"] = abs(fd - float(np.sum(g * v))) / abs(fd)
    # 4 the Green function: the free-space limit and the sphere estimate
    Gd, _, _ = green_box(8.0, 48, 48.0)
    Gd2, _, _ = green_box(8.0, 64, 64.0)
    out["G_D_d8_L48"], out["G_D_d8_L64"], out["free_space_d8"] = Gd, Gd2, 1.0 / 8.0
    # 5 degrees of the seeds
    for kind in ("ref", "like", "unlike"):
        n, pin, geo = build(kind, 8.0, 48, 48.0)
        rd = reads(n, pin, geo, kind, 8.0, 48, 48.0)
        out[f"seed_{kind}"] = {
            k2: rd[k2] for k2 in rd if k2.startswith("degree") or k2 in ("E", "E_central")
        }
    ok = (
        out["E_rel_vs_production"] < 1e-10
        and out["grad_rel_vs_production"] < 1e-8
        and out["grad_fd_rel"] < 1e-5
    )
    out["PASS"] = bool(ok)
    print(json.dumps(out, indent=1))
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)


def collect():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    J = load_json()
    rows = J["rows"]
    k = 8.0 * (1.0 - DELTA) ** 4
    out = {"k": k, "slope_free_space": 8 * np.pi * k, "boxes": {}}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for (n_, L), col in (
        ((48, 48.0), "C0"),
        ((32, 48.0), "C1"),
        ((64, 64.0), "C2"),
        ((64, 48.0), "C3"),
    ):
        tab = []
        for d in (6.0, 8.0, 10.0, 12.0):
            t = {
                kind: rows.get(job_tag({"kind": kind, "d": d, "n": n_, "L": L}))
                for kind in ("ref", "like", "unlike")
            }
            if any(v is None for v in t.values()):
                continue
            Gd = green_box(d, n_, L)[0]
            U = {kd: t[kd]["E"] - 2 * t["ref"]["E"] for kd in ("like", "unlike")}
            Uc = {
                kd: t[kd]["end_reads"]["E_central"] - 2 * t["ref"]["end_reads"]["E_central"]
                for kd in ("like", "unlike")
            }
            tab.append(
                {
                    "d": d,
                    "G_D": Gd,
                    "U_like": U["like"],
                    "U_unlike": U["unlike"],
                    "U_grav": 0.5 * (U["like"] + U["unlike"]),
                    "U_odd": 0.5 * (U["like"] - U["unlike"]),
                    "U_box_pred": 8 * np.pi * k * Gd,
                    "U_free_pred": 8 * np.pi * k / d,
                    "pin_estimate_even": -2.0
                    * k
                    * sum(
                        (l / (l + 1.0))
                        * (4 * np.pi / (2 * l + 1))
                        * R_C ** (2 * l + 1)
                        / d ** (2 * l + 2)
                        for l in range(1, 40)
                    ),
                    "lattice_scale_h_over_d2": 8 * np.pi * k * (L / n_) / d**2,
                    "U_like_central": Uc["like"],
                    "U_unlike_central": Uc["unlike"],
                    "residual": {kd: t[kd]["residual_end"] for kd in t},
                    "iters": {kd: t[kd]["iters"] for kd in t},
                    "degrees": {
                        kd: [
                            t[kd]["end_reads"].get("degree_core1"),
                            t[kd]["end_reads"].get("degree_core2"),
                            t[kd]["end_reads"]["degree_box"],
                        ]
                        for kd in t
                    },
                    "E_transverse": {
                        kd: t[kd]["end_reads"]["E_transverse_vs_box_coulomb"] for kd in t
                    },
                    "E_coulomb_same_mask": {
                        kd: t[kd]["end_reads"]["E_coulomb_pred_same_mask"] for kd in t
                    },
                    "seed_U": {
                        kd: t[kd]["seed_reads"]["E"] - 2 * t["ref"]["seed_reads"]["E"]
                        for kd in ("like", "unlike")
                    },
                }
            )
        if not tab:
            continue
        box = {"rows": tab}
        if len(tab) >= 2:
            gd = np.array([r["G_D"] for r in tab])
            for key in ("U_like", "U_unlike", "U_odd"):
                u = np.array([r[key] for r in tab])
                sl, ic = np.polyfit(gd, u, 1)
                sl1, ic1 = np.polyfit(1.0 / np.array([r["d"] for r in tab]), u, 1)
                box[f"fit_{key}"] = {
                    "slope_vs_G_D": float(sl),
                    "intercept_vs_G_D": float(ic),
                    "slope_vs_inv_d": float(sl1),
                    "intercept_vs_inv_d": float(ic1),
                    "slope_over_8pik_vs_G_D": float(sl / (8 * np.pi * k)),
                    "slope_over_8pik_vs_inv_d": float(sl1 / (8 * np.pi * k)),
                }
        out["boxes"][f"n{n_}_L{L:g}"] = box
        ds = [r["d"] for r in tab]
        axs[0].plot(ds, [r["U_like"] for r in tab], "o-", color=col, label=f"U(++) n{n_} L{L:g}")
        axs[0].plot(ds, [r["U_unlike"] for r in tab], "s-", color=col, label=f"U(+-) n{n_} L{L:g}")
        axs[0].plot(ds, [r["U_box_pred"] for r in tab], ":", color=col)
        axs[0].plot(ds, [-r["U_box_pred"] for r in tab], ":", color=col)
        axs[1].plot(ds, [r["U_grav"] for r in tab], "o-", color=col, label=f"U_grav n{n_} L{L:g}")
        axs[1].plot(ds, [r["pin_estimate_even"] for r in tab], "--", color=col, alpha=0.5)
    dd = np.linspace(6, 12, 50)
    axs[0].plot(dd, 8 * np.pi * k / dd, "k--", lw=0.8, label="free space 8 pi k / d")
    axs[0].plot(dd, -8 * np.pi * k / dd, "k--", lw=0.8)
    axs[0].set_xlabel("d"), axs[0].set_ylabel("U(d)"), axs[0].legend(fontsize=6)
    axs[0].set_title("pair energy (dotted: the grounded-box Green function)")
    axs[1].axhline(0, color="k", lw=0.5)
    axs[1].set_xlabel("d"), axs[1].set_title(
        "charge-even half-sum (dashed: the pin estimate)"
    ), axs[1].legend(fontsize=6)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "..", "plots"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "..", "plots", "m5_32_r22_2_pair.png"), dpi=130)
    with open(OUT_JSON.replace(".json", "_collect.json"), "w") as f:
        json.dump(out, f, indent=1)
    for name, box in out["boxes"].items():
        print(name)
        for r in box["rows"]:
            print(
                f"  d {r['d']:4.0f}  U++ {r['U_like']:+8.4f}  U+- {r['U_unlike']:+8.4f}  box {r['U_box_pred']:7.4f}  free {r['U_free_pred']:7.4f}"
                f"  U_grav {r['U_grav']:+8.4f}  pin {r['pin_estimate_even']:+7.4f}  h/d2 {r['lattice_scale_h_over_d2']:6.3f}"
                f"  res {max(r['residual'].values()):.1e}  ET like {r['E_transverse']['like']:.3f} unlike {r['E_transverse']['unlike']:.3f}"
            )
        for key in ("fit_U_like", "fit_U_unlike", "fit_U_odd"):
            if key in box:
                print(
                    f"  {key}: slope/8pik vs G_D {box[key]['slope_over_8pik_vs_G_D']:+.4f} (intercept {box[key]['intercept_vs_G_D']:+.4f}); "
                    f"vs 1/d {box[key]['slope_over_8pik_vs_inv_d']:+.4f} (intercept {box[key]['intercept_vs_inv_d']:+.4f})"
                )


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    elif mode == "refine":
        # added at EXECUTE (2026-09-19 12:55 EDT): the central-stencil read of the like pair sits 8 to 25 percent
        # under the stack's, the unlike pair's agrees; the d 8 triple at h 0.75 in the SAME box gives the h sequence
        run_pool(
            [{"kind": kd, "d": 8.0, "n": 64, "L": 48.0} for kd in ("ref", "unlike", "like")], 3
        )
    elif mode == "collect":
        collect()
    elif mode == "one":
        print(
            json.dumps(
                {k: v for k, v in run_job(json.loads(sys.argv[2])).items() if k != "trace"},
                indent=1,
            )
        )


if __name__ == "__main__":
    main()
