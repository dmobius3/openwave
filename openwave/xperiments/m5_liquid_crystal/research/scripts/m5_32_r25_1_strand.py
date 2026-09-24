"""M5.32 R25-1: the straight transverse strand on the certified biaxial vacuum,
the delta and w ladders against the Bogomolny tension (the author's Packet A,
amended on the 2026-09-22 19:24 UTC reply).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta; the code
branch s = -1: M_vac = diag(8, 1, delta, 0). The certified action per cell
    E = 4 sum_{i<j} <F_ij, F_ij>_eta + w sum_{p<=4} (tr N^p - C_p)^2,
    F_ij = [d_i M, d_j M]_eta,  C_p = sum of the vacuum eigenvalues' p-th powers,
w = W1 x {6.25, 25, 100} (W1 the certified weight). In the BLOCK SECTOR only
the transverse pair moves: with the pair block S = s_0 I + b(rho) (cos(m phi)
sigma_z + sin(m phi) sigma_x), s_0 = b_0 = delta / 2 (the pair (delta, 0)),
    u = 32 m^2 (b b' / rho)^2  (the winding costs nothing, the core everything),
    V4 = w (f - f_0)^2 [4 + 36 s_0^2 + (12 s_0^2 + 2 (f + f_0))^2],  f = b^2,
so the tension of a straight line obeys the Bogomolny bound (R25-0 a)
    T >= T_BPS = pi sqrt(32 w K) b_0^4 m,  K = 4 + 36 s_0^2 + 144 s_0^4,
attained at m = 1 by f = f_0 (1 - exp(-kappa rho^2)), kappa = sqrt(w K / 32),
the core radius (f = f_0 / 2) growing as w^(-1/4). T / delta^4 = 0.330, 0.302,
0.299, 0.299 at W1 x 25 for delta 0.3, 0.1, 0.03, 0.01 (the K(s_0) correction).

THE MEASUREMENT
---------------
A z-invariant slab n x n x NZ at spacing h, the x and y faces pinned (depth
1.6) at the BPS profile, z free (the stack's derivative is one-sided at the
faces). Seed = the BPS profile (the exact block minimizer, an admissible
configuration no block-sector state beats) plus a z-DEPENDENT random kick of
KICK x delta in the 6 spatial entries of the free cells (M_00 is stiff under
V4 and slaved by it; M_0i stays 0, the static sector, see static_sector): the
channels that leave the block (the director row, M_33, the z variation) are
all excited. RUN-TIME DEVIATION 2026-09-23: the first pool kicked M_0i too and
every row ran away (the fullkick witness files).
FIRE (FIRE_ITERS) then L-BFGS on the 7 static entries to the gate
    fmax < max(GATE_REL x T_BPS / h, GATE_FLOOR)
(the delta 0.01 rows sit at the round-off floor and say so). Reads: T =
(E_u + V4) / (NZ h), T / T_BPS, T / T_1D (the 1D pinned-axis minimizer with the
exact block potential at the same delta, w), the block departure per channel
relative to delta, the z variation, the core radius (interpolated), the L-BFGS
drop after FIRE (a FIRE stall if over 1e-4 relative). A row above the bound at
the gate is re-seeded once unkicked; still above by 2 percent, it is reported
as an instrument fault, never a label.

Ladder (22 rows): delta {0.3, 0.1, 0.03, 0.01} x w1s {6.25, 25, 100} at n 48
L 48 (h 1); plus h {1.5, 1.0, 0.75, 0.5} at L 48 (n 32, 48, 64, 96) and
L {48, 72, 96} at h 1 (n 48, 72, 96) for delta 0.3 and 0.03 at w1s 25.

Pre-registered labels (amended 2026-09-23 before any row is relaxed):
    STRAND_IN_BLOCK      every row returns to the block (departure under 1e-6 of
                         delta) with T within 2 percent of T_BPS at h 1 or finer
    STRAND_LEAVES_BLOCK  a row at the gate ends over 2 percent BELOW T_BPS with a
                         departure over 1e-3 of delta, at both h 1 and h 0.5 for
                         that (delta, w)
    UNRESOLVED           otherwise
The delta^4 and sqrt(w) laws are identities at the bound: REPORTED with the
K(s_0) correction, never labeled.

Modes: smoke | run [workers] | collect | plot | jobs. Output
data/m5_32_r25_1_strand.json, arrays data/m5_32_r25_1/ (local).
Regenerate: run 12 about 30 min (the n 96 rows the long pole); collect seconds.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r25_1_strand.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r25_1")
NZ = 4
SG = -8.0  # vac4 = diag(-SG, 1, delta, 0) = diag(8, 1, delta, 0)
PIN_DEPTH = 1.6
KICK = 0.02
KICK_EXCESS = 0.25  # the kicked seed's excess over the BPS seed, in units of T_bps
FIRE_ITERS = 2000
LBFGS_ITERS = 4000
LBFGS_CHUNK = 500
GATE_REL, GATE_FLOOR = 1e-6, 1e-13
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
R0 = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")
W1 = B3.W1
IU4 = np.triu_indices(4)
OFF4 = np.where(IU4[0] == IU4[1], 1.0, 2.0)


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the bound and the 1D minimizer =================
def K_lead(s0):
    return 4.0 + 36.0 * s0**2 + 144.0 * s0**4


def T_bps(delta, w, m=1):
    s0 = b0 = delta / 2.0
    return float(np.pi * np.sqrt(32.0 * w * K_lead(s0)) * b0**4 * m)


def kappa_bps(delta, w, m=1):
    return float(np.sqrt(w * K_lead(delta / 2.0) / (32.0 * m * m)))


def V_block(f, f0, s0, w):
    d = f - f0
    return w * ((2 * d) ** 2 + (6 * s0 * d) ** 2 + (d * (12 * s0**2 + 2 * (f + f0))) ** 2)


def T_1d_of(f, rho, f0, s0, w, m=1):
    """the 1D tension by the midpoint scheme (immune to the central-difference zigzag)."""
    fp = np.diff(f) / np.diff(rho)
    rm = 0.5 * (rho[1:] + rho[:-1])
    fm = 0.5 * (f[1:] + f[:-1])
    u = 8.0 * m * m * (fp / rm) ** 2
    return float(2 * np.pi * np.sum((u + V_block(fm, f0, s0, w)) * rm * np.diff(rho)))


def V_slaved(f, s0, w, C):
    """RUN-TIME READ (2026-09-23): the block potential minimized over (M_00, M_33) at fixed pair
    weight f = b^2: the bound freezes M_00 = 8 and M_33 = 1, the stack lets them relax per cell,
    and the radial gradients of the diagonal slots commute with the pair's winding, so the 1D
    kinetic term is unchanged and the Bogomolny argument goes through with V_slaved in place of
    the block V (T_slaved below). Returns (V_min, m00, m33)."""
    from scipy.optimize import minimize

    b = np.sqrt(max(f, 0.0))

    def V(x):
        ev = np.array([-x[0], s0 + b, s0 - b, x[1]])
        return w * sum((np.sum(ev**q) - C[q - 1]) ** 2 for q in (1, 2, 3, 4))

    r = minimize(
        V,
        [8.0, 1.0],
        method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-30, "maxiter": 4000},
    )
    return float(r.fun), float(r.x[0]), float(r.x[1])


def T_slaved(delta, w, m=1):
    """the Bogomolny value with M_00 and M_33 slaved: 4 pi sqrt(8) m int_0^f0 sqrt(V_slaved) df
    (the frozen version of this integral with the exact block V is the pinned 1D minimizer)."""
    from scipy.integrate import quad

    s0 = delta / 2.0
    f0 = s0**2
    C = [(-8.0) ** q + 1.0 + delta**q for q in (1, 2, 3, 4)]
    val = quad(lambda f: np.sqrt(V_slaved(f, s0, w, C)[0]), 0.0, f0, limit=100)[0]
    return 4.0 * np.pi * np.sqrt(8.0) * m * val


def T_1d_min(delta, w, m=1, rho_max=40.0, npts=801):
    """the 1D minimizer with the axis PINNED at f(0) = 0 (the continuum's core) and the exact potential."""
    from scipy.optimize import minimize

    s0 = b0 = delta / 2.0
    f0 = b0**2
    rc = np.linspace(1e-3, rho_max, npts)
    f_init = f0 * (1 - np.exp(-kappa_bps(delta, w, m) * rc**2))

    def obj(x):
        return T_1d_of(np.concatenate([[0.0], x, [f0]]), rc, f0, s0, w, m)

    r = minimize(obj, f_init[1:-1], method="L-BFGS-B", options={"maxiter": 3000, "ftol": 1e-14})
    return float(r.fun)


# ================= the slab =================
def slab_cfg(n, L, delta):
    return B3.base_cfg(s=-1.0, n=n, L=float(L), delta=delta)


def slab_coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    return X, Y, np.sqrt(X * X + Y * Y), np.arctan2(Y, X)


def bps_field(n, nz, h, delta, w, m=1, kappa=None):
    """the block-sector field with the BPS profile: M = diag(8, 1) + the pair block, z-invariant."""
    X, Y, rho, phi = slab_coords(n, h)
    k = kappa_bps(delta, w, m) if kappa is None else kappa
    s0 = b0 = delta / 2.0
    b = b0 * np.sqrt(1.0 - np.exp(-k * rho**2))
    M = np.zeros((n, n, nz, 4, 4))
    M[..., 0, 0] = -SG
    M[..., 3, 3] = 1.0  # the director eigenvalue
    M[..., 1, 1] = (s0 + b * np.cos(m * phi))[:, :, None]
    M[..., 2, 2] = (s0 - b * np.cos(m * phi))[:, :, None]
    M[..., 1, 2] = M[..., 2, 1] = (b * np.sin(m * phi))[:, :, None]
    return M


def free_mask_slab(n, nz, h):
    wc = max(1, int(np.ceil(PIN_DEPTH / h)))
    free = np.ones((n, n, nz), dtype=bool)
    free[:wc] = free[-wc:] = False
    free[:, :wc] = free[:, -wc:] = False
    return free


def density_u(M, cfg):
    """the curvature density per cell, 4 sum <F, F>_eta (the stack's e_parts without the sum)."""
    e = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e = e + wt * 4.0 * B3.inner_eta(F, F)
    return cfg["h"] ** 3 * e


def energy_parts(M, cfg, w):
    eu, _ = B3.e_parts(M, cfg)
    ev, _ = R0.v4_energy_grad(M, cfg, R0.roots_of(cfg), w, need_grad=False)
    return float(eu), float(ev)


def grad_w(M, cfg, w):
    """the stack's gradient with V4 at W1 replaced by V4 at w (the R20 hook, by subtraction)."""
    G = B3.grad(M, cfg)
    q = R0.roots_of(cfg)
    _, g1 = R0.v4_energy_grad(M, cfg, q, W1, True)
    _, gw = R0.v4_energy_grad(M, cfg, q, w, True)
    return static_sector(G - g1 + gw)


def static_sector(G):
    """RUN-TIME DEVIATION (2026-09-23 15:47 UTC): the descent stays in the block-diagonal
    static sector M_0i = 0, the sector every rung's static measurement lives in (R20 § 6.10:
    the static gradient is block-diagonal; R24 § 6.15: the time-space components lower the
    certified static action without bound and are not the Hamiltonian's sector). The first
    R25-1 pool kicked M_0i and every row ran away through a boost texture (E to -1e15, kept
    as the witness data/m5_32_r25_1_strand_fullkick.json). The gradient's time row is zeroed
    here; M_00 stays free (unkicked)."""
    G = G.copy()
    G[..., 0, 1:] = 0.0
    G[..., 1:, 0] = 0.0
    return G


def fire_slab(M0, cfg, w, free, iters, dt0=0.02, dt_max=0.2, log_every=500, tag=""):
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    fr = free[..., None, None].astype(float)
    F = -grad_w(M, cfg, w) * fr
    E0 = sum(energy_parts(M, cfg, w))
    for it in range(1, iters + 1):
        P = float(np.sum(F * v))
        if P > 0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * F * (vn / max(fn, 1e-300))
            if n_up > 5:
                dt, alpha = min(dt * 1.1, dt_max), alpha * 0.99
        else:
            v[:] = 0.0
            dt, alpha, n_up = dt * 0.5, 0.1, 0
        v = v + dt * F
        M = B3.sym4(M + dt * v * fr)
        F = -grad_w(M, cfg, w) * fr
        if it % log_every == 0 or it == iters:
            E = sum(energy_parts(M, cfg, w))
            log(f"{tag} FIRE it {it:5d} E {E:.9e} fmax {np.max(np.abs(F)):.2e} dt {dt:.3f}")
    E1 = sum(energy_parts(M, cfg, w))
    return M, {"E_start": E0, "E_end": E1, "fmax_end": float(np.max(np.abs(F))), "iters": iters}


def pack(M, free):
    return M[free][:, IU4[0], IU4[1]].ravel().copy()


def unpack(x, M_base, free):
    M = M_base.copy()
    blk = np.zeros((int(free.sum()), 4, 4))
    blk[:, IU4[0], IU4[1]] = x.reshape(-1, 10)
    blk = blk + blk.swapaxes(-1, -2) - np.einsum("...ii->...i", blk)[..., None] * np.eye(4)
    M[free] = blk
    return M


def lbfgs_slab(M0, cfg, w, free, gate, max_iter, chunk=LBFGS_CHUNK, tag=""):
    from scipy.optimize import minimize

    M = M0.copy()
    done, chunks, verdict = 0, [], "FALLING"
    E_start = sum(energy_parts(M, cfg, w))
    while done < max_iter:
        base = M.copy()

        def fun(x, base=base):
            Mx = unpack(x, base, free)
            eu, ev = energy_parts(Mx, cfg, w)
            G = grad_w(Mx, cfg, w)[free][:, IU4[0], IU4[1]] * OFF4
            return eu + ev, G.ravel()

        k = min(chunk, max_iter - done)
        res = minimize(
            fun,
            pack(M, free),
            jac=True,
            method="L-BFGS-B",
            options={"maxcor": 20, "maxiter": k, "maxfun": 3 * k, "gtol": 1e-16, "ftol": 1e-18},
        )
        M = unpack(np.asarray(res.x), base, free)
        its = int(res.nit)
        done += max(1, its)
        G = grad_w(M, cfg, w)
        fmax = float(np.max(np.abs(G[free])))
        E = sum(energy_parts(M, cfg, w))
        chunks.append(
            {"iters": done, "E": E, "fmax": fmax, "chunk_iters": its, "scipy": str(res.message)}
        )
        log(f"{tag} LBFGS its {done} E {E:.9e} fmax {fmax:.2e} (gate {gate:.1e})")
        if fmax < gate:
            verdict = "AT_GATE"
            break
        if its < 3:
            verdict = "LINE_SEARCH_STALL"
            break
    return M, {
        "E_start": E_start,
        "E_end": E,
        "fmax_end": fmax,
        "iters": done,
        "verdict": verdict,
        "chunks": chunks,
    }


# ================= reads =================
def departure(M, delta):
    """the block departure per channel, relative to delta, and the z variation."""
    return {
        "M0i": float(np.max(np.abs(M[..., 0, 1:4]))) / delta,
        "director_row": float(np.max(np.abs(M[..., 3, 1:3]))) / delta,
        "M33_minus_1": float(np.max(np.abs(M[..., 3, 3] - 1.0))) / delta,
        "M00_minus_8": float(np.max(np.abs(M[..., 0, 0] + SG))) / delta,
        "z_variation": float(np.max(np.abs(M - M[:, :, :1]))) / delta,
    }


def dep_max(d):
    return max(d[k] for k in ("M0i", "director_row", "M33_minus_1", "M00_minus_8"))


def core_radius(M, cfg, delta):
    """where the pair gap is half the vacuum's, interpolated along the +x axis at z = 0 (rows y nearest 0)."""
    n, h = cfg["n"], cfg["h"]
    lam = np.linalg.eigvalsh(M[..., 1:4, 1:4])
    gap = lam[..., 1] - lam[..., 0]  # the pair (delta, 0)
    jy = n // 2  # the row nearest y = +h/2
    prof = 0.5 * (gap[:, jy, 0] + gap[:, jy - 1, 0])  # the average of the two rows about y = 0
    x = (np.arange(n) - (n - 1) / 2.0) * h
    half = n // 2
    xs, gs = x[half:], prof[half:]
    idx = np.where(gs >= 0.5 * delta)[0]
    if len(idx) == 0 or idx[0] == 0:
        return None
    i = idx[0]
    x0, x1, g0, g1 = xs[i - 1], xs[i], gs[i - 1], gs[i]
    return float(x0 + (0.5 * delta - g0) / max(g1 - g0, 1e-300) * (x1 - x0))


# ================= jobs =================
def job_tag(j):
    return f"d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}"


def jobs_all():
    jobs = []
    for delta in (0.3, 0.1, 0.03, 0.01):
        for w1s in (6.25, 25.0, 100.0):
            jobs.append({"delta": delta, "w1s": w1s, "n": 48, "L": 48.0, "arm": "ladder"})
    for delta in (0.3, 0.03):
        for n in (32, 64, 96):  # h 1.5, 0.75, 0.5 (h 1.0 is the ladder row)
            jobs.append({"delta": delta, "w1s": 25.0, "n": n, "L": 48.0, "arm": "h"})
        for n, L in ((72, 72.0), (96, 96.0)):  # h 1.0, L 72, 96 (L 48 is the ladder row)
            jobs.append({"delta": delta, "w1s": 25.0, "n": n, "L": L, "arm": "L"})
    return jobs


def run_job(j, kicked=True, fire_iters=FIRE_ITERS, lbfgs_iters=LBFGS_ITERS):
    t0 = time.time()
    tag = job_tag(j) + ("" if kicked else "_unkicked")
    n, L, delta, w = j["n"], j["L"], j["delta"], W1 * j["w1s"]
    cfg = slab_cfg(n, L, delta)
    h = cfg["h"]
    row = dict(j, tag=tag, h=h, w=w, kicked=kicked)
    try:
        free = free_mask_slab(n, NZ, h)
        Tb = T_bps(delta, w)
        gate = max(GATE_REL * Tb / h, GATE_FLOOR)
        row.update(
            T_bps=Tb,
            kappa_bps=kappa_bps(delta, w),
            gate=gate,
            gate_at_floor=bool(gate == GATE_FLOOR),
        )
        M0 = bps_field(n, NZ, h, delta, w)
        eu0, ev0 = energy_parts(M0, cfg, w)
        row["T_seed_bps"] = (eu0 + ev0) / (NZ * h)
        if kicked:
            rng = np.random.default_rng(zlib.crc32(tag.encode()))
            kick = KICK * delta * rng.standard_normal(M0.shape)
            kick[..., 0, 0] = 0.0  # M_00 is stiff under V4 (the R20 trap): not a departure channel
            kick[..., 0, 1:] = 0.0  # the static sector (see static_sector): no time-space kick
            kick[..., 1:, 0] = 0.0
            kick = B3.sym4(kick) * free[..., None, None]
            # RUN-TIME DEVIATION 2 (2026-09-23): the kick is sized against the BOUND, an excess
            # of KICK_EXCESS x T_bps above the seed, not against delta: at 0.02 delta the kicked
            # seed sat 21x the bound at delta 0.3 and 15000x at delta 0.01 (the kick energy
            # scales as delta^2, the bound as delta^4), beyond what FIRE + L-BFGS relax
            for _ in range(3):
                eu1, ev1 = energy_parts(M0 + kick, cfg, w)
                exc = (eu1 + ev1) / (NZ * h) - row["T_seed_bps"]
                if exc <= 0:
                    break
                kick = kick * np.sqrt(KICK_EXCESS * row["T_bps"] / exc)
            M0 = M0 + kick
        eu1, ev1 = energy_parts(M0, cfg, w)
        row["T_seed_kicked"] = (eu1 + ev1) / (NZ * h)
        row["departure_seed"] = departure(M0, delta)
        os.makedirs(OUT_NPZ, exist_ok=True)
        M, fr = fire_slab(M0, cfg, w, free, fire_iters, tag=tag)
        row["fire"] = fr
        M, lb = lbfgs_slab(M, cfg, w, free, gate, lbfgs_iters, tag=tag)
        row["lbfgs"] = lb
        eu, ev = energy_parts(M, cfg, w)
        T = (eu + ev) / (NZ * h)
        row.update(
            E_u_per_len=eu / (NZ * h),
            V4_per_len=ev / (NZ * h),
            T=T,
            T_over_delta4=T / delta**4,
            T_over_T_bps=T / Tb,
            T_1d=T_1d_min(delta, w),
            departure=departure(M, delta),
            core_radius=core_radius(M, cfg, delta),
            core_radius_bps=float(np.sqrt(np.log(2.0) / kappa_bps(delta, w))),
            fire_stall=bool(abs(lb["E_end"] - fr["E_end"]) > 1e-4 * max(abs(fr["E_end"]), 1e-300)),
            lbfgs_drop_rel=float((fr["E_end"] - lb["E_end"]) / max(abs(fr["E_end"]), 1e-300)),
        )
        row["T_over_T_1d"] = T / row["T_1d"]
        row["departure_max"] = dep_max(row["departure"])
        row["status"] = "OK" if np.all(np.isfinite(M)) else "NONFINITE"
        row["verdict"] = lb["verdict"]
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M.astype(np.float64))
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc())
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} {row.get('status')} {row.get('verdict')} T/T_bps {row.get('T_over_T_bps')} "
        f"dep {row.get('departure_max')} wall {row['wall_s']}"
    )
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R25-1", "rows": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def _run_one(j):
    row = run_job(j, kicked=True)
    if row.get("status") == "OK" and row["verdict"] == "AT_GATE" and row["T_over_T_bps"] > 1.02:
        # above the bound at the gate: the optimizer, not the strand; once more unkicked
        row2 = run_job(j, kicked=False)
        row2["reseeded_unkicked"] = True
        row2["kicked_row"] = {
            k: row.get(k) for k in ("T", "T_over_T_bps", "departure_max", "verdict")
        }
        row2["tag"] = row["tag"]
        row2["instrument_fault"] = bool(row2.get("T_over_T_bps", 9) > 1.02)
        return row2
    return row


def run_pool(workers, jobs=None):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in (jobs or jobs_all()) if rows.get(job_tag(j), {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(_run_one, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            J = load_json()
            J["rows"][row["tag"]] = row
            save_json(J)
    log("pool done")


# ================= collect =================
def classify(rows):
    """the pre-registered labels on the rows (dicts with T_over_T_bps, departure_max, verdict, h, delta, w1s)."""
    ok = [r for r in rows.values() if r.get("status") == "OK" and r.get("verdict") == "AT_GATE"]
    in_block = bool(ok) and all(
        r["departure_max"] < 1e-6 and (r["h"] > 1.0 or abs(r["T_over_T_bps"] - 1.0) <= 0.02)
        for r in ok
    )
    leaves = []
    for r in ok:
        if r["T_over_T_bps"] < 0.98 and r["departure_max"] > 1e-3 and r["h"] == 1.0:
            twin = [
                q
                for q in ok
                if q["delta"] == r["delta"]
                and q["w1s"] == r["w1s"]
                and q["h"] == 0.5
                and q["T_over_T_bps"] < 0.98
                and q["departure_max"] > 1e-3
            ]
            if twin:
                leaves.append(r["tag"])
    n_rows = len(rows)
    if leaves:
        label = "STRAND_LEAVES_BLOCK"
    elif in_block and len(ok) == n_rows and n_rows > 0:
        label = "STRAND_IN_BLOCK"
    else:
        label = "UNRESOLVED"
    return {
        "label": label,
        "rows_at_gate": len(ok),
        "rows": n_rows,
        "leaves_rows": leaves,
        "faults": [r["tag"] for r in rows.values() if r.get("instrument_fault")],
    }


def collect():
    J = load_json()
    rows = J["rows"]
    out = {"classify": classify(rows)}
    lad = {}
    for r in rows.values():
        if r.get("status") != "OK":
            continue
        key = f"d{r['delta']:g}_w{r['w1s']:g}"
        lad.setdefault(key, []).append(
            {
                k: r.get(k)
                for k in (
                    "tag",
                    "h",
                    "L",
                    "T",
                    "T_over_delta4",
                    "T_over_T_bps",
                    "T_over_T_1d",
                    "departure_max",
                    "core_radius",
                    "core_radius_bps",
                    "verdict",
                    "fire_stall",
                )
            }
        )
    out["ladder"] = lad
    # RUN-TIME READ: the slaved bound per (delta, w) and each row against it
    slaved = {}
    for r in rows.values():
        if r.get("status") != "OK":
            continue
        key = (r["delta"], r["w1s"])
        if key not in slaved:
            slaved[key] = T_slaved(r["delta"], W1 * r["w1s"])
        r["T_slaved"] = slaved[key]
        r["T_over_T_slaved"] = r["T"] / slaved[key]
    for key, lst in lad.items():
        for e in lst:
            e["T_over_T_slaved"] = rows[e["tag"]].get("T_over_T_slaved")
    out["slaved_bound"] = {
        f"d{d:g}_w{ws:g}": {
            "T_slaved": t,
            "T_bps_lead": T_bps(d, W1 * ws),
            "T_slaved_over_T_bps": t / T_bps(d, W1 * ws),
            "m33_at_f0_over_4": V_slaved(
                (d / 2) ** 2 / 4, d / 2, W1 * ws, [(-8.0) ** q + 1.0 + d**q for q in (1, 2, 3, 4)]
            )[2],
        }
        for (d, ws), t in sorted(slaved.items())
    }
    out["slaved_bound"]["note"] = (
        "K_eff(s0) = 4 (175936 s0^4 + 1262016 s0^3 + 1677108 s0^2 - 2064132 s0 + 642249) / 845261 at "
        "quadratic order (sympy, the two shifts solved), limit 3.039 = 0.760 of the frozen 4 as delta -> 0, "
        "so T_slaved / T_bps -> 0.872; the shift is the director eigenvalue (M_33 - 1 = -0.25 (f - f0) at "
        "small delta), M_00 moves by 4e-4 (f - f0)"
    )
    # the laws, reported: T / delta^4 across delta and T / sqrt(w) across w at n 48 L 48
    base = {
        (r["delta"], r["w1s"]): r
        for r in rows.values()
        if r.get("status") == "OK" and r["n"] == 48 and r["L"] == 48.0
    }
    out["delta_law_at_w25"] = {
        f"{d:g}": {
            "T_over_delta4": base[(d, 25.0)]["T_over_delta4"],
            "predicted": T_bps(d, W1 * 25.0) / d**4,
        }
        for d in (0.3, 0.1, 0.03, 0.01)
        if (d, 25.0) in base
    }
    out["w_law_at_delta"] = {
        f"{d:g}": {
            f"{ws:g}": {
                "T_over_sqrt_w": base[(d, ws)]["T"] / np.sqrt(W1 * ws),
                "predicted": T_bps(d, W1 * ws) / np.sqrt(W1 * ws),
            }
            for ws in (6.25, 25.0, 100.0)
            if (d, ws) in base
        }
        for d in (0.3, 0.03)
    }
    # the h ladder: the discretization order of T against T_1d
    hl = {}
    for d in (0.3, 0.03):
        pts = sorted(
            [
                (r["h"], r["T"], r["T_1d"])
                for r in rows.values()
                if r.get("status") == "OK"
                and r["delta"] == d
                and r["w1s"] == 25.0
                and r["L"] == 48.0
            ]
        )
        if len(pts) >= 3:
            hs = np.array([p[0] for p in pts])
            err = np.array([abs(p[1] - p[2]) / p[2] for p in pts])
            good = err > 0
            order = (
                float(np.polyfit(np.log(hs[good]), np.log(err[good]), 1)[0])
                if good.sum() >= 2
                else None
            )
            hl[f"{d:g}"] = {"points_h_T_T1d": pts, "order": order}
    out["h_ladder"] = hl
    J["collect"] = out
    save_json(J)
    print(json.dumps(out, indent=1))
    return out


def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    J = load_json()
    rows = [r for r in J["rows"].values() if r.get("status") == "OK"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for ws, mk in ((6.25, "o"), (25.0, "s"), (100.0, "^")):
        pts = sorted(
            [
                (r["delta"], r["T_over_T_bps"])
                for r in rows
                if r["w1s"] == ws and r["n"] == 48 and r["L"] == 48.0
            ]
        )
        if pts:
            ax[0].plot([p[0] for p in pts], [p[1] for p in pts], mk + "-", label=f"W1 x {ws:g}")
    ax[0].axhline(1.0, color="k", lw=0.8)
    # the slaved bound (M_00 and M_33 relaxed per radius), from the collect
    sl = (J.get("collect") or {}).get("slaved_bound") or {}
    spts = sorted(
        (v["T_bps_lead"] and float(k.split("_")[0][1:]), v["T_slaved_over_T_bps"])
        for k, v in sl.items()
        if k != "note" and k.endswith("_w25")
    )
    if spts:
        ax[0].plot(
            [p[0] for p in spts],
            [p[1] for p in spts],
            "k--",
            lw=1.2,
            label="slaved bound (M_33, M_00 relaxed)",
        )
    ax[0].set_xscale("log")
    ax[0].set_xlabel("delta")
    ax[0].set_ylabel("T / T_BPS")
    ax[0].legend()
    for d, mk in ((0.3, "o"), (0.03, "s")):
        pts = sorted(
            [
                (r["h"], r.get("T_over_T_slaved") or r["T"] / r["T_1d"])
                for r in rows
                if r["delta"] == d and r["w1s"] == 25.0 and r["L"] == 48.0
            ]
        )
        if pts:
            ax[1].plot([p[0] for p in pts], [p[1] for p in pts], mk + "-", label=f"delta {d:g}")
    ax[1].axhline(1.0, color="k", lw=0.8)
    ax[1].set_xlabel("h")
    ax[1].set_ylabel("T / T_slaved")
    ax[1].legend()
    fig.suptitle("R25-1: the straight strand against the Bogomolny tension")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "m5_32_r25_1_strand.png"), dpi=130)
    log("plot written")


# ================= smoke =================
def smoke():
    out = {}
    # (1) the bound against the 1D minimizer and the stack's BPS energy on a small slab
    d, w = 0.3, W1 * 25.0
    Tb, T1 = T_bps(d, w), T_1d_min(d, w)
    n, L = 24, 24.0
    cfg = slab_cfg(n, L, d)
    M = bps_field(n, NZ, 1.0, d, w)
    eu, ev = energy_parts(M, cfg, w)
    Ts = (eu + ev) / NZ
    out["bound"] = {
        "T_bps": Tb,
        "T_1d": T1,
        "T_stack_h1_L24": Ts,
        "T_1d_over_T_bps": T1 / Tb,
        "T_stack_over_T_1d": Ts / T1,
    }
    ok1 = 0.99 < T1 / Tb < 1.01 and 0.97 < Ts / T1 < 1.01
    # (2) the winding additivity: the m 2 BPS field costs twice the m 1 one
    M2 = bps_field(n, NZ, 1.0, d, w, m=2)
    eu2, ev2 = energy_parts(M2, cfg, w)
    out["additivity"] = {"T_m2_over_T_m1": (eu2 + ev2) / (eu + ev), "expected": 2.0}
    ok2 = abs(out["additivity"]["T_m2_over_T_m1"] - 2.0) < 0.05
    # (3) the pure-winding field has no exterior density: E_u vanishes with b constant
    # (the four cells about the axis carry the unresolved sign jump on the lattice; beyond 4 h the
    # density is a lattice residue under 1e-4 of theirs, and on the BPS field exponentially localized)
    Mw = bps_field(n, NZ, 1.0, d, w, kappa=1e6)
    ew = density_u(Mw, cfg)
    eb = density_u(M, cfg)
    _, _, rho, _ = slab_coords(n, 1.0)
    far = (rho > 4.0)[:, :, None] & np.ones((1, 1, NZ), dtype=bool)
    out["exterior_density"] = {
        "pure_winding_max_beyond_4h": float(np.max(np.abs(ew[far]))),
        "pure_winding_max_axis_cells": float(np.max(np.abs(ew))),
        "bps_max_beyond_3_core_radii": float(
            np.max(
                np.abs(eb[(rho > 3.0 * np.sqrt(np.log(2.0) / kappa_bps(d, w)))[:, :, None] & far])
            )
        ),
        "bps_max": float(np.max(np.abs(eb))),
    }
    # (the lattice residue of the pure winding is O(h^4) per cell and falls as rho^-6: R25-0 b reads it)
    ok3 = out["exterior_density"]["pure_winding_max_beyond_4h"] < 1e-4 * max(
        out["exterior_density"]["pure_winding_max_axis_cells"], 1e-300
    )
    # (4) the wiring: a kicked seed through 3 FIRE and 3 L-BFGS iterations, the reads
    j = {"delta": d, "w1s": 25.0, "n": 16, "L": 16.0, "arm": "smoke"}
    global OUT_NPZ
    keep = OUT_NPZ
    OUT_NPZ = os.path.join(DATA, "m5_32_r25_1_smoke_tmp")
    try:
        row = run_job(j, kicked=True, fire_iters=3, lbfgs_iters=3)
    finally:
        import shutil

        shutil.rmtree(OUT_NPZ, ignore_errors=True)
        OUT_NPZ = keep
    out["wiring"] = {
        k: row.get(k)
        for k in (
            "status",
            "verdict",
            "T_seed_bps",
            "T_seed_kicked",
            "T",
            "T_over_T_bps",
            "departure_max",
            "core_radius",
            "gate",
        )
    }
    out["wiring"]["departure_seed_max"] = dep_max(row["departure_seed"])
    ok4 = (
        row.get("status") == "OK"
        and out["wiring"]["departure_seed_max"] > 1e-3
        and row["T_seed_kicked"] > row["T_seed_bps"]
    )

    # (5) classify on synthetic rows
    def syn(tag, ratio, dep, h=1.0, delta=0.3, w1s=25.0, verdict="AT_GATE"):
        return {
            "tag": tag,
            "status": "OK",
            "verdict": verdict,
            "T_over_T_bps": ratio,
            "departure_max": dep,
            "h": h,
            "delta": delta,
            "w1s": w1s,
        }

    in_block = {t: syn(t, 1.003, 1e-8) for t in ("a", "b")}
    leaves = {"a": syn("a", 0.95, 1e-2, h=1.0), "b": syn("b", 0.95, 1e-2, h=0.5)}
    unres = {"a": syn("a", 1.003, 1e-8), "b": syn("b", 1.05, 1e-8, verdict="FALLING")}
    out["classify"] = {
        "in_block": classify(in_block)["label"],
        "leaves": classify(leaves)["label"],
        "unresolved": classify(unres)["label"],
    }
    ok5 = out["classify"] == {
        "in_block": "STRAND_IN_BLOCK",
        "leaves": "STRAND_LEAVES_BLOCK",
        "unresolved": "UNRESOLVED",
    }
    tags = [job_tag(j) for j in jobs_all()]
    out["jobs"] = {"count": len(tags), "unique": len(set(tags)) == len(tags)}
    out["PASS"] = bool(
        ok1 and ok2 and ok3 and ok4 and ok5 and out["jobs"]["unique"] and len(tags) == 22
    )
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12)
    elif mode == "collect":
        collect()
    elif mode == "plot":
        plot()
    elif mode == "jobs":
        for j in jobs_all():
            print(job_tag(j), j)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
