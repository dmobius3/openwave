"""M5.21.1 phase P1: deep 3D hedgehog statics vs HIS three structural
claims + the stability bar (Duda 2026-07-15, tasks/m5_20_convo.md).

THE OBJECT + THE FUNCTIONAL (all audited upstream, m5_20_2_a_eom):
    E_static = INT w [ u_eta + V4 ],
    u_eta = 4 SUM_{i<j} <[A_i, A_j]_eta, [A_i, A_j]_eta>_eta   (axisym
    channels A_rho, A_phi = [J, M]/rho, A_z),
    V4 = wscale SUM_{p=1..4} (Tr_eta(M^p) - C_p)^2,
    C_p = g^p + 1 + delta^p  (s = +1; P0 measured the statics of
    block-diagonal seeds sign-mirror-equivalent, so P1 runs the
    verified branch), FIRE descent on the free DOFs (outer pin).

HIS THREE STRUCTURAL CLAIMS (2026-07-15, measured, not imposed)
    C1  "vortex in z axis - requires regularization with two equal
        eigenvalues": the transverse eigenvalue pair on the axis
        equalizes (vacuum far-field anisotropy is delta - 0 = delta)
    C2  "in the center all three spatial eigenvalues have to be
        equal": the 3-equal core; deep relax decides whether it
        SURVIVES (M5.21 found a seed-adjacent transient + a Q8 slide)
    C3  "511keV mass mainly in the center, vortex should be extremely
        light": E(<r) containment (r50/r90), the core-ball energy
        fraction, and the vortex-line energy per unit length

THE STABILITY BAR (his own: the electron "has to be stable")
    SB1  FIRE convergence: force drop + E tail slope + center drift
         (the slide read: drift rate ~0 = pinned, else the Q8 slide
         survives his regularization language)
    SB2  Hessian lowest modes at the endpoint (Lanczos, HV by complex
         step, FD-gated): sector-split into (a) the time-BLOCK-
         DIAGONAL sector (the statics stability read) and (b) the
         full 10-dim per-cell space incl. time-mixing (the 4D-sector
         hazard context, feeds P2); negative modes counted per sector
    SB3  3D spot-check (the m5_21_d stack, STATIC descent): embed the
         axisym endpoint in a 48^3 box, add the l=2 non-axisym bump,
         FIRE on the 3D gradient: a_break(it) decaying / flat = no
         lower non-axisym basin at this resolution; growing = the
         axisym reduction froze a real descent direction

FILM  the relaxation sequence renders through m5_film.film_strip in
      BOTH templates (basic + thermal), first row t = 0 (the seed).

Run:  python m5_21_1_b_statics.py b1 [iters] | b2 | b3 [N iters] | all
Out:  ../data/m5_21_1_b_statics.json, ../data/m5_21_1_b_endpoint.npz,
      ../plots/m5_21_1_b_film_basic.png / _film_thermal.png /
      _profiles.png
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ARGV = sys.argv[1:]
sys.argv = sys.argv[:1]

from m5_16_axisym import pin_mask                                   # noqa: E402
from m5_17_energy import cell_weights, grid_coords                  # noqa: E402
from m5_20_2_a_eom import (G_T, H, NR, NZ, WSCALE, grad_static_4,   # noqa: E402
                           total_energy_4, u_eta_density, v4_density)
from m5_21_b_electron import (core_spectrum, defect_center,         # noqa: E402
                              electron_seed, meridional_charge)
from m5_film import film_strip                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
DELTA = 0.3
A_CORE = (1.0 + DELTA) / 3.0
VAC_SPEC = (G_T, 1.0, DELTA, 0.0)
SNAPS = (0, 2000, 6000, 12000, 24000, 48000)

W4 = cell_weights(NR, NZ, H)[..., None, None]
RHO4 = ((np.arange(NR - 1) + 0.5) * H)[:, None, None, None]
PIN = pin_mask(NR, NZ)
FREE4 = (~PIN)[..., None, None].astype(float)


def e_fn(M):
    return total_energy_4(M, WSCALE, DELTA)


def g_fn(M):
    return grad_static_4(M, WSCALE, DELTA, w4=W4, rho4=RHO4)


# ================= observables =================
def axis_two_equalness(M, n_rows=2):
    """C1: transverse eigenvalue split on the z axis. For the first
    n_rows rho rows: eigen-decompose the spatial block, take the
    eigenvector most aligned with z as the axis eigenvalue, report the
    TRANSVERSE pair difference d(z) = lam_t1 - lam_t2 (vacuum combed
    frame would give delta - 0 = delta; his regularization = 0)."""
    S = M[:n_rows, :, 1:4, 1:4]
    lam, V = np.linalg.eigh(0.5 * (S + np.swapaxes(S, -1, -2)))
    zdot = np.abs(V[..., 2, :])                     # |v . zhat| per mode
    ax = np.argmax(zdot, axis=-1)
    idx = np.arange(3)[None, None, :]
    tr_mask = idx != ax[..., None]
    lam_t = lam[tr_mask].reshape(n_rows, NZ, 2)
    d = np.abs(lam_t[..., 1] - lam_t[..., 0]).mean(axis=0)
    z = (np.arange(NZ) - NZ / 2 + 0.5) * H
    far = np.abs(z) > 20.0
    mid = (np.abs(z) > 5.0) & (np.abs(z) < 60.0)
    return {"z": z.tolist(), "d_transverse": d.tolist(),
            "d_far_median": float(np.median(d[far])),
            "d_mid_median": float(np.median(d[mid])),
            "vacuum_reference": DELTA}


def containment(M, r_tube=4.0):
    """C3: spherical containment + the vortex-line energy per length."""
    R, Z = grid_coords(NR, NZ, H)
    w = cell_weights(NR, NZ, H)
    dens = (u_eta_density(M, H) + v4_density(M, WSCALE, DELTA)) * w
    Rc, Zc = R[: NR - 1, 1:-1], Z[: NR - 1, 1:-1]
    r = np.sqrt(Rc ** 2 + Zc ** 2)
    E_tot = float(dens.sum())
    rb = np.arange(1.0, 100.0, 1.0)
    csum = np.array([dens[r <= x].sum() for x in rb])
    def r_at(frac):
        i = int(np.searchsorted(csum, frac * csum[-1]))
        return float(rb[min(i, len(rb) - 1)])
    # line energy per unit z inside the rho-tube, away from the core
    tube = Rc <= r_tube
    e_line = np.array([dens[:, j][tube[:, j]].sum()
                       for j in range(dens.shape[1])])
    z = Zc[0]
    far = (np.abs(z) > 20.0) & (np.abs(z) < 80.0)
    lam_line = float(np.median(e_line[far]))
    core_ball = float(dens[r <= 8.0].sum())
    tube_far = float(e_line[np.abs(z) > 8.0].sum())
    return {"E_total": E_tot, "r50": r_at(0.5), "r90": r_at(0.9),
            "core_ball_r8_frac": core_ball / E_tot,
            "vortex_line_energy_per_length": lam_line,
            "vortex_tube_beyond_r8_frac": tube_far / E_tot,
            "e_line_z": e_line.tolist(), "z": z.tolist(),
            "cum_r": rb.tolist(), "cum_E": csum.tolist()}


def snapshot_obs(M, it):
    return {"it": it, "E": e_fn(M),
            "q_meridional": meridional_charge(M),
            "core": core_spectrum(M),
            "center": defect_center(M),
            "axis": {k: v for k, v in axis_two_equalness(M).items()
                     if not isinstance(v, list)},
            "containment": {k: v for k, v in containment(M).items()
                            if not isinstance(v, list)}}


# ================= FIRE with snapshots =================
def fire_relax_snap(M0, snaps=SNAPS, dt0=0.02, dt_max=0.2,
                    log_every=1000):
    """the m5_21_b FIRE (same algorithm, same knobs) with state capture
    at the pre-registered snapshot iterations."""
    w = np.ones((NR, NZ))
    w[: NR - 1, 1:-1] = cell_weights(NR, NZ, H)
    precond = (1.0 / w)[..., None, None]
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha = dt0, 0.1
    n_up = 0
    hist = []
    states = [{"it": 0, "t": 0.0, "M": M0.copy(), "vac_spec": VAC_SPEC}]
    F = -g_fn(M) * precond * FREE4
    f0 = float(np.max(np.abs(F)))
    max_iter = max(snaps)
    t_start = time.time()
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
        F = -g_fn(M) * precond * FREE4
        if it % log_every == 0 or it == max_iter:
            E = e_fn(M)
            fmax = float(np.max(np.abs(F)))
            hist.append({"it": it, "E": E, "fmax": fmax, "dt": dt})
            print(f"  it {it:6d} E {E:12.6f} fmax {fmax:.3e} "
                  f"dt {dt:.3f} [{time.time() - t_start:.0f}s]",
                  flush=True)
        if it in snaps:
            states.append({"it": it, "t": float(it), "M": M.copy(),
                           "vac_spec": VAC_SPEC})
    return M, states, {"fmax_seed": f0, "trace": hist}


def phase_b1(iters=None):
    snaps = SNAPS if iters is None else tuple(
        int(x) for x in np.unique(np.round(
            np.array(SNAPS) * iters / max(SNAPS))))
    M0 = electron_seed()
    print(f"[B1] deep relax, snaps {snaps}", flush=True)
    M, states, rel = fire_relax_snap(M0, snaps=snaps)
    obs = [snapshot_obs(st["M"], st["it"]) for st in states]
    # the slide read: center drift + E tail slope over the last segment
    c = np.array([o["center"] for o in obs])
    its = np.array([o["it"] for o in obs], dtype=float)
    drift = (np.linalg.norm(c[-1] - c[-2])
             / max(its[-1] - its[-2], 1.0))
    tr = rel["trace"]
    tail = [r for r in tr if r["it"] >= 0.5 * its[-1]]
    dE_dit = ((tail[-1]["E"] - tail[0]["E"])
              / max(tail[-1]["it"] - tail[0]["it"], 1.0))
    out = {"snapshots": obs,
           "slide": {"center_drift_per_iter_tail": float(drift),
                     "dE_per_iter_tail": float(dE_dit),
                     "force_drop": rel["fmax_seed"]
                     / max(tr[-1]["fmax"], 1e-300)},
           "fire": rel}
    np.savez_compressed(os.path.join(DATA, "m5_21_1_b_endpoint.npz"), M=M)
    # films: BOTH templates (the film standard)
    film_strip(states, os.path.join(PLOTS, "m5_21_1_b_film_basic.png"),
               template="basic", delta=DELTA, h=H, g=G_T, wscale=WSCALE,
               suptitle="M5.21.1-P1 electron hedgehog deep statics "
                        "(FIRE iterations as frames)")
    film_strip(states, os.path.join(PLOTS, "m5_21_1_b_film_thermal.png"),
               template="thermal", delta=DELTA, h=H, g=G_T,
               suptitle="M5.21.1-P1 electron hedgehog deep statics "
                        "(thermal panel set)", step=8)
    # profiles figure
    end = states[-1]["M"]
    ax_prof = axis_two_equalness(end)
    cont = containment(end)
    fig, axes = plt.subplots(1, 4, figsize=(18, 3.8))
    axes[0].plot(cont["cum_r"], np.array(cont["cum_E"])
                 / cont["cum_E"][-1], lw=1.2)
    for frac, rr in (("r50", cont["r50"]), ("r90", cont["r90"])):
        axes[0].axvline(rr, ls="--", lw=0.8, color="gray")
    axes[0].set_title(f"E(<r)/E: r50 {cont['r50']:.0f}, "
                      f"r90 {cont['r90']:.0f}", fontsize=9)
    axes[0].set_xlabel("r")
    axes[1].plot(ax_prof["z"], ax_prof["d_transverse"], lw=1.0)
    axes[1].axhline(DELTA, ls="--", lw=0.8, color="r",
                    label="combed vacuum (delta)")
    axes[1].set_title(f"C1 axis transverse split: mid-median "
                      f"{ax_prof['d_mid_median']:.3f}", fontsize=9)
    axes[1].set_xlabel("z")
    axes[1].legend(fontsize=7)
    axes[2].semilogy(cont["z"], np.maximum(cont["e_line_z"], 1e-12),
                     lw=1.0)
    axes[2].set_title(f"C3 tube (rho<4) energy per z: line "
                      f"{cont['vortex_line_energy_per_length']:.2e}",
                      fontsize=9)
    axes[2].set_xlabel("z")
    core_a = [o["core"]["core_mean"] for o in obs]
    core_s = [o["core"]["core_spread"] for o in obs]
    axes[3].plot(its, core_a, "o-", label="core mean a")
    axes[3].plot(its, core_s, "s-", label="core spread")
    axes[3].axhline(A_CORE, ls="--", lw=0.8, color="gray",
                    label="(1+delta)/3")
    axes[3].set_title("C2 core spectrum vs iteration", fontsize=9)
    axes[3].set_xlabel("FIRE iteration")
    axes[3].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_1_b_profiles.png"), dpi=130)
    plt.close(fig)
    out["profiles"] = {"axis": ax_prof, "containment": cont}
    print(f"[B1] E {obs[-1]['E']:.5f} q {obs[-1]['q_meridional']:.4f} "
          f"core a {obs[-1]['core']['core_mean']:.4f} spread "
          f"{obs[-1]['core']['core_spread']:.4f} r50/r90 "
          f"{cont['r50']:.0f}/{cont['r90']:.0f} axis-split(mid) "
          f"{ax_prof['d_mid_median']:.4f} drift "
          f"{out['slide']['center_drift_per_iter_tail']:.2e}/it",
          flush=True)
    return out


# ================= SB2: Hessian lowest modes =================
IU = np.triu_indices(4)
OFFD = (IU[0] != IU[1]).astype(float)
SC = 1.0 + OFFD * (np.sqrt(2.0) - 1.0)      # Frobenius-metric packing
TIME_MIX = np.array([i for i in range(10)
                     if IU[0][i] == 0 and IU[1][i] > 0])


def _pack(M):
    return (M[..., IU[0], IU[1]] * SC)[~PIN].ravel()


def _unpack(x):
    M = np.zeros((NR, NZ, 4, 4))
    blk = np.zeros(((~PIN).sum(), 10))
    blk[:] = x.reshape(-1, 10) / SC
    full = np.zeros(((~PIN).sum(), 4, 4))
    full[:, IU[0], IU[1]] = blk
    full[:, IU[1], IU[0]] = blk
    M[~PIN] = full
    return M


def hv_cs(M0, V, eps=1e-30):
    """Hessian-vector product by complex step on the (polynomial)
    gradient."""
    G = grad_static_4(M0.astype(complex) + 1j * eps * V, WSCALE, DELTA,
                      w4=W4, rho4=RHO4)
    return np.imag(G) / eps * FREE4


def phase_b2(k=3):
    from scipy.sparse.linalg import LinearOperator, eigsh
    M0 = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))["M"]
    rng = np.random.default_rng(11)
    # gate: complex-step HV vs FD directional derivative of the gradient
    worst = 0.0
    for _ in range(2):
        D = np.zeros((NR, NZ, 4, 4))
        Dc = rng.normal(size=((~PIN).sum(), 4, 4))
        Dc = 0.5 * (Dc + np.swapaxes(Dc, -1, -2))
        D[~PIN] = Dc
        hv = hv_cs(M0, D)
        hh = 1e-6
        fd = (grad_static_4(M0 + hh * D, WSCALE, DELTA, w4=W4, rho4=RHO4)
              - grad_static_4(M0 - hh * D, WSCALE, DELTA, w4=W4,
                              rho4=RHO4)) / (2 * hh) * FREE4
        worst = max(worst, float(np.max(np.abs(hv - fd))
                                 / max(np.max(np.abs(fd)), 1e-300)))
    print(f"[SB2 gate] complex-step HV vs FD rel {worst:.2e}", flush=True)
    out = {"hv_gate_rel": worst, "sectors": {}}
    n = int((~PIN).sum() * 10)

    def spectrum(sector):
        mask = np.ones(10)
        if sector == "block_diag":
            mask[TIME_MIX] = 0.0

        def mv(x):
            xm = (x.reshape(-1, 10) * mask).ravel()
            y = _pack(hv_cs(M0, _unpack(xm)))
            return (y.reshape(-1, 10) * mask).ravel()

        L = LinearOperator((n, n), matvec=mv)
        t0 = time.time()
        try:
            vals = eigsh(L, k=k, which="SA", tol=1e-3, maxiter=250,
                         return_eigenvectors=False)
            note = "converged"
        except Exception as ex:                              # noqa: BLE001
            vals = getattr(ex, "eigenvalues", np.array([]))
            note = f"partial: {type(ex).__name__}"
        return {"lowest": sorted(np.asarray(vals).tolist()),
                "note": note, "wall_s": round(time.time() - t0, 1)}

    for sector in ("block_diag", "full"):
        out["sectors"][sector] = spectrum(sector)
        print(f"[SB2 {sector}] lowest {out['sectors'][sector]['lowest']} "
              f"({out['sectors'][sector]['note']}, "
              f"{out['sectors'][sector]['wall_s']}s)", flush=True)
    return out


def phase_b2x():
    """SB2x: the convergent replacement for the eigsh probe (the first
    attempt hit ArpackNoConvergence in both sectors: SA on a 327k-dim
    indefinite operator without shift-invert; kept in the JSON as the
    negative record). Two reads:
      (i)  directional curvatures D2E(v)/||v||_F^2 on the 128x256
           endpoint along INFORMATIVE directions (the slide force,
           a twist, boost textures, randoms per sector);
      (ii) extreme-eigenvalue ESTIMATES per sector by shifted power
           iteration on the 64x128 c0 state (4x cheaper matvec),
           residual-qualified."""
    from m5_21_1_c_4d import FREE42 as F42
    from m5_21_1_c_4d import RHO42 as R42
    from m5_21_1_c_4d import W42 as Wt42
    from m5_21_1_c_4d import NR2, NZ2, PIN2, hv_cs2
    from m5_21_c_clockrun import local_gen
    from m5_21_a_snap import eig_fields
    out = {}
    # ---- (i) directional probes on the endpoint ----
    M0 = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))["M"]
    R, Z = grid_coords(NR, NZ, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    env = np.exp(-((r / 20.0) ** 4))[..., None, None]
    rng = np.random.default_rng(23)
    F = -g_fn(M0) * FREE4
    _, V = eig_fields(M0)
    W = local_gen(V[..., :, 0])
    bump = np.exp(-(((R - 2.0) ** 2 + Z ** 2)
                    / 18.0))[..., None, None]
    dirs = {"slide_force": F,
            "twist_local": env * (W @ M0 - M0 @ W)}
    B = np.zeros_like(M0)
    B[..., 0, 1] = bump[..., 0, 0]
    B[..., 1, 0] = bump[..., 0, 0]
    dirs["boost01_core"] = B
    for i in range(2):
        D = rng.normal(size=M0.shape)
        D = 0.5 * (D + np.swapaxes(D, -1, -2))
        Db = D.copy()
        Db[..., 0, 1:] = 0.0
        Db[..., 1:, 0] = 0.0
        dirs[f"rand_blockdiag_{i}"] = Db
        Dt = np.zeros_like(D)
        Dt[..., 0, 1:] = D[..., 0, 1:]
        Dt[..., 1:, 0] = D[..., 1:, 0]
        dirs[f"rand_timemix_{i}"] = Dt
    rows = {}
    for name, v in dirs.items():
        v = v * FREE4
        n2 = float(np.sum(v * v))
        if n2 <= 0:
            continue
        hv = hv_cs(M0, v)
        rows[name] = {"curvature": float(np.sum(v * hv)) / n2}
        print(f"[SB2x dir {name}] D2E/|v|^2 = "
              f"{rows[name]['curvature']:+.4e}", flush=True)
    out["directional_endpoint"] = rows
    # ---- (ii) shifted power iteration on the 64x128 c0 state ----
    Mc = np.load(os.path.join(DATA, "m5_21_1_c_state_sp.npz"))["M"]
    iu = np.triu_indices(4)
    offd = (iu[0] != iu[1]).astype(float)
    sc = 1.0 + offd * (np.sqrt(2.0) - 1.0)
    tmix = np.array([i for i in range(10)
                     if iu[0][i] == 0 and iu[1][i] > 0])

    def pack2(Mx):
        return (Mx[..., iu[0], iu[1]] * sc)[~PIN2].ravel()

    def unpack2(x):
        Mx = np.zeros((NR2, NZ2, 4, 4))
        blk = x.reshape(-1, 10) / sc
        full = np.zeros(((~PIN2).sum(), 4, 4))
        full[:, iu[0], iu[1]] = blk
        full[:, iu[1], iu[0]] = blk
        Mx[~PIN2] = full
        return Mx

    def spectrum_est(sector, n_top=30, n_bot=140):
        mask = np.ones(10)
        if sector == "block_diag":
            mask[tmix] = 0.0

        def mv(x):
            xm = (x.reshape(-1, 10) * mask).ravel()
            y = pack2(hv_cs2(Mc, unpack2(xm), G_T, Wt42, R42, F42))
            return (y.reshape(-1, 10) * mask).ravel()

        rng2 = np.random.default_rng(5)
        x = rng2.normal(size=int((~PIN2).sum() * 10))
        x = (x.reshape(-1, 10) * mask).ravel()
        x /= np.linalg.norm(x)
        lam_top = 0.0
        for _ in range(n_top):
            y = mv(x)
            lam_top = float(x @ y)
            x = y / max(np.linalg.norm(y), 1e-300)
        shift = 1.05 * abs(lam_top)
        x = rng2.normal(size=x.size)
        x = (x.reshape(-1, 10) * mask).ravel()
        x /= np.linalg.norm(x)
        lam_b = 0.0
        for _ in range(n_bot):
            y = shift * x - mv(x)
            lam_b = float(x @ y)
            x = y / max(np.linalg.norm(y), 1e-300)
        y = mv(x)
        lam_min = float(x @ y)
        resid = float(np.linalg.norm(y - lam_min * x))
        return {"lam_top_est": lam_top, "lam_min_est": lam_min,
                "power_residual": resid,
                "n_matvec": n_top + n_bot + 1}

    for sector in ("block_diag", "full"):
        t0 = time.time()
        est = spectrum_est(sector)
        est["wall_s"] = round(time.time() - t0, 1)
        out[f"power_{sector}_64grid"] = est
        print(f"[SB2x power {sector}] lam_top ~ {est['lam_top_est']:.3e}"
              f" lam_min ~ {est['lam_min_est']:+.3e} (resid "
              f"{est['power_residual']:.2e}, {est['wall_s']}s)",
              flush=True)
    return out


# ================= SB3: 3D spot-check (static descent) =================
def phase_b3(N=48, iters=1500):
    from m5_21_d_3dcheck import (a_break, embed_axisym, gate_gd3a,
                                 grad_static_3d, perturb_nonaxisym,
                                 total_energy_3d)
    ok3, det3 = gate_gd3a(N=24)
    print(f"[SB3 gate GD3a] {'PASS' if ok3 else 'FAIL'} "
          + json.dumps(det3, default=float)[:160], flush=True)
    Mt = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))["M"]
    M0 = embed_axisym(Mt, N)
    Mp = perturb_nonaxisym(M0, amp=0.02, sig=4.0)
    free = np.zeros((N, N, N, 1, 1))
    free[2:-2, 2:-2, 2:-2] = 1.0
    E_axi = total_energy_3d(M0, delta=DELTA)
    M = Mp.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = 0.02, 0.1, 0
    F = -grad_static_3d(M, delta=DELTA) * free
    trace = [{"it": 0, "E": total_energy_3d(M, delta=DELTA),
              "a_break": float(a_break(M))}]
    t0 = time.time()
    for it in range(1, iters + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, 0.2)
                alpha *= 0.99
        else:
            v[:] = 0.0
            dt *= 0.5
            alpha = 0.1
            n_up = 0
        v += dt * F
        M += dt * v
        F = -grad_static_3d(M, delta=DELTA) * free
        if it % 150 == 0 or it == iters:
            trace.append({"it": it,
                          "E": total_energy_3d(M, delta=DELTA),
                          "a_break": float(a_break(M))})
            print(f"  3D it {it:5d} E {trace[-1]['E']:.4f} a_break "
                  f"{trace[-1]['a_break']:.4e} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    ab = [r["a_break"] for r in trace]
    verdict = ("decaying" if ab[-1] < 0.5 * ab[0]
               else "flat" if ab[-1] < 2.0 * ab[0] else "GROWING")
    out = {"gate_gd3a": {"ok": bool(ok3), **det3}, "N": N,
           "E_axisym_embedded": E_axi, "trace": trace,
           "a_break_first_last": [ab[0], ab[-1]],
           "E_drop_vs_axisym": trace[-1]["E"] - E_axi,
           "verdict_nonaxisym": verdict}
    print(f"[SB3] a_break {ab[0]:.3e} -> {ab[-1]:.3e} ({verdict}); "
          f"E {trace[-1]['E']:.4f} vs embedded axisym {E_axi:.4f}",
          flush=True)
    return out


def main(which="all", arg1=None, arg2=None):
    os.makedirs(DATA, exist_ok=True)
    res = {"task": "M5.21.1", "phase": "P1", "delta": DELTA, "g": G_T,
           "wscale": WSCALE, "grid": [NR, NZ, H], "sign_note":
           "s = +1 branch (P0 measured block-diagonal statics "
           "sign-mirror-equivalent)"}
    fn = os.path.join(DATA, "m5_21_1_b_statics.json")
    if os.path.exists(fn):
        with open(fn) as f:
            res.update(json.load(f))
    if which in ("b1", "all"):
        res["B1"] = phase_b1(int(arg1) if arg1 else None)
    if which in ("b2", "all"):
        res["SB2"] = phase_b2()
    if which in ("b2x", "all"):
        res["SB2x"] = phase_b2x()
    if which in ("b3", "all"):
        res["SB3"] = phase_b3(int(arg1) if which == "b3" and arg1
                              else 48,
                              int(arg2) if which == "b3" and arg2
                              else 1500)
    with open(fn, "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("wrote data/m5_21_1_b_statics.json")
    return res


if __name__ == "__main__":
    main(ARGV[0] if ARGV else "all",
         ARGV[1] if len(ARGV) > 1 else None,
         ARGV[2] if len(ARGV) > 2 else None)
