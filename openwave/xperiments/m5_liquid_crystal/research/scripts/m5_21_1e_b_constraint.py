"""M5.21.1e phase B: the constraint hypothesis on the hedgehog slide.

Question (from the spec review, tasks/m5_21_1e_task_details.md): is the
M5.21.1 P1 spreading dilution the paper-anticipated escape channel that
FRAME (arXiv:2108.07896v7) flags TWICE as needing a volume-preserving
constraint (p. 8: "det(M) = const ... to prevent using only long axes
which allow for low curvature (hence energy)"; p. 12 again for 4D), and
that Faber (hep-th/9910221) closes with his hard S^3 norm constraint?

Three machine-checkable arms on the 64x128 axisym stack (s = +1 branch,
the M5.21.1 P1 working branch):

  d  DIAGNOSIS   re-run the unconstrained slide with snapshots and read,
                 per snapshot: the Derrick virial ratio u/(3 V4), the
                 fraction of the static force that lies IN the
                 isospectral (orbit) tangent space, and the spectral
                 drift away from the vacuum targets. Amplitude-mode
                 hypothesis: the force is mostly tangent-ORTHOGONAL
                 (spectral) and the drift grows monotonically.
                 + a one-shot read of the frozen 128x256 P1 endpoint.
  w  LADDER      the same relax at wscale x {10, 100, 1000}: the
                 penalty-strength bridge to the constraint limit (the
                 physical regime: the P4 law stiff-mode ~ g^3 makes the
                 potential effectively hard at g ~ 1e10). Gate: does the
                 spread (r90 growth, dE tail) arrest with stiffness?
  i  ISOSPECTRAL exact orbit-class descent: project force AND velocity
                 onto the tangent space of {Lam M Lam^T : Lam in
                 SO(1,3)} per cell, so the eta-spectrum field of the
                 SEED (his 2-equal vortex + 3-equal center profile) is
                 an exact invariant. Runs: rotations-only (J12,J13,J23)
                 and full so(1,3) (+boosts). Gate GT1: tangent
                 directions are exactly V4-flat (cyclicity). Gate GT2:
                 invariant drift stays small (monitored + corrected by a
                 pin-to-seed retraction, itself complex-step gated,
                 GT0). Verdict: does a LOCALIZED stationary state exist
                 in the orbit class (projected residual falls, r90
                 bounded) where the unconstrained slide spreads?

Outputs: data/m5_21_1e_constraint.json, plots/m5_21_1e_constraint.png,
films plots/m5_21_1e_iso_film_{basic,thermal}.png (film standard).
Modes: gates | d | w | i | plots | all
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

from m5_17_energy import cell_weights, grid_coords                  # noqa: E402
from m5_20_2_a_eom import (ETA, G_T, H, WSCALE, c4_of,              # noqa: E402
                           grad_static_4, u_eta_density, v4_density)
from m5_21_1_c_4d import (FREE42, NR2, NZ2, PIN2, RHO42, W42,       # noqa: E402
                          core64, q_meridional64, seed64)
from m5_film import film_strip                                      # noqa: E402

DELTA = 0.3
SG = +1.0 * G_T                    # the M5.21.1 P1 working branch (s = +1)
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
PLOTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plots")
W2 = cell_weights(NR2, NZ2, H)
SNAPS = (0, 1000, 2000, 4000, 8000)


# ---------------- so(1,3) generators + tangent machinery ----------------
def so13_gens():
    J, K = [], []
    for (a, b) in ((1, 2), (1, 3), (2, 3)):
        W = np.zeros((4, 4))
        W[a, b], W[b, a] = -1.0, 1.0
        J.append(W)
    for i in (1, 2, 3):
        W = np.zeros((4, 4))
        W[0, i], W[i, 0] = 1.0, 1.0
        K.append(W)
    return J, K


def tangent_basis(M, gens):
    """t_k = W_k M + M W_k^T per cell -> (nr, nz, k, 4, 4)."""
    return np.stack([W @ M + M @ W.T for W in gens], axis=2)


def project_tangent(G, T):
    """Frobenius projection of G onto span{t_k} per cell."""
    nr, nz, k = T.shape[0], T.shape[1], T.shape[2]
    Tf = T.reshape(nr, nz, k, 16)
    Gf = G.reshape(nr, nz, 16)
    gram = np.einsum("...ka,...la->...kl", Tf, Tf)
    rhs = np.einsum("...ka,...a->...k", Tf, Gf)
    tr = np.einsum("...kk->...", gram)[..., None, None]
    reg = gram + (1e-12 * tr + 1e-300) * np.eye(k)
    coef = np.linalg.solve(reg, rhs[..., None])[..., 0]
    return np.einsum("...k,...ka->...a", coef, Tf).reshape(nr, nz, 4, 4)


def invariants(M):
    """I_p = tr((eta M)^p), p = 1..4 -> (nr, nz, 4)."""
    EM = np.broadcast_to(ETA, M.shape) @ M
    P, out = EM, []
    for p in range(4):
        out.append(np.einsum("...aa->...", P))
        P = P @ EM
    return np.stack(out, axis=-1)


def grad_pin(M, I_tgt, wp=1.0):
    """d/dM of P_pin = wp sum_p (I_p - I_tgt_p)^2 (symmetrized)."""
    EM = np.broadcast_to(ETA, M.shape) @ M
    powers = [np.broadcast_to(np.eye(4), M.shape), EM, EM @ EM]
    I = invariants(M)
    dv = np.zeros_like(M)
    for p in range(1, 5):
        coef = (2.0 * wp * (I[..., p - 1] - I_tgt[..., p - 1])
                * p)[..., None, None]
        dv = dv + coef * (powers[p - 1] @ ETA if p <= 3
                          else (powers[2] @ EM) @ ETA)
    return 0.5 * (dv + np.swapaxes(dv, -1, -2))


def p_pin_c(M, I_tgt, wp=1.0):
    """complex-safe P_pin (no float cast): used by the GT0 complex step."""
    I = invariants(M)
    return wp * np.sum((I - I_tgt) ** 2)


def p_pin(M, I_tgt, wp=1.0):
    return float(np.real(p_pin_c(M, I_tgt, wp)))


# ---------------- shared observables ----------------
def e_split(M, wscale=WSCALE):
    u = float(np.sum(u_eta_density(M, H) * W2))
    v = float(np.sum(v4_density(M, wscale, DELTA, g=SG) * W2))
    return u, v


def containment(M, wscale=WSCALE):
    """spherical cumulative-energy radii on the total static density."""
    dens = (u_eta_density(M, H) + v4_density(M, wscale, DELTA, g=SG)) * W2
    R, Z = grid_coords(NR2, NZ2, H)
    r = np.sqrt(R ** 2 + Z ** 2)[: NR2 - 1, 1:-1].ravel()
    d = dens.ravel()
    o = np.argsort(r)
    cs = np.cumsum(d[o])
    tot = max(cs[-1], 1e-300)
    r50 = float(r[o][np.searchsorted(cs, 0.5 * tot)])
    r90 = float(r[o][np.searchsorted(cs, 0.9 * tot)])
    ball = float(np.sum(d[r <= 8.0]) / tot)
    return r50, r90, ball


def diag_row(M, it, wscale=WSCALE, gens=None, I_ref=None):
    u, v = e_split(M, wscale)
    r50, r90, ball = containment(M, wscale)
    G = grad_static_4(M, wscale, DELTA, g=SG, w4=W42, rho4=RHO42) * FREE42
    row = {"it": it, "E_u": u, "E_v4": v, "E": u + v,
           "virial_u_over_3v": u / max(3.0 * v, 1e-300),
           "r50": r50, "r90": r90, "coreball": ball,
           "q": q_meridional64(M), **core64(M),
           "gnorm": float(np.sqrt(np.sum(G * G)))}
    if gens is not None:
        T = tangent_basis(M, gens)
        Gt = project_tangent(G, T)
        row["tan_frac"] = float(np.sum(Gt * Gt) / max(np.sum(G * G), 1e-300))
    if I_ref is not None:
        I = invariants(M)
        row["drift_max"] = float(np.max(np.abs(I - I_ref)
                                        / (1.0 + np.abs(I_ref))))
        row["drift_l2"] = float(np.sqrt(np.sum((I - I_ref) ** 2)))
    return row


# ---------------- gates ----------------
def gates():
    rng = np.random.default_rng(7)
    M = seed64(DELTA, SG) + 0.05 * _sym_noise(rng)
    J, K = so13_gens()
    out = {}
    # GT0: grad_pin vs complex step
    I_tgt = invariants(seed64(DELTA, SG))
    V = _sym_noise(rng)
    g_an = float(np.sum(grad_pin(M, I_tgt) * V))
    eps = 1e-30
    g_cs = float(np.imag(p_pin_c(M + 1j * eps * V, I_tgt)) / eps)
    out["GT0_pin_grad_rel"] = abs(g_an - g_cs) / max(abs(g_cs), 1e-300)
    # GT1: tangent V4-flatness (all 6 generators, complex step on V4)
    worst = 0.0
    for W in J + K:
        t = W @ M + M @ W.T
        v0 = float(np.sum(v4_density(M, WSCALE, DELTA, g=SG) * W2))
        v1 = float(np.imag(np.sum(v4_density(M + 1j * eps * t, WSCALE,
                                             DELTA, g=SG) * W2)) / eps)
        worst = max(worst, abs(v1) / max(abs(v0), 1e-300))
    out["GT1_tangent_v4_flat_rel"] = worst
    # GT2: projector idempotence + symmetry of tangent
    T = tangent_basis(M, J + K)
    G = grad_static_4(M, WSCALE, DELTA, g=SG, w4=W42, rho4=RHO42)
    Gt = project_tangent(G, T)
    Gtt = project_tangent(Gt, T)
    out["GT2_proj_idem_rel"] = float(np.sqrt(np.sum((Gtt - Gt) ** 2))
                                     / max(np.sqrt(np.sum(Gt * Gt)), 1e-300))
    out["GT2_tangent_sym"] = float(np.max(np.abs(T - np.swapaxes(T, -1, -2))))
    print(json.dumps(out, indent=1))
    ok = (out["GT0_pin_grad_rel"] < 1e-10
          and out["GT1_tangent_v4_flat_rel"] < 1e-10
          and out["GT2_proj_idem_rel"] < 1e-6)   # regularized projector
    out["PASS"] = bool(ok)
    return out


def _sym_noise(rng):
    A = rng.normal(size=(NR2, NZ2, 4, 4))
    return 0.5 * (A + np.swapaxes(A, -1, -2))


# ---------------- FIRE (unconstrained; fire64 clone with wscale knob) ----
def fire_w(M0, wscale, max_iter, dt0, dt_max, snaps=SNAPS, tag="",
           gens=None, I_ref=None, log_every=1000):
    w = np.ones((NR2, NZ2))
    w[: NR2 - 1, 1:-1] = W2
    precond = (1.0 / w)[..., None, None]
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    rows = [diag_row(M, 0, wscale, gens, I_ref)]
    states = [{"it": 0, "t": 0.0, "M": M0.copy()}]
    F = -grad_static_4(M, wscale, DELTA, g=SG, w4=W42, rho4=RHO42) \
        * precond * FREE42
    t0 = time.time()
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
        if not np.all(np.isfinite(M)):
            rows.append({"it": it, "stop": "non-finite"})
            break
        F = -grad_static_4(M, wscale, DELTA, g=SG, w4=W42, rho4=RHO42) \
            * precond * FREE42
        if it in snaps or it == max_iter:
            rows.append(diag_row(M, it, wscale, gens, I_ref))
            states.append({"it": it, "t": float(it), "M": M.copy()})
        if it % log_every == 0:
            print(f"  {tag} it {it:6d} E {rows[-1].get('E', float('nan')):12.6f}"
                  f" [{time.time() - t0:.0f}s]", flush=True)
    return M, rows, states


# ---------------- projected isospectral descent ----------------
def fire_iso(M0, gens, max_iter, dt0=0.02, dt_max=0.2, snaps=SNAPS,
             tag="iso", corr_every=200, corr_steps=20, log_every=1000):
    w = np.ones((NR2, NZ2))
    w[: NR2 - 1, 1:-1] = W2
    precond = (1.0 / w)[..., None, None]
    M = M0.copy()
    I_ref = invariants(M0)
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    rows = [diag_row(M, 0, WSCALE, gens, I_ref)]
    rows[0]["res_rel"] = 1.0
    states = [{"it": 0, "t": 0.0, "M": M0.copy()}]

    def proj_force(Mx):
        G = grad_static_4(Mx, WSCALE, DELTA, g=SG, w4=W42, rho4=RHO42)
        F0 = -G * precond * FREE42
        T = tangent_basis(Mx, gens)
        return project_tangent(F0, T), T, G

    F, T, G0 = proj_force(M)
    f_seed = np.sqrt(np.sum(F * F))
    t0 = time.time()
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
        v = project_tangent(v, T)
        M += dt * v
        if not np.all(np.isfinite(M)):
            rows.append({"it": it, "stop": "non-finite"})
            break
        if it % corr_every == 0:
            # pin-to-seed retraction (drift correction), logged honestly
            p0 = p_pin(M, I_ref)
            for _ in range(corr_steps):
                gp = grad_pin(M, I_ref) * FREE42
                gn = float(np.max(np.abs(gp)))
                if gn < 1e-14:
                    break
                step = min(0.2 / gn, 1.0)
                M2 = M - step * gp
                if p_pin(M2, I_ref) < p_pin(M, I_ref):
                    M = M2
                else:
                    break
            p1 = p_pin(M, I_ref)
            if it % log_every == 0:
                print(f"  {tag} retraction P_pin {p0:.3e} -> {p1:.3e}",
                      flush=True)
        F, T, G0 = proj_force(M)
        if it in snaps or it == max_iter:
            row = diag_row(M, it, WSCALE, gens, I_ref)
            row["res_rel"] = float(np.sqrt(np.sum(F * F)) / max(f_seed,
                                                                1e-300))
            rows.append(row)
            states.append({"it": it, "t": float(it), "M": M.copy()})
        if it % log_every == 0:
            print(f"  {tag} it {it:6d} E {rows[-1].get('E', float('nan')):12.6f}"
                  f" res {rows[-1].get('res_rel', float('nan')):.3e}"
                  f" drift {rows[-1].get('drift_max', float('nan')):.2e}"
                  f" [{time.time() - t0:.0f}s]", flush=True)
    return M, rows, states


# ---------------- phases ----------------
def _merge(update):
    path = os.path.join(DATA, "m5_21_1e_constraint.json")
    cur = {}
    if os.path.exists(path):
        with open(path) as f:
            cur = json.load(f)
    cur.update(update)
    with open(path, "w") as f:
        json.dump(cur, f, indent=1)
    print(f"wrote {path}", flush=True)


def phase_d(iters=8000):
    J, K = so13_gens()
    gens = J + K
    M0 = seed64(DELTA, SG)
    I_vac = invariants(M0 * 0 + np.diag([-SG, 1.0, DELTA, 0.0]))
    M, rows, _ = fire_w(M0, WSCALE, iters, 0.02, 0.2, tag="d-base",
                        gens=gens, I_ref=I_vac)
    out = {"baseline": rows}
    # one-shot on the frozen 128x256 P1 endpoint
    p = os.path.join(DATA, "m5_21_1_b_endpoint.npz")
    if os.path.exists(p):
        z = np.load(p)
        Mep = z[z.files[0]]
        w2 = cell_weights(Mep.shape[0], Mep.shape[1], H)
        u = float(np.sum(u_eta_density(Mep, H) * w2))
        v = float(np.sum(v4_density(Mep, WSCALE, DELTA, g=SG) * w2))
        G = grad_static_4(Mep, WSCALE, DELTA, g=SG)
        T = tangent_basis(Mep, gens)
        Gt = project_tangent(G, T)
        out["endpoint128"] = {
            "E_u": u, "E_v4": v, "virial_u_over_3v": u / max(3 * v, 1e-300),
            "tan_frac": float(np.sum(Gt * Gt) / max(np.sum(G * G), 1e-300))}
    _merge({"phase_d": out})


def phase_w(iters=8000):
    out = {}
    for mult in (10.0, 100.0, 1000.0):
        dt0 = 0.02 / np.sqrt(mult)
        M, rows, states = fire_w(seed64(DELTA, SG), WSCALE * mult, iters,
                                 dt0, 10 * dt0, tag=f"w-x{mult:g}")
        out[f"x{mult:g}"] = rows
        if mult == 1000.0:
            np.savez_compressed(os.path.join(
                DATA, "m5_21_1e_w1000_endpoint.npz"), M=M)
    _merge({"phase_w": out})


def phase_i(iters=8000, only=None):
    J, K = so13_gens()
    path = os.path.join(DATA, "m5_21_1e_constraint.json")
    out = {}
    if os.path.exists(path):
        with open(path) as f:
            out = json.load(f).get("phase_i", {})
    # rot3/full6: the orbit class of the STANDARD seed (r_c = 4 core
    # blend). smallcore: r_c = 1, testing whether the orientation escape
    # threads through the fat degenerate-core region (the hatch-size
    # hypothesis; added after the first rot3/full6 read).
    for tag, gens, rc in (("rot3", J, 4.0), ("full6", J + K, 4.0),
                          ("smallcore", J, 1.0)):
        if only is not None and tag not in only:
            continue
        M, rows, states = fire_iso(seed64(DELTA, SG, r_c=rc), gens, iters,
                                   tag=f"i-{tag}")
        out[tag] = rows
        _merge({"phase_i": out})     # per-run merge: a film crash cannot
        # lose a finished descent (lesson from the first phase_i attempt)
        if tag == "rot3":
            for tpl in ("basic", "thermal"):
                film_strip(states, os.path.join(
                    PLOTS, f"m5_21_1e_iso_film_{tpl}.png"), template=tpl,
                    delta=DELTA, g=SG)
            np.savez_compressed(os.path.join(
                DATA, "m5_21_1e_iso_endpoint.npz"), M=M)
    _merge({"phase_i": out})


def phase_plots():
    with open(os.path.join(DATA, "m5_21_1e_constraint.json")) as f:
        d = json.load(f)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    base = [r for r in d["phase_d"]["baseline"] if "E" in r]
    its = [r["it"] for r in base]
    ax[0, 0].plot(its, [r["E"] for r in base], "o-", label="wscale x1")
    for k, rows in d.get("phase_w", {}).items():
        rr = [r for r in rows if "E" in r]
        ax[0, 0].plot([r["it"] for r in rr], [r["E"] for r in rr], "o-",
                      label=f"wscale {k}")
    ax[0, 0].set_title("E(it): does stiffness arrest the slide?")
    ax[0, 0].set_yscale("log")
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].plot(its, [r["r90"] for r in base], "o-", label="x1 r90")
    for k, rows in d.get("phase_w", {}).items():
        rr = [r for r in rows if "r90" in r]
        ax[0, 1].plot([r["it"] for r in rr], [r["r90"] for r in rr], "o-",
                      label=f"{k} r90")
    for k, rows in d.get("phase_i", {}).items():
        rr = [r for r in rows if "r90" in r]
        ax[0, 1].plot([r["it"] for r in rr], [r["r90"] for r in rr], "s--",
                      label=f"iso {k} r90")
    ax[0, 1].set_title("r90(it): containment")
    ax[0, 1].legend(fontsize=7)
    ax[0, 2].plot(its, [r["virial_u_over_3v"] for r in base], "o-")
    ax[0, 2].axhline(1.0, color="k", lw=0.5)
    ax[0, 2].set_title("Derrick virial u/(3 V4), baseline")
    ax[1, 0].plot(its, [r.get("tan_frac", np.nan) for r in base], "o-")
    ax[1, 0].set_title("force fraction INSIDE the orbit tangent (baseline)")
    ax[1, 0].set_ylim(0, 1)
    ax[1, 1].plot(its, [r.get("drift_l2", np.nan) for r in base], "o-",
                  label="baseline vs vacuum")
    for k, rows in d.get("phase_i", {}).items():
        rr = [r for r in rows if "drift_l2" in r]
        ax[1, 1].plot([r["it"] for r in rr], [r["drift_l2"] for r in rr],
                      "s--", label=f"iso {k} vs seed")
    ax[1, 1].set_title("spectral drift (L2 of invariants)")
    ax[1, 1].set_yscale("log")
    ax[1, 1].legend(fontsize=7)
    for k, rows in d.get("phase_i", {}).items():
        rr = [r for r in rows if "res_rel" in r]
        ax[1, 2].plot([r["it"] for r in rr], [r["res_rel"] for r in rr],
                      "s--", label=f"iso {k}")
    ax[1, 2].set_title("isospectral projected residual (rel to seed)")
    ax[1, 2].set_yscale("log")
    ax[1, 2].legend(fontsize=7)
    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(PLOTS, "m5_21_1e_constraint.png")
    fig.savefig(p, dpi=110)
    print(f"wrote {p}", flush=True)


def main(which="all"):
    if which in ("gates", "all"):
        g = gates()
        _merge({"gates": g})
        if not g["PASS"]:
            print("GATES FAILED, stopping")
            return
    if which in ("d", "all"):
        phase_d()
    if which in ("w", "all"):
        phase_w()
    if which in ("i", "all"):
        phase_i()
    if which == "i2":
        phase_i(only=("smallcore",))
    if which in ("plots", "all"):
        phase_plots()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
