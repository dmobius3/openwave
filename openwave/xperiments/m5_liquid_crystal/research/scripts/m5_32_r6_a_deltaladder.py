"""M5.32 R6.a: the DELTA LADDER of the clock inertia on the R2 candidate
(the free rung triggered by R4: is the IR-extensive clock inertia a
toy-point artifact of delta = 0.3?).

EQUATIONS FIRST
---------------
Candidate (R2, every symbol as in m5_32_r2_b_bounded.py, imported as RB,
never modified; the R4 producer m5_32_r4_clock.py and its audit are the
reference instruments, read-only):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    I1 : q_eta(F) = tr(eta F eta F^T),   I1_h : q_h(F) = tr(h F h F^T)
    h(M) = eta + 2 (eta u)(eta u)^T,  u = timelike unit eigenvector of M eta
    E(M; omega) = E_stat_lambda(M) + omega^2 kin_lambda(M)     (exact in omega)
    E_stat = 4 h^3 sum_x sum_{i<j} q_lambda(F_ij) + V4
    kin    = 4 h^3 sum_x sum_i     q_lambda(comm_eta(a0, A_i))
The base object at every delta (the R4 base, m5_21_8 dressed(m = 0)):
    Mb = Qh d4 Qh^T,  d4 = diag(-s g, 1, delta, 0) = diag(g, 1, delta, 0) at s = -1
    Qh = R3(phi) R2(theta) : the twisting hedgehog (direction-only, V4 = 0)
    a0_unit = d/dt [Qh(t) d4 Qh(t)^T]_{t=0} = Qh (G1 d4 - d4 G1) Qh^T
    G1 = the (2, 3) rotation generator  =>  (G1 d4 - d4 G1) = delta (E_23 + E_32)
so the record clock generator is EXACTLY proportional to delta:
    a0_unit(delta) = delta * Qh (E_23 + E_32) Qh^T                    (identity A)
and on the whole boost-dressed family (M = Qb Mb Qb^T, a0 = Qb a0_unit Qb^T)
    kin(amp, R, L; delta) = delta^2 * kappa(amp, R, L; delta),
    kappa smooth at delta -> 0 (delta enters kappa only through the d4 entry
    in the jets A_i, which has a finite limit d4 -> diag(g, 1, 0, 0)).
Identity A is the analytic content of the R4 measurement kin ~ delta^2.03
at fixed texture; this rung measures whether the L-EXTENSIVE part and the
CORE part of kin scale differently in delta (the pre-registered hypothesis
kin(L) = kin_core + c delta^2 L^p needs kin_core to NOT scale as delta^2).
The T2 relaxation (existence probe): FIRE descent of the certified 4D
stack (m5_21_3_a_4d.py: sym stencil, eta-commutator curvature, V4 = the
trace-target eigenvalue penalty, pinned shell) from Mb at each delta, the
m5_21_8 relax budget (maxit 3000, dt0 0.02, dt_max 0.2); the hedgehog
EXISTS at delta if the descent stays finite, the far-field degree stays
-1 on every read surface, and E_end <= E_seed; it DEGENERATES if the
degree is lost or the field runs away.
Fixed-J electron (I1 frame, omega the multiplier):
    E_J(amp, R) = E_stat(amp, R) + J^2 / (4 kin(amp, R)),  omega* = J / (2 kin*)
    J = 2 kin_core omega_t,  kin_core = kin(L 48) frac(r < 12) (the inertia
    inside r < 12 on the n 32 box; the fitted intercept a of kin = a + b L^p
    is reported but is NEGATIVE on the R4 record: 426.5 / 658.3 / 890.4 at
    L 48 / 72 / 96 gives a = -37, so it is no positive core),
    omega_t in {0.1, 0.3};  J = 200 = the R4 record charge (control)
    min over amp on the grid AMPS (parabolic refinement), R* = argmin over
    the R grid (interior / wall test), the drift omega*(L 72) / omega*(L 48)
    - 1, the G3 gate |drift| <= 10 %.
Static dressing gain: G(R) = min_amp E_stat(amp, R) - E_stat(0) at lambda 1
(lambda 0 control), R in {6, 9, 12, 18, 24} (+ 36 from the L 72 box).
Radial reads: e_u(x), k(x) = the per-cell curvature and inertia densities
(h^3-weighted, branch-averaged); r_half = the radius enclosing half of
E_u; frac(r < 12) = the fraction of kin inside r < 12.
kin(L) fit: three boxes (n, L) = (32, 48), (48, 72), (64, 96) at h = 1.5,
kin = a + b L^p solved exactly (p by root of the difference ratio), with the
p = 1 least-squares fit as the second reading; b(delta) exponent = the
log-log slope over the ladder.

STAGES (python3 m5_32_r6_a_deltaladder.py delta <delta> [--maxit N] [--quick]
        | collect):
    delta   the whole per-delta program, partial -> ../checkpoints/m5_32_r6/
    collect merge partials -> ../data/m5_32_r6_deltaladder.json + plot
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r6")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RB = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
AUD22 = _load("m5_22_e_audit", "m5_22_e_audit.py")
B3, B8 = RB.B3, RB.B8
ETA = B3.ETA
G_MAIN, S_MAIN = RB.G_MAIN, RB.S_MAIN
DELTAS = (0.3, 0.1, 0.03, 0.01)
LAMBDAS = (0.0, 0.75, 1.0)
OMEGA_T = (0.1, 0.3)
J_RECORD = 200.0
BOXES = ((32, 48.0), (48, 72.0), (64, 96.0))
FIXJ_BOXES = ((32, 48.0), (48, 72.0))
R_LADDER = (6.0, 9.0, 12.0, 18.0)
AMPS = np.round(np.concatenate([np.arange(0.0, 0.0601, 0.0025),
                                [0.07, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3]]), 6)
DEGREE_WIDTHS = (6.0, 12.0, 18.0)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def ck(name):
    os.makedirs(CKPT, exist_ok=True)
    return os.path.join(CKPT, name)


def dump(name, obj):
    with open(ck(name), "w") as f:
        json.dump(obj, f, indent=1)


def cfg_of(n, L, delta):
    return B3.base_cfg(s=S_MAIN, g=G_MAIN, n=n, L=float(L), delta=float(delta))


# ================= per-cell densities =================
def densities(M, a0, cfg):
    """(e_u(x), k(x)): per-cell curvature and inertia densities, h^3-weighted,
    branch-averaged (the lambda-blind I1 reads; at amp = 0 I1 == I1_h)."""
    h3, h = cfg["h"] ** 3, cfg["h"]
    eu = np.zeros(M.shape[:3])
    kk = np.zeros(M.shape[:3])
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                eu += wt * 4.0 * B3.inner_eta(F, F)
            if a0 is not None:
                F = B3.comm_eta(a0, A[i])
                kk += wt * 4.0 * B3.inner_eta(F, F)
    return h3 * eu, h3 * kk


def radial_profile(dens, r, radii=(3.0, 6.0, 9.0, 12.0, 18.0)):
    tot = float(np.sum(dens))
    order = np.argsort(r.ravel())
    cum = np.cumsum(dens.ravel()[order])
    rs = r.ravel()[order]
    out = {"total": tot, "frac_inside": {}}
    for R in radii:
        out["frac_inside"][f"{R:g}"] = float(np.sum(dens[r < R]) / tot) if tot != 0 else None
    if tot > 0:
        out["r_half"] = float(rs[int(np.searchsorted(cum, 0.5 * tot))])
        out["r_90"] = float(rs[int(np.searchsorted(cum, 0.9 * tot))])
    return out


def degree_of(M, cfg, widths=DEGREE_WIDTHS):
    """the Q37 instrument (m5_22_e_audit.read_charge_from_M) on the spatial
    3x3 block, centered cube surfaces of half-width w."""
    n, h = cfg["n"], cfg["h"]
    out = {}
    M3 = np.ascontiguousarray(M[..., 1:4, 1:4])
    for w in widths:
        k = int(np.floor(w / h + (n - 1) / 2.0 + 1e-9))
        ilo, ihi = n - 1 - k, k
        if ilo < 0 or ihi > n - 1:
            continue        # surface outside the grid (small quick boxes)
        Q, cf = AUD22.read_charge_from_M(M3, ilo, ihi)
        out[f"{w:g}"] = {"Q": float(Q), "conflicts": int(cf), "ilo": ilo, "ihi": ihi}
    return out


def spectrum_report(M, cfg, delta):
    """the M eta spectrum: the (delta, 0) pair gap over the free cells and at
    the cell nearest the origin."""
    n = cfg["n"]
    lam = np.sort(np.linalg.eigvals(M @ ETA).real, axis=-1)   # ascending
    free = ~B3.pin_shell(n, cfg["h"])
    # M eta spectrum ascending: [s g, 0, delta, 1] = [-32, 0, delta, 1] at s = -1;
    # the (0, delta) pair sits at indices 1, 2
    gap_d0 = lam[..., 2] - lam[..., 1]
    c = n // 2
    return {"gap_delta0_min_free": float(np.min(gap_d0[free])),
            "gap_delta0_min_free_over_delta": float(np.min(gap_d0[free]) / delta),
            "gap_delta0_median_free": float(np.median(gap_d0[free])),
            "spectrum_near_origin": lam[c, c, c].tolist(),
            "spectrum_vacuum": sorted([S_MAIN * G_MAIN, 1.0, delta, 0.0])}


# ================= the dressing family =================
class Family:
    def __init__(self, n, L, delta):
        self.n, self.L, self.delta = n, float(L), float(delta)
        self.cfg = cfg_of(n, L, delta)
        self.R, self.K, self.K2 = RB.boost_geom(self.cfg)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0b = B8.a0_unit(self.cfg, 0.0)
        self.cache = {}
        self.n_reads = 0

    def b_of(self, amp, R, kind):
        r = self.R
        if kind == "rigid":
            return amp * np.tanh(r / 2.0)
        return amp * np.tanh(r / 2.0) * np.exp(-(r / R) ** 2)

    def read(self, amp, R, kind="exp2"):
        key = (round(float(amp), 10), round(float(R), 8), kind)
        if key not in self.cache:
            Md, a0d = RB.dress(self.Mb, self.a0b, self.b_of(*key), self.K, self.K2)
            rd = RB.reads(Md, self.cfg, a0d)
            rd.update(amp=key[0], R=key[1], kind=kind)
            self.cache[key] = {k: rd[k] for k in ("amp", "R", "kind", "A_I1", "A_I1h", "V4",
                                                  "kin_I1", "kin_I1h", "min_gap",
                                                  "timelike_eigvec_everywhere")}
            self.n_reads += 1
        return self.cache[key]


def es_kin(rd, lam):
    return RB.e_stat_lam(rd, lam), RB.kin_lam(rd, lam)


def min_over_amp(rows, J, lam):
    """rows: the amp grid at one R. E_J on the grid, parabolic refinement at
    an interior argmin; kin <= 0 excluded (the lambda 0 guard edge)."""
    amps = np.array([r["amp"] for r in rows])
    es = np.array([es_kin(r, lam)[0] for r in rows])
    kk = np.array([es_kin(r, lam)[1] for r in rows])
    ok = kk > 0
    E = np.where(ok, es + J * J / (4.0 * np.where(ok, kk, 1.0)), np.inf)
    i = int(np.argmin(E))
    if not np.isfinite(E[i]):
        return None
    a_star, E_star = amps[i], E[i]
    interior = 0 < i < len(rows) - 1 and np.isfinite(E[i - 1]) and np.isfinite(E[i + 1])
    refined = False
    if interior:
        a, b, c = E[i - 1], E[i], E[i + 1]
        da = amps[i + 1] - amps[i]
        if c - 2 * b + a > 1e-14 and abs(amps[i] - amps[i - 1] - da) < 1e-12:
            a_star = amps[i] - 0.5 * da * (c - a) / (c - 2 * b + a)
            E_star = b - 0.125 * (c - a) ** 2 / (c - 2 * b + a)
            refined = True
    k_star = float(np.interp(a_star, amps, kk))
    es_star = float(np.interp(a_star, amps, es))
    edge = "interior" if interior else ("amp_min" if i == 0 else "amp_max")
    if ok.sum() and i == int(np.where(ok)[0][-1]) and not ok.all():
        edge = "kin_guard_edge"
    return {"amp_star": float(a_star), "E_J": float(E_star), "E_stat": es_star, "kin": k_star,
            "J2_over_4kin": float(J * J / (4.0 * k_star)), "omega_star": float(J / (2.0 * k_star)),
            "amp_edge": edge, "parabolic_refined": refined, "n_kin_positive": int(ok.sum())}


def interior_test(Rs, EJ):
    E = np.array([np.inf if e is None else e for e in EJ])
    i = int(np.argmin(E))
    fin = E[np.isfinite(E)]
    return {"R_star": float(Rs[i]), "argmin_index": i, "interior": bool(0 < i < len(Rs) - 1),
            "at_wall": bool(i == len(Rs) - 1),
            "monotone_decreasing_in_R": bool(len(fin) > 1 and np.all(np.diff(fin) < 0)),
            "E_J_of_R": [None if not np.isfinite(e) else float(e) for e in E]}


# ================= fits =================
def fit_aLp(Ls, ks):
    """kin = a + b L^p through three points (exact): p from the difference ratio."""
    from scipy.optimize import brentq
    L1, L2, L3 = Ls
    k1, k2, k3 = ks
    out = {"L": list(Ls), "kin": list(ks)}
    d21, d32 = k2 - k1, k3 - k2
    ok = d21 != 0 and (d32 / d21) > 0
    if ok:
        ratio = d32 / d21

        def f(p):
            return (L3 ** p - L2 ** p) / (L2 ** p - L1 ** p) - ratio
        try:
            p = brentq(f, 0.02, 6.0)
            b = d21 / (L2 ** p - L1 ** p)
            a = k1 - b * L1 ** p
            out.update(p=float(p), b=float(b), a=float(a), exact3=True)
        except ValueError:
            ok = False
    if not ok:
        out.update(p=None, b=None, a=None, exact3=False,
                   note="difference ratio outside the reachable range (non-monotone or flat)")
    # the p = 1 least-squares reading
    A = np.vstack([np.ones(3), np.array(Ls)]).T
    co, *_ = np.linalg.lstsq(A, np.array(ks), rcond=None)
    res = np.array(ks) - A @ co
    out["linear"] = {"a": float(co[0]), "b": float(co[1]),
                     "max_abs_residual": float(np.max(np.abs(res))),
                     "max_rel_residual": float(np.max(np.abs(res)) / np.max(np.abs(ks)))}
    return out


# ================= the relaxation =================
def backtracking_fire(M, cfg, free_mask, it_cap=3000, dt0=0.02, dt_max=0.2,
                      f_tol=1e-8, log_every=250, tag=""):
    """energy-monotone FIRE on E = B3.e_total (the certified action, V4 from
    cfg): a step that raises E is rejected and dt halves (the R4 stage_relax
    control); returns (M, info)."""
    free = free_mask[..., None, None].astype(float)
    E_prev = float(B3.e_total(M, cfg))
    F = -B3.grad(M, cfg) * free
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    n_acc = n_rej = it = 0
    stop = "it_cap"
    trace = [{"it": 0, "acc": 0, "E": E_prev, "fmax": float(np.max(np.abs(F))), "dt": dt}]
    while it < it_cap:
        it += 1
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max); alpha *= 0.99
        else:
            v[:] = 0.0; alpha = 0.1; n_up = 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E = float(B3.e_total(M_try, cfg))
        if not np.isfinite(E) or E > E_prev + 1e-12 * max(abs(E_prev), 1.0):
            n_rej += 1; dt *= 0.5; v[:] = 0.0; alpha, n_up = 0.1, 0
            if dt < 1e-7:
                stop = "dt_collapsed"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -B3.grad(M, cfg) * free
        fmax = float(np.max(np.abs(F)))
        if n_acc % log_every == 0:
            trace.append({"it": it, "acc": n_acc, "E": E, "fmax": fmax, "dt": dt, "rej": n_rej})
            log(f"{tag} RELAX acc {n_acc:5d} it {it:5d} E {E:.6f} fmax {fmax:.3e} dt {dt:.2e} rej {n_rej}")
        if fmax < f_tol:
            stop = "f_tol"
            break
    trace.append({"it": it, "acc": n_acc, "E": E_prev, "fmax": float(np.max(np.abs(F))), "dt": dt, "rej": n_rej})
    return M, {"stop": stop, "trace": trace, "accepted": n_acc, "rejected": n_rej, "iterations": it}


# ================= the per-delta program =================
def stage_delta(delta, maxit=3000, quick=False):
    boxes = ((16, 24.0), (24, 36.0), (32, 48.0)) if quick else BOXES
    fboxes = ((16, 24.0), (24, 36.0)) if quick else FIXJ_BOXES
    amps = AMPS[::4] if quick else AMPS
    tag = f"d{delta:g}"
    out = {"delta": delta, "g": G_MAIN, "s": S_MAIN, "stencil": "sym",
           "identity_A": "a0_unit = delta * Qh (E_23 + E_32) Qh^T exactly (G1 d4 - d4 G1 = delta (E_23 + E_32))",
           "amps": amps.tolist(), "omega_t": list(OMEGA_T), "J_record": J_RECORD,
           "hedgehog": {}, "kin_ladder": {}, "fixedj": {}, "static_gain": {}, "rigid": {}}
    os.makedirs(CKPT, exist_ok=True)

    # ---------- 1. the hedgehog at delta (n 32, L 48) + the T2 relaxation ----------
    n, L = (16, 24.0) if quick else (32, 48.0)
    cfg = cfg_of(n, L, delta)
    Mb = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    X, Y, Z = B3.coords(n, cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    eu, ev = B3.e_parts(Mb, cfg)
    deu, dk = densities(Mb, a0, cfg)
    a0_over_delta = float(np.max(np.abs(a0)) / delta)
    hh = {"n": n, "L": L, "h": cfg["h"], "delta": delta,
          "ansatz": {"E_stat": float(eu + ev), "E_u": float(eu), "V4": float(ev),
                     "kin": float(np.sum(dk)), "kin_over_delta2": float(np.sum(dk) / delta ** 2),
                     "a0_max_over_delta": a0_over_delta,
                     "E_u_radial": radial_profile(deu, r), "kin_radial": radial_profile(dk, r),
                     "degree": degree_of(Mb, cfg), "spectrum": spectrum_report(Mb, cfg, delta)}}
    log(f"{tag} HEDGEHOG ansatz (n {n}, L {L:g}, h {cfg['h']:g}): E_stat {eu + ev:.5f} kin {np.sum(dk):.5f} "
        f"kin/d^2 {np.sum(dk) / delta ** 2:.3f} r_half(E_u) {hh['ansatz']['E_u_radial']['r_half']:.2f} "
        f"kin frac<12 {hh['ansatz']['kin_radial']['frac_inside']['12']:.4f} "
        f"degree {[v['Q'] for v in hh['ansatz']['degree'].values()]}")
    out["hedgehog"] = hh
    dump(f"delta_{tag}.json", out)
    # the T2 relaxation (energy-monotone backtracking FIRE on the certified
    # action B3.e_total / B3.grad, both delta-consistent through cfg; the
    # m5_21_8 relax iteration budget, the R4 stage_relax step control)
    free = ~B3.pin_shell(n, cfg["h"])
    t1 = time.time()
    Mr, info = backtracking_fire(Mb.copy(), cfg, free, it_cap=maxit, tag=tag)
    finite = bool(np.all(np.isfinite(Mr)))
    rel = {"maxit": maxit, "stop": info["stop"], "wall_s": round(time.time() - t1, 1),
           "trace": info["trace"], "finite": finite, "accepted": info["accepted"],
           "rejected": info["rejected"], "iterations": info["iterations"],
           "method": "backtracking FIRE (energy-monotone, dt halves on an uphill "
                     "step; stop on dt < 1e-7, f_tol 1e-8, or the iteration cap) on "
                     "E = E_u + V4 of the certified 4D stack, pinned shell depth 1.6"}
    if finite:
        eu1, ev1 = B3.e_parts(Mr, cfg)
        deu1, dk1 = densities(Mr, a0, cfg)
        dv = float(np.sqrt(np.sum((Mr - Mb) ** 2)) / max(np.sqrt(np.sum((Mb - B3.vac4(cfg)) ** 2)), 1e-300))
        deg = degree_of(Mr, cfg)
        rel.update({"E_stat": float(eu1 + ev1), "E_u": float(eu1), "V4": float(ev1),
                    "E_seed": float(eu + ev), "rel_move": dv,
                    "kin_frozen_a0": float(np.sum(dk1)), "kin_frozen_a0_over_delta2": float(np.sum(dk1) / delta ** 2),
                    "E_u_radial": radial_profile(deu1, r), "kin_radial": radial_profile(dk1, r),
                    "degree": deg, "spectrum": spectrum_report(Mr, cfg, delta),
                    "max_abs_M0i": float(np.max(np.abs(Mr[..., 0, 1:4])))})
        degrees_ok = all(abs(abs(v["Q"]) - 1.0) < 1e-6 for v in deg.values())
        rel["exists"] = bool(degrees_ok and eu1 + ev1 <= eu + ev + 1e-9)
        rel["degenerates"] = not rel["exists"]
        np.savez_compressed(ck(f"relaxed_{tag}_n{n}_L{L:g}.npz"), M=Mr.astype(np.float32),
                            delta=delta, h=cfg["h"])
    else:
        rel.update(exists=False, degenerates=True)
    out["hedgehog"]["relaxed"] = rel
    log(f"{tag} RELAX {info['stop']} in {rel['wall_s']} s: " + (
        f"E {rel['E_stat']:.5f} (seed {rel['E_seed']:.5f}) move {rel['rel_move']:.3e} "
        f"kin(frozen a0) {rel['kin_frozen_a0']:.4f} degree {[v['Q'] for v in rel['degree'].values()]} "
        f"gap(delta,0)/delta {rel['spectrum']['gap_delta0_min_free_over_delta']:.3f} exists {rel['exists']}"
        if finite else "NON-FINITE"))
    dump(f"delta_{tag}.json", out)

    # ---------- 2. the undressed clock inertia across the L ladder ----------
    kl = {"boxes": [], "note": "amp = 0: lambda-blind (F_0i has no time row on the undressed hedgehog)"}
    for (nn, LL) in boxes:
        c = cfg_of(nn, LL, delta)
        Mx = B8.dressed(c, 0.0)
        ax = B8.a0_unit(c, 0.0)
        Xx, Yx, Zx = B3.coords(nn, c["h"])
        rx = np.sqrt(Xx * Xx + Yx * Yx + Zx * Zx)
        deux, dkx = densities(Mx, ax, c)
        eux, evx = B3.e_parts(Mx, c)
        row = {"n": nn, "L": LL, "h": c["h"], "E_stat": float(eux + evx), "kin": float(np.sum(dkx)),
               "kin_over_delta2": float(np.sum(dkx) / delta ** 2),
               "kin_radial": radial_profile(dkx, rx), "E_u_radial": radial_profile(deux, rx)}
        kl["boxes"].append(row)
        log(f"{tag} KIN (n {nn}, L {LL:g}, h {c['h']:g}): E_stat {row['E_stat']:.5f} kin {row['kin']:.5f} "
            f"kin/d^2 {row['kin_over_delta2']:.3f} frac<12 {row['kin_radial']['frac_inside']['12']:.4f}")
        del Mx, ax, deux, dkx
    Ls = [b["L"] for b in kl["boxes"]]
    ks = [b["kin"] for b in kl["boxes"]]
    kl["fit"] = fit_aLp(Ls, ks)
    kl["fit_over_delta2"] = fit_aLp(Ls, [k / delta ** 2 for k in ks])
    fit = kl["fit"]
    # the operational core: the inertia inside r < 12 on the smallest box (the
    # fitted intercept a is reported but is not a positive core at delta 0.3:
    # kin(48/72/96) = 426.5/658.3/890.4 gives a = -37 on the R4 record)
    kin_core = float(ks[0] * kl["boxes"][0]["kin_radial"]["frac_inside"]["12"])
    kl["kin_core"] = kin_core
    kl["kin_core_definition"] = "kin(L 48) * frac(r < 12): the inertia inside r < 12 on the n 32, L 48 box"
    kl["kin_core_fit_intercept"] = fit["a"] if fit["a"] is not None else fit["linear"]["a"]
    out["kin_ladder"] = kl
    log(f"{tag} FIT kin = a + b L^p: a {fit['a']} b {fit['b']} p {fit['p']}; linear a {fit['linear']['a']:.3f} "
        f"b {fit['linear']['b']:.4f} maxrel {fit['linear']['max_rel_residual']:.2e}; kin_core {kin_core:.4f}")
    dump(f"delta_{tag}.json", out)

    # ---------- 3. + 4. + 5. the family grids ----------
    Js = {f"Jcore_{w:g}": 2.0 * kin_core * w for w in OMEGA_T}
    Js["J_200"] = J_RECORD
    out["J_values"] = Js
    for (nn, LL) in fboxes:
        fam = Family(nn, LL, delta)
        Rs = sorted(set(list(R_LADDER) + [LL / 2.0]))
        btag = f"n{nn}_L{LL:g}"
        grid = {}
        for R in Rs:
            grid[R] = [fam.read(a, R, "exp2") for a in amps]
            log(f"{tag} GRID {btag} R {R:g} done ({fam.n_reads} reads): E1(amp .02) "
                f"{es_kin(fam.read(0.02, R), 1.0)[0]:.4f} kin1 {es_kin(fam.read(0.02, R), 1.0)[1]:.4f}")
        fx = {"n": nn, "L": LL, "h": fam.cfg["h"], "R_grid": Rs,
              "grid": {f"R_{R:g}": grid[R] for R in Rs}, "scan": {}}
        for lam in LAMBDAS:
            for Jk, J in Js.items():
                byR = {f"R_{R:g}": min_over_amp(grid[R], J, lam) for R in Rs}
                it = interior_test(Rs, [byR[f"R_{R:g}"]["E_J"] if byR[f"R_{R:g}"] else None for R in Rs])
                star = byR[f"R_{it['R_star']:g}"]
                fx["scan"][f"lam_{lam:g}_{Jk}"] = {"lam": lam, "J": J, "byR": byR, "R_test": it,
                                                   "star": star}
                if lam > 0 or Jk == "J_200":
                    log(f"{tag} FIXJ {btag} lam {lam:g} {Jk} (J {J:.4g}): R* {it['R_star']:g} "
                        f"({'interior' if it['interior'] else 'WALL' if it['at_wall'] else 'R_min'}) "
                        + (f"amp* {star['amp_star']:.4f} w* {star['omega_star']:.5f} E {star['E_J']:.4f} "
                           f"kin* {star['kin']:.4f} edge {star['amp_edge']}" if star else "UNBOUNDED"))
        # the static dressing gain (J = 0)
        sg = {}
        for lam in LAMBDAS:
            byR = {}
            for R in Rs:
                es = np.array([es_kin(rd, lam)[0] for rd in grid[R]])
                i = int(np.argmin(es))
                byR[f"R_{R:g}"] = {"amp_star": float(amps[i]), "E_stat": float(es[i]),
                                   "E_stat_amp0": float(es[0]), "gain": float(es[i] - es[0]),
                                   "amp_edge": "interior" if 0 < i < len(es) - 1 else ("amp_min" if i == 0 else "amp_max")}
            g6 = byR["R_6"]["gain"]
            for R in Rs:
                byR[f"R_{R:g}"]["gain_over_G6"] = float(byR[f"R_{R:g}"]["gain"] / g6) if g6 != 0 else None
            gains = [byR[f"R_{R:g}"]["gain"] for R in Rs]
            incr = [gains[i + 1] - gains[i] for i in range(len(gains) - 1)]
            sg[f"lam_{lam:g}"] = {"byR": byR, "gains": gains, "increments": incr,
                                  "saturating": bool(all(abs(incr[i + 1]) < abs(incr[i]) for i in range(len(incr) - 1))),
                                  "last_increment_over_G6": float(incr[-1] / g6) if g6 != 0 else None}
            log(f"{tag} GAIN {btag} lam {lam:g}: G(R) " + " ".join(f"{g:.3f}" for g in gains)
                + " ; G/G6 " + " ".join(f"{g / g6:.3f}" if g6 != 0 else "-" for g in gains))
        fx["static_gain"] = sg
        out["fixedj"][btag] = fx
        out["static_gain"][btag] = sg
        dump(f"delta_{tag}.json", out)
        del fam, grid
    # the rigid family (the R2 / R4 control) on the three boxes
    for (nn, LL) in boxes:
        fam = Family(nn, LL, delta)
        rows = [fam.read(a, 1e9, "rigid") for a in amps]
        btag = f"n{nn}_L{LL:g}"
        rg = {"n": nn, "L": LL, "h": fam.cfg["h"], "rows": rows, "scan": {}}
        for lam in LAMBDAS:
            for Jk, J in Js.items():
                rg["scan"][f"lam_{lam:g}_{Jk}"] = min_over_amp(rows, J, lam)
        out["rigid"][btag] = rg
        s = rg["scan"]["lam_1_J_200"]
        log(f"{tag} RIGID {btag} lam 1 J 200: " + (f"amp* {s['amp_star']:.4f} w* {s['omega_star']:.5f} "
                                                  f"kin* {s['kin']:.4f} E {s['E_J']:.4f}" if s else "UNBOUNDED"))
        dump(f"delta_{tag}.json", out)
        del fam, rows
    # the drift gate per (lam, J) between the two fixed-J boxes
    drift = {}
    b1, b2 = [f"n{nn}_L{LL:g}" for (nn, LL) in fboxes]
    for key in out["fixedj"][b1]["scan"]:
        s1, s2 = out["fixedj"][b1]["scan"][key]["star"], out["fixedj"][b2]["scan"][key]["star"]
        if s1 and s2:
            d = (s2["omega_star"] - s1["omega_star"]) / s1["omega_star"]
            r12 = (out["fixedj"][b2]["scan"][key]["byR"]["R_12"]["omega_star"]
                   / out["fixedj"][b1]["scan"][key]["byR"]["R_12"]["omega_star"] - 1.0
                   if out["fixedj"][b1]["scan"][key]["byR"]["R_12"] and out["fixedj"][b2]["scan"][key]["byR"]["R_12"] else None)
            drift[key] = {"omega_L1": s1["omega_star"], "omega_L2": s2["omega_star"], "drift": float(d),
                          "gate_le_10pct": bool(abs(d) <= 0.10),
                          "R_star_L1": out["fixedj"][b1]["scan"][key]["R_test"]["R_star"],
                          "R_star_L2": out["fixedj"][b2]["scan"][key]["R_test"]["R_star"],
                          "drift_at_R12": r12,
                          "E_L1": s1["E_J"], "E_L2": s2["E_J"], "kin_L1": s1["kin"], "kin_L2": s2["kin"]}
    rig = {}
    tags = [f"n{nn}_L{LL:g}" for (nn, LL) in boxes]
    for key in out["rigid"][tags[0]]["scan"]:
        ws = [out["rigid"][t]["scan"][key]["omega_star"] if out["rigid"][t]["scan"][key] else None for t in tags]
        rig[key] = {"omega_star": ws, "boxes": tags,
                    "omega_times_L": [w * out["rigid"][t]["L"] if w else None for w, t in zip(ws, tags)]}
    out["drift_gate"] = drift
    out["rigid_drift"] = rig
    out["wall_s"] = round(time.time() - T0, 1)
    dump(f"delta_{tag}.json", out)
    log(f"{tag} DONE in {out['wall_s']} s")
    return out


# ================= collect =================
def first_box(p):
    """the smallest fixed-J box tag of a per-delta partial (n32_L48 in production)."""
    return sorted(p["fixedj"], key=lambda t: int(t[1:].split("_")[0]))[0]


def stage_collect():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    parts = {}
    for p in sorted(glob.glob(ck("delta_d*.json"))):
        d = json.load(open(p))
        parts[f"{d['delta']:g}"] = d
    deltas = sorted((float(k) for k in parts), reverse=True)
    out = {"task": "M5.32 R6.a: the delta ladder of the clock inertia (free rung triggered by R4)",
           "candidate": "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4, h = eta + 2(eta u)(eta u)^T",
           "point": {"g": G_MAIN, "s": S_MAIN, "deltas": deltas, "stencil": "sym"},
           "hypothesis": "kin(L) = kin_core + c delta^2 L^p with kin_core NOT ~ delta^2: the clock localizes as delta -> 0",
           "identity_A": parts[f"{deltas[0]:g}"]["identity_A"],
           "per_delta": parts}
    # hedgehog table
    out["hedgehog_table"] = {}
    for d in deltas:
        h = parts[f"{d:g}"]["hedgehog"]
        r = h.get("relaxed", {})
        out["hedgehog_table"][f"{d:g}"] = {
            "n": h["n"], "L": h["L"], "h": h["h"],
            "E_stat_ansatz": h["ansatz"]["E_stat"], "r_half_ansatz": h["ansatz"]["E_u_radial"]["r_half"],
            "degree_ansatz": {k: v["Q"] for k, v in h["ansatz"]["degree"].items()},
            "kin_ansatz": h["ansatz"]["kin"], "kin_frac_r12_ansatz": h["ansatz"]["kin_radial"]["frac_inside"]["12"],
            "E_stat_relaxed": r.get("E_stat"), "r_half_relaxed": (r.get("E_u_radial") or {}).get("r_half"),
            "degree_relaxed": {k: v["Q"] for k, v in r.get("degree", {}).items()},
            "rel_move": r.get("rel_move"), "stop": r.get("stop"), "exists": r.get("exists"),
            "kin_relaxed_frozen_a0": r.get("kin_frozen_a0"),
            "gap_delta0_over_delta_relaxed": (r.get("spectrum") or {}).get("gap_delta0_min_free_over_delta")}
    # kin(L) fits and the b(delta) exponent
    kt = {}
    for d in deltas:
        kl = parts[f"{d:g}"]["kin_ladder"]
        kt[f"{d:g}"] = {"L": [b["L"] for b in kl["boxes"]], "n": [b["n"] for b in kl["boxes"]],
                        "kin": [b["kin"] for b in kl["boxes"]],
                        "kin_over_delta2": [b["kin_over_delta2"] for b in kl["boxes"]],
                        "frac_r12": [b["kin_radial"]["frac_inside"]["12"] for b in kl["boxes"]],
                        "frac_r6": [b["kin_radial"]["frac_inside"]["6"] for b in kl["boxes"]],
                        "r_half_kin": [b["kin_radial"]["r_half"] for b in kl["boxes"]],
                        "fit": kl["fit"], "kin_core": kl["kin_core"]}
    out["kin_table"] = kt

    def slope(xs, ys):
        xs, ys = np.log(np.array(xs)), np.log(np.array(ys))
        A = np.vstack([np.ones_like(xs), xs]).T
        co, *_ = np.linalg.lstsq(A, ys, rcond=None)
        return float(co[1])
    ok = [d for d in deltas if kt[f"{d:g}"]["fit"]["b"] is not None and kt[f"{d:g}"]["fit"]["b"] > 0]
    exps = {}
    if len(ok) >= 2:
        exps["b_exact3_exponent"] = slope(ok, [kt[f"{d:g}"]["fit"]["b"] for d in ok])
        exps["a_exact3_exponent"] = (slope(ok, [kt[f"{d:g}"]["fit"]["a"] for d in ok])
                                     if all(kt[f"{d:g}"]["fit"]["a"] > 0 for d in ok) else None)
    exps["b_linear_exponent"] = slope(deltas, [kt[f"{d:g}"]["fit"]["linear"]["b"] for d in deltas])
    exps["a_linear_exponent"] = (slope(deltas, [kt[f"{d:g}"]["fit"]["linear"]["a"] for d in deltas])
                                 if all(kt[f"{d:g}"]["fit"]["linear"]["a"] > 0 for d in deltas) else None)
    exps["kin_L48_exponent"] = slope(deltas, [kt[f"{d:g}"]["kin"][0] for d in deltas])
    exps["kin_L96_exponent"] = slope(deltas, [kt[f"{d:g}"]["kin"][-1] for d in deltas])
    exps["pairwise_b_linear"] = {f"{deltas[i]:g}->{deltas[i + 1]:g}": slope(deltas[i:i + 2], [kt[f"{d:g}"]["fit"]["linear"]["b"] for d in deltas[i:i + 2]])
                                 for i in range(len(deltas) - 1)}
    out["delta_exponents"] = exps
    # the fixed-J drift table and the gain table
    out["drift_table"] = {f"{d:g}": parts[f"{d:g}"]["drift_gate"] for d in deltas}
    out["rigid_table"] = {f"{d:g}": parts[f"{d:g}"]["rigid_drift"] for d in deltas}
    out["gain_table"] = {f"{d:g}": {bt: {lk: {"gains": v["gains"], "G_over_G6": [b["gain_over_G6"] for b in v["byR"].values()],
                                             "saturating": v["saturating"], "last_increment_over_G6": v["last_increment_over_G6"]}
                                        for lk, v in sg.items()}
                                   for bt, sg in parts[f"{d:g}"]["static_gain"].items()}
                         for d in deltas}
    # the gate verdict
    gate = {}
    for d in deltas:
        dg = parts[f"{d:g}"]["drift_gate"]
        rows = {k: v for k, v in dg.items() if not k.startswith("lam_0_")}
        gate[f"{d:g}"] = {"n_cases": len(rows), "n_pass": int(sum(v["gate_le_10pct"] for v in rows.values())),
                          "n_wall_L48": int(sum(v["R_star_L1"] >= 24.0 - 1e-9 for v in rows.values())),
                          "drifts": {k: v["drift"] for k, v in rows.items()}}
    out["G3_gate_per_delta"] = gate
    # the R4 reproduction control at delta 0.3
    ctrl = {}
    if "0.3" in parts:
        p = parts["0.3"]
        bt0 = first_box(p)
        ref = {"rigid_omega_J200_lam1_R4audit": [0.1484154170952781, 0.09816222004771805, 0.0733],
               "localized_wall_omega_J200_lam1_R4audit": {"n32_L48": 0.18798143340570123, "n48_L72": 0.12230730035311002},
               "static_gain_lam1_R4audit_n32_L48": [-0.741, -2.387, -3.716, -5.424, -6.672],
               "base_E_stat_kin_R4": [62.851744331478166, 426.5070121483972]}
        try:
            a = json.load(open(os.path.join(DATA, "m5_32_r4_audit_clock.json")))
            ref["rigid_omega_J200_lam1_R4audit"] = [a["K3_rigid_R2_audit"][t]["J_200"]["omega_star"] for t in ("n32_L48", "n48_L72", "n64_L96")]
        except Exception:
            pass
        got = {"rigid_omega_J200_lam1": p["rigid_drift"]["lam_1_J_200"]["omega_star"],
               "localized_wall_omega_J200_lam1": {bt: p["fixedj"][bt]["scan"]["lam_1_J_200"]["star"]["omega_star"] for bt in p["fixedj"]},
               "static_gain_lam1_n32_L48": p["static_gain"][bt0]["lam_1"]["gains"],
               "base_E_stat_kin": [p["hedgehog"]["ansatz"]["E_stat"], p["hedgehog"]["ansatz"]["kin"]]}
        ctrl = {"reference": ref, "measured": got,
                "rel_dev_rigid": [abs(x / y - 1) for x, y in zip(got["rigid_omega_J200_lam1"], ref["rigid_omega_J200_lam1_R4audit"]) if x],
                "rel_dev_wall": {bt: abs(got["localized_wall_omega_J200_lam1"][bt] / ref["localized_wall_omega_J200_lam1_R4audit"][bt] - 1)
                                 for bt in ref["localized_wall_omega_J200_lam1_R4audit"] if bt in got["localized_wall_omega_J200_lam1"]},
                "rel_dev_base": [abs(x / y - 1) for x, y in zip(got["base_E_stat_kin"], ref["base_E_stat_kin_R4"])]}
    out["control_delta_0.3_vs_R4"] = ctrl
    # -------- plot --------
    os.makedirs(PLOTS, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(deltas)))
    for c, d in zip(cols, deltas):
        k = kt[f"{d:g}"]
        ax[0, 0].plot(k["L"], k["kin_over_delta2"], "o-", color=c, label=f"delta {d:g}")
        ax[0, 2].plot(k["L"], k["frac_r12"], "o-", color=c, label=f"delta {d:g} (r<12)")
        ax[0, 2].plot(k["L"], k["frac_r6"], "s--", color=c, label=f"delta {d:g} (r<6)")
    ax[0, 0].set_xlabel("L (h = 1.5)"); ax[0, 0].set_ylabel("kin / delta^2 (undressed)")
    ax[0, 0].set_title("undressed clock inertia across the L ladder"); ax[0, 0].legend(fontsize=7)
    ax[0, 1].loglog(deltas, [kt[f"{d:g}"]["fit"]["linear"]["b"] for d in deltas], "o-", label=f"b (p = 1 fit), slope {exps['b_linear_exponent']:.3f}")
    if all(kt[f"{d:g}"]["fit"]["linear"]["a"] > 0 for d in deltas):
        ax[0, 1].loglog(deltas, [kt[f"{d:g}"]["fit"]["linear"]["a"] for d in deltas], "s-", label=f"a = kin_core (p = 1 fit), slope {exps['a_linear_exponent']:.3f}")
    ax[0, 1].loglog(deltas, [kt[f"{d:g}"]["kin"][0] for d in deltas], "^:", label=f"kin(L 48), slope {exps['kin_L48_exponent']:.3f}")
    ax[0, 1].set_xlabel("delta"); ax[0, 1].set_title("delta scaling of the core and the L-extensive parts"); ax[0, 1].legend(fontsize=7)
    ax[0, 2].set_xlabel("L"); ax[0, 2].set_ylabel("fraction of kin inside r"); ax[0, 2].set_title("radial localization of kin (undressed)"); ax[0, 2].legend(fontsize=6)
    for c, d in zip(cols, deltas):
        bt0 = first_box(parts[f"{d:g}"])
        dg = parts[f"{d:g}"]["drift_gate"]
        for key, mk in (("lam_1_Jcore_0.1", "o-"), ("lam_1_Jcore_0.3", "s--"), ("lam_1_J_200", "^:")):
            if key in dg:
                v = dg[key]
                ax[1, 0].plot([48, 72], [v["omega_L1"], v["omega_L2"]], mk, color=c, label=f"delta {d:g} {key[6:]}")
        for lk, mk in (("lam_1", "o-"), ("lam_0", "x--")):
            sg = parts[f"{d:g}"]["static_gain"][bt0][lk]
            Rs = parts[f"{d:g}"]["fixedj"][bt0]["R_grid"]
            ax[1, 1].plot(Rs, [b["gain_over_G6"] for b in sg["byR"].values()], mk, color=c, label=f"delta {d:g} {lk}")
        sc = parts[f"{d:g}"]["fixedj"][bt0]["scan"].get("lam_1_Jcore_0.3")
        if sc:
            ax[1, 2].plot(Rs, [e if e is not None else np.nan for e in sc["R_test"]["E_J_of_R"]], "o-", color=c, label=f"delta {d:g} lam 1 Jcore_0.3")
    ax[1, 0].set_xlabel("L (h = 1.5)"); ax[1, 0].set_ylabel("omega*"); ax[1, 0].set_title("fixed-J omega* across the box ladder (gate <= 10 %)"); ax[1, 0].legend(fontsize=6)
    ax[1, 1].set_xlabel("R"); ax[1, 1].set_ylabel("G(R) / G(6)"); ax[1, 1].set_title("static dressing gain vs radius (n 32, L 48)"); ax[1, 1].legend(fontsize=6)
    ax[1, 2].set_xlabel("R"); ax[1, 2].set_ylabel("min_amp E_J(R)"); ax[1, 2].set_title("fixed-J energy along R (n 32, L 48, lambda 1)"); ax[1, 2].legend(fontsize=6)
    fig.suptitle("M5.32 R6.a: the delta ladder of the clock inertia (g = 32, s = -1, h = 1.5)")
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r6_deltaladder.png"), dpi=110)
    out["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(os.path.join(DATA, "m5_32_r6_deltaladder.json"), "w") as f:
        json.dump(out, f, indent=1)
    log("COLLECT written")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["delta", "collect"])
    ap.add_argument("delta", nargs="?", type=float, default=0.3)
    ap.add_argument("--maxit", type=int, default=3000)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    if a.stage == "delta":
        stage_delta(a.delta, maxit=a.maxit, quick=a.quick)
    else:
        stage_collect()
    log(f"done {a.stage}")


if __name__ == "__main__":
    main()
