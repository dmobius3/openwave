"""M5.32 R7 arm (b), class C4: the lower-order time-row term K_T as the
localizer of the boost dressing, on the R4 fixed-J instrument.

THE CANDIDATE
    L = -4 [(1 - lambda) I1 + lambda I1_h] - c2 K_T - V4
    K_T = 1/2 sum_mu eta^{mu mu} [ tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu) ]
    A_mu = d_mu M,  h = h_cov = eta + 2 (eta u)(eta u)^T,
    u = the timelike unit eigenvector of M eta (u^T eta u = -1).

EQUATIONS FIRST
---------------
(1) The u-frame form. Both traces are invariant under the internal Lorentz
    map A -> Lam A Lam^T, h -> Lam^-T h Lam^-1 (Lam in O(1,3)); in the frame
    Lam u = e0 one has h = 1 (eta e0 = -e0, so h = eta + 2 e0 e0^T = 1) and
        tr(A A) - tr(eta A eta A) = sum_ab (1 - eta_a eta_b) A_ab^2
                                  = 4 sum_j A_0j^2,
    so   K_T = sum_mu eta^{mu mu} k_mu,   k_mu = 2 sum_j (A_mu)_{0j}^2 >= 0
    (the time-row gradient squared, in the u-frame). On a block-diagonal
    static field (u = e0, A_i block-diagonal) K_T = 0: the Coulomb sector is
    untouched. On the family M = Qb Mb Qb^T (u = Qb e0 exactly) the u-frame
    jet is Qb^-1 A_i Qb^-T = G_i Mb + Mb G_i^T + d_i Mb with G_i = Qb^-1 d_i Qb
    in so(1,3); d_i Mb has no time row, so (A_i)_0j = sum_k (G_i)_0k Mb_kj
    + Mb_00 (G_i)_0j: the dressing's boost generator times ~g.
(2) Sign convention and the Legendre energy. With A_0 = omega a0,
        -c2 K_T = c2 k_0 - c2 sum_i k_i = c2 omega^2 kappa(a0) - c2 sum_i k_i,
    kappa(a0) = 2 sum_j (a0)_{0j}^2 (u-frame) >= 0. The Lagrangian piece
    c2 omega^2 kappa is a standard positive kinetic term (T = c2 omega^2
    kappa) and c2 sum_i k_i a positive static cost, so with H = omega dL/domega
    - L (quadratic in omega: H_kin = L_kin):
        E(omega) = E_stat + omega^2 kin,
        E_stat   = h^3 sum_x [ 4 sum_{i<j} q_lam(F_ij) + c2 sum_i k_i ] + E_V
        kin      = h^3 sum_x [ 4 sum_i q_lam(comm_eta(a0, A_i)) + c2 kappa(a0) ]
    (same structure as the curvature part: L = -4 I1 gives T = 4 omega^2
    sum_i q(F_0i) and V = 4 sum_{i<j} q(F_ij)). The omega^2 floor from K_T is
    c2 kappa >= 0 for c2 > 0.
(3) The fixed-J electron (the R4 instrument, unchanged):
        E_J(amp, R) = E_stat(amp, R) + J^2 / (4 kin(amp, R)),
        omega* = J / (2 kin*),  R* = argmin_R min_amp E_J.
    Family: b(r) = amp tanh(r/2) exp(-(r/R)^2) about the time axis on the
    certified m5_21_8 hedgehog Mb (m = 0), a0 = Qb a0_unit Qb^T, h = Qb^-T
    Qb^-1 analytic (checked against the registry eigensolver).
(4) Derrick for the dressing (R6.b): b_mu(r) = b(r/mu): E_KT ~ mu^1 (two
    derivatives, three volume powers), so c2 > 0 penalizes a WIDE dressing
    while staying bounded; at fixed amp E_KT(R) ~ R for R >> the tanh core.
    Continuum estimate for the exp2 window (tanh ~ 1, Mb_00 = g dominant):
    E_KT ~ 2 g^2 4 pi amp^2 R [3 sqrt(pi)/(8 sqrt 2) + sqrt(pi/2)]
         ~ 4.4e4 amp^2 R at g = 32.
(5) The static dressing gain with the term:
        G_c2(R) = min_amp [E_stat(amp, R; lam, c2)] - E_stat(0)
    (E_stat(0) = the undressed hedgehog, independent of c2 since K_T = 0 there).
(6) Gradient of E_KT wrt M (relax stage, a0 held): through A_i,
    dK/dA_i = h A_i h - eta A_i eta (chained with d1_adj); through h(M) at
    fixed A: d tr(h A h A) = 2 tr(dh S), S = A h A, dh = 2[(eta du)(eta u)^T
    + (eta u)(eta du)^T], so d(1/2 tr(hAhA)) = 4 u^T eta S eta du, and du by
    the eigenvector-perturbation formula of the R2 instrument. FD-gated.

Every number carries (n, L, h, family). Nothing is tuned: the c2 ladder is
the pre-registered geometric grid {0.03, 0.1, 0.3, 1, 3} extended DOWNWARD
by {0.001, 0.003, 0.01} (the continuum estimate (4) puts the crossover with
the R4 static gain near c2 ~ 0.03, so the pre-registered grid alone would
have no resolution below it); the R grid {6, 9, 12, 18, 24, L/2} is extended
downward by {3, 4.5} so that a minimizer running to SMALL R is seen as an
edge rather than a false interior point.

Stages (each writes partials to ../checkpoints/m5_32_r7/):
    gate                       gates 1 + controls (b), (c)
    grid n=32 L=48             the (R x amp) parts grid + the (lam, c2, J) scan
    grid n=48 L=72
    morse / relax / garm       conditional on an interior R* (task item 5)
    collect                    drifts, G7 range, controls (a), (d), plots, JSON
Out: ../data/m5_32_r7_kt_lattice.json, ../plots/m5_32_r7_kt_*.png
"""
from __future__ import annotations

import glob
import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
CK = os.path.join(HERE, "..", "checkpoints", "m5_32_r7")
os.makedirs(CK, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r7_kt_lattice.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R2 = _load("m5_32_r2_audit_lattice", "m5_32_r2_audit_lattice.py")
TX = R2.TX
B3 = R2.B3
B8 = R2.B8
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA, S = 32.0, 0.3, -1.0
KIN_BASE_REF = 426.5070121483972      # n 32, L 48, h 1.5 undressed (R2 / R4 record)
LAMS = (0.75, 1.0)
C2S = (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
C2_PREREG = (0.03, 0.1, 0.3, 1.0, 3.0)
C2_MUT = -0.1
J_RECORD = (50.0, 200.0, 800.0)
J_OMEGA_T = {0.1: 2.0 * KIN_BASE_REF * 0.1, 0.3: 2.0 * KIN_BASE_REF * 0.3}
JS = J_RECORD + tuple(J_OMEGA_T.values())
AMPS_AUDIT = np.round(np.concatenate([np.arange(0.0, 0.0601, 0.0025),
                                      [0.07, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3]]), 6)
AMPS = np.unique(np.round(np.concatenate([AMPS_AUDIT, [0.001, 0.002, 0.003, 0.004,
                                                       0.00625, 0.00875, 0.01125]]), 6))
KIND = "exp2"
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def dump(name, obj):
    with open(os.path.join(CK, name), "w") as f:
        json.dump(obj, f, indent=1)


def rd(name):
    p = os.path.join(CK, name)
    return json.load(open(p)) if os.path.exists(p) else None


# ================= the K_T density =================
def tr_hAhA(A, h):
    return np.einsum("...ab,...bc,...cd,...da->...", h, A, h, A, optimize=True)


def tr_eAeA(A):
    return np.einsum("a,...ab,b,...ba->...", np.diag(ETA), A, np.diag(ETA), A, optimize=True)


def kt_density(A, h):
    """k(A) = 1/2 [tr(h A h A) - tr(eta A eta A)] per cell (= 2 sum_j A_0j^2
    in the u-frame)."""
    return 0.5 * (tr_hAhA(A, h) - tr_eAeA(A))


def kt_parts(M, cfg, a0=None, h=None):
    """(E_KT_spatial, kappa) h^3-weighted on the sym stencil per branch;
    E_stat gets + c2 E_KT_spatial, kin gets + c2 kappa."""
    h3 = cfg["h"] ** 3
    if h is None:
        h = R2.h_of(M)
    e = k = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            e += wt * np.sum(kt_density(A[i], h))
    if a0 is not None:
        k = np.sum(kt_density(a0, h))
    return h3 * e, h3 * k


def uframe_timerow(A, Qb):
    """2 sum_j (Qb^-1 A Qb^-T)_{0j}^2 per cell: the u-frame formula, direct."""
    Qi = np.linalg.inv(Qb)
    Au = np.einsum("...ab,...bc,...dc->...ad", Qi, A, Qi)
    return 2.0 * np.sum(Au[..., 0, 1:] ** 2, axis=-1)


# ================= the localized family (R4 audit LocFamily + K_T) =================
class Fam:
    def __init__(self, n, L):
        self.n, self.L = n, float(L)
        self.cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
        self.Mb = B8.dressed(self.cfg, 0.0)
        self.a0 = B8.a0_unit(self.cfg, 0.0)
        X, Y, Z = B3.coords(n, self.cfg["h"])
        self.R = np.sqrt(X * X + Y * Y + Z * Z)
        self.cache = {}

    def b_of(self, amp, R, mu=1.0):
        r = self.R / mu
        return amp * np.tanh(r / 2.0) * np.exp(-(r / R) ** 2)

    def field(self, amp, R, mu=1.0):
        Qb = R2.qb_field(self.cfg, self.b_of(amp, R, mu))
        Md = B3.sym4(R2.conj(Qb, self.Mb))
        a0d = B3.sym4(R2.conj(Qb, self.a0))
        return Qb, Md, a0d

    def parts(self, amp, R, mu=1.0, registry_u=False):
        """(E_u0, E_u1, E_V, kin0, kin1, E_KT, kappa), all h^3-weighted."""
        key = (round(float(amp), 12), round(float(R), 9), round(mu, 9), registry_u)
        if key in self.cache:
            return self.cache[key]
        Qb, Md, a0d = self.field(amp, R, mu)
        h = R2.h_of(Md) if registry_u else R2.h_from_Q(Qb)
        p = tuple(float(v) for v in R2.lattice_parts(Md, self.cfg, a0d, h=h))
        e, k = kt_parts(Md, self.cfg, a0d, h=h)
        out = p + (float(e), float(k))
        self.cache[key] = out
        return out


def es_kin(p, lam, c2):
    eu0, eu1, ev, k0, k1, ekt, kkt = p
    return (R2.mix(eu0, eu1, lam) + ev + c2 * ekt,
            R2.mix(k0, k1, lam) + c2 * kkt)


# ================= minimization helpers =================
def parabola3(x, y):
    """vertex of the parabola through three points (x nonuniform OK);
    None if not convex."""
    (x0, x1, x2), (y0, y1, y2) = x, y
    d0, d1, d2 = (x0 - x1) * (x0 - x2), (x1 - x0) * (x1 - x2), (x2 - x0) * (x2 - x1)
    a = y0 / d0 + y1 / d1 + y2 / d2
    if a <= 1e-14:
        return None
    b = -(y0 * (x1 + x2) / d0 + y1 * (x0 + x2) / d1 + y2 * (x0 + x1) / d2)
    xs = -b / (2 * a)
    if not (min(x) < xs < max(x)):
        return None
    c = y0 * x1 * x2 / d0 + y1 * x0 * x2 / d1 + y2 * x0 * x1 / d2
    return float(xs), float(a * xs * xs + b * xs + c)


def min_over_amp(rows, J, lam, c2):
    """rows: [(amp, parts)]; E_J on the amp grid, 3-point parabolic
    refinement around an interior argmin; kin <= 0 excluded."""
    amps = np.array([a for a, _ in rows])
    es = np.array([es_kin(p, lam, c2)[0] for _, p in rows])
    kk = np.array([es_kin(p, lam, c2)[1] for _, p in rows])
    ok = kk > 0
    E = np.where(ok, es + J * J / (4.0 * np.where(ok, kk, 1.0)), np.inf)
    i = int(np.argmin(E))
    if not np.isfinite(E[i]):
        return None
    a_s, E_s, refined = float(amps[i]), float(E[i]), False
    interior = 0 < i < len(rows) - 1 and np.isfinite(E[i - 1]) and np.isfinite(E[i + 1])
    if interior:
        v = parabola3(amps[i - 1:i + 2], E[i - 1:i + 2])
        if v is not None:
            a_s, E_s, refined = v[0], v[1], True
    k_s = float(np.interp(a_s, amps, kk))
    es_s = float(np.interp(a_s, amps, es))
    ekt_s = float(np.interp(a_s, amps, [p[5] for _, p in rows]))
    edge = "interior" if interior else ("amp_min" if i == 0 else "amp_max")
    if ok.sum() and i == int(np.where(ok)[0][-1]) and not ok.all():
        edge = "kin_guard_edge"
    return {"amp_star": a_s, "E_J": E_s, "E_stat": es_s, "kin": k_s, "E_KT": ekt_s,
            "c2_E_KT": c2 * ekt_s, "J2_over_4kin": float(J * J / (4.0 * k_s)),
            "omega_star": float(J / (2.0 * k_s)), "grid_index": i, "amp_edge": edge,
            "parabolic_refined": refined, "n_kin_positive": int(ok.sum())}


def scan_R(Rs, byR):
    EJ = np.array([np.inf if byR[R] is None else byR[R]["E_J"] for R in Rs])
    i = int(np.argmin(EJ))
    fin = EJ[np.isfinite(EJ)]
    mono = bool(np.all(np.diff(fin) < 0))
    out = {"R_star_grid": float(Rs[i]), "argmin_index": i,
           "interior_grid": bool(0 < i < len(Rs) - 1),
           "at_wall": bool(i == len(Rs) - 1), "at_R_min": bool(i == 0),
           "monotone_decreasing_in_R": mono,
           "E_J_of_R": [None if not np.isfinite(e) else float(e) for e in EJ]}
    if out["interior_grid"]:
        v = parabola3(Rs[i - 1:i + 2], EJ[i - 1:i + 2])
        out["R_star_parabolic"] = None if v is None else v[0]
    return out


# ================= stage gate =================
def stage_gate(n=32, L=48.0):
    fam = Fam(n, L)
    cfg = fam.cfg
    out = {"n": n, "L": L, "h": cfg["h"], "family": KIND, "stencil": cfg["stencil"]}
    # 1a: K_T = 0 on the undressed hedgehog (analytic h at Q = 1, and the registry h)
    Qb, Md, a0d = fam.field(0.0, 6.0)
    e_an, k_an = kt_parts(Md, cfg, a0d, h=R2.h_from_Q(Qb))
    e_reg, k_reg = kt_parts(Md, cfg, a0d, h=R2.h_of(Md))
    # the scale: the same integrals of tr(h A h A) alone
    h3 = cfg["h"] ** 3
    scale = sum(wt * float(np.sum(tr_hAhA(A[i], R2.h_from_Q(Qb))))
                for br, (A, wt) in B3.a_fields(Md, cfg).items() for i in range(3)) * h3
    out["gate_undressed_zero"] = {"E_KT_analytic_h": e_an, "E_KT_registry_h": e_reg,
                                  "kappa_analytic_h": k_an, "kappa_registry_h": k_reg,
                                  "scale_sum_tr_hAhA": scale,
                                  "rel_analytic": abs(e_an) / scale, "rel_registry": abs(e_reg) / scale,
                                  "pass_1e-14": bool(abs(e_an) / scale <= 1e-14 and abs(e_reg) / scale <= 1e-14)}
    log(f"gate undressed: E_KT {e_an:.3e} (registry {e_reg:.3e}) / scale {scale:.3e}; kappa {k_an:.3e}")
    # 1b: K_T > 0 on a dressed member; analytic vs registry h; the u-frame identity
    amp, R = 0.025, 12.0
    Qb, Md, a0d = fam.field(amp, R)
    hA, hR = R2.h_from_Q(Qb), R2.h_of(Md)
    e_an, k_an = kt_parts(Md, cfg, a0d, h=hA)
    e_reg, k_reg = kt_parts(Md, cfg, a0d, h=hR)
    e_uf = 0.0
    dev = 0.0
    for br, (A, wt) in B3.a_fields(Md, cfg).items():
        for i in range(3):
            d_lab = kt_density(A[i], hA)
            d_uf = uframe_timerow(A[i], Qb)
            e_uf += wt * float(np.sum(d_uf))
            dev = max(dev, float(np.max(np.abs(d_lab - d_uf)) / (np.max(np.abs(d_uf)) + 1e-300)))
    e_uf *= h3
    kappa_uf = h3 * float(np.sum(uframe_timerow(a0d, Qb)))
    # the record clock's time row in the u-frame: max |(Qb^-1 a0 Qb^-T)_0j|
    Qi = np.linalg.inv(Qb)
    a0u = np.einsum("...ab,...bc,...dc->...ad", Qi, a0d, Qi)
    out["gate_dressed"] = {"amp": amp, "R": R, "E_KT_analytic_h": e_an, "E_KT_registry_h": e_reg,
                           "rel_dev_analytic_vs_registry": abs(e_an - e_reg) / abs(e_reg),
                           "E_KT_uframe_direct": e_uf, "max_rel_dev_density_lab_vs_uframe": dev,
                           "positive": bool(e_an > 0),
                           "kappa_record_clock": k_an, "kappa_uframe_direct": kappa_uf,
                           "max_abs_a0_timerow_uframe": float(np.max(np.abs(a0u[..., 0, 1:]))),
                           "max_abs_a0_spatial_uframe": float(np.max(np.abs(a0u[..., 1:, 1:]))),
                           "continuum_estimate_4.4e4_amp2_R": 4.4e4 * amp * amp * R}
    log(f"gate dressed (amp {amp}, R {R}): E_KT {e_an:.4f} (registry {e_reg:.4f}, u-frame {e_uf:.4f}, "
        f"density dev {dev:.1e}); kappa {k_an:.3e}; a0 time row {out['gate_dressed']['max_abs_a0_timerow_uframe']:.2e} "
        f"vs spatial {out['gate_dressed']['max_abs_a0_spatial_uframe']:.2e}; estimate {4.4e4 * amp * amp * R:.1f}")
    # 1c: E_KT(R) at fixed amp, and the dilation exponent
    Rs = [3.0, 4.5, 6.0, 9.0, 12.0, 18.0, 24.0]
    ekt = [fam.parts(amp, Rq)[5] for Rq in Rs]
    sl = np.polyfit(np.log(Rs[2:]), np.log(ekt[2:]), 1)[0]
    pair = [float(np.log(ekt[i + 1] / ekt[i]) / np.log(Rs[i + 1] / Rs[i])) for i in range(len(Rs) - 1)]
    mus = (0.8, 1.0, 1.25)
    edil = [fam.parts(amp, R, mu=m)[5] for m in mus]
    sl_mu = float(np.polyfit(np.log(mus), np.log(edil), 1)[0])
    out["gate_E_KT_vs_R"] = {"amp": amp, "R": Rs, "E_KT": ekt, "loglog_slope_R6_to_24": float(sl),
                             "pair_slopes": pair, "expected": "1 (Derrick, R >> core 2)",
                             "dilation": {"R": R, "mu": list(mus), "E_KT": edil, "loglog_slope": sl_mu,
                                          "expected": "1 for the dressing alone; the hedgehog background is not dilated"}}
    log(f"gate E_KT(R) slope {sl:.3f} pairs {np.round(pair, 3)}; dilation slope {sl_mu:.3f}")
    # controls (b), (c)
    p0 = fam.parts(0.0, 6.0)
    out["control_b_undressed_c2_independent"] = {
        "E_KT": p0[5], "E_stat_lam1_c2": {f"{c2:g}": es_kin(p0, 1.0, c2)[0] for c2 in (0.0, 0.03, 3.0)},
        "identity": "E_stat(0; c2) = E_stat(0; 0) + c2 * E_KT(0), E_KT(0) = 0 by gate 1a"}
    vac = np.broadcast_to(B3.vac4(cfg), Md.shape).copy()
    ev_, kv_ = kt_parts(vac, cfg, None, h=R2.h_of(vac))
    eu, ev = B3.e_parts(vac, cfg)
    out["control_c_vacuum_null"] = {"E_KT": ev_, "E_u_cert": float(eu), "E_V_cert": float(ev)}
    log(f"control c vacuum: E_KT {ev_:.2e} E_u {eu:.2e} E_V {ev:.2e}")
    out["runtime_s"] = round(time.time() - T0, 1)
    dump("gate.json", out)
    return out


# ================= stage grid =================
def stage_grid(n, L):
    fam = Fam(n, L)
    cfg = fam.cfg
    tag = f"n{n}_L{L:g}"
    Rs = [3.0, 4.5, 6.0, 9.0, 12.0, 18.0, 24.0]
    if L / 2.0 not in Rs:
        Rs.append(L / 2.0)
    Rs = sorted(set(Rs))
    grid_name = f"grid_{tag}.json"
    saved = rd(grid_name)
    grid = {}
    if saved is not None:
        grid = {float(R): [(a, tuple(p)) for a, p in rows] for R, rows in saved["grid"].items()}
        log(f"{tag}: grid restored ({len(grid)} R columns)")
    for R in Rs:
        if R in grid and len(grid[R]) == len(AMPS):
            continue
        grid[R] = [(float(a), fam.parts(a, R)) for a in AMPS]
        p = grid[R][AMPS.tolist().index(0.025)][1]
        log(f"{tag} R={R:g}: E_stat(lam1,c2=0) at amp 0.025 {es_kin(p, 1.0, 0.0)[0]:.3f} "
            f"kin {es_kin(p, 1.0, 0.0)[1]:.2f} E_KT {p[5]:.3f} kappa {p[6]:.2e}")
        dump(grid_name, {"n": n, "L": L, "h": cfg["h"], "family": KIND, "amps": AMPS.tolist(),
                         "R_grid": Rs, "grid": {f"{R:g}": grid[R] for R in grid}})
    p00 = grid[Rs[0]][0][1]
    E0 = {f"lam_{lam:g}": es_kin(p00, lam, 0.0)[0] for lam in LAMS}
    kin0 = {f"lam_{lam:g}": es_kin(p00, lam, 0.0)[1] for lam in LAMS}
    out = {"n": n, "L": L, "h": cfg["h"], "family": KIND, "stencil": cfg["stencil"],
           "amps": AMPS.tolist(), "R_grid": Rs, "E_stat_undressed": E0, "kin_undressed": kin0,
           "E_KT_at_amp_0.025": {f"{R:g}": grid[R][AMPS.tolist().index(0.025)][1][5] for R in Rs},
           "scan": {}, "gain": {}}
    c2_list = (0.0,) + C2S + (C2_MUT,)
    for lam in LAMS:
        for c2 in c2_list:
            key = f"lam_{lam:g}_c2_{c2:g}"
            part = rd(f"scan_{tag}_{key}.json")
            if part is not None:
                out["scan"][key] = part["scan"]
                out["gain"][key] = part["gain"]
                continue
            sc = {}
            for J in JS:
                byR = {R: min_over_amp(grid[R], J, lam, c2) for R in Rs}
                st = scan_R(Rs, byR)
                best = byR[st["R_star_grid"]]
                rec = {"J": J, "byR": {f"{R:g}": byR[R] for R in Rs}, "R_test": st}
                # refinement: a short amp column at the parabolic R*
                R_ref = st.get("R_star_parabolic")
                if R_ref is not None and best is not None and best["amp_edge"] == "interior":
                    i = best["grid_index"]
                    sub = [(float(a), fam.parts(a, R_ref)) for a in AMPS[max(0, i - 2):i + 3]]
                    m = min_over_amp(sub, J, lam, c2)
                    rec["refined"] = dict(m, R=R_ref)
                    use = dict(m, R=R_ref) if (m is not None and m["E_J"] < best["E_J"]) else dict(best, R=st["R_star_grid"])
                else:
                    use = None if best is None else dict(best, R=st["R_star_grid"])
                dressed = use is not None and use["amp_star"] > 0.0 and use["amp_edge"] != "amp_min"
                interior = bool(st["interior_grid"] and dressed)
                rec["result"] = None if use is None else {
                    "R_star": use["R"], "interior": interior, "dressed": bool(dressed),
                    "R_edge": "interior" if st["interior_grid"] else ("wall" if st["at_wall"] else "R_min"),
                    "amp_star": use["amp_star"], "amp_edge": use["amp_edge"], "omega_star": use["omega_star"],
                    "E_total": use["E_J"], "E_stat": use["E_stat"], "kin": use["kin"],
                    "c2_E_KT": use["c2_E_KT"], "J2_over_4kin": use["J2_over_4kin"],
                    "kin_dressing": use["kin"] - kin0[f"lam_{lam:g}"]}
                sc[f"J_{J:.6g}"] = rec
            # the static gain G_c2(R)
            gain = {}
            for R in Rs:
                amps = np.array([a for a, _ in grid[R]])
                es = np.array([es_kin(p, lam, c2)[0] for _, p in grid[R]])
                i = int(np.argmin(es))
                a_s, e_s = float(amps[i]), float(es[i])
                if 0 < i < len(amps) - 1:
                    v = parabola3(amps[i - 1:i + 2], es[i - 1:i + 2])
                    if v is not None:
                        a_s, e_s = v
                gain[f"{R:g}"] = {"amp_star": a_s, "E_stat": e_s, "gain": e_s - E0[f"lam_{lam:g}"],
                                  "amp_edge": "interior" if 0 < i < len(amps) - 1 else ("amp_min" if i == 0 else "amp_max")}
            gs = np.array([gain[f"{R:g}"]["gain"] for R in Rs])
            ig = int(np.argmin(gs))
            gain["summary"] = {"R_star_static": float(Rs[ig]), "gain_min": float(gs[ig]),
                               "interior": bool(0 < ig < len(Rs) - 1),
                               "saturates": bool(len(gs) > 2 and abs(gs[-1] - gs[-2]) < 0.05 * abs(gs[-1]) + 1e-12),
                               "gain_of_R": gs.tolist()}
            out["scan"][key] = sc
            out["gain"][key] = gain
            dump(f"scan_{tag}_{key}.json", {"scan": sc, "gain": gain})
            r = sc["J_200"]["result"]
            log(f"{tag} {key}: J200 R* {r['R_star']:g} [{r['R_edge']}{'' if r['dressed'] else ', undressed'}] "
                f"amp* {r['amp_star']:.4f} omega* {r['omega_star']:.4f} E {r['E_total']:.2f} c2E_KT {r['c2_E_KT']:.2f}; "
                f"static G* {gain['summary']['gain_min']:.2f} at R {gain['summary']['R_star_static']:g}"
                f"{' (int)' if gain['summary']['interior'] else ''}")
    # control (a): c2 = 0 on the audit's amp subset
    ctrl = {}
    sub_idx = [AMPS.tolist().index(a) for a in AMPS_AUDIT.tolist()]
    for lam in LAMS:
        for J in (200.0, 800.0):
            byR = {R: min_over_amp([grid[R][i] for i in sub_idx], J, lam, 0.0) for R in Rs if R >= 6.0}
            Rq = sorted(byR)
            st = scan_R(Rq, byR)
            ctrl[f"lam_{lam:g}_J_{J:g}"] = {"R_star": st["R_star_grid"], "at_wall": st["at_wall"],
                                            "omega_star_wall": byR[L / 2.0]["omega_star"],
                                            "omega_star_R12": byR[12.0]["omega_star"],
                                            "E_J_wall": byR[L / 2.0]["E_J"], "amp_star_wall": byR[L / 2.0]["amp_star"]}
    out["control_a_c2_0_audit_subset"] = ctrl
    out["runtime_s"] = round(time.time() - T0, 1)
    dump(f"box_{tag}.json", out)
    log(f"{tag} done")
    return out


# ================= stage morse (2x2 Hessian over (amp, R) at an interior point) =================
def stage_morse(n, L, lam, c2, J, amp, R, da=0.002, dR=0.5):
    fam = Fam(n, L)

    def EJ(a, r):
        es, kk = es_kin(fam.parts(a, r), lam, c2)
        return es + J * J / (4.0 * kk)

    f0 = EJ(amp, R)
    fa = (EJ(amp + da, R) - 2 * f0 + EJ(amp - da, R)) / da ** 2
    fR = (EJ(amp, R + dR) - 2 * f0 + EJ(amp, R - dR)) / dR ** 2
    faR = (EJ(amp + da, R + dR) - EJ(amp + da, R - dR) - EJ(amp - da, R + dR) + EJ(amp - da, R - dR)) / (4 * da * dR)
    ga = (EJ(amp + da, R) - EJ(amp - da, R)) / (2 * da)
    gR = (EJ(amp, R + dR) - EJ(amp, R - dR)) / (2 * dR)
    H = np.array([[fa, faR], [faR, fR]])
    ev = np.linalg.eigvalsh(H)
    out = {"n": n, "L": L, "h": fam.cfg["h"], "family": KIND, "lam": lam, "c2": c2, "J": J,
           "amp": amp, "R": R, "E_J": f0, "grad": [ga, gR], "hessian": H.tolist(),
           "eigen": ev.tolist(), "morse_index": int(np.sum(ev < 0)), "da": da, "dR": dR}
    dump(f"morse_n{n}_L{L:g}_lam{lam:g}_c2{c2:g}_J{J:g}.json", out)
    log(f"morse ({lam}, {c2}, J {J}) at amp {amp:.4f} R {R:.2f}: grad {ga:+.3e} {gR:+.3e} eig {ev} index {out['morse_index']}")
    return out


# ================= stage relax (true fixed-J descent with K_T) =================
def kt_grad(M, cfg, a0, c2):
    """c2 [E_KT_spatial + kappa(a0)]-gradient wrt symmetric M (a0 held)."""
    RB = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
    h3, h = cfg["h"] ** 3, cfg["h"]
    u0, V, lamv, sig, k0, ok, gap = RB.tl_eig(M)
    if not np.all(ok):
        return None
    hh = RB.h_of(u0)
    Gd = np.zeros_like(M)
    Ssum = np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        for i in range(3):
            dA = (hh @ A[i] @ hh) - (ETA @ A[i] @ ETA)          # dK/dA_i (symmetric A)
            Gd += wt * B3.d1_adj(dA, i, h, br)
            Ssum += wt * (A[i] @ hh @ A[i])
    if a0 is not None:
        Ssum += a0 @ hh @ a0
    v = 4.0 * np.einsum("ab,...bc,cd,...d->...a", ETA, Ssum, ETA, u0)
    vu = np.einsum("...a,...ak->...k", v, V)
    l0 = np.take_along_axis(lamv, k0[..., None], axis=-1)[..., 0]
    den = l0[..., None] - lamv
    mask = (np.arange(4)[None, :] != k0.reshape(-1, 1)).reshape(den.shape)
    c = np.where(mask, sig * vu / np.where(mask, den, 1.0), 0.0)
    w = np.einsum("...ak,...k->...a", V, c)
    Gd += (w @ ETA)[..., :, None] * (u0 @ ETA)[..., None, :]
    return c2 * h3 * B3.sym4(Gd)


def stage_relax(n, L, lam, c2, J, amp, R, steps_acc=300, it_cap=900, dt0=0.02, dt_max=0.2):
    RB = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
    R4 = _load("m5_32_r4_clock", "m5_32_r4_clock.py")
    fam = Fam(n, L)
    cfg = fam.cfg
    Qb, M, a0 = fam.field(amp, R)

    def EJ_grad(Mq):
        E, Gs, info = RB.energy_grad(Mq, cfg, lam)
        if Gs is None:
            return np.nan, None, info, None
        rdd = RB.reads(Mq, cfg, a0)
        ekt, kkt = kt_parts(Mq, cfg, a0)
        k = RB.kin_lam(rdd, lam) + c2 * kkt
        if k <= 0:
            return np.nan, None, dict(info, ok=False), None
        Gk = R4.kin_grad_lam(Mq, a0, cfg, lam)
        Gkt = kt_grad(Mq, cfg, a0, c2)
        if Gk is None or Gkt is None:
            return np.nan, None, dict(info, ok=False), None
        # kt_grad carries both the spatial and the kappa pieces; split: the
        # kappa piece enters through kin, the spatial through E_stat
        Gkt_sp = kt_grad(Mq, cfg, None, c2)
        Gkap = Gkt - Gkt_sp
        Es = E + c2 * ekt
        return Es + J * J / (4 * k), Gs + Gkt_sp - J * J / (4 * k * k) * (Gk + Gkap), info, (Es, k)

    # FD gate of the K_T gradient (Richardson pair) on a random smooth direction
    rng = np.random.default_rng(7)
    D = B3.sym4(rng.normal(size=M.shape)) * np.exp(-(fam.R / 8.0) ** 2)[..., None, None]
    Gkt = kt_grad(M, cfg, a0, c2)
    lin = float(np.sum(Gkt * D))
    fdv = []
    for eps in (1e-3, 5e-4):
        ep, _ = kt_parts(M + eps * D, cfg, a0); em, _ = kt_parts(M - eps * D, cfg, a0)
        kp = kt_parts(M + eps * D, cfg, a0)[1]; km = kt_parts(M - eps * D, cfg, a0)[1]
        fdv.append(c2 * ((ep + kp) - (em + km)) / (2 * eps))
    fd = (4 * fdv[1] - fdv[0]) / 3.0
    fd_gate = {"linear": lin, "fd_richardson": float(fd), "rel_dev": abs(lin - fd) / max(abs(fd), 1e-300)}
    log(f"relax FD gate of the K_T gradient: lin {lin:.6e} fd {fd:.6e} rel {fd_gate['rel_dev']:.2e}")
    free = (~B3.pin_shell(n, cfg["h"]))[..., None, None].astype(float)
    EJ0, Gr, info, parts = EJ_grad(M)
    out = {"n": n, "L": L, "h": cfg["h"], "family": KIND, "lam": lam, "c2": c2, "J": J,
           "start": {"amp": amp, "R": R, "E_J": EJ0, "E_stat": parts[0], "kin": parts[1], "omega": J / (2 * parts[1])},
           "fd_gate_kt_grad": fd_gate, "method": "R4 relax pattern: E_J[M] = E_stat[M] + J^2/(4 kin[M]), a0 held, "
           "exact gradients (RB.energy_grad + R4.kin_grad_lam + kt_grad), backtracking FIRE, pinned shell 1.6",
           "steps_acc_budget": steps_acc, "trace": []}
    log(f"RELAX start E_J {EJ0:.5f} E_stat {parts[0]:.4f} kin {parts[1]:.3f}")
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    F = -Gr * free
    E_prev, n_acc, n_rej, it = EJ0, 0, 0, 0
    stop = "budget"
    while it < it_cap and n_acc < steps_acc:
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
        E, Gr, info, parts = EJ_grad(M_try)
        reject = (Gr is None) or not np.isfinite(E) or E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1; dt *= 0.5; v[:] = 0.0; alpha, n_up = 0.1, 0
            if dt < 1e-7:
                stop = "STALLED (dt collapsed)" if Gr is not None else "LOCUS-HIT"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -Gr * free
        if n_acc % 25 == 0 or n_acc == steps_acc:
            row = {"it": it, "acc": n_acc, "E_J": float(E), "E_stat": parts[0], "kin": parts[1],
                   "omega": J / (2 * parts[1]), "fmax": float(np.max(np.abs(F))), "dt": dt, "rej": n_rej}
            out["trace"].append(row)
            log(f"RELAX acc {n_acc:4d} it {it:4d} E_J {E:.5f} E_stat {parts[0]:.4f} kin {parts[1]:.3f} "
                f"omega {row['omega']:.5f} fmax {row['fmax']:.2e} dt {dt:.1e} rej {n_rej}")
    end = {"E_J": float(E_prev), "E_stat": float(parts[0]) if parts else None,
           "kin": float(parts[1]) if parts else None, "omega": float(J / (2 * parts[1])) if parts else None}
    st = out["start"]
    out["end"] = end
    out["drift"] = {q: float((end[q] - st[q]) / abs(st[q])) for q in ("E_J", "E_stat", "kin", "omega") if end[q] is not None}
    out["stop"], out["accepted"], out["rejected"], out["iterations"] = stop, n_acc, n_rej, it
    tr = out["trace"]
    if len(tr) >= 4:
        q = [r for r in tr if r["acc"] >= 0.75 * n_acc]
        out["last_quarter_rel"] = float(abs(q[-1]["E_J"] - q[0]["E_J"]) / max(abs(q[-1]["E_J"]), 1.0))
        out["verdict"] = "PLATEAU" if out["last_quarter_rel"] <= 1e-3 else "FALLING (still descending at the budget)"
    # the radial profile of the relaxed dressing: the u-frame boost rapidity read off u(M)
    u0 = RB.tl_eig(M)[0]
    rap = np.arctanh(np.clip(np.linalg.norm(u0[..., 1:], axis=-1) / np.abs(u0[..., 0]), 0, 1 - 1e-12))
    rb = fam.R.ravel(); rp = rap.ravel()
    bins = np.arange(0, L / 2 + 3, 3.0)
    prof = [float(np.mean(rp[(rb >= bins[i]) & (rb < bins[i + 1])])) for i in range(len(bins) - 1)]
    prof0 = [float(np.mean(fam.b_of(amp, R).ravel()[(rb >= bins[i]) & (rb < bins[i + 1])])) for i in range(len(bins) - 1)]
    out["rapidity_profile"] = {"r_bins": bins.tolist(), "relaxed": prof, "start_family": prof0}
    out["wall_s"] = round(time.time() - T0, 1)
    np.savez_compressed(os.path.join(CK, f"relax_n{n}_L{L:g}_lam{lam:g}_c2{c2:g}_J{J:g}.npz"), M=M, a0=a0)
    dump(f"relax_n{n}_L{L:g}_lam{lam:g}_c2{c2:g}_J{J:g}.json", out)
    log(f"RELAX end {stop}: drift " + ", ".join(f"{k} {v:+.4%}" for k, v in out["drift"].items()))
    return out


# ================= stage garm (rigid boost dressing gain vs g with K_T) =================
def stage_garm(lam, c2, gs=(8.0, 32.0), branches=("A",), nm=121, span=3.0):
    GM = _load("m5_21_11_d_garm", "m5_21_11_d_garm.py")
    out = {"lam": lam, "c2": c2, "arms": []}
    for br in branches:
        M3 = np.load(os.path.join(DATA, f"m5_21_11_end_t11lad_{br}_n{GM.N}_d{GM.DELTA:g}.npz"))["M"].astype(np.float64)
        for g in gs:
            cfg = B3.base_cfg(s=-1.0, g=g, n=GM.N, L=48.0, delta=GM.DELTA)
            M4 = B3.embed34(M3, cfg)
            m_his = float(np.arctanh(1.0 / g))
            ms = np.linspace(-span * m_his, span * m_his, nm)
            Es = []
            for m in ms:
                Qb = GM.qb_field(cfg, m)
                Md = B3.sym4(R2.conj(Qb, M4))
                h = R2.h_from_Q(Qb)
                p = R2.lattice_parts(Md, cfg, None, h=h)
                ekt, _ = kt_parts(Md, cfg, None, h=h)
                Es.append(R2.mix(p[0], p[1], lam) + p[2] + c2 * ekt)
            Es = np.array(Es)
            i = int(np.argmin(Es))
            E0 = float(Es[nm // 2])
            m_s, E_s, edge = float(ms[i]), float(Es[i]), not (0 < i < nm - 1)
            if not edge:
                v = parabola3(ms[i - 1:i + 2], Es[i - 1:i + 2])
                if v is not None:
                    m_s, E_s = v
            out["arms"].append({"branch": br, "g": g, "n": GM.N, "L": 48.0, "h": cfg["h"], "m_his": m_his,
                                "E0": E0, "m_star": m_s, "E_star": E_s, "gain": E_s - E0, "edge_minimum": edge,
                                "curve": [{"m": float(a), "E": float(b)} for a, b in zip(ms[::6], Es[::6])]})
            log(f"garm {br} g {g:g}: gain {E_s - E0:+.4f} at m* {m_s:+.4f} (m_his {m_his:.4f}){' EDGE' if edge else ''}")
    dump(f"garm_lam{lam:g}_c2{c2:g}.json", out)
    return out


# ================= stage collect =================
def stage_collect():
    gate = rd("gate.json")
    boxes = {t: rd(f"box_{t}.json") for t in ("n32_L48", "n48_L72")}
    boxes = {t: b for t, b in boxes.items() if b is not None}
    out = {"task": "M5.32 R7.b: class C4, K_T on the R4 fixed-J instrument",
           "candidate": "L = -4[(1-lambda) I1 + lambda I1_h] - c2 K_T - V4, K_T = 1/2 sum_mu eta^mumu [tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu)]",
           "conventions": {"energy": "E = E_stat + omega^2 kin; E_stat = h^3 sum [4 sum_{i<j} q_lam(F_ij) + c2 sum_i k_i] + E_V; kin = h^3 sum [4 sum_i q_lam(comm_eta(a0, A_i)) + c2 kappa(a0)]; k_i = 1/2[tr(h A_i h A_i) - tr(eta A_i eta A_i)] = 2 sum_j (A_i)_0j^2 in the u-frame; kappa the same on a0; the omega^2 floor c2 kappa >= 0 for c2 > 0",
                           "sign": "-c2 K_T in L = +c2 omega^2 kappa (kinetic, positive) - c2 sum_i k_i (static cost); H = omega dL/domega - L",
                           "fixed_J": "E_J = E_stat + J^2/(4 kin); omega* = J/(2 kin*); min over the amp grid (3-point parabola) per R; R* = argmin on the R grid, 3-point parabola in R, re-minimized on a 5-amp column at the parabolic R",
                           "family": "b(r) = amp tanh(r/2) exp(-(r/R)^2) about the time axis on the m5_21_8 hedgehog (m = 0); a0 = Qb a0_unit Qb^T; h = Qb^-T Qb^-1 analytic",
                           "toy": {"g": G, "delta": DELTA, "s": S}, "amps": AMPS.tolist(),
                           "c2_ladder": list(C2S), "c2_preregistered": list(C2_PREREG), "c2_mutation": C2_MUT,
                           "J_ladder": {"record": list(J_RECORD), "omega_t": {f"{w:g}": J for w, J in J_OMEGA_T.items()}, "KIN_BASE_REF": KIN_BASE_REF},
                           "interior": "R* strictly inside the R grid AND amp* > 0 (a dressing exists)",
                           "drift_gate": "|omega*(L72)/omega*(L48) - 1| <= 0.10 at fixed h = 1.5"},
           "gate": gate, "boxes": boxes, "localization": [], "drift": {}, "G7_range": {}, "controls": {}}
    c2_list = (0.0,) + C2S + (C2_MUT,)
    # the localization table + drifts
    for lam in LAMS:
        for c2 in c2_list:
            key = f"lam_{lam:g}_c2_{c2:g}"
            for J in JS:
                Jk = f"J_{J:.6g}"
                row = {"lam": lam, "c2": c2, "J": J}
                for t, b in boxes.items():
                    r = b["scan"][key][Jk]["result"]
                    row[t] = r
                if all(t in row and row[t] is not None for t in ("n32_L48", "n48_L72")):
                    a, bb = row["n32_L48"], row["n48_L72"]
                    row["omega_drift"] = bb["omega_star"] / a["omega_star"] - 1.0
                    row["R_drift"] = bb["R_star"] / a["R_star"] - 1.0
                    row["amp_drift"] = (bb["amp_star"] / a["amp_star"] - 1.0) if a["amp_star"] > 0 else None
                    row["kin_dressing_drift"] = ((bb["kin_dressing"] / a["kin_dressing"] - 1.0)
                                                 if abs(a["kin_dressing"]) > 1e-9 else None)
                    row["E_total_drift"] = bb["E_total"] / a["E_total"] - 1.0
                    row["interior_both"] = bool(a["interior"] and bb["interior"])
                    row["drift_gate_met"] = bool(abs(row["omega_drift"]) <= 0.10)
                out["localization"].append(row)
    # the G7 range per (lam, J): the c2 interval with interior R* in both boxes (+ the drift gate)
    for lam in LAMS:
        for J in JS:
            rows = [r for r in out["localization"] if r["lam"] == lam and r["J"] == J and r["c2"] > 0]
            c_int = [r["c2"] for r in rows if r.get("interior_both")]
            c_int48 = [r["c2"] for r in rows if r["n32_L48"] and r["n32_L48"]["interior"]]
            c_gate = [r["c2"] for r in rows if r.get("interior_both") and r.get("drift_gate_met")]
            out["G7_range"][f"lam_{lam:g}_J_{J:.6g}"] = {
                "c2_interior_L48": c_int48, "c2_interior_both_boxes": c_int, "c2_interior_and_drift_gate": c_gate,
                "factor_interior_both": (max(c_int) / min(c_int)) if c_int else None,
                "factor_with_drift_gate": (max(c_gate) / min(c_gate)) if c_gate else None,
                "G7_met": bool(c_gate and max(c_gate) / min(c_gate) >= 4.0)}
    # controls (a): against the R4 audit record
    ref = json.load(open(os.path.join(DATA, "m5_32_r4_audit_clock.json")))
    ca = {}
    for t, b in boxes.items():
        rk = "box_" + t
        if rk in ref["parts"]:
            for k, v in b["control_a_c2_0_audit_subset"].items():
                sc = ref["parts"][rk]["kinds"]["exp2"]["scan"][k]
                Rw = f"R_{b['L'] / 2:g}"
                ca[f"{t}_{k}"] = {"omega_star_wall_here": v["omega_star_wall"], "omega_star_wall_audit": sc["byR"][Rw]["omega_star"],
                                  "rel_dev": abs(v["omega_star_wall"] / sc["byR"][Rw]["omega_star"] - 1),
                                  "R_star_here": v["R_star"], "R_star_audit": sc["R_test"]["R_star"],
                                  "E_J_wall_here": v["E_J_wall"], "E_J_wall_audit": sc["byR"][Rw]["E_J"]}
    out["controls"]["a_c2_0_vs_R4_audit"] = ca
    out["controls"]["b_undressed_c2_independent"] = gate["control_b_undressed_c2_independent"] if gate else None
    out["controls"]["c_vacuum_null"] = gate["control_c_vacuum_null"] if gate else None
    mut = {}
    for lam in LAMS:
        for t, b in boxes.items():
            r = b["scan"][f"lam_{lam:g}_c2_{C2_MUT:g}"]["J_200"]["result"]
            g = b["gain"][f"lam_{lam:g}_c2_{C2_MUT:g}"]["summary"]
            mut[f"lam_{lam:g}_{t}"] = {"R_star": r["R_star"], "R_edge": r["R_edge"], "amp_edge": r["amp_edge"],
                                       "amp_star": r["amp_star"], "E_total": r["E_total"],
                                       "static_gain_min": g["gain_min"], "static_R_star": g["R_star_static"],
                                       "unbounded_signature": bool(r["amp_edge"] == "amp_max" or r["R_edge"] == "wall")}
    out["controls"]["d_mutation_c2_-0.1"] = mut
    # conditional stages
    out["morse"] = {os.path.basename(p)[6:-5]: json.load(open(p)) for p in sorted(glob.glob(os.path.join(CK, "morse_*.json")))}
    out["relax"] = {os.path.basename(p)[6:-5]: json.load(open(p)) for p in sorted(glob.glob(os.path.join(CK, "relax_*.json")))}
    out["garm"] = {os.path.basename(p)[5:-5]: json.load(open(p)) for p in sorted(glob.glob(os.path.join(CK, "garm_*.json")))}
    out["runtime_s_collect"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    plots(out)
    log("collected")
    return out


def plots(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    boxes = out["boxes"]
    # plot 1: E_J(R) per c2 at J = 200, both lambdas, both boxes
    fig, ax = plt.subplots(2, len(boxes), figsize=(6 * len(boxes), 8), squeeze=False)
    for j, (t, b) in enumerate(boxes.items()):
        Rs = b["R_grid"]
        for i, lam in enumerate(LAMS):
            for c2 in (0.0,) + C2S:
                sc = b["scan"][f"lam_{lam:g}_c2_{c2:g}"]["J_200"]
                E = sc["R_test"]["E_J_of_R"]
                ax[i, j].plot(Rs, E, "o-", ms=3, label=f"c2 {c2:g}")
                r = sc["result"]
                if r and r["interior"]:
                    ax[i, j].plot([r["R_star"]], [r["E_total"]], "k*", ms=9)
            ax[i, j].set_xscale("log"); ax[i, j].set_yscale("symlog", linthresh=10)
            ax[i, j].set_xlabel("R (dressing radius)"); ax[i, j].set_ylabel("min_amp E_J")
            ax[i, j].set_title(f"lambda {lam:g}, J 200, n {b['n']} L {b['L']:g} h {b['h']:g} ({KIND})")
            ax[i, j].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r7_kt_localization.png"), dpi=110); plt.close(fig)
    # plot 2: the static gain G_c2(R)
    fig, ax = plt.subplots(1, len(boxes), figsize=(6 * len(boxes), 4.5), squeeze=False)
    for j, (t, b) in enumerate(boxes.items()):
        Rs = b["R_grid"]
        for lam, ls in zip(LAMS, ("-", "--")):
            for c2 in (0.0,) + C2S:
                g = b["gain"][f"lam_{lam:g}_c2_{c2:g}"]["summary"]["gain_of_R"]
                ax[0, j].plot(Rs, g, ls, marker="o", ms=3, label=f"lam {lam:g} c2 {c2:g}")
        ax[0, j].set_xscale("log"); ax[0, j].set_xlabel("R"); ax[0, j].set_ylabel("G_c2(R) = min_amp E_stat - E_stat(0)")
        ax[0, j].set_title(f"static dressing gain, n {b['n']} L {b['L']:g} h {b['h']:g}"); ax[0, j].legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r7_kt_gain.png"), dpi=110); plt.close(fig)
    # plot 3: R*, amp*, omega* vs c2 (J = 200) per box
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for lam, ls in zip(LAMS, ("-", "--")):
        for t, b in boxes.items():
            rows = [r for r in out["localization"] if r["lam"] == lam and r["J"] == 200.0 and r["c2"] > 0 and r.get(t)]
            c2s = [r["c2"] for r in rows]
            ax[0].plot(c2s, [r[t]["R_star"] for r in rows], ls, marker="o", ms=4, label=f"lam {lam:g} {t}")
            ax[1].plot(c2s, [r[t]["amp_star"] for r in rows], ls, marker="o", ms=4, label=f"lam {lam:g} {t}")
            ax[2].plot(c2s, [r[t]["omega_star"] for r in rows], ls, marker="o", ms=4, label=f"lam {lam:g} {t}")
    for a, yl in zip(ax, ("R* (wall = L/2, R_min = 3)", "amp*", "omega* (J 200)")):
        a.set_xscale("log"); a.set_xlabel("c2"); a.set_ylabel(yl); a.legend(fontsize=7)
    ax[0].set_yscale("log")
    fig.suptitle("fixed-J minimizer vs c2 (h 1.5, exp2 family)")
    fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r7_kt_range.png"), dpi=110); plt.close(fig)


if __name__ == "__main__":
    mode = sys.argv[1]
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "gate":
        stage_gate()
    elif mode == "grid":
        stage_grid(int(kw.get("n", 32)), float(kw.get("L", 48)))
    elif mode == "morse":
        stage_morse(int(kw["n"]), float(kw["L"]), float(kw["lam"]), float(kw["c2"]), float(kw["J"]),
                    float(kw["amp"]), float(kw["R"]))
    elif mode == "relax":
        stage_relax(int(kw["n"]), float(kw["L"]), float(kw["lam"]), float(kw["c2"]), float(kw["J"]),
                    float(kw["amp"]), float(kw["R"]), steps_acc=int(kw.get("steps", 300)))
    elif mode == "garm":
        stage_garm(float(kw["lam"]), float(kw["c2"]))
    elif mode == "collect":
        stage_collect()
