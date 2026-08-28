"""M5.32 R2 arm (c): the LATTICE SCREEN of the C2 sign map (the lambda-family,
J2, Pgrad) at the certified toy point, on the Coulomb pair instrument, on
the M5.21.11 g-arm, and under the full lattice-jet Lorentz action.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), jets A_i = d_i M
(i = 1..3, the certified sym stencil: fwd / bwd branches, the density
evaluated per branch and averaged), A_0 = omega a0 (a0 a channel).
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu
Terms (registries m5_32_lagrangian.py + m5_32_terms_ext.py, imported,
never modified; every definition string hashed there):
    I1      = (1/2) F_{mu nu a b} F^{mu nu a b}         (certified)
    I1_frob = I1 with the Frobenius internal metric      (control)
    I1_h    = sum_{mu<nu} eta^mu eta^nu tr(h F h F^T),
              h = eta + 2 (eta u)(eta u)^T, u the timelike unit eigenvector
              of M eta (u^T eta u = -1); h = 1 in the vacuum eigenframe
    J1      = sum_{mu<nu} eta^mu eta^nu tr(F eta M eta F eta M eta)
    J2      = sum_{mu<nu} eta^mu eta^nu tr(F eta F eta M eta M eta)
    Pgrad   = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t), P_t = u u^T eta,
              d_mu u by first-order eigenvector perturbation of M eta
The lambda-family (the R1 close-out's live C2 member):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    E_lambda = +4 [(1 - lambda) I1 + lambda I1_h]_static-Hamiltonian + V4
Every term is exactly quadratic in omega, I(omega) = A + B omega
+ C omega^2 (3-point exact), Hamiltonian H = C omega^2 - A. Energy-read
kinetic coefficient C_E(T) := -4 C_T (so C_E(I1) = the certified
B3.kin_of, checked); for the family K(lambda) = (1 - lambda) C_E(I1)
+ lambda C_E(I1_h) per channel; "stable" = K >= 0 on every channel,
"boost reversed" = K > 0 on the boost channels (I1's is < 0).

Channels (the R1 audit construction, m5_32_r1_audit_screen.channels,
imported): a0 = w(x)(G M + M G^T)/||.|| the TANGENT (physical) form,
a0 = w(x)(G M - M G^T)/||.|| the record PROBE form (antisymmetric), G in
{clock_local, plane_1d, rot_z, rot_x, boost_z, boost_x} plus the EXTRA
channels (radial / azimuthal boost twists, boost + rotation mixes, a
random smooth symmetric a0, a random antisymmetric probe).

Coulomb pair instrument (m5_32_r1_b_pair.py: the certified m5_21_4
seeds, FIRE loop, pinned shell, charge suite, imported): static 3x3 M on
n = 32, L = 48, it = 120, embedded as M4 = diag(g_tt) (+) M3 with
g_tt = 32 (the toy point, s = -1: M_vac = diag(g, 1, delta, 0)),
A_0 = 0, A_i = 0 (+) d_i M3. Energy relaxed:
    E[M] = 4 h^3 sum_br wt sum_cells D(F, M4) + h^3 sum V_T2(M)
    D = (1 - lambda) I1 + lambda I1_h           (identity run)
    D = I1 + c J2                               (footprint run)
Gradients (exact, chained through the certified stencil adjoints):
    I1_h: W_{mu nu} = eta_mu eta_nu h F_{mu nu} h (h frozen; ordered pairs,
          the 1/2 of the ordered-pair sum times the 2 of d/dF), plus the
          LOCAL M-dependence through u: with D_h = 2 sum_{mu<nu} eta_mu
          eta_nu F h F^T, d dens = 4 (eta du)^T D_h (eta u), and from
          (N - l0) du = -(dN - dl0) u, dl0 = -(w^T dN u), N = M eta,
          w = eta u: du = -(I + u w^T)(N - l0)^+ (I + u w^T) dM eta u
          (Moore-Penrose pseudo-inverse of the singular resolvent, the
          u-component of du fixed by u^T eta du = 0), so
          grad_M = -4 sym[(R^T y) w^T], R = (I + u w^T)(N - l0)^+(I + u w^T),
          y = eta D_h eta u.
    J2:   W_{mu nu} = -(1/2) eta_mu eta_nu (P F E + E F P), E = eta,
          P = E M E M E (ordered pairs); local: with Z = F E F E,
          Q = E M E: grad_M = sym[Z^T Q + E M Z^T E] per pair (eta-weighted).
    Both gradients are checked against central finite differences of the
    lattice energy along random directions (gradient_check in the JSON).
Static-sector facts used (and checked numerically, not assumed):
    on block-diagonal M4 with F in the spatial block, u = e0 exactly, so
    h = 1 on the block and I1_h = I1 (density AND gradient: the local
    piece is proportional to F^{0b} = 0); Pgrad is quadratic in du and du
    is linear in the off-block jets, so Pgrad = 0 with zero gradient on
    the whole 3x3 sector (checked as density = 0 and a directional
    derivative = 0 on the relaxed control end states).
Reads: E_int(d) = E_pair(d) - 2 E_single, fit A + B/d over d = 12, 18,
24 (B > 0 repulsive), same / anti; sign pattern relative to the lambda =
0 control; the (a, b, c_) content (I1 = 2a, I2 = 4b, I6 = 4c_) of the
relaxed single; degrees from the Mermin-Ho flux suite.

g-arm (import I19, m5_21_11_d_garm.py construction): the N = 48,
delta = 0.3 production endpoints (branches A, C, B) embedded at g and
rigidly dressed by the radial boost hedgehog Qb(m) = 1 + sinh(m) K +
(cosh(m) - 1) K2; per g the curve E_lambda(m) on 121 points of
[-3 m_his, 3 m_his], m_his = artanh(1/g), with
    E_lambda(m) = (1 - lambda) E_u(m) + lambda 4 h^3 sum I1_h(m) + E_V(m)
(E_u the certified e_parts curvature energy = 4 h^3 sum I1);
GAIN(g) := E_lambda(m*) - E_lambda(0), m* the argmin (parabolic refine,
the record's recipe); q := slope of log(-GAIN) vs log(artanh(1/g)) over
the record set g = 8, 16, 32 (the certified action: q ~ 0, FLAT; the
M5.21.11 F4 bar needed q >= 2). Controls under every lambda: C1 vacuum
null (the constant vacuum dressed: minimum at m = 0, zero gain), C2 field
identity (the Qb pipeline vs m5_21_8_b_lattice.dressed), C3 record match
(lambda = 0 gains vs m5_21_11_garm.json to 1e-10).

Covariance (G4, lattice level): the R1.c stage_inv (b) construction on a
random smooth field (n = 14, L = 21, rng 3203): A'_mu = (Lambda^{-T})_mu^nu
Lambda A_nu Lambda^T, M' = Lambda M Lambda^T on the certified stencil
jets; drift = |E' - E| / |E| of the h^3 density sum; the no-eta controls
(the ext registry's eta_time_row mutant of each new term, plus I1 with a
delta derivative pair and I2 on the plain commutator) MUST fail under
boosts.

Stages: chan, coulomb, garm, inv (or all). Every stage merges into the
JSON so a stage can be re-run alone.
Usage: python3 m5_32_r2_c_screen.py [--stage all] [--workers 12] [--plot]
Out: ../data/m5_32_r2_screen.json, ../plots/m5_32_r2_screen.png,
     ../data/m5_32_r2_b/*.npz (local, gitignored)
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r2_screen.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r2_b")
PLOT = os.path.join(PLOTS, "m5_32_r2_screen.png")
T_START = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
AUD = _load("m5_32_r1_audit_screen", "m5_32_r1_audit_screen.py")
PB = _load("m5_32_r1_b_pair", "m5_32_r1_b_pair.py")
GM = _load("m5_21_11_d_garm", "m5_21_11_d_garm.py")
L8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")
B3 = LAG.B3
INS4 = GM.INS4
PAIR, INS = PB.PAIR, PB.INS
ALLREG = EXT.all_terms()

ETA = LAG.ETA
ETAD = np.diag(ETA)
TERMS = ["I1", "I1_frob", "I1_h", "J1", "J2", "Pgrad"]
LAMBDAS = (0.0, 0.5, 0.75, 1.0)
G_TOY, DELTA_TOY, S_TOY = 32.0, 0.3, -1.0
N_PAIR, L_PAIR, IT_PAIR = 32, 48.0, 120
DS = (12.0, 18.0, 24.0)
J2_C = (-0.25, -0.1, 0.1, 0.25)
PGRAD_C = (-1.0, -0.1, 0.1, 1.0)
TOL = 1e-12


def log(msg):
    print(f"[{time.time() - T_START:8.1f}s] {msg}", flush=True)


def toy_params():
    return LAG.default_params(s=S_TOY, g=G_TOY, delta=DELTA_TOY)


def merge_json(key, payload):
    old = {}
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            old = json.load(f)
    old[key] = payload
    old["elapsed_s_total"] = time.time() - T_START
    with open(OUT_JSON, "w") as f:
        json.dump(old, f, indent=1, default=float)
    return old


# ================= stage CHAN =================
def c_e_of(term, M, cfg, p, a0):
    A_, B_, C_ = LAG.omega_decompose(ALLREG[term], M, cfg, p, a0)
    return {"A": float(A_), "B": float(B_), "C": float(C_),
            "C_E": float(-4.0 * C_)}


def lambda_interval(kin1, kinh, chan_names):
    """{lambda : (1-l) kin1 + l kinh >= 0 (> 0 on boosts)} over channels."""
    lo, hi, empty = -np.inf, np.inf, False
    per = {}
    for ch in chan_names:
        k1, kh = kin1[ch], kinh[ch]
        boost = ch.split("|")[0].startswith("boost") or "boost" in ch
        slope = kh - k1                       # K(l) = k1 + l (kh - k1)
        if abs(slope) <= TOL:
            ok = (k1 > TOL) if boost else (k1 >= -TOL)
            per[ch] = "unconstrained" if ok else "EMPTY (flat, K <= 0)"
            empty = empty or (not ok)
        elif slope > 0:
            l0 = -k1 / slope
            per[ch] = f"lambda {'>' if boost else '>='} {l0:.6g}"
            lo = max(lo, l0)
        else:
            l0 = -k1 / slope
            per[ch] = f"lambda {'<' if boost else '<='} {l0:.6g}"
            hi = min(hi, l0)
    return {"per_channel": per, "lo": float(lo), "hi": float(hi),
            "nonempty": bool((not empty) and lo < hi)}


def chan_table(cfg, p, M, chans, terms, label):
    rows = {}
    for nm, a0 in chans.items():
        rows[nm] = {t: c_e_of(t, M, cfg, p, a0) for t in terms}
        rows[nm]["kin_of_B3"] = float(B3.kin_of(M, a0, cfg))
        log(f"{label} {nm}: " + ", ".join(
            f"{t} {rows[nm][t]['C_E']:+.5f}" for t in terms)
            + f", B3 {rows[nm]['kin_of_B3']:+.5f}")
    return rows


def stage_chan():
    rng = np.random.default_rng(3220)
    cfg = B3.base_cfg(s=S_TOY, g=G_TOY, n=32, L=48.0, delta=DELTA_TOY)
    p = toy_params()
    M = AUD.hedgehog_M(cfg, G_TOY)
    chans = AUD.channels(cfg, M, rng, extras=True)
    log(f"CHAN toy point: n={cfg['n']} L={cfg['L']} h={cfg['h']} "
        f"g={cfg['g']} delta={cfg['delta']} s={cfg['s']} "
        f"stencil={cfg['stencil']}; {len(chans)} channels")
    rows = chan_table(cfg, p, M, chans, TERMS, "CHAN")
    # block-diagonality of the background (decides I1_h = I1_frob here)
    offblock = float(np.max(np.abs(M[..., 0, 1:])))
    # controls: C_E(I1) = kin_of; vs the R1 audit table on the
    # deterministic channels
    with open(os.path.join(DATA, "m5_32_r1_audit_screen.json")) as f:
        aud = json.load(f)["S1_S2"]["C_E_table"]
    ctrl = {}
    for nm in rows:
        r = {"I1_vs_kin_of_rel": LAG._rel(rows[nm]["I1"]["C_E"],
                                          rows[nm]["kin_of_B3"])}
        if nm in aud and not nm.startswith("x_random"):
            r["I1_vs_R1_audit_rel"] = LAG._rel(rows[nm]["I1"]["C_E"],
                                               aud[nm]["I1"])
            r["I1_frob_vs_R1_audit_rel"] = LAG._rel(
                rows[nm]["I1_frob"]["C_E"], aud[nm]["I1_frob"])
        r["I1_h_minus_I1_frob"] = rows[nm]["I1_h"]["C_E"] \
            - rows[nm]["I1_frob"]["C_E"]
        ctrl[nm] = r
    worst_ctrl = max(max(v for k, v in r.items() if k.endswith("_rel"))
                     for r in ctrl.values())
    # the lambda family
    kin1 = {nm: rows[nm]["I1"]["C_E"] for nm in rows}
    kinh = {nm: rows[nm]["I1_h"]["C_E"] for nm in rows}
    fam = {}
    for lam in LAMBDAS:
        fam[f"lambda_{lam:g}"] = {nm: (1 - lam) * kin1[nm] + lam * kinh[nm]
                                  for nm in rows}
    tan_all = [c for c in rows if c.endswith("|tan")]
    tan_cat = [c for c in tan_all if not c.startswith("x_")]
    ivs = {"tangent_catalog": lambda_interval(kin1, kinh, tan_cat),
           "tangent_plus_extras": lambda_interval(kin1, kinh, tan_all),
           "everything_incl_probes": lambda_interval(kin1, kinh,
                                                     list(rows))}
    for k, v in ivs.items():
        log(f"CHAN lambda interval [{k}]: lo {v['lo']:.6g} hi {v['hi']:.6g}"
            f" nonempty={v['nonempty']}")
    # the stored non-vacuum 3D field (M5.21.11 end state), embedded
    cfg3 = B3.base_cfg(s=1.0)
    p3 = LAG.default_params(s=1.0, g=cfg3["g"], delta=cfg3["delta"])
    M3 = B3.embed34(LAG.load_stored3(), cfg3)
    ch3 = AUD.channels(cfg3, M3, rng, extras=False)
    log(f"CHAN stored3: {LAG.STORED3_NPZ}, cfg s={cfg3['s']} g={cfg3['g']} "
        f"delta={cfg3['delta']} n={cfg3['n']} L={cfg3['L']}")
    rows3 = chan_table(cfg3, p3, M3, ch3, TERMS, "CHAN stored3")
    kin1_3 = {nm: rows3[nm]["I1"]["C_E"] for nm in rows3}
    kinh_3 = {nm: rows3[nm]["I1_h"]["C_E"] for nm in rows3}
    fam3 = {f"lambda_{lam:g}": {nm: (1 - lam) * kin1_3[nm]
                                + lam * kinh_3[nm] for nm in rows3}
            for lam in LAMBDAS}
    ivs3 = {"tangent": lambda_interval(kin1_3, kinh_3,
                                       [c for c in rows3
                                        if c.endswith("|tan")]),
            "everything": lambda_interval(kin1_3, kinh_3, list(rows3))}
    out = {"toy_point": {k: cfg[k] for k in ("s", "g", "delta", "n", "L",
                                               "h", "stencil")},
           "background": "analytic hedgehog C14.m4h_batch(P, g) (Qh d4 "
                         "Qh^T, spatial rotation of the vacuum)",
           "background_max_offblock_time_row": offblock,
           "channels": list(rows),
           "rows": rows,
           "C_E_table": {t: {nm: rows[nm][t]["C_E"] for nm in rows}
                         for t in TERMS},
           "B_table": {t: {nm: rows[nm][t]["B"] for nm in rows}
                       for t in TERMS},
           "controls": ctrl, "controls_worst_rel": float(worst_ctrl),
           "controls_pass_1e-10": bool(worst_ctrl <= 1e-10),
           "lambda_family_K": fam, "lambda_intervals": ivs,
           "stored3": {"file": LAG.STORED3_NPZ,
                       "cfg": {k: cfg3[k] for k in ("s", "g", "delta", "n",
                                                    "L", "h", "stencil")},
                       "rows": rows3,
                       "C_E_table": {t: {nm: rows3[nm][t]["C_E"]
                                         for nm in rows3} for t in TERMS},
                       "lambda_family_K": fam3, "lambda_intervals": ivs3},
           "reading": ("K(lambda) = (1 - lambda) C_E(I1) + lambda C_E(I1_h)"
                       " per channel; C_E = -4 C, C the omega^2 "
                       "coefficient of the Lagrangian read; interval = "
                       "intersection of the per-channel half-lines "
                       "(boosts strict)")}
    merge_json("chan", out)
    return out


# ================= stage COULOMB =================
def embed_pair(M3, gtt):
    M4 = np.zeros(M3.shape[:3] + (4, 4), dtype=M3.dtype)
    M4[..., 1:, 1:] = M3
    M4[..., 0, 0] = gtt
    return M4


def i1h_pieces(A, M4):
    """(density, W over ordered pairs, local grad wrt M4) of I1_h."""
    F = LAG.F_of_A(A)
    h = EXT.h_cov_np(M4)                              # (..., 4, 4)
    dens = 0.0
    W = np.zeros_like(F)
    Dh = np.zeros(M4.shape)
    for mu in range(4):
        for nu in range(4):
            if mu == nu:
                continue
            e = ETAD[mu] * ETAD[nu]
            Fm = F[..., mu, nu, :, :]
            hFh = h @ Fm @ h
            W[..., mu, nu, :, :] = e * hFh
            if mu < nu:
                dens = dens + e * np.einsum("...ab,...ab->...", hFh, Fm)
                Dh += 2.0 * e * (Fm @ h @ Fm.swapaxes(-1, -2))
    # local piece through u(M): grad_M = -4 sym[(R^T y) w^T]
    u0 = EXT.timelike_eig_np(M4)[0]
    N = M4 @ ETA
    w = u0 @ ETA                                      # (eta u)^T
    # l0 = (u^T eta N u) / (u^T eta u) with u^T eta u = -1
    l0 = -np.einsum("...a,...a->...", w, (N @ u0[..., None])[..., 0])
    # (N - l0) du = -(dN - dl0) u, dl0 = -(w^T dN u): RHS = -(I + u w^T) dN u;
    # the u-component of du fixed by u^T eta du = 0: du = (I + u w^T) du_pinv
    Rp = np.linalg.pinv(N - l0[..., None, None] * np.eye(4))
    Pu = np.eye(4) + u0[..., :, None] * w[..., None, :]
    R = Pu @ Rp @ Pu
    y = np.einsum("ab,...bc,cd,...d->...a", ETA, Dh, ETA, u0)
    Ry = np.einsum("...ba,...b->...a", R, y)             # R^T y
    Gloc = -4.0 * Ry[..., :, None] * w[..., None, :]
    Gloc = 0.5 * (Gloc + Gloc.swapaxes(-1, -2))
    return dens, W, Gloc


def j2_pieces(A, M4):
    F = LAG.F_of_A(A)
    E = ETA
    P = E @ M4 @ E @ M4 @ E
    Q = E @ M4 @ E
    dens = 0.0
    W = np.zeros_like(F)
    Gloc = np.zeros(M4.shape)
    for mu in range(4):
        for nu in range(4):
            if mu == nu:
                continue
            e = ETAD[mu] * ETAD[nu]
            Fm = F[..., mu, nu, :, :]
            W[..., mu, nu, :, :] = -0.5 * e * (P @ Fm @ E + E @ Fm @ P)
            if mu < nu:
                Z = Fm @ E @ Fm @ E
                dens = dens + e * np.einsum("...ab,...ba->...", Z, M4 @ Q)
                ZT = Z.swapaxes(-1, -2)
                Gloc += e * (ZT @ Q + E @ M4 @ ZT @ E)
    Gloc = 0.5 * (Gloc + Gloc.swapaxes(-1, -2))
    return dens, W, Gloc


class ExtAction:
    """E = 4 h^3 sum_br wt sum [(1-lam) I1 + lam I1_h + c J2] + V_T2 on the
    3x3 pair instrument (M4 = diag(gtt) (+) M3)."""

    def __init__(self, lam=0.0, c_j2=0.0, gtt=G_TOY):
        self.lam, self.c, self.gtt = float(lam), float(c_j2), float(gtt)
        self.K1 = PB.K_of("I1")

    def _dens_and_grad(self, A, M4, need_grad):
        F = LAG.F_of_A(A)
        dens = (1.0 - self.lam) * LAG.density_from_K(F, self.K1)
        W = (1.0 - self.lam) * LAG.dW_from_K(F, self.K1) if need_grad else 0.0
        Gloc = 0.0
        if self.lam != 0.0:
            d, Wh, Gh = i1h_pieces(A, M4)
            dens = dens + self.lam * d
            if need_grad:
                W = W + self.lam * Wh
                Gloc = Gloc + self.lam * Gh
        if self.c != 0.0:
            d, Wj, Gj = j2_pieces(A, M4)
            dens = dens + self.c * d
            if need_grad:
                W = W + self.c * Wj
                Gloc = Gloc + self.c * Gj
        return dens, W, Gloc

    def e_parts(self, M, cfg, st=None):
        st = st or cfg["stencil"]
        h = cfg["h"]
        M4 = embed_pair(M, self.gtt)
        e_u = 0.0
        for br, wt in INS.branches(st):
            A = PB.jets4(M, h, br)
            dens, _, _ = self._dens_and_grad(A, M4, False)
            e_u += wt * 4.0 * float(np.sum(dens))
        e_v = float(np.sum(INS.v_density(M, cfg)))
        return h ** 3 * e_u, 0.0, h ** 3 * e_v

    def grad(self, M, cfg):
        h = cfg["h"]
        M4 = embed_pair(M, self.gtt)
        G = np.zeros_like(M)
        for br, wt in INS.branches(cfg["stencil"]):
            A = PB.jets4(M, h, br)
            _, W, Gloc = self._dens_and_grad(A, M4, True)
            dA = LAG.dA_from_W(W, A)
            for ax in range(3):
                G += wt * 4.0 * INS.d1_adj(dA[1 + ax][..., 1:, 1:],
                                           ax, h, br)
            if not np.isscalar(Gloc):
                G += wt * 4.0 * Gloc[..., 1:, 1:]
        G = 0.5 * (G + G.swapaxes(-1, -2))
        G += INS.v_grad(M, cfg)
        return (h ** 3) * G

    def install(self):
        INS.grad = self.grad
        INS.e_parts = self.e_parts


def fd_check(act, M, cfg, rng, eps=1e-4, ndir=2):
    """central-difference directional derivative vs <grad, D>."""
    g = act.grad(M, cfg)
    out = []
    for _ in range(ndir):
        D = rng.normal(size=M.shape)
        D = 0.5 * (D + D.swapaxes(-1, -2))
        D /= np.sqrt(np.sum(D * D))
        ep = sum(act.e_parts(M + eps * D, cfg))
        em = sum(act.e_parts(M - eps * D, cfg))
        fd = (ep - em) / (2 * eps)
        an = float(np.sum(g * D))
        out.append({"fd": fd, "analytic": an,
                    "rel": abs(fd - an) / max(abs(an), 1e-300)})
    return out


def fd_check_4x4(lam, c_j2, rng, n=6, eps=1e-5):
    """the same gradient formulas on a GENERIC (non-block) 4x4 lattice
    field: the local u-chain and the J2 local piece are exercised."""
    cfg = B3.base_cfg(s=S_TOY, g=G_TOY, n=n, L=9.0, delta=DELTA_TOY)
    h = cfg["h"]
    M = B3.vac4(cfg)[None, None, None] + 0.4 * B3.sym4(
        rng.normal(size=(n, n, n, 4, 4)))
    act = ExtAction(lam, c_j2)

    def energy(M):
        tot = 0.0
        for br, wt in B3.branches("sym"):
            A = np.zeros((4,) + M.shape)
            for ax in range(3):
                A[1 + ax] = B3.d1(M, ax, h, br)
            dens, _, _ = act._dens_and_grad(A, M, False)
            tot += wt * float(np.sum(dens))
        return h ** 3 * tot

    G = np.zeros_like(M)
    for br, wt in B3.branches("sym"):
        A = np.zeros((4,) + M.shape)
        for ax in range(3):
            A[1 + ax] = B3.d1(M, ax, h, br)
        _, W, Gloc = act._dens_and_grad(A, M, True)
        dA = LAG.dA_from_W(W, A)
        for ax in range(3):
            G += wt * B3.d1_adj(dA[1 + ax], ax, h, br)
        if not np.isscalar(Gloc):
            G += wt * Gloc
    G = h ** 3 * B3.sym4(G)
    out = []
    for _ in range(2):
        D = B3.sym4(rng.normal(size=M.shape))
        D /= np.sqrt(np.sum(D * D))
        fd = (energy(M + eps * D) - energy(M - eps * D)) / (2 * eps)
        an = float(np.sum(G * D))
        out.append({"fd": fd, "analytic": an,
                    "rel": abs(fd - an) / max(abs(an), 1e-300)})
    return out


def run_ext_point(lam, c_j2, kind, d, tag_prefix, save=True):
    t0 = time.time()
    cfg = PB.cfg_of()
    act = ExtAction(lam, c_j2)
    act.install()
    tag = f"{tag_prefix}_{kind}_d{d:g}"
    base = {"lambda": lam, "c_j2": c_j2, "kind": kind, "d": d,
            "it": IT_PAIR, "n": N_PAIR, "L": L_PAIR, "tag": tag}
    try:
        M0, M, info = PAIR.heal(cfg, kind, d, IT_PAIR)
        e_u, _, e_v = act.e_parts(M, cfg)
    except Exception as e:                       # noqa: BLE001
        return dict(base, status="DIVERGED", stop=f"exception: {e!r}",
                    wall_s=time.time() - t0)
    E = e_u + e_v
    row = dict(base, E=float(E), E_curv=float(e_u), E_v=float(e_v),
               stop=info["stop"], wall_s=time.time() - t0)
    finite = bool(np.isfinite(E)) and bool(np.all(np.isfinite(M))) \
        and info["stop"] != "non-finite"
    row["status"] = "OK" if finite else "DIVERGED"
    if finite:
        Es = [t["E"] for t in info["trace"]]
        row["E_trace_min"] = float(min(Es))
        ints = PB.curv_integrals(M, cfg, ("I1", "I2", "I6"))
        row["integrals"] = ints
        row["abc"] = {"a": ints["I1"] / 2.0, "b": ints["I2"] / 4.0,
                      "c_": ints["I6"] / 4.0}
        # the J2 and Pgrad integrals on the end state (footprint content)
        h = cfg["h"]
        M4 = embed_pair(M, G_TOY)
        j2 = pg = 0.0
        for br, wt in INS.branches(cfg["stencil"]):
            A = PB.jets4(M, h, br)
            j2 += wt * float(np.sum(EXT.J2_np(A, M4, None)))
            pg += wt * float(np.sum(np.abs(EXT.Pgrad_np(A, M4, None))))
        row["integrals"]["J2"] = h ** 3 * j2
        row["integrals"]["Pgrad_abs"] = h ** 3 * pg
        dq = d if d > 0 else 18.0
        row["charge"] = PAIR.charge_suite(M, cfg, dq)
        if kind == "single":
            row["r_half"] = float(INS.r_half(M, cfg))
        else:
            zs, gv = PAIR.core_zs(M, cfg, with_gaps=True)
            row["cores"], row["core_gaps"] = zs, gv
        if save:
            os.makedirs(OUT_NPZ, exist_ok=True)
            np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"),
                                M=M.astype(np.float64))
    print(json.dumps({k: v for k, v in row.items()}, default=float),
          flush=True)
    return row


def _worker(args):
    return run_ext_point(*args)


def fit_rows(rows, lam, c):
    """E_int(d) fits for one action (lambda, c)."""
    sel = [r for r in rows if r["lambda"] == lam and r["c_j2"] == c]
    single = [r for r in sel if r["kind"] == "single"]
    out = {"lambda": lam, "c_j2": c, "n_rows": len(sel),
           "diverged": [r["tag"] for r in sel if r["status"] != "OK"]}
    if not single or single[0]["status"] != "OK":
        out["status"] = "NO SINGLE"
        return out
    s = single[0]
    out["E_single"] = s["E"]
    out["abc_single"] = s["abc"]
    out["integrals_single"] = s["integrals"]
    out["r_half"] = s.get("r_half")
    out["degree_single_far"] = s["charge"]["far"]
    for kind in ("same", "anti"):
        pts = sorted([r for r in sel if r["kind"] == kind
                      and r["status"] == "OK"], key=lambda r: r["d"])
        ds = [r["d"] for r in pts]
        eint = [r["E"] - 2 * s["E"] for r in pts]
        rec = {"d": ds, "E_int": eint,
               "far_flux": [r["charge"]["far"] for r in pts]}
        if len(ds) >= 3:
            A_, B_, r2 = PB.fit_inv_d(ds, eint)
            rec.update({"A": A_, "B": B_, "R2": r2,
                        "trend": "E_int falls with d (repulsive)"
                        if eint[-1] < eint[0] else
                        "E_int rises with d (attractive)"})
        out[kind] = rec
    out["status"] = "OK"
    return out


def stage_coulomb(workers=12):
    t0 = time.time()
    rng = np.random.default_rng(3221)
    cfg = PB.cfg_of()
    out = {"instrument": {"n": N_PAIR, "L": L_PAIR, "h": cfg["h"],
                          "it": IT_PAIR, "stencil": cfg["stencil"],
                          "bc": cfg["bc"], "term": cfg["term"],
                          "g_tt_embed": G_TOY, "delta": cfg["delta"]}}
    # (0) gradient checks: generic 4x4 field, then the 3x3 seeds
    gc = {"generic_4x4_n6_L9": {
        "lambda_1": fd_check_4x4(1.0, 0.0, rng),
        "lambda_0.5": fd_check_4x4(0.5, 0.0, rng),
        "J2_c0.25": fd_check_4x4(0.0, 0.25, rng),
        "J2_c-0.25_lambda_0.5": fd_check_4x4(0.5, -0.25, rng)}}
    gc["pair_seeds_3x3"] = {}
    for kind, d in (("single", 0.0), ("same", 18.0)):
        M0 = PAIR.seed_pair(cfg, kind, d)
        for lab, act in (("lambda_1", ExtAction(1.0, 0.0)),
                         ("J2_c0.25", ExtAction(0.0, 0.25)),
                         ("J2_c-0.1", ExtAction(0.0, -0.1))):
            gc["pair_seeds_3x3"][f"{kind}_d{d:g}|{lab}"] = fd_check(
                act, M0, cfg, rng)
        # the identity at the seed: lambda = 1 vs lambda = 0 energy + grad
        a0, a1 = ExtAction(0.0), ExtAction(1.0)
        g0, g1 = a0.grad(M0, cfg), a1.grad(M0, cfg)
        e0, e1 = a0.e_parts(M0, cfg), a1.e_parts(M0, cfg)
        gcert = PB.ModifiedAction("I2", 0.0).grad(M0, cfg)
        gc["pair_seeds_3x3"][f"{kind}_d{d:g}|identity_seed"] = {
            "E_u_lambda0": e0[0], "E_u_lambda1": e1[0],
            "abs_diff_E_u": abs(e0[0] - e1[0]),
            "grad_max_abs_diff_lambda1_vs_0": float(np.max(np.abs(g1 - g0))),
            "grad_max_abs_diff_lambda0_vs_r1b": float(
                np.max(np.abs(g0 - gcert))),
            "grad_max_abs": float(np.max(np.abs(g0)))}
    worst_fd = max(x["rel"] for grp in gc.values() for v in grp.values()
                   if isinstance(v, list) for x in v)
    gc["worst_fd_rel"] = float(worst_fd)
    out["gradient_check"] = gc
    log(f"COULOMB gradient checks: worst fd rel {worst_fd:.2e}")
    # (1) the jobs: control lambda = 0, identity lambda = 1, J2 footprint
    configs = [("single", 0.0)] + [(k, d) for k in ("same", "anti")
                                   for d in DS]
    jobs = []
    for kind, d in configs:
        jobs.append((0.0, 0.0, kind, d, "lam0"))
        jobs.append((1.0, 0.0, kind, d, "lam1"))
        for c in J2_C:
            jobs.append((0.0, c, kind, d, f"J2_c{c:+g}"))
    log(f"COULOMB {len(jobs)} relaxations on {workers} workers")
    rows = []
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        futs = {ex.submit(_worker, j): j for j in jobs}
        done = 0
        for fut in as_completed(futs):
            j = futs[fut]
            try:
                rows.append(fut.result())
            except Exception as e:                      # noqa: BLE001
                rows.append({"lambda": j[0], "c_j2": j[1], "kind": j[2],
                             "d": j[3], "tag": f"{j[4]}_{j[2]}_d{j[3]:g}",
                             "status": "DIVERGED",
                             "stop": f"exception: {e!r}"})
            done += 1
            if done % 10 == 0 or done == len(jobs):
                log(f"COULOMB [{done}/{len(jobs)}]")
    out["rows"] = rows
    # (2) the identity lambda = 1 vs lambda = 0, and vs the record
    with open(os.path.join(DATA, "m5_21_4_ladder_it120.json")) as f:
        rec = json.load(f)
    ident = {}
    for kind, d in configs:
        r0 = [r for r in rows if r["lambda"] == 0.0 and r["c_j2"] == 0.0
              and r["kind"] == kind and r["d"] == d][0]
        r1 = [r for r in rows if r["lambda"] == 1.0 and r["c_j2"] == 0.0
              and r["kind"] == kind and r["d"] == d][0]
        key = f"{kind}_d{d:g}"
        e = {"E_lambda0": r0.get("E"), "E_lambda1": r1.get("E"),
             "abs_diff_E": abs(r0["E"] - r1["E"])
             if r0["status"] == "OK" and r1["status"] == "OK" else None}
        try:
            Ma = np.load(os.path.join(OUT_NPZ, r0["tag"] + ".npz"))["M"]
            Mb = np.load(os.path.join(OUT_NPZ, r1["tag"] + ".npz"))["M"]
            e["field_max_abs_diff"] = float(np.max(np.abs(Ma - Mb)))
        except Exception:                                # noqa: BLE001
            e["field_max_abs_diff"] = None
        if kind == "single":
            e["E_record"] = rec["E_single"]
        else:
            rr = [x for x in rec["rows"] if x["kind"] == kind
                  and x["d"] == d]
            e["E_record"] = rr[0]["E"] if rr else None
        if e["E_record"] is not None and r0["status"] == "OK":
            e["abs_diff_lambda0_vs_record"] = abs(r0["E"] - e["E_record"])
        ident[key] = e
    ident["max_abs_diff_E_lambda1_vs_0"] = max(
        v["abs_diff_E"] for k, v in ident.items()
        if isinstance(v, dict) and v.get("abs_diff_E") is not None)
    ident["max_field_diff"] = max(
        v["field_max_abs_diff"] for k, v in ident.items()
        if isinstance(v, dict) and v.get("field_max_abs_diff") is not None)
    ident["max_abs_diff_lambda0_vs_record"] = max(
        v["abs_diff_lambda0_vs_record"] for k, v in ident.items()
        if isinstance(v, dict) and "abs_diff_lambda0_vs_record" in v)
    ident["pass_1e-10"] = bool(ident["max_abs_diff_E_lambda1_vs_0"] <= 1e-10)
    out["identity_lambda1"] = ident
    log(f"COULOMB identity: max|E(1) - E(0)| = "
        f"{ident['max_abs_diff_E_lambda1_vs_0']:.3e}, field "
        f"{ident['max_field_diff']:.3e}, lambda0 vs record "
        f"{ident['max_abs_diff_lambda0_vs_record']:.3e}")
    # (3) footprint fits
    fits = {"control_lambda0": fit_rows(rows, 0.0, 0.0),
            "lambda1": fit_rows(rows, 1.0, 0.0)}
    for c in J2_C:
        fits[f"J2_c{c:+g}"] = fit_rows(rows, 0.0, c)
    ctrl = fits["control_lambda0"]
    for k, f_ in fits.items():
        if f_.get("status") != "OK":
            continue
        rel = {}
        for kind in ("same", "anti"):
            if "B" in f_.get(kind, {}) and "B" in ctrl.get(kind, {}):
                rel[kind] = {
                    "B": f_[kind]["B"], "B_control": ctrl[kind]["B"],
                    "same_sign_as_control": bool(
                        np.sign(f_[kind]["B"]) == np.sign(ctrl[kind]["B"])),
                    "E_int_ratio_to_control": [
                        a / b if abs(b) > 1e-300 else None
                        for a, b in zip(f_[kind]["E_int"],
                                        ctrl[kind]["E_int"])]}
        f_["relative_to_control"] = rel
        log(f"COULOMB fit {k}: same B {f_['same'].get('B')}, anti B "
            f"{f_['anti'].get('B')}, E_single {f_['E_single']:.6f}, "
            f"diverged {f_['diverged']}")
    out["fits"] = fits
    # (4) Pgrad on the static sector: density and directional derivative
    pg = {}
    for kind, d in (("single", 0.0), ("same", 18.0), ("anti", 24.0)):
        r0 = [r for r in rows if r["lambda"] == 0.0 and r["c_j2"] == 0.0
              and r["kind"] == kind and r["d"] == d][0]
        M = np.load(os.path.join(OUT_NPZ, r0["tag"] + ".npz"))["M"]
        h = cfg["h"]
        M4 = embed_pair(M, G_TOY)
        dmax = 0.0
        for br, wt in INS.branches(cfg["stencil"]):
            A = PB.jets4(M, h, br)
            dmax = max(dmax, float(np.max(np.abs(EXT.Pgrad_np(A, M4, None)))))
        # directional derivative of the Pgrad lattice integral (3x3 dirs)
        D = rng.normal(size=M.shape)
        D = 0.5 * (D + D.swapaxes(-1, -2))
        D /= np.sqrt(np.sum(D * D))

        def pg_int(Mx):
            tot = 0.0
            for br, wt in INS.branches(cfg["stencil"]):
                A = PB.jets4(Mx, h, br)
                tot += wt * float(np.sum(EXT.Pgrad_np(A, embed_pair(Mx, G_TOY),
                                                      None)))
            return h ** 3 * tot

        eps = 1e-4
        dd = (pg_int(M + eps * D) - pg_int(M - eps * D)) / (2 * eps)
        pg[f"{kind}_d{d:g}"] = {"max_abs_density": dmax,
                                "directional_derivative": dd,
                                "integral": pg_int(M)}
    pg["reading"] = ("Pgrad = sum eta^mu q(d_mu P_t, d_mu P_t) is quadratic "
                     "in d_mu u; on the block-diagonal static sector d_mu u "
                     "= 0 identically (the off-block jets vanish), so the "
                     "density, the lattice integral and every first "
                     "variation within the 3x3 sector vanish: the Pgrad "
                     "footprint on this instrument is IDENTICAL to the "
                     "control at every c (no relaxation can move), and "
                     "the term is UNDECIDABLE here (needs 4x4 dof)")
    pg["c_grid_declared"] = list(PGRAD_C)
    out["pgrad_static_sector"] = pg
    out["elapsed_s"] = time.time() - t0
    merge_json("coulomb", out)
    return out


# ================= stage GARM =================
def i1h_lattice(M4, cfg):
    """4 h^3 sum_br wt sum I1_h (the static curvature energy of I1_h)."""
    h = cfg["h"]
    tot = 0.0
    for br, wt in INS4.branches(cfg["stencil"]):
        A = np.zeros((4,) + M4.shape)
        for ax in range(3):
            A[1 + ax] = INS4.d1(M4, ax, h, br)
        tot += wt * float(np.sum(EXT.I1_h_np(A, M4, None)))
    return 4.0 * h ** 3 * tot


def e_pieces(M4, cfg, m):
    Qb = GM.qb_field(cfg, m)
    Md = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, M4, Qb))
    e_u, e_v = INS4.e_parts(Md, cfg)
    return float(e_u), float(e_v), i1h_lattice(Md, cfg)


def garm_arm(branch, g, nm=121, span=3.0, lambdas=(0.0, 0.5, 1.0)):
    t0 = time.time()
    tag = f"t11lad_{branch}_n{GM.N}_d{GM.DELTA:g}"
    M3 = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))["M"] \
        .astype(np.float64)
    cfg = INS4.base_cfg(s=-1.0, g=g, n=GM.N, L=48.0, delta=GM.DELTA)
    M4 = INS4.embed34(M3, cfg)
    m_his = float(np.arctanh(1.0 / g))
    ms = np.linspace(-span * m_his, span * m_his, nm)
    pieces = np.array([e_pieces(M4, cfg, m) for m in ms])   # (nm, 3)
    p0 = np.array(e_pieces(M4, cfg, 0.0))
    out = {"branch": branch, "g": g, "s": -1.0, "n": GM.N, "L": 48.0,
           "h": cfg["h"], "delta": GM.DELTA, "stencil": cfg["stencil"],
           "nm": nm, "span": span, "m_his": m_his,
           "E0_pieces": {"E_u": p0[0], "E_V": p0[1], "E_I1h": p0[2]},
           "per_lambda": {}, "wall_s": None}
    for lam in lambdas:
        Es = (1 - lam) * pieces[:, 0] + lam * pieces[:, 2] + pieces[:, 1]
        E0 = (1 - lam) * p0[0] + lam * p0[2] + p0[1]
        i = int(np.argmin(Es))
        m_star, E_star, edge = float(ms[i]), float(Es[i]), False
        if 0 < i < nm - 1:
            a, b, c = Es[i - 1], Es[i], Es[i + 1]
            dm = ms[1] - ms[0]
            den = c - 2 * b + a
            if den > 0:
                m_star = float(ms[i] - 0.5 * dm * (c - a) / den)
                E_star = float(b - 0.125 * (c - a) ** 2 / den)
        else:
            edge = True
        out["per_lambda"][f"lambda_{lam:g}"] = {
            "E0": float(E0), "m_star": m_star, "E_star": E_star,
            "gain": float(E_star - E0), "edge_minimum": edge,
            "E_edges": [float(Es[0]), float(Es[-1])],
            "curve": [{"m": float(m), "E": float(E)}
                      for m, E in zip(ms[::6], Es[::6])]}
    out["wall_s"] = time.time() - t0
    print(json.dumps({k: out[k] for k in ("branch", "g", "m_his")}
                     | {k: (v["gain"], v["m_star"], v["edge_minimum"])
                        for k, v in out["per_lambda"].items()}),
          flush=True)
    return out


def _garm_worker(args):
    return garm_arm(*args)


def q_fit(arms, lam, gs):
    rows = [a for a in arms if a["g"] in gs]
    rows.sort(key=lambda a: a["g"])
    gains = [a["per_lambda"][f"lambda_{lam:g}"]["gain"] for a in rows]
    edge = [a["per_lambda"][f"lambda_{lam:g}"]["edge_minimum"] for a in rows]
    rec = {"g": [a["g"] for a in rows], "gain": gains, "edge_minimum": edge}
    if all(gn < 0 for gn in gains):
        x = np.log([np.arctanh(1.0 / a["g"]) for a in rows])
        y = np.log([-gn for gn in gains])
        slope, icpt = np.polyfit(x, y, 1)
        rec["q_lsq"] = float(slope)
        rec["pair_slopes"] = [float((y[i + 1] - y[i]) / (x[i + 1] - x[i]))
                              for i in range(len(x) - 1)]
        rec["status"] = "FIT" + (" (edge minima present)" if any(edge)
                                 else "")
    else:
        rec["q_lsq"] = None
        rec["status"] = "NO GAIN on some g (gain >= 0: the dressing does " \
                        "not lower the energy; q undefined)"
    return rec


def stage_garm(workers=12):
    t0 = time.time()
    gs = (8.0, 16.0, 32.0, 128.0)
    jobs = [(br, g) for br in GM.BRANCHES for g in gs]
    log(f"GARM {len(jobs)} arms on {workers} workers (n = {GM.N}, "
        f"121 m-points each, lambda = 0 / 0.5 / 1)")
    arms = []
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs)),
                             mp_context=ctx) as ex:
        for fut in as_completed([ex.submit(_garm_worker, j) for j in jobs]):
            arms.append(fut.result())
    arms.sort(key=lambda a: (GM.BRANCHES.index(a["branch"]), a["g"]))
    fits = {}
    for br in GM.BRANCHES:
        fits[br] = {}
        sub = [a for a in arms if a["branch"] == br]
        for lam in (0.0, 0.5, 1.0):
            fits[br][f"lambda_{lam:g}"] = {
                "record_set_8_16_32": q_fit(sub, lam, (8.0, 16.0, 32.0)),
                "full_set_8_16_32_128": q_fit(sub, lam, gs)}
    # C3: lambda = 0 vs the record garm.json
    with open(os.path.join(DATA, "m5_21_11_garm.json")) as f:
        rec = json.load(f)
    c3 = []
    for a in arms:
        rr = [r for r in rec["arms"] if r["branch"] == a["branch"]
              and r["g"] == a["g"] and r["s"] == -1.0]
        if rr:
            l0 = a["per_lambda"]["lambda_0"]
            c3.append({"branch": a["branch"], "g": a["g"],
                       "gain_here": l0["gain"], "gain_record": rr[0]["gain"],
                       "abs_diff_gain": abs(l0["gain"] - rr[0]["gain"]),
                       "abs_diff_E0": abs(l0["E0"] - rr[0]["E0"]),
                       "m_star_here": l0["m_star"],
                       "m_star_record": rr[0]["m_star"]})
    c3_worst = max(x["abs_diff_gain"] for x in c3)
    rec_q = {br: rec["fits"][br]["q_lsq"] for br in GM.BRANCHES}
    # C1 vacuum null under every lambda
    c1 = []
    for g in (8.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0, delta=0.3)
        n = cfg["n"]
        M4 = np.zeros((n, n, n, 4, 4))
        M4[:] = INS4.vac4(cfg)
        mh = float(np.arctanh(1.0 / g))
        ms = np.linspace(-2.5 * mh, 2.5 * mh, 21)
        pieces = np.array([e_pieces(M4, cfg, m) for m in ms])
        for lam in (0.0, 0.5, 1.0):
            Es = (1 - lam) * pieces[:, 0] + lam * pieces[:, 2] + pieces[:, 1]
            i = int(np.argmin(Es))
            c1.append({"g": g, "lambda": lam, "n": n, "L": 48.0,
                       "m_min": float(ms[i]), "E_min": float(Es[i]),
                       "E_edge": float(Es[0]),
                       "E_at_m0": float(Es[len(ms) // 2]),
                       "pass": bool(abs(ms[i]) < 1e-12 and Es[i] < 1e-10
                                    and Es[0] > 1.0)})
            log(f"GARM C1 g={g} lambda={lam}: m_min {ms[i]:+.4f} "
                f"E_min {Es[i]:.4e} E_edge {Es[0]:.4e} pass={c1[-1]['pass']}")
    # C2 field identity (lambda-independent)
    c2 = []
    for g in (8.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0)
        base = L8.dressed(cfg, 0.0)
        Mt = L8.dressed(cfg, -0.025)
        Qb = GM.qb_field(cfg, -0.025)
        Mm = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, base, Qb))
        c2.append({"g": g, "m": -0.025,
                   "field_diff_max": float(np.abs(Mt - Mm).max())})
    out = {"definition": ("GAIN(g) = E_lambda(m*) - E_lambda(0), m* = argmin "
                          "over 121 points of [-3, 3] artanh(1/g) with the "
                          "record's parabolic refinement; q = slope of "
                          "log(-GAIN) vs log(artanh(1/g)); E_lambda(m) = "
                          "(1-lambda) E_u + lambda 4 h^3 sum I1_h + E_V on "
                          "the rigid Qb(m) dressing of the embedded N = 48 "
                          "delta = 0.3 endpoint"),
           "arms": arms, "fits": fits,
           "record_q_lsq": rec_q, "record_f4_pass_all": rec["f4_pass_all"],
           "C3_record_match": {"rows": c3, "worst_abs_diff_gain": c3_worst,
                               "pass_1e-10": bool(c3_worst <= 1e-10)},
           "C1_vacuum_null": c1,
           "C2_field_identity": c2,
           "elapsed_s": time.time() - t0}
    for br in GM.BRANCHES:
        for lam in (0.0, 0.5, 1.0):
            f_ = fits[br][f"lambda_{lam:g}"]["record_set_8_16_32"]
            log(f"GARM {br} lambda={lam}: gains {f_['gain']} q {f_['q_lsq']} "
                f"{f_['status']}")
    merge_json("garm", out)
    return out


# ================= stage INV =================
def stage_inv():
    from scipy.ndimage import gaussian_filter
    t0 = time.time()
    rng2 = np.random.default_rng(3203)
    p2 = LAG.default_params(s=S_TOY, g=G_TOY, delta=DELTA_TOY)
    cfg2 = B3.base_cfg(s=S_TOY, g=G_TOY, n=14, L=21.0, delta=DELTA_TOY)
    MB2 = np.stack([[gaussian_filter(rng2.normal(size=(14,) * 3), 2.0)
                     for _ in range(4)] for _ in range(4)], axis=-1)
    MB2 = MB2.reshape(14, 14, 14, 4, 4)
    MB2 = B3.vac4(cfg2)[None, None, None] + 0.5 * B3.sym4(MB2)
    jets = LAG.lattice_jets(MB2, cfg2)
    K_I1_noeta = LAG._K_from_pattern("mnab", "mnab", 0.5, deriv_metric="delta")
    K_I2 = LAG.REGISTRY["I2"]._K()
    ext_names = ["I1_h", "J1", "J2", "Pgrad"]

    def dens_sum(t, A, M):
        if t == "I1_deriv_delta":
            return float(np.sum(LAG.density_from_K(LAG.F_of_A(A), K_I1_noeta)))
        if t == "I2_plain_bracket":
            P = np.einsum("m...ab,n...bc->...mnac", A, A, optimize=True)
            F = P - P.swapaxes(-4, -3)
            return float(np.sum(LAG.density_from_K(F, K_I2)))
        if t.endswith("|noeta"):
            return float(np.sum(EXT.REGISTRY_EXT[t[:-6]].density(A, M, p2)))
        return float(np.sum(ALLREG[t].density(A, M, p2)))

    names = ["I1"] + ext_names + ["I1_frob", "I1_deriv_delta",
                                  "I2_plain_bracket"] \
        + [t + "|noeta" for t in ext_names]
    Ls = []
    for kind in ("boost", "rotation"):
        for k in range(2):
            Ls.append((f"{kind}_{k}", LAG._lorentz(rng2, kind)))
    Ls.append(("mixed_0", LAG._lorentz(rng2, "boost", 0.4)
               @ LAG._lorentz(rng2, "rotation", 0.4)))
    trans = {}
    for key, L in Ls:
        Linv_T = np.linalg.inv(L).T
        Mp = np.einsum("ab,...bc,dc->...ad", L, MB2, L)
        for t in names:
            if t.endswith("|noeta"):
                EXT.set_mutant("eta_time_row")
            e0, e1 = 0.0, 0.0
            for A, wt in jets:
                Ap = np.einsum("mn,n...ab->m...ab", Linv_T,
                               np.einsum("ab,n...bc,dc->n...ad", L, A, L))
                e0 += wt * dens_sum(t, A, MB2)
                e1 += wt * dens_sum(t, Ap, Mp)
            if t.endswith("|noeta"):
                EXT.ETA_X = np.diag([-1.0, 1.0, 1.0, 1.0])
                EXT.XI_X = LAG.sp.diag(-1, 1, 1, 1)
                EXT.MUTANT = None
            trans.setdefault(t, {})[key] = {
                "drift": float(abs(e1 - e0) / max(abs(e0), 1e-300)),
                "E0": e0, "E1": e1}
        log(f"INV {key}: " + ", ".join(
            f"{t} {trans[t][key]['drift']:.1e}" for t in names))
    cov = ["I1"] + ext_names
    ctrl = [t for t in names if t not in cov]
    worst_cov = max(max(v["drift"] for v in trans[t].values()) for t in cov)
    min_ctrl_boost = min(min(v["drift"] for k, v in trans[t].items()
                             if k.startswith("boost") or k.startswith("mixed"))
                         for t in ctrl)
    out = {"drifts": trans, "covariant_terms": cov, "control_terms": ctrl,
           "worst_covariant": float(worst_cov),
           "min_control_boost_drift": float(min_ctrl_boost),
           "pass_covariant_1e-10": bool(worst_cov <= 1e-10),
           "controls_fail_under_boost": bool(min_ctrl_boost > 1e-3),
           "field": "random smooth 4x4 field, n = 14, L = 21, h = 1.5, "
                    "rng 3203, vacuum(s=-1, g=32, delta=0.3) + 0.5 "
                    "sym4(gaussian_filter(normal, 2)); sym stencil",
           "transform": "A'_mu = (Lambda^{-T})_mu^nu Lambda A_nu Lambda^T, "
                        "M' = Lambda M Lambda^T (u' = Lambda u follows "
                        "from N' = Lambda N Lambda^{-1})",
           "controls": "the ext registry's eta_time_row mutant of each new "
                       "term (contraction metric -> identity), I1 with a "
                       "delta derivative pair, I2 on the plain commutator",
           "elapsed_s": time.time() - t0}
    log(f"INV worst covariant {worst_cov:.2e}, min control boost drift "
        f"{min_ctrl_boost:.2e}")
    merge_json("inv", out)
    return out


# ================= plot =================
def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with open(OUT_JSON) as f:
        d = json.load(f)
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    # (a) K(lambda) per tangent + extra channel
    ax = axs[0]
    if "chan" in d:
        c = d["chan"]
        lams = np.linspace(-0.25, 1.5, 8)
        for nm in c["channels"]:
            if not nm.endswith("|tan"):
                continue
            k1 = c["C_E_table"]["I1"][nm]
            kh = c["C_E_table"]["I1_h"][nm]
            K = (1 - lams) * k1 + lams * kh
            boost = "boost" in nm
            ax.plot(lams, K, "-" if boost else "--", lw=2 if boost else 1,
                    label=nm.split("|")[0])
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(0.5, color="gray", ls=":", lw=0.8)
        iv = c["lambda_intervals"]["tangent_plus_extras"]
        ax.set_title(f"CHAN: K(lambda) tangent+extra, toy point\n"
                     f"stable+reversed for lambda in ({iv['lo']:.3g}, "
                     f"{iv['hi']:.3g})", fontsize=9)
        ax.set_xlabel("lambda"); ax.set_ylabel("K = C_E (energy kin coeff)")
        ax.legend(fontsize=6, ncol=2)
    # (b) footprint E_int(d)
    ax = axs[1]
    if "coulomb" in d:
        fits = d["coulomb"]["fits"]
        for k, f_ in fits.items():
            if f_.get("status") != "OK":
                continue
            for kind, mk in (("same", "o-"), ("anti", "s--")):
                r = f_.get(kind, {})
                if "E_int" in r:
                    ax.plot(r["d"], r["E_int"], mk, label=f"{k} {kind}",
                            lw=2 if k.startswith("control") else 1)
        ax.set_xlabel("d"); ax.set_ylabel("E_int = E_pair - 2 E_single")
        ax.set_title("COULOMB pair instrument (n=32, L=48, it=120)\n"
                     "lambda=0 control, lambda=1, J2 footprint", fontsize=9)
        ax.legend(fontsize=6, ncol=2)
    # (c) g-arm gain vs g
    ax = axs[2]
    if "garm" in d:
        for br in ("A", "C", "B"):
            for lam, mk in ((0.0, "o-"), (0.5, "s--"), (1.0, "^:")):
                sub = [a for a in d["garm"]["arms"] if a["branch"] == br]
                gs = [a["g"] for a in sub]
                gn = [a["per_lambda"][f"lambda_{lam:g}"]["gain"] for a in sub]
                ax.plot(gs, gn, mk, label=f"{br} lambda={lam:g}")
        ax.set_xscale("log"); ax.set_xlabel("g")
        ax.set_ylabel("GAIN = E(m*) - E(0)")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title("g-arm (n=48, delta=0.3): boost-dressing gain\n"
                     "lambda = 0 (record) / 0.5 / 1", fontsize=9)
        ax.legend(fontsize=6, ncol=3)
    fig.suptitle("M5.32 R2.c lattice screen of the C2 sign map", fontsize=11)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(PLOT, dpi=130)
    log(f"plot -> {PLOT}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "chan", "coulomb", "garm", "inv", "plot"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    if a.stage in ("all", "chan"):
        stage_chan()
    if a.stage in ("all", "inv"):
        stage_inv()
    if a.stage in ("all", "coulomb"):
        stage_coulomb(a.workers)
    if a.stage in ("all", "garm"):
        stage_garm(a.workers)
    if a.stage in ("all", "plot") or a.plot:
        plot()
    log("done")


if __name__ == "__main__":
    main()
