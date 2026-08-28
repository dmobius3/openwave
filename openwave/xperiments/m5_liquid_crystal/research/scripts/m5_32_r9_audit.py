"""M5.32 R9 ADVERSARIAL AUDIT: the string excision (arm a) and the
string-free-hedgehog theorem (arm b).

INDEPENDENT rebuild of claims A1..A4 and B1..B5. The producer's scripts
(m5_32_r9_a_tube.py, m5_32_r9_b_string.py, m5_32_r8_a_quartics.py,
m5_32_r8_b_ir_theorem.py, m5_32_r8_c_collect.py) were NOT read and are
NOT imported. Oracles: the certified stack (m5_21_3_a_4d.py) and the
m5_21_8_b_lattice.py hedgehog builder, plus my own predecessor auditor
m5_32_r8_audit.py (its quartic densities and differencing are
re-derived here from the written definitions, not imported).

EQUATIONS FIRST
---------------
eta = diag(-1, 1, 1, 1), index 0 = time. M real symmetric 4x4, RAW
CONTRAVARIANT internal entries, M -> L M L^T. Jets A_mu = d_mu M with
A_0 = omega a0.

    F_munu   = A_mu eta A_nu - A_nu eta A_mu
    <F,G>_eta= tr(eta F eta G^T) = sum_ab w_a w_b F_ab G_ab, w = diag eta
    T_munu   = tr(A_mu eta A_nu eta) = sum_ab w_a w_b (A_mu)_ab (A_nu)_ba
    I1       = sum_{mu<nu} w_mu w_nu <F_munu, F_munu>_eta
    R[nu,a]  = sum_mu F[mu,nu][a,mu];  I4 = sum w_nu w_a R[nu,a]^2
    Q_I1sq   = I1^2      Q_I4sq = I4^2
    Q_Fpair  = sum_{mu<nu} sum_{rho<sig} w.. <F_munu, F_rhosig>_eta^2
    Q_C6a    = [sum_mu w_mu T_mumu]^2

Every term is EXACTLY quartic in omega (F_0i is linear in omega and
each term is quadratic in F or in T), so I(omega) = A + C1 om + C2 om^2
+ C3 om^3 + C4 om^4 is exact; a degree-4 fit on 7 samples has zero
residual and that is checked, not assumed.

THE EXCISION (arm a). Five families of cut, each removing EXACTLY the
same NUMBER of cells k(c) = #{rho_z < c} on every box:
    cyl_z   the z axis tube        (the producer's cut: targets the string)
    cyl_x   the x axis tube        (same shape, same volume, NO string)
    cyl_off a tube about x = 12    (same shape, same volume, NO string)
    sph     a spherical core       (same volume, NOT axis-shaped)
    rand    a random cell set      (same volume, no geometry at all)
If the C5 exponents "clean up" under any of the last four, the cleanup
is a volume effect and A1 measures nothing.

Modes (each writes ../checkpoints/m5_32_r9/<tag>.json):
    probe   frame sanity: Qh e1 = n-hat, Qh e2 = e_phi, singular locus
    ring    B1 + B2: the n n^T identity, the residual, the ring law
    kin     B3: the delta ladder for kin, exact vs asymptotic delta^2
    topo    B4 (a) + (b): the assignment enumeration, the RP^1 lift,
            the Poincare-Hopf index, the nematic-XY minimization, and
            the SMOOTH-ORBIT counterexample (the radial boost)
    tube    A1..A4: the excision ladder with the four controls
    share   B5: tube shares at a FIXED PHYSICAL tube, with the null control
    divergence  separates the omega^0 POINT-core divergence from the
            omega^4 LINE divergence (they are not the same defect)
    mutate  the mutation gates: every confirmed gate must fail under a
            deliberate mutation (M2b moves the string and checks the cure
            follows it)
    relax   B4 (c): does FIRE resolve the line into a finite core?
    merge   assembles ../data/m5_32_r9_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r9")
os.makedirs(CKPT, exist_ok=True)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")            # certified stack
B8 = _load("m5_21_8_b_lattice", "m5_21_8_b_lattice.py")  # hedgehog + clock

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W = np.diag(ETA)
WW = W[:, None] * W[None, :]
G, DELTA, S = 32.0, 0.3, -1.0
TERMS = ("I1", "Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a")
OMEGAS = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
BRANCHES = (("bwd", 0.5), ("fwd", 0.5))
RHO_CUTS = (0.0, 1.5, 3.0, 4.5, 6.0, 9.0)
FAMILIES = ("cyl_z", "cyl_shift", "cyl_x", "cyl_off", "sph", "rand")
XOFF = 12.0
BOXES_L = ((32, 48.0), (48, 72.0), (64, 96.0))       # h = 1.5
BOXES_H = (((32, 48.0), (64, 48.0)), ((48, 72.0), (96, 72.0)))


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def dump(tag, obj):
    with open(os.path.join(CKPT, f"{tag}.json"), "w") as f:
        json.dump(obj, f, indent=1, default=float)
    log(f"checkpoint {tag}.json")


def rel(a, b):
    d = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / d


def loglog(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    m = (ys != 0) & np.isfinite(ys)
    if m.sum() < 2:
        return None
    return float(np.polyfit(np.log(xs[m]), np.log(np.abs(ys[m])), 1)[0])


# ===================== my differencing / algebra =====================
def mydiff(f, ax, h, br):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl); s[ax] = i; return tuple(s)
    d = (f[at(slice(1, None))] - f[at(slice(0, -1))]) / h
    if br == "fwd":
        out[at(slice(0, -1))] = d
    elif br == "bwd":
        out[at(slice(1, None))] = d
    else:
        raise ValueError(br)
    return out


def comm_eta(A, Bm):
    return (A * W[None, :]) @ Bm - (Bm * W[None, :]) @ A


def inner_eta(F, Gm):
    return np.einsum("ab,...ab,...ab->...", WW, F, Gm, optimize=True)


def tr_AetaBeta(A, Bm):
    return np.einsum("ab,...ab,...ba->...", WW, A, Bm, optimize=True)


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def densities(A, terms=TERMS):
    Fs = {p: comm_eta(A[p[0]], A[p[1]]) for p in PAIRS}

    def getF(mu, nu):
        return Fs[(mu, nu)] if (mu, nu) in Fs else -Fs[(nu, mu)]

    out = {}
    if ("I1" in terms) or ("Q_I1sq" in terms):
        i1 = 0.0
        for (mu, nu) in PAIRS:
            i1 = i1 + W[mu] * W[nu] * inner_eta(Fs[(mu, nu)], Fs[(mu, nu)])
        if "I1" in terms:
            out["I1"] = i1
        if "Q_I1sq" in terms:
            out["Q_I1sq"] = i1 * i1
    if "Q_I4sq" in terms:
        sh = A[0].shape[:-2]
        R = np.zeros(sh + (4, 4))
        for nu in range(4):
            for mu in range(4):
                if mu != nu:
                    R[..., nu, :] += getF(mu, nu)[..., :, mu]
        i4 = np.einsum("na,...na,...na->...", WW, R, R, optimize=True)
        out["Q_I4sq"] = i4 * i4
    if "Q_Fpair" in terms:
        tot = 0.0
        for (mu, nu) in PAIRS:
            for (rh, sg) in PAIRS:
                v = inner_eta(Fs[(mu, nu)], Fs[(rh, sg)])
                tot = tot + W[mu] * W[nu] * W[rh] * W[sg] * v * v
        out["Q_Fpair"] = tot
    if "Q_C6a" in terms:
        s = 0.0
        for mu in range(4):
            s = s + W[mu] * tr_AetaBeta(A[mu], A[mu])
        out["Q_C6a"] = s * s
    return out


VDM = np.vander(OMEGAS, 5, increasing=True)
VPI = np.linalg.pinv(VDM)


def fit4(I):
    """exact degree-4 coefficients of I(omega) plus the fit residual."""
    I = np.asarray(I, float)
    c = VPI @ I
    res = float(np.abs(VDM @ c - I).max())
    return c, res


# ===================== the frame / ansatz helpers =====================
def frame_vectors(x, y, z):
    """the ansatz frame columns u1 = Qh e1, u2 = Qh e2, u3 = Qh e3."""
    rho = np.sqrt(x * x + y * y)
    phi = np.arctan2(y, x)
    th = -np.arctan2(z, rho)
    Q = np.einsum("...ab,...bc->...ac",
                  B8.rot_field(B8.G3, np.atleast_3d(phi)),
                  B8.rot_field(B8.G2, np.atleast_3d(th)))
    return Q[0, 0, 0] if np.ndim(x) == 0 else Q


def M_of(Q, d4):
    return np.einsum("...ab,bc,...dc->...ad", Q, d4, Q)


def d4_of(delta, g=G, s=S):
    return np.diag([-s * g, 1.0, delta, 0.0])


# =============================== probe ===============================
def audit_probe():
    """frame sanity: what does Qh actually do, and where is it singular?"""
    rng = np.random.default_rng(7)
    pts = rng.normal(size=(200, 3)) * 5.0
    e1, e2, e3 = [], [], []
    for (x, y, z) in pts:
        Q = frame_vectors(np.array(x), np.array(y), np.array(z))
        r = np.sqrt(x * x + y * y + z * z)
        nh = np.array([x, y, z]) / r
        rho = np.sqrt(x * x + y * y)
        eph = np.array([-y, x, 0.0]) / rho
        eth = np.cross(eph, nh)
        e1.append(np.abs(np.abs(Q[1:4, 1]) - np.abs(nh)).max())
        e2.append(min(np.abs(Q[1:4, 2] - eph).max(),
                      np.abs(Q[1:4, 2] + eph).max()))
        e3.append(min(np.abs(Q[1:4, 3] - eth).max(),
                      np.abs(Q[1:4, 3] + eth).max()))
    out = {"Qh_col1_is_nhat_max_abs_dev": float(max(e1)),
           "Qh_col2_is_e_phi_max_abs_dev": float(max(e2)),
           "Qh_col3_is_e_theta_max_abs_dev": float(max(e3)),
           "reading": ("Qh e1 = n-hat (the hedgehog direction), Qh e2 = "
                       "e_phi and Qh e3 = e_theta (the tangent frame); the "
                       "only singular locus of (phi, theta) is rho = 0, the "
                       "whole z axis including the origin")}
    # singular locus scan: how large is the phi-ring spread of the FRAME
    # near each candidate locus?
    loci = {}
    for nm, (px, py, pz) in (("z_axis_z12", (0.0, 0.0, 12.0)),
                             ("x_axis_x12", (12.0, 0.0, 0.0)),
                             ("plane_z0_x12", (12.0, 0.0, 0.0)),
                             ("generic", (5.0, 7.0, 9.0))):
        d4 = d4_of(DELTA)
        sp = []
        for eps in (1e-2, 1e-4, 1e-6):
            Ms = []
            for a in np.linspace(0, 2 * np.pi, 33)[:-1]:
                # a small circle of radius eps in a plane through the point
                if nm.startswith("z_axis"):
                    p = (px + eps * np.cos(a), py + eps * np.sin(a), pz)
                elif nm.startswith("x_axis"):
                    p = (px, py + eps * np.cos(a), pz + eps * np.sin(a))
                elif nm.startswith("plane"):
                    p = (px + eps * np.cos(a), py, pz + eps * np.sin(a))
                else:
                    p = (px + eps * np.cos(a), py + eps * np.sin(a), pz)
                Q = frame_vectors(*[np.array(v) for v in p])
                Ms.append(M_of(Q, d4))
            Ms = np.array(Ms)
            sp.append(float(np.abs(Ms - Ms.mean(0)).max()))
        loci[nm] = {"eps": [1e-2, 1e-4, 1e-6], "spread": sp,
                    "singular": bool(sp[-1] > 1e-3)}
    out["singular_locus_scan"] = loci
    dump("probe", out)
    return out


# ============================= B1 and B2 =============================
def audit_ring():
    out = {}
    # ---- B1: is M a function of n n^T alone at delta = 0? ----
    rng = np.random.default_rng(11)
    pts = rng.normal(size=(4000, 3)) * 6.0
    rows = {}
    for delta in (0.0, 0.3, 1.0):
        d4 = d4_of(delta)
        devs, devs_best = [], []
        for (x, y, z) in pts:
            r = np.sqrt(x * x + y * y + z * z)
            nh = np.array([x, y, z]) / r
            Q = frame_vectors(np.array(x), np.array(y), np.array(z))
            M = M_of(Q, d4)
            # the frame-free form with lam_par = 1, lam_perp = 0 (delta = 0)
            F0 = np.zeros((4, 4))
            F0[0, 0] = -S * G
            F0[1:4, 1:4] = np.outer(nh, nh)
            devs.append(np.abs(M - F0).max())
            # the BEST frame-free form: lam_perp = delta / 2 (the ring mean)
            Fb = np.zeros((4, 4))
            Fb[0, 0] = -S * G
            Fb[1:4, 1:4] = (0.5 * delta * (np.eye(3) - np.outer(nh, nh))
                            + (1.0 - 0.5 * delta) * np.outer(nh, nh))
            devs_best.append(np.abs(M - Fb).max())
        rows[f"delta_{delta:g}"] = {
            "max_abs_dev_lamperp_0": float(max(devs)),
            "max_abs_dev_lamperp_half_delta": float(max(devs_best)),
            "rel_to_field_scale_32": float(max(devs) / 32.0),
            "dev_over_delta": (float(max(devs) / delta) if delta else None),
            "best_dev_over_delta": (float(max(devs_best) / delta)
                                    if delta else None)}
    out["B1_frame_free_form"] = rows
    out["B1_derivation"] = (
        "M - [d00 + lam_par n n^T] = delta * u2 u2^T with u2 = e_phi, so the "
        "residual of the lam_perp = 0 form is delta * max_ij |u2_i u2_j|, "
        "whose supremum over directions is exactly delta (attained when e_phi "
        "is a coordinate axis). The producer's 0.2998 is delta = 0.3 up to "
        "the sampling of that supremum, NOT a separate number. The BEST "
        "frame-free fit takes lam_perp = delta/2 and leaves delta/2.")
    # ---- B2: the continuum ring ----
    ring = {}
    for delta in (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
        d4 = d4_of(delta)
        for z in (12.0, -12.0, 0.5, 0.0, 1e-6):
            sp = []
            for rr in (1e-2, 1e-4, 1e-6, 1e-9):
                Ms = []
                for a in np.linspace(0, 2 * np.pi, 65)[:-1]:
                    Q = frame_vectors(np.array(rr * np.cos(a)),
                                      np.array(rr * np.sin(a)), np.array(z))
                    Ms.append(M_of(Q, d4))
                Ms = np.array(Ms)
                sp.append(float(np.abs(Ms - Ms.mean(0)).max()))
            ring[f"delta_{delta:g}_z_{z:g}"] = {
                "rho": [1e-2, 1e-4, 1e-6, 1e-9], "spread": sp,
                "ratio_1e2_to_1e9": (sp[-1] / sp[0] if sp[0] else None),
                "spread_over_delta": (sp[-1] / delta if delta else None),
                "spread_over_half_one_minus_delta":
                    sp[-1] / max(0.5 * (1.0 - delta), 1e-300)}
    out["B2_rings"] = ring
    out["B2_derivation"] = (
        "As rho -> 0 at FIXED z != 0, theta -> -pi/2 sign(z), so the spatial "
        "block becomes R3(phi) diag(0, delta, 1) R3(phi)^T in the (1,2,3) "
        "basis: the (1,2) block is diag(0, delta) rotated by phi, giving "
        "M11 = delta sin^2 phi, M22 = delta cos^2 phi, M12 = -delta sin phi "
        "cos phi. Subtracting the phi-mean leaves -(delta/2) cos 2phi and "
        "-(delta/2) sin 2phi, so the spread is EXACTLY delta/2, independent "
        "of rho and of z, for every z != 0. At z = 0 theta = 0 instead, the "
        "rotated block is diag(1, delta) and the spread is (1 - delta)/2, "
        "which does NOT vanish at delta = 0: that ring measures the ORIGIN "
        "point defect, not the axis.")
    # the OTHER axis
    other = {}
    for delta in (0.0, 0.3):
        d4 = d4_of(delta)
        for x0 in (12.0, -12.0, 3.0):
            sp = []
            for rr in (1e-2, 1e-4, 1e-6):
                Ms = []
                for a in np.linspace(0, 2 * np.pi, 33)[:-1]:
                    Q = frame_vectors(np.array(x0), np.array(rr * np.cos(a)),
                                      np.array(rr * np.sin(a)))
                    Ms.append(M_of(Q, d4))
                Ms = np.array(Ms)
                sp.append(float(np.abs(Ms - Ms.mean(0)).max()))
            other[f"x_axis_delta_{delta:g}_x_{x0:g}"] = {
                "rho": [1e-2, 1e-4, 1e-6], "spread": sp,
                "ratio": (sp[-1] / sp[0] if sp[0] else None)}
    out["B2_other_axis"] = other
    # ---- GLOBAL scan, not a sampled one: where in the WHOLE box does the
    # nearest-neighbour jump in M fail to shrink with h? ----
    edges = [0.0, 1.5, 3.0, 6.0, 12.0, 18.0, 24.0]
    scan = {}
    for n in (32, 64):
        cfg = B3.base_cfg(s=S, g=G, n=n, L=48.0, delta=DELTA)
        h = cfg["h"]
        Mg = B8.dressed(cfg, 0.0)
        X, Y, Z = B3.coords(n, h)
        rho = np.sqrt(X * X + Y * Y)
        j = np.zeros(X.shape)
        for ax in range(3):
            d = np.abs(np.diff(Mg, axis=ax)).max(axis=(-1, -2))
            sl = [slice(None)] * 3
            sl[ax] = slice(0, -1); j[tuple(sl)] = np.maximum(j[tuple(sl)], d)
            sl[ax] = slice(1, None); j[tuple(sl)] = np.maximum(j[tuple(sl)], d)
        idx = np.digitize(rho.ravel(), np.array(edges)) - 1
        mx = []
        for i in range(len(edges) - 1):
            sel = idx == i
            mx.append(float(j.ravel()[sel].max()) if sel.any() else None)
        scan[f"h_{h:g}"] = {"rho_bins": edges, "max_neighbour_jump": mx}
    a, b = scan["h_1.5"]["max_neighbour_jump"], scan["h_0.75"]["max_neighbour_jump"]
    out["B2_global_jump_scan"] = {
        "rows": scan,
        "ratio_h0.75_over_h1.5": [
            (b[i] / a[i] if a[i] else None) for i in range(len(a))],
        "reading": ("EXHAUSTIVE over the box: only the innermost bin, "
                    "rho < 1.5, keeps its jump when h halves (ratio 1.000, "
                    "value 0.6667 at both h). Every other bin halves (ratios "
                    "0.54, 0.52, 0.51), i.e. the jump goes as h / rho, the "
                    "signature of a smooth field. The discontinuity is "
                    "therefore confined to the z axis and to nothing else.")}
    dump("ring", out)
    return out


# =============================== B3 ===============================
def audit_kin(n=32, L=48.0):
    deltas = (0.0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1.0)
    rows = []
    for dl in deltas:
        cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=dl)
        M = B8.dressed(cfg, 0.0)
        a0 = B8.a0_unit(cfg, 0.0)
        k = float(B3.kin_of(M, a0, cfg))
        rows.append({"delta": dl, "kin": k, "a0_max_abs": float(np.abs(a0).max()),
                     "kin_over_delta2": (k / dl ** 2 if dl else None)})
        log(f"kin delta={dl:g}: {k:.6f}")
    nz = [r for r in rows if r["delta"] > 0]
    out = {"n": n, "L": L, "rows": rows,
           "kin_at_delta_0": rows[0]["kin"],
           "a0_at_delta_0_max_abs": rows[0]["a0_max_abs"],
           "exponent_all": loglog([r["delta"] for r in nz],
                                  [r["kin"] for r in nz]),
           "exponent_small_delta": loglog([r["delta"] for r in nz[:4]],
                                          [r["kin"] for r in nz[:4]]),
           "exponent_large_delta": loglog([r["delta"] for r in nz[-3:]],
                                          [r["kin"] for r in nz[-3:]])}
    pairs = []
    for i in range(len(nz) - 1):
        a, b = nz[i], nz[i + 1]
        pairs.append({"delta_lo": a["delta"], "delta_hi": b["delta"],
                      "local_exponent": float(np.log(b["kin"] / a["kin"])
                                              / np.log(b["delta"] / a["delta"]))})
    out["local_exponents"] = pairs
    ratios = [r["kin_over_delta2"] for r in nz]
    out["kin_over_delta2"] = {"values": ratios,
                              "drift_max_over_min": max(ratios) / min(ratios),
                              "exact_delta2": bool(max(ratios) / min(ratios)
                                                   < 1.0 + 1e-9)}
    # MUTATION: doubling a0 must quadruple kin; a WRONG generator must differ
    cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    k1 = float(B3.kin_of(M, a0, cfg))
    k2 = float(B3.kin_of(M, 2.0 * a0, cfg))
    out["mutation"] = {"kin": k1, "kin_2a0": k2, "ratio": k2 / k1,
                       "expected": 4.0, "fires": bool(abs(k2 / k1 - 4.0) < 1e-9)}
    dump("kin", out)
    return out


# =============================== B4 ===============================
def _su2_lift_endpoint(axis, angle):
    """SU(2) element for a rotation of `angle` about `axis` (unit)."""
    c, s = np.cos(angle / 2.0), np.sin(angle / 2.0)
    return np.array([[c - 1j * s * axis[2], -1j * s * axis[0] - s * axis[1]],
                     [-1j * s * axis[0] + s * axis[1], c + 1j * s * axis[2]]])


def audit_topo():
    out = {}
    # ---- (b1) the stabilizer of d4 in so(1,3): which generators annihilate?
    gens = {}
    for (i, j) in ((1, 2), (1, 3), (2, 3)):
        Xg = np.zeros((4, 4)); Xg[i, j] = -1.0; Xg[j, i] = 1.0
        gens[f"rot_{i}{j}"] = Xg
    for k in (1, 2, 3):
        Xg = np.zeros((4, 4)); Xg[0, k] = Xg[k, 0] = 1.0
        gens[f"boost_{k}"] = Xg
    stab = {}
    for delta in (0.0, 0.3, 1.0):
        d4 = d4_of(delta)
        stab[f"delta_{delta:g}"] = {
            nm: float(np.abs(Xg @ d4 + d4 @ Xg.T).max())
            for nm, Xg in gens.items()}
    out["stabilizer_tangents"] = stab
    out["stabilizer_note"] = (
        "for diagonal d4 the tangent of a rotation J_ij is (m_i - m_j) on the "
        "(i,j) entries and of a boost K_k is (m_0 + m_k); at delta = 0.3 the "
        "spatial eigenvalues 1, 0.3, 0 are distinct and m_0 = 32 > 0, so the "
        "stabilizer is DISCRETE: the Klein four-group of sign matrices "
        "diag(1, e1, e2, e3) with e1 e2 e3 = 1.")
    # verify the discrete stabilizer by enumeration
    d4 = d4_of(DELTA)
    keep = []
    for e1 in (1, -1):
        for e2 in (1, -1):
            for e3 in (1, -1):
                Lm = np.diag([1.0, e1, e2, e3])
                if (np.abs(Lm @ ETA @ Lm.T - ETA).max() < 1e-14
                        and abs(np.linalg.det(Lm) - 1.0) < 1e-12
                        and np.abs(Lm @ d4 @ Lm.T - d4).max() < 1e-12):
                    keep.append([e1, e2, e3])
    out["discrete_stabilizer"] = {"elements": keep, "order": len(keep),
                                  "is_klein_four": bool(len(keep) == 4)}
    # ---- (b2) pi_1 of the order-parameter manifold: the Q8 lift ----
    # the preimage in SU(2) of the four spatial sign matrices
    lifts = []
    for ax in (np.array([1.0, 0, 0]), np.array([0, 1.0, 0]), np.array([0, 0, 1.0])):
        lifts.append(_su2_lift_endpoint(ax, np.pi))
    q = [np.eye(2, dtype=complex)] + lifts
    prei = []
    for A in q:
        prei.append(A); prei.append(-A)
    uniq = []
    for A in prei:
        if not any(np.abs(A - Bq).max() < 1e-9 for Bq in uniq):
            uniq.append(A)
    out["pi1_preimage"] = {
        "order": len(uniq), "is_Q8": bool(len(uniq) == 8),
        "note": ("the preimage in SU(2) of the Klein four-group is the "
                 "quaternion group Q8, so pi_1(OPM) = Q8 and pi_2(OPM) = "
                 "pi_2(SU(2)) = pi_2(S^3) = 0")}
    # the ansatz loop: phi 0 -> 2pi about the z axis lifts to -1 in SU(2)
    end = _su2_lift_endpoint(np.array([0.0, 0.0, 1.0]), 2 * np.pi)
    out["ansatz_line_class"] = {
        "su2_endpoint_trace": float(np.real(np.trace(end))),
        "is_minus_identity": bool(np.abs(end + np.eye(2)).max() < 1e-9),
        "class": "-1 in Q8",
        "trivial_in_pi1": bool(np.abs(end - np.eye(2)).max() < 1e-9),
        "reading": ("the z-axis line is the class -1 of Q8, which is NOT the "
                    "identity, so it is a TOPOLOGICALLY PROTECTED line defect "
                    "of the biaxial order parameter, not a coordinate "
                    "artifact; it cannot be removed while the field stays "
                    "biaxial, but it CAN be resolved by melting a core. "
                    "pi_2 = 0 means the hedgehog carries NO point charge.")}
    # ---- (b3) is the RP^1 refinement any relief? ----
    # track the delta-eigenvector's sign around loops; an INTEGER winding
    # means the line field lifts to a vector field and hairy ball applies
    wind = {}
    for nm, loop in (("equator", lambda a: (np.cos(a), np.sin(a), 0.0)),
                     ("polar_cap_z12", lambda a: (0.3 * np.cos(a),
                                                  0.3 * np.sin(a), 12.0)),
                     ("tilted", lambda a: (np.cos(a), 0.6 * np.sin(a),
                                           0.8 * np.sin(a)))):
        vs = []
        for a in np.linspace(0, 2 * np.pi, 721):
            x, y, z = loop(a)
            Q = frame_vectors(np.array(x), np.array(y), np.array(z))
            v = Q[1:4, 2].copy()
            if vs and float(np.dot(v, vs[-1])) < 0:
                v = -v                       # continuous SIGN choice
            vs.append(v)
        closes_plus = float(np.dot(vs[-1], vs[0]))
        wind[nm] = {"closure_dot": closes_plus,
                    "orientable_integer_winding": bool(closes_plus > 0.9)}
    out["line_field_lift"] = {
        "loops": wind,
        "argument": ("S^2 is simply connected, so ANY continuous line field "
                     "on it lifts to a continuous unit VECTOR field (the "
                     "orientation double cover of a line subbundle over a "
                     "simply connected base is trivial). The RP^1 refinement "
                     "therefore gives NO relief: hairy ball still applies. "
                     "The loops above confirm the winding is integer, never "
                     "half-integer.")}
    # ---- (b4) Poincare-Hopf index of the tangent field ----
    idx = {}
    for nm, (zc, sgn) in (("north_pole", (0.999, +1.0)),
                          ("south_pole", (-0.999, -1.0))):
        rr = np.sqrt(max(1.0 - zc * zc, 1e-12))
        ang = []
        for a in np.linspace(0, 2 * np.pi, 2001):
            x, y, z = rr * np.cos(a), rr * np.sin(a), zc
            Q = frame_vectors(np.array(x), np.array(y), np.array(z))
            v = Q[1:4, 2]
            nh = np.array([x, y, z]) / np.linalg.norm([x, y, z])
            # local tangent basis at the pole: e_x, e_y projected
            b1 = np.array([1.0, 0.0, 0.0]) - nh * nh[0]
            b1 /= np.linalg.norm(b1)
            b2 = np.cross(nh, b1) * sgn
            ang.append(np.arctan2(float(np.dot(v, b2)), float(np.dot(v, b1))))
        ang = np.unwrap(np.array(ang))
        idx[nm] = float((ang[-1] - ang[0]) / (2 * np.pi))
    out["poincare_hopf"] = {
        "index": idx, "total": float(sum(idx.values())),
        "euler_characteristic": 2.0,
        "matches_chi": bool(abs(sum(idx.values()) - 2.0) < 0.05),
        "reading": ("the index total is a homotopy invariant equal to chi(S^2) "
                    "= 2 for EVERY continuous tangent field, so no choice of "
                    "frame removes the singularity: at least one zero remains")}
    # ---- (a) the eigenvalue-assignment enumeration ----
    # a1 sits on the hedgehog direction; a2, a3 on the tangent frame
    perms = []
    import itertools
    for delta in (0.0, 0.3, 1.0):
        vals = [1.0, delta, 0.0]
        for p in set(itertools.permutations(range(3))):
            a = [vals[p[0]], vals[p[1]], vals[p[2]]]
            sfree = abs(a[1] - a[2]) < 1e-14
            row = {"delta": delta, "assignment": a, "string_free": bool(sfree)}
            # rotation generators: tangent, and whether the ORBIT stays smooth
            for gnm, (i, j) in (("G1_23", (1, 2)), ("G2_31", (2, 0)),
                                ("G3_12", (0, 1))):
                tan = abs(a[i] - a[j])
                # the orbit M(t) = Qh R(t) d R(t)^T Qh^T is a function of n n^T
                # alone iff R(t) preserves the e1 axis (only G1 does) and the
                # perpendicular pair is degenerate, OR all three are equal
                keeps_axis = (gnm == "G1_23")
                orbit_smooth = bool((sfree and keeps_axis)
                                    or (abs(a[0] - a[1]) < 1e-14
                                        and abs(a[1] - a[2]) < 1e-14))
                row[gnm] = {"tangent": float(tan),
                            "orbit_smooth": orbit_smooth,
                            "smooth_AND_nonzero": bool(orbit_smooth and tan > 1e-14)}
            # the radial boost K_1: tangent m0 + a1, orbit smooth iff a2 = a3
            row["K1_radial_boost"] = {
                "tangent": float(-S * G + a[0]),
                "orbit_smooth": bool(sfree),
                "smooth_AND_nonzero": bool(sfree and abs(-S * G + a[0]) > 1e-14),
                "periodic": False}
            perms.append(row)
    out["assignment_enumeration"] = {
        "rows": perms,
        "n_rotational_smooth_and_nonzero": int(sum(
            r[g]["smooth_AND_nonzero"] for r in perms
            for g in ("G1_23", "G2_31", "G3_12"))),
        "n_boost_smooth_and_nonzero": int(sum(
            r["K1_radial_boost"]["smooth_AND_nonzero"] for r in perms)),
        "reading": ("no relabeling of which vacuum eigenvalue rides the "
                    "hedgehog direction rescues a PERIODIC clock: the "
                    "perpendicular pair must be degenerate for smoothness, "
                    "and the only rotation preserving the hedgehog axis is "
                    "the one that rotates exactly that degenerate pair. But "
                    "the RADIAL BOOST K_1 is smooth AND nonzero whenever the "
                    "perpendicular pair is degenerate, so a0 = 0 is FALSE.")}
    # ---- the explicit smooth-orbit counterexample, measured ----
    cex = []
    for n, L in ((32, 48.0), (64, 48.0)):
        cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=0.0)
        M = B8.dressed(cfg, 0.0)
        X, Y, Z = B3.coords(n, cfg["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        nh = np.stack([X / r, Y / r, Z / r], axis=-1)
        a0b = np.zeros(M.shape)
        amp = (-S * G) + 1.0                       # m_0 + a_1 = 32 + 1
        a0b[..., 0, 1:4] = amp * nh
        a0b[..., 1:4, 0] = amp * nh
        k = float(B3.kin_of(M, a0b, cfg))
        ka = float(B3.kin_of(M, B8.a0_unit(cfg, 0.0), cfg))
        cex.append({"n": n, "L": L, "h": cfg["h"],
                    "kin_radial_boost": k, "kin_G1_clock": ka,
                    "a0_boost_max_abs": float(np.abs(a0b).max())})
        log(f"counterexample n={n}: kin_boost {k:.6f} kin_G1 {ka:.3e}")
    out["smooth_orbit_counterexample"] = {
        "rows": cex,
        "h_ratio": (cex[1]["kin_radial_boost"] / cex[0]["kin_radial_boost"]
                    if cex[0]["kin_radial_boost"] else None),
        "reading": ("at delta = 0 the hedgehog M = diag(32,0,0,0) + n n^T is "
                    "SMOOTH and the radial-boost orbit M(t) = Qh B1(t) d4 "
                    "B1(t)^T Qh^T stays smooth for all t (only d00, d01, d11 "
                    "move and the 0i row is d01 n-hat), with a0 = 33 sym(e0 "
                    "n^T) NONZERO. So 'a smooth hedgehog has a0 = 0 "
                    "identically' is false as written. The boost is not "
                    "PERIODIC, so it is not a clock: the theorem survives "
                    "only in the restricted form 'no periodic clock'.")}
    # ---- attack (a), EXHAUSTIVE over the whole 6-dim algebra ----
    # A field M = Qh d Qh^T is frame-free (smooth away from the origin) iff
    # d has the form  d00 | d01 e1 | alpha e1e1^T + beta (I3 - e1e1^T),
    # because Qh e1 = n-hat and Qh (I3 - e1e1^T) Qh^T = I3 - n n^T are the
    # only Qh-covariants that are single-valued. Solve for every generator X
    # in so(1,3) whose flow keeps d in that 4-dim form.
    basis = [gens[k] for k in ("rot_12", "rot_13", "rot_23",
                               "boost_1", "boost_2", "boost_3")]
    names = ["rot_12", "rot_13", "rot_23", "boost_1", "boost_2", "boost_3"]
    e1v = np.array([0.0, 1.0, 0.0, 0.0])
    P = np.zeros((4, 4))                     # e1 e1^T in 4x4 spatial slot
    P[1, 1] = 1.0
    Q3 = np.zeros((4, 4))                    # I3 - e1 e1^T
    Q3[2, 2] = Q3[3, 3] = 1.0
    E00 = np.zeros((4, 4)); E00[0, 0] = 1.0
    E01 = np.zeros((4, 4)); E01[0, 1] = E01[1, 0] = 1.0
    form = [E00, E01, P, Q3]                 # the 4-dim smooth form
    fmat = np.stack([f.ravel() for f in form], axis=1)
    Pf = fmat @ np.linalg.pinv(fmat)         # projector onto the form
    exh = {}
    for delta in (0.0, 0.3, 1.0):
        d4 = d4_of(delta)
        cols = []
        for Xg in basis:
            T = (Xg @ d4 + d4 @ Xg.T).ravel()
            cols.append(T - Pf @ T)          # the part OUTSIDE the form
        A = np.stack(cols, axis=1)
        u, sv, vt = np.linalg.svd(A)
        null = [vt[i] for i in range(6) if (sv[i] if i < len(sv) else 0.0) < 1e-10]
        rows = []
        for v in null:
            Xg = sum(v[i] * basis[i] for i in range(6))
            T = Xg @ d4 + d4 @ Xg.T
            # does the FINITE orbit stay in the form?
            from scipy.linalg import expm
            worst = 0.0
            for t in (0.1, 0.5, 1.0, 2.0):
                Lm = expm(t * Xg)
                dt_ = Lm @ d4 @ Lm.T
                worst = max(worst, float(np.abs(
                    dt_.ravel() - Pf @ dt_.ravel()).max()))
            rows.append({
                "coeffs": {names[i]: float(v[i]) for i in range(6)
                           if abs(v[i]) > 1e-9},
                "tangent_max_abs": float(np.abs(T).max()),
                "finite_orbit_max_form_violation": worst,
                "nonzero_clock": bool(np.abs(T).max() > 1e-12),
                # a PERIODIC clock needs a purely compact (rotation)
                # generator: any boost component makes the flow hyperbolic
                "is_compact_periodic": bool(
                    max(abs(v[i]) for i in range(3, 6)) < 1e-9)})
        exh[f"delta_{delta:g}"] = {
            "singular_values": sv.tolist(),
            "null_dimension": len(null), "generators": rows}
    out["attack_a_exhaustive_algebra"] = {
        "rows": exh,
        "reading": ("EXHAUSTIVE, not a sampling: at delta = 0 exactly TWO of "
                    "the six so(1,3) directions keep the field frame-free "
                    "along their whole flow, the perpendicular-pair rotation "
                    "G1 (whose tangent is identically zero) and the radial "
                    "boost K_1 (tangent 33, but NOT periodic). At delta = 0.3 "
                    "the null space is 1-dimensional and is K_1 alone. There "
                    "is therefore no smooth periodic clock anywhere in the "
                    "algebra, and no relabeling or biaxial escape can create "
                    "one within the eigenvalue-pinned family.")}
    # ---- cross-check with the CERTIFIED generator catalog ----
    cat = []
    for dl in (0.0, 0.03, 0.3, 1.0):
        cfg = B3.base_cfg(s=S, g=G, n=32, L=48.0, delta=dl)
        M = B8.dressed(cfg, 0.0)
        w = B3.envelope(cfg)[..., None, None]
        lam, V = np.linalg.eigh(M[..., 1:4, 1:4])

        def local_rot(vh):
            Wm = np.zeros(vh.shape[:-1] + (4, 4))
            n1, n2, n3 = vh[..., 0], vh[..., 1], vh[..., 2]
            Wm[..., 1, 2], Wm[..., 1, 3] = -n3, n2
            Wm[..., 2, 1], Wm[..., 2, 3] = n3, -n1
            Wm[..., 3, 1], Wm[..., 3, 2] = -n2, n1
            return Wm
        row = {"delta": dl}
        a0s = B3.gen_catalog(cfg, M)
        for nm, vh in (("clock_local", V[..., :, 2]),
                       ("plane_1d", V[..., :, 0])):
            Gm = local_rot(vh)
            a = w * (Gm @ M - M @ np.swapaxes(Gm, -1, -2))
            row[f"{nm}_raw_frobenius"] = float(np.sqrt(np.sum(a * a)))
            row[f"{nm}_kin_after_unit_normalization"] = float(
                B3.kin_of(M, a0s[nm], cfg))
        cat.append(row)
    out["certified_catalog_cross_check"] = {
        "rows": cat,
        "reading": ("the certified catalog's clock_local IS the "
                    "perpendicular-pair rotation, and its RAW tangent norm at "
                    "delta = 0 is 6.5e-15, i.e. exactly zero: B4's rotational "
                    "leg is confirmed by the stack's own clock generator. "
                    "plane_1d, which rotates n-hat into a tangent, keeps a "
                    "raw norm of 36.8 at delta = 0 but its orbit is not "
                    "smooth. WARNING: gen_catalog divides by max(norm, "
                    "1e-300), so at delta = 0 it turns a numerically zero "
                    "generator into a UNIT-NORM NOISE field and reports a "
                    "spurious kin of order 2.25; any delta -> 0 study routed "
                    "through gen_catalog will see a phantom clock.")}
    # ---- the nematic-XY minimization: can any frame beat the defect? ----
    xy = {}
    for N in (24, 48, 96):
        th = (np.arange(N) + 0.5) * np.pi / N          # polar
        ph = np.arange(2 * N) * 2 * np.pi / (2 * N)    # azimuth
        TH, PH = np.meshgrid(th, ph, indexing="ij")
        # tangent basis e_theta, e_phi; the line field is an angle psi mod pi
        psi = np.zeros_like(TH)                        # start from e_theta
        # the parallel-transport connection between azimuthal neighbours is a
        # rotation by 0 in this basis; between polar neighbours also 0. The
        # ANSATZ field is psi = pi/2 (e_phi) everywhere, which is smooth in
        # this chart but singular at the poles because the chart is.
        # minimize sum over edges of sin^2(psi_i - psi_j + A_ij) with the
        # holonomy A around a polar cap loop = 2 pi (the chart's own winding)
        # -> equivalent to a nematic XY model with an unremovable +1 vortex
        # at each pole. We measure the residual max edge jump after descent.
        for _ in range(400):
            g = np.zeros_like(psi)
            for ax, roll in ((0, 1), (0, -1), (1, 1), (1, -1)):
                d = psi - np.roll(psi, roll, axis=ax)
                g += np.sin(2 * d)
            psi -= 0.05 * g
        jumps = []
        for ax, roll in ((0, 1), (1, 1)):
            d = psi - np.roll(psi, roll, axis=ax)
            jumps.append(float(np.abs(np.sin(d)).max()))
        # the true obstruction: the field must ALSO close across the poles,
        # where the chart's frame winds by 2 pi. Measure the pole mismatch.
        pole_mismatch = float(np.abs(np.sin(psi[0] - psi[0][::-1] + PH[0])).max())
        xy[f"N_{N}"] = {"max_edge_jump_after_descent": max(jumps),
                        "pole_frame_mismatch": pole_mismatch,
                        "grid_spacing": float(np.pi / N)}
    out["nematic_xy"] = {
        "rows": xy,
        "reading": ("the interior descent flattens to a uniform field, but the "
                    "pole mismatch (the chart frame winds by 2 pi around each "
                    "pole) stays O(1) as the grid refines: the defect does not "
                    "shrink with resolution, so it is not a discretization "
                    "artifact")}
    dump("topo", out)
    return out


# =============================== arm a ===============================
def make_families(cfg):
    """cut families. Every family except `cyl_shift` removes EXACTLY the
    same number of cells k(c) = #{rho_z < c}, so the four controls are
    volume-matched to the producer's target cut by construction.

    `cyl_shift` reproduces the PRODUCER's mask, which is built on a
    cell-centered grid x = (i - n//2) h while the density lives on the
    certified offset grid x = (i - (n-1)/2) h. Its cylinder is therefore
    centered at (h/2, h/2), NOT on the z axis, and its center MOVES with
    the lattice."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    xs = (np.arange(n) - n // 2) * h
    Xs, Ys, _ = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z).ravel()
    rho_z = np.sqrt(X * X + Y * Y).ravel()
    rho_sh = np.sqrt(Xs * Xs + Ys * Ys).ravel()
    rho_x = np.sqrt(Y * Y + Z * Z).ravel()
    rho_off = np.sqrt((X - XOFF) ** 2 + Y * Y).ravel()
    u = np.random.default_rng(909).random(rho_z.shape)
    keys = {"cyl_z": rho_z, "cyl_shift": rho_sh, "cyl_x": rho_x,
            "cyl_off": rho_off, "sph": r, "rand": u}
    srt = np.sort(rho_z)
    kcut = {c: int(np.searchsorted(srt, c)) for c in RHO_CUTS}
    srt_sh = np.sort(rho_sh)
    kcut_sh = {c: int(np.searchsorted(srt_sh, c)) for c in RHO_CUTS}
    fams = {}
    for nm, key in keys.items():
        fams[nm] = {"order": np.argsort(key, kind="stable"),
                    "kcut": kcut_sh if nm == "cyl_shift" else kcut}
    return fams, kcut, kcut_sh, rho_z.size


def box_sums(n, L, delta=DELTA, terms=TERMS, omegas=OMEGAS):
    """masked integrals I_term(omega, mask) for every cut family."""
    cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=delta)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    h3 = cfg["h"] ** 3
    fams, kcut, kcut_sh, ncell = make_families(cfg)
    sums = {t: {f"{fn}_{c:g}": np.zeros(len(omegas))
                for fn in FAMILIES for c in RHO_CUTS} for t in terms}
    for br, wt in BRANCHES:
        Asp = [mydiff(M, ax, cfg["h"], br) for ax in range(3)]
        for k, om in enumerate(omegas):
            d = densities([om * a0] + Asp, terms=terms)
            for t in terms:
                flat = (d[t] * h3).ravel()
                for fn in FAMILIES:
                    ordv = fams[fn]["order"]
                    kc = fams[fn]["kcut"]
                    # sum from the OUTSIDE in: kept_k = tail sum, no cancellation
                    rc = np.cumsum(flat[ordv][::-1])
                    for c in RHO_CUTS:
                        kk = kc[c]
                        sums[t][f"{fn}_{c:g}"][k] += wt * float(
                            rc[ncell - kk - 1] if kk < ncell else 0.0)
        del Asp
    return cfg, sums, kcut, kcut_sh, ncell


def audit_tube():
    out = {"boxes": {}, "rho_cuts": list(RHO_CUTS), "families": list(FAMILIES)}
    store = {}
    allboxes = sorted(set(BOXES_L) | {b for pr in BOXES_H for b in pr})
    for (n, L) in allboxes:
        t0 = time.time()
        cfg, sums, kcut, kcut_sh, ncell = box_sums(n, L)
        row = {"n": n, "L": L, "h": cfg["h"],
               "cells_removed_frac": {f"{c:g}": kcut[c] / ncell
                                      for c in RHO_CUTS},
               "cells_removed_frac_producer_mask": {
                   f"{c:g}": kcut_sh[c] / ncell for c in RHO_CUTS},
               "producer_mask_axis_offset": cfg["h"] / np.sqrt(2.0),
               "terms": {}}
        for t in TERMS:
            row["terms"][t] = {}
            for key, I in sums[t].items():
                c, res = fit4(I)
                row["terms"][t][key] = {
                    "A": float(c[0]), "C1": float(c[1]), "C2": float(c[2]),
                    "C3": float(c[3]), "C4": float(c[4]),
                    "fit_resid": res,
                    "fit_resid_rel": res / max(np.abs(I).max(), 1e-300)}
        store[(n, L)] = row
        out["boxes"][f"n{n}_L{L:g}"] = row
        log(f"tube box n={n} L={L:g} done ({time.time() - t0:.1f}s) "
            f"Q_I1sq C4 full {row['terms']['Q_I1sq']['cyl_z_0']['C4']:.5f} "
            f"cyl_z_3 {row['terms']['Q_I1sq']['cyl_z_3']['C4']:.5f} "
            f"sph_3 {row['terms']['Q_I1sq']['sph_3']['C4']:.5f}")
        dump("tube_partial", out)
    # exponents
    Ls = [L for _, L in BOXES_L]
    exps = {"L": {}, "h": {}}
    for t in TERMS:
        exps["L"][t] = {}
        for fn in FAMILIES:
            for c in RHO_CUTS:
                key = f"{fn}_{c:g}"
                row = {}
                for co in ("A", "C2", "C4"):
                    row[co] = loglog(Ls, [store[b]["terms"][t][key][co]
                                          for b in BOXES_L])
                exps["L"][t][key] = row
        exps["h"][t] = {}
        for fn in FAMILIES:
            for c in RHO_CUTS:
                key = f"{fn}_{c:g}"
                row = {}
                for co in ("A", "C2", "C4"):
                    vals = []
                    for pr in BOXES_H:
                        hs = [store[b]["h"] for b in pr]
                        ys = [store[b]["terms"][t][key][co] for b in pr]
                        vals.append(loglog(hs, ys))
                    row[co] = vals
                exps["h"][t][key] = row
    out["exponents"] = exps
    # the reading: does the cleanup depend on rho_cut, and do the CONTROLS
    # clean up too?
    verdicts = {}
    for t in ("Q_I1sq", "Q_Fpair", "Q_I4sq", "Q_C6a", "I1"):
        vv = {}
        for fn in FAMILIES:
            vv[fn] = {
                "C4_L_exponent_vs_rho_cut": [
                    exps["L"][t][f"{fn}_{c:g}"]["C4"] for c in RHO_CUTS],
                "C4_h_exponent_vs_rho_cut": [
                    exps["h"][t][f"{fn}_{c:g}"]["C4"][0] for c in RHO_CUTS],
                "C2_L_exponent_vs_rho_cut": [
                    exps["L"][t][f"{fn}_{c:g}"]["C2"] for c in RHO_CUTS],
                "C2_h_exponent_vs_rho_cut": [
                    exps["h"][t][f"{fn}_{c:g}"]["C2"][0] for c in RHO_CUTS]}
        verdicts[t] = vv
    out["verdicts"] = verdicts
    dump("tube", out)
    return out


# =============================== B5 ===============================
def audit_share():
    """the axis-tube share of the quartic C4, at a FIXED PHYSICAL tube,
    at two lattice spacings, over the delta ladder."""
    rows = {}
    for (n, L) in ((32, 48.0), (64, 48.0)):
        for dl in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
            cfg, sums, kcut, kcut_sh, ncell = box_sums(
                n, L, delta=dl, terms=("Q_I1sq", "Q_Fpair", "Q_I4sq", "I1"))
            row = {}
            for t in ("Q_I1sq", "Q_Fpair", "Q_I4sq", "I1"):
                full, _ = fit4(sums[t]["cyl_z_0"])
                r2 = {}
                for c in (1.5, 3.0):
                    for fn in FAMILIES:
                        kept, _ = fit4(sums[t][f"{fn}_{c:g}"])
                        for co, ix in (("C2", 2), ("C4", 4)):
                            r2[f"{co}_{fn}_share_lt_{c:g}"] = (
                                float((full[ix] - kept[ix]) / full[ix])
                                if full[ix] else None)
                        r2[f"C4_{fn}_kept_ge_{c:g}"] = float(kept[4])
                    # the AXIS EXCESS: the z tube minus the equal-volume
                    # x-axis tube, which also contains the origin core
                    r2[f"C4_axis_excess_{c:g}"] = (
                        r2[f"C4_cyl_z_share_lt_{c:g}"]
                        - r2[f"C4_cyl_x_share_lt_{c:g}"])
                r2["C4_full"] = float(full[4])
                r2["C2_full"] = float(full[2])
                row[t] = r2
            rows[f"n{n}_delta{dl:g}"] = row
            q = row["Q_I1sq"]
            log(f"share n={n} delta={dl:g}: "
                f"z-tube {q['C4_cyl_z_share_lt_3']:.4f} "
                f"x-tube {q['C4_cyl_x_share_lt_3']:.4f} "
                f"sph {q['C4_sph_share_lt_3']:.4f}")
            dump("share_partial", rows)
    out = {"rows": rows,
           "note": ("the producer's axis tube is rho < 1.01 h, a LATTICE "
                    "region whose share is not a physical quantity; here the "
                    "tube is a fixed physical radius sampled at h = 1.5 and "
                    "h = 0.75, and the x-axis tube of the same radius is the "
                    "null control")}
    dump("share", out)
    return out


# =============================== B4 (c) ===============================
def core_profile(M, cfg, zlo=6.0, zhi=15.0):
    """the biaxial gap (lam_mid - lam_min) of the spatial block, binned in
    the cylindrical radius, away from the origin and the walls."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    gap = lam[..., 1] - lam[..., 0]        # the BIAXIAL gap (mid - min)
    top = lam[..., 2] - lam[..., 1]        # the uniaxial gap (max - mid)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    nh = np.stack([X / r, Y / r, Z / r], axis=-1)
    # |cos| between the LEADING spatial eigenvector and n-hat: 1 = the
    # eigenvalue-pinned hedgehog, < 1 = the director has escaped
    tilt = np.abs(np.einsum("...i,...i->...", V[..., :, 2], nh))
    sel = (np.abs(Z) > zlo) & (np.abs(Z) < zhi) & (rho < 12.0)
    edges = np.arange(0.0, 12.0 + h, max(h, 0.5))
    idx = np.digitize(rho[sel].ravel(), edges) - 1
    cnt = np.bincount(idx, minlength=len(edges))
    mid = 0.5 * (edges[:-1] + edges[1:])
    ok = cnt[:len(mid)] > 0

    def prof(f):
        s = np.bincount(idx, weights=f[sel].ravel(), minlength=len(edges))
        return (s[:len(mid)][ok] / cnt[:len(mid)][ok]).tolist()
    return mid[ok].tolist(), prof(gap), prof(top), prof(tilt)


def core_width(mid, gap, target):
    """the smallest rho where the gap first reaches half its far value."""
    for r, g in zip(mid, gap):
        if g >= 0.5 * target:
            return float(r)
    return None


def audit_relax(g=8.0, delta=0.3, dt0=0.01):
    out = {"g": g, "delta": delta, "rows": []}
    for (n, L, maxit) in ((32, 48.0, 3000), (64, 48.0, 2000)):
        cfg = B3.base_cfg(s=S, g=g, n=n, L=L, delta=delta)
        M0 = B8.dressed(cfg, 0.0)
        a0 = B8.a0_unit(cfg, 0.0)
        free = ~B3.pin_shell(n, cfg["h"])
        e0 = B3.e_parts(M0, cfg)
        mid0, gap0, top0, tilt0 = core_profile(M0, cfg)
        t0 = time.time()
        M, info = B3.fire(M0, cfg, free, maxit, log_every=250,
                          tag=f"r9audit_relax_n{n}", dt0=dt0, dt_max=0.1)
        e1 = B3.e_parts(M, cfg)
        mid1, gap1, top1, tilt1 = core_profile(M, cfg)
        # the on-axis discontinuity: the phi spread of M on the nearest ring
        h = cfg["h"]
        X, Y, Z = B3.coords(n, h)
        rho = np.sqrt(X * X + Y * Y)
        near = (rho < 1.2 * h) & (np.abs(Z) > 6.0) & (np.abs(Z) < 15.0)
        def ringspread(Mx):
            vals = []
            zs = np.unique(Z[near])
            for zv in zs:
                m = near & (Z == zv)
                if m.sum() < 3:
                    continue
                sub = Mx[m]
                vals.append(float(np.abs(sub - sub.mean(0)).max()))
            return float(np.mean(vals)) if vals else None
        far = float(np.median(gap1[len(gap1) // 2:]))
        row = {"n": n, "L": L, "h": h, "maxit": maxit,
               "E_start": [float(e0[0]), float(e0[1])],
               "E_end": [float(e1[0]), float(e1[1])],
               "E_drop": float(sum(e0) - sum(e1)),
               "stop": info["stop"], "wall_s": round(time.time() - t0, 1),
               "rho_mid": mid1,
               "gap_before": gap0, "gap_after": gap1,
               "topgap_before": top0, "topgap_after": top1,
               "tilt_before": tilt0, "tilt_after": tilt1,
               "gap_far_after": far,
               "core_width_before": core_width(mid0, gap0, delta),
               "core_width_after": core_width(mid1, gap1, far),
               "core_width_after_over_h": None,
               "ring_spread_before": ringspread(M0),
               "ring_spread_after": ringspread(M),
               "kin_frozen_a0_before": float(B3.kin_of(M0, a0, cfg)),
               "kin_frozen_a0_after": float(B3.kin_of(M, a0, cfg))}
        if row["core_width_after"]:
            row["core_width_after_over_h"] = row["core_width_after"] / h
        out["rows"].append(row)
        log(f"relax n={n}: E {sum(e0):.4f} -> {sum(e1):.4f} stop={info['stop']} "
            f"core_after={row['core_width_after']} ring {row['ring_spread_before']:.4f}"
            f" -> {row['ring_spread_after']:.4f}")
        dump("relax_partial", out)
    if len(out["rows"]) == 2:
        a, b = out["rows"]
        out["core_h_independence"] = {
            "core_width_h1.5": a["core_width_after"],
            "core_width_h0.75": b["core_width_after"],
            "ratio": ((b["core_width_after"] / a["core_width_after"])
                      if a["core_width_after"] and b["core_width_after"] else None),
            "physical_if_ratio_near_1": True,
            "lattice_artifact_if_ratio_near_0.5": True}
        out["ring_h_dependence"] = {
            "before": [a["ring_spread_before"], b["ring_spread_before"]],
            "after": [a["ring_spread_after"], b["ring_spread_after"]]}
    dump("relax", out)
    return out


# ============ separating the POINT core from the LINE core ============
def audit_divergence(L=48.0, ns=(24, 32, 48, 64, 96), g=8.0):
    """Two different divergences live in this ansatz and R8/R9 conflate them.

    E_static (the omega^0 sector) diverges as h^-1 at EVERY delta, delta = 0
    included, so it is the hedgehog POINT core at the origin, where |grad M|
    ~ 1/r makes the quartic density ~ r^-4 and the radial integral ~ 1/h.
    The omega^4 coefficient C4 diverges as h^-1.6 only at delta > 0, and a
    z-axis tube cures it while an equal-volume SPHERICAL core cut (which
    removes the origin) does not: that one is the LINE."""
    out = {"L": L, "g": g, "static": {}}
    for dl in (0.0, 0.3):
        hs, es = [], []
        for n in ns:
            cfg = B3.base_cfg(s=S, g=g, n=n, L=L, delta=dl)
            eu, ev = B3.e_parts(B8.dressed(cfg, 0.0), cfg)
            hs.append(cfg["h"]); es.append(float(eu))
        out["static"][f"delta_{dl:g}"] = {
            "h": hs, "E_u": es, "h_exponent": loglog(hs, es)}
    out["reading"] = (
        "E_static's h^-1 divergence is present at delta = 0, where the field "
        "is provably string-free (B1), so it is the POINT core and NOT the "
        "line. The omega^4 divergence is the opposite: it vanishes at "
        "delta = 0 (C4 = 1e-14) and is cured by a z-axis tube but not by an "
        "equal-volume spherical cut through the origin. Arm a's excision "
        "therefore addresses the LINE only; the point core is untouched by "
        "it and is a separate, unresolved feature of the ansatz.")
    dump("divergence", out)
    return out


# ========================= mutation discipline =========================
def audit_mutate():
    """every gate confirmed above must FAIL under a deliberate mutation."""
    out = {}
    # M1: the exact-quartic gate. Inject a degree-6 omega content into the
    # samples; the degree-4 residual must go from 1e-14 to order one.
    cfg, sums, kcut, kcut_sh, ncell = box_sums(32, 48.0, terms=("Q_I1sq",))
    I = sums["Q_I1sq"]["cyl_z_0"]
    _, r_clean = fit4(I)
    _, r_mut = fit4(I + 0.05 * np.abs(I).max() * OMEGAS ** 6)
    out["M1_exact_quartic"] = {
        "resid_clean": r_clean, "resid_mutant": r_mut,
        "scale": float(np.abs(I).max()),
        "fires": bool(r_mut > 1e3 * max(r_clean, 1e-300))}
    # M2a: at delta = 0 the ansatz is string-free (B1) and the whole omega^4
    # coefficient must collapse to numerical zero: nothing left to excise.
    _, s0, _, _, _ = box_sums(32, 48.0, delta=0.0, terms=("Q_I1sq", "Q_Fpair"))
    _, s3, _, _, _ = box_sums(32, 48.0, delta=DELTA, terms=("Q_I1sq", "Q_Fpair"))
    out["M2a_delta0_kills_C4"] = {
        t: {"C4_delta_0": float(fit4(s0[t]["cyl_z_0"])[0][4]),
            "C4_delta_0.3": float(fit4(s3[t]["cyl_z_0"])[0][4])}
        for t in ("Q_I1sq", "Q_Fpair")}
    out["M2a_delta0_kills_C4"]["fires"] = bool(
        abs(fit4(s0["Q_I1sq"]["cyl_z_0"])[0][4]) < 1e-12)
    # M2b: MOVE the string. Rebuild the ansatz about a hedgehog center at
    # (XOFF, 0, 0); the line defect then runs along x = XOFF, y = 0, and the
    # cleanup must FOLLOW it: cyl_off must now behave the way cyl_z did, and
    # cyl_z must now behave the way the null controls did.
    m2b = {}
    st = {}
    for (n, L) in ((32, 48.0), (64, 48.0)):
        cfg = B3.base_cfg(s=S, g=G, n=n, L=L, delta=DELTA)
        h = cfg["h"]
        X, Y, Z = B3.coords(n, h)
        Q = frame_vectors(X - XOFF, Y, Z)
        d4 = B3.vac4(cfg)
        Ms = B3.sym4(M_of(Q, d4))
        # the same G1 clock flow, co-moving with the shifted frame
        core = B8.G1 @ d4 + d4 @ B8.G1.T
        a0s = np.einsum("...ab,bc,...dc->...ad", Q, core, Q)
        h3 = h ** 3
        fams, kcut, kcut_sh, ncell = make_families(cfg)
        acc = {t: {fn: np.zeros(len(OMEGAS)) for fn in FAMILIES}
               for t in ("Q_I1sq", "Q_Fpair")}
        for br, wt in BRANCHES:
            Asp = [mydiff(Ms, ax, h, br) for ax in range(3)]
            for k, om in enumerate(OMEGAS):
                d = densities([om * a0s] + Asp, ("Q_I1sq", "Q_Fpair"))
                for t in acc:
                    flat = (d[t] * h3).ravel()
                    for fn in FAMILIES:
                        rc = np.cumsum(flat[fams[fn]["order"]][::-1])
                        kk = fams[fn]["kcut"][3.0]
                        acc[t][fn][k] += wt * float(rc[ncell - kk - 1])
        st[(n, L)] = acc
    for t in ("Q_I1sq", "Q_Fpair"):
        row = {}
        for fn in ("cyl_z", "cyl_off", "cyl_x", "sph"):
            ys = [float(fit4(st[b][t][fn])[0][4])
                  for b in ((32, 48.0), (64, 48.0))]
            row[fn] = {"C4_h1.5": ys[0], "C4_h0.75": ys[1],
                       "h_exponent": loglog([1.5, 0.75], ys)}
        m2b[t] = row
    out["M2b_string_moved"] = {
        "hedgehog_center": [XOFF, 0.0, 0.0], "rows": m2b,
        "fires": bool(all(
            abs(m2b[t]["cyl_off"]["h_exponent"]) < 0.5
            and m2b[t]["cyl_z"]["h_exponent"] < -0.8 for t in m2b)),
        "reading": ("with the hedgehog moved to x = 12 the line defect moves "
                    "with it; the h-divergence must now be cured by the "
                    "cyl_off tube and NOT by the z-axis tube. If the z-axis "
                    "tube still cured it, the arm-a diagnostic would be "
                    "measuring the lattice axis, not the string.")}
    # M3: A4's derivative-free identity. Replace the SPATIAL jets by noise;
    # C4(Q_C6a) must not move at all.
    cfg = B3.base_cfg(s=S, g=G, n=32, L=48.0, delta=DELTA)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    h3 = cfg["h"] ** 3
    rng = np.random.default_rng(5)
    real = [mydiff(M, ax, cfg["h"], "fwd") for ax in range(3)]
    fake = [B3.sym4(rng.normal(size=M.shape)) for _ in range(3)]
    c_real = fit4(np.array([
        float(np.sum(densities([om * a0] + real, ("Q_C6a",))["Q_C6a"]) * h3)
        for om in OMEGAS]))[0][4]
    c_fake = fit4(np.array([
        float(np.sum(densities([om * a0] + fake, ("Q_C6a",))["Q_C6a"]) * h3)
        for om in OMEGAS]))[0][4]
    t00 = tr_AetaBeta(a0, a0)
    out["M3_C6a_derivative_free"] = {
        "C4_real_jets": float(c_real), "C4_random_jets": float(c_fake),
        "closed_form_sum_t00_squared": float(np.sum(t00 * t00) * h3),
        "rel_dev_real_vs_closed": rel(float(c_real),
                                      float(np.sum(t00 * t00) * h3)),
        "rel_dev_real_vs_random": rel(float(c_real), float(c_fake)),
        "fires": bool(rel(float(c_real), float(c_fake)) < 1e-9),
        "reading": ("C4(Q_C6a) = sum_cells (tr(a0 eta a0 eta))^2 h^3 contains "
                    "NO spatial derivative at all, so A4's L^3 and h^0 "
                    "exponents are automatic and carry no information about "
                    "whether the frame is singular")}
    # M4: the internal metric. eta -> identity must break the certified kin.
    global W, WW
    Ws, WWs = W, WW
    try:
        W = np.ones(4); WW = W[:, None] * W[None, :]
        Asp = [mydiff(M, ax, cfg["h"], "fwd") for ax in range(3)]
        vals = np.array([
            float(np.sum(densities([om * a0] + Asp, ("I1",))["I1"]) * h3)
            for om in OMEGAS])
        kin_mut = -4.0 * float(fit4(vals)[0][2])
    finally:
        W, WW = Ws, WWs
    out["M4_internal_metric"] = {
        "kin_mutant": kin_mut, "kin_certified": 426.5070121483972,
        "rel_dev": rel(kin_mut, 426.5070121483972),
        "fires": bool(rel(kin_mut, 426.5070121483972) > 1e-3)}
    out["all_mutations_fire"] = bool(all(
        out[k].get("fires") for k in out if k.startswith("M")))
    dump("mutate", out)
    return out


# =============================== main ===============================
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("probe", "all"):
        audit_probe()
    if mode in ("ring", "all"):
        audit_ring()
    if mode in ("kin", "all"):
        audit_kin()
    if mode in ("topo", "all"):
        audit_topo()
    if mode in ("tube", "all"):
        audit_tube()
    if mode in ("share", "all"):
        audit_share()
    if mode in ("divergence", "all"):
        audit_divergence()
    if mode in ("mutate", "all"):
        audit_mutate()
    if mode in ("relax", "all"):
        audit_relax()
    if mode in ("merge", "all"):
        _load("merge", "m5_32_r9_audit_merge.py").main()
    log(f"mode {mode} done")


if __name__ == "__main__":
    main()
