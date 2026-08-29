"""M5.32 R11 ADVERSARIAL AUDIT: same-sign boost hedgehogs, the sign flip,
and (F_abcd F^abcd)^2 as the omega-divergence cure.

Written blind: the producer's script and JSON are NOT read before this
runs.  Every number below comes from an implementation written here.

EQUATIONS FIRST
---------------
(1) The notebook density (both PDFs, y = 0 plane, g = m = 1 after the
    /(m^4 g^4) normalization).  With
        A = d^2 + x^2 - 2 d z + z^2 = x^2 + (z - d)^2
        B = d^2 + x^2 + 2 d z + z^2 = x^2 + (z + d)^2
        S = sqrt(A) sqrt(B)
        N = d^8
          + d^6 (2 x^2 - 2 z^2 + S)
          + (x^2 + z^2)^3 (x^2 + z^2 + S)
          + d^2 (x^4 - z^4) (2 x^2 + 2 z^2 + S)
          + d^4 (4 x^4 + 2 z^4 - z^2 S + x^2 (4 z^2 + S))
    the 2026-08-17 notebook's Out[13] is
        Hs_1708 = +8 g^4 m^4 N / (A^3 B^3)
    and the 2026-08-29 notebook's Out[12] is
        Hs_2908 = -8 g^4 m^4 N / (A^3 B^3)      (the ONLY change)
    with the same energy functional
        en(d) = Int_0^inf dz Int_0^inf dx  4 pi x Hs [A > 0.001] / (m^4 g^4)
    and the same least-squares fit en(d) ~ a + b/d on d = 0.1 .. 3.0.

(2) The certified static action (registry m5_32_lagrangian.py):
        F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu,  A_i = d_i M
        <F,G>_eta = sum_{ab} eta_a eta_b F_ab G_ab,  eta = diag(-1,1,1,1)
        I1_static  = sum_{i<j} <F_ij, F_ij>_eta
        E_cert     = 4 h^3 sum_cells I1_static + V4
        V4         = w h^3 sum_cells sum_{p=1..4} (tr((M eta)^p) - C_p)^2,
                     C_p = (s g)^p + 1 + delta^p,  w = W1
        E_flip     = -4 h^3 sum_cells I1_static + V4
    Stencil: sym = 1/2 (fwd + bwd), density per branch then averaged.

(3) The boost hedgehog:
        theta(r) = m r / sqrt(r^2 + r_c^2),  nhat = x/r
        K[0,i] = K[i,0] = nhat_i,  K^2[0,0] = 1, K^2[i,j] = nhat_i nhat_j
        o      = I + sinh(theta) K + (cosh(theta) - 1) K^2   (in SO(1,3))
        M      = o M0 o^T,  M0 = diag(-s g, 1, delta, 0) = diag(32,1,0.3,0)
    o Lorentz => (M eta) = o (M0 eta) o^{-1} => V4 = 0 pointwise.
    Pair: o = o1 o2, o_k the same builder centered at z = -d and z = +d.

(4) Exact dilation identity (the no-floor argument).  With M_s(x) = M(s x),
    I1 is quartic in first derivatives, V4 is derivative-free, so
        E_u[M_s ; box L/s, n]  =  s * E_u[M ; box L, n]     (EXACT on the
                                   lattice: identical samples, h -> h/s)
        V4 [M_s ; box L/s, n]  =  s^{-3} * V4[M ; box L, n]
    hence E_flip = -E_u + V4 -> -inf like -s.  Holding L FIXED instead
    gives E_u[M_s ; L, n] = s * E_u[M ; sL, n], which is contaminated by
    the coarser grid, so the fixed-box squeeze exponent is NOT the
    dilation exponent.

(5) Legendre read of a quartic term.  L_extra = sigma c5 (s + kappa w^2)^2,
    k = kappa w^2, H = w dL/dw - L = sigma c5 (-s^2 + 2 s k + 3 k^2).

WHAT IS CHECKED (claims C1..C4 of the R11 producer)
---------------------------------------------------
C1  the two PDFs differ only by the global sign of Hs; my own polar
    quadrature of (1) reproduces 863.733 + 167.668/d, so the 08-29 fit is
    the negation.  Plus a cutoff-sensitivity probe not in the claim.
C2  the flipped static sector has no floor: exact dilation identity (4)
    measured on the lattice, the m ladder, the fixed-box squeeze ladder,
    V4 on the Lorentz orbit, h-convergence at m = 0.4, and a V4-rescue
    test on an OFF-orbit field.
C3  the same-sign pair under the certified action, my own field builder
    and my own lattice sum, at n = 48 L = 48 and a n = 64 grid probe.
C4  the R8 JSON arithmetic (H2 = -4 C2_I1, thresholds, omega*, drift,
    statics ratio), the sign convention verified by an independent
    kinetic-energy identity, V4 on the R8 field, A_static reproduced,
    and the Legendre algebra by sympy.

OUTPUTS
-------
../data/m5_32_r11_audit.json   (every number, per claim)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r11_audit.json")
T0 = time.time()

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = 0.000724023879
G_MAIN, S_MAIN, DELTA = 32.0, -1.0, 0.3
RC = 0.5


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# =====================================================================
# C1: my own quadrature of the notebook density
# =====================================================================
def hs_num(x, z, d):
    """N / (A^3 B^3) with the 8 g^4 m^4 prefactor divided out later."""
    A = x * x + (z - d) ** 2
    B = x * x + (z + d) ** 2
    S = np.sqrt(A) * np.sqrt(B)
    x2, z2 = x * x, z * z
    r2 = x2 + z2
    N = (d ** 8
         + d ** 6 * (2 * x2 - 2 * z2 + S)
         + r2 ** 3 * (r2 + S)
         + d ** 2 * (x2 * x2 - z2 * z2) * (2 * x2 + 2 * z2 + S)
         + d ** 4 * (4 * x2 * x2 + 2 * z2 * z2 - z2 * S
                     + x2 * (4 * z2 + S)))
    return N / (A ** 3 * B ** 3)


def _gauss_panels(edges, npt):
    """nodes/weights of npt-point Gauss-Legendre on each panel."""
    gx, gw = np.polynomial.legendre.leggauss(npt)
    a = edges[..., :-1]
    b = edges[..., 1:]
    mid = 0.5 * (a + b)
    half = 0.5 * (b - a)
    nodes = mid[..., None] + half[..., None] * gx
    wts = half[..., None] * gw
    sh = nodes.shape[:-2] + (nodes.shape[-2] * nodes.shape[-1],)
    return nodes.reshape(sh), wts.reshape(sh)


def en_of_d(d, cutoff=0.001, nphi_pan=16, nphi=20, nt_pan=24, nt=20,
            tfloor=1e-9):
    """en(d) for the POSITIVE (08-17) density, by polar quadrature about
    the singular center (x, z) = (0, d).

    x = rho sin(phi), z = d + rho cos(phi), phi in [0, pi] covers x >= 0;
    z > 0 caps rho at -d/cos(phi) for phi > pi/2; the Boole cutoff is
    rho > rho_c = sqrt(cutoff).  Substituting rho = rho_c / t maps the
    unbounded arm to t in (0, 1] and the integrand to
        4 pi x Hs rho_c^2 / t^3.
    """
    rc = np.sqrt(cutoff)
    # phi nodes: two regions split at pi/2 (t_min has a kink there)
    e1 = np.linspace(0.0, np.pi / 2, nphi_pan + 1)
    e2 = np.linspace(np.pi / 2, np.pi, nphi_pan + 1)
    p1, w1 = _gauss_panels(e1, nphi)
    p2, w2 = _gauss_panels(e2, nphi)
    phi = np.concatenate([p1, p2])
    wphi = np.concatenate([w1, w2])
    c = np.cos(phi)
    tmin = np.where(c >= 0.0, tfloor, np.maximum(rc * (-c) / d, tfloor))
    # geometric panels in t from tmin to 1 (per phi)
    lo = np.log(tmin)
    edges = np.exp(lo[:, None] + (0.0 - lo)[:, None]
                   * np.linspace(0.0, 1.0, nt_pan + 1)[None, :])
    tn, tw = _gauss_panels(edges, nt)          # (nphi, nt_pan*nt)
    rho = rc / tn
    xx = rho * np.sin(phi)[:, None]
    zz = d + rho * c[:, None]
    val = 4.0 * np.pi * xx * hs_num(xx, zz, d) * rc ** 2 / tn ** 3
    inner = np.sum(val * tw, axis=1)
    return 8.0 * float(np.sum(inner * wphi))


def fit_a_b(ds, en):
    Amat = np.stack([np.ones_like(ds), 1.0 / ds], axis=1)
    coef, *_ = np.linalg.lstsq(Amat, en, rcond=None)
    return float(coef[0]), float(coef[1])


def claim_c1():
    log("C1: notebook quadrature")
    res = {}
    ds = np.arange(1, 31) * 0.1
    en = np.array([en_of_d(float(d)) for d in ds])
    a, b = fit_a_b(ds, en)
    # refinement gate
    en_hi = np.array([en_of_d(float(d), nphi_pan=24, nphi=28,
                              nt_pan=36, nt=28, tfloor=1e-11)
                      for d in ds])
    a2, b2 = fit_a_b(ds, en_hi)
    res["fit_1708_coarse"] = {"a": a, "b": b}
    res["fit_1708_refined"] = {"a": a2, "b": b2}
    res["quadrature_rel_drift_a"] = abs(a2 - a) / abs(a2)
    res["quadrature_rel_drift_b"] = abs(b2 - b) / abs(b2)
    res["pdf_1708_stated"] = {"a": 863.733, "b": 167.668}
    res["pdf_2908_stated"] = {"a": -863.733, "b": -167.668}
    res["rel_vs_1708"] = {"a": abs(a2 - 863.733) / 863.733,
                          "b": abs(b2 - 167.668) / 167.668}
    res["fit_2908_predicted_by_sign_flip"] = {"a": -a2, "b": -b2}
    res["rel_vs_2908"] = {"a": abs(-a2 + 863.733) / 863.733,
                          "b": abs(-b2 + 167.668) / 167.668}
    res["en_table"] = [[float(d), float(e)] for d, e in zip(ds, en_hi)]
    # independent cross-check of one d with scipy nested quad
    from scipy.integrate import quad
    d0 = 1.0
    rc = np.sqrt(0.001)

    def inner(phi):
        cc = np.cos(phi)
        rmax = np.inf if cc >= 0 else -d0 / cc

        def f(rho):
            x = rho * np.sin(phi)
            z = d0 + rho * cc
            return 4.0 * np.pi * x * hs_num(x, z, d0) * rho
        pts = [rc, 10 * rc, 1.0, 10.0, 100.0]
        pts = [p for p in pts if p < rmax] + ([rmax] if np.isfinite(rmax)
                                              else [np.inf])
        tot = 0.0
        for aa, bb in zip(pts[:-1], pts[1:]):
            tot += quad(f, aa, bb, limit=200)[0]
        return tot
    sq = (quad(inner, 0.0, np.pi / 2, limit=200)[0]
          + quad(inner, np.pi / 2, np.pi, limit=200)[0]) * 8.0
    res["scipy_crosscheck_d1"] = {"scipy": float(sq),
                                  "mine": float(en_of_d(1.0)),
                                  "rel": abs(sq - en_of_d(1.0)) / abs(sq)}
    # singularity exponent + cutoff sensitivity (NOT part of the claim)
    rr = np.array([1e-3, 3e-3, 1e-2, 3e-2])
    hv = np.array([hs_num(r / np.sqrt(2.0), 1.0 + r / np.sqrt(2.0), 1.0)
                   for r in rr])
    res["core_exponent_of_Hs_at_d1"] = float(
        np.polyfit(np.log(rr), np.log(np.abs(hv)), 1)[0])
    cut = {}
    for cc in (0.0005, 0.001, 0.002):
        e = np.array([en_of_d(float(d), cutoff=cc) for d in ds])
        aa, bb = fit_a_b(ds, e)
        cut[str(cc)] = {"a": aa, "b": bb}
    res["cutoff_sensitivity"] = cut
    return res


# =====================================================================
# my own lattice layer (independent of the registry)
# =====================================================================
def coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def d1_mine(f, ax, h, br):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim

    def at(i):
        s = list(sl)
        s[ax] = i
        return tuple(s)
    if br == "fwd":
        out[at(slice(0, -1))] = (f[at(slice(1, None))]
                                 - f[at(slice(0, -1))]) / h
    elif br == "bwd":
        out[at(slice(1, None))] = (f[at(slice(1, None))]
                                   - f[at(slice(0, -1))]) / h
    else:
        raise ValueError(br)
    return out


def i1_static_sum(M, h):
    """sum_cells sum_{i<j} <F_ij, F_ij>_eta, sym stencil (branch average)."""
    ee = np.diag(ETA)
    w = ee[:, None] * ee[None, :]
    tot = 0.0
    for br in ("fwd", "bwd"):
        A = [d1_mine(M, ax, h, br) for ax in range(3)]
        Ae = [a @ ETA for a in A]
        s = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                F = Ae[i] @ A[j] - Ae[j] @ A[i]
                s += np.sum(w * F * F)
        tot += 0.5 * s
    return float(tot)


def v4_sum(M, s=S_MAIN, g=G_MAIN, delta=DELTA, w=W1):
    Me = M @ ETA
    P = Me.copy()
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    sg = s * g
    cp = [sg ** k + 1.0 + delta ** k for k in range(1, 5)]
    vd = sum((t[k] - cp[k]) ** 2 for k in range(4))
    return float(w * np.sum(vd)), float(w * np.max(vd))


def e_cert(M, h, sign=+1.0):
    """sign = +1 certified, sign = -1 the flipped curvature sector."""
    u = i1_static_sum(M, h)
    v, vmax = v4_sum(M)
    return (sign * 4.0 * h ** 3 * u + h ** 3 * v, 4.0 * h ** 3 * u,
            h ** 3 * v, h ** 3 * vmax)


def boost_o(X, Y, Z, m, rc, cx=0.0, cy=0.0, cz=0.0):
    x, y, z = X - cx, Y - cy, Z - cz
    r = np.sqrt(x * x + y * y + z * z)
    safe = np.where(r > 0, r, 1.0)
    nx, ny, nz = x / safe, y / safe, z / safe
    th = m * r / np.sqrt(r * r + rc * rc)
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = np.where(r > 0, 1.0, 0.0)
    for i, a in enumerate((nx, ny, nz)):
        for j, b in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * b
    return (np.eye(4) + np.sinh(th)[..., None, None] * K
            + (np.cosh(th) - 1.0)[..., None, None] * K2)


def M0_vac():
    return np.diag([-S_MAIN * G_MAIN, 1.0, DELTA, 0.0])


def hedgehog(n, L, m, rc=RC, centers=((0.0, 0.0, 0.0),), scale=1.0):
    """M = o M0 o^T with o the ordered product over centers; `scale`
    multiplies the coordinates (the squeeze)."""
    h = L / n
    X, Y, Z = coords(n, h)
    X, Y, Z = scale * X, scale * Y, scale * Z
    o = None
    for (cx, cy, cz) in centers:
        ok = boost_o(X, Y, Z, m, rc, cx, cy, cz)
        o = ok if o is None else o @ ok
    M = o @ M0_vac() @ o.swapaxes(-1, -2)
    return 0.5 * (M + M.swapaxes(-1, -2)), h


def orbit_defect(M):
    """max |spectrum(M eta) - spectrum(M0 eta)| over cells (Lorentz-orbit
    membership, the thing that forces V4 = 0)."""
    Me = M @ ETA
    ev = np.sort_complex(np.linalg.eigvals(Me.reshape(-1, 4, 4)))
    ref = np.sort_complex(np.linalg.eigvals(M0_vac() @ ETA))
    return float(np.max(np.abs(ev - ref[None, :])))


# =====================================================================
# C2
# =====================================================================
def claim_c2(LAG):
    log("C2: flipped sector no-floor")
    res = {}
    # --- gate: my I1 == the registry I1 on a random field
    rng = np.random.default_rng(1132)
    B3 = LAG.B3
    cfg = B3.base_cfg(n=10, L=15.0, s=S_MAIN, g=G_MAIN)
    Mr = B3.vac4(cfg) + 0.4 * B3.sym4(rng.normal(size=(10, 10, 10, 4, 4)))
    mine = 4.0 * cfg["h"] ** 3 * i1_static_sum(Mr, cfg["h"])
    reg = float(LAG.term_energy(LAG.REGISTRY["I1"], Mr, cfg,
                                LAG.default_params(s=S_MAIN, g=G_MAIN),
                                -4.0))
    eu_cert, ev_cert = B3.e_parts(Mr, cfg)
    v_mine = cfg["h"] ** 3 * v4_sum(Mr)[0]
    res["gate_I1_mine_vs_registry_rel"] = abs(mine - reg) / abs(reg)
    res["gate_I1_mine_vs_certified_e_parts_rel"] = \
        abs(mine - float(eu_cert)) / abs(float(eu_cert))
    res["gate_V4_mine_vs_certified_rel"] = \
        abs(v_mine - float(ev_cert)) / abs(float(ev_cert))

    # --- the m ladder, both grids, holding L = 32
    ms = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6]
    lad = {}
    for (n, L) in ((32, 32.0), (48, 32.0)):
        rows = []
        for m in ms:
            M, h = hedgehog(n, L, m)
            ef, eu, v4, v4max = e_cert(M, h, sign=-1.0)
            rows.append({"m": m, "E_flip": ef, "E_u_certified": eu,
                         "V4": v4, "V4_max_cell": v4max,
                         "orbit_defect": orbit_defect(M) if n == 32 and
                         m in (0.05, 1.6) else None})
            log(f"  n={n} L={L:g} m={m:g} E_flip {ef:.6g} V4 {v4:.3e}")
        sl = np.polyfit(np.log([r["m"] for r in rows[:3]]),
                        np.log([-r["E_flip"] for r in rows[:3]]), 1)[0]
        lad[f"n{n}_L{int(L)}"] = {"rows": rows,
                                  "small_m_exponent_first3": float(sl)}
    res["m_ladder"] = lad

    # --- h convergence at m = 0.4, fixed L = 32
    hconv = []
    for n in (32, 48, 64, 96, 128):
        M, h = hedgehog(n, 32.0, 0.4)
        ef, eu, v4, _ = e_cert(M, h, sign=-1.0)
        hconv.append({"n": n, "h": h, "E_flip": ef, "V4": v4})
        log(f"  hconv n={n} h={h:.4g} E_flip {ef:.6g}")
    res["h_convergence_m0.4_L32"] = hconv
    res["h_convergence_note"] = (
        "r_c = 0.5 is below h at n = 32 (h = 1.0) and n = 48 (h = 0.667); "
        "the core is only resolved from n = 96 (h = 0.333) upward")

    # --- WHERE the curvature energy lives: grow the box at FIXED h
    #     (far-field weight) vs refine h at FIXED box (core weight) vs
    #     vary r_c at a resolved h (the core scaling law)
    box_fixed_h = []
    for (n, L) in ((32, 32.0), (48, 48.0), (64, 64.0), (96, 96.0)):
        M, h = hedgehog(n, L, 0.4)
        box_fixed_h.append({"n": n, "L": L, "h": h,
                            "E_u": e_cert(M, h)[1]})
    rc_scan = []
    for rc in (0.5, 1.0, 2.0, 4.0):
        M, h = hedgehog(128, 32.0, 0.4, rc=rc)
        rc_scan.append({"r_c": rc, "h": h, "E_u": e_cert(M, h)[1]})
    res["dominance"] = {
        "box_at_fixed_h": box_fixed_h,
        "box_growth_rel_over_3x_L": abs(box_fixed_h[-1]["E_u"]
                                        - box_fixed_h[0]["E_u"])
        / box_fixed_h[0]["E_u"],
        "rc_scan_at_h0.25": rc_scan,
        "rc_exponent": float(np.polyfit(
            np.log([r["r_c"] for r in rc_scan]),
            np.log([r["E_u"] for r in rc_scan]), 1)[0]),
        "reading": ("tripling L at fixed h moves E_u by 2.6 % while "
                    "refining h at fixed L moves it by 34 %, and E_u ~ "
                    "r_c^-1: the curvature energy is CORE dominated, not "
                    "far-field dominated")}

    # --- squeeze ladder: fixed box (the producer's protocol)
    ss = [1.0, 2.0, 3.0, 4.0]
    fixed, ex_fixed = {}, {}
    for n in (32, 48):
        rows = []
        for s in ss:
            M, h = hedgehog(n, 32.0, 0.4, scale=s)
            ef, eu, v4, _ = e_cert(M, h, sign=-1.0)
            rows.append({"s": s, "E_flip": ef})
        fixed[f"n{n}"] = rows
        ex_fixed[f"n{n}"] = float(np.polyfit(
            np.log(ss), np.log([-r["E_flip"] for r in rows]), 1)[0])
    # --- squeeze ladder: box scaled with the squeeze (exact dilation)
    scaled = []
    for s in ss:
        M, h = hedgehog(32, 32.0 / s, 0.4, scale=s)
        ef, eu, v4, _ = e_cert(M, h, sign=-1.0)
        scaled.append({"s": s, "E_flip": ef, "E_u": eu, "V4": v4})
    ex_scaled = float(np.polyfit(np.log(ss),
                                 np.log([-r["E_flip"] for r in scaled]),
                                 1)[0])
    res["squeeze_fixed_box"] = {"rows_by_grid": fixed,
                                "exponent_by_grid": ex_fixed}
    res["squeeze_scaled_box"] = {"rows": scaled, "exponent": ex_scaled}

    # --- V4 rescue test on an OFF-orbit field (V4 != 0)
    n, L = 32, 32.0
    M, h = hedgehog(n, L, 0.4)
    Moff = M.copy()
    Moff[..., 1, 1] += 0.7          # break the spectrum -> V4 > 0
    off = []
    for s in ss:
        # dilate the off-orbit field exactly: same samples, h -> h/s
        ef, eu, v4, _ = e_cert(Moff, h / s, sign=-1.0)
        off.append({"s": s, "E_flip": ef, "E_u": eu, "V4": v4})
    res["v4_rescue_offorbit_dilation"] = {
        "rows": off,
        "E_u_exponent": float(np.polyfit(np.log(ss),
                                         np.log([r["E_u"] for r in off]),
                                         1)[0]),
        "V4_exponent": float(np.polyfit(np.log(ss),
                                        np.log([r["V4"] for r in off]),
                                        1)[0])}
    return res


# =====================================================================
# C3
# =====================================================================
def claim_c3():
    log("C3: same-sign pair under the certified action")
    res = {}
    for (n, L) in ((48, 48.0), (64, 48.0)):
        rows = []
        Ms, h = hedgehog(n, L, 0.1)
        e_single = e_cert(Ms, h, sign=+1.0)
        for d in (8.0, 10.0, 12.0):
            Mp, h = hedgehog(n, L, 0.1,
                             centers=((0.0, 0.0, -d), (0.0, 0.0, d)))
            ec, eu, v4, v4max = e_cert(Mp, h, sign=+1.0)
            ef = -eu + v4
            rows.append({"d": d, "E_pair_certified": ec,
                         "E_int_certified": ec - 2.0 * e_single[0],
                         "E_pair_flipped": ef,
                         "E_int_flipped": ef - 2.0 * (-e_single[1]
                                                      + e_single[2]),
                         "V4": v4, "V4_max_cell": v4max,
                         "orbit_defect": orbit_defect(Mp) if n == 48 and
                         d == 8.0 else None})
            log(f"  n={n} d={d:g} E_pair {ec:.6g} "
                f"E_int {rows[-1]['E_int_certified']:.6g} V4 {v4:.3e}")
        dd = np.array([r["d"] for r in rows])
        ee = np.array([r["E_int_certified"] for r in rows])
        A = np.stack([np.ones_like(dd), 1.0 / dd], axis=1)
        coef, *_ = np.linalg.lstsq(A, ee, rcond=None)
        res[f"n{n}_L{int(L)}"] = {
            "E_single_certified": e_single[0],
            "E_single_u": e_single[1], "E_single_V4": e_single[2],
            "rows": rows,
            "fit_a_plus_b_over_d": {"a": float(coef[0]), "b": float(coef[1])},
            "dE_dd_sign": "positive (attractive)"
            if ee[-1] > ee[0] else "negative (repulsive)"}
    return res


# =====================================================================
# C4
# =====================================================================
def claim_c4(LAG):
    log("C4: (F F)^2 as the omega cure")
    res = {}
    with open(os.path.join(DATA, "m5_32_r8_quartics.json")) as f:
        r8 = json.load(f)
    q = r8["scaling"]["Q_I1sq"]
    i1 = r8["scaling"]["I1"]
    lad = r8["ladder"]["per_term"]["Q_I1sq"]
    Ls = [b["L"] for b in r8["boxes"]]
    res["r8_ladder_note"] = r8["ladder"]["note"]
    res["Q_I1sq_definition"] = q["definition"]

    # (a) H2 = -4 C2_I1
    h2_mine = [-4.0 * c for c in i1["C2"]]
    res["H2_check"] = {"mine": h2_mine, "stored": lad["certified_inertia_H2"],
                       "max_rel": max(abs(a - b) / abs(b) for a, b in
                                      zip(h2_mine, lad["certified_inertia_H2"]))}
    # (b) c5 threshold
    c5_mine = [h / abs(c) for h, c in zip(h2_mine, q["C2"])]
    res["c5_threshold_check"] = {
        "mine": c5_mine, "stored": lad["c5_required_to_flip_H2"],
        "max_rel": max(abs(a - b) / abs(b) for a, b in
                       zip(c5_mine, lad["c5_required_to_flip_H2"])),
        "exponent_in_L": float(np.polyfit(np.log(Ls), np.log(c5_mine), 1)[0])}
    # (c) omega* at c5 = 1.5 * c5_threshold(L = 96), ENERGY reading
    c5 = 1.5 * c5_mine[-1]
    om_energy, om_legendre = [], []
    for k in range(3):
        num = -(h2_mine[k] + c5 * q["C2"][k])
        om_energy.append(float(np.sqrt(num / (2.0 * c5 * q["C4"][k]))))
        om_legendre.append(float(np.sqrt(num / (2.0 * 3.0 * c5
                                                * q["C4"][k]))))
    res["omega_star"] = {
        "c5_used": c5,
        "energy_reading": om_energy,
        "energy_reading_drift":
            (om_energy[0] - om_energy[-1]) / om_energy[0],
        "legendre_reading_w4_times_3": om_legendre,
        "legendre_reading_drift":
            (om_legendre[0] - om_legendre[-1]) / om_legendre[0]}
    # (d) statics deformation at threshold
    dfm = [c5_mine[k] * q["A_static"][k] / (4.0 * i1["A_static"][k])
           for k in range(3)]
    res["statics_deformation_ratio"] = {
        "mine": dfm,
        "sign_note": ("with E_contrib = c (C w^2 - A) the c5 term's STATIC "
                      "energy contribution is -c5 A_static, i.e. it does "
                      "not deform the hedgehog's static energy, it drives "
                      "the total NEGATIVE (62.85 -> -470.8 at L = 48)"),
        "total_static_energy_at_threshold":
            [4.0 * i1["A_static"][k] - c5_mine[k] * q["A_static"][k]
             for k in range(3)]}
    # (e) the sign convention, verified independently: for the certified
    #     I1 the energy's w^2 coefficient must equal 4 * kin_of
    B3 = LAG.B3
    rng = np.random.default_rng(1134)
    cfg = B3.base_cfg(n=12, L=18.0, s=S_MAIN, g=G_MAIN)
    Mr = B3.vac4(cfg) + 0.4 * B3.sym4(rng.normal(size=(12, 12, 12, 4, 4)))
    a0 = B3.sym4(rng.normal(size=Mr.shape))
    p = LAG.default_params(s=S_MAIN, g=G_MAIN)
    A0, B0, C0 = LAG.omega_decompose(LAG.REGISTRY["I1"], Mr, cfg, p, a0)
    kin = float(B3.kin_of(Mr, a0, cfg))
    # B3.kin_of already carries the factor 4 (E = 4(U + w^2 T) + V4), so
    # the identity to test is  H2 = -4 C2_I1 == kin_of  (NOT 4 kin_of)
    res["sign_convention_check"] = {
        "C2_lagrangian": float(C0), "minus4_C2": float(-4.0 * C0),
        "certified_kin_of": kin,
        "rel": abs(-4.0 * C0 - kin) / abs(kin),
        "static_A": float(A0), "4A_vs_certified_Eu_rel":
            abs(4.0 * float(A0) - float(B3.e_parts(Mr, cfg)[0]))
            / abs(float(B3.e_parts(Mr, cfg)[0]))}
    # (f) the R8 field: is it a Lorentz orbit of the vacuum (V4 = 0)?
    #     rebuilt here independently (rotation hedgehog, m = 0)
    n, L = 32, 48.0
    h = L / n
    X, Y, Z = coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    th = -np.arctan2(Z, rho)

    def rot(gen, ang):
        G2m = gen @ gen
        return (np.eye(4) + np.sin(ang)[..., None, None] * gen
                + (1 - np.cos(ang))[..., None, None] * G2m)
    G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0
    G2 = np.zeros((4, 4)); G2[1, 3] = 1.0; G2[3, 1] = -1.0
    G3 = np.zeros((4, 4)); G3[1, 2] = -1.0; G3[2, 1] = 1.0
    Q = rot(G3, phi) @ rot(G2, th) @ rot(G1, np.zeros_like(phi))
    Mr8 = Q @ M0_vac() @ Q.swapaxes(-1, -2)
    Mr8 = 0.5 * (Mr8 + Mr8.swapaxes(-1, -2))
    v4tot, v4max = v4_sum(Mr8)
    a_static_mine = i1_static_sum(Mr8, h) * h ** 3
    res["r8_field_check"] = {
        "field": "dressed(cfg, m = 0): a pure SO(3) rotation hedgehog, "
                 "NOT a relaxed/dressed electron",
        "V4_total": h ** 3 * v4tot, "V4_max_cell": h ** 3 * v4max,
        "orbit_defect": orbit_defect(Mr8),
        "A_static_I1_mine": a_static_mine,
        "A_static_I1_stored": i1["A_static"][0],
        "rel": abs(a_static_mine - i1["A_static"][0]) / i1["A_static"][0]}
    # (g) the Legendre algebra, sympy
    c5s, ss_, ks, w, sig = sp.symbols("c5 s kappa omega sigma")
    Lx = sig * c5s * (ss_ + ks * w ** 2) ** 2
    H = sp.expand(w * sp.diff(Lx, w) - Lx)
    k = ks * w ** 2
    claim = sp.expand(sig * c5s * (-ss_ ** 2 + 2 * ss_ * k + 3 * k ** 2))
    res["legendre_algebra"] = {
        "H": str(sp.simplify(H)),
        "matches_producer_form": bool(sp.simplify(H - claim) == 0),
        "sigma_plus1_static_coeff": str(sp.expand(H.subs(sig, 1)).coeff(w, 0)),
        "sigma_minus1_w4_coeff": str(sp.expand(H.subs(sig, -1)).coeff(w, 4))}
    return res


# =====================================================================
def main():
    LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
    out = {"task": "M5.32 R11 adversarial audit (independent)",
           "python": sys.executable,
           "conventions": {"eta": "diag(-1,1,1,1)", "s": S_MAIN,
                           "g": G_MAIN, "delta": DELTA, "w": W1,
                           "r_c": RC, "stencil": "sym (fwd/bwd average)"}}
    out["C1"] = claim_c1()
    out["C2"] = claim_c2(LAG)
    out["C3"] = claim_c3()
    out["C4"] = claim_c4(LAG)
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    log(f"wrote {OUT}  ({out['runtime_s']} s)")


if __name__ == "__main__":
    main()
