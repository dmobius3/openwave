"""M5.32 R1 arm (c): the lattice TIME-SECTOR screen of the covariant term
basis at the certified toy point (g = 32, delta = 0.3, s = -1, n = 32,
L = 48, sym stencil).

EQUATIONS FIRST
---------------
Registry (m5_32_lagrangian.py, imported, never modified): per term T the
Lagrangian density on jets A_mu, F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu,
    I1 = (1/2) F_{mu nu a b} F^{mu nu a b}   (the certified action)
    I1_frob = I1 with the Frobenius internal metric (variant A; control)
    I2 = F_{mu nu a b} F^{a b mu nu},  I3 = F_{mu nu a b} F^{mu a nu b}
    I4 = R_{nu a} R^{nu a},  I5 = R_{nu a} R^{a nu},  I6 = R^2,
    R_{nu a} = F_{mu nu a}^{mu}.
The lattice read (h^3-weighted sum over cells, both sym-stencil branches):
    L_T(omega) = h^3 sum_cells T(A_i = d_i M, A_0 = omega a0)
               = A_T + B_T omega + C_T omega^2       (exact, 3-point)
    A_T = L_T(0), B_T = (L_T(1) - L_T(-1))/2, C_T = (L_T(1)+L_T(-1))/2 - A_T
(the module's omega_decompose). Legendre: H_T = C_T omega^2 - A_T, the
omega-linear B_T drops out of the Hamiltonian.

Certified normalization (the FACTOR-4 BRIDGE): E_cert = -4 H_I1 + V4, so
the certified per-channel kin of M5.21.16 CHAN is
    kin_I1(channel) = -4 C_I1(channel)          (E_kin = kin omega^2).
Every term is reported in the SAME normalization,
    C_E(T) := -4 C_T   (energy-read kinetic coefficient; C_E >= 0 stable)
and the modified action is
    L'(c) = -4 (I1 + c T) - V4,   E'(c) = -4 (H_I1 + c H_T) + V4
so the kinetic coefficient of E'(c) on a channel is LINEAR in c,
    K(c) = C_E(I1) + c C_E(T).
Stability interval of T: {c : K(c) >= 0 on EVERY channel}; each channel
contributes one half-line (c >= -C_E(I1)/C_E(T) if C_E(T) > 0, c <= ... if
C_E(T) < 0, unconstrained / empty if C_E(T) = 0 with C_E(I1) >= 0 / < 0).
Boost-sign reversal: K(c) > 0 on boost_z and boost_x (I1's is < 0).

Static energy read of the dressed / bump / family configurations under
L'(c):  E'_static(c) = +4 (A_I1 + c A_T) + V4  (H_T = -A_T, times the
coefficient -4; V4 = 0 on every boost conjugation of the family; checked).

Kernel order (M5.21.16 s 3.5 construction, m5_21_16_d_twobody.py): the
vacuum dressed by one radial boost bump of amplitude eps,
    b(x) = eps (r/rho) exp(-(r/rho)^2), rho = 3, M = Qb M_vac Qb^T,
    ratio(T) = A_T(2 eps) / A_T(eps)   (4 = quadratic kernel, 16 = quartic)
    order(T) = log2 ratio / 1.
Two bumps at +-d/2 z-hat:  E_int(d) = E(pair) - E(single at +d/2)
- E(single at -d/2) on the same box; convention F = -dE_int/dd, dE_int/dd
< 0 = REPULSIVE. The omega-carrying bump: A_0 = omega a0 with a0 = Jz M -
M Jz (the rot_z clock generator acting on the dressed vacuum; C_T is the
omega^2 coefficient of the bump, its interaction sign read the same way).

Covariance (item 5): (a) the certified INV stage of m5_21_16_b_field.py
(internal conjugation only, M -> L M L^T, derivative index NOT rotated):
exact for I1 (derivative indices contract among themselves) and NOT a
symmetry of the mixed terms, so the certified test is reported as-is and
(b) the full action of the module's selftest is applied to the LATTICE
jets, A'_mu = (Lambda^{-T})_mu^nu Lambda A_nu Lambda^T, and the drift of
the h^3 sum of densities is the lattice covariance number; (c) the no-eta
controls (must FAIL under boosts): I1 with the derivative pair contracted
by delta, I2 on the plain-commutator bracket F = A_mu A_nu - A_nu A_mu,
and the registry controls I1_frob, I3_mixed_eta.

Stages: CHAN, DRESSED, KERNEL, BOUNDED, LADDER (the h-ladder control of
BOUNDED, I1 only), INV (all, or `--stage NAME`).
Out: ../data/m5_32_r1_screen.json, ../plots/m5_32_r1_screen_chan.png
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
B3 = LAG.B3
B16 = _load("m5_21_16_b_field", "m5_21_16_b_field.py")
C14 = B16.C14
B8 = C14.B8
D16 = _load("m5_21_16_d_twobody", "m5_21_16_d_twobody.py")
REG = LAG.REGISTRY

TERMS = ["I1", "I1_frob", "I2", "I3", "I4", "I5", "I6"]
COEF = -4.0                  # the certified normalization of I1
G_MAIN, S_MAIN = 32.0, -1.0
TOL = 1e-12
T_START = time.time()


def log(msg):
    print(f"[{time.time() - T_START:8.1f}s] {msg}", flush=True)


def abc(term, M, cfg, p, a0):
    A, B, C = LAG.omega_decompose(REG[term], M, cfg, p, a0)
    return {"A": float(A), "B": float(B), "C": float(C),
            "C_E": float(COEF * C)}


def static(term, M, cfg, p):
    return float(LAG.term_lagrangian(REG[term], M, cfg, p))


def interval(kin1, kinT, tol=TOL):
    """{c : kin1[ch] + c kinT[ch] >= 0 for all ch}; per-channel half-lines."""
    lo, hi, empty = -np.inf, np.inf, False
    per = {}
    for ch in kin1:
        k1, kT = kin1[ch], kinT[ch]
        if abs(kT) <= tol:
            per[ch] = "unconstrained" if k1 >= -tol else "EMPTY (kT = 0, k1 < 0)"
            if k1 < -tol:
                empty = True
        elif kT > 0:
            per[ch] = f"c >= {-k1 / kT:.6g}"
            lo = max(lo, -k1 / kT)
        else:
            per[ch] = f"c <= {-k1 / kT:.6g}"
            hi = min(hi, -k1 / kT)
    ok = (not empty) and lo <= hi
    out = {"per_channel": per, "lo": float(lo), "hi": float(hi),
           "nonempty": bool(ok)}
    if ok:
        # an interior test point
        if np.isfinite(lo) and np.isfinite(hi):
            ct = 0.5 * (lo + hi)
        elif np.isfinite(lo):
            ct = lo + 1.0
        elif np.isfinite(hi):
            ct = hi - 1.0
        else:
            ct = 0.0
        kin_ct = {ch: kin1[ch] + ct * kinT[ch] for ch in kin1}
        boosts = [ch for ch in kin1 if ch.startswith("boost")]
        out["test_c"] = float(ct)
        out["kin_at_test_c"] = kin_ct
        out["boost_reversed_at_test_c"] = bool(
            all(kin_ct[ch] > tol for ch in boosts)) if boosts else None
        out["all_stable_at_test_c"] = bool(
            all(kin_ct[ch] >= -tol for ch in kin1))
    return out


def endpoints_for_probe(iv):
    """the two c values of the boundedness probe (item 4)."""
    if not iv["nonempty"]:
        return [-0.5, 0.5], "interval empty: +-0.5"
    lo, hi = iv["lo"], iv["hi"]
    if np.isfinite(lo) and np.isfinite(hi):
        return [lo, hi], "interval endpoints"
    if np.isfinite(lo):
        return [lo, lo + 4.0], "lower endpoint and lower + 4 (upper infinite)"
    if np.isfinite(hi):
        return [hi - 4.0, hi], "upper - 4 and upper endpoint (lower infinite)"
    return [-0.5, 0.5], "unconstrained: +-0.5"


# ---------------- CHAN ----------------
def stage_chan():
    cfg = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=32, L=48.0)
    p = LAG.default_params(s=S_MAIN, g=G_MAIN)
    M = B16.lattice_family_M(cfg, G_MAIN)
    a0s = B3.gen_catalog(cfg, M)
    chans = list(a0s.keys())
    with open(os.path.join(DATA, "m5_21_16_field.json")) as f:
        stored = json.load(f)["CHAN"]["rows"]
    rows = {}
    for t in TERMS:
        rows[t] = {}
        for ch in chans:
            rows[t][ch] = abc(t, M, cfg, p, a0s[ch])
        log(f"CHAN {t}: C_E = " + ", ".join(
            f"{ch} {rows[t][ch]['C_E']:+.4f}" for ch in chans))
    # pointwise omega-linear content: is B zero per cell or only after
    # the integral? max |B density| vs max |C density| per term, channel
    pw = {}
    for t in TERMS:
        pw[t] = {}
        for ch in chans:
            dp = dm = d0 = 0.0
            for (A, wt), (Ap, _), (Am, _) in zip(
                    LAG.lattice_jets(M, cfg),
                    LAG.lattice_jets(M, cfg, a0s[ch], 1.0),
                    LAG.lattice_jets(M, cfg, a0s[ch], -1.0)):
                d0 = d0 + wt * REG[t].density(A, M, p)
                dp = dp + wt * REG[t].density(Ap, M, p)
                dm = dm + wt * REG[t].density(Am, M, p)
            Bd = 0.5 * (dp - dm)
            Cd = 0.5 * (dp + dm) - d0
            pw[t][ch] = {"max_abs_B_density": float(np.max(np.abs(Bd))),
                         "max_abs_C_density": float(np.max(np.abs(Cd))),
                         "sum_B_density_h3": float(cfg["h"] ** 3 * Bd.sum())}
        log(f"CHAN pointwise {t}: max|B| = " + ", ".join(
            f"{ch} {pw[t][ch]['max_abs_B_density']:.2e}" for ch in chans))
    # controls vs the stored certified CHAN table
    ctrl = {}
    for ch in chans:
        ctrl[ch] = {
            "I1_C_E": rows["I1"][ch]["C_E"],
            "stored_kin_eta": stored[ch]["kin_eta"],
            "rel_I1": LAG._rel(rows["I1"][ch]["C_E"], stored[ch]["kin_eta"]),
            "I1_frob_C_E": rows["I1_frob"][ch]["C_E"],
            "stored_kin_flip": stored[ch]["kin_flip"],
            "rel_frob": LAG._rel(rows["I1_frob"][ch]["C_E"],
                                 stored[ch]["kin_flip"])}
    worst = max(max(v["rel_I1"], v["rel_frob"]) for v in ctrl.values())
    kin1 = {ch: rows["I1"][ch]["C_E"] for ch in chans}
    ivs = {}
    for t in TERMS:
        kinT = {ch: rows[t][ch]["C_E"] for ch in chans}
        ivs[t] = interval(kin1, kinT)
        log(f"CHAN interval {t}: [{ivs[t]['lo']:.4g}, {ivs[t]['hi']:.4g}] "
            f"nonempty={ivs[t]['nonempty']} "
            f"boost_reversed={ivs[t].get('boost_reversed_at_test_c')}")
    return {"cfg": {k: cfg[k] for k in ("n", "L", "h", "g", "s", "delta",
                                        "stencil")},
            "channels": chans,
            "rows": rows,
            "control_vs_stored_CHAN": ctrl,
            "control_worst_rel": float(worst),
            "control_pass": bool(worst <= 1e-12),
            "C_E_table": {t: {ch: rows[t][ch]["C_E"] for ch in chans}
                          for t in TERMS},
            "B_table": {t: {ch: rows[t][ch]["B"] for ch in chans}
                        for t in TERMS},
            "B_pointwise": pw,
            "stability_intervals": ivs,
            "reading": ("K(c) = C_E(I1) + c C_E(T) per channel; the "
                        "interval is the intersection of the six "
                        "half-lines; boost_reversed = both boost "
                        "channels strictly positive at the interior "
                        "test point")}


# ---------------- DRESSED ----------------
def dressed_pair(cfg, scale, rs, bstar):
    """(M, a0) of the b*-dressed electron on the lattice, b -> scale*b*
    (the M5.21.14 lattice cross-check construction, verbatim)."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = scale * np.interp(R.ravel(), rs, bstar).reshape(R.shape)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None]
          * K + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    Md = B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
    a0d = B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, a0b, Qb))
    return Md, a0d


def stage_dressed(chan_ivs=None):
    cfg = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=32, L=48.0)
    p = LAG.default_params(s=S_MAIN, g=G_MAIN)
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        rec = json.load(f)
    rs, bstar = np.array(rec["rs"]), np.array(rec["b_star"])
    out = {"cfg": {k: cfg[k] for k in ("n", "L", "h", "g", "s", "delta")},
           "b_star_source": "m5_21_14_minimize.json (rs, b_star)",
           "scales": {}}
    for scale in (0.0, 1.0, 2.0):
        Md, a0d = dressed_pair(cfg, scale, rs, bstar)
        row = {"V4": float(LAG.term_lagrangian(REG["V4"], Md, cfg, p))}
        for t in TERMS:
            row[t] = abc(t, Md, cfg, p, a0d)
        out["scales"][f"b_{scale:g}bstar"] = row
        log(f"DRESSED scale {scale:g}: " + ", ".join(
            f"{t} C_E {row[t]['C_E']:+.3f}" for t in TERMS)
            + f", V4 {row['V4']:.2e}")
    # controls vs the M5.21.14 lattice cross-check record
    lat = rec["lattice_crosscheck"]
    r0, r1 = out["scales"]["b_0bstar"], out["scales"]["b_1bstar"]
    kin_base = r0["I1"]["C_E"]
    kin_corr = r1["I1"]["C_E"] - r0["I1"]["C_E"]
    e_corr = -COEF * (r1["I1"]["A"] - r0["I1"]["A"]) + r1["V4"] - r0["V4"]
    out["control_vs_m5_21_14_lattice"] = {
        "kin_base": kin_base, "stored_kin_base": lat["kin_base_lattice"],
        "rel_kin_base": LAG._rel(kin_base, lat["kin_base_lattice"]),
        "kin_corr": kin_corr, "stored_kin_corr": lat["kin_corr_lattice"],
        "rel_kin_corr": LAG._rel(kin_corr, lat["kin_corr_lattice"]),
        "E_corr": e_corr, "stored_E_corr": lat["E_corr_lattice"],
        "rel_E_corr": LAG._rel(e_corr, lat["E_corr_lattice"])}
    out["control_pass"] = bool(max(
        out["control_vs_m5_21_14_lattice"][k]
        for k in ("rel_kin_base", "rel_kin_corr", "rel_E_corr")) <= 1e-9)
    # per-term sign of C_E and the c-interval on the dressed electron
    ivs = {}
    for key, row in out["scales"].items():
        kin1 = {"dressed": row["I1"]["C_E"]}
        ivs[key] = {}
        for t in TERMS:
            iv = interval(kin1, {"dressed": row[t]["C_E"]})
            iv["C_E_sign"] = ("runaway (C_E < 0)" if row[t]["C_E"] < 0
                              else "stable (C_E >= 0)")
            ivs[key][t] = iv
    out["stability_intervals"] = ivs
    # table: value (A), A/B/C per term at each scale
    out["table"] = {key: {t: {k: row[t][k] for k in ("A", "B", "C", "C_E")}
                          for t in TERMS}
                    for key, row in out["scales"].items()}
    return out


# ---------------- KERNEL ----------------
def stage_kernel():
    cfg = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=48, L=72.0)
    p = LAG.default_params(s=S_MAIN, g=G_MAIN)
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0

    def a0_of(M):
        return Jz @ M - M @ Jz

    out = {"cfg": {k: cfg[k] for k in ("n", "L", "h", "g", "s", "delta")},
           "rho": D16.RHO, "amps": [0.05, 0.1, 0.2],
           "single_bump": {}, "two_bump": {}}
    for amp in out["amps"]:
        M1 = D16.dressed_vac(cfg, [(0.0, 0.0, 0.0)], amp=amp)
        a01 = a0_of(M1)
        row = {}
        for t in TERMS:
            row[t] = abc(t, M1, cfg, p, a01)
        row["V4"] = float(LAG.term_lagrangian(REG["V4"], M1, cfg, p))
        out["single_bump"][f"amp_{amp:g}"] = row
        log(f"KERNEL single bump amp {amp:g} done")
    ratios = {}
    for t in TERMS:
        a = [out["single_bump"][f"amp_{x:g}"][t]["A"] for x in out["amps"]]
        c = [out["single_bump"][f"amp_{x:g}"][t]["C"] for x in out["amps"]]
        b = [out["single_bump"][f"amp_{x:g}"][t]["B"] for x in out["amps"]]

        def rat(v):
            r = []
            for i in range(2):
                if abs(v[i]) > 1e-14 * max(1.0, abs(v[i + 1])):
                    r.append(float(v[i + 1] / v[i]))
                else:
                    r.append(None)
            return r
        ra, rc, rb = rat(a), rat(c), rat(b)
        ratios[t] = {
            "A_values": a, "A_ratio_2eps_over_eps": ra,
            "A_order_log2": [None if r is None or r <= 0 else float(np.log2(r))
                             for r in ra],
            "C_values": c, "C_ratio_2eps_over_eps": rc,
            "C_order_log2": [None if r is None or r <= 0 else float(np.log2(r))
                             for r in rc],
            "B_values": b, "B_ratio_2eps_over_eps": rb}
        log(f"KERNEL {t}: A ratios {ra}, C ratios {rc}")
    out["ratios"] = ratios
    ctrl = ratios["I1"]["A_ratio_2eps_over_eps"][0]
    with open(os.path.join(DATA, "m5_21_16_twobody.json")) as f:
        tb = json.load(f)
    out["control_I1_ratio"] = {
        "measured": ctrl, "stored": tb["kernel_amp_scaling_ratio"]["eta"],
        "rel": LAG._rel(ctrl, tb["kernel_amp_scaling_ratio"]["eta"])}
    # two-bump overlap at the record separations
    ds = [10.0, 14.0]
    two = {}
    for d in ds:
        c1, c2 = (0, 0, -d / 2), (0, 0, d / 2)
        fields = {"pair": D16.dressed_vac(cfg, [c1, c2]),
                  "a": D16.dressed_vac(cfg, [c1]),
                  "b": D16.dressed_vac(cfg, [c2])}
        vals = {k: {t: abc(t, Mf, cfg, p, a0_of(Mf)) for t in TERMS}
                for k, Mf in fields.items()}
        row = {}
        for t in TERMS:
            row[t] = {q: float(vals["pair"][t][q] - vals["a"][t][q]
                               - vals["b"][t][q]) for q in ("A", "B", "C")}
            row[t]["E_int_static"] = -COEF * row[t]["A"]      # E = +4 A
            row[t]["E_int_omega2_coeff"] = COEF * row[t]["C"]  # E_kin = -4 C
            row[t]["single_A"] = vals["a"][t]["A"]
            row[t]["single_C"] = vals["a"][t]["C"]
        two[f"d_{d:g}"] = row
        log(f"KERNEL two-bump d {d:g} done")
    # control I1 vs the stored twobody rows
    st = {r["d"]: r["E_int_eta"] for r in tb["rows"]}
    out["control_I1_twobump"] = {
        f"d_{d:g}": {"measured": two[f"d_{d:g}"]["I1"]["E_int_static"],
                     "stored": st[d],
                     "rel": LAG._rel(two[f"d_{d:g}"]["I1"]["E_int_static"],
                                     st[d])} for d in ds}
    signs = {}
    for t in TERMS:
        for kind, key in (("static", "E_int_static"),
                          ("omega", "E_int_omega2_coeff")):
            e10, e14 = two["d_10"][t][key], two["d_14"][t][key]
            ref = abs(two["d_10"][t]["single_A" if kind == "static"
                                     else "single_C"]) * abs(COEF)
            if max(abs(e10), abs(e14)) < 1e-9 * max(abs(ref), 1e-300):
                verdict = "null (below 1e-9 of the single-bump value)"
            else:
                dEdd = (e14 - e10) / 4.0
                verdict = ("repulsive (dE_int/dd < 0)" if dEdd < 0
                           else "attractive (dE_int/dd > 0)")
            signs.setdefault(t, {})[kind] = {
                "E_int_d10": e10, "E_int_d14": e14, "verdict": verdict}
    out["two_bump"] = two
    out["two_bump_signs"] = signs
    out["convention"] = ("F = -dE_int/dd; E_int in the certified "
                         "normalization (E = -4 H_T = 4 A_T - 4 C_T omega^2); "
                         "omega bump: a0 = Jz M - M Jz")
    return out


# ---------------- BOUNDED ----------------
def stage_bounded(chan):
    cfg = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=64, L=24.0)
    p = LAG.default_params(s=S_MAIN, g=G_MAIN)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    Mb = B8.dressed(cfg, 0.0)

    def dress(bl):
        Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None]
              * K + (np.cosh(bl) - 1.0)[..., None, None] * K2)
        return B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))

    family = [("base", np.zeros_like(R))]
    for amp in (0.02, 0.05, 0.1, 0.2, 0.4):
        family.append((f"wide_amp_{amp:g}", amp * np.tanh(R / 2.0)))
    for lam in (2.0, 1.0, 0.5):
        for A in (0.05, 0.1, 0.2, 0.4):
            family.append((f"saw_lam_{lam:g}_A_{A:g}", C14.saw_of(A, lam, R)))
    vals = {}
    for name, bl in family:
        Md = dress(bl)
        vals[name] = {t: static(t, Md, cfg, p) for t in TERMS}
        vals[name]["V4"] = float(LAG.term_lagrangian(REG["V4"], Md, cfg, p))
        log(f"BOUNDED {name}: I1 A = {vals[name]['I1']:.6g}, "
            f"V4 = {vals[name]['V4']:.3g}")
    out = {"cfg": {k: cfg[k] for k in ("n", "L", "h", "g", "s", "delta")},
           "family_members": [nm for nm, _ in family],
           "A_values": vals,
           "resolution_note": ("h = 0.375: the lam = 0.5 sawtooth has 1.3 "
                               "cells per half period (unresolved, "
                               "reported as-is); lam = 1 has 2.7, lam = 2 "
                               "has 5.3"),
           "per_term": {}}
    base = vals["base"]
    for t in TERMS:
        iv = chan["stability_intervals"][t]
        cs, how = endpoints_for_probe(iv)
        if t == "I1":
            cs, how = [0.0], "I1 alone (the reference action)"
        rows = {}
        for c in cs:
            ec = {}
            for name, _ in family[1:]:
                ec[name] = float(-COEF * ((vals[name]["I1"] - base["I1"])
                                          + c * (vals[name][t] - base[t]))
                                 + vals[name]["V4"] - base["V4"])
            for fam, pref in (("wide", "wide_"), ("saw2", "saw_lam_2_"),
                              ("saw1", "saw_lam_1_"), ("saw0.5", "saw_lam_0.5_")):
                keys = [k for k in ec if k.startswith(pref)]
                es = [ec[k] for k in keys]
                imin = int(np.argmin(es))
                rows.setdefault(f"c_{c:g}", {})[fam] = {
                    "E_corr": dict(zip(keys, es)),
                    "min": float(es[imin]), "argmin": keys[imin],
                    "interior": bool(imin < len(es) - 1),
                    "goes_negative": bool(min(es) < 0)}
        out["per_term"][t] = {"c_values": cs, "c_choice": how,
                              "probe": rows}
        log(f"BOUNDED {t}: c = {cs} ({how})")
    return out


# ---------------- LADDER (resolution control of the BOUNDED lattice) ----
def stage_ladder():
    """I1 only (the certified e_parts, cheap): E_corr of four family
    members on the h ladder 0.375 / 0.25 / 0.1875 at L = 24. Reads
    whether the n = 64 sawtooth numbers of BOUNDED are resolved."""
    out = {"L": 24.0, "rows": {}}
    for n in (64, 96, 128):
        cfg = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=n, L=24.0)
        X, Y, Z = B3.coords(n, cfg["h"])
        R = np.sqrt(X * X + Y * Y + Z * Z)
        nx, ny, nz = X / R, Y / R, Z / R
        K = np.zeros(X.shape + (4, 4))
        K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
        K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
        K2 = np.zeros_like(K)
        K2[..., 0, 0] = 1.0
        for i, a in enumerate((nx, ny, nz)):
            for j, bb in enumerate((nx, ny, nz)):
                K2[..., 1 + i, 1 + j] = a * bb
        Mb = B8.dressed(cfg, 0.0)
        e0 = sum(B3.e_parts(Mb, cfg))
        row = {}
        for name, bl in (("wide_amp_0.05", 0.05 * np.tanh(R / 2.0)),
                         ("saw_lam_2_A_0.05", C14.saw_of(0.05, 2.0, R)),
                         ("saw_lam_2_A_0.4", C14.saw_of(0.4, 2.0, R)),
                         ("saw_lam_1_A_0.4", C14.saw_of(0.4, 1.0, R))):
            Qb = (np.eye(4)[None, None, None]
                  + np.sinh(bl)[..., None, None] * K
                  + (np.cosh(bl) - 1.0)[..., None, None] * K2)
            Md = B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
            row[name] = float(sum(B3.e_parts(Md, cfg)) - e0)
        out["rows"][f"n{n}_h{cfg['h']:g}"] = row
        log(f"LADDER n {n}: {row}")
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        bnd = json.load(f)["BND"]
    out["quadrature_record_saw_lam_2"] = {
        f"A_{a:g}": e for a, e in zip(bnd[0]["amps"], bnd[0]["E"])}
    out["reading"] = ("the wide family and the small-amplitude sawtooth "
                      "converge to a NEGATIVE well (the certified eta "
                      "runaway); the large-amplitude sawtooth is NOT "
                      "resolved at h = 0.375 (sign flips along the "
                      "ladder toward the positive quadrature record), so "
                      "the BOUNDED large-A sawtooth entries are lattice "
                      "artifacts and only the wide + A <= 0.05 rows are "
                      "read")
    return out


# ---------------- INV ----------------
def stage_inv():
    out = {}
    # (a) the certified internal-conjugation INV stage, reproduced exactly
    rng = np.random.default_rng(21160)
    B16.stage_ident(rng)                       # consume the IDENT draws
    from scipy.ndimage import gaussian_filter
    cfg = B3.base_cfg(n=14, L=21.0)
    p = LAG.default_params()
    MB = np.stack([[gaussian_filter(rng.normal(size=(14,) * 3), 2.0)
                    for _ in range(4)] for _ in range(4)], axis=-1)
    MB = MB.reshape(14, 14, 14, 4, 4)
    MB = B3.vac4(cfg)[None, None, None] + 0.5 * B3.sym4(MB)
    gens = {}
    Grot = np.zeros((4, 4)); Grot[1, 2], Grot[2, 1] = -0.3, 0.3
    gens["so3_rot"] = Grot
    Gboo = np.zeros((4, 4)); Gboo[0, 1] = Gboo[1, 0] = 0.25
    gens["so13_boost"] = Gboo
    Gcmp = np.zeros((4, 4)); Gcmp[0, 1], Gcmp[1, 0] = -0.25, 0.25
    gens["so4_compact"] = Gcmp
    conj = {}
    for t in TERMS + ["I3_mixed_eta"]:
        E0 = static(t, MB, cfg, p)
        row = {"E0": E0}
        for nm, Gm in gens.items():
            L = expm(Gm)
            ML = np.einsum("ab,...bc,dc->...ad", L, MB, L)
            row[nm] = float(abs(static(t, ML, cfg, p) - E0)
                            / max(abs(E0), 1e-300))
        conj[t] = row
        log(f"INV conj {t}: " + ", ".join(f"{k} {row[k]:.2e}" for k in gens))
    with open(os.path.join(DATA, "m5_21_16_field.json")) as f:
        st = json.load(f)["INV"]
    out["internal_conjugation"] = {
        "drifts": conj,
        "control_I1_vs_stored_eta": {k: {"measured": conj["I1"][k],
                                        "stored": st["eta"][k]}
                                     for k in gens},
        "control_I1_frob_vs_stored_flip": {
            k: {"measured": conj["I1_frob"][k], "stored": st["flip"][k]}
            for k in gens},
        "reading": ("conjugation alone rotates the internal indices and "
                    "leaves the derivative index fixed: a symmetry of I1 "
                    "(derivative pair self-contracted) but NOT of the mixed "
                    "terms I2..I6, whose derivative-internal contractions "
                    "need the full Lorentz action (b)")}

    # (b) the full action on the LATTICE jets, and (c) the no-eta controls
    rng2 = np.random.default_rng(3203)
    p2 = LAG.default_params(s=S_MAIN, g=G_MAIN)
    cfg2 = B3.base_cfg(s=S_MAIN, g=G_MAIN, n=14, L=21.0)
    MB2 = np.stack([[gaussian_filter(rng2.normal(size=(14,) * 3), 2.0)
                     for _ in range(4)] for _ in range(4)], axis=-1)
    MB2 = MB2.reshape(14, 14, 14, 4, 4)
    MB2 = B3.vac4(cfg2)[None, None, None] + 0.5 * B3.sym4(MB2)
    jets = LAG.lattice_jets(MB2, cfg2)
    K_I1_noeta = LAG._K_from_pattern("mnab", "mnab", 0.5, deriv_metric="delta")
    K_I2 = REG["I2"]._K()

    def dens_sum(t, A, M):
        if t == "I1_deriv_delta":
            return float(np.sum(LAG.density_from_K(LAG.F_of_A(A), K_I1_noeta)))
        if t == "I2_plain_bracket":
            P = np.einsum("m...ab,n...bc->...mnac", A, A, optimize=True)
            F = P - P.swapaxes(-4, -3)
            return float(np.sum(LAG.density_from_K(F, K_I2)))
        return float(np.sum(REG[t].density(A, M, p2)))

    names = TERMS + ["I3_mixed_eta", "I1_deriv_delta", "I2_plain_bracket"]
    trans = {}
    for kind in ("boost", "rotation"):
        for k in range(2):
            L = LAG._lorentz(rng2, kind)
            Linv_T = np.linalg.inv(L).T
            key = f"{kind}_{k}"
            for t in names:
                e0, e1 = 0.0, 0.0
                for A, wt in jets:
                    Mp = np.einsum("ab,...bc,dc->...ad", L, MB2, L)
                    Ap = np.einsum("mn,n...ab->m...ab", Linv_T,
                                   np.einsum("ab,n...bc,dc->n...ad", L, A, L))
                    e0 += wt * dens_sum(t, A, MB2)
                    e1 += wt * dens_sum(t, Ap, Mp)
                trans.setdefault(t, {})[key] = float(
                    abs(e1 - e0) / max(abs(e0), 1e-300))
            log(f"INV full-action {key} done")
    cov_terms = ["I1", "I2", "I3", "I4", "I5", "I6"]
    ctrl_terms = ["I1_frob", "I3_mixed_eta", "I1_deriv_delta",
                  "I2_plain_bracket"]
    worst_cov = max(max(trans[t].values()) for t in cov_terms)
    min_ctrl_boost = min(min(v for k, v in trans[t].items()
                             if k.startswith("boost")) for t in ctrl_terms)
    out["full_action_lattice"] = {
        "drifts": trans, "worst_covariant": float(worst_cov),
        "min_control_boost_drift": float(min_ctrl_boost),
        "pass_covariant_1e-10": bool(worst_cov <= 1e-10),
        "controls_fail_under_boost": bool(min_ctrl_boost > 1e-3),
        "field": "random smooth 4x4 field, n = 14, L = 21, rng 3203, "
                 "vacuum + 0.5 sym4(gaussian_filter(normal, 2))",
        "transform": "A'_mu = (Lambda^{-T})_mu^nu Lambda A_nu Lambda^T, "
                     "M' = Lambda M Lambda^T; jets = the certified sym "
                     "stencil branches; static field (A_0 = 0 before the "
                     "boost, nonzero after)"}
    return out


# ---------------- plot ----------------
def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
            "#008300", "#4a3aa7"]
    chan = out["CHAN"]
    chans = chan["channels"]
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    nt = len(TERMS)
    wdt = 0.8 / nt
    a = ax[0]
    for i, t in enumerate(TERMS):
        ys = [chan["C_E_table"][t][ch] for ch in chans]
        a.bar(np.arange(len(chans)) + (i - nt / 2 + 0.5) * wdt, ys, wdt,
              color=cols[i], label=t, linewidth=0)
    a.axhline(0, color="k", lw=0.6)
    a.set_xticks(np.arange(len(chans)))
    a.set_xticklabels(chans, rotation=20)
    a.set_yscale("symlog", linthresh=1.0)
    a.set_ylabel("C_E = -4 C (energy kinetic coefficient; >= 0 stable)")
    a.set_title("CHAN: omega^2 coefficient per channel per term\n"
                "(analytic boost-hedgehog family, g = 32, n = 32, L = 48)")
    a.legend(frameon=False, ncol=2, fontsize=8)
    a.grid(axis="y", alpha=0.25)
    a = ax[1]
    dr = out["DRESSED"]
    keys = list(dr["scales"].keys())
    for i, t in enumerate(TERMS):
        ys = [dr["scales"][k][t]["C_E"] for k in keys]
        a.bar(np.arange(len(keys)) + (i - nt / 2 + 0.5) * wdt, ys, wdt,
              color=cols[i], label=t, linewidth=0)
    a.axhline(0, color="k", lw=0.6)
    a.set_xticks(np.arange(len(keys)))
    a.set_xticklabels(["b = 0", "b = b*", "b = 2 b*"])
    a.set_yscale("symlog", linthresh=10.0)
    a.set_ylabel("C_E = -4 C on the dressed electron")
    a.set_title("DRESSED: omega^2 coefficient on the M5.21.14 boost-dressed\n"
                "electron (clock a0 = the family flow), per term")
    a.grid(axis="y", alpha=0.25)
    fig.suptitle("M5.32 R1 (c): lattice time-sector screen of the term basis")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "m5_32_r1_screen_chan.png"), dpi=140)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all")
    args = ap.parse_args()
    path = os.path.join(DATA, "m5_32_r1_screen.json")
    out = {"task": "M5.32 R1 arm (c): lattice time-sector screen",
           "terms": TERMS, "normalization": "L' = -4 (I1 + c T) - V4; "
           "C_E = -4 C_T; kin(c) = C_E(I1) + c C_E(T)",
           "python": sys.version.split()[0], "numpy": np.__version__}
    if os.path.exists(path) and args.stage != "all":
        with open(path) as f:
            out.update(json.load(f))
    stages = (["CHAN", "DRESSED", "KERNEL", "BOUNDED", "LADDER", "INV"]
              if args.stage == "all" else [args.stage])
    timings = out.get("runtime_s", {})
    for st in stages:
        t0 = time.time()
        if st == "CHAN":
            out["CHAN"] = stage_chan()
        elif st == "DRESSED":
            out["DRESSED"] = stage_dressed()
        elif st == "KERNEL":
            out["KERNEL"] = stage_kernel()
        elif st == "BOUNDED":
            out["BOUNDED"] = stage_bounded(out["CHAN"])
        elif st == "LADDER":
            out["LADDER"] = stage_ladder()
        elif st == "INV":
            out["INV"] = stage_inv()
        timings[st] = round(time.time() - t0, 1)
        out["runtime_s"] = timings
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        log(f"stage {st} written ({timings[st]} s)")
    if "CHAN" in out and "DRESSED" in out:
        plot(out)
    log("done")


if __name__ == "__main__":
    main()
