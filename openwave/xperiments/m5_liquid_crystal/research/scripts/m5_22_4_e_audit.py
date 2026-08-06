"""M5.22.4 independent adversarial audit (auditor-owned implementations).

Attacks the five claims of findings/m5_22_4_note.md with independent
code. The only borrowed machinery is PAIR.orient_v1 (the audited
eigenvector-orientation continuity fix), used for orientation ONLY;
every derivative, Gamma, curvature, flux, kin, and energy computation
in this file is auditor-owned and written from the stated definitions,
not copied from the scripts under test.

C1  algebraic identity (Gvec_i x Gvec_j)_a = e_a.(d_i e_a x d_j e_a)
    for a = 1, 2, 3, tested on random smooth synthetic rotation-field
    frames (Rodrigues), two resolutions, convergence-order read.
C2  full-F comp3 flux + basic mermin flux at half 18 on the four
    target states, own integrator, vs m5_22_4_fullf_all.json.
C3  kin(M; a0) for clock_local and boost_z on the prot and deut P1
    endpoints, own a0 build + own sym-stencil kin, vs the p2 rows.
C4  ladder arithmetic from the row JSONs (every rung above ctrl,
    offset = omega^2 kin at 0.8, static parts vs ctrl) + own static
    functional (u + trace-target V4) on the prot omega=0.8 endpoint.
C5  far-field oriented Mermin-Ho charge of the omega=0.8 endpoints
    (prot, deut) vs the P1 charge, own integrator.

Run: python3 m5_22_4_e_audit.py     (~2 min)
Out: ../data/m5_22_4_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


PAIR = _load("pair_audit", "m5_21_4_a_pair.py")   # orient_v1 ONLY

# stack constants, restated from the definitions under audit
W1 = 0.000724023879
SG = 8.0
DELTA = 0.3
H = 1.5
N = 32
RENV = 18.0
ETA_SIGNS = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(ETA_SIGNS)


# ================= auditor-owned numerics =================
def cdiff(f, ax, h):
    """central difference, boundary planes zeroed (interior use only)."""
    out = np.zeros_like(f)
    sl_p = [slice(None)] * f.ndim
    sl_m = [slice(None)] * f.ndim
    sl_c = [slice(None)] * f.ndim
    sl_p[ax], sl_m[ax], sl_c[ax] = \
        slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(sl_c)] = (f[tuple(sl_p)] - f[tuple(sl_m)]) / (2.0 * h)
    return out


def one_sided(f, ax, h, kind):
    """fwd/bwd one-sided difference, unfilled plane left zero
    (the m5_21_3_a_4d branch semantics, restated)."""
    out = np.zeros_like(f)
    lo = [slice(None)] * f.ndim
    hi = [slice(None)] * f.ndim
    lo[ax], hi[ax] = slice(0, -1), slice(1, None)
    d = (f[tuple(hi)] - f[tuple(lo)]) / h
    tgt = lo if kind == "fwd" else hi
    out[tuple(tgt)] = d
    return out


def box_flux(B, i0, i1, z0, z1, h):
    """lattice box flux: sum of the outward normal component over the
    six faces (each face spans the full inclusive transverse range),
    times h^2 / 4pi."""
    s = B[i1, i0:i1 + 1, z0:z1 + 1, 0].sum() \
        - B[i0, i0:i1 + 1, z0:z1 + 1, 0].sum()
    s += B[i0:i1 + 1, i1, z0:z1 + 1, 1].sum() \
        - B[i0:i1 + 1, i0, z0:z1 + 1, 1].sum()
    s += B[i0:i1 + 1, i0:i1 + 1, z1, 2].sum() \
        - B[i0:i1 + 1, i0:i1 + 1, z0, 2].sum()
    return float(s * h * h / (4.0 * np.pi))


def flux_half(B, half, h=H, n=N, placement="authors"):
    """faces for the physical half-width: authors snap the center to
    round((n-1)/2) = 16 and go +-k; the symmetric alternative uses
    (4, 27) style faces. Both are closed boxes; plateau => agreement."""
    k = int(round(half / h))
    c = int(round((n - 1) / 2.0))
    if placement == "authors":
        i0, i1 = c - k, c + k
    else:
        i0, i1 = c - k, c + k - 1
    return box_flux(B, i0, i1, i0, i1, h)


def mermin_density(nhat, h):
    """B = (F23, -F13, F12) with F_ij = n.(d_i n x d_j n), own diffs."""
    dn = [cdiff(nhat, ax, h) for ax in range(3)]
    f12 = np.einsum("...a,...a->...", nhat, np.cross(dn[0], dn[1]))
    f13 = np.einsum("...a,...a->...", nhat, np.cross(dn[0], dn[2]))
    f23 = np.einsum("...a,...a->...", nhat, np.cross(dn[1], dn[2]))
    return np.stack([f23, -f13, f12], axis=-1)


def inner_eta(F, G):
    """<F,G>_eta = tr(eta F eta G^T) = sum_ab s_a s_b F_ab G_ab."""
    S = ETA_SIGNS[:, None] * ETA_SIGNS[None, :]
    return np.einsum("...ab,...ab->...", S * F, G)


def comm_eta(A, B):
    return A @ ETA @ B - B @ ETA @ A


def coords3(n=N, h=H):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def interior_mask(n=N, margin=2):
    m = np.zeros((n, n, n), dtype=bool)
    m[margin:n - margin, margin:n - margin, margin:n - margin] = True
    return m


# ================= C1: the algebraic identity =================
def rodrigues(theta):
    """R = exp(hat(theta)) per voxel, theta shape (..., 3)."""
    ang = np.linalg.norm(theta, axis=-1)
    ang_s = np.where(ang < 1e-12, 1e-12, ang)
    k = theta / ang_s[..., None]
    K = np.zeros(theta.shape[:-1] + (3, 3))
    K[..., 0, 1], K[..., 0, 2] = -k[..., 2], k[..., 1]
    K[..., 1, 0], K[..., 1, 2] = k[..., 2], -k[..., 0]
    K[..., 2, 0], K[..., 2, 1] = -k[..., 1], k[..., 0]
    I = np.broadcast_to(np.eye(3), K.shape)
    s = np.sin(ang)[..., None, None]
    c = np.cos(ang)[..., None, None]
    R = I + s * K + (1.0 - c) * (K @ K)
    R = np.where((ang < 1e-12)[..., None, None], I, R)
    return R


def identity_errors(O, h):
    """max |(Gvec_i x Gvec_j)_a - e_a.(d_i e_a x d_j e_a)| over the
    interior, per internal component a, joined over the 3 (i,j) pairs;
    plus the RHS scale for the relative read."""
    n = O.shape[0]
    G = []
    for ax in range(3):
        dO = cdiff(O, ax, h)
        Gm = np.einsum("...ka,...kb->...ab", O, dO)
        G.append(np.stack([Gm[..., 2, 1], Gm[..., 0, 2],
                           Gm[..., 1, 0]], axis=-1))
    e = [O[..., :, a] for a in range(3)]
    de = [[cdiff(e[a], ax, h) for ax in range(3)] for a in range(3)]
    itr = interior_mask(n, margin=2)
    errs, scales = [], []
    for a in range(3):
        e_max, s_max = 0.0, 0.0
        for (i, j) in ((0, 1), (0, 2), (1, 2)):
            lhs = np.cross(G[i], G[j])[..., a]
            rhs = np.einsum("...a,...a->...", e[a],
                            np.cross(de[a][i], de[a][j]))
            e_max = max(e_max, float(np.abs((lhs - rhs)[itr]).max()))
            s_max = max(s_max, float(np.abs(rhs[itr]).max()))
        errs.append(e_max)
        scales.append(s_max)
    return errs, scales


def c1():
    rng = np.random.default_rng(2264)
    out = {"trials": []}
    worst_rel, worst_ratio_lo, worst_ratio_hi = 0.0, np.inf, 0.0
    for trial in range(3):
        n1 = 25
        n2 = 2 * n1 - 1                      # same box, h -> h/2
        # one analytic gaussian-sum parameter draw, evaluated on BOTH
        # grids (the field is the same smooth function of x)
        params = [(rng.uniform(-0.6, 0.6, size=3),
                   rng.uniform(0.35, 0.7),
                   rng.uniform(-1.2, 1.2, size=3)) for _ in range(5)]

        def theta_on(n):
            x = np.linspace(-1.0, 1.0, n)
            X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
            th = np.zeros((n, n, n, 3))
            for c, s, a in params:
                g = np.exp(-((X - c[0]) ** 2 + (Y - c[1]) ** 2
                             + (Z - c[2]) ** 2) / (2 * s * s))
                th += g[..., None] * a
            return th
        row = {}
        for lab, n in (("coarse", n1), ("fine", n2)):
            h = 2.0 / (n - 1)
            O = rodrigues(theta_on(n))
            errs, scales = identity_errors(O, h)
            row[lab] = {"h": h, "err_max_per_comp": errs,
                        "rhs_scale_per_comp": scales}
        ratios = [row["coarse"]["err_max_per_comp"][a]
                  / max(row["fine"]["err_max_per_comp"][a], 1e-300)
                  for a in range(3)]
        rels = [row["fine"]["err_max_per_comp"][a]
                / max(row["fine"]["rhs_scale_per_comp"][a], 1e-300)
                for a in range(3)]
        row["ratio_coarse_over_fine"] = ratios
        row["fine_rel_err"] = rels
        out["trials"].append(row)
        worst_rel = max(worst_rel, max(rels))
        worst_ratio_lo = min(worst_ratio_lo, min(ratios))
        worst_ratio_hi = max(worst_ratio_hi, max(ratios))
    out["worst_fine_rel_err"] = worst_rel
    out["ratio_range"] = [worst_ratio_lo, worst_ratio_hi]
    # second-order stencils on a product of first derivatives: expect
    # error ~ h^2, i.e. ratio ~ 4; accept [2.5, 6]
    out["verdict"] = "PASS" if (worst_rel < 5e-2
                                and worst_ratio_lo > 2.5
                                and worst_ratio_hi < 6.0) else "REFUTED"
    return out


# ================= shared frame build (C2) =================
def short_axis_oriented(M):
    """own eigenvalue remap so the SHORT axis leads, then the allowed
    orient_v1 continuity machinery on the remapped tensor."""
    lam, vec = np.linalg.eigh(M)
    w = np.zeros_like(lam)
    w[..., 0] = 2.0                          # short axis -> leading
    w[..., 1] = 0.0
    w[..., 2] = 1.0
    M2 = np.einsum("...ai,...i,...bi->...ab", vec, w, vec)
    return PAIR.orient_v1(M2)


def full_frame_own(M):
    e3, ncf3 = PAIR.orient_v1(M)
    e1, ncf1 = short_axis_oriented(M)
    e1 = e1 - np.einsum("...a,...a->...", e1, e3)[..., None] * e3
    e1 = e1 / np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(e3, e1)
    return np.stack([e1, e2, e3], axis=-1), int(ncf3), int(ncf1)


def comp3_field_own(M, h):
    """the long-axis internal component of the full-F curvature,
    assembled (R23_3, -R13_3, R12_3), all own numerics."""
    O, ncf3, ncf1 = full_frame_own(M)
    G = []
    for ax in range(3):
        dO = cdiff(O, ax, h)
        Gm = np.einsum("...ka,...kb->...ab", O, dO)
        G.append(np.stack([Gm[..., 2, 1], Gm[..., 0, 2],
                           Gm[..., 1, 0]], axis=-1))
    r23 = np.cross(G[1], G[2])[..., 2]
    r13 = np.cross(G[0], G[2])[..., 2]
    r12 = np.cross(G[0], G[1])[..., 2]
    return np.stack([r23, -r13, r12], axis=-1), ncf3, ncf1


C2_TARGETS = [
    ("m5_22", "P-0.5_plane_sc6_n32_pinned_d0.3"),
    ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3"),
    ("m5_22_1", "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3"),
    ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3"),
]


def c2():
    with open(os.path.join(DATA, "m5_22_4_fullf_all.json")) as f:
        ref = json.load(f)
    out = {"states": {}}
    ok = True
    worst_dq, worst_gap = 0.0, 0.0
    for src, tag in C2_TARGETS:
        z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
        M = z["M"].astype(np.float64)
        e3, _ = PAIR.orient_v1(M)
        basic = mermin_density(e3, H)
        comp3, ncf3, ncf1 = comp3_field_own(M, H)
        row = {
            "basic_half18": flux_half(basic, 18.0),
            "comp3_half18": flux_half(comp3, 18.0),
            "basic_half18_symbox": flux_half(basic, 18.0,
                                             placement="sym"),
            "comp3_half18_symbox": flux_half(comp3, 18.0,
                                             placement="sym"),
            "orient_conflicts": [ncf3, ncf1],
        }
        itr = interior_mask()
        d = (comp3 - basic)[itr]
        row["identity_diff_max"] = float(np.abs(d).max())
        row["identity_diff_rms"] = float(np.sqrt((d * d).mean()))
        row["ref_comp3_half18"] = ref[tag]["comp3"]["half18"]
        row["ref_basic_half18"] = ref[tag]["basic"]["half18"]
        row["d_comp3_vs_ref"] = abs(row["comp3_half18"]
                                    - row["ref_comp3_half18"])
        row["d_basic_vs_ref"] = abs(row["basic_half18"]
                                    - row["ref_basic_half18"])
        row["comp3_minus_basic"] = abs(row["comp3_half18"]
                                       - row["basic_half18"])
        worst_dq = max(worst_dq, row["d_comp3_vs_ref"],
                       row["d_basic_vs_ref"])
        worst_gap = max(worst_gap, row["comp3_minus_basic"])
        if row["d_comp3_vs_ref"] > 1e-3 or row["d_basic_vs_ref"] > 1e-3:
            ok = False
        out["states"][tag] = row
    out["worst_abs_dev_vs_ref"] = worst_dq
    out["worst_comp3_minus_basic"] = worst_gap
    out["claim_basic_matches_comp3_6e-4"] = bool(worst_gap <= 6.5e-4)
    out["verdict"] = "PASS" if (ok and worst_gap <= 1e-3) else "REFUTED"
    return out


# ================= C3: kin, own build =================
def envelope_own(n=N, h=H, renv=RENV):
    X, Y, Z = coords3(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return np.exp(-((r / renv) ** 4))


def a0_clock_local(M4):
    """rotation about the local leading spatial eigenvector, envelope
    weighted, unit global Frobenius norm."""
    lam, V = np.linalg.eigh(M4[..., 1:4, 1:4])
    nh = V[..., :, 2]
    W = np.zeros(M4.shape)
    n1, n2, n3 = nh[..., 0], nh[..., 1], nh[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    return _finish_a0(W, M4)


def a0_boost_z(M4):
    K = np.zeros((4, 4))
    K[0, 3] = K[3, 0] = 1.0
    W = np.broadcast_to(K, M4.shape)
    return _finish_a0(W, M4)


def _finish_a0(G, M4):
    w = envelope_own()[..., None, None]
    a0 = w * (G @ M4 - M4 @ G.swapaxes(-1, -2))
    nrm = np.sqrt(np.sum(a0 * a0))
    return a0 / max(nrm, 1e-300)


def kin_own(M4, a0, h=H):
    """kin = 4 sum_i <comm_eta(a0, A_i), comm_eta(a0, A_i)>_eta h^3,
    sym stencil = mean of fwd/bwd one-sided branches."""
    k = 0.0
    for kind in ("fwd", "bwd"):
        for ax in range(3):
            A = one_sided(M4, ax, h, kind)
            F = comm_eta(a0, A)
            k += 0.5 * 4.0 * float(np.sum(inner_eta(F, F)))
    return k * h ** 3


def c3():
    out = {"states": {}}
    ok = True
    for key in ("prot", "deut"):
        M4 = np.load(os.path.join(
            DATA, f"m5_22_4_p1_{key}.npz"))["M"].astype(np.float64)
        with open(os.path.join(DATA,
                               f"m5_22_4_row_p2_{key}.json")) as f:
            ref = json.load(f)["rows"]
        row = {}
        for nm, build in (("clock_local", a0_clock_local),
                          ("boost_z", a0_boost_z)):
            k = kin_own(M4, build(M4))
            kr = ref[nm]["kin"]
            row[nm] = {"kin_own": k, "kin_ref": kr,
                       "rel_dev": abs(k - kr) / max(abs(kr), 1e-300),
                       "sign_match": bool(np.sign(k) == np.sign(kr))}
            if not row[nm]["sign_match"] or row[nm]["rel_dev"] > 1e-3:
                ok = False
        # the note's sign claim: rotations > 0, boosts < 0
        if not (row["clock_local"]["kin_own"] > 0
                and row["boost_z"]["kin_own"] < 0):
            ok = False
        out["states"][key] = row
    out["verdict"] = "PASS" if ok else "REFUTED"
    return out


# ================= C4: ladder arithmetic + own static functional ====
def e_static_own(M4, h=H, sg=SG, delta=DELTA):
    """(E_u, E_v): u = 4 sum_{i<j} <comm_eta(A_i, A_j)>^2_eta, sym
    stencil; V4 = W1 sum_{p=1..4} (tr((M eta)^p) - C_p)^2; h^3."""
    e_u = 0.0
    for kind in ("fwd", "bwd"):
        A = [one_sided(M4, ax, h, kind) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                e_u += 0.5 * 4.0 * float(np.sum(inner_eta(F, F)))
    Me = M4 @ ETA
    P = Me.copy()
    vd = 0.0
    for p in range(1, 5):
        if p > 1:
            P = P @ Me
        tr = np.einsum("...kk->...", P)
        cp = sg ** p + 1.0 + delta ** p
        vd = vd + (tr - cp) ** 2
    return e_u * h ** 3, float(np.sum(vd)) * W1 * h ** 3


def c4():
    keys = ("prot", "neut", "d2", "deut")
    out = {"states": {}}
    ok = True
    for key in keys:
        with open(os.path.join(
                DATA, f"m5_22_4_row_p3_{key}_clock_local.json")) as f:
            p3 = json.load(f)
        with open(os.path.join(DATA,
                               f"m5_22_4_row_ctrl_{key}.json")) as f:
            ct = json.load(f)
        ctrl_e = ct["E_end"]
        rungs = [r for r in p3["ladder"] if r["omega"] > 0]
        row = {"ctrl_E": ctrl_e,
               "rungs": [], "all_rungs_above_ctrl": True,
               "max_E_decomp_err": 0.0}
        for r in rungs:
            dec = abs(r["E"] - (r["E_u"] + r["E_v"]
                                + r["omega"] ** 2 * r["kin"]))
            row["max_E_decomp_err"] = max(row["max_E_decomp_err"], dec)
            above = r["E"] > ctrl_e
            if not above:
                row["all_rungs_above_ctrl"] = False
            row["rungs"].append({"omega": r["omega"],
                                 "E_minus_ctrl": r["E"] - ctrl_e,
                                 "above": above})
        r8 = rungs[-1]
        assert abs(r8["omega"] - 0.8) < 1e-12
        stat8 = r8["E_u"] + r8["E_v"]
        row["static_at_0.8"] = stat8
        row["static_minus_ctrl"] = stat8 - ctrl_e
        row["omega2kin_at_0.8"] = 0.8 ** 2 * r8["kin"]
        row["offset_minus_omega2kin"] = (r8["E"] - ctrl_e) \
            - row["omega2kin_at_0.8"]
        if not row["all_rungs_above_ctrl"]:
            ok = False
        if abs(row["static_minus_ctrl"]) > 2e-3:
            ok = False
        if row["max_E_decomp_err"] > 1e-9:
            ok = False
        out["states"][key] = row
    # own recompute of the prot omega=0.8 endpoint static parts
    M4 = np.load(os.path.join(
        DATA, "m5_22_4_p3_prot_clock_local.npz"))["M"].astype(np.float64)
    eu, ev = e_static_own(M4)
    with open(os.path.join(
            DATA, "m5_22_4_row_p3_prot_clock_local.json")) as f:
        r8 = json.load(f)["ladder"][-1]
    rec = {"E_u_own": eu, "E_v_own": ev,
           "E_u_ref": r8["E_u"], "E_v_ref": r8["E_v"],
           "sum_own": eu + ev, "sum_ref": r8["E_u"] + r8["E_v"],
           "d_u": abs(eu - r8["E_u"]), "d_v": abs(ev - r8["E_v"]),
           "note": "npz stores float32; sub-1e-3 deviations expected"}
    out["prot_p3_recompute"] = rec
    if abs(eu + ev - rec["sum_ref"]) / rec["sum_ref"] > 1e-3:
        ok = False
    out["verdict"] = "PASS" if ok else "REFUTED"
    return out


# ================= C5: charge preservation =================
def c5():
    out = {"states": {}}
    ok = True
    for key, p1_q in (("prot", None), ("deut", None)):
        q = {}
        for lab, fn in (
                ("p1", f"m5_22_4_p1_{key}.npz"),
                ("p3_w0.8", f"m5_22_4_p3_{key}_clock_local.npz")):
            M4 = np.load(os.path.join(DATA, fn))["M"].astype(np.float64)
            M3 = M4[..., 1:4, 1:4]
            nhat, ncf = PAIR.orient_v1(M3)
            B = mermin_density(nhat, H)
            f165 = flux_half(B, 16.5)
            f18 = flux_half(B, 18.0)
            q[lab] = {"flux_half16.5": f165, "flux_half18": f18,
                      "q_far": 0.5 * (f165 + f18),
                      "n_conflicts": int(ncf)}
        with open(os.path.join(DATA,
                               f"m5_22_4_row_p1_{key}.json")) as f:
            qref = json.load(f)["ring_end"]["q_far"]
        row = {"own": q, "q_far_p1_ref": qref,
               "abs_q_change_own": abs(abs(q["p3_w0.8"]["q_far"])
                                       - abs(q["p1"]["q_far"])),
               "d_p1_own_vs_ref": abs(q["p1"]["q_far"] - qref)}
        if row["abs_q_change_own"] > 5e-4 or \
                row["d_p1_own_vs_ref"] > 1e-3:
            ok = False
        out["states"][key] = row
    out["verdict"] = "PASS" if ok else "REFUTED"
    return out


# ================= main =================
def main():
    res = {}
    for nm, fn in (("C1_identity", c1), ("C2_flux_table", c2),
                   ("C3_kin_signs", c3), ("C4_ladder", c4),
                   ("C5_charge", c5)):
        print(f"== {nm} ==", flush=True)
        res[nm] = fn()
        print(json.dumps({"verdict": res[nm]["verdict"]}), flush=True)
    with open(os.path.join(DATA, "m5_22_4_audit.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: v["verdict"] for k, v in res.items()},
                     indent=1))
    return res


if __name__ == "__main__":
    main()
