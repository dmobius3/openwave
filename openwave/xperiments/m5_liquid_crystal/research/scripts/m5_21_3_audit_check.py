"""M5.21.3 INDEPENDENT ADVERSARIAL AUDIT (own re-implementation).

Tries to REFUTE claims C1-C7 of findings/m5_21_3_note.md with code
written from scratch: own difference stencils, own eta inner product in
diagonal form (sum_{a,c} eta_a eta_c F_ac G_ac instead of the task's
tr(eta F eta G^T) matmul route), own V4 both via an eigenvalue route
and via an own power route, own Taylor expm for the SO(1,3) check,
fresh RNG seeds for the curvature directions.

The ONLY task-code reuse is the authorized gates re-run (C1a), executed
with the task module's DATA redirected to a temp dir so no repo file is
touched. Output: ../data/m5_21_3_audit.json (new file) + stdout table.

Run: python3 m5_21_3_audit_check.py
"""
from __future__ import annotations

import importlib
import json
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(ETA_D)
W1 = 0.000724023879
RENV = 10.0
G_COUPL = 8.0


# ---------------- my own numerics (independent implementations) -----
def a_d1(f, ax, h, st):
    """First differences, own slicing code. fwd leaves the last plane
    zero, bwd the first plane zero, cen = central with one-sided edges
    (matches the conventions under audit)."""
    out = np.zeros_like(f)
    lo = [slice(None)] * f.ndim
    lo[ax] = slice(0, -1)
    lo = tuple(lo)
    hi = [slice(None)] * f.ndim
    hi[ax] = slice(1, None)
    hi = tuple(hi)
    if st == "fwd":
        out[lo] = (f[hi] - f[lo]) / h
    elif st == "bwd":
        out[hi] = (f[hi] - f[lo]) / h
    else:
        mid = [slice(None)] * f.ndim
        mid[ax] = slice(1, -1)
        lo2 = [slice(None)] * f.ndim
        lo2[ax] = slice(0, -2)
        hi2 = [slice(None)] * f.ndim
        hi2[ax] = slice(2, None)
        out[tuple(mid)] = (f[tuple(hi2)] - f[tuple(lo2)]) / (2 * h)
        for edge, nbr in ((0, 1), (-1, -2)):
            e = [slice(None)] * f.ndim
            e[ax] = edge
            nb = [slice(None)] * f.ndim
            nb[ax] = nbr
            out[tuple(e)] = (f[tuple(e)] - f[tuple(nb)]) / h \
                if edge == -1 else (f[tuple(nb)] - f[tuple(e)]) / h
    return out


def m_eta(A, B):
    """A eta B using the diagonal of eta directly (own route)."""
    return np.einsum("...ac,c,...cb->...ab", A, ETA_D, B, optimize=True)


def comm_eta(A, B):
    return m_eta(A, B) - m_eta(B, A)


def inner_sum(F, G):
    """SUM over cells of tr(eta F eta G^T), diagonal form:
    sum_{a,c} eta_a eta_c F_ac G_ac (no matmuls)."""
    return float(np.einsum("...ac,...ac,a,c->", F, G, ETA_D, ETA_D,
                           optimize=True))


def e_u(M, h, st="sym", az_extra=None):
    """h^3 * 4 * sum_{cells, i<j} <F_ij, F_ij>_eta; sym = mean of the
    fwd and bwd builds. az_extra = the twist channel A_z shift."""
    sts = [("fwd", 0.5), ("bwd", 0.5)] if st == "sym" else [(st, 1.0)]
    tot = 0.0
    for br, wt in sts:
        A = [a_d1(M, ax, h, br) for ax in range(3)]
        if az_extra is not None:
            A[2] = A[2] + az_extra
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(A[i], A[j])
                tot += wt * 4.0 * inner_sum(F, F)
    return h ** 3 * tot


def v4(M, h, sg, delta, method="pow"):
    """W1 * sum_p (tr((M eta)^p) - C_p)^2, h^3-weighted. method='eig'
    computes the traces from the eigenvalues of M eta (a structurally
    different route); method='pow' via own einsum powers."""
    Me = M * ETA_D                        # (M eta)_ab = M_ab eta_b
    if method == "eig":
        lam = np.linalg.eigvals(Me)
        tr = [np.sum(lam ** p, axis=-1).real for p in range(1, 5)]
    else:
        P = Me
        tr = [np.einsum("...aa->...", P)]
        for _ in range(3):
            P = np.einsum("...ab,...bc->...ac", P, Me, optimize=True)
            tr.append(np.einsum("...aa->...", P))
    tot = 0.0
    for p in range(1, 5):
        cp = sg ** p + 1.0 + delta ** p
        tot += np.sum((tr[p - 1] - cp) ** 2)
    return float(h ** 3 * W1 * tot)


def e_tot(M, h, sg, delta, method="pow"):
    return e_u(M, h) + v4(M, h, sg, delta, method)


def kin(M, a0, h, st="sym"):
    """kin = h^3 * 4 * sum_i <comm_eta(a0, A_i), comm_eta(a0, A_i)>."""
    sts = [("fwd", 0.5), ("bwd", 0.5)] if st == "sym" else [(st, 1.0)]
    tot = 0.0
    for br, wt in sts:
        for ax in range(3):
            F = comm_eta(a0, a_d1(M, ax, h, br))
            tot += wt * 4.0 * inner_sum(F, F)
    return h ** 3 * tot


def envelope(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return np.exp(-((r / RENV) ** 4))


def a0_global(M, W, w, sign=-1.0):
    """a0 = w * (W M + sign * M W^T), unit Frobenius norm over the
    grid. sign=-1 is the task convention (G M - M G^T); sign=+1 is the
    tangent of the conjugation orbit M -> L M L^T used by gate G1 and
    the films (d/dt = G M + M G^T)."""
    raw = w[..., None, None] * (W @ M + sign * (M @ W.T))
    return raw / np.sqrt(np.sum(raw * raw)), raw


def a0_clock(M, w):
    """clock_local, task convention: rotation generator about the local
    leading eigenvector of the spatial block, enveloped, normalized."""
    _, V = np.linalg.eigh(M[..., 1:4, 1:4])
    n1, n2, n3 = V[..., 0, 2], V[..., 1, 2], V[..., 2, 2]
    W = np.zeros(M.shape)
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    raw = w[..., None, None] * (
        np.einsum("...ab,...bc->...ac", W, M, optimize=True)
        - np.einsum("...ab,...cb->...ac", M, W, optimize=True))
    return raw / np.sqrt(np.sum(raw * raw))


def my_expm(A, terms=48):
    out = np.eye(A.shape[0])
    term = np.eye(A.shape[0])
    for k in range(1, terms):
        term = term @ A / k
        out = out + term
    return out


def load_npz(name):
    Z = np.load(os.path.join(DATA, name))
    M = Z["M"].astype(np.float64)
    return M, float(Z["h"]), float(Z["s"]), float(Z["delta"])


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


# ---------------- C1: gates ----------------
def audit_c1():
    out = {}
    with open(os.path.join(DATA, "m5_21_3_gates.json")) as f:
        rec = json.load(f)
    sys.path.insert(0, HERE)
    task = importlib.import_module("m5_21_3_a_4d")
    tmp = tempfile.mkdtemp(prefix="m5_21_3_audit_gates_")
    task.DATA = tmp                      # no repo file gets touched
    rerun = task.gates()
    # roundoff-scale diagnostics (g0*, g1_so13) are not expected to
    # reproduce bit-exactly across BLAS paths; the audit criterion is
    # that the RERUN passes the same recorded bars
    bars = {"g0_static": ("<", 5e-9), "g0_kin": ("<", 5e-9),
            "g1_so13": ("<", 1e-9), "g1_negctrl": (">", 1e-6),
            "g3_regression": ("<", 1e-12)}
    cmp = {}
    ok = True
    for k, (op, bar) in bars.items():
        rv = float(rerun[k])
        same = bool((rv < bar and rec[k] < bar) if op == "<"
                    else (rv > bar and rec[k] > bar))
        cmp[k] = {"recorded": rec[k], "rerun": rv, "bar": f"{op} {bar}",
                  "both_pass_bar": same}
        ok = ok and same
    for k in ("g2_vac_s+1", "g2_vac_s-1"):
        rv = [float(x) for x in rerun[k]]
        same = bool(all(abs(x) < 1e-16 for x in rv)
                    and all(abs(x) < 1e-16 for x in rec[k]))
        cmp[k] = {"recorded": rec[k], "rerun": rv, "bar": "< 1e-16",
                  "both_pass_bar": same}
        ok = ok and same
    ok = ok and bool(rerun["pass"]) and bool(rec["pass"])
    out["gates_rerun_matches_recorded"] = ok
    out["gates_compare"] = cmp

    # my own SO(1,3) check: own field, own generator, own expm, own E
    rng = np.random.default_rng(20260718)
    n6, h6 = 6, 1.5
    sg, delta = 8.0, 0.3
    Mr = rng.normal(size=(n6, n6, n6, 4, 4))
    Mr = 0.5 * (Mr + Mr.swapaxes(-1, -2))
    Mr += np.diag([-sg, 1.0, delta, 0.0])
    K1 = np.zeros((4, 4)); K1[0, 1] = K1[1, 0] = 1.0
    K3 = np.zeros((4, 4)); K3[0, 3] = K3[3, 0] = 1.0
    J3 = np.zeros((4, 4)); J3[1, 2], J3[2, 1] = -1.0, 1.0
    J1 = np.zeros((4, 4)); J1[2, 3], J1[3, 2] = -1.0, 1.0
    Gm = 0.07 * K1 - 0.12 * K3 + 0.21 * J3 - 0.09 * J1
    out["gen_eta_antisym_dev"] = float(np.max(np.abs(Gm.T @ ETA
                                                     + ETA @ Gm)))
    L = my_expm(Gm)
    out["L_preserves_eta_dev"] = float(np.max(np.abs(L.T @ ETA @ L
                                                     - ETA)))
    E0 = e_tot(Mr, h6, sg, delta)
    ML = np.einsum("ab,...bc,dc->...ad", L, Mr, L)
    out["own_so13_rel"] = rel(e_tot(ML, h6, sg, delta), E0)
    Lb = L + 0.03 * rng.normal(size=(4, 4))
    MLb = np.einsum("ab,...bc,dc->...ad", Lb, Mr, Lb)
    out["own_negctrl_rel"] = rel(e_tot(MLb, h6, sg, delta), E0)
    out["verdict"] = ("CONFIRMED" if ok and out["own_so13_rel"] < 1e-9
                      and out["own_negctrl_rel"] > 1e-6
                      and out["L_preserves_eta_dev"] < 1e-12
                      else "REFUTED")
    return out


# ---------------- C2: the saddle ----------------
def audit_c2(M, h, sg, delta, w):
    out = {}
    E0 = e_tot(M, h, sg, delta, "pow")
    out["E0_mine_pow"] = E0
    out["E0_mine_eig"] = e_tot(M, h, sg, delta, "eig")
    out["E0_task_hess_row"] = 6.239697665771402
    out["E0_rel_vs_task"] = rel(E0, out["E0_task_hess_row"])
    out["Eu_mine"] = e_u(M, h)
    out["Eu_task_p2_row"] = 4.832936949526051
    out["Eu_rel_vs_task"] = rel(out["Eu_mine"], out["Eu_task_p2_row"])
    rng = np.random.default_rng(12345)     # NOT the task's seed 31
    t = 1e-3
    curvs = []
    P_first = None
    for k in range(6):
        P = np.zeros_like(M)
        for i in range(1, 4):
            f = rng.normal(size=M.shape[:3]) * w
            P[..., 0, i] = f
            P[..., i, 0] = f
        P /= np.sqrt(np.sum(P * P))
        if P_first is None:
            P_first = P
        c = (e_tot(M + t * P, h, sg, delta) - 2 * E0
             + e_tot(M - t * P, h, sg, delta)) / t ** 2
        curvs.append(float(c))
    out["curv_timemix_fresh"] = curvs
    out["n_neg"] = int(sum(c < 0 for c in curvs))
    # t-robustness on the first direction
    rob = {}
    for tt in (5e-4, 2e-3):
        c = (e_tot(M + tt * P_first, h, sg, delta) - 2 * E0
             + e_tot(M - tt * P_first, h, sg, delta)) / tt ** 2
        rob[f"t={tt:g}"] = float(c)
    out["t_robustness_dir0"] = rob
    # block-diagonal control directions (should be positive)
    ctrl = []
    for _ in range(3):
        P = np.zeros_like(M)
        B = rng.normal(size=M.shape[:3] + (3, 3))
        B = 0.5 * (B + B.swapaxes(-1, -2))
        P[..., 1:4, 1:4] = B * w[..., None, None]
        P /= np.sqrt(np.sum(P * P))
        c = (e_tot(M + t * P, h, sg, delta) - 2 * E0
             + e_tot(M - t * P, h, sg, delta)) / t ** 2
        ctrl.append(float(c))
    out["curv_ctrl_fresh"] = ctrl
    all_neg = out["n_neg"] == 6
    in_range = all(-0.30 < c < -0.20 for c in curvs)
    ctrl_pos = all(c > 0 for c in ctrl)
    out["verdict"] = ("CONFIRMED" if all_neg and in_range and ctrl_pos
                      and out["E0_rel_vs_task"] < 1e-8
                      else "REFUTED")
    return out


# ---------------- C3: kin sign table ----------------
def audit_c3(M, h, w):
    out = {}
    Wrz = np.zeros((4, 4)); Wrz[1, 2], Wrz[2, 1] = -1.0, 1.0
    Wbz = np.zeros((4, 4)); Wbz[0, 3] = Wbz[3, 0] = 1.0
    a_rz, raw_rz = a0_global(M, Wrz, w, sign=-1.0)
    a_bz, raw_bz = a0_global(M, Wbz, w, sign=-1.0)
    out["kin_rot_z_mine"] = kin(M, a_rz, h)
    out["kin_boost_z_mine"] = kin(M, a_bz, h)
    out["kin_rot_z_task"] = 0.25750939814530854
    out["kin_boost_z_task"] = -0.08344602895705688
    out["rot_z_rel"] = rel(out["kin_rot_z_mine"], out["kin_rot_z_task"])
    out["boost_z_rel"] = rel(out["kin_boost_z_mine"],
                             out["kin_boost_z_task"])
    # cross-stencil (the row_confirm block, independently)
    out["kin_boost_z_fwd_mine"] = kin(M, a_bz, h, "fwd")
    out["kin_boost_z_bwd_mine"] = kin(M, a_bz, h, "bwd")
    out["kin_boost_z_cen_mine"] = kin(M, a_bz, h, "cen")
    # STRUCTURAL: the task a0 = G M - M G^T is ANTISYMMETRIC for
    # symmetric M; the conjugation-orbit tangent (G1 gate / films
    # convention M -> L M L^T) is G M + M G^T, symmetric.
    out["a0_rot_z_antisym_resid"] = float(
        np.max(np.abs(raw_rz + raw_rz.swapaxes(-1, -2)))
        / np.max(np.abs(raw_rz)))
    out["a0_boost_z_antisym_resid"] = float(
        np.max(np.abs(raw_bz + raw_bz.swapaxes(-1, -2)))
        / np.max(np.abs(raw_bz)))
    a_rz_c, _ = a0_global(M, Wrz, w, sign=+1.0)
    a_bz_c, _ = a0_global(M, Wbz, w, sign=+1.0)
    out["kin_rot_z_conjtangent"] = kin(M, a_rz_c, h)
    out["kin_boost_z_conjtangent"] = kin(M, a_bz_c, h)
    ok = (out["rot_z_rel"] < 1e-3 and out["boost_z_rel"] < 1e-3
          and out["kin_rot_z_mine"] > 0 and out["kin_boost_z_mine"] < 0)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    return out


# ---------------- C4: the decoupling ----------------
def audit_c4():
    out = {}
    with open(os.path.join(DATA, "m5_21_3_all.json")) as f:
        allr = json.load(f)
    lad = allr["p3_p1_s+1_boost_z"]["ladder"]
    r08 = [r for r in lad if r["omega"] == 0.8][0]
    ctrl = allr["ctrl_p1_s+1"]
    diff = r08["E"] - ctrl["E_end"]
    w2k = 0.8 ** 2 * r08["kin"]
    out["E_star_0p8"] = r08["E"]
    out["E_ctrl"] = ctrl["E_end"]
    out["diff"] = diff
    out["omega2_kin_end"] = w2k
    out["ratio"] = diff / w2k
    out["pct_dev"] = abs(diff - w2k) / abs(w2k) * 100.0
    # internal consistency of the recorded row
    out["row_internal_resid"] = abs(
        r08["E_u"] + r08["E_v"] + 0.8 ** 2 * r08["kin"] - r08["E"])
    # static parts vs control (the "4 decimals" wording check)
    out["Eu_gap"] = abs(r08["E_u"] - ctrl["E_u"])
    out["Ev_gap"] = abs(r08["E_v"] - ctrl["E_v"])
    # soft cross-check: kin re-derived from the saved ladder endpoint
    Mp, hp, sp, dp = load_npz("m5_21_3_p3_p1_s+1_boost_z.npz")
    wp = envelope(Mp.shape[0], hp)
    Wbz = np.zeros((4, 4)); Wbz[0, 3] = Wbz[3, 0] = 1.0
    a_bz, _ = a0_global(Mp, Wbz, wp, sign=-1.0)
    out["kin_from_saved_endpoint"] = kin(Mp, a_bz, hp)
    out["kin_endpoint_note"] = ("recorded kin_end used the rung's "
                                "frozen a0 from the warm start, so "
                                "only approximate re-derivation is "
                                "possible from the npz")
    out["verdict"] = "CONFIRMED" if out["pct_dev"] < 2.0 else "REFUTED"
    return out


# ---------------- C5: the twist null ----------------
def audit_c5(M, h, w):
    out = {}
    a_cl = a0_clock(M, w)
    ks = [-0.03, -0.015, 0.0, 0.015, 0.03]   # my own k grid
    es = [e_u(M, h, az_extra=k * a_cl) for k in ks]
    co = np.polyfit(np.array(ks), np.array(es), 2)
    d, b, e0 = float(co[0]), float(co[1]), float(co[2])
    kstar = -b / (2.0 * d)
    out["ks"] = ks
    out["E_at_k"] = es
    out["b_polyfit"] = b
    out["d_polyfit"] = d
    out["k_star"] = kstar
    out["dE_at_kstar"] = b * kstar + d * kstar ** 2
    out["b_direct_oddpart"] = (es[3] - es[1]) / (2 * 0.015)
    out["palindrome_dev"] = max(abs(es[0] - es[4]), abs(es[1] - es[3]))
    out["E_scale"] = es[2]
    ok = (abs(b) < 1e-9 and abs(out["dE_at_kstar"]) < 1e-18
          and abs(out["b_direct_oddpart"]) < 1e-9)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    return out


# ---------------- C6: block diagonality ----------------
def audit_c6():
    out = {}
    ok = True
    for name in ("m5_21_3_p1_s+1.npz", "m5_21_3_p1_s-1.npz"):
        M, _, _, _ = load_npz(name)
        off = float(np.max(np.abs(M[..., 0, 1:4])))
        out[name] = off
        ok = ok and off == 0.0
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    return out


def main():
    res = {"audit": "M5.21.3 independent adversarial audit",
           "auditor": "second-agent re-implementation, fresh seeds",
           "date": "2026-07-18"}

    print("== C1: gates ==", flush=True)
    res["C1_gates"] = audit_c1()
    print(json.dumps({k: v for k, v in res["C1_gates"].items()
                      if k != "gates_compare"}, indent=1), flush=True)

    M, h, s, delta = load_npz("m5_21_3_p1_s+1.npz")
    sg = s * G_COUPL
    w = envelope(M.shape[0], h)

    print("== C2: saddle curvatures (fresh directions) ==", flush=True)
    res["C2_saddle"] = audit_c2(M, h, sg, delta, w)
    print(json.dumps(res["C2_saddle"], indent=1), flush=True)

    print("== C3: kin sign table ==", flush=True)
    res["C3_kin"] = audit_c3(M, h, w)
    print(json.dumps(res["C3_kin"], indent=1), flush=True)

    print("== C4: decoupling ==", flush=True)
    res["C4_decoupling"] = audit_c4()
    print(json.dumps(res["C4_decoupling"], indent=1), flush=True)

    print("== C5: twist null ==", flush=True)
    res["C5_twist"] = audit_c5(M, h, w)
    print(json.dumps(res["C5_twist"], indent=1), flush=True)

    print("== C6: block diagonality ==", flush=True)
    res["C6_blockdiag"] = audit_c6()
    print(json.dumps(res["C6_blockdiag"], indent=1), flush=True)

    # C7 wording findings are assembled by the auditor from the numbers
    res["C7_wording"] = {
        "verdict": "SEE REPORT",
        "flags": [
            "Sec 6 'static parts match the control to 4 decimals': "
            "measured gaps are E_u 2.1e-4, E_V 6.8e-4 (3 decimals for "
            "E_u, between 2 and 3 for E_V); '4 decimals' overstates",
            "Sec 6 verdict sentence 'on NEITHER reading does free 4x4 "
            "energy minimization select nonzero time derivatives': on "
            "the eta-reading E*(omega) is monotone DECREASING (every "
            "probed omega beats the equal-depth static control), so "
            "descent does prefer omega != 0; supported claim is 'no "
            "FINITE stationary omega* / no stable dynamical state is "
            "landed', which the same paragraph states correctly",
            "Sec 8 points the discretization burden at a 'Sec 9 "
            "confirm block' that is not written in the note (the "
            "confirm DATA row exists: data/m5_21_3_row_confirm.json)",
            "a0 convention: task a0 = G M - M G^T is exactly "
            "ANTISYMMETRIC for symmetric M (see C3 antisym resid), "
            "so dM/dt = omega a0 leaves the symmetric configuration "
            "space; the films and gate G1 use M -> L M L^T whose "
            "tangent is G M + M G^T (symmetric); C3 also reports kin "
            "under the conjugation-tangent variant for comparison",
        ]}
    with open(os.path.join(DATA, "m5_21_3_audit.json"), "w") as f:
        json.dump(res, f, indent=1)
    print("\n== VERDICTS ==")
    for k in ("C1_gates", "C2_saddle", "C3_kin", "C4_decoupling",
              "C5_twist", "C6_blockdiag"):
        print(f"  {k:16s} {res[k]['verdict']}")
    print("wrote data/m5_21_3_audit.json")


if __name__ == "__main__":
    main()
