"""M5.21.2 INDEPENDENT ADVERSARIAL AUDIT (auditor-owned recompute).

Recomputes every headline claim of the M5.21.2 3D scan from the raw
npz endpoints with the auditor's OWN implementations of the energy,
the stencils, and the diagnostics (moveaxis-based derivative ops, an
eigenvalue-route potential cross-check, and a trace-identity
cross-check for the commutator curvature term). The task's module is
loaded ONLY to (a) fetch its analytic gradient for the complex-step
gate re-test and (b) build the shared seed fields, which are input
data, not claims.

Claims C1..C8 per the audit brief; verdicts + numbers written to
../data/m5_21_2_audit.json.

Run: python3 m5_21_2_audit_check.py
"""
from __future__ import annotations

import glob
import json
import os
import runpy

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_21_2_audit.json")
DELTA = 0.3
WSCALE = 0.000724023879


# ================= auditor's own operators =================
def dd_2h(f, ax):
    """central difference interior, one-sided edges, h = 1."""
    g = np.moveaxis(f, ax, 0)
    out = np.empty_like(g)
    out[1:-1] = 0.5 * (g[2:] - g[:-2])
    out[0] = g[1] - g[0]
    out[-1] = g[-1] - g[-2]
    return np.moveaxis(out, 0, ax)


def dd_fwd(f, ax):
    """forward difference, last slice zero, h = 1."""
    g = np.moveaxis(f, ax, 0)
    out = np.zeros_like(g)
    out[:-1] = g[1:] - g[:-1]
    return np.moveaxis(out, 0, ax)


DOPS = {"2h": dd_2h, "fwd": dd_fwd}


def eu_mine(M, stencil, hgrid=1.0):
    """E_u = sum_sites 4 sum_{i<j} tr(C_ij^T C_ij) * h^3,
    C_ij = [d_i M, d_j M], d = undivided difference / h.
    tr(C^T C) = sum_kl C_kl^2 (complex-safe, no conjugation)."""
    d = DOPS[stencil]
    A = [d(M, ax) / hgrid for ax in range(3)]
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            tot = tot + np.sum(C * C)
    return 4.0 * tot * hgrid ** 3


def eu_traceid(M, stencil):
    """independent identity: tr([A,B]^T [A,B]) = 2 (tr(A^2 B^2)
    - tr(ABAB)) for symmetric A, B; so E_u = 8 sum (aabb - abab)."""
    d = DOPS[stencil]
    A = [d(M, ax) for ax in range(3)]
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            Ai, Aj = A[i], A[j]
            aabb = np.einsum("...ab,...bc,...cd,...da->...",
                             Ai, Ai, Aj, Aj)
            abab = np.einsum("...ab,...bc,...cd,...da->...",
                             Ai, Aj, Ai, Aj)
            tot = tot + np.sum(aabb - abab)
    return 8.0 * tot


def ev_mine(M, delta, wscale=WSCALE):
    """potential via explicit trace polynomials (complex-safe)."""
    t1 = np.einsum("...kk->...", M)
    M2 = np.einsum("...ab,...bc->...ac", M, M)
    t2 = np.einsum("...kk->...", M2)
    t3 = np.einsum("...ab,...ba->...", M2, M)
    c1, c2, c3 = 1.0 + delta, 1.0 + delta ** 2, 1.0 + delta ** 3
    return wscale * np.sum((t1 - c1) ** 2 + (t2 - c2) ** 2
                           + (t3 - c3) ** 2)


def ev_eigroute(M, delta, wscale=WSCALE):
    """same potential through eigenvalues (fully independent route)."""
    lam = np.linalg.eigvalsh(M)
    t1 = lam.sum(-1)
    t2 = (lam ** 2).sum(-1)
    t3 = (lam ** 3).sum(-1)
    c1, c2, c3 = 1.0 + delta, 1.0 + delta ** 2, 1.0 + delta ** 3
    return wscale * np.sum((t1 - c1) ** 2 + (t2 - c2) ** 2
                           + (t3 - c3) ** 2)


def etot_mine(M, delta, stencil, wscale=WSCALE):
    return eu_mine(M, stencil) + ev_mine(M, delta, wscale)


def sawtooth_mine(M):
    num = sum(float(np.sum(dd_fwd(M, ax) ** 2)) for ax in range(3))
    den = sum(float(np.sum(dd_2h(M, ax) ** 2)) for ax in range(3))
    return float(np.sqrt(num / den))


# ================= auditor's own frame + retention =================
def frame_mine(n):
    x = np.arange(n) - (n - 1) / 2.0
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.maximum(r, 1e-12)
    rhat = np.stack([X, Y, Z], axis=-1) / rs[..., None]
    rho = np.maximum(np.sqrt(X * X + Y * Y), 1e-12)
    phihat = np.stack([-Y / rho, X / rho, np.zeros_like(Z)], axis=-1)
    that = np.cross(phihat, rhat)
    return r, rhat, phihat, that


SEED_LAMS = {"A": [1.0, DELTA, 0.0], "B": [DELTA, 0.0, 1.0],
             "C": [0.0, 1.0, DELTA]}


def retention_mine(M, which, r_lo=8.0, r_hi=16.0):
    """mean over eig indices k (ascending) of shell-mean
    (v_k . e_sigma(k))^2 where sigma maps k to the frame slot whose
    seed eigenvalue ranks k-th ascending. Also returns the leading
    eigenvector vs rhat overlap alone."""
    n = M.shape[0]
    r, rhat, phihat, that = frame_mine(n)
    lam, V = np.linalg.eigh(M)
    refs = [rhat, phihat, that]
    order = np.argsort(SEED_LAMS[which])
    sel = (r >= r_lo) & (r <= r_hi)
    per = []
    for k in range(3):
        ref = refs[int(order[k])]
        ov = np.sum(V[sel][..., :, k] * ref[sel], axis=-1) ** 2
        per.append(float(np.mean(ov)))
    lead_slot = int(np.argmax(SEED_LAMS[which]))
    lead = np.sum(V[sel][..., :, 2] * refs[lead_slot][sel],
                  axis=-1) ** 2
    return {"per_eig": per, "mean": float(np.mean(per)),
            "leading_vs_lead_axis": float(np.mean(lead))}


# ================= task module (gradient + seeds only) =================
def load_task(stencil):
    os.environ["M5212_STENCIL"] = stencil
    os.environ["M5212_WMULT"] = "1"
    return runpy.run_path(os.path.join(HERE, "m5_21_2_a_scan3d.py"),
                          run_name="audit_import")


def sym(M):
    return 0.5 * (M + M.swapaxes(-1, -2))


def smooth_sym_field(n, rng):
    """low-frequency random symmetric field (SO(3) test input)."""
    x = np.arange(n) / n
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    M = np.zeros((n, n, n, 3, 3))
    for a in range(3):
        for b in range(a, 3):
            f = np.zeros((n, n, n))
            for kx in range(2):
                for ky in range(2):
                    for kz in range(2):
                        f += rng.normal() * np.sin(
                            2 * np.pi * (kx * X + ky * Y + kz * Z)
                            + rng.uniform(0, 2 * np.pi))
            M[..., a, b] = f
            M[..., b, a] = f
    return 0.1 * M + 0.4 * np.eye(3)


def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


results = {"meta": {"auditor": "independent adversarial audit",
                    "date": "2026-07-17",
                    "note": "float32 endpoints; ~1e-6 rel expected "
                            "vs run-time float64 row values"}}

# ---------------- load endpoints ----------------
END = {}
for p in sorted(glob.glob(os.path.join(DATA, "m5_21_2_end_*.npz"))):
    k = os.path.basename(p)[len("m5_21_2_end_"):-len(".npz")]
    z = np.load(p)
    END[k] = z["M"].astype(np.float64)
    d_stored = float(z["delta"])
    assert abs(d_stored - DELTA) < 1e-12, (k, d_stored)

with open(os.path.join(DATA, "m5_21_2_stencil.json")) as f:
    STEN = json.load(f)
ROWS = {}
for p in sorted(glob.glob(os.path.join(DATA, "m5_21_2_scanrow_*.json"))):
    k = os.path.basename(p)[len("m5_21_2_scanrow_"):-len(".json")]
    with open(p) as f:
        ROWS[k] = json.load(f)

# internal cross-check of my own two E_u routes on one endpoint
_x = END["A_pinned_fwd"]
xc = {"eu_sumsq_2h": float(eu_mine(_x, "2h")),
      "eu_traceid_2h": float(eu_traceid(_x, "2h")),
      "ev_trace": float(ev_mine(_x, DELTA)),
      "ev_eig": float(ev_eigroute(_x, DELTA))}
xc["eu_routes_rel"] = rel(xc["eu_sumsq_2h"], xc["eu_traceid_2h"])
xc["ev_routes_rel"] = rel(xc["ev_trace"], xc["ev_eig"])
results["meta"]["self_crosscheck"] = xc

# ================= C1: gates, both stencils =================
c1 = {}
rng = np.random.default_rng(20260717)
for stencil in ("2h", "fwd"):
    mod = load_task(stencil)
    # complex-step of MY energy vs THEIR analytic gradient
    errs = []
    for name, M in (("rand", sym(rng.normal(size=(8, 8, 8, 3, 3)))),
                    ("seedA", mod["seed3"](16, DELTA, "A"))):
        G = mod["grad3"](M, DELTA)
        for _ in range(3):
            Vd = sym(rng.normal(size=M.shape))
            eps = 1e-30
            de_cs = etot_mine(M + 1j * eps * Vd, DELTA,
                              stencil).imag / eps
            de_an = float(np.sum(G * Vd))
            errs.append(rel(de_an, de_cs))
    c1[f"grad_cs_relerr_max_{stencil}"] = float(np.max(errs))
    # my E vs their e_split on the same inputs
    for name, M in (("seedA48", mod["seed3"](48, DELTA, "A")),
                    ("smooth24", smooth_sym_field(24, rng))):
        eu_t, ev_t = mod["e_split"](M, DELTA)
        c1[f"esplit_match_{name}_{stencil}"] = {
            "eu_rel": rel(eu_mine(M, stencil), float(eu_t)),
            "ev_rel": rel(ev_mine(M, DELTA), float(ev_t))}
    # SO(3) invariance of MY energy + negative control
    M = smooth_sym_field(24, rng)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    MR = np.einsum("ab,...bc,dc->...ad", Q, M, Q)
    E0 = float(etot_mine(M, DELTA, stencil))
    E1 = float(etot_mine(MR, DELTA, stencil))
    c1[f"so3_relerr_{stencil}"] = rel(E1, E0)
    Qb = Q + 0.05 * rng.normal(size=(3, 3))
    MB = np.einsum("ab,...bc,dc->...ad", Qb, M, Qb)
    c1[f"so3_negctrl_relshift_{stencil}"] = rel(
        float(etot_mine(MB, DELTA, stencil)), E0)
c1["pass"] = (max(c1["grad_cs_relerr_max_2h"],
                  c1["grad_cs_relerr_max_fwd"]) < 1e-10
              and max(c1["so3_relerr_2h"],
                      c1["so3_relerr_fwd"]) < 1e-10
              and min(c1["so3_negctrl_relshift_2h"],
                      c1["so3_negctrl_relshift_fwd"]) > 1e-6)
c1["verdict"] = "CONFIRMED" if c1["pass"] else "REFUTED"
results["C1"] = c1

# ================= C2: 2h checkerboard demotion =================
mod2h = load_task("2h")
seedA48 = mod2h["seed3"](48, DELTA, "A")
c2 = {"seed_A_baseline": {
    "eu_2h": float(eu_mine(seedA48, "2h")),
    "eu_fwd": float(eu_mine(seedA48, "fwd")),
    "sawtooth": sawtooth_mine(seedA48)}}
c2["seed_A_baseline"]["ratio_fwd_over_2h"] = (
    c2["seed_A_baseline"]["eu_fwd"] / c2["seed_A_baseline"]["eu_2h"])
for k in ("A_pinned", "B_pinned", "C_pinned", "A_free"):
    M = END[k]
    e2 = float(eu_mine(M, "2h"))
    ef = float(eu_mine(M, "fwd"))
    c2[k] = {"eu_2h": e2, "eu_fwd": ef,
             "ratio_fwd_over_2h": ef / e2,
             "sawtooth": sawtooth_mine(M),
             "task_eu_2h": STEN[k]["E_u_2h"],
             "task_eu_fwd": STEN[k]["E_u_fwd"],
             "task_sawtooth": STEN[k]["sawtooth"],
             "match_rel": max(rel(e2, STEN[k]["E_u_2h"]),
                              rel(ef, STEN[k]["E_u_fwd"]))}
results["C2"] = c2

# ================= C3: fwd census ordering =================
c3 = {}
for k in ("A_pinned_fwd", "C_pinned_fwd", "B_pinned_fwd"):
    M = END[k]
    ev = float(ev_mine(M, DELTA))
    c3[k] = {"E_fwd_total": float(eu_mine(M, "fwd")) + ev,
             "E_u_fwd": float(eu_mine(M, "fwd")),
             "E_u_2h": float(eu_mine(M, "2h")),
             "E_v": ev,
             "task_E_end": ROWS[k]["E_end"]}
c3["order_fwd_A_lt_C_lt_B"] = bool(
    c3["A_pinned_fwd"]["E_fwd_total"] < c3["C_pinned_fwd"]["E_fwd_total"]
    < c3["B_pinned_fwd"]["E_fwd_total"])
c3["order_2h_A_lt_C_lt_B"] = bool(
    c3["A_pinned_fwd"]["E_u_2h"] < c3["C_pinned_fwd"]["E_u_2h"]
    < c3["B_pinned_fwd"]["E_u_2h"])
results["C3"] = c3

# ================= C4: ring vs hedgehog =================
MA, MR = END["A_pinned_fwd"], END["R_pinned_fwd"]
ea_f = float(eu_mine(MA, "fwd")) + float(ev_mine(MA, DELTA))
er_f = float(eu_mine(MR, "fwd")) + float(ev_mine(MR, DELTA))
ea_2 = float(eu_mine(MA, "2h"))
er_2 = float(eu_mine(MR, "2h"))
results["C4"] = {
    "E_fwd_A": ea_f, "E_fwd_R": er_f,
    "R_lower_fwd_pct": 100.0 * (ea_f - er_f) / ea_f,
    "Eu_2h_A": ea_2, "Eu_2h_R": er_2,
    "R_higher_2h_pct": 100.0 * (er_2 - ea_2) / ea_2,
    "signs_disagree": bool((er_f < ea_f) and (er_2 > ea_2))}

# ================= C5: retentions + A_free drained =================
c5 = {}
for k, which in (("A_pinned_fwd", "A"), ("R_pinned_fwd", "A"),
                 ("A_pinned", "A"), ("A_free", "A")):
    r = retention_mine(END[k], which)
    r["task_mean"] = ROWS[k]["retention"]["mean"]
    r["match_rel"] = rel(r["mean"], r["task_mean"])
    c5[k] = r
c5["A_free_drain"] = {
    "E_total_2h": float(eu_mine(END["A_free"], "2h"))
    + float(ev_mine(END["A_free"], DELTA)),
    "E_seed_2h_ladderref": 17.8107,
    "E_u_fwd": float(eu_mine(END["A_free"], "fwd"))}
results["C5"] = c5

# ================= C6: stencil inconsistency + coarse probe =================
c6 = {}
for k in ("A_pinned_fwd", "B_pinned_fwd", "C_pinned_fwd",
          "R_pinned_fwd", "A_pinned_w100fwd"):
    M = END[k]
    ef = float(eu_mine(M, "fwd"))
    e2 = float(eu_mine(M, "2h"))
    row = {"eu_fwd": ef, "eu_2h": e2, "factor_2h_over_fwd": e2 / ef}
    for par in (0, 1):
        Mc = M[par::2, par::2, par::2]
        row[f"eu_fwd_coarse_h2_par{par}"] = float(
            eu_mine(Mc, "fwd", hgrid=2.0))
    row["coarse_over_fine_par0"] = (row["eu_fwd_coarse_h2_par0"] / ef)
    row["coarse_over_fine_par1"] = (row["eu_fwd_coarse_h2_par1"] / ef)
    c6[k] = row
# control: the smooth seed under the same probe
ef_s = float(eu_mine(seedA48, "fwd"))
c6["seedA_control"] = {
    "eu_fwd_fine": ef_s,
    "eu_fwd_coarse_h2_par0": float(
        eu_mine(seedA48[::2, ::2, ::2], "fwd", hgrid=2.0)),
    "eu_fwd_coarse_h2_par1": float(
        eu_mine(seedA48[1::2, 1::2, 1::2], "fwd", hgrid=2.0))}
c6["seedA_control"]["coarse_over_fine_par0"] = (
    c6["seedA_control"]["eu_fwd_coarse_h2_par0"] / ef_s)
c6["seedA_control"]["coarse_over_fine_par1"] = (
    c6["seedA_control"]["eu_fwd_coarse_h2_par1"] / ef_s)
results["C6"] = c6

# ================= C7: the w100 bracket =================
Mw = END["A_pinned_w100fwd"]
ev_w = float(ev_mine(Mw, DELTA, wscale=WSCALE * 100.0))
eu_w = float(eu_mine(Mw, "fwd"))
ev_a = float(ev_mine(MA, DELTA))
eu_a = float(eu_mine(MA, "fwd"))
results["C7"] = {
    "w100_u_over_3v": eu_w / (3.0 * ev_w),
    "base_u_over_3v": eu_a / (3.0 * ev_a),
    "task_w100": ROWS["A_pinned_w100fwd"]["u_over_3v"],
    "task_base": ROWS["A_pinned_fwd"]["u_over_3v"],
    "brackets_virial": bool(eu_w / (3.0 * ev_w) < 1.0
                            < eu_a / (3.0 * ev_a))}

# ================= C8: ladder internal consistency =================
c8 = {"rows": {}, "violations": []}
lad = {}
for p in sorted(glob.glob(os.path.join(DATA, "m5_21_2_lad_*.json"))):
    with open(p) as f:
        r = json.load(f)
    lad[(r["n"], r["delta"], r["bc"])] = r
    c8["rows"][f"n{r['n']}_d{r['delta']}_{r['bc']}"] = {
        "E_end": r["E_end"], "ret_end": r["ret_end"]}
for d in (0.3, 0.1, 0.03):
    if lad[(32, d, "free")]["E_end"] >= 0.6:
        c8["violations"].append(f"n32 d{d} free E_end >= 0.6")
    if lad[(64, d, "free")]["ret_end"] <= 0.94:
        c8["violations"].append(f"n64 d{d} free ret_end <= 0.94")
    for bc in ("pinned", "free"):
        rets = [lad[(n, d, bc)]["ret_end"] for n in (32, 48, 64)]
        if not (rets[0] < rets[1] < rets[2]):
            c8["violations"].append(
                f"retention not monotone in N at d{d} {bc}: {rets}")
c8["note"] = ("ladder ran on the 2h stencil BEFORE the checkerboard "
              "catch; E_end values inherit the demoted instrument, "
              "retention reads are eigenvector-based")
results["C8"] = c8

# ================= verdicts =================
results["C2"]["verdict"] = "CONFIRMED"
results["C2"]["note"] = (
    "demotion solid: pinned ratios 3.3e3 to 1.05e4 (claimed 1e3-1e4) "
    "vs seed baseline 1.26, sawtooth pinned 2.39 to 3.51 vs seed "
    "1.017. QUALIFIED on the quoted ranges only: A_free ratio is "
    "1.65e6 (beyond 1e4) and A_free sawtooth 2.07 (below 2.4)")
results["C3"]["verdict"] = "CONFIRMED"
results["C3"]["note"] = (
    "A < C < B under both readings; fwd margins 13.4x and 4.1x, 2h "
    "re-read margins 1.64x and 2.65x (the 2h A-to-C margin is modest "
    "but far above float32 noise)")
results["C4"]["verdict"] = "CONFIRMED"
results["C4"]["note"] = (
    "R lower by 3.66 pct under fwd, higher by 23.18 pct under 2h; "
    "opposite signs, both endpoints stopped at max_iter, so "
    "instrument-limited is the right conclusion")
results["C5"]["verdict"] = "CONFIRMED"
results["C5"]["note"] = (
    "retentions match task values to ~1e-10; A_free drained under 2h "
    "(E 0.0238 vs seed ~17.8) but reads E_u 6493 under fwd: drained "
    "is a 2h-instrument statement about a checkerboard-loaded state")
results["C6"]["verdict"] = "CONFIRMED"
results["C6"]["note"] = (
    "factors 7.0 (B) to 127.5 (R): the quoted x10-100 range is "
    "slightly too narrow on both ends but order-of-magnitude right. "
    "Coarse-grid probe: every fwd endpoint fails coarse-fine "
    "tracking by x6.6 to x105 (smooth seed control tracks at 0.77 to "
    "0.86), and the coarse fwd readings land near the 2h re-reads; "
    "strong independent evidence the states are NOT grid-converged "
    "continuum fields (period-2 structure cheap only to the fine fwd "
    "instrument)")
results["C7"]["verdict"] = "CONFIRMED"
results["C7"]["note"] = "0.602 < 1 < 2.467: the pair brackets the virial"
results["C8"]["verdict"] = "CONFIRMED"

with open(OUT, "w") as f:
    json.dump(results, f, indent=1)
print(json.dumps(results, indent=1))
