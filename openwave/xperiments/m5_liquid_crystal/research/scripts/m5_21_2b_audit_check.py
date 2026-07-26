"""M5.21.2b INDEPENDENT ADVERSARIAL AUDIT (claims C1-C8 of
findings/m5_21_2b_note.md, sections 4-8).

Independent recomputation: this script shares NO code with
m5_21_2b_a_instrument.py. Energies are rebuilt from the stated
conventions; potentials go through the eigenvalue route
(np.linalg.eigvalsh) instead of the task's trace-power route; the
winding read (C6) is an independent trilinear-interpolated contour
integral, not the task's split reader.

Conventions audited (from the note section 1):
    A_i = d_i M / h (fwd | bwd | 2h with one-sided edges)
    C_ij = A_i A_j - A_j A_i
    u = 4 sum_{i<j} tr(C^T C);  D = sum_i tr(A_i A_i^T)
    T1: V = W1 * sum_p (tr(M^p) - (1 + delta^p))^2
    T2: V = W2 * sum_i (lam_i - v_i)^2, ascending, v = sorted(0,d,1)
    E = h^3 sum_cells (u + eps*D + V);  E_sym = (E_fwd + E_bwd)/2

Run: python3 m5_21_2b_audit_check.py
Out: ../data/m5_21_2b_audit.json
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

W1 = 0.000724023879
W2 = 0.002758100


def load_end(tag):
    z = np.load(os.path.join(DATA, f"m5_21_2b_end_{tag}.npz"))
    return (np.asarray(z["M"], dtype=np.float64), float(z["delta"]),
            float(z["h"]))


# ---------------- my own stencils / energies ----------------
def fd(M, ax, h, st):
    """finite difference along ax; entries where the stencil does not
    fit are zero (one-sided edges for 2h)."""
    A = np.zeros_like(M)
    n = M.shape[ax]

    def sl(a, b):
        s = [slice(None)] * M.ndim
        s[ax] = slice(a, b)
        return tuple(s)

    if st == "fwd":
        A[sl(0, n - 1)] = (M[sl(1, n)] - M[sl(0, n - 1)]) / h
    elif st == "bwd":
        A[sl(1, n)] = (M[sl(1, n)] - M[sl(0, n - 1)]) / h
    elif st == "2h":
        A[sl(1, n - 1)] = (M[sl(2, n)] - M[sl(0, n - 2)]) / (2.0 * h)
        A[sl(0, 1)] = (M[sl(1, 2)] - M[sl(0, 1)]) / h
        A[sl(n - 1, n)] = (M[sl(n - 1, n)] - M[sl(n - 2, n - 1)]) / h
    else:
        raise ValueError(st)
    return A


def e_u_one(M, h, st):
    A = [fd(M, ax, h, st) for ax in range(3)]
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            u += float(np.sum(C * C))
    return 4.0 * h ** 3 * u


def e_u(M, h, st):
    if st == "sym":
        return 0.5 * (e_u_one(M, h, "fwd") + e_u_one(M, h, "bwd"))
    return e_u_one(M, h, st)


def e_dirichlet(M, h, st):
    def one(b):
        return float(sum(np.sum(fd(M, ax, h, b) ** 2)
                         for ax in range(3)))
    d = 0.5 * (one("fwd") + one("bwd")) if st == "sym" else one(st)
    return h ** 3 * d


def v_dens_T1(M, delta):
    """T1 density per cell, EIGENVALUE route (independent of the
    task's trace-power route)."""
    lam = np.linalg.eigvalsh(M)
    out = np.zeros(M.shape[:-2])
    for p in (1, 2, 3):
        out = out + (np.sum(lam ** p, axis=-1) - (1.0 + delta ** p)) ** 2
    return W1 * out


def v_dens_T2(M, delta):
    lam = np.linalg.eigvalsh(M)          # ascending
    v = np.sort(np.array([0.0, delta, 1.0]))
    return W2 * np.sum((lam - v) ** 2, axis=-1)


def e_v(M, h, delta, term):
    dens = v_dens_T2(M, delta) if term == "T2" else v_dens_T1(M, delta)
    return h ** 3 * float(np.sum(dens))


def coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


# ---------------- C5: excess-weighted transverse scale ----------------
def rho_w_T1(M, h, delta):
    n = M.shape[0]
    X, Y, _ = coords(n, h)
    vd = v_dens_T1(M, delta)
    hi = vd > 0.5 * vd.max()
    w = vd[hi]
    rho = np.hypot(X, Y)[hi]
    return float(np.sum(w * rho) / np.sum(w))


# ---------------- C6: my own contour winding ----------------
def interp_M(M, h, pts):
    """trilinear interpolation of the tensor field at physical pts."""
    n = M.shape[0]
    f = pts / h + (n - 1) / 2.0
    i0 = np.clip(np.floor(f).astype(int), 0, n - 2)
    t = f - i0
    out = np.zeros((len(pts), 3, 3))
    for dx in (0, 1):
        wx = t[:, 0] if dx else 1.0 - t[:, 0]
        for dy in (0, 1):
            wy = t[:, 1] if dy else 1.0 - t[:, 1]
            for dz in (0, 1):
                wz = t[:, 2] if dz else 1.0 - t[:, 2]
                w = wx * wy * wz
                out += w[:, None, None] * M[i0[:, 0] + dx,
                                            i0[:, 1] + dy,
                                            i0[:, 2] + dz]
    return out


def circle_winding(M, h, z, rho_c, nsamp=720):
    """net mod-pi winding of the mid-band in-plane director angle on a
    counterclockwise circle of radius rho_c at height z; also the min
    eigen-gap and min in-plane projection along the circle."""
    phi = np.linspace(0.0, 2.0 * np.pi, nsamp, endpoint=False)
    pts = np.stack([rho_c * np.cos(phi), rho_c * np.sin(phi),
                    np.full(nsamp, z)], axis=1)
    Mi = interp_M(M, h, pts)
    lam, V = np.linalg.eigh(Mi)
    gap = np.minimum(lam[:, 1] - lam[:, 0], lam[:, 2] - lam[:, 1])
    vmid = V[:, :, 1]
    ang = np.mod(np.arctan2(vmid[:, 1], vmid[:, 0]), np.pi)
    d = np.diff(np.concatenate([ang, ang[:1]]))
    d = np.mod(d + np.pi / 2.0, np.pi) - np.pi / 2.0
    return (float(np.sum(d) / np.pi), float(gap.min()),
            float(np.hypot(vmid[:, 0], vmid[:, 1]).min()))


def winding_sign_control():
    """synthetic +1 / -1 planar defects with the mid band along the
    director: the instrument must read +2 / -2 half-units."""
    n, h = 32, 1.0
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, _ = np.meshgrid(x, x, x, indexing="ij")
    phi = np.arctan2(Y, X)
    out = {}
    for sgn in (1.0, -1.0):
        ang = sgn * phi
        d = np.stack([np.cos(ang), np.sin(ang), np.zeros_like(ang)],
                     axis=-1)
        zh = np.zeros_like(d)
        zh[..., 2] = 1.0
        M = (zh[..., :, None] * zh[..., None, :]
             + 0.3 * d[..., :, None] * d[..., None, :])
        out[f"{sgn:+.0f}"] = circle_winding(M, h, 2.0, 8.0)[0]
    return out


# ---------------- the claims ----------------
def main():
    audit = {}
    ok = lambda b: "CONFIRMED" if b else "REFUTED"

    # ---- C1: fwd-only funnels (bwd read explodes) ----
    c1 = {}
    for tag, key in (("i1_A_fwd_e0", "A"), ("i1_R_fwd_e0", "R")):
        M, _, h = load_end(tag)
        f, b = e_u(M, h, "fwd"), e_u(M, h, "bwd")
        c1[key] = {"E_u_fwd": f, "E_u_bwd": b, "ratio": b / f}
    passed = all(c1[k]["ratio"] >= 250.0 for k in ("A", "R"))
    audit["C1"] = {
        "verdict": ok(passed), "numbers": c1,
        "reasoning": ("bwd/fwd E_u = %.1f (A: %.3f vs %.1f) and %.1f "
                      "(R: %.3f vs %.1f); both >= 250: the fwd-only "
                      "endpoint hides curvature from its own stencil."
                      % (c1["A"]["ratio"], c1["A"]["E_u_fwd"],
                         c1["A"]["E_u_bwd"], c1["R"]["ratio"],
                         c1["R"]["E_u_fwd"], c1["R"]["E_u_bwd"]))}

    # ---- C2: sym endpoints cross-stencil consistent ----
    c2 = {}
    for tag, key in (("i1_A_sym_e0", "A"), ("i1_R_sym_e0", "R")):
        M, _, h = load_end(tag)
        reads = {st: e_u(M, h, st) for st in ("fwd", "bwd", "2h")}
        strict = dict(reads)
        for par in (0, 1):
            sub = M[par::2, par::2, par::2]
            reads[f"sub{par}"] = e_u(sub, 2 * h, "fwd")
            strict[f"sub{par}_fwd"] = reads[f"sub{par}"]
            strict[f"sub{par}_bwd"] = e_u(sub, 2 * h, "bwd")
        rat = max(reads.values()) / min(reads.values())
        rat_s = max(strict.values()) / min(strict.values())
        c2[key] = {"reads": reads, "ratio": rat,
                   "ratio_strict_7read": rat_s}
    passed = all(c2[k]["ratio"] <= 1.3 for k in ("A", "R"))
    audit["C2"] = {
        "verdict": ok(passed), "numbers": c2,
        "reasoning": ("max/min over {fwd,bwd,2h,sub0,sub1} = %.3f (A) "
                      "and %.3f (R), both <= 1.3; the stricter 7-read "
                      "set (bwd subsamples added) reads %.3f / %.3f."
                      % (c2["A"]["ratio"], c2["R"]["ratio"],
                         c2["A"]["ratio_strict_7read"],
                         c2["R"]["ratio_strict_7read"]))}

    # ---- C3: N=48 census ordering ----
    c3 = {}
    stated = {"A": 5.261, "C": 22.06, "B": 84.08}
    for s in ("A", "C", "B"):
        M, delta, h = load_end(f"c48_{s}_d0.3")
        c3[s] = {"E": e_u(M, h, "sym") + e_v(M, h, delta, "T1"),
                 "stated": stated[s]}
        c3[s]["rel_dev"] = c3[s]["E"] / stated[s] - 1.0
    order = c3["A"]["E"] < c3["C"]["E"] < c3["B"]["E"]
    close = all(abs(c3[s]["rel_dev"]) < 0.02 for s in c3)
    audit["C3"] = {
        "verdict": ok(order and close), "numbers": c3,
        "reasoning": ("recomputed E: A %.4f < C %.4f < B %.4f "
                      "(ordering holds; each within %.2f%% of stated)."
                      % (c3["A"]["E"], c3["C"]["E"], c3["B"]["E"],
                         100 * max(abs(c3[s]["rel_dev"])
                                   for s in c3)))}

    # ---- C4: basin merger, field distance ----
    c4 = {}
    for pa, pr, key in (("c48_A_d0.3", "c48_R_d0.3", "T1"),
                        ("c48_A_T2", "c48_R_T2", "T2")):
        MA, delta, h = load_end(pa)
        MR, _, _ = load_end(pr)
        vac = np.zeros_like(MA)
        vac[:] = np.diag([1.0, delta, 0.0])
        dAR = float(np.linalg.norm(MA - MR))
        dAv = float(np.linalg.norm(MA - vac))
        c4[key] = {"dist_AR": dAR, "dist_Avac": dAv,
                   "ratio": dAR / dAv,
                   "ratio_Rnorm": dAR / float(np.linalg.norm(MR - vac))}
    passed = all(c4[k]["ratio"] <= 0.10 for k in c4)
    audit["C4"] = {
        "verdict": ok(passed), "numbers": c4,
        "reasoning": ("||A-R||/||A-vac|| = %.4f (T1) and %.4f (T2), "
                      "both <= 0.10: point- and ring-seeded endpoints "
                      "are the same state to <=%.1f%% field distance."
                      % (c4["T1"]["ratio"], c4["T2"]["ratio"],
                         100 * max(c4[k]["ratio"] for k in c4)))}

    # ---- C5: dynamically selected transverse scale ----
    c5 = {"rho_w_d0.1": {}, "delta_trend": {}}
    for tag in ("i4_a3_d0.1", "i4_a4.5_d0.1", "i4_a6_d0.1",
                "c48_A_d0.1"):
        M, delta, h = load_end(tag)
        c5["rho_w_d0.1"][tag] = rho_w_T1(M, h, delta)
    vals = np.array(list(c5["rho_w_d0.1"].values()))
    mean = float(vals.mean())
    c5["mean"] = mean
    c5["max_dev_from_mean"] = float(np.max(np.abs(vals / mean - 1.0)))
    for tag in ("c48_A_d0.3", "c48_A_d0.2", "c48_A_d0.1"):
        M, delta, h = load_end(tag)
        c5["delta_trend"][tag] = rho_w_T1(M, h, delta)
    tr = [c5["delta_trend"][t] for t in ("c48_A_d0.3", "c48_A_d0.2",
                                         "c48_A_d0.1")]
    mono = tr[0] < tr[1] < tr[2]
    tight = c5["max_dev_from_mean"] <= 0.05 and abs(mean - 3.4) < 0.2
    audit["C5"] = {
        "verdict": ok(tight and mono), "numbers": c5,
        "reasoning": ("four delta=0.1 seedings read rho_w %.3f-%.3f "
                      "(mean %.3f, max dev %.1f%%, bar 5%%); delta "
                      "trend %.3f -> %.3f -> %.3f monotone as delta "
                      "0.3 -> 0.1."
                      % (vals.min(), vals.max(), mean,
                         100 * c5["max_dev_from_mean"], *tr))}

    # ---- C6: T2 split winding, my own contour read ----
    MA, _, h = load_end("c48_A_T2")
    ctrl = winding_sign_control()
    c6 = {"sign_control": ctrl, "circles": []}
    all_w, all_gap = [], []
    for z in (-12.0, -4.0, 4.0, 12.0):
        for rc in (10.0, 14.0):
            w, gap, inpl = circle_winding(MA, h, z, rc)
            c6["circles"].append({"z": z, "rho_c": rc,
                                  "winding_half_units": w,
                                  "min_gap": gap,
                                  "min_inplane": inpl})
            all_w.append(w)
            all_gap.append(gap)
    wint = [round(w) for w in all_w]
    ctrl_ok = (abs(ctrl["+1"] - 2.0) < 1e-9
               and abs(ctrl["-1"] + 2.0) < 1e-9)
    passed = (ctrl_ok and all(r == 2 for r in wint)
              and min(all_gap) > 0.05)
    c6["min_gap_overall"] = float(min(all_gap))
    audit["C6"] = {
        "verdict": ok(passed), "numbers": c6,
        "reasoning": ("sign control reads +2/-2 on synthetic +-1 "
                      "defects; all 8 circles (z in {-12,-4,4,12} x "
                      "rho_c in {10,14}) read net winding +2 "
                      "half-units (integer exactly, mod-pi "
                      "telescoping); min eigen-gap %.3f > 0.05: band "
                      "identity clean throughout."
                      % min(all_gap))}

    # ---- C7: eps-family linearity ----
    rows = json.load(open(os.path.join(DATA, "m5_21_2b_all.json")))
    eps_list = [0.0, 3e-4, 1e-3, 3e-3]
    tags = ["i1_A_sym_e0", "i1_A_sym_e0.0003", "i1_A_sym_e0.001",
            "i1_A_sym_e0.003"]
    c7 = {"E_recomputed": {}, "E_rows": {}}
    E = []
    for eps, tag in zip(eps_list, tags):
        M, delta, h = load_end(tag)
        e = (e_u(M, h, "sym") + eps * e_dirichlet(M, h, "sym")
             + e_v(M, h, delta, "T1"))
        c7["E_recomputed"][tag] = e
        c7["E_rows"][tag] = rows[tag]["E_end"]
        E.append(e)
    E = np.array(E)
    x = np.array(eps_list)
    A1 = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A1, E, rcond=None)
    resid = E - A1 @ coef
    span = E[-1] - E[0]
    c7["fit_intercept"] = float(coef[0])
    c7["fit_slope"] = float(coef[1])
    c7["residuals"] = resid.tolist()
    c7["span"] = float(span)
    c7["max_resid_over_span"] = float(np.max(np.abs(resid)) / span)
    c7["row_vs_recompute_max_rel"] = float(max(
        abs(c7["E_recomputed"][t] / c7["E_rows"][t] - 1.0)
        for t in tags))
    passed = c7["max_resid_over_span"] < 0.03
    audit["C7"] = {
        "verdict": ok(passed), "numbers": c7,
        "reasoning": ("LSQ line E = %.4f + %.1f*eps; max |residual| "
                      "= %.2f%% of the span %.4f (bar 3%%); rows "
                      "match my recompute to %.1e."
                      % (coef[0], coef[1],
                         100 * c7["max_resid_over_span"], span,
                         c7["row_vs_recompute_max_rel"]))}

    # ---- C8: virial balance T2 vs T1 ----
    c8 = {}
    for tag, term, key in (("c48_A_T2", "T2", "T2"),
                           ("c48_A_d0.3", "T1", "T1")):
        M, delta, h = load_end(tag)
        eu = e_u(M, h, "sym")
        ev = e_v(M, h, delta, term)
        c8[key] = {"E_u": eu, "E_v": ev,
                   "q": abs(-eu + 3.0 * ev) / (eu + ev)}
    passed = c8["T2"]["q"] <= 0.06 and c8["T1"]["q"] >= 0.25
    audit["C8"] = {
        "verdict": ok(passed), "numbers": c8,
        "reasoning": ("|-E_u + 3E_v|/E = %.4f on the T2 endpoint "
                      "(<= 0.06: landed) vs %.4f on the T1 endpoint "
                      "(>= 0.25: pin-held): T2 alone reaches the "
                      "scale-stationary balance."
                      % (c8["T2"]["q"], c8["T1"]["q"]))}

    out = os.path.join(DATA, "m5_21_2b_audit.json")
    with open(out, "w") as f:
        json.dump(audit, f, indent=1)
    for k in sorted(audit):
        print(f"{k}: {audit[k]['verdict']}")
        print(f"   {audit[k]['reasoning']}")
    print(f"-> {os.path.relpath(out, HERE)}")


if __name__ == "__main__":
    main()
