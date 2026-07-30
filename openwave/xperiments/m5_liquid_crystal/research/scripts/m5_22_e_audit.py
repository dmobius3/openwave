#!/usr/bin/env python3
"""M5.22 ADVERSARIAL AUDIT (independent re-implementation, AI_HYGIENE.md par. 1).

Self-contained: numpy only (trapezoid integration; no scipy needed). Does NOT import
m5_22_a/b/c or any m5_21_* script. Energy and charge instruments re-derived from the
equations in findings/m5_22_note.md par. 1 and findings/m5_21_2b_note.md par. 1:

  E   = h^3 sum_cells(u + V_T2)
  u   = 4 sum_{i<j} tr([A_i,A_j]^T [A_i,A_j]),  A_i = d_i M / h
  sym = (E_fwd + E_bwd)/2, fwd/bwd one-sided differences
  V_T2= w2 sum_i (lambda_i - v_i)^2, lambda ascending, v = sorted(1, delta, 0)

  charge: degree of the oriented-lift leading eigenvector over a closed centered
  lattice cube, computed here NOT as the note's finite-difference Mermin-Ho flux but
  as the exact sum of spherical solid angles of the image quads (van Oosterom &
  Strackee), an independent method with the same continuum meaning.

2D seed winding (C1): numerical integration of the ANALYTIC d(ang)/dphi along the
semicircle (x,y) = (-d cos phi, d sin phi), d = 10 (not sampling + unwrap).

Verdicts -> ../data/m5_22_audit.json
Run: python3 m5_22_e_audit.py            (from the scripts directory)
"""

import json
import heapq
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

TRAP = getattr(np, "trapezoid", None) or np.trapz


# ---------------------------------------------------------------- energy ----
def _sl(ax, a, b, nd=5):
    s = [slice(None)] * nd
    s[ax] = slice(a, b)
    return tuple(s)


def _u_of(A):
    u = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            u = u + np.einsum("...ij,...ij->...", C, C)
    return 4.0 * u


def eu_side(Mf, h, side, conv):
    """E_u for one one-sided stencil under an explicit edge convention."""
    n = Mf.shape[0]
    d = [(Mf[_sl(ax, 1, n)] - Mf[_sl(ax, 0, n - 1)]) / h for ax in range(3)]
    if conv == "crop":
        lo, hi = (0, n - 1) if side == "fwd" else (1, n)
        A = []
        for ax in range(3):
            x = d[ax]
            for ax2 in range(3):
                if ax2 != ax:
                    x = x[_sl(ax2, lo, hi)]
            A.append(x)
        return (h ** 3) * float(_u_of(A).sum())
    A = []
    for ax in range(3):
        D = np.zeros_like(Mf)
        if side == "fwd":
            D[_sl(ax, 0, n - 1)] = d[ax]
            if conv == "clamp":
                D[_sl(ax, n - 1, n)] = d[ax][_sl(ax, n - 2, n - 1)]
            elif conv == "wrap":
                D[_sl(ax, n - 1, n)] = (Mf[_sl(ax, 0, 1)] - Mf[_sl(ax, n - 1, n)]) / h
        else:
            D[_sl(ax, 1, n)] = d[ax]
            if conv == "clamp":
                D[_sl(ax, 0, 1)] = d[ax][_sl(ax, 0, 1)]
            elif conv == "wrap":
                D[_sl(ax, 0, 1)] = (Mf[_sl(ax, 0, 1)] - Mf[_sl(ax, n - 1, n)]) / h
        A.append(D)
    return (h ** 3) * float(_u_of(A).sum())


def ev_of(Mf, h, delta, w2):
    lam = np.linalg.eigvalsh(Mf)
    v = np.array(sorted((1.0, delta, 0.0)))
    return (h ** 3) * w2 * float(((lam - v) ** 2).sum())


def energy(Mf, h, delta, w2, conv):
    ef = eu_side(Mf, h, "fwd", conv)
    eb = eu_side(Mf, h, "bwd", conv)
    ev = ev_of(Mf, h, delta, w2)
    return {"E_u_fwd": ef, "E_u_bwd": eb, "E_u": 0.5 * (ef + eb),
            "E_v": ev, "E_end": 0.5 * (ef + eb) + ev}


def load_end(tag):
    z = np.load(DATA / f"m5_22_end_{tag}.npz")
    return (np.asarray(z["M"], dtype=np.float64), float(z["h"]),
            float(z["delta"]), z)


# ---------------------------------------------------------------- charge ----
def surface_points(ilo, ihi):
    pts = set()
    rng = range(ilo, ihi + 1)
    for a in range(3):
        bc = [x for x in range(3) if x != a]
        for pos in (ilo, ihi):
            for u in rng:
                for v in rng:
                    idx = [0, 0, 0]
                    idx[a] = pos
                    idx[bc[0]] = u
                    idx[bc[1]] = v
                    pts.add(tuple(idx))
    return sorted(pts)


def lift_signs(pts, vecs):
    """Greedy max-|dot|-first (Prim) orientation lift of a director field."""
    index = {p: k for k, p in enumerate(pts)}
    N = len(pts)
    adj = [[] for _ in range(N)]
    for k, p in enumerate(pts):
        for a in range(3):
            q = list(p)
            q[a] += 1
            j = index.get(tuple(q))
            if j is not None:
                adj[k].append(j)
                adj[j].append(k)
    signs = np.zeros(N, dtype=np.int64)
    signs[0] = 1
    heap = []

    def push(k):
        for j in adj[k]:
            if signs[j] == 0:
                dot = float(vecs[k] @ vecs[j])
                heapq.heappush(heap, (-abs(dot), k, j, dot))

    push(0)
    done = 1
    while done < N and heap:
        _, k, j, dot = heapq.heappop(heap)
        if signs[j] != 0:
            continue
        signs[j] = signs[k] * (1 if dot >= 0 else -1)
        done += 1
        push(j)
    ncf = 0
    for k in range(N):
        for j in adj[k]:
            if j > k and signs[k] * signs[j] * float(vecs[k] @ vecs[j]) < 0:
                ncf += 1
    return signs, ncf


def _tri_omega(a, b, c):
    num = float(a @ np.cross(b, c))
    den = 1.0 + float(a @ b) + float(b @ c) + float(c @ a)
    return 2.0 * np.arctan2(num, den)


def surface_degree(pts, vecs_signed, ilo, ihi):
    """Signed degree = sum of solid angles of image quads / 4pi (outward oriented)."""
    index = {p: k for k, p in enumerate(pts)}
    ident = np.eye(3)
    total = 0.0
    for a in range(3):
        b, c = [x for x in range(3) if x != a]
        par = float(np.cross(ident[b], ident[c]) @ ident[a])
        for pos, sgn in ((ihi, 1.0), (ilo, -1.0)):
            orient = par * sgn
            for u in range(ilo, ihi):
                for v in range(ilo, ihi):
                    q = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        idx = [0, 0, 0]
                        idx[a] = pos
                        idx[b] = u + du
                        idx[c] = v + dv
                        q.append(vecs_signed[index[tuple(idx)]])
                    om = _tri_omega(q[0], q[1], q[2]) + _tri_omega(q[0], q[2], q[3])
                    total += orient * om
    return total / (4.0 * np.pi)


def read_charge_from_M(Mf, ilo, ihi):
    pts = surface_points(ilo, ihi)
    mats = np.array([Mf[p] for p in pts])
    _, V = np.linalg.eigh(mats)
    vecs = V[..., -1]
    vecs = vecs / np.linalg.norm(vecs, axis=-1, keepdims=True)
    signs, ncf = lift_signs(pts, vecs)
    Q = surface_degree(pts, signs[:, None] * vecs, ilo, ihi)
    return Q, ncf


def read_charge_from_vecs(pts, vecs, ilo, ihi):
    signs, ncf = lift_signs(pts, vecs)
    Q = surface_degree(pts, signs[:, None] * vecs, ilo, ihi)
    return Q, ncf


# ----------------------------------------------------- analytic seed math ----
def step_f(t):
    return np.arctan(t) / np.pi + 0.5


def ang_family(fam, x, y, s):
    if fam == "E":
        R = 2.0
        return -s * (np.arctan2(-x, R - y) - np.arctan2(-x, R))
    R = 1.0
    base = s * (np.arctan2(y - R, x) - np.arctan2(-R, x))
    if fam == "N":
        return base - np.pi * step_f(5 * x)
    if fam == "P":
        return (base - (2 * np.pi / 3) * step_f(5 * x)
                + (np.pi / 3) * (step_f(5 * (x - 1)) + step_f(5 * (x + 1))))
    raise ValueError(fam)


def analytic_dirs(pts, coords, fam, s, c):
    zhat = np.array([0.0, 0.0, 1.0])
    out = np.zeros((len(pts), 3))
    for k, p in enumerate(pts):
        X, Y, Z = coords[p[0]], coords[p[1]], coords[p[2]]
        rho = float(np.hypot(X, Y))
        a = ang_family(fam, Z / c, rho / c, s)
        rhohat = np.array([X / rho, Y / rho, 0.0])
        v = np.cos(a) * zhat + np.sin(a) * rhohat
        out[k] = v / np.linalg.norm(v)
    return out


def q2d_analytic(fam, s, d=10.0, npts=800001):
    """(1/pi) int_0^pi d(ang)/dphi dphi along (x,y) = (-d cos phi, d sin phi).

    d/dphi atan2(b, a) = (a b' - b a') / (a^2 + b^2): smooth, no unwrapping.
    """
    phi = np.linspace(0.0, np.pi, npts)
    x = -d * np.cos(phi)
    y = d * np.sin(phi)
    xp = d * np.sin(phi)
    yp = d * np.cos(phi)

    def datan2(b, a, bp, ap):
        return (a * bp - b * ap) / (a * a + b * b)

    def dstep(shift):
        return 5.0 * xp / (np.pi * (1.0 + 25.0 * (x - shift) ** 2))

    if fam == "E":
        R = 2.0
        d1 = datan2(-x, R - y, -xp, -yp)
        d2 = datan2(-x, np.full_like(x, R), -xp, np.zeros_like(x))
        dang = -s * (d1 - d2)
    else:
        R = 1.0
        d1 = datan2(y - R, x, yp, xp)
        d2 = datan2(np.full_like(x, -R), x, np.zeros_like(x), xp)
        dang = s * (d1 - d2)
        if fam == "N":
            dang = dang - np.pi * dstep(0.0)
        else:
            dang = dang - (2 * np.pi / 3) * dstep(0.0) \
                + (np.pi / 3) * (dstep(1.0) + dstep(-1.0))
    return float(TRAP(dang, x=phi) / np.pi)


# ------------------------------------------------------- stationarity probe --
def smooth_direction(Mf, shell=2):
    n = Mf.shape[0]
    S = np.zeros_like(Mf)
    for ax in range(3):
        S[_sl(ax, 1, n)] += Mf[_sl(ax, 0, n - 1)]
        S[_sl(ax, 0, n - 1)] += Mf[_sl(ax, 1, n)]
    dvec = S / 6.0 - Mf
    for ax in range(3):
        dvec[_sl(ax, 0, shell)] = 0.0
        dvec[_sl(ax, n - shell, n)] = 0.0
    return dvec


def stationarity_probe(Mf, h, delta, w2, conv, t=0.05):
    dvec = smooth_direction(Mf)
    dn = float(np.linalg.norm(dvec))
    e0 = energy(Mf, h, delta, w2, conv)["E_end"]
    ep = energy(Mf + t * dvec, h, delta, w2, conv)["E_end"]
    em = energy(Mf - t * dvec, h, delta, w2, conv)["E_end"]
    dEdt = (ep - em) / (2 * t)
    curv = (ep + em - 2 * e0) / (t * t)
    return {"E0": e0, "dE_dt": dEdt, "dE_dt_per_unit": dEdt / dn,
            "curvature": curv, "dir_norm": dn}


# -------------------------------------------------- analytic gradient (own) --
def grad_energy(Mf, h, delta, w2):
    """Analytic gradient of E_sym under the pad0 edge convention (own derivation).

    d(tr(C^T C))/dA_i = 2[C, A_j] with C = [A_i, A_j] antisymmetric, A symmetric;
    times the u-prefactor 4 -> 8[C, A_j]. Stencil adjoints for pad0 one-sided
    differences; V_T2 gradient = V diag(2 w2 (lambda - v)) V^T per cell.
    """
    n = Mf.shape[0]
    G = np.zeros_like(Mf)
    for side in ("fwd", "bwd"):
        A = oneside_A_list(Mf, h, side)
        gA = [np.zeros_like(Mf) for _ in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                gA[i] += 8.0 * (C @ A[j] - A[j] @ C)
                gA[j] -= 8.0 * (C @ A[i] - A[i] @ C)
        for ax in range(3):
            g = gA[ax]
            D = np.zeros_like(Mf)
            if side == "fwd":
                gg = g.copy()
                gg[_sl(ax, n - 1, n)] = 0.0
                D[_sl(ax, 1, n)] += gg[_sl(ax, 0, n - 1)]
                D -= gg
            else:
                gg = g.copy()
                gg[_sl(ax, 0, 1)] = 0.0
                D += gg
                D[_sl(ax, 0, n - 1)] -= gg[_sl(ax, 1, n)]
            G += D / h
    G *= 0.5
    lam, V = np.linalg.eigh(Mf)
    vtar = np.array(sorted((1.0, delta, 0.0)))
    S = 2.0 * w2 * (lam - vtar)
    G = G + np.einsum("...ik,...k,...jk->...ij", V, S, V)
    return G * h ** 3


def oneside_A_list(Mf, h, side):
    n = Mf.shape[0]
    out = []
    for ax in range(3):
        d = (Mf[_sl(ax, 1, n)] - Mf[_sl(ax, 0, n - 1)]) / h
        D = np.zeros_like(Mf)
        if side == "fwd":
            D[_sl(ax, 0, n - 1)] = d
        else:
            D[_sl(ax, 1, n)] = d
        out.append(D)
    return out


def grad_gate(rng):
    """FD gate for grad_energy on a small random symmetric field."""
    n = 8
    h = 1.0
    delta, w2 = 0.3, 0.0027581
    M = rng.standard_normal((n, n, n, 3, 3))
    M = 0.5 * (M + np.swapaxes(M, -1, -2)) * 0.4
    G = grad_energy(M, h, delta, w2)
    worst = 0.0
    for _ in range(4):
        D = rng.standard_normal((n, n, n, 3, 3))
        D = 0.5 * (D + np.swapaxes(D, -1, -2))
        t = 1e-6
        num = (energy(M + t * D, h, delta, w2, "pad0")["E_end"]
               - energy(M - t * D, h, delta, w2, "pad0")["E_end"]) / (2 * t)
        ana = float((G * D).sum())
        worst = max(worst, abs(num - ana) / max(abs(num), 1e-12))
    return worst


def descent_probe(Mf, h, delta, w2, steps=150, alpha=0.02, shell=2):
    """Own plain-GD descent with the analytic gradient, pinned shell held.

    A genuinely stationary endpoint sheds only the float32-rounding roughness
    (tiny dE, end slope ~ 0); a sliding state keeps producing energy drop.
    """
    n = Mf.shape[0]
    M = Mf.copy()
    e0 = energy(M, h, delta, w2, "pad0")["E_end"]
    ehist = [e0]
    fmax0 = None
    for k in range(steps):
        G = grad_energy(M, h, delta, w2)
        for ax in range(3):
            G[_sl(ax, 0, shell)] = 0.0
            G[_sl(ax, n - shell, n)] = 0.0
        if fmax0 is None:
            fmax0 = float(np.abs(G).max())
        M -= alpha * G
        if (k + 1) % 25 == 0:
            e = energy(M, h, delta, w2, "pad0")["E_end"]
            if e > ehist[-1]:
                alpha *= 0.5
            ehist.append(e)
    e1 = energy(M, h, delta, w2, "pad0")["E_end"]
    Gend = grad_energy(M, h, delta, w2)
    for ax in range(3):
        Gend[_sl(ax, 0, shell)] = 0.0
        Gend[_sl(ax, n - shell, n)] = 0.0
    fmax1 = float(np.abs(Gend).max())
    tail = ehist[-2] - e1 if len(ehist) >= 2 else 0.0
    return {"E_start": e0, "E_end_probe": e1, "dE_total": e0 - e1,
            "dE_last25": tail, "fmax_start": fmax0, "fmax_end": fmax1,
            "steps": steps, "alpha_final": alpha}


# ---------------------------------------------------------------- helpers ----
def rowjson(name):
    with open(DATA / name) as f:
        return json.load(f)


def coords_of(n, h):
    return (np.arange(n) - (n - 1) / 2.0) * h


def fmt(x, nd=6):
    return float(np.round(x, nd))


# ======================================================================= main
def main():
    audit = {}
    print("=== M5.22 adversarial audit (independent instruments) ===")

    # ---------------- C1: seed 2D windings, analytic-derivative integration --
    print("\n-- C1: seed-charge transcription (GATE A) --")
    targets = {
        ("N", -1.0): 1.01273, ("N", -0.5): 0.0127307, ("N", 0.0): -0.987269,
        ("N", 0.5): -1.98727, ("N", 1.0): -2.98727,
        ("P", -1.0): 1.99991, ("P", -0.5): 0.999914, ("P", 0.0): -0.0000857,
        ("P", 0.5): -1.00009, ("P", 1.0): -2.00009,
    }
    c1 = {}
    worst = 0.0
    for (fam, s), tgt in targets.items():
        q = q2d_analytic(fam, s)
        err = abs(q - tgt)
        worst = max(worst, err)
        c1[f"{fam}_s{s:g}"] = {"mine": q, "author": tgt, "abs_err": err}
        print(f"  {fam} s={s:+.1f}: mine {q:+.7f} author {tgt:+.7f} err {err:.2e}")
    for s in (-0.5, 0.5):
        q = q2d_analytic("E", s)
        c1[f"E_s{s:g}"] = {"mine": q, "expected": -np.sign(s)}
        print(f"  E s={s:+.1f}: mine {q:+.7f} (expected {-np.sign(s):+.0f})")
    c1_verdict = "CONFIRMED" if worst < 5e-5 else "REFUTED"
    audit["C1"] = {"verdict": c1_verdict, "worst_abs_err": worst, "detail": c1,
                   "method": "analytic d(ang)/dphi integrated (trapezoid, 8e5 pts)"}
    print(f"  C1 {c1_verdict}: worst |err| = {worst:.2e}")

    # ---------------- energy edge-convention scan on P-0.5 n48 ---------------
    print("\n-- energy edge-convention scan (P-0.5 n48) --")
    Mf, h, delta, _ = load_end("P-0.5_plane_sc6_n48_pinned_d0.3")
    row = rowjson("m5_22_row_P-0.5_plane_sc6_n48_pinned_d0.3.json")
    w2 = float(row["w2"])
    scan = {}
    for conv in ("pad0", "clamp", "crop", "wrap"):
        ef = eu_side(Mf, h, "fwd", conv)
        eb = eu_side(Mf, h, "bwd", conv)
        scan[conv] = {"E_u_fwd": ef, "E_u_bwd": eb,
                      "err_vs_row": abs(0.5 * (ef + eb) - row["E_u"])}
        print(f"  {conv:6s}: E_u_fwd {ef:.9f} E_u_bwd {eb:.9f} "
              f"|sym-row| {scan[conv]['err_vs_row']:.3e}")
    conv = min(scan, key=lambda k: scan[k]["err_vs_row"])
    ev = ev_of(Mf, h, delta, w2)
    print(f"  chosen convention: {conv}; E_v mine {ev:.9f} row {row['E_v']:.9f} "
          f"(diff {abs(ev - row['E_v']):.3e})")
    audit["instrument"] = {
        "edge_convention_scan": {k: {kk: fmt(vv, 9) for kk, vv in v.items()}
                                 for k, v in scan.items()},
        "chosen_convention": conv,
        "E_v_check": {"mine": ev, "row": row["E_v"]},
        "charge_method": "orientation lift (max-|dot| Prim) + exact quad solid "
                         "angles (van Oosterom-Strackee), NOT the note's FD flux",
        "caveats": [
            "Edge convention: the note's equations do not pin the one-sided "
            "stencils' boundary handling; pad0 (missing edge derivative = 0) is "
            "the implemented convention. clamp/crop shift E_u by 0.07-0.41 "
            "(~1-6%), so the convention is load-bearing and should be stated.",
            "Charge reads 0.995/1.002 are FD-flux discretization error, not "
            "physics: the solid-angle degree of the same lifted fields is "
            "EXACTLY -1.000000 at every half-width for both states. The note's "
            "'0.7% apart' (read vi) quantifies the instrument, not the states: "
            "proton and lepton far-field degrees are exactly equal integers.",
            "Row n_conflicts (2016/1918) is not a far-surface property: the "
            "far-cube lift is cleanly orientable (0 conflicts at w = 13.5, "
            "16.5, 19.5 with a max-|dot|-first lift).",
            "float32 endpoint storage sets a stationarity re-verification "
            "floor: a true f_tol state reads max|force| ~ 2e-7 from the npz "
            "(claimed 1e-8); smoothing-type probes bottom out at ~5e-4/unit.",
        ],
    }

    # ---------------- C2: energy recompute on 7 endpoints --------------------
    print("\n-- C2: energy values --")
    cases = [
        ("P-0.5_plane_sc6_n48_pinned_d0.3", 8.252),
        ("P-1_plane_sc6_n48_pinned_d0.3", 12.719),
        ("E-0.5_plane_sc6_n48_pinned_d0.3", 6.253),
        ("N-0.5_plane_sc6_n48_pinned_d0.3", 2.058),
        ("P-0.5_plane_sc6_n32_pinned_d0.3", 8.496),
        ("P-1_plane_sc6_n32_pinned_d0.3", 12.730),
        ("E+0.5_plane_sc6_n32_pinned_d0.3", 6.029),
    ]
    c2 = {}
    myE = {}
    worst_rel = 0.0
    for tag, noteval in cases:
        Mf, h, delta, _ = load_end(tag)
        r = rowjson(f"m5_22_row_{tag}.json")
        e = energy(Mf, h, delta, float(r["w2"]), conv)
        myE[tag] = e
        rel = abs(e["E_end"] - r["E_end"]) / abs(r["E_end"])
        worst_rel = max(worst_rel, rel)
        c2[tag] = {"mine": e["E_end"], "row": r["E_end"], "note": noteval,
                   "rel_err": rel}
        print(f"  {tag}: mine {e['E_end']:.6f} row {r['E_end']:.6f} "
              f"note {noteval} rel {rel:.2e}")
    c2_verdict = "CONFIRMED" if worst_rel < 1e-4 else (
        "PARTIAL" if worst_rel < 1e-2 else "REFUTED")
    audit["C2"] = {"verdict": c2_verdict, "worst_rel_err": worst_rel,
                   "edge_convention": conv, "detail": c2}
    print(f"  C2 {c2_verdict}: worst rel err {worst_rel:.2e}")

    # ---------------- C3: charge instrument, controls + endpoint reads -------
    print("\n-- C3: charge controls + endpoint charges --")
    n32c = coords_of(32, 1.5)
    ilo, ihi = 4, 27  # |coord| <= 17.25 on the n32 grid
    pts = surface_points(ilo, ihi)
    # (a) radial hedgehog, degree +1
    hh = np.zeros((len(pts), 3))
    for k, p in enumerate(pts):
        r = np.array([n32c[p[0]], n32c[p[1]], n32c[p[2]]])
        hh[k] = r / np.linalg.norm(r)
    q_h, cf_h = read_charge_from_vecs(pts, hh, ilo, ihi)
    # (b) degree-2: azimuth doubled on the position sphere
    d2 = np.zeros((len(pts), 3))
    for k, p in enumerate(pts):
        X, Y, Z = n32c[p[0]], n32c[p[1]], n32c[p[2]]
        r = np.sqrt(X * X + Y * Y + Z * Z)
        st = np.hypot(X, Y) / r
        ph = np.arctan2(Y, X)
        d2[k] = [st * np.cos(2 * ph), st * np.sin(2 * ph), Z / r]
    q_2, cf_2 = read_charge_from_vecs(pts, d2, ilo, ihi)
    print(f"  hedgehog control: Q = {q_h:+.6f} (conflicts {cf_h}); expect +1")
    print(f"  degree-2 control: Q = {q_2:+.6f} (conflicts {cf_2}); expect +2")

    ends = {}
    for tag, claimed in (("P-0.5_plane_sc6_n48_pinned_d0.3", 0.995),
                         ("E-0.5_plane_sc6_n48_pinned_d0.3", 1.002)):
        Mf, h, delta, _ = load_end(tag)
        prof = {}
        for w, (a, b) in (("13.5", (10, 37)), ("16.5", (7, 40)),
                          ("19.5", (4, 43))):
            Q, cf = read_charge_from_M(Mf, a, b)
            prof[w] = {"Q": Q, "conflicts": cf}
            print(f"  {tag} w={w}: Q = {Q:+.4f} (conflicts {cf})")
        ends[tag] = {"claimed_absQ": claimed, "profile": prof,
                     "mine_absQ_far": abs(prof["19.5"]["Q"])}
    ok_ctrl = abs(q_h - 1) < 0.02 and abs(q_2 - 2) < 0.02
    ok_ends = all(abs(v["mine_absQ_far"] - v["claimed_absQ"]) < 0.05
                  for v in ends.values())
    c3_verdict = "CONFIRMED" if (ok_ctrl and ok_ends) else (
        "PARTIAL" if ok_ctrl else "REFUTED")
    audit["C3"] = {"verdict": c3_verdict,
                   "hedgehog": {"Q": q_h, "conflicts": cf_h},
                   "degree2": {"Q": q_2, "conflicts": cf_2},
                   "endpoints": ends}
    print(f"  C3 {c3_verdict}")

    # ---------------- C4: the P-1 escape ------------------------------------
    print("\n-- C4: P-1 analytic director, 3D degree vs 2D winding --")
    n48c = coords_of(48, 1.0)
    c4 = {"q2d_P_s-1": audit["C1"]["detail"]["P_s-1"]["mine"]}
    for w, (a, b) in (("13.5", (10, 37)), ("16.5", (7, 40)), ("19.5", (4, 43))):
        pts48 = surface_points(a, b)
        vec = analytic_dirs(pts48, n48c, "P", -1.0, 6.0)
        Q, cf = read_charge_from_vecs(pts48, vec, a, b)
        c4[f"deg3d_w{w}"] = {"Q": Q, "conflicts": cf}
        print(f"  P s=-1 analytic, w={w}: Q3d = {Q:+.5f} (conflicts {cf})")
    pts48 = surface_points(4, 43)
    vec = analytic_dirs(pts48, n48c, "P", -0.5, 6.0)
    Qhalf, cfh = read_charge_from_vecs(pts48, vec, 4, 43)
    c4["deg3d_P_s-0.5_w19.5"] = {"Q": Qhalf, "conflicts": cfh}
    print(f"  cross-check P s=-1/2 analytic, w=19.5: Q3d = {Qhalf:+.5f}")
    esc = all(abs(c4[f"deg3d_w{w}"]["Q"]) < 0.1 for w in ("13.5", "16.5", "19.5"))
    wind2 = abs(c4["q2d_P_s-1"] - 2.0) < 1e-3
    c4_verdict = "CONFIRMED" if (esc and wind2 and abs(abs(Qhalf) - 1) < 0.1) \
        else "REFUTED"
    audit["C4"] = {"verdict": c4_verdict, "detail": c4,
                   "note": "2D winding +2 (1.99991) does not lift: 3D degree ~ 0; "
                           "the s=-1/2 seed lifts to |Q3d| ~ 1"}
    print(f"  C4 {c4_verdict}")

    # ---------------- C5: mass ordering -------------------------------------
    print("\n-- C5: mass ordering P-1 > P-0.5 --")
    o48 = myE["P-1_plane_sc6_n48_pinned_d0.3"]["E_end"] \
        - myE["P-0.5_plane_sc6_n48_pinned_d0.3"]["E_end"]
    o32 = myE["P-1_plane_sc6_n32_pinned_d0.3"]["E_end"] \
        - myE["P-0.5_plane_sc6_n32_pinned_d0.3"]["E_end"]
    c5_verdict = "CONFIRMED" if (o48 > 0 and o32 > 0) else "REFUTED"
    audit["C5"] = {"verdict": c5_verdict,
                   "E_P1_minus_EP05_n48": o48, "E_P1_minus_EP05_n32": o32}
    print(f"  n48: E(P-1)-E(P-0.5) = {o48:+.4f}; n32: {o32:+.4f} -> {c5_verdict}")

    # ---------------- C6: dissolution / stationarity ------------------------
    print("\n-- C6: stationarity (own analytic gradient + descent probe) --")
    rng = np.random.default_rng(5522)
    ggate = grad_gate(rng)
    print(f"  gradient FD gate (8^3 random field, 4 dirs): worst rel {ggate:.2e}")
    Mn, hn, dn, _ = load_end("N-0.5_plane_sc6_n48_pinned_d0.3")
    Mp, hp, dp, _ = load_end("P-1_plane_sc6_n48_pinned_d0.3")
    pn_s = stationarity_probe(Mn, hn, dn, w2, conv)
    pp_s = stationarity_probe(Mp, hp, dp, w2, conv)
    print(f"  smoothing probe (float32-floor-limited): N {pn_s['dE_dt_per_unit']:+.2e}"
          f" vs P {pp_s['dE_dt_per_unit']:+.2e} per unit: NOT discriminating")
    pn = descent_probe(Mn, hn, dn, w2)
    pp = descent_probe(Mp, hp, dp, w2)
    print(f"  N-0.5 n48 descent: dE_total {pn['dE_total']:+.3e} "
          f"last25 {pn['dE_last25']:+.3e} fmax {pn['fmax_start']:.2e}"
          f" -> {pn['fmax_end']:.2e}")
    print(f"  P-1   n48 descent: dE_total {pp['dE_total']:+.3e} "
          f"last25 {pp['dE_last25']:+.3e} fmax {pp['fmax_start']:.2e}"
          f" -> {pp['fmax_end']:.2e}")
    Mn32, h32, d32, _ = load_end("N-0.5_plane_sc6_n32_pinned_d0.3")
    e_n48 = myE["N-0.5_plane_sc6_n48_pinned_d0.3"]["E_end"]
    e_n32 = energy(Mn32, h32, d32, w2, conv)["E_end"]
    r48 = rowjson("m5_22_row_N-0.5_plane_sc6_n48_pinned_d0.3.json")
    r32 = rowjson("m5_22_row_N-0.5_plane_sc6_n32_pinned_d0.3.json")
    sliding = e_n48 < e_n32
    ratio = pn["dE_total"] / max(pp["dE_total"], 1e-300)
    c6_verdict = "CONFIRMED" if (ggate < 1e-6 and pn["dE_total"] > 10 * pp["dE_total"]
                                 and pn["dE_last25"] > 10 * abs(pp["dE_last25"])
                                 and sliding) else "PARTIAL"
    audit["C6"] = {"verdict": c6_verdict, "grad_fd_gate": ggate,
                   "descent_N05_n48": pn, "descent_P1_n48": pp,
                   "dE_ratio_N_over_P": ratio,
                   "smoothing_probe_N": pn_s, "smoothing_probe_P": pp_s,
                   "smoothing_probe_caveat": "float32 storage floor ~5e-4/unit "
                   "dominates both states; not discriminating",
                   "E_N05_n48_mine": e_n48, "E_N05_n32_mine": e_n32,
                   "E_N05_n48_row": r48["E_end"], "E_N05_n32_row": r32["E_end"],
                   "still_sliding": bool(sliding)}
    print(f"  E(N-0.5): n48 {e_n48:.4f} < n32 {e_n32:.4f} : {sliding}; "
          f"descent dE ratio N/P = {ratio:.1f}x -> {c6_verdict}")

    # ---------------- C7: E-family below the certified 2b state --------------
    print("\n-- C7: lepton endpoint vs certified 2b value --")
    r2b = rowjson("m5_21_2b_row_c48_R_T2.json")
    allj = rowjson("m5_21_2b_all.json")
    e2b_all = allj.get("c48_R_T2", {}).get("E_end")
    mineE = myE["E-0.5_plane_sc6_n48_pinned_d0.3"]["E_end"]
    below = (r2b["E_end"] - mineE) / r2b["E_end"]
    w2match = abs(float(r2b["w2"]) - 0.0027581) < 1e-9
    c7_verdict = "CONFIRMED" if (mineE < r2b["E_end"] and w2match
                                 and abs(below - 0.086) < 0.01) else "PARTIAL"
    audit["C7"] = {"verdict": c7_verdict, "mine_E_lepton_n48": mineE,
                   "certified_2b_row": r2b["E_end"], "certified_2b_all": e2b_all,
                   "w2_2b": float(r2b["w2"]), "fraction_below": below}
    print(f"  mine {mineE:.4f} vs 2b {r2b['E_end']:.4f} "
          f"(all.json {e2b_all}); {100*below:.2f}% below -> {c7_verdict}")

    # ---------------- C8: the Q38 scan refit ---------------------------------
    print("\n-- C8: Q38 quark-shift scan refit --")
    dxs = [0.0, 0.25, 0.5, 0.75, 1.0]
    names = ["dx0", "dx0.25", "dx0.5", "dx0.75", "dx1.0"]
    Ejson, Emine = [], []
    for nm in names:
        r = rowjson(f"m5_22_row_q38_{nm}.json")
        Ejson.append(float(r["E_end"]))
        Mq, hq, dq, _ = load_end(f"q38_{nm}")
        Emine.append(energy(Mq, hq, dq, float(r["w2"]), conv)["E_end"])
    sj = np.polyfit(dxs, Ejson, 1)[0]
    sm = np.polyfit(dxs, Emine, 1)[0]
    lat_m = sm / 6.0
    ratio_lo, ratio_hi = 6.2 / lat_m, 7.0 / lat_m
    print(f"  E(dx) mine: {[f'{e:.4f}' for e in Emine]}")
    print(f"  slope mine {sm:.4f}/author-unit = {lat_m:.5f}/lattice "
          f"(json fit {sj:.4f}); M5.21.4 tension 6.2-7.0 -> ratio "
          f"{ratio_lo:.0f}-{ratio_hi:.0f}x")
    c8_verdict = "CONFIRMED" if (abs(sm - 0.738) < 0.02 and 40 < ratio_lo < 60) \
        else "PARTIAL"
    audit["C8"] = {"verdict": c8_verdict, "E_mine": Emine, "E_json": Ejson,
                   "slope_mine_author_unit": sm, "slope_mine_lattice": lat_m,
                   "slope_json": sj, "ratio_vs_m5214_tension":
                   [ratio_lo, ratio_hi]}
    print(f"  C8 {c8_verdict}")

    # ---------------- C9: lepton mirror degeneracy ---------------------------
    print("\n-- C9: E+0.5 vs E-0.5 n32 mirror relation --")
    Ma, ha, da, za = load_end("E+0.5_plane_sc6_n32_pinned_d0.3")
    Mb, hb, db, zb = load_end("E-0.5_plane_sc6_n32_pinned_d0.3")
    M0a = np.asarray(za["M0"], dtype=np.float64)
    M0b = np.asarray(zb["M0"], dtype=np.float64)
    ea = myE["E+0.5_plane_sc6_n32_pinned_d0.3"]["E_end"]
    eb2 = energy(Mb, hb, db, w2, conv)["E_end"]
    dE_rel = abs(ea - eb2) / ea

    def scan_transforms(A, B):
        best = None
        na = float(np.linalg.norm(A))
        for flip in (None, 0, 1, 2):
            T0 = B if flip is None else np.flip(B, axis=flip)
            for conj in (None, 0, 1, 2):
                T = T0
                if conj is not None:
                    sgn = np.ones(3)
                    sgn[conj] = -1.0
                    T = T * (sgn[:, None] * sgn[None, :])
                rel = float(np.linalg.norm(A - T)) / na
                if best is None or rel < best[2]:
                    best = (flip, conj, rel)
        ident = float(np.linalg.norm(A - B)) / na
        return best, ident

    bseed, iseed = scan_transforms(M0a, M0b)
    bend, iend = scan_transforms(Ma, Mb)
    print(f"  seeds: best transform flip={bseed[0]} conj={bseed[1]} "
          f"rel {bseed[2]:.3e} (identity rel {iseed:.3e})")
    print(f"  ends : best transform flip={bend[0]} conj={bend[1]} "
          f"rel {bend[2]:.3e} (identity rel {iend:.3e})")
    print(f"  E equality: {ea:.8f} vs {eb2:.8f} (rel {dE_rel:.2e})")
    c9_verdict = "CONFIRMED" if (dE_rel < 1e-5 and bend[2] < 1e-3) else (
        "PARTIAL" if dE_rel < 1e-5 else "REFUTED")
    audit["C9"] = {"verdict": c9_verdict, "E_plus": ea, "E_minus": eb2,
                   "E_rel_diff": dE_rel,
                   "seed_best": {"flip_axis": bseed[0], "conj_axis": bseed[1],
                                 "rel_dist": bseed[2], "identity_rel": iseed},
                   "end_best": {"flip_axis": bend[0], "conj_axis": bend[1],
                                "rel_dist": bend[2], "identity_rel": iend}}
    print(f"  C9 {c9_verdict}")

    # ---------------- save ----------------------------------------------------
    def tofloat(o):
        if isinstance(o, dict):
            return {k: tofloat(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [tofloat(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        return o

    out = DATA / "m5_22_audit.json"
    with open(out, "w") as f:
        json.dump(tofloat(audit), f, indent=1)
    print(f"\nverdicts -> {out}")
    print("\n=== verdict summary ===")
    for k in ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"):
        print(f"  {k}: {audit[k]['verdict']}")


if __name__ == "__main__":
    main()
