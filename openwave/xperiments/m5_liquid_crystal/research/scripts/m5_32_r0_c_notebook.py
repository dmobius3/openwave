"""M5.32 R0 arm (c): sympy reproduction of the author's 2026-08-17 notebook
"newton for boost hedgehogs.nb"
(../../theory/duda_2026-08-17_newton_for_boost_hedgehogs.pdf, 2 pages).

THE NOTEBOOK, TRANSCRIBED (index conventions: internal indices 1..4 in the
notebook = 0..3 here, index 0 = time; the field is the symmetric 4x4 M):

  eta = diag(-1, 1, 1, 1)                     (defined, NEVER used below)
  M0  = diag(g, 0, 0, 0)                      ("lengths of axes")
  Gamma_mu = [[0,    tG1,  tG2,  tG3],
              [tG1,  0,   -G3,   G2],
              [tG2,  G3,   0,   -G1],
              [tG3, -G2,   G1,   0]]           (3D rot + boost)
  G4 = coefficients of Gamma_mu wrt (G1, G2, G3, tG1, tG2, tG3)
  Gb = G4[4..6]  = the three boost generators: (Gb_j)_{0j} = (Gb_j)_{j0} = 1
  f[r_] := 1/Sqrt[r]                           ("distance dependence to check")
  o1 = Normal[Series[MatrixExp[m f[x^2+y^2+(z+d)^2] ((x,y,z+d).Gb)], {m,0,1}]]
  o2 = Normal[Series[MatrixExp[m f[x^2+y^2+(z-d)^2] ((x,y,z-d).Gb)], {m,0,1}]]
     NOTE: f is applied to the SQUARED distance, so the effective profile is
     f(r^2) = 1/r, i.e. each hedgehog is m * (unit vector) . Gb.  A literal
     1/sqrt(r) profile would make the energy density fall only as 1/r and the
     half-plane integral over (0, inf)^2 diverge linearly; the notebook's
     finite numbers are only consistent with the 1/r reading.
  o  = Normal[Series[o1.o2, {m,0,1}]]
  M  = Normal[Series[o.M0.Transpose[o], {m,0,1}]]
  dM = Normal[Series[{D[M,x], D[M,y], D[M,z]}, {m,0,1}]]
  coms = { dM_i.dM_j - dM_j.dM_i : (i,j) in {(1,2),(2,3),(3,1)} }   (PLAIN
     matrix commutator; no eta, no xi inserted)
  Hs = Normal[Series[Sum_com (com_23^2 + com_34^2 + com_42^2), {m,0,4}]] /. y->0
     (1-based entries of the spatial block; "space-space curvature")
  en[d] = NIntegrate[4 Pi x Hs Boole[x^2+(z-d)^2 > 0.001] / m^4 / g^4,
                     {x,0,Inf}, {z,0,Inf}],  d = 0.1, 0.2, ..., 3.0
     (half-plane (x,z) with the axial weight 2 pi x, doubled by the z -> -z
     symmetry; the cutoff is on the SQUARED distance from the +d center, the
     only center inside z > 0; r_c = sqrt(0.001) = 0.0316)
  Fit[en, {1, 1/d}, d]  ->  863.733 + 167.668/d      (B > 0: like charges)

Reproduction here:
  * sympy builds Gb, o1, o2, o, M, dM exactly as above (series truncation by
    polynomial coefficient extraction in m; dM is EXACTLY linear in m at
    first order, so com is exactly m^2 and Hs exactly m^4: the notebook's
    Series[..., {m,0,4}] is the identity there).
  * The m-coefficient matrices of dM (small expressions) are lambdified;
    the commutators, the squares and the m-order bookkeeping are done as
    vectorized numpy 4x4 algebra at the quadrature nodes (the fully expanded
    symbolic Hs is a 40 kB rational expression whose 1/r^4 tail arises from
    cancellation; the matrix route has none).
  * Hs(y = 0) is checked against the notebook's printed closed form Out[13]
    at random points (relative difference reported).
  * Half-plane integral: polar chart about the +d center,
    x = s sin(th), z = d + s cos(th), th in (0, pi), s in (r_c, s_max(th)),
    s_max = S for cos(th) >= 0 and min(S, d/(-cos(th))) otherwise (the
    z > 0 wall); S = 1e7 stands for the notebook's infinity (the tail
    beyond S is O(1/S)).  Gauss-Legendre product rule in (ln s, th), the th
    range split at pi/2, node counts doubled until the energies agree to
    1e-7 relative (the convergence table is written to the JSON).
  * Fit = ordinary least squares on {1, 1/d} over the notebook's 30 d values.
  * Sensitivity ladder: cutoffs r_c^2 in {1e-4, 1e-3, 1e-2} x domains
    {S = 1e7 ("inf"), S = 20}.

Pre-registered extensions (separate functions, not the notebook):
  (i)   FULL eta-contraction: F_ij = dM_i xi dM_j - dM_j xi dM_i,
        <F,F>_eta = tr(eta F eta F^T), summed over i<j (static: F_0i = 0),
        xi = eta = diag(-1,1,1,1); and the plain-commutator variant.
  (ii)  full vacuum spectrum M0 = diag(g, 1, delta, 0), g = 32, delta = 0.3.
  (iii) second order in m (o1, o2, o, M, dM to O(m^2)); H then carries m^4,
        m^5 and (incomplete without O(m^3) in M) m^6 pieces, fitted per order.

Out: ../data/m5_32_r0_notebook.json, ../plots/m5_32_r0_notebook_fit.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

NB_A, NB_B = 863.733, 167.668
NB_D = [round(0.1 * k, 10) for k in range(1, 31)]
NB_CUT2 = 0.001
S_INF = 1.0e7

x, y, z, d = sp.symbols("x y z d", positive=True)
m, g, delta = sp.symbols("m g delta", positive=True)
XI = np.diag([-1.0, 1.0, 1.0, 1.0])


# ----------------------------------------------------------------- sympy side
def boost_generators():
    Gb = []
    for j in (1, 2, 3):
        G = sp.zeros(4, 4)
        G[0, j] = 1
        G[j, 0] = 1
        Gb.append(G)
    return Gb


def trunc(expr, order):
    """keep the terms of a polynomial-in-m expression through m**order."""
    e = sp.expand(expr)
    return sum(e.coeff(m, k) * m**k for k in range(order + 1))


def trunc_mat(Mx, order):
    return Mx.applyfunc(lambda e: trunc(e, order))


def hedgehog(sign, order):
    """o_k = exp(m f(r^2) (x, y, z + sign d) . Gb) through m**order, f = 1/sqrt."""
    Gb = boost_generators()
    v = (x, y, z + sign * d)
    r2 = x**2 + y**2 + (z + sign * d) ** 2
    X = m * (1 / sp.sqrt(r2)) * sum((vi * Gi for vi, Gi in zip(v, Gb)),
                                    sp.zeros(4, 4))
    o = sp.eye(4)
    P = sp.eye(4)
    for k in range(1, order + 1):
        P = P * X
        o = o + P / sp.factorial(k)
    return o


def build_dM(M0, order=1):
    o1 = hedgehog(+1, order)
    o2 = hedgehog(-1, order)
    o = trunc_mat(o1 * o2, order)
    M = trunc_mat(o * M0 * o.T, order)
    return [trunc_mat(M.diff(v), order) for v in (x, y, z)]


class DMNumeric:
    """dM_i = sum_k m^k D_i^(k); each D^(k) lambdified at y = 0."""

    def __init__(self, M0, order, subs=None):
        dM = build_dM(M0, order)
        subs = subs or {}
        self.order = order
        self.fns = {}
        for i in range(3):
            for k in range(1, order + 1):
                Dk = dM[i].applyfunc(
                    lambda e: sp.expand(e).coeff(m, k).subs(y, 0).subs(subs))
                self.fns[(i, k)] = [[sp.lambdify((x, z, d), Dk[a, b], "numpy")
                                     for b in range(4)] for a in range(4)]

    def __call__(self, xx, zz, dd):
        """returns D[i][k] as arrays of shape (4, 4, N)."""
        N = xx.shape[0]
        out = {}
        for (i, k), rows in self.fns.items():
            A = np.zeros((4, 4, N))
            for a in range(4):
                for b in range(4):
                    A[a, b] = np.broadcast_to(rows[a][b](xx, zz, dd), (N,))
            out[(i, k)] = A
        return out


def notebook_out13():
    """the printed Out[13] (Hs at y = 0), for the closed-form check."""
    rm = sp.sqrt(d**2 + x**2 - 2 * d * z + z**2)
    rp = sp.sqrt(d**2 + x**2 + 2 * d * z + z**2)
    num = (d**8
           + d**6 * (2 * x**2 - 2 * z**2 + rm * rp)
           + (x**2 + z**2) ** 3 * (x**2 + z**2 + rm * rp)
           + d**2 * (x**4 - z**4) * (2 * x**2 + 2 * z**2 + rm * rp)
           + d**4 * (4 * x**4 + 2 * z**4 - z**2 * rm * rp
                     + x**2 * (4 * z**2 + rm * rp)))
    return 8 * g**4 * m**4 * num / (rm**6 * rp**6)


# ---------------------------------------------------------- numpy curvature
def mm(A, B):
    return np.einsum("abn,bcn->acn", A, B)


def com(A, B, kind):
    if kind == "plain":
        return mm(A, B) - mm(B, A)
    e = np.diag(XI)[None, :, None]          # A xi B = (A with columns scaled) B
    return mm(A * e, B) - mm(B * e, A)


def h_spatial(F):
    """the notebook Hs: com_23^2 + com_34^2 + com_42^2 (1-based)."""
    return F[1, 2] ** 2 + F[2, 3] ** 2 + F[3, 1] ** 2


def h_eta(F):
    """extension (i): tr(eta F eta F^T) = sum_ab eta_a eta_b F_ab^2."""
    e = np.diag(XI)
    return np.einsum("a,b,abn->n", e, e, F**2)


PAIRS = [(0, 1), (1, 2), (2, 0)]


def density_by_order(D, kind, hfun, order):
    """H = sum_k m^k H_k from dM_i = sum_k m^k D_i^(k); returns {k: H_k}."""
    # com_ij = sum_{k,l} m^(k+l) [D_i^(k), D_j^(l)]_kind
    comorders = {}
    for i, j in PAIRS:
        for k in range(1, order + 1):
            for l in range(1, order + 1):
                c = com(D[(i, k)], D[(j, l)], kind)
                comorders.setdefault((i, j), {}).setdefault(k + l, []).append(c)
    Hk = {}
    for ij, byorder in comorders.items():
        summed = {p: sum(v) for p, v in byorder.items()}
        for p, Fp in summed.items():
            for q, Fq in summed.items():
                Hk[p + q] = Hk.get(p + q, 0) + hprod(Fp, Fq, hfun)
    return Hk


def hprod(F, G, hfun):
    """bilinear form behind hfun: h(F+G) = h(F) + h(G) + 2 hprod(F, G)."""
    return 0.25 * (hfun(F + G) - hfun(F - G))


# ---------------------------------------------------------------- quadrature
def gl(n, a, b):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * t + 0.5 * (b + a), 0.5 * (b - a) * w


def nodes(dd, cut2, S, ns, nth):
    """polar product rule; returns xx, zz, weight (weight includes 4 pi x s ds dth)."""
    rc = np.sqrt(cut2)
    xs, zs, ws = [], [], []
    for a, b in ((0.0, np.pi / 2), (np.pi / 2, np.pi)):
        th, wth = gl(nth, a, b)
        for t, wt in zip(th, wth):
            c = np.cos(t)
            smax = S if c >= 0 else min(S, dd / (-c))
            if smax <= rc:
                continue
            u, wu = gl(ns, np.log(rc), np.log(smax))
            s = np.exp(u)
            xs.append(s * np.sin(t))
            zs.append(dd + s * c)
            ws.append(4.0 * np.pi * (s * np.sin(t)) * s * s * wu * wt)
    return np.concatenate(xs), np.concatenate(zs), np.concatenate(ws)


def energy_orders(dmn, kind, hfun, dd, cut2, S, ns, nth, norm):
    xx, zz, ww = nodes(dd, cut2, S, ns, nth)
    D = dmn(xx, zz, dd)
    Hk = density_by_order(D, kind, hfun, dmn.order)
    return {k: float(np.dot(ww, Hk[k]) / norm) for k in Hk}


def fit_1_over_d(ds, es):
    ds = np.asarray(ds, float)
    es = np.asarray(es, float)
    Amat = np.column_stack([np.ones_like(ds), 1.0 / ds])
    coef, *_ = np.linalg.lstsq(Amat, es, rcond=None)
    res = es - Amat @ coef
    return float(coef[0]), float(coef[1]), res


def scan(dmn, kind, hfun, ds, cut2, S, ns, nth, norm, orders):
    t0 = time.time()
    per = {k: [] for k in orders}
    for dd in ds:
        e = energy_orders(dmn, kind, hfun, dd, cut2, S, ns, nth, norm)
        for k in orders:
            per[k].append(e.get(k, 0.0))
    rows = {}
    for k in orders:
        A, B, res = fit_1_over_d(ds, per[k])
        rows[k] = {"A": A, "B": B, "sign_B": "+" if B > 0 else ("-" if B < 0 else "0"),
                   "E": per[k], "residuals": res.tolist(),
                   "rms_residual": float(np.sqrt(np.mean(res**2)))}
    return rows, time.time() - t0


def compare(row):
    row["A_rel_diff_pct"] = 100.0 * (row["A"] - NB_A) / NB_A
    row["B_rel_diff_pct"] = 100.0 * (row["B"] - NB_B) / NB_B
    row["gate_1pct_and_sign"] = bool(abs(row["A_rel_diff_pct"]) < 1
                                     and abs(row["B_rel_diff_pct"]) < 1
                                     and row["B"] > 0)
    return row


# ------------------------------------------------------------------ the runs
def run_notebook(out):
    t0 = time.time()
    dmn = DMNumeric(sp.diag(g, 0, 0, 0), 1, {g: 1})
    out["symbolic_runtime_s"] = time.time() - t0

    # closed-form check vs Out[13] (g = 1, m = 1; the density is m^4 g^4 x ...)
    ref = sp.lambdify((x, z, d), notebook_out13().subs({g: 1, m: 1}), "numpy")
    rng = np.random.default_rng(0)
    P = rng.uniform(0.2, 3.0, size=(400, 3))
    D = dmn(P[:, 0], P[:, 1], P[:, 2])
    H4 = density_by_order(D, "plain", h_spatial, 1)[4]
    R = ref(P[:, 0], P[:, 1], P[:, 2])
    rel = float(np.max(np.abs(H4 - R) / np.abs(R)))
    out["out13_closed_form_max_rel_diff"] = rel
    print(f"[notebook] Hs vs printed Out[13]: max rel diff {rel:.2e} over 400 pts")

    # quadrature convergence at the notebook configuration
    conv = []
    prev = None
    for ns, nth in ((100, 40), (200, 80), (400, 160), (800, 320)):
        rows, dt = scan(dmn, "plain", h_spatial, [0.3, 1.0, 3.0], NB_CUT2, S_INF,
                        ns, nth, 1.0, [4])
        E = rows[4]["E"]
        entry = {"ns": ns, "nth": nth, "E_d0.3_1_3": E, "runtime_s": dt}
        if prev is not None:
            entry["max_rel_change"] = float(max(abs(a - b) / abs(b)
                                                for a, b in zip(E, prev)))
        conv.append(entry)
        print(f"[conv] ns={ns} nth={nth} E={E} "
              f"chg={entry.get('max_rel_change', float('nan')):.2e} ({dt:.1f}s)")
        prev = E
    out["quadrature_convergence"] = conv
    NS, NTH = 400, 160

    rows, dt = scan(dmn, "plain", h_spatial, NB_D, NB_CUT2, S_INF, NS, NTH, 1.0, [4])
    row = compare(rows[4])
    row.update({"construction": "notebook: Hs spatial block, plain commutator, "
                                "M0=diag(g,0,0,0), O(m^1), cut2=0.001, S=inf",
                "d_values": NB_D, "cut2": NB_CUT2, "S": "inf", "runtime_s": dt,
                "ns": NS, "nth": NTH})
    out["notebook_reproduction"] = row
    print(f"[notebook] A = {row['A']:.3f} ({row['A_rel_diff_pct']:+.3f}%)  "
          f"B = {row['B']:.3f} ({row['B_rel_diff_pct']:+.3f}%)  "
          f"rms res {row['rms_residual']:.2f}  gate={row['gate_1pct_and_sign']}  "
          f"{dt:.0f}s")

    ladder = []
    for cut2 in (1e-4, 1e-3, 1e-2):
        for S in (S_INF, 20.0):
            rows, dt = scan(dmn, "plain", h_spatial, NB_D, cut2, S, NS, NTH, 1.0, [4])
            r = {k: rows[4][k] for k in ("A", "B", "sign_B", "rms_residual")}
            r.update({"cut2": cut2, "S": "inf" if S == S_INF else S, "runtime_s": dt})
            ladder.append(r)
            print(f"[ladder] cut2={cut2:g} S={r['S']}: A={r['A']:.3f} "
                  f"B={r['B']:.3f} rms={r['rms_residual']:.2f} ({dt:.0f}s)")
    out["sensitivity_ladder"] = ladder
    return NS, NTH


def ext_row(name, dmn, kind, hfun, norm, orders, rows_out, NS, NTH, note=""):
    rows, dt = scan(dmn, kind, hfun, NB_D, NB_CUT2, S_INF, NS, NTH, norm, orders)
    for k in orders:
        r = rows[k]
        r = {kk: r[kk] for kk in ("A", "B", "sign_B", "rms_residual", "E")}
        r["construction"] = f"{name}, m^{k} coefficient{note if k == 6 else ''}"
        r["runtime_s"] = dt / len(orders)
        if abs(r["B"]) < 1e-9 and abs(r["A"]) < 1e-9:
            r["verdict"] = "vanishes identically"
        elif r["B"] > 0:
            r["verdict"] = ("B > 0: energy falls with separation, the like-charge "
                            "(wrong Newton) sign persists")
        else:
            r["verdict"] = "B < 0: energy rises with separation, attractive"
        rows_out.append(r)
        print(f"[ext] {r['construction']}: A={r['A']:.6g} B={r['B']:.6g} "
              f"sign {r['sign_B']} rms={r['rms_residual']:.3g} ({dt:.0f}s)")


def run_extensions(out, NS, NTH):
    rows = []
    # (i) full eta contraction, M0 = diag(g,0,0,0), first order
    dmn = DMNumeric(sp.diag(g, 0, 0, 0), 1, {g: 1})
    ext_row("(i) <F,F>_eta, xi-commutator, diag(g,0,0,0), O(m^1)",
            dmn, "xi", h_eta, 1.0, [4], rows, NS, NTH)
    ext_row("(i') <F,F>_eta, plain commutator, diag(g,0,0,0), O(m^1)",
            dmn, "plain", h_eta, 1.0, [4], rows, NS, NTH)
    # diagnostic: time-row content of F (xi-commutator) at random points
    rng = np.random.default_rng(1)
    P = rng.uniform(0.2, 3.0, size=(50, 3))
    D = dmn(P[:, 0], P[:, 1], P[:, 2])
    trow = max(float(np.max(np.abs(com(D[(i, 1)], D[(j, 1)], "xi")[0])))
               for i, j in PAIRS)
    out["ext_i_time_row_of_F_max_abs"] = trow

    # (ii) full vacuum spectrum g = 32, delta = 0.3 (normalized by g^4)
    gv, dv = 32.0, 0.3
    dmn2 = DMNumeric(sp.diag(g, 1, delta, 0), 1, {g: gv, delta: dv})
    ext_row("(ii) Hs spatial, plain commutator, diag(32,1,0.3,0), O(m^1)",
            dmn2, "plain", h_spatial, gv**4, [4], rows, NS, NTH)
    ext_row("(ii) <F,F>_eta, xi-commutator, diag(32,1,0.3,0), O(m^1)",
            dmn2, "xi", h_eta, gv**4, [4], rows, NS, NTH)
    D2 = dmn2(P[:, 0], P[:, 1], P[:, 2])
    out["ext_ii_time_row_of_F_max_abs_over_g4"] = max(
        float(np.max(np.abs(com(D2[(i, 1)], D2[(j, 1)], "xi")[0]))) / gv**4
        for i, j in PAIRS)

    # (iii) second order in m, diag(g,0,0,0) and diag(32,1,0.3,0)
    t0 = time.time()
    dmn3 = DMNumeric(sp.diag(g, 0, 0, 0), 2, {g: 1})
    dmn4 = DMNumeric(sp.diag(g, 1, delta, 0), 2, {g: gv, delta: dv})
    out["ext_iii_symbolic_runtime_s"] = time.time() - t0
    note = " [incomplete: needs O(m^3) in M]"
    ext_row("(iii) Hs spatial, plain comm, diag(g,0,0,0), O(m^2) in M",
            dmn3, "plain", h_spatial, 1.0, [4, 5, 6], rows, NS, NTH, note)
    ext_row("(iii) <F,F>_eta, xi-comm, diag(g,0,0,0), O(m^2) in M",
            dmn3, "xi", h_eta, 1.0, [4, 5, 6], rows, NS, NTH, note)
    ext_row("(iii)+(ii) <F,F>_eta, xi-comm, diag(32,1,0.3,0), O(m^2) in M",
            dmn4, "xi", h_eta, gv**4, [4, 5, 6], rows, NS, NTH, note)
    D4 = dmn4(P[:, 0], P[:, 1], P[:, 2])
    out["ext_iii_ii_time_row_of_F_m3_max_abs_over_g4"] = max(
        float(np.max(np.abs((com(D4[(i, 1)], D4[(j, 2)], "xi")
                             + com(D4[(i, 2)], D4[(j, 1)], "xi"))[0]))) / gv**4
        for i, j in PAIRS)
    out["extensions"] = rows


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    row = out["notebook_reproduction"]
    ds = np.array(row["d_values"])
    es = np.array(row["E"])
    dd = np.linspace(0.2, 3.0, 300)
    fig, ax = plt.subplots(figsize=(7, 4.4), dpi=130)
    ax.plot(dd, row["A"] + row["B"] / dd, color="#2563EB", lw=2,
            label=f"this fit: {row['A']:.1f} + {row['B']:.1f}/d")
    ax.plot(dd, NB_A + NB_B / dd, color="#D97706", lw=2, ls="--",
            label=f"notebook fit: {NB_A} + {NB_B}/d")
    ax.plot(ds, es, "o", ms=5, color="#374151", mfc="white", mew=1.5,
            label="E(d), this quadrature")
    ax.set_xlim(0.2, 3.05)
    ax.set_ylim(880, 1320)
    ax.set_xlabel("d (half separation; centers at z = ±d)")
    ax.set_ylabel("E(d) = ∫ 4π x Hs / (m⁴ g⁴)")
    ax.set_title("M5.32 R0 (c): the 2026-08-17 notebook reproduced")
    ax.grid(color="#E5E7EB", lw=0.8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    path = os.path.join(PLOTS, "m5_32_r0_notebook_fit.png")
    fig.savefig(path)
    print("[plot]", path)


def main():
    T0 = time.time()
    out = {"notebook_target": {"A": NB_A, "B": NB_B, "d_values": NB_D,
                               "cut2": NB_CUT2},
           "quadrature": "Gauss-Legendre product rule in (ln s, theta) on the "
                         "polar chart about the +d center, theta split at pi/2, "
                         "s in (sqrt(cut2), s_max(theta)), S_inf = 1e7"}
    NS, NTH = run_notebook(out)
    run_extensions(out, NS, NTH)
    out["total_runtime_s"] = time.time() - T0
    os.makedirs(DATA, exist_ok=True)
    os.makedirs(PLOTS, exist_ok=True)
    path = os.path.join(DATA, "m5_32_r0_notebook.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("[json]", path)
    plot(out)
    print(f"[done] {out['total_runtime_s']:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
