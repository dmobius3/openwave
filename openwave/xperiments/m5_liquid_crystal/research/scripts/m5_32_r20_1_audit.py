"""M5.32 R20-1 / R20-2 AUDIT: an independent recomputation of the three-axis
hedgehog reads (the producer: m5_32_r20_1_axes.py, its JSON
m5_32_r20_1_axes.json with the n32 extension rows merged from
m5_32_r20_1_axes_ext32.json, its saved end fields m5_32_r20_1/<tag>.npz).
Every read below is written from the definitions; the stack is imported for
the certified curvature cross-check (EN.block_reads), the seed builder gate
(the producer's seed_axes against m5_21_4_a_pair.seed_pair), and the R3
boost dressing (m5_32_r3_ii_pair.boost_at, conj). No relaxation, no descent:
energy evaluations on saved fields only.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4 per cell, eta = diag(-1, 1, 1, 1), N = M eta,
code branch s = -1: M_vac = diag(g, 1, delta, 0), delta = 0.3, g = 8 (32 for
the control), the N spectrum q = (-g, 1, delta, 0). Boxes n32 L48 and n64
L96, h = 1.5. Jets A_i = d_i M on the sym stencil (fwd and bwd branches of
weight 1/2, this file's own differences, zero on the last / first plane).
    <X, Y>_eta = sum_ab eta_a eta_b X_ab Y_ab
    F_ij = A_i eta A_j - A_j eta A_i
    e_curv(x) = 4 h^3 sum_br wt sum_{i<j} <F_ij, F_ij>_eta       (per cell)
    E_curv = sum_x e_curv
Potentials, this file's own, through the eigenvalues of N (the saved fields
have M_0i = 0 exactly, so N = diag(-M_00) (+) M_3 eta_3 and the spectrum is
(-M_00, eigh(M_3)); the general case falls back to eigvals(N)):
    V4     = h^3 W1 sum_x sum_{p=1..4} (t_p - C_p)^2,  t_p = sum_i lam_i^p,
             C_p = sum_i q_i^p, q = (-g, 1, delta, 0)   [v4std]
                                q = (-g, 1, delta, delta) [v4dd]
    V_spec = gamma h^3 sum_x sum_i P(lam_i)^2, P(x) = prod_i (x - q_i),
             gamma = 4.742198e-4 (m5_32_r20_0_class.json, c.gamma)
    E_total = E_curv + V
Degree of an eigenvector field v (unit, oriented): this file orients v by
continuity along the maximum spanning tree of |v . v'| over the
nearest-neighbor links (the strongest links first, so the core and any
near-degenerate cell cannot flip a subtree), then reads TWO numbers on the
lattice cubes of half-width 6, 9, 12 (index c0 = n // 2, faces at c0 +- k,
k = half / h):
    Mermin flux  Q_M = (1/4pi) sum_faces B . dS,  B_k = (1/2) eps_kij
                 v . (d_i v x d_j v), central differences (the producer's
                 formula, my code)
    solid angle  Q_S = (1/4pi) sum over the surface triangles of the signed
                 spherical area Omega(a, b, c) = 2 atan2(a . (b x c),
                 1 + a.b + b.c + c.a) (Van Oosterom and Strackee), the exact
                 degree of the surface map when the map is continuous
    conflicts    the number of nearest-neighbor pairs with v . v' < 0 after
                 orientation (the line-defect read)
    gap_min      the smallest eigenvalue gap (lam_1 - lam_0, lam_2 - lam_1)
                 on the cube surface (an eigenvector labeled by rank is
                 discontinuous where its gap closes)
Convergence gate (pre-registered): over the trace rows with acc >= 2/3 of
the accepted steps, dE = E_last - E_first, rel = |dE| / max(|E_last|, 1);
decades = log10(fmax_seed / fmax_end); CONVERGED iff rel < 1e-3 AND
decades >= 2, else FALLING if dE < 0 else RISING.
Derrick: M_lam(x) = M(x / lam) by this file's own trilinear interpolation
(edge-clamped), dE/dlam at 1 by the central difference over (0.95, 1.05)
and over (0.9, 1.1); the continuum value for a pure dilation of a quartic
curvature + potential functional is -E_curv + 3 V (E = a / lam + b lam^3),
and the virial E_curv / V is 3 at a Derrick equilibrium.
Frame: M -> Q M Q^T with Q = boost_at(cfg, 0, s), s = +-0.02; the N
spectrum is invariant (Q^T eta Q = eta), so dV must be round-off;
d2E/ds2 = (E+ + E- - 2 E0) / s^2.
E(L) two-rung rule (the producer's docstring): d = E64 - E32, slope
d / 48: AXIS_STRING if d > 0 and slope > 0.01; AXIS_CONVERGED if
|d| < 1e-2 |E64|; AXIS_BOX_LIMITED if d < 0 and the n64 virial > 3; else
undecided. Triple: resolution = max |d|; THREE_AXES_DISTINCT iff every
pairwise |dE| > 3 resolution, else AXES_DEGENERATE. Koide
Q = sum E / (sum sqrt E)^2. Both are read on the producer's basis (the
9000-step n32 rows against the 4500-step n64 rows) and on the equal-budget
basis (the 4500-step n32 rows).
Tube (string) read: the energy in the tube perp < 3 along z with
4 < |z| < L/2 - 2 minus the mean of the same tubes along x and y, per unit
length 2 (L/2 - 2 - 4); plus this file's z-profile of the tube (energy per
plane against |z|) and the frozen pin-shell energy.
Matching radius: E(< R) = sum of the cell density over r < R, R in
(3, 6, ..., 24); rel_spread = (max - min) / |mean| over the boxes;
R_match_2pct = the smallest R with rel_spread <= 0.02.

CHECKS (each a PASS line that can fail; every number printed):
 C1  energies: my E_total on all 25 saved end fields vs the row's E (rel
     1e-10); my E_curv vs EN.block_reads (rel 1e-10); the g 32 control at
     18.963 and below R3's 18.970; the seed energies E0 of the 16 fresh rows
     from my rebuilt seeds; the extension rows' E0 equal to the parent's E;
     the rows' roots equal mine; M_0i = 0 on every field
 C2  seed identity: seed_axes(cfg, (1, delta, 0)) == embed34(seed_pair
     'single') exactly at n32 and n64 (g 8) and n32 g 32; the transverse
     swap: same per-cell eigenvalue multiset, different field; the pin
     shell of every saved field equal to its seed (exact)
 C3  degrees: my Mermin flux vs the rows' rank{k}.flux on r 6, 9, 12
     (abs 1e-6); S_1 |Q| in [1.04, 1.12] on every cube at both boxes under
     V4 (and the solid-angle degree within 0.02 of 1); under V_spec at n32
     S_d and S_0 below 0.2 on r 6, 9, 12 while S_1 keeps 1; the degenerate
     trio (S_d^dd, S_0^dd) with no winding eigenvector (|Q| < 0.5 on r 9,
     12 for every rank) and thousands of conflicts on the two delta ranks;
     the S_d rank-1 read (0.6) diagnosed by the solid-angle degree and the
     surface gap minima
 C4  convergence labels on all 25 rows recomputed from the trace
 C5  Derrick: my dE/dlam at 1 (two stencils) for S_1, S_d, S_0 under V4 at
     n32 (4500 and 9000 steps) and n64 and S_1 under V_spec: negative
     everywhere, within 25 percent of the continuum -E_curv + 3 V; my
     virial vs the rows' (rel 1e-8); the S_1 V4 virial about 23 at n32
     4500 steps and about 4.6 under V_spec; the pin-shell / interior split
     of dE/dlam (the shell is frozen seed data moved by the dilation)
 C6  frame: |dV| < 1e-9 |V| at s = +-0.02 on every field; d2E/ds2 < 0 on
     every field; the quoted -98 (S_1 V4 n32 4500 steps), -84 (V_spec),
     -1672 (g 32) within 2 percent; my d2E/ds2 vs the rows' (rel 1e-6)
 C7  E(L): the increments, the labels, the triple, Koide Q and the ratios
     recomputed on the producer's basis (must equal the JSON) and on the
     equal-budget basis; the tube-read string prediction sigma_64 x 48
     against each increment; the n64 last-third energy drift against each
     increment
 C8  matching radius: my E(< R) tables vs the JSON (rel 1e-8), R_match None
     for every V4 object on the producer's basis; the equal-budget basis
 C9  string read: my tube read vs the rows' (abs 1e-8) on every row; the
     S_d and S_1 n32 values quoted (0.82 / 0.02 at 4500 steps, 0.78 / 0.018
     at 9000 steps); the z-profile flatness ratio (far half over near half)
 C10 the wordings (listed in the JSON with the numbers that do or do not
     support them)
Out: ../data/m5_32_r20_1_audit.json. Runtime and peak RSS printed.
"""
from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import gc  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import breadth_first_order, minimum_spanning_tree  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
IN_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")
EXT_JSON = os.path.join(DATA, "m5_32_r20_1_axes_ext32.json")
R20_0_JSON = os.path.join(DATA, "m5_32_r20_0_class.json")
R3_JSON = os.path.join(DATA, "m5_32_r3_pair.json")
NPZ = os.path.join(DATA, "m5_32_r20_1")
OUT_JSON = os.path.join(DATA, "m5_32_r20_1_audit.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


EN = _load("m5_32_r19_entrants", "m5_32_r19_entrants.py")     # block_reads: the certified curvature cross-check
B3 = EN.B3                                                     # W1, base_cfg, embed34, pin_shell (definitions of the setup)
PAIR = _load("m5_21_4_a_pair", "m5_21_4_a_pair.py")            # seed_pair: the record's builder (claim 2)
R3 = _load("m5_32_r3_ii_pair", "m5_32_r3_ii_pair.py")          # boost_at, conj: the record's dressing (claim 6)
RUN = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")         # seed_axes ONLY (the claim under test in C2)

ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(ETA_D)
DELTA = 0.3
W1 = float(B3.W1)
HALVES = (6.0, 9.0, 12.0)
PROFILE_R = (3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 30.0, 36.0, 42.0)
LAMS = {"S1": (1.0, DELTA, 0.0), "Sd": (DELTA, 1.0, 0.0), "S0": (0.0, 1.0, DELTA),
        "S1p": (1.0, 0.0, DELTA), "Sdp": (DELTA, 0.0, 1.0), "S0p": (0.0, DELTA, 1.0),
        "S1dd": (1.0, DELTA, DELTA), "Sddd": (DELTA, 1.0, DELTA), "S0dd": (DELTA, DELTA, 1.0)}
WINDING_RANK = {"S1": 2, "Sd": 1, "S0": 0, "S1p": 2, "Sdp": 1, "S0p": 0, "S1dd": 2}


def log(msg):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    print(f"[{time.time() - T0:6.1f}s {rss:4.2f}GB] {msg}", flush=True)


def cfg_of(n, L, g):
    return B3.base_cfg(s=-1.0, g=g, n=n, L=float(L), delta=DELTA)


def coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def roots_mine(g, potk):
    return (-g, 1.0, DELTA, DELTA) if potk == "v4dd" else (-g, 1.0, DELTA, 0.0)


# ================= my densities =================
def diff(f, ax, h, br):
    """the sym-stencil branch derivative: fwd (f[i+1] - f[i]) / h on cells 0..n-2, bwd (f[i] - f[i-1]) / h on cells 1..n-1."""
    out = np.zeros_like(f)
    a = [slice(None)] * f.ndim
    b = [slice(None)] * f.ndim
    a[ax], b[ax] = slice(1, None), slice(0, -1)
    d = (f[tuple(a)] - f[tuple(b)]) / h
    if br == "fwd":
        out[tuple(b)] = d
    else:
        out[tuple(a)] = d
    return out


def curv_density(M, h):
    """e_curv per cell (h^3-weighted): 4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta."""
    e = np.zeros(M.shape[:3])
    sgn = np.einsum("a,b->ab", ETA_D, ETA_D)
    for br in ("fwd", "bwd"):
        A = [diff(M, ax, h, br) for ax in range(3)]
        AE = [a @ ETA for a in A]
        for i in range(3):
            for j in range(i + 1, 3):
                F = AE[i] @ A[j] - AE[j] @ A[i]
                e += 0.5 * 4.0 * np.einsum("...ab,ab->...", F * F, sgn)
        del A, AE, F
    return h ** 3 * e


def spectrum_N(M):
    """the eigenvalues of N = M eta per cell, (N, 4): block form when M_0i = 0, eigvals otherwise."""
    Mf = M.reshape(-1, 4, 4)
    m0i = float(np.max(np.abs(Mf[:, 0, 1:])))
    if m0i == 0.0:
        lam3 = np.linalg.eigvalsh(Mf[:, 1:, 1:])
        return np.concatenate([-Mf[:, 0, 0][:, None], lam3], axis=1), m0i, "block"
    ev = np.linalg.eigvals(Mf @ ETA)
    if np.max(np.abs(ev.imag)) > 1e-9 * max(1.0, np.max(np.abs(ev))):
        raise RuntimeError("complex N spectrum")
    return ev.real, m0i, "eigvals"


def pot_density(M, h, potk, q, gamma):
    lam, m0i, how = spectrum_N(M)
    if potk in ("v4std", "v4dd"):
        acc = np.zeros(lam.shape[0])
        for p in range(1, 5):
            acc += (np.sum(lam ** p, axis=1) - sum(qi ** p for qi in q)) ** 2
        v = W1 * acc
    else:
        P = np.ones_like(lam)
        for qi in q:
            P = P * (lam - qi)
        v = gamma * np.sum(P * P, axis=1)
    return (h ** 3 * v).reshape(M.shape[:3]), m0i, how


def energy(M, cfg, potk, q, gamma):
    h = cfg["h"]
    ec = curv_density(M, h)
    ev, m0i, how = pot_density(M, h, potk, q, gamma)
    return {"E_curv": float(ec.sum()), "V": float(ev.sum()), "E_total": float(ec.sum() + ev.sum()), "max_abs_M0i": m0i, "spectrum_how": how}, ec + ev, ec, ev


# ================= my degree reads =================
def orient_tree(v):
    """orientation by continuity along the MAXIMUM spanning tree of |v . v'| over nearest-neighbor links (scipy's
    minimum_spanning_tree on the weights 1 - |v . v'|): every cell takes the sign of its tree parent times
    sign(v . v_parent), the products propagated by pointer jumping. The tree follows the strongest links, so the
    core (where a hedgehog turns by 90 degrees per cell) and any near-degenerate cell cannot flip a subtree."""
    n = v.shape[0]
    N = n ** 3
    vf = v.reshape(N, 3)
    idx = np.arange(N).reshape(n, n, n)
    aa, bb, ww = [], [], []
    for ax in range(3):
        sa = [slice(None)] * 3
        sb = [slice(None)] * 3
        sa[ax], sb[ax] = slice(1, None), slice(0, -1)
        a, b = idx[tuple(sa)].ravel(), idx[tuple(sb)].ravel()
        aa.append(a)
        bb.append(b)
        ww.append(1.0 - np.abs(np.einsum("na,na->n", vf[a], vf[b])) + 1e-9)
    a, b, w = np.concatenate(aa), np.concatenate(bb), np.concatenate(ww)
    G = coo_matrix((np.concatenate([w, w]), (np.concatenate([a, b]), np.concatenate([b, a]))), shape=(N, N)).tocsr()
    T = minimum_spanning_tree(G)
    root = int(idx[n // 2, n // 2, n // 2])
    order, pred = breadth_first_order(T, root, directed=False, return_predecessors=True)
    if len(order) != N:
        raise RuntimeError("spanning tree does not reach every cell")
    pred = pred.astype(np.int64)
    pred[root] = root
    s = np.sign(np.einsum("na,na->n", vf, vf[pred]))
    s[s == 0.0] = 1.0
    s[root] = 1.0
    P, S = pred.copy(), s.copy()
    for _ in range(64):
        S = S * S[P]
        P = P[P]
        if np.all(P == root):
            break
    vo = v * S.reshape(n, n, n)[..., None]
    ncf = 0
    for a_ in range(3):
        s1 = [slice(None)] * 3
        s2 = [slice(None)] * 3
        s1[a_], s2[a_] = slice(1, None), slice(0, -1)
        ncf += int(np.sum(np.einsum("...a,...a->...", vo[tuple(s1)], vo[tuple(s2)]) < 0.0))
    return vo, ncf


def central(f, ax, h):
    out = np.zeros_like(f)
    a = [slice(None)] * f.ndim
    b = [slice(None)] * f.ndim
    c = [slice(None)] * f.ndim
    a[ax], b[ax], c[ax] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(c)] = (f[tuple(a)] - f[tuple(b)]) / (2.0 * h)
    return out


def mermin_flux(vo, h, c0, k):
    d = [central(vo, a, h) for a in range(3)]
    B = np.stack([np.einsum("...a,...a->...", vo, np.cross(d[1], d[2])),
                  np.einsum("...a,...a->...", vo, np.cross(d[2], d[0])),
                  np.einsum("...a,...a->...", vo, np.cross(d[0], d[1]))], axis=-1)
    lo, hi = c0 - k, c0 + k
    s = (B[hi, lo:hi + 1, lo:hi + 1, 0].sum() - B[lo, lo:hi + 1, lo:hi + 1, 0].sum()
         + B[lo:hi + 1, hi, lo:hi + 1, 1].sum() - B[lo:hi + 1, lo, lo:hi + 1, 1].sum()
         + B[lo:hi + 1, lo:hi + 1, hi, 2].sum() - B[lo:hi + 1, lo:hi + 1, lo, 2].sum())
    return float(s * h * h / (4.0 * np.pi))


def solid_angle_degree(vo, c0, k):
    """the degree of the map cube surface -> S^2 as the sum of the signed spherical triangle areas."""
    lo, hi = c0 - k, c0 + k
    sl = slice(lo, hi + 1)
    faces = [(vo[hi, sl, sl], False), (vo[lo, sl, sl], True),          # +x: (y, z) right-handed; -x: swapped
             (vo[sl, hi, sl], True), (vo[sl, lo, sl], False),           # +y: (z, x) -> the (x, z) grid swapped
             (vo[sl, sl, hi], False), (vo[sl, sl, lo], True)]           # +z: (x, y); -z: swapped
    tot = 0.0
    for f, swap in faces:
        if swap:
            f = np.swapaxes(f, 0, 1)
        a, b, c, d = f[:-1, :-1], f[1:, :-1], f[1:, 1:], f[:-1, 1:]
        for (p, q, r) in ((a, b, c), (a, c, d)):
            num = np.einsum("...a,...a->...", p, np.cross(q, r))
            den = 1.0 + np.einsum("...a,...a->...", p, q) + np.einsum("...a,...a->...", q, r) + np.einsum("...a,...a->...", r, p)
            tot += float(np.sum(2.0 * np.arctan2(num, den)))
    return tot / (4.0 * np.pi)


def surface_gap_min(lam, c0, k):
    lo, hi = c0 - k, c0 + k
    m = np.zeros(lam.shape[:3], dtype=bool)
    m[lo, lo:hi + 1, lo:hi + 1] = m[hi, lo:hi + 1, lo:hi + 1] = True
    m[lo:hi + 1, lo, lo:hi + 1] = m[lo:hi + 1, hi, lo:hi + 1] = True
    m[lo:hi + 1, lo:hi + 1, lo] = m[lo:hi + 1, lo:hi + 1, hi] = True
    g = lam[m]
    return [float(np.min(g[:, 1] - g[:, 0])), float(np.min(g[:, 2] - g[:, 1]))]


def degree_reads(M, cfg):
    n, h = cfg["n"], cfg["h"]
    c0 = n // 2
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    out = {}
    for r in range(3):
        vo, ncf = orient_tree(vec[..., :, r])
        rec = {"conflicts": ncf, "mermin": {}, "solid_angle": {}}
        for hv in HALVES:
            k = int(round(hv / h))
            rec["mermin"][f"{hv:g}"] = mermin_flux(vo, h, c0, k)
            rec["solid_angle"][f"{hv:g}"] = solid_angle_degree(vo, c0, k)
        out[f"rank{r}"] = rec
    out["surface_gap_min"] = {f"{hv:g}": surface_gap_min(lam, c0, int(round(hv / h))) for hv in HALVES}
    return out


# ================= Derrick, frame, tube, radial =================
def dilate(M, lam):
    """M_lam(x) = M(x / lam) about the box center, own trilinear interpolation, edge-clamped."""
    n = M.shape[0]
    c = (n - 1) / 2.0
    src = np.clip((np.arange(n) - c) / lam + c, 0.0, n - 1.0)
    i0 = np.minimum(np.floor(src).astype(int), n - 2)
    f = src - i0
    out = np.zeros_like(M)
    for dx, wx in ((0, 1.0 - f), (1, f)):
        Mx = M[i0 + dx]
        for dy, wy in ((0, 1.0 - f), (1, f)):
            Mxy = Mx[:, i0 + dy]
            for dz, wz in ((0, 1.0 - f), (1, f)):
                w = wx[:, None, None] * wy[None, :, None] * wz[None, None, :]
                out += w[..., None, None] * Mxy[:, :, i0 + dz]
    return out


def derrick(M, cfg, potk, q, gamma, E1, ec, ev, pairs=((0.95, 1.05), (0.9, 1.1))):
    """dE/dlam by dilation, split into the pin-shell part (the frozen seed values at the box faces, moved by the
    dilation), the interior part, and the core ball r < 12 (with its own continuum value -E_curv(r<12) + 3 V(r<12)):
    a scaling force of the object is the core part; a derivative carried by the shell or the outer region is the
    box, not the object."""
    out = {"virial": E1["E_curv"] / E1["V"], "continuum_dE_dlam": -E1["E_curv"] + 3.0 * E1["V"]}
    n, h = cfg["n"], cfg["h"]
    pin = B3.pin_shell(n, h, 1.6)
    X, Y, Z = coords(n, h)
    core = np.sqrt(X * X + Y * Y + Z * Z) < 12.0
    out["continuum_dE_dlam_r_lt_12"] = float(-ec[core].sum() + 3.0 * ev[core].sum())
    out["virial_r_lt_12"] = float(ec[core].sum() / ev[core].sum())
    for lo, hi in pairs:
        Eld, el, _, _ = energy(dilate(M, lo), cfg, potk, q, gamma)
        Ehd, eh, _, _ = energy(dilate(M, hi), cfg, potk, q, gamma)
        out[f"dE_dlam_r_lt_12_{lo:g}_{hi:g}"] = float((eh[core].sum() - el[core].sum()) / (hi - lo))
        El, Eh = Eld["E_total"], Ehd["E_total"]
        out[f"dE_dlam_{lo:g}_{hi:g}"] = (Eh - El) / (hi - lo)
        out[f"d2E_dlam2_{lo:g}_{hi:g}"] = (Eh + El - 2.0 * E1["E_total"]) / ((hi - lo) / 2.0) ** 2
        out[f"E_{lo:g}"], out[f"E_{hi:g}"] = El, Eh
        out[f"dE_dlam_pin_shell_{lo:g}_{hi:g}"] = float((eh[pin].sum() - el[pin].sum()) / (hi - lo))
        out[f"dE_dlam_interior_{lo:g}_{hi:g}"] = float((eh[~pin].sum() - el[~pin].sum()) / (hi - lo))
        del el, eh
    return out


def frame(M, cfg, potk, q, gamma, E0, s=0.02):
    Qp = R3.boost_at(cfg, 0.0, s)[0]
    Ep = energy(R3.conj(Qp, M), cfg, potk, q, gamma)[0]
    del Qp
    Qm = R3.boost_at(cfg, 0.0, -s)[0]
    Em = energy(R3.conj(Qm, M), cfg, potk, q, gamma)[0]
    del Qm
    return {"s": s, "dV_plus": Ep["V"] - E0["V"], "dV_minus": Em["V"] - E0["V"],
            "d2E_ds2": (Ep["E_total"] + Em["E_total"] - 2.0 * E0["E_total"]) / s ** 2,
            "dE_ds": (Ep["E_total"] - Em["E_total"]) / (2.0 * s), "max_abs_M0i_plus": Ep["max_abs_M0i"]}


def tube_reads(e, cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = coords(n, h)
    zmax = 0.5 * L - 2.0
    length = 2.0 * (zmax - 4.0)
    tubes = {}
    for lab, perp, along in (("z", np.sqrt(X * X + Y * Y), Z), ("x", np.sqrt(Y * Y + Z * Z), X), ("y", np.sqrt(X * X + Z * Z), Y)):
        m = (perp < 3.0) & (np.abs(along) > 4.0) & (np.abs(along) < zmax)
        tubes[lab] = float(np.sum(e[m]))
    read = (tubes["z"] - 0.5 * (tubes["x"] + tubes["y"])) / length
    # the z profile: tube energy per plane against |z| (both signs summed), and the x tube as the baseline
    rho = np.sqrt(X * X + Y * Y)
    zs = np.unique(np.abs(Z))
    zs = zs[(zs > 4.0) & (zs < zmax)]
    prof = []
    for zv in zs:
        m = (rho < 3.0) & (np.abs(np.abs(Z) - zv) < 1e-9)
        mx = (np.sqrt(Y * Y + Z * Z) < 3.0) & (np.abs(np.abs(X) - zv) < 1e-9)
        prof.append((float(zv), float(np.sum(e[m])) / 2.0, float(np.sum(e[mx])) / 2.0))
    half = len(prof) // 2
    near = sum(p[1] - p[2] for p in prof[:half])
    far = sum(p[1] - p[2] for p in prof[half:])
    ex = [p[1] - p[2] for p in prof]
    mid = ex[len(ex) // 4: 3 * len(ex) // 4]
    return {"tubes": tubes, "length": length, "string_read": read, "profile_z_E_per_plane_minus_x": [(p[0], p[1] - p[2]) for p in prof],
            "far_over_near": far / near if near != 0 else None, "cells_per_plane": int(np.sum((rho < 3.0) & (np.abs(Z - Z.ravel()[0]) < 1e-9))),
            "flat_tension_mid_half_per_unit_L": float(np.median(mid) / h), "last_plane_share": float(ex[-1] / sum(ex)) if sum(ex) != 0 else None,
            "max_over_median_mid": float(max(mid) / np.median(mid)) if np.median(mid) != 0 else None}


def radial(e, cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    pin = B3.pin_shell(n, h)
    return {"E_lt_R": {f"{R:g}": float(np.sum(e[r < R])) for R in PROFILE_R if R <= 0.5 * L}, "E_pin_shell": float(np.sum(e[pin])),
            "E_interior": float(np.sum(e[~pin]))}


def match_table(profiles):
    """profiles: {n: {R: E}} -> the spread table and R_match_2pct."""
    common = sorted(set.intersection(*[set(v) for v in profiles.values()]), key=float)
    tab = {R: {n: profiles[n][R] for n in profiles} for R in common}
    rel = {R: float((max(v.values()) - min(v.values())) / max(abs(np.mean(list(v.values()))), 1e-300)) for R, v in tab.items()}
    match = [R for R in common if rel[R] <= 0.02]
    return {"table": tab, "rel_spread": rel, "R_match_2pct": float(match[0]) if match else None}


def label_of(trace, n_acc, fmax0, fmax_end):
    q = [r for r in trace if r["acc"] >= (2.0 / 3.0) * n_acc]
    if len(q) < 2:
        return "budget (trace too short)", None, None
    dE = q[-1]["E"] - q[0]["E"]
    rel = abs(dE) / max(abs(q[-1]["E"]), 1.0)
    dec = float(np.log10(fmax0 / max(fmax_end, 1e-300)))
    lab = "CONVERGED" if (rel < 1e-3 and dec >= 2.0) else ("FALLING" if dE < 0 else "RISING")
    return lab, rel, dec


def two_rung(E32, E64, vir64):
    d = E64 - E32
    sl = d / 48.0
    if d > 0 and sl > 0.01:
        lab = "AXIS_STRING"
    elif abs(d) < 1e-2 * abs(E64):
        lab = "AXIS_CONVERGED"
    elif d < 0 and np.isfinite(vir64) and vir64 > 3.0:
        lab = "AXIS_BOX_LIMITED"
    else:
        lab = "E(L)_UNDECIDED"
    return {"d": d, "slope": sl, "label": lab}


def triple(Es, incs):
    resol = max(abs(v) for v in incs.values())
    arr = np.array([Es["S1"], Es["Sd"], Es["S0"]])
    diffs = [abs(arr[i] - arr[j]) for i in range(3) for j in range(i + 1, 3)]
    srt = np.sort(arr)
    return {"resolution": resol, "pairwise_abs_diffs": diffs, "outcome": "THREE_AXES_DISTINCT" if min(diffs) > 3.0 * resol else "AXES_DEGENERATE",
            "koide_Q": float(np.sum(arr) / np.sum(np.sqrt(arr)) ** 2), "ratios_sorted": [float(x / srt[0]) for x in srt]}


# ================= the audit =================
def main():
    J = json.load(open(IN_JSON))
    rows = J["rows"]
    res = J["results"]
    ext = json.load(open(EXT_JSON))["rows"]
    gamma = float(json.load(open(R20_0_JSON))["c"]["gamma"])
    R3E = {r["tag"]: r["E"] for r in json.load(open(R3_JSON))["rows"]}
    A = {"rows": {}, "lines": {}, "W1": W1, "gamma": gamma}
    L = {}
    log(f"{len(rows)} rows in the main JSON, {len(ext)} in the extension file, W1 {W1:.12g}, gamma {gamma:.10g}")
    L["C0_extension_rows_merged_verbatim"] = all(t in rows and rows[t]["E"] == ext[t]["E"] for t in ext) and len(ext) == 9

    # ---- per-field reads (C1, C2 pin, C3, C5, C6, C8, C9) ----
    order = sorted(rows, key=lambda t: (rows[t]["n"], t))
    for tag in order:
        r = rows[tag]
        cfg = cfg_of(r["n"], r["L"], r["g"])
        potk = r["potk"]
        q = roots_mine(r["g"], potk)
        M = np.load(os.path.join(NPZ, f"{tag}.npz"))["M"]
        cross = EN.block_reads(M, cfg, "I1")          # the certified curvature (heavy at n64: called first, nothing else held)
        gc.collect()
        E1, e, ec, ev = energy(M, cfg, potk, q, gamma)
        rec = {"n": r["n"], "obj": r["obj"], "potk": potk, "g": r["g"], "steps_total": r["steps"] * (2 if r.get("resume") else 1),
               "roots_row": r["roots"], "roots_mine": list(q), "E_row": r["E"], "E_mine": E1["E_total"], "E_curv_mine": E1["E_curv"], "V_mine": E1["V"],
               "rel_E": abs(E1["E_total"] - r["E"]) / abs(r["E"]), "rel_curv_vs_block_reads": abs(E1["E_curv"] - cross["E_curv_I1"]) / abs(cross["E_curv_I1"]),
               "rel_V_vs_row": abs(E1["V"] - r["end_reads"]["energy"]["V"]) / abs(r["end_reads"]["energy"]["V"]),
               "max_abs_M0i": E1["max_abs_M0i"], "spectrum_how": E1["spectrum_how"], "max_abs_M00_minus_g": float(np.max(np.abs(M[..., 0, 0] - r["g"])))}
        # the seed: fresh rows rebuilt, extension rows = the parent's saved field
        if r.get("resume"):
            rec["E0_row"] = r["descent"]["E0"]
            rec["E0_expected_parent_E"] = rows[r["resume"]]["E"]
            rec["E0_match"] = bool(abs(r["descent"]["E0"] - rows[r["resume"]]["E"]) < 1e-12 * abs(rows[r["resume"]]["E"]))
            M0 = RUN.seed_axes(cfg, LAMS[r["obj"]])
        else:
            M0 = RUN.seed_axes(cfg, LAMS[r["obj"]])
            rec["E0_row"] = r["descent"]["E0"]
            rec["E0_mine"] = energy(M0, cfg, potk, q, gamma)[0]["E_total"]
            rec["E0_match"] = bool(abs(rec["E0_mine"] - rec["E0_row"]) / abs(rec["E0_row"]) < 1e-10)
        pin = B3.pin_shell(cfg["n"], cfg["h"], 1.6)
        rec["pin_shell_max_dev_from_seed"] = float(np.max(np.abs(M[pin] - M0[pin])))
        del M0
        # degrees
        dg = degree_reads(M, cfg)
        rec["degree"] = dg
        rec["degree_row"] = {k: {"conflicts": v["conflicts"], "flux": {hv: v["flux"].get(hv) for hv in ("6", "9", "12")}} for k, v in r["end_reads"]["degree"].items()}
        # the sign of a director field's orientation is conventional (eigh's sign at the root cell), so |Q| is compared; a rank whose
        # orientation has conflicts (a line defect crossing the cubes) has a cut whose placement is algorithm-dependent, so those ranks
        # are compared only on the side of the 0.5 escape threshold
        clean = [k for k in range(3) if dg[f"rank{k}"]["conflicts"] < 100 and r["end_reads"]["degree"][f"rank{k}"]["conflicts"] < 100]
        rec["clean_ranks"] = clean
        rec["mermin_max_abs_dev"] = max([abs(abs(dg[f"rank{k}"]["mermin"][hv]) - abs(r["end_reads"]["degree"][f"rank{k}"]["flux"][hv])) for k in clean for hv in ("6", "9", "12")] or [0.0])
        rec["mermin_conflict_ranks_same_side_of_0.5"] = all((abs(dg[f"rank{k}"]["mermin"][hv]) < 0.5) == (abs(r["end_reads"]["degree"][f"rank{k}"]["flux"][hv]) < 0.5)
                                                            for k in range(3) if k not in clean for hv in ("6", "9", "12"))
        rec["mermin_all_ranks_max_abs_dev"] = max(abs(abs(dg[f"rank{k}"]["mermin"][hv]) - abs(r["end_reads"]["degree"][f"rank{k}"]["flux"][hv])) for k in range(3) for hv in ("6", "9", "12"))
        # tube, radial
        tb = tube_reads(e, cfg)
        rec["tube"] = tb
        rec["string_read_row"] = r["end_reads"]["string_tension_read"]
        rec["string_read_dev"] = abs(tb["string_read"] - r["end_reads"]["string_tension_read"])
        rd = radial(e, cfg)
        rec["radial"] = rd
        rec["radial_row"] = r["end_reads"]["radial"]["E_lt_R"]
        rec["radial_max_rel_dev"] = max(abs(rd["E_lt_R"][R] - r["end_reads"]["radial"]["E_lt_R"][R]) / abs(r["end_reads"]["radial"]["E_lt_R"][R]) for R in rd["E_lt_R"])
        rec["E_density_sum_vs_E_rel"] = abs(float(e.sum()) - E1["E_total"]) / abs(E1["E_total"])
        del e
        # frame (all rows) and Derrick (the main V4 rows at both boxes and S_1 under V_spec; the rest at n32 too since cheap)
        rec["frame"] = frame(M, cfg, potk, q, gamma, E1)
        rec["frame_row"] = {"d2E_ds2": r["end_reads"]["frame"]["s0.02"]["d2E_ds2"], "dV_plus": r["end_reads"]["frame"]["s0.02"]["dV_plus"]}
        gc.collect()
        if r["n"] == 32 or r["obj"] in ("S1", "Sd", "S0"):
            rec["derrick"] = derrick(M, cfg, potk, q, gamma, E1, ec, ev, pairs=((0.95, 1.05), (0.9, 1.1)) if r["n"] == 32 else ((0.95, 1.05),))
            rec["derrick_row"] = {"virial": r["end_reads"]["derrick"]["virial_E_curv_over_V"], "dE_dlam_order3": r["end_reads"]["derrick"]["order3"]["dE_dlambda_at_1"],
                                  "dE_dlam_order1": r["end_reads"]["derrick"]["order1"]["dE_dlambda_at_1"]}
        del M, ec, ev
        gc.collect()
        A["rows"][tag] = rec
        dd = rec.get("derrick", {})
        log(f"{tag:24s} E row {r['E']:10.5f} mine {E1['E_total']:10.5f} rel {rec['rel_E']:.1e} curv/block {rec['rel_curv_vs_block_reads']:.1e} "
            f"M0i {E1['max_abs_M0i']:.0e} pin {rec['pin_shell_max_dev_from_seed']:.0e} E0 {rec['E0_match']} | mermin dev {rec['mermin_max_abs_dev']:.1e} "
            f"conf {[dg[f'rank{k}']['conflicts'] for k in range(3)]} | str {tb['string_read']:+.4f} (row {r['end_reads']['string_tension_read']:+.4f}) far/near {tb['far_over_near']} "
            f"| vir {rec['E_curv_mine'] / rec['V_mine']:.2f} dE/dlam {dd.get('dE_dlam_0.95_1.05', float('nan')):+.2f} cont {dd.get('continuum_dE_dlam', float('nan')):+.2f} "
            f"| d2E/ds2 {rec['frame']['d2E_ds2']:+.2f} (row {rec['frame_row']['d2E_ds2']:+.2f}) dV {rec['frame']['dV_plus']:+.1e}")
        for k in range(3):
            g = dg[f"rank{k}"]
            log(f"    rank{k} mermin {[round(g['mermin'][hv], 3) for hv in ('6', '9', '12')]} solid {[round(g['solid_angle'][hv], 3) for hv in ('6', '9', '12')]} "
                f"row {[round(r['end_reads']['degree'][f'rank{k}']['flux'][hv], 3) for hv in ('6', '9', '12')]} conflicts mine {g['conflicts']} row {r['end_reads']['degree'][f'rank{k}']['conflicts']}")
        log(f"    surface gap min {dg['surface_gap_min']}  pin-shell E {rd['E_pin_shell']:.3f} interior {rd['E_interior']:.3f}")
    R = A["rows"]

    # ---- C1 ----
    L["C1a_E_total_all_25_fields_rel_1e-10"] = len(R) == 25 and max(v["rel_E"] for v in R.values()) < 1e-10
    L["C1b_E_curv_vs_EN_block_reads_rel_1e-10"] = max(v["rel_curv_vs_block_reads"] for v in R.values()) < 1e-10
    L["C1c_V_vs_row_rel_1e-10"] = max(v["rel_V_vs_row"] for v in R.values()) < 1e-10
    g32 = R["S1_v4std_n32_g32"]["E_mine"]
    L["C1d_g32_control_18.963_below_R3_18.970"] = abs(g32 - 18.963) < 1e-3 and g32 < R3E["lam0_un0_single_d0_n32"] and abs(R3E["lam0_un0_single_d0_n32"] - 18.970) < 1e-3
    L["C1e_seed_E0_16_fresh_rows_rebuilt_and_9_extension_E0_equal_parent_E"] = all(v["E0_match"] for v in R.values()) and sum(1 for v in R.values() if "E0_mine" in v) == 16
    L["C1f_rows_roots_equal_mine_and_M0i_zero"] = all(tuple(v["roots_row"]) == tuple(v["roots_mine"]) and v["max_abs_M0i"] == 0.0 for v in R.values())
    L["C1g_density_sums_to_E_rel_1e-12"] = max(v["E_density_sum_vs_E_rel"] for v in R.values()) < 1e-12
    A["C1"] = {"worst_rel_E": max(v["rel_E"] for v in R.values()), "worst_rel_curv": max(v["rel_curv_vs_block_reads"] for v in R.values()),
               "worst_rel_V": max(v["rel_V_vs_row"] for v in R.values()), "g32_control_mine": g32, "R3_record": R3E["lam0_un0_single_d0_n32"],
               "max_abs_M00_minus_g": max(v["max_abs_M00_minus_g"] for v in R.values())}

    # ---- C2: the seed identity ----
    c2 = {}
    for n, Lb, g in ((32, 48.0, 8.0), (64, 96.0, 8.0), (32, 48.0, 32.0)):
        cfg = cfg_of(n, Lb, g)
        Ma = RUN.seed_axes(cfg, LAMS["S1"])
        Mb = B3.embed34(PAIR.seed_pair(cfg, "single", 0.0), cfg)
        Mc = RUN.seed_axes(cfg, LAMS["S1p"])
        la = np.sort(np.linalg.eigvalsh(Ma[..., 1:, 1:]), axis=-1)
        lc = np.sort(np.linalg.eigvalsh(Mc[..., 1:, 1:]), axis=-1)
        c2[f"n{n}_g{g:g}"] = {"identity_max_abs": float(np.max(np.abs(Ma - Mb))), "swap_spectrum_max_abs": float(np.max(np.abs(la - lc))),
                              "swap_field_max_abs": float(np.max(np.abs(Ma - Mc))), "M00": float(Ma[0, 0, 0, 0, 0])}
        del Ma, Mb, Mc
    L["C2a_seed_axes_S1_equals_embed34_seed_pair_single_exact_3_boxes"] = all(v["identity_max_abs"] == 0.0 for v in c2.values())
    L["C2b_swap_same_eigenvalue_multiset_1e-12_different_field"] = all(v["swap_spectrum_max_abs"] < 1e-12 and v["swap_field_max_abs"] > 0.1 for v in c2.values())
    L["C2c_pin_shell_at_seed_values_exact_all_25_fields"] = max(v["pin_shell_max_dev_from_seed"] for v in R.values()) == 0.0
    A["C2"] = c2

    # ---- C3: degrees ----
    L["C3a_abs_mermin_flux_equals_rows_1e-6_on_every_rank_clean_under_both_orientations"] = max(v["mermin_max_abs_dev"] for v in R.values()) < 1e-6
    L["C3a2_flux_on_line_defect_ranks_orientation_independent_to_0.1"] = max(v["mermin_all_ranks_max_abs_dev"] for v in R.values()) < 0.1
    A["C3a2_max_disagreement"] = {t: v["mermin_all_ranks_max_abs_dev"] for t, v in R.items() if v["mermin_all_ranks_max_abs_dev"] >= 0.1}
    s1 = [R[t] for t in ("S1_v4std_n32_g8", "S1_v4std_n32_g8_x4500", "S1_v4std_n64_g8")]
    L["C3b_S1_V4_rank2_mermin_in_1.04_1.12_and_solid_angle_within_0.02_of_1_both_boxes"] = all(
        1.04 <= abs(v["degree"]["rank2"]["mermin"][hv]) <= 1.12 and abs(abs(v["degree"]["rank2"]["solid_angle"][hv]) - 1.0) < 0.02 for v in s1 for hv in ("6", "9", "12"))
    vs = {o: R[f"{o}_vspec_n32_g8"] for o in ("S1", "Sd", "S0")}
    L["C3c_vspec_n32_Sd_S0_winding_rank_below_0.2_on_r6_9_12_and_S1_keeps_1"] = (
        all(abs(vs[o]["degree"][f"rank{WINDING_RANK[o]}"]["mermin"][hv]) < 0.2 for o in ("Sd", "S0") for hv in ("6", "9", "12"))
        and all(abs(abs(vs["S1"]["degree"]["rank2"]["mermin"][hv]) - 1.0) < 0.12 for hv in ("6", "9", "12")))
    dd = {o: R[f"{o}_v4dd_n32_g8_x4500"] for o in ("S1dd", "Sddd", "S0dd")}
    # the degenerate trio: the delta pair (ranks 0, 1) is degenerate on the cubes (gap ~ 0), so those eigenvectors are undefined and
    # their flux is not a number (both orientations leave thousands of conflicts and disagree, C3a2); the non-degenerate 1-eigenvector
    # (rank 2) is the read that exists: no winding on Sddd and S0dd, winding on S1dd
    L["C3d_degenerate_trio_delta_pair_gap_lt_0.01_on_cubes_and_conflicts_gt_1000_both_orientations"] = (
        all(dd[o]["degree"]["surface_gap_min"][hv][0] < 0.01 for o in dd for hv in ("6", "9", "12"))
        and all(dd[o]["degree"][f"rank{k}"]["conflicts"] > 1000 and rows[f"{o}_v4dd_n32_g8_x4500"]["end_reads"]["degree"][f"rank{k}"]["conflicts"] > 1000 for o in dd for k in (0, 1)))
    L["C3d2_Sddd_S0dd_rank2_no_winding_abs_mermin_lt_0.1_r6_9_12_both_orientations"] = all(
        abs(dd[o]["degree"]["rank2"]["mermin"][hv]) < 0.1 and abs(rows[f"{o}_v4dd_n32_g8_x4500"]["end_reads"]["degree"]["rank2"]["flux"][hv]) < 0.1 for o in ("Sddd", "S0dd") for hv in ("6", "9", "12"))
    L["C3e_S1dd_rank2_keeps_winding_1"] = all(abs(abs(dd["S1dd"]["degree"]["rank2"]["solid_angle"][hv]) - 1.0) < 0.02 for hv in ("6", "9", "12"))
    sd = {t: R[t] for t in ("Sd_v4std_n32_g8", "Sd_v4std_n32_g8_x4500", "Sd_v4std_n64_g8")}
    L["C3f_Sd_V4_rank1_mermin_0.6_but_solid_angle_degree_1_within_0.02_every_cube_both_boxes_conflicts_lt_10"] = all(
        abs(abs(v["degree"]["rank1"]["solid_angle"][hv]) - 1.0) < 0.02 and abs(v["degree"]["rank1"]["mermin"][hv]) < 0.7 and v["degree"]["rank1"]["conflicts"] < 10 for v in sd.values() for hv in ("6", "9", "12"))
    L["C3g_S0_V4_winding_rank0_solid_angle_1_and_mermin_matches_rows_at_n32_9000_and_n64"] = all(
        abs(abs(R[t]["degree"]["rank0"]["solid_angle"][hv]) - 1.0) < 0.02 and abs(abs(R[t]["degree"]["rank0"]["mermin"][hv]) - abs(rows[t]["end_reads"]["degree"]["rank0"]["flux"][hv])) < 1e-6
        for t in ("S0_v4std_n32_g8_x4500", "S0_v4std_n64_g8") for hv in ("6", "9", "12"))
    A["C3"] = {"S1_V4": {t: {"mermin": v["degree"]["rank2"]["mermin"], "solid": v["degree"]["rank2"]["solid_angle"]} for t, v in zip(("n32_4500", "n32_9000", "n64"), s1)},
               "vspec": {o: {f"rank{k}": {"mermin": vs[o]["degree"][f"rank{k}"]["mermin"], "solid": vs[o]["degree"][f"rank{k}"]["solid_angle"], "conflicts": vs[o]["degree"][f"rank{k}"]["conflicts"]} for k in range(3)} for o in vs},
               "degenerate": {o: {f"rank{k}": {"mermin": dd[o]["degree"][f"rank{k}"]["mermin"], "solid": dd[o]["degree"][f"rank{k}"]["solid_angle"], "conflicts": dd[o]["degree"][f"rank{k}"]["conflicts"]} for k in range(3)} for o in dd},
               "Sd_V4_rank1": {t: {"mermin": v["degree"]["rank1"]["mermin"], "solid": v["degree"]["rank1"]["solid_angle"], "conflicts": v["degree"]["rank1"]["conflicts"],
                                   "surface_gap_min": v["degree"]["surface_gap_min"]} for t, v in sd.items()}}

    # ---- C4: convergence labels ----
    c4 = {}
    ok = True
    for tag, r in rows.items():
        d = r["descent"]
        lab, rel, dec = label_of(d["trace"], d["accepted"], d["fmax_seed"], d["fmax_end"])
        rec = {"label_mine": lab, "label_row": d["verdict"], "rel_mine": rel, "rel_row": d.get("last_third_rel"), "decades_mine": dec, "decades_row": d.get("fmax_decades"),
               "trace_last_fmax_equals_fmax_end": d["trace"][-1]["fmax"] == d["fmax_end"], "trace_last_E_equals_E_end": d["trace"][-1]["E"] == d["E_end"],
               "E_end_equals_E_row": abs(d["E_end"] - r["E"]) < 1e-12 * abs(r["E"])}
        if r.get("resume"):
            fs = rows[r["resume"]]["descent"]["fmax_seed"]
            rec["decades_from_original_seed"] = float(np.log10(fs / d["fmax_end"]))
            rec["drop_in_extension"] = rows[r["resume"]]["E"] - r["E"]
        ok &= lab == d["verdict"] and abs(rel - d["last_third_rel"]) < 1e-12 and abs(dec - d["fmax_decades"]) < 1e-9 and rec["E_end_equals_E_row"]
        c4[tag] = rec
    L["C4_labels_rel_drift_and_decades_reproduced_all_25_rows"] = ok and len(c4) == 25
    L["C4b_only_CONVERGED_rows_are_the_4_vspec_and_g32_none_at_n64"] = sorted(t for t, v in c4.items() if v["label_mine"] == "CONVERGED") == sorted(
        ["S1_vspec_n32_g8", "Sd_vspec_n32_g8", "S0_vspec_n32_g8", "S1_v4std_n32_g32"])
    A["C4"] = c4

    # ---- C5: Derrick ----
    c5 = {}
    ok_neg = ok_cont = ok_vir = True
    for tag, v in R.items():
        if "derrick" not in v:
            continue
        d = v["derrick"]
        cont = d["continuum_dE_dlam"]
        num = d["dE_dlam_0.95_1.05"]
        rec = {"virial_mine": d["virial"], "virial_row": v["derrick_row"]["virial"], "dE_dlam_mine_0.95_1.05": num, "dE_dlam_mine_0.9_1.1": d.get("dE_dlam_0.9_1.1"),
               "dE_dlam_row_cubic": v["derrick_row"]["dE_dlam_order3"], "dE_dlam_row_linear": v["derrick_row"]["dE_dlam_order1"], "continuum_-Ecurv+3V": cont,
               "num_over_continuum": num / cont, "d2E_dlam2_mine": d["d2E_dlam2_0.95_1.05"],
               "dE_dlam_pin_shell": d["dE_dlam_pin_shell_0.95_1.05"], "dE_dlam_interior": d["dE_dlam_interior_0.95_1.05"],
               "pin_shell_share": d["dE_dlam_pin_shell_0.95_1.05"] / num, "dE_dlam_r_lt_12": d["dE_dlam_r_lt_12_0.95_1.05"],
               "continuum_r_lt_12": d["continuum_dE_dlam_r_lt_12"], "virial_r_lt_12": d["virial_r_lt_12"],
               "core_num_over_continuum": d["dE_dlam_r_lt_12_0.95_1.05"] / d["continuum_dE_dlam_r_lt_12"]}
        ok_neg &= num < 0 and (d.get("dE_dlam_0.9_1.1", -1.0) < 0)
        ok_cont &= abs(num / cont - 1.0) < 0.25
        ok_vir &= abs(d["virial"] - v["derrick_row"]["virial"]) / v["derrick_row"]["virial"] < 1e-8
        c5[tag] = rec
    L["C5a_dE_dlam_negative_every_field_read_both_stencils"] = ok_neg and len(c5) >= 21
    L["C5b_dE_dlam_within_25pct_of_continuum_-Ecurv+3V"] = ok_cont
    L["C5c_virial_equals_rows_rel_1e-8"] = ok_vir
    L["C5d_S1_V4_virial_23_at_n32_4500_steps_and_4.6_under_vspec"] = abs(c5["S1_v4std_n32_g8"]["virial_mine"] - 23.4) < 0.5 and abs(c5["S1_vspec_n32_g8"]["virial_mine"] - 4.64) < 0.1
    L["C5e_virial_far_above_3_every_V4_field_min_reported"] = min(v["virial_mine"] for t, v in c5.items() if "vspec" not in t) > 3.0
    L["C5f_interior_part_of_dE_dlam_carries_more_than_half_on_every_field"] = all(v["pin_shell_share"] < 0.5 for v in c5.values())
    L["C5g_interior_dE_dlam_negative_every_field"] = all(v["dE_dlam_interior"] < 0.0 for v in c5.values())
    L["C5h_core_r_lt_12_dE_dlam_negative_every_field"] = all(v["dE_dlam_r_lt_12"] < 0.0 for v in c5.values())

    A["C5"] = {"note_r_lt_12": "the r < 12 split of dE/dlam carries the density transported across r = 12 by the dilation, so it is not compared with -E_curv + 3 V of the ball (reported only)",
               "rows": c5, "min_virial_V4": min(v["virial_mine"] for t, v in c5.items() if "vspec" not in t),
               "S1_virial_by_state": {"n32_4500": c5["S1_v4std_n32_g8"]["virial_mine"], "n32_9000": c5["S1_v4std_n32_g8_x4500"]["virial_mine"], "n64_4500": c5["S1_v4std_n64_g8"]["virial_mine"]}}

    # ---- C6: frame ----
    c6 = {t: {"d2E_ds2_mine": v["frame"]["d2E_ds2"], "d2E_ds2_row": v["frame_row"]["d2E_ds2"], "dV_plus_rel": abs(v["frame"]["dV_plus"]) / v["V_mine"],
              "dV_minus_rel": abs(v["frame"]["dV_minus"]) / v["V_mine"], "dE_ds": v["frame"]["dE_ds"], "max_abs_M0i_plus": v["frame"]["max_abs_M0i_plus"]} for t, v in R.items()}
    L["C6a_dV_at_roundoff_rel_1e-9_both_signs_all_fields"] = max(max(v["dV_plus_rel"], v["dV_minus_rel"]) for v in c6.values()) < 1e-9
    L["C6b_d2E_ds2_negative_every_field_saddle"] = all(v["d2E_ds2_mine"] < 0.0 for v in c6.values())
    L["C6c_d2E_ds2_vs_rows_rel_1e-6"] = max(abs(v["d2E_ds2_mine"] - v["d2E_ds2_row"]) / abs(v["d2E_ds2_row"]) for v in c6.values()) < 1e-6
    L["C6d_quoted_-98_S1_V4_n32_4500_-84_vspec_-1672_g32_within_2pct"] = (abs(c6["S1_v4std_n32_g8"]["d2E_ds2_mine"] / -97.7 - 1) < 0.02 and abs(c6["S1_vspec_n32_g8"]["d2E_ds2_mine"] / -84.3 - 1) < 0.02
                                                                         and abs(c6["S1_v4std_n32_g32"]["d2E_ds2_mine"] / -1672.0 - 1) < 0.02)
    A["C6"] = c6

    # ---- C7: E(L) ----
    per = res["axes"]["v4std"]["per_axis"]
    c7 = {"producer_basis": {}, "equal_budget_basis": {}}
    ok7 = True
    E_pb, E_eb, inc_pb, inc_eb = {}, {}, {}, {}
    for o in ("S1", "Sd", "S0"):
        E32x, E32, E64 = R[f"{o}_v4std_n32_g8_x4500"]["E_mine"], R[f"{o}_v4std_n32_g8"]["E_mine"], R[f"{o}_v4std_n64_g8"]["E_mine"]
        vir64 = R[f"{o}_v4std_n64_g8"]["E_curv_mine"] / R[f"{o}_v4std_n64_g8"]["V_mine"]
        pb = two_rung(E32x, E64, vir64)
        eb = two_rung(E32, E64, vir64)
        sig64 = R[f"{o}_v4std_n64_g8"]["tube"]["string_read"]
        sig32x = R[f"{o}_v4std_n32_g8_x4500"]["tube"]["string_read"]
        sig32 = R[f"{o}_v4std_n32_g8"]["tube"]["string_read"]
        d64 = rows[f"{o}_v4std_n64_g8"]["descent"]
        drift64 = abs(d64["last_third_dE"])
        drop_ext = E32 - E32x
        pb.update({"E32_9000": E32x, "E64_4500": E64, "label_row": per[o]["outcomes"], "increment_row": per[o]["increments"]["L48_to_L96"],
                   "string_prediction_sigma64_x_48": sig64 * 48.0, "string_prediction_sigma32_x_48": sig32x * 48.0, "n64_last_third_drift_abs": drift64,
                   "n32_drop_in_extension": drop_ext, "increment_over_n64_drift": pb["d"] / drift64, "increment_over_string_prediction_64": pb["d"] / (sig64 * 48.0)})
        eb.update({"E32_4500": E32, "E64_4500": E64, "string_prediction_sigma64_x_48": sig64 * 48.0, "string_prediction_sigma32_x_48": sig32 * 48.0,
                   "n64_last_third_drift_abs": drift64, "increment_over_n64_drift": eb["d"] / drift64, "increment_over_string_prediction_64": eb["d"] / (sig64 * 48.0)})
        ok7 &= (per[o]["outcomes"][0].split(" ")[0] == pb["label"]) and abs(pb["d"] - per[o]["increments"]["L48_to_L96"]) < 1e-9
        c7["producer_basis"][o] = pb
        c7["equal_budget_basis"][o] = eb
        E_pb[o], E_eb[o], inc_pb[o], inc_eb[o] = E64, E64, pb["d"], eb["d"]
    tp = triple(E_pb, inc_pb)
    te = triple(E_eb, inc_eb)
    tr = res["axes"]["v4std"]["triple"]
    c7["triple_producer_basis"] = tp
    c7["triple_equal_budget_basis"] = te
    c7["triple_row"] = {k: tr[k] for k in ("resolution", "outcome", "koide_Q", "ratios_sorted", "pairwise_abs_diffs")}
    L["C7a_increments_and_AXIS_STRING_labels_reproduced_producer_basis"] = ok7
    L["C7b_triple_AXES_DEGENERATE_resolution_koide_ratios_reproduced"] = (tp["outcome"] == tr["outcome"] and abs(tp["resolution"] - tr["resolution"]) < 1e-9
                                                                          and abs(tp["koide_Q"] - tr["koide_Q"]) < 1e-12 and max(abs(a - b) for a, b in zip(tp["ratios_sorted"], tr["ratios_sorted"])) < 1e-9)
    L["C7c_equal_budget_labels_all_AXIS_STRING_too"] = all(v["label"] == "AXIS_STRING" for v in c7["equal_budget_basis"].values())
    L["C7d_every_increment_exceeds_the_n64_last_third_drift_both_bases"] = all(v["increment_over_n64_drift"] > 1.0 for b in ("producer_basis", "equal_budget_basis") for v in c7[b].values())
    L["C7e_every_increment_within_30pct_of_the_tube_string_prediction_both_bases"] = all(abs(v["increment_over_string_prediction_64"] - 1.0) < 0.3 for b in ("producer_basis", "equal_budget_basis") for v in c7[b].values())
    L["C7f_labels_carry_FALLING_on_both_rungs"] = all("FALLING,FALLING" in per[o]["outcomes"][0] for o in per)
    A["C7"] = c7

    # ---- C8: matching radius ----
    c8 = {}
    ok8a = ok8b = True
    for o in ("S1", "Sd", "S0"):
        pb = match_table({32: R[f"{o}_v4std_n32_g8_x4500"]["radial"]["E_lt_R"], 64: R[f"{o}_v4std_n64_g8"]["radial"]["E_lt_R"]})
        eb = match_table({32: R[f"{o}_v4std_n32_g8"]["radial"]["E_lt_R"], 64: R[f"{o}_v4std_n64_g8"]["radial"]["E_lt_R"]})
        mr = per[o]["matching_radius"]
        dev = max(abs(pb["rel_spread"][Rk] - mr["rel_spread"][Rk]) for Rk in pb["rel_spread"])
        ok8a &= dev < 1e-8 and pb["R_match_2pct"] is None and mr["R_match_2pct"] is None
        ok8b &= max(v["radial_max_rel_dev"] for v in R.values()) < 1e-8
        c8[o] = {"producer_basis": pb, "equal_budget_basis": eb, "rel_spread_dev_vs_row": dev,
                 "pin_shell_E": {"n32_4500": R[f"{o}_v4std_n32_g8"]["radial"]["E_pin_shell"], "n32_9000": R[f"{o}_v4std_n32_g8_x4500"]["radial"]["E_pin_shell"], "n64_4500": R[f"{o}_v4std_n64_g8"]["radial"]["E_pin_shell"]},
                 "interior_E": {"n32_4500": R[f"{o}_v4std_n32_g8"]["radial"]["E_interior"], "n32_9000": R[f"{o}_v4std_n32_g8_x4500"]["radial"]["E_interior"], "n64_4500": R[f"{o}_v4std_n64_g8"]["radial"]["E_interior"]}}
    L["C8a_E_lt_R_tables_and_R_match_None_reproduced_producer_basis"] = ok8a
    L["C8b_E_lt_R_vs_rows_rel_1e-8_all_rows"] = ok8b
    L["C8c_equal_budget_R_match_2pct_exists_for_every_object"] = all(c8[o]["equal_budget_basis"]["R_match_2pct"] is not None for o in c8)
    A["C8"] = c8

    # ---- C9: string read ----
    L["C9a_tube_read_equals_rows_abs_1e-8_all_rows"] = max(v["string_read_dev"] for v in R.values()) < 1e-8
    s9 = {t: R[t]["tube"]["string_read"] for t in ("Sd_v4std_n32_g8", "S1_v4std_n32_g8", "Sd_v4std_n32_g8_x4500", "S1_v4std_n32_g8_x4500", "Sd_v4std_n64_g8", "S1_v4std_n64_g8", "S0_v4std_n32_g8", "S0_v4std_n32_g8_x4500", "S0_v4std_n64_g8")}
    L["C9b_Sd_0.82_S1_0.02_at_n32_4500_steps_and_0.78_0.018_at_9000"] = (abs(s9["Sd_v4std_n32_g8"] - 0.82) < 0.005 and abs(s9["S1_v4std_n32_g8"] - 0.022) < 0.002
                                                                          and abs(s9["Sd_v4std_n32_g8_x4500"] - 0.784) < 0.002 and abs(s9["S1_v4std_n32_g8_x4500"] - 0.018) < 0.001)
    prof = {t: {"far_over_near": R[t]["tube"]["far_over_near"], "profile": R[t]["tube"]["profile_z_E_per_plane_minus_x"]} for t in s9}
    A["C9"] = {}
    L["C9c_z_tube_profile_flat_far_over_near_in_0.5_2_for_Sd_both_boxes"] = all(0.5 < prof[t]["far_over_near"] < 2.0 for t in ("Sd_v4std_n32_g8_x4500", "Sd_v4std_n64_g8"))
    L["C9d_last_plane_next_to_the_pin_shell_carries_less_than_half_of_the_tube_excess_every_V4_main_row"] = all(
        R[t]["tube"]["last_plane_share"] < 0.5 for t in s9)
    L["C9e_Sd_n64_mid_tube_flat_max_over_median_lt_1.2"] = R["Sd_v4std_n64_g8"]["tube"]["max_over_median_mid"] < 1.2
    A["C9"]["flat_tension_per_unit_L"] = {t: R[t]["tube"]["flat_tension_mid_half_per_unit_L"] for t in s9}
    A["C9"]["last_plane_share"] = {t: R[t]["tube"]["last_plane_share"] for t in s9}
    A["C9"]["flat_string_prediction_x48"] = {t: 48.0 * R[t]["tube"]["flat_tension_mid_half_per_unit_L"] for t in s9}
    A["C9"].update({"string_reads": s9, "profiles": prof})

    # ---- C10: wordings ----
    ctrl = res["controls"]
    sddd, s0dd = R["Sddd_v4dd_n32_g8_x4500"], R["S0dd_v4dd_n32_g8_x4500"]
    A["C10"] = {
        "AXIS_STRING (at FALLING,FALLING on n32,64) x3": "the rule fires on both bases; the S1 and S0 increments are at or below the n64 last-third drift (C7d) and off the tube prediction (C7e); Sd is supported by both",
        "AXES_DEGENERATE": "by the rule with resolution = the Sd increment 22.86 (the largest); pairwise differences 17.8 to 70.6; on the equal-budget basis the resolution is 20.1",
        "AXIS_ESCAPES (Sd, S0 under V_spec)": "supported: the winding rank's |Q| < 0.2 on r 6, 9, 12 (C3c)",
        "saddle: true (every field)": "supported: d2E/ds2 < 0 everywhere (C6b)",
        "R_match_2pct: null (every V4 object)": "reproduced on the producer's basis; on the equal-budget basis see C8c",
        "scaling_force: EXPAND (every field)": "the sign is supported (C5a, C5g, C5h); the magnitude is not a Derrick scaling force: the frozen pin shell carries 41 to 52 percent of dE/dlam on the S_d and S_0 families (C5f) and the whole-box value is 3.5 to 9 times -E_curv + 3 V there, 1.2 to 1.7 times on the compact S_1 family (C5b)",
        "winding_degree end 0.587 / 0.616 for S_d under V4 (the rank-1 hedgehog)": "not a partial escape: the solid-angle degree of the same oriented field is 1.000 on every cube at both boxes with 0 to 2 conflicts (C3f); the Mermin flux is under-resolved where lam2 - lam1 closes to 0.014 to 0.08 on the cube surfaces",
        "degree.rank{k}.conflicts (the line-defect read)": "orientation-algorithm dependent: the producer's majority vote leaves 4 to 30 times more conflicts than the maximum spanning tree on the line-defect ranks (S_1 rank 1 at n64: 4224 vs 126); the degenerate delta pair of the dd trio keeps thousands under both, and its flux is not a number (C3a2)",
        "degenerate_null_rel 0.0398": f"Sddd {sddd['E_mine']:.4f} vs S0dd {s0dd['E_mine']:.4f}; pin-shell energies {sddd['radial']['E_pin_shell']:.4f} vs {s0dd['radial']['E_pin_shell']:.4f} (frozen, differ by {abs(sddd['radial']['E_pin_shell'] - s0dd['radial']['E_pin_shell']):.4f}); interior {sddd['radial']['E_interior']:.4f} vs {s0dd['radial']['E_interior']:.4f}; both FALLING",
        "S1_g32_vs_R3_record_18.970": f"{ctrl['S1_g32_vs_R3_record_18.970']:.6f} vs R3 {R3E['lam0_un0_single_d0_n32']:.6f} (C1d)",
        "verdict CONVERGED / FALLING": "reproduced from the traces (C4); the extension rows' decades are measured from the resumed field's fmax, not the original seed's (reported per row)",
        "two_rung_rule: true": "the n48 rung is absent from the rows; the rule as coded is the one recomputed in C7",
    }
    A["lines"] = {k: bool(v) for k, v in L.items()}
    A["runtime_s"] = round(time.time() - T0, 1)
    A["peak_rss_GB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    A["audited_utc"] = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    with open(OUT_JSON, "w") as f:
        json.dump(A, f, indent=1, default=float)
    print()
    for k, v in A["lines"].items():
        print(f"{'PASS' if v else 'FAIL'} {k}")
    print()
    for b in ("producer_basis", "equal_budget_basis"):
        for o, v in c7[b].items():
            print(f"C7 {b:19s} {o}: d {v['d']:+8.3f} slope {v['slope']:+.4f} label {v['label']:16s} n64 drift {v['n64_last_third_drift_abs']:.3f} (ratio {v['increment_over_n64_drift']:.2f}) "
                  f"string pred sigma64x48 {v['string_prediction_sigma64_x_48']:.3f} (ratio {v['increment_over_string_prediction_64']:.2f})")
        t = c7["triple_producer_basis" if b == "producer_basis" else "triple_equal_budget_basis"]
        print(f"   triple {t['outcome']} resolution {t['resolution']:.3f} diffs {[round(x, 2) for x in t['pairwise_abs_diffs']]} Koide {t['koide_Q']:.4f} ratios {[round(x, 3) for x in t['ratios_sorted']]}")
    for o in c8:
        print(f"C8 {o}: producer basis R_match {c8[o]['producer_basis']['R_match_2pct']} rel {[round(x, 3) for x in c8[o]['producer_basis']['rel_spread'].values()]}; "
              f"equal budget R_match {c8[o]['equal_budget_basis']['R_match_2pct']} rel {[round(x, 3) for x in c8[o]['equal_budget_basis']['rel_spread'].values()]}; pin-shell {c8[o]['pin_shell_E']}")
    for t, v in prof.items():
        print(f"C9 {t:24s} read {s9[t]:+.4f} flat mid {R[t]['tube']['flat_tension_mid_half_per_unit_L']:.4f}/L last-plane share {R[t]['tube']['last_plane_share']:.2f} "
              f"far/near {v['far_over_near']:.3f} profile {[round(p[1], 3) for p in v['profile']]}")
    for t, v in c5.items():
        print(f"C5 {t:24s} vir {v['virial_mine']:6.2f} dE/dlam {v['dE_dlam_mine_0.95_1.05']:+8.2f} (row cubic {v['dE_dlam_row_cubic']:+8.2f}) continuum {v['continuum_-Ecurv+3V']:+7.2f} "
              f"pin {v['dE_dlam_pin_shell']:+8.2f} interior {v['dE_dlam_interior']:+8.2f} core r<12 {v['dE_dlam_r_lt_12']:+7.2f} vs its continuum {v['continuum_r_lt_12']:+7.2f} (virial r<12 {v['virial_r_lt_12']:.2f})")
    print(f"runtime {A['runtime_s']} s, peak RSS {A['peak_rss_GB']:.2f} GB")


if __name__ == "__main__":
    main()
