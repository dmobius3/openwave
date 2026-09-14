"""M5.32 R21-1 to R21-4 ADVERSARIAL AUDIT: an independent attempt to refute the claims of the
run script m5_32_r21_1_runs.py and its collected results data/m5_32_r21_1_runs.json (key
"results", the per-row "rows"; the optional n64 row from data/m5_32_r21_1_runs_n64.json, merged at
collect). Every READ below is this file's own (its own one-sided stencils and curvature density,
its own V4 through the N spectrum, its own pin mask, its own tube / shell / radial sums, its own
dilation by map_coordinates at a different spline order, its own boost curvature and even
extrapolation, its own surface orientation and preimage-counting degree, its own norm profile, its
own root-law and product fits). The stack is imported ONLY for the certified objects the claims
are about: m5_32_r20_1_axes.energy_grad (the residual max |G| the gate is defined on) and
seed_axes (the seeds the pins hold), m5_21_3_a_4d as B3 (base_cfg, pin_shell, W1), and
m5_32_r3_ii_pair boost_at / conj (the dressing the frame curvature is defined along). No descent,
no L-BFGS: energy evaluations on the saved fields data/m5_32_r21_1/<tag>.npz (polished),
<tag>_fire.npz (the FIRE end) and data/m5_32_r20_1/*.npz (the R20 seeds) only.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4 per cell, eta = diag(-1, 1, 1, 1), N = M eta, code branch s = -1:
M_vac = diag(g, 1, delta, 0), the N spectrum q = (-g, 1, delta, 0). Box n32 L48 (h 1.5); the n64
row L96 (h 1.5). Coordinates x_i = (i - (n - 1) / 2) h, so the central 32^3 of the n64 box sits on
the n32 cell centers exactly (x_{i+16}(64) = x_i(32)).
    A_i = d_i M on the two one-sided branches (fwd: (M[i+1] - M[i]) / h on cells 0..n-2, zero on
          the last plane; bwd: (M[i] - M[i-1]) / h on cells 1..n-1, zero on the first plane)
    F_ij = A_i eta A_j - A_j eta A_i,  <X, Y>_eta = sum_ab eta_a eta_b X_ab Y_ab
    e_curv(x) = 4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta          E_curv = sum_x e_curv
    e_V(x)    = h^3 w sum_{p=1..4} (t_p - C_p)^2, t_p = sum_i lam_i(N)^p, C_p = sum_i q_i^p,
                w = W1 x w1s (W1 = B3.W1, w1s in 1, 5, 25); lam(N) = (-M_00, eigh(M_3x3)) when
                M_0i = 0 exactly (checked), eigvals(N) otherwise
    E = E_curv + V
The residual: G = the stack's gradient (energy_grad); max |G| over the FREE cells, this file's own
mask (all cells under B_free, the cells outside the outer ceil(1.6 / h) = 2 planes per face under
B_seed and B_far, 21952 of 32768 at n32); POLISHED iff max |G| < 1e-3.
The pin: under B_seed the shell cells of the polished and FIRE fields equal the seed's bitwise (the
seed rebuilt by seed_axes(cfg, lam) at the row's g and delta, or the saved R20 field the row started
from); under B_far the shell equals the central 32^3 of the R20 n64 end field (this file's own
slice [16:48]^3).
Degree of the rank-k eigenvector line field on the cube surface c0 +- k cells (c0 = n // 2, k =
half / h, half 6 / 9 / 12): the surface cells are oriented by a greedy Prim walk along the strongest
links |v . v'| over the surface graph; a conflict = a surface link whose oriented vectors disagree
(a line defect through the surface, no integer degree); the degree by PREIMAGE COUNTING: for 24
random unit directions u, the signed number of oriented spherical triangles (the two triangles of
each surface quad, outward-oriented) that contain u (all three of det(p, q, u), det(q, r, u),
det(r, p, u) of one sign, on u's hemisphere), the median over u and the spread max - min (zero
for a well-defined map); the Van Oosterom / Strackee solid-angle sum over the same triangles as
this file's second number. |degree| is the invariant (the sign is the start cell's convention).
Tube (string) read: E in the tube rho < 3 along z with 4 < |z| < L/2 - 2 minus the mean of the same
tubes along x and y, per unit length 2 (L/2 - 2 - 4) = 36; the flat-middle tension = the median
over the middle half of the per-plane excess (z-tube plane minus x-tube plane, both signs summed
and halved), per unit length h.
Shells: E(< R) = sum e over r < R; the face layer = sum e over the pin cells (pinned or not).
Norm: n3 = tr(M_3x3^2) on the axis cells rho < 0.75 h, 4 < |z| < L/2 - 2, against the vacuum's
1 + delta^2; crossings = axis cells with min(lam_1 - lam_0, lam_2 - lam_1) < 0.02.
Derrick: M_lam(x) = M(x / lam) about the box center by scipy map_coordinates at spline ORDER 2
(the producer used 3, with 1 as its check; this file also runs 1), edge-clamped; dE/dlam at 1 by
the central difference (E(1.05) - E(0.95)) / 0.1 and by the 4-point stencil over 0.9 / 0.95 /
1.05 / 1.1; the continuum value for a quartic-curvature + potential functional is -E_curv + 3 V;
virial = E_curv / V (3 at a Derrick equilibrium); r_half = the radius enclosing half the density;
R_* = r_half (E_curv / 3 V)^(1/4).
Frame: M -> Q M Q^T, Q = boost_at(cfg, 0, s) (eta-orthogonal, so the N spectrum and both
potentials are invariant: dV at round-off); c(s) = (E(+s) + E(-s) - 2 E(0)) / s^2 at s = 0.005,
0.01, 0.02, 0.05; the even extrapolation c0 = the least squares of c0 + c1 s^2 on 0.01 / 0.02 /
0.05 (the producer's three points, this file's own solve) and the Richardson pair
(4 c(0.005) - c(0.01)) / 3 as the check.
Fits: log |c| = log A + 2 log |x - x0| by least squares over (log A, x0) (a grid over x0 then
scipy least_squares); the producer minimized the max relative miss instead. The miss of a row =
|pred - |c|| / max(|c|, 0.2 median |c|). Product: |c| = C (g - g0)^2 (delta - delta0)^2 over the
eight rows with C the geometric-mean fit. Labels per the producer's docstring:
FRAME_NOT_PRODUCT iff the product max miss > 0.2; FRAME_TWO_ROOTS iff |g0 - 1| <= 0.3 and
|delta0 - 1| <= 0.15; FRAME_ONE_ROOT iff one of those or |g0| <= 0.3 (the pure g^2 law); else
FRAME_ROOTS_ELSEWHERE. The g 16 separation: (g - 1)^2 predicts c(32)/c(8) = (31/7)^2 = 19.61 and
c(16)/c(8) = (15/7)^2 = 4.59; g^2 predicts 16.0 and 4.0.
Outcomes (the producer's docstring, re-derived on this file's numbers):
    per axis  INTERIOR_SET iff E(< 12) agrees to 2 percent over the three boundaries and
              |degree| within 0.1 of 1 on r 9 under all three; BOUNDARY_SET iff the polished
              total under B_free or B_far differs from B_seed's by more than 5 percent or the
              winding is lost (|degree| < 0.5) under any; else BOUNDARY_UNDECIDED
    string    STRING_ESCAPES iff under B_free the tube read < 0.02 or the degree is lost;
              STRING_HELD iff within 30 percent of B_seed's with the degree kept; else UNDECIDED;
              the route where a winding is lost: EXCHANGE iff the axis norm min >= 0.9 vac and a
              crossing (gap < 0.02) exists; MELTING iff the norm min < 0.7 vac; else UNDECIDED
    box       AXIS_BOX_LIMITED iff dE/dlam < 0 and R_* > 12; AXIS_LOCALIZED iff R_* < 8 and
              |virial - 3| <= 0.9; else BOX_UNDECIDED
    triple    THREE_AXES_DISTINCT iff every pairwise |dE| of the three B_free energies exceeds
              three times the largest boundary difference of any axis; else AXES_DEGENERATE;
              Koide Q = sum E / (sum sqrt E)^2
    ladder    SCALE_IS_THE_POTENTIALS iff the S_1 w25 B_seed row has R_* < 8, |virial - 3| <= 0.9
              and |dE/dlam| <= 0.1 |dE/dlam of the S_1 w1 B_seed row|; SCALE_IS_THE_BOXS iff
              R_* > 12 or virial > 6; else SCALE_UNDECIDED
The results block was re-collected by the producer at 2026-09-14 22:31 UTC, while this audit was
running, with a surface-orientability reader (end_reads.degree_surface: the field oriented on the
cube surface by a maximum spanning tree, "undefined" where a surface link conflicts) and a box
outcome that reads the scan against the identity -E_curv + 3 V; the numeric reads did not change.
Every label comparison below is against the re-collected file; the new reader is itself checked
against this file's surface orientation (the conflict counts must agree exactly, C3g).
Convergence gate on a FIRE trace (the R20 gate): over the trace rows with acc >= 2/3 of the
accepted steps, rel = |E_last - E_first| / max(|E_last|, 1), decades = log10(fmax_seed / fmax_end);
CONVERGED iff rel < 1e-3 and decades >= 2, else FALLING if dE < 0 else RISING.

CHECKS (each a PASS line that can fail; every number printed):
 C1  energies: my E_total on every polished field vs row["E"] (rel 1e-10), on every FIRE field vs
     row["E_fire"]; my E_curv and V vs the row's parts; the seed energies vs descent.E0 / polish.E_in;
     M_0i = 0 on every field; the rows' roots and W1_eff equal mine; the n64 row merged verbatim
 C2  the residual gate: the stack's max |G| on my free mask vs polish.fmax_end (rel 1e-6) and the
     label POLISHED iff < 1e-3 on every row; the FIRE-end max |G| vs descent.fmax_end; the shell of
     every pinned field bitwise equal to its seed's; the B_far shell equal to the restricted n64
     field; the B_free shells moved; the residual spread across the boundary triple; the rows on
     which the polish RAISED max |G|
 C3  the boundary triple: my E, E(< 12), |degree| on r 9 (preimage counting) per boundary; the
     labels re-derived; my degrees vs the producer's solid-angle degrees (abs 0.02 on the cubes
     where my surface map is conflict-free); on every cube with conflicts, the solid-angle integer
     under TWO surface orientations (greedy Prim, plain BFS from another start): a closed
     triangulated map has an integer solid-angle sum for any vertex signs, so an integer that
     changes with the orientation is a sign choice, not a degree (C3d2, C3d3); the S_d and S_0
     winding ranks on r 9 (C3e), the S_0 B_free rank-1 carrier (C3f)
 C4  the string: my tube read on S_d under B_free / B_seed (abs 1e-8 vs the rows), the flat-middle
     tension, the face layer under B_free; the label under the tube read and under the flat-middle
     read; the degree-kept clause tested on its own (C4b2); the norm-rule escape routes on the rows
     the producer flagged and the rank-transfer carrier (C4e, C4e2); the axis-norm caveat
 C5  the box: my dE/dlam (order 2, both stencils; order 1 as the check) vs the row's order-3 value
     (rel 0.15), the sign, the continuum -E_curv + 3 V and its sign; the slope split over the face
     layer / the interior / the core r < 12 (C5g, C5h); my virial (rel 1e-8), r_half (rel 0.01),
     R_* (rel 0.01); the labels
 C6  the triple: pairwise differences, the resolution, the label, Koide Q and the ratios
 C7  the W1 ladder: my R_*, virial, dE/dlam, shell share per rung; the label; the S_1 w25 seed vs
     free comparison
 C8  the frame: my c(s) at 0.01 / 0.02 / 0.05 vs the row's (rel 1e-6), dV at round-off, my c0 vs
     the row's c0_extrapolated (rel 1e-4) and vs my Richardson pair; my root fits (g0, delta0), the
     product fit's max miss, the label; the g 16 separation; the fit on the POLISHED g rows alone
     and on the FALLING g rows alone
 C9  the traces: every polish trace monotone non-increasing, the iterations and evaluations; every
     FIRE verdict reproduced from its trace by the R20 gate
 C10 the wordings in results (listed with the numbers that do or do not support them)
Out: ../data/m5_32_r21_1_audit.json. Runtime and peak RSS printed.
Usage: python3 m5_32_r21_1_audit.py   (a few minutes single-threaded, under 3 GB)
"""

from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import gc  # noqa: E402
import heapq  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
RUNS_JSON = os.path.join(DATA, "m5_32_r21_1_runs.json")
N64_JSON = os.path.join(DATA, "m5_32_r21_1_runs_n64.json")
R20_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")
NPZ = os.path.join(DATA, "m5_32_r21_1")
R20_NPZ = os.path.join(DATA, "m5_32_r20_1")
OUT_JSON = os.path.join(DATA, "m5_32_r21_1_audit.json")
T0 = time.time()
LINES = {}


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


R20 = _load(
    "m5_32_r20_1_axes", "m5_32_r20_1_axes.py"
)  # energy_grad (the residual's definition), seed_axes (the pinned seed)
B3 = R20.B3  # base_cfg, pin_shell, W1
R3 = R20.R3  # boost_at, conj (the dressing)
LAG = R20.LAG  # default_params (the stack's parameter dict for energy_grad)

ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(ETA_D)
SGN = np.einsum("a,b->ab", ETA_D, ETA_D)
W1 = float(B3.W1)
HALVES = (6.0, 9.0, 12.0)
PROFILE_R = (3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0)
LAM_OF = {
    "S1": lambda d: (1.0, d, 0.0),
    "Sd": lambda d: (d, 1.0, 0.0),
    "S0": lambda d: (0.0, 1.0, d),
}
WRANK = {"S1": 2, "Sd": 1, "S0": 0}
FRAME_S = (0.005, 0.01, 0.02, 0.05)


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def log(msg):
    print(f"[{time.time() - T0:6.1f}s {rss_gb():4.2f}GB] {msg}", flush=True)


def line(key, ok, detail):
    LINES[key] = {"pass": bool(ok), "detail": detail}
    print(f"{'PASS' if ok else 'FAIL'} {key}: {detail}", flush=True)


def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


# ================= own energy layer =================
def own_coords(n, h):
    x = (np.arange(n) - 0.5 * (n - 1)) * h
    return np.meshgrid(x, x, x, indexing="ij")


def own_pin(n, h, depth=1.6):
    w = int(np.ceil(depth / h))
    P = np.ones((n, n, n), dtype=bool)
    P[w : n - w, w : n - w, w : n - w] = False
    return P


def own_d1(f, ax, h, side):
    out = np.zeros_like(f)
    hi = [slice(None)] * f.ndim
    lo = [slice(None)] * f.ndim
    hi[ax], lo[ax] = slice(1, None), slice(0, -1)
    d = (f[tuple(hi)] - f[tuple(lo)]) / h
    tgt = [slice(None)] * f.ndim
    tgt[ax] = slice(0, -1) if side == "fwd" else slice(1, None)
    out[tuple(tgt)] = d
    return out


def own_curv_density(M, h):
    """4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta per cell, F_ij = A_i eta A_j - A_j eta A_i."""
    e = np.zeros(M.shape[:3])
    for side in ("fwd", "bwd"):
        A = [own_d1(M, ax, h, side) for ax in range(3)]
        AE = [
            a * ETA_D[None, None, None, None, :] for a in A
        ]  # A eta (eta diagonal: scale the columns)
        for i in range(3):
            for j in range(i + 1, 3):
                F = AE[i] @ A[j] - AE[j] @ A[i]
                e += 0.5 * 4.0 * np.einsum("...ab,ab->...", F * F, SGN)
        del A, AE, F
    return h**3 * e


def own_spectrum(M):
    Mf = M.reshape(-1, 4, 4)
    m0i = float(np.max(np.abs(Mf[:, 0, 1:])))
    if m0i == 0.0:
        lam3 = np.linalg.eigvalsh(Mf[:, 1:, 1:])
        return np.concatenate([-Mf[:, 0, 0][:, None], lam3], axis=1), m0i, "block"
    ev = np.linalg.eigvals(Mf @ ETA)
    if np.max(np.abs(ev.imag)) > 1e-9 * max(1.0, float(np.max(np.abs(ev)))):
        raise RuntimeError("complex N spectrum")
    return ev.real, m0i, "eigvals"


def own_v4_density(M, h, q, w):
    lam, m0i, how = own_spectrum(M)
    acc = np.zeros(lam.shape[0])
    for p in range(1, 5):
        acc += (np.sum(lam**p, axis=1) - sum(qi**p for qi in q)) ** 2
    return (h**3 * w * acc).reshape(M.shape[:3]), m0i, how


def own_energy(M, cfg, q, w):
    h = cfg["h"]
    ec = own_curv_density(M, h)
    ev, m0i, how = own_v4_density(M, h, q, w)
    E = {
        "E_curv": float(ec.sum()),
        "V": float(ev.sum()),
        "E_total": float(ec.sum() + ev.sum()),
        "max_abs_M0i": m0i,
        "spectrum_how": how,
    }
    return E, ec + ev, ec, ev


# ================= own reads =================
def surface_index(n, c0, k):
    lo, hi = c0 - k, c0 + k
    m = np.zeros((n, n, n), dtype=bool)
    sl = slice(lo, hi + 1)
    m[lo, sl, sl] = m[hi, sl, sl] = True
    m[sl, lo, sl] = m[sl, hi, sl] = True
    m[sl, sl, lo] = m[sl, sl, hi] = True
    cells = np.argwhere(m)
    return cells, lo, hi


def own_degree(vec, gap, c0, k, rng, ndir=24):
    """(median preimage count, spread, conflicts, solid-angle degree, cells, min gap on the surface, min gap at the
    conflict links) of the rank eigenvector on the cube surface c0 +- k; gap = the rank's smaller eigenvalue gap per
    cell (a conflict link whose ends sit at a closed gap is a discontinuity of the rank-labeled eigenvector, an
    exchange with the neighboring rank, not a jump of a continuous field)."""
    n = vec.shape[0]
    cells, lo, hi = surface_index(n, c0, k)
    key = {tuple(c): i for i, c in enumerate(cells)}
    V = vec[cells[:, 0], cells[:, 1], cells[:, 2]]
    V = V / np.linalg.norm(V, axis=1)[:, None]
    Gs = gap[cells[:, 0], cells[:, 1], cells[:, 2]]
    nb = [[] for _ in range(len(cells))]
    for i, c in enumerate(cells):
        for ax in range(3):
            for dd in (-1, 1):
                t = list(c)
                t[ax] += dd
                j = key.get(tuple(t))
                if j is not None:
                    nb[i].append(j)
    sign = np.zeros(len(cells))
    start = 0
    sign[start] = 1.0
    heap = [(-abs(float(V[start] @ V[j])), start, j) for j in nb[start]]
    heapq.heapify(heap)
    while heap:
        _, i, j = heapq.heappop(heap)
        if sign[j] != 0.0:
            continue
        sign[j] = sign[i] * (1.0 if float(V[i] @ V[j]) >= 0.0 else -1.0)
        for j2 in nb[j]:
            if sign[j2] == 0.0:
                heapq.heappush(heap, (-abs(float(V[j] @ V[j2])), j, j2))
    if np.any(sign == 0.0):
        raise RuntimeError("surface graph not connected")
    VO = V * sign[:, None]
    conflicts = 0
    gap_conf = []
    for i in range(len(cells)):
        for j in nb[i]:
            if j > i and float(VO[i] @ VO[j]) < 0.0:
                conflicts += 1
                gap_conf.append(min(Gs[i], Gs[j]))
    # a second orientation: a plain breadth-first walk from the LAST cell (weakest-link-blind); on an orientable surface field
    # it reproduces the first up to a global sign, on a torn one the two integers below can differ
    sign2 = np.zeros(len(cells))
    start2 = len(cells) - 1
    sign2[start2] = 1.0
    queue = [start2]
    while queue:
        i = queue.pop(0)
        for j in nb[i]:
            if sign2[j] == 0.0:
                sign2[j] = sign2[i] * (1.0 if float(V[i] @ V[j]) >= 0.0 else -1.0)
                queue.append(j)
    VO2 = V * sign2[:, None]
    tris = []
    for ax in range(3):
        o1, o2 = [a for a in range(3) if a != ax]
        for plane, is_lo in ((hi, False), (lo, True)):
            flip = is_lo ^ (ax == 1)
            for a in range(lo, hi):
                for b in range(lo, hi):

                    def at(u, v_):
                        idx = [0, 0, 0]
                        idx[ax], idx[o1], idx[o2] = plane, u, v_
                        return key[tuple(idx)]

                    p, q_, r_, s_ = at(a, b), at(a + 1, b), at(a + 1, b + 1), at(a, b + 1)
                    for t in ((p, q_, r_), (p, r_, s_)):
                        tris.append((t[0], t[2], t[1]) if flip else t)
    T = np.array(tris)
    P, Q, R = VO[T[:, 0]], VO[T[:, 1]], VO[T[:, 2]]
    cPQ, cQR, cRP = np.cross(P, Q), np.cross(Q, R), np.cross(R, P)
    degs = []
    for _ in range(ndir):
        u = rng.standard_normal(3)
        u /= np.linalg.norm(u)
        s1, s2, s3 = np.sign(cPQ @ u), np.sign(cQR @ u), np.sign(cRP @ u)
        inside = (s1 == s2) & (s2 == s3) & ((P + Q + R) @ u > 0.0)
        degs.append(int(np.sum(s1[inside])))

    def solid_of(W):
        P_, Q_, R_ = W[T[:, 0]], W[T[:, 1]], W[T[:, 2]]
        num = np.einsum("na,na->n", P_, np.cross(Q_, R_))
        den = (
            1.0
            + np.einsum("na,na->n", P_, Q_)
            + np.einsum("na,na->n", Q_, R_)
            + np.einsum("na,na->n", R_, P_)
        )
        return float(np.sum(2.0 * np.arctan2(num, den)) / (4.0 * np.pi))

    solid, solid2 = solid_of(VO), solid_of(VO2)
    return (
        float(np.median(degs)),
        int(max(degs) - min(degs)),
        conflicts,
        solid,
        len(cells),
        float(np.min(Gs)),
        float(min(gap_conf)) if gap_conf else None,
        solid2,
    )


def degree_reads(M, cfg, rng):
    n, h = cfg["n"], cfg["h"]
    c0 = n // 2
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    g01, g12 = lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1]
    gaps = {0: g01, 1: np.minimum(g01, g12), 2: g12}
    out = {}
    for k in range(3):
        rec = {}
        for hv in HALVES:
            d, sp, cf, sa, nc, gmin, gconf, sa2 = own_degree(
                vec[..., :, k], gaps[k], c0, int(round(hv / h)), rng
            )
            rec[f"{hv:g}"] = {
                "preimage": d,
                "spread": sp,
                "conflicts": cf,
                "solid_angle": sa,
                "solid_angle_bfs_orientation": sa2,
                "cells": nc,
                "surface_gap_min": gmin,
                "gap_min_at_conflicts": gconf,
                "orientation_dependent": bool(abs(abs(sa) - abs(sa2)) > 0.5),
            }
        out[f"rank{k}"] = rec
    return out, lam


def tube_reads(e, cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = own_coords(n, h)
    zmax = 0.5 * L - 2.0
    length = 2.0 * (zmax - 4.0)
    tubes = {}
    for lab, perp, along in (
        ("z", np.sqrt(X * X + Y * Y), Z),
        ("x", np.sqrt(Y * Y + Z * Z), X),
        ("y", np.sqrt(X * X + Z * Z), Y),
    ):
        m = (perp < 3.0) & (np.abs(along) > 4.0) & (np.abs(along) < zmax)
        tubes[lab] = float(np.sum(e[m]))
    read = (tubes["z"] - 0.5 * (tubes["x"] + tubes["y"])) / length
    rho = np.sqrt(X * X + Y * Y)
    zs = np.unique(np.abs(Z))
    zs = zs[(zs > 4.0) & (zs < zmax)]
    ex = []
    for zv in zs:
        mz = (rho < 3.0) & (np.abs(np.abs(Z) - zv) < 1e-9)
        mx = (np.sqrt(Y * Y + Z * Z) < 3.0) & (np.abs(np.abs(X) - zv) < 1e-9)
        ex.append(0.5 * (float(np.sum(e[mz])) - float(np.sum(e[mx]))))
    q1, q3 = len(ex) // 4, 3 * len(ex) // 4
    mid = ex[q1:q3]
    return {
        "tubes": tubes,
        "length": length,
        "string_read": float(read),
        "per_plane_excess": [float(v) for v in ex],
        "planes_abs_z": [float(v) for v in zs],
        "flat_mid_tension": float(np.median(mid) / h),
        "near_quarter_tension": float(np.median(ex[:q1]) / h),
        "far_quarter_tension": float(np.median(ex[q3:]) / h),
    }


def radial_reads(e, cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = own_coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    pin = own_pin(n, h)
    outer = np.zeros((n, n, n), dtype=bool)
    outer[0], outer[-1], outer[:, 0], outer[:, -1], outer[:, :, 0], outer[:, :, -1] = (
        True,
        True,
        True,
        True,
        True,
        True,
    )
    o = np.argsort(r.ravel())
    cum = np.cumsum(e.ravel()[o])
    r_half = float(r.ravel()[o][np.searchsorted(cum, 0.5 * cum[-1])])
    edges = np.arange(8.0, 16.0 + 0.5 * h, h)
    rc, dens = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        if m.sum():
            rc.append(0.5 * (a + b))
            dens.append(float(np.mean(e[m])))
    return {
        "E_lt_R": {f"{R:g}": float(np.sum(e[r < R])) for R in PROFILE_R if R <= 0.5 * L},
        "E_pin_shell": float(np.sum(e[pin])),
        "E_outer_plane": float(np.sum(e[outer])),
        "pin_cells": int(pin.sum()),
        "r_half": r_half,
        "tail_slope": float(np.polyfit(rc, dens, 1)[0]),
    }


def norm_reads(M, lam, cfg):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = own_coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    ax = (rho < 0.75 * h) & (np.abs(Z) > 4.0) & (np.abs(Z) < 0.5 * L - 2.0)
    nrm = np.einsum("...ij,...ji->...", M[..., 1:, 1:], M[..., 1:, 1:])
    vac = 1.0 + cfg["delta"] ** 2
    gmin = np.minimum(lam[ax][:, 1] - lam[ax][:, 0], lam[ax][:, 2] - lam[ax][:, 1])
    return {
        "vacuum": vac,
        "axis_cells": int(ax.sum()),
        "axis_norm_min_rel": float(np.min(nrm[ax]) / vac),
        "axis_norm_mean_rel": float(np.mean(nrm[ax]) / vac),
        "axis_crossings_gap_lt_0.02": int(np.sum(gmin < 0.02)),
        "axis_gap_min": float(np.min(gmin)),
    }


def dilate(M, lam, order):
    n = M.shape[0]
    c = (n - 1) / 2.0
    src = (np.arange(n) - c) / lam + c
    X, Y, Z = np.meshgrid(src, src, src, indexing="ij")
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()])
    out = np.empty_like(M)
    for a in range(4):
        for b in range(a, 4):
            v = map_coordinates(M[..., a, b], coords, order=order, mode="nearest").reshape(n, n, n)
            out[..., a, b] = v
            out[..., b, a] = v
    return out


def derrick_reads(M, cfg, q, w, E1, r_half, ec, ev):
    """dE/dlam by dilation, and its split over the face layer (the pin cells), the interior, and the core r < 12: a dilation
    on a truncated box carries energy across the faces, so the split says whether the sign is the object's or the box's.
    """
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = own_coords(n, h)
    core = np.sqrt(X * X + Y * Y + Z * Z) < 12.0
    pin = own_pin(n, h)
    out = {
        "virial": E1["E_curv"] / E1["V"],
        "continuum_dE_dlam": -E1["E_curv"] + 3.0 * E1["V"],
        "r_half": r_half,
        "R_star": float(r_half * (E1["E_curv"] / (3.0 * E1["V"])) ** 0.25),
        "continuum_dE_dlam_r_lt_12": float(-ec[core].sum() + 3.0 * ev[core].sum()),
        "virial_r_lt_12": float(ec[core].sum() / ev[core].sum()),
    }
    for order in (2, 1):
        Es, dens = {}, {}
        for lam in (0.9, 0.95, 1.05, 1.1):
            Ed, ed, _, _ = own_energy(dilate(M, lam, order), cfg, q, w)
            Es[lam] = Ed["E_total"]
            dens[lam] = ed
        out[f"order{order}"] = {
            "dE_dlam_central_0.05": (Es[1.05] - Es[0.95]) / 0.1,
            "dE_dlam_4pt": (-Es[1.1] + 8.0 * Es[1.05] - 8.0 * Es[0.95] + Es[0.9]) / (12.0 * 0.05),
            "E_of_lambda": {f"{k:g}": v for k, v in Es.items()},
            "dE_dlam_face_layer": float((dens[1.05][pin].sum() - dens[0.95][pin].sum()) / 0.1),
            "dE_dlam_interior": float((dens[1.05][~pin].sum() - dens[0.95][~pin].sum()) / 0.1),
            "dE_dlam_r_lt_12": float((dens[1.05][core].sum() - dens[0.95][core].sum()) / 0.1),
            "dE_dlam_r_ge_12": float((dens[1.05][~core].sum() - dens[0.95][~core].sum()) / 0.1),
        }
        del dens
    return out


def frame_reads(M, cfg, q, w, E0):
    out = {"c": {}, "dV_rel": {}}
    for sv in FRAME_S:
        Ep = own_energy(R3.conj(R3.boost_at(cfg, 0.0, sv)[0], M), cfg, q, w)[0]
        Em = own_energy(R3.conj(R3.boost_at(cfg, 0.0, -sv)[0], M), cfg, q, w)[0]
        out["c"][f"{sv:g}"] = float((Ep["E_total"] + Em["E_total"] - 2.0 * E0["E_total"]) / sv**2)
        out["dV_rel"][f"{sv:g}"] = float(
            max(abs(Ep["V"] - E0["V"]), abs(Em["V"] - E0["V"])) / abs(E0["V"])
        )
        out.setdefault("M0i_plus", {})[f"{sv:g}"] = Ep["max_abs_M0i"]
    s3 = np.array([0.01, 0.02, 0.05])
    c3 = np.array([out["c"]["0.01"], out["c"]["0.02"], out["c"]["0.05"]])
    A = np.stack([np.ones(3), s3**2], axis=1)
    coef = np.linalg.solve(A.T @ A, A.T @ c3)
    out["c0_lsq_3pt"] = float(coef[0])
    out["c1_lsq_3pt"] = float(coef[1])
    out["c0_richardson_0.005_0.01"] = float((4.0 * out["c"]["0.005"] - out["c"]["0.01"]) / 3.0)
    return out


def fit_root_law(xs, cs):
    """log |c| = log A + 2 log |x - x0|, least squares over (log A, x0): a grid over x0 then a least_squares polish."""
    xs, ac = np.asarray(xs, float), np.abs(np.asarray(cs, float))
    lc = np.log(ac)

    def sse(x0):
        w = 2.0 * np.log(np.maximum(np.abs(xs - x0), 1e-9))
        logA = float(np.mean(lc - w))
        return float(np.sum((lc - logA - w) ** 2)), logA

    grid = np.linspace(-2.0 * xs.max(), 2.0 * xs.max(), 40001)
    grid = grid[np.min(np.abs(grid[:, None] - xs[None, :]), axis=1) > 1e-3]
    best = min(grid, key=lambda x0: sse(x0)[0])
    res = least_squares(
        lambda p: lc - p[0] - 2.0 * np.log(np.maximum(np.abs(xs - p[1]), 1e-9)),
        x0=[sse(best)[1], best],
    )
    logA, x0 = (
        (float(res.x[0]), float(res.x[1]))
        if res.success and res.cost * 2 <= sse(best)[0] + 1e-12
        else (sse(best)[1], float(best))
    )
    pred = np.exp(logA) * (xs - x0) ** 2
    miss = np.abs(pred - ac) / np.maximum(ac, 0.2 * np.median(ac))
    return {
        "A": float(np.exp(logA)),
        "x0": x0,
        "max_rel_miss": float(np.max(miss)),
        "rms_log_resid": float(np.sqrt(sse(x0)[0] / len(xs))),
        "rel_miss_by_row": {f"{x:g}": float(m) for x, m in zip(xs, miss)},
    }


def product_fit(rows, g0, d0):
    ac = np.array([abs(c) for _, _, c in rows])
    basis = np.array([(g - g0) ** 2 * (d - d0) ** 2 for g, d, _ in rows])
    C = float(np.exp(np.mean(np.log(ac) - np.log(np.maximum(basis, 1e-300)))))
    miss = np.abs(C * basis - ac) / np.maximum(ac, 0.2 * np.median(ac))
    return {
        "C": C,
        "g0": g0,
        "delta0": d0,
        "max_rel_miss": float(np.max(miss)),
        "rel_miss_by_row": {f"g{g:g}_d{d:g}": float(m) for (g, d, _), m in zip(rows, miss)},
    }


def frame_label(g0, d0, prod_miss):
    if prod_miss > 0.2:
        return "FRAME_NOT_PRODUCT"
    g_at_1, d_at_1, g_at_0 = abs(g0 - 1.0) <= 0.3, abs(d0 - 1.0) <= 0.15, abs(g0) <= 0.3
    if g_at_1 and d_at_1:
        return "FRAME_TWO_ROOTS"
    if g_at_1 or d_at_1 or g_at_0:
        return "FRAME_ONE_ROOT" + (" (the pure g^2 law)" if g_at_0 else "")
    return "FRAME_ROOTS_ELSEWHERE"


def fire_gate(des):
    tr = des.get("trace", [])
    n_acc = des.get("accepted", 0)
    q = [r for r in tr if r["acc"] >= (2.0 / 3.0) * n_acc] if n_acc else []
    if len(q) < 2:
        return "budget (trace too short)", None, None
    dE = q[-1]["E"] - q[0]["E"]
    rl = abs(dE) / max(abs(q[-1]["E"]), 1.0)
    dec = float(np.log10(des["fmax_seed"] / max(des["fmax_end"], 1e-300)))
    lab = "CONVERGED" if (rl < 1e-3 and dec >= 2.0) else ("FALLING" if dE < 0 else "RISING")
    return lab, float(rl), dec


def seed_field(r, cfg):
    sf = r.get("seed_from")
    if sf == "n64restrict":
        M64 = np.load(os.path.join(R20_NPZ, f"{r['obj']}_v4std_n64_g8.npz"))["M"]
        return M64[16:48, 16:48, 16:48].copy(), "restricted n64 end field"
    if sf == "r20_x4500":
        return (
            np.load(os.path.join(R20_NPZ, f"{r['obj']}_v4std_n32_g8_x4500.npz"))["M"],
            "R20 x4500 field",
        )
    if sf == "r20_g32":
        return np.load(os.path.join(R20_NPZ, "S1_v4std_n32_g32.npz"))["M"], "R20 g32 field"
    if sf == "r20_n64":
        return (
            np.load(os.path.join(R20_NPZ, f"{r['obj']}_v4std_n64_g8.npz"))["M"],
            "R20 n64 end field",
        )
    return R20.seed_axes(cfg, LAM_OF[r["obj"]](r["delta"])), "seed_axes rebuilt"


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return str(o)
    return o


# ================= the audit =================
def main():
    J = json.load(open(RUNS_JSON))
    rows = J["rows"]
    res = J["results"]
    J64 = json.load(open(N64_JSON))["rows"]
    R20J = json.load(open(R20_JSON))["rows"]
    A = {"task": "M5.32 R21-1 to R21-4 audit", "W1": W1, "rows": {}}
    rng = np.random.default_rng(2121)
    order = sorted(rows, key=lambda t: (rows[t]["n"], t))
    log(
        f"{len(rows)} rows (n32 {sum(1 for t in rows if rows[t]['n'] == 32)}, n64 {sum(1 for t in rows if rows[t]['n'] == 64)}), W1 {W1:.12g}"
    )

    # C1 preliminaries: the n64 row merged verbatim, the roots and weights
    n64_ok = (
        all(
            t in rows
            and rows[t]["E"] == r["E"]
            and rows[t]["polish"]["fmax_end"] == r["polish"]["fmax_end"]
            for t, r in J64.items()
        )
        and len(J64) == 1
    )
    line(
        "C1g_n64_row_merged_verbatim_from_the_side_pool_file",
        n64_ok,
        f"{list(J64)} E {[r['E'] for r in J64.values()]}",
    )

    frame_rows = {}
    for tag in order:
        r = rows[tag]
        n, L, g, d, w1s = r["n"], r["L"], r["g"], r["delta"], r["w1s"]
        cfg = B3.base_cfg(s=-1.0, g=g, n=n, L=float(L), delta=d)
        p = LAG.default_params(s=-1.0, g=g, delta=d)
        q = (-g, 1.0, d, 0.0)
        w = W1 * w1s
        pot = None if w1s == 1.0 else ("v4", q, w)
        pinned = r["bnd"] != "free"
        pin = own_pin(n, cfg["h"])
        free = ~pin if pinned else np.ones((n, n, n), dtype=bool)
        rec = {
            "n": n,
            "obj": r["obj"],
            "bnd": r["bnd"],
            "g": g,
            "delta": d,
            "w1s": w1s,
            "steps": r["steps"],
            "label_row": r["label"],
            "roots_row": r["roots"],
            "roots_mine": list(q),
            "W1_eff_row": r["W1_eff"],
            "W1_eff_mine": w,
            "pin_mask_equals_B3": bool(np.array_equal(pin, B3.pin_shell(n, cfg["h"]))),
            "free_cells_mine": int(free.sum()),
            "free_cells_row": r["polish"]["free_cells"],
        }
        Mp = np.load(os.path.join(NPZ, f"{tag}.npz"))["M"]
        E1, e, ec, ev = own_energy(Mp, cfg, q, w)
        rec.update(
            {
                "E_row": r["E"],
                "E_mine": E1["E_total"],
                "rel_E": rel(E1["E_total"], r["E"]),
                "rel_curv": rel(E1["E_curv"], r["end_reads"]["energy"]["E_curv"]),
                "rel_V": rel(E1["V"], r["end_reads"]["energy"]["V"]),
                "max_abs_M0i": E1["max_abs_M0i"],
                "max_abs_M00_minus_g": float(np.max(np.abs(Mp[..., 0, 0] - g))),
                "max_abs_M00_minus_g_on_shell": float(np.max(np.abs(Mp[pin][:, 0, 0] - g))),
                "E_curv_mine": E1["E_curv"],
                "V_mine": E1["V"],
            }
        )
        # the seed and the pin
        M0, seed_kind = seed_field(r, cfg)
        E0 = own_energy(M0, cfg, q, w)[0]
        rec["seed_kind"] = seed_kind
        rec["E0_mine"] = E0["E_total"]
        rec["E0_row"] = r["descent"]["E0"] if r["steps"] > 0 else r["polish"]["E_in"]
        rec["E0_seed_reads_row"] = r["seed_reads"]["energy"]["E_total"]
        rec["rel_E0"] = max(
            rel(E0["E_total"], rec["E0_row"]), rel(E0["E_total"], rec["E0_seed_reads_row"])
        )
        rec["shell_bitwise_equal_seed_polished"] = bool(np.array_equal(Mp[pin], M0[pin]))
        rec["shell_max_dev_from_seed_polished"] = float(np.max(np.abs(Mp[pin] - M0[pin])))
        if r.get("seed_from") in ("r20_x4500", "r20_g32"):
            Ms = R20.seed_axes(cfg, LAM_OF[r["obj"]](d))
            rec["shell_bitwise_equal_seed_axes_rebuilt"] = bool(np.array_equal(Mp[pin], Ms[pin]))
            del Ms
        del M0
        # the FIRE end
        fpath = os.path.join(NPZ, f"{tag}_fire.npz")
        if os.path.exists(fpath):
            Mf = np.load(fpath)["M"]
            Ef = own_energy(Mf, cfg, q, w)[0]
            rec["E_fire_row"] = r["E_fire"]
            rec["E_fire_mine"] = Ef["E_total"]
            rec["rel_E_fire"] = rel(Ef["E_total"], r["E_fire"])
            rec["fire_max_abs_M0i"] = Ef["max_abs_M0i"]
            Mfs, _ = seed_field(r, cfg)
            rec["shell_bitwise_equal_seed_fire"] = bool(np.array_equal(Mf[pin], Mfs[pin]))
            rec["shell_max_dev_from_seed_fire"] = float(np.max(np.abs(Mf[pin] - Mfs[pin])))
            del Mfs
            if n == 32:
                _, Gf, _ = R20.energy_grad(Mf, cfg, p, pot)
                rec["fire_fmax_free_mine"] = float(np.max(np.abs(Gf[free])))
                rec["fire_fmax_row"] = r["descent"]["fmax_end"]
                rec["rel_fire_fmax"] = rel(rec["fire_fmax_free_mine"], rec["fire_fmax_row"])
                del Gf
            del Mf
            gc.collect()
        # the residual gate on the polished field (the stack's gradient, my mask); the n64 gradient exceeds the memory budget and is not run
        if n == 32:
            Eg, G, info = R20.energy_grad(Mp, cfg, p, pot)
            rec["E_stack"] = float(Eg)
            rec["rel_E_stack_vs_mine"] = rel(E1["E_total"], float(Eg))
            rec["fmax_free_mine"] = float(np.max(np.abs(G[free])))
            rec["fmax_pinned_cells"] = float(np.max(np.abs(G[pin]))) if pinned else None
            rec["fmax_row"] = r["polish"]["fmax_end"]
            rec["rel_fmax"] = rel(rec["fmax_free_mine"], rec["fmax_row"])
            rec["label_mine"] = (
                "POLISHED"
                if rec["fmax_free_mine"] < 1e-3
                else f"FALLING (at max |G| = {rec['fmax_free_mine']:.2e})"
            )
            rec["label_agrees"] = rec["label_mine"].split(" ")[0] == r["label"].split(" ")[0]
            rec["polish_raised_fmax"] = bool(r["polish"]["fmax_end"] > r["polish"]["fmax_in"])
            rec["polish_fmax_in"] = r["polish"]["fmax_in"]
            del G
            gc.collect()
        # the reads
        dg, lam = degree_reads(Mp, cfg, rng)
        rec["degree"] = dg
        rec["degree_row_solid_angle"] = {
            k: v["degree"] for k, v in r["end_reads"]["solid_angle"].items()
        }
        devs = []
        for k in range(3):
            for hv in HALVES:
                mine = dg[f"rank{k}"][f"{hv:g}"]
                if (
                    mine["conflicts"] == 0
                    and mine["spread"] == 0
                    and f"{hv:g}" in rec["degree_row_solid_angle"][f"rank{k}"]
                ):
                    devs.append(
                        abs(
                            abs(mine["preimage"])
                            - abs(rec["degree_row_solid_angle"][f"rank{k}"][f"{hv:g}"])
                        )
                    )
        rec["degree_abs_dev_vs_row_orientable_cubes"] = max(devs) if devs else None
        rec["degree_orientable_cubes"] = len(devs)
        wr = WRANK[r["obj"]]
        rec["winding_rank"] = wr
        rec["deg9_winding_mine"] = dg[f"rank{wr}"]["9"]["preimage"]
        rec["deg9_winding_row"] = r["end_reads"]["solid_angle"][f"rank{wr}"]["degree"]["9"]
        rec["deg9_by_rank_mine"] = {f"rank{k}": dg[f"rank{k}"]["9"]["preimage"] for k in range(3)}

        # the producer's re-collected surface reader (added to the runs JSON at 22:31 UTC, mid-audit): its surface conflict counts
        # and its kept / lost / undefined readings against mine, rank by rank and cube by cube
        def my_reading(d):
            if d["conflicts"] > 0:
                return "undefined"
            return (
                "kept"
                if abs(abs(d["preimage"]) - 1.0) < 0.1
                else ("lost" if abs(d["preimage"]) < 0.5 else "fractional")
            )

        ds = r["end_reads"].get("degree_surface")
        if ds is not None:
            cmp = []
            for k in range(3):
                for hv in ("6", "9", "12"):
                    if hv in ds[f"rank{k}"]:
                        mine, theirs = dg[f"rank{k}"][hv], ds[f"rank{k}"][hv]
                        cmp.append(
                            {
                                "rank": k,
                                "cube": hv,
                                "conflicts_mine": mine["conflicts"],
                                "conflicts_row": theirs["surface_conflicts"],
                                "reading_mine": my_reading(mine),
                                "reading_row": theirs["reading"].split(" ")[0],
                            }
                        )
            rec["surface_reader_vs_mine"] = cmp
            rec["surface_conflicts_all_equal"] = all(
                c["conflicts_mine"] == c["conflicts_row"] for c in cmp
            )
            rec["surface_readings_all_equal"] = all(
                c["reading_mine"] == c["reading_row"] for c in cmp
            )
            rec["surface_reader_cubes"] = len(cmp)
        tb = tube_reads(e, cfg)
        rec["tube"] = tb
        rec["string_read_row"] = r["end_reads"]["string_tension_read"]
        rec["string_read_dev"] = abs(tb["string_read"] - r["end_reads"]["string_tension_read"])
        rd = radial_reads(e, cfg)
        rec["radial"] = rd
        rec["radial_max_rel_dev"] = max(
            rel(rd["E_lt_R"][R], r["end_reads"]["radial"]["E_lt_R"][R])
            for R in rd["E_lt_R"]
            if R in r["end_reads"]["radial"]["E_lt_R"]
        )
        rec["face_layer_row"] = r["end_reads"]["face_layer"]["E"]
        rec["face_layer_dev"] = abs(rd["E_pin_shell"] - r["end_reads"]["face_layer"]["E"])
        rec["face_layer_share"] = rd["E_pin_shell"] / E1["E_total"]
        rec["density_sum_rel_dev"] = rel(float(e.sum()), E1["E_total"])
        rec["norm"] = norm_reads(Mp, lam, cfg)
        rec["norm_row"] = {
            k: r["end_reads"]["norm3"][k]
            for k in ("axis_norm_min_rel", "axis_crossings_gap_lt_0.02", "axis_gap_min")
        }
        del e, lam
        if n == 32:
            rec["derrick"] = derrick_reads(Mp, cfg, q, w, E1, rd["r_half"], ec, ev)
            dr = r["end_reads"]["derrick"]
            rec["derrick_row"] = {
                "virial": dr["virial_E_curv_over_V"],
                "dE_dlam_order3": dr["order3"]["dE_dlambda_at_1"],
                "dE_dlam_order1": dr["order1"]["dE_dlambda_at_1"],
                "r_half": dr["r_half_energy"],
                "R_star": dr["R_star_from_virial"],
            }
            if (
                r["obj"] == "S1"
                and r["bnd"] == "seed"
                and w1s == 1.0
                and "frame3" in r["end_reads"]
            ):
                rec["frame"] = frame_reads(Mp, cfg, q, w, E1)
                f3 = r["end_reads"]["frame3"]
                rec["frame_row"] = {
                    "c": f3["d2E_ds2"],
                    "c0_extrapolated": f3["c0_extrapolated"],
                    "dV": f3["dV"],
                }
                frame_rows[tag] = (g, d, rec["frame"], rec["frame_row"], r["label"])
        del Mp, ec, ev
        gc.collect()
        A["rows"][tag] = rec
        log(
            f"{tag:26s} E row {r['E']:9.5f} mine {E1['E_total']:9.5f} rel {rec['rel_E']:.1e} | E0 rel {rec['rel_E0']:.1e} ({seed_kind}) | shell eq {rec['shell_bitwise_equal_seed_polished']} "
            f"dev {rec['shell_max_dev_from_seed_polished']:.1e} | fmax mine {rec.get('fmax_free_mine', float('nan')):.3e} row {r['polish']['fmax_end']:.3e} {r['label'][:8]} "
            f"| deg9 {rec['deg9_by_rank_mine']} (row winding {rec['deg9_winding_row']:+.3f}) | str {tb['string_read']:+.4f} flat {tb['flat_mid_tension']:+.4f} | face {rd['E_pin_shell']:.3f}"
        )
        for k in range(3):
            v = dg[f"rank{k}"]
            log(
                f"    rank{k} preimage {[v[hv]['preimage'] for hv in ('6', '9', '12')]} spread {[v[hv]['spread'] for hv in ('6', '9', '12')]} conflicts {[v[hv]['conflicts'] for hv in ('6', '9', '12')]} "
                f"solid {[round(v[hv]['solid_angle'], 3) for hv in ('6', '9', '12')]} row solid {[round(rec['degree_row_solid_angle'][f'rank{k}'].get(hv, float('nan')), 3) for hv in ('6', '9', '12')]} "
                f"bfs solid {[round(v[hv]['solid_angle_bfs_orientation'], 3) for hv in ('6', '9', '12')]} "
                f"gap min surf {[round(v[hv]['surface_gap_min'], 4) for hv in ('6', '9', '12')]} at conflicts {[None if v[hv]['gap_min_at_conflicts'] is None else round(v[hv]['gap_min_at_conflicts'], 4) for hv in ('6', '9', '12')]}"
            )
        if "derrick" in rec:
            dd = rec["derrick"]
            log(
                f"    derrick vir {dd['virial']:.3f} (row {rec['derrick_row']['virial']:.3f}) dE/dlam o2 {dd['order2']['dE_dlam_4pt']:+.3f} o1 {dd['order1']['dE_dlam_4pt']:+.3f} (row o3 {rec['derrick_row']['dE_dlam_order3']:+.3f}) "
                f"cont {dd['continuum_dE_dlam']:+.3f} r_half {dd['r_half']:.3f} (row {rec['derrick_row']['r_half']:.3f}) R* {dd['R_star']:.3f} (row {rec['derrick_row']['R_star']:.3f})"
            )
        if "frame" in rec:
            fr = rec["frame"]
            log(
                f"    frame c {[round(fr['c'][s], 4) for s in ('0.005', '0.01', '0.02', '0.05')]} c0 lsq {fr['c0_lsq_3pt']:.5f} rich {fr['c0_richardson_0.005_0.01']:.5f} (row c0 {rec['frame_row']['c0_extrapolated']:.5f}) dV rel max {max(fr['dV_rel'].values()):.1e}"
            )
    R = A["rows"]
    n32 = {t: v for t, v in R.items() if v["n"] == 32}

    # ---- C1 ----
    line(
        "C1a_E_total_polished_all_21_fields_rel_1e-10",
        len(R) == 21 and max(v["rel_E"] for v in R.values()) < 1e-10,
        f"worst rel {max(v['rel_E'] for v in R.values()):.2e} over {len(R)} fields",
    )
    line(
        "C1b_E_curv_and_V_parts_rel_1e-10",
        max(max(v["rel_curv"], v["rel_V"]) for v in R.values()) < 1e-10,
        f"worst curv {max(v['rel_curv'] for v in R.values()):.2e} V {max(v['rel_V'] for v in R.values()):.2e}",
    )
    ff = {t: v for t, v in R.items() if "E_fire_mine" in v}
    line(
        "C1c_E_fire_all_17_FIRE_fields_rel_1e-10",
        len(ff) == 17 and max(v["rel_E_fire"] for v in ff.values()) < 1e-10,
        f"{len(ff)} FIRE fields, worst rel {max(v['rel_E_fire'] for v in ff.values()):.2e}",
    )
    line(
        "C1d_seed_E0_every_row_rel_1e-10_my_seed_vs_descent_E0_or_polish_E_in",
        max(v["rel_E0"] for v in R.values()) < 1e-10,
        f"worst rel {max(v['rel_E0'] for v in R.values()):.2e}; kinds {sorted(set(v['seed_kind'] for v in R.values()))}",
    )
    line(
        "C1e_M0i_zero_every_polished_and_FIRE_field",
        all(v["max_abs_M0i"] == 0.0 and v.get("fire_max_abs_M0i", 0.0) == 0.0 for v in R.values()),
        f"max |M_0i| {max(v['max_abs_M0i'] for v in R.values()):.1e}",
    )
    line(
        "C1f_roots_and_W1_eff_equal_mine_pin_mask_equals_B3_free_cells_equal",
        all(
            tuple(v["roots_row"]) == tuple(v["roots_mine"])
            and abs(v["W1_eff_row"] - v["W1_eff_mine"]) < 1e-15
            and v["pin_mask_equals_B3"]
            and v["free_cells_mine"] == v["free_cells_row"]
            for v in R.values()
        ),
        f"free cells {sorted(set(v['free_cells_mine'] for v in R.values()))}",
    )
    line(
        "C1h_stack_energy_equals_mine_rel_1e-10_n32",
        max(v["rel_E_stack_vs_mine"] for v in n32.values()) < 1e-10,
        f"worst {max(v['rel_E_stack_vs_mine'] for v in n32.values()):.2e}",
    )
    line(
        "C1i_density_sums_to_E_rel_1e-12",
        max(v["density_sum_rel_dev"] for v in R.values()) < 1e-12,
        f"worst {max(v['density_sum_rel_dev'] for v in R.values()):.2e}",
    )
    m00 = {t: v["max_abs_M00_minus_g"] for t, v in R.items()}
    line(
        "C1j_M00_moved_off_g_in_the_interior_of_every_row_but_stays_g_exactly_on_every_B_seed_shell",
        all(v > 0 for v in m00.values())
        and all(
            v["max_abs_M00_minus_g_on_shell"] == 0.0 for v in R.values() if v["bnd"] == "seed"
        ),
        f"max |M_00 - g| interior {max(m00.values()):.3f} min {min(m00.values()):.1e}; on B_seed shells {max(v['max_abs_M00_minus_g_on_shell'] for v in R.values() if v['bnd'] == 'seed'):.1e}; "
        f"on B_far shells (relaxed n64 values) {max(v['max_abs_M00_minus_g_on_shell'] for v in R.values() if v['bnd'] == 'far'):.1e}",
    )

    # ---- C2 ----
    line(
        "C2a_max_abs_G_on_my_free_mask_equals_polish_fmax_end_rel_1e-6_all_20_n32_rows",
        len(n32) == 20 and max(v["rel_fmax"] for v in n32.values()) < 1e-6,
        f"worst rel {max(v['rel_fmax'] for v in n32.values()):.2e}",
    )
    line(
        "C2b_label_POLISHED_iff_max_abs_G_lt_1e-3_reproduced_all_n32_rows",
        all(v["label_agrees"] for v in n32.values()),
        f"POLISHED: {sorted(t for t, v in n32.items() if v['label_mine'] == 'POLISHED')}",
    )
    npol = sum(1 for v in R.values() if v["label_row"] == "POLISHED")
    line(
        "C2c_only_3_of_21_rows_reach_the_gate_g1_g2_before_the_polish_g4_by_the_callback",
        npol == 3
        and sorted(t for t, v in R.items() if v["label_row"] == "POLISHED")
        == ["S1_Bseed_g1_d0.3_w1_n32", "S1_Bseed_g2_d0.3_w1_n32", "S1_Bseed_g4_d0.3_w1_n32"],
        f"{npol} POLISHED; the rest FALLING at max |G| from {min(v['fmax_free_mine'] for v in n32.values() if v['label_mine'] != 'POLISHED'):.2e} to {max(v['fmax_free_mine'] for v in n32.values() if v['label_mine'] != 'POLISHED'):.2e}",
    )
    line(
        "C2d_FIRE_end_max_abs_G_on_my_free_mask_equals_descent_fmax_end_rel_1e-6",
        all(v["rel_fire_fmax"] < 1e-6 for v in n32.values() if "rel_fire_fmax" in v),
        f"worst rel {max(v['rel_fire_fmax'] for v in n32.values() if 'rel_fire_fmax' in v):.2e} over {sum(1 for v in n32.values() if 'rel_fire_fmax' in v)} rows",
    )
    pinned_rows = {t: v for t, v in R.items() if v["bnd"] != "free"}
    line(
        "C2e_pinned_shell_bitwise_equal_seed_on_polished_and_FIRE_fields_all_pinned_rows",
        all(
            v["shell_bitwise_equal_seed_polished"] and v.get("shell_bitwise_equal_seed_fire", True)
            for v in pinned_rows.values()
        ),
        f"{len(pinned_rows)} pinned rows, max shell dev {max(v['shell_max_dev_from_seed_polished'] for v in pinned_rows.values()):.1e}",
    )
    r20p = {
        t: v
        for t, v in R.items()
        if "shell_bitwise_equal_seed_axes_rebuilt" in v and v["bnd"] == "seed"
    }
    line(
        "C2f_polish_only_R20_rows_shell_equals_seed_axes_rebuilt_bitwise_4_rows",
        len(r20p) == 4 and all(v["shell_bitwise_equal_seed_axes_rebuilt"] for v in r20p.values()),
        f"{sorted(r20p)}",
    )
    far = {t: v for t, v in R.items() if v["bnd"] == "far"}
    line(
        "C2g_B_far_shell_equals_my_restriction_of_the_R20_n64_field_bitwise_3_rows",
        len(far) == 3 and all(v["shell_bitwise_equal_seed_polished"] for v in far.values()),
        f"{sorted(far)}",
    )
    frees = {t: v for t, v in R.items() if v["bnd"] == "free"}
    line(
        "C2h_B_free_shells_moved_off_the_seed_every_free_row",
        all(v["shell_max_dev_from_seed_polished"] > 1e-3 for v in frees.values()),
        f"shell max dev {{{', '.join(f'{t}: {v['shell_max_dev_from_seed_polished']:.3f}' for t, v in frees.items())}}}",
    )
    raised = sorted(t for t, v in n32.items() if v.get("polish_raised_fmax"))
    line(
        "C2i_polish_raised_max_abs_G_on_no_row",
        len(raised) == 0,
        f"{len(raised)} rows raised: "
        + "; ".join(
            f"{t} {n32[t]['polish_fmax_in']:.2e} -> {n32[t]['fmax_free_mine']:.2e}" for t in raised
        ),
    )
    spread = {}
    for obj in ("S1", "Sd", "S0"):
        fm = {
            b: n32[f"{obj}_B{b}_g8_d0.3_w1_n32"]["fmax_free_mine"] for b in ("seed", "far", "free")
        }
        spread[obj] = {"fmax": fm, "max_over_min": max(fm.values()) / min(fm.values())}
    A["C2_residual_spread"] = spread
    line(
        "C2j_boundary_triple_at_equal_residual_max_over_min_fmax_lt_2_every_axis",
        all(v["max_over_min"] < 2.0 for v in spread.values()),
        "; ".join(
            f"{o}: "
            + ", ".join(f"{b} {x:.2e}" for b, x in v["fmax"].items())
            + f" (x{v['max_over_min']:.1f})"
            for o, v in spread.items()
        ),
    )
    line(
        "C2k_re_collected_residual_spread_max_over_min_equals_mine_rel_1e-8_and_the_per_axis_wording_says_not_one_residual",
        all(
            rel(
                spread[o]["max_over_min"],
                res["per_axis"][o].get("residual_spread_max_over_min", float("nan")),
            )
            < 1e-8
            and "not one residual" in res["per_axis"][o]["outcomes"][0]
            for o in spread
        ),
        f"row spreads {[round(res['per_axis'][o].get('residual_spread_max_over_min', float('nan')), 3) for o in spread]}",
    )

    # ---- C3 ----
    per = {}
    for obj in ("S1", "Sd", "S0"):
        tri = {b: n32[f"{obj}_B{b}_g8_d0.3_w1_n32"] for b in ("seed", "far", "free")}
        Es = {b: v["E_mine"] for b, v in tri.items()}
        e12 = {b: v["radial"]["E_lt_R"]["12"] for b, v in tri.items()}
        dg9 = {b: v["deg9_winding_mine"] for b, v in tri.items()}
        rel12 = (max(e12.values()) - min(e12.values())) / abs(np.mean(list(e12.values())))
        reldiff = max(abs(Es["far"] - Es["seed"]), abs(Es["free"] - Es["seed"])) / abs(Es["seed"])
        wind_all = all(abs(abs(v) - 1.0) < 0.1 for v in dg9.values())
        wind_lost = any(abs(v) < 0.5 for v in dg9.values())
        outc = []
        if rel12 <= 0.02 and wind_all:
            outc.append("INTERIOR_SET")
        if reldiff > 0.05 or wind_lost:
            outc.append("BOUNDARY_SET")
        if not outc:
            outc.append("BOUNDARY_UNDECIDED")
        per[obj] = {
            "E": Es,
            "E_lt_12": e12,
            "deg9": dg9,
            "rel_spread_E_lt_12": float(rel12),
            "rel_boundary_diff_total": float(reldiff),
            "max_boundary_abs_diff": float(
                max(
                    abs(Es["far"] - Es["seed"]),
                    abs(Es["free"] - Es["seed"]),
                    abs(Es["free"] - Es["far"]),
                )
            ),
            "outcome_mine": outc,
            "outcome_row": res["per_axis"][obj]["outcomes"],
            "E_lt_12_row": res["per_axis"][obj]["E_lt_12"],
            "deg9_row": res["per_axis"][obj]["degree_r9"],
        }
    A["C3_per_axis"] = per
    line(
        "C3a_per_axis_labels_BOUNDARY_SET_reproduced_on_my_numbers_all_three",
        all(
            per[o]["outcome_mine"] == ["BOUNDARY_SET"]
            and res["per_axis"][o]["outcomes"][0].startswith("BOUNDARY_SET")
            for o in per
        ),
        "; ".join(
            f"{o}: E(<12) spread {v['rel_spread_E_lt_12']:.3f}, total diff {v['rel_boundary_diff_total']:.3f}, deg9 {[round(x, 2) for x in v['deg9'].values()]}"
            for o, v in per.items()
        ),
    )
    line(
        "C3b_per_axis_numbers_equal_rows_E_lt_12_rel_1e-8_spread_and_diff_abs_1e-8",
        all(
            max(
                rel(per[o]["E_lt_12"][b], per[o]["E_lt_12_row"][b])
                for b in ("seed", "far", "free")
            )
            < 1e-8
            and abs(per[o]["rel_spread_E_lt_12"] - res["per_axis"][o]["rel_spread_E_lt_12"]) < 1e-8
            and abs(
                per[o]["rel_boundary_diff_total"] - res["per_axis"][o]["rel_boundary_diff_total"]
            )
            < 1e-8
            for o in per
        ),
        "",
    )
    line(
        "C3c_my_preimage_degree_equals_row_solid_angle_abs_0.02_on_every_orientable_cube",
        all(
            v["degree_abs_dev_vs_row_orientable_cubes"] is not None
            and v["degree_abs_dev_vs_row_orientable_cubes"] < 0.02
            for v in R.values()
        ),
        f"worst {max(v['degree_abs_dev_vs_row_orientable_cubes'] for v in R.values() if v['degree_abs_dev_vs_row_orientable_cubes'] is not None):.1e}; orientable cubes per row {[v['degree_orientable_cubes'] for v in R.values()]}",
    )
    w9 = {t: R[t]["degree"][f"rank{R[t]['winding_rank']}"]["9"] for t in R}
    reads1 = [t for t in R if abs(abs(R[t]["deg9_winding_row"]) - 1) < 0.1]
    torn1 = [t for t in reads1 if w9[t]["conflicts"] > 0]
    line(
        "C3d_winding_rank_r9_surface_orientable_no_conflicts_on_every_row_the_producer_reads_degree_1",
        len(torn1) == 0,
        f"{len(reads1)} rows read 1 by the producer; {len(torn1)} of them have conflicts on my r 9 surface: "
        + "; ".join(
            f"{t[:22]} {w9[t]['conflicts']} links, gap at conflicts {w9[t]['gap_min_at_conflicts']:.3f}"
            for t in torn1
        ),
    )
    orient_ok = all(
        not v["degree"][f"rank{k}"][hv]["orientation_dependent"]
        for v in R.values()
        for k in range(3)
        for hv in ("6", "9", "12")
        if v["degree"][f"rank{k}"][hv]["conflicts"] == 0
    )
    dep = [
        (t, k, hv)
        for t, v in R.items()
        for k in range(3)
        for hv in ("6", "9", "12")
        if v["degree"][f"rank{k}"][hv]["conflicts"] > 0
        and v["degree"][f"rank{k}"][hv]["orientation_dependent"]
    ]
    conf_cubes = sum(
        1
        for v in R.values()
        for k in range(3)
        for hv in ("6", "9", "12")
        if v["degree"][f"rank{k}"][hv]["conflicts"] > 0
    )
    line(
        "C3d2_solid_angle_integer_never_changes_with_the_orientation_on_conflict_free_cubes_and_does_change_on_some_conflicted_cube",
        orient_ok and len(dep) > 0,
        f"orientation-dependent integers on {len(dep)} of {conf_cubes} conflicted cubes (e.g. {dep[:4]}), on 0 conflict-free cubes: {orient_ok}",
    )
    wdep = [t for t in reads1 if w9[t]["orientation_dependent"]]
    line(
        "C3d3_winding_rank_r9_integer_is_orientation_INdependent_on_every_row_the_producer_reads_1",
        len(wdep) == 0,
        f"orientation-dependent on {len(wdep)} rows: "
        + ", ".join(
            f"{t[:22]} (greedy {w9[t]['solid_angle']:+.0f}, bfs {w9[t]['solid_angle_bfs_orientation']:+.0f}, producer {R[t]['deg9_winding_row']:+.0f})"
            for t in wdep
        ),
    )
    s1c = {b: w9[f"S1_B{b}_g8_d0.3_w1_n32"]["conflicts"] for b in ("seed", "far", "free")}
    line(
        "C3e_S1_winding_read_1_by_preimage_under_all_three_boundaries_Sd_and_S0_winding_rank_torn_on_r9_under_all_three",
        all(abs(abs(per["S1"]["deg9"][b]) - 1) < 0.1 for b in ("seed", "far", "free"))
        and all(
            w9[f"{o}_B{b}_g8_d0.3_w1_n32"]["conflicts"] > 0
            for o in ("Sd", "S0")
            for b in ("seed", "far", "free")
        ),
        f"S1 preimage {per['S1']['deg9']} (conflicts {s1c}); Sd conflicts {[w9[f'Sd_B{b}_g8_d0.3_w1_n32']['conflicts'] for b in ('seed', 'far', 'free')]} gap {[round(w9[f'Sd_B{b}_g8_d0.3_w1_n32']['gap_min_at_conflicts'], 3) for b in ('seed', 'far', 'free')]}; "
        f"S0 conflicts {[w9[f'S0_B{b}_g8_d0.3_w1_n32']['conflicts'] for b in ('seed', 'far', 'free')]}",
    )
    sr = {t: v for t, v in R.items() if "surface_reader_vs_mine" in v}
    line(
        "C3g_re_collected_surface_reader_conflict_counts_equal_mine_exactly_every_rank_and_cube_all_21_rows",
        len(sr) == 21 and all(v["surface_conflicts_all_equal"] for v in sr.values()),
        f"{len(sr)} rows with the degree_surface read, {sum(v['surface_reader_cubes'] for v in sr.values())} rank-cubes; mismatches on {[t for t, v in sr.items() if not v['surface_conflicts_all_equal']]}",
    )
    mism = [
        (t, c)
        for t, v in sr.items()
        for c in v["surface_reader_vs_mine"]
        if c["reading_mine"] != c["reading_row"]
    ]
    line(
        "C3g2_re_collected_kept_lost_undefined_readings_equal_mine_every_rank_and_cube",
        len(mism) == 0,
        f"{len(mism)} reading mismatches: {[(t[:20], c['rank'], c['cube'], c['reading_mine'], c['reading_row']) for t, c in mism[:6]]}",
    )
    s0f = n32["S0_Bfree_g8_d0.3_w1_n32"]["degree"]["rank1"]["9"]
    line(
        "C3f_S0_B_free_rank1_carrier_at_r9_established_conflict_free_with_an_orientation_independent_integer_1",
        s0f["conflicts"] == 0
        and not s0f["orientation_dependent"]
        and abs(abs(s0f["preimage"]) - 1) < 0.1,
        f"S0 free rank1 r9: preimage {s0f['preimage']}, conflicts {s0f['conflicts']}, gap at conflicts {s0f['gap_min_at_conflicts']}, greedy {s0f['solid_angle']:+.0f} bfs {s0f['solid_angle_bfs_orientation']:+.0f} producer +1",
    )

    # ---- C4 ----
    sdf, sds = n32["Sd_Bfree_g8_d0.3_w1_n32"], n32["Sd_Bseed_g8_d0.3_w1_n32"]
    t_free, t_seed = sdf["tube"]["string_read"], sds["tube"]["string_read"]
    flat_free, flat_seed = sdf["tube"]["flat_mid_tension"], sds["tube"]["flat_mid_tension"]
    deg = sdf["deg9_winding_mine"]

    def string_label(tf, ts, dg):
        lost = abs(dg) < 0.5
        if tf < 0.02 or lost:
            return "STRING_ESCAPES"
        if abs(tf - ts) <= 0.3 * abs(ts) and abs(abs(dg) - 1) < 0.1:
            return "STRING_HELD"
        return "STRING_UNDECIDED"

    st = {
        "tension_free": t_free,
        "tension_seed": t_seed,
        "flat_mid_free": flat_free,
        "flat_mid_seed": flat_seed,
        "deg9_free": deg,
        "face_layer_free": sdf["radial"]["E_pin_shell"],
        "label_on_tube_read": string_label(t_free, t_seed, deg),
        "label_on_flat_middle": string_label(flat_free, flat_seed, deg),
        "label_row": res["string"]["outcomes"],
        "per_plane_excess_free": sdf["tube"]["per_plane_excess"],
        "per_plane_excess_seed": sds["tube"]["per_plane_excess"],
        "planes": sdf["tube"]["planes_abs_z"],
    }
    A["C4_string"] = st
    line(
        "C4a_tube_read_equals_rows_abs_1e-8_all_21_rows",
        max(v["string_read_dev"] for v in R.values()) < 1e-8,
        f"worst {max(v['string_read_dev'] for v in R.values()):.1e}",
    )
    line(
        "C4b_Sd_B_free_tube_read_below_0.02_STRING_ESCAPES_reproduced_by_the_tension_clause",
        t_free < 0.02
        and st["label_on_tube_read"] == "STRING_ESCAPES"
        and len(res["string"]["outcomes"]) == 1
        and res["string"]["outcomes"][0].startswith("STRING_ESCAPES"),
        f"tube read free {t_free:+.5f} seed {t_seed:+.4f}; the re-collected label: {res['string']['outcomes'][0]}",
    )
    sdw = sdf["degree"]["rank1"]["9"]
    line(
        "C4b2_Sd_B_free_degree_KEPT_clause_established_rank1_r9_conflict_free_orientation_independent_1",
        sdw["conflicts"] == 0
        and not sdw["orientation_dependent"]
        and abs(abs(sdw["preimage"]) - 1) < 0.1,
        f"rank1 r9: preimage {sdw['preimage']}, conflicts {sdw['conflicts']}, gap at conflicts {sdw['gap_min_at_conflicts']:.4f}, greedy {sdw['solid_angle']:+.0f} bfs {sdw['solid_angle_bfs_orientation']:+.0f} producer -1",
    )
    line(
        "C4c_Sd_B_free_flat_middle_tension_also_below_0.02_label_unchanged_on_the_flat_middle_read",
        flat_free < 0.02 and st["label_on_flat_middle"] == "STRING_ESCAPES",
        f"flat mid free {flat_free:+.5f} (near quarter {sdf['tube']['near_quarter_tension']:+.5f}, far quarter {sdf['tube']['far_quarter_tension']:+.5f}); seed flat {flat_seed:+.4f}",
    )
    line(
        "C4d_Sd_B_free_face_layer_below_0.3_and_below_a_fifth_of_the_total",
        sdf["radial"]["E_pin_shell"] < 0.3 and sdf["face_layer_share"] < 0.2,
        f"face layer {sdf['radial']['E_pin_shell']:.4f} share {sdf['face_layer_share']:.3f}; seed face layer {sds['radial']['E_pin_shell']:.3f} share {sds['face_layer_share']:.3f}",
    )
    routes = {}
    for t, v in n32.items():
        if v["w1s"] != 1.0 or v["g"] != 8.0 or v["delta"] != 0.3:
            continue
        nr = v["norm"]
        route = (
            "ESCAPE_BY_EXCHANGE"
            if (nr["axis_norm_min_rel"] >= 0.9 and nr["axis_crossings_gap_lt_0.02"] > 0)
            else (
                "ESCAPE_BY_MELTING" if nr["axis_norm_min_rel"] < 0.7 else "ESCAPE_ROUTE_UNDECIDED"
            )
        )
        carrier = [
            k
            for k, x in v["deg9_by_rank_mine"].items()
            if abs(abs(x) - 1) < 0.1
            and k != f"rank{v['winding_rank']}"
            and v["degree"][k]["9"]["conflicts"] == 0
        ]
        routes[f"{v['obj']}_B{v['bnd']}"] = {
            "winding_lost_by_producer": f"{v['obj']}_B{v['bnd']}" in res["escape_routes"],
            "winding_r9_torn_mine": v["degree"][f"rank{v['winding_rank']}"]["9"]["conflicts"] > 0,
            "preimage_winding_r9": v["deg9_winding_mine"],
            "route_by_norm_rule_mine": route,
            "route_row": res["escape_routes"].get(f"{v['obj']}_B{v['bnd']}", {}).get("route"),
            "axis_norm_min_rel": nr["axis_norm_min_rel"],
            "crossings": nr["axis_crossings_gap_lt_0.02"],
            "conflict_free_carrier_mine": carrier,
            "carrier_row": res["escape_routes"]
            .get(f"{v['obj']}_B{v['bnd']}", {})
            .get("winding_transferred_to"),
        }
    A["C4_routes"] = routes
    for k, v in routes.items():
        v["route_row"] = (
            res["escape_routes"]
            .get(k, {})
            .get("route_by_norm_rule", res["escape_routes"].get(k, {}).get("route"))
        )
    flagged = {k: v for k, v in routes.items() if v["winding_lost_by_producer"]}
    mine_not_kept = sorted(
        k
        for k, v in routes.items()
        if v["winding_r9_torn_mine"] or abs(v["preimage_winding_r9"]) < 0.9
    )
    line(
        "C4e_re_collected_escape_route_rows_equal_the_rows_where_my_r9_winding_read_is_not_kept_and_the_norm_rule_routes_agree",
        sorted(flagged) == mine_not_kept
        and all(v["route_by_norm_rule_mine"] == v["route_row"] for v in flagged.values()),
        f"rows {sorted(flagged)}; "
        + "; ".join(
            f"{k}: {v['route_by_norm_rule_mine']} (row {v['route_row']}) norm {v['axis_norm_min_rel']:.3f} crossings {v['crossings']}"
            for k, v in flagged.items()
        ),
    )
    line(
        "C4e2_no_conflict_free_rank_carries_a_transferred_winding_on_any_flagged_row_mine_and_the_re_collected_reader_agree",
        all(
            v["conflict_free_carrier_mine"] == [] and v["carrier_row"] == []
            for v in flagged.values()
        ),
        f"carriers mine {[v['conflict_free_carrier_mine'] for v in flagged.values()]} row {[v['carrier_row'] for v in flagged.values()]}",
    )
    line(
        "C4f_norm_reads_equal_rows_axis_norm_min_rel_1e-8_crossings_exact",
        all(
            rel(v["norm"]["axis_norm_min_rel"], v["norm_row"]["axis_norm_min_rel"]) < 1e-8
            and v["norm"]["axis_crossings_gap_lt_0.02"]
            == v["norm_row"]["axis_crossings_gap_lt_0.02"]
            for v in R.values()
        ),
        "",
    )
    kept = [t for t, v in n32.items() if abs(abs(v["deg9_winding_mine"]) - 1) < 0.1]
    lostr = [t for t, v in n32.items() if abs(v["deg9_winding_mine"]) < 0.5]
    line(
        "C4g_axis_norm_alone_does_not_discriminate_norm_min_rel_within_10pct_of_vacuum_on_a_kept_row_AND_on_a_lost_row",
        any(abs(n32[t]["norm"]["axis_norm_min_rel"] - 1) <= 0.1 for t in kept)
        and any(abs(n32[t]["norm"]["axis_norm_min_rel"] - 1) <= 0.1 for t in lostr),
        f"kept rows norm min rel {[round(n32[t]['norm']['axis_norm_min_rel'], 3) for t in kept]}; lost rows {[round(n32[t]['norm']['axis_norm_min_rel'], 3) for t in lostr]}",
    )

    # ---- C5 ----
    box = {}
    for obj in ("S1", "Sd", "S0"):
        v = n32[f"{obj}_Bfree_g8_d0.3_w1_n32"]
        dd = v["derrick"]
        de = dd["order2"]["dE_dlam_4pt"]
        lab = (
            "AXIS_BOX_LIMITED"
            if (de < 0 and dd["R_star"] > 12.0)
            else (
                "AXIS_LOCALIZED"
                if (dd["R_star"] < 8.0 and abs(dd["virial"] - 3.0) <= 0.9)
                else "BOX_UNDECIDED"
            )
        )
        box[obj] = {
            "dE_dlam_mine_o2": de,
            "dE_dlam_mine_o1": dd["order1"]["dE_dlam_4pt"],
            "dE_dlam_row_o3": v["derrick_row"]["dE_dlam_order3"],
            "continuum": dd["continuum_dE_dlam"],
            "R_star_mine": dd["R_star"],
            "R_star_row": v["derrick_row"]["R_star"],
            "virial_mine": dd["virial"],
            "virial_row": v["derrick_row"]["virial"],
            "label_mine": lab,
            "label_row": res["box"][obj]["outcome"],
        }
    A["C5_box"] = box
    dk = {t: v for t, v in n32.items() if "derrick" in v}
    line(
        "C5a_dE_dlam_order2_within_15pct_of_row_order3_and_same_sign_all_n32_rows",
        all(
            rel(v["derrick"]["order2"]["dE_dlam_4pt"], v["derrick_row"]["dE_dlam_order3"]) < 0.15
            and np.sign(v["derrick"]["order2"]["dE_dlam_4pt"])
            == np.sign(v["derrick_row"]["dE_dlam_order3"])
            for v in dk.values()
        ),
        f"worst rel {max(rel(v['derrick']['order2']['dE_dlam_4pt'], v['derrick_row']['dE_dlam_order3']) for v in dk.values()):.3f}",
    )
    line(
        "C5b_dE_dlam_negative_every_n32_row_both_orders_both_stencils",
        all(
            v["derrick"][o][s] < 0
            for v in dk.values()
            for o in ("order2", "order1")
            for s in ("dE_dlam_4pt", "dE_dlam_central_0.05")
        ),
        "",
    )
    line(
        "C5c_virial_equals_rows_rel_1e-8_r_half_and_R_star_rel_1e-2",
        all(
            rel(v["derrick"]["virial"], v["derrick_row"]["virial"]) < 1e-8
            and rel(v["derrick"]["r_half"], v["derrick_row"]["r_half"]) < 1e-2
            and rel(v["derrick"]["R_star"], v["derrick_row"]["R_star"]) < 1e-2
            for v in dk.values()
        ),
        f"worst r_half rel {max(rel(v['derrick']['r_half'], v['derrick_row']['r_half']) for v in dk.values()):.1e}",
    )
    for o in box:
        box[o]["identity_row"] = res["box"][o].get("dE_dlambda_identity")
    line(
        "C5d_pre_registered_box_rule_gives_AXIS_BOX_LIMITED_on_my_numbers_all_three_axes_the_re_collected_outcome_is_BOX_UNDECIDED_with_the_identity_equal_to_my_continuum_rel_1e-8",
        all(
            box[o]["label_mine"] == "AXIS_BOX_LIMITED"
            and box[o]["label_row"].startswith("BOX_UNDECIDED")
            and box[o]["identity_row"] is not None
            and rel(box[o]["identity_row"], box[o]["continuum"]) < 1e-8
            for o in box
        ),
        "; ".join(
            f"{o}: dE/dlam {v['dE_dlam_mine_o2']:+.3f} (row scan {v['dE_dlam_row_o3']:+.3f}, row identity {v['identity_row']}, my continuum {v['continuum']:+.3f}) R* {v['R_star_mine']:.2f} vir {v['virial_mine']:.2f}"
            for o, v in box.items()
        ),
    )
    line(
        "C5e_B_free_dilation_slope_has_the_sign_of_the_continuum_-Ecurv+3V_every_axis",
        all(box[o]["dE_dlam_mine_o2"] * box[o]["continuum"] > 0 for o in box),
        f"dilation / continuum ratios {[round(box[o]['dE_dlam_mine_o2'] / box[o]['continuum'], 2) for o in box]} (continuum {[round(box[o]['continuum'], 3) for o in box]}: the Derrick balance says CONTRACT, the dilation says EXPAND)",
    )
    line(
        "C5f_B_free_virial_below_3_every_axis_the_potential_dominates_the_scaling_force",
        all(box[o]["virial_mine"] < 3.0 for o in box),
        f"{[round(box[o]['virial_mine'], 2) for o in box]}",
    )
    split = {}
    for o in box:
        d2 = n32[f"{o}_Bfree_g8_d0.3_w1_n32"]["derrick"]["order2"]
        split[o] = {
            "total": d2["dE_dlam_central_0.05"],
            "face_layer": d2["dE_dlam_face_layer"],
            "interior": d2["dE_dlam_interior"],
            "r_lt_12": d2["dE_dlam_r_lt_12"],
            "r_ge_12": d2["dE_dlam_r_ge_12"],
            "continuum_r_lt_12": n32[f"{o}_Bfree_g8_d0.3_w1_n32"]["derrick"][
                "continuum_dE_dlam_r_lt_12"
            ],
        }
        box[o]["split"] = split[o]
    line(
        "C5g_B_free_dilation_slope_carried_by_the_interior_not_the_face_layer_interior_part_negative_and_larger_than_the_face_part_every_axis",
        all(
            abs(split[o]["interior"]) > abs(split[o]["face_layer"]) and split[o]["interior"] < 0
            for o in box
        ),
        "; ".join(
            f"{o}: total {v['total']:+.3f} = face {v['face_layer']:+.3f} + interior {v['interior']:+.3f}; r<12 {v['r_lt_12']:+.3f} (continuum r<12 {v['continuum_r_lt_12']:+.3f}), r>=12 {v['r_ge_12']:+.3f}"
            for o, v in split.items()
        ),
    )
    line(
        "C5h_B_free_core_r_lt_12_dilation_slope_negative_every_axis_the_object_itself_wants_to_expand",
        all(split[o]["r_lt_12"] < 0 for o in box),
        f"r<12 slopes {[round(split[o]['r_lt_12'], 3) for o in box]} vs continuum r<12 {[round(split[o]['continuum_r_lt_12'], 3) for o in box]}",
    )

    # ---- C6 ----
    Ef = {o: per[o]["E"]["free"] for o in ("S1", "Sd", "S0")}
    diffs = [abs(Ef["S1"] - Ef["Sd"]), abs(Ef["S1"] - Ef["S0"]), abs(Ef["Sd"] - Ef["S0"])]
    resol = max(per[o]["max_boundary_abs_diff"] for o in per)
    arr = np.array(list(Ef.values()))
    srt = np.sort(arr)
    trip = {
        "E_free": Ef,
        "pairwise_abs_diffs": diffs,
        "resolution_3x": 3.0 * resol,
        "outcome_mine": "THREE_AXES_DISTINCT" if min(diffs) > 3.0 * resol else "AXES_DEGENERATE",
        "outcome_row": res["triple"]["outcome"],
        "koide_Q": float(np.sum(arr) / np.sum(np.sqrt(arr)) ** 2),
        "ratios_sorted": [float(x / srt[0]) for x in srt],
    }
    A["C6_triple"] = trip
    line(
        "C6a_triple_AXES_DEGENERATE_reproduced_resolution_and_koide_rel_1e-8",
        trip["outcome_mine"] == trip["outcome_row"] == "AXES_DEGENERATE"
        and rel(trip["resolution_3x"], res["triple"]["resolution_3x_max_boundary_diff"]) < 1e-8
        and rel(trip["koide_Q"], res["triple"]["koide_Q"]) < 1e-8,
        f"diffs {[round(x, 3) for x in diffs]} vs 3x resolution {3 * resol:.1f}; Koide {trip['koide_Q']:.4f}; ratios {[round(x, 3) for x in trip['ratios_sorted']]}",
    )
    line(
        "C6b_resolution_is_the_Sd_seed_minus_free_gap_and_exceeds_every_pairwise_difference_50_fold",
        resol == per["Sd"]["max_boundary_abs_diff"] and 3 * resol > 50 * max(diffs),
        f"resolution {resol:.2f} = Sd seed {per['Sd']['E']['seed']:.2f} minus free {per['Sd']['E']['free']:.2f}; largest pairwise {max(diffs):.2f}",
    )

    # ---- C7 ----
    lad = {}
    for obj, wv, b in (
        ("S1", 1.0, "seed"),
        ("S1", 5.0, "seed"),
        ("S1", 25.0, "seed"),
        ("S1", 25.0, "free"),
        ("Sd", 1.0, "seed"),
        ("Sd", 25.0, "seed"),
    ):
        v = n32[f"{obj}_B{b}_g8_d0.3_w{wv:g}_n32"]
        dd = v["derrick"]
        lad[f"{obj}_w{wv:g}_B{b}"] = {
            "E": v["E_mine"],
            "label": v["label_row"],
            "fmax": v["fmax_free_mine"],
            "R_star": dd["R_star"],
            "virial": dd["virial"],
            "dE_dlam": dd["order2"]["dE_dlam_4pt"],
            "shell_share": v["radial"]["E_pin_shell"] / v["E_mine"],
            "E_lt_12": v["radial"]["E_lt_R"]["12"],
            "deg9": v["deg9_winding_mine"],
            "r_half": dd["r_half"],
        }
    q25 = lad["S1_w25_Bseed"]
    base = lad["S1_w1_Bseed"]["dE_dlam"]
    if (
        q25["R_star"] < 8.0
        and abs(q25["virial"] - 3.0) <= 0.9
        and abs(q25["dE_dlam"]) <= 0.1 * abs(base)
    ):
        lab = "SCALE_IS_THE_POTENTIALS"
    elif q25["R_star"] > 12.0 or q25["virial"] > 6.0:
        lab = "SCALE_IS_THE_BOXS"
    else:
        lab = "SCALE_UNDECIDED"
    lad["outcome_mine"] = lab
    lad["outcome_row"] = res["w1_ladder"]["outcome"]
    lad["dE_dlam_reference_used_by_the_producer"] = res["w1_ladder"]["R20_S1_dE_dlambda_reference"]
    lad["dE_dlam_R20_S1_x4500_row"] = R20J["S1_v4std_n32_g8_x4500"]["end_reads"]["derrick"][
        "order3"
    ]["dE_dlambda_at_1"]
    A["C7_ladder"] = lad
    line(
        "C7a_ladder_label_SCALE_UNDECIDED_reproduced_on_my_numbers",
        lab == "SCALE_UNDECIDED" and res["w1_ladder"]["outcome"] == "SCALE_UNDECIDED",
        f"S1 w25 seed: R* {q25['R_star']:.2f} (gate < 8 / > 12), virial {q25['virial']:.2f}, dE/dlam {q25['dE_dlam']:+.2f} vs 0.1 x {abs(base):.2f}",
    )
    rs = [lad[k]["R_star"] for k in ("S1_w1_Bseed", "S1_w5_Bseed", "S1_w25_Bseed")]
    ss = [lad[k]["shell_share"] for k in ("S1_w1_Bseed", "S1_w5_Bseed", "S1_w25_Bseed")]
    line(
        "C7b_S1_R_star_and_shell_share_fall_monotonically_up_the_ladder_w1_w5_w25",
        rs[0] > rs[1] > rs[2] and ss[0] > ss[1] > ss[2],
        f"R* {[round(x, 2) for x in rs]} shell share {[round(x, 3) for x in ss]} r_half {[round(lad[k]['r_half'], 2) for k in ('S1_w1_Bseed', 'S1_w5_Bseed', 'S1_w25_Bseed')]}",
    )
    line(
        "C7c_S1_w25_dE_dlam_does_NOT_shrink_with_W1_it_grows_in_magnitude",
        abs(lad["S1_w25_Bseed"]["dE_dlam"]) > abs(lad["S1_w1_Bseed"]["dE_dlam"]),
        f"dE/dlam w1 {lad['S1_w1_Bseed']['dE_dlam']:+.2f} w5 {lad['S1_w5_Bseed']['dE_dlam']:+.2f} w25 {lad['S1_w25_Bseed']['dE_dlam']:+.2f}",
    )
    s25, f25 = lad["S1_w25_Bseed"], lad["S1_w25_Bfree"]
    line(
        "C7d_S1_w25_seed_vs_free_energies_within_3pct_R_star_within_10pct_both_winding_kept",
        rel(f25["E"], s25["E"]) < 0.03
        and rel(f25["R_star"], s25["R_star"]) < 0.1
        and abs(abs(s25["deg9"]) - 1) < 0.1
        and abs(abs(f25["deg9"]) - 1) < 0.1,
        f"E seed {s25['E']:.3f} free {f25['E']:.3f}; R* {s25['R_star']:.2f} / {f25['R_star']:.2f}; virial {s25['virial']:.2f} / {f25['virial']:.2f}; fmax {s25['fmax']:.2e} / {f25['fmax']:.2e}",
    )
    line(
        "C7e_Sd_w25_seed_shell_share_above_0.7_the_frozen_string_on_the_shell_dominates_both_Sd_seed_rows",
        lad["Sd_w1_Bseed"]["shell_share"] > 0.7 and lad["Sd_w25_Bseed"]["shell_share"] > 0.7,
        f"Sd shell share w1 {lad['Sd_w1_Bseed']['shell_share']:.3f} w25 {lad['Sd_w25_Bseed']['shell_share']:.3f}; E {lad['Sd_w1_Bseed']['E']:.2f} / {lad['Sd_w25_Bseed']['E']:.2f}",
    )

    # ---- C8 ----
    fr = {}
    for tag, (g, d, mine, rowf, lab) in frame_rows.items():
        fr[tag] = {
            "g": g,
            "delta": d,
            "label": lab,
            "c_mine": mine["c"],
            "c_row": rowf["c"],
            "c0_lsq_mine": mine["c0_lsq_3pt"],
            "c0_rich_mine": mine["c0_richardson_0.005_0.01"],
            "c0_row": rowf["c0_extrapolated"],
            "dV_rel_max": max(mine["dV_rel"].values()),
            "rel_c_vs_row": max(rel(mine["c"][s], rowf["c"][s]) for s in ("0.01", "0.02", "0.05")),
            "rel_c0_vs_row": rel(mine["c0_lsq_3pt"], rowf["c0_extrapolated"]),
            "rel_c0_lsq_vs_rich": rel(mine["c0_lsq_3pt"], mine["c0_richardson_0.005_0.01"]),
        }
    A["C8_frame_rows"] = fr
    line(
        "C8a_frame_curvatures_equal_rows_at_s_0.01_0.02_0.05_rel_1e-6_all_8_rows",
        len(fr) == 8 and max(v["rel_c_vs_row"] for v in fr.values()) < 1e-6,
        f"worst {max(v['rel_c_vs_row'] for v in fr.values()):.1e}",
    )
    line(
        "C8b_both_potentials_silent_along_the_boost_dV_rel_below_1e-9_all_s_all_rows",
        max(v["dV_rel_max"] for v in fr.values()) < 1e-9,
        f"worst dV rel {max(v['dV_rel_max'] for v in fr.values()):.1e}",
    )
    line(
        "C8c_c0_extrapolation_equals_rows_rel_1e-4_and_my_Richardson_pair_rel_1e-4",
        max(v["rel_c0_vs_row"] for v in fr.values()) < 1e-4
        and max(v["rel_c0_lsq_vs_rich"] for v in fr.values()) < 1e-4,
        f"worst vs row {max(v['rel_c0_vs_row'] for v in fr.values()):.1e}, lsq vs Richardson {max(v['rel_c0_lsq_vs_rich'] for v in fr.values()):.1e}",
    )
    line(
        "C8d_saddle_c0_negative_every_row",
        all(v["c0_lsq_mine"] < 0 for v in fr.values()),
        f"c0 {[round(v['c0_lsq_mine'], 2) for v in fr.values()]}",
    )
    gp = {v["g"]: v for v in fr.values() if v["delta"] == 0.3}
    dp = {v["delta"]: v for v in fr.values() if v["g"] == 8.0}
    gs = sorted(gp)
    ds = sorted(dp)
    fg = fit_root_law(gs, [gp[g]["c0_lsq_mine"] for g in gs])
    fd = fit_root_law(ds, [dp[d]["c0_lsq_mine"] for d in ds])
    rows8 = [(g, 0.3, gp[g]["c0_lsq_mine"]) for g in gs] + [
        (8.0, d, dp[d]["c0_lsq_mine"]) for d in ds if d != 0.3
    ]
    pf = product_fit(rows8, fg["x0"], fd["x0"])
    lab_mine = frame_label(fg["x0"], fd["x0"], pf["max_rel_miss"])
    pf_prod = product_fit(rows8, res["frame"]["fit_g"]["g0"], res["frame"]["fit_delta"]["delta0"])
    pol = [g for g in gs if gp[g]["label"] == "POLISHED"]
    fal = [g for g in gs if gp[g]["label"] != "POLISHED"]
    fg_pol = fit_root_law(pol, [gp[g]["c0_lsq_mine"] for g in pol]) if len(pol) >= 3 else None
    fg_fal = fit_root_law(fal, [gp[g]["c0_lsq_mine"] for g in fal]) if len(fal) >= 3 else None
    r32, r16 = (
        gp[32.0]["c0_lsq_mine"] / gp[8.0]["c0_lsq_mine"],
        gp[16.0]["c0_lsq_mine"] / gp[8.0]["c0_lsq_mine"],
    )
    frame = {
        "fit_g_mine": fg,
        "fit_delta_mine": fd,
        "product_fit_mine": pf,
        "label_mine": lab_mine,
        "label_row": res["frame"]["outcome"],
        "fit_g_row": res["frame"]["fit_g"],
        "fit_delta_row": res["frame"]["fit_delta"],
        "product_fit_on_producer_roots": pf_prod,
        "product_max_miss_row": res["frame"]["product_fit"]["max_rel_miss"],
        "fit_g_polished_rows_only": fg_pol,
        "fit_g_falling_rows_only": fg_fal,
        "polished_g": pol,
        "falling_g": fal,
        "ratio_c32_over_c8": r32,
        "ratio_c16_over_c8": r16,
        "g_minus_1_law": {"32/8": (31 / 7) ** 2, "16/8": (15 / 7) ** 2},
        "g2_law": {"32/8": 16.0, "16/8": 4.0},
        "ratios_all": {
            f"{gs[i + 1]:g}/{gs[i]:g}": gp[gs[i + 1]]["c0_lsq_mine"] / gp[gs[i]]["c0_lsq_mine"]
            for i in range(len(gs) - 1)
        },
        "E_by_g": {f"{g:g}": rows[f"S1_Bseed_g{g:g}_d0.3_w1_n32"]["E"] for g in gs},
    }
    A["C8_frame"] = frame
    line(
        "C8e_product_fit_misses_20pct_on_my_roots_AND_on_the_producers_roots_FRAME_NOT_PRODUCT_reproduced",
        pf["max_rel_miss"] > 0.2
        and pf_prod["max_rel_miss"] > 0.2
        and lab_mine == "FRAME_NOT_PRODUCT"
        and res["frame"]["outcome"].startswith("FRAME_NOT_PRODUCT"),
        f"my roots g0 {fg['x0']:+.3f} delta0 {fd['x0']:+.3f} product max miss {pf['max_rel_miss']:.3f}; producer's g0 {res['frame']['fit_g']['g0']:+.3f} delta0 {res['frame']['fit_delta']['delta0']:+.3f} miss {res['frame']['product_fit']['max_rel_miss']:.3f} (mine on their roots {pf_prod['max_rel_miss']:.3f})",
    )
    line(
        "C8f_g_root_from_my_log_least_squares_within_0.3_of_the_producers_minimax_root",
        abs(fg["x0"] - res["frame"]["fit_g"]["g0"]) < 0.3,
        f"g0 mine {fg['x0']:+.3f} (row {res['frame']['fit_g']['g0']:+.3f}), max miss {fg['max_rel_miss']:.3f} (row {res['frame']['fit_g']['max_rel_miss']:.3f}); neither at 1 (within 0.3): {abs(fg['x0'] - 1) > 0.3 and abs(res['frame']['fit_g']['g0'] - 1) > 0.3}",
    )
    line(
        "C8f2_delta_root_within_0.15_of_the_producers_and_far_from_1_both_fits",
        abs(fd["x0"] - res["frame"]["fit_delta"]["delta0"]) < 0.15
        and abs(fd["x0"] - 1) > 0.15
        and abs(res["frame"]["fit_delta"]["delta0"] - 1) > 0.15,
        f"delta0 mine {fd['x0']:+.3f} (row {res['frame']['fit_delta']['delta0']:+.3f}), miss {fd['max_rel_miss']:.3f} (row {res['frame']['fit_delta']['max_rel_miss']:.3f})",
    )
    fm_high = [
        round(n32[f"S1_Bseed_g{g:g}_d0.3_w1_n32"]["fmax_free_mine"], 3) for g in (8.0, 16.0, 32.0)
    ]
    line(
        "C8g_g16_row_separates_the_laws_c32_over_c8_and_c16_over_c8_closer_to_g2_16_4_than_to_g_minus_1_squared_19.6_4.59",
        abs(r32 - 16.0) < abs(r32 - (31 / 7) ** 2) and abs(r16 - 4.0) < abs(r16 - (15 / 7) ** 2),
        f"c32/c8 {r32:.3f} (g^2 16.0, (g-1)^2 19.61), c16/c8 {r16:.3f} (4.0, 4.59); the three rows FALLING at max |G| {fm_high}",
    )
    line(
        "C8h_POLISHED_g_rows_1_2_4_and_FALLING_g_rows_8_16_32_agree_on_g0_within_0.3",
        fg_pol is not None and fg_fal is not None and abs(fg_pol["x0"] - fg_fal["x0"]) < 0.3,
        f"g0 polished rows {fg_pol['x0'] if fg_pol else None:+.3f} (miss {fg_pol['max_rel_miss'] if fg_pol else float('nan'):.3f}); g0 FALLING rows {fg_fal['x0'] if fg_fal else None:+.3f} (miss {fg_fal['max_rel_miss'] if fg_fal else float('nan'):.3f}); successive ratios {[round(x, 2) for x in frame['ratios_all'].values()]}",
    )
    Eg = frame["E_by_g"]
    line(
        "C8i_polished_E_monotone_in_g_over_the_six_g_rows",
        all(Eg[f"{gs[i + 1]:g}"] > Eg[f"{gs[i]:g}"] for i in range(len(gs) - 1))
        or all(Eg[f"{gs[i + 1]:g}"] < Eg[f"{gs[i]:g}"] for i in range(len(gs) - 1)),
        f"E(g) {[round(Eg[f'{g:g}'], 3) for g in gs]} labels {[gp[g]['label'][:8] for g in gs]}",
    )

    # ---- C9 ----
    tr9 = {}
    for t, r in rows.items():
        pt = r["polish"]["trace"]
        Es_ = [x["E"] for x in pt]
        mono = (
            all(Es_[i + 1] <= Es_[i] + 1e-12 * max(abs(Es_[i]), 1.0) for i in range(len(Es_) - 1))
            if len(Es_) > 1
            else True
        )
        mono = mono and (r["polish"]["E_end"] <= r["polish"]["E_in"] + 1e-12)
        tr9[t] = {
            "polish_monotone": mono,
            "iters": r["polish"]["iters"],
            "n_eval": r["polish"]["n_eval"],
            "E_in": r["polish"]["E_in"],
            "E_end": r["polish"]["E_end"],
            "trace_len": len(pt),
        }
        if r["steps"] > 0:
            lab, rl, dec = fire_gate(r["descent"])
            tr9[t].update(
                {
                    "fire_verdict_mine": lab,
                    "fire_verdict_row": r["descent"]["verdict"],
                    "last_third_rel": rl,
                    "decades": dec,
                }
            )
    A["C9_traces"] = tr9
    line(
        "C9a_polish_energy_monotone_non_increasing_every_row",
        all(v["polish_monotone"] for v in tr9.values()),
        f"iters {[v['iters'] for v in tr9.values()]}; evals {[v['n_eval'] for v in tr9.values()]}",
    )
    line(
        "C9b_FIRE_verdicts_reproduced_from_the_traces_all_17_FIRE_rows_all_FALLING",
        sum(1 for v in tr9.values() if "fire_verdict_mine" in v) == 17
        and all(
            v["fire_verdict_mine"] == v["fire_verdict_row"] == "FALLING"
            for v in tr9.values()
            if "fire_verdict_mine" in v
        ),
        "; ".join(
            f"{t[:14]} rel {v['last_third_rel']:.1e} dec {v['decades']:.2f}"
            for t, v in tr9.items()
            if "fire_verdict_mine" in v
        ),
    )
    line(
        "C9c_polish_hit_the_3000_iteration_cap_on_every_row_that_did_not_reach_the_gate_n64_at_200",
        all(
            v["iters"] == (200 if rows[t]["n"] == 64 else 3000)
            for t, v in tr9.items()
            if rows[t]["label"] != "POLISHED"
        ),
        "",
    )

    # ---- the n64 row ----
    v64 = R["Sd_Bfree_g8_d0.3_w1_n64"]
    r64 = rows["Sd_Bfree_g8_d0.3_w1_n64"]
    A["n64"] = {
        "E_mine": v64["E_mine"],
        "E_row": r64["E"],
        "E_fire_mine": v64["E_fire_mine"],
        "E0_mine": v64["E0_mine"],
        "E0_R20_row": R20J["Sd_v4std_n64_g8"]["E"],
        "tube": v64["tube"]["string_read"],
        "flat_mid": v64["tube"]["flat_mid_tension"],
        "deg9": v64["deg9_winding_mine"],
        "face_layer": v64["radial"]["E_pin_shell"],
        "norm_min_rel": v64["norm"]["axis_norm_min_rel"],
        "shell_dev": v64["shell_max_dev_from_seed_polished"],
        "row_summary": res["n64_string"],
    }
    line(
        "C10n_n64_row_E_E_fire_E0_reproduced_rel_1e-10_E0_equals_the_R20_n64_Sd_row",
        v64["rel_E"] < 1e-10
        and v64["rel_E_fire"] < 1e-10
        and rel(v64["E0_mine"], R20J["Sd_v4std_n64_g8"]["E"]) < 1e-10,
        f"E {v64['E_mine']:.4f} fire {v64['E_fire_mine']:.4f} E0 {v64['E0_mine']:.4f}; tube {v64['tube']['string_read']:+.4f} flat {v64['tube']['flat_mid_tension']:+.4f} |deg9| {abs(v64['deg9_winding_mine']):.2f} face {v64['radial']['E_pin_shell']:.3f} norm min {v64['norm']['axis_norm_min_rel']:.3f}",
    )
    line(
        "C10o_n64_S_d_B_free_tube_read_below_0.02_at_200_polish_iterations_but_above_the_n32_free_value_the_string_not_yet_gone_at_n64",
        v64["tube"]["string_read"] < 0.02
        and v64["tube"]["string_read"] > n32["Sd_Bfree_g8_d0.3_w1_n32"]["tube"]["string_read"],
        f"n64 {v64['tube']['string_read']:+.5f} vs n32 {n32['Sd_Bfree_g8_d0.3_w1_n32']['tube']['string_read']:+.5f}; n64 max |G| {r64['polish']['fmax_end']:.2f} (row), FIRE 2000 steps",
    )

    # ---- C10 wordings ----
    W = []
    W.append(
        {
            "wording": f"results re-collected at {res.get('collected_utc')} (mid-audit) with a surface-orientability reader and a box note citing this audit's C5e; the numeric reads did not change, the labels did",
            "numbers": "C3g / C3g2 compare the new reader with mine; C2k, C4b, C4e, C5d compare against the re-collected labels",
            "supported": True,
        }
    )
    W.append(
        {
            "wording": "docstring: 'the three-axis hedgehog under three boundaries at equal residual'",
            "numbers": f"max |G| across the triple: "
            + "; ".join(f"{o} x{v['max_over_min']:.0f}" for o, v in spread.items())
            + f"; 18 of 21 rows FALLING (1.7e-2 to 7.9), 3 POLISHED; the re-collected per-axis labels now say 'not one residual'",
            "supported": False,
        }
    )
    W.append(
        {
            "wording": "docstring: 'POLISHED iff max |G| < 1e-3 ... else the residual is recorded and the row is labeled FALLING'",
            "numbers": "C2a / C2b: reproduced on every n32 row",
            "supported": True,
        }
    )
    W.append(
        {
            "wording": "results.per_axis: BOUNDARY_SET on S1, Sd, S0 (with the residual note)",
            "numbers": "; ".join(
                f"{o}: total diff {v['rel_boundary_diff_total']:.2f}, E(<12) spread {v['rel_spread_E_lt_12']:.2f}"
                for o, v in per.items()
            ),
            "supported": True,
            "caveat": "the seed rows carry the frozen R20 shell (Sd seed shell 46.0 of 55.1, S0 11.5 of 15.5): the 'total' compared is dominated by pinned seed data, not by the relaxed interior; E(<12) differs by 36 to 46 percent, the interior itself is boundary-set",
        }
    )
    W.append(
        {
            "wording": "results.per_axis degree_r9: S_d winding (rank1) kept at -1.0 under all three boundaries; S_1 kept; S_0 lost",
            "numbers": f"my r 9 surface: S1 preimage {per['S1']['deg9']} ; Sd rank1 conflicts {[w9[f'Sd_B{b}_g8_d0.3_w1_n32']['conflicts'] for b in ('seed', 'far', 'free')]} at gaps {[round(w9[f'Sd_B{b}_g8_d0.3_w1_n32']['gap_min_at_conflicts'], 3) for b in ('seed', 'far', 'free')]}, orientation-dependent integers {[w9[f'Sd_B{b}_g8_d0.3_w1_n32']['orientation_dependent'] for b in ('seed', 'far', 'free')]}",
            "supported": False,
            "caveat": "a solid-angle sum over a closed triangulated surface is an integer for ANY vertex signs (the geodesic interpolation is continuous), so the producer's exact -1.0 on a surface where the rank-1 / rank-2 gap closes to 0.01 is a sign-choice integer, not a degree; the S_d 'winding kept' read is not established",
        }
    )
    W.append(
        {
            "wording": f"results.string (re-collected): {res['string']['outcomes'][0]}",
            "numbers": f"tube read free {t_free:+.5f} < 0.02; flat middle {flat_free:+.5f}; rank1 r9 conflicts mine {sdw['conflicts']} (row {res['string'].get('winding_r9_free', [None, None])[1]}), preimage {sdw['preimage']}",
            "supported": True,
            "caveat": "the label fires on the tension clause and the re-collected wording now carries the undefined winding; the pinned string was the frozen seed on the shell (46.0 of 55.1 under B_seed)",
        }
    )
    W.append(
        {
            "wording": "results.escape_routes (re-collected): eight rows, 'a disclination pierces the r 9 cube ... the biaxial-ring route' wherever the surface is non-orientable; the S_0 rows 'no eigenvector carries a winding'",
            "numbers": "; ".join(
                f"{k}: {v['route_by_norm_rule_mine']} norm {v['axis_norm_min_rel']:.3f} crossings {v['crossings']}"
                for k, v in flagged.items()
            ),
            "supported": True,
            "caveat": "the norm-rule routes and the conflict counts reproduce; 'the biaxial-ring route' names the non-orientability, no biaxiality read enters the label; on S1_Bseed g8 and S1_Bfree the non-orientable read (4 and 6 links at gaps 0.009 / 0.012) sits on rows at max |G| 0.12 and 0.42, a residual-state read",
        }
    )
    W.append(
        {
            "wording": f"results.box (re-collected): {res['box']['S1']['outcome']}",
            "numbers": "; ".join(
                f"{o}: scan {v['dE_dlam_mine_o2']:+.2f} vs identity {v['continuum']:+.2f}; face layer {v['split']['face_layer']:+.2f}, interior {v['split']['interior']:+.2f}, r<12 {v['split']['r_lt_12']:+.2f} (identity r<12 {v['split']['continuum_r_lt_12']:+.2f}), R* {v['R_star_mine']:.1f}"
                for o, v in box.items()
            ),
            "supported": True,
            "caveat": "the pre-registered rule gives AXIS_BOX_LIMITED on the same numbers (C5d); the re-collected note attributes the scan's sign to face padding, but my split (C5g / C5h) puts 82 to 97 percent of the negative slope in the interior and the r < 12 core slope is negative on all three axes while the r < 12 identity is positive: not a face-padding artefact alone; the identity assumes the vacuum outside the region and the field is not at the vacuum inside r 12 either",
        }
    )
    W.append(
        {
            "wording": "results.triple: AXES_DEGENERATE (resolution 161.6)",
            "numbers": f"pairwise diffs {[round(x, 2) for x in diffs]} vs 3 x {resol:.1f}",
            "supported": True,
            "caveat": "the resolution is the Sd seed-minus-free gap (53.9), itself the frozen shell energy of the seed row; the label cannot come out DISTINCT with that resolution whatever the free energies are",
        }
    )
    W.append(
        {
            "wording": "results.w1_ladder: SCALE_UNDECIDED; R20_S1_dE_dlambda_reference -5.61",
            "numbers": f"S1 w25 seed R* {q25['R_star']:.2f}, virial {q25['virial']:.2f}, dE/dlam {q25['dE_dlam']:+.2f}",
            "supported": True,
            "caveat": f"the reference used is the polished R20 S1 row (-5.61), not the R20 x4500 row's own read ({lad['dE_dlam_R20_S1_x4500_row']:+.2f}); the dE/dlam clause fails either way",
        }
    )
    W.append(
        {
            "wording": "results.frame: FRAME_NOT_PRODUCT (rows not POLISHED: g8,g16,g32,d0.3,d0.6,d1)",
            "numbers": f"product max miss mine {pf['max_rel_miss']:.3f}; g0 polished-only {fg_pol['x0'] if fg_pol else float('nan'):+.2f} vs FALLING-only {fg_fal['x0'] if fg_fal else float('nan'):+.2f}",
            "supported": True,
            "caveat": "the six g rows mix three POLISHED rows (g 1, 2, 4) with three rows at max |G| 0.12 to 1.0; the two subsets do not share a root (C8h), so the root fit is over rows in different states",
        }
    )
    W.append(
        {
            "wording": "docstring: 'the g 16 row separates the (g - 1)^2 law from the g^2 law'",
            "numbers": f"c32/c8 {r32:.3f}, c16/c8 {r16:.3f} (g^2: 16, 4; (g-1)^2: 19.6, 4.59)",
            "supported": bool(abs(r16 - 4.0) < abs(r16 - (15 / 7) ** 2)),
            "caveat": "on the three FALLING rows; g0 from those rows alone is near 0 (the pure g^2 law), from the polished rows near -0.9",
        }
    )
    W.append(
        {
            "wording": "runs.py FIRE traces: every FIRE row FALLING",
            "numbers": "C9b reproduced",
            "supported": True,
        }
    )
    W.append(
        {
            "wording": "results.n64_string: tension 0.0156, degree -1, FALLING at 1.85",
            "numbers": f"mine {v64['tube']['string_read']:+.4f}, |deg9| {abs(v64['deg9_winding_mine']):.2f}",
            "supported": True,
            "caveat": "200 polish iterations only; the stack's residual at n64 not recomputed here (memory budget)",
        }
    )
    A["C10_wordings"] = W
    line(
        "C10_wordings_listed_and_the_equal_residual_wording_unsupported",
        any(not w["supported"] for w in W),
        f"{len(W)} wordings, {sum(1 for w in W if not w['supported'])} unsupported",
    )

    A["lines"] = LINES
    A["runtime_s"] = round(time.time() - T0, 1)
    A["peak_rss_GB"] = round(rss_gb(), 3)
    with open(OUT_JSON, "w") as f:
        json.dump(_jsonable(A), f, indent=1)
    npass = sum(1 for v in LINES.values() if v["pass"])
    print()
    for k, v in LINES.items():
        print(f"{'PASS' if v['pass'] else 'FAIL'} {k}")
    print(
        f"\n{npass}/{len(LINES)} PASS, runtime {A['runtime_s']} s, peak RSS {A['peak_rss_GB']:.2f} GB, wrote {OUT_JSON}"
    )


if __name__ == "__main__":
    main()
