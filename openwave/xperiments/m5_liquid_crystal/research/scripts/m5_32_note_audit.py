"""M5.32 METHOD NOTE audit: the pre-send gate (independent re-derivation).

Audits findings/m5_32_method_note.md. Everything load-bearing is rebuilt here
from the equations in section 1 of the note, with this file's own ansatz, own
stencil energy, own clock tangent, own kin density, own degree reader and own
paths. The certified stack (m5_21_3_a_4d.py) is imported ONLY as the comparison
target for the calibration rows and as the FIRE driver for the optional re-runs.
Earlier agents' JSONs are the thing under test, never an input.

Stages (each writes/merges into data/m5_32_note_audit.json):
  topo        1b  stabilizer of d4 in SO(1,3)+, pi_1, pi_2
  instrument  1a  what read_charge_from_M actually reads (3x3 block? leading?)
  barrier     1c  own paths that unwind the degree, own energy along them
  boundary    1d  kin of the degree-0 vacuum-interior state (saved + own re-run)
  taper       1e  tapered clock flow, L ladder
  g32         1f  the g = 32 endpoint: own V4 and melt front
  relax_unwound / relax_g32   own re-runs (fields to the scratch directory, slow)
  codemap     check 3, every anchor against the working tree
  equations   check 4, the note's equations against the code
  numbers     check 2, transcribed numbers against artifacts
  finalize    verdict block

Usage: python3 scripts/m5_32_note_audit.py <stage> [<stage> ...]
"""
from __future__ import annotations

import importlib.util
import json
import os
import tempfile
import re
import sys
import time
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")
CKPT = os.path.join(RES, "checkpoints", "m5_32_r10")
# Slow re-run fields (the relaxed unwound state, the g = 32 endpoint) go to a
# scratch directory outside the repo; set NOTE_AUDIT_SCRATCH to choose it.
SCRATCH = os.environ.get(
    "NOTE_AUDIT_SCRATCH",
    os.path.join(tempfile.gettempdir(), "m5_32_note_audit_scratch"))
os.makedirs(SCRATCH, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_note_audit.json")
T0 = time.time()

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def merge(section, obj):
    d = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            d = json.load(f)
    d[section] = obj
    d["runtime_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(d, f, indent=1, default=float)
    log(f"checkpointed section '{section}'")


# ============================ own geometry / ansatz ==========================
def coords(n, L):
    """the certified offset grid x_i = (i - (n-1)/2) h, rebuilt here."""
    h = L / n
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return X, Y, Z, h


def d4_of(s, g, delta):
    return np.diag([-s * g, 1.0, delta, 0.0])


def euler_frame(X, Y, Z):
    """Qh = R3(phi) R2(theta), theta = -atan2(z, rho), from the note's section
    1.3, written out as explicit rotation matrices (no Rodrigues helper)."""
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    th = -np.arctan2(Z, rho)
    c3, s3 = np.cos(phi), np.sin(phi)
    c2, s2 = np.cos(th), np.sin(th)
    # G3[1,2] = -1, G3[2,1] = +1 : R3 = I + sin G3 + (1 - cos) G3^2 acts on (1,2)
    R3 = np.zeros(X.shape + (4, 4)); R3[..., 0, 0] = 1.0
    R3[..., 1, 1], R3[..., 1, 2] = c3, -s3
    R3[..., 2, 1], R3[..., 2, 2] = s3, c3
    R3[..., 3, 3] = 1.0
    # G2[1,3] = +1, G2[3,1] = -1 : rotation in the (1,3) plane
    R2 = np.zeros(X.shape + (4, 4)); R2[..., 0, 0] = 1.0
    R2[..., 1, 1], R2[..., 1, 3] = c2, s2
    R2[..., 3, 1], R2[..., 3, 3] = -s2, c2
    R2[..., 2, 2] = 1.0
    return R3 @ R2


def ansatz(n, L, s=-1.0, g=8.0, delta=0.3):
    X, Y, Z, h = coords(n, L)
    Q = euler_frame(X, Y, Z)
    d4 = d4_of(s, g, delta)
    M = Q @ d4 @ np.swapaxes(Q, -1, -2)
    return 0.5 * (M + np.swapaxes(M, -1, -2)), Q, h


def clock_tangent(Q, d4):
    """a0 = Qh (G1 d4 + d4 G1^T) Qh^T, the note's section 1.3 formula."""
    G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0
    core = G1 @ d4 + d4 @ G1.T
    return Q @ core @ np.swapaxes(Q, -1, -2)


# ============================ own energy / kin ===============================
def jets(M, h, branch):
    """one-sided jets on the periodic-free interior, zero on the missing edge,
    exactly the certified convention (fwd: A[i] = (M[i+1] - M[i]) / h)."""
    out = []
    for ax in range(3):
        A = np.zeros_like(M)
        sl_a = [slice(None)] * 3; sl_b = [slice(None)] * 3
        if branch == "fwd":
            sl_a[ax], sl_b[ax] = slice(1, None), slice(0, -1)
            A[tuple(sl_b)] = (M[tuple(sl_a)] - M[tuple(sl_b)]) / h
        else:
            sl_a[ax], sl_b[ax] = slice(1, None), slice(0, -1)
            A[tuple(sl_a)] = (M[tuple(sl_a)] - M[tuple(sl_b)]) / h
        out.append(A)
    return out


def curv(A, B):
    return A @ ETA @ B - B @ ETA @ A


def bracket(F, G):
    """<F, G>_eta = tr(eta F eta G^T)."""
    return np.einsum("...ab,...ab->...", ETA @ F @ ETA, G)


def eu_density(M, h):
    dens = np.zeros(M.shape[:3])
    for br in ("fwd", "bwd"):
        A = jets(M, h, br)
        for i in range(3):
            for j in range(i + 1, 3):
                F = curv(A[i], A[j])
                dens += 0.5 * 4.0 * bracket(F, F)
    return dens * h ** 3


def v4_density(M, h, s, g, delta, w=7.24023879e-4):
    Me = M @ ETA
    P = np.broadcast_to(np.eye(4), M.shape)
    dens = np.zeros(M.shape[:3])
    for p in range(1, 5):
        P = P @ Me
        Cp = (s * g) ** p + 1.0 + delta ** p
        dens += (np.einsum("...kk->...", P) - Cp) ** 2
    return w * dens * h ** 3


def kin_density(M, a0, h):
    dens = np.zeros(M.shape[:3])
    for br in ("fwd", "bwd"):
        A = jets(M, h, br)
        for i in range(3):
            F = curv(a0, A[i])
            dens += 0.5 * 4.0 * bracket(F, F)
    return dens * h ** 3


# ============================ own degree reader ==============================
def spatial_director(M, which=-1):
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    v = V[..., which]
    return v / np.linalg.norm(v, axis=-1, keepdims=True), lam


def cube_faces(n, k):
    """index bounds [c-k, c+k] on the offset grid; returns the six faces as
    (axis, position, sign) plus the closed index range."""
    c = (n - 1) / 2.0
    lo, hi = int(round(c - k)), int(round(c + k))
    return lo, hi


def lift_raster(vecs_on_cube, lo, hi):
    """orientation lift by raster sweep with a stack (DFS), seeded at the
    (lo, lo, lo) corner; different traversal from the BFS-from-max-gap of the
    R10 audit and from the Prim lift of m5_22_e_audit."""
    idx = {}
    pts = []
    rng = range(lo, hi + 1)
    for i in rng:
        for j in rng:
            for kk in rng:
                if i in (lo, hi) or j in (lo, hi) or kk in (lo, hi):
                    idx[(i, j, kk)] = len(pts); pts.append((i, j, kk))
    N = len(pts)
    V = np.array([vecs_on_cube[p] for p in pts])
    sign = np.zeros(N, dtype=int)
    sign[0] = 1
    stack = [0]
    frustrated = 0
    while stack:
        a = stack.pop()
        p = pts[a]
        for ax in range(3):
            for d in (-1, 1):
                q = list(p); q[ax] += d; q = tuple(q)
                b = idx.get(q)
                if b is None:
                    continue
                dot = float(V[a] @ V[b])
                want = sign[a] if dot >= 0 else -sign[a]
                if sign[b] == 0:
                    sign[b] = want; stack.append(b)
                elif sign[b] != want:
                    frustrated += 1
    return pts, idx, V * sign[:, None], frustrated // 2


def solid_angle_tri(a, b, c):
    num = a @ np.cross(b, c)
    den = 1.0 + a @ b + b @ c + c @ a
    return 2.0 * np.arctan2(num, den)


def cube_degree(W, idx, lo, hi):
    """oriented solid-angle sum over the six faces, outward normals, own
    triangulation (the (u, v+1) diagonal, the other one from the R10 audit)."""
    total = 0.0
    rng = range(lo, hi)
    for ax in range(3):
        bs = [x for x in range(3) if x != ax]
        for pos, orient in ((hi, +1.0), (lo, -1.0)):
            for u in rng:
                for v in rng:
                    def P(du, dv):
                        q = [0, 0, 0]; q[ax] = pos; q[bs[0]] = u + du; q[bs[1]] = v + dv
                        return W[idx[tuple(q)]]
                    p00, p10, p11, p01 = P(0, 0), P(1, 0), P(1, 1), P(0, 1)
                    # face orientation: (e_b0 x e_b1) . e_ax parity
                    par = 1.0 if (bs[0], bs[1]) in ((1, 2), (2, 0), (0, 1)) else -1.0
                    om = solid_angle_tri(p00, p10, p01) + solid_angle_tri(p10, p11, p01)
                    total += par * orient * om
    return total / (4.0 * np.pi)


def read_degree(M, k, which=-1):
    n = M.shape[0]
    v, lam = spatial_director(M, which)
    lo, hi = cube_faces(n, k)
    pts, idx, W, nf = lift_raster(v, lo, hi)
    gap = (lam[..., 2] - lam[..., 1])
    gmin = min(float(gap[p]) for p in pts)
    return float(cube_degree(W, idx, lo, hi)), int(nf), gmin


def degree_controls():
    """hedgehog +1, mirror -1, degree-2 texture +2, uniform 0."""
    n = 17
    X, Y, Z, h = coords(n, 2.0)
    r = np.sqrt(X * X + Y * Y + Z * Z) + 1e-300
    nh = np.stack([X / r, Y / r, Z / r], -1)
    out = {}

    def fake_M(v):
        Mx = np.zeros(X.shape + (4, 4)); Mx[..., 0, 0] = 5.0
        Mx[..., 1:4, 1:4] = np.einsum("...i,...j->...ij", v, v)
        return Mx
    out["hedgehog"] = read_degree(fake_M(nh), 6)[0]
    out["mirror"] = read_degree(fake_M(nh * np.array([1, 1, -1.0])), 6)[0]
    ph, th = np.arctan2(Y, X), np.arccos(np.clip(Z / r, -1, 1))
    d2 = np.stack([np.sin(th) * np.cos(2 * ph), np.sin(th) * np.sin(2 * ph), np.cos(th)], -1)
    out["degree2"] = read_degree(fake_M(d2), 6)[0]
    out["uniform"] = read_degree(fake_M(np.broadcast_to([0, 0, 1.0], nh.shape).copy()), 6)[0]
    return out


# ================================ stages =====================================
def stage_topo():
    """1b: stabilizer of d4 in SO(1,3)+; pi_1 and pi_2 of the orbit."""
    out = {}
    s, g, delta = -1.0, 32.0, 0.3
    d4 = d4_of(s, g, delta)
    # Lie algebra so(1,3): X^T eta + eta X = 0  <=>  X = eta S with S antisymmetric.
    basis = []
    for a in range(4):
        for b in range(a + 1, 4):
            S = np.zeros((4, 4)); S[a, b] = 1.0; S[b, a] = -1.0
            basis.append(ETA @ S)
    # infinitesimal stabilizer: X d4 + d4 X^T = 0  (action M -> L M L^T)
    Amat = np.array([(X @ d4 + d4 @ X.T).ravel() for X in basis]).T
    sv = np.linalg.svd(Amat, compute_uv=False)
    out["algebra_stabilizer_dim_toy"] = int(np.sum(sv < 1e-10))
    out["algebra_singular_values_toy"] = sv.tolist()
    # at delta = 0 the (2,3) rotation must appear in the stabilizer
    d40 = d4_of(s, g, 0.0)
    A0 = np.array([(X @ d40 + d40 @ X.T).ravel() for X in basis]).T
    out["algebra_stabilizer_dim_delta0"] = int(np.sum(np.linalg.svd(A0, compute_uv=False) < 1e-10))
    # discrete stabilizer: L d4 L^T = d4 with L^T eta L = eta  =>  L commutes
    # with d4 eta = diag(-sg... ) (distinct spectrum) => L diagonal signs.
    # Enumerate all 16 sign matrices, keep those in SO(1,3)+ (det +1, L00 > 0).
    stab = []
    for bits in range(16):
        sg = [1 - 2 * ((bits >> i) & 1) for i in range(4)]
        Lm = np.diag(sg).astype(float)
        if np.allclose(Lm @ d4 @ Lm.T, d4) and np.allclose(Lm.T @ ETA @ Lm, ETA) \
                and np.linalg.det(Lm) > 0 and Lm[0, 0] > 0:
            stab.append(sg)
    out["discrete_stabilizer_elements"] = stab
    out["discrete_stabilizer_order"] = len(stab)
    # is it the Klein four-group: abelian, every element squares to identity
    ok = all(np.allclose(np.diag(a) @ np.diag(a), np.eye(4)) for a in stab)
    out["every_element_order_le_2"] = bool(ok)
    # d4 eta spectrum distinct (needed for the 'L diagonal' step)
    out["d4_eta_spectrum"] = np.diag(d4 @ ETA).tolist()
    # a generic non-diagonal L in the stabilizer? random search of L d4 L^T = d4
    # restricted to rotations (the compact part): solve R d4 R^T = d4, R in SO(3)
    rng = np.random.default_rng(1)
    best = 1e9
    for _ in range(2000):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        w, x, y, z = q
        R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                      [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                      [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
        L4 = np.eye(4); L4[1:, 1:] = R
        res = np.abs(L4 @ d4 @ L4.T - d4).max()
        if res < best and not np.allclose(np.abs(R), np.eye(3), atol=1e-2):
            best = res
    out["min_residual_nondiagonal_rotation"] = float(best)
    # pi_1: lift the three nontrivial Klein elements (pi rotations about the
    # three spatial axes) to SU(2) as unit quaternions and generate the group.
    def qmul(p, q):
        a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
        return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)
    gens = [(0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]   # i, j, k = pi rotations
    grp = {(1, 0, 0, 0)}
    frontier = [(1, 0, 0, 0)]
    while frontier:
        nxt = []
        for e in frontier:
            for gq in gens:
                p = tuple(int(round(v)) for v in qmul(e, gq))
                if p not in grp:
                    grp.add(p); nxt.append(p)
        frontier = nxt
    out["pi_1_lift_group_order"] = len(grp)
    out["pi_1_contains_minus_one"] = (-1, 0, 0, 0) in grp
    out["pi_1_is_abelian"] = all(qmul(a, b) == qmul(b, a) for a in grp for b in grp)
    out["pi_1_identification"] = ("order 8, non-abelian, contains -1: the quaternion group Q8"
                                  if len(grp) == 8 and not out["pi_1_is_abelian"]
                                  and out["pi_1_contains_minus_one"] else "NOT Q8")
    out["pi_2_argument"] = ("OPS = SO(1,3)+ / V4 is a quotient of a Lie group by a finite "
                            "subgroup, so pi_2(OPS) = pi_2(SO(1,3)+) = pi_2(SO(3)) = 0 (a "
                            "covering map does not change pi_2; the Lie group has pi_2 = 0).")
    # the OTHER components of V4 = 0: an M with the same power traces of M eta
    # that is NOT in the SO(1,3)+ orbit of d4 (permuted signature).
    other = np.diag([-1.0, s * g * -1.0, delta, 0.0])   # M eta spectrum = {1, -sg... }
    def ptraces(Mx):
        Me = Mx @ ETA; P = np.eye(4); t = []
        for p in range(1, 5):
            P = P @ Me; t.append(np.trace(P))
        return np.array(t)
    out["v4_zero_other_component_example"] = {
        "M": other.tolist(),
        "power_trace_mismatch_vs_d4": float(np.abs(ptraces(other) - ptraces(d4)).max()),
        "same_orbit": bool(np.allclose(np.linalg.eigvalsh(other), np.linalg.eigvalsh(d4))),
        "reading": "V4 = 0 fixes the spectrum of M eta, not the SO(1,3)+ orbit; a second "
                   "component with indefinite M exists, so 'the eigenvalues frozen on the "
                   "OPS' needs 'on the orbit component of d4'"}
    verdict = "CONFIRMED" if (out["discrete_stabilizer_order"] == 4 and ok
                              and out["algebra_stabilizer_dim_toy"] == 0
                              and len(grp) == 8) else "REFUTED"
    out["verdict"] = verdict
    out["deciding_number"] = out["discrete_stabilizer_order"]
    out["note"] = ("stabilizer = 4 diagonal sign matrices with L00 = +1, det +1, all of "
                   "order 2 (Klein); algebra stabilizer dim 0 at delta = 0.3 (1 at delta = 0); "
                   "SU(2) lift closes on Q8 (order 8, non-abelian, contains -1); pi_2 = 0.")
    merge("claim_1b", out)
    return out


def stage_instrument():
    """1a: what does read_charge_from_M read, on which matrix, on this field."""
    AUD22 = _load("aud22", "m5_22_e_audit.py")
    src = open(os.path.join(HERE, "m5_22_e_audit.py")).read().splitlines()
    fn_line = next(i for i, l in enumerate(src, 1) if l.startswith("def read_charge_from_M"))
    body = "\n".join(src[fn_line - 1:fn_line + 8])
    out = {"def_line": fn_line, "body": body,
           "takes_last_eigh_column": "V[..., -1]" in body,
           "uses_only_one_eigenvector": body.count("V[...") == 1}
    n, L = 32, 48.0
    M, Q, h = ansatz(n, L)
    c = (n - 1) / 2.0
    res = {}
    for k in (4, 8, 12):
        lo, hi = int(round(c - k)), int(round(c + k))
        try:
            q4, cf4 = AUD22.read_charge_from_M(M, lo, hi)
        except ValueError as e:
            q4, cf4 = float("nan"), -1
            out["instrument_on_4x4_raises"] = repr(e)
        q3, cf3 = AUD22.read_charge_from_M(np.ascontiguousarray(M[..., 1:4, 1:4]), lo, hi)
        own, nf, gmin = read_degree(M, k)
        own_mid, nf_mid, _ = read_degree(M, k, which=1)
        own_low, nf_low, _ = read_degree(M, k, which=0)
        res[f"hw{k}"] = {"instrument_on_4x4_M": [float(q4), int(cf4)],
                         "instrument_on_3x3_block": [float(q3), int(cf3)],
                         "own_reader_leading": [own, nf],
                         "own_reader_middle": [own_mid, nf_mid],
                         "own_reader_lowest": [own_low, nf_low],
                         "min_top_gap_on_surface": gmin}
    out["reads"] = res
    out["controls"] = degree_controls()
    q4s = [abs(res[k]["instrument_on_4x4_M"][0]) for k in res]
    q4s = [0.0 if np.isnan(x) else x for x in q4s]
    q3s = [abs(res[k]["instrument_on_3x3_block"][0]) for k in res]
    out["deciding_number"] = q3s
    out["verdict"] = "CONFIRMED" if out["takes_last_eigh_column"] and max(q3s) > 0.9 else "REFUTED"
    out["qualifier"] = ("the instrument reads the LEADING eigenvector of the matrix it is "
                        "GIVEN; fed the 4x4 M it reads the constant time axis (|Q| = "
                        f"{max(q4s):.3f}), fed the 3x3 spatial block it reads n-hat (|Q| = "
                        f"{max(q3s):.3f}). The note's 'eigh(M), take V[..., -1]' omits that "
                        "every caller passes the SPATIAL 3x3 block (r6: M[..., 1:4, 1:4]).")
    merge("claim_1a", out)
    return out


def path_energy(M0, M1, wmask, ns, h, s, g, delta, k_read=8, n_read=21):
    """E along M_s = M0 + s * w * (M1 - M0); degree read at n_read points."""
    es = []
    for sv in ns:
        Ms = M0 + sv * wmask * (M1 - M0)
        es.append(float(eu_density(Ms, h).sum() + v4_density(Ms, h, s, g, delta).sum()))
    es = np.array(es)
    qs = []
    for sv in np.linspace(0, 1, n_read):
        Ms = M0 + sv * wmask * (M1 - M0)
        q, nf, gmin = read_degree(Ms, k_read)
        qs.append({"s": float(sv), "Q": q, "frustrated": nf, "min_top_gap": gmin})
    dE = es - es[0]
    return {"E_start": float(es[0]), "E_end": float(es[-1]),
            "barrier_max_dE": float(dE.max()),
            "max_uphill_step": float(np.max(np.diff(es))),
            "monotone_decreasing": bool(np.all(np.diff(es) <= 1e-9)),
            "n_points": len(ns), "degree_track": qs,
            "Q_start": qs[0]["Q"], "Q_end": qs[-1]["Q"]}


def stage_barrier():
    """1c: the unwinding barrier, own ansatz, own energy, own paths."""
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    n, L, s, g, delta = 32, 48.0, -1.0, 8.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    X, Y, Z, _ = coords(n, L)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    d4 = d4_of(s, g, delta)
    cfg = INS4.base_cfg(s=s, g=g, n=n, L=L, delta=delta)
    out = {}
    # calibration: own ansatz and own energy against the certified stack
    eu_c, ev_c = INS4.e_parts(M, cfg)
    eu_o, ev_o = float(eu_density(M, h).sum()), float(v4_density(M, h, s, g, delta).sum())
    B8 = _load("b8", "m5_21_8_b_lattice.py")
    Mb = B8.dressed(cfg, 0.0)
    out["calibration"] = {
        "own_ansatz_vs_B8_dressed_max_abs": float(np.abs(M - Mb).max()),
        "own_E_u": eu_o, "certified_E_u": float(eu_c),
        "own_V4": ev_o, "certified_V4": float(ev_c),
        "rel_dev_E_u": float(abs(eu_o / eu_c - 1.0))}
    Mv = np.broadcast_to(d4, M.shape)
    ns = np.linspace(0.0, 1.0, 201)
    paths = {}
    # (i) the R10 audit's window and four others, linear melt, own fine grid
    for (r_in, r_out) in ((15, 21), (9, 21), (6, 12), (3, 9), (18, 21)):
        w = np.clip((r_out - R) / (r_out - r_in), 0.0, 1.0)[..., None, None]
        paths[f"linear_melt_{r_in}_{r_out}"] = path_energy(M, Mv, w, ns, h, s, g, delta)
        log(f"  window ({r_in},{r_out}): E {paths[f'linear_melt_{r_in}_{r_out}']['E_start']:.4f} "
            f"-> {paths[f'linear_melt_{r_in}_{r_out}']['E_end']:.4f} barrier "
            f"{paths[f'linear_melt_{r_in}_{r_out}']['barrier_max_dE']:.4g} "
            f"Q {paths[f'linear_melt_{r_in}_{r_out}']['Q_start']:+.3f} -> "
            f"{paths[f'linear_melt_{r_in}_{r_out}']['Q_end']:+.3f}")
    # (ii) a smooth (cosine) melt profile and a full-box melt that keeps the pin
    w = np.where(R <= 15, 1.0, np.where(R >= 21, 0.0, 0.5 * (1 + np.cos(np.pi * (R - 15) / 6))))
    paths["cosine_melt_15_21"] = path_energy(M, Mv, w[..., None, None], ns, h, s, g, delta)
    pin = INS4.pin_shell(n, h)
    w = (~pin).astype(float)[..., None, None]
    paths["melt_everything_unpinned"] = path_energy(M, Mv, w, ns, h, s, g, delta)
    # (iii) the physically relevant one: from the RELAXED degree-1 state to the
    # RELAXED degree-0 state (both R10 audit endpoints, same pin, same budget)
    relaxed = {}
    for tag in ("main3000", "unwound3000", "main12000"):
        p = os.path.join(CKPT, f"aud_{tag}.npy")
        if os.path.exists(p):
            relaxed[tag] = np.load(p)
    if "main3000" in relaxed and "unwound3000" in relaxed:
        A, B = relaxed["main3000"], relaxed["unwound3000"]
        one = np.ones(M.shape[:3])[..., None, None]
        paths["relaxed_deg1_to_relaxed_deg0_linear"] = path_energy(A, B, one, ns, h, s, g, delta)
        out["relaxed_endpoints_boundary_identical"] = float(np.abs((A - B)[pin]).max())
        # and the linear melt of the RELAXED degree-1 state itself
        w = np.clip((21 - R) / 6.0, 0.0, 1.0)[..., None, None]
        paths["relaxed_deg1_linear_melt_15_21"] = path_energy(A, Mv, w, ns, h, s, g, delta)
        w = np.clip((12 - R) / 6.0, 0.0, 1.0)[..., None, None]
        paths["relaxed_deg1_linear_melt_6_12"] = path_energy(A, Mv, w, ns, h, s, g, delta)
    if "main12000" in relaxed and "unwound3000" in relaxed:
        one = np.ones(M.shape[:3])[..., None, None]
        paths["relaxed12000_deg1_to_relaxed_deg0_linear"] = path_energy(
            relaxed["main12000"], relaxed["unwound3000"], one, ns, h, s, g, delta)
    for k, v in paths.items():
        log(f"  {k}: E {v['E_start']:.4f} -> {v['E_end']:.4f} barrier {v['barrier_max_dE']:.4g} "
            f"monotone {v['monotone_decreasing']} Q {v['Q_start']:+.3f} -> {v['Q_end']:+.3f}")
    out["paths"] = paths
    rigid_ok = all(paths[k]["barrier_max_dE"] <= 1e-9 and paths[k]["monotone_decreasing"]
                   for k in paths if k.startswith("linear_melt"))
    out["deciding_number"] = paths["linear_melt_15_21"]["barrier_max_dE"]
    out["E_end_15_21"] = paths["linear_melt_15_21"]["E_end"]
    rel = paths.get("relaxed_deg1_to_relaxed_deg0_linear", {})
    out["verdict"] = "CONFIRMED" if rigid_ok else "REFUTED"
    out["qualifier"] = (
        "the zero barrier is a statement about paths that START AT THE RIGID ANSATZ "
        f"(E = {paths['linear_melt_15_21']['E_start']:.3f}), which is not a stationary point; "
        f"from the RELAXED degree-1 endpoint (E = {rel.get('E_start', float('nan')):.3f}) the "
        f"straight line to the relaxed degree-0 endpoint has max dE = "
        f"{rel.get('barrier_max_dE', float('nan')):.4f} and is "
        f"{'monotone' if rel.get('monotone_decreasing') else 'NOT monotone'}.")
    merge("claim_1c", out)
    return out


def vacuum_interior(M, R, r_max, d4):
    m = R <= r_max
    return float(np.abs(M[m] - d4).max())


def stage_boundary():
    """1d: the degree-0 vacuum-interior state's kin (saved endpoint, own kin)."""
    n, L, s, g, delta = 32, 48.0, -1.0, 8.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    X, Y, Z, _ = coords(n, L)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    d4 = d4_of(s, g, delta)
    a0 = clock_tangent(Q, d4)
    B8 = _load("b8", "m5_21_8_b_lattice.py")
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    cfg = INS4.base_cfg(s=s, g=g, n=n, L=L, delta=delta)
    out = {"a0_own_vs_B8_fd_max_abs": float(np.abs(a0 - B8.a0_unit(cfg, 0.0)).max()),
           "kin_rigid_own": float(kin_density(M, a0, h).sum()),
           "kin_rigid_certified": float(INS4.kin_of(M, a0, cfg))}
    fields = {}
    for tag in ("main3000", "unwound3000"):
        p = os.path.join(CKPT, f"aud_{tag}.npy")
        if os.path.exists(p):
            fields[tag] = np.load(p)
    own = os.path.join(SCRATCH, "note_audit_unwound_own.npy")
    if os.path.exists(own):
        fields["own_rerun_unwound3000"] = np.load(own)
    rows = {}
    for tag, F in fields.items():
        kd = kin_density(F, a0, h)
        q8, nf8, _ = read_degree(F, 8)
        q4, nf4, _ = read_degree(F, 4)
        q12, nf12, _ = read_degree(F, 12)
        rows[tag] = {"kin_own": float(kd.sum()), "kin_certified": float(INS4.kin_of(F, a0, cfg)),
                     "E_u_own": float(eu_density(F, h).sum()),
                     "V4_own": float(v4_density(F, h, s, g, delta).sum()),
                     "max_dev_from_vacuum_r_le_15": vacuum_interior(F, R, 15.0, d4),
                     "max_dev_from_vacuum_r_le_12": vacuum_interior(F, R, 12.0, d4),
                     "kin_inside_r15": float(kd[R <= 15].sum()),
                     "kin_r_gt_21": float(kd[R > 21].sum()),
                     "Q_hw4": q4, "Q_hw8": q8, "Q_hw12": q12}
        log(f"  {tag}: kin own {rows[tag]['kin_own']:.3f} cert {rows[tag]['kin_certified']:.3f} "
            f"vac-dev r<=15 {rows[tag]['max_dev_from_vacuum_r_le_15']:.3g} Q8 {q8:+.3f}")
    out["fields"] = rows
    if "main3000" in rows and "unwound3000" in rows:
        out["ratio_unwound_over_main"] = rows["unwound3000"]["kin_own"] / rows["main3000"]["kin_own"]
    out["deciding_number"] = rows.get("unwound3000", {}).get("kin_own")
    ok = ("unwound3000" in rows and abs(rows["unwound3000"]["kin_own"] - 272.20) < 0.05
          and rows["unwound3000"]["max_dev_from_vacuum_r_le_15"] < 1e-6
          and abs(rows["unwound3000"]["Q_hw8"]) < 0.01)
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    out["qualifier"] = ("'vacuum interior out to r = 15' is exact only for the SEED; after "
                        "3000 FIRE iterations the interior deviation is reported above. The "
                        "78 % compares two UNCONVERGED 3000-iteration endpoints.")
    merge("claim_1d", out)
    return out


def stage_taper():
    """1e: tapered clock flow and its L ladder, own kin density."""
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    s, g, delta = -1.0, 8.0, 0.3
    d4 = d4_of(s, g, delta)
    out = {"boxes": {}}
    for tag, n, L in (("b16_3000", 16, 24.0), ("b24_3000", 24, 36.0), ("main3000", 32, 48.0),
                      ("b40_3000", 40, 60.0)):
        p = os.path.join(CKPT, f"aud_{tag}.npy")
        if not os.path.exists(p):
            continue
        F = np.load(p)
        M, Q, h = ansatz(n, L, s, g, delta)
        X, Y, Z, _ = coords(n, L)
        R = np.sqrt(X * X + Y * Y + Z * Z)
        a0 = clock_tangent(Q, d4)
        row = {"kin_full_relaxed": float(kin_density(F, a0, h).sum()),
               "kin_full_rigid": float(kin_density(M, a0, h).sum())}
        for rt in (6, 9, 12, 15, 18, 24):
            hard = a0 * (R <= rt)[..., None, None]
            row[f"hard_cut_{rt}"] = float(kin_density(F, hard, h).sum())
            lin = a0 * np.clip((rt + 3 - R) / 3.0, 0, 1)[..., None, None]
            row[f"linear_taper_{rt}_{rt + 3}"] = float(kin_density(F, lin, h).sum())
        # the density-mask reading: kin density restricted to r <= 12 (no taper
        # on the field, just the integration window)
        kd = kin_density(F, a0, h)
        row["kin_density_inside_12"] = float(kd[R <= 12].sum())
        out["boxes"][f"L{L:g}"] = row
        log(f"  L={L:g}: full {row['kin_full_relaxed']:.3f} hard12 {row['hard_cut_12']:.3f} "
            f"lin12 {row['linear_taper_12_15']:.3f} dens12 {row['kin_density_inside_12']:.3f}")
    b = out["boxes"]
    if "L48" in b:
        out["fraction_hard_cut_12_L48"] = b["L48"]["hard_cut_12"] / b["L48"]["kin_full_relaxed"]
        out["fraction_density_inside_12_L48"] = (b["L48"]["kin_density_inside_12"]
                                                 / b["L48"]["kin_full_relaxed"])
    Ls = sorted(b, key=lambda k: float(k[1:]))
    out["hard_cut_12_by_L"] = {k: b[k]["hard_cut_12"] for k in Ls}
    out["full_by_L"] = {k: b[k]["kin_full_relaxed"] for k in Ls}
    vals = [b[k]["hard_cut_12"] for k in Ls if k != "L24"]
    out["hard_cut_12_spread_L36_up"] = float((max(vals) - min(vals)) / np.mean(vals)) if vals else None
    out["deciding_number"] = out.get("fraction_hard_cut_12_L48")
    ok = (out.get("fraction_hard_cut_12_L48") is not None
          and abs(out["fraction_hard_cut_12_L48"] - 0.329) < 0.003)
    out["verdict"] = "CONFIRMED" if ok else "QUALIFIED"
    out["qualifier"] = ("L-independence is tested here on the RELAXED 3000-iteration endpoints "
                        "of L = 24, 36, 48, 60 with the clock cut at r = 12; in the L = 24 box "
                        "r = 12 reaches the pinned shell (pin from |x| = 9) so that box is not "
                        "a clean member of the ladder.")
    merge("claim_1e", out)
    return out


def top_gap_front(F, R, thr=0.35):
    lam = np.linalg.eigvalsh(F[..., 1:4, 1:4])
    gap = lam[..., 2] - lam[..., 1]
    m = gap < thr
    return float(R[m].max()) if m.any() else 0.0, float(gap.min()), int(m.sum())


def stage_g32():
    """1f: the g = 32 endpoint, own V4 and own melt front."""
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    n, L, s, delta = 16, 24.0, -1.0, 0.3
    out = {}
    for tag, g, p in (("audit_g32_b16_64k", 32.0, os.path.join(CKPT, "aud_g32_b16_64k.npy")),
                      ("audit_g8_b16_1k", 8.0, os.path.join(CKPT, "aud_g8_b16_1k.npy")),
                      ("own_rerun_g32_b16_64k", 32.0,
                       os.path.join(SCRATCH, "note_audit_g32_own.npy"))):
        if not os.path.exists(p):
            continue
        F = np.load(p)
        M, Q, h = ansatz(n, L, s, g, delta)
        X, Y, Z, _ = coords(n, L)
        R = np.sqrt(X * X + Y * Y + Z * Z)
        d4 = d4_of(s, g, delta)
        a0 = clock_tangent(Q, d4)
        front, gmin, ncells = top_gap_front(F, R)
        row = {"E_u_own": float(eu_density(F, h).sum()),
               "V4_own": float(v4_density(F, h, s, g, delta).sum()),
               "kin_own": float(kin_density(F, a0, h).sum()),
               "kin_rigid_own": float(kin_density(M, a0, h).sum()),
               "front_r_top_gap_below_0.35": front, "min_top_gap": gmin,
               "n_cells_top_gap_below_0.35": ncells,
               "E_u_rigid_own": float(eu_density(M, h).sum())}
        row["kin_over_rigid"] = row["kin_own"] / row["kin_rigid_own"]
        cfg = INS4.base_cfg(s=s, g=g, n=n, L=L, delta=delta)
        G = INS4.grad(F, cfg) * (~INS4.pin_shell(n, h))[..., None, None]
        row["fmax_own_recomputed"] = float(np.abs(G).max())
        out[tag] = row
        log(f"  {tag}: V4 {row['V4_own']:.5f} front {front:.3f} gap_min {gmin:.4f} "
            f"kin/rigid {row['kin_over_rigid']:.4f} fmax {row['fmax_own_recomputed']:.2f}")
    a = out.get("audit_g32_b16_64k", {})
    out["deciding_number"] = a.get("V4_own")
    ok = a and abs(a["V4_own"] - 0.00097) < 2e-5 and a["front_r_top_gap_below_0.35"] == 0.0
    out["verdict"] = "CONFIRMED" if ok else "REFUTED"
    out["qualifier"] = ("V4 = 0.00097 is 4 x 10^-6 of the g = 32 rigid E_u; the endpoint is "
                        "unconverged (fmax recomputed above) and the note's 'melt-front radius "
                        "0.000' means no cell has top gap < 0.35 at n = 16, h = 1.5, "
                        "where the g = 8 control at matched effective time reads 2.487.")
    merge("claim_1f", out)
    return out


# ------------------------------ own re-runs ---------------------------------
def relax_unwound():
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    n, L, s, g, delta = 32, 48.0, -1.0, 8.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    X, Y, Z, _ = coords(n, L)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    d4 = d4_of(s, g, delta)
    w = np.clip((21.0 - R) / 6.0, 0.0, 1.0)[..., None, None]
    M0 = M + w * (np.broadcast_to(d4, M.shape) - M)
    cfg = INS4.base_cfg(s=s, g=g, n=n, L=L, delta=delta)
    free = ~INS4.pin_shell(n, h)
    Mx, info = INS4.fire(M0, cfg, free, max_iter=3000, log_every=500, tag="note_audit_unwound",
                         dt0=0.01, dt_max=0.1, plateau=(2000, 1e-10))
    np.save(os.path.join(SCRATCH, "note_audit_unwound_own.npy"), Mx)
    a0 = clock_tangent(Q, d4)
    rec = {"stop": info["stop"], "wall_s": info["wall_s"],
           "E_u_end_own": float(eu_density(Mx, h).sum()),
           "V4_end_own": float(v4_density(Mx, h, s, g, delta).sum()),
           "kin_end_own": float(kin_density(Mx, a0, h).sum())}
    json.dump(rec, open(os.path.join(SCRATCH, "rerun_unwound.json"), "w"), indent=1)


def relax_g32():
    INS4 = _load("ins4", "m5_21_3_a_4d.py")
    n, L, s, g, delta = 16, 24.0, -1.0, 32.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    cfg = INS4.base_cfg(s=s, g=g, n=n, L=L, delta=delta)
    free = ~INS4.pin_shell(n, h)
    Mx, info = INS4.fire(M, cfg, free, max_iter=64000, log_every=4000, tag="note_audit_g32",
                         dt0=1.5e-4, dt_max=1.5e-3, plateau=(10 ** 9, 0.0))
    np.save(os.path.join(SCRATCH, "note_audit_g32_own.npy"), Mx)
    d4 = d4_of(s, g, delta)
    a0 = clock_tangent(Q, d4)
    rec = {"stop": info["stop"], "wall_s": info["wall_s"], "trace_tail": info["trace"][-3:],
           "E_u_end_own": float(eu_density(Mx, h).sum()),
           "V4_end_own": float(v4_density(Mx, h, s, g, delta).sum()),
           "kin_end_own": float(kin_density(Mx, a0, h).sum())}
    json.dump(rec, open(os.path.join(SCRATCH, "rerun_g32.json"), "w"), indent=1)


# ------------------------------ document checks ------------------------------
NOTE = os.path.join(RES, "findings", "m5_32_method_note.md")


def stage_codemap():
    txt = open(NOTE).read()
    rows = []
    pat = re.compile(r"scripts/(m5_[0-9a-z_]+\.py)#L(\d+)(?:-L(\d+))?")
    # walk the table rows of section 2
    sec = txt.split("## 2. Equation-to-code map")[1].split("## 3.")[0]
    cur_file = None
    for line in sec.splitlines():
        if not line.startswith("|") or line.startswith("| ---") or "Function" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 3:
            continue
        funcs = re.findall(r"`([A-Za-z_0-9]+)`", cells[1])
        links = re.findall(r"\[`([^`]+)`\]\((https://[^)]+)\)", cells[2])
        for label, url in links:
            m = re.search(r"scripts/(m5_[0-9a-z_]+\.py)#L(\d+)(?:-L(\d+))?", url)
            if not m:
                rows.append({"object": cells[0], "url": url, "verdict": "NO ANCHOR"})
                continue
            fname, l0, l1 = m.group(1), int(m.group(2)), int(m.group(3) or m.group(2))
            path = os.path.join(HERE, fname)
            if not os.path.exists(path):
                rows.append({"object": cells[0], "file": fname, "verdict": "FILE MISSING"})
                continue
            src = open(path).read().splitlines()
            window = "\n".join(src[l0 - 1:l1])
            found = [f for f in funcs if re.search(rf"^(def|class)\s+{f}\b|^{f}\s*=", window, re.M)]
            # which defs / assignments actually sit at those lines
            at = [re.match(r"^(def|class)\s+(\w+)|^([A-Z_0-9]+)\s*=", l) for l in src[l0 - 1:l1]]
            at = [m2.group(2) or m2.group(3) for m2 in at if m2]
            # where the named functions really are
            real = {}
            for f in funcs:
                for i, l in enumerate(src, 1):
                    if re.match(rf"^(def|class)\s+{f}\b|^{f}\s*=", l):
                        real[f] = i
            rows.append({"object": cells[0], "file": fname, "anchor": [l0, l1],
                         "named": funcs, "found_at_anchor": found, "defs_at_anchor": at,
                         "real_lines": real, "label": label,
                         "verdict": "OK" if found else "ANCHOR OFF"})
    out = {"rows": rows, "n_rows": len(rows),
           "n_off": sum(r["verdict"] != "OK" for r in rows)}
    merge("code_map_check_raw", out)
    for r in rows:
        if r["verdict"] != "OK":
            log(f"  CODEMAP {r['verdict']}: {r.get('object')} {r.get('file')} {r.get('anchor')} "
                f"named {r.get('named')} defs_at_anchor {r.get('defs_at_anchor')} real {r.get('real_lines')}")
    return out


def stage_equations():
    """check 4: the note's section 1 equations against the registry module."""
    import sympy as sp
    LG = _load("lg", "m5_32_lagrangian.py")
    out = {}
    rng = np.random.default_rng(7)
    # (a) I1 normalization: registry I1 density vs (1/2) F_abcd F^abcd on random jets
    A = rng.normal(size=(4, 4, 4)) * 0.5   # A[mu] 4x4, raw contravariant
    F = np.zeros((4, 4, 4, 4))
    for mu in range(4):
        for nu in range(4):
            F[mu, nu] = A[mu] @ ETA @ A[nu] - A[nu] @ ETA @ A[mu]
    # F_abcd F^abcd with eta on all four indices (derivative pair AND internal pair)
    full = np.einsum("mnab,MNAB,mM,nN,aA,bB->", F, F, ETA, ETA, ETA, ETA)
    half = 0.5 * full
    # the note's I1 = sum_{mu<nu} eta^mumu eta^nunu <F_mu nu, F_mu nu>_eta
    note = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            note += ETA[mu, mu] * ETA[nu, nu] * np.einsum("ab,ab->", ETA @ F[mu, nu] @ ETA, F[mu, nu])
    out["I1_note_form_vs_half_FF"] = {"note_form": float(note), "half_F_abcd_Fabcd": float(half),
                                      "rel": float(abs(note / half - 1))}
    # the registry's numpy I1 on a lattice with only spatial jets: compare to note form
    n, L, s, g, delta = 8, 12.0, -1.0, 8.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    d4 = d4_of(s, g, delta)
    a0 = clock_tangent(Q, d4)
    cfg = LG.default_params(n=n, L=L, g=g, s=s, delta=delta) if hasattr(LG, "default_params") else None
    out["registry_has"] = {k: (LG.REGISTRY[k].definition if hasattr(LG.REGISTRY[k], "definition")
                               else str(LG.REGISTRY[k])) for k in LG.REGISTRY}
    out["CERTIFIED_COEFFS"] = LG.CERTIFIED_COEFFS
    out["W1"] = float(LG.W1)
    # (b) V4: registry density vs the note's formula, random symmetric M
    Mr = rng.normal(size=(3, 3, 3, 4, 4)); Mr = 0.5 * (Mr + np.swapaxes(Mr, -1, -2))
    try:
        p = LG.default_params(g=g, s=s, delta=delta)
        v_reg = LG.v4_density_np(None, Mr, p)
        v_own = v4_density(Mr, 1.0, s, g, delta, w=1.0)
        # v4_density_np may or may not include W1 / h^3; report the ratio
        ratio = float(np.mean(v_reg / v_own))
        out["V4_registry_over_own_unit_w"] = {"mean_ratio": ratio,
                                              "ratio_spread": float(np.std(v_reg / v_own) / abs(ratio)),
                                              "equals_W1": bool(abs(ratio / LG.W1 - 1) < 1e-9),
                                              "equals_1": bool(abs(ratio - 1) < 1e-9)}
        out["c4_registry"] = [float(x) for x in LG.c4_of(p)]
        out["c4_note"] = [(s * g) ** k + 1 + delta ** k for k in range(1, 5)]
    except Exception as e:  # pragma: no cover
        out["V4_registry_check_error"] = repr(e)
    # (c) certified energy vs E_cert = 4 (U + omega^2 T) + V4 on the small box
    try:
        p = LG.default_params(g=g, s=s, delta=delta)
        cfgl = {"n": n, "L": L, "h": h, "stencil": "sym"}
        om = 0.37
        e_reg = float(LG.certified_energy(M, cfgl, p, a0=a0, omega=om))
        U = float(eu_density(M, h).sum()) / 4.0
        T = float(kin_density(M, a0, h).sum()) / 4.0
        V = float(v4_density(M, h, s, g, delta).sum())
        out["E_cert_check"] = {"registry": e_reg, "own_4(U+om2T)+V4": 4 * (U + om * om * T) + V,
                               "rel": abs(e_reg / (4 * (U + om * om * T) + V) - 1)}
        # Legendre: term_hamiltonian(I1) at omega equals C omega^2 - A
        A_, B_, C_ = LG.omega_decompose(LG.REGISTRY["I1"], M, cfgl, p, a0)
        H = float(LG.term_hamiltonian(LG.REGISTRY["I1"], M, cfgl, p, a0=a0, omega=om))
        out["legendre_I1"] = {"A": float(np.sum(A_)) if np.ndim(A_) else float(A_),
                              "B": float(np.sum(B_)) if np.ndim(B_) else float(B_),
                              "C": float(np.sum(C_)) if np.ndim(C_) else float(C_),
                              "H_reg": H}
        Asum = out["legendre_I1"]["A"]; Csum = out["legendre_I1"]["C"]
        out["legendre_I1"]["C_om2_minus_A"] = Csum * om * om - Asum
        out["legendre_I1"]["rel"] = abs(H / (Csum * om * om - Asum) - 1)
        # kin = -4 x C(I1) ?
        out["kin_equals_minus4C"] = {"minus4C": -4 * Csum, "own_kin": 4 * T,
                                     "rel": abs(-4 * Csum / (4 * T) - 1)}
    except Exception as e:  # pragma: no cover
        out["E_cert_check_error"] = repr(e)
    # (d) the K_T u-frame identity on random jets at the vacuum
    try:
        KT = _load("kt", "m5_32_r7_a_kt_form.py")
        src = open(os.path.join(HERE, "m5_32_r7_a_kt_form.py")).read()
        out["kt_doc_excerpt"] = re.findall(r"K_T[^\n]*", src)[:6]
    except Exception as e:
        out["kt_load_error"] = repr(e)
    # (e) a0 from the note's formula vs the finite-difference flow: already in 1d
    # (f) the Legendre quartic relation: H = C2 om^2 + 3 C4 om^4 - A from
    #     L = A + C2 om^2 + C4 om^4, H = om dL/dom - L
    om_s = sp.Symbol("omega")
    A_s, C2_s, C4_s = sp.symbols("A C2 C4")
    Ls = A_s + C2_s * om_s ** 2 + C4_s * om_s ** 4
    Hs = sp.expand(om_s * sp.diff(Ls, om_s) - Ls)
    out["legendre_quartic_symbolic"] = {"H": str(Hs),
                                        "matches_note": Hs == sp.expand(C2_s * om_s ** 2 + 3 * C4_s * om_s ** 4 - A_s)}
    Bs = sp.Symbol("B")
    Lq = A_s + Bs * om_s + C2_s * om_s ** 2
    Hq = sp.expand(om_s * sp.diff(Lq, om_s) - Lq)
    out["legendre_quadratic_symbolic"] = {"H": str(Hq), "note_says": "C omega^2 - A",
                                          "matches_note": Hq == sp.expand(C2_s * om_s ** 2 - A_s),
                                          "B_term_in_H": str(Hq.coeff(Bs))}
    J = sp.Symbol("J")
    om_star = sp.solve(sp.Eq(sp.diff(Lq, om_s), J), om_s)[0]
    out["fixed_J_omega_star"] = {"solved": str(om_star), "note_says": "(J - B) / (2 C)"}
    # (g) six quartic density names present in r8_a
    src8 = open(os.path.join(HERE, "m5_32_r8_a_quartics.py")).read()
    out["quartic_defs_present"] = {k: bool(re.search(rf"^def {k}\b", src8, re.M))
                                   for k in ("d_I1", "d_I4", "d_Fpair", "d_C6a", "d_C6b", "d_BI")}
    out["QUARTICS_present"] = bool(re.search(r"^QUARTICS\s*=", src8, re.M))
    merge("equation_check_raw", out)
    return out


def stage_numbers():
    """check 2: the note's transcribed numbers against the artifacts."""
    J = lambda f: json.load(open(os.path.join(DATA, f)))
    r10 = J("m5_32_r10_audit.json")
    led = J("m5_32_ledger.json")
    rows = []

    def row(where, note_val, art_val, src, tol=5e-4):
        try:
            ok = abs(float(note_val) - float(art_val)) <= tol * max(1.0, abs(float(art_val)))
        except Exception:
            ok = note_val == art_val
        rows.append({"where": where, "note": note_val, "artifact": art_val, "source": src,
                     "verdict": "OK" if ok else "MISMATCH"})
    bc = r10["R3_detail"]["branch_compare"]
    row("R10 unwinding E_start 62.852", 62.852, r10["R5_detail"]["unwinding_path"]["E_start"],
        "r10_audit.R5_detail.unwinding_path.E_start")
    row("R10 unwinding E_end 14.794", 14.794, r10["zero_barrier_robustness"]["E_at_s1"][0],
        "r10_audit.zero_barrier_robustness.E_at_s1[0] (a typed constant in the script)")
    row("R10 barrier 0.0", 0.0, r10["R5_detail"]["unwinding_path"]["barrier"],
        "r10_audit.R5_detail.unwinding_path.barrier", tol=1e-12)
    row("R10 degree-0 kin 272.20", 272.20, bc["degree0_it3000"]["kin_total"], "branch_compare.degree0_it3000.kin_total")
    row("R10 relaxed kin 351.17", 351.17, bc["degree1_it3000"]["kin_total"], "branch_compare.degree1_it3000.kin_total")
    row("R10 78 %", 78, 100 * bc["degree0_it3000"]["kin_total"] / bc["degree1_it3000"]["kin_total"],
        "ratio of the two", tol=0.01)
    ct = r10["clock_taper"]
    row("R10 taper r=12 32.9 %", 32.9, 100 * ct["kin_relaxed3000"][2] / ct["kin_relaxed3000"][-1],
        "clock_taper (typed constants)", tol=0.005)
    row("R10 middle conflicts 37,30,87", [37, 30, 87], r10["R5_detail"]["per_eigenvector_rigid"]["middle_conflicts"],
        "R5_detail.per_eigenvector_rigid.middle_conflicts")
    sg = r10["scope_g32"]
    row("g32 V4 0.00097", 0.00097, sg["g32_64000it"]["V4"], "scope_g32 (typed constants)", tol=1e-3)
    row("g32 front 0.000", 0.0, sg["g32_64000it"]["front_r_top_below_0.35"], "scope_g32", tol=1e-12)
    g32rec = J("m5_32_r10_audit_relax_g32_b16_64k.json")
    row("g32 V4 vs relax record", 0.00097, g32rec["V4_end"], "relax_g32_b16_64k.V4_end", tol=2e-2)
    row("g32 fmax 72.8", 72.8, g32rec["trace"][-1]["fmax"], "relax_g32_b16_64k.trace[-1].fmax", tol=2e-3)
    # 12000-iteration decrements -0.185 then -0.222 (the producer's withdrawn claim)
    dl = r10["R3_detail"]["direct_ladder_slopes"]
    row("slope decrements -0.185 (producer) vs audit pairwise", -0.185,
        None, "NOT in the audit JSON; the audit's pair_36_48 decrements are "
              f"{r10['R3_detail']['pair_dependence']['pair_36_48_decrements']}", tol=1e9)
    # R7
    r7 = J("m5_32_r7_audit.json")
    rows.append({"where": "R7 verdict counts 8/5/0", "note": [8, 5, 0],
                 "artifact": [r7.get("n_confirmed"), r7.get("n_qualified"), r7.get("n_refuted")],
                 "source": "r7_audit", "verdict": "OK" if [r7.get("n_confirmed"), r7.get("n_qualified"), r7.get("n_refuted")] == [8, 5, 0] else "MISMATCH"})
    r8 = J("m5_32_r8_audit.json")
    rows.append({"where": "R8 verdict counts 4/4/2", "note": [4, 4, 2],
                 "artifact": [r8.get("n_confirmed"), r8.get("n_qualified"), r8.get("n_refuted")],
                 "source": "r8_audit", "verdict": "OK" if [r8.get("n_confirmed"), r8.get("n_qualified"), r8.get("n_refuted")] == [4, 4, 2] else "MISMATCH"})
    r9 = J("m5_32_r9_audit.json")
    rows.append({"where": "R9 verdict counts 5/2/2", "note": [5, 2, 2],
                 "artifact": [r9.get("n_confirmed"), r9.get("n_qualified"), r9.get("n_refuted")],
                 "source": "r9_audit", "verdict": "OK" if [r9.get("n_confirmed"), r9.get("n_qualified"), r9.get("n_refuted")] == [5, 2, 2] else "MISMATCH"})
    rows.append({"where": "R10 verdict counts 2/3", "note": [2, 0, 3],
                 "artifact": [r10.get("n_confirmed"), r10.get("n_qualified"), r10.get("n_refuted")],
                 "source": "r10_audit", "verdict": "OK" if [r10.get("n_confirmed"), r10.get("n_qualified"), r10.get("n_refuted")] == [2, 0, 3] else "MISMATCH"})
    out = {"rows": rows, "n_mismatch": sum(r["verdict"] != "OK" for r in rows),
           "r7_keys": list(r7.keys())[:30], "r8_keys": list(r8.keys())[:30]}
    merge("number_check_auto", out)
    return out


STAGES = {"topo": stage_topo, "instrument": stage_instrument, "barrier": stage_barrier,
          "boundary": stage_boundary, "taper": stage_taper, "g32": stage_g32,
          "relax_unwound": relax_unwound, "relax_g32": relax_g32,
          "codemap": stage_codemap, "equations": stage_equations, "numbers": stage_numbers}



def stage_barrier_detail():
    """1c continued: where the rigid-start paths go uphill, and the degree read
    on a surface INSIDE each melt window (the R10 audit only read hw = 8)."""
    n, L, s, g, delta = 32, 48.0, -1.0, 8.0, 0.3
    M, Q, h = ansatz(n, L, s, g, delta)
    X, Y, Z, _ = coords(n, L)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    d4 = d4_of(s, g, delta)
    Mv = np.broadcast_to(d4, M.shape)
    ns = np.linspace(0.0, 1.0, 201)
    out = {}
    for (r_in, r_out) in ((15, 21), (9, 21), (6, 12), (3, 9), (18, 21)):
        w = np.clip((r_out - R) / (r_out - r_in), 0.0, 1.0)[..., None, None]
        es = np.array([float(eu_density(M + sv * w * (Mv - M), h).sum()
                             + v4_density(M + sv * w * (Mv - M), h, s, g, delta).sum()) for sv in ns])
        de = np.diff(es)
        up = np.where(de > 0)[0]
        # the inner surface: half-width in cells inside r_in (strictly melted)
        k_in = max(1, int(np.floor(r_in / h)) - 1)
        reads = {}
        for k in sorted({k_in, 8}):
            q0 = read_degree(M, k)[0]
            q1 = read_degree(M + 1.0 * w * (Mv - M), k)[0]
            reads[f"hw{k}"] = {"r_units": k * h, "absQ_start": abs(q0), "absQ_end": abs(q1)}
        out[f"window_{r_in}_{r_out}"] = {
            "E_start": float(es[0]), "E_end": float(es[-1]), "E_min": float(es.min()),
            "s_at_min": float(ns[int(es.argmin())]),
            "n_uphill_steps": int(len(up)), "max_uphill_step": float(de.max()),
            "total_uphill": float(de[de > 0].sum()),
            "first_uphill_s": float(ns[up[0] + 1]) if len(up) else None,
            "reads": reads}
        log(f"  window ({r_in},{r_out}): min {es.min():.4f} at s={ns[es.argmin()]:.3f}, "
            f"uphill steps {len(up)} total {de[de > 0].sum():.4f}; reads {reads}")
    merge("claim_1c_detail", out)
    return out


STAGES["barrier_detail"] = stage_barrier_detail




def stage_finalize():
    """assemble the verdict block from the checkpointed sections."""
    d = json.load(open(OUT))
    c1a, c1b, c1c, c1d, c1e, c1f = (d.get(k, {}) for k in
                                    ("claim_1a", "claim_1b", "claim_1c", "claim_1d", "claim_1e", "claim_1f"))
    rr_g32 = {}
    pg = os.path.join(SCRATCH, "rerun_g32.json")
    if os.path.exists(pg):
        rr_g32 = json.load(open(pg))
    own = c1f.get("own_rerun_g32_b16_64k")
    if own and not rr_g32:
        rr_g32 = {"V4_end_own": own["V4_own"], "E_u_end_own": own["E_u_own"], "kin_end_own": own["kin_own"],
                  "source": "field saved by relax_g32 (post-save summary crashed on a since-fixed einsum)"}
    P = c1c.get("paths", {})
    rel = P.get("relaxed_deg1_to_relaxed_deg0_linear", {})
    relm = P.get("relaxed_deg1_linear_melt_15_21", {})
    rel12 = P.get("relaxed12000_deg1_to_relaxed_deg0_linear", {})
    g32_line = ("still running at finalize" if not rr_g32 else
                f"V4 {rr_g32.get('V4_end_own', float('nan')):.6f}, E_u {rr_g32.get('E_u_end_own', float('nan')):.4f}, "
                f"kin {rr_g32.get('kin_end_own', float('nan')):.3f}")
    claims = {
        "1a": {"verdict": "CONFIRMED", "deciding_number": c1a.get("deciding_number"),
               "note": ("read_charge_from_M takes eigh(...)[..., -1] of the matrix it is GIVEN; on the 4x4 M "
                        "it raises (3-vector cross product); every caller passes the 3x3 spatial block, whose "
                        "leading eigenvector is n-hat. Own reader |Q| = 1 on hw 4, 8, 12; the middle eigenvector "
                        "admits no lift (4, 71, 99 frustrated with my lift vs 37, 30, 87: lift-dependent).")},
        "1b": {"verdict": "CONFIRMED", "deciding_number": c1b.get("discrete_stabilizer_order"),
               "note": ("L d4 L^T = d4, L in SO(1,3)+, forces L to commute with d4 eta (distinct spectrum), so "
                        "L is a diagonal sign matrix; 4 survive det +1, L00 +1, all order 2: Klein. Algebra "
                        "stabilizer dim 0 (1 at delta = 0). SU(2) lift closes on Q8. pi_2 of a Lie group mod a "
                        "finite subgroup is 0. Qualifier: V4 = 0 has an indefinite second component; say 'orbit of d4'.")},
        "1c": {"verdict": "QUALIFIED", "deciding_number": c1c.get("deciding_number"),
               "note": ("NUMBER CONFIRMED (own ansatz, energy, reader; 201-point paths): no melt window rises "
                        "above the rigid-ansatz start 62.8517; endpoint 14.7940; |Q| 1 -> 0 inside each window. "
                        "CLAIM QUALIFIED: the path starts at the UNRELAXED ansatz, 49 units above the relaxed "
                        f"state; from the relaxed 3000-iteration state (E 13.800) the straight line to the relaxed "
                        f"degree-0 state rises {rel.get('barrier_max_dE', float('nan')):.3f}, the melt of the relaxed "
                        f"state rises {relm.get('barrier_max_dE', float('nan')):.3f}, from the 12000-iteration state "
                        f"{rel12.get('barrier_max_dE', float('nan')):.3f}; FIRE keeps |Q| = 1 on every surface through "
                        "12000 iterations. 'Monotone' fails in detail (4 uphill steps, 0.037, on (15,21); the "
                        "(6,12) and (3,9) windows end 1.2 and 0.9 above their minimum). 'Taper windows' are melt windows.")},
        "1d": {"verdict": "CONFIRMED", "deciding_number": c1d.get("deciding_number"),
               "note": ("own seed (linear melt 15 to 21 of my ansatz), certified FIRE 3000 its, own kin density: "
                        "272.204 = 77.5 % of 351.170, own rerun and saved endpoint identical; interior deviation "
                        "from vacuum at r <= 15 is 2.5e-3 after relaxation (exact only for the seed); |Q| = 0 on hw 4, 8, 12.")},
        "1e": {"verdict": "CONFIRMED", "deciding_number": c1e.get("deciding_number"),
               "note": ("'taper at r = 12' is a LINEAR ramp to zero over 12 to 15 (unstated): 115.385 = 32.86 % "
                        "on L = 48; 115.407 / 115.385 / 115.385 on L = 36 / 48 / 60 (L = 24: 118.567, +2.8 %). "
                        "A hard cut at 12 gives 29.0 %. L-independence is near-tautological for a compactly "
                        "supported clock on a box-independent interior.")},
        "1f": {"verdict": "CONFIRMED", "deciding_number": c1f.get("deciding_number"),
               "note": ("own V4 on the saved endpoint 0.000970, no cell with top gap < 0.35 (min 0.616), "
                        f"kin/rigid 0.9840, own E_u 48.72; own 64000-iteration rerun: {g32_line}. "
                        "But 'fmax 72.8' is the 63500-iteration log row; the endpoint fmax is 5.14.")},
    }
    number_check = d.get("number_check_auto", {}).get("rows", []) + [
        {"where": "section 1.1 'measured boost drift 32.7'", "note": 32.7, "artifact": 32.16,
         "source": "r0_audit C1 drift_eta_all_rule_under_contra_law I3 = 32.1636; no artifact carries 32.7",
         "verdict": "MISMATCH"},
        {"where": "section 6 'fmax 72.8'", "note": 72.8, "artifact": 5.141,
         "source": "relax_g32_b16_64k.trace: 72.797 at it 63500, 5.141 at it 64000 (endpoint)", "verdict": "MISMATCH"},
        {"where": "R10 row 'over five taper windows'", "note": "taper", "artifact": "melt windows (r_in, r_out)",
         "source": "r10_audit.zero_barrier_robustness.windows_r_in_r_out", "verdict": "MISMATCH (wording)"},
        {"where": "R10 row 'energy monotone 62.852 -> 14.794'", "note": "monotone",
         "artifact": "4 uphill steps (0.037 total) on (15,21); 62 steps (1.22 total) on (6,12)",
         "source": "own 201-point paths, claim_1c_detail", "verdict": "MISMATCH (wording)"},
        {"where": "R6 row 'variation <= 2e-7 on 50 dressed points'", "note": "2e-7 / 50",
         "artifact": "producer V4 max rel variation 1.86e-7 (r6_orbitblind N1); the audit measured 20 dressings at 8.5e-10",
         "source": "r6_orbitblind.json, r6_audit.json", "verdict": "OK (producer figure)"},
        {"where": "section 6 'decrements -0.185 then -0.222'", "note": "-0.185, -0.222",
         "artifact": "the producer's withdrawn numbers; the audit's pairwise decrements are -0.267/-0.403 (24,36) and -0.159/-0.130 (36,48)",
         "source": "r10_audit.R3_detail.pair_dependence", "verdict": "QUALIFIED"},
        {"where": "R0 '17/17', R3 '348 reads', '44 heals', R4 '96/96'", "note": "as written",
         "artifact": "found only in tasks/m5_32_task_details.md", "source": "task record", "verdict": "OK"},
        {"where": "five-window barrier table, clock_taper table, scope_g32 block", "note": "as written",
         "artifact": "TYPED CONSTANTS in m5_32_r10_audit.py lines 1036-1092, no producing code",
         "source": "scripts/m5_32_r10_audit.py",
         "verdict": "OK by re-derivation here (all five E_at_s1, the 115.385 taper, the g32 endpoint reproduce); not clickable"},
    ]
    cm = d.get("code_map_check_raw", {})
    code_map = [
        {"rows": cm.get("n_rows"), "anchors_off": cm.get("n_off"), "verdict": "CONFIRMED",
         "note": "all 31 anchors resolve to the named function on the working tree"},
        {"verdict": "REFUTED (send order)",
         "note": "every link is blob/main on an unmerged branch: at send time every code link is a 404; merge first or pin to the commit"},
        {"verdict": "QUALIFIED",
         "note": "the instrument row links m5_22_e_audit.read_charge_from_M, but the R10 numbers came from the re-implementation directors()/read_surface() in m5_32_r10_audit.py; the linked function cannot take the 4x4 M"},
        {"verdict": "QUALIFIED",
         "note": "m5_32_r10_audit.py, which produced the barrier, taper, degree-0 and g = 32 numbers, is in no code-map row and not in the inspection set"}]
    eq = d.get("equation_check_raw", {})
    equation_check = [
        {"item": "I1 = (1/2) F_abcd F^abcd", "verdict": "CONFIRMED", "number": eq.get("I1_note_form_vs_half_FF", {}).get("rel")},
        {"item": "V4 with C_p and w = 7.24023879e-4", "verdict": "CONFIRMED", "number": eq.get("W1")},
        {"item": "E_cert = 4 (U + omega^2 T) + V4", "verdict": "CONFIRMED", "number": eq.get("E_cert_check", {}).get("rel")},
        {"item": "quadratic Legendre H = C omega^2 - A, B drops", "verdict": "CONFIRMED", "number": eq.get("legendre_I1", {}).get("rel")},
        {"item": "quartic Legendre H = C2 omega^2 + 3 C4 omega^4 - A", "verdict": "CONFIRMED", "number": 0},
        {"item": "omega* = (J - B) / (2 C)", "verdict": "CONFIRMED", "number": 0},
        {"item": "kin = -4 C(I1) = kin_of", "verdict": "CONFIRMED", "number": eq.get("kin_equals_minus4C", {}).get("rel")},
        {"item": "a0 = Qh (G1 d4 + d4 G1^T) Qh^T vs the FD flow", "verdict": "CONFIRMED", "number": c1d.get("a0_own_vs_B8_fd_max_abs")},
        {"item": "ansatz M = Qh d4 Qh^T with R3(phi) R2(theta)", "verdict": "CONFIRMED", "number": c1c.get("calibration", {}).get("own_ansatz_vs_B8_dressed_max_abs")},
        {"item": "K_T u-frame form", "verdict": "CONFIRMED by the r7_a docstring; numerics audited at R7 A2/A3", "number": None},
        {"item": "six quartic densities present", "verdict": "CONFIRMED", "number": 6},
        {"item": "L_cert = -4 I1 - V4, E_cert = +4 U + V4", "verdict": "CONFIRMED (CERTIFIED_COEFFS I1 -4, V4 -1)", "number": None},
    ]
    overclaims = [
        {"where": "R10 row, 4.1, section 7, ledger U2: 'barrier exactly 0.0'", "verdict": "QUALIFIED",
         "note": "a property of paths from the RIGID ansatz; from the relaxed state the probes rise 0.73 / 3.04 / 4.49 and FIRE keeps |Q| = 1 through 12000 iterations; the audit's own 'metastability at most' was dropped"},
        {"where": "R10 row 'energy monotone'", "verdict": "REFUTED in detail",
         "note": "not monotone on any window at 201 points; barrier relative to the start is 0 on all five"},
        {"where": "4.1 and R10 row 'makes it L-INDEPENDENT'", "verdict": "QUALIFIED",
         "note": "tautological for a clock supported inside r = 15 on a box-independent interior; present as the bound argument"},
        {"where": "R9 row 'the clock vanishes identically at delta = 0'", "verdict": "QUALIFIED",
         "note": "true for the periodic (2,3) clock; the R9 audit (B4 REFUTED) found the radial boost K_1 smooth with a0 nonzero, kin -6.1e6"},
        {"where": "R9 row", "verdict": "UNDERSOLD",
         "note": "the audit's headline (the line resolves into a finite core, radius 3.6 to 4.0, the clock survives at 351.17 / 351.14, h-convergent) is absent"},
        {"where": "R7 row 'K_T localizes the dressing'", "verdict": "QUALIFIED",
         "note": "audit B2: the interior R* is the dressing switching off (amp* 0.0251 -> 0, worth 23 % -> 0.96 % of E_J, minimum 0.43 % deep)"},
        {"where": "R2 row 'bounded below by a pointwise theorem'", "verdict": "QUALIFIED",
         "note": "holds wherever M eta has a real timelike eigenvector; past t* = (g+1)/2 the family is undefined and R1.4 recorded I1_h unbounded below along a boost-mixing path"},
        {"where": "4.1 'eigenvalues frozen on the order-parameter space'", "verdict": "QUALIFIED",
         "note": "V4 = 0 has an indefinite second component; the pi_2 = 0 statement is about the SO(1,3)+ orbit of d4"},
        {"where": "4.1 '37, 30 and 87 frustrated edges'", "verdict": "QUALIFIED",
         "note": "lift-dependent counts (my lift 4, 71, 99); the invariant statement is 'admits no consistent lift'"},
        {"where": "section 6 'decrements NOT shrinking'", "verdict": "QUALIFIED",
         "note": "producer numbers the audit superseded; pair-dependent, neither convergence nor decay established"},
        {"where": "section 7 R8 '73 to 98 %', R0 '10/10 at <= 2.2e-15'", "verdict": "CONFIRMED", "note": "r8_audit summary; r0_baseline max_rel_diff 2.233e-15"},
    ]
    omissions = [
        "The barrier from the RELAXED degree-1 state (minimum-energy path or saddle) was never computed; only linear interpolations from the rigid ansatz. A physicist reading 'barrier exactly 0.0' will assume the relaxed object was tested.",
        "FIRE never unwound the degree: |Q| = 1 on hw 1 to 12 at 3000, 6000 and 12000 iterations (r10_audit melt_front); the note reports the collapsing top gap as 'the escape in progress' without saying the degree never moved under descent.",
        "The degree instrument reads the 3x3 spatial block; after relaxation the time row M_0i is nonzero and invisible to it.",
        "The five-window barrier, the taper ladder and the g = 32 comparison exist only as typed constants in m5_32_r10_audit.py (no producing code).",
        "The 78 % and 32.9 % rows compare 3000-iteration UNCONVERGED endpoints at equal budget; the rows do not carry that flag.",
        "No g = 32 counterpart of the degree-0 state, the taper or the barrier.",
        "The sign of the RP^2 degree is a lift convention; 4.1 says 'degree-1 reading' without 'up to sign'.",
        "The gates G1 to G7 and classes C0 to C8 are used throughout and never defined in the note.",
    ]
    reader_test = {
        "functional_in_one_click": "yes once merged (I1, V4, CERTIFIED_COEFFS rows -> m5_32_lagrangian.py, docstring carries the equations); NO at send time: unmerged branch, blob/main links 404",
        "inspection_set_physics_first": "yes; but the R10 audit script that produced the headline numbers is missing from it",
        "headline_numbers_with_gate_and_convergence": "the barrier 0.0 carries no convergence evidence and no gate of its own; the R10 gate column is the producer's prediction; 'fmax 72.8' is the wrong row",
        "vocabulary_the_author_would_not_recognize": ["Q37 as the name of the degree", "G1..G7, C0..C8 (undefined)", "Q49 to Q59", "CERTIFIED_COEFFS, INS4.kin_of, 'certified stencil', 'the record', 'rung', 'CHAN'", "'taper windows'", "'melt-front radius' (undefined: r_max of cells with top gap < 0.35)"],
        "module_size": "989 lines against the standard's ~200-line target; the docstring carries the equations",
        "section_1_5_instrument_description": "'eigh(M), take V[..., -1]' is wrong for the 4x4 M; must say the spatial 3x3 block",
    }
    blocking = [
        "B1 Re-scope the lead claim 'unwinding barrier exactly 0.0' to what was measured: a straight-line melt from the UNRELAXED ansatz (E 62.85, 49 units above the relaxed state) never rises above its start. From the relaxed state the same probes rise 0.73 to 4.49 and FIRE holds |Q| = 1 through 12000 iterations. Say: no barrier protects the ansatz reading; the relaxed state's barrier was not computed; pi_2 = 0 is the protection argument. Same fix in section 7, ledger U2 and the outbound message.",
        "B2 Sections 1.5 and 4.1 describe the instrument as 'eigh(M), take V[..., -1]'. On the 4x4 M that is the constant time axis (degree 0) and the linked function raises; every caller passes the 3x3 SPATIAL block. Write 'eigh of the spatial 3x3 block'.",
        "B3 'energy monotone 62.852 -> 14.794' is false at 201 points (4 uphill steps on that window, 62 on (6,12)); write 'never above the start'. 'Five taper windows' are melt windows.",
        "B4 All code links are blob/main on an unmerged branch: every link is dead at the reader's first click. Merge before sending, or pin the links to the commit.",
        "B5 The gates G1 to G7 and classes C0 to C8 are never defined in the note; add a one-line legend or a link to the task record section 8 table.",
        "B6 'fmax 72.8' (section 6) is the 63500-iteration row, the endpoint is 5.14; 'boost drift 32.7' (section 1.1) matches no artifact, the R0 audit carries 32.16.",
    ]
    minor = [
        "M1 R9 row: 'the clock vanishes identically' -> 'the periodic clock' (the audit found a nonzero smooth radial-boost flow, kin -6.1e6), and add the audit's headline (finite core, clock survives at 351.17 / 351.14).",
        "M2 R7 row: 'K_T localizes the dressing' is the audit's QUALIFIED B2 (the dressing switches off; 0.43 % deep minimum).",
        "M3 R2 row: add 'wherever M eta has a real timelike eigenvector'.",
        "M4 4.1: 'eigenvalues frozen on the order-parameter space' -> 'on the SO(1,3)+ orbit of d4'.",
        "M5 4.1: frustrated-edge counts are lift-dependent (4, 71, 99 with another lift); say 'admits no consistent lift'.",
        "M6 'makes it L-INDEPENDENT' is tautological for a compactly supported clock; state that 'taper at r = 12' is a linear ramp to zero over 12 to 15.",
        "M7 Section 6 'decrements -0.185 then -0.222' are the producer's withdrawn numbers; cite the audit's pair-dependent decrements.",
        "M8 Section 6 should add: the relaxed-state barrier, FIRE never moving the degree, the 3x3-only instrument, the equal-budget unconverged comparisons.",
        "M9 The five-window barrier table, the taper table and the g = 32 block are typed constants in m5_32_r10_audit.py with no producing code; add the script to the code map and the inspection set.",
        "M10 The instrument row links m5_22_e_audit.read_charge_from_M while the R10 numbers came from m5_32_r10_audit.directors/read_surface.",
        "M11 Vocabulary: Q37, Q49 to Q59, CHAN, 'the record', INS4.kin_of, 'melt-front radius' (define it).",
        "M12 The registry module is 989 lines against the standard's ~200-line target.",
        "M13 The R10 gate column holds the producer's prediction, not a gate for the unwinding or boundary claims; those carry no convergence evidence.",
        "M14 The RP^2 degree sign is a lift convention; write 'degree +-1'.",
    ]
    summary = ("Five of six load-bearing claims re-derive with independent code (1a, 1b, 1d, 1e, 1f). The sixth, "
               "the unwinding barrier, reproduces as a NUMBER (0.0 relative to the rigid-ansatz start on all five "
               "melt windows, endpoint 14.7940) but not as the CLAIM the note and the outbound message lead with: "
               "the path starts 49 energy units above the relaxed state, from the relaxed state the same probes "
               "rise 0.73 to 4.49, FIRE never moved the degree in 12000 iterations, and the energy is not monotone. "
               "Equations, 31 code anchors and the transcribed numbers hold except 'fmax 72.8' (mid-run row; "
               "endpoint 5.14) and 'boost drift 32.7' (artifact 32.16). Two legibility failures at the first click: "
               "blob/main links on an unmerged branch, and gates G1-G7 used undefined. Six blocking, fourteen minor.")
    out = {"claims": claims, "number_check": number_check, "code_map_check": code_map,
           "equation_check": equation_check, "overclaims": overclaims, "omissions": omissions,
           "reader_test": reader_test, "blocking": blocking, "minor": minor, "summary": summary,
           "n_blocking": len(blocking), "n_minor": len(minor), "own_rerun_g32": rr_g32,
           "runtime_s": d.get("runtime_s")}
    for k, v in out.items():
        d[k] = v
    with open(OUT, "w") as f:
        json.dump(d, f, indent=1, default=float)
    log("finalized")
    return out


STAGES["finalize"] = stage_finalize


if __name__ == "__main__":
    for st in sys.argv[1:]:
        log(f"=== stage {st}")
        STAGES[st]()
    log("done")
