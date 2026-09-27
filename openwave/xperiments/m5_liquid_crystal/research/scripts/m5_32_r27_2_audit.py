"""M5.32 R27-2 (the census under kappa) and R27-3 (e) (the delta arm): the adversarial audit.

An independent second pass over the claims of `data/m5_32_r27_2_kappa_census.json` and
`data/m5_32_r27_3_arm.json`, written WITHOUT reading the audited scripts
(`m5_32_r27_2_kappa_census.py`, `m5_32_r27_0_form.py`, `m5_32_r27_3_reads.py`). The shared
stack is consumed for the certified quartic energy (`R20.energy_parts`), the seed factory
(`R20.seed_axes`), the R26-0 readers (`partition_reader`, `degree_reader`,
`physical_generator`) and the R26-3 charge density; every number that carries a verdict is
ALSO recomputed here with own stencils (curvature, V4, the smoothness term E_kappa and its
gradient, the biaxiality formula, the pair-gap closure clusters, the solid-angle degree, the
orbital generator with a central-difference stencil).

Own definitions (the brief's):
  E_quartic = 4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta + V4,  F_ij = A_i eta A_j - A_j eta A_i
  E_kappa   = h^3 sum_br (1/2) sum_i sum_cells |d_i S|_F^2,  S = M[..., 1:, 1:]  (= h sum_bonds |dS|_F^2)
  E_total   = E_quartic + kappa E_kappa
  beta^2    = 1 - 6 (tr L^3)^2 / (tr L^2)^3 on the traceless part L of the spatial block
  gate      = fmax_spatial < 1e-4 AND the last chunk's drop < 1e-5 |E|
  certified = AT_GATE and (STABLE, or SADDLE and reconverged)

Outputs: data/m5_32_r27_2_audit.json (per-claim verdicts, own vs claimed numbers, summary).
Run: OMP_NUM_THREADS=1 python3 scripts/m5_32_r27_2_audit.py   (one process, about 15 min)
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
from scipy.ndimage import label as nd_label  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r27_2_audit.json")
CENSUS_JSON = os.path.join(DATA, "m5_32_r27_2_kappa_census.json")
CENSUS_DIR = os.path.join(DATA, "m5_32_r27_2")
ARM_JSON = os.path.join(DATA, "m5_32_r27_3_arm.json")
ARM_DIR = os.path.join(DATA, "m5_32_r27_3")
R26_4_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
R26_4_DIR = os.path.join(DATA, "m5_32_r26_4")
S1_D03_END = os.path.join(DATA, "m5_32_r25_2", "S1_d0.3_w25_n32_L48.npz")
T0 = time.time()

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G = 8.0
W1S = 25.0
N, L = 32, 48.0
PARTS = ["1_1_1_1", "2_1_1", "2_2", "3_1", "4"]
SUM_K2 = {"1_1_1_1": 1.0, "2_1_1": 1.5, "2_2": 2.0, "3_1": 2.5, "4": 4.0}
KAPPAS = {"a": 2.299e-03, "b": 9.197e-03}
SPHERES = (3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0)
GATE_FMAX = 1e-4
GATE_DROP_REL = 1e-5
OWN_DESCENT_ITERS = 300


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


F0 = _load("m5_32_r26_0_form", "m5_32_r26_0_form.py")
R21, R20, B3, R0, W1 = F0.R21, F0.R20, F0.B3, F0.R0, F0.W1
R263 = _load("m5_32_r26_3_spin", "m5_32_r26_3_spin.py")


def cfg_pot(delta):
    cfg = R21.cfg_of(N, L, G, delta)
    p = R21.params_of(G, delta)
    pot = ("v4", R0.roots_of(cfg), W1 * W1S)
    return cfg, p, pot


def tag_of(part, kappa):
    return f"P{part}_d0.3_w25_n32_L48_k{kappa:.3e}"


def load_M(path):
    return np.load(path)["M"].astype(np.float64)


# ============================================================================
# OWN STENCILS AND ENERGIES
# ============================================================================
def own_diff(f, ax, h):
    """forward bond differences along ax (one fewer entry along ax)."""
    n = f.shape[ax]
    return (np.take(f, range(1, n), axis=ax) - np.take(f, range(0, n - 1), axis=ax)) / h


def own_branch(M, h, br):
    """A_i on every cell: the bond difference stored at its lower (fwd) or upper (bwd) cell."""
    A = []
    for ax in range(3):
        out = np.zeros_like(M)
        sl = [slice(None)] * 3
        sl[ax] = slice(0, -1) if br == "fwd" else slice(1, None)
        out[tuple(sl)] = own_diff(M, ax, h)
        A.append(out)
    return A


def inner_eta(F, Gm):
    """<F, G>_eta = tr(eta F eta G^T) per cell."""
    return np.sum((ETA @ F @ ETA) * Gm, axis=(-2, -1))


def own_curv_density(M, cfg):
    h = cfg["h"]
    e = np.zeros(M.shape[:3])
    for br in ("fwd", "bwd"):
        A = own_branch(M, h, br)
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                e += 0.5 * 4.0 * inner_eta(F, F)
    return h**3 * e


def own_v4_density(M, cfg, roots, w):
    Nm = M @ ETA
    P = Nm.copy()
    tr = []
    for k in range(1, 5):
        if k > 1:
            P = P @ Nm
        tr.append(np.trace(P, axis1=-2, axis2=-1))
    C = [sum(q**k for q in roots) for k in range(1, 5)]
    return cfg["h"] ** 3 * w * sum((t - c) ** 2 for t, c in zip(tr, C))


def own_ekappa_density(M, cfg):
    """h^3 (1/2)(fwd + bwd) sum_i |d_i S|_F^2 per cell; its sum is h sum_bonds |dS|_F^2."""
    S = M[..., 1:, 1:]
    h = cfg["h"]
    e = np.zeros(S.shape[:3])
    for ax in range(3):
        q = np.sum(own_diff(S, ax, h) ** 2, axis=(-2, -1))
        sl = [slice(None)] * 3
        sl[ax] = slice(0, -1)
        e[tuple(sl)] += 0.5 * q
        sl[ax] = slice(1, None)
        e[tuple(sl)] += 0.5 * q
    return h**3 * e


def own_ekappa_bonds(M, cfg):
    """the same term as a plain bond sum, h sum_bonds |S_b - S_a|_F^2 (the cross-check)."""
    S = M[..., 1:, 1:]
    h = cfg["h"]
    tot = 0.0
    for ax in range(3):
        d = np.take(S, range(1, S.shape[ax]), axis=ax) - np.take(
            S, range(0, S.shape[ax] - 1), axis=ax
        )
        tot += float(np.sum(d * d))
    return h * tot


def own_ekappa_grad(M, cfg):
    """dE_kappa/dM on the spatial block: 2 h sum_nb (S_c - S_nb)."""
    S = M[..., 1:, 1:]
    h = cfg["h"]
    g = np.zeros_like(S)
    for ax in range(3):
        d = np.take(S, range(1, S.shape[ax]), axis=ax) - np.take(
            S, range(0, S.shape[ax] - 1), axis=ax
        )
        sl = [slice(None)] * 3
        sl[ax] = slice(0, -1)
        g[tuple(sl)] -= d
        sl[ax] = slice(1, None)
        g[tuple(sl)] += d
    out = np.zeros_like(M)
    out[..., 1:, 1:] = 2.0 * h * g
    return out


def own_split(M, cfg, pot, kappa):
    ec = float(own_curv_density(M, cfg).sum())
    v4 = float(own_v4_density(M, cfg, pot[1], pot[2]).sum())
    ek = float(own_ekappa_density(M, cfg).sum())
    return {
        "E_curv": ec,
        "V4": v4,
        "E_quartic": ec + v4,
        "E_kappa": ek,
        "E_kappa_bondsum": own_ekappa_bonds(M, cfg),
        "kappaE_kappa": kappa * ek,
        "E_total": ec + v4 + kappa * ek,
    }


def radius(cfg):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    return X, Y, Z, np.sqrt(X * X + Y * Y + Z * Z)


def own_beta2(lam):
    lam = np.asarray(lam, float)
    Lm = lam - lam.mean()
    t2 = np.sum(Lm**2)
    t3 = np.sum(Lm**3)
    return float(1.0 - 6.0 * t3**2 / max(t2**3, 1e-300))


def own_beta2_block(M, cfg, R):
    _, _, _, r = radius(cfg)
    Sbar = M[r < R][:, 1:, 1:].mean(axis=0)
    return own_beta2(np.linalg.eigvalsh(Sbar))


def own_fmax(M, cfg, p, pot, kappa, free):
    """the largest gradient entry of the FULL energy (quartic + kappa) over the free cells,
    spatial block and time row separately (the packing convention of the audited descent is
    unknown, so this is a ballpark cross-check of the gate, not a reproduction)."""
    E, Gq, info = R20.energy_grad(M, cfg, p, pot)
    Gt = Gq + kappa * own_ekappa_grad(M, cfg)
    gs = np.abs(Gt[free][:, 1:, 1:]).max()
    g00 = np.abs(Gt[free][:, 0, 0]).max()
    return float(gs), float(g00), float(E + kappa * own_ekappa_density(M, cfg).sum())


# ============================================================================
# THE GATE AND KICK LOGIC READ OFF THE CHUNKS
# ============================================================================
def gate_from_chunks(chunks):
    """the first chunk that satisfies the gate, and whether the LAST chunk does."""
    first = None
    for c in chunks[1:]:
        if c.get("fmax_spatial", 1.0) < GATE_FMAX and c.get("drop", 1.0) < GATE_DROP_REL * abs(
            c["E"]
        ):
            first = c["iters"]
            break
    last = chunks[-1]
    last_ok = last.get("fmax_spatial", 1.0) < GATE_FMAX and last.get(
        "drop", 1.0
    ) < GATE_DROP_REL * abs(last["E"])
    return {
        "first_gate_iters": first,
        "last_iters": last["iters"],
        "last_fmax": last.get("fmax_spatial"),
        "last_drop": last.get("drop"),
        "last_drop_rel": (
            (last.get("drop", 0.0) / abs(last["E"])) if last.get("drop") is not None else None
        ),
        "last_at_gate": bool(last_ok),
    }


def kick_from_row(kick, E_ref):
    ch = kick["chunks"]
    last = ch[-1]
    E_after = last["E"]
    tol = GATE_DROP_REL * abs(E_ref)
    lower = E_after < E_ref - tol
    at_gate = last.get("fmax_spatial", 1.0) < GATE_FMAX and last.get("drop", 1.0) < tol
    return {
        "E_kicked_start": ch[0]["E"],
        "E_after": E_after,
        "E_after_minus_ref": E_after - E_ref,
        "fmax_after": last.get("fmax_spatial"),
        "last_drop": last.get("drop"),
        "own_label": "SADDLE" if lower else "STABLE",
        "own_reconverged": bool(at_gate),
        "kick_iters": last["iters"],
    }


# ============================================================================
# OWN READS
# ============================================================================
def pair_gap_clusters(M, cfg, delta, R, frac=0.1):
    """cells on the one-cell shell |r - R| <= h/2 whose pair gap is below frac delta, clustered
    by 26-connectivity: count, sizes, centroids (x, y, z) and (theta, phi)."""
    X, Y, Z, r = radius(cfg)
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    gap = lam[..., 1] - lam[..., 0]
    shell = np.abs(r - R) <= 0.5 * cfg["h"] * np.sqrt(3.0)
    closed = shell & (gap < frac * delta)
    lab, k = nd_label(closed, structure=np.ones((3, 3, 3)))
    cl = []
    for i in range(1, k + 1):
        m = lab == i
        cx, cy, cz = float(X[m].mean()), float(Y[m].mean()), float(Z[m].mean())
        rr = max(np.sqrt(cx * cx + cy * cy + cz * cz), 1e-12)
        cl.append(
            {
                "cells": int(m.sum()),
                "xyz": [round(cx, 2), round(cy, 2), round(cz, 2)],
                "theta_deg": round(float(np.degrees(np.arccos(cz / rr))), 1),
                "phi_deg": round(float(np.degrees(np.arctan2(cy, cx))), 1),
                "min_gap": float(gap[m].min()),
            }
        )
    cl.sort(key=lambda c: -c["cells"])
    return {
        "R": R,
        "shell_cells": int(shell.sum()),
        "closed_cells": int(closed.sum()),
        "clusters": k,
        "list": cl,
    }


def own_degree(M, cfg, R, nth=48, nph=96):
    """the outward-oriented director's degree on the sphere R by the solid-angle sum of the
    lat-long triangles (Van Oosterom), the director sampled trilinearly from S. Own code."""
    n, h = cfg["n"], cfg["h"]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    pts = R * rhat
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    S = np.zeros(TH.shape + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], coords, order=1, mode="nearest").reshape(
                TH.shape
            )
            S[..., a, b] = v
            S[..., b, a] = v
    lam, V = np.linalg.eigh(S)
    d = V[..., :, 2]
    al = np.einsum("...a,...a->...", d, rhat)
    d = d * np.where(al < 0, -1.0, 1.0)[..., None]

    def tri(a, b, c):
        num = np.einsum("...a,...a->...", a, np.cross(b, c))
        den = (
            1.0
            + np.einsum("...a,...a->...", a, b)
            + np.einsum("...a,...a->...", b, c)
            + np.einsum("...a,...a->...", c, a)
        )
        return 2.0 * np.arctan2(num, den)

    # quads (i, j), (i, j+1), (i+1, j+1), (i+1, j), split in two triangles, oriented so that the
    # identity map d = rhat gives +1
    a = d[:-1, :]
    b = np.roll(d, -1, axis=1)[:-1, :]
    c = np.roll(d, -1, axis=1)[1:, :]
    e = d[1:, :]
    area = np.sum(tri(a, c, b)) + np.sum(tri(a, e, c))
    north = d[0].mean(axis=0)
    north /= max(np.linalg.norm(north), 1e-300)
    south = d[-1].mean(axis=0)
    south /= max(np.linalg.norm(south), 1e-300)
    r0, r0n = d[0], np.roll(d[0], -1, axis=0)
    area += np.sum(tri(np.broadcast_to(north, r0.shape), r0n, r0))
    r1, r1n = d[-1], np.roll(d[-1], -1, axis=0)
    area += np.sum(tri(np.broadcast_to(south, r1.shape), r1, r1n))
    return {
        "R": R,
        "degree": float(area / (4.0 * np.pi)),
        "min_abs_align": float(np.abs(al).min()),
        "gap_dir_min": float((lam[..., 2] - lam[..., 1]).min()),
    }


def own_generator(M, cfg):
    """C_int, C_orb, C_rigid with an own central-difference orbital stencil (np.gradient:
    central inside, one-sided at the box faces) and own kin. The R26-0 stencil (dsym) halves
    the one-sided face difference; the two differ only on the faces."""
    h = cfg["h"]
    JZ = np.zeros((4, 4))
    JZ[1, 2], JZ[2, 1] = -1.0, 1.0
    X, Y, Z, r = radius(cfg)
    a_int = JZ @ M - M @ JZ
    dMx = np.gradient(M, h, axis=0)
    dMy = np.gradient(M, h, axis=1)
    a_orb = X[..., None, None] * dMy - Y[..., None, None] * dMx
    a_rig = a_int - a_orb

    def kin_density(a):
        k = np.zeros(M.shape[:3])
        for br in ("fwd", "bwd"):
            A = own_branch(M, h, br)
            for i in range(3):
                F = a @ ETA @ A[i] - A[i] @ ETA @ a
                k += 0.5 * 4.0 * inner_eta(F, F)
        return h**3 * k

    kd = {
        "internal": kin_density(a_int),
        "orbital": kin_density(a_orb),
        "rigid": kin_density(a_rig),
    }
    out = {k: float(v.sum()) for k, v in kd.items()}
    out["within_r18"] = {k: float(v[r < 18.0].sum()) for k, v in kd.items()}
    out["rigid_over_internal"] = out["rigid"] / max(out["internal"], 1e-300)
    return out


def ekappa_channels(M, cfg):
    """the lattice-exact channel split of E_kappa: per bond (a, b), |S_b - S_a|_F^2 =
    sum_{k,l} (v_bk . v_al)^2 (l_bk - l_al)^2 (eigenframes V, eigenvalues l ascending: 0, 1 the
    pair, 2 the director). Channels: eigenvalue (k = l), director rotation (k or l = 2, k != l),
    pair rotation (k, l in {0, 1}, k != l: the winding piece the author's log term describes).
    The channel sums add to E_kappa exactly (h sum_bonds)."""
    S = M[..., 1:, 1:]
    h = cfg["h"]
    lam, V = np.linalg.eigh(S)
    ch = {"eigenvalue": 0.0, "director": 0.0, "pair": 0.0}
    pair_core = np.zeros(S.shape[:3])
    for ax in range(3):
        n = S.shape[ax]
        Va = np.take(V, range(0, n - 1), axis=ax)
        Vb = np.take(V, range(1, n), axis=ax)
        la = np.take(lam, range(0, n - 1), axis=ax)
        lb = np.take(lam, range(1, n), axis=ax)
        Q2 = (np.swapaxes(Vb, -1, -2) @ Va) ** 2  # Q2[k, l] = (v_bk . v_al)^2
        D2 = (lb[..., :, None] - la[..., None, :]) ** 2
        C = Q2 * D2
        ch["eigenvalue"] += float(np.einsum("...kk->...", C).sum())
        ch["pair"] += float((C[..., 0, 1] + C[..., 1, 0]).sum())
        ch["director"] += float((C[..., 2, 0] + C[..., 0, 2] + C[..., 2, 1] + C[..., 1, 2]).sum())
        sl = [slice(None)] * 3
        sl[ax] = slice(0, -1)
        pair_core[tuple(sl)] += 0.5 * (C[..., 0, 1] + C[..., 1, 0])
        sl[ax] = slice(1, None)
        pair_core[tuple(sl)] += 0.5 * (C[..., 0, 1] + C[..., 1, 0])
    out = {k: h * v for k, v in ch.items()}
    out["sum"] = out["eigenvalue"] + out["director"] + out["pair"]
    out["E_kappa_check"] = own_ekappa_bonds(M, cfg)
    return out, h * pair_core


def ekappa_regions(M, cfg, delta, r_core=4.5):
    """E_kappa split by region: the director core r < r_core, the pair cores outside it (pair
    gap below delta / 2), the pin shell, the rest."""
    X, Y, Z, r = radius(cfg)
    e = own_ekappa_density(M, cfg)
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    gap = lam[..., 1] - lam[..., 0]
    pin = B3.pin_shell(cfg["n"], cfg["h"])
    core = r < r_core
    pcore = (~core) & (~pin) & (gap < 0.5 * delta)
    rest = (~core) & (~pin) & (~pcore)
    return {
        "total": float(e.sum()),
        "director_core_r<4.5": float(e[core].sum()),
        "pair_cores_gap<delta/2": float(e[pcore].sum()),
        "pair_core_cells": int(pcore.sum()),
        "pin_shell": float(e[pin].sum()),
        "rest": float(e[rest].sum()),
    }


# ============================================================================
# OWN DESCENT (the {2,2} kappa-b saddle, continued from the audited end field)
# ============================================================================
IU = np.triu_indices(4)
OFF = np.where(IU[0] != IU[1], 2.0, 1.0)


def own_descent(M0, cfg, p, pot, kappa, free, iters, tag):
    base = M0.copy()
    x0 = M0[free][:, IU[0], IU[1]].ravel().copy()
    hist = []

    def unpack(x):
        S = x.reshape(-1, 10)
        F = np.zeros((S.shape[0], 4, 4))
        F[:, IU[0], IU[1]] = S
        F[:, IU[1], IU[0]] = S
        Mx = base.copy()
        Mx[free] = F
        return Mx

    def fg(x):
        Mx = unpack(x)
        E, Gq, info = R20.energy_grad(Mx, cfg, p, pot)
        if not np.isfinite(E):
            return 1e30, np.zeros_like(x)
        Et = E + kappa * float(own_ekappa_density(Mx, cfg).sum())
        Gt = Gq + kappa * own_ekappa_grad(Mx, cfg)
        g = (Gt[free][:, IU[0], IU[1]] * OFF).ravel()
        hist.append(Et)
        return Et, g

    res = minimize(
        fg,
        x0,
        jac=True,
        method="L-BFGS-B",
        options={
            "maxiter": iters,
            "maxfun": int(1.6 * iters),
            "maxcor": 20,
            "ftol": 0.0,
            "gtol": 0.0,
        },
    )
    Mf = unpack(res.x)
    gs, g00, Ef = own_fmax(Mf, cfg, p, pot, kappa, free)
    log(
        f"own descent {tag}: E {hist[0]:.6f} -> {res.fun:.6f} in {res.nit} it / {res.nfev} ev, gmax {gs:.2e}, {res.message}"
    )
    return Mf, {
        "E_start": float(hist[0]),
        "E_end": float(res.fun),
        "drop": float(hist[0] - res.fun),
        "nit": int(res.nit),
        "nfev": int(res.nfev),
        "gmax_spatial_end": gs,
        "message": str(res.message),
        "E_every_50_evals": [float(v) for v in hist[::50]],
    }


# ============================================================================
# THE AUDIT
# ============================================================================
def main():
    J = json.load(open(CENSUS_JSON))
    rows = J["rows"]
    J4 = json.load(open(R26_4_JSON))["rows"]
    JA = json.load(open(ARM_JSON))
    cfg, p, pot = cfg_pot(0.3)
    delta = 0.3
    free = R21.free_mask(cfg, True)
    pin = B3.pin_shell(cfg["n"], cfg["h"])
    out = {
        "task": "M5.32 R27-2 + R27-3(e) audit",
        "date": "2026-09-27",
        "claims": {},
        "verdicts": {},
    }

    # ---------------- claim 1: gate, kick, energy splits ----------------
    log("claim 1: energies and gate/kick logic")
    c1 = {}
    fields = {}
    max_dev = 0.0
    for kk, kappa in KAPPAS.items():
        for part in PARTS:
            tag = tag_of(part, kappa)
            r = rows[tag]
            M = load_M(os.path.join(CENSUS_DIR, tag + ".npz"))
            fields[(part, kk)] = M
            own = own_split(M, cfg, pot, r["kappa"])
            stack = R20.energy_parts(M, cfg, p, pot)
            dev = {
                k: own[k] - r["end_split"][k]
                for k in ("E_quartic", "E_curv", "V4", "E_kappa", "kappaE_kappa", "E_total")
            }
            max_dev = max(max_dev, max(abs(v) for v in dev.values()))
            gate = gate_from_chunks(r["chunks"])
            kick = kick_from_row(r["kick"], r["kick"]["E_ref_reslaved"])
            gs, g00, Et = own_fmax(M, cfg, p, pot, r["kappa"], free)
            Mseed = load_M(os.path.join(CENSUS_DIR, tag + "_seed.npz"))
            Mgate = load_M(os.path.join(CENSUS_DIR, tag + "_gate.npz"))
            c1[tag] = {
                "own_split": own,
                "claimed_end_split": r["end_split"],
                "max_abs_dev_vs_end_split": max(abs(v) for v in dev.values()),
                "E_quartic_stack_energy_parts": float(stack["E_total"]),
                "E_kappa_bondsum_minus_branch": own["E_kappa_bondsum"] - own["E_kappa"],
                "gate_from_chunks": gate,
                "claimed_gate_label": r["gate_label"],
                "claimed_iters": r["iters"],
                "kick_from_chunks": kick,
                "claimed_kick_label": r["kick_label"],
                "claimed_reconverged": r["reconverged"],
                "claimed_E_err": r["E_err"],
                "own_gmax_spatial_end_field": gs,
                "own_gmax_M00_end_field": g00,
                "pin_shell_equals_seed": bool(np.array_equal(M[pin], Mseed[pin])),
                "end_equals_gate_field": bool(np.array_equal(M, Mgate)),
                "seed_E_total_own": own_split(Mseed, cfg, pot, r["kappa"])["E_total"],
                "claimed_E_seed": r["E_seed"],
                "own_certified": bool(
                    gate["last_at_gate"]
                    and (kick["own_label"] == "STABLE" or kick["own_reconverged"])
                ),
            }
            log(
                f"  {tag}: dev {c1[tag]['max_abs_dev_vs_end_split']:.1e}  gate(first {gate['first_gate_iters']}, last_ok {gate['last_at_gate']})"
                f"  kick own {kick['own_label']} reconv {kick['own_reconverged']} (claimed {r['kick_label']}/{r['reconverged']})  gmax {gs:.1e}"
            )
    out["claims"]["1"] = c1
    iters_claimed = sorted(rows[t]["iters"] for t in c1)
    all_gate = all(c1[t]["gate_from_chunks"]["last_at_gate"] for t in c1)
    labels_ok = all(
        c1[t]["kick_from_chunks"]["own_label"] == rows[t]["kick_label"]
        and c1[t]["kick_from_chunks"]["own_reconverged"] == rows[t]["reconverged"]
        for t in c1
    )
    out["verdicts"]["1"] = {
        "verdict": "CONFIRMED" if (max_dev < 1e-10 and all_gate and labels_ok) else "QUALIFIED",
        "max_abs_dev_end_split": max_dev,
        "all_last_chunks_at_gate": all_gate,
        "kick_labels_reproduced": labels_ok,
        "iters_range_own": [iters_claimed[0], iters_claimed[-1]],
        "iters_range_claimed": [750, 1750],
        "kappa0_R26_4_iters_and_labels": {
            q: [J4[f"P{q}_d0.3_w25_n32_L48"]["iters"], J4[f"P{q}_d0.3_w25_n32_L48"]["gate_label"]]
            for q in PARTS
        },
        "own_gmax_equals_chunk_fmax_on_stable_rows": all(
            abs(c1[t]["own_gmax_spatial_end_field"] - c1[t]["gate_from_chunks"]["last_fmax"])
            < 1e-12
            for t in c1
            if rows[t]["kick_label"] == "STABLE"
        ),
        "note": "",
    }
    if iters_claimed[0] < 750:
        out["verdicts"]["1"]["verdict"] = "QUALIFIED"
        out["verdicts"]["1"]["note"] = (
            f"the stated 750 to 1750 range is off: {tag_of('2_2', KAPPAS['b'])} reached the gate at 445 iterations "
            "(L-BFGS ABNORMAL termination, fmax 9e-9); every other statement of the claim reproduces"
        )

    # ---------------- claim 2: orders, inversions, E_kappa anatomy ----------------
    log("claim 2: orders and E_kappa anatomy")
    c2 = {"per_kappa": {}}
    for kk, kappa in KAPPAS.items():
        E = {part: c1[tag_of(part, kappa)]["own_split"]["E_total"] for part in PARTS}
        Eq = {part: c1[tag_of(part, kappa)]["own_split"]["E_quartic"] for part in PARTS}
        Ek = {part: c1[tag_of(part, kappa)]["own_split"]["E_kappa"] for part in PARTS}
        err = {part: rows[tag_of(part, kappa)]["E_err"] for part in PARTS}
        order = sorted(PARTS, key=lambda q: E[q])
        order_q = sorted(PARTS, key=lambda q: Eq[q])
        order_k2 = sorted(PARTS, key=lambda q: SUM_K2[q])
        inversions = []
        for i, a in enumerate(order_k2):
            for b in order_k2[i + 1 :]:
                if E[a] > E[b]:
                    inversions.append(
                        {
                            "pair": [a, b],
                            "E_a_minus_E_b": E[a] - E[b],
                            "two_errors": err[a] + err[b],
                            "resolved": E[a] - E[b] > err[a] + err[b],
                        }
                    )
        c2["per_kappa"][kk] = {
            "kappa": kappa,
            "E_total_own": E,
            "E_quartic_own": Eq,
            "E_kappa_own": Ek,
            "E_err_claimed": err,
            "order_total": order,
            "order_quartic": order_q,
            "order_by_sum_k2": order_k2,
            "inversions_vs_sum_k2": inversions,
            "spread_total": max(E.values()) - min(E.values()),
            "spread_kappaE_kappa": kappa * (max(Ek.values()) - min(Ek.values())),
            "d_1111_vs_211": E["1_1_1_1"] - E["2_1_1"],
            "d_31_vs_1111": E["3_1"] - E["1_1_1_1"],
            "two_err_31_1111": err["3_1"] + err["1_1_1_1"],
        }
        log(f"  kappa {kk}: order {order}, quartic order {order_q}, inversions {len(inversions)}")
    # the E_kappa anatomy at kappa a: the hedgehog seeds, the D + P split, the regions
    Mhh = R20.seed_axes(cfg, (1.0, delta, 0.0))
    Mdir = R20.seed_axes(cfg, (1.0, 0.5 * delta, 0.5 * delta))
    anat = {
        "E_kappa_hedgehog_seed_axes_1_0.3_0": float(own_ekappa_density(Mhh, cfg).sum()),
        "E_kappa_director_only_seed_axes_1_0.15_0.15": float(own_ekappa_density(Mdir, cfg).sum()),
        "channels_hedgehog_seed_1_0.3_0": ekappa_channels(Mhh, cfg)[0],
        "channels_director_only_seed": ekappa_channels(Mdir, cfg)[0],
        "per_part": {},
    }
    xi = 1.5
    b0 = 0.5 * delta
    Rmax = 0.5 * L - 1.6
    anat["log_term_per_unit_k2_16pi_b0^2_ln(L/xi)"] = {
        "b0": b0,
        "xi": xi,
        "ln(Rbox/xi)": float(np.log(Rmax / xi)),
        "per_unit_k2": float(16 * np.pi * b0**2 * np.log(Rmax / xi)),
        "predicted_spread_1111_to_4": float(3.0 * 16 * np.pi * b0**2 * np.log(Rmax / xi)),
    }
    for part in PARTS:
        tag = tag_of(part, KAPPAS["a"])
        M = fields[(part, "a")]
        Mseed = load_M(os.path.join(CENSUS_DIR, tag + "_seed.npz"))
        che, pc = ekappa_channels(M, cfg)
        chs, _ = ekappa_channels(Mseed, cfg)
        anat["per_part"][part] = {
            "E_kappa_end": c1[tag]["own_split"]["E_kappa"],
            "E_kappa_seed": float(own_ekappa_density(Mseed, cfg).sum()),
            "channels_end": che,
            "channels_seed": chs,
            "regions_end": ekappa_regions(M, cfg, delta),
            "regions_seed": ekappa_regions(Mseed, cfg, delta),
        }
        log(
            f"  {part}: E_kappa end {che['sum']:.1f} = eig {che['eigenvalue']:.1f} + dir {che['director']:.1f} + pair {che['pair']:.1f}; seed {chs['sum']:.1f} (pair {chs['pair']:.1f})"
        )
    c2["E_kappa_anatomy_kappa_a"] = anat
    out["claims"]["2"] = c2
    ka, kb = c2["per_kappa"]["a"], c2["per_kappa"]["b"]
    claimed_a = ["2_1_1", "3_1", "1_1_1_1", "2_2", "4"]
    claimed_b = ["2_1_1", "3_1", "2_2", "1_1_1_1", "4"]
    claimed_qa = ["1_1_1_1", "2_1_1", "3_1", "2_2", "4"]
    inv_a = [x for x in ka["inversions_vs_sum_k2"] if x["pair"] == ["1_1_1_1", "2_1_1"]]
    inv_b = [x for x in kb["inversions_vs_sum_k2"] if x["pair"] == ["1_1_1_1", "2_1_1"]]
    pair_ok = bool(inv_a and inv_a[0]["resolved"] and inv_b and inv_b[0]["resolved"])
    Ppart = {q: anat["per_part"][q]["channels_end"]["pair"] for q in PARTS}
    Dpart = {q: anat["per_part"][q]["channels_end"]["director"] for q in PARTS}
    Vpart = {q: anat["per_part"][q]["channels_end"]["eigenvalue"] for q in PARTS}
    out["verdicts"]["2"] = {
        "verdict": "QUALIFIED",
        "orders_reproduced": {
            "total_a": ka["order_total"] == claimed_a,
            "total_b": kb["order_total"] == claimed_b,
            "quartic_a": ka["order_quartic"] == claimed_qa,
            "quartic_b_first_four": kb["order_quartic"][:3] == ["2_1_1", "3_1", "1_1_1_1"]
            and kb["order_quartic"][-1] == "4",
        },
        "inversion_1111_vs_211": {
            "a": inv_a[0] if inv_a else None,
            "b": inv_b[0] if inv_b else None,
            "both_resolved": pair_ok,
        },
        "inversion_31_vs_1111_kappa_a": {
            "d": ka["d_31_vs_1111"],
            "two_errors": ka["two_err_31_1111"],
            "resolved": ka["d_31_vs_1111"] < -ka["two_err_31_1111"],
            "ratio_to_two_errors": abs(ka["d_31_vs_1111"]) / ka["two_err_31_1111"],
        },
        "E_kappa_director_channel_kappa_a": Dpart,
        "E_kappa_pair_channel_kappa_a": Ppart,
        "E_kappa_eigenvalue_channel_kappa_a": Vpart,
        "E_kappa_pair_channel_spread": max(Ppart.values()) - min(Ppart.values()),
        "E_kappa_director_channel_spread": max(Dpart.values()) - min(Dpart.values()),
        "E_kappa_eigenvalue_channel_spread": max(Vpart.values()) - min(Vpart.values()),
        "E_kappa_director_channel_hedgehog_seed": anat["channels_hedgehog_seed_1_0.3_0"][
            "director"
        ],
        "log_term_predicted_spread": anat["log_term_per_unit_k2_16pi_b0^2_ln(L/xi)"][
            "predicted_spread_1111_to_4"
        ],
        "E_kappa_pair_channel_seeds": {
            q: anat["per_part"][q]["channels_seed"]["pair"] for q in PARTS
        },
        "note": (
            "the orders, both inversions of {1,1,1,1} against {2,1,1} (0.089 and 0.504 against two-error bounds "
            "1.4e-5 and 2.8e-5) and the resolved {3,1} < {1,1,1,1} at kappa a (0.00127, 40 times the two errors) "
            "all reproduce; REFUTED is fair for the prediction AS STATED (the total order by sum k^2). But the "
            "prediction was never testable on this box: the lattice-exact channel split of E_kappa gives the "
            "director-rotation channel 658 to 739 (spread 81), the eigenvalue channel 68 to 76 (spread 8) and "
            "the pair-rotation channel (the only piece the log term describes) 98 to 124 (spread 25), against a "
            "log-term spread of 9.2 for ln(22.4/1.5); the seeds' pair channel IS ordered by sum k^2 (30, 52, 57, "
            "75, 131) and the relaxation erases that order (103, 98, 99, 102, 124). The winding log term cannot "
            "be isolated on one box; it needs a box ladder at fixed partition"
        ),
    }

    # ---------------- claim 3: kappa a order vs R26-4 kappa 0 ----------------
    log("claim 3: the kappa-0 order")
    E0 = {}
    E0_own = {}
    for part in PARTS:
        t4 = f"P{part}_d0.3_w25_n32_L48"
        E0[part] = J4[t4]["E"]
        M4 = load_M(os.path.join(R26_4_DIR, t4 + ".npz"))
        E0_own[part] = {
            "E_quartic_own": own_split(M4, cfg, pot, 0.0)["E_quartic"],
            "E_kappa_own": own_split(M4, cfg, pot, 0.0)["E_kappa"],
            "gate_label": J4[t4]["gate_label"],
            "E_err": J4[t4]["E_err"],
        }
    order0 = sorted(PARTS, key=lambda q: E0[q])
    c3 = {
        "E_kappa0_R26_4_claimed": E0,
        "own_reads_of_R26_4_end_fields": E0_own,
        "order_kappa0": order0,
        "order_kappa_a": ka["order_total"],
        "order_kappa_b": kb["order_total"],
        "spread_kappa0": max(E0.values()) - min(E0.values()),
        "spread_kappa_a": ka["spread_total"],
        "spread_kappa_b": kb["spread_total"],
        "E_kappa_of_kappa0_end_fields": {q: E0_own[q]["E_kappa_own"] for q in PARTS},
        "kappa0_uncertified_rows": [q for q in PARTS if E0_own[q]["gate_label"] != "AT_GATE"],
    }
    out["claims"]["3"] = c3
    same = order0 == ka["order_total"]
    out["verdicts"]["3"] = {
        "verdict": "QUALIFIED" if same else "REFUTED",
        "order_equal_kappa_a": same,
        "order_equal_kappa_b": order0 == kb["order_total"],
        "spreads_kappa0_a_b": [c3["spread_kappa0"], c3["spread_kappa_a"], c3["spread_kappa_b"]],
        "note": (
            "the kappa-a order equals the kappa-0 order and the total spread narrows 0.439 -> 0.385, but two of the "
            "five kappa-0 rows ({4}, {2,1,1}) were FALLING (E_err 1e-3, never certified) and at kappa b the order "
            "changes ({2,2} drops below {1,1,1,1}) while the spread WIDENS to 0.647: 'converged them' holds at kappa a only"
        ),
    }

    # ---------------- claim 4: interiors under kappa ----------------
    log("claim 4: partitions on spheres and pair-gap closure clusters")
    c4 = {}
    claimed_int = {
        "1_1_1_1": [1, 1, 1, 1],
        "2_1_1": [2, 1, 1],
        "2_2": [2, 2],
        "3_1": [3, 1, 1, -1],
        "4": [4, 1, 1, 1, 1, -1, -1, -1, -1],
    }
    claimed_r9 = {
        "1_1_1_1": [1, 1, 1, 1],
        "2_1_1": [3, 1],
        "2_2": [4],
        "3_1": [1, 1, 1, 1],
        "4": [2, 1, 1],
    }
    ok4 = True
    for part in PARTS:
        M = fields[(part, "a")]
        reads = {}
        for R in SPHERES:
            rd = F0.partition_reader(M, cfg, R, delta)
            reads[f"R{R:g}"] = {
                "partition": rd["partition"],
                "total": rd["total_half_units"],
                "unreadable": rd["unreadable_plaquettes"],
                "flagged": rd.get("cap_ring_flagged"),
            }
        smallest = reads["R3"]["partition"]
        r9 = reads["R9"]["partition"]
        clusters = {f"R{R:g}": pair_gap_clusters(M, cfg, delta, R) for R in (3.0, 6.0, 9.0)}
        clusters_loose = {
            f"frac{fr:g}": pair_gap_clusters(M, cfg, delta, 9.0, frac=fr)["clusters"]
            for fr in (0.1, 0.25, 0.5)
        }
        largest = reads["R21"]["partition"]
        c4[part] = {
            "reader_by_sphere": reads,
            "claimed_interior": claimed_int[part],
            "claimed_r9": claimed_r9[part],
            "own_R3_matches_claimed_interior": smallest == claimed_int[part],
            "own_R9_matches_claimed_r9": r9 == claimed_r9[part],
            "pair_gap_closure_clusters": clusters,
            "own_R21_matches_claimed_interior": largest == claimed_int[part],
            "r9_cluster_count_by_gap_threshold_x_delta": clusters_loose,
            "r9_carriers_claimed": len(claimed_r9[part]),
        }
        ok4 = ok4 and (largest == claimed_int[part]) and (r9 == claimed_r9[part])
        log(
            f"  {part}: R3 {smallest} R9 {r9}; r9 closure clusters {clusters['R9']['clusters']} ({[c['cells'] for c in clusters['R9']['list']]})"
        )
    out["claims"]["4"] = c4
    out["verdicts"]["4"] = {
        "verdict": "QUALIFIED" if ok4 else "REFUTED",
        "reader_reproduces_r9_and_the_quoted_interior_as_R21": ok4,
        "R3_reads": {q: c4[q]["reader_by_sphere"]["R3"]["partition"] for q in PARTS},
        "R21_reads": {q: c4[q]["reader_by_sphere"]["R21"]["partition"] for q in PARTS},
        "flagged_cap_rings_R3_to_R12": {
            q: sorted(
                {
                    R
                    for R, v in c4[q]["reader_by_sphere"].items()
                    if v["flagged"] and float(R[1:]) <= 12
                }
            )
            for q in PARTS
        },
        "r9_closure_cluster_counts_by_threshold": {
            q: c4[q]["r9_cluster_count_by_gap_threshold_x_delta"] for q in PARTS
        },
        "r9_carriers_claimed": {q: len(claimed_r9[q]) for q in PARTS},
        "note": (
            "the r 9 partitions reproduce exactly; the quoted 'interior partitions' are the R 21 reads (the "
            "largest readable sphere, 1.4 h inside the pinned faces), NOT the inner spheres: at R 3 the reader "
            "gives [1,1,1,1], [4], [4], [1,1,1,1], [4] and the R 3 to R 12 reads of {2,1,1}, {2,2}, {4} rest on "
            "FLAGGED polar cap rings (fallback rings with up to 30 percent unreadable samples); the own pair-gap "
            "closure count on the r 9 shell finds fewer closed cores than carriers at gap < 0.1 delta, so the r 9 "
            "reorganization is a reader statement, not an independent core count"
        ),
    }

    # ---------------- claim 5: beta^2 and the virial ----------------
    log("claim 5: beta^2 and the virial")
    claimed_b2a = {"1_1_1_1": 0.029, "2_1_1": 0.019, "2_2": 0.111, "3_1": 0.014, "4": 0.139}
    c5 = {}
    ok5 = True
    for kk, kappa in KAPPAS.items():
        for part in PARTS:
            tag = tag_of(part, kappa)
            M = fields[(part, kk)]
            b2 = own_beta2_block(M, cfg, 6.0)
            b2_stack = F0.biaxiality_reads(M, cfg, delta)["block_mean_r6"]
            own = c1[tag]["own_split"]
            vir = own["E_curv"] / (3.0 * own["V4"])
            claimed = rows[tag]["own_reads"]
            c5[tag] = {
                "beta2_r6_own": b2,
                "beta2_r6_stack": b2_stack,
                "beta2_r6_claimed": claimed["biaxiality"]["block_mean_r6"],
                "virial_own": vir,
                "virial_claimed": claimed["virial"]["E_u_over_3V"],
            }
            ok5 = (
                ok5
                and abs(b2 - claimed["biaxiality"]["block_mean_r6"]) < 1e-9
                and abs(vir - claimed["virial"]["E_u_over_3V"]) < 1e-9
            )
            if kk == "a":
                ok5 = ok5 and abs(round(b2, 3) - claimed_b2a[part]) < 1.5e-3
    vir_a = [c5[tag_of(q, KAPPAS["a"])]["virial_own"] for q in PARTS]
    vir_b = [c5[tag_of(q, KAPPAS["b"])]["virial_own"] for q in PARTS]
    b2_b = [c5[tag_of(q, KAPPAS["b"])]["beta2_r6_own"] for q in PARTS]
    out["claims"]["5"] = c5
    out["verdicts"]["5"] = {
        "verdict": (
            "CONFIRMED"
            if ok5
            and 0.85 < min(vir_a)
            and max(vir_a) < 0.98
            and 0.73 < min(vir_b)
            and max(vir_b) < 0.81
            else "QUALIFIED"
        ),
        "virial_range_a": [min(vir_a), max(vir_a)],
        "virial_range_b": [min(vir_b), max(vir_b)],
        "beta2_r6_range_b": [min(b2_b), max(b2_b)],
        "note": "",
    }

    # ---------------- claim 6: the free-shell twins ----------------
    log("claim 6: the twins")
    c6 = {}
    for part in ("3_1", "1_1_1_1"):
        base = os.path.join(CENSUS_DIR, f"P{part}_d0.3_w25_n32_L48_k0.000e+00_free")
        if os.path.exists(base + ".npz"):
            M = load_M(base + ".npz")
            src, done, chunks = "end", None, None
        else:
            Z = np.load(base + "_stage.npz", allow_pickle=True)
            M = Z["M"].astype(np.float64)
            done = int(Z["done"])
            chunks = json.loads(str(Z["chunks"]))
            src = "stage"
        Mseed = load_M(base + "_seed.npz")
        own = own_split(M, cfg, pot, 0.0)
        rho, _ = R263.charge_density(M, cfg)
        _, _, _, r = radius(cfg)
        deg = {}
        for R in (9.0, 18.0):
            deg[f"R{R:g}"] = {
                "own_solid_angle": own_degree(M, cfg, R),
                "stack_degree_reader": F0.degree_reader(M, cfg, R)["degree"],
                "charge_inside_R26_3_density": float(rho[r < R].sum()),
            }
        # where is the hedgehog now: the cell with the smallest director gap, and the charge inside r < 6
        lam = np.linalg.eigvalsh(M[..., 1:, 1:])
        gdir = lam[..., 2] - lam[..., 1]
        X, Y, Z_, _ = radius(cfg)
        inner = ~pin
        idx = np.argmin(np.where(inner, gdir, 1e9))
        ii = np.unravel_index(idx, gdir.shape)
        t4 = f"P{part}_d0.3_w25_n32_L48"
        M4 = load_M(os.path.join(R26_4_DIR, t4 + ".npz"))
        lam4 = np.linalg.eigvalsh(M4[..., 1:, 1:])
        align = {}
        for R in (9.0, 18.0):
            shell = np.abs(r - R) <= 0.5 * cfg["h"] * np.sqrt(3.0)
            _, Vv = np.linalg.eigh(M[shell][:, 1:, 1:])
            dd = Vv[:, :, 2]
            rh = np.stack([X[shell], Y[shell], Z_[shell]], -1) / R
            align[f"R{R:g}"] = float(np.abs(np.einsum("ca,ca->c", dd, rh)).mean())
        c6[part] = {
            "mean_abs_director_dot_rhat_on_shell": align,
            "min_director_gap_pinned_kappa0_R26_4_field": float(
                (lam4[..., 2] - lam4[..., 1])[~pin].min()
            ),
            "source": src,
            "iters_done": done,
            "E_chunks_tail": [(c["iters"], c["E"]) for c in chunks[-6:]] if chunks else None,
            "E_quartic_own_latest": own["E_quartic"],
            "E_kappa_own_latest": own["E_kappa"],
            "E_seed_own": own_split(Mseed, cfg, pot, 0.0)["E_quartic"],
            "pin_shell_moved_maxabs": float(np.abs(M[pin] - Mseed[pin]).max()),
            "degree": deg,
            "charge_inside_r6": float(rho[r < 6.0].sum()),
            "charge_abs_total_inside_pin": float(np.abs(rho[~pin]).sum()),
            "min_director_gap_cell_xyz": [float(X[ii]), float(Y[ii]), float(Z_[ii])],
            "min_director_gap": float(gdir[ii]),
            "min_pair_gap": float((lam[..., 1] - lam[..., 0])[~pin].min()),
            "mean_block_eigs_r<6": [
                float(v) for v in np.linalg.eigvalsh(M[r < 6.0][:, 1:, 1:].mean(0))
            ],
            "reference_pinned_E_quartic_kappa0_R26_4": E0[part],
        }
        log(
            f"  twin {part} ({src}, {done}): E_quartic {own['E_quartic']:.4f}; degree r9 own {deg['R9']['own_solid_angle']['degree']:.3f} stack {deg['R9']['stack_degree_reader']:.3f} charge {deg['R9']['charge_inside_R26_3_density']:.3f}; r18 own {deg['R18']['own_solid_angle']['degree']:.3f} charge {deg['R18']['charge_inside_R26_3_density']:.3f}"
        )
    out["claims"]["6"] = c6
    out["verdicts"]["6"] = {
        "verdict": "CONFIRMED",
        "source": {q: [c6[q]["source"], c6[q]["iters_done"]] for q in c6},
        "E_quartic_latest": {q: c6[q]["E_quartic_own_latest"] for q in c6},
        "charge_inside_r9_r18": {
            q: [
                c6[q]["degree"]["R9"]["charge_inside_R26_3_density"],
                c6[q]["degree"]["R18"]["charge_inside_R26_3_density"],
            ]
            for q in c6
        },
        "min_director_gap_free_vs_pinned": {
            q: [c6[q]["min_director_gap"], c6[q]["min_director_gap_pinned_kappa0_R26_4_field"]]
            for q in c6
        },
        "mean_abs_align_r9_r18": {q: c6[q]["mean_abs_director_dot_rhat_on_shell"] for q in c6},
        "note": (
            "read at the latest stage files (the twins were still running): the energies fall as claimed and the "
            "unit charge is gone: no director core survives anywhere in the box (min director gap far above the "
            "pinned kappa-0 fields' core values), the charge integral inside r 9 is at the percent level, the "
            "mean |d . rhat| on r 9 is far from the hedgehog's 1. The sphere degree reads (own and R26-0) return "
            "integers that mean nothing here: the outward lift is undefined where |d . rhat| ~ 1e-4, as it is on "
            "every sphere of a near-uniform director"
        ),
    }

    # ---------------- claim 7: the delta arm ----------------
    log("claim 7: the delta arm")
    c7 = {}
    claimed7 = {
        "S1_d0.1_w25_n32_L48": (17467, 10369, 20981),
        "S1_d0.03_w25_n32_L48": (13569, 3135, 13430),
    }
    for tag, rr in JA["rows"].items():
        d = rr["delta"]
        cfgd, pd, potd = cfg_pot(d)
        M = load_M(os.path.join(ARM_DIR, tag + ".npz"))
        pg = F0.physical_generator(M, cfgd)
        og = own_generator(M, cfgd)
        own = own_split(M, cfgd, potd, 0.0)
        gate = gate_from_chunks(rr["chunks"])
        kick = kick_from_row(rr["kick"], rr["kick"]["E_ref_reslaved"])
        c7[tag] = {
            "delta": d,
            "E_quartic_own": own["E_quartic"],
            "E_claimed": rr["E"],
            "stack_generator": {k: pg[k] for k in ("internal", "orbital", "rigid")},
            "stack_generator_within_r18": pg["within_r"]["18"],
            "own_generator": og,
            "claimed_C_int_orb_rigid": claimed7[tag],
            "gate_from_chunks": gate,
            "kick_from_chunks": kick,
            "claimed_gate_label": rr["gate_label"],
            "claimed_kick_label": rr["kick_label"],
            "partition_r9_stack": F0.partition_reader(M, cfgd, 9.0, d)["partition"],
            "iters": rr["iters"],
        }
        log(
            f"  {tag}: E {own['E_quartic']:.6f} (claimed {rr['E']:.6f}); C stack {pg['internal']:.0f}/{pg['orbital']:.0f}/{pg['rigid']:.0f} own {og['internal']:.0f}/{og['orbital']:.0f}/{og['rigid']:.0f}; gate last_ok {gate['last_at_gate']}; kick own {kick['own_label']} reconv {kick['own_reconverged']} fmax_after {kick['fmax_after']:.1e}"
        )
    Mref = load_M(S1_D03_END)
    pgr = F0.physical_generator(Mref, cfg)
    ogr = own_generator(Mref, cfg)
    c7["S1_d0.3_reference_R26_3"] = {
        "stack_generator": {k: pgr[k] for k in ("internal", "orbital", "rigid")},
        "own_generator": ogr,
        "claimed": (14071, 1403, 14031),
    }
    out["claims"]["7"] = c7
    crig = {
        0.03: c7["S1_d0.03_w25_n32_L48"]["stack_generator"]["rigid"],
        0.1: c7["S1_d0.1_w25_n32_L48"]["stack_generator"]["rigid"],
        0.3: pgr["rigid"],
    }
    crig_own = {
        0.03: c7["S1_d0.03_w25_n32_L48"]["own_generator"]["rigid"],
        0.1: c7["S1_d0.1_w25_n32_L48"]["own_generator"]["rigid"],
        0.3: ogr["rigid"],
    }
    nonmono = (crig[0.1] > crig[0.03]) and (crig[0.1] > crig[0.3])
    nonmono_own = (crig_own[0.1] > crig_own[0.03]) and (crig_own[0.1] > crig_own[0.3])
    k03 = c7["S1_d0.03_w25_n32_L48"]["kick_from_chunks"]
    out["verdicts"]["7"] = {
        "verdict": "QUALIFIED",
        "C_rigid_by_delta_stack": crig,
        "C_rigid_by_delta_own_stencil": crig_own,
        "non_monotone_stack": nonmono,
        "non_monotone_own": nonmono_own,
        "d0.03_kick_fmax_after": k03["fmax_after"],
        "d0.03_kick_reconverged_own": k03["own_reconverged"],
        "d0.03_C_rigid_at_8000_vs_end": [
            JA["collect"]["rows"]["S1_d0.03_w25_n32_L48"]["C_rigid_at_8000"],
            crig[0.03],
        ],
        "note": (
            "the three inertias reproduce with the stack and with the own face-stencil (C_orb differs by 0.1 "
            "percent, faces only) and both descents pass the gate from their chunks; the non-monotonicity is "
            "robust (C_rigid at delta 0.1 is 50 percent above both neighbors) but the 0.03 < 0.3 side (13430 vs "
            "14031, 4 percent) is inside the delta-0.03 row's own stage-to-end drift (12903 at 8000 iterations "
            "-> 13430 at 13750); the delta-0.03 kick returned to the energy (4e-6 above the reference) with fmax "
            "1.1e-4, a 10 percent miss of the gate, so STABLE rests on the energy, not on a reconverged kicked field"
        ),
    }

    # ---------------- the own descent on the {2,2} kappa-b saddle ----------------
    log("own descent: continue the {2,2} kappa-b kicked end field")
    tag = tag_of("2_2", KAPPAS["b"])
    M = fields[("2_2", "b")]
    Mf, desc = own_descent(M, cfg, p, pot, rows[tag]["kappa"], free, OWN_DESCENT_ITERS, tag)
    desc["partition_R3_after"] = F0.partition_reader(Mf, cfg, 3.0, delta)["partition"]
    desc["partition_R9_after"] = F0.partition_reader(Mf, cfg, 9.0, delta)["partition"]
    desc["split_after"] = own_split(Mf, cfg, pot, rows[tag]["kappa"])
    desc["E_ref_gate"] = rows[tag]["E_gate"]
    desc["E_end_claimed"] = rows[tag]["E"]
    desc["E_211_kappa_b"] = c1[tag_of("2_1_1", KAPPAS["b"])]["own_split"]["E_total"]
    desc["E_31_kappa_b"] = c1[tag_of("3_1", KAPPAS["b"])]["own_split"]["E_total"]
    out["claims"]["1"]["own_descent_2_2_kappa_b"] = desc
    out["verdicts"]["1"]["own_descent_2_2_kappa_b"] = {
        "E_start": desc["E_start"],
        "E_end": desc["E_end"],
        "drop": desc["drop"],
        "gmax_end": desc["gmax_spatial_end"],
        "partition_R3_R9_after": [desc["partition_R3_after"], desc["partition_R9_after"]],
    }
    out["verdicts"]["1"]["note"] += (
        f"; own 300-iteration L-BFGS continuation of the {{2,2}} kappa-b kicked end field: E "
        f"{desc['E_start']:.6f} -> {desc['E_end']:.6f} (drop {desc['drop']:.1e}), gmax still "
        f"{desc['gmax_spatial_end']:.1e}, partition [4] at R 3 and R 9 unchanged: the saddle escape is a slow "
        "creep, the row's E 16.4522 is an unconverged upper bound, SADDLE / not reconverged is the right label"
    )

    # ---------------- the summary ----------------
    verdict_counts = {}
    for k, v in out["verdicts"].items():
        verdict_counts[v["verdict"]] = verdict_counts.get(v["verdict"], 0) + 1
    out["summary"] = {
        "counts": verdict_counts,
        "refuted": [k for k, v in out["verdicts"].items() if v["verdict"] == "REFUTED"],
        "qualified": {
            k: v["note"] for k, v in out["verdicts"].items() if v["verdict"] == "QUALIFIED"
        },
        "blind_spots": [
            "E_err is the last chunk's drop as the rows define it; no independent Richardson-style estimate of the remaining tail was made beyond the geometric drop ratio",
            f"the twins were read at their {c6['3_1']['iters_done']}-iteration stage files (still running); a later field can differ",
            "the own descent ran 300 L-BFGS iterations from the audited {2,2} kappa-b end field, not from an independent seed",
            "the channel split of E_kappa assigns each bond's mixing by the eigenframes; where two eigenvalues are nearly degenerate the frame is ill-conditioned but the channel weight (l_k - l_l)^2 vanishes there",
        ],
        "wall_s": time.time() - T0,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    log(f"wrote {OUT_JSON}")
    return out


if __name__ == "__main__":
    main()
