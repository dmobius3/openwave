"""M5.32 R22-1 adversarial audit: an independent attempt to refute the core-seed claims.

Audited: m5_32_r22_1_cores.py (rows data/m5_32_r22_1_cores.json and _ext.json, collect
data/m5_32_r22_1_cores_collect.json, end fields data/m5_32_r22_1/<tag>.npz) and
m5_32_r22_1b_split_halo.py (data/m5_32_r22_1b_split_halo.json).

The production stack is consumed read-only (R21 cfg_of / params_of, R20 energy_grad /
energy_parts, R0 roots_of, B3 coords / pin_shell / W1) and the seeds come from the audited
seed_core. Own pieces:
    energy      own one-sided stencils, own commutator sum, own V4 (Q1)
    reduced     the L-BFGS descent on the 6 spatial entries per free cell with M_00 slaved,
    polish      per cell, to its exact V4 minimizer (M_00 does not enter the curvature when
                M_0i = 0, so the reduced function has the same minima; its gradient is the
                spatial block of the production gradient by the envelope theorem)
    degree      Van Oosterom-Strackee solid angles on the triangulated cube surface, the
                eigenvector oriented on the surface graph itself (breadth-first), the
                conflicting links counted; plus a signed preimage count of probe directions
    split test  E(+a) + E(-a) - 2 E(0) of a pair-split perturbation on the exact uniaxial
                hedgehog exterior, against the second variation derived here:
                    e_2 = 8 (1 - delta)^2 / r^2 [ T + 2 eps_r^2 + 8 eps^2 omega_r^2 ]
                          - 16 (1 - delta)^2 eps^2 / r^4,
                T the tangential spin-weight-2 form with angular eigenvalue lam / r^2 on the
                l = 2 profile sin^2(theta), lam = 0 or 4 by the sign the geometry fixes

Checks print as  Q<id>.<n> PASS|FAIL ...  and land in data/m5_32_r22_1_audit.json.
A FAIL means the claim as worded was broken; the value field carries the counter-number.

Run: python m5_32_r22_1_audit.py   (two worker processes at most, about 15 minutes)
     python m5_32_r22_1_audit.py static   (the checks without the polish pool, no JSON)
"""

import importlib.util
import itertools
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import curve_fit, least_squares, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
NPZ = os.path.join(DATA, "m5_32_r22_1")
ROWS_JSON = os.path.join(DATA, "m5_32_r22_1_cores.json")
EXT_JSON = os.path.join(DATA, "m5_32_r22_1_cores_ext.json")
COLLECT_JSON = os.path.join(DATA, "m5_32_r22_1_cores_collect.json")
HALO_JSON = os.path.join(DATA, "m5_32_r22_1b_split_halo.json")
OUT_JSON = os.path.join(DATA, "m5_32_r22_1_audit.json")

G_T = 8.0
DELTA = 0.3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
IU3 = np.triu_indices(3)
OFF3 = np.where(IU3[0] != IU3[1], 2.0, 1.0)
IU4 = np.triu_indices(4)
OFF4 = np.where(IU4[0] != IU4[1], 2.0, 1.0)
T0 = time.time()
CHECKS = []


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CORES = _load("m5_32_r22_1_cores", "m5_32_r22_1_cores.py")
R21, R20, B3, R0 = CORES.R21, CORES.R20, CORES.B3, CORES.R0
W1 = B3.W1


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def check(cid, claim, method, value, expected, ok):
    verdict = "PASS" if ok else "FAIL"
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": value,
            "expected": expected,
            "verdict": verdict,
        }
    )
    print(f"{cid} {verdict} {claim} | value {json.dumps(value, default=str)[:420]}", flush=True)


def setup(n, L, delta, w1s):
    cfg = R21.cfg_of(n, L, G_T, delta)
    p = R21.params_of(G_T, delta)
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * w1s)
    return cfg, p, pot


def field(tag):
    return np.load(os.path.join(NPZ, tag + ".npz"))["M"]


def radius(n, h):
    X, Y, Z = B3.coords(n, h)
    return X, Y, Z, np.sqrt(X * X + Y * Y + Z * Z)


# ================= own energy (Q1) =================
def own_energy(M, h, roots, w):
    """own stencils: forward (zero at the last cell) and backward (zero at the first), averaged."""

    def one_sided(f, ax, fwd):
        d = np.diff(f, axis=ax) / h
        pad = [(0, 0)] * f.ndim
        pad[ax] = (0, 1) if fwd else (1, 0)
        return np.pad(d, pad)

    e_curv = 0.0
    for fwd in (True, False):
        A = [one_sided(M, ax, fwd) for ax in range(3)]
        for i, j in ((0, 1), (0, 2), (1, 2)):
            F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
            e_curv += 0.5 * 4.0 * np.sum((ETA @ F @ ETA) * F)
    N = M @ ETA
    P = N.copy()
    v = 0.0
    for k in range(1, 5):
        if k > 1:
            P = P @ N
        v = v + (np.trace(P, axis1=-2, axis2=-1) - sum(q**k for q in roots)) ** 2
    return h**3 * e_curv, h**3 * w * float(np.sum(v))


def v4_cell(x, roots, w):
    """own V4 of one cell from the 10 upper-triangle entries (no volume factor)."""
    Mc = np.zeros((4, 4))
    Mc[IU4] = x
    Mc = Mc + Mc.T - np.diag(np.diag(Mc))
    N = Mc @ ETA
    return w * sum(
        (np.trace(np.linalg.matrix_power(N, k)) - sum(q**k for q in roots)) ** 2
        for k in range(1, 5)
    )


# ================= the reduced polish (own instrument) =================
def solve_m00(S, roots, m_start):
    """per cell argmin over m of sum_p ((-m)^p + tr S^p - C_p)^2 (Newton, vectorized)."""
    cp = [sum(q**k for q in roots) for k in range(1, 5)]
    s, P = [], S
    for k in range(1, 5):
        if k > 1:
            P = P @ S
        s.append(np.trace(P, axis1=-2, axis2=-1))
    x = -m_start
    for _ in range(8):
        r = [x**k + s[k - 1] - cp[k - 1] for k in range(1, 5)]
        d = [k * x ** (k - 1) for k in range(1, 5)]
        dd = [k * (k - 1) * x ** (k - 2) if k > 1 else 0.0 * x for k in range(1, 5)]
        g = sum(2.0 * r[k] * d[k] for k in range(4))
        H = sum(2.0 * d[k] ** 2 + 2.0 * r[k] * dd[k] for k in range(4))
        x = x - g / H
    return -x


def polish_job(job):
    """worker: L-BFGS from an end field; mode 'reduced' (M_00 slaved) or 'plain' (10 entries)."""
    cfg, p, pot = setup(32, 48.0, DELTA, job["w1s"])
    M = field(job["tag"])
    mask = R21.free_mask(cfg, True)
    st = {"n_eval": 0, "it": 0, "trace": [], "last": None}
    reduced = job["mode"] == "reduced"

    def build(x):
        Mx = M.copy()
        sub = Mx[mask]
        if reduced:
            F = np.zeros((sub.shape[0], 3, 3))
            F[:, IU3[0], IU3[1]] = x.reshape(-1, 6)
            F[:, IU3[1], IU3[0]] = x.reshape(-1, 6)
            sub[:, 1:, 1:] = F
            sub[:, 0, 0] = solve_m00(F, pot[1], sub[:, 0, 0])
        else:
            F = np.zeros((sub.shape[0], 4, 4))
            F[:, IU4[0], IU4[1]] = x.reshape(-1, 10)
            F[:, IU4[1], IU4[0]] = x.reshape(-1, 10)
            sub = F
        Mx[mask] = sub
        return Mx

    def fun(x):
        Mx = build(x)
        E, G, _ = R20.energy_grad(Mx, cfg, p, pot)
        st["n_eval"] += 1
        Gm = G[mask]
        Gs = Gm.copy()
        Gs[:, 0, 0] = 0.0
        st["last"] = (float(E), float(np.abs(Gm).max()), float(np.abs(Gs).max()))
        if reduced:
            return float(E), (Gm[:, 1:, 1:][:, IU3[0], IU3[1]] * OFF3).ravel()
        return float(E), (Gm[:, IU4[0], IU4[1]] * OFF4).ravel()

    def cb(_xk):
        st["it"] += 1
        st["trace"].append((st["it"], st["n_eval"]) + st["last"])

    sub0 = M[mask]
    x0 = (sub0[:, 1:, 1:][:, IU3[0], IU3[1]] if reduced else sub0[:, IU4[0], IU4[1]]).ravel()
    E_in = fun(x0.copy())[0]
    first = st["last"]
    res = minimize(
        fun,
        x0.copy(),
        jac=True,
        method="L-BFGS-B",
        callback=cb,
        options={
            "maxcor": 20,
            "maxiter": job["iters"],
            "maxfun": 3 * job["iters"],
            "gtol": 1e-14,
            "ftol": 1e-16,
        },
    )
    Mp = build(np.asarray(res.x))
    E_end, G_end, _ = R20.energy_grad(Mp, cfg, p, pot)
    Gs = G_end[mask].copy()
    fmax_all = float(np.abs(Gs).max())
    Gs[:, 0, 0] = 0.0
    return {
        "job": job,
        "E_in": float(E_in),
        "in_E_fmax_fmaxsp": first,
        "E_end": float(E_end),
        "fmax_end_all": fmax_all,
        "fmax_end_spatial": float(np.abs(Gs).max()),
        "iters": st["it"],
        "n_eval": st["n_eval"],
        "trace": st["trace"],
        "M": Mp,
        "wall_s": round(time.time() - T0, 1),
    }


# ================= own symmetry images, degree reader =================
def images12(S):
    """the 6 axis permutations times the full inversion, (R S)(x) = R S(R^T x) R^T."""
    out = []
    for perm in itertools.permutations(range(3)):
        F = np.transpose(S, axes=(perm[0], perm[1], perm[2], 3, 4))
        F = F[..., list(perm), :][..., :, list(perm)]
        out.append((perm, +1, F))
        out.append((perm, -1, F[::-1, ::-1, ::-1]))
    return out


def rms12(Ma, Mb, mask):
    A = Ma[..., 1:, 1:]
    best = None
    for perm, sgn, T in images12(Mb[..., 1:, 1:]):
        v = float(np.sqrt(np.mean(np.sum((A - T) ** 2, axis=(-1, -2))[mask])))
        if best is None or v < best[0]:
            best = (v, list(perm), sgn)
    ident = float(np.sqrt(np.mean(np.sum((A - Mb[..., 1:, 1:]) ** 2, axis=(-1, -2))[mask])))
    return {"rms_min12": best[0], "perm": best[1], "inversion": best[2], "rms_identity": ident}


def surface_degree(v, lo, hi, probes=7, seed=3):
    """degree of the unit field v through the closed surface of the cell cube [lo, hi]^3."""
    nodes, tris = {}, []

    def nid(c):
        return nodes.setdefault(c, len(nodes))

    m = hi - lo
    for ax in range(3):
        u_ax, v_ax = (ax + 1) % 3, (ax + 2) % 3  # u x v = +ax
        for side in (hi, lo):
            for a in range(m):
                for b in range(m):
                    q = []
                    for da, db in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        c = [0, 0, 0]
                        c[ax], c[u_ax], c[v_ax] = side, lo + a + da, lo + b + db
                        q.append(nid(tuple(c)))
                    if side == lo:
                        q = q[::-1]
                    tris.append((q[0], q[1], q[2]))
                    tris.append((q[0], q[2], q[3]))
    cells = np.array(list(nodes.keys()))
    vs = v[cells[:, 0], cells[:, 1], cells[:, 2]].copy()
    links = set()
    for a, b, c in tris:
        for e in ((a, b), (b, c), (c, a)):
            links.add((min(e), max(e)))
    adj = {}
    for a, b in links:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    sign = np.zeros(len(nodes))
    sign[0] = 1.0
    queue = [0]
    while queue:
        a = queue.pop(0)
        for b in adj[a]:
            if sign[b] == 0.0:
                sign[b] = sign[a] * (1.0 if np.dot(vs[a], vs[b]) >= 0.0 else -1.0)
                queue.append(b)
    vo = vs * sign[:, None]
    conflicts = sum(1 for a, b in links if np.dot(vo[a], vo[b]) < 0.0)
    # global sign: outward at the first node
    ctr = 0.5 * (lo + hi)
    outward = np.sum(vo * (cells - ctr)) >= 0.0
    if not outward:
        vo = -vo
    T = np.array(tris)
    p, q, r = vo[T[:, 0]], vo[T[:, 1]], vo[T[:, 2]]
    num = np.einsum("ia,ia->i", p, np.cross(q, r))
    den = (
        1.0
        + np.einsum("ia,ia->i", p, q)
        + np.einsum("ia,ia->i", q, r)
        + np.einsum("ia,ia->i", r, p)
    )
    omega = float(np.sum(2.0 * np.arctan2(num, den))) / (4.0 * np.pi)
    # signed preimage count of probe directions (an integer reader, independent of the angles)
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(probes):
        t = rng.normal(size=3)
        t /= np.linalg.norm(t)
        Mx = np.stack([p, q, r], axis=-1)  # columns
        det = np.linalg.det(Mx)
        ok = np.abs(det) > 1e-14
        co = np.linalg.solve(Mx[ok], np.broadcast_to(t, (int(ok.sum()), 3))[..., None])[..., 0]
        inside = np.all(co > 0.0, axis=1)
        counts.append(int(np.sum(np.sign(det[ok][inside]))))
    return {
        "solid_angle_degree": omega,
        "conflicting_links": int(conflicts),
        "links": len(links),
        "probe_counts": counts,
        "degree": (int(round(omega)) if conflicts == 0 else "undefined"),
    }


def degree_reads(M):
    """degree, conflicting links and the smallest (top - mid) gap on the r 6, 9, 12 cell cubes."""
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    rec = {}
    for nm, lo, hi in (("r6", 12, 19), ("r9", 10, 21), ("r12", 8, 23)):
        rec[nm] = surface_degree(vec[..., :, 2], lo, hi)
        box = np.zeros(lam.shape[:3], dtype=bool)
        box[lo : hi + 1, lo : hi + 1, lo : hi + 1] = True
        box[lo + 1 : hi, lo + 1 : hi, lo + 1 : hi] = False
        rec[nm]["min_gap_top_mid_on_surface"] = float((lam[..., 2] - lam[..., 1])[box].min())
    return rec


# ================= Q1 =================
def all_rows():
    rows = {}
    for path in (ROWS_JSON, EXT_JSON):
        with open(path) as f:
            rows.update(json.load(f)["rows"])
    return {t: r for t, r in rows.items() if os.path.exists(os.path.join(NPZ, t + ".npz"))}


def q1(rows):
    dev, pins, m0i, table = [], {}, 0.0, {}
    for tag, r in rows.items():
        cfg, p, pot = setup(r["n"], r["L"], r["delta"], r["w1s"])
        M = field(tag)
        parts = R20.energy_parts(M, cfg, p, pot)
        Eg = R20.energy_grad(M, cfg, p, pot, need_grad=False)[0]
        table[tag] = {k: parts[k] for k in ("E_total", "E_curv", "V")}
        dev.append(
            max(
                abs(parts["E_total"] - r["E"]),
                abs(parts["E_curv"] - r["end_E"]["E_curv"]),
                abs(parts["V"] - r["end_E"]["V"]),
                abs(Eg - r["E"]),
            )
        )
        m0i = max(m0i, float(np.abs(M[..., 0, 1:]).max()))
        if r["bnd"] == "pin":
            seed = CORES.seed_core(cfg, r["seed"], r["delta"])
            ps = B3.pin_shell(cfg["n"], cfg["h"])
            pins[tag] = bool(np.array_equal(M[ps], seed[ps]))
    check(
        "Q1.1",
        "every end field's E, E_curv, V reproduce the JSON (production energy_parts and energy_grad)",
        f"{len(rows)} end fields recomputed",
        {"max_abs_dev": max(dev), "fields": len(rows)},
        "< 1e-9",
        max(dev) < 1e-9,
    )
    quoted = {
        "bia_pin_d0.3_w25_n32_L48": 5.227,
        "perm_pin_d0.3_w25_n32_L48": 5.323,
        "obl_pin_d0.3_w25_n32_L48": 5.503,
        "rad_pin_d0.3_w25_n32_L48": 5.709,
        "big_pin_d0.3_w25_n32_L48": 6.001,
        "perm_pin_d0.3_w25_n32_L48_x9000": 4.380,
        "bia_pin_d0.3_w25_n32_L48_x9000": 4.409,
        "obl_pin_d0.3_w25_n32_L48_x9000": 4.445,
        "rad_pin_d0.3_w25_n32_L48_x9000": 5.044,
        "bia_pin_d0.3_w1_n32_L48": 2.27921,
        "rad_pin_d0.3_w1_n32_L48": 2.27999,
        "bia_free_d0.3_w1_n32_L48": 2.2155,
        "rad_free_d0.3_w1_n32_L48": 2.2188,
        "rad_pin_d0.89_w25_n32_L48": 0.00407,
    }
    bad = {}
    for t, q in quoted.items():
        digits = len(str(q).split(".")[1])
        if t not in table or round(table[t]["E_total"], digits) != q:
            bad[t] = table.get(t, {}).get("E_total")
    check(
        "Q1.2",
        "the quoted energies match the recomputed ones to the quoted digits",
        "rounding of the recomputed E_total",
        {"mismatches": bad},
        "none",
        not bad,
    )
    own = {}
    for tag in ("bia_pin_d0.3_w25_n32_L48_x9000", "rad_free_d0.3_w1_n32_L48"):
        r = rows[tag]
        cfg, p, pot = setup(r["n"], r["L"], r["delta"], r["w1s"])
        ec, v = own_energy(field(tag), cfg["h"], pot[1], pot[2])
        own[tag] = {
            "own_E_curv": ec,
            "own_V": v,
            "dev_curv": abs(ec - table[tag]["E_curv"]),
            "dev_V": abs(v - table[tag]["V"]),
        }
    ok = all(o["dev_curv"] < 1e-9 and o["dev_V"] < 1e-9 for o in own.values())
    check(
        "Q1.3",
        "an independent energy (own stencils, own commutators, own V4) agrees",
        "own_energy on one pinned W1 x 25 and one free W1 field",
        own,
        "< 1e-9",
        ok,
    )
    check(
        "Q1.4",
        "the pinned shell is bitwise the seed's on every 'pin' row",
        "np.array_equal on B3.pin_shell",
        {"rows": len(pins), "not_equal": [t for t, b in pins.items() if not b]},
        "all equal",
        all(pins.values()),
    )
    check("Q1.5", "M_0i = 0 on every end field", "max |M_0i|", m0i, "0", m0i == 0.0)
    return table


# ================= Q2 =================
def q2(rows):
    drops = {}
    for tag, r in rows.items():
        if r["w1s"] == 25.0 and r["delta"] == 0.3:
            ch = r["polish_chunks"]
            drops[tag] = ch[-2]["E_end"] - ch[-1]["E_end"]
    lo, hi = min(drops.values()), max(drops.values())
    out_of_range = {t: round(d, 4) for t, d in drops.items() if not 0.03 <= d <= 0.245}
    check(
        "Q2.1",
        "every W1 x 25 row still falls by 0.03 to 0.24 per 500 iterations at its end",
        "last-chunk energy drops from polish_chunks",
        {"min": lo, "max": hi, "outside_0.03_0.24": out_of_range},
        "all inside; all positive",
        not out_of_range,
    )
    check(
        "Q2.1b",
        "no W1 x 25 row is converged (every last-chunk drop is positive and far above 1e-6)",
        "the same drops",
        {"min_drop": lo},
        "> 1e-3",
        lo > 1e-3,
    )

    def f_exp(k, e, a, tau):
        return e + a * np.exp(-k / tau)

    def f_pow(k, e, a, pw):
        return e + a * k ** (-pw)

    def f_hyp(k, e, a, k0):
        return e + a / (k + k0)

    fits = {}
    for tag, r in rows.items():
        if not tag.endswith("_x9000"):
            continue
        E = np.array([c["E_end"] for c in r["polish_chunks"]])
        k = np.arange(1, len(E) + 1, dtype=float)
        rec = {}
        for start in (6, 9):
            kk, ee = k[start - 1 :], E[start - 1 :]
            for name, fn, p0, bnd in (
                ("exp", f_exp, (ee[-1] - 0.5, 2.0, 8.0), ((-50, 0, 0.5), (ee[-1], 1e3, 1e4))),
                ("power", f_pow, (ee[-1] - 1.0, 5.0, 0.5), ((-50, 0, 0.01), (ee[-1], 1e4, 10))),
                (
                    "hyperbola",
                    f_hyp,
                    (ee[-1] - 1.0, 20.0, 5.0),
                    ((-50, 0, -0.9 * start), (ee[-1], 1e5, 1e4)),
                ),
            ):
                try:
                    popt, _ = curve_fit(fn, kk, ee, p0=p0, bounds=bnd, maxfev=20000)
                    res = float(np.sqrt(np.mean((fn(kk, *popt) - ee) ** 2)))
                    rec[f"{name}_from{start}"] = {"E_inf": float(popt[0]), "rms_resid": res}
                except Exception as e:  # noqa: BLE001
                    rec[f"{name}_from{start}"] = {"error": repr(e)[:80]}
        vals = [v["E_inf"] for v in rec.values() if "E_inf" in v]
        rec["E_inf_min"], rec["E_inf_max"], rec["E_9000"] = min(vals), max(vals), float(E[-1])
        fits[tag.split("_")[0]] = rec
    spread = {s: round(v["E_inf_max"] - v["E_inf_min"], 3) for s, v in fits.items()}
    overlap_lo = max(v["E_inf_min"] for v in fits.values())
    overlap_hi = min(v["E_inf_max"] for v in fits.values())
    check(
        "Q2.2",
        "the 9000-iteration traces do not fix a floor (extrapolations by decay law disagree)",
        "three decay laws (exponential, power, hyperbola) from chunk 6 and from chunk 9, scipy curve_fit",
        {
            "fits": fits,
            "spread_of_E_inf_per_seed": spread,
            "common_window_all_seeds": [overlap_lo, overlap_hi],
        },
        "spread per seed above the 0.03 to 0.66 gaps between the seeds",
        min(spread.values()) > 0.3,
    )
    return fits


# ================= Q3 =================
def q3(rows):
    out = {}
    for tag in [t for t in rows if t.endswith("_x9000")] + ["big_pin_d0.3_w25_n32_L48"]:
        r = rows[tag]
        cfg, p, pot = setup(r["n"], r["L"], r["delta"], r["w1s"])
        M = field(tag)
        _, G, _ = R20.energy_grad(M, cfg, p, pot)
        mask = R21.free_mask(cfg, True)
        Ga = np.abs(G) * mask[..., None, None]
        i = np.unravel_index(np.argmax(Ga), Ga.shape)
        Gs = Ga.copy()
        Gs[..., 0, 0] = 0.0
        H00 = cfg["h"] ** 3 * pot[2] * 2.0 * sum((k * G_T ** (k - 1)) ** 2 for k in range(1, 5))
        out[tag] = {
            "fmax_all": float(Ga.max()),
            "fmax_json": r["fmax_end"],
            "argmax_entry": [int(i[3]), int(i[4])],
            "argmax_cell": [int(i[0]), int(i[1]), int(i[2])],
            "fmax_excluding_00": float(Gs.max()),
            "rms_G00": float(np.sqrt(np.mean(G[..., 0, 0][mask] ** 2))),
            "rms_G_spatial": float(np.sqrt(np.mean(G[..., 1:, 1:][mask] ** 2))),
            "predicted_E_decrease_from_M00": float(np.sum(G[..., 0, 0][mask] ** 2) / (2.0 * H00)),
            "max_M00_error": float(np.abs(G[..., 0, 0][mask]).max() / H00),
        }
    ok = all(
        v["argmax_entry"] == [0, 0] and abs(v["fmax_all"] - v["fmax_json"]) < 1e-9
        for v in out.values()
    )
    check(
        "Q3.1",
        "the end-of-chunk max |G| is carried by the (0,0) entry (the stiff M_00 direction)",
        "production gradient on the end fields, arg max over entries and free cells",
        out,
        "argmax entry (0,0); fmax reproduces the JSON",
        ok,
    )
    sp = [v["fmax_excluding_00"] for v in out.values()]
    check(
        "Q3.2",
        "with the (0,0) entry excluded the residual is far below the recorded fmax_end, yet still above the 1e-3 gate",
        "max |G| over the 9 other entries on the free cells; the M_00 part as a predicted energy decrease G_00^2 / (2 H_00)",
        {
            "fmax_excluding_00_min": min(sp),
            "fmax_excluding_00_max": max(sp),
            "fmax_recorded_min": min(v["fmax_all"] for v in out.values()),
            "fmax_recorded_max": max(v["fmax_all"] for v in out.values()),
            "predicted_E_decrease_from_M00_max": max(
                v["predicted_E_decrease_from_M00"] for v in out.values()
            ),
        },
        "1e-3 < spatial residual < 0.1 on every row",
        1e-3 < min(sp) and max(sp) < 0.1,
    )
    # V4 Hessian at the vacuum: own finite differences on one cell
    roots = (-G_T, 1.0, DELTA, DELTA)
    x0 = np.zeros(10)
    x0[[0, 4, 7, 9]] = [G_T, 1.0, DELTA, DELTA]
    t = 2e-5
    Hm = np.zeros((10, 10))
    for a in range(10):
        for b in range(10):
            ea, eb = np.eye(10)[a] * t, np.eye(10)[b] * t
            Hm[a, b] = (
                v4_cell(x0 + ea + eb, roots, 1.0)
                - v4_cell(x0 + ea - eb, roots, 1.0)
                - v4_cell(x0 - ea + eb, roots, 1.0)
                + v4_cell(x0 - ea - eb, roots, 1.0)
            ) / (4 * t * t)
    ev = np.linalg.eigvalsh(Hm)
    zero = int(np.sum(np.abs(ev) < 1e-2))  # w units; the massive three are 3.09, 35.0, 8.46e6
    check(
        "Q3.3",
        "V4 Hessian along M_00 is 8.46e6 w; 7 of the 10 directions are V4-flat at the uniaxial vacuum",
        "own V4, central second differences on one cell",
        {
            "H_00_over_w": float(Hm[0, 0]),
            "closed_form": 2.0 * sum((k * G_T ** (k - 1)) ** 2 for k in range(1, 5)),
            "zero_modes": zero,
            "eigenvalues_over_w": [float(x) for x in ev],
        },
        "8.46e6; 7",
        abs(Hm[0, 0] / 8.46285e6 - 1) < 1e-3 and zero == 7,
    )
    return out


# ================= Q4 (static part) =================
def q4_static(rows):
    cfg, p, pot = setup(32, 48.0, DELTA, 1.0)
    ta, tb = "rad_pin_d0.3_w1_n32_L48", "bia_pin_d0.3_w1_n32_L48"
    A, B = field(ta), field(tb)
    X, Y, Z, r = radius(32, cfg["h"])
    mask = r < 12.0
    fm = R21.free_mask(cfg, True)
    res = {}
    for nm, F in (("rad", A), ("bia", B)):
        E, G, _ = R20.energy_grad(F, cfg, p, pot)
        res[nm] = {"E": float(E), "fmax": float((np.abs(G) * fm[..., None, None]).max())}
    rel = abs(res["rad"]["E"] - res["bia"]["E"]) / res["bia"]["E"]
    check(
        "Q4.1",
        "both W1 pinned rows are under the 1e-3 max |G| gate and agree in energy to 0.03 percent",
        "production gradient on the two end fields",
        {"rows": res, "rel_dE": rel},
        "fmax < 1e-3; rel dE about 3.4e-4",
        res["rad"]["fmax"] < 1e-3 and res["bia"]["fmax"] < 1e-3 and rel < 4e-4,
    )
    hist = {}
    for nm, t in (("rad", ta), ("bia", tb)):
        ch = rows[t]["polish_chunks"]
        hist[nm] = {
            "chunk_iters": [c["iters"] for c in ch],
            "chunk_E": [c["E_end"] for c in ch],
            "chunk_fmax": [c["fmax_end"] for c in ch],
            "drop_in_last_full_chunk": ch[-3]["E_end"] - ch[-2]["E_end"],
            "collect_last_chunk_drop_is_over_iters": ch[-1]["iters"],
        }
    ok = all(
        h["chunk_iters"][-1] > 50
        and h["drop_in_last_full_chunk"] < abs(res["rad"]["E"] - res["bia"]["E"])
        for h in hist.values()
    )
    check(
        "Q4.2",
        "the W1 pinned rows ARE converged (the gate was met by a settled descent)",
        "polish_chunks: how the gate was met, and the energy drop of the last full 500-iteration chunk against dE",
        {"history": hist, "dE_between_rows": abs(res["rad"]["E"] - res["bia"]["E"])},
        "gate met inside a chunk, last full-chunk drop below dE",
        ok,
    )
    own = rms12(A, B, mask)
    sa, sb = CORES.seed_core(cfg, "rad", DELTA), CORES.seed_core(cfg, "bia", DELTA)
    own_seed = rms12(sa, sb, mask)
    check(
        "Q4.3",
        "end RMS 0.044 over r < 12 (min over the 12 exact symmetries), seed RMS 0.093, ratio above 0.25",
        "own 12 images (axis permutations times the inversion)",
        {"end": own, "seed": own_seed, "ratio": own["rms_min12"] / own_seed["rms_min12"]},
        "0.0438, 0.0932",
        abs(own["rms_min12"] - 0.0438) < 2e-4 and abs(own_seed["rms_min12"] - 0.0932) < 2e-4,
    )
    # (a) channels of the difference in the local eigenframe of the mean field
    D = B - A
    lam, vec = np.linalg.eigh(0.5 * (A + B)[..., 1:, 1:])
    Df = np.einsum("...ia,...ij,...jb->...ab", vec, D[..., 1:, 1:], vec)
    tot = np.sum(D**2, axis=(-1, -2))
    chan = {
        "top_diag (V4-massive)": Df[..., 2, 2] ** 2,
        "pair_trace (V4-massive)": 0.5 * (Df[..., 0, 0] + Df[..., 1, 1]) ** 2,
        "M00 (V4-massive)": D[..., 0, 0] ** 2,
        "pair_split (V4-flat at the vacuum)": 0.5 * (Df[..., 1, 1] - Df[..., 0, 0]) ** 2,
        "pair_offdiag (V4-flat at the vacuum)": 2.0 * Df[..., 0, 1] ** 2,
        "top_axis_rotation (V4-flat)": 2.0 * (Df[..., 0, 2] ** 2 + Df[..., 1, 2] ** 2),
    }
    frac = {k: float(v[mask].sum() / tot[mask].sum()) for k, v in chan.items()}
    massive = sum(v for k, v in frac.items() if "massive" in k)
    gap = {
        f"r{R:g}": [
            float(x)
            for x in np.linalg.eigvalsh(0.5 * (A + B)[..., 1:, 1:])[np.abs(r - R) < 1.125].mean(0)
        ]
        for R in (3.0, 6.0, 9.0, 12.0)
    }
    check(
        "Q4.4",
        "the field difference lies along the V4 zero modes of the local vacuum frame",
        "cell-by-cell projection of D on the eigenframe of the mean field, fractions of sum |D|^2 over r < 12",
        {
            "fractions": frac,
            "massive_total": massive,
            "mean_field_eigenvalues_on_shells": gap,
            "note": "the W1 fields sit far from the (1, 0.3, 0.3) vacuum inside r 12, so 'flat at the vacuum' is nominal there",
        },
        "massive fraction under 0.05",
        massive < 0.05,
    )
    # (b) rotation-invariant reads
    la, lb = np.linalg.eigvalsh(A[..., 1:, 1:]), np.linalg.eigvalsh(B[..., 1:, 1:])
    prof = {}
    for R in (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0):
        sh = np.abs(r - R) < 1.125
        prof[f"{R:g}"] = {
            "rad_lam": [float(x) for x in la[sh].mean(0)],
            "bia_lam": [float(x) for x in lb[sh].mean(0)],
            "rad_biax": float((la[sh][:, 1] - la[sh][:, 0]).mean()),
            "bia_biax": float((lb[sh][:, 1] - lb[sh][:, 0]).mean()),
        }
    dmax = max(max(abs(a - b) for a, b in zip(v["rad_lam"], v["bia_lam"])) for v in prof.values())
    bmax = max(abs(v["rad_biax"] - v["bia_biax"]) for v in prof.values())
    inv_rms = float(np.sqrt(np.mean(np.sum((la - lb) ** 2, axis=-1)[mask])))
    check(
        "Q4.5",
        "the two endpoints carry the same rotation-invariant reads (a frame-only difference)",
        "sorted eigenvalue and biaxiality shell means; pointwise sorted-eigenvalue RMS over r < 12",
        {
            "max_shell_eigenvalue_diff": dmax,
            "max_shell_biaxiality_diff": bmax,
            "pointwise_eigenvalue_rms": inv_rms,
            "field_rms": own["rms_identity"],
            "profiles": prof,
        },
        "differences under 0.005",
        dmax < 0.005 and bmax < 0.005,
    )
    # line scan rad -> bia (both share the pinned shell, so the segment is admissible)
    scan = {}
    for t in (-0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25):
        e = R20.energy_parts(A + t * D, cfg, p, pot)
        scan[f"{t:g}"] = {"E": e["E_total"], "E_curv": e["E_curv"], "V": e["V"]}
    ts = np.array([float(k) for k in scan])
    es = np.array([v["E"] for v in scan.values()])
    c = np.polyfit(ts, es, 2)
    t_min = float(-c[1] / (2 * c[0]))
    barrier = float(max(es[(ts >= 0) & (ts <= 1)]) - max(es[ts == 0][0], es[ts == 1][0]))
    slope0 = float(2 * c[0] * 0.0 + c[1])
    check(
        "Q4.6",
        "rad and bia are two equal-energy minima (each stationary along the segment joining them)",
        "E on the straight segment A + t (B - A), t in [-0.25, 1.25], quadratic fit",
        {
            "scan": scan,
            "quadratic_min_at_t": t_min,
            "barrier": barrier,
            "dE_dt_at_rad": slope0,
            "curvature_d2E_dt2": float(2 * c[0]),
            "norm_D": float(np.sqrt(np.sum(D**2))),
        },
        "slope at t = 0 consistent with 0 and an interior maximum",
        abs(slope0) < 1e-4 and barrier > 0,
    )
    return {"A": A, "B": B, "mask": mask, "E": res}


# ================= Q5 =================
def q5(rows):
    import sympy as sp

    d, e = sp.symbols("delta epsilon", positive=True)
    lam = (1, d + e, d - e)
    Gs = [sp.expand(sp.prod([lam[a] - lam[b] for b in range(3) if b != a])) for a in range(3)]
    ok = sp.simplify(Gs[0] ** 2 - (((1 - d) ** 2 - e**2) ** 2)) == 0
    ok &= sp.simplify(Gs[1] ** 2 - 4 * e**2 * (1 - d - e) ** 2) == 0
    ok &= sp.simplify(Gs[2] ** 2 - 4 * e**2 * (1 - d + e) ** 2) == 0
    g = sp.symbols("g", positive=True)
    spec = (-g, 1, d + e, d - e)
    base = (-g, 1, d, d)
    V = sum((sum(q**k for q in spec) - sum(q**k for q in base)) ** 2 for k in range(1, 5))
    A_sym = sp.simplify(sp.series(sp.expand(V), e, 0, 5).removeO() / e**4)
    ok &= sp.simplify(A_sym - (4 + 36 * d**2 + 144 * d**4)) == 0
    w = W1 * 25.0
    A_num = 4 + 36 * DELTA**2 + 144 * DELTA**4
    c0 = float(np.sqrt(8 * (1 - DELTA) ** 2 / (w * A_num)))
    check(
        "Q5.1",
        "G_1^2 = ((1-delta)^2 - eps^2)^2, G_2,3^2 = 4 eps^2 (1-delta -/+ eps)^2, V4 = w A eps^4, c_0 = 5.08",
        "own sympy expansion",
        {"A_delta0.3": A_num, "c_0": c0},
        "identities hold; 5.08",
        bool(ok) and abs(c0 - 5.08) < 0.01,
    )
    # the quartic resistance with the three massive modes (top, pair mean, M_00) relaxed at fixed eps
    roots = (-G_T, 1.0, DELTA, DELTA)

    def vmin(eps):
        def res(x):
            t, m, s = x
            spec_ = (-G_T + s, 1.0 + t, DELTA + m + eps, DELTA + m - eps)
            return [sum(q**k for q in spec_) - sum(q**k for q in roots) for k in range(1, 5)]

        sol = least_squares(res, np.zeros(3), xtol=1e-15, ftol=1e-15, gtol=1e-15)
        return float(np.sum(np.asarray(sol.fun) ** 2)), sol.x

    A_eff = [float(vmin(x)[0] / x**4) for x in (0.02, 0.05, 0.1)]
    shifts = vmin(0.02)[1] / 0.02**2  # (top, pair mean, N_00) shifts per eps^2
    quad = 2.0 * (1 - DELTA) ** 2 - 4.0 * (1 - DELTA) ** 3 * (
        shifts[0] - shifts[1]
    )  # -d(G_1^2)/d(eps^2)
    c0_rel = float(np.sqrt(4.0 * quad / (w * A_eff[0]))) if quad > 0 else None
    check(
        "Q5.2",
        "the local V4 resistance to the split is w A eps^4 with A = 8.41 (delta 0.3)",
        "own V4 on one cell, minimized over the three V4-massive eigenvalue shifts at fixed eps",
        {
            "A_fixed_top_and_pair_mean": A_num,
            "A_relaxed": A_eff,
            "relaxed_shifts_per_eps2_(top, pair_mean, N00)": [float(x) for x in shifts],
            "G1_sq_loss_per_eps2_fixed": 2.0 * (1 - DELTA) ** 2,
            "G1_sq_loss_per_eps2_relaxed": float(quad),
            "c_0_fixed": c0,
            "c_0_same_balance_with_relaxed_V4": c0_rel,
        },
        "A_relaxed within 10 percent of A",
        abs(A_eff[0] / A_num - 1) < 0.1,
    )

    # (b) frame-current cost against the Coulomb gain on the end fields (central differences)
    cfg, p, pot = setup(32, 48.0, DELTA, 25.0)
    h = cfg["h"]
    X, Y, Z, r = radius(32, h)
    frame_rows = {}
    for tag in ("bia_pin_d0.3_w25_n32_L48_x9000", "rad_pin_d0.3_w25_n32_L48_x9000"):
        S = field(tag)[..., 1:, 1:]
        lam_, vec = np.linalg.eigh(S)
        top, mid, sm = lam_[..., 2], lam_[..., 1], lam_[..., 0]
        eps, pm = 0.5 * (mid - sm), 0.5 * (mid + sm)
        Gw = {
            2: ((top - mid) * (top - sm)) ** 2,
            1: ((mid - top) * (mid - sm)) ** 2,
            0: ((sm - top) * (sm - mid)) ** 2,
        }

        def rho2(T):
            dT = [np.gradient(T, h, axis=ax) for ax in range(3)]
            out = 0.0
            for i, j in ((0, 1), (0, 2), (1, 2)):
                F = dT[i] @ dT[j] - dT[j] @ dT[i]
                out = out + 0.5 * np.sum(F * F, axis=(-1, -2))
            return out

        rho = {a: rho2(vec[..., :, a, None] * vec[..., None, :, a]) for a in (0, 1, 2)}
        actual = 8.0 * rho2(S)  # 4 sum_{i<j} tr(F F^T) = 8 * (sum tr(F F^T) / 2)
        rec = {}
        for R in (6.0, 9.0, 12.0, 15.0):
            sh = np.abs(r - R) < 1.125
            gain = 8.0 * np.sum(((top - pm) ** 4 - Gw[2])[sh] * rho[2][sh])
            cost = 8.0 * np.sum((Gw[1] * rho[1] + Gw[0] * rho[0])[sh])
            framesum = 8.0 * np.sum(sum(Gw[a] * rho[a] for a in (0, 1, 2))[sh])
            rec[f"{R:g}"] = {
                "gain_top_axis": float(gain),
                "cost_tangent_currents": float(cost),
                "cost_over_gain": float(cost / gain),
                "frame_formula_total": float(framesum),
                "actual_curvature_total": float(np.sum(actual[sh])),
                "eps_mean": float(eps[sh].mean()),
                "rho1_r4_mean": float((rho[2] * r**4)[sh].mean()),
            }
        frame_rows[tag.split("_")[0] + "_9000"] = rec
    cg = [v["cost_over_gain"] for rec in frame_rows.values() for v in rec.values()]
    fr = [
        v["frame_formula_total"] / v["actual_curvature_total"]
        for rec in frame_rows.values()
        for v in rec.values()
    ]
    check(
        "Q5.3",
        "on the end fields the tangent-frame current cost stays below the top-axis gain, read with the frame formula "
        "e = 8 sum_a G_a^2 rho_a^2 that the mechanism rests on",
        "sign-free projector currents rho_a^2 = sum_{i<j} tr(F F^T)/2, F = [d_i P_a, d_j P_a], central differences, shell sums; "
        "the frame formula total against the actual curvature (same stencil) as the validity control",
        {
            "shells": frame_rows,
            "cost_over_gain_min_max": [min(cg), max(cg)],
            "frame_formula_over_actual_min_max": [min(fr), max(fr)],
            "note": "the fixed-eigenvalue frame formula misses the actual curvature by up to the factor shown (eigenvalue gradients, "
            "eigenvector exchanges), so this budget cannot decide the sign; the direct test Q5.4 does",
        },
        "cost / gain < 1 on every shell and frame formula within 25 percent of the actual curvature",
        max(cg) < 1.0 and 0.75 < min(fr) and max(fr) < 1.25,
    )

    # (c) the direct test on the exact uniaxial exterior
    M0 = CORES.seed_core(cfg, "rad", DELTA)
    e0 = R20.energy_parts(M0, cfg, p, pot)
    free = ~B3.pin_shell(32, h)
    rh = np.stack([X, Y, Z], axis=-1) / r[..., None]

    def frame_of(axis):
        ax = np.asarray(axis, float) / np.linalg.norm(axis)
        c = np.einsum("...a,a->...", rh, ax)
        s = np.sqrt(np.maximum(1.0 - c * c, 1e-300))
        ph = np.cross(ax, rh)
        ph = ph / np.maximum(np.linalg.norm(ph, axis=-1), 1e-300)[..., None]
        return np.cross(ph, rh), ph, s

    def g_shell(r0, sig):
        return lambda q: np.where(
            np.abs(q - r0) < sig, np.cos(np.pi * (q - r0) / (2 * sig)) ** 2, 0.0
        )

    def g_broad(q):
        taper = np.where(
            q < 15.0, 1.0, np.where(q < 21.0, np.cos(np.pi * (q - 15.0) / 12.0) ** 2, 0.0)
        )
        return np.sqrt(q / 10.0) * taper * (1.0 - np.exp(-((q / 5.0) ** 2)))

    trials = [
        ("z", (0, 0, 1), "shell r 13 +/- 6", g_shell(13.0, 6.0), 2, "split", (0.05, 0.1)),
        ("z", (0, 0, 1), "shell r 13 +/- 6", g_shell(13.0, 6.0), 2, "partner", (0.05,)),
        ("111", (1, 1, 1), "shell r 13 +/- 6", g_shell(13.0, 6.0), 2, "split", (0.05,)),
        ("z", (0, 0, 1), "shell r 10 +/- 4", g_shell(10.0, 4.0), 2, "split", (0.05,)),
        ("z", (0, 0, 1), "shell r 10 +/- 4", g_shell(10.0, 4.0), 4, "split", (0.05,)),
        ("z", (0, 0, 1), "shell r 15 +/- 5", g_shell(15.0, 5.0), 2, "split", (0.05,)),
        ("z", (0, 0, 1), "broad sqrt(r), r 4 to 21", g_broad, 2, "split", (0.05, 0.02)),
        ("z", (0, 0, 1), "broad sqrt(r), r 4 to 21", g_broad, 2, "partner", (0.05,)),
        ("x", (1, 0, 0), "broad sqrt(r), r 4 to 21", g_broad, 2, "split", (0.05,)),
    ]
    out = []
    for name, axis, shape, gfun, k, kind, amps in trials:
        th, ph, s = frame_of(axis)
        gr = gfun(r) * free
        f = gr * s**k
        if kind == "split":
            T = th[..., :, None] * th[..., None, :] - ph[..., :, None] * ph[..., None, :]
        else:
            T = th[..., :, None] * ph[..., None, :] + ph[..., :, None] * th[..., None, :]
        Dm = np.zeros_like(M0)
        Dm[..., 1:, 1:] = f[..., None, None] * T
        gp = (gfun(r + 1e-4) - gfun(r - 1e-4)) / 2e-4
        pred = {}
        if k == 2:
            for lam_t in (0.0, 4.0):
                dens = s**4 * (2.0 * gp**2 + (lam_t - 2.0) * gr**2 / r**2) / r**2
                pred[f"lam{lam_t:g}"] = float(8.0 * (1 - DELTA) ** 2 * np.sum(dens * free) * h**3)
        for a in amps:
            ep = R20.energy_parts(M0 + a * Dm, cfg, p, pot)
            em = R20.energy_parts(M0 - a * Dm, cfg, p, pot)
            out.append(
                {
                    "frame_axis": name,
                    "radial_shape": shape,
                    "angular": f"sin^{k}",
                    "tensor": kind,
                    "a": a,
                    "c2_curv": (ep["E_curv"] + em["E_curv"] - 2 * e0["E_curv"]) / (2 * a * a),
                    "c2_total": (ep["E_total"] + em["E_total"] - 2 * e0["E_total"]) / (2 * a * a),
                    "dV_mean": 0.5 * (ep["V"] + em["V"]) - e0["V"],
                    "continuum_c2_if_tangential_eigenvalue_0": pred.get("lam0"),
                    "continuum_c2_if_tangential_eigenvalue_4": pred.get("lam4"),
                    "no_gradient_prediction_of_1b": (
                        float(
                            -16.0
                            * (1 - DELTA) ** 2
                            * np.sum(f**2 * (1.0 - 2.0 * (1 - s * s) / (s * s)) / r**4)
                            * h**3
                        )
                        if kind == "split"
                        else None
                    ),
                }
            )
    neg = [o for o in out if o["c2_curv"] < 0]
    check(
        "Q5.4",
        "the uniaxial Coulomb tail is unstable to the split at every radius (d2E/da2 < 0 for a split perturbation)",
        "E_curv(+a) + E_curv(-a) - 2 E_curv(0) on the exact 'rad' seed exterior, production energy_parts; "
        "3 frames, 2 tensors, shell and box-filling radial shapes",
        {
            "trials": out,
            "negative_curvature_trials": len(neg),
            "c2_curv_min": min(o["c2_curv"] for o in out),
            "c2_curv_max": max(o["c2_curv"] for o in out),
            "note": "the second variation of a split eps(r) sin^2(theta) is 8 (1-delta)^2 / r^2 [2 eps_r^2 + (lam - 2) eps^2 / r^2]; "
            "lam = 4 reproduces the clean shells (r 10 +/- 4, r 15 +/- 5) to 0.3 percent, so the form is positive for "
            "every radial profile r^s (2 s^2 + 2 > 0); the -16 (1-delta)^2 eps^2 / r^4 of the audited note is there, "
            "outweighed by gradient terms of the same order that it omits",
        },
        "at least one trial with c2_curv < 0",
        len(neg) > 0,
    )

    # (c2) every split texture at once: the lowest Hessian eigenvalue in the split sector on a shell
    region = free & (r > 8.0) & (r < 16.0)
    nc = int(region.sum())
    th, ph, _ = frame_of((0, 0, 1))
    th, ph = th[region], ph[region]
    B1 = (th[:, :, None] * th[:, None, :] - ph[:, :, None] * ph[:, None, :]) / np.sqrt(2.0)
    B2 = (th[:, :, None] * ph[:, None, :] + ph[:, :, None] * th[:, None, :]) / np.sqrt(2.0)

    def shifted(v, t):
        c = v.reshape(-1, 2)
        Mx = M0.copy()
        sub = Mx[region]
        sub[:, 1:, 1:] += t * (c[:, 0, None, None] * B1 + c[:, 1, None, None] * B2)
        Mx[region] = sub
        return Mx

    def gproj(Mx):
        Gs = R20.energy_grad(Mx, cfg, p, pot)[1][region][:, 1:, 1:]
        return np.stack(
            [np.sum(Gs * B1, axis=(1, 2)), np.sum(Gs * B2, axis=(1, 2))], axis=1
        ).ravel()

    def hv(v):
        t = 2e-3 / np.abs(v).max()
        return (gproj(shifted(v, t)) - gproj(shifted(v, -t))) / (2.0 * t)

    rng = np.random.default_rng(1)
    qv = rng.normal(size=2 * nc)
    Qs, al, be = [qv / np.linalg.norm(qv)], [], []
    for _ in range(80):
        wv = hv(Qs[-1])
        al.append(float(wv @ Qs[-1]))
        for qq in Qs:
            wv = wv - (wv @ qq) * qq
        be.append(float(np.linalg.norm(wv)))
        Qs.append(wv / be[-1])
    Tm = np.diag(al) + np.diag(be[:-1], 1) + np.diag(be[:-1], -1)
    ev, U = np.linalg.eigh(Tm)
    yv = np.array(Qs[:-1]).T @ U[:, 0]
    ep, em = R20.energy_parts(shifted(yv, 0.05), cfg, p, pot), R20.energy_parts(
        shifted(yv, -0.05), cfg, p, pot
    )
    d2 = (ep["E_total"] + em["E_total"] - 2.0 * e0["E_total"]) / 0.05**2
    wgt = np.sum(yv.reshape(-1, 2) ** 2, axis=1)
    lz = {
        "cells": nc,
        "lanczos_steps": 80,
        "lowest_ritz_values": [float(x) for x in ev[:3]],
        "largest": float(ev[-1]),
        "residual_of_lowest": float(abs(be[-1] * U[-1, 0])),
        "direct_d2E_da2_along_the_ritz_vector": float(d2),
        "mean_r_of_the_mode": float(np.sum(wgt * r[region]) / wgt.sum()),
        "local_value_claimed_by_1b_at_r9_r12_r15": [
            -16 * (1 - DELTA) ** 2 * h**3 / R**4 for R in (9.0, 12.0, 15.0)
        ],
        "note": "finite-difference Hessian-vector products of the production gradient (step 2e-3 in sup norm; a step of 0.5 "
        "in 2-norm gives spurious negative Ritz values through the quartic nonlinearity, caught by the direct check)",
    }
    check(
        "Q5.8",
        "on (certified quartic + V4) the uniaxial exterior is a saddle in the split sector (some split texture lowers E)",
        "Lanczos (full reorthogonalization) on the Hessian restricted to the 2 split coordinates per cell of the shell 8 < r < 16, "
        "'rad' seed, all other entries clamped; the Ritz vector re-tested by a direct energy second difference",
        lz,
        "lowest eigenvalue < 0",
        ev[0] < 0.0 and d2 < 0.0,
    )

    # (d) the radial law of the halo
    halo = {}
    for nm, tag in (
        ("rad_3000", "rad_pin_d0.3_w25_n32_L48"),
        ("rad_9000", "rad_pin_d0.3_w25_n32_L48_x9000"),
        ("bia_3000", "bia_pin_d0.3_w25_n32_L48"),
        ("bia_9000", "bia_pin_d0.3_w25_n32_L48_x9000"),
        ("perm_9000", "perm_pin_d0.3_w25_n32_L48_x9000"),
        ("obl_9000", "obl_pin_d0.3_w25_n32_L48_x9000"),
    ):
        lam_ = np.linalg.eigvalsh(field(tag)[..., 1:, 1:])
        eps = 0.5 * (lam_[..., 1] - lam_[..., 0])
        Rs = np.arange(6.0, 18.1, 1.5)
        em = np.array([eps[np.abs(r - R) < 0.75].mean() for R in Rs])
        sl_in = np.polyfit(np.log(Rs[:5]), np.log(em[:5]), 1)[0]
        sl_out = np.polyfit(np.log(Rs[4:]), np.log(em[4:]), 1)[0]
        halo[nm] = {
            "r": [float(x) for x in Rs],
            "eps": [float(x) for x in em],
            "eps_r2": [float(x) for x in em * Rs**2],
            "loglog_slope_r6_to_12": float(sl_in),
            "loglog_slope_r12_to_18": float(sl_out),
        }
    sl = [v["loglog_slope_r6_to_12"] for k_, v in halo.items() if "9000" in k_]
    ratio = [max(v["eps_r2"]) / c0 for k_, v in halo.items() if "9000" in k_]
    check(
        "Q5.5",
        "the halo is the local-balance power law eps = c_0 / r^2 with c_0 = 5.08",
        "shell means of eps every 1.5 from r 6 to 18, log-log slopes, eps r^2 against c_0",
        {"halo": halo, "slopes_9000_r6_to_12": sl, "max_eps_r2_over_c0_at_9000": ratio},
        "slope near -2 and eps r^2 near 5.08",
        all(abs(s_ + 2.0) < 0.3 for s_ in sl) and max(ratio) < 1.5,
    )

    # (e) the running weight, and what carries it
    run = {}
    for nm, tag in (
        ("rad", "rad_pin_d0.3_w25_n32_L48_x9000"),
        ("bia", "bia_pin_d0.3_w25_n32_L48_x9000"),
        ("perm", "perm_pin_d0.3_w25_n32_L48_x9000"),
        ("obl", "obl_pin_d0.3_w25_n32_L48_x9000"),
    ):
        lam_ = np.linalg.eigvalsh(field(tag)[..., 1:, 1:])
        sh = np.abs(r - 9.0) < 1.125
        top, mid, sm = lam_[..., 2][sh], lam_[..., 1][sh], lam_[..., 0][sh]
        g1 = float((((top - mid) * (top - sm)) ** 2).mean())
        uni = float(((top - 0.5 * (mid + sm)) ** 4).mean())
        run[nm] = {
            "G1_sq_r9": g1,
            "uniaxial_part_(top - pair mean)^4": uni,
            "share_of_the_drop_from_the_split": (uni - g1) / (0.7**4 - g1),
        }
    vals = [v["G1_sq_r9"] for v in run.values()]
    check(
        "Q5.6",
        "G_1^2 is 0.12 to 0.15 at r 9 against the uniaxial 0.2401",
        "shell mean of ((top - mid)(top - small))^2 at r 9",
        run,
        "0.12 to 0.155",
        0.12 <= min(vals) and max(vals) <= 0.155,
    )
    shares = [v["share_of_the_drop_from_the_split"] for v in run.values()]
    check(
        "Q5.7",
        "the charge weight runs BECAUSE of the split (the split carries most of the drop of G_1^2 at r 9)",
        "G_1^2 = (top - pair mean)^4 - 2 (top - pair mean)^2 eps^2 + eps^4, shell means",
        {"split_share": shares},
        "share above 0.5",
        min(shares) > 0.5,
    )
    return {"c0": c0}


# ================= Q6, Q7, Q8 =================
def q6(table):
    rel = {}
    for s in ("rad", "bia", "obl", "perm", "big"):
        a, b = (
            table[f"{s}_pin_d0.3_w25_n32_L48"]["E_total"],
            table[f"{s}_free_d0.3_w25_n32_L48"]["E_total"],
        )
        rel[s] = (b - a) / a
    check(
        "Q6.1",
        "pin and free agree within 1 percent at 3000 iterations for every seed at W1 x 25",
        "recomputed E_total, (free - pin) / pin",
        rel,
        "all |rel| < 0.01",
        all(abs(v) < 0.01 for v in rel.values()),
    )
    return rel


def q7(rows):
    tag = "rad_pin_d0.89_w25_n32_L48"
    cfg, p, pot = setup(32, 48.0, 0.89, 25.0)
    fm = R21.free_mask(cfg, True)
    M, Ms = field(tag), CORES.seed_core(cfg, "rad", 0.89)

    def fmax(F, c, pp, po):
        _, G, _ = R20.energy_grad(F, c, pp, po)
        Ga = np.abs(G) * fm[..., None, None]
        Gs = Ga.copy()
        Gs[..., 0, 0] = 0.0
        return float(Ga.max()), float(Gs.max())

    end, seed = fmax(M, cfg, p, pot), fmax(Ms, cfg, p, pot)
    cfg3, p3, pot3 = setup(32, 48.0, DELTA, 25.0)
    seed3 = fmax(CORES.seed_core(cfg3, "rad", DELTA), cfg3, p3, pot3)
    X, Y, Z, r = radius(32, cfg["h"])
    moved = float(np.sqrt(np.mean(np.sum((M - Ms) ** 2, axis=(-1, -2))[r < 12.0])))
    scale = ((1 - 0.89) / (1 - DELTA)) ** 4
    val = {
        "polish_iters": rows[tag]["polish_iters"],
        "fmax_end_all_and_spatial": end,
        "fmax_seed_all_and_spatial": seed,
        "fmax_seed_delta0.3_all_and_spatial": seed3,
        "(1-delta)^4": 0.11**4,
        "gradient_scale_ratio_to_delta0.3": scale,
        "gate_scaled_by_(1-delta)^4": 1e-3 * scale,
        "end_spatial_fmax_over_scaled_gate": end[1] / (1e-3 * scale),
        "rms_moved_from_seed_r_lt_12": moved,
        "E_seed": rows[tag]["seed_E"]["E_total"],
        "E_end": rows[tag]["E"],
    }
    check(
        "Q7.1",
        "the delta 0.89 row is uninformative: the absolute 1e-3 gate is met trivially at gradient scale (1-delta)^4",
        "production gradient on the seed and on the end field; the gate rescaled by ((1-0.89)/(1-0.3))^4",
        val,
        "end residual far above the rescaled gate",
        end[1] > 10 * 1e-3 * scale,
    )


def q8(rows):
    out = {}
    for tag, r in rows.items():
        if r["n"] != 32 or r["bnd"] != "pin" or r["delta"] != 0.3:
            continue
        rec = {"json_flux_degree_r9": r["end_shells"]["degree_top_r9"]}
        rec.update(degree_reads(field(tag)))
        out[tag] = rec
    ok = all(
        v["r9"]["degree"] == 1 and all(c == 1 for c in v["r9"]["probe_counts"])
        for v in out.values()
    )
    check(
        "Q8.1",
        "the one-charge sector is kept: the degree of the top eigenvector through the r 9 cube is 1",
        "own reader: surface-oriented eigenvector, solid-angle sum on the triangulated cell cube 10..21, signed probe preimages",
        out,
        "degree 1, zero conflicting links, on every pinned n 32 row",
        ok,
    )
    bad6 = {
        t: {
            "conflicting_links": v["r6"]["conflicting_links"],
            "min_gap_top_mid": v["r6"]["min_gap_top_mid_on_surface"],
        }
        for t, v in out.items()
        if v["r6"]["degree"] != 1
    }
    flux = {t: v["json_flux_degree_r9"] for t, v in out.items() if t.endswith("_x9000")}
    check(
        "Q8.2",
        "the charge stays in a compact core: the top eigenvector is orientable with degree 1 on the r 6 cube too",
        "the same reader on the cell cube 12..19 (half-width 5.25)",
        {
            "rows_not_degree_1_at_r6": bad6,
            "note": "the 0.93 to 0.98 of the audited flux reader at r 9 is reader error: the solid-angle degree there is exactly 1",
            "audited_flux_reads_r9_at_9000": flux,
        },
        "degree 1 on every row",
        not bad6,
    )
    return out


def q_n48(rows):
    ta, tb = "rad_pin_d0.3_w25_n48_L48", "bia_pin_d0.3_w25_n48_L48"
    have = [t for t in (ta, tb) if t in rows]
    val = {"final_fields_present": have}
    for t in have:
        ch = rows[t]["polish_chunks"]
        val[t] = {
            "E": rows[t]["E"],
            "polish_iters": rows[t]["polish_iters"],
            "label": rows[t]["label"],
            "last_chunk_drop": ch[-2]["E_end"] - ch[-1]["E_end"] if len(ch) > 1 else None,
            "E_n32_same_seed_at_3000": rows[t.replace("n48", "n32")]["E"],
        }
    ok = False
    if len(have) == 2:
        cfg, _, _ = setup(48, 48.0, DELTA, 25.0)
        X, Y, Z, r = radius(48, cfg["h"])
        end = rms12(field(ta), field(tb), r < 12.0)
        sd = rms12(
            CORES.seed_core(cfg, "rad", DELTA), CORES.seed_core(cfg, "bia", DELTA), r < 12.0
        )
        tol = max(
            0.01 * abs(val[ta]["E"]),
            3 * max(val[ta]["last_chunk_drop"], val[tb]["last_chunk_drop"]),
        )
        val.update(
            {
                "rms_end": end["rms_min12"],
                "rms_seed": sd["rms_min12"],
                "dE": val[ta]["E"] - val[tb]["E"],
                "tol_E": tol,
            }
        )
        ok = abs(val["dE"]) < tol and end["rms_min12"] < 0.25 * sd["rms_min12"]
    converged = all("FALLING" not in str(val[t]["label"]) for t in have) and len(have) == 2
    check(
        "Q9.2",
        "the n 48 refinement rows can decide BRANCH_LATTICE_ONLY (both present, both at the gate, the rule applied)",
        "the pre-registered rule with own 12 images, if both final fields exist",
        val,
        "both rows present and converged",
        converged and ok,
    )


# ================= the polish pool and its reads =================
def run_polish_pool():
    jobs = [
        {"tag": "rad_pin_d0.3_w1_n32_L48", "w1s": 1.0, "mode": "reduced", "iters": 300},
        {"tag": "bia_pin_d0.3_w1_n32_L48", "w1s": 1.0, "mode": "reduced", "iters": 300},
        {"tag": "bia_pin_d0.3_w25_n32_L48_x9000", "w1s": 25.0, "mode": "reduced", "iters": 150},
        {"tag": "rad_pin_d0.3_w25_n32_L48_x9000", "w1s": 25.0, "mode": "reduced", "iters": 150},
        {"tag": "bia_pin_d0.3_w25_n32_L48_x9000", "w1s": 25.0, "mode": "plain", "iters": 40},
    ]
    with ProcessPoolExecutor(max_workers=2, mp_context=mp.get_context("spawn")) as ex:
        res = list(ex.map(polish_job, jobs))
    return {(r["job"]["tag"], r["job"]["mode"]): r for r in res}


def q_polish(rows, q4, pol, fits):
    cfg, p, pot = setup(32, 48.0, DELTA, 1.0)
    A, B, mask = q4["A"], q4["B"], q4["mask"]
    ra, rb = (
        pol[("rad_pin_d0.3_w1_n32_L48", "reduced")],
        pol[("bia_pin_d0.3_w1_n32_L48", "reduced")],
    )
    before, after = rms12(A, B, mask), rms12(ra["M"], rb["M"], mask)
    la, lb = np.linalg.eigvalsh(ra["M"][..., 1:, 1:]), np.linalg.eigvalsh(rb["M"][..., 1:, 1:])
    val = {
        "rad": {k: ra[k] for k in ("E_in", "E_end", "iters", "n_eval", "fmax_end_spatial")},
        "bia": {k: rb[k] for k in ("E_in", "E_end", "iters", "n_eval", "fmax_end_spatial")},
        "dE_before": ra["E_in"] - rb["E_in"],
        "dE_after": ra["E_end"] - rb["E_end"],
        "rms_before": before["rms_min12"],
        "rms_after": after["rms_min12"],
        "rms_rad_moved": rms12(A, ra["M"], mask)["rms_identity"],
        "rms_bia_moved": rms12(B, rb["M"], mask)["rms_identity"],
        "eigenvalue_rms_after": float(np.sqrt(np.mean(np.sum((la - lb) ** 2, axis=-1)[mask]))),
        "E_trace_rad_every_50": [t[2] for t in ra["trace"][49::50]],
        "E_trace_bia_every_50": [t[2] for t in rb["trace"][49::50]],
        "degree_reads_rad_after": {
            k: (v["degree"], v["conflicting_links"]) for k, v in degree_reads(ra["M"]).items()
        },
        "degree_reads_bia_after": {
            k: (v["degree"], v["conflicting_links"]) for k, v in degree_reads(rb["M"]).items()
        },
    }
    gate_dE = abs(val["dE_before"])
    still = max(ra["E_in"] - ra["E_end"], rb["E_in"] - rb["E_end"])
    check(
        "Q4.7",
        "a further polish leaves the two W1 fields where they are (converged endpoints, apart at one energy)",
        "300 L-BFGS iterations on copies of each field, pinned shell held, M_00 slaved (reduced polish), production energy_grad",
        val,
        "each field's further energy drop below the 7.8e-4 gap between them",
        still < gate_dE,
    )
    closer = val["rms_after"] < 0.5 * val["rms_before"]
    check(
        "Q4.8",
        "decision: same branch, two distinct minima, or undecidable at this gate",
        "Q4.2, Q4.5, Q4.6, Q4.7 together",
        {
            "rms_ratio_after_over_before": val["rms_after"] / val["rms_before"],
            "closer_by_half": closer,
            "reading": "neither row was converged at the gate; see the final report for the decision",
        },
        "informational (PASS iff the claim 'two converged equal-energy endpoints' survives)",
        still < gate_dE and not closer,
    )

    # W1 x 25: the instrument
    pb, pr = (
        pol[("bia_pin_d0.3_w25_n32_L48_x9000", "reduced")],
        pol[("rad_pin_d0.3_w25_n32_L48_x9000", "reduced")],
    )
    pp = pol[("bia_pin_d0.3_w25_n32_L48_x9000", "plain")]
    n_pl = pp["iters"]
    prod_rate = {
        s: (
            rows[f"{s}_pin_d0.3_w25_n32_L48_x9000"]["polish_chunks"][-2]["E_end"]
            - rows[f"{s}_pin_d0.3_w25_n32_L48_x9000"]["polish_chunks"][-1]["E_end"]
        )
        / 500.0
        for s in ("bia", "rad")
    }
    val = {
        "production_drop_per_iteration_last_chunk": prod_rate,
        "plain_restart": {
            "iters": n_pl,
            "E_in": pp["E_in"],
            "E_end": pp["E_end"],
            "drop_per_it": (pp["E_in"] - pp["E_end"]) / n_pl,
            "fmax_all_in": pp["in_E_fmax_fmaxsp"][1],
            "fmax_all_after_1_it": pp["trace"][0][3],
            "fmax_spatial_after_1_it": pp["trace"][0][4],
            "fmax_all_end": pp["fmax_end_all"],
        },
        "reduced_bia": {
            "iters": pb["iters"],
            "E_in": pb["E_in"],
            "E_end": pb["E_end"],
            "drop_per_it": (pb["E_in"] - pb["E_end"]) / pb["iters"],
            "E_at_it_40": pb["trace"][min(39, len(pb["trace"]) - 1)][2],
            "E_every_25": [t[2] for t in pb["trace"][24::25]],
            "fmax_end_spatial": pb["fmax_end_spatial"],
        },
        "reduced_rad": {
            "iters": pr["iters"],
            "E_in": pr["E_in"],
            "E_end": pr["E_end"],
            "drop_per_it": (pr["E_in"] - pr["E_end"]) / pr["iters"],
            "E_every_25": [t[2] for t in pr["trace"][24::25]],
            "fmax_end_spatial": pr["fmax_end_spatial"],
        },
        "extrapolated_floor_window_bia": [fits["bia"]["E_inf_min"], fits["bia"]["E_inf_max"]],
        "extrapolated_floor_window_rad": [fits["rad"]["E_inf_min"], fits["rad"]["E_inf_max"]],
    }
    speed = val["reduced_bia"]["drop_per_it"] / prod_rate["bia"]
    val["speedup_reduced_over_production_bia"] = speed
    check(
        "Q2.3",
        "the slow fall at W1 x 25 is the landscape's (the 10-entry L-BFGS polish is an adequate instrument there)",
        "the same L-BFGS with M_00 slaved per cell to its V4 minimizer, from the 9000-iteration fields, against the production rate",
        val,
        "energy drop per iteration within 3 x of the production polish",
        speed < 3.0,
    )
    check(
        "Q3.4",
        "max |G| falls from its end-of-chunk spike to 0.02 to 0.07 after ONE iteration of a fresh L-BFGS",
        "plain 10-entry L-BFGS restarted from the bia 9000 field",
        {"fmax_in": pp["in_E_fmax_fmaxsp"][1], "fmax_after_1_it": pp["trace"][0][3]},
        "in: 1.2; after one iteration: under 0.1",
        pp["in_E_fmax_fmaxsp"][1] > 1.0 and pp["trace"][0][3] < 0.1,
    )
    # where the reduced polish takes the two W1 x 25 fields
    cfg25, p25, pot25 = setup(32, 48.0, DELTA, 25.0)
    X, Y, Z, r = radius(32, cfg25["h"])
    reads = {}
    for nm, res_ in (("bia", pb), ("rad", pr)):
        lam, vec = np.linalg.eigh(res_["M"][..., 1:, 1:])
        eps = 0.5 * (lam[..., 1] - lam[..., 0])
        parts = R20.energy_parts(res_["M"], cfg25, p25, pot25)
        reads[nm] = {
            "E_curv": parts["E_curv"],
            "V": parts["V"],
            "eps_r2": {
                f"{R:g}": float(eps[np.abs(r - R) < 1.125].mean() * R * R)
                for R in (6.0, 9.0, 12.0, 15.0, 18.0)
            },
            "degree_reads": {
                k: (v["degree"], v["conflicting_links"], v["min_gap_top_mid_on_surface"])
                for k, v in degree_reads(res_["M"]).items()
            },
            "pin_bitwise": bool(
                np.array_equal(
                    res_["M"][~R21.free_mask(cfg25, True)],
                    field(res_["job"]["tag"])[~R21.free_mask(cfg25, True)],
                )
            ),
        }
    reads["rms_bia_vs_rad_before"] = rms12(
        field("bia_pin_d0.3_w25_n32_L48_x9000"), field("rad_pin_d0.3_w25_n32_L48_x9000"), r < 12.0
    )
    reads["rms_bia_vs_rad_after"] = rms12(pb["M"], pr["M"], r < 12.0)
    ordered = pb["E_end"] < pr["E_end"]
    check(
        "Q9.1",
        "one descent path ordered by how biaxial the seed's core was; no compact endpoint within 9000 iterations",
        "reads of the two reduced-polish continuations (150 iterations each)",
        reads,
        "informational: PASS iff the ordering survives and the r 9 cube keeps degree 1 on both",
        ordered
        and reads["bia"]["degree_reads"]["r9"][0] == 1
        and reads["rad"]["degree_reads"]["r9"][0] == 1,
    )
    return val, reads


def write_json(rows, pending, partial):
    n_pass = sum(c["verdict"] == "PASS" for c in CHECKS)
    out = {
        "task": "M5.32 R22-1 adversarial audit",
        "partial_static_only": partial,
        "rows_audited": sorted(rows),
        "rows_still_computing": pending,
        "checks": CHECKS,
        "totals": {"PASS": n_pass, "FAIL": len(CHECKS) - n_pass, "checks": len(CHECKS)},
        "wall_s": round(time.time() - T0, 1),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(
        f"TOTAL PASS {n_pass} FAIL {len(CHECKS) - n_pass} of {len(CHECKS)}; wall {out['wall_s']} s"
    )


def main():
    rows = all_rows()
    log(f"{len(rows)} rows with a final field")
    pending = [
        t
        for t in (
            "bia_pin_d0.3_w25_n48_L48",
            "rad_pin_d0.3_w25_n48_L48",
            "rad_pin_d0.3_w25_n64_L96",
        )
        if not os.path.exists(os.path.join(NPZ, t + ".npz"))
    ]
    table = q1(rows)
    fits = q2(rows)
    q3(rows)
    q4 = q4_static(rows)
    q5(rows)
    q6(table)
    q7(rows)
    q8(rows)
    if len(sys.argv) > 1 and sys.argv[1] == "static":
        log("static mode: the polish pool is skipped")
        write_json(rows, pending, partial=True)
        return
    q_n48(rows)
    log("static checks done; the polish pool (2 workers) starts")
    pol = run_polish_pool()
    log("polish pool done")
    q_polish(rows, q4, pol, fits)
    write_json(rows, pending, partial=False)


if __name__ == "__main__":
    main()
