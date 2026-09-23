"""M5.32 R25-2: a charge's strand on the certified biaxial vacuum, the box ladder
of the compact electron (the first OpenWave measurement of strand confinement).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. Code branch
s = -1: M_vac = diag(8, 1, delta, 0), the certified BIAXIAL vacuum (the roots
R0.roots_of(cfg), not the degenerate pair of R22 to R24). The static energy is
the R20 instrument's (m5_32_r20_1_axes.energy_grad, consumed read-only):
    E[M] = 4 h^3 sum_br wt sum_cells sum_{i<j} <F_ij, F_ij>_eta + V4[M],
    F_ij = [A_i, A_j]_eta, A_i = d_i M (sym stencil),
    V4 = w sum_{p<=4} (tr N^p - C_p)^2, C_p = sum_i q_i^p, w = W1 x 25.
The object: the director hedgehog with the transverse pair (delta, 0) laid in
the (phi-hat, theta-hat) frame, R20.seed_axes(cfg, (1, delta, 0)). That frame
is singular on the whole polar axis, so the seed carries the strand as two
half-lines from the core to the walls, each a full turn of the transverse line
field (2 half-turns); the far sphere carries transverse index 2, the Euler
number of n*TS^2 for a unit charge.
The strand's tension per half-turn in the block sector is the Bogomolny value
(the author's 2026-09-22 19:24 UTC reply, reproduced at PLAN by our own script)
    T_half = pi sqrt(32 w K) (delta / 2)^4,  K = 4 + 36 s0^2 + 144 s0^4,  s0 = delta / 2
(2.675e-3 at delta 0.3, W1 x 25), additive in the index. A centered charge whose
two strands end on the faces gains dl = dL / 2 per strand under a box change dL,
each strand carrying 2 half-turns, so the predicted slope of the box ladder is
    dE / dL = 2 T_half        (the general form T_half sum_s m_s dl_s / dL
                               is also computed from the read windings).

THE INSTRUMENT
--------------
The R23-1 reduced descent at c = 0, consumed read-only (m5_32_r23_1_cscan.Reduced:
the six spatial entries per free cell with M_00 slaved per cell, scipy L-BFGS-B
in chunks of CS.CHUNK = 250, memory 20), the pinned shell (depth 1.6), resumable
stage files. The gate is a decade under R23's: fmax_spatial < 1e-4 read inside a
chunk and the last-chunk drop under 1e-5 abs(E), because the slope signal is
5e-3 per box unit. After the gate, the kick-and-continue clause of R25-K: 0.02
Gaussian on the six free spatial entries, M_00 re-slaved, two more chunks; the
quoted E is the kicked end if it is lower by over 1e-4 (kick_label SADDLE), else
the gate energy (STABLE); a row that never reaches the gate is FALLING.

THE ROWS
--------
run:          delta 0.3, w1s 25: h 1.5 at (n 32, L 48), (n 48, L 72), (n 64, L 96);
              h 1.0 at (n 48, L 48), (n 64, L 64)
run_stretch:  delta 0.1 at (n 32, L 48) and (n 48, L 72), the iteration cap x 3
run_pair:     a stretch arm behind the puncture gate of R25-0 f (not run by
              default): the +- pair at d in {9, 12, 15} on n 48 L 72, the director
              texture n = stereographic inverse of w_1 conj(w_2), w_k = (x + i y) /
              (r_k + z_k) about z = +-d/2, the transverse frame pulled back from
              the plane (e1 = normalize(dn / dRe zeta)), cores melted to the
              isotropic mean over r_k / 4

THE READS (own code below)
--------------------------
Energy parts; the winding of the transverse LINE field (the middle eigenvector
of the spatial block) on circles of radius 2 h around the z, x and y axes at
offsets +-6, +-9, ... up to L/2 - 3 h: a loop is read only where the top
eigenvector has abs(n_axis) > 0.95 on every sample, the winding is the sum of
the angle increments mod pi over pi (2 = one full turn = index 1); the tube read
T_read(z) = the density summed over rho < 6 on each z plane with abs(z) > 9,
per unit length; the spin-gate numbers on the rigid clock about the strand axis
(B3.gen_catalog(cfg, M)["rot_z"], the M5.21.3 envelope renv 10): C = kin,
omega_* = sqrt(E / (3 C)), J_* = sqrt(4 C E / 3), the gate 2 J omega / E = 1 on
E_J = E + J^2 / (4 C), E_rot / E = 1 / 4 by construction.

PRE-REGISTERED LABELS
---------------------
    CHARGE_STRAND_MATCHES   the fitted dE/dL within 25 percent of the prediction on
                            BOTH spacings, every row STABLE or SADDLE-then-reconverged
    CHARGE_STRAND_DIFFERS   outside 25 percent on both spacings with the rows at the gate
    INSUFFICIENT            otherwise
    SPIN_GATE_BOX_DEPENDENT omega_*(L 96) / omega_*(L 48) at h 1.5 off 1 by over 10 percent
    SPIN_GATE_PINNED        within 10 percent on both spacings; else UNRESOLVED
    CHARGE_PAIR_LINEAR      (pair arm) the two consecutive slopes agree within 25 percent
                            and match 4 T_half within 25 percent
    CHARGE_PAIR_COULOMB     their ratio within 25 percent of (9 x 12) / (12 x 15) = 0.6;
                            else UNRESOLVED
The delta 0.1 stretch reports sigma(0.1) / sigma(0.3) against (1 / 81) K(0.05) / K(0.15).

Modes: smoke | run [workers] | run_stretch [workers] | run_pair [workers] |
collect | plot | jobs. Output: data/m5_32_r25_2_charge.json (rows + collect),
arrays in data/m5_32_r25_2/ (local), plots/m5_32_r25_2_charge.png.
Regenerate: run 12 about 3 to 4 h (the n 64 L 96 row the long pole);
run_stretch about the same; collect and plot seconds; smoke about a minute.
"""

import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r25_2_charge.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r25_2")
R25_1_JSON = os.path.join(DATA, "m5_32_r25_1_strand.json")
PLOT_PNG = os.path.join(PLOTS, "m5_32_r25_2_charge.png")

G = 8.0
W1S = 25.0
GATE, GATE_DROP = 1e-4, 1e-5
MAX_ITER = {16: 4, 32: 8000, 48: 8000, 64: 6000}
STRETCH_FACTOR = 3
KICK_AMP = 0.02
KICK_CHUNKS = 2
KICK_LOWER = 1e-4
LOOP_RADIUS_H = 2.0
LOOP_OFFSETS = (6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0)
LOOP_SAMPLES = 96
AXIAL_TOL = 0.95
TUBE_RHO, TUBE_ZMIN = 6.0, 9.0
TOL_SLOPE, TOL_SPIN = 0.25, 0.10
PAIR_DS = (9.0, 12.0, 15.0)
EXPECTED_L = {"1.5": (48.0, 72.0, 96.0), "1": (48.0, 64.0)}
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, B3, R0, W1 = CS.R21, CS.R20, CS.B3, CS.R0, CS.W1
IU3 = CS.IU3


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the Bogomolny tension and the prediction =================
def K_of(s0):
    return 4.0 + 36.0 * s0**2 + 144.0 * s0**4


def T_half_bps(delta, w):
    """pi sqrt(32 w K) b0^4 with b0 = s0 = delta / 2 (the block sector's exact minimum)."""
    s0 = delta / 2.0
    return float(np.pi * np.sqrt(32.0 * w * K_of(s0)) * s0**4)


def t_half_of(delta, w1s, h):
    """T_half from the R25-1 slab rows at the same delta, w1s, h if present, else the Bogomolny value."""
    fallback = T_half_bps(delta, W1 * w1s)
    if not os.path.exists(R25_1_JSON):
        return fallback, "bogomolny (R25-1 JSON absent)"
    try:
        with open(R25_1_JSON) as f:
            J = json.load(f)
        cands = []
        for tag, r in J.get("rows", {}).items():
            if not isinstance(r, dict) or r.get("T") is None:
                continue
            if abs(float(r.get("delta", -1)) - delta) > 1e-12:
                continue
            if abs(float(r.get("w1s", -1)) - w1s) > 1e-12:
                continue
            rh = r.get("h")
            if rh is None and r.get("n") and r.get("L"):
                rh = float(r["L"]) / float(r["n"])
            if rh is None or abs(float(rh) - h) > 1e-9:
                continue
            cands.append((float(r["T"]), tag))
        if cands:
            cands.sort()
            return cands[0][0], f"R25-1 row {cands[0][1]}"
    except Exception as e:  # noqa: BLE001
        return fallback, f"bogomolny (R25-1 JSON unreadable: {e!r})"
    return fallback, "bogomolny (no R25-1 row at this delta, w1s, h)"


# ================= jobs =================
def job_tag(j):
    base = f"{j['kind']}_d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}"
    if j["kind"] == "pair":
        base += f"_dd{j['d']:g}"
    return base


def jobs_main():
    rows = [(32, 48.0), (48, 72.0), (64, 96.0), (48, 48.0), (64, 64.0)]
    return [dict(kind="S1", delta=0.3, w1s=W1S, n=n, L=L, cap=MAX_ITER[n]) for n, L in rows]


def jobs_stretch():
    rows = [(32, 48.0), (48, 72.0)]
    return [
        dict(kind="S1", delta=0.1, w1s=W1S, n=n, L=L, cap=STRETCH_FACTOR * MAX_ITER[n])
        for n, L in rows
    ]


def jobs_pair():
    return [
        dict(kind="pair", delta=0.3, w1s=W1S, n=48, L=72.0, d=d, cap=MAX_ITER[48]) for d in PAIR_DS
    ]


# ================= seeds =================
def frame_pullback(nv):
    """e1 = normalize(dn / dRe zeta), zeta = (n_x + i n_y) / (1 + n_z): the constant frame of the
    plane pulled back through the stereographic chart, one index-2 zero at n = -z."""
    a = nv[..., 0] / (1.0 + nv[..., 2])
    b = nv[..., 1] / (1.0 + nv[..., 2])
    D = 1.0 + a * a + b * b
    e = np.stack([2.0 * (1.0 - a * a + b * b), -4.0 * a * b, -4.0 * a], -1) / (D * D)[..., None]
    nrm = np.linalg.norm(e, axis=-1)
    e = e / np.maximum(nrm, 1e-300)[..., None]
    return e, nrm


def seed_pair(cfg, d, delta):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)

    def w_of(zc):
        zk = Z - zc
        rk = np.sqrt(X * X + Y * Y + zk * zk)
        return (X + 1j * Y) / (rk + zk + 1e-300), rk

    w1, r1 = w_of(+d / 2.0)
    w2, r2 = w_of(-d / 2.0)
    w = w1 * np.conj(w2)
    den = 1.0 + np.abs(w) ** 2
    nv = np.stack([2.0 * w.real, 2.0 * w.imag, 1.0 - np.abs(w) ** 2], -1) / den[..., None]
    # the chart is singular at n = -z: regularize the frame there by a tiny tilt
    nz = nv[..., 2]
    bad = nz < -1.0 + 1e-9
    nv[bad] = np.array([1e-6, 0.0, -1.0]) / np.sqrt(1.0 + 1e-12)
    e1, _ = frame_pullback(nv)
    S = nv[..., :, None] * nv[..., None, :] + delta * e1[..., :, None] * e1[..., None, :]
    a = (1.0 + delta) / 3.0
    u = (1.0 - np.exp(-((r1 / 4.0) ** 2))) * (1.0 - np.exp(-((r2 / 4.0) ** 2)))
    M3 = u[..., None, None] * S + (1.0 - u[..., None, None]) * a * np.eye(3)
    return B3.embed34(M3, cfg)


def seed_of(j, cfg):
    if j["kind"] == "S1":
        return R20.seed_axes(cfg, (1.0, j["delta"], 0.0)), "seed_axes (1, delta, 0)"
    if j["kind"] == "pair":
        return seed_pair(cfg, j["d"], j["delta"]), f"pair texture w1 conj(w2), d {j['d']:g}"
    raise ValueError(j["kind"])


# ================= the descent (the R23-1 reduced instrument at c = 0) =================
def spatial_fmax(M, cfg, p, pot, mask):
    E, Gm, info = CS.energy_grad(M, cfg, p, pot, 0.0)
    Gf = Gm[mask]
    return float(E), float(np.max(np.abs(Gf[:, 1:, 1:]))), float(np.max(np.abs(Gf[:, 0, 0])))


def chunk_light(M, cfg, p, pot, mask, done):
    E, fsp, f00 = spatial_fmax(M, cfg, p, pot, mask)
    return {
        "E": E,
        "iters": done,
        "fmax_spatial": fsp,
        "fmax_M00": f00,
        "M00_range": [float(M[..., 0, 0].min()), float(M[..., 0, 0].max())],
    }


def descend_chunks(M, cfg, p, pot, mask, tag, cap, chunks, done, stage, n_chunks=None):
    """L-BFGS-B chunks on the reduced field until the gate, the cap, or n_chunks; resumable."""
    from scipy.optimize import minimize

    verdict = "FALLING"
    k_run = 0
    while done < cap:
        red = CS.Reduced(M, cfg, p, pot, 0.0, True)
        st = {"it": 0, "gate_hit": None}

        def cb(xk, red=red, st=st):
            st["it"] += 1
            if st["it"] >= 5 and red.last is not None and red.last[1] < GATE:
                st["gate_hit"] = st["it"]

        res = minimize(
            red.fun,
            red.pack(M),
            jac=True,
            method="L-BFGS-B",
            callback=cb,
            options={
                "maxcor": 20,
                "maxiter": CS.CHUNK,
                "maxfun": 3 * CS.CHUNK,
                "gtol": 1e-14,
                "ftol": 1e-16,
            },
        )
        M = red.build(np.asarray(res.x))
        done += max(1, st["it"])
        k_run += 1
        rec = chunk_light(M, cfg, p, pot, mask, done)
        rec["gate_hit_inside_chunk"] = st["gate_hit"]
        rec["chunk_iters"] = st["it"]
        rec["scipy"] = str(res.message)
        rec["drop"] = chunks[-1]["E"] - rec["E"]
        chunks.append(rec)
        np.savez_compressed(stage, M=M, done=done, chunks=json.dumps(chunks))
        log(
            f"{tag} its {done} E {rec['E']:.7g} drop {rec['drop']:.2e} "
            f"fmax_sp {rec['fmax_spatial']:.2e} (gate {GATE:.0e})"
        )
        if (
            rec["fmax_spatial"] < GATE
            and 0 <= rec["drop"] < GATE_DROP * abs(rec["E"])
            and st["it"] >= 5
        ):
            verdict = "AT_GATE"
            break
        if st["it"] < 5 and rec["drop"] <= 0:
            verdict = "LINE_SEARCH_STALL"
            break
        if n_chunks is not None and k_run >= n_chunks:
            verdict = "CHUNKS_DONE"
            break
    return M, done, chunks, verdict


def kick_field(M, mask, pot, tag, amp=KICK_AMP):
    """0.02 Gaussian on the six spatial upper-triangle entries of the free cells, M_00 re-slaved."""
    seed = int(hashlib.sha256(tag.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    Mk = M.copy()
    S = Mk[..., 1:, 1:].copy()
    nf = int(mask.sum())
    blk = np.zeros((nf, 3, 3))
    blk[:, IU3[0], IU3[1]] = amp * rng.standard_normal((nf, 6))
    blk = blk + blk.swapaxes(-1, -2) - np.einsum("...ii->...i", blk)[..., None] * np.eye(3)
    S[mask] = S[mask] + blk
    Mk[..., 1:, 1:] = S
    m00 = Mk[..., 0, 0].copy()
    m00[mask] = CS.solve_m00(S[mask], pot[1], pot[2], 0.0, m00[mask])
    Mk[..., 0, 0] = m00
    return Mk, seed


# ================= reads =================
def _cell_index(x, n, h):
    return int(np.clip(np.rint(x / h + (n - 1) / 2.0), 0, n - 1))


def winding_reads(M, cfg):
    """the transverse line field (middle eigenvector) on circles of radius 2 h around each axis."""
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    n_top, e_mid = V[..., :, 2], V[..., :, 1]
    R = LOOP_RADIUS_H * h
    ph = np.linspace(0.0, 2.0 * np.pi, LOOP_SAMPLES, endpoint=False)
    out = []
    zmax = 0.5 * L - 3.0 * h
    for axis in (2, 0, 1):
        u, v = [a for a in (0, 1, 2) if a != axis]
        for off in LOOP_OFFSETS:
            if off > zmax:
                break
            for sgn in (+1.0, -1.0):
                pts = np.zeros((LOOP_SAMPLES, 3))
                pts[:, axis] = sgn * off
                pts[:, u] = R * np.cos(ph)
                pts[:, v] = R * np.sin(ph)
                idx = tuple(
                    np.array([_cell_index(pts[k, a], n, h) for k in range(LOOP_SAMPLES)])
                    for a in range(3)
                )
                nt = n_top[idx]
                em = e_mid[idx]
                axial = float(np.min(np.abs(nt[:, axis])))
                read = axial > AXIAL_TOL
                wnd = None
                if read:
                    th = np.arctan2(em[:, v], em[:, u])
                    dth = np.diff(np.concatenate([th, th[:1]]))
                    dth = (dth + np.pi / 2.0) % np.pi - np.pi / 2.0
                    wnd = float(np.sum(dth) / np.pi)
                out.append(
                    {
                        "axis": "xyz"[axis],
                        "offset": sgn * off,
                        "read": bool(read),
                        "min_abs_n_axis": axial,
                        "winding_half_turns": wnd,
                    }
                )
    return out


def tube_reads(M, cfg, pot):
    """T_read(z): the density summed over rho < TUBE_RHO on each z plane with abs(z) > TUBE_ZMIN."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    e = R20.density(M, cfg, pot)
    rho = np.sqrt(X * X + Y * Y)
    zs = Z[0, 0, :]
    out = []
    for k, z in enumerate(zs):
        if abs(z) <= TUBE_ZMIN or abs(z) > 0.5 * cfg["L"] - 2.0 * h:
            continue
        disk = rho[:, :, k] < TUBE_RHO
        out.append([float(z), float(e[:, :, k][disk].sum() / h)])
    return out


def strand_geometry(windings, cfg):
    """the read strands: per axis and sign, the winding on the outermost read loop; dl/dL = 1/2 for a
    strand from a centered core to a face."""
    strands = []
    for axis in "zxy":
        for sgn in (+1.0, -1.0):
            loops = [w for w in windings if w["axis"] == axis and np.sign(w["offset"]) == sgn]
            rd = [w for w in loops if w["read"] and w["winding_half_turns"] is not None]
            if not rd:
                continue
            outer = max(rd, key=lambda w: abs(w["offset"]))
            m = float(np.rint(abs(outer["winding_half_turns"])))
            if m > 0:
                strands.append({"axis": axis, "sign": sgn, "m_half_turns": m, "dl_dL": 0.5})
    return strands


def spin_gate_reads(M, cfg, E):
    a0 = B3.gen_catalog(cfg, M)["rot_z"]
    C = float(B3.kin_of(M, a0, cfg))
    out = {"C_rot_z": C, "renv": cfg["renv"], "E_stat": float(E)}
    if C > 0 and E > 0:
        out["omega_star"] = float(np.sqrt(E / (3.0 * C)))
        out["J_star"] = float(np.sqrt(4.0 * C * E / 3.0))
        out["E_rot_over_E"] = 0.25
    else:
        out["omega_star"] = None
        out["J_star"] = None
        out["note"] = "C or E not positive: the gate needs a positive rigid inertia"
    return out


def reads_row(M, cfg, p, pot):
    parts = R20.energy_parts(M, cfg, p, pot)
    wnd = winding_reads(M, cfg)
    out = {
        "energy": parts,
        "windings": wnd,
        "strands": strand_geometry(wnd, cfg),
        "tube_z_T": tube_reads(M, cfg, pot),
        "spin_gate": spin_gate_reads(M, cfg, parts["E_total"]),
        "loop_rule": f"loops read only where abs(n_axis) > {AXIAL_TOL} on every sample",
    }
    try:
        out["chunk_reads"] = CS.chunk_reads(M, cfg, p, pot, 0.0)
    except Exception as e:  # noqa: BLE001
        out["chunk_reads"] = {"error": repr(e)}
    return out


# ================= the job =================
def run_job(j):
    t0 = time.time()
    tag = job_tag(j)
    ug = CS.SL.ups_guard()
    cfg = R21.cfg_of(j["n"], j["L"], G, j["delta"])
    p = R21.params_of(G, j["delta"])
    pot = ("v4", R0.roots_of(cfg), W1 * j["w1s"])
    mask = R21.free_mask(cfg, True)
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    row = dict(j, tag=tag, h=cfg["h"], roots=list(pot[1]), W1_eff=pot[2], gate=GATE)
    try:
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done, chunks = Zs["M"], int(Zs["done"]), json.loads(str(Zs["chunks"]))
            row["resumed_at"] = done
        else:
            M, note = seed_of(j, cfg)
            row["seed_note"] = note
            if np.abs(M[..., 0, 1:]).max() > 0:
                raise RuntimeError("the start field is not block-diagonal")
            done = 0
            chunks = [dict(chunk_light(M, cfg, p, pot, mask, 0), note="the start field")]
            np.savez_compressed(stage, M=M, done=0, chunks=json.dumps(chunks))
        M, done, chunks, verdict = descend_chunks(
            M, cfg, p, pot, mask, tag, j["cap"], chunks, done, stage
        )
        row.update(chunks=chunks, iters=done, gate_label=verdict)
        E_gate = chunks[-1]["E"]
        row["E_gate"] = E_gate
        np.savez_compressed(os.path.join(OUT_NPZ, tag + "_gate.npz"), M=M.astype(np.float64))
        if ug is not None and ug.wrap_up():
            row["stop"] = "UPS wrap-up before the kick (resumable)"
        # the kick-and-continue clause
        if verdict == "AT_GATE":
            Mk, seed = kick_field(M, mask, pot, tag)
            E_k0 = spatial_fmax(Mk, cfg, p, pot, mask)[0]
            kchunks = [dict(chunk_light(Mk, cfg, p, pot, mask, 0), note="the kicked field")]
            kstage = os.path.join(OUT_NPZ, tag + "_kick_stage.npz")
            Mk, kdone, kchunks, kverdict = descend_chunks(
                Mk, cfg, p, pot, mask, tag + "_kick", 10**9, kchunks, 0, kstage, KICK_CHUNKS
            )
            E_k = kchunks[-1]["E"]
            row["kick"] = {
                "amp": KICK_AMP,
                "seed": seed,
                "E_kicked_start": E_k0,
                "E_after": E_k,
                "iters": kdone,
                "chunks": kchunks,
                "verdict": kverdict,
                "fmax_after": kchunks[-1]["fmax_spatial"],
            }
            if E_k < E_gate - KICK_LOWER:
                row["kick_label"] = "SADDLE"
                row["E"] = E_k
                M = Mk
                row["reconverged"] = bool(kchunks[-1]["fmax_spatial"] < GATE)
            else:
                row["kick_label"] = "STABLE"
                row["E"] = E_gate
                row["reconverged"] = True
        else:
            row["kick_label"] = "FALLING"
            row["E"] = E_gate
            row["reconverged"] = False
        row["end_reads"] = reads_row(M, cfg, p, pot)
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M.astype(np.float64))
        row["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc(), E=None)
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} gate {row.get('gate_label')} "
        f"kick {row.get('kick_label')} E {row.get('E')} wall {row['wall_s']}"
    )
    return row


# ================= JSON + pool =================
def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R25-2", "rows": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in jobs if rows.get(job_tag(j), {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


# ================= collect =================
def _fit_slope(Ls, Es):
    Ls, Es = np.asarray(Ls, float), np.asarray(Es, float)
    if len(Ls) < 2:
        return None
    A = np.stack([Ls, np.ones_like(Ls)], 1)
    sol, *_ = np.linalg.lstsq(A, Es, rcond=None)
    return float(sol[0])


def classify_ladder(rows, delta=0.3, w1s=W1S):
    """rows: {tag: row} of the S1 kind at this delta; returns the ladder read and the labels."""
    by_h = {}
    for tag, r in rows.items():
        if r.get("kind") != "S1" or abs(r.get("delta", -1) - delta) > 1e-12:
            continue
        if abs(r.get("w1s", -1) - w1s) > 1e-12:
            continue
        by_h.setdefault(f"{r['h']:g}", []).append(r)
    out = {"per_h": {}, "delta": delta, "w1s": w1s}
    ok_all, within, outside, n_h = True, 0, 0, 0
    for hk, rs in sorted(by_h.items(), key=lambda kv: -float(kv[0])):
        rs = sorted(rs, key=lambda r: r["L"])
        h = float(hk)
        t_half, src = t_half_of(delta, w1s, h)
        pts = [(r["L"], r.get("E"), r.get("kick_label"), r.get("reconverged")) for r in rs]
        good = [pt for pt in pts if pt[1] is not None and pt[2] in ("STABLE", "SADDLE")]
        have = {float(pt[0]) for pt in pts}
        complete = all(L in have for L in EXPECTED_L.get(hk, ()))
        rows_ok = (
            complete
            and len(good) == len(pts)
            and all((pt[2] == "STABLE") or (pt[2] == "SADDLE" and pt[3]) for pt in pts)
        )
        slope = _fit_slope([g[0] for g in good], [g[1] for g in good]) if len(good) >= 2 else None
        pred_centered = 2.0 * t_half
        # the general form from the read windings, on the largest box
        pred_read = None
        if rs and rs[-1].get("end_reads"):
            strands = rs[-1]["end_reads"].get("strands") or []
            if strands:
                pred_read = t_half * sum(s["m_half_turns"] * s["dl_dL"] for s in strands)
        rec = {
            "t_half": t_half,
            "t_half_source": src,
            "points_L_E_kick_reconv": pts,
            "slope": slope,
            "slope_pred_centered": pred_centered,
            "slope_pred_from_windings": pred_read,
            "rows_ok": rows_ok,
            "ladder_complete": complete,
        }
        if slope is not None:
            rec["ratio_to_pred"] = slope / pred_centered
            n_h += 1
            if abs(slope / pred_centered - 1.0) <= TOL_SLOPE:
                within += 1
            else:
                outside += 1
        else:
            ok_all = False
        ok_all = ok_all and rows_ok
        # the spin gate on this spacing
        om = {
            f"{r['L']:g}": (r.get("end_reads") or {}).get("spin_gate", {}).get("omega_star")
            for r in rs
        }
        rec["omega_star_by_L"] = om
        out["per_h"][hk] = rec
    if n_h >= 2 and ok_all and within == n_h:
        label = "CHARGE_STRAND_MATCHES"
    elif n_h >= 2 and ok_all and outside == n_h:
        label = "CHARGE_STRAND_DIFFERS"
    else:
        label = "INSUFFICIENT"
    out["label"] = label
    # the spin gate
    spin = "UNRESOLVED"
    r15 = out["per_h"].get("1.5", {}).get("omega_star_by_L", {})
    r10 = out["per_h"].get("1", {}).get("omega_star_by_L", {})
    ratios = {}
    if r15.get("96") and r15.get("48"):
        ratios["h1.5_96_over_48"] = r15["96"] / r15["48"]
    if r10.get("64") and r10.get("48"):
        ratios["h1_64_over_48"] = r10["64"] / r10["48"]
    if "h1.5_96_over_48" in ratios and abs(ratios["h1.5_96_over_48"] - 1.0) > TOL_SPIN:
        spin = "SPIN_GATE_BOX_DEPENDENT"
    elif len(ratios) == 2 and all(abs(v - 1.0) <= TOL_SPIN for v in ratios.values()):
        spin = "SPIN_GATE_PINNED"
    out["spin_gate"] = {"ratios": ratios, "label": spin}
    return out


def classify_pair(rows, delta=0.3, w1s=W1S):
    Es = {}
    for tag, r in rows.items():
        if r.get("kind") == "pair" and r.get("E") is not None:
            Es[float(r["d"])] = r["E"]
    out = {"E_by_d": {f"{d:g}": e for d, e in sorted(Es.items())}}
    if all(d in Es for d in PAIR_DS):
        s1 = (Es[12.0] - Es[9.0]) / 3.0
        s2 = (Es[15.0] - Es[12.0]) / 3.0
        t_half, src = t_half_of(delta, w1s, 1.5)
        out.update(slope_9_12=s1, slope_12_15=s2, t_half=t_half, t_half_source=src)
        ratio = s2 / s1 if s1 != 0 else None
        out["slope_ratio"] = ratio
        lin = (
            ratio is not None
            and abs(ratio - 1.0) <= TOL_SLOPE
            and abs(0.5 * (s1 + s2) / (4.0 * t_half) - 1.0) <= TOL_SLOPE
        )
        coul = ratio is not None and abs(ratio / 0.6 - 1.0) <= TOL_SLOPE
        out["label"] = (
            "CHARGE_PAIR_LINEAR" if lin else ("CHARGE_PAIR_COULOMB" if coul else "UNRESOLVED")
        )
    else:
        out["label"] = "INSUFFICIENT"
    return out


def collect():
    J = load_json()
    rows = J["rows"]
    C = {"main": classify_ladder(rows, 0.3), "pair": classify_pair(rows)}
    st = classify_ladder(rows, 0.1)
    C["stretch_delta_0.1"] = st
    s03 = C["main"]["per_h"].get("1.5", {}).get("slope")
    s01 = st["per_h"].get("1.5", {}).get("slope")
    if s03 and s01:
        C["stretch_delta_0.1"]["sigma_ratio_0.1_over_0.3"] = s01 / s03
    C["stretch_delta_0.1"]["sigma_ratio_expected"] = (1.0 / 81.0) * K_of(0.05) / K_of(0.15)
    J["collect"] = C
    save_json(J)
    print(json.dumps(C, indent=1))
    return C


def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    J = load_json()
    rows = J["rows"]
    C = J.get("collect") or collect()
    os.makedirs(PLOTS, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    for hk, rec in C["main"]["per_h"].items():
        pts = [(L, E) for L, E, kl, rc in rec["points_L_E_kick_reconv"] if E is not None]
        if not pts:
            continue
        Ls, Es = zip(*pts)
        axs[0].plot(Ls, Es, "o-", label=f"h {hk}: slope {rec['slope']}")
        L0 = np.array([min(Ls), max(Ls)])
        axs[0].plot(L0, Es[0] + rec["slope_pred_centered"] * (L0 - Ls[0]), "--", color="gray")
    axs[0].set_xlabel("L")
    axs[0].set_ylabel("E")
    axs[0].set_title("E(L); dashed: 2 T_half")
    axs[0].legend(fontsize=8)
    for tag, r in rows.items():
        tz = (r.get("end_reads") or {}).get("tube_z_T")
        if tz:
            z, t = zip(*tz)
            axs[1].plot(z, t, ".-", label=tag, lw=0.8)
    axs[1].set_xlabel("z")
    axs[1].set_ylabel("T_read(z)")
    axs[1].set_title("the tube read (rho < 6)")
    axs[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(PLOT_PNG, dpi=130)
    print("plot", PLOT_PNG)


# ================= smoke =================
def smoke():
    out = {}
    # (1) the winding reader on the analytic seed at n 32 L 48 (no descent)
    cfg = R21.cfg_of(32, 48.0, G, 0.3)
    M = R20.seed_axes(cfg, (1.0, 0.3, 0.0))
    wnd = winding_reads(M, cfg)
    zread = [w for w in wnd if w["axis"] == "z" and w["read"]]
    xy_read = [w for w in wnd if w["axis"] in "xy" and w["read"]]
    z_w = [w["winding_half_turns"] for w in zread]
    xy_w = [w["winding_half_turns"] for w in xy_read]
    # on the hedgehog n = r-hat is axial on every axis, so the x and y loops pass the rule and
    # must read ZERO winding (no strand there); the z loops read 2 half-turns (one full turn)
    out["winding_seed"] = {
        "z_loops_read": len(zread),
        "z_windings": z_w,
        "xy_loops_read": len(xy_read),
        "xy_windings": xy_w,
        "ok": bool(zread)
        and all(abs(abs(v) - 2.0) < 0.05 for v in z_w)
        and all(abs(v) < 0.05 for v in xy_w),
    }
    pot = ("v4", R0.roots_of(cfg), W1 * W1S)
    p = R21.params_of(G, 0.3)
    tz = tube_reads(M, cfg, pot)
    out["tube_seed"] = {
        "n_planes": len(tz),
        "T_range": [min(t for _, t in tz), max(t for _, t in tz)],
    }
    out["strands_seed"] = strand_geometry(wnd, cfg)
    out["spin_gate_seed"] = spin_gate_reads(M, cfg, R20.energy_parts(M, cfg, p, pot)["E_total"])
    # (2) the descent and kick wiring on n 16 L 24, two iterations per chunk
    import shutil
    import tempfile

    global OUT_NPZ
    keep = OUT_NPZ
    OUT_NPZ = tempfile.mkdtemp(prefix="r25_2_smoke_", dir=DATA)
    chunk_keep = CS.CHUNK
    CS.CHUNK = 2
    try:
        j = dict(kind="S1", delta=0.3, w1s=W1S, n=16, L=24.0, cap=MAX_ITER[16])
        row = run_job(j)
        out["wire"] = {
            "status": row["status"],
            "gate_label": row.get("gate_label"),
            "kick_label": row.get("kick_label"),
            "iters": row.get("iters"),
            "E_start": row["chunks"][0]["E"] if row.get("chunks") else None,
            "E_end": row.get("E"),
            "stop": row.get("stop"),
        }
        # the kick wiring itself (the n 16 row never reaches the gate inside 4 iterations)
        cfg16 = R21.cfg_of(16, 24.0, G, 0.3)
        pot16 = ("v4", R0.roots_of(cfg16), W1 * W1S)
        mask16 = R21.free_mask(cfg16, True)
        M16 = np.load(os.path.join(OUT_NPZ, job_tag(j) + ".npz"))["M"]
        Mk, seed = kick_field(M16, mask16, pot16, "smoke")
        out["kick_wire"] = {
            "seed": seed,
            "rms_change_free": float(np.sqrt(np.mean((Mk - M16)[mask16] ** 2))),
            "pinned_untouched": bool(np.allclose(Mk[~mask16], M16[~mask16])),
            "block_diagonal": bool(np.abs(Mk[..., 0, 1:]).max() == 0.0),
        }
        out["wire"]["ok"] = row["status"] == "OK" and out["kick_wire"]["pinned_untouched"]
    finally:
        CS.CHUNK = chunk_keep
        shutil.rmtree(OUT_NPZ)
        OUT_NPZ = keep
    # (3) the classify function on synthetic ladders
    t_half = T_half_bps(0.3, W1 * W1S)

    def synth(fac, drop=None):
        rows = {}
        for n, L in [(32, 48.0), (48, 72.0), (64, 96.0), (48, 48.0), (64, 64.0)]:
            if drop == (n, L):
                continue
            h = L / n
            rows[f"S1_d0.3_w25_n{n}_L{L:g}"] = dict(
                kind="S1",
                delta=0.3,
                w1s=W1S,
                n=n,
                L=L,
                h=h,
                E=10.0 + fac * 2.0 * t_half * L,
                kick_label="STABLE",
                reconverged=True,
                end_reads={"spin_gate": {"omega_star": 1.0}, "strands": []},
            )
        return rows

    lab = {
        "match": classify_ladder(synth(1.0))["label"],
        "twice": classify_ladder(synth(2.0))["label"],
        "missing": classify_ladder(synth(1.0, drop=(64, 96.0)))["label"],
        "spin_pinned": classify_ladder(synth(1.0))["spin_gate"]["label"],
    }
    out["classify"] = lab
    out["classify"]["ok"] = (
        lab["match"] == "CHARGE_STRAND_MATCHES"
        and lab["twice"] == "CHARGE_STRAND_DIFFERS"
        and lab["missing"] == "INSUFFICIENT"
        and lab["spin_pinned"] == "SPIN_GATE_PINNED"
    )
    pr = {
        f"pair_d{d:g}": dict(kind="pair", d=d, E=1.0 + 4.0 * t_half * d, delta=0.3, w1s=W1S)
        for d in PAIR_DS
    }
    pc = {
        f"pair_d{d:g}": dict(kind="pair", d=d, E=1.0 - 0.5 / d, delta=0.3, w1s=W1S)
        for d in PAIR_DS
    }
    out["classify_pair"] = {
        "linear": classify_pair(pr)["label"],
        "coulomb": classify_pair(pc)["label"],
    }
    out["classify_pair"]["ok"] = (
        out["classify_pair"]["linear"] == "CHARGE_PAIR_LINEAR"
        and out["classify_pair"]["coulomb"] == "CHARGE_PAIR_COULOMB"
    )
    # (4) wiring
    tags = [job_tag(j) for j in jobs_main() + jobs_stretch() + jobs_pair()]
    out["wiring"] = {
        "jobs": len(tags),
        "tags_unique": len(set(tags)) == len(tags),
        "r25_1_json_present": os.path.exists(R25_1_JSON),
        "t_half_bogomolny_delta_0.3_w25": t_half,
        "t_half_fallback_ok": abs(t_half - 2.675e-3) < 2e-5,
        "t_half_source_now": t_half_of(0.3, W1S, 1.5)[1],
        "slope_pred_centered": 2.0 * t_half,
    }
    out["PASS"] = bool(
        out["winding_seed"]["ok"]
        and out["wire"]["ok"]
        and out["classify"]["ok"]
        and out["classify_pair"]["ok"]
        and out["wiring"]["tags_unique"]
        and out["wiring"]["t_half_fallback_ok"]
    )
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_main(), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "run_stretch":
        run_pool(jobs_stretch(), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "run_pair":
        run_pool(jobs_pair(), sys.argv[2] if len(sys.argv) > 2 else 12)
    elif mode == "collect":
        collect()
    elif mode == "plot":
        plot()
    elif mode == "jobs":
        for j in jobs_main() + jobs_stretch() + jobs_pair():
            print(job_tag(j), j)
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
