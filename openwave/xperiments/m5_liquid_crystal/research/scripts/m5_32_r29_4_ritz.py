"""M5.32 R29-4: the corrected Ritz matrix across the stored K-positive modes of the R29-1 pencil, per rung, with the
reduced-pencil residual (the author's instrument clarification of 2026-10-29).

EQUATIONS FIRST
---------------
R29-1 reads the lowest omega^2 of the pencil (U+^T H U+, B+) about a relaxed field M*: H the Hessian of the static
energy (central differences of the gradient), B+ the K weight on the K-positive coordinates, the generator span
projected out. Its "schur-post" step relaxed the K-null coordinates by conjugate gradients but held each saved bare
mode fixed and read mode 0: an auxiliary-relaxed Rayleigh quotient, not an eigenvalue of the eliminated pencil.

R29-4 forms, on every stored field, the Ritz matrix across the saved K-positive modes and diagonalises it:
    Z          the union of every saved K-positive mode file of the rung, restricted to the K-positive coordinates
               (lossless: the saved modes live there), the generator span removed by the B-orthogonal projection
               z -> z - Y (Y^T B Y)^-1 Y^T B z, then B-orthonormalised through the eigen-decomposition of Z B Z^T
               (near-null directions dropped at 1e-10 of the largest; the rank and the spectrum are reported);
    H_eff Z_j = H (Z_j + U0 z_j),  z_j the solution of  H00 z = -(H Z_j)_0  by R29-1's preconditioned conjugate
               gradients (the per-cell block inverse as the preconditioner; 1e-6 relative residual, 30 steps cap),
               the FIELD returned rather than the scalar D . H_eff D of R29-1;
    A_ij      = Z_i . H_eff Z_j   (the Frobenius inner product equals the dot product of the K-eigen coordinates),
               symmetrised, the asymmetry  max |A_ij - A_ji| / max |A|  reported as the measure of the CG truncation;
    B_ij      = K(Z_i, Z_j) = sum_x kappa_x z_i z_j  (the polarisation of K_of; exact, K is quadratic; the identity
               on the orthonormalised basis, the deviation reported);
    A c = lambda B c,   omega^2_Ritz = lambda_min / 2;
    residual  = || P (H_eff y - lambda B y) || / (|lambda| ||B y||),   y = sum_j c_j Z_j,  P the Euclidean projector
               onto the K-positive coordinates orthogonal to the span of B Y (the Lagrange multipliers of the generator
               constraint live in span(B Y): a constrained eigenpair has its residual there), assembled from the stored
               H_eff Z_j coordinates (linear in them, no new gradients); the unprojected norm is reported beside it;
    the same Ritz matrix without the elimination (H Z_j in place of H_eff Z_j, the bare pencil) beside it, with its own
    residual: it also reads the n 72 rung, whose LOBPCG stopped at its iteration cap.

CONTROL (fails_if): the bare Ritz matrix restricted to the members of one rung's FINAL mode file must return that
rung's recorded bare LOBPCG value (the collect record's omega2_bare); if it does not, the basis reader is wrong and
nothing else is read. Then the 1/L and 1/L^2 fits of R29-1's collect are rerun on the Ritz-corrected values beside
the bare and the auxiliary-relaxed fits (fit-only errors).

Pre-registered words (the task doc, § R29-4): CONSISTENT / MOVED / NOT CONVERGED.

Modes (python3 m5_32_r29_4_ritz.py MODE ...):
    prep     [TAG ...]     the Ritz basis of each rung (reads the modes, builds Z, writes *_ritz_Z.npy + *_ritz_meta.npz)
    apply    [--workers W] the pool over the (rung, basis vector) jobs: H Z_j and H_eff Z_j, one npz per job (resumable)
    collect                the Ritz matrices, the control, the residuals, the fits, the words, the plot
    run      [--workers W] prep + apply + collect
Set OMP_NUM_THREADS=1 per worker (the script sets the four thread variables to 1 before NumPy is imported; the pool
size is the machine's, never written here).
Out: ../data/m5_32_r29_1/*_ritz_*  (the basis and the H fields, local), ../data/m5_32_r29_4_ritz.json (tracked),
     ../plots/m5_32_r29_4_ritz.png.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import glob
import importlib.util
import json
import multiprocessing as mp
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
IN_DIR = os.path.join(DATA, "m5_32_r29_1")
OUT_JSON = os.path.join(DATA, "m5_32_r29_4_ritz.json")
PLOT = os.path.join(HERE, "..", "plots", "m5_32_r29_4_ritz.png")
COLLECT_R29_1 = os.path.join(DATA, "m5_32_r29_1_gap.json")
T0 = time.time()

# the five stored fields (tag, n, L, h) of the R29 record; the h 1.0 box is the mesh check
RUNGS = [("rad_pin_d0.3_w25_c0_n72_L108", 72, 108.0, 1.5),
         ("rad_pin_d0.3_w25_c0_n48_L72", 48, 72.0, 1.5),
         ("rad_pin_d0.3_w25_c0_n48_L48", 48, 48.0, 1.0),
         ("rad_pin_d0.3_w25_c0_n36_L54", 36, 54.0, 1.5),
         ("rad_pin_d0.3_w25_c0_n32_L48", 32, 48.0, 1.5)]
EPS = 1e-3            # the finite-difference step of R29-1's S2_apply
CG_TOL, CG_ITERS = 1e-6, 30
RANK_TOL = 1e-10      # the B-Gram cut, relative to the largest Gram eigenvalue


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R29 = _load("m5_32_r29_1_gap", "m5_32_r29_1_gap.py")
C23 = R29.C23


# ================= the rung, cached per process =================
_RUNG = {}


def get_rung(tag, n, L):
    if tag in _RUNG:
        return _RUNG[tag]
    M = np.load(os.path.join(C23.OUT_NPZ, tag + ".npz"))["M"]
    assert M.shape[0] == n
    R = R29.Rung(n, L, M)
    R.build_basis()
    _RUNG.clear()
    _RUNG[tag] = R
    return R


def generator_coords(R):
    """the generator fields in the K-positive coordinates, the independent ones (as R29-1's pencil keeps them)."""
    gens = R29.generators(R)
    Yc = np.array([R.field_to_coords(D, R.pos) for D in gens.values()]).T       # (dim, 12)
    u, s, _ = np.linalg.svd(Yc, full_matrices=False)
    keep = s > 1e-8 * s[0]
    return u[:, keep], list(gens.keys()), s.tolist()


# ================= prep: the Ritz basis of a rung =================
def mode_files(tag):
    return sorted(glob.glob(os.path.join(IN_DIR, f"{tag}_Kpos_modes*.npz")))


def prep_rung(tag, n, L):
    R = get_rung(tag, n, L)
    b = R.kev[R.pos]                                                          # the K weight on K+, positive
    Y, gnames, gsv = generator_coords(R)                                        # Euclidean-orthonormal columns
    YtBY = Y.T @ (b[:, None] * Y)
    members, Zr = [], []
    for f in mode_files(tag):
        Zf = np.load(f)["modes"]
        for i in range(Zf.shape[0]):
            z = R.field_to_coords(Zf[i], R.pos)
            back = R.coords_to_field(z, R.pos)
            loss = float(np.sqrt(np.sum((back - Zf[i]) ** 2)) / np.sqrt(np.sum(Zf[i] ** 2)))   # the K+ restriction loss (0 expected)
            zg = z - Y @ np.linalg.solve(YtBY, Y.T @ (b * z))                 # the B-orthogonal projection out of the generator span
            gen_part = float(np.sqrt(np.sum(b * (z - zg) ** 2)) / np.sqrt(np.sum(b * z * z)))
            members.append({"file": os.path.basename(f), "index": i, "restriction_loss": loss, "generator_component_B_norm_ratio": gen_part,
                            "final_file": f.endswith(f"{tag}_Kpos_modes.npz")})
            Zr.append(zg)
    Zr = np.array(Zr)                                                           # (m, dim)
    Gm = Zr @ (b[:, None] * Zr.T)                                               # the B Gram
    w, V = np.linalg.eigh(Gm)
    keep = w > RANK_TOL * w.max()
    Zo = (V[:, keep] / np.sqrt(w[keep])).T @ Zr                                 # (r, dim), B-orthonormal
    C = Zr @ (b[:, None] * Zo.T)                                                # the raw members in the orthonormal basis: Zr ~ C Zo
    recon = float(np.max(np.sqrt(np.sum(b * (Zr - C @ Zo) ** 2, axis=1)) / np.sqrt(np.sum(b * Zr * Zr, axis=1))))
    np.save(os.path.join(IN_DIR, f"{tag}_ritz_Z.npy"), Zo)
    np.savez(os.path.join(IN_DIR, f"{tag}_ritz_meta.npz"), b=b, Y=Y, C=C, gram_eigs=w, members=json.dumps(members),
             gnames=json.dumps(gnames), gsv=np.array(gsv), n=n, L=L, rank=int(keep.sum()), recon=recon)
    rec = {"tag": tag, "n": n, "L": L, "dim_Kpos": int(R.pos.sum()), "members": len(members), "rank": int(keep.sum()),
           "gram_eigenvalues": w.tolist(), "generator_constraints": int(Y.shape[1]), "max_member_reconstruction_B_err": recon,
           "max_restriction_loss": max(m["restriction_loss"] for m in members),
           "max_generator_component": max(m["generator_component_B_norm_ratio"] for m in members), "files": [os.path.basename(f) for f in mode_files(tag)]}
    log(f"prep {tag}: {len(members)} members, rank {rec['rank']}, gram {w.min():.2e}..{w.max():.2e}, gen comp max {rec['max_generator_component']:.2e}")
    return rec


def _prep_worker(args):
    return prep_rung(*args)


# ================= apply: H Z_j and H_eff Z_j for one basis vector =================
def heff_field(R, D, eps=EPS, iters=CG_ITERS, tol=CG_TOL):
    """R29-1's S2_schur_cg, returning the fields: (H D, H (D + U0 z), cg residual ratio, steps, |U0 z| / |D|)."""
    HD = R.S2_apply(D, eps)
    inv = R.null_block_inverse()
    BX, BW = R29.BX, R29.BW

    def to_null(F):
        x = np.einsum("fab,cab->fc", F[R.free], BX) / BW
        return np.einsum("fij,fi->fj", R.kU, x * R.sw)[R.nul]

    def from_null(zn):
        return R.coords_to_field(zn, R.nul)

    def H00(zn):
        return to_null(R.S2_apply(from_null(zn), eps))

    def Pinv(rn):
        zz = np.zeros(R.kU.shape[:2]); zz[R.nul] = rn
        return np.einsum("fij,fj->fi", inv, zz)[R.nul]

    bvec = -to_null(HD)
    z = Pinv(bvec); r = bvec - H00(z); y = Pinv(r); pdir = y.copy(); ry = float(r @ y); b2 = float(np.sqrt(bvec @ bvec)); steps = 0
    for steps in range(1, iters + 1):
        Ap = H00(pdir); alpha = ry / float(pdir @ Ap)
        z = z + alpha * pdir; r = r - alpha * Ap
        if np.sqrt(float(r @ r)) < tol * b2:
            break
        y = Pinv(r); ry_new = float(r @ y); pdir = y + (ry_new / ry) * pdir; ry = ry_new
    Dz = from_null(z)
    HeffD = R.S2_apply(D + Dz, eps)
    return HD, HeffD, float(np.sqrt(float(r @ r)) / b2), steps, float(np.sqrt(np.sum(Dz * Dz)) / np.sqrt(np.sum(D * D)))


def job_path(tag, j):
    return os.path.join(IN_DIR, f"{tag}_ritz_H_{j:02d}.npz")


def _apply_worker(args):
    tag, n, L, j = args
    out = job_path(tag, j)
    if os.path.exists(out):
        return {"tag": tag, "j": j, "skipped": True}
    t0 = time.time()
    R = get_rung(tag, n, L)
    Zo = np.load(os.path.join(IN_DIR, f"{tag}_ritz_Z.npy"), mmap_mode="r")
    D = R.coords_to_field(np.array(Zo[j]), R.pos)
    R.ngrad = 0
    HD, HeffD, res, steps, zr = heff_field(R, D)
    hd = R.field_to_coords(HD, R.pos); he = R.field_to_coords(HeffD, R.pos)
    # the parts of H D outside K+ (held at zero by the pencil): reported, not used
    hd_out = float(np.sqrt(np.sum(HD * HD) - np.sum(hd * hd)) / np.sqrt(np.sum(HD * HD)))
    np.savez(out, hd=hd, he=he, cg_res=res, cg_steps=steps, z_ratio=zr, gradients=R.ngrad, wall=time.time() - t0, hd_outside_Kpos=hd_out)
    return {"tag": tag, "j": j, "cg_res": res, "cg_steps": steps, "z_ratio": zr, "gradients": R.ngrad, "wall_s": round(time.time() - t0, 1), "hd_outside_Kpos": hd_out}


def apply_all(workers, tags=None):
    jobs = []
    for tag, n, L, h in RUNGS:
        if tags and tag not in tags:
            continue
        meta = np.load(os.path.join(IN_DIR, f"{tag}_ritz_meta.npz"))
        for j in range(int(meta["rank"])):
            jobs.append((tag, n, L, j))                                       # RUNGS is ordered largest first
    log(f"apply: {len(jobs)} jobs on {workers} workers")
    try:
        sys.path.insert(0, os.path.expanduser("~/Library/Application Support/UPSSentinel")); import ups_guard; ups_guard.register(); guard = ups_guard
    except Exception:
        guard = None
    rows = []
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        for row in pool.imap_unordered(_apply_worker, jobs):
            rows.append(row)
            if not row.get("skipped"):
                log(f"{row['tag']} j{row['j']}: cg res {row['cg_res']:.1e} in {row['cg_steps']} steps, |z|/|D| {row['z_ratio']:.3f}, {row['gradients']} grads, {row['wall_s']} s")
            if guard is not None and guard.wrap_up():
                log("UPS wrap-up flag: stopping at a job boundary (finished jobs are on disk; rerun `apply` to resume)")
                pool.terminate(); break
    return rows


# ================= collect =================
def ols(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
    nres = len(x) - 2
    resid = float(np.sqrt(res[0] / max(1, nres))) if len(res) else 0.0
    cov = np.linalg.inv(A.T @ A) * (res[0] / nres if (len(res) and nres > 0) else 0.0)
    return {"intercept": float(coef[0]), "slope": float(coef[1]), "residual": resid, "intercept_se_fit": float(np.sqrt(cov[0, 0])), "n": int(len(x))}


def ritz(A, Bm):
    A = 0.5 * (A + A.T)
    w, Vb = np.linalg.eigh(Bm)
    Wh = (Vb / np.sqrt(w)) @ Vb.T
    lam, U = np.linalg.eigh(Wh @ A @ Wh)
    return lam, Wh @ U


def collect_rung(tag, n, L, h, rec29):
    meta = np.load(os.path.join(IN_DIR, f"{tag}_ritz_meta.npz"))
    Zo = np.load(os.path.join(IN_DIR, f"{tag}_ritz_Z.npy"))
    b, Y, C = meta["b"], meta["Y"], meta["C"]
    members = json.loads(str(meta["members"])); r = int(meta["rank"])
    HD = np.zeros_like(Zo); HE = np.zeros_like(Zo); stats = []
    for j in range(r):
        z = np.load(job_path(tag, j))
        HD[j], HE[j] = z["hd"], z["he"]
        stats.append({"j": j, "cg_res": float(z["cg_res"]), "cg_steps": int(z["cg_steps"]), "z_ratio": float(z["z_ratio"]), "gradients": int(z["gradients"]), "wall_s": float(z["wall"]), "hd_outside_Kpos": float(z["hd_outside_Kpos"])})
    Bm = Zo @ (b[:, None] * Zo.T)
    B_dev = float(np.max(np.abs(Bm - np.eye(r))))
    out = {"tag": tag, "n": n, "L": L, "h": h, "rank": r, "members": len(members), "B_identity_deviation": B_dev, "jobs": stats}
    BY = b[:, None] * Y
    Q, _ = np.linalg.qr(BY)                                                    # the span of B Y (the multiplier space)
    Qy, _ = np.linalg.qr(Y)
    for name, Hm in (("bare", HD), ("corrected", HE)):
        A = Zo @ Hm.T                                                          # A_ij = Z_i . H Z_j
        asym = float(np.max(np.abs(A - A.T)) / np.max(np.abs(A)))
        lam, Cv = ritz(A, Bm)
        c = Cv[:, 0]
        y = c @ Zo; Hy = c @ Hm; By = b * y
        res = Hy - lam[0] * By
        nrm = abs(lam[0]) * float(np.sqrt(By @ By))
        rp = res - Q @ (Q.T @ res)
        ry = res - Qy @ (Qy.T @ res)
        out[name] = {"ritz_values_lambda": lam.tolist(), "omega2_lowest": float(lam[0] / 2), "omega2_second": float(lam[1] / 2) if r > 1 else None,
                     "asymmetry": asym, "residual_projected_out_of_BY_span": float(np.sqrt(rp @ rp) / nrm), "residual_projected_out_of_Y_span": float(np.sqrt(ry @ ry) / nrm),
                     "residual_unprojected": float(np.sqrt(res @ res) / nrm),
                     "residual_abs_B_normalised": float(np.sqrt(rp @ rp) / np.sqrt(By @ y)),      # LOBPCG's own convention: ||A x - lambda B x||, x B-normalised
                     "generator_component_of_ritz_vector_B": float(np.sqrt(np.sum(b * (Y @ (Y.T @ (b * y))) ** 2)) / np.sqrt(By @ y)),
                     "ritz_vector_weights_on_members": None}
        # the control: every stored member's bare Rayleigh quotient (through C and A, no new gradients) against the
        # LOBPCG runs recorded for this rung (the gap_*.json files): a stored FILE is validated when its mode 0 returns
        # one recorded run's lowest value, and the record's ladder run is flagged when it is the one on disk. R29-1's
        # mode_gap wrote every run's modes to the same `{tag}_Kpos_modes.npz`, so a later diagnostic run can have
        # overwritten the ladder run's modes (the n 32 rung: the pass-2 modes were overwritten by the shell-warm run).
        if name == "bare":
            runs = {}
            for f in sorted(glob.glob(os.path.join(IN_DIR, f"gap_{tag}*.json"))):
                d = json.load(open(f)); kp = d.get("K_positive", {})
                if kp.get("modes"):
                    runs[os.path.basename(f)] = [mm["omega2_S2_over_2K"] for mm in kp["modes"]]
            member_rq = {}
            for k, mm in enumerate(members):
                rq = float(C[k] @ A @ C[k]) / float(C[k] @ Bm @ C[k]) / 2
                member_rq[f"{mm['file']}[{mm['index']}]"] = rq
            files = {}
            for mm in members:
                files.setdefault(mm["file"], []).append(mm)
            ctrl = {}
            for fname, mems in files.items():
                m0 = [mm for mm in mems if mm["index"] == 0][0]
                rq0 = member_rq[f"{fname}[0]"]
                best = min(((abs(rq0 / v[0] - 1.0), rn) for rn, v in runs.items()), default=(None, None))
                ctrl[fname] = {"mode0_bare_omega2": rq0, "matches_run": best[1], "rel_diff": best[0], "is_the_ladder_run": best[1] == rec29["gap"],
                               "ritz_over_this_file": None}
                idx = [k for k, mm in enumerate(members) if mm["file"] == fname]
                Cf = C[idx]
                lam_f, _ = ritz(Cf @ A @ Cf.T, Cf @ Bm @ Cf.T)
                ctrl[fname]["ritz_over_this_file"] = float(lam_f[0] / 2)
            out["control_stored_files_vs_recorded_runs"] = ctrl
            out["member_bare_omega2"] = member_rq
            out["recorded_runs_omega2"] = runs
            out["ladder_run_modes_on_disk"] = any(c["is_the_ladder_run"] and c["rel_diff"] is not None and c["rel_diff"] < 0.01 for c in ctrl.values())
            out["control_ok"] = all(c["rel_diff"] is not None and c["rel_diff"] < 0.01 for c in ctrl.values())
    out["recorded"] = {"omega2_bare": rec29["omega2_bare"], "omega2_auxiliary_relaxed": rec29["omega2_cg"], "omega2_auxiliary_relaxed_second": rec29.get("omega2_cg_second"),
                       "lobpcg_iterations": rec29.get("iterations")}
    out["drop_vs_auxiliary_relaxed"] = float(1.0 - out["corrected"]["omega2_lowest"] / rec29["omega2_cg"])
    out["drop_vs_bare_recorded"] = float(1.0 - out["corrected"]["omega2_lowest"] / rec29["omega2_bare"])
    out["bare_ritz_drop_vs_bare_recorded"] = float(1.0 - out["bare"]["omega2_lowest"] / rec29["omega2_bare"])
    out["max_cg_residual"] = max(s["cg_res"] for s in stats); out["max_cg_steps"] = max(s["cg_steps"] for s in stats)
    out["gradients_total"] = sum(s["gradients"] for s in stats); out["wall_total_s"] = round(sum(s["wall_s"] for s in stats), 1)
    return out


def collect():
    rec29 = {r["tag"]: r for r in json.load(open(COLLECT_R29_1))["rungs"]}
    rungs = [collect_rung(tag, n, L, h, rec29[tag]) for tag, n, L, h in RUNGS]
    rungs.sort(key=lambda r: (r["h"] != 1.5, r["L"]))
    ladder = [r for r in rungs if r["h"] == 1.5]
    mesh = [r for r in rungs if r["h"] == 1.0][0]
    same_box = [r for r in ladder if r["L"] == mesh["L"]][0]
    x = np.array([1.0 / r["L"] for r in ladder])
    series = {"bare_lobpcg": [r["recorded"]["omega2_bare"] for r in ladder], "auxiliary_relaxed": [r["recorded"]["omega2_auxiliary_relaxed"] for r in ladder],
              "ritz_bare": [r["bare"]["omega2_lowest"] for r in ladder], "ritz_corrected": [r["corrected"]["omega2_lowest"] for r in ladder]}
    mesh_vals = {"bare_lobpcg": (same_box["recorded"]["omega2_bare"], mesh["recorded"]["omega2_bare"]), "auxiliary_relaxed": (same_box["recorded"]["omega2_auxiliary_relaxed"], mesh["recorded"]["omega2_auxiliary_relaxed"]),
                 "ritz_bare": (same_box["bare"]["omega2_lowest"], mesh["bare"]["omega2_lowest"]), "ritz_corrected": (same_box["corrected"]["omega2_lowest"], mesh["corrected"]["omega2_lowest"])}
    fits = {}
    for k, ys in series.items():
        y = np.array(ys)
        f1, f2 = ols(x, y), ols(x * x, y)
        md = float(mesh_vals[k][0] - mesh_vals[k][1])
        fits[k] = {"points_(L, omega2)": [[r["L"], v] for r, v in zip(ladder, ys)], "fit_in_1_over_L": f1, "fit_in_1_over_L2": f2,
                   "mesh_difference_(h 1.5 minus h 1.0, same box)": md, "intercept_error_quadrature": float(np.sqrt(f1["intercept_se_fit"] ** 2 + md ** 2)),
                   "intercept_1_over_L_within_2se_of_zero": bool(abs(f1["intercept"]) < 2 * f1["intercept_se_fit"]), "intercept_1_over_L2_sign": int(np.sign(f2["intercept"])),
                   "mesh_term_largest": bool(abs(md) > max(abs(f1["intercept"]), f1["intercept_se_fit"]))}
    # the check that this fit routine reproduces the R29-1 record's bare fit
    rec_fit = json.load(open(COLLECT_R29_1))["fit_omega2_bare"]["fit_in_1_over_L"]
    fits["reproduces_r29_1_bare_fit"] = {"intercept_rel_diff": float(fits["bare_lobpcg"]["fit_in_1_over_L"]["intercept"] / rec_fit["intercept"] - 1.0),
                                         "se_rel_diff": float(fits["bare_lobpcg"]["fit_in_1_over_L"]["intercept_se_fit"] / rec_fit["intercept_se_fit"] - 1.0)}
    # the words
    drops = {r["tag"]: r["drop_vs_auxiliary_relaxed"] for r in rungs}
    resid = {r["tag"]: r["corrected"]["residual_projected_out_of_BY_span"] for r in rungs}
    control_ok = all(r["control_ok"] for r in rungs)
    fc = fits["ritz_corrected"]; fb = fits["bare_lobpcg"]
    not_conv = [t for t, v in resid.items() if v > 1e-2]
    moved = [t for t, v in drops.items() if v > 0.10]
    verdict_change = (fc["intercept_1_over_L_within_2se_of_zero"] != fb["intercept_1_over_L_within_2se_of_zero"]) or (fc["intercept_1_over_L2_sign"] != fb["intercept_1_over_L2_sign"]) or (fc["mesh_term_largest"] != fb["mesh_term_largest"])
    out_ctrl = {r["tag"]: {"control_ok": r["control_ok"], "ladder_run_modes_on_disk": r["ladder_run_modes_on_disk"]} for r in rungs}
    if not control_ok:
        word = "CONTROL FAILED: a stored mode file does not return any recorded LOBPCG value; nothing else is read"
    else:
        word = ("NOT CONVERGED at " + ", ".join(not_conv) + "; " if not_conv else "") + ("MOVED" if (moved or verdict_change) else "CONSISTENT")
    out = {"rungs": rungs, "fits": fits, "drops_vs_auxiliary_relaxed": drops, "corrected_residuals": resid, "control_ok": control_ok, "controls": out_ctrl,
           "lobpcg_style_abs_residuals": {r["tag"]: {"bare": r["bare"]["residual_abs_B_normalised"], "corrected": r["corrected"]["residual_abs_B_normalised"]} for r in rungs},
           "rungs_not_converged": not_conv, "rungs_moved_above_10_percent": moved, "fit_verdict_changed": verdict_change, "word": word,
           "pre_registered_read_r29_1": json.load(open(COLLECT_R29_1))["fit_omega2_bare"].get("pre_registered_read"),
           "eps": EPS, "cg_tol": CG_TOL, "cg_iters": CG_ITERS, "rank_tol": RANK_TOL}
    json.dump(out, open(OUT_JSON, "w"), indent=1, default=float)
    plot(out, ladder, mesh)
    return out


def plot(out, ladder, mesh):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.array([1.0 / r["L"] for r in ladder]); xx = np.linspace(0, x.max() * 1.05, 50)
    for k, lab, mk in (("bare_lobpcg", "bare LOBPCG (R29-1)", "o"), ("auxiliary_relaxed", "auxiliary-relaxed Rayleigh quotient (R29-1 schur-post)", "s"),
                       ("ritz_bare", "Ritz, bare pencil (R29-4)", "^"), ("ritz_corrected", "Ritz, eliminated pencil (R29-4)", "D")):
        f = out["fits"][k]; y = [p[1] for p in f["points_(L, omega2)"]]
        ax[0].plot(x, y, mk, label=lab); ax[0].plot(xx, f["fit_in_1_over_L"]["intercept"] + f["fit_in_1_over_L"]["slope"] * xx, "-", lw=0.8, alpha=0.6)
    ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xlabel("1 / L"); ax[0].set_ylabel("lowest omega^2 on K+"); ax[0].set_title("the ladder at h 1.5, fits in 1/L"); ax[0].legend(fontsize=7)
    ax[0].plot([1.0 / mesh["L"]], [mesh["corrected"]["omega2_lowest"]], "kx", label="h 1.0 mesh check (corrected)")
    tags = [r["n"] for r in ladder] + [f"{mesh['n']} (h1)"]
    rr = ladder + [mesh]
    ax[1].bar(np.arange(len(rr)) - 0.2, [r["bare"]["residual_projected_out_of_BY_span"] for r in rr], 0.4, label="bare")
    ax[1].bar(np.arange(len(rr)) + 0.2, [r["corrected"]["residual_projected_out_of_BY_span"] for r in rr], 0.4, label="corrected")
    ax[1].axhline(1e-2, color="r", lw=0.8, ls="--", label="NOT CONVERGED bar"); ax[1].set_yscale("log"); ax[1].set_xticks(range(len(rr))); ax[1].set_xticklabels([f"n {t}" for t in tags])
    ax[1].set_ylabel("reduced-pencil residual of the Ritz vector"); ax[1].set_title("residuals"); ax[1].legend(fontsize=7)
    short = ("NOT CONVERGED at every rung" if len(out["rungs_not_converged"]) == len(out["rungs"]) else f"NOT CONVERGED at {len(out['rungs_not_converged'])} rungs") + ("; MOVED at " + ", ".join(f"n {r['n']}" for r in out["rungs"] if r["tag"] in out["rungs_moved_above_10_percent"]) if out["rungs_moved_above_10_percent"] else "; CONSISTENT")
    fig.suptitle(f"M5.32 R29-4: the Ritz matrix across the stored K+ modes; word: {short}; the R29 read unchanged (NOT RESOLVED at h 1.5)", fontsize=9)
    fig.tight_layout(); os.makedirs(os.path.dirname(PLOT), exist_ok=True); fig.savefig(PLOT, dpi=130); plt.close(fig)


def regen(tag, n, L, nev, warm, maxiter):
    """a fresh K-positive LOBPCG run with R29-1's pencil (the generator span projected out), saved as
    `{tag}_Kpos_modes_r29_4regen.npz`, for a rung whose ladder-run modes were overwritten on disk (the n 32 rung)."""
    R = get_rung(tag, n, L)
    gens = R29.generators(R)
    res = R.pencil(R.pos, nev, Y_fields=list(gens.values()), label="gap", warm=warm, maxiter=maxiter)
    np.savez_compressed(os.path.join(IN_DIR, f"{tag}_Kpos_modes_r29_4regen.npz"), modes=R.last_mode_fields)
    rec = {"field": tag, "n": n, "L": L, "warm": warm, "K_positive": res}
    json.dump(rec, open(os.path.join(IN_DIR, f"gap_{tag}_r29_4regen.json"), "w"), indent=1, default=float)
    log(f"regen {tag}: lowest omega2 {res['modes'][0]['omega2_S2_over_2K']}, {res['iterations']} iterations, {res['gradients']} gradients, {res['wall_s']} s")
    return rec


def verify(tag, n, L):
    """the Ritz vector of the corrected pencil re-read with R29-1's own instrument: D . H_eff D by `S2_schur_cg` and
    K by `K_of` on the field of the Ritz vector, against lambda_min / 2 of the Ritz matrix (a check that the matrix
    assembly, the coordinates and the CG fields agree with the scalar instrument of R29-1's schur-post); also every
    stored member's own corrected Rayleigh quotient, so a drop can be attributed to one member or to the coupling."""
    meta = np.load(os.path.join(IN_DIR, f"{tag}_ritz_meta.npz")); Zo = np.load(os.path.join(IN_DIR, f"{tag}_ritz_Z.npy"))
    b, C = meta["b"], meta["C"]; members = json.loads(str(meta["members"])); r = int(meta["rank"])
    HE = np.array([np.load(job_path(tag, j))["he"] for j in range(r)]); HD = np.array([np.load(job_path(tag, j))["hd"] for j in range(r)])
    Bm = Zo @ (b[:, None] * Zo.T)
    A = Zo @ HE.T; lam, Cv = ritz(A, Bm); c = Cv[:, 0]; y = c @ Zo
    Ab = Zo @ HD.T
    R = get_rung(tag, n, L)
    D = R.coords_to_field(y, R.pos); D /= np.sqrt(np.sum(D * D))
    s2_bare = float(np.sum(D * R.S2_apply(D))); k = R.K_of(D)
    s2_cg, res, steps, zr = R.S2_schur_cg(D)
    out = {"tag": tag, "omega2_ritz_corrected": float(lam[0] / 2), "omega2_direct_S2_schur_cg_on_ritz_vector": s2_cg / (2 * k), "omega2_direct_bare_on_ritz_vector": s2_bare / (2 * k),
           "cg_residual": res, "cg_steps": steps, "null_response_ratio": zr,
           "rel_diff_corrected": float((s2_cg / (2 * k)) / (lam[0] / 2) - 1.0),
           "member_corrected_omega2": {f"{m['file']}[{m['index']}]": float(C[i] @ A @ C[i]) / float(C[i] @ Bm @ C[i]) / 2 for i, m in enumerate(members)},
           "member_bare_omega2": {f"{m['file']}[{m['index']}]": float(C[i] @ Ab @ C[i]) / float(C[i] @ Bm @ C[i]) / 2 for i, m in enumerate(members)},
           "ritz_vector_weights_on_orthonormal_basis": c.tolist()}
    json.dump(out, open(os.path.join(IN_DIR, f"ritz_verify_{tag}.json"), "w"), indent=1)
    log(f"verify {tag}: Ritz {lam[0] / 2:.6e} direct {s2_cg / (2 * k):.6e} (rel {out['rel_diff_corrected']:.1e}); lowest member corrected {min(out['member_corrected_omega2'].values()):.6e}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("tags", nargs="*")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--nev", type=int, default=6); ap.add_argument("--warm", default=None); ap.add_argument("--maxiter", type=int, default=400)
    a = ap.parse_args()
    sel = [r for r in RUNGS if (not a.tags or r[0] in a.tags)]
    if a.mode == "regen":
        t, n, L, h = sel[0]
        regen(t, n, L, a.nev, a.warm, a.maxiter)
        return
    if a.mode == "verify":
        for t, n, L, h in sel:
            verify(t, n, L)
        return
    if a.mode in ("prep", "run"):
        ctx = mp.get_context("spawn")
        with ctx.Pool(min(len(sel), a.workers)) as pool:
            recs = pool.map(_prep_worker, [(t, n, L) for t, n, L, h in sel])
        pj = os.path.join(IN_DIR, "ritz_prep.json")
        allrec = {r["tag"]: r for r in (json.load(open(pj)) if os.path.exists(pj) else [])}
        allrec.update({r["tag"]: r for r in recs})
        json.dump(list(allrec.values()), open(pj, "w"), indent=1)
    if a.mode in ("apply", "run"):
        apply_all(a.workers, tags=[r[0] for r in sel] if a.tags else None)
    if a.mode in ("collect", "run"):
        out = collect()
        print(json.dumps({"word": out["word"], "drops": out["drops_vs_auxiliary_relaxed"], "residuals": out["corrected_residuals"], "control_ok": out["control_ok"],
                          "fits": {k: {"intercept_1_over_L": v["fit_in_1_over_L"]["intercept"], "se": v["fit_in_1_over_L"]["intercept_se_fit"], "mesh": v["mesh_difference_(h 1.5 minus h 1.0, same box)"]}
                                   for k, v in out["fits"].items() if "fit_in_1_over_L" in v}}, indent=1))


if __name__ == "__main__":
    main()
