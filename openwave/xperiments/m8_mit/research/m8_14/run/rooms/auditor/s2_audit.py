"""Stage 2 audit of the two returns' T4 numerics (R = 1).

Run as `./py s2_audit.py`.  Writes s2_audit.json.  Relative paths only.

The rooms' T4 scripts are copied to lib_a/ and lib_b/ and their assembly routines are
loaded from the copies (Room A's top-level driver is cut off before its main loop).
Comparison values are my stage-1 values from stage1.json, unchanged.

Checks (each prints PASS/FAIL; each one is also run against a planted defect, see PLANTS):
  K1  mesh: level-0 y-nodes equal pi/2 -+ (pi/2)(k/16)^2, 8 w-cells, levels nested by
      midpoint insertion, ndof = (2N-1)(M-1) at every level (fiber row, y = pi row, w = +-W removed)
  K2  seam read from the assembled code path: the nodal field rebuilt from the room's dof map
      satisfies phi(pi,-w_j) = -phi(0,w_j), vanishes on w = +-W and on the fiber; and the
      bottom eigenvector correlates with the closed-form bottom eigenfunction's interpolant
  K3  eigen-residual |K x - lam M x| / (lam |M x|) < 1e-8 for all six pairs at every level
  K4  no eigenvalue missed: Sylvester inertia of K - sigma M (dense LDL^T, levels 0-2) gives
      0 eigenvalues below 0.999999*lam1 and exactly 6 below (lam6 + lam7)/2
  K5  rerun eigenvalues equal the room's JSON eigenvalues (rel 1e-9)
  K6  p, extrapolation and error recomputed from the room's JSON levels by the spec formulas
      equal the room's JSON p, extrap, err (rel 1e-9)
  K7  V3: |room extrap - my stage-1 bottom| / my bottom < 1e-3
"""
import importlib.util
import json
import math
import os
import shutil
import sys

import numpy as np
import scipy.linalg as sla
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
WIDTHS = {"1/4": 0.25, "1/2": 0.5, "1": 1.0, "7/5": 1.4, "3/2": 1.5, "pi/2": math.pi / 2, "17/10": 1.7}
QUICK = "--quick" in sys.argv      # levels 0-2 only (for development)
NLEV = 3 if QUICK else 5

lines = []


def report(name, ok, planted=False):
    tag = "PASS" if ok else "FAIL"
    s = f"{'[planted] ' if planted else ''}{tag}  {name}"
    print(s, flush=True)
    lines.append(s)
    return ok


# ---------------------------------------------------------------- load the rooms' code from copies
def load_rooms():
    for r in ("a", "b"):
        d = os.path.join(HERE, f"lib_{r}")
        os.makedirs(d, exist_ok=True)
        shutil.copy(os.path.join(HERE, f"room_{r}", "t4_fem.py"), os.path.join(d, "t4_fem.py"))
    src = open(os.path.join(HERE, "lib_a", "t4_fem.py")).read()
    src = src[:src.index("results = {}")]          # definitions only, no driver
    A = {"__name__": "room_a_t4", "__file__": os.path.join(HERE, "lib_a", "t4_fem.py")}
    saved = sys.argv
    sys.argv = ["t4_fem.py"]
    exec(compile(src, "lib_a/t4_fem.py", "exec"), A)
    spec = importlib.util.spec_from_file_location("room_b_t4", os.path.join(HERE, "lib_b", "t4_fem.py"))
    B = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(B)
    sys.argv = saved
    return A, B


A, B = load_rooms()


_B_dof_map = B.dof_map


def _B_dof_map_periodic(ys, ws):
    """Room B's dof map with its plant (a) only: periodic seam sign. Its --plant flag would also
    switch the weight to signed cos, which makes M indefinite, so only the seam sign is planted."""
    idx, sgn, n = _B_dof_map(ys, ws)
    sgn[-1, :] = 1.0
    return idx, sgn, n


def set_plant(room, on):
    if room == "a":
        A["plant"] = "seam" if on else None
    else:
        B.dof_map = _B_dof_map_periodic if on else _B_dof_map


def mesh(room, W, lev):
    """y-nodes, w-nodes of the room's level `lev`."""
    if room == "a":
        yy = A["ynodes"](lev)
        ww = np.linspace(-W, W, A["M0"] * 2**lev + 1)
        return yy, ww
    ys, nw = B.base_mesh()
    ws = np.linspace(-W, W, nw + 1)
    for _ in range(lev):
        ys, ws = B.refine(ys), B.refine(ws)
    return ys, ws


def system(room, W, lev):
    """K, M and a function (i, j) -> (dof, sign) for the room's level `lev`."""
    yy, ww = mesh(room, W, lev)
    if room == "a":
        K, M, yy2, ww2, dof, _ = A["assemble"](W, lev)
        assert np.array_equal(yy, yy2) and np.allclose(ww, ww2)
        return K, M, yy, ww, dof
    K, M, idx, sgn = B.assemble(yy, ww)
    return K, M, yy, ww, (lambda i, j: (int(idx[i, j]), float(sgn[i, j]) if idx[i, j] >= 0 else 0.0))


def nodal(vec, yy, ww, dof):
    U = np.zeros((len(yy), len(ww)))
    for i in range(len(yy)):
        for j in range(len(ww)):
            d, s = dof(i, j)
            U[i, j] = s * vec[d] if d >= 0 else 0.0
    return U


def closed_bottom(yv, wv, W):
    nu = math.pi / (2 * W)
    up = yv <= math.pi / 2
    return np.where(up, 1.0, -1.0) * np.abs(np.cos(yv)) ** nu * np.cos(math.pi * wv / (2 * W))


# ---------------------------------------------------------------- inputs
stage1 = json.load(open(os.path.join(HERE, "stage1.json")))["v2"]
res = {"a": json.load(open(os.path.join(HERE, "room_a", "results.json")))["t4"],
       "b": json.load(open(os.path.join(HERE, "room_b", "results.json")))["t4"]}
out = {"checks": {}, "planted": {}, "v3": {}, "residuals": {}, "inertia": {}, "a_vs_b": {}}


# ---------------------------------------------------------------- K1 mesh
def k1(room, W, grading=2.0, planted=False):
    ok = True
    y0, w0 = mesh(room, W, 0)
    k = np.arange(17)
    side = (math.pi / 2) * (k / 16) ** grading
    want = np.concatenate([math.pi / 2 - side[::-1], math.pi / 2 + side[1:]])
    ok &= np.allclose(y0, want, rtol=0, atol=1e-15) and len(w0) == 9 and np.allclose(w0, np.linspace(-W, W, 9))
    prev = y0
    for lev in range(1, 5):
        yy, ww = mesh(room, W, lev)
        ok &= np.array_equal(yy[0::2], prev) and np.allclose(yy[1::2], 0.5 * (prev[:-1] + prev[1:]))
        ok &= len(ww) == 8 * 2**lev + 1
        prev = yy
    nd = [(32 * 2**l - 1) * (8 * 2**l - 1) for l in range(5)]
    ok &= nd == [217, 945, 3937, 16065, 64897]
    if room == "a":
        ok &= res["a"][[l for l, v in WIDTHS.items() if v == W][0]]["ndof"] == nd
    return report(f"K1 room {room} W={W:.6f} mesh: grading (k/16)^2, 8 cells, nested midpoint levels, ndof {nd}",
                  ok, planted)


for room in ("a", "b"):
    for lab, W in WIDTHS.items():
        k1(room, W)
# planted: a cubic grading as the expectation (equivalently, a room that graded cubically) must FAIL
out["planted"]["K1_cubic_grading"] = not k1("a", 1.4, grading=3.0, planted=True)


# ---------------------------------------------------------------- K2 seam from the code path
def k2(room, W, lev=2, planted=False):
    set_plant(room, planted)
    K, M, yy, ww, dof = system(room, W, lev)
    vals, vecs = spla.eigsh(K, k=1, M=M, sigma=0.0, which="LM")
    U = nodal(vecs[:, 0], yy, ww, dof)
    U /= np.max(np.abs(U))
    ifib = (len(yy) - 1) // 2
    seam = np.max(np.abs(U[-1, ::-1] + U[0, :]))            # phi(pi, -w_j) + phi(0, w_j)
    zero = max(np.max(np.abs(U[:, 0])), np.max(np.abs(U[:, -1])), np.max(np.abs(U[ifib, :])))
    Y, Wg = np.meshgrid(yy, ww, indexing="ij")
    C = closed_bottom(Y, Wg, W)
    corr = abs(np.sum(U * C)) / math.sqrt(np.sum(U * U) * np.sum(C * C))
    set_plant(room, False)
    ok = seam < 1e-12 and zero == 0.0 and corr > 0.999
    report(f"K2 room {room} W={W:.6f} level {lev}: seam residual {seam:.1e}, data on arcs/fiber {zero:.1e}, "
           f"corr. with closed-form bottom {corr:.6f}", ok, planted)
    return ok, {"seam_residual": float(seam), "zero_data": float(zero), "corr_closed_form": float(corr)}


for room in ("a", "b"):
    for lab, W in WIDTHS.items():
        ok, info = k2(room, W)
        out["checks"][f"K2_{room}_{lab}"] = info
    okp, info = k2(room, 1.4, planted=True)     # the room's own periodic-seam plant
    out["planted"][f"K2_{room}_periodic_seam"] = {"fired": not okp, **info}

# ---------------------------------------------------------------- K3, K4, K5 per level
for room in ("a", "b"):
    for lab, W in WIDTHS.items():
        rr, ii = [], []
        for lev in range(NLEV):
            K, M, yy, ww, dof = system(room, W, lev)
            vals, vecs = spla.eigsh(K, k=7, M=M, sigma=0.0, which="LM")
            o = np.argsort(vals)
            vals, vecs = vals[o], vecs[:, o]
            r = [float(np.linalg.norm(K @ vecs[:, q] - vals[q] * (M @ vecs[:, q]))
                       / (vals[q] * np.linalg.norm(M @ vecs[:, q]))) for q in range(6)]
            rr.append(max(r))
            jl = res[room][lab]["levels"][lev]
            rel = max(abs(vals[q] - jl[q]) / jl[q] for q in range(6))
            report(f"K3 room {room} W={lab} level {lev}: max eigen-residual {max(r):.1e}", max(r) < 1e-8)
            report(f"K5 room {room} W={lab} level {lev}: rerun vs JSON eigenvalues, max rel {rel:.1e}", rel < 1e-9)
            entry = {"level": lev, "max_residual": max(r), "rerun_vs_json_rel": float(rel),
                     "lam7": float(vals[6])}
            if lev <= 2:
                Kd, Md = K.toarray(), M.toarray()
                def below(sig):
                    _, D, _ = sla.ldl(Kd - sig * Md)   # D is block diagonal (1x1, 2x2), hence tridiagonal
                    ev = sla.eigvalsh_tridiagonal(np.diag(D).copy(), np.diag(D, 1).copy())
                    return int(np.sum(ev < 0))
                n0 = below(0.999999 * vals[0])
                n6 = below(0.5 * (vals[5] + vals[6]))
                ok = n0 == 0 and n6 == 6
                report(f"K4 room {room} W={lab} level {lev}: inertia below 0.999999*lam1: {n0}, "
                       f"below (lam6+lam7)/2: {n6}", ok)
                entry.update({"count_below_lam1": n0, "count_below_mid67": n6})
                if room == "a" and lab == "7/5" and lev == 1:
                    # planted: a solver that missed lam1 (reports lam2..lam7 as the six lowest)
                    m6 = below(0.5 * (vals[6] + vals[6] * 1.01))
                    out["planted"]["K4_missed_lam1"] = not report(
                        f"K4 room a W=7/5 level 1 with lam1 missed: count {m6} vs 6", m6 == 6, True)
                    # planted: residual of a pair whose eigenvalue is shifted by 1e-6 relative
                    x = vecs[:, 0]
                    rp = np.linalg.norm(K @ x - vals[0] * (1 + 1e-6) * (M @ x)) / (vals[0] * np.linalg.norm(M @ x))
                    out["planted"]["K3_shifted_eigenvalue"] = not report(
                        f"K3 room a W=7/5 level 1, eigenvalue shifted 1e-6: residual {rp:.1e}", rp < 1e-8, True)
                    # planted: JSON value altered by 1e-7 relative
                    rel_p = abs(vals[0] - jl[0] * (1 + 1e-7)) / jl[0]
                    out["planted"]["K5_altered_json"] = not report(
                        f"K5 room a W=7/5 level 1, JSON lam1 altered 1e-7: rel {rel_p:.1e}", rel_p < 1e-9, True)
            ii.append(entry)
        out["residuals"][f"{room}_{lab}"] = ii

# ---------------------------------------------------------------- K6 estimates from the JSON levels
def estimates(lv, idx=(2, 3, 4)):
    a, b, c = (lv[i][0] for i in idx)
    p = math.log2((a - b) / (b - c))
    ex = c + (c - b) / (2**p - 1)
    return p, ex, abs(ex - c)


for room in ("a", "b"):
    for lab in WIDTHS:
        d = res[room][lab]
        p, ex, er = estimates(d["levels"])
        rel = max(abs(p - d["p"]) / abs(p), abs(ex - d["extrap"]) / ex, abs(er - d["err"]) / er)
        report(f"K6 room {room} W={lab}: p {p:.6f} extrap {ex:.12f} err {er:.3e} vs JSON, max rel {rel:.1e}",
               rel < 1e-9)
d = res["a"]["7/5"]
p, ex, er = estimates(d["levels"], idx=(1, 2, 3))     # planted: wrong triple of levels
out["planted"]["K6_wrong_levels"] = not report(
    f"K6 room a W=7/5 with levels 1,2,3: extrap {ex:.12f} vs JSON {d['extrap']:.12f}",
    abs(ex - d["extrap"]) / ex < 1e-9, True)

# ---------------------------------------------------------------- K7 / V3
for room in ("a", "b"):
    for lab in WIDTHS:
        mine = stage1[lab]["bottom"]
        theirs = res[room][lab]["extrap"]
        rel = (theirs - mine) / mine
        ok = report(f"K7 V3 room {room} W={lab}: extrap {theirs:.12f} vs mine {mine:.12f}, rel {rel:+.3e}",
                    abs(rel) < 1e-3)
        fin = res[room][lab]["levels"][-1]
        out["v3"][f"{room}_{lab}"] = {
            "room_bottom": theirs, "room_err": res[room][lab]["err"], "room_p": res[room][lab]["p"],
            "mine": mine, "mine_err": stage1[lab]["err"], "rel_diff": rel, "within_1e-3": ok,
            "abs_diff_over_room_err": abs(theirs - mine) / res[room][lab]["err"],
            "finest_rel_vs_mine_lowest6": [(f - m) / m for f, m in zip(fin, stage1[lab]["lowest"])]}
mine = stage1["17/10"]["bottom"]
bad = mine * (1 + 2e-3)                                 # planted: a room bottom off by 2e-3
out["planted"]["K7_offset_2e-3"] = not report(
    f"K7 V3 planted room bottom {bad:.12f} vs mine {mine:.12f}", abs(bad - mine) / mine < 1e-3, True)

# ---------------------------------------------------------------- A vs B
for lab in WIDTHS:
    la, lb = np.array(res["a"][lab]["levels"]), np.array(res["b"][lab]["levels"])
    out["a_vs_b"][lab] = {"max_rel_levels": float(np.max(np.abs(la - lb) / la)),
                          "extrap_rel": (res["a"][lab]["extrap"] - res["b"][lab]["extrap"]) / res["b"][lab]["extrap"]}
    print(f"A vs B W={lab}: max rel over levels {out['a_vs_b'][lab]['max_rel_levels']:.1e}, "
          f"extrap rel {out['a_vs_b'][lab]['extrap_rel']:+.1e}")

nfail = sum(1 for s in lines if s.startswith("FAIL"))
npass = sum(1 for s in lines if s.startswith("PASS"))
fired = {k: (v["fired"] if isinstance(v, dict) else v) for k, v in out["planted"].items()}
print(f"checks: {npass} PASS, {nfail} FAIL; planted defects fired: {fired}")
out["summary"] = {"pass": npass, "fail": nfail, "planted_fired": fired, "lines": lines, "quick": QUICK}
with open(os.path.join(HERE, "s2_audit_quick.json" if QUICK else "s2_audit.json"), "w") as fh:
    json.dump(out, fh, indent=1)
