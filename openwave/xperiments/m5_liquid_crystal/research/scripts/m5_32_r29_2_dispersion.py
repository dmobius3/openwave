"""M5.32 R29-2: the S2 dispersion readout about the vacuum (the author's search kit, S2), on the certified stencils.

EQUATIONS FIRST
---------------
The quadratic pencils of R29-0 (c) (same functions, imported from m5_32_r29_0_form.py, never modified), read
over k ON THE LATTICE: the certified stencil set "sym" averages the energies of the forward and backward one-sided
differences, whose plane-wave symbol is |k_lat|^2 = (4/h^2) sin^2(k h / 2) per axis at h = 1.5, so the continuum
limit is read as O(k^2 h^2). For a term L = -(kappa/2) Q(A) - V4 about M_vac the pencil is
    (k_lat^2 kappa Q + H) X = omega^2 kappa Q X,      H = Hess V4 (exact, rank 4),
and the map per term records, over kappa in {+1, -1} (c2 for K_T, kappa_p for K_spec) and k along the three axes
and the (1,1,1) diagonal: the norm sign of every propagating mode, its dispersion omega^2(k) and mass gap, the
helicity content, and whether a pair with helicity-2 content above 0.9 propagates with c_T = c (omega^2 = k_lat^2).
The torsion row (locked variant, R29-0's torsion_qep, validated on Maxwell): the (a1, a2, a3) scan on a grid that
contains the TEGR ratio at both overall signs, mu_K in {0, 30}, both solderings, the same k set.
The record's own -4 I1 and I1_h give the vacuum no quadratic norm (M5.34 R1v): reported, not scanned.

Gates (pre-registered in the task doc): SUCCESS (the author's) = a region of coefficients with no negative norm on
any propagating mode, a massless spin-2 pair with c_T = c in the continuum limit, every other mode gapped or pure
gauge; NULL (the author's) = a map, per term, of which ghost or extra massless mode blocks it; OUR READ = the yes or
no on "the ghost-free region exists only for the TEGR-like torsion combination".

Run from this directory:  python3 m5_32_r29_2_dispersion.py  ->  ../data/m5_32_r29_2_dispersion.json  (minutes)
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r29_2_dispersion.json")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# import the R29-0 module WITHOUT running its main body: exec only the definitions (everything above "(a) sign convention")
src = open(os.path.join(HERE, "m5_32_r29_0_form.py")).read()
head = src[: src.index('log("(a) sign convention")')]
# the torsion solver and the helicity helpers live further down; pull their definitions too
def _grab(name, src):
    i = src.index(f"def {name}(")
    j = src.find("\n\n\n", i)
    return src[i:j] + "\n"
head += "\n" + _grab("helicity_projectors", src) + _grab("norm_vec10", src) + _grab("content", src) + _grab("pencil_rows", src)
head += "\n" + _grab("torsion_qep", src) + _grab("hel16", src)
F = {"__file__": os.path.join(HERE, "m5_32_r29_0_form.py"), "__name__": "m5_32_r29_0_defs"}
exec(compile(head, "m5_32_r29_0_form.py(defs)", "exec"), F)
ETA, BX, BN, W1 = F["ETA"], F["BX"], F["BN"], F["W1"]
quad_matrix, q_T, q_cov, q_h, q_spec, hess_v4, h_of_vac, torsion_qep = (F[k] for k in ("quad_matrix", "q_T", "q_cov", "q_h", "q_spec", "hess_v4", "h_of_vac", "torsion_qep"))
sla = F["sla"]

H_LAT = 1.5
VACUA = {
    "toy (32, 1, 0.3, 0), w W1": dict(roots=(-32.0, 1.0, 0.3, 0.0), M=np.diag([32.0, 1.0, 0.3, 0.0]), w=W1),
    "R23 (8, 1, 0.3, 0.3), w 25 W1": dict(roots=(-8.0, 1.0, 0.3, 0.3), M=np.diag([8.0, 1.0, 0.3, 0.3]), w=25.0 * W1),
}
KHATS = {"x1": (1.0, 0.0, 0.0), "x2": (0.0, 1.0, 0.0), "x3": (0.0, 0.0, 1.0), "diag": tuple(np.ones(3) / np.sqrt(3))}
K_CONT = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]            # continuum k; the lattice symbol is applied per axis


def k_lat2(k, khat):
    return float(sum((2.0 / H_LAT * np.sin(k * c * H_LAT / 2.0)) ** 2 for c in khat))


def helicity_projectors(khat):
    return F["helicity_projectors"](khat)


def mfield_map(vname, v):
    M, roots, w, hm = v["M"], v["roots"], v["w"], h_of_vac(v["M"])
    H = hess_v4(M, roots, w)
    terms = {
        "K_T": quad_matrix(lambda X: q_T(X, hm)),
        "K_cov": quad_matrix(q_cov),
        "K_h": quad_matrix(lambda X: q_h(X, hm)),
        "K_spec": 2.0 * quad_matrix(lambda X: q_spec(X, M)),
    }
    out = {"HessV4_eigenvalues": np.linalg.eigvalsh(H).tolist()}
    for lab, Q in terms.items():
        for kappa in (1.0, -1.0):
            rec = {"coefficient": kappa, "per_k": {}}
            n_static = 10 - int(np.sum(np.abs(np.linalg.eigvalsh(Q)) > 1e-12))
            rec["static_directions"] = n_static
            worst_neg, massless_h2 = 0, []
            for kn, khat in KHATS.items():
                P = helicity_projectors(khat)
                rows = []
                for k in K_CONT:
                    kl2 = k_lat2(k, khat)
                    w_, V = sla.eig(kl2 * kappa * Q + H, kappa * Q)
                    modes = []
                    for i in range(10):
                        if not np.isfinite(w_[i]) or abs(w_[i]) > 1e12:
                            continue
                        x = V[:, i]
                        nrm = float(np.real(np.conj(x) @ (kappa * Q) @ x))
                        c = F["content"](x, P)
                        modes.append({"omega2": float(np.real(w_[i])), "omega2_imag": float(np.imag(w_[i])),
                                      "norm_sign": int(np.sign(nrm)) if abs(nrm) > 1e-12 else 0, "h": c, "dominant": BN[int(np.argmax(np.abs(x)))]})
                    rows.append({"k": k, "k_lat2": kl2, "modes": modes})
                    worst_neg = max(worst_neg, sum(1 for m in modes if m["norm_sign"] < 0))
                    for m in modes:
                        if m["h"]["|h|=2"] > 0.9 and abs(m["omega2"] - kl2) < 1e-6 * (1 + kl2):
                            massless_h2.append((kn, k, m["norm_sign"]))
                rec["per_k"][kn] = rows
            rec["max_negative_norm_modes_over_k"] = worst_neg
            rec["massless_helicity2_modes_(khat, k, sign)"] = massless_h2[:12]
            rec["complex_omega2_anywhere"] = bool(any(abs(m["omega2_imag"]) > 1e-9 for kn in rec["per_k"] for r in rec["per_k"][kn] for m in r["modes"]))
            out[f"{lab} (coef {kappa:+g})"] = rec
    out["-4 I1, I1_h"] = "no quadratic norm about the uniform vacuum (M5.34 R1v); not scanned"
    return out


def torsion_map():
    g_t, d_t = 32.0, 0.3
    H_M = hess_v4(np.diag([g_t, 1.0, d_t, 0.0]), (-g_t, 1.0, d_t, 0.0), W1)
    pts = {}
    base = [(0.25, 0.5, -1.0), (0.25, 0.5, 0.0), (0.25, 0.0, -1.0), (0.25, 0.0, 0.0), (0.25, 0.25, -1.0), (0.25, 1.0, -1.0), (0.25, 0.5, -2.0), (0.25, 0.5, -0.5), (0.0, 0.5, -1.0)]
    for s in (+1.0, -1.0):
        for (a1, a2, a3) in base:
            pts[f"({s * a1:+g}, {s * a2:+g}, {s * a3:+g})"] = (s * a1, s * a2, s * a3)
    out = {"points": {}, "setup": "E0 = diag(sqrt 32, 1, sqrt 0.3, 0), V4 on, k continuum values (the torsion term is written in the continuum at the form level; the lattice symbol applies as for the M-field terms, O(k^2 h^2))"}
    for pname, (a1, a2, a3) in pts.items():
        rec = {}
        for mu_K in (0.0, 30.0):
            for sold in ("rank3_pseudoinverse", "null_leg_dropped"):
                worst_neg, n_cplx, h2 = 0, 0, []
                n_unst, max_growth = 0, 0.0
                per = {}
                for kn, khat in KHATS.items():
                    for k in (0.25, 0.5, 1.0):
                        r = torsion_qep(g_t, d_t, a1, a2, a3, mu_K, sold, khat, k, H_M)
                        real = [m for m in r["modes"] if abs(m["omega"][1]) < 1e-9]
                        cplx = [m for m in r["modes"] if abs(m["omega"][1]) >= 1e-9]
                        neg = sum(1 for m in real if m["norm_sign"] < 0)
                        unstable = [m for m in cplx if abs(m["omega"][1]) > 1e-2]          # a growth rate above 1e-2: an instability; below: a marginal near-zero mode
                        worst_neg = max(worst_neg, neg); n_cplx += len(cplx); n_unst += len(unstable)
                        max_growth = max(max_growth, max((abs(m["omega"][1]) for m in cplx), default=0.0))
                        for m in real:
                            if m["helicity_content"]["|h|=2"] > 0.9:
                                h2.append((kn, k, round(m["omega2"], 5), m["norm_sign"]))
                        per[f"{kn} k {k:g}"] = {"real": len(real), "complex": len(cplx), "unstable_(growth>1e-2)": len(unstable), "negative_norm": neg,
                                                 "omega2_sign": [(round(m["omega2"], 5), m["norm_sign"]) for m in real],
                                                 "complex_omega": [(round(m["omega"][0], 5), round(m["omega"][1], 5)) for m in cplx[:6]]}
                rec[f"mu_K {mu_K:g}, {sold}"] = {"max_negative_norm_over_k": worst_neg, "complex_modes_total": n_cplx, "unstable_modes_total": n_unst, "max_growth_rate": max_growth,
                                                   "helicity2_modes_(khat,k,omega2,sign)": h2[:10], "per_k": per}
        out["points"][pname] = rec
        log(f"torsion {pname} done")
    return out


out = {"versions": {"python": sys.version.split()[0], "numpy": np.__version__}}
log("M-field maps")
out["M_field"] = {vn: mfield_map(vn, v) for vn, v in VACUA.items()}
log("torsion map")
out["torsion_locked_variant"] = torsion_map()
# the summary map and the pre-registered read
summ = {}
for vn, rows in out["M_field"].items():
    summ[vn] = {k: {"max_negative_norm_modes": r["max_negative_norm_modes_over_k"], "static": r["static_directions"],
                    "massless_helicity2": len(r["massless_helicity2_modes_(khat, k, sign)"]), "complex": r["complex_omega2_anywhere"]}
                for k, r in rows.items() if isinstance(r, dict)}
tsum = {}
for pname, rec in out["torsion_locked_variant"]["points"].items():
    tsum[pname] = {k: {"max_negative_norm": r["max_negative_norm_over_k"], "complex": r["complex_modes_total"], "unstable": r["unstable_modes_total"], "max_growth": r["max_growth_rate"], "helicity2_modes": len(r["helicity2_modes_(khat,k,omega2,sign)"])}
                   for k, r in rec.items()}
ghost_free = [p for p, rec in tsum.items() if any(v["max_negative_norm"] == 0 and v["unstable"] == 0 for v in rec.values())]
ghost_free_strict = [p for p, rec in tsum.items() if any(v["max_negative_norm"] == 0 and v["complex"] == 0 for v in rec.values())]
out["summary"] = {"M_field": summ, "torsion": tsum, "torsion_points_with_a_ghost_free_reading_(no_negative_norm,_no_unstable_mode)": ghost_free,
                  "torsion_points_strict_(no_complex_mode_at_all)": ghost_free_strict,
                  "read": ("the author's SUCCESS needs no negative norm on any propagating mode AND a massless helicity-2 pair at c_T = c; "
                           "OUR READ is the yes or no on 'the ghost-free region exists only for the TEGR-like torsion combination': "
                           "yes if a ghost-free torsion point exists and every such point is TEGR-like, no otherwise; written at the results post from this map")}
out["wall_s"] = round(time.time() - T0, 1)
json.dump(out, open(OUT, "w"), indent=1, default=float)
log(f"wrote {os.path.relpath(OUT, HERE)}")
