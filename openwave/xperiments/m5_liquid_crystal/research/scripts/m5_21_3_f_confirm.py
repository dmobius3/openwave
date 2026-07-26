"""M5.21.3 confirm block (reproducible; supersedes the inline run):
(1) cross-stencil kin re-reads on the s = +1 endpoint;
(2) a fresh n = 24 (h = 2.0) static lift + its kin sign pattern;
(3) THE AUDIT-ADOPTED VARIANT: the conjugation-orbit tangent
    a0_conj = w * (G M + M G^T) (the symmetric, physical velocity of
    M -> L M L^T; the instrument's a0 = w * (G M - M G^T) is the
    antisymmetric variant, audit § C7.4) for the FULL generator
    catalog, both branches: signs must survive, magnitudes are the
    quotable ones.

Run: python3 m5_21_3_f_confirm.py
Out: ../data/m5_21_3_row_confirm.json (extended)
"""
import json
import os
import runpy

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
g4 = runpy.run_path(os.path.join(HERE, "m5_21_3_a_4d.py"),
                    run_name="not_main")

out = {}
path = os.path.join(DATA, "m5_21_3_row_confirm.json")
if os.path.exists(path):
    with open(path) as f:
        out = json.load(f)


def conj_catalog(cfg, M):
    """the conjugation-orbit tangents a0 = w (G M + M G^T), unit
    Frobenius norm (same generators as gen_catalog)."""
    w = g4["envelope"](cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    cat = {}

    def local_rot(vhat):
        W = np.zeros(vhat.shape[:-1] + (4, 4))
        n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
        W[..., 1, 2], W[..., 1, 3] = -n3, n2
        W[..., 2, 1], W[..., 2, 3] = n3, -n1
        W[..., 3, 1], W[..., 3, 2] = -n2, n1
        return W
    cat["clock_local"] = local_rot(V[..., :, 2])
    cat["plane_1d"] = local_rot(V[..., :, 0])
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    Jx = np.zeros((4, 4)); Jx[2, 3], Jx[3, 2] = -1.0, 1.0
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0
    for nm, Gm in (("rot_z", Jz), ("rot_x", Jx),
                   ("boost_z", Kz), ("boost_x", Kx)):
        cat[nm] = np.broadcast_to(Gm, M.shape)
    a0s = {}
    for nm, Gm in cat.items():
        a0 = w * (Gm @ M + M @ Gm.swapaxes(-1, -2))
        nrm = np.sqrt(np.sum(a0 * a0))
        a0s[nm] = a0 / max(nrm, 1e-300)
    return a0s


for tb in ("p1_s+1", "p1_s-1"):
    M, cfg = g4["load_p1"](tb)
    for nm, a0 in conj_catalog(cfg, M).items():
        out[f"kinconj_{tb}_{nm}"] = float(g4["kin_of"](M, a0, cfg))

# (1) cross-stencil re-reads (idempotent refresh)
M, cfg = g4["load_p1"]("p1_s+1")
a0s = g4["gen_catalog"](cfg, M)
for nm in ("clock_local", "boost_z"):
    for st in ("sym", "fwd", "bwd", "2h"):
        out[f"kin_{nm}_{st}"] = float(g4["kin_of"](M, a0s[nm], cfg, st))

# (2) fresh n = 24 lift (skipped if already recorded)
if "kin24_boost_z" not in out:
    cfg24 = g4["base_cfg"](n=24, maxit=4000, tag="p1n24_s+1")
    Z = np.load(g4["SEED3_NPZ"])
    M3 = Z["M"].astype(np.float64)
    idx = np.linspace(0, 31, 24).round().astype(int)
    M3r = M3[np.ix_(idx, idx, idx)]
    M0 = np.zeros((24, 24, 24, 4, 4))
    M0[..., 1:4, 1:4] = M3r
    M0[..., 0, 0] = -cfg24["sg"]
    free = ~g4["pin_shell"](24, cfg24["h"])
    M24, info = g4["fire"](M0, cfg24, free, 4000, tag="p1n24")
    a24 = g4["gen_catalog"](cfg24, M24)
    for nm in ("clock_local", "plane_1d", "rot_z", "boost_z",
               "boost_x"):
        out[f"kin24_{nm}"] = float(g4["kin_of"](M24, a24[nm], cfg24))
    eu, ev = g4["e_parts"](M24, cfg24)
    out["n24_E"] = float(eu + ev)
    out["n24_stop"] = info["stop"]

with open(path, "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps({k: v for k, v in out.items()
                  if k.startswith("kinconj_")}, indent=1))
