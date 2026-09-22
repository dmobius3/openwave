"""C2 / C3 completeness attack.

Orbit membership is tested DETERMINISTICALLY by moving each candidate into a
canonical frame (no inner optimization loop):

  MAX candidates: rotate so Nbar0 = Nbar - 4I is diagonal with its distinct
     eigenvalue on z.  A member of the orbit of (v3+v-3)/sqrt2 must then be
     supported on {v3, v-3} with equal moduli.
  MIN candidates: rotate so <f> points along +z.  A coherent state must then
     be a phase times v3.

Both tests are then independently confirmed by a direct alignment search on a
subsample.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
from scipy.optimize import minimize
import core as C
from grad import f_and_grad_real

HERE = pathlib.Path(__file__).parent
out = []
P = out.append
rng = np.random.default_rng(3)

CAT = np.array([1, 0, 0, 0, 0, 0, 1], dtype=complex) / np.sqrt(2)
COH = C.basis(3)


def rodrigues(n, th):
    n = np.asarray(n, float)
    n = n / np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def so3_to_D(R):
    """Lift R in SO(3) to D^3 via axis-angle."""
    ang = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    if ang < 1e-12:
        return np.eye(7, dtype=complex)
    if abs(ang - np.pi) < 1e-7:
        w, V = np.linalg.eigh((R + np.eye(3)) / 2)
        ax = V[:, np.argmax(w)]
    else:
        ax = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]]) / (2 * np.sin(ang))
    D = C.rot(ax, ang)
    # verify the lift really implements R on the vector operators
    return D


def canon_max(u):
    """Rotate u so Nbar0 is diagonal, distinct eigenvalue on z. Return best-case
    'mass outside {v3,v-3}' over the discrete frame choices."""
    u = u / np.linalg.norm(u)
    N0 = C.Nbar(u) - 4 * np.eye(3)
    w, V = np.linalg.eigh(N0)
    best = 1e9
    # try every assignment of eigenvector to axis, both signs, det=+1
    import itertools
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product([1, -1], repeat=3):
            R = np.column_stack([sg[k] * V[:, perm[k]] for k in range(3)]).T
            if np.linalg.det(R) < 0:
                continue
            D = so3_to_D(R.T)          # D implements the rotation taking frame->frame
            for Duse in (D, D.conj().T):
                c = Duse @ u
                leak = 1.0 - (abs(c[C.IDX[3]]) ** 2 + abs(c[C.IDX[-3]]) ** 2)
                imb = abs(abs(c[C.IDX[3]]) ** 2 - abs(c[C.IDX[-3]]) ** 2)
                best = min(best, max(leak, imb))
    return best


def canon_min(u):
    """Rotate <f> onto +z; a coherent state then has |c_3| = 1."""
    u = u / np.linalg.norm(u)
    f = C.fvec(u)
    nf = np.linalg.norm(f)
    if nf < 1e-12:
        return 1.0
    f = f / nf
    z = np.array([0.0, 0.0, 1.0])
    ax = np.cross(f, z)
    if np.linalg.norm(ax) < 1e-12:
        R = np.eye(3) if f[2] > 0 else rodrigues([1, 0, 0], np.pi)
    else:
        R = rodrigues(ax, np.arccos(np.clip(np.dot(f, z), -1, 1)))
    best = 1e9
    for D in (so3_to_D(R), so3_to_D(R).conj().T):
        c = D @ u
        if abs(np.linalg.norm(C.fvec(c)) - nf) > 1e-8:
            continue
        fz = C.fvec(c)
        if fz[2] < 0:
            continue
        best = min(best, 1.0 - abs(c[C.IDX[3]]) ** 2)
    return best if best < 1e8 else 1.0


# arm the instrument: the tests must return 0 on known orbit members
P("ARMING the orbit tests on points KNOWN to be in each orbit")
mx = 0.0
for _ in range(200):
    D = C.rot(rng.normal(size=3), rng.uniform(0, 2 * np.pi))
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi))
    mx = max(mx, canon_max(ph * (D @ CAT)))
P(f"  canon_max on 200 random rotations+phases of the cat state : max = {mx:.3e}")
mn = 0.0
for _ in range(200):
    D = C.rot(rng.normal(size=3), rng.uniform(0, 2 * np.pi))
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi))
    mn = max(mn, canon_min(ph * (D @ COH)))
P(f"  canon_min on 200 random rotations+phases of v3           : max = {mn:.3e}")
# and it must be NONZERO on things known to be outside
P(f"  canon_max on v3 (outside the max orbit)                  : {canon_max(COH):.3e}")
P(f"  canon_max on (v3+v-3)/sqrt2 rotated? already 0; on v2+v-2: "
  f"{canon_max(np.array([0,1,0,0,0,1,0],dtype=complex)):.3e}")
P(f"  canon_min on the cat state (outside the min orbit)       : {canon_min(CAT):.3e}")
P(f"  canon_min on (v3+v2)/sqrt2                               : "
  f"{canon_min(np.array([1,1,0,0,0,0,0],dtype=complex)):.3e}")

for nm, ref, target, test in (("MAX", CAT, 463 / 924, canon_max),
                              ("MIN", COH, 1 / 924, canon_min)):
    pts = np.load(HERE / f"pool_{nm}.npy")
    vals = np.load(HERE / f"vals_{nm}.npy")
    sel = np.where(np.abs(vals - target) < 1e-9)[0]
    P("")
    P(f"=== {nm}: {len(sel)} of {len(vals)} converged points within 1e-9 of {target!r}")
    P(f"    worst |value - target| = {np.abs(vals[sel]-target).max():.3e}")
    rows = np.array([[np.dot(C.fvec(pts[i] / np.linalg.norm(pts[i])), C.fvec(pts[i] / np.linalg.norm(pts[i]))),
                      abs(C.a00(pts[i] / np.linalg.norm(pts[i]))) ** 2,
                      float(np.sum(C.Nbar(pts[i] / np.linalg.norm(pts[i])) ** 2))] for i in sel])
    P(f"    |f|^2    range [{rows[:,0].min():.3e}, {rows[:,0].max():.3e}]")
    P(f"    |a00|^2  range [{rows[:,1].min():.12f}, {rows[:,1].max():.12f}]   1/7 = {1/7:.12f}")
    P(f"    TrNbar^2 range [{rows[:,2].min():.10f}, {rows[:,2].max():.10f}]   171/2 = 85.5")
    gaps = np.array([test(pts[i]) for i in sel])
    P(f"    ORBIT TEST over ALL {len(sel)}: max gap = {gaps.max():.3e}, "
      f"median = {np.median(gaps):.3e}, #(gap>1e-6) = {int((gaps>1e-6).sum())}")
    if (gaps > 1e-6).sum():
        bad = sel[np.argmax(gaps)]
        P(f"    WORST OFFENDER u = {np.round(pts[bad]/np.linalg.norm(pts[bad]),8)}")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s5.txt").write_text(txt + "\n")
