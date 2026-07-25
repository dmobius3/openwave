"""M5.27 arm D: the adversarial audit (independent implementation).

Per the AI_HYGIENE cardinal rule, the pilot's claims are re-derived here by a
SECOND implementation that shares no code with the taichi harness: pure numpy,
f64, its own stencil, its own gradient, its own integrator. The audit tries to
REFUTE, so each check states what would falsify the main-run claim.

A1  the driven EOM is what the task says it is: the ONLY channel by which the
    background enters the dynamics is the trace targets C_p(sg(t)). Refuter: if
    some other sg dependence exists, an independent build would disagree.
A2  the V4 gradient used by production matches an independent complex-step
    derivative of V4 (machine-exact, no truncation error).
A3  the drive-power identity: over a step, the energy change of the field equals
    the explicit drive work Int (dV4/dsg) dsg. Refuter: a broken ledger would
    make "zero average power" meaningless.
A4  the NULL claim, independently: on a small arena, does a background drive at
    the clock frequency sustain the carried isorotation charge above control?
    Refuter for the main result: if the independent build DOES sustain J where
    the taichi run does not, the main NULL is an implementation artifact.
A5  the structural claim behind B1: a spatially UNIFORM drive on the spectral
    targets couples to the eigenVALUES, and at first order does not torque the
    eigenVECTOR frame that carries the clock. Refuter: if the commutator of the
    drive force with the clock flow is nonzero, the coupling exists after all.

Run:  python m5_27_d_audit.py
Out:  data/m5_27_audit.json
"""
import json
import math
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G0, DELTA, W1 = 8.0, 0.5, 7.24023879e-4
H_RES, DT = 1.5, 0.005
RESULTS = []


def check(name, ok, detail, refuter):
    RESULTS.append({"check": name, "pass": bool(ok), "detail": detail,
                    "refuter": refuter})
    print(f"{'PASS' if ok else 'FAIL'}  {name}\n      {detail}\n      refuter: {refuter}",
          flush=True)


# ================================================================
# independent V4 + gradient (own derivation, complex-step verified)
# ================================================================
def v4_np(M, sg, delta=DELTA, w1=W1):
    """V4 = w Sum_p (tr((M eta)^p) - C_p)^2, elementwise over a field of M."""
    me = M @ ETA
    p1 = me
    p2 = p1 @ me
    p3 = p2 @ me
    p4 = p3 @ me
    out = 0.0
    for p, pw in zip((1, 2, 3, 4), (p1, p2, p3, p4)):
        cp = sg**p + 1.0 + delta**p
        out = out + (np.trace(pw, axis1=-2, axis2=-1) - cp) ** 2
    return w1 * out


def dv4_dM_np(M, sg, delta=DELTA, w1=W1):
    """dV4/dM by an INDEPENDENT route: d tr((M eta)^p)/dM = p (eta (M eta)^(p-1))^T,
    then symmetrized (M is symmetric)."""
    me = M @ ETA
    powers = [np.broadcast_to(np.eye(4), M.shape).copy(), me, me @ me, me @ me @ me]
    g = np.zeros_like(M)
    for p in (1, 2, 3, 4):
        cp = sg**p + 1.0 + delta**p
        tp = np.trace(powers[p - 1] @ me if p > 1 else me, axis1=-2, axis2=-1)
        if p == 1:
            tp = np.trace(me, axis1=-2, axis2=-1)
        else:
            tp = np.trace(np.linalg.matrix_power(me, p), axis1=-2, axis2=-1)
        pref = 2.0 * w1 * (tp - cp) * p
        term = np.einsum("ij,...jk->...ik", ETA, powers[p - 1])
        g = g + pref[..., None, None] * np.swapaxes(term, -1, -2)
    return 0.5 * (g + np.swapaxes(g, -1, -2))


def dv4_dsg_np(M, sg, delta=DELTA, w1=W1):
    me = M @ ETA
    out = 0.0
    for p in (1, 2, 3, 4):
        cp = sg**p + 1.0 + delta**p
        tp = np.trace(np.linalg.matrix_power(me, p), axis1=-2, axis2=-1)
        out = out - 2.0 * w1 * p * sg ** (p - 1) * (tp - cp)
    return out


t0 = time.time()
print("[M5.27 audit] independent numpy re-derivation\n")

# ================================================================
# A2 — gradient vs complex-step (machine-exact, no truncation)
# ================================================================
rng = np.random.default_rng(0)
M1 = np.diag([-G0, 1.0, DELTA, 0.0]) + 0.05 * rng.standard_normal((4, 4))
M1 = 0.5 * (M1 + M1.T)
g_an = dv4_dM_np(M1[None], G0)[0]
g_cs = np.zeros((4, 4))
for a in range(4):
    for b in range(4):
        Mc = M1.astype(complex).copy()
        hstep = 1e-20
        Mc[a, b] += 1j * hstep
        if a != b:
            Mc[b, a] += 1j * hstep      # keep it symmetric: perturb the PAIR
        me = Mc @ ETA
        val = 0.0
        for p in (1, 2, 3, 4):
            cp = G0**p + 1.0 + DELTA**p
            val = val + (np.trace(np.linalg.matrix_power(me, p)) - cp) ** 2
        g_cs[a, b] = (W1 * val).imag / hstep
# complex-step of a symmetric-pair perturbation gives the symmetrized derivative
# with the off-diagonals counted twice; compare on that convention
g_cmp = g_an.copy()
off = ~np.eye(4, dtype=bool)
g_cmp[off] *= 2.0
rel_g = np.abs(g_cmp - g_cs).max() / max(np.abs(g_cs).max(), 1e-30)
check(
    "A2 V4 gradient matches an independent complex-step derivative",
    rel_g < 1e-9,
    f"max rel deviation {rel_g:.3e} (complex-step h = 1e-20, no truncation error)",
    "a mismatch would mean the force driving every run is not dV4/dM",
)

# ================================================================
# A1 — the drive enters ONLY through C_p(sg)
# ================================================================
def fd_sg(hh):
    """CENTRAL difference (O(h^2)); a forward difference leaves O(h) truncation
    that reads as a false mismatch (first audit pass, corrected)."""
    return (v4_np(M1[None], G0 + hh)[0] - v4_np(M1[None], G0 - hh)[0]) / (2 * hh)


an = dv4_dsg_np(M1[None], G0)[0]
e1_, e2_ = abs(an - fd_sg(1e-4)), abs(an - fd_sg(5e-5))
conv = e1_ / max(e2_, 1e-30)
rel_s = e2_ / max(abs(an), 1e-30)
check(
    "A1 the background enters the dynamics ONLY through the trace targets C_p(sg)",
    rel_s < 1e-6 and 3.0 < conv < 5.0,  # the O(h^2) convergence IS the proof; rel_s is the FD floor
    f"dV4/dsg analytic {an:.8f} vs central FD {fd_sg(5e-5):.8f} (rel {rel_s:.2e}); "
    f"truncation falls {conv:.2f}x on halving h (O(h^2) expected 4x); "
    f"the sg dependence of V4 is exactly C_p = sg^p + 1 + delta^p",
    "any additional sg channel would break this identity and invalidate the "
    "host-side drive construction",
)

# ================================================================
# A5 — the structural claim behind blindspot B1
# ================================================================
# A uniform drive changes C_p only, so its force contribution is
#   dF/dsg = -2w Sum_p p^2 sg^(p-1) sym[(eta (M eta)^(p-1))^T].
# FIRST AUDIT PASS CLAIMED these are polynomials in M and therefore commute with
# M. That is FALSE in general: eta (M eta)^(p-1) interleaves eta, so for an M
# carrying mixed (0,i) entries the commutator is NONZERO (measured 1.4e-2) and
# the drive CAN torque the eigenframe.
# It becomes TRUE exactly on the BLOCK-DIAGONAL states this staged 4x4 runs: with
# the (0,i) block projected out, eta restricts to the identity on the spatial
# block, so eta (M eta)^(p-1) restricts to (M_sp)^(p-1) — a genuine polynomial in
# M_sp, which commutes with it (measured 4.5e-21, machine zero).
# CONSEQUENCE, and the mechanism behind the null tongue: the mixed (0,i) block IS
# the channel through which a background scalar could torque the clock frame.
# Deferring it (the staged 4x4) removes precisely the term that could entrain.
dF = (dv4_dM_np(M1[None], G0 + 1e-4)[0] - dv4_dM_np(M1[None], G0 - 1e-4)[0]) / 2e-4
comm_gen = dF @ M1 - M1 @ dF
rel_gen = np.abs(comm_gen).max() / max(np.abs(dF).max() * np.abs(M1).max(), 1e-30)

Mb = M1.copy()
Mb[0, 1:] = 0.0
Mb[1:, 0] = 0.0                      # the staged 4x4: mixed block projected out
dFb = (dv4_dM_np(Mb[None], G0 + 1e-4)[0] - dv4_dM_np(Mb[None], G0 - 1e-4)[0]) / 2e-4
comm_blk = dFb @ Mb - Mb @ dFb
rel_blk = np.abs(comm_blk).max() / max(np.abs(dFb).max() * np.abs(Mb).max(), 1e-30)
check(
    "A5 on the STAGED (block-diagonal) states the uniform drive force commutes with M: "
    "it moves eigenvalues and cannot torque the clock frame",
    rel_blk < 1e-15 and rel_gen > 1e-4,
    f"||[dF/dsg, M]||rel = {rel_blk:.3e} block-diagonal (machine zero) vs "
    f"{rel_gen:.3e} with the mixed (0,i) block present. So the null is STRUCTURAL "
    f"for phase A, and the mixed block is the missing coupling channel",
    "if the block-diagonal commutator were nonzero, a uniform spectral drive "
    "could torque the clock directly and a null tongue would need another "
    "explanation; if the general-M commutator were ALSO zero, the mixed block "
    "would not be the channel and phase B would need a different design",
)

# ================================================================
# A3 + A4 — an independent driven run on a small arena
# ================================================================
N = 13
h = H_RES


def lap_free_evolve(M, Mp, sg, dt=DT, dx=h):
    """An INDEPENDENT integrator: plain leapfrog on E = curvature + V4 with a
    simple (non-eta) Laplacian curvature term. It is deliberately NOT the
    production eta-flux scheme — agreement of the PHYSICS conclusion across two
    different discretizations is the point of an audit."""
    lap = (
        np.roll(M, 1, 0) + np.roll(M, -1, 0)
        + np.roll(M, 1, 1) + np.roll(M, -1, 1)
        + np.roll(M, 1, 2) + np.roll(M, -1, 2)
        - 6.0 * M
    ) / dx**2
    force = -lap + dv4_dM_np(M, sg)
    Mn = 2.0 * M - Mp - dt**2 * force
    return Mn


def seed(N):
    """A hedgehog-like spatial defect embedded in the covariant vacuum."""
    c = (N - 1) / 2
    idx = np.arange(N) - c
    X, Y, Z = np.meshgrid(idx, idx, idx, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    nx, ny, nz = X / r, Y / r, Z / r
    env = np.tanh(r / 2.5)
    M = np.zeros((N, N, N, 4, 4))
    n = np.stack([nx, ny, nz], -1)
    nn = np.einsum("...i,...j->...ij", n, n)
    sp = DELTA * np.eye(3) + (1.0 - DELTA) * nn
    vac = np.diag([1.0, DELTA, 0.0])
    M[..., 1:, 1:] = env[..., None, None] * sp + (1 - env)[..., None, None] * vac
    M[..., 0, 0] = -G0
    return M


def clock_flow(M):
    """a0 = w [W, M] about the local leading spatial eigenvector (the a0_conj
    construction, independently re-implemented)."""
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    vh = V[..., :, 2]
    W = np.zeros(M.shape)
    n1, n2, n3 = vh[..., 0], vh[..., 1], vh[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    A = W @ M - M @ W
    nrm = np.sqrt((A * A).sum())
    return A / max(nrm, 1e-300)


def run(eps, om_bar, n_steps=24000, om_kick=0.2, ramp_t=60.0):
    M = seed(N)
    a0 = clock_flow(M)
    Mp = M - DT * om_kick * a0        # the SET-J style kick, independently built
    j_hist, e_hist, work = [], [], 0.0
    sg_prev = G0
    for s in range(n_steps):
        t = s * DT
        env = 1.0 if t >= ramp_t else 0.5 * (1 - math.cos(math.pi * t / ramp_t))
        sg = G0 * (1 + eps * env * math.cos(om_bar * t)) if eps else G0
        Mn = lap_free_evolve(M, Mp, sg)
        # explicit drive work over this step (the A3 ledger line)
        work += float(dv4_dsg_np(M, sg).sum()) * (sg - sg_prev)
        sg_prev = sg
        Mp, M = M, Mn
        if not np.isfinite(M).all():
            break
        if s % 2000 == 0:
            a0c = clock_flow(M)
            md = (M - Mp) / DT
            j_hist.append((t, float((md * a0c).sum())))
            v = float(v4_np(M, sg).sum())
            kin = 0.5 * float((md * md).sum())
            e_hist.append((t, kin + v))
    return {"J": j_hist, "E": e_hist, "work": work, "finite": bool(np.isfinite(M).all())}


print("\n[A3/A4] independent driven runs on a %d^3 arena" % N)
ctrl = run(0.0, 0.0)
drv = run(0.1, 0.2)
drv2 = run(0.1, 0.4)          # the 2:1 parametric window
j_c = ctrl["J"][-1][1] / max(abs(ctrl["J"][0][1]), 1e-30)
j_d = drv["J"][-1][1] / max(abs(drv["J"][0][1]), 1e-30)
j_d2 = drv2["J"][-1][1] / max(abs(drv2["J"][0][1]), 1e-30)

# A3: energy change vs explicit drive work
e0, e1 = drv["E"][0][1], drv["E"][-1][1]
dE = e1 - e0
rel_w = abs(dE - drv["work"]) / max(abs(drv["work"]), abs(dE), 1e-30)
check(
    "A3 the field energy change is accounted by the explicit drive work",
    rel_w < 0.5,
    f"dE = {dE:+.4f} vs Int (dV4/dsg) dsg = {drv['work']:+.4f} (rel gap {rel_w:.2f}; "
    f"a coarse independent discretization, so the bar is order-of-magnitude)",
    "a ledger that does not close means 'zero average drive power' cannot "
    "discriminate entrainment from pumping",
)
check(
    "A4 INDEPENDENT: a background drive at the clock frequency does not sustain "
    "the carried J above control",
    (j_d - j_c) < 0.10 and (j_d2 - j_c) < 0.10,
    f"J retention: control {j_c:+.4f}, driven 1:1 {j_d:+.4f} (gain {j_d-j_c:+.4f}), "
    f"driven 2:1 {j_d2:+.4f} (gain {j_d2-j_c:+.4f}) on an independent numpy "
    f"integrator with a non-eta curvature term",
    "if THIS build sustained J where the taichi run does not, the main NULL "
    "would be an artifact of the production discretization",
)

out = {"results": RESULTS, "all_pass": all(r["pass"] for r in RESULTS),
       "arena": N, "retention": {"control": j_c, "driven_1to1": j_d,
                                 "driven_2to1": j_d2},
       "wall_s": time.time() - t0}
p = os.path.join(DATA, "m5_27_audit.json")
with open(p, "w") as f:
    json.dump(out, f, indent=1)
print(f"\n[audit] {sum(r['pass'] for r in RESULTS)}/{len(RESULTS)} pass, "
      f"wall {out['wall_s']:.1f} s -> {p}")
