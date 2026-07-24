"""
M5.8.2c-2 — headless repro of the _topo_dressed4d GUI explosion (64³, the exact
launcher parameter path, CPU f32). Three variants isolate the culprit:

  (a) V-on + signed flux  — the GUI config (expected to reproduce the explosion)
  (b) V-OFF + signed flux — isolates the signed-flux/static-mask instability
  (c) V-on + pure-Euclid flux (mask forced 0) — isolates the V-force on the
      OFF-MINIMUM dressed seed (the well pins Tr(M_sp²)→1+δ² but the dressing's
      g-leak puts the seed OFF the minimum; at 64³ ldg_c = K·c²/dx⁴ is ~50× the
      24³ headless value — the static V-force scales violently with resolution)

Tracks max|M|, max|M−M_seed|, and the blow-up LOCATION (radius + mask value at
argmax) — the geography tells whether it's the mask boundary, the dressed
shell, or the core.

USAGE:  python m5_8_2c2_gui_repro.py
"""
import os
import sys
import time

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import taichi as ti

ti.init(arch=ti.cpu, default_fp=ti.f32)

from openwave.xperiments.m5_liquid_crystal import medium
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde
from openwave.common import constants

EDGE = 1e-15
N_TARGET = 64**3
STEPS = 200


def smooth_mask(m, passes=3):
    """Separable 6-neighbor box blur — kills the hard signed/Euclid constitutive
    seam (the step-120 core pump: the 0/1 jump at the melt zone drives the core
    through the flux divergence of its signed neighbors)."""
    out = m.astype(np.float32).copy()
    for _ in range(passes):
        for ax in range(3):
            p = np.swapaxes(out, 0, ax)
            q = p.copy()
            q[1:-1] = 0.25 * (p[2:] + p[:-2] + 2.0 * p[1:-1])
            out = np.swapaxes(q, 0, ax)
    return out


def run_variant(tag, ldg_k, force_euclid, smooth=False, steps=STEPS, km=0.0):
    tf = medium.TensorField([EDGE, EDGE, EDGE], N_TARGET, [0.5, 0.5, 0.5], viz_stride=2)
    N = tf.nx
    c = N // 2
    # the exact launcher parameter path
    sim_speed = 1.0
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND * sim_speed
    dt_rs = tf.dx_am * 0.95 / (c_amrs / sim_speed * 3**0.5)
    dt_rs *= 0.5                                          # DT_SCALE_4D
    delta = 0.30
    ldg_c = ldg_k * c_amrs**2 / tf.dx_am**4
    ldg_a = -2.0 * ldg_c * (1.0 + delta * delta)
    r0_vox = 0.06 * N
    rhoc_vox = 3.0
    rw_vox = 0.29 * N
    seeds.seed_dressed_hedgehog_M(tf, c, c, c, r0_vox, rhoc_vox, delta, 0.13, rw_vox, 0.0)
    pde.compute_stable_mask(tf)
    if force_euclid:
        tf.stable_mask.fill(0.0)
    if smooth:
        tf.stable_mask.from_numpy(smooth_mask(tf.stable_mask.to_numpy()))
    pde.compute_tstar(tf)                       # V pins to the dressed seed (2c-2 fix)
    mask = tf.stable_mask.to_numpy()
    M_seed = tf.M_am.to_numpy()
    print(f"\n[{tag}] N={N} dx={tf.dx_am:.4f} dt={dt_rs:.4f} ldg_c={ldg_c:.3e} "
          f"stable={100 * mask.mean():.1f}%")
    n_pl = float(tf.nx * tf.ny + tf.ny * tf.nz + tf.nx * tf.nz)
    t0 = time.time()
    for n in range(steps):
        pde.compute_curvature_flux_4d(tf)
        pde.sample_v03_drift(tf, dt_rs)
        vm = [tf.v03_sums[a_] / n_pl for a_ in range(3)]
        pde.evolve_M_4d(tf, c_amrs, dt_rs, ldg_c, vm[0], vm[1], vm[2], km)
        tf.swap_matrix_buffers()
        if n % 50 == 49 or n < 3:
            M = tf.M_am.to_numpy()
            d = np.abs(M - M_seed).max(axis=(-1, -2))
            mx = float(np.abs(M).max())
            am = np.unravel_index(np.argmax(d), d.shape)
            r_am = np.sqrt(sum((am[i] - c) ** 2 for i in range(3)))
            print(f"   step {n + 1:4d} [{time.time() - t0:5.0f}s] max|M|={mx:10.3e}  "
                  f"max|ΔM|={d.max():10.3e} at r={r_am:5.1f} vox (mask={mask[am]:.0f})")
            if not np.isfinite(mx) or mx > 1e6:
                print(f"   → EXPLODED by step {n + 1}")
                return False, n + 1
    print("   → bounded through the run")
    return True, steps


def main():
    print("=" * 76)
    print("M5.8.2c-2 — GUI explosion repro at 64³ (the exact launcher path, CPU f32)")
    print("=" * 76)
    # Fix-stack history: (a/b/c) culprit = signed flux, core seam → (d/e) smoothing
    # delays only → (f) diagonal inertia km=30 bounds V-OFF 500 ✓ but (g) V-on
    # still pumped slowly @750 (the dressed seed sat OFF the V minimum — constant
    # static push). FINAL STACK: smooth mask + km inertia + V pinned to the
    # DRESSED t*(x) (compute_tstar — now in run_variant for ALL variants):
    # HISTORY: signed flux + cheap inertia grows on every config (h/i @700-750 even
    # with the full heuristic stack) — structural (2b-2/2c-1: only the spectral-
    # projection constrained kernel is stable; Taichi port = follow-up). The GUI
    # v1 SAFE config = SIGNED_FLUX_4D off ⇒ Euclidean flux, time axis live:
    j_ok, j_n = run_variant("j: SAFE GUI v1 (V-on, Euclid flux, t*, km=30), 1500",
                            1.0, True, smooth=True, steps=1500, km=30.0)
    print("\n" + "=" * 76)
    print(f"  j (SAFE GUI v1, 1500): {'bounded' if j_ok else f'EXPLODED @ {j_n}'}")
    print("  → " + ("SAFE config validated — maintainer GUI re-test" if j_ok
          else "even the safe config fails — full stop, deeper diagnosis"))
    print("=" * 76)


if __name__ == "__main__":
    main()
