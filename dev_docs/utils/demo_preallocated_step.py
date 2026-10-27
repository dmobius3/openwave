"""Demonstrate a preallocated NumPy lattice step against the naive one (T8).

The naive step builds every intermediate as a new array (`np.roll`, each operator); the
preallocated step writes the same math into buffers allocated once, with periodic shifts as
slice adds. Both run from the same random field and must agree to 1e-12. The two timings print
for the machine at hand only: see PERFORMANCE_GUIDELINES.md § Write Into Fixed Buffers.

Usage: python dev_docs/utils/demo_preallocated_step.py [N] [steps]   (default 64 50)
"""

import sys
import time

import numpy as np


def lap_naive(f):
    return -2.0 * f.ndim * f + sum(np.roll(f, 1, a) + np.roll(f, -1, a) for a in range(f.ndim))


def lap_into(f, out):
    np.multiply(f, -2.0 * f.ndim, out=out)
    for ax in range(f.ndim):

        def s(sl):
            return tuple(sl if a == ax else slice(None) for a in range(f.ndim))

        lo, hi, first, last = (
            s(slice(None, -1)),
            s(slice(1, None)),
            s(slice(0, 1)),
            s(slice(-1, None)),
        )
        np.add(out[hi], f[lo], out=out[hi])  # out[i] += f[i-1]
        np.add(out[first], f[last], out=out[first])  # wrap-around
        np.add(out[lo], f[hi], out=out[lo])  # out[i] += f[i+1]
        np.add(out[last], f[first], out=out[last])  # wrap-around


def run(u, v, steps, prealloc, dt=0.05):
    lap, tmp = np.empty_like(u), np.empty_like(u)
    for _ in range(steps):
        if prealloc:  # v += dt * (lap(u) - u**3); u += dt * v, with no temporaries
            lap_into(u, lap)
            np.multiply(u, u, out=tmp)
            np.multiply(tmp, u, out=tmp)
            np.subtract(lap, tmp, out=tmp)
            np.multiply(tmp, dt, out=tmp)
            np.add(v, tmp, out=v)
            np.multiply(v, dt, out=tmp)
            np.add(u, tmp, out=u)
        else:
            v += dt * (lap_naive(u) - u**3)
            u += dt * v
    return u, v


n = int(sys.argv[1]) if len(sys.argv) > 1 else 64
steps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
u0 = 0.1 * np.random.default_rng(0).standard_normal((n, n, n))
results = {}
for prealloc in (False, True):
    t0 = time.perf_counter()
    results[prealloc] = run(u0.copy(), np.zeros_like(u0), steps, prealloc)
    print(
        f"{'preallocated' if prealloc else 'naive':>12}: {time.perf_counter() - t0:.2f} s on this machine"
    )
err = max(np.abs(a - b).max() for a, b in zip(results[False], results[True]))
assert err < 1e-12, f"preallocated step disagrees with the naive one: {err:.2e}"
print(f"agreement: max |naive - preallocated| over u and v = {err:.1e} (< 1e-12)")
