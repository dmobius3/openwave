#!/usr/bin/env python3
"""Pre-compute the n=20 SVD and cache it. Run standalone with nohup."""
import sys, os, json, time, hashlib
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'm8_5b'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'm8_5b', 'pilot'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'm8_5b', 'production'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'm8_5b', 'gates'))

RESULTS = os.path.join(os.path.dirname(__file__), 'results')
CACHE_BASIS = os.path.join(RESULTS, 'C_svd_n20.npy')
CACHE_META = os.path.join(RESULTS, 'svd_n20.json')

if os.path.exists(CACHE_BASIS) and os.path.exists(CACHE_META):
    print("n=20 SVD already cached, nothing to do.")
    sys.exit(0)

from p0.group import build_icosians
from route_a_twosided import pairs_left
from route_a_repn import invariant_dim_and_basis

elems = build_icosians()
pairs = pairs_left(elems)

print(f"Starting n=20 SVD at {time.strftime('%H:%M:%S')}...")
print(f"Matrix will be {120*441}x{441} with full_matrices=True")
print(f"Expected U allocation: ~{120*441*120*441*16/1e9:.1f} GB")
sys.stdout.flush()

t0 = time.time()
k, basis, info = invariant_dim_and_basis(pairs, 20)
elapsed = time.time() - t0

print(f"n=20 SVD done in {elapsed:.1f}s: dim={k}, gap={info['gap']}")

h = hashlib.sha256(np.ascontiguousarray(basis).tobytes()).hexdigest()
np.save(CACHE_BASIS, basis)
with open(CACHE_META, 'w') as f:
    json.dump({
        "n": 20, "dim": k,
        "gap": {"gap": info["gap"], "state": info["state"]},
        "basis_shape": list(basis.shape),
        "basis_hash": h,
        "time_s": elapsed,
    }, f, indent=2)

print(f"Cached to {CACHE_BASIS} and {CACHE_META}")
print(f"Hash: {h}")
print("DONE")
