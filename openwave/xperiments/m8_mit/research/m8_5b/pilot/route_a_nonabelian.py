#!/usr/bin/env python3
"""
M8.5-B engineering pilot, route (a) prototype 3: nonabelian identification map.

NON-EVIDENTIARY pilot code (frozen pre-registration section 6). Tuning set only.
Deck group here is 2T, the binary tetrahedral group, order 24: nonabelian, and NOT
the target group. No 2I-specific quantity is computed anywhere in this file.

WHY THIS RUN EXISTS
  Prototype 2 showed the cyclic identification map collapses to a Fourier selection
  rule, because Hopf coordinates diagonalize the maximal-torus action and a cyclic
  deck group sits inside that torus. A nonabelian group does not, so the constraint
  couples sectors and must be solved rather than selected. This measures that case:
  is it practical, and is it well conditioned?

THE MEASUREMENT, AND WHY IT IS CHARACTER-FREE
  L^2(S^3) at level n carries V_n tensor V_n, dimension (n+1)^2, with the deck group
  acting through the LEFT factor only. So

      dim (level n)^Gamma  =  (n+1) * dim (V_n)^Gamma

  and the pilot obtains dim (V_n)^Gamma as the NULLITY of the stacked operator

      [ rho_n(g_1) - I ;  rho_n(g_2) - I ]

  for generators g_i, computed by singular values. That is a transformation-matrix
  rank on the numerical group action, which is what the frozen section 2 rule permits
  for route (a). No character table, no irreducible labels, no Molien series enters
  the primary measurement.

  A character count is computed SEPARATELY and only as a cross-check, since that is
  route (b)'s method and must never be route (a)'s primary path.
"""

import argparse
import json
import time

import numpy as np

# --------------------------------------------------------------------------
# group construction: explicit unit quaternions, closed numerically
# --------------------------------------------------------------------------

def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def close_group(gens, tol=1e-10, cap=500):
    elems = [np.array([1.0, 0, 0, 0])]

    def known(q):
        return any(np.linalg.norm(q - e) < tol for e in elems)

    frontier = list(elems)
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = qmul(a, g)
                if not known(b):
                    elems.append(b)
                    nxt.append(b)
                    if len(elems) > cap:
                        raise ValueError("closure exceeded cap")
        frontier = nxt
    return elems


def quat_to_su2(q):
    w, x, y, z = q
    return np.array([[w + 1j * x, y + 1j * z],
                     [-y + 1j * z, w - 1j * x]], dtype=complex)


def sym_power(M, n):
    """Matrix of Sym^n(M) on the monomial basis u^(n-k) v^k, k = 0..n."""
    dim = n + 1
    out = np.zeros((dim, dim), dtype=complex)
    # column k: image of u^(n-k) v^k  =  (M00 u + M10 v)^(n-k) (M01 u + M11 v)^k
    for k in range(dim):
        polyA = np.array([1.0 + 0j])
        for _ in range(n - k):
            polyA = np.convolve(polyA, np.array([M[0, 0], M[1, 0]]))
        polyB = np.array([1.0 + 0j])
        for _ in range(k):
            polyB = np.convolve(polyB, np.array([M[0, 1], M[1, 1]]))
        coeffs = np.convolve(polyA, polyB)  # in u^(n) ... v^(n) order
        out[:, k] = coeffs
    return out


# --------------------------------------------------------------------------
# the two measurements
# --------------------------------------------------------------------------

def invariant_dim_by_rank(gens_su2, n, tol_rel=1e-8):
    """Route (a) style: nullity of the stacked (rho_n(g) - I), from singular values."""
    dim = n + 1
    blocks = []
    for M in gens_su2:
        blocks.append(sym_power(M, n) - np.eye(dim))
    A = np.vstack(blocks)
    s = np.linalg.svd(A, compute_uv=False)
    smax = s[0] if s.size else 0.0
    cutoff = max(tol_rel * smax, 1e-12)
    nullity = int(np.sum(s < cutoff))
    # conditioning: gap between the largest "zero" and the smallest "nonzero"
    below = s[s < cutoff]
    above = s[s >= cutoff]
    gap = (above.min() / below.max()) if below.size and above.size else np.inf
    return nullity, float(gap), float(smax)


def invariant_dim_by_character(elems, n):
    """Cross-check ONLY. Route (b)'s method; never route (a)'s primary path."""
    total = 0.0 + 0j
    for q in elems:
        M = quat_to_su2(q)
        total += np.trace(sym_power(M, n))
    return total.real / len(elems)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=14)
    ap.add_argument("--json")
    args = ap.parse_args()

    # 2T generators: a unit and a Hurwitz unit. Order 24, nonabelian, not the target.
    g1 = np.array([0.0, 1.0, 0.0, 0.0])                 # i
    g2 = np.array([0.5, 0.5, 0.5, 0.5])                 # (1+i+j+k)/2
    elems = close_group([g1, g2])
    print(f"deck group: 2T, closed order = {len(elems)} (24 expected)")
    # nonabelian check, done rather than asserted
    A, B = quat_to_su2(g1), quat_to_su2(g2)
    print(f"nonabelian: {not np.allclose(A @ B, B @ A)}")
    print()

    gens_su2 = [quat_to_su2(g1), quat_to_su2(g2)]

    print(f"  {'n':>3} {'dim V_n':>8} {'inv (rank)':>11} {'inv (char)':>11} "
          f"{'agree':>6} {'sv gap':>10} {'sec':>8}")
    rows = []
    t_total = 0.0
    for n in range(args.n_max + 1):
        t0 = time.time()
        nullity, gap, smax = invariant_dim_by_rank(gens_su2, n)
        dt = time.time() - t0
        t_total += dt
        chi = invariant_dim_by_character(elems, n)
        agree = abs(nullity - chi) < 1e-6
        rows.append({"n": n, "dim": n + 1, "invariant_by_rank": nullity,
                     "invariant_by_character": round(chi, 9),
                     "agree": bool(agree), "sv_gap": gap, "seconds": round(dt, 4)})
        print(f"  {n:>3} {n+1:>8} {nullity:>11} {chi:>11.4f} "
              f"{'yes' if agree else 'NO':>6} {gap:>10.2e} {dt:>8.4f}")

    allagree = all(r["agree"] for r in rows)
    worst_gap = min(r["sv_gap"] for r in rows if np.isfinite(r["sv_gap"]))
    print()
    print(f"total time for n = 0..{args.n_max}: {t_total:.3f}s")
    print(f"rank vs character agreement: {'ALL AGREE' if allagree else 'DISAGREEMENT'}")
    print(f"worst singular-value gap (kernel separation): {worst_gap:.2e}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"group": "2T", "order": len(elems),
                       "note": "pilot, non-evidentiary; not the target group",
                       "primary_method": "stacked (rho_n(g)-I) nullity by SVD, character-free",
                       "crosscheck_method": "character average, route (b) style, cross-check only",
                       "total_seconds": round(t_total, 3),
                       "all_agree": bool(allagree),
                       "worst_sv_gap": worst_gap,
                       "rows": rows}, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
