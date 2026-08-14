"""Route (b) production core: two-sided character averaging on S^3.

ROUTE-B ONLY.  This is route (b)'s owned method (section 2 method disjointness:
"Route (b) owns the character-averaging prediction").  No route (a) module may
import it and it imports nothing from route (a).

It supersedes the pilot's LEFT-ACTION-ONLY implementation.  The pilot
`route_b_spectral.py` / `route_b_oneform.py` computed
`(n+1) * (1/|Gamma|) sum_g chi_n(g)`, which cannot express an inhomogeneous
action at all: section 2 shows that writing the action on S^3 in C^2 as
`(z1, z2) -> (zeta^s1 z1, zeta^s2 z2)`, left multiplication alone forces
`s1 = s2` and right multiplication alone forces `s2 = -s1`, so any case with
`s2 != +-s1 (mod q)` needs both factors.

FROZEN FORMULAS, section 6.1, transcribed rather than rederived.

  scalars     L^2(S^3) = sum_n V_n(u) tensor V_n(v)
              invariant dimension  (1/|Gamma|) sum_gamma chi_n(u_gamma) chi_n(v_gamma)

  one-forms   Omega^1 = sum_n V_n(u) tensor [V_n tensor V_2](v)
              V_n tensor V_2 = V_{n+2} + V_n + V_{n-2}
              invariant dimension  (1/|Gamma|) sum_gamma chi_n(u_gamma) chi_m(v_gamma)

  SECTOR, EIGENVALUE AND RANGE (frozen explicitly so the towers cannot be
  interchanged while still totalling correctly).  `n` labels the SCALAR factor
  V_n(u) throughout, never the resulting representation and never the eigenvalue
  index:

      m = n     exact    lambda = n(n+2)    for n >= 1
      m = n+2   coexact  lambda = (n+2)^2   for n >= 0
      m = n-2   coexact  lambda = n^2       for n >= 2

  A summand outside its range is a COMPUTED ZERO CELL wherever the schema
  requires that sector, never a silent omission.

WHY THE LEFT FACTOR CARRIES u AND THE RIGHT CARRIES v.  Section 6.1: "forced by
the section 2 pullback law since the frame is carried by Ad(v) alone".  This
assignment is NOT interchangeable, and no character check can detect swapping it,
which is why section 2 makes the manufactured pointwise pullback test a required
gate rather than a development convenience.

CENTRAL KERNEL.  `[u, v]` and `[-u, -v]` induce the same isometry, and
`chi_n(-u) chi_n(-v) = (-1)^n (-1)^n chi_n(u) chi_n(v)`, so both representatives
contribute the same term.  Averaging must therefore run over the EFFECTIVE group,
one representative per isometry, or the normalization is wrong.  Effective
representatives are selected here by the 4x4 rotation, which is the same surface
G1b canonicalizes on.
"""

import numpy as np

__all__ = ["qmul", "su2_character", "close_pairs_effective",
           "scalar_invariant_dims", "scalar_levels",
           "oneform_levels", "ROUTE", "METHOD"]

ROUTE = "b"
METHOD = "two-sided character averaging over action pairs (section 6.1)"

TOL = 1e-10


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def su2_character(g, n, tol=1e-12):
    """chi_n(g) = sin((n+1) alpha) / sin(alpha) for a unit quaternion of
    half-angle alpha, with the degenerate branches at 0 and pi."""
    a = np.arccos(np.clip(g[0], -1.0, 1.0))
    s = np.sin(a)
    if abs(s) < tol:
        return (n + 1) * (1.0 if abs(a) < tol else (-1.0) ** n)
    return np.sin((n + 1) * a) / s


def _isometry(u, v):
    """The 4x4 rotation of x -> u x v, used only to identify effective elements."""
    w, x, y, z = u
    L = np.array([[w, -x, -y, -z], [x, w, -z, y], [y, z, w, -x], [z, -y, x, w]])
    w, x, y, z = v
    R = np.array([[w, -x, -y, -z], [x, w, z, -y], [y, -z, w, x], [z, y, -x, w]])
    return L @ R


def close_pairs_effective(gen_pairs, cap=2000, tol=TOL):
    """Close the group from RAW ACTION PAIRS, returning one pair per effective element.

    Section 6.1 requires "the group closed independently from the raw pairs".
    Closure is on the pairs; de-duplication is on the induced isometry, so
    `[u, v]` and `[-u, -v]` collapse to a single effective element and the
    average below is normalized by the effective order.
    """
    ident = (np.array([1.0, 0, 0, 0]), np.array([1.0, 0, 0, 0]))
    reps, mats = [ident], [np.eye(4)]

    def seen(M):
        return any(np.abs(M - E).max() < tol for E in mats)

    frontier = [ident]
    while frontier:
        nxt = []
        for (a, b) in frontier:
            for (u, v) in gen_pairs:
                c = (qmul(a, u), qmul(b, v))
                M = _isometry(*c)
                if not seen(M):
                    reps.append(c)
                    mats.append(M)
                    nxt.append(c)
                    if len(reps) > cap:
                        raise ValueError("closure exceeded cap")
        frontier = nxt
    return reps


def scalar_invariant_dims(pairs, nmax):
    """(1/|Gamma|) sum_gamma chi_n(u) chi_n(v), for n = 0..nmax."""
    order = len(pairs)
    return [sum(su2_character(u, n) * su2_character(v, n) for u, v in pairs) / order
            for n in range(nmax + 1)]


def scalar_levels(pairs, nmax):
    """Scalar spectrum: eigenvalue n(n+2), multiplicity dim (V_n x V_n)^Gamma.

    Note the multiplicity is the invariant dimension DIRECTLY.  The pilot's
    left-only code multiplied dim (V_n)^Gamma by (n+1) to account for the second
    factor; the two-sided average already spans both, and multiplying again
    would double-count.
    """
    dims = scalar_invariant_dims(pairs, nmax)
    return [{"n": n, "eigenvalue": n * (n + 2), "multiplicity_raw": d}
            for n, d in enumerate(dims)]


def oneform_levels(pairs, nmax):
    """One-form spectrum under the frozen section 6.1 index map.

    Emits a record for every (n, sector) the schema requires, including summands
    outside their range, which are COMPUTED ZERO CELLS carrying
    `in_range = False` rather than being omitted.
    """
    order = len(pairs)

    def d(n, m):
        return sum(su2_character(u, n) * su2_character(v, m) for u, v in pairs) / order

    rows = []
    for n in range(nmax + 1):
        rows.append({"n": n, "m": n, "sector": "exact",
                     "eigenvalue": n * (n + 2), "in_range": n >= 1,
                     "multiplicity_raw": d(n, n) if n >= 1 else 0.0})
        rows.append({"n": n, "m": n + 2, "sector": "coexact",
                     "eigenvalue": (n + 2) ** 2, "in_range": True,
                     "multiplicity_raw": d(n, n + 2)})
        rows.append({"n": n, "m": n - 2, "sector": "coexact",
                     "eigenvalue": n * n, "in_range": n >= 2,
                     "multiplicity_raw": d(n, n - 2) if n >= 2 else 0.0})
    return rows


# --- binding the averaged group to the group G1b certified -------------------

def averaged_group_digest(pairs, places=10):
    """Canonical digest of the effective group route (b) ACTUALLY averages over.

    Closure here is route (b)'s own (`close_pairs_effective`); only the
    CANONICALIZATION is shared with the gate, so that "identical" means the same
    thing on both sides.  Two independent closures agreeing on this digest is
    element-identity, not merely equal order.
    """
    import hashlib

    def norm(M):
        return tuple(round(float(x), places) + 0.0 for x in np.asarray(M).ravel())

    rows = sorted(norm(_isometry(u, v)) for u, v in pairs)
    payload = ";".join(",".join(f"{x:.{places}f}" for x in r) for r in rows)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def certified_average_group(gen_pairs, g1b_digest):
    """Close, then REFUSE to proceed unless the group matches G1b's certification.

    Without this, route (b) could average over a group that G1b never inspected:
    the gate would certify one canonicalized action while the spectral core
    silently deduplicated a slightly different one.  That is displaced
    verification one layer down, so the binding is structural rather than
    documented.
    """
    reps = close_pairs_effective(gen_pairs)
    digest = averaged_group_digest(reps)
    if digest != g1b_digest:
        raise ValueError(
            "route (b) closed a group that is not the one G1b certified: "
            f"averaged digest {digest[:16]}..., certified {g1b_digest[:16]}...")
    return reps, digest
