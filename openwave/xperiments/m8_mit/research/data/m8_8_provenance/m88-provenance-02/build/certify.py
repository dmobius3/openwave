"""DETERMINISTIC, SELF-CONTAINED certificate that  im(d3) = ker(d2)  integrally.

Reproducibility: no seed, no random choice. The column set is the pivot set of a
fixed elimination at a fixed prime; the dropped rows are 0, 1, 2 in order.
Re-running reproduces identical output and identical JSON.

The argument, with every premise established in this file:

  (a) CONTAINMENT.  d3.d2 = 0 exactly over Z[2I], so im(d3) is inside ker(d2).
      d2.d1 = 0 and eps.d1 = 0 are checked too, since the rank ceilings below
      lean on them.

  (b) EXACT RANKS, no Smith form at size 240.  For any prime p, rank_p(M) is a
      LOWER bound on rank_Q(M). The chain relations give matching UPPER bounds:
        rank d1 <= 119   because eps.d1 = 0 puts im(d1) inside ker(eps), rank 119
        rank d2 <= 240 - rank d1   because d2.d1 = 0
        rank d3 <= 240 - rank d2   because d3.d2 = 0
      Measuring rank_p >= 119, 121, 119 pins all three to equality, hence
        rank ker(d2) = 240 - rank d2 = 119.

  (c) SATURATION.  For an integer matrix of rank r, the gcd of ALL r x r minors
      is the product of the elementary divisors, and it divides the gcd of any
      SUBSET. Exhibiting one maximal minor of determinant +-1 forces every
      elementary divisor to 1, so im(d3) is primitive in C_2.

  (a)+(b)+(c): im(d3) and ker(d2) are both saturated rank-119 sublattices of
  C_2 and one contains the other, so they are equal. A saturated sublattice
  cannot properly contain another of the same rank, since the quotient would be
  torsion inside a torsion-free group.

No prime list is load-bearing: p enters only as a lower bound, and the bound is
closed from above by the chain relations.
"""
import sys, json, hashlib, pathlib, functools; sys.path.insert(0, '.')
from math import gcd
import numpy as np
from complex import *
import sat

# The guard below protects THIS FILE'S OWN numpy elimination (rank_p, further down),
# which works in int64: its update forms factor * pivot_row with both entries below p,
# so it is exact only while (p-1)^2 < 2^63. Past that the products wrap SILENTLY and
# the routine returns a WRONG rank instead of raising. Refuse rather than absorb.
# Bound adopted from the maintainer audit (af178091).
#
# For the record, since an earlier note here claimed otherwise: sat.py needs no guard.
# It is pure-Python arbitrary-precision arithmetic and cannot overflow at any prime.
# The int64 exposure in this codebase is certify.py's rank_p and nothing else.
MAX_INT64_PRIME = 3_037_000_499
P = 10**9 + 7
if P > MAX_INT64_PRIME:
    raise ValueError(f"p={P} exceeds the int64-safe bound {MAX_INT64_PRIME}; "
                     "the elimination would overflow and return a wrong rank")
F = pathlib.Path('m8_8_construction_packet.json')
PK = json.loads(F.read_text())
dec = lambda terms: {g: c for c, g in terms}
D = {n: [[dec(e) for e in row] for row in PK["boundary_maps"][n]] for n in ("d1","d2","d3")}
fail = []
measured = {}          # every check's ACTUAL outcome, keyed for the emitted certificate
def chk(label, cond, detail="", key=None):
    cond = bool(cond)
    print(f"  [{'PASS' if cond else '**FAIL**'}] {label}" + (f"   {detail}" if detail else ""))
    if not cond: fail.append(label)
    if key: measured[key] = cond
    return cond

print(f"packet   {F.name}")
print(f"sha256   {hashlib.sha256(F.read_bytes()).hexdigest()}\n")

print("(a) containment and the chain relations, exact over Z[2I]")
def mm(A, B):
    return [[functools.reduce(radd, (rmul(A[i][k], B[k][j]) for k in range(len(B))), ZR)
             for j in range(len(B[0]))] for i in range(len(A))]
chk("d3 . d2 = 0, so im(d3) is contained in ker(d2)",
    all(not e for r in mm(D["d3"], D["d2"]) for e in r), key="d3_d2_zero")
chk("d2 . d1 = 0", all(not e for r in mm(D["d2"], D["d1"]) for e in r), key="d2_d1_zero")
chk("eps . d1 = 0", not any(sum(e.values()) for r in D["d1"] for e in r), key="eps_d1_zero")

print("\n(b) exact ranks, from modular lower bounds against chain-complex ceilings")
def zmat(M):
    rows = []
    for row in M:
        for g in range(120):
            v = []
            for e in row:
                col = [0]*120
                for h, c in e.items(): col[MUL[g][h]] += c
                v += col
            rows.append(v)
    return rows
def rank_p(rows, p=P):
    if p > MAX_INT64_PRIME:
        raise ValueError(f"p={p} exceeds the int64-safe bound {MAX_INT64_PRIME}")
    A = np.array(rows, dtype=np.int64) % p; R, C = A.shape; rk = 0
    for c in range(C):
        nz = np.nonzero(A[rk:, c])[0]
        if not nz.size: continue
        r = rk + nz[0]; A[[rk, r]] = A[[r, rk]]
        A[rk] = (A[rk] * pow(int(A[rk, c]), p-2, p)) % p
        col = A[:, c].copy(); col[rk] = 0
        A = (A - np.outer(col, A[rk])) % p; rk += 1
        if rk == R: break
    return rk
Z1, Z2, Z3 = zmat(D["d1"]), zmat(D["d2"]), zmat(D["d3"])
lo1, lo2, lo3 = rank_p(Z1), rank_p(Z2), rank_p(Z3)
hi1 = 119                       # eps.d1 = 0 puts im(d1) in ker(eps), which has rank 119
hi2 = 240 - lo1                 # d2.d1 = 0
hi3 = 240 - lo2                 # d3.d2 = 0
print(f"      rank d1: lower {lo1} (mod {P}), upper {hi1} (eps.d1 = 0)")
print(f"      rank d2: lower {lo2} (mod {P}), upper {hi2} (d2.d1 = 0)")
print(f"      rank d3: lower {lo3} (mod {P}), upper {hi3} (d3.d2 = 0)")
chk("rank d1 = 119 exactly", lo1 == hi1 == 119, key="rank_d1_exact")
chk("rank d2 = 121 exactly", lo2 == hi2 == 121, key="rank_d2_exact")
chk("rank d3 = 119 exactly", lo3 == hi3 == 119, key="rank_d3_exact")
chk("therefore rank ker(d2) = 240 - 121 = 119", (240 - lo2) == 119, key="ker_d2_rank_119")

print("\n(c) saturation of im(d3): gcd of maximal minors")
cols = sat.pivot_cols(Z3, P)   # same prime everywhere so the certificate records ONE p; sat.py itself is pure-Python and cannot overflow
chk(f"deterministic pivot columns number 119", len(cols) == 119)
g, minors = 0, []
for drop in (0, 1, 2):
    sub = [[Z3[i][c] for c in cols] for i in range(120) if i != drop]
    d = sat.bareiss(sub)
    g = gcd(g, abs(d)); minors.append({"dropped_row": drop, "det": str(d), "running_gcd": g})
    print(f"      minor with row {drop} dropped: det = {d}, running gcd = {g}")
    if g == 1: break
chk("gcd of exhibited maximal minors is 1, so every elementary divisor is 1", g == 1,
    key="all_elementary_divisors_1")

concl = not fail
print(f"\n=> im(d3) and ker(d2) are both saturated rank-119 sublattices of C_2,")
print(f"   and (a) puts one inside the other.")
print(f"=> im(d3) = ker(d2) as integral lattices: {'CERTIFIED' if concl else 'NOT CERTIFIED'}")

cert = {
  "packet_sha256": hashlib.sha256(F.read_bytes()).hexdigest(),
  # Every flag below is a MEASURED outcome. An earlier version wrote these as Python
  # literals, so the emitted file asserted its own premises even on a run where they
  # had failed, which is exactly the shape the maintainer's D11 gate names. The file
  # must never claim more than the run established.
  "containment": {k: measured[k] for k in ("d3_d2_zero", "d2_d1_zero", "eps_d1_zero")},
  "ranks": {"prime_for_lower_bounds": P,
            "d1": {"lower": lo1, "upper": hi1, "exact_established": measured["rank_d1_exact"]},
            "d2": {"lower": lo2, "upper": hi2, "exact_established": measured["rank_d2_exact"]},
            "d3": {"lower": lo3, "upper": hi3, "exact_established": measured["rank_d3_exact"]},
            "ker_d2_rank_119": measured["ker_d2_rank_119"]},
  "saturation": {"column_selection": f"pivot columns of deterministic elimination at p = {P}",
                 "columns": cols, "row_selection": "rows 0, 1, 2 dropped in order",
                 "minors": minors, "gcd": g,
                 "all_elementary_divisors_1": measured["all_elementary_divisors_1"]},
  "checks_failed": fail,
  "conclusion": "im(d3) = ker(d2) as integral lattices" if concl else "NOT CERTIFIED",
}
if "--emit" in sys.argv:
    out = pathlib.Path("saturation_certificate.json")
    out.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    print(f"\nwritten to {out}")
else:
    # Deliberately does NOT write by default. The certificate is a manifest-listed file,
    # so writing on every run means an auditor who reruns before hashing sees a manifest
    # failure caused by their own audit. Pass --emit to regenerate it on purpose.
    print("\n(not written; pass --emit to regenerate saturation_certificate.json)")
sys.exit(0 if concl else 1)
