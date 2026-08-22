"""Exact saturation test for a candidate d3, with both directions certified.

im(d3) = ker(d2) as INTEGRAL lattices, equivalently H_2(C_*) = 0 integrally,
equivalently every elementary divisor of the Z-matrix of d3 equals 1.

  NOT saturated : rank_p(M3) < 119 for some prime p.   Definitive.
  saturated     : gcd of sampled maximal (119x119) minors is 1.  Definitive,
                  because the gcd of ALL maximal minors is the product of the
                  elementary divisors, and it divides the gcd of any subset.

Neither direction is a finite-prime heuristic; the accept side is exact.
"""
import random
from math import gcd
N = 120

def zmat(M3, MUL):
    rows = []
    for g in range(N):
        v = []
        for e in M3:
            col = [0]*N
            for h, c in e.items(): col[MUL[g][h]] += c
            v += col
        rows.append(v)
    return rows

def rank_p(rows, p):
    A = [[x % p for x in r] for r in rows]; R, C = len(A), len(A[0]); rk = 0
    for c in range(C):
        piv = next((r for r in range(rk, R) if A[r][c]), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p-2, p)
        A[rk] = [(x*inv) % p for x in A[rk]]
        for r in range(R):
            if r != rk and A[r][c]:
                f = A[r][c]; A[r] = [(a-f*b) % p for a, b in zip(A[r], A[rk])]
        rk += 1
        if rk == R: break
    return rk

def bareiss(M):
    """Exact integer determinant, fraction-free."""
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            sw = next((i for i in range(k+1, n) if M[i][k]), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def pivot_cols(rows, p=10**9+7):
    """119 columns that stay independent mod a large prime."""
    A = [[x % p for x in r] for r in rows]; R, C = len(A), len(A[0])
    rk, cols = 0, []
    for c in range(C):
        piv = next((r for r in range(rk, R) if A[r][c]), None)
        if piv is None: continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p-2, p)
        A[rk] = [(x*inv) % p for x in A[rk]]
        for r in range(R):
            if r != rk and A[r][c]:
                f = A[r][c]; A[r] = [(a-f*b) % p for a, b in zip(A[r], A[rk])]
        cols.append(c); rk += 1
        if rk == R: break
    return cols

def saturated(rows, rng, screen=(2,3,5,7,11,13,17,19,23), minors=3, verbose=False):
    """Returns (verdict, detail). verdict True/False are both certified."""
    for p in screen:
        if rank_p(rows, p) < 119:
            return False, f"rank drops mod {p}: {rank_p(rows,p)} < 119"
    cols = pivot_cols(rows)
    if len(cols) != 119: return False, f"rational rank {len(cols)} != 119"
    g = 0
    for k in range(minors):
        drop = rng.randrange(N)                        # delete one row
        sub = [[rows[i][c] for c in cols] for i in range(N) if i != drop]
        d = bareiss(sub)
        g = gcd(g, abs(d))
        if verbose: print(f"      minor {k} (row {drop} dropped): |det| digits {len(str(abs(d)))}, gcd so far {g}")
        if g == 1:
            return True, f"gcd of {k+1} maximal minors is 1: every elementary divisor is 1"
    return None, f"inconclusive from {minors} minors, gcd = {g}"
