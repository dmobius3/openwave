"""Item 2 supplement: the seven maps M_K(v) = [[v (x) Theta v]_K (x) v]_3 span only a 4-dimensional space.
Basis used: B_K'(v) = [[v (x) v]_K' (x) Theta v]_3, K' in {0, 2, 4, 6} ([v (x) v]_K' = 0 for odd K' by CG exchange symmetry).
Solve M_K = sum_K' W_{K K'} B_K' by least squares on random v at two precisions, identify W_{KK'}^2 as rationals,
verify the identified exact W reproduces M_K (residual), and print the relations among the M_K.
Usage: ./py relations.py DPS"""
import sys
import json
import random
from fractions import Fraction as Fr
import mpmath
from common import HERE, cg_exact, identify_rational

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp = mpmath.mp
mp.dps = DPS
ZERO = mp.mpc(0)


def cg(j1, m1, j2, m2, J, M):
    S, A = cg_exact(j1, m1, j2, m2, J, M)
    return mp.mpf(0) if S == 0 else mp.mpf(S.numerator) / S.denominator * mp.sqrt(mp.mpf(A.numerator) / A.denominator)


def theta(u):
    return [(1 if (i - 3) % 2 == 0 else -1) * mp.conj(u[6 - i]) for i in range(7)]


def couple(x, j1, y, j2, J):
    return [sum((cg(j1, m1, j2, M - m1, J, M) * x[m1 + j1] * y[M - m1 + j2]
                 for m1 in range(-j1, j1 + 1) if abs(M - m1) <= j2), ZERO) for M in range(-J, J + 1)]


def Mmap(v, K):
    return couple(couple(v, 3, theta(v), 3, K), K, v, 3, 3)


def Bmap(v, K):
    return couple(couple(v, 3, v, 3, K), K, theta(v), 3, 3)


rng = random.Random(99)
vs = [[mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(7)] for _ in range(5)]
Kp = [0, 2, 4, 6]
# odd K' vanish
oddmax = max(abs(x) for v in vs for K in (1, 3, 5) for x in couple(v, 3, v, 3, K))
print("max |[v x v]_K'| for odd K':", mpmath.nstr(oddmax, 3))
Bm = mp.matrix(len(vs) * 7 * 2, 4)
rowsM = {K: [] for K in range(7)}
r = 0
Bvals = {K: [Bmap(v, K) for v in vs] for K in Kp}
for iv, v in enumerate(vs):
    for i in range(7):
        for part in (0, 1):
            for c, K in enumerate(Kp):
                x = Bvals[K][iv][i]
                Bm[r, c] = mp.re(x) if part == 0 else mp.im(x)
            r += 1
W = {}
tol = mp.mpf(10) ** (-(DPS - 15))
out = {"dps": DPS, "W": {}}
for K in range(7):
    rhs = mp.matrix(len(vs) * 14, 1)
    r = 0
    for v in vs:
        mk = Mmap(v, K)
        for i in range(7):
            rhs[r] = mp.re(mk[i]); r += 1
            rhs[r] = mp.im(mk[i]); r += 1
    sol, res = mp.qr_solve(Bm, rhs)
    row = []
    for c in range(4):
        x = sol[c]
        fr, err = identify_rational(x * x, mp, 10 ** 12, tol)
        assert fr is not None, (K, c, x)
        sign = 0 if fr == 0 else (1 if x > 0 else -1)
        row.append((sign, fr))
    W[K] = row
    out["W"][str(K)] = [("-" if s < 0 else "+" if s > 0 else "0") + f"sqrt({f})" for s, f in row]
    print(f"M_{K} = " + " + ".join(f"({'-' if s < 0 else ''}sqrt({f})) B_{Kp[c]}" for c, (s, f) in enumerate(row)) +
          f"   [lstsq residual {mpmath.nstr(res, 3)}]")
# verify the exact W on fresh vectors (can fail)
dev = 0
for _ in range(3):
    v = [mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(7)]
    for K in range(7):
        mk = Mmap(v, K)
        rec = [ZERO] * 7
        for c, (s, f) in enumerate(W[K]):
            coef = s * mp.sqrt(mp.mpf(f.numerator) / f.denominator)
            b = Bmap(v, Kp[c])
            rec = [x + coef * y for x, y in zip(rec, b)]
        dev = max(dev, max(abs(x - y) for x, y in zip(mk, rec)))
print("max |M_K - sum W B| on fresh v with exact W:", mpmath.nstr(dev, 3))
out["verify_dev"] = mpmath.nstr(dev, 5)
assert dev < tol
(HERE / f"out_relations_{DPS}.json").write_text(json.dumps(out, indent=1))
