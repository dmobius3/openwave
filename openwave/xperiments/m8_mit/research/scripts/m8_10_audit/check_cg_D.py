"""Cross-checks of the CG routine and the D^j construction (no tabulated CG used).
1. Racah formula vs highest-weight/lowering construction, exact, all triples j1,j2 <= 6 and the
   triples used later (listed below).
2. CS phase: J+ matrix elements positive in the D construction (derivative of D at identity).
3. D^j is a homomorphism, unitary; conj(D_mk) = (-1)^(m-k) D_{-m,-k}.
4. Product formula D^j1_{m1k1} D^j2_{m2k2} = sum_J CG CG D^J_{MK} at random g (ties CG to D).
"""
import itertools
import json

import numpy as np
import sympy as sp

import lib

C = lib.Checker("cg_D")
rng = np.random.default_rng(12345)

# 1. Racah vs lowering
triples = set()
for j1 in range(0, 7):
    for j2 in range(0, 7):
        for J in range(abs(j1 - j2), j1 + j2 + 1):
            triples.add((j1, j2, J))
# triples used by the algebra (3,3,K), (K,3,J), (3,J,K2), (J,3,K2), (K2,3,3), (K2,J,3) with J<=9, K2<=12
for J in range(0, 10):
    for K2 in range(abs(J - 3), J + 4):
        triples |= {(3, J, K2), (J, 3, K2), (K2, 3, 3), (K2, J, 3), (6, 3, J)}
triples = sorted(t for t in triples if abs(t[0] - t[1]) <= t[2] <= t[0] + t[1])
maxdiff = 0.0
n_exact_mismatch = 0
for (j1, j2, J) in triples:
    tab = lib.cg_lowering_table(j1, j2, J)
    for m1 in range(-j1, j1 + 1):
        for m2 in range(-j2, j2 + 1):
            M = m1 + m2
            if abs(M) > J:
                continue
            a = lib.cg(j1, m1, j2, m2, J, M)
            b = tab.get((m1, m2, M), 0)
            d = sp.N(a - b, 30)
            maxdiff = max(maxdiff, abs(float(d)))
            if abs(d) > 1e-25:
                n_exact_mismatch += 1
print("triples compared:", len(triples), " max |Racah - lowering| =", maxdiff)
C.check("Racah CG == lowering CG on all compared triples (30-digit)",
        lambda z: z[0] == 0 and z[1] < 1e-25, (n_exact_mismatch, maxdiff),
        (n_exact_mismatch + 1, maxdiff), "count one extra mismatch")

# symbolic exact equality on a sub-sample (all (3,3,K))
bad = 0
for K in range(7):
    tab = lib.cg_lowering_table(3, 3, K)
    for m1 in range(-3, 4):
        for m2 in range(-3, 4):
            if abs(m1 + m2) <= K:
                if sp.simplify(lib.cg(3, m1, 3, m2, K, m1 + m2) - tab[(m1, m2, m1 + m2)]) != 0:
                    bad += 1
C.check("Racah == lowering symbolically for all (3,3,K)", lambda z: z == 0, bad, bad + 1,
        "pretend one mismatch")

# orthogonality of Racah table
for (j1, j2) in [(3, 3), (6, 3), (3, 9)]:
    Js = range(abs(j1 - j2), j1 + j2 + 1)
    rows = [(J, M) for J in Js for M in range(-J, J + 1)]
    cols = [(m1, m2) for m1 in range(-j1, j1 + 1) for m2 in range(-j2, j2 + 1)]
    Mx = np.array([[float(lib.cg(j1, m1, j2, m2, J, M)) for (m1, m2) in cols] for (J, M) in rows])
    err = np.abs(Mx @ Mx.T - np.eye(len(rows))).max()
    C.check(f"CG matrix ({j1}x{j2}) orthogonal, err={err:.1e}", lambda e: e < 1e-12, err, err + 1e-3,
            "add 1e-3 to the error")

# 2. CS phase of the D construction: J+ from derivative of D at identity
for n in (1, 2, 6, 12, 18):
    eps = 1e-6
    # U = exp(eps E) with E = [[0,1],[0,0]] -> U = [[1, eps],[0,1]] (not unitary but D is polynomial)
    Dp = lib.D_num(n, 1.0, eps, 0.0, 1.0)
    Jp_num = (Dp - np.eye(n + 1)) / eps
    Jz, Jp, Jm = lib.spin_mats_sym(n)
    Jp_ex = np.array(Jp.evalf(), dtype=complex)
    err = np.abs(Jp_num - Jp_ex).max()
    C.check(f"D^({n}/2) raising generator == CS J+ (err {err:.1e})", lambda e: e < 1e-4, err, err + 1,
            "add 1 to error")

# 3. homomorphism, unitarity, conjugation identity
def rand_su2(k):
    q = rng.normal(size=(k, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    return lib.su2_from_quat_float(*q.T)

for n in (1, 6, 13, 18):
    a1 = rand_su2(3)
    a2 = rand_su2(3)
    U1 = np.array([[a1[0], a1[1]], [a1[2], a1[3]]]).transpose(2, 0, 1)
    U2 = np.array([[a2[0], a2[1]], [a2[2], a2[3]]]).transpose(2, 0, 1)
    U12 = U1 @ U2
    D1 = lib.D_num(n, *a1)
    D2 = lib.D_num(n, *a2)
    D12 = lib.D_num(n, U12[:, 0, 0], U12[:, 0, 1], U12[:, 1, 0], U12[:, 1, 1])
    e_hom = np.abs(D1 @ D2 - D12).max()
    e_uni = np.abs(D1 @ D1.conj().transpose(0, 2, 1) - np.eye(n + 1)).max()
    sgn = np.array([[(-1) ** (i - k) for k in range(n + 1)] for i in range(n + 1)])
    e_conj = np.abs(D1.conj() - sgn * D1[:, ::-1, ::-1]).max()
    C.check(f"D^({n}/2) homomorphism err {e_hom:.1e}", lambda e: e < 1e-10, e_hom, e_hom + 1, "add 1")
    C.check(f"D^({n}/2) unitary err {e_uni:.1e}", lambda e: e < 1e-10, e_uni, e_uni + 1, "add 1")
    C.check(f"conj D = (-1)^(m-k) D_(-m,-k) for n={n}, err {e_conj:.1e}", lambda e: e < 1e-10,
            e_conj, np.abs(D1.conj() - D1[:, ::-1, ::-1]).max() if n % 2 == 0 and n > 0 else e_conj + 1,
            "drop the (-1)^(m-k) sign")

# 4. product formula
a = rand_su2(2)
for (j1, j2) in [(3, 3), (6, 3), (3, 9)]:
    D1 = lib.D_num(2 * j1, *a)
    D2 = lib.D_num(2 * j2, *a)
    DJ = {J: lib.D_num(2 * J, *a) for J in range(abs(j1 - j2), j1 + j2 + 1)}
    err = 0.0
    err_mut = 0.0
    for m1, k1, m2, k2 in itertools.product(range(-j1, j1 + 1), range(-j1, j1 + 1),
                                            range(-j2, j2 + 1), range(-j2, j2 + 1)):
        lhs = D1[:, m1 + j1, k1 + j1] * D2[:, m2 + j2, k2 + j2]
        rhs = 0
        rhs_m = 0
        for J in DJ:
            M, K = m1 + m2, k1 + k2
            if abs(M) > J or abs(K) > J:
                continue
            c1 = float(lib.cg(j1, m1, j2, m2, J, M))
            c2 = float(lib.cg(j1, k1, j2, k2, J, K))
            rhs = rhs + c1 * c2 * DJ[J][:, M + J, K + J]
            rhs_m = rhs_m + abs(c1) * c2 * DJ[J][:, M + J, K + J]
        err = max(err, np.abs(lhs - rhs).max())
        err_mut = max(err_mut, np.abs(lhs - rhs_m).max())
    C.check(f"product formula D^{j1} D^{j2} = sum CG CG D^J, err {err:.1e}", lambda e: e < 1e-10,
            err, err_mut, "replace one CG by its absolute value (phase convention broken)")

json.dump({"records": C.records, "all_ok": C.all_ok()},
          open(lib.ROOM / "check_cg_D.json", "w"), indent=1)
print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
