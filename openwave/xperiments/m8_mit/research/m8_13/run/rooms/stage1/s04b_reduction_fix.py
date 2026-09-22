"""s04b_reduction_fix.py -- corrected build of N, then the four constants n_J.

CORRECTION to s04_reduction.py (recorded, not smoothed over):
  In s04 I paired the indices as N_{(i,l),(j,k)} = sum_Q (T_Q)_{ij} conj((T_Q)_{lk})
  and evaluated it as A^dag . (T B T^dagger).  That is one transpose too many.
  The correct expansion of |<u|T_Q|u>|^2 is

    |<u|T|u>|^2 = sum_{i,j,k,l} conj(u_i u_k) (T)_{ij} conj((T)_{lk}) (u_j u_l),

  i.e. the BRA pair is (i,k) and the KET pair is (j,l), and the contraction is
  A^dag . (T B Tbar) with Tbar the ENTRYWISE conjugate (no transpose).
  With the wrong pairing the matrix was not block-scalar (n_2 came out
  double-valued, {125/924, -25/462}), which is what exposed the error.

Stage 1: numeric build + check of block-scalar structure.
Stage 2: exact sympy evaluation of the four constants on highest-weight vectors.
"""

from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG
import json

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}


def cgx(j1, m1, j2, m2, Jt, Q):
    if m1 + m2 != Q or abs(m1) > j1 or abs(m2) > j2 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(j1, m1, j2, m2, Jt, Q).doit()))


def T6_exact():
    T = []
    for Q in range(-6, 7):
        M = sp.zeros(7, 7)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                M[i, k] = cgx(3, m1, 3, -mm, 6, Q) * (-1) ** (3 + mm)
        T.append(M)
    return T


def sym_blocks_exact():
    basis = {}
    for Jt in [0, 2, 4, 6]:
        vecs = []
        for Q in range(-Jt, Jt + 1):
            S = sp.zeros(7, 7)
            for m1 in range(-3, 4):
                m2 = Q - m1
                if abs(m2) <= 3:
                    S[IDX[m1], IDX[m2]] = cgx(3, m1, 3, m2, Jt, Q)
            vecs.append(S)
        basis[Jt] = vecs
    return basis


def main():
    T6 = T6_exact()
    blocks = sym_blocks_exact()
    T6n = [np.array(sp.matrix2numpy(t, dtype=complex)) for t in T6]

    labels, mats = [], []
    for J in [0, 2, 4, 6]:
        for q, v in enumerate(blocks[J]):
            labels.append((J, q))
            mats.append(np.array(sp.matrix2numpy(v, dtype=complex)))

    def Nentry_num(A, B):
        tot = 0j
        for T in T6n:
            M = T @ B @ np.conjugate(T)
            tot += np.sum(np.conjugate(A) * M)
        return tot

    Nn = np.zeros((28, 28), dtype=complex)
    for a in range(28):
        for b in range(28):
            Nn[a, b] = Nentry_num(mats[a], mats[b])

    print("Hermiticity of N|_Sym2 : max|N - N^H| = %.3e" % np.max(np.abs(Nn - Nn.conj().T)))
    # block-scalar test
    worst_off = 0.0
    diag = {}
    for a in range(28):
        for b in range(28):
            if labels[a][0] == labels[b][0] and a == b:
                diag.setdefault(labels[a][0], []).append(Nn[a, a].real)
            else:
                worst_off = max(worst_off, abs(Nn[a, b]))
    print("max |off-block-diagonal or off-diagonal entry| = %.3e" % worst_off)
    for J in sorted(diag):
        arr = np.array(diag[J])
        print("  block J=%d  dim %2d  diagonal entries: min %.15f  max %.15f  spread %.2e"
              % (J, len(arr), arr.min(), arr.max(), arr.max() - arr.min()))

    # ---------- exact evaluation of n_J ----------
    print("\nexact n_J via <psi_J| N |psi_J> on the HIGHEST-WEIGHT vector of each block:")
    res = {}
    for J in [0, 2, 4, 6]:
        psi = blocks[J][-1]  # Q = +J, highest weight, sparse
        tot = 0
        for T in T6:
            M = T * psi * T.applyfunc(sp.conjugate)
            tot += sum(sp.conjugate(psi[i, k]) * M[i, k] for i in range(7) for k in range(7))
        val = sp.nsimplify(sp.simplify(sp.expand(tot)))
        res[J] = val
        print("   n_%d = %s = %s" % (J, val, sp.N(val, 30)))

    # cross-check: a second, independent exact evaluation on the Q=0 vector
    print("\n cross-check on the Q=0 vector of each block (must agree):")
    for J in [0, 2, 4, 6]:
        psi = blocks[J][J]  # Q = 0
        tot = 0
        for T in T6:
            M = T * psi * T.applyfunc(sp.conjugate)
            tot += sum(sp.conjugate(psi[i, k]) * M[i, k] for i in range(7) for k in range(7))
        val = sp.nsimplify(sp.simplify(sp.expand(tot)))
        print("   n_%d(Q=0) = %s   agrees: %s" % (J, val, sp.simplify(val - res[J]) == 0))

    with open(HERE / "s04b_nJ.json", "w") as fh:
        json.dump({str(J): str(res[J]) for J in res}, fh, indent=2)
    print("\nwrote s04b_nJ.json")

    # ---------- verify the reduction identity on random states ----------
    print("\nverifying  rhat_6(u) = sum_J n_J p_J(u)  on random complex u:")
    rng = np.random.default_rng(7)
    nJf = {J: float(res[J]) for J in res}
    blocks_n = {J: [np.array(sp.matrix2numpy(v, dtype=complex)) for v in blocks[J]] for J in blocks}
    for trial in range(6):
        u = rng.normal(size=7) + 1j * rng.normal(size=7)
        u /= np.linalg.norm(u)
        S = np.outer(u, u)
        r6 = sum(abs(np.einsum("i,ij,j->", u.conj(), T, u)) ** 2 for T in T6n).real
        p = {J: sum(abs(np.sum(np.conjugate(v) * S)) ** 2 for v in blocks_n[J]) for J in blocks_n}
        print("   rhat6 = %.15f   sum n_J p_J = %.15f   sum p_J = %.15f"
              % (r6, sum(nJf[J] * p[J] for J in p), sum(p.values())))


if __name__ == "__main__":
    main()
