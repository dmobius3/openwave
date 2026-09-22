"""s05_recoupling.py -- the full linear map p -> rhat.

Same Schur argument as s04b, applied to EVERY multipole weight rhat_K, K = 0..6:
each rhat_K(u) is a rotation-invariant (2,2)-form in u, hence a Hermitian form on
Sym^2(V3) that is scalar on each of the four multiplicity-free blocks.  So

    rhat_K(u) = sum_{J in {0,2,4,6}} A[K,J] * p_J(u),     sum_J p_J(u) = 1.

This script computes A exactly, checks the two identities it must satisfy
(row K=0 constant 1/7, columns summing to 1), and sets up the linear programme
      maximise A[6,.] p   s.t.  p >= 0, sum p = 1, A[K,.] p >= 0 for all K,
whose value is a rigorous UPPER bound for max rhat_6 (every achievable p obeys
all those constraints, because rhat_K >= 0 is a sum of squares).
"""

from pathlib import Path
import sympy as sp
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_6j
import json

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}


def cgx(j1, m1, j2, m2, Jt, Q):
    if m1 + m2 != Q or abs(m1) > j1 or abs(m2) > j2 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(j1, m1, j2, m2, Jt, Q).doit()))


def TK_exact(K):
    T = []
    for Q in range(-K, K + 1):
        M = sp.zeros(7, 7)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                M[i, k] = cgx(3, m1, 3, -mm, K, Q) * (-1) ** (3 + mm)
        T.append(M)
    return T


def sym_block_hw(J):
    """highest-weight vector (Q = +J) of the spin-J block of Sym^2, as a 7x7 matrix"""
    S = sp.zeros(7, 7)
    for m1 in range(-3, 4):
        m2 = J - m1
        if abs(m2) <= 3:
            S[IDX[m1], IDX[m2]] = cgx(3, m1, 3, m2, J, J)
    return S


def main():
    Js = [0, 2, 4, 6]
    hw = {J: sym_block_hw(J) for J in Js}
    A = sp.zeros(7, 4)
    for K in range(7):
        TK = TK_exact(K)
        for jj, J in enumerate(Js):
            psi = hw[J]
            tot = 0
            for T in TK:
                M = T * psi * T.applyfunc(sp.conjugate)
                tot += sum(sp.conjugate(psi[i, k]) * M[i, k] for i in range(7) for k in range(7))
            A[K, jj] = sp.nsimplify(sp.simplify(sp.expand(tot)))
        print("row K=%d :" % K, [A[K, j] for j in range(4)])

    print("\nchecks on A:")
    print("  row K=0 all equal 1/7 :", all(sp.simplify(A[0, j] - sp.Rational(1, 7)) == 0 for j in range(4)))
    colsums = [sp.simplify(sum(A[K, j] for K in range(7))) for j in range(4)]
    print("  column sums (must all be 1):", colsums)

    # denominators as a table
    print("\n  A * 924 :")
    for K in range(7):
        print("    K=%d :" % K, [sp.nsimplify(A[K, j] * 924) for j in range(4)])

    # ---- independent cross-check against the Racah 6j crossing matrix ----
    print("\n independent 6j cross-check   A[K,J] =? (2J+1) {3 3 K ; 3 3 J}^2 * (2K+1)/(something)")
    for K in range(7):
        row = []
        for J in Js:
            s = wigner_6j(3, 3, K, 3, 3, J)
            row.append(sp.nsimplify(sp.simplify((2 * J + 1) * (2 * K + 1) * s ** 2)))
        print("    K=%d  (2J+1)(2K+1){6j}^2 =" % K, row, "   A row =", [A[K, j] for j in range(4)])
        print("           ratio A/6jexpr =",
              [sp.simplify(A[K, j] / row[j]) if row[j] != 0 else None for j in range(4)])

    with open(HERE / "s05_A.json", "w") as fh:
        json.dump({"A": [[str(A[K, j]) for j in range(4)] for K in range(7)], "Js": Js}, fh, indent=2)
    print("\nwrote s05_A.json")

    # ---------------- the LP ----------------
    print("\nLP relaxation: maximise A[6,.]p  s.t. p>=0, sum p =1, A[K,.]p >= 0 all K")
    p = sp.symbols("p0 p2 p4 p6", nonnegative=True)
    obj = sum(A[6, j] * p[j] for j in range(4))
    cons = [sum(p) - 1]
    print("  objective:", obj)
    for K in range(7):
        print("  constraint K=%d : " % K, sum(A[K, j] * p[j] for j in range(4)), ">= 0")

    # enumerate vertices of the polytope {p>=0, sum p=1, A[K,.]p>=0} in 3 free dims
    import itertools
    import numpy as np
    Af = np.array([[float(A[K, j]) for j in range(4)] for K in range(7)])
    # constraint rows: p_j >= 0  (4), A[K,.]p >= 0 (7), and equality sum p = 1
    rows = []
    for j in range(4):
        e = np.zeros(4)
        e[j] = 1.0
        rows.append(e)
    for K in range(7):
        rows.append(Af[K])
    rows = np.array(rows)
    verts = []
    for comb in itertools.combinations(range(len(rows)), 3):
        Mx = np.vstack([rows[list(comb)], np.ones(4)])
        rhs = np.array([0.0, 0.0, 0.0, 1.0])
        try:
            sol = np.linalg.solve(Mx, rhs)
        except np.linalg.LinAlgError:
            continue
        if np.all(sol >= -1e-10) and np.all(rows @ sol >= -1e-10):
            verts.append((float(Af[6] @ sol), sol, comb))
    verts.sort(key=lambda t: -t[0])
    print("\n  polytope vertices, best first (objective, p, active set):")
    seen = set()
    for v, sol, comb in verts:
        key = tuple(np.round(sol, 9))
        if key in seen:
            continue
        seen.add(key)
        print("    %.15f   p = %s   active %s" % (v, np.round(sol, 12), comb))
    print("\n  LP optimum = %.15f ;  463/924 = %.15f ; 6/7 = %.15f"
          % (verts[0][0] if verts else float("nan"), 463 / 924, 6 / 7))


if __name__ == "__main__":
    main()
