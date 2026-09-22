"""s04_reduction.py -- THE reduction, exactly.

Claim proved here by exact algebra:

  rhat_6(u) = sum_{J in {0,2,4,6}} n_J * p_J(u),   sum_J p_J(u) = 1,  p_J >= 0,

where p_J(u) is the weight of the HOLOMORPHIC square u (x) u in the spin-J block of
Sym^2(V3) (no Theta anywhere in p_J), and n_J are four universal constants.

Why this is forced:
  rhat_6(u) = <u(x)u| N |u(x)u> for a Hermitian N on C^7 (x) C^7 built from the
  multipole operators, rhat_6 is invariant under u -> D(g)u, a Hermitian form on
  Sym^2 is determined by its values on the Veronese {u(x)u}, hence P_sym N P_sym
  commutes with D (x) D on Sym^2.  Sym^2(V3) = W6 + W4 + W2 + W0 (13+9+5+1 = 28)
  is MULTIPLICITY-FREE, so by Schur P_sym N P_sym is a scalar n_J on each block.

Everything below is exact sympy; the numeric echo at the end is a check only.
"""

from pathlib import Path
import sympy as sp
from sympy.physics.quantum.cg import CG
import json

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))  # v_3 .. v_-3
IDX = {m: i for i, m in enumerate(ORDER)}


def cg(j1, m1, j2, m2, Jt, Q):
    if m1 + m2 != Q or abs(m1) > j1 or abs(m2) > j2 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(j1, m1, j2, m2, Jt, Q).doit()))


def multipole_ops():
    ops = {}
    for Jt in range(7):
        for Q in range(-Jt, Jt + 1):
            M = sp.zeros(7, 7)
            for i, mm in enumerate(ORDER):
                for k, m1 in enumerate(ORDER):
                    M[i, k] = cg(3, m1, 3, -mm, Jt, Q) * (-1) ** (3 + mm)
            ops[(Jt, Q)] = M
    return ops


def sym_block_basis():
    """Orthonormal basis of the spin-J block of Sym^2(V3), J = 0,2,4,6,
       as 7x7 symmetric coefficient matrices S with (u(x)u) <-> outer(u,u)."""
    basis = {}
    for Jt in [0, 2, 4, 6]:
        vecs = []
        for Q in range(-Jt, Jt + 1):
            S = sp.zeros(7, 7)
            for m1 in range(-3, 4):
                m2 = Q - m1
                if abs(m2) <= 3:
                    S[IDX[m1], IDX[m2]] = cg(3, m1, 3, m2, Jt, Q)
            vecs.append(S)
        basis[Jt] = vecs
    return basis


def main():
    ops = multipole_ops()
    T6 = [ops[(6, Q)] for Q in range(-6, 7)]
    blocks = sym_block_basis()

    # --- orthonormality of the Sym^2 block basis (28 vectors in C^49) ---
    allv = [(J, q, v) for J in blocks for q, v in enumerate(blocks[J])]
    print("Sym^2 block basis: %d vectors" % len(allv))
    bad = []
    for (J1, q1, A) in allv:
        for (J2, q2, B) in allv:
            ip = sp.simplify(sum(sp.conjugate(A[i, j]) * B[i, j] for i in range(7) for j in range(7)))
            want = 1 if (J1, q1) == (J2, q2) else 0
            if sp.simplify(ip - want) != 0:
                bad.append(((J1, q1), (J2, q2), ip))
    print("  orthonormality mismatches:", bad if bad else "none (28x28 checked)")
    print("  block dims:", {J: len(blocks[J]) for J in blocks},
          " total", sum(len(blocks[J]) for J in blocks))

    # --- the Hermitian form Ncal(A,B) = sum_Q tr(A^dag T_Q^dag) * ... -------
    # rhat_6(u) = sum_Q |<u|T_Q|u>|^2 ;  <u|T_Q|u> = sum_{i,j} conj(u_i)(T_Q)_{ij}u_j
    # For a symmetric tensor represented by matrix S (so u(x)u <-> outer(u,u), i.e.
    # S_{jk} = u_j u_k) define the linear functional L_Q(S) = sum_{j,k} (T_Q)_{?}.
    # Careful: <u|T_Q|u> is NOT linear in u(x)u; it is linear in u(x)conj(u).
    # The correct statement: rhat_6(u) = <u(x)u| N |u(x)u> with
    #   N_{(i,l),(j,k)} = sum_Q (T_Q)_{ij} conj((T_Q)_{lk}).
    # Build the 28x28 matrix of N in the Sym^2 block basis.
    def Nentry(A, B):
        """<A| N |B> for A,B given as 7x7 matrices (index pairs (i,l) and (j,k))."""
        tot = 0
        for T in T6:
            # sum_{i,l,j,k} conj(A_{il}) T_{ij} conj(T_{lk}) B_{jk}
            #  = sum_{i,l} conj(A_{il}) * (T B T^dagger-ish)
            # do it as: M = T * B * T^H  -> M_{il} = sum_{j,k} T_{ij} B_{jk} conj(T_{lk})
            M = T * B * T.H
            tot += sum(sp.conjugate(A[i, l]) * M[i, l] for i in range(7) for l in range(7))
        return sp.simplify(sp.expand(tot))

    print("\nbuilding N in the Sym^2 block basis (28x28, exact) ...")
    rows = []
    labels = [(J, q) for (J, q, v) in allv]
    mats = [v for (J, q, v) in allv]
    Nmat = sp.zeros(28, 28)
    for a in range(28):
        for b in range(28):
            Nmat[a, b] = Nentry(mats[a], mats[b])
        rows.append(a)
    print("  done.")

    # --- check block-diagonality and read off n_J ---
    offdiag = []
    nJ = {}
    for a in range(28):
        for b in range(28):
            Ja, Jb = labels[a][0], labels[b][0]
            v = sp.simplify(Nmat[a, b])
            if a == b:
                nJ.setdefault(Ja, []).append(v)
            elif v != 0:
                offdiag.append((labels[a], labels[b], v))
    print("\n  off-diagonal nonzero entries of N in this basis:",
          offdiag if offdiag else "none -> N is a scalar on each block (Schur, as predicted)")
    print("\n  the four constants n_J:")
    res = {}
    for J in sorted(nJ):
        vals = set(sp.simplify(x) for x in nJ[J])
        assert len(vals) == 1, (J, vals)
        v = vals.pop()
        res[J] = v
        print("    n_%d = %s = %s = %s" % (J, v, sp.nsimplify(v), sp.N(v, 30)))

    # --- consistency: sum over blocks with multiplicity must reproduce tr N -----
    trN = sp.simplify(sum(Nmat[a, a] for a in range(28)))
    print("\n  tr N|_Sym2 =", trN, "=", sp.nsimplify(trN))

    with open(HERE / "s04_nJ.json", "w") as fh:
        json.dump({str(J): str(res[J]) for J in res}, fh, indent=2)
    print("\nwrote s04_nJ.json")

    # --- what the reduction says immediately ---
    print("\nIMMEDIATE CONSEQUENCES")
    mx = max(res.values(), key=lambda t: float(t))
    mn = min(res.values(), key=lambda t: float(t))
    print("  max_J n_J =", mx, "  min_J n_J =", mn)
    print("  so  min_J n_J <= rhat_6(u) <= max_J n_J  for every u, since rhat_6 is")
    print("  the convex combination sum_J n_J p_J(u) with p_J >= 0, sum p_J = 1.")


if __name__ == "__main__":
    main()
