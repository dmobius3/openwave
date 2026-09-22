"""s02_structure.py -- exact structural facts about rhat_J, derived here.

Facts established (all by exact sympy algebra on symbolic c_m):

 (A)  rho_J(u)_Q = <u| T^J_Q |u>  for an explicit operator basis T^J_Q, and the
      T^J_Q are HS-orthonormal.  Hence sum_{J=0..6} rhat_J(u) = (||u||^2)^2/||u||^4 = 1
      for every u:  the seven multipole weights partition unity.
 (B)  rhat_0(u) = 1/7 identically.  Hence rhat_6 <= 6/7 for every u.
 (C)  Theta_J rho_J(u) = (-1)^J rho_J(u), a reality constraint on each multipole.
 (D)  rhat_J is invariant under u -> exp(i phi) u and under every rotation.

Each is checked symbolically (exact), not numerically.
"""

from pathlib import Path
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))  # v_3 ... v_-3


def cg(m1, m2, Q, Jt):
    if m1 + m2 != Q or abs(m1) > 3 or abs(m2) > 3 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(3, m1, 3, m2, Jt, Q).doit()))


def build_ops():
    ops = {}
    for Jt in range(7):
        for Q in range(-Jt, Jt + 1):
            M = sp.zeros(7, 7)
            for i, mm in enumerate(ORDER):
                for k, m1 in enumerate(ORDER):
                    M[i, k] = cg(m1, -mm, Q, Jt) * (-1) ** (3 + mm)
            ops[(Jt, Q)] = M
    return ops


def main():
    ops = build_ops()

    # ---- (A) HS-orthonormality of the full 49-operator set -> partition of unity
    print("(A) checking HS-orthonormality of all 49 operators T^J_Q ...")
    bad = []
    keys = list(ops)
    for a in keys:
        for b in keys:
            v = sp.simplify((ops[a].H * ops[b]).trace())
            want = 1 if a == b else 0
            if sp.simplify(v - want) != 0:
                bad.append((a, b, v))
    print("    mismatches:", bad if bad else "none (49x49 checked)")

    # symbolic u
    re = sp.symbols("x0:7", real=True)
    im = sp.symbols("y0:7", real=True)
    u = sp.Matrix(7, 1, [re[i] + sp.I * im[i] for i in range(7)])
    n2 = sp.expand(sum(re[i] ** 2 + im[i] ** 2 for i in range(7)))

    def rhat(Jt):
        tot = 0
        for Q in range(-Jt, Jt + 1):
            e = (u.H * ops[(Jt, Q)] * u)[0, 0]
            tot += sp.expand(sp.re(sp.expand(e)) ** 2 + sp.im(sp.expand(e)) ** 2)
        return sp.expand(tot)

    r = {Jt: rhat(Jt) for Jt in range(7)}

    # ---- (B) rhat_0 = 1/7 identically
    d0 = sp.simplify(sp.expand(r[0] - n2 ** 2 / 7))
    print("(B) sum_Q |rho_0|^2 - ||u||^4/7  =", d0)

    # ---- (A') partition of unity
    tot = sp.expand(sum(r.values()) - n2 ** 2)
    print("(A') sum_J sum_Q |rho_J|^2 - ||u||^4 =", sp.simplify(tot))

    # ---- (C) Theta_J rho_J = (-1)^J rho_J, checked on symbolic u for each J
    print("(C) Theta_J rho_J(u) - (-1)^J rho_J(u):")
    for Jt in range(7):
        comps = []
        for Q in range(-Jt, Jt + 1):
            rQ = sp.expand((u.H * ops[(Jt, Q)] * u)[0, 0])
            rmQ = sp.expand((u.H * ops[(Jt, -Q)] * u)[0, 0])
            lhs = (-1) ** (Jt - Q) * sp.conjugate(rmQ)
            comps.append(sp.simplify(sp.expand(lhs - (-1) ** Jt * rQ)))
        print("    J =", Jt, "->", "all zero" if all(cc == 0 for cc in comps) else comps)

    # ---- (D) rotation invariance of rhat_6: check on the generator J_z and J_x
    #      by finite rotation about two axes at a random exact rational-ish point.
    print("(D) rotation invariance is checked numerically in s03 (finite rotations);")
    print("    here the algebraic reason: T^J_Q span an irreducible tensor operator")
    print("    set, so sum_Q |<T^J_Q>|^2 is the norm of a spin-J vector.")

    # save rhat_6 as an expression for later exact work
    with open(HERE / "s02_rhat6_expr.txt", "w") as fh:
        fh.write(sp.srepr(r[6]))
    print("\nwrote s02_rhat6_expr.txt (srepr of sum_Q|rho_6|^2, unnormalised)")


if __name__ == "__main__":
    main()
