"""s01_setup.py -- build the objects of worklist section 1 exactly, and check item 0.

Everything here is exact (sympy Rational / sqrt).  Nothing numerical.

Conventions follow worklist.md section 1 literally:
  V3 = C^7, basis v_3 ... v_-3, u = sum_m c_m v_m
  J_z v_m = m v_m,  J_+- v_m = sqrt(12 - m(m+-1)) v_{m+-1}
  (Theta u)_m = (-1)^(3-m) conj(c_{-m})
  rho6(u)_Q = sum_{m1} <3 m1; 3 (Q-m1) | 6 Q> c_{m1} (Theta u)_{Q-m1}
  rhat6(u)  = sum_Q |rho6_Q|^2 / ||u||^4
"""

from pathlib import Path
import json
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
J = 3  # spin of V3


def cg(m1, m2, Q, Jtot=6):
    """<3 m1; 3 m2 | Jtot Q>, Condon-Shortley, via sympy."""
    if m1 + m2 != Q:
        return sp.Integer(0)
    if abs(m1) > 3 or abs(m2) > 3 or abs(Q) > Jtot:
        return sp.Integer(0)
    return sp.simplify(CG(3, m1, 3, m2, Jtot, Q).doit())


def cg_stretched_closed_form(m1, m2):
    """Independent route: the closed form for the stretched coupling
       <j1 m1; j2 m2 | (j1+j2) M>
         = sqrt( (2j1)!(2j2)!(J+M)!(J-M)!
                 / [ (2J)! (j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)! ] )
       with j1 = j2 = 3, J = 6, M = m1+m2.  Derived, not looked up: it is the
       normalisation of the highest-weight component of the symmetric product,
       checked below against sympy for every (m1,m2)."""
    M = m1 + m2
    if abs(M) > 6:
        return sp.Integer(0)
    f = sp.factorial
    num = f(6) * f(6) * f(6 + M) * f(6 - M)
    den = f(12) * f(3 + m1) * f(3 - m1) * f(3 + m2) * f(3 - m2)
    return sp.sqrt(sp.Rational(num, den))


def main():
    out = {}

    # ---------- item 0, part a:  <3 3; 3 -3 | 6 0> ----------
    a_sympy = cg(3, -3, 0)
    a_closed = cg_stretched_closed_form(3, -3)
    print("item 0a: <3 3; 3 -3 | 6 0>")
    print("   sympy CG            :", a_sympy, "=", sp.nsimplify(a_sympy))
    print("   closed form          :", a_closed)
    print("   difference simplifies:", sp.simplify(a_sympy - a_closed))
    print("   squared              :", sp.simplify(a_sympy**2))
    print("   50-digit value       :", sp.N(a_sympy, 50))
    out["cg_3_3_3_m3_6_0"] = str(a_sympy)
    out["cg_3_3_3_m3_6_0_sq"] = str(sp.simplify(a_sympy**2))

    # cross-check the two routes on ALL (m1,m2)
    bad = []
    for m1 in range(-3, 4):
        for m2 in range(-3, 4):
            d = sp.simplify(cg(m1, m2, m1 + m2) - cg_stretched_closed_form(m1, m2))
            if d != 0:
                bad.append((m1, m2, d))
    print("   stretched closed form vs sympy, mismatches over all 49 pairs:", bad)
    out["cg_routes_agree"] = (len(bad) == 0)

    # ---------- item 0, part b:  Theta(Theta u) ----------
    c = sp.symbols("c_m3 c_m2 c_m1 c_0 c_p1 c_p2 c_p3")  # m = -3..3
    idx = {m: i for i, m in enumerate(range(-3, 4))}
    cs = {m: c[idx[m]] for m in range(-3, 4)}

    def Theta(vec):
        return {m: (-1) ** (3 - m) * sp.conjugate(vec[-m]) for m in range(-3, 4)}

    tt = Theta(Theta(cs))
    diffs = [sp.simplify(tt[m] - cs[m]) for m in range(-3, 4)]
    print("item 0b: Theta(Theta u) - u componentwise:", diffs)
    out["theta_squared_is_identity"] = all(d == 0 for d in diffs)
    # symbolic reason, independent of sympy's conjugate():
    print("   reason: (-1)^(3-m) * (-1)^(3+m) = (-1)^6 = +1 for every integer m")

    # ---------- the multipole operators T^J_Q with rho_J(u)_Q = <u|T^J_Q|u> ----
    # rho_J(u)_Q = sum_{m1,m2} <3m1;3m2|JQ> c_{m1} (Theta u)_{m2}
    #            = sum_{m1,m}  <3m1;3,-m|JQ> (-1)^(3+m) c_{m1} conj(c_m)
    # so (T^J_Q)_{m,m1} = <3 m1; 3 -m | J Q> (-1)^(3+m)   (row m = bra, col m1 = ket)
    ops = {}
    order = list(range(3, -4, -1))  # v_3 ... v_-3 as in the worklist
    for Jt in range(0, 7):
        for Q in range(-Jt, Jt + 1):
            M = sp.zeros(7, 7)
            for i, mm in enumerate(order):        # bra index
                for k, m1 in enumerate(order):    # ket index
                    M[i, k] = cg(m1, -mm, Q, Jt) * (-1) ** (3 + mm)
            ops[(Jt, Q)] = M
    sp.pickling = None

    # sanity: tr(T^J_Q^dagger T^J'_Q') = delta  (the CG matrix is orthogonal)
    checks = []
    for (J1, Q1) in [(0, 0), (1, 0), (6, 0), (6, 3), (4, -2)]:
        for (J2, Q2) in [(0, 0), (1, 0), (6, 0), (6, 3), (4, -2)]:
            v = sp.simplify((ops[(J1, Q1)].H * ops[(J2, Q2)]).trace())
            checks.append(((J1, Q1), (J2, Q2), v))
    print("item 0 extra: orthonormality tr(T+ T) spot checks")
    for ch in checks:
        print("   ", ch)

    with open(HERE / "s01_out.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote", HERE / "s01_out.json")


if __name__ == "__main__":
    main()
