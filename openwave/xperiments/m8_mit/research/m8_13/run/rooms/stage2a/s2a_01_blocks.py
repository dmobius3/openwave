"""s2a_01_blocks.py -- rebuild the objects of step 3 from scratch, exactly.

Nothing here is taken on trust from step3_text.md.  Spin-3 operators are built
from the stage-1 conventions (J_z v_m = m v_m, J_+- v_m = sqrt(12 - m(m+-1)) v_{m+-1}),
and every claim of the text's setup is re-derived:

  (1)  f^2 = 12 I
  (2)  e1 fx^2 + e2 fy^2 + e3 fz^2 = a(fz^2 - 4) + b(f+^2 + f-^2)
       under e3 = 2a/3, e1 - e2 = 4b, e1+e2+e3 = 0
  (3)  |e| = 1  <=>  (2/3)a^2 + 8b^2 = 1
  (4)  the parity/flip decomposition: one zero eigenvalue and three 2x2 blocks,
       and the blocks M1, M2, M3 as the text writes them
"""

from pathlib import Path
import sympy as sp

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))          # m = 3 .. -3
IDX = {m: i for i, m in enumerate(ORDER)}
a, b, L = sp.symbols("a b L", real=True)


def spin3_ops():
    Jz = sp.diag(*[sp.Integer(m) for m in ORDER])
    Jp = sp.zeros(7, 7)
    Jm = sp.zeros(7, 7)
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = sp.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = sp.sqrt(12 - m * (m - 1))
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2 * sp.I)
    return Jx, Jy, Jz, Jp, Jm


def main():
    fx, fy, fz, fp, fm = spin3_ops()

    print("(1) f^2 = fx^2+fy^2+fz^2 :")
    f2 = sp.simplify(fx * fx + fy * fy + fz * fz)
    print("    equals 12 I :", sp.simplify(f2 - 12 * sp.eye(7)) == sp.zeros(7, 7))

    # (2)+(3) the substitution
    e1 = -a / 3 + 2 * b
    e2 = -a / 3 - 2 * b
    e3 = 2 * a / 3
    print("\n(2) trace of e :", sp.simplify(e1 + e2 + e3))
    print("    e3 - 2a/3 =", sp.simplify(e3 - 2 * a / 3), "   (e1-e2) - 4b =", sp.simplify(e1 - e2 - 4 * b))
    lhs = sp.expand(e1 * fx * fx + e2 * fy * fy + e3 * fz * fz)
    rhs = sp.expand(a * (fz * fz - 4 * sp.eye(7)) + b * (fp * fp + fm * fm))
    print("    e1 fx^2 + e2 fy^2 + e3 fz^2  ==  a(fz^2-4) + b(f+^2+f-^2) :",
          sp.simplify(lhs - rhs) == sp.zeros(7, 7))

    print("\n(3) |e|^2 =", sp.simplify(sp.expand(e1**2 + e2**2 + e3**2)),
          "  == (2/3)a^2 + 8b^2 :",
          sp.simplify(e1**2 + e2**2 + e3**2 - (sp.Rational(2, 3) * a**2 + 8 * b**2)) == 0)

    # (4) the parity / flip adapted basis
    s = sp.sqrt(2)
    def vec(d):
        v = sp.zeros(7, 1)
        for m, co in d.items():
            v[IDX[m]] = co
        return v

    basis = {
        "s3 = (|3>+|-3>)/sqrt2": vec({3: 1 / s, -3: 1 / s}),
        "s1 = (|1>+|-1>)/sqrt2": vec({1: 1 / s, -1: 1 / s}),
        "t3 = (|3>-|-3>)/sqrt2": vec({3: 1 / s, -3: -1 / s}),
        "t1 = (|1>-|-1>)/sqrt2": vec({1: 1 / s, -1: -1 / s}),
        "p2 = (|2>+|-2>)/sqrt2": vec({2: 1 / s, -2: 1 / s}),
        "z0 = |0>": vec({0: 1}),
        "n2 = (|2>-|-2>)/sqrt2": vec({2: 1 / s, -2: -1 / s}),
    }
    names = list(basis)
    B = sp.Matrix.hstack(*[basis[n] for n in names])
    print("\n(4) adapted basis orthonormal :", sp.simplify(B.H * B - sp.eye(7)) == sp.zeros(7, 7))

    A = rhs
    Ab = sp.simplify(sp.expand(B.H * A * B))
    print("    A in the adapted basis (rows/cols =", names, "):")
    sp.pprint(Ab)

    M1 = Ab[0:2, 0:2]
    M2 = Ab[2:4, 2:4]
    M3 = Ab[4:6, 4:6]
    zero = Ab[6, 6]
    off = sp.zeros(7, 7)
    for i in range(7):
        for j in range(7):
            blk = lambda k: 0 if k < 2 else (1 if k < 4 else (2 if k < 6 else 3))
            if blk(i) != blk(j):
                off[i, j] = Ab[i, j]
    print("\n    block-diagonal (all inter-block entries zero) :", sp.simplify(off) == sp.zeros(7, 7))
    print("    M1 =", M1.tolist())
    print("    M2 =", M2.tolist())
    print("    M3 =", M3.tolist())
    print("    seventh diagonal entry (the 'zero eigenvalue') =", sp.simplify(zero))

    M1_text = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a + 12 * b]])
    M3_text = sp.Matrix([[0, sp.sqrt(240) * b], [sp.sqrt(240) * b, -4 * a]])
    print("\n    M1 matches the text's M1 :", sp.simplify(M1 - M1_text) == sp.zeros(2, 2))
    print("    M3 matches the text's M3 :", sp.simplify(M3 - M3_text) == sp.zeros(2, 2))
    print("    M2 char poly == M1(-b) char poly :",
          sp.simplify(sp.expand(M2.charpoly(L).as_expr()
                                - M1_text.subs(b, -b).charpoly(L).as_expr())) == 0)
    print("    (M2 as built =", M2.tolist(), ", i.e. M1(-b) conjugated by diag(1,-1))")

    with open(HERE / "s2a_01_out.txt", "w") as fh:
        fh.write("M1 = %s\nM2 = %s\nM3 = %s\nzero = %s\n"
                 % (M1.tolist(), M2.tolist(), M3.tolist(), sp.simplify(zero)))


if __name__ == "__main__":
    main()
