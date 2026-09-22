"""s2a_04_multiplicity.py -- multiplicity of L0 at each equality point, and the
(a,b) <-> e bijection.

Two things matter for grading 3c:
 (i) the map (a,b) -> e is a bijection from the ellipse onto the unit traceless
     diagonal e, so enumerating points of the ellipse really does enumerate all e;
 (ii) at each equality point, HOW MANY blocks reach L0.  The text's closing line
     asserts the top eigenspace is span{|3,3>_n, |3,-3>_n}, i.e. 2-dimensional, so
     the multiplicity of L0 must be 2 at every equality point.  Which block supplies
     the second copy is exactly the question the text's tracing skips.
"""

from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}
a, b = sp.symbols("a b", real=True)
S6 = sp.sqrt(6)
L0 = 15 / S6


def ops():
    Jz = np.diag([float(m) for m in ORDER])
    Jp = np.zeros((7, 7))
    Jm = np.zeros((7, 7))
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = np.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = np.sqrt(12 - m * (m - 1))
    return Jz, Jp, Jm


def main():
    # (i) bijection
    e1 = -a / 3 + 2 * b
    e2 = -a / 3 - 2 * b
    e3 = 2 * a / 3
    print("(i) (a,b) -> e = (%s, %s, %s)" % (e1, e2, e3))
    print("    inverse: a = 3 e3/2, b = (e1-e2)/4 ; check:",
          sp.simplify(3 * e3 / 2 - a) == 0, sp.simplify((e1 - e2) / 4 - b) == 0)
    print("    so the map is a linear BIJECTION between the ellipse and the unit")
    print("    traceless diagonal e; enumerating ellipse points enumerates all such e.\n")

    Jz, Jp, Jm = ops()
    FZ2 = Jz @ Jz
    PP = Jp @ Jp + Jm @ Jm
    L0f = float(L0)

    pts = {"(u,v)=(1/2,0)    e=(-1,-1,2)/sqrt6": (np.sqrt(6) / 2, 0.0),
           "(u,v)=(-1/4,1/8) e=(2,-1,-1)/sqrt6": (-np.sqrt(6) / 4, np.sqrt(6) / 8),
           "(u,v)=(-1/4,-1/8) e=(-1,2,-1)/sqrt6": (-np.sqrt(6) / 4, -np.sqrt(6) / 8)}

    print("(ii) full 7x7 operator A = a(fz^2-4) + b(f+^2+f-^2) at each equality point:")
    for name, (av, bv) in pts.items():
        A = av * (FZ2 - 4 * np.eye(7)) + bv * PP
        w = np.linalg.eigvalsh(A)
        mult = int(np.sum(np.abs(w - L0f) < 1e-10))
        print("  %-38s spectrum = %s" % (name, np.array2string(w, precision=9)))
        print("  %-38s lam_max = %.15f   L0 = %.15f   multiplicity of L0 = %d"
              % ("", w.max(), L0f, mult))

    print("\n  Every equality point carries L0 with multiplicity 2, matching the text's")
    print("  closing claim that the top eigenspace is 2-dimensional.")
    print("  From s2a_03, the two blocks supplying those copies are:")
    print("     (1/2,0)      : M1 and M2")
    print("     (-1/4, 1/8)  : M1 and M3")
    print("     (-1/4,-1/8)  : M2 and M3")
    print("  So at TWO of the three points M3 supplies one of the two copies, and at")
    print("  ONE of them (-1/4,-1/8) BOTH copies come from blocks (M2, M3) whose")
    print("  equality loci the text never traces.")

    # what a reader who used only the two named loci would conclude
    print("\n(iii) what the two NAMED loci alone give:")
    named = {"(1/2,0) from b=0,a>0": (np.sqrt(6) / 2, 0.0),
             "(-1/4,1/8) from a=-2b,b>0": (-np.sqrt(6) / 4, np.sqrt(6) / 8)}
    for name, (av, bv) in named.items():
        A = av * (FZ2 - 4 * np.eye(7)) + bv * PP
        w = np.linalg.eigvalsh(A)
        print("   %-28s lam_max = %.15f  multiplicity %d"
              % (name, w.max(), int(np.sum(np.abs(w - L0f) < 1e-10))))
    print("   count of e-vectors obtained: 2   (text asserts 3)")

    # the permutation NOT obtained
    av, bv = -np.sqrt(6) / 4, -np.sqrt(6) / 8
    A = av * (FZ2 - 4 * np.eye(7)) + bv * PP
    w = np.linalg.eigvalsh(A)
    print("   the missing one, e=(-1,2,-1)/sqrt6 at (a,b)=(%.12f,%.12f):" % (av, bv))
    print("     lam_max = %.15f = L0 : %s, multiplicity %d"
          % (w.max(), abs(w.max() - L0f) < 1e-12, int(np.sum(np.abs(w - L0f) < 1e-10))))


if __name__ == "__main__":
    main()
