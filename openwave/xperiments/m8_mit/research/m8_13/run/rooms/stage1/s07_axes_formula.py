"""s07_axes_formula.py -- the three-axes formula for rhat_6 on the real form.

Derivation (checked numerically here, written out in the return):

  A Theta-invariant state has an antipodally symmetric Majorana constellation,
  i.e. 3 AXES +-n_1, +-n_2, +-n_3.  Its spinor polynomial is P = q_1 q_2 q_3 with
  q_k the quadratic whose roots are the pair +-n_k, and for a unit spinor z with
  Bloch vector m one has |q_k(z)|^2 = (1 - (m.n_k)^2)/4.  With the Fock-norm
  identity ||P||_F^2 = n!(n+1) <|P|^2>_{S^3} and the Hopf pushforward of the
  normalised measure of S^3 to that of S^2,

      rhat_6 = p_6 = (13/49) * <F^2> / <F>^2 ,      F(m) = prod_k (1 - (m.n_k)^2),

  the averages taken over the unit sphere S^2.

Exact sphere moments use the Wick/pairing rule
      < (m.v_1) ... (m.v_{2k}) > = [ sum over perfect matchings of prod v_i.v_j ] / (2k+1)!!
which is checked against direct numerical integration below.
"""

from pathlib import Path
import itertools
import numpy as np
import sympy as sp

HERE = Path(__file__).parent


# ---------------------------------------------------------------- numeric side
def state_from_axes(ns):
    """spin-3 coefficient vector c_m (m = 3..-3) of the Theta-invariant state whose
       Majorana constellation is {+-n_1, +-n_2, +-n_3}."""
    # q_k(z) = (w/2) z1^2 - n_z z1 z2 - (wbar/2) z2^2 ,  w = n_x + i n_y
    poly = np.array([1.0 + 0j])  # coefficients in z1^j z2^(deg-j), index = power of z1 ascending
    for n in ns:
        w = n[0] + 1j * n[1]
        q = np.array([-np.conjugate(w) / 2, -n[2], w / 2])  # z2^2, z1z2, z1^2
        poly = np.convolve(poly, q)
    # poly[j] = coefficient of z1^j z2^(6-j)
    from math import factorial
    c = np.array([poly[j] * np.sqrt(factorial(j) * factorial(6 - j)) for j in range(7)])
    # c[j] corresponds to m = j-3 ; reorder to m = 3..-3
    c = c[::-1]
    return c / np.linalg.norm(c)


def rand_axes(rng):
    ns = []
    for _ in range(3):
        v = rng.normal(size=3)
        ns.append(v / np.linalg.norm(v))
    return ns


def F_average_numeric(ns, power, nsamp=400000, rng=None):
    rng = rng or np.random.default_rng(3)
    m = rng.normal(size=(nsamp, 3))
    m /= np.linalg.norm(m, axis=1)[:, None]
    F = np.ones(nsamp)
    for n in ns:
        F *= 1 - (m @ n) ** 2
    return float(np.mean(F ** power))


# ---------------------------------------------------------------- exact side
def wick(vlist, c):
    """< prod_i (m . v_i) > over S^2, v_i given as labels 0,1,2 with Gram entries c[(i,j)]."""
    L = len(vlist)
    if L % 2 == 1:
        return sp.Integer(0)
    k = L // 2
    dd = sp.Integer(1)
    for t in range(1, 2 * k + 2, 2):
        dd *= t  # (2k+1)!!
    tot = sp.Integer(0)
    idx = list(range(L))

    def matchings(items):
        if not items:
            yield []
            return
        a = items[0]
        for i in range(1, len(items)):
            b = items[i]
            rest = items[1:i] + items[i + 1:]
            for mm in matchings(rest):
                yield [(a, b)] + mm

    for mm in matchings(idx):
        term = sp.Integer(1)
        for (i, j) in mm:
            term *= c[(vlist[i], vlist[j])]
        tot += term
    return sp.expand(tot / dd)


def main():
    c12, c13, c23 = sp.symbols("c12 c13 c23", real=True)
    C = {(0, 0): sp.Integer(1), (1, 1): sp.Integer(1), (2, 2): sp.Integer(1),
         (0, 1): c12, (1, 0): c12, (0, 2): c13, (2, 0): c13, (1, 2): c23, (2, 1): c23}

    # <F> : F = prod (1 - x_k^2)
    print("computing <F> exactly ...")
    Fav = sp.Integer(0)
    for a in (0, 1):
        for b in (0, 1):
            for d in (0, 1):
                coef = (-1) ** (a + b + d)
                vl = [0] * (2 * a) + [1] * (2 * b) + [2] * (2 * d)
                Fav += coef * wick(vl, C)
    Fav = sp.simplify(sp.expand(Fav))
    print("  <F> =", sp.factor(Fav))

    # <F^2> : (1-x^2)^2 = 1 - 2x^2 + x^4
    print("computing <F^2> exactly (degree 12 moments) ...")
    coefs = {0: sp.Integer(1), 1: sp.Integer(-2), 2: sp.Integer(1)}
    F2av = sp.Integer(0)
    for a in (0, 1, 2):
        for b in (0, 1, 2):
            for d in (0, 1, 2):
                coef = coefs[a] * coefs[b] * coefs[d]
                vl = [0] * (2 * a) + [1] * (2 * b) + [2] * (2 * d)
                F2av += coef * wick(vl, C)
    F2av = sp.simplify(sp.expand(F2av))
    print("  <F^2> =", sp.factor(F2av))

    R = sp.simplify(sp.Rational(13, 49) * F2av / Fav ** 2)
    print("\n  rhat_6(c12,c13,c23) = (13/49) <F^2>/<F>^2")
    print("   numerator   :", sp.factor(sp.numer(sp.together(R))))
    print("   denominator :", sp.factor(sp.denom(sp.together(R))))

    with open(HERE / "s07_R.txt", "w") as fh:
        fh.write("Fav = " + sp.srepr(Fav) + "\n")
        fh.write("F2av = " + sp.srepr(F2av) + "\n")

    # ---------- checks ----------
    print("\nCHECK 1: three coincident axes (state v_0), expect rhat_6 = 100/231")
    v = R.subs({c12: 1, c13: 1, c23: 1})
    print("   R(1,1,1) =", sp.nsimplify(v), "=", sp.N(v, 25), "   100/231 =", sp.N(sp.Rational(100, 231), 25))

    print("\nCHECK 2: coplanar hexagon axes at 0, 60, 120 deg, expect 463/924")
    v = R.subs({c12: sp.Rational(1, 2), c13: sp.Rational(-1, 2), c23: sp.Rational(1, 2)})
    print("   R =", sp.nsimplify(v), "=", sp.N(v, 25), "   463/924 =", sp.N(sp.Rational(463, 924), 25))

    print("\nCHECK 3: three orthogonal axes (octahedron state)")
    v = R.subs({c12: 0, c13: 0, c23: 0})
    print("   R(0,0,0) =", sp.nsimplify(v), "=", sp.N(v, 25))

    print("\nCHECK 4: formula vs direct state computation, random axes")
    rng = np.random.default_rng(5)
    import sympy.physics.quantum.cg as _cg
    from sympy.physics.quantum.cg import CG
    ORDER = list(range(3, -4, -1))

    def cgx(m1, m2, Q, Jt):
        if m1 + m2 != Q or abs(m1) > 3 or abs(m2) > 3 or abs(Q) > Jt:
            return 0.0
        return float(sp.N(CG(3, m1, 3, -m2, Jt, Q).doit(), 25))

    T6 = []
    for Q in range(-6, 7):
        M = np.zeros((7, 7), dtype=complex)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                M[i, k] = cgx(m1, mm, Q, 6) * (-1.0) ** (3 + mm)
        T6.append(M)

    def rh6(u):
        return float(sum(abs(np.einsum("i,ij,j->", u.conj(), M, u)) ** 2 for M in T6).real)

    Rf = sp.lambdify((c12, c13, c23), R, "numpy")
    for _ in range(6):
        ns = rand_axes(rng)
        u = state_from_axes(ns)
        a = rh6(u)
        b = float(Rf(ns[0] @ ns[1], ns[0] @ ns[2], ns[1] @ ns[2]))
        print("   direct %.14f   formula %.14f   diff %.2e" % (a, b, abs(a - b)))

    print("\nCHECK 5: exact <F>,<F^2> vs Monte-Carlo (4e5 samples, ~1e-3 accuracy)")
    Fav_f = sp.lambdify((c12, c13, c23), Fav, "numpy")
    F2av_f = sp.lambdify((c12, c13, c23), F2av, "numpy")
    for _ in range(3):
        ns = rand_axes(rng)
        g = (ns[0] @ ns[1], ns[0] @ ns[2], ns[1] @ ns[2])
        print("   <F>  exact %.8f   MC %.8f" % (Fav_f(*g), F_average_numeric(ns, 1, rng=rng)))
        print("   <F^2> exact %.8f   MC %.8f" % (F2av_f(*g), F_average_numeric(ns, 2, rng=rng)))


if __name__ == "__main__":
    main()
