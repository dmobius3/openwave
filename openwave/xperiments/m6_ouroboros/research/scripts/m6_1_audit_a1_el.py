"""AUDIT A1: independent EL derivation for L = -a F^2 - a G^2 + c J.A - q (J.J)^2.

Method: full 4D sympy Euler operator on 8 generic component functions of (t,x,y,z),
no gauge assumption; verify against the hand-derived covariant identity

  deltaS/deltaA_nu = eta^{nu nu} [ c J_nu + 4a ( s*box A_nu - d_nu (divA) ) ]
  deltaS/deltaJ_nu = eta^{nu nu} [ c A_nu - 4q (J.J) J_nu + 4a ( s*box J_nu - d_nu (divJ) ) ]

with box = d_t^2 - lap (the paper's Table 1 definition), s = eta_00 (i.e. d^mu d_mu = s*box),
divV = eta^{mm} d_m V_m.  Then solve for (a,c,q) that reproduce the PRINTED pair
(2.1) box A = J, (2.2) box J = A - 4g (J.J) J in Lorenz gauge, per signature.
"""
import sympy as sp

t, x, y, z = sp.symbols("t x y z")
X = [t, x, y, z]
a, c, q, g = sp.symbols("a c q g", positive=True)  # positivity irrelevant, just symbols
A = [sp.Function(f"A{i}")(*X) for i in range(4)]
J = [sp.Function(f"J{i}")(*X) for i in range(4)]


def box(f):
    return sp.diff(f, t, 2) - sp.diff(f, x, 2) - sp.diff(f, y, 2) - sp.diff(f, z, 2)


def run(name, eta_diag):
    eta = eta_diag  # list of +-1, diagonal metric, eta^ = eta_ (entries +-1)
    s = eta[0]      # d^mu d_mu = s * box  (s=+1 mostly-minus, s=-1 mostly-plus)

    def Flow(V):
        return [[sp.diff(V[n], X[m]) - sp.diff(V[m], X[n]) for n in range(4)]
                for m in range(4)]

    def sq(F):  # F_{mn} F^{mn}
        return sum(F[m][n] ** 2 * eta[m] * eta[n] for m in range(4) for n in range(4))

    FA, FJ = Flow(A), Flow(J)
    JdotA = sum(eta[m] * J[m] * A[m] for m in range(4))
    JdotJ = sum(eta[m] * J[m] ** 2 for m in range(4))
    L = -a * sq(FA) - a * sq(FJ) + c * JdotA - q * JdotJ ** 2

    def euler(fields):
        out = []
        for nu in range(4):
            e = sp.diff(L, fields[nu])
            for mu in range(4):
                e -= sp.diff(sp.diff(L, sp.Derivative(fields[nu], X[mu])), X[mu])
            out.append(sp.expand(e))
        return out

    ELA, ELJ = euler(A), euler(J)
    divA = sum(eta[m] * sp.diff(A[m], X[m]) for m in range(4))
    divJ = sum(eta[m] * sp.diff(J[m], X[m]) for m in range(4))

    okA = all(
        sp.simplify(ELA[nu] - eta[nu] * (c * J[nu] + 4 * a * (s * box(A[nu]) - sp.diff(divA, X[nu])))) == 0
        for nu in range(4))
    okJ = all(
        sp.simplify(ELJ[nu] - eta[nu] * (c * A[nu] - 4 * q * JdotJ * J[nu]
                                         + 4 * a * (s * box(J[nu]) - sp.diff(divJ, X[nu])))) == 0
        for nu in range(4))
    print(f"[{name}] covariant identity for deltaS/deltaA verified: {okA}")
    print(f"[{name}] covariant identity for deltaS/deltaJ verified: {okJ}")

    # In Lorenz gauge the EL system reads (from the verified identity, setting div=0):
    #   box A_nu = -(s c/4a) J_nu
    #   box J_nu = -(s c/4a) A_nu + (s q/a) (J.J) J_nu
    coefA = sp.simplify(-s * c / (4 * a))
    coefJq = sp.simplify(s * q / a)
    print(f"[{name}] Lorenz-gauge EL:  box A = ({coefA}) J ;  box J = ({coefA}) A + ({coefJq}) (J.J) J")

    # printed pair demands: coefA = 1 and coefJq = -4g
    sol = sp.solve([sp.Eq(coefA, 1), sp.Eq(coefJq, -4 * g), sp.Eq(a, sp.Rational(1, 4))],
                   [a, c, q], dict=True)
    print(f"[{name}] (a,c,q) reproducing printed (2.1)/(2.2) with a=1/4 kinetic norm: {sol}")

    # printed Lagrangian coefficients (a,c,q) = (1,1,g): what EL pair results?
    subs = {a: 1, c: 1, q: g}
    print(f"[{name}] EL of PRINTED L (-F^2-G^2+J.A-g(J.J)^2):  "
          f"box A = ({coefA.subs(subs)}) J ;  box J = ({coefA.subs(subs)}) A + ({coefJq.subs(subs)}) (J.J) J")
    # overall rescaling L -> k L : (a,c,q) -> (k, k, k g); show ratios invariant
    k = sp.symbols("k", positive=True)
    subsk = {a: k, c: k, q: k * g}
    print(f"[{name}] after L -> k*L:  box A = ({sp.simplify(coefA.subs(subsk))}) J ; "
          f"quartic coef = ({sp.simplify(coefJq.subs(subsk))})  -> rescaling cannot fix\n")


run("mostly-minus (+,-,-,-)", [1, -1, -1, -1])
run("mostly-plus  (-,+,+,+)", [-1, 1, 1, 1])
