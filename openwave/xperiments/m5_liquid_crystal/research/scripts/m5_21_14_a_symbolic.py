"""M5.21.14 S1: the first nontrivial (1/g) term of the boost-hedgehog
dressing, derived symbolically and gated.

The author's recipe (m5_21_convo.md 2026-08-09): work symbolically in
the 4x4 case with the spherical boost hedgehog
Qb = MatrixExp[b(r) rhat . boostgenerators], find the first nontrivial
term in the (1/g) expansion, to be included in the 3x3 case; the
pre-stated criterion: the term needs a NEGATIVE Hamiltonian
contribution to get oscillations.

Conventions = the certified 4D instrument (m5_21_3_a_4d.py):
eta = diag(-1,1,1,1); E_u density = 4 sum_{i<j} <F_ij,F_ij>_eta with
F_ij = dM eta dM commutator; kin density = 4 sum_i <F_0i,F_0i>_eta
with a0 = dM/dt; embed34 puts sigma = -s*g in the (0,0) slot.
Boost: Qb = I + sinh(b) K + (cosh(b)-1) K^2, K = rhat.(sym boosts).

Stages (each a machine gate, JSON verdicts):
  V0  lemmas: Qb eta Qb^T = eta (exact, symbolic) and
      Qb(b) Qb(-b) = I  =>  E_V invariant under ANY b(r) (similarity);
      spot-checked exactly on the brute instance at full b.
  T1  structured block expansion at leading order in eps = 1/g with
      b = eps*beta(r) (the m* ~ 1/g regime): claimed forms
        A_i   = blockdiag(0, G_i) + st*(e0 v_i^T + v_i e0^T) + O(eps)
        v_i   = d_i(beta nhat),  G_i = d_i M3,  st = -s (st^2 = 1)
        T1_u  = 4 sum_{i<j} [ 2<[G_i,G_j], v_j v_i^T - v_i v_j^T>
                + |v_j v_i^T - v_i v_j^T|^2 - 2|G_i v_j - G_j v_i|^2 ]
        T1_k  = -8 sum_i |Mdot3 v_i|^2        (omega^2 multiplier)
      verified as EXACT matrix-algebra identities with FREE v_i.
  BR  brute-force check: concrete random-rational M3(x), beta(r),
      exact sympy differentiation of the FULL dressed field, eps-series
      at a random off-axis point, both s signs: the eps^0 coefficient
      minus the undressed 3x3 density == T1 evaluated there (50-digit
      arithmetic, tol 1e-30 relative).
  ORD the ordering question: at FIXED b the density is polynomial in
      sigma (deg <= 4 possible); report the measured leading power ->
      the (1/g) expansion is only nontrivial in the scaled variable
      beta = g*b (the regime the author's m* ~ 1/g law selects).

Out: ../data/m5_21_14_symbolic.json
"""
from __future__ import annotations

import json
import os
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

ETA = sp.diag(-1, 1, 1, 1)
E0 = sp.Matrix([1, 0, 0, 0])


def frob2(A):
    """|A|_F^2."""
    return sp.trace(A * A.T)


def inner_eta(F, G):
    """<F,G>_eta = tr(eta F eta G^T)."""
    return sp.trace(ETA * F * ETA * G.T)


def comm_eta(A, B):
    return A * ETA * B - B * ETA * A


def boost_k(n1, n2, n3):
    K = sp.zeros(4, 4)
    K[0, 1], K[0, 2], K[0, 3] = n1, n2, n3
    K[1, 0], K[2, 0], K[3, 0] = n1, n2, n3
    return K


def qb_of(b, n1, n2, n3):
    K = boost_k(n1, n2, n3)
    K2 = K * K
    return sp.eye(4) + sp.sinh(b) * K + (sp.cosh(b) - 1) * K2


def sym3(name):
    a = sp.symbols(f"{name}11 {name}12 {name}13 {name}22 {name}23 {name}33")
    return sp.Matrix([[a[0], a[1], a[2]],
                      [a[1], a[3], a[4]],
                      [a[2], a[4], a[5]]])


def blockdiag03(G3):
    M = sp.zeros(4, 4)
    M[1:4, 1:4] = G3
    return M


def timespace(v):
    """e0 v^T + v e0^T with v a spatial 3-vector."""
    v4 = sp.Matrix([0, v[0], v[1], v[2]])
    return E0 * v4.T + v4 * E0.T


# ================= V0: the lemmas =================
def stage_v0():
    b, n1, n2 = sp.symbols("b n1 n2", real=True)
    n3 = sp.sqrt(1 - n1 ** 2 - n2 ** 2)
    Qb = qb_of(b, n1, n2, n3)
    lem1 = sp.simplify(Qb * ETA * Qb.T - ETA)
    lem2 = sp.simplify(Qb * qb_of(-b, n1, n2, n3) - sp.eye(4))
    return {"eta_preserved": lem1 == sp.zeros(4, 4),
            "inverse_is_minus_b": lem2 == sp.zeros(4, 4)}


# ================= T1: structured identity, free v_i =================
def stage_t1():
    st = sp.symbols("st")  # st = -s, st^2 = 1
    G = [sym3(f"g{i}") for i in (1, 2, 3)]
    V = [sp.Matrix(sp.symbols(f"v{i}1 v{i}2 v{i}3")) for i in (1, 2, 3)]
    Md3 = sym3("md")

    A = [blockdiag03(G[i]) + st * timespace(V[i]) for i in range(3)]
    A0 = blockdiag03(Md3)

    dens_u = sp.Integer(0)
    for i in range(3):
        for j in range(i + 1, 3):
            F = comm_eta(A[i], A[j])
            dens_u += 4 * inner_eta(F, F)
    dens_u = sp.expand(dens_u).subs(st ** 2, 1)

    dens_k = sp.Integer(0)
    for i in range(3):
        F = comm_eta(A0, A[i])
        dens_k += 4 * inner_eta(F, F)
    dens_k = sp.expand(dens_k).subs(st ** 2, 1)

    # the undressed 3x3 densities (plain commutators, eta trivial there)
    base_u = sp.expand(sum(
        4 * frob2(G[i] * G[j] - G[j] * G[i])
        for i in range(3) for j in range(i + 1, 3)))
    base_k = sp.expand(sum(
        4 * frob2(Md3 * G[i] - G[i] * Md3) for i in range(3)))

    # claimed first nontrivial terms
    t1_u = sp.Integer(0)
    for i in range(3):
        for j in range(i + 1, 3):
            C = G[i] * G[j] - G[j] * G[i]
            W = V[j] * V[i].T - V[i] * V[j].T
            w = G[i] * V[j] - G[j] * V[i]
            t1_u += 4 * (2 * sp.trace(C * W.T) + frob2(W)
                         - 2 * (w.T * w)[0, 0])
    t1_u = sp.expand(t1_u)

    t1_k = sp.expand(-8 * sum(((Md3 * V[i]).T * (Md3 * V[i]))[0, 0]
                              for i in range(3)))

    ok_u = sp.expand(dens_u - base_u - t1_u) == 0
    ok_k = sp.expand(dens_k - base_k - t1_k) == 0
    # st (= sign-knob) independence at leading order
    ok_st = (sp.expand(dens_u).has(st) is False
             and sp.expand(dens_k).has(st) is False)
    return {"t1_static_identity": bool(ok_u),
            "t1_kin_identity": bool(ok_k),
            "sign_knob_free_at_leading_order": bool(ok_st),
            "t1_static_str": ("4*sum_{i<j}[ 2<[G_i,G_j], v_j v_i^T - "
                              "v_i v_j^T> + |v_j v_i^T - v_i v_j^T|^2 "
                              "- 2|G_i v_j - G_j v_i|^2 ]"),
            "t1_kin_str": "-8*sum_i |Mdot3 v_i|^2  (omega^2 multiplier)"}


# ================= BR: brute-force eps-series check =================
def _brute_instance(sign_s, seed_shift=0):
    """exact dressed construction at a concrete random-rational
    instance; returns dict of exact eps^0 residuals (50-digit)."""
    x, y, z, eps = sp.symbols("x y z eps", positive=False)
    r = sp.sqrt(x * x + y * y + z * z)
    n1, n2, n3 = x / r, y / r, z / r

    # concrete rational coefficient pools (deterministic, seed-shifted)
    def q(k):
        vals = [sp.Rational(p, 7 + ((p + k + seed_shift) % 5))
                for p in (2, -3, 1, 4, -1, 3, -2, 5, 1, -4)]
        return vals[(k + seed_shift) % len(vals)]

    xs = (x, y, z)

    def rand_sym3(off):
        M = sp.zeros(3, 3)
        c = 0
        for a in range(3):
            for bb in range(a, 3):
                e = q(off + c)
                c += 1
                for k in range(3):
                    e += q(off + 10 + c + k) * xs[k]
                    c += 1
                M[a, bb] = e
                M[bb, a] = e
        return M

    M3 = rand_sym3(0)
    Md3 = rand_sym3(97)
    beta = q(51) + q(52) * (x * x + y * y + z * z)  # beta(r), poly in r^2

    st = -sp.Integer(sign_s)          # sigma = -s*g = st/eps
    sigma = st / eps
    b = eps * beta
    Qb = qb_of(b, n1, n2, n3)
    M4 = blockdiag03(M3)
    M4[0, 0] = sigma
    Mdr = Qb * M4 * Qb.T
    A = [sp.diff(Mdr, w) for w in xs]
    a0 = Qb * blockdiag03(Md3) * Qb.T

    pt = {x: sp.Rational(2, 3), y: sp.Rational(-1, 2),
          z: sp.Rational(3, 4)}

    # E_V invariance spot-check at FULL b (eps = 1 here, any value ok)
    ev_res = []
    Me_d = (Mdr * ETA).subs(pt).subs(eps, 1)
    Me_0 = (M4 * ETA).subs(pt).subs(eps, 1)
    P_d, P_0 = sp.eye(4), sp.eye(4)
    for p in range(1, 5):
        P_d, P_0 = P_d * Me_d, P_0 * Me_0
        ev_res.append(abs(sp.N(sp.trace(P_d) - sp.trace(P_0), 50)))

    # the exact densities at the point, as series in eps
    Apt = [sp.N(Ai.subs(pt), 50) for Ai in A]
    a0pt = sp.N(a0.subs(pt), 50)

    def dens_u_of(Alist):
        tot = sp.Integer(0)
        for i in range(3):
            for j in range(i + 1, 3):
                F = comm_eta(Alist[i], Alist[j])
                tot += 4 * inner_eta(F, F)
        return tot

    def dens_k_of(a0m, Alist):
        tot = sp.Integer(0)
        for i in range(3):
            F = comm_eta(a0m, Alist[i])
            tot += 4 * inner_eta(F, F)
        return tot

    # series in eps: entries of Apt are numeric-coefficient functions
    # of eps only; series each density scalar directly
    du = sp.series(sp.expand(dens_u_of(Apt)), eps, 0, 2).removeO()
    dk = sp.series(sp.expand(dens_k_of(a0pt, Apt)), eps, 0, 2).removeO()
    du0 = du.coeff(eps, 0)
    dk0 = dk.coeff(eps, 0)
    neg_pows = [du.coeff(eps, p) for p in (-4, -3, -2, -1)] + \
               [dk.coeff(eps, p) for p in (-4, -3, -2, -1)]
    neg_ok = all(sp.N(abs(c), 50) < sp.Float("1e-30") for c in neg_pows
                 if c is not None)

    # the claimed T1 at the point
    G_pt = [sp.N(sp.diff(M3, w).subs(pt), 50) for w in xs]
    bn = sp.Matrix([beta * n1, beta * n2, beta * n3])
    V_pt = [sp.N(sp.diff(bn, w).subs(pt), 50) for w in xs]
    M3pt = sp.N(M3.subs(pt), 50)
    Md3pt = sp.N(Md3.subs(pt), 50)

    base_u = sum(4 * frob2(G_pt[i] * G_pt[j] - G_pt[j] * G_pt[i])
                 for i in range(3) for j in range(i + 1, 3))
    base_k = sum(4 * frob2(Md3pt * G_pt[i] - G_pt[i] * Md3pt)
                 for i in range(3))
    t1_u = sp.Integer(0)
    for i in range(3):
        for j in range(i + 1, 3):
            C = G_pt[i] * G_pt[j] - G_pt[j] * G_pt[i]
            W = V_pt[j] * V_pt[i].T - V_pt[i] * V_pt[j].T
            w = G_pt[i] * V_pt[j] - G_pt[j] * V_pt[i]
            t1_u += 4 * (2 * sp.trace(C * W.T) + frob2(W)
                         - 2 * (w.T * w)[0, 0])
    t1_k = -8 * sum(((Md3pt * V_pt[i]).T * (Md3pt * V_pt[i]))[0, 0]
                    for i in range(3))

    res_u = sp.N(du0 - base_u - t1_u, 50)
    res_k = sp.N(dk0 - base_k - t1_k, 50)
    scale_u = max(abs(sp.N(base_u, 50)), abs(sp.N(t1_u, 50)), 1)
    scale_k = max(abs(sp.N(base_k, 50)), abs(sp.N(t1_k, 50)), 1)
    return {"s": sign_s,
            "ev_invariance_max_abs": float(max(ev_res)),
            "no_negative_eps_powers": bool(neg_ok),
            "t1_u_rel_residual": float(abs(res_u) / scale_u),
            "t1_k_rel_residual": float(abs(res_k) / scale_k),
            "t1_u_value": float(sp.N(t1_u, 30)),
            "t1_k_value": float(sp.N(t1_k, 30))}


def stage_brute():
    return [_brute_instance(-1), _brute_instance(+1, seed_shift=3)]


# ================= ORD: fixed-b sigma-polynomial =================
def stage_order():
    """at FIXED b, collect the density's polynomial degree in sigma."""
    x, y, z, sg = sp.symbols("x y z sigma")
    r = sp.sqrt(x * x + y * y + z * z)
    n1, n2, n3 = x / r, y / r, z / r
    b = sp.Rational(1, 3) + sp.Rational(1, 5) * (x * x + y * y + z * z)
    Qb = qb_of(b, n1, n2, n3)
    M3 = sp.diag(1, sp.Rational(3, 10), 0)  # vacuum-like spatial block
    M4 = blockdiag03(M3)
    M4[0, 0] = sg
    Mdr = Qb * M4 * Qb.T
    xs = (x, y, z)
    A = [sp.diff(Mdr, w) for w in xs]
    pt = {x: sp.Rational(2, 3), y: sp.Rational(-1, 2),
          z: sp.Rational(3, 4)}
    Apt = [sp.N(Ai.subs(pt), 30) for Ai in A]
    dens = sp.Integer(0)
    for i in range(3):
        for j in range(i + 1, 3):
            F = comm_eta(Apt[i], Apt[j])
            dens += 4 * inner_eta(F, F)
    poly = sp.Poly(sp.expand(dens), sg)
    coeffs = {f"sigma^{m[0]}": float(abs(sp.N(c, 30)))
              for m, c in zip(poly.monoms(), poly.coeffs())}
    deg = poly.degree()
    return {"sigma_degree_at_fixed_b": int(deg),
            "abs_coeffs": coeffs,
            "reading": ("at fixed b the static density GROWS as "
                        f"sigma^{deg}: the (1/g) expansion is only "
                        "nontrivial in the scaled variable "
                        "beta = g*b, the regime the author's "
                        "m* ~ 1/g law selects")}


def main():
    t0 = time.time()
    out = {}
    out["V0"] = stage_v0()
    print(json.dumps({"V0": out["V0"]}), flush=True)
    out["T1"] = stage_t1()
    print(json.dumps({"T1": {k: v for k, v in out["T1"].items()
                             if not k.endswith("_str")}}), flush=True)
    out["BR"] = stage_brute()
    print(json.dumps({"BR": out["BR"]}), flush=True)
    out["ORD"] = stage_order()
    print(json.dumps({"ORD": {"sigma_degree_at_fixed_b":
                              out["ORD"]["sigma_degree_at_fixed_b"]}}),
          flush=True)
    gates = {
        "V0": all(out["V0"].values()),
        "T1": (out["T1"]["t1_static_identity"]
               and out["T1"]["t1_kin_identity"]
               and out["T1"]["sign_knob_free_at_leading_order"]),
        "BR": all(b["no_negative_eps_powers"]
                  and b["ev_invariance_max_abs"] < 1e-30
                  and b["t1_u_rel_residual"] < 1e-30
                  and b["t1_k_rel_residual"] < 1e-30
                  for b in out["BR"]),
    }
    out["gates"] = gates
    out["all_green"] = all(gates.values())
    out["runtime_s"] = round(time.time() - t0, 1)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_21_14_symbolic.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"gates": gates, "all_green": out["all_green"],
                      "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
