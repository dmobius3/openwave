"""Exact structural checks used by the arguments (items 2, 8, 9, 10, 11), each with a recorded mutation."""
import json, pickle
from fractions import Fraction as F
import sympy as sp
from mq import MQ, ZERO, ONE, mq
import su2
from su2 import couple, theta, cg, cg_table_hw, D_exact
from engine import block, N, DN, compress, block_fibre, inner, norm2, vdot, vadd, vscale, P
from linalg_mq import rank
import s_main
from s_main import U, analyze, rho, mu
from checks import check, save_log

R = {}


def to_sp(x):
    return sum((sp.Rational(a.numerator, a.denominator) + sp.I * sp.Rational(b.numerator, b.denominator)) * sp.sqrt(s)
               for s, (a, b) in x.t.items()) if not x.is_zero() else sp.Integer(0)


# ---------------------------------------------------------------- rerun analysis (populates the CG cache)
ST = {}
for name in U:
    for d in (3, 4):
        ST[(name, d)] = analyze(name, d)
# ---------------------------------------------------------------- 1. every CG coefficient used, cross-checked
# instrument the engine's coupling routine to record every (j1, j2, L) coupling it performs
triples = set()
import engine
engine.couple = lambda x, j1, y, j2, L: _couple_instr(x, j1, y, j2, L)
def _couple_instr(x, j1, y, j2, L):
    triples.add((j1, j2, L))
    return couple(x, j1, y, j2, L)
for name in ("U1", "U5", "U6"):
    for d in (3, 4):
        analyze(name, d)
mism = 0
count = 0
for (j1, j2, J) in sorted(triples):
    tab = cg_table_hw(j1, j2, J)
    for m1 in range(-j1, j1 + 1):
        for m2 in range(-j2, j2 + 1):
            M = m1 + m2
            if abs(M) > J:
                continue
            a = cg(j1, m1, j2, m2, J, M)
            b = tab.get((m1, m2, M), ZERO)
            count += 1
            if not (a - b).is_zero():
                mism += 1
print(f"CG cross-check: {len(triples)} coupling triples, {count} coefficients, mismatches: {mism}")
R["cg_crosscheck"] = {"triples": len(triples), "coefficients": count, "mismatches": mism, "max_spins": max(max(t) for t in triples)}
t0 = sorted(triples)[-1]
tab0 = cg_table_hw(*t0)
kk = next(iter(tab0))
check("Racah CG == highest-weight CG for every coefficient of every coupling used", lambda m: m == 0, mism,
      sum(1 for k in [kk] if not (cg(t0[0], kk[0], t0[1], kk[1], t0[2], kk[2]) - tab0[kk] * (-1)).is_zero()),
      "negate one highest-weight coefficient")
engine.couple = couple

# ---------------------------------------------------------------- 2. quartic identity int|Phi|^4 = 1 + beta_sigma r6hat
pts = [n for n in U]
rows = []
for n in pts:
    uh = ST[(n, 3)][1]["uhat"]
    rows.append([vdot(couple(uh, 3, uh, 3, K), couple(uh, 3, uh, 3, K)) for K in (0, 2, 4, 6)])
rk = rank(rows)
print("rank of [ ||[u(x)u]_K||^2 ]_{points x K=0,2,4,6} =", rk)
check("invariant quartic basis evaluated at the 7 points has rank 4 (evaluation injective on the 4-dim space)",
      lambda r: r == 4, rk, rank([r[:3] + [r[0]] for r in rows]), "replace K=6 column by the K=0 column")
beta = {3: F(28, 39), 4: F(21, 52)}
for d in (3, 4):
    devs = [ST[(n, d)][1]["Q"] - (ONE + ST[(n, d)][1]["r6"] * beta[d]) for n in pts]
    ok = all(x.is_zero() for x in devs)
    print(f"sector {d}: int|Phi|^4 - (1 + {beta[d]} r6hat) at 7 points:", [x.pretty() for x in devs])
    check(f"sector {d}: int|Phi|^4 == 1 + {beta[d]} r6hat at all 7 points (hence identically on the block sphere)",
          lambda b: all((ST[(n, d)][1]["Q"] - (ONE + ST[(n, d)][1]["r6"] * b)).is_zero() for n in pts), beta[d],
          beta[d] + F(1, 1000), "beta + 1/1000")
R["beta"] = {d: str(b) for d, b in beta.items()}

# ---------------------------------------------------------------- 3. r6hat along the U5 curve and the U6 family (sympy, exact)
t, x, y = sp.symbols("t x y", real=True)
def e(m):
    v = [ZERO] * 7; v[m + 3] = ONE; return v
def rr(a, b):
    return [to_sp(c) for c in rho(a, b)]
basis_pairs = {}
def rho_sym(coeffs):
    # coeffs: dict m -> sympy coefficient; rho(u) = sum_{m,n} c_m conj(c_n) [v_m (x) Theta v_n]_6
    out = [0] * 13
    for m, cm in coeffs.items():
        for n, cn in coeffs.items():
            v = rr(e(m), e(n))
            out = [o + cm * sp.conjugate(cn) * vv for o, vv in zip(out, v)]
    return out
c5 = {2: sp.cos(t), -3: sp.sin(t)}
r5 = rho_sym(c5)
r6_t = sp.simplify(sum(sp.expand(a * sp.conjugate(a)) for a in r5))
print("r6hat(u(t)) =", sp.factor(sp.simplify(sp.expand(r6_t.rewrite(sp.cos)))))
S2 = sp.symbols("S2", positive=True)  # sin^2 t
r6_S = sp.simplify(sp.expand(sp.expand_trig(r6_t)).subs(sp.cos(t) ** 2, 1 - S2).subs(sp.sin(t) ** 2, S2))
r6_S = sp.simplify(r6_S.subs(sp.cos(t), sp.sqrt(1 - S2)).subs(sp.sin(t), sp.sqrt(S2)))
print("r6hat along U5 curve as function of S = sin^2 t:", sp.factor(r6_S))
crit5 = sp.solve(sp.diff(r6_S, S2), S2)
print("critical S in (0,1):", crit5)
R["r6hat_U5curve"] = str(sp.factor(r6_S)); R["r6hat_U5curve_critical_S"] = [str(c) for c in crit5]
check("d r6hat/d(sin^2 t) vanishes at sin^2 t = 12/25", lambda s_: sp.simplify(sp.diff(r6_S, S2).subs(S2, s_)) == 0,
      sp.Rational(12, 25), sp.Rational(1, 4), "evaluate at sin^2 t = 1/4 instead")
check("r6hat(sin^2 t = 12/25) == 9/35 (matches exact engine)", lambda v: sp.simplify(v - sp.Rational(9, 35)) == 0,
      r6_S.subs(S2, sp.Rational(12, 25)), r6_S.subs(S2, sp.Rational(1, 4)), "value at 1/4")
# U6 family u(z) = v3 + z v0 + v-3, z = x + i y
z = x + sp.I * y
c6 = {3: sp.Integer(1), 0: z, -3: sp.Integer(1)}
r6v = rho_sym(c6)
num = sp.expand(sum(sp.expand(a * sp.conjugate(a)) for a in r6v))
den = sp.expand((2 + x ** 2 + y ** 2) ** 2)
f6 = sp.simplify(num / den)
print("r6hat(u(z)) =", sp.factor(f6))
fx = sp.simplify(f6.subs(y, 0))
crit6 = sp.solve(sp.diff(fx, x), x)
print("critical real x:", crit6)
R["r6hat_U6family"] = str(sp.factor(f6)); R["r6hat_U6family_critical_x"] = [str(c) for c in crit6]
z0 = sp.sqrt(sp.Rational(23, 10))
check("d r6hat/dx vanishes at z = z0 = sqrt(23/10)", lambda v: sp.simplify(sp.diff(f6, x).subs({x: v, y: 0})) == 0, z0, z0 + sp.Rational(1, 10), "z0 + 1/10")
check("d r6hat/dy vanishes at z = z0 (y = 0)", lambda v: sp.simplify(sp.diff(f6, y).subs({x: z0, y: v})) == 0, 0, sp.Rational(1, 10), "y = 1/10")
check("r6hat(z0) == 200/903", lambda v: sp.simplify(v - sp.Rational(200, 903)) == 0, f6.subs({x: z0, y: 0}), f6.subs({x: 1, y: 0}), "value at z = 1")

# ---------------------------------------------------------------- 4. lambda_4 identity and sign
for (n, d), (Rr, S_) in ST.items():
    lam = S_["Gf"]
    along = vdot(S_["u"], S_["Gf"]) * F(d, 7)
    levels = [J for J in S_["Nf"] if J != 3]
    ident = sum((norm2({J: S_["Nf"][J]}) * F(-3, mu(J)) for J in levels), ZERO)
    wrong = sum((norm2({J: S_["Nf"][J]}) * F(-3, mu(J)) for J in levels[1:]), ZERO)
    check(f"{n} s{d}: <Phi, DN xi>/g == -3 sum_n ||Pi_n N||^2/(n(n+2)-48)", lambda v: (along - v).is_zero(), ident, wrong, "drop the lowest level from the sum")
    if n != "U5q":
        assert all(mu(J) > 0 for J in levels)

# ---------------------------------------------------------------- 5. antiunitary symmetry
# (a) Theta_col(conj P) has columns in range P  (=> T(Phi_u) := conj(Phi_u) W = Phi_{Theta u} with W in Hom(sigma, conj sigma))
for d in (3, 4):
    Pm = P[d]
    Pp = [[Pm[-k + 3][a].conj() * ((-1) ** (k % 2)) for a in range(7)] for k in range(-3, 4)]  # Pp[k][a] = (-1)^k conj P[-k][a]
    PPp = [[sum((Pm[i][k] * Pp[k][a] for k in range(7)), ZERO) for a in range(7)] for i in range(7)]
    ok = all((PPp[i][a] - Pp[i][a]).is_zero() for i in range(7) for a in range(7))
    other = P[7 - d]
    Po = [[sum((other[i][k] * Pp[k][a] for k in range(7)), ZERO) for a in range(7)] for i in range(7)]
    check(f"sector {d}: P Theta_col(conj P) == Theta_col(conj P) (conjugate sector coincides with the sector)",
          lambda M: all((M[i][a] - Pp[i][a]).is_zero() for i in range(7) for a in range(7)), PPp, Po, "use the other sector's projector")
# (b) r0 = exp(-i pi J_y) (quaternion j -> [[0,-1],[1,0]]): D^J(r0) Theta x == (-1)^J conj(x)  for J = 0..9
okall = True
for J in range(0, 10):
    Dr = D_exact(J, [[ZERO, -ONE], [ONE, ZERO]])
    for m in range(-J, J + 1):
        for ph in (ONE, MQ.I):
            xv = [ZERO] * (2 * J + 1); xv[m + J] = ph
            lhs = [sum((Dr[i][k] * theta(xv, J)[k] for k in range(2 * J + 1)), ZERO) for i in range(2 * J + 1)]
            rhs = [c.conj() * ((-1) ** J) for c in xv]
            okall &= all((a - b).is_zero() for a, b in zip(lhs, rhs))
check("D^J(r0) Theta x == (-1)^J conj(x) for all J <= 9 (exact)", lambda v: v, okall,
      all((a - b).is_zero() for a, b in zip([sum((D_exact(3, [[ZERO, -ONE], [ONE, ZERO]])[i][k] * theta(e(2), 3)[k] for k in range(7)), ZERO) for i in range(7)], e(2))),
      "compare with +conj(x) at J=3")
# (c) S' = -(r0 o T) acts on the block by u -> conj(u): check that pi_6 DN xi fibres (Gf), kappa are real at U1..U6
for (n, d), (Rr, S_) in ST.items():
    if n == "U5q":
        continue
    realG = all(c.is_real() for c in S_["Gf"]) and all(c.is_real() for c in S_["kappa"]) and all(c.is_real() for c in S_["u"])
    check(f"{n} s{d}: fibres of Phi, Pi_6 DN_Phi[xi], kappa are real (S'-invariant)", lambda v: v, realG,
          all(c.is_real() for c in vscale(S_["Gf"], MQ.I)) and any(not c.is_zero() for c in S_["Gf"]), "multiply the DN xi fibre by i")

# ---------------------------------------------------------------- 6. Sym^3 mechanism for the level-16 zeros at U2, U3, U4
for n in ("U2", "U3", "U4"):
    uh = ST[(n, 3)][1]["uhat"]
    th_ = theta(uh, 3)
    c = vdot(uh, th_)
    prop = all((a - c * b).is_zero() for a, b in zip(th_, uh))
    check(f"{n}: Theta u == c u with |c| = 1 (c = {c.pretty()})", lambda v: v, prop,
          all((a - c * b).is_zero() for a, b in zip(theta(ST[("U1", 3)][1]["uhat"], 3), ST[("U1", 3)][1]["uhat"])), "use U1 instead")
# multiplicity of V_8 in Sym^3 V_3 by weight counting
from itertools import combinations_with_replacement as cwr
cnt = {}
for tpl in cwr(range(-3, 4), 3):
    cnt[sum(tpl)] = cnt.get(sum(tpl), 0) + 1
mult = {J: cnt.get(J, 0) - cnt.get(J + 1, 0) for J in range(0, 10)}
print("multiplicities of V_J in Sym^3 V_3:", mult)
R["Sym3_V3_multiplicities"] = mult
check("V_8 does not occur in Sym^3 V_3", lambda m: m[8] == 0, mult, {**mult, 8: 1}, "set mult 1")

# ---------------------------------------------------------------- 7. item 11: Re<Phi, DN_Phi[e_t]> at sin^2 t = 1/4 (and at U5)
for n in ("U5q", "U5"):
    for d in (3, 4):
        Rr, S_ = ST[(n, d)]
        et = S_["tang"]["e_t"]
        Phi = block(S_["u"], d)
        Ef = block(et, d)
        val = inner(Phi, compress(DN(Phi, Ef, levels={3}))).re()
        alt = inner(block(S_["Fu"], d), Ef).re() * 3        # = 3 <N(Phi), e_t>_R  (symmetry of DN)
        R[f"item11_Re<Phi,DN e_t>_{n}_s{d}"] = val.pretty()
        print(f"{n} sector {d}: Re<Phi, DN_Phi[e_t]> =", val.pretty(), " (3 Re<Pi_6 N(Phi), e_t> =", alt.pretty(), ")")
        check(f"{n} s{d}: Re<Phi,DN e_t> == 3 Re<N(Phi), e_t> (two computations)", lambda v: (val - v).is_zero(), alt, alt + F(1, 10**9), "+1e-9")

with open("out/verify.json", "w") as f:
    json.dump(R, f, indent=1, default=str)
save_log("verify")
