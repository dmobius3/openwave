"""Exact (sympy) SU(2) part of items 2, 4, 6, 8, 10, 11.

All block quantities are SU(2)-covariant functions of the fibre vector; Gamma enters only through
the constants of s01 (|r_K|^2 and tau_J).  Structure used (proved in RETURN.md, item 2):
  |Phi|^2 = 1 + (level-12 part rho_6(u') (x) r_6),   u' = sqrt(7/d) u_hat,  |r_6|^2 = d(7-d)/7,
  int |Phi|^4 = 1 + (7(7-d)/(13 d)) rhat6(u_hat),
  ||Pi_{2J} N(Phi)||^2 = (7/d)^3 A_J(u_hat) tau_J   (J != 3),   A_J(u) = ||[rho_6(u) (x) u]_J||^2.
Here we compute rhat6, A_J, gradients and Hessians of f(u) = ||rho_6(u)||^2 exactly."""
from common import cg_exact, save_results, load_results, check
import sympy as sp
from sympy import sqrt, Rational as R, I

J3 = 3


def cgs(j1, m1, j2, m2, J, M):
    s, v2 = cg_exact(j1, m1, j2, m2, J, M)
    return 0 if s == 0 else s * sqrt(R(v2.numerator, v2.denominator))


def theta_s(x):
    j = (len(x) - 1) // 2
    return [(-1) ** (m % 2) * sp.conjugate(x[-m + j]) for m in range(-j, j + 1)]


def couple_s(x, j1, y, j2, J):
    out = []
    for M in range(-J, J + 1):
        s = 0
        for m1 in range(-j1, j1 + 1):
            m2 = M - m1
            if abs(m2) <= j2:
                c = cgs(j1, m1, j2, m2, J, M)
                if c != 0:
                    s += c * x[m1 + j1] * y[m2 + j2]
        out.append(sp.expand(s))
    return out


def norm2(x):
    return sp.expand(sum(sp.expand(a * sp.conjugate(a)) for a in x))


def rho(u, K):
    return couple_s(u, 3, theta_s(u), 3, K)


def vec(d):
    v = [0] * 7
    for m, c in d.items():
        v[m + 3] = c
    return v


def simp(x):
    return sp.radsimp(sp.simplify(x))


def ex(x):
    """exact canonical form (no numerical guessing)."""
    return sp.radsimp(sp.expand(sp.sympify(x)))


cT, sT = sqrt(13) / 5, 2 * sqrt(3) / 5          # U5: sin^2 t = 12/25
z0 = sqrt(R(23, 10))
cQ, sQ = sqrt(3) / 2, R(1, 2)                   # item 11: sin^2 t = 1/4
U = {
    'U1': vec({3: 1}),
    'U2': vec({0: 1}),
    'U3': vec({2: 1, -2: 1}),
    'U4': vec({3: 1, -3: 1}),
    'U5': vec({2: cT, -3: sT}),
    'U6': vec({3: 1, 0: z0, -3: 1}),
    'U5_item11': vec({2: cQ, -3: sQ}),
}
# unit tangent fibre vectors (section 1.6)
W = {
    'e_t': vec({2: -sT, -3: cT}),
    'e_t_item11': vec({2: -sQ, -3: cQ}),
}
u6n2 = R(43, 10)
tx = [(1 if m == 0 else 0) - z0 / u6n2 * U['U6'][m + 3] for m in range(-3, 4)]
tx = [sp.radsimp(t / sqrt(norm2(tx))) for t in tx]
W['tau_x'] = tx

# symbolic real coordinates
a = sp.symbols('a0:7', real=True)
b = sp.symbols('b0:7', real=True)
usym = [a[i] + I * b[i] for i in range(7)]
rho6 = rho(usym, 6)
f = sp.expand(sum(sp.expand(r * sp.conjugate(r)) for r in rho6))
coords = list(a) + list(b)
grad = [sp.diff(f, c) for c in coords]
hess = [[sp.diff(gc, c) for c in coords] for gc in grad]


def realvec(u):
    return [sp.re(x) for x in u] + [sp.im(x) for x in u]


def subs_at(expr, u):
    rv = realvec(u)
    return sp.radsimp(sp.expand(expr.subs(dict(zip(coords, rv)))))


out = {}
results = {}
gc = load_results('gamma_constants')
for k, u in U.items():
    n2 = norm2(u)
    uh = [sp.radsimp(x / sqrt(n2)) for x in u]
    fval = simp(norm2(rho(uh, 6)))
    rec = {'rhat6': str(fval)}
    # stationarity on the unit sphere: grad f(u_hat) = 4 f(u_hat) u_hat (f homogeneous of degree 4)
    g = [subs_at(gc_, uh) for gc_ in grad]
    rv = realvec(uh)
    tang = [sp.radsimp(gi - 4 * fval * ri) for gi, ri in zip(g, rv)]
    stationary = all(sp.simplify(t) == 0 for t in tang)
    rec['stationary'] = stationary
    rec['riemannian_grad_norm2'] = str(simp(sum(t * t for t in tang)))
    # all rho_K norms (for reference)
    rec['rhoK_norm2'] = {K: str(simp(norm2(rho(uh, K)))) for K in range(7)}
    # A_J = ||[rho_6(u_hat) (x) u_hat]_J||^2
    r6 = rho(uh, 6)
    rec['A_J'] = {J: str(simp(norm2(couple_s(r6, 6, uh, 3, J)))) for J in range(3, 10)}
    results[k] = (uh, fval, g)
    out[k] = rec
    print(k, 'rhat6 =', fval, ' stationary:', stationary, ' |Riemannian grad|^2 =', rec['riemannian_grad_norm2'])
    print('    A_J (J=3..9):', rec['A_J'])


# Hessian-based quantities
def hvec(uh):
    return [[subs_at(h, uh) for h in row] for row in hess]


def h_along(H, fval, w):
    """d^2/ds^2 f(cos s u + sin s w) at s = 0, for unit u, unit w real-orthogonal to u."""
    wv = realvec(w)
    return ex(sp.radsimp(sp.expand(
        sum(wv[i] * H[i][j] * wv[j] for i in range(14) for j in range(14)) - 4 * fval)))


def d_along(g, w):
    wv = realvec(w)
    return ex(sp.radsimp(sp.expand(sum(gi * wi for gi, wi in zip(g, wv)))))


sec_d = {'3': 3, '4': 4}
kappa = {d: R(7 * (7 - d), 52 * d) for d in (3, 4)}   # form = kappa_d * r6'' ; Q = 1 + 4 kappa_d rhat6
item6 = {}
for k, wname in (('U5', 'e_t'), ('U6', 'tau_x')):
    uh, fval, g = results[k]
    H = hvec(uh)
    w1 = W[wname]
    w2 = [I * x for x in w1]
    h1 = h_along(H, fval, w1)
    h2 = h_along(H, fval, w2)
    rec = {'r6pp_' + wname: str(h1), 'r6pp_i' + wname: str(h2)}
    if k == 'U6':
        wp = [sp.radsimp((x + y) / sqrt(2)) for x, y in zip(w1, w2)]
        wm = [sp.radsimp((x - y) / sqrt(2)) for x, y in zip(w1, w2)]
        cross = ex(sp.radsimp((h_along(H, fval, wp) - h_along(H, fval, wm)) / 2))
        rec['r6pp_cross'] = str(cross)
    for name, d in sec_d.items():
        rec['form_%s_sector%s' % (wname, name)] = str(ex(kappa[d] * h1))
        rec['form_i%s_sector%s' % (wname, name)] = str(ex(kappa[d] * h2))
        if k == 'U6':
            rec['form_cross_sector%s' % name] = str(ex(kappa[d] * cross))
    item6[k] = rec
    print(k, rec)
    # nondegeneracy: Riemannian Hessian on W = (span_R{u, iu})^perp, exact rank
    rv = realvec(uh)
    iv = realvec([I * x for x in uh])
out['item6'] = item6

# exact Riemannian Hessian ranks on W at all six points (item 10)
ranks = {}
for k in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    uh, fval, g = results[k]
    H = sp.Matrix(hvec(uh)) - 4 * fval * sp.eye(14)
    rv = sp.Matrix(realvec(uh))
    iv = sp.Matrix(realvec([I * x for x in uh]))
    P = sp.eye(14) - rv * rv.T - iv * iv.T
    M = (P * H * P).applyfunc(lambda x: ex(sp.radsimp(sp.expand(x))))
    rk = M.rank(simplify=True)
    ranks[k] = rk
    print(k, 'exact rank of the Riemannian Hessian of rhat6 on W (real dim 12):', rk)
out['hessian_rank_on_W'] = ranks

# item 11
uh, fval, g = results['U5_item11']
w = W['e_t_item11']
d1 = d_along(g, w)
out['item11'] = {'rhat6': str(fval), 'stationary': out['U5_item11']['stationary'],
                 'd_rhat6_ds_along_e_t': str(d1)}
for name, d in sec_d.items():
    # Re<Phi, DN_Phi[e_t]> = 3 Re<N(Phi), e_t> = 3 d/ds (1/4) int|Phi_s|^4 = 3 kappa_d d/ds rhat6
    out['item11']['Re_Phi_DN_et_sector' + name] = str(ex(3 * kappa[d] * d1))
print('item 11:', out['item11'])

# combine with Gamma constants: exact Q, level norms, xi norms, lambda_4
comb = {}
for k in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    rec = out[k]
    fval = ex(rec['rhat6'])
    for name, d in sec_d.items():
        tau = {J: ex(gc['tau_%d_%s' % (J, name)]) for J in range(3, 10) if 'tau_%d_%s' % (J, name) in gc}
        Q = ex(1 + 4 * kappa[d] * fval)
        lev = {}
        levxi = {}
        lam4 = 0
        for J, t in tau.items():
            if J == 3:
                continue
            nN = ex(R(7, d) ** 3 * ex(rec['A_J'][J]) * t)
            mu = 4 * J * (J + 1) - 48
            lev[2 * J] = str(nN)
            levxi[2 * J] = str(ex(nN / mu ** 2))
            lam4 += -3 * nN / mu
        comb['%s_sector%s' % (k, name)] = {'Q': str(Q), 'normPiN2': lev, 'normPixi2_over_g2': levxi,
                                            'lambda4_over_g2': str(ex(lam4))}
        print(k, 'sector', name, 'Q =', Q, ' lambda4/g^2 =', ex(lam4))
        print('     ||Pi_n xi||^2/g^2:', levxi)
out['combined'] = comb
save_results('exact_su2', out)
