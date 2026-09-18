"""Item 1: dim Hom_Gamma(sigma, V_{n/2}) for n <= 18, derived from the generators.

Method (DERIVED):
  * Gamma is enumerated from q1, q2 (s00); characters of V_j are exact: chi_j(h) = U_{2j}(Re h).
  * The two constituents of V_3 are found as the eigenspaces of a generic element of the
    commutant (1/|G|) sum_h D^3(h) M D^3(h)^dagger; their characters chi_sigma(h) = tr(P D^3(h))
    are computed numerically and identified in Z[phi] (exact), class by class.
  * dim Hom = (1/|G|) sum_h chi_sigma(h) chi_j(h) exactly in Q(sqrt5); independently, the number of
    orthonormal intertwiners found by group averaging (Reynolds operator) at 50 digits.
Also: residuals of the intertwiners at q1, q2, and the Gamma-dependent constants used later
(norms of r_K = sum_a [Theta eta_a (x) eta_a]_K and tau_J), identified at 80 and 120 digits."""
from common import *
import time

out = {}
G = group()
n = len(G)
t0 = time.time()
mp.dps = 50
S = Setup('A')
print('setup at 50 digits: %.1fs' % (time.time() - t0))
print('commutant eigenvalue clusters (value, multiplicity):',
      [(nstr(v, 8), k) for v, k in S.commutant_eigen_groups])
check('V_3 restricted to Gamma has exactly two constituents, of dimensions 3 and 4',
      sorted(k for _, k in S.commutant_eigen_groups) == [3, 4])

# class labels
cls_of = {}
reps = []
for a in G:
    if a in cls_of:
        continue
    cl = {qmul(qmul(b, a), qinv(b)) for b in G}
    for c in cl:
        cls_of[c] = len(reps)
    reps.append(a)

chars = {}
for name in ('3', '4'):
    P = S.projectors[name]
    vals = {}
    worst_classfn = mpf(0)
    for i, h in enumerate(G):
        D = S.D(3, i)
        tr = mp.re(mp.fsum(P[r][k] * D[k][r] for r in range(7) for k in range(7)))
        c = cls_of[h]
        if c in vals:
            worst_classfn = max(worst_classfn, abs(vals[c] - tr))
        else:
            vals[c] = tr
    check('sector %s: tr(P D^3(h)) is a class function' % name, worst_classfn < mpf(10) ** -40,
          nstr(worst_classfn, 3))
    exact = {}
    for c, v in vals.items():
        s, r = identify(v)
        # parse into Q5
        if abs(v) < mpf(10) ** -40:
            q = Q5(0)
        else:
            rel = mp.pslq([v, 1, mp.sqrt(5)], maxcoeff=1000, maxsteps=10000, tol=mpf(10) ** -40)
            q = Q5(Fraction(-rel[1], rel[0]), Fraction(-rel[2], rel[0]))
        check('sector %s class %d: character identified as %s' % (name, c, q),
              abs(q.to_mp() - v) < mpf(10) ** -40, nstr(abs(q.to_mp() - v), 3))
        exact[c] = q
    chars[name] = exact

# exact checks on the characters
for name, d in (('3', 3), ('4', 4)):
    ch = chars[name]
    norm = sum((ch[cls_of[h]] * ch[cls_of[h]] for h in G), Q5(0)) / n
    check('sector %s: (1/|G|) sum chi^2 = 1 exactly (irreducible)' % name, norm == Q5(1), str(norm))
    check('sector %s: chi(1) = %d' % (name, d), ch[cls_of[ONE]] == Q5(d))
tot = all(chars['3'][cls_of[h]] + chars['4'][cls_of[h]] == chebU(6, h[0]) for h in G)
check('chi_3 + chi_4 = chi_{V_3} exactly on all 120 elements', tot)
def elt_order(q):
    k, p = 1, q
    while p != ONE:
        p = qmul(p, q)
        k += 1
    return k


out['characters'] = {name: {'class_%d(order %d, Re=%s)' % (c, elt_order(reps[c]), str(reps[c][0])): str(v)
                            for c, v in chars[name].items()} for name in chars}

# exact dimensions for n = 0..18
dims = {}
for name in ('3', '4'):
    ch = chars[name]
    row = {}
    for nn in range(0, 19):
        s = sum((ch[cls_of[h]] * chebU(nn, h[0]) for h in G), Q5(0)) / n
        assert s.b == 0 and s.a.denominator == 1
        row[nn] = int(s.a)
    dims[name] = row
    print('sector %s  dim Hom(sigma, V_{n/2}), n = 0..18:' % name, [row[k] for k in range(19)])
triv = {K: int((sum((chebU(2 * K, h[0]) for h in G), Q5(0)) / n).a) for K in range(0, 13)}
print('dim V_K^Gamma for K = 0..12:', [triv[K] for K in range(13)])
out['dims'] = dims
out['invariant_dims_VK'] = triv

# numeric cross-check: number of orthonormal intertwiners found by group averaging
for name in ('3', '4'):
    sec = S.sectors[name]
    numdims = [sec.dim(J) for J in range(10)]
    check('sector %s: numerical intertwiner counts equal exact dims at even n' % name,
          numdims == [dims[name][2 * J] for J in range(10)], str(numdims))
    check('sector %s: exact dims vanish at every odd n' % name,
          all(dims[name][k] == 0 for k in range(1, 19, 2)))

# residuals of every intertwiner used later
req = mpf(10) ** -40
res = {}
for name in ('3', '4'):
    sec = S.sectors[name]
    worst_int = mpf(0)
    worst_orth = mpf(0)
    worst_cross = mpf(0)
    for J in range(10):
        for gi, eta in enumerate(sec.eta.get(J, [])):
            for idx in (S.idx_q1, S.idx_q2):
                lhs = mat_mul(S.D(J, idx), eta)
                rhs = mat_mul(eta, sec.sigma[idx])
                worst_int = max(worst_int, mat_maxdiff(lhs, rhs))
            EtE = mat_mul(mat_adj(eta), eta)
            worst_orth = max(worst_orth, mat_maxdiff(EtE, [[1 if a == b else 0 for b in range(sec.d)]
                                                            for a in range(sec.d)]))
            for gj, eta2 in enumerate(sec.eta[J]):
                if gj > gi:
                    worst_cross = max(worst_cross, max(abs(x) for r in mat_mul(mat_adj(eta), eta2) for x in r))
    # sigma real (Theta-real basis) and a homomorphism
    worst_real = max(abs(mp.im(x)) for sg in sec.sigma for r in sg for x in r)
    i12 = G.index(qmul(Q1, Q2))
    hom = mat_maxdiff(mat_mul(sec.sigma[S.idx_q1], sec.sigma[S.idx_q2]), sec.sigma[i12])
    check('sector %s: max |D^j(h) eta - eta sigma(h)|, h = q1, q2, all levels' % name, worst_int < req,
          nstr(worst_int, 3))
    check('sector %s: max |eta^dag eta - I|' % name, worst_orth < req, nstr(worst_orth, 3))
    if any(sec.dim(J) == 2 for J in range(10)):
        check('sector %s: max |eta^dag eta\'| at dimension two' % name, worst_cross < req, nstr(worst_cross, 3))
    else:
        print('sector %s: no level n <= 18 with dim Hom = 2 (no eta\' to check)' % name)
    check('sector %s: sigma(h) real (Theta eta = eta) and sigma(q1)sigma(q2) = sigma(q1 q2)' % name,
          worst_real < req and hom < req, '%s, %s' % (nstr(worst_real, 3), nstr(hom, 3)))
    res[name] = {'intertwining': nstr(worst_int, 3), 'orthonormal': nstr(worst_orth, 3),
                 'cross': nstr(worst_cross, 3) if any(sec.dim(J) == 2 for J in range(10)) else 'n/a (no level of dim 2)',
                 'sigma_imag': nstr(worst_real, 3), 'required': '1e-40', 'digits': 50}
out['residuals'] = res


# Gamma constants at two precisions
def gamma_constants():
    SS = Setup('A')
    vals = {}
    for name in ('3', '4'):
        sec = SS.sectors[name]
        for K in range(7):
            vals['r%d_norm2_%s' % (K, name)] = vnorm2(sec.rs(3, 0, 3, 0, K))
        r6 = sec.rs(3, 0, 3, 0, 6)
        for J in range(3, 10):
            if sec.dim(J) == 0:
                continue
            c = sec.mulcoef(('r6',), r6, 6, 3, 0, J)
            vals['tau_%d_%s' % (J, name)] = mpf(sec.d) / (2 * J + 1) * mp.fsum(abs(x) ** 2 for x in c)
    return vals


ids = {}
allv = {}
for p in PRECISIONS:
    with mp.workdps(p):
        t0 = time.time()
        v = gamma_constants()
        print('Gamma constants at %d digits: %.1fs' % (p, time.time() - t0))
        for k, x in v.items():
            s, r = identify(x)
            ids.setdefault(k, []).append(s)
            allv.setdefault(k, []).append(x)
gc = {}
for k, lst in ids.items():
    ok = lst[0] is not None and lst[0] == lst[1]
    check('Gamma constant %s identified identically at %d and %d digits: %s' % ((k,) + PRECISIONS + (lst[0],)), ok)
    gc[k] = lst[0]
out['gamma_constants'] = gc
save_results('item1', out)
save_results('gamma_constants', gc)
