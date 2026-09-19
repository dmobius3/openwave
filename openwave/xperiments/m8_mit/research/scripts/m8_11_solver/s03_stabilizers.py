"""Item 3: rotations fixing Phi up to a phase, their character, and the dimension of the
character subspace of the block (and of every V_J, J <= 9, used in item 9).

Rotations act on block sections through the fibre vector: r.Phi_{sigma,u} = Phi_{sigma, rho(r) u},
rho(r) = conj(D^3(r)) (s00 checks conj D^j = Theta-conjugate of D^j).  So the stabilizer does not
depend on the sector.  Method:
  * Lie algebra of the stabilizer: real solutions c of (sum_a c_a X_a) u in C u, X_a = i conj(J_a);
  * finite part: Majorana roots.  With u <-> P_u(x,y) = sum u_m x^(3+m) y^(3-m)/sqrt((3+m)!(3-m)!),
    D(g) acts by P -> P((x,y) g), so the root multiset of D(g)u is (roots of u) g^-1.  A rotation
    fixing [u] must permute the roots; we enumerate every g in SU(2) (mod +-1) sending a fixed pair
    of roots to a pair of roots, and keep those that fix [u].  This search is complete.
  * the closed-form generators of common.stabilizer_generators() are checked to generate exactly
    the group found by the search."""
from common import *

mp.dps = 30
TOL = mpf(10) ** -20
U = fibre_vectors()
out = {}


def conj_vec(u):
    return [mp.conj(x) for x in u]


def roots(u):
    """Majorana roots of u under D, as unit row vectors (x, y)."""
    j = 3
    c = [u[m + j] / mp.sqrt(fac(j + m) * fac(j - m)) for m in range(-j, j + 1)]
    nz = [i for i in range(7) if abs(c[i]) > TOL]
    deg, low = max(nz), min(nz)
    rts = []
    if deg > low:
        rts = mp.polyroots([c[i] for i in range(deg, low - 1, -1)], maxsteps=400, extraprec=200)
    pts = [[r, mpc(1)] for r in rts] + [[mpc(0), mpc(1)]] * low + [[mpc(1), mpc(0)]] * (6 - deg)
    return [vscale(p, 1 / mp.sqrt(vnorm2(p))) for p in pts]


def same_multiset(A, B):
    B = list(B)
    for p in A:
        for i, q in enumerate(B):
            if abs(abs(p[0] * mp.conj(q[0]) + p[1] * mp.conj(q[1])) - 1) < TOL:
                B.pop(i)
                break
        else:
            return False
    return True


def rowmul(p, W):
    return [p[0] * W[0][0] + p[1] * W[1][0], p[0] * W[0][1] + p[1] * W[1][1]]


def inv2(g):
    return [[g[1][1], -g[0][1]], [-g[1][0], g[0][0]]]


def fixes(u, g):
    v = mat_vec(rho_left(3, g), u)
    c = vdot(u, v) / vnorm2(u)
    return mp.sqrt(vnorm2(vadd(v, u, -c))) < TOL, c


# Majorana equivariance check at a generic element (this can fail if the convention were wrong)
gr = rot((1, 2, 3), mpf('0.7'))
ut = [mpc(k + 1, 3 - k) for k in range(7)]
check('Majorana roots are equivariant: roots(D(g)u) = roots(u) g^-1',
      same_multiset(roots(mat_vec(wigner_D(6, gr), ut)), [rowmul(p, inv2(gr)) for p in roots(ut)]))

Jx, Jy, Jz = spin_matrices(3)
gens_cf = stabilizer_generators()


def lie_stabilizer_dim(u):
    """dim of {c in R^3 : (sum c_a i conj(J_a)) u in C u}."""
    un = vscale(u, 1 / mp.sqrt(vnorm2(u)))
    cols = []
    for Ja in (Jx, Jy, Jz):
        v = mat_vec(mat_conj(Ja), un)
        v = vscale(v, mpc(0, 1))
        v = vadd(v, un, -vdot(un, v))          # remove the C u component
        cols.append([mp.re(x) for x in v] + [mp.im(x) for x in v])
    A = mp.matrix(14, 3)
    for k in range(3):
        for i in range(14):
            A[i, k] = cols[k][i]
    sv = mp.svd_r(A, compute_uv=False)
    return sum(1 for i in range(3) if sv[i] < TOL)


def close_group(gens):
    els = [[[mpc(1), mpc(0)], [mpc(0), mpc(1)]]]
    frontier = list(els)
    while frontier:
        new = []
        for e in frontier:
            for g in gens:
                h = mat_mul(e, g)
                if not any(mat_maxdiff(h, x) < TOL for x in els):
                    els.append(h)
                    new.append(h)
                    if len(els) > 500:
                        return els
        frontier = new
    return els


def isotypic_dim(J, gens):
    """dim {x in V_J : rho(g) x = chi(g) x for all generators}."""
    n = 2 * J + 1
    rows = []
    for g, chi in gens:
        Rm = rho_left(J, g)
        for i in range(n):
            rows.append([Rm[i][k] - (chi if i == k else 0) for k in range(n)])
    A = mp.matrix(rows)
    sv = mp.svd_c(A, compute_uv=False)
    return sum(1 for i in range(n) if sv[i] < TOL)


for k in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    u = U[k]
    rec = {}
    ldim = lie_stabilizer_dim(u)
    rec['lie_algebra_dim'] = ldim
    R_ = roots(conj_vec(u))          # stabilizer of [u] under rho = stabilizer of [conj u] under D
    distinct = []
    for p in R_:
        if not any(abs(abs(p[0] * mp.conj(q[0]) + p[1] * mp.conj(q[1])) - 1) < TOL for q in distinct):
            distinct.append(p)
    rec['distinct_roots'] = len(distinct)
    # generators check
    for g, chi in gens_cf[k]:
        ok, c = fixes(u, g)
        check('%s: closed-form generator fixes [u] with character %s' % (k, nstr(chi, 8)),
              ok and abs(c - chi) < TOL, nstr(c, 8))
    if ldim == 0:
        # complete enumeration through a non-orthogonal pair of distinct roots
        pair = None
        for a_ in range(6):
            for b_ in range(6):
                ip = abs(R_[a_][0] * mp.conj(R_[b_][0]) + R_[a_][1] * mp.conj(R_[b_][1]))
                if TOL < ip < 1 - TOL:
                    pair = (a_, b_)
                    break
            if pair:
                break
        p1, p2 = R_[pair[0]], R_[pair[1]]
        ip = p1[0] * mp.conj(p2[0]) + p1[1] * mp.conj(p2[1])
        found = []
        for a_ in range(6):
            for b_ in range(6):
                q1, q2 = R_[a_], R_[b_]
                iq = q1[0] * mp.conj(q2[0]) + q1[1] * mp.conj(q2[1])
                if abs(abs(ip) - abs(iq)) > TOL:
                    continue
                ph = ip / iq
                q2p = vscale(q2, mp.conj(ph) / abs(ph))
                P = mp.matrix([p1, p2])
                Qm = mp.matrix([q1, q2p])
                Wm = mp.inverse(P) * Qm
                Wm = Wm / mp.sqrt(mp.det(Wm))
                g = inv2([[Wm[0, 0], Wm[0, 1]], [Wm[1, 0], Wm[1, 1]]])
                ok, c = fixes(u, g)
                if ok and not any(mat_maxdiff(g, h) < TOL or mat_maxdiff(g, [[-x for x in r] for r in h]) < TOL
                                  for h, _ in found):
                    found.append((g, c))
        rec['order_SO3'] = len(found)
        grp = close_group([g for g, _ in gens_cf[k]])
        rec['order_SU2_from_generators'] = len(grp)
        check('%s: closed-form generators generate a group of order %d in SU(2) = 2 x (search count %d)'
              % (k, len(grp), len(found)), len(grp) == 2 * len(found))
        inside = all(any(mat_maxdiff(g, h) < TOL for h in grp) for g, _ in found)
        check('%s: every rotation found by the complete search lies in the generated group' % k, inside)
        # character is a homomorphism on the generated group (spot check: chi(g) from action)
        chis = []
        for h in grp:
            ok, c = fixes(u, h)
            chis.append(c)
        check('%s: every element of the generated group fixes [u] (character values of modulus 1)' % k,
              all(abs(abs(c) - 1) < TOL for c in chis))
        rec['character_values'] = sorted({nstr(mp.chop(c, tol=TOL), 6) for c in chis})
    else:
        rec['order_SO3'] = 'infinite'
        # roots: U1 six-fold at one point; U2 three-fold at two antipodal points
        rec['root_multiplicities'] = [sum(1 for p in R_ if abs(abs(p[0] * mp.conj(q[0]) + p[1] * mp.conj(q[1])) - 1) < TOL)
                                      for q in distinct]
    rec['isotypic_dim_by_J'] = {J: isotypic_dim(J, gens_cf[k]) for J in range(0, 10)}
    out[k] = rec
    print(k, rec)

save_results('stabilizers', out)
