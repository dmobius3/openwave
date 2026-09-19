"""Items 2, 4, 5, 6, 7, 8, 11: the full Gamma-dependent computation in the section algebra of
common.py (Clebsch-Gordan products of sections, intertwiners derived from the generators),
run from scratch at 80 and at 120 digits.  Every reported number is identified (identify():
0, Q, +-sqrt(Q), Q(sqrt5), +-sqrt(Q(sqrt5)); see common.identify) at both precisions and the two
identifications must agree.  Values also obtained exactly in s02 are compared with these."""
from common import *
import sympy
import time

UKEYS = ('U1', 'U2', 'U3', 'U4', 'U5', 'U6')


def mu(J):
    return 4 * J * (J + 1) - 48


def W_basis(sec, Phi):
    """Orthonormal real basis of W = {kappa in block : <Phi, kappa> = 0}."""
    def P(h):
        return sec_add(h, sec_scale(Phi, sec.inner(Phi, h)), -1)
    Wb = []
    for m in range(7):
        for c in (mpc(1), mpc(0, 1)):
            v = [mpc(0)] * 7
            v[m] = c
            h = P({3: [v]})
            for w in Wb:
                h = sec_add(h, w, -mp.re(sec.inner(w, h)))
            nrm = mp.sqrt(sec.norm2(h))
            if nrm > mpf(10) ** (-(mp.dps // 2)):
                Wb.append(sec_scale(h, 1 / nrm))
    return Wb


def left_act(sec_vecs, g):
    """Left translation of a section given as {J: [fibre vectors]}."""
    return {J: [mat_vec(rho_left(J, g), u) for u in us] for J, us in sec_vecs.items()}


def antiunitary(sec_vecs):
    """A' = -rho(R_y(pi)) o Theta on every level (conj of a section; eta's are Theta-real)."""
    g = rot((0, 1, 0), mp.pi)
    return {J: [vscale(mat_vec(rho_left(J, g), theta(u)), -1) for u in us] for J, us in sec_vecs.items()}


def sec_maxdiff(A, B):
    m = mpf(0)
    for J in set(A) | set(B):
        for u, v in zip(A.get(J, []), B.get(J, [])):
            m = max(m, max(abs(x - y) for x, y in zip(u, v)))
    return m


def analyze(S, name, key, u, et, tx):
    sec = S.sectors[name]
    r = {}
    Phi = sec.block(u)
    up = Phi[3][0]
    f = sec.scal(Phi, Phi)
    N = sec.mul(f, Phi, range(0, 10))
    Qc = sec.inner(Phi, N)
    Q = mp.re(Qc)
    r['Q'] = Q
    r['Q_imag'] = mp.im(Qc)
    r['Pi6N_minus_QPhi_norm'] = mp.sqrt(sec.level_norm2(sec_add({3: N[3]}, sec_scale(Phi, Qc), -1), 3))
    levels = [J for J in range(10) if sec.dim(J) > 0]
    r['normPiN2'] = {2 * J: sec.level_norm2(N, J) for J in levels}
    xi = {J: [vscale(v, -1 / mpf(mu(J))) for v in N[J]] for J in levels if J != 3}
    r['normPixi2_over_g2'] = {2 * J: (sec.level_norm2(xi, J) if J != 3 else mpf(0)) for J in levels}
    # xi transforms like Phi under the stabilizer (item 4)
    worst = mpf(0)
    for g, chi in stabilizer_generators()[key]:
        worst = max(worst, sec_maxdiff(left_act(xi, g), sec_scale(xi, chi)))
        worst = max(worst, sec_maxdiff(left_act(Phi, g), sec_scale(Phi, chi)))
    r['xi_covariance_residual'] = worst
    # item 5
    DNxi = sec.DN(Phi, xi)
    along = sec.inner(Phi, DNxi)
    r['DNxi_along_Phi_re'] = mp.re(along)
    r['DNxi_along_Phi_im'] = mp.im(along)
    perp = sec_add(DNxi, sec_scale(Phi, along), -1)
    r['DNxi_perp_norm'] = mp.sqrt(sec.norm2(perp))
    tangents = None
    if key == 'U5':
        E1 = sec.block(et)
        tangents = (('e_t', E1), ('ie_t', sec_scale(E1, mpc(0, 1))))
    if key == 'U6':
        E1 = sec.block(tx)
        tangents = (('tau_x', E1), ('tau_y', sec_scale(E1, mpc(0, 1))))
    if tangents:
        rem = perp
        for tn, E in tangents:
            c = mp.re(sec.inner(E, perp))
            r['DNxi_perp_comp_' + tn] = c
            rem = sec_add(rem, E, -c)
        r['DNxi_perp_remainder_norm'] = mp.sqrt(sec.norm2(rem))
        # item 6
        for tn, E in tangents:
            r['form_' + tn] = mp.re(sec.inner(E, sec.DN(Phi, E))) - Q
            r['tangent_check_' + tn] = abs(sec.inner(Phi, E)) + abs(sec.norm2(E) - 1)
        if key == 'U6':
            r['form_cross'] = mp.re(sec.inner(tangents[0][1], sec.DN(Phi, tangents[1][1])))
        # antiunitary symmetry A' (item 9): fixes Phi, xi, first tangent; flips the second
        r['antiunitary_residual'] = max(sec_maxdiff(antiunitary(Phi), Phi), sec_maxdiff(antiunitary(xi), xi),
                                        sec_maxdiff(antiunitary(tangents[0][1]), tangents[0][1]),
                                        sec_maxdiff(antiunitary(tangents[1][1]), sec_scale(tangents[1][1], -1)))
    # item 7: L on W
    Wb = W_basis(sec, Phi)
    n = len(Wb)
    DNW = [sec.DN(Phi, w) for w in Wb]
    L = mp.matrix(n)
    for j in range(n):
        for i in range(n):
            L[i, j] = mp.re(sec.inner(Wb[i], DNW[j])) - (Q if i == j else 0)
    sym = max(abs(L[i, j] - L[j, i]) for i in range(n) for j in range(n))
    r['L_asymmetry'] = sym
    ev, V = mp.eigsy(L)
    evs = sorted([ev[i] for i in range(n)])
    tolz = mpf(10) ** (-(mp.dps // 2))
    zero_idx = [a for a in range(n) if abs(ev[a]) < tolz]
    r['L_eigenvalues'] = evs
    r['L_kernel_dim'] = len(zero_idx)
    r['L_min_abs_nonzero_eigenvalue'] = min(abs(ev[a]) for a in range(n) if a not in zero_idx)
    # orbit tangents
    Jx, Jy, Jz = spin_matrices(3)
    T = []
    for Ja in (Jx, Jy, Jz):
        t = {3: [vscale(mat_vec(mat_conj(Ja), up), mpc(0, 1))]}
        t = sec_add(t, sec_scale(Phi, sec.inner(Phi, t)), -1)
        T.append([mp.re(sec.inner(w, t)) for w in Wb])
    Tm = mp.matrix(n, 3)
    for a in range(3):
        for i in range(n):
            Tm[i, a] = T[a][i]
    svT = mp.svd_r(Tm, compute_uv=False)
    orbit_dim = sum(1 for a in range(3) if svT[a] > tolz)
    r['orbit_tangent_dim'] = orbit_dim
    LT = L * Tm
    r['L_on_orbit_tangents'] = max(abs(LT[i, a]) for i in range(n) for a in range(3))
    b = [-mp.re(sec.inner(w, perp)) for w in Wb]
    r['rhs_kernel_component'] = mp.sqrt(mp.fsum(mp.fsum(V[i, a] * b[i] for i in range(n)) ** 2 for a in zero_idx))
    kap = [mp.fsum(V[i, a] * mp.fsum(V[k, a] * b[k] for k in range(n)) / ev[a]
                   for a in range(n) if a not in zero_idx) for i in range(n)]
    kappa = {3: [[mpc(0)] * 7]}
    for c, w in zip(kap, Wb):
        kappa = sec_add(kappa, w, c)
    r['kappa_norm2_over_g2'] = sec.norm2(kappa)
    r['kappa_orth_kernel'] = max([abs(mp.fsum(V[i, a] * kap[i] for i in range(n))) for a in zero_idx] + [mpf(0)])
    if tangents:
        rem = kappa
        for tn, E in tangents:
            c = mp.re(sec.inner(E, kappa))
            r['kappa_comp_' + tn] = c
            rem = sec_add(rem, E, -c)
        r['kappa_remainder_norm'] = mp.sqrt(sec.norm2(rem))
        if key == 'U5':
            # i e_t is an orbit direction (free): its component in ker L
            t = [mp.re(sec.inner(w, tangents[1][1])) for w in Wb]
            Lt = [mp.fsum(L[i, j] * t[j] for j in range(n)) for i in range(n)]
            r['L_on_ie_t'] = max(abs(x) for x in Lt)
    eqk = sec_add(sec_add(sec.DN(Phi, kappa), kappa, -Q), DNxi)
    eqk = sec_add(eqk, sec_scale(Phi, sec.inner(Phi, eqk)), -1)
    r['eq_perp_norm_with_kappa'] = mp.sqrt(sec.norm2(eqk))
    r['eq_perp_norm_without_kappa'] = r['DNxi_perp_norm']
    # item 8
    lk = sec.inner(Phi, sec.DN(Phi, kappa))
    r['lambda4_over_g2_without_kappa'] = mp.re(along)
    r['lambda4_over_g2_with_kappa'] = mp.re(along) + mp.re(lk)
    r['Phi_DN_kappa_abs'] = abs(lk)
    r['lambda4_formula'] = -3 * mp.fsum(r['normPiN2'][2 * J] / mu(J) for J in levels if J != 3)
    return r


def item11(S, name):
    sec = S.sectors[name]
    c, s = mp.sqrt(3) / 2, mpf(1) / 2
    u = vadd(vscale(basis_vec(3, 2), c), vscale(basis_vec(3, -3), s))
    e = vadd(vscale(basis_vec(3, 2), -s), vscale(basis_vec(3, -3), c))
    Phi = sec.block(u)
    E = sec.block(e)
    N = sec.N(Phi, (3,))
    Qc = sec.inner(Phi, N)
    return {'Re_Phi_DN_et': mp.re(sec.inner(Phi, sec.DN(Phi, E))),
            'Pi6N_minus_QPhi_norm': mp.sqrt(sec.level_norm2(sec_add({3: N[3]}, sec_scale(Phi, Qc), -1), 3)),
            'Pi6N_along_e_t_re': mp.re(sec.inner(E, N))}


def run(dps, ident='A'):
    with mp.workdps(dps):
        t0 = time.time()
        S = Setup(ident)
        U = fibre_vectors()
        et, tx = tangent_vectors(U)
        res = {}
        for name in ('3', '4'):
            sec = S.sectors[name]
            th = max(sec_maxdiff({J: [theta(col(e, a)) for a in range(sec.d)]},
                                 {J: [col(e, a) for a in range(sec.d)]})
                     for J in sec.eta for e in sec.eta[J])
            res['theta_real_eta_' + name] = th
            for key in UKEYS:
                res[(key, name)] = analyze(S, name, key, U[key], et, tx)
            res[('item11', name)] = item11(S, name)
            res['expand_residual_' + name] = sec.max_expand_residual
            res['dropped_r_' + name] = sec.max_dropped_r
        print('  run at %d digits (identification %s): %.1fs' % (dps, ident, time.time() - t0))
        return res


# quantities that are identified exactly (the rest are residuals / diagnostics)
EXACT_KEYS = ['Q', 'Q_imag', 'Pi6N_minus_QPhi_norm', 'DNxi_along_Phi_re', 'DNxi_along_Phi_im', 'DNxi_perp_norm',
              'DNxi_perp_comp_e_t', 'DNxi_perp_comp_ie_t', 'DNxi_perp_comp_tau_x', 'DNxi_perp_comp_tau_y',
              'DNxi_perp_remainder_norm', 'form_e_t', 'form_ie_t', 'form_tau_x', 'form_tau_y', 'form_cross',
              'kappa_norm2_over_g2', 'kappa_comp_e_t', 'kappa_comp_ie_t', 'kappa_comp_tau_x', 'kappa_comp_tau_y',
              'kappa_remainder_norm', 'eq_perp_norm_without_kappa', 'lambda4_over_g2_without_kappa',
              'lambda4_over_g2_with_kappa', 'L_kernel_dim', 'orbit_tangent_dim']
RESIDUAL_KEYS = ['xi_covariance_residual', 'antiunitary_residual', 'L_asymmetry', 'L_on_orbit_tangents',
                 'rhs_kernel_component', 'kappa_orth_kernel', 'eq_perp_norm_with_kappa', 'Phi_DN_kappa_abs',
                 'L_on_ie_t', 'tangent_check_e_t', 'tangent_check_ie_t', 'tangent_check_tau_x', 'tangent_check_tau_y']

if __name__ == '__main__':
    LO, HI = PRECISIONS
    runs = {p: run(p) for p in PRECISIONS}
    exact = load_results('exact_su2')
    out = {}
    allok = True

    def ident2(label, v50, v80):
        global allok
        if isinstance(v50, int):
            ok = v50 == v80
            if not ok:
                allok = False
            return v50, ok
        with mp.workdps(LO):
            s50, _ = identify(v50)
        with mp.workdps(HI):
            s80, _ = identify(v80)
        ok = s50 is not None and s50 == s80
        if not ok:
            allok = False
            print('   identification mismatch/failure for %s: %s | %s  (%s)' % (label, s50, s80, nstr(v80, 30)))
        return (s50 if ok else 'UNIDENTIFIED ' + nstr(v80, 40)), ok

    def exact_cmp(label, expr_str, v80):
        with mp.workdps(HI + 10):
            ev = mpf(str(sympy.N(sympy.sympify(expr_str), HI + 15)))
            d = abs(ev - v80)
        return check('%s agrees with the exact s02 value %s' % (label, expr_str), d < mpf(10) ** -70, nstr(d, 3))

    for name in ('3', '4'):
        for tag in ('theta_real_eta_', 'expand_residual_', 'dropped_r_'):
            check('sector %s: %s (%d digits) below 1e-70' % (name, tag.strip('_'), HI),
                  runs[HI][tag + name] < mpf(10) ** -70, nstr(runs[HI][tag + name], 3))
    for key in UKEYS + ('item11',):
        for name in ('3', '4'):
            r50, r80 = runs[LO][(key, name)], runs[HI][(key, name)]
            rec = {}
            for k in r80:
                if k in ('normPiN2', 'normPixi2_over_g2'):
                    rec[k] = {}
                    for n_ in r80[k]:
                        rec[k][n_] = ident2('%s %s %s n=%d' % (key, name, k, n_), r50[k][n_], r80[k][n_])[0]
                elif k == 'L_eigenvalues':
                    # reported numerically only (not asked for exactly); agreement of the two runs checked
                    rec[k] = [nstr(b, 20) for b in r80[k]]
                    check('%s s%s eigenvalues of L agree at both precisions' % (key, name),
                          max(abs(a - b) for a, b in zip(r50[k], r80[k])) < mpf(10) ** -60)
                elif k in RESIDUAL_KEYS or k in ('L_min_abs_nonzero_eigenvalue', 'lambda4_formula'):
                    rec[k] = nstr(r80[k], 5)
                else:
                    rec[k] = ident2('%s %s %s' % (key, name, k), r50[k], r80[k])[0]
            out['%s_sector%s' % (key, name)] = rec
            print('%s sector %s:' % (key, name))
            for k, v in rec.items():
                print('    %-34s %s' % (k, v))
            # comparisons with the exact SU(2) x Gamma-constant values of s02
            if key in UKEYS:
                cmb = exact['combined']['%s_sector%s' % (key, name)]
                exact_cmp('%s s%s Q' % (key, name), cmb['Q'], r80['Q'])
                exact_cmp('%s s%s lambda4/g^2 (with kappa)' % (key, name), cmb['lambda4_over_g2'],
                          r80['lambda4_over_g2_with_kappa'])
                for n_, v in cmb['normPixi2_over_g2'].items():
                    exact_cmp('%s s%s ||Pi_%s xi||^2/g^2' % (key, name, n_), v, r80['normPixi2_over_g2'][int(n_)])
                check('%s s%s lambda4 equals -3 sum ||Pi_n N||^2/(n(n+2)-48)' % (key, name),
                      abs(r80['lambda4_formula'] - r80['lambda4_over_g2_with_kappa']) < mpf(10) ** -70)
                check('%s s%s xi (and Phi) transform by the item-3 character' % (key, name),
                      r80['xi_covariance_residual'] < mpf(10) ** -70, nstr(r80['xi_covariance_residual'], 3))
                check('%s s%s ker L = orbit tangents (dims %d = %d, L T = 0)' %
                      (key, name, r80['L_kernel_dim'], r80['orbit_tangent_dim']),
                      r80['L_kernel_dim'] == r80['orbit_tangent_dim'] and r80['L_on_orbit_tangents'] < mpf(10) ** -70)
                check('%s s%s order-a^5 block equation solvable (rhs has no ker L component)' % (key, name),
                      r80['rhs_kernel_component'] < mpf(10) ** -70, nstr(r80['rhs_kernel_component'], 3))
                check('%s s%s with kappa the equation has no component orthogonal to Phi' % (key, name),
                      r80['eq_perp_norm_with_kappa'] < mpf(10) ** -70, nstr(r80['eq_perp_norm_with_kappa'], 3))
                if key in ('U5', 'U6'):
                    it6 = exact['item6'][key]
                    pairs = ((('e_t', 'e_t'), ('ie_t', 'ie_t')) if key == 'U5' else
                             (('tau_x', 'tau_x'), ('tau_y', 'itau_x'), ('cross', 'cross')))
                    for mine, theirs in pairs:
                        kk = 'form_%s_sector%s' % (theirs, name) if theirs != 'cross' else 'form_cross_sector%s' % name
                        exact_cmp('%s s%s form on %s' % (key, name, mine), it6[kk], r80['form_' + mine])
                    check('%s s%s antiunitary A\' fixes Phi, xi, first tangent and flips the second' % (key, name),
                          r80['antiunitary_residual'] < mpf(10) ** -70, nstr(r80['antiunitary_residual'], 3))
            else:
                exact_cmp('item11 s%s Re<Phi, DN e_t>' % name, exact['item11']['Re_Phi_DN_et_sector' + name],
                          r80['Re_Phi_DN_et'])
    check('all exact values identified identically at both precisions', allok)
    save_results('numeric', out)
