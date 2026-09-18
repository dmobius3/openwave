"""Item 9: every vanishing quantity of items 2, 4-8, classified by the exact reason that makes it
vanish.  Reasons (proved in RETURN.md, item 9):
  CONSTRUCTION  zero by definition (xi has no level-6 part; kappa solves the equation);
  REALITY       <N(Phi), xi> and <N(Phi), kappa> are real / zero (exact identities);
  ISOTYPIC      the rotation stabilizer (item 3) acts on Phi by chi; the quantity lies in the
                chi-isotypic subspace of V_J (x) Hom, and that subspace is 0 (or is C Phi, or C Phi + C E);
  SYM3          Theta u = c u, so rho_6(u) (x) u lies in Sym^3(V_3), which has no spin-J constituent;
  ANTIUNITARY   A' = -rho(R_y(pi)) o complex conjugation fixes Phi, xi and E_1 and sends E_2 -> -E_2;
  ORBIT         the direction is tangent to the SU(2) orbit of Phi at a critical point.
A zero with no applicable reason is reported UNRESOLVED; this script prints each check."""
from common import *
import sympy as sp

num = load_results('numeric')
stab = load_results('stabilizers')
ex = load_results('exact_su2')

# Sym^3(V_3): multiplicity of spin J by weight counting (exact integers)
cnt = {}
for m1 in range(-3, 4):
    for m2 in range(m1, 4):
        for m3 in range(m2, 4):
            cnt[m1 + m2 + m3] = cnt.get(m1 + m2 + m3, 0) + 1
sym3 = {J: cnt.get(J, 0) - cnt.get(J + 1, 0) for J in range(0, 10)}
print('multiplicities of spin J in Sym^3(V_3), J = 0..9:', [sym3[J] for J in range(10)])

# Theta u proportional to u? (exact, sympy)
from sympy import sqrt, Rational as R
Uex = {'U1': {3: 1}, 'U2': {0: 1}, 'U3': {2: 1, -2: 1}, 'U4': {3: 1, -3: 1},
       'U5': {2: sqrt(13) / 5, -3: 2 * sqrt(3) / 5}, 'U6': {3: 1, 0: sqrt(R(23, 10)), -3: 1}}
theta_prop = {}
for k, d in Uex.items():
    u = [d.get(m, 0) for m in range(-3, 4)]
    tu = [(-1) ** (m % 2) * sp.conjugate(u[-m + 3]) for m in range(-3, 4)]
    # proportional iff all 2x2 minors vanish
    prop = all(sp.simplify(u[i] * tu[j] - u[j] * tu[i]) == 0 for i in range(7) for j in range(7))
    theta_prop[k] = prop
print('Theta u proportional to u:', theta_prop)

out = {'zeros': [], 'unresolved': []}


def record(where, qty, reasons):
    ok = len(reasons) > 0
    out['zeros'].append({'where': where, 'quantity': qty, 'reasons': reasons})
    if not ok:
        out['unresolved'].append({'where': where, 'quantity': qty})
    check('%s: %s = 0 (computed) explained by %s' % (where, qty, ', '.join(reasons) if reasons else 'NOTHING'), ok)


for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    iso = {int(J): v for J, v in stab[key]['isotypic_dim_by_J'].items()}
    for name in ('3', '4'):
        r = num['%s_sector%s' % (key, name)]
        w = '%s sector %s' % (key, name)
        # item 2
        if r['Pi6N_minus_QPhi_norm'] == '0':
            reasons = []
            if iso[3] == 1:
                reasons.append('ISOTYPIC (block chi-subspace = C Phi)')
            if ex[key]['stationary']:
                reasons.append('exact stationarity of rhat6 (s02) + Pi6 N = gradient')
            record(w, 'Pi_6 N(Phi) - Q Phi', reasons)
        if r['Q_imag'] == '0':
            record(w, 'Im Q', ['REALITY (Q = int |Phi|^4)'])
        # item 4
        for n_, v in r['normPixi2_over_g2'].items():
            if v != '0':
                continue
            J = int(n_) // 2
            reasons = []
            if J == 3:
                reasons.append('CONSTRUCTION (xi orthogonal to the block)')
            else:
                if iso[J] == 0:
                    reasons.append('ISOTYPIC (chi-subspace of V_%d is 0)' % J)
                if theta_prop[key] and sym3[J] == 0:
                    reasons.append('SYM3 (Theta u = c u; Sym^3 V_3 has no spin %d)' % J)
                if ex[key]['A_J'][str(J)] == '0':
                    reasons.append('exact A_%d = 0 in s02' % J)
            record(w, '||Pi_%s xi||^2' % n_, reasons)
        # item 5
        if r['DNxi_along_Phi_im'] == '0':
            record(w, 'Im <Phi, DN_Phi[xi]>', ['REALITY (= Im 3<N, xi>, and <N, xi> = -g sum ||Pi_n N||^2/mu_n)'])
        if r['DNxi_perp_norm'] == '0':
            record(w, '||P_perp Pi_6 DN_Phi[xi]||', ['ISOTYPIC (block chi-subspace = C Phi)'] if iso[3] == 1 else [])
        for comp in ('DNxi_perp_comp_ie_t', 'DNxi_perp_comp_tau_y', 'kappa_comp_ie_t', 'kappa_comp_tau_y'):
            if r.get(comp) == '0':
                reasons = ['ANTIUNITARY (A\' fixes Phi, xi, E1; A\'E2 = -E2)']
                if comp == 'kappa_comp_ie_t':
                    reasons.append('ORBIT (i e_t in ker L = free direction; kappa chosen orthogonal)')
                record(w, comp, reasons)
        for comp in ('DNxi_perp_remainder_norm', 'kappa_remainder_norm'):
            if r.get(comp) == '0':
                record(w, comp, ['ISOTYPIC (block chi-subspace = C Phi + C E1, dim %d)' % iso[3]] if iso[3] == 2 else [])
        # item 6
        if r.get('form_ie_t') == '0':
            record(w, 'form on i e_t (and rhat6\'\')', ['ORBIT (i e_t is tangent to the R_z orbit; Hessian of an '
                                                        'invariant function vanishes on orbit tangents at a critical point)'])
        if r.get('form_cross') == '0':
            record(w, 'form cross term (tau_x, tau_y)', ['ANTIUNITARY'])
        # item 7
        if r['kappa_norm2_over_g2'] == '0':
            record(w, 'kappa', ['ISOTYPIC (right-hand side P_perp Pi_6 DN[xi] = 0; kappa orthogonal to ker L)']
                   if iso[3] == 1 else [])
        if r['eq_perp_norm_without_kappa'] == '0':
            record(w, 'equation perp component without kappa', ['ISOTYPIC'] if iso[3] == 1 else [])
        # item 8
        if r['lambda4_over_g2_with_kappa'] == r['lambda4_over_g2_without_kappa']:
            record(w, 'lambda4(with kappa) - lambda4(without kappa)',
                   ['REALITY (<Phi, DN_Phi kappa> = 2<N,kappa> + conj<N,kappa>, <N,kappa> = Q<Phi,kappa> = 0)'])
    # item 6 exact SU(2) zeros (sector independent)
for key, keyn in (('U5', 'r6pp_ie_t'), ('U6', 'r6pp_cross')):
    if ex['item6'][key][keyn] == '0':
        record(key, 'rhat6 second derivative ' + keyn, ['ORBIT' if key == 'U5' else 'ANTIUNITARY'])

check('no unresolved zero', len(out['unresolved']) == 0, '%d unresolved' % len(out['unresolved']))
out['sym3_multiplicities'] = sym3
out['theta_proportional'] = theta_prop
save_results('item9', out)
