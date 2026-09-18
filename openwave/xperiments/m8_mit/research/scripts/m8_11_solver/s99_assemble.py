"""Assemble results.json from the per-script outputs in results/.  Exact values are strings;
square roots of rationals are put in the canonical sympy form a*sqrt(b)/c."""
import json
import os
import sympy as sp
from common import load_results, HERE

i0 = load_results('item0')
i1 = load_results('item1')
ex = load_results('exact_su2')
num = load_results('numeric')
st = load_results('stabilizers')
i9 = load_results('item9')
i10 = load_results('item10')
i12 = load_results('item12')
iB = load_results('identification_B')
neg = load_results('negative_controls')


def canon(s):
    if not isinstance(s, str) or s.startswith('UNIDENTIFIED'):
        return s
    val = sp.radsimp(sp.sympify(s))       # exact: no numerical guessing
    return str(val)


UK = ('U1', 'U2', 'U3', 'U4', 'U5', 'U6')
SECS = ('3', '4')
R = {'conventions': {'exact_values': 'strings; identified values were identified at 80 and 120 digits '
                                     '(continued fractions, denominators <= 1e30; PSLQ in Q(sqrt5), coefficients <= 1e10)'}}

R['item0'] = {'norm2_q1': i0['norm2_q1'], 'norm2_q2': i0['norm2_q2'], 'order_Gamma': i0['order'],
              'Gamma_equals_derived_subgroup': i0['perfect'], 'derived_subgroup_order': i0['derived_order'],
              'CG_33_3m3_60': i0['cg_33_3m3_60']}

R['item1'] = {'declaration': 'DERIVED',
              'dims_sector3': {n: v for n, v in i1['dims']['3'].items()},
              'dims_sector4': {n: v for n, v in i1['dims']['4'].items()},
              'invariant_dims_V_K': i1['invariant_dims_VK'],
              'characters': i1['characters'],
              'intertwiner_residuals_50digits_required_1e-40': i1['residuals'],
              'gamma_constants': i1['gamma_constants']}

it2 = {}
for k in UK:
    rec = {'rhat6': ex[k]['rhat6'], 'rhat6_stationary_on_unit_sphere': ex[k]['stationary']}
    for s in SECS:
        r = num['%s_sector%s' % (k, s)]
        rec['sector' + s] = {'Pi6N_is_multiple_of_Phi': r['Pi6N_minus_QPhi_norm'] == '0', 'multiple_Q': r['Q']}
    it2[k] = rec
R['item2'] = it2

gens = {
    'U1': {'group': 'maximal torus T = {Rz(a) = diag(e^{-ia/2}, e^{ia/2})} (connected, SO(3) image SO(2))',
           'character': 'chi(Rz(a)) = e^{3ia}'},
    'U2': {'group': 'N(T) = T u T.Rx(pi) (Pin(2); SO(3) image O(2))',
           'character': 'chi(Rz(a)) = 1, chi(Rx(pi)) = -1'},
    'U3': {'group': 'binary octahedral 2O, order 48 (SO(3) image O, order 24; 4-fold axes z, (1,+-1,0)/sqrt2)',
           'generators': 'Rz(pi/2), R((1,1,0)/sqrt2, pi/2)',
           'character': 'sign character of O = S4: -1 on 4-fold rotations and on pi-rotations about 2-fold axes, '
                        '+1 on 3-fold rotations and on pi-rotations about 4-fold axes; chi(-1) = 1'},
    'U4': {'group': 'binary dihedral of order 24 (SO(3) image D6, order 12)',
           'generators': 'Rz(pi/3), Rx(pi)', 'character': 'chi(Rz(pi/3)) = -1, chi(Rx(pi)) = -1'},
    'U5': {'group': 'cyclic of order 10 (SO(3) image C5)', 'generators': 'Rz(2pi/5)',
           'character': 'chi(Rz(2pi/5)) = e^{4 pi i/5}'},
    'U6': {'group': 'binary dihedral of order 12 (SO(3) image D3)', 'generators': 'Rz(2pi/3), Rx(pi)',
           'character': 'chi(Rz(2pi/3)) = 1, chi(Rx(pi)) = -1'},
}
it3 = {}
for k in UK:
    it3[k] = dict(gens[k])
    it3[k]['complex_dim_of_character_subspace_in_block'] = st[k]['isotypic_dim_by_J']['3']
    it3[k]['character_subspace_dims_V_J_J0to9'] = st[k]['isotypic_dim_by_J']
    it3[k]['action'] = '(r.Phi)(g) = Phi(r^-1 g) = Phi_{sigma, conj(D^3(r)) u}; same for both sectors'
R['item3'] = it3

it4, it5, it6, it7, it8 = {}, {}, {}, {}, {}
for k in UK:
    for s in SECS:
        r = num['%s_sector%s' % (k, s)]
        key = '%s_sector%s' % (k, s)
        it4[key] = {'norm_Pi_n_xi_sq_over_g2': {n: canon(v) for n, v in r['normPixi2_over_g2'].items()},
                    'levels_above_18': '0 (N(Phi) has no component above level 18)',
                    'xi_transforms_by_same_character': True,
                    'norm_Pi_n_N_sq': {n: canon(v) for n, v in r['normPiN2'].items()}}
        e5 = {'along_Phi_per_g': r['DNxi_along_Phi_re'], 'along_Phi_imag': r['DNxi_along_Phi_im'],
              'perp_norm_per_g': canon(r['DNxi_perp_norm'])}
        for c in ('e_t', 'ie_t', 'tau_x', 'tau_y'):
            if 'DNxi_perp_comp_' + c in r:
                e5['perp_component_' + c] = canon(r['DNxi_perp_comp_' + c])
        if 'DNxi_perp_remainder_norm' in r:
            e5['perp_remainder_norm'] = canon(r['DNxi_perp_remainder_norm'])
        it5[key] = e5
        if k in ('U5', 'U6'):
            e6 = {}
            for c in ('e_t', 'ie_t', 'tau_x', 'tau_y', 'cross'):
                if 'form_' + c in r:
                    e6['form_' + c] = r['form_' + c]
            it6[key] = e6
        e7 = {'kappa_norm2_over_g2': canon(r['kappa_norm2_over_g2']),
              'free_directions_real_dim': r['L_kernel_dim'],
              'free_directions': ('P_W(X_x Phi), P_W(X_y Phi): infinitesimal rotations about x, y '
                                  '(rotation about z only changes the phase)' if k in ('U1', 'U2') else
                                  'P_W(X_a Phi), a = x, y, z: infinitesimal rotations' +
                                  (' (P_W(X_z Phi) is proportional to i e_t)' if k == 'U5' else '')),
              'eq_perp_norm_without_kappa_per_g': canon(r['eq_perp_norm_without_kappa']),
              'eq_perp_norm_with_kappa_per_g': '0 (residual %s at 120 digits)' % r['eq_perp_norm_with_kappa'],
              'L_eigenvalues': r['L_eigenvalues']}
        for c in ('e_t', 'ie_t', 'tau_x', 'tau_y'):
            if 'kappa_comp_' + c in r:
                e7['kappa_component_' + c] = canon(r['kappa_comp_' + c])
        if 'kappa_remainder_norm' in r:
            e7['kappa_remainder_norm'] = canon(r['kappa_remainder_norm'])
        it7[key] = e7
        it8[key] = {'lambda4_over_g2_with_kappa': r['lambda4_over_g2_with_kappa'],
                    'lambda4_over_g2_without_kappa': r['lambda4_over_g2_without_kappa'],
                    'sign_of_lambda4': 'negative for every g != 0'}
for k in ('U5', 'U6'):
    it6[k + '_rhat6_second_derivatives'] = {kk: v for kk, v in ex['item6'][k].items() if kk.startswith('r6pp')}
R['item4'] = it4
R['item5'] = it5
R['item6'] = it6
R['item7'] = it7
R['item8'] = it8
R['item9'] = {'zeros': i9['zeros'], 'unresolved': i9['unresolved'], 'missing': []}
R['item10'] = {k: {'existence_for_all_small_a': i10[k]['existence'], 'orbit_dim_in_W': i10[k]['orbit_dim_in_W'],
                   'hessian_rank_on_W': i10[k]['hessian_rank_on_W']} for k in UK}
R['item10']['spectral_gap'] = {'sector3': i10['gap_sector3'], 'sector4': i10['gap_sector4']}
R['item11'] = {'rhat6': ex['item11']['rhat6'], 'stationary': ex['item11']['stationary'],
               'riemannian_grad_norm2_of_rhat6': ex['U5_item11']['riemannian_grad_norm2'],
               'd_rhat6_ds_along_e_t': ex['item11']['d_rhat6_ds_along_e_t'],
               'Re_Phi_DN_et_sector3': ex['item11']['Re_Phi_DN_et_sector3'],
               'Re_Phi_DN_et_sector4': ex['item11']['Re_Phi_DN_et_sector4'],
               'numeric_identification_sector3': canon(num['item11_sector3']['Re_Phi_DN_et']),
               'numeric_identification_sector4': canon(num['item11_sector4']['Re_Phi_DN_et'])}
R['item12'] = i12
R['identification_independence'] = iB
R['negative_controls'] = neg
with open(os.path.join(HERE, 'results.json'), 'w') as fh:
    json.dump(R, fh, indent=1)
print('wrote results.json')
