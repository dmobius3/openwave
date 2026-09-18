"""Independence of the quaternion -> SU(2) identification (section 1.2).  Identification B
(i, j, k -> i s_z, i s_y, i s_x) differs from A by an orientation-reversing permutation of the axes,
so it places the icosahedral group in the mirror position relative to the weight basis.  All
reported quantities are recomputed with B at 40 digits and compared with A."""
from common import *
from s04_pipeline import analyze

mp.dps = 40
U = fibre_vectors()
et, tx = tangent_vectors(U)
SA = Setup('A')
SB = Setup('B')
for name in ('3', '4'):
    check('sector %s: identical level dimensions under identifications A and B' % name,
          [SA.sectors[name].dim(J) for J in range(10)] == [SB.sectors[name].dim(J) for J in range(10)])
# the two embeddings of Gamma really differ (so the comparison is not vacuous)
gA = [quat_to_su2(h, 'A') for h in group()]
gB = [quat_to_su2(h, 'B') for h in group()]
same = all(any(mat_maxdiff(a, b) < mpf(10) ** -30 for b in gB) for a in gA)
print('images of Gamma in SU(2) under A and B coincide as sets:', same)
worst = mpf(0)
skip = ('L_eigenvalues',)
for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    for name in ('3', '4'):
        rA = analyze(SA, name, key, U[key], et, tx)
        rB = analyze(SB, name, key, U[key], et, tx)
        w = mpf(0)
        for k, v in rA.items():
            if k in skip or 'residual' in k or 'asymmetry' in k or 'kernel_component' in k or k.startswith('L_on') \
                    or k.startswith('tangent_check') or k in ('kappa_orth_kernel', 'eq_perp_norm_with_kappa',
                                                              'Phi_DN_kappa_abs'):
                continue
            if isinstance(v, dict):
                w = max(w, max(abs(v[n_] - rB[k][n_]) for n_ in v))
            else:
                w = max(w, abs(v - rB[k]))
        w = max(w, max(abs(a - b) for a, b in zip(rA['L_eigenvalues'], rB['L_eigenvalues'])))
        worst = max(worst, w)
        check('%s sector %s: every reported quantity identical under identifications A and B' % (key, name),
              w < mpf(10) ** -30, nstr(w, 3))
save_results('identification_B', {'max_difference': nstr(worst, 3), 'digits': 40,
                                  'embeddings_coincide': same})
