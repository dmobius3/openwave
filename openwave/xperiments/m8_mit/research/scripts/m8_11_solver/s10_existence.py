"""Item 10: collect, for each of U1..U6, the hypotheses of the Lyapunov-Schmidt / equivariant
implicit-function argument of RETURN.md and where each is verified.
  H0  spectral gap: 48 = n(n+2) only at n = 6 among the levels of the sector (item 1);
  H1  Pi_6 N(Phi) = Q Phi (item 2: exact stationarity of rhat6, numerical Pi_6 N - Q Phi = 0);
  H2  nondegeneracy modulo symmetry: ker L = tangent space of the orbit inside W
      (exact rank of the Riemannian Hessian of rhat6 on W, s02; numerical ker L = orbit, s04);
  H3  symmetry: left rotations and phases preserve sector and equation (proof in RETURN.md);
  H4  g != 0 (given)."""
from common import *

d1 = load_results('item1')
ex = load_results('exact_su2')
num = load_results('numeric')
stab = load_results('stabilizers')
out = {}
for name in ('3', '4'):
    levels = [int(n) for n, v in d1['dims'][name].items() if v > 0]
    gaps = {n: n * (n + 2) - 48 for n in levels}
    zero_at = [n for n, gp in gaps.items() if gp == 0]
    mingap = min(abs(gp) for n, gp in gaps.items() if gp != 0)
    # beyond n = 18: n(n+2) - 48 >= 20*22 - 48 = 392 > mingap
    check('sector %s: 48 is an eigenvalue only at level 6 (levels %s)' % (name, levels), zero_at == [6])
    out['gap_sector' + name] = {'levels_le_18': levels, 'min_gap': mingap}
    print('   sector %s: distance from 48 to the rest of the spectrum = %d' % (name, mingap))
for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
    orbit = 3 - stab[key]['lie_algebra_dim']
    rank = ex['hessian_rank_on_W'][key]
    rec = {'stationary_exact': ex[key]['stationary'], 'orbit_dim_in_W': orbit, 'hessian_rank_on_W': rank}
    h1 = ex[key]['stationary']
    h2 = rank + orbit == 12
    for name in ('3', '4'):
        r = num['%s_sector%s' % (key, name)]
        h1 = h1 and r['Pi6N_minus_QPhi_norm'] == '0'
        h2 = h2 and int(r['L_kernel_dim']) == orbit == int(r['orbit_tangent_dim'])
        rec['Q_sector' + name] = r['Q']
    check('%s: H1 critical point (exact) in both sectors' % key, h1)
    check('%s: H2 nondegenerate modulo symmetry: rank %d + orbit dim %d = 12, ker L = orbit tangents' %
          (key, rank, orbit), h2)
    rec['existence'] = bool(h1 and h2)
    out[key] = rec
save_results('item10', out)
