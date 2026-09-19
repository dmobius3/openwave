"""Item 12: check (-Delta - 48) xi + g (N(Phi) - Pi_6 N(Phi)) = 0 as a vector identity whose
two sides come from different operations:
  * left side: -Delta xi obtained by applying 4 * Casimir, built from explicit spin matrices
    J_x, J_y, J_z (Condon-Shortley J_+), to the fibre vectors of the xi of item 4 (s04 algebra);
  * right side: N(Phi) = |Phi|^2 Phi evaluated pointwise on a product grid of Euler angles
    g = R_z(alpha) R_y(beta) R_z(gamma) (Phi(g) = u^T D^3(g) eta from the 2x2 matrix g), and
    projected onto every level J <= 11 by quadrature, keeping the full coefficient tensor
    F^J[M,K,a] = (2J+1) int conj(D^J_MK) N_a  (all copies, and any part outside Hom).
The quadrature (24 x 14 x 24 points: trapezoid in alpha, gamma; Gauss-Legendre in cos beta) is
exact for integrands of total spin <= 22, i.e. for projections of N (spin <= 9) up to J = 13.
-Delta = 4 * Casimir: on the unit S^3 = SU(2) the right-invariant fields Y_q (q = i, j, k unit
quaternions) are orthonormal and Delta = sum Y_q^2; under the identification q -> -i sigma_a,
Y_q acts on u^T D(g) eta by u^T -> u^T (-2i J_a), so -Delta acts on u by 4 sum (J_a^T)^2.
Working precision 30 digits; required residual 1e-25."""
from common import *
from s04_pipeline import mu
import time

mp.dps = 30
REQ = mpf(10) ** -25
NA = 24
NB = 14
JMAX = 11

def make_grid(na, nb):
    """trapezoid grid in alpha, gamma (na points) and Gauss-Legendre in cos(beta) (nb nodes)."""
    Pn = mp.taylor(lambda x: mp.legendre(nb, x), 0, nb)
    xs = sorted(mp.re(r) for r in mp.polyroots(Pn[::-1], maxsteps=500, extraprec=300))

    def dP(x):
        return nb * (x * mp.legendre(nb, x) - mp.legendre(nb - 1, x)) / (x * x - 1)

    ws = [2 / ((1 - x * x) * dP(x) ** 2) for x in xs]
    betas = [mp.acos(x) for x in xs]
    alphas = [2 * mp.pi * k / na for k in range(na)]
    dsmall = {(i, J): wigner_D(2 * J, rot((0, 1, 0), betas[i])) for i in range(nb) for J in range(JMAX + 1)}
    return {'na': na, 'nb': nb, 'xs': xs, 'ws': ws, 'betas': betas, 'alphas': alphas, 'dsmall': dsmall}


Jmats = {J: spin_matrices(J) for J in range(JMAX + 1)}


def minus_laplacian(u, J):
    out = [mpc(0)] * len(u)
    for Ja in Jmats[J]:
        JT = [[Ja[k][i] for k in range(len(u))] for i in range(len(u))]
        out = vadd(out, mat_vec(JT, mat_vec(JT, u)))
    return vscale(out, 4)


def project_quadrature(sec, up, grid):
    NA, NB = grid['na'], grid['nb']
    alphas, betas, ws, dsmall = grid['alphas'], grid['betas'], grid['ws'], grid['dsmall']
    d = sec.d
    eta = sec.eta[3][0]
    # pointwise N on the grid, then Fourier in alpha, gamma for each beta
    F = {J: [[[mpc(0)] * d for _ in range(2 * J + 1)] for _ in range(2 * J + 1)] for J in range(JMAX + 1)}
    ex = {(M, k): mp.expj(M * alphas[k]) for M in range(-JMAX, JMAX + 1) for k in range(NA)}
    for i in range(NB):
        Nv = {}
        for ka in range(NA):
            for kg in range(NA):
                g = mat_mul(mat_mul(rot((0, 0, 1), alphas[ka]), rot((0, 1, 0), betas[i])), rot((0, 0, 1), alphas[kg]))
                D = wigner_D(6, g)
                row = [mp.fsum(up[m] * D[m][k] for m in range(7)) for k in range(7)]
                phi = [mp.fsum(row[k] * eta[k][a] for k in range(7)) for a in range(d)]
                n2 = mp.fsum(abs(x) ** 2 for x in phi)
                Nv[(ka, kg)] = [n2 * x for x in phi]
        # G[M][K][a] = mean over alpha, gamma of N e^{iM alpha} e^{iK gamma}
        H = {}
        for ka in range(NA):
            for K in range(-JMAX, JMAX + 1):
                H[(ka, K)] = [mp.fsum(Nv[(ka, kg)][a] * ex[(K, kg)] for kg in range(NA)) / NA for a in range(d)]
        G = {}
        for M in range(-JMAX, JMAX + 1):
            for K in range(-JMAX, JMAX + 1):
                G[(M, K)] = [mp.fsum(H[(ka, K)][a] * ex[(M, ka)] for ka in range(NA)) / NA for a in range(d)]
        for J in range(JMAX + 1):
            dj = dsmall[(i, J)]
            for M in range(-J, J + 1):
                for K in range(-J, J + 1):
                    c = (2 * J + 1) * ws[i] / 2 * dj[M + J][K + J]
                    for a in range(d):
                        F[J][M + J][K + J][a] += c * G[(M, K)][a]
    return F


def tensor_of(sec, vecs, J):
    """sum_gamma u^gamma (x) eta^gamma_J as a (2J+1) x (2J+1) x d tensor."""
    d = sec.d
    T = [[[mpc(0)] * d for _ in range(2 * J + 1)] for _ in range(2 * J + 1)]
    for u, eta in zip(vecs, sec.eta.get(J, [])):
        for M in range(2 * J + 1):
            for K in range(2 * J + 1):
                for a in range(d):
                    T[M][K][a] += u[M] * eta[K][a]
    return T


def tmax(T):
    return max(abs(x) for P in T for r in P for x in r)


def tsub(A, B):
    return [[[x - y for x, y in zip(ra, rb)] for ra, rb in zip(PA, PB)] for PA, PB in zip(A, B)]


if __name__ == '__main__':
    t0 = time.time()
    S = Setup('A')
    U = fibre_vectors()
    grid = make_grid(NA, NB)
    xs, ws = grid['xs'], grid['ws']
    check('Gauss-Legendre weights sum to 2 and integrate x^26 exactly',
          abs(mp.fsum(ws) - 2) < REQ and abs(mp.fsum(w * x ** 26 for w, x in zip(ws, xs)) - mpf(2) / 27) < REQ)
    # Euler factorisation used for the projection (checked, can fail)
    g_e = mat_mul(mat_mul(rot((0, 0, 1), grid['alphas'][5]), rot((0, 1, 0), grid['betas'][3])), rot((0, 0, 1), grid['alphas'][7]))
    De = wigner_D(14, g_e)
    err = max(abs(De[M + 7][K + 7] - mp.expj(-M * grid['alphas'][5]) * grid['dsmall'][(3, 7)][M + 7][K + 7] * mp.expj(-K * grid['alphas'][7]))
              for M in range(-7, 8) for K in range(-7, 8))
    check('D^J(Rz(a)Ry(b)Rz(c)) = e^{-iMa} d^J_MK(b) e^{-iKc}', err < REQ, nstr(err, 3))

    out = {}
    worst_all = mpf(0)
    for name in ('3', '4'):
        sec = S.sectors[name]
        for key in ('U1', 'U2', 'U3', 'U4', 'U5', 'U6'):
            Phi = sec.block(U[key])
            N = sec.N(Phi)                                         # algebraic (for comparison only)
            xi = {J: [vscale(v, -1 / mpf(mu(J))) for v in N[J]] for J in N if J != 3}
            F = project_quadrature(sec, Phi[3][0], grid)
            res_eq = mpf(0)
            res_alg = mpf(0)
            res_out = mpf(0)
            above = mpf(0)
            for J in range(JMAX + 1):
                FJ = F[J]
                if J >= 10:
                    above = max(above, tmax(FJ))
                    continue
                # left side: (-Delta - 48) xi via the Casimir
                lhs_vecs = [vadd(minus_laplacian(u, J), u, -48) for u in xi.get(J, [])]
                lhs = tensor_of(sec, lhs_vecs, J) if lhs_vecs else [[[mpc(0)] * sec.d for _ in range(2 * J + 1)]
                                                                   for _ in range(2 * J + 1)]
                rhs = FJ if J != 3 else [[[mpc(0)] * sec.d for _ in range(7)] for _ in range(7)]  # N - Pi_6 N
                res_eq = max(res_eq, tmax([[[x + y for x, y in zip(ra, rb)] for ra, rb in zip(PA, PB)]
                                           for PA, PB in zip(lhs, rhs)]))
                # quadrature N versus the algebraic N (every copy), and the part of F outside Hom
                res_alg = max(res_alg, tmax(tsub(FJ, tensor_of(sec, N.get(J, []), J))))
                proj = []
                for eta in sec.eta.get(J, []):
                    proj.append([mp.fsum(FJ[M][K][a] * mp.conj(eta[K][a]) for K in range(2 * J + 1)
                                         for a in range(sec.d)) / sec.d for M in range(2 * J + 1)])
                res_out = max(res_out, tmax(tsub(FJ, tensor_of(sec, proj, J))) if proj else tmax(FJ))
            worst_all = max(worst_all, res_eq)
            out['%s_sector%s' % (key, name)] = {'residual_equation': nstr(res_eq, 3),
                                                'quadrature_vs_algebra': nstr(res_alg, 3),
                                                'component_outside_Hom': nstr(res_out, 3),
                                                'max_component_levels_20_22': nstr(above, 3)}
            check('%s sector %s: |(-Delta-48)xi + N - Pi6 N| over all levels <= 18' % (key, name), res_eq < REQ,
                  nstr(res_eq, 3))
            check('%s sector %s: quadrature N equals algebraic N at every level (all copies)' % (key, name),
                  res_alg < REQ, nstr(res_alg, 3))
            check('%s sector %s: quadrature N has no part outside V_J (x) Hom(sigma, V_J)' % (key, name),
                  res_out < REQ, nstr(res_out, 3))
            check('%s sector %s: quadrature N has no component at levels 20, 22' % (key, name), above < REQ,
                  nstr(above, 3))
    print('elapsed %.0fs' % (time.time() - t0))
    out['required'] = '1e-25 at 30 digits'
    out['worst_residual'] = nstr(worst_all, 3)
    save_results('item12', out)
