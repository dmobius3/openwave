"""M8.11 fresh derivation, step 3: the author's independent high-precision route (algebraic identification).

The values are identified by limit_denominator(10^25) with a residual below 10^-(dps - 12): very strong candidate exact
values, not symbolic evaluations (F1). The file name is historical.

No code is shared with m811_ops.py (no product_level, no DN). The route uses the paper's density structure:
  * for block sections the mixed density Phi_a^dagger Phi_b is c0 + d12, with c0 = (d/7)<a, b> and d12 a single
    level-12 matrix coefficient X_d R^6(g) y_d(a, b); that the density has levels 0 and 12 only (Lemma 4.1,
    M8.1.2 C1) is checked here, not assumed;
  * a level-12 density times a block section has, at level 2J, the factorized form X_J R^J(g) Y_J, with
    X_J = (X_d (x) P) C_J^T fixed by the sector and Y_J = C_J (y_d (x) c) on the left index;
  * so every quantity reduces to left-index Clebsch-Gordan algebra times the constants C_J = tr(X_J^+ X_J)/(2J+1).
Shared, and disclosed: m810_exact's mpmath 2I group, D^3, isotypic projectors and exact Racah coefficients (the
pinned M8.10 primitives). Everything else is new.

Field prediction, stated before the run: a rotation reverses each slice coordinate (R_z(pi/5) at the pyramid,
R_z(pi/3) at the prism), so the forcing and the tilt are odd in it. They are a rational times sqrt(39) at the
pyramid (sin t cos t = 2 sqrt(39)/25) and a rational times sqrt(230) at the prism (x0 = sqrt(230)/10). The prism
half failed on the first run and is corrected, with the reason, in the prism block: the field there is sqrt(115). Q, the
second variation, the level norms, lambda_4 and ||v||^2 are rational. Precision M811_DPS (default 80). g = 1.
"""
import json
import os
import sys
from fractions import Fraction

os.environ['M810_DPS'] = os.environ.get('M811_DPS', '80')
import mpmath as mp
import m810_exact as X

DPS = mp.mp.dps
TOL = mp.mpf(10) ** -(DPS - 20)
TOL_ID, MAXDEN = mp.mpf(10) ** -(DPS - 12), 10 ** 25
LEVELS = {'3p': {6, 10, 14, 16, 18}, '4': {6, 8, 12, 14, 16, 18}}
bad, count, OUT = [], [0], {}


def check(name, ok, detail=''):
    count[0] += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''))
    if not ok:
        bad.append(name)


def ms(j):
    return list(range(j, -j - 1, -1))


def Wmat(j):
    W, m = mp.matrix(2 * j + 1, 2 * j + 1), ms(j)
    for i, mi in enumerate(m):
        W[i, m.index(-mi)] = (-1) ** (j - mi)
    return W


def coupling(j1, j2, J):
    C, m2s = mp.matrix(2 * J + 1, (2 * j1 + 1) * (2 * j2 + 1)), ms(j2)
    for a, M in enumerate(ms(J)):
        for i1, m1 in enumerate(ms(j1)):
            if abs(M - m1) <= j2:
                C[a, i1 * (2 * j2 + 1) + m2s.index(M - m1)] = X.cg(j1, m1, j2, M - m1, J, M)
    return C


def conj(M):
    R = mp.matrix(M.rows, M.cols)
    for i in range(M.rows):
        for j in range(M.cols):
            R[i, j] = mp.conj(M[i, j])
    return R


def dag(M):
    return conj(M).T


def kron(A, B):
    R = mp.matrix(A.rows * B.rows, A.cols * B.cols)
    for i in range(A.rows):
        for j in range(A.cols):
            if A[i, j] != 0:
                for k in range(B.rows):
                    for l in range(B.cols):
                        R[i * B.rows + k, j * B.cols + l] = A[i, j] * B[k, l]
    return R


def inner(a, b):                               # Hermitian, conjugate-linear in a
    return mp.fsum(mp.conj(a[i]) * b[i] for i in range(a.rows))


def maxabs(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def spin(J):
    n, m = 2 * J + 1, ms(J)
    jp = mp.matrix(n, n)
    for c in range(1, n):
        jp[c - 1, c] = mp.sqrt(J * (J + 1) - m[c] * (m[c] + 1))
    return (jp + jp.T) / 2, (jp - jp.T) * mp.mpc(0, -0.5), mp.diag([mp.mpf(x) for x in m])


SPIN = {J: spin(J) for J in range(10)}


def DJ(J, q):
    s = mp.sqrt(q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    Jx, Jy, Jz = SPIN[J]
    return mp.expm(mp.mpc(0, -1) * 2 * mp.atan2(s, q[0]) * (q[1] / s * Jx + q[2] / s * Jy + q[3] / s * Jz))


def e(m):
    v = mp.matrix(7, 1)
    v[3 - m] = 1
    return v


def den(n):
    return n * (n + 2) - 48


def ident(x, scale=1):
    y = mp.re(x) / scale
    f = Fraction(mp.nstr(y, DPS - 5)).limit_denominator(MAXDEN)
    return f, abs(y - mp.mpf(f.numerator) / f.denominator)


W3, W6 = Wmat(3), Wmat(6)
C33 = {K: coupling(3, 3, K) for K in range(7)}
C63 = {J: coupling(6, 3, J) for J in range(3, 10)}
IW = kron(mp.eye(7), W3.T)


class Sector:
    def __init__(self, s):
        self.s, self.d, self.P = s, {'3p': 3, '4': 4}[s], X.PROJ[s]
        Pp = self.P - mp.mpf(self.d) / 7 * mp.eye(7)
        b = IW * mp.matrix([Pp[k, l] for k in range(7) for l in range(7)])
        self.q = {K: C33[K] * b for K in range(7)}
        self.xd = W6.T * self.q[6]                              # X_d = xd^T
        self.L12 = mp.re(inner(self.xd, self.xd)) / 13
        self.XJ = {J: kron(self.xd.T, self.P) * C63[J].T for J in C63}
        self.CJ = {J: mp.re(sum(inner(self.XJ[J][:, c], self.XJ[J][:, c]) for c in range(2 * J + 1))) / (2 * J + 1) for J in C63}
        self.sup = [J for J in C63 if J != 3 and 2 * J in LEVELS[s]]
        self.k = mp.sqrt(mp.mpf(7) / self.d)

    def yd(self, a, b, twist=True):
        return W6.T * (C33[6] * ((IW if twist else mp.eye(49)) * kron(conj(a), b)))

    def Y(self, J, a, b, c):
        return C63[J] * kron(self.yd(a, b), c)

    def c0(self, a, b):
        return mp.mpf(self.d) / 7 * inner(a, b)

    def dd(self, x, y, z, w):                                  # int d(x, y) d(z, w)
        return self.c0(x, y) * self.c0(z, w) + self.L12 * inner(self.yd(y, x), self.yd(z, w))

    def Q(self, u):
        return mp.re(self.dd(u, u, u, u))

    def grad(self, u, a):                                      # <Phi_a, N(Phi_u)>
        return self.dd(u, u, a, u)

    def H(self, u, a, b):                                      # Re<a, DN_Phi[b]> - Q Re c0(a, b)
        return mp.re(self.dd(u, b, a, u) + self.dd(u, u, a, b) + self.dd(b, u, a, u)) - self.Q(u) * mp.re(self.c0(a, b))

    def lemma(self, u, a):                                     # <Phi, DN_Phi[a]>
        return self.dd(u, a, u, u) + self.dd(u, u, u, a) + self.dd(a, u, u, u)

    def forcing(self, u, t):                                   # <t, DN_Phi[xi]>, xi = -A^-1 Pi_perp N(Phi)
        Yu = {J: self.Y(J, u, u, u) for J in self.sup}
        s1 = sum(-self.CJ[J] / den(2 * J) * inner(self.Y(J, u, t, u), Yu[J]) for J in self.sup)
        s2 = sum(-self.CJ[J] / den(2 * J) * inner(self.Y(J, u, u, t), Yu[J]) for J in self.sup)
        s3 = sum(-self.CJ[J] / den(2 * J) * inner(self.Y(J, t, u, u), Yu[J]) for J in self.sup)
        return s1 + s2 + mp.conj(s3)

    def levels(self, u):
        return {2 * J: self.CJ[J] * mp.re(inner(self.Y(J, u, u, u), self.Y(J, u, u, u))) for J in self.sup}

    def lam4(self, u):
        return sum(-3 * v / den(n) for n, v in self.levels(u).items())

    def tangential_grad(self, u):
        return max(abs(self.grad(u, self.k * e(m) - self.c0(u, self.k * e(m)) / self.c0(u, u) * u)) for m in ms(3))


def pyramid(S, sgn=1, s2=mp.mpf(12) / 25):
    c, s = mp.sqrt(1 - s2), sgn * mp.sqrt(s2)
    return S.k * (c * e(2) + s * e(-3)), S.k * (-s * e(2) + c * e(-3))


def prism(S, x):
    return S.k * (e(3) + x * e(0) + e(-3)) / mp.sqrt(2 + x * x), S.k * (2 * e(0) - x * (e(3) + e(-3))) / mp.sqrt(4 + 2 * x * x)


X0 = mp.sqrt(mp.mpf(23) / 10)
NOTE = {'3p': {5: 1.318681e-1, 7: 8.144796e-2, 8: 5.042017e-2, 9: 5.715647e-2},
        '4': {4: 1.318681e-1, 6: 1.318681e-1, 7: 5.042017e-2, 8: 8.144796e-2, 9: 7.471167e-2}}
PAR = {'3p': dict(Qpy=Fraction(77, 65), Hpy=Fraction(-56, 165), Qpr=Fraction(5831, 5031), Hpr=(Fraction(6440, 18447), Fraction(56, 429)),
                  f39=Fraction(7188839, 81061695000), beta=Fraction(1026977, 3930264000), lam=Fraction(-267786421, 58544557500),
                  xi={10: Fraction(77, 5475600), 14: Fraction(281211, 189112352000), 16: Fraction(553, 1591200000), 18: Fraction(30233, 56458242360)}),
       '4': dict(Qpy=Fraction(287, 260), Hpy=Fraction(-21, 110), Qpr=Fraction(609, 559), Hpr=(Fraction(2415, 12298), Fraction(21, 286)),
                 f39=Fraction(-19565553, 192146240000), beta=Fraction(-931693, 1746784000), lam=Fraction(-4797453339, 2497901120000),
                 xi={8: Fraction(3591, 346112000), 12: Fraction(63, 83200000), 14: Fraction(361557, 931014656000), 16: Fraction(34839, 147097600000),
                     18: Fraction(2371131, 8029616691200)})}
FLOAT = {n: json.load(open(f'out/{n}.json')) for n in ('pyramid', 'prism')}
feq = lambda x, f: abs(x - mp.mpf(f.numerator) / f.denominator) < TOL
rel = lambda a, b: abs(a / b - 1) if b != 0 else abs(a)
rng = mp.mpf(0)
qs = [[mp.mpf(x) for x in v] for v in ((0.3, 0.5, -0.7, 0.2), (-0.6, 0.1, 0.4, 0.55))]
qs = [[x / mp.sqrt(sum(y * y for y in v)) for x in v] for v in qs]
S0 = {}

for s in ('3p', '4'):
    S, P = Sector(s), X.PROJ[s]
    S0[s] = S
    print(f'\n=== sector {s} (d = {S.d}) ===')
    o = OUT[s] = {}

    # The density structure, checked.
    check('Lemma 4.1 / C1: the density has no multipole K = 1..5', max(maxabs(S.q[K]) for K in range(1, 6)) < TOL)
    a, b = S.k * (e(2) + 2 * e(-1)) / mp.sqrt(5), S.k * (e(3) - e(0)) / mp.sqrt(2)
    p0 = C33[0] * (IW * kron(conj(a), b))
    check('the K = 0 part is c0 = (d/7)<a, b>', abs((p0.T * S.q[0])[0, 0] - S.c0(a, b)) < TOL)
    check("D^3 here equals m810_exact's D^3", max(maxabs(DJ(3, q) - X.D3(q)) for q in qs) < TOL)
    u, _ = pyramid(S)
    worst_d, worst_n, worst_bad = mp.mpf(0), mp.mpf(0), mp.mpf(0)
    for q in qs:
        R = {J: dag(DJ(J, q)) for J in range(10)}
        Phi = P * R[3] * u
        dens = inner(Phi, Phi)
        worst_d = max(worst_d, abs(dens - S.c0(u, u) - (S.xd.T * R[6] * S.yd(u, u))[0, 0]))
        worst_bad = max(worst_bad, abs(dens - S.c0(u, u) - (S.xd.T * R[6] * S.yd(u, u, twist=False))[0, 0]))
        rhs = S.c0(u, u) * Phi
        for J in C63:
            rhs += S.XJ[J] * R[J] * S.Y(J, u, u, u)
        worst_n = max(worst_n, maxabs(dens * Phi - rhs))
    check('pointwise at two group elements: |Phi|^2 = c0 + X_d R^6 y_d', worst_d < TOL, mp.nstr(worst_d, 3))
    check('pointwise: N(Phi) = c0 Phi + sum_J X_J R^J Y_J, the factorized cubic', worst_n < TOL, mp.nstr(worst_n, 3))
    check('arm: without the W twist in y_d the density identity fails', worst_bad > mp.mpf('1e-3'), mp.nstr(worst_bad, 3))
    kap = sum(S.XJ[3][i, i] for i in range(7)) / S.d
    check('the level-6 part is a block section: X_3 = kappa P', maxabs(S.XJ[3] - kap * P) < TOL)

    # The sector constants C_J: forbidden levels vanish, the rest are the post-derivation note's C_{sigma,n}.
    forb = [J for J in C63 if J != 3 and 2 * J not in LEVELS[s]]
    check(f'forbidden levels {[2 * J for J in forb]}: C_J = 0 exactly', all(S.CJ[J] < TOL for J in forb),
          ', '.join(mp.nstr(S.CJ[J], 3) for J in forb))
    CJx = {}
    for J in S.sup:
        f, r = ident(S.CJ[J])
        CJx[2 * J] = f
        check(f'C at level {2 * J}: rational {f}, matches the note (7 digits)', r < TOL_ID and rel(S.CJ[J], NOTE[s][J]) < 1e-6,
              f'note {NOTE[s][J]}')
    o['C'] = {n: str(f) for n, f in CJx.items()}

    # The pyramid.
    u, t = pyramid(S)
    Q, H = S.Q(u), S.H(u, t, t)
    check('pyramid: Q = parent, exactly', feq(Q, PAR[s]['Qpy']), str(PAR[s]['Qpy']))
    check('pyramid: critical in the whole block', S.tangential_grad(u) < TOL)
    check('pyramid: second variation along e_t = parent, exactly', feq(H, PAR[s]['Hpy']), str(PAR[s]['Hpy']))
    check('pyramid: the direction i e_t is null (rotation orbit)', abs(S.H(u, 1j * t, 1j * t)) < TOL)
    F = S.forcing(u, t)
    f39, r39 = ident(F, mp.sqrt(39))
    check('pyramid: forcing has no orbit (imaginary) component', abs(mp.im(F)) < TOL)
    check('pyramid: forcing = rational x sqrt(39), as predicted, and equals the audited value', r39 < TOL_ID and f39 == PAR[s]['f39'], f'{f39} sqrt(39)')
    check('arm: the forcing is not a small-denominator rational x sqrt(230)', ident(F, mp.sqrt(230))[1] > TOL_ID)
    Fm = S.forcing(*pyramid(S, -1))
    check('pyramid parity: F(-t) = -F(t)', abs(Fm + F) < TOL)
    beta = -mp.re(F) / H
    b39, rb = ident(beta, mp.sqrt(39))
    check('pyramid: tilt beta = rational x sqrt(39), equal to the inherited form', rb < TOL_ID and b39 == PAR[s]['beta'], f'{b39} sqrt(39)')
    lv = S.levels(u)
    xi = {n: ident(v / den(n) ** 2) for n, v in lv.items()}
    check("pyramid: xi level norms rational and equal to #547's blind diagnostics", all(r < TOL_ID for _, r in xi.values())
          and {n: f for n, (f, _) in xi.items()} == PAR[s]['xi'])
    lam, rl = ident(S.lam4(u))
    check("pyramid: lambda_4 rational and equal to #547's blind diagnostic", rl < TOL_ID and lam == PAR[s]['lam'], str(lam))
    check('pyramid: lambda_4 lemma, <Phi, DN_Phi[e_t]> = 0 at the critical ray', abs(S.lemma(u, t)) < TOL)
    uc, tc = pyramid(S, 1, mp.mpf(1) / 4)
    check('arm: at sin^2 t = 1/4 the gradient and the lemma pairing are nonzero', S.tangential_grad(uc) > mp.mpf('1e-3') and abs(S.lemma(uc, tc)) > mp.mpf('1e-3'))
    fl = FLOAT['pyramid'][s]
    check('pyramid: the float pipeline agrees (tilt, lambda_4, forcing)', rel(beta, fl['tilt_beta']) < 1e-12 and rel(S.lam4(u), fl['lambda4']) < 1e-12
          and rel(mp.re(F), fl['forcing_along_et_section']) < 1e-12)
    o['pyramid'] = dict(Q=str(PAR[s]['Qpy']), H=str(PAR[s]['Hpy']), forcing=f'{f39}*sqrt(39)', tilt=f'{b39}*sqrt(39)', tilt_norm2=str(b39 ** 2 * 39),
                        lambda4=str(lam), xi={n: str(f) for n, (f, _) in xi.items()})

    # The prism.
    u, tx = prism(S, X0)
    ty = 1j * tx
    Q, Hx, Hy, Hc = S.Q(u), S.H(u, tx, tx), S.H(u, ty, ty), S.H(u, tx, ty)
    check('prism: Q = parent, exactly', feq(Q, PAR[s]['Qpr']), str(PAR[s]['Qpr']))
    check('prism: critical in the whole block', S.tangential_grad(u) < TOL)
    check('prism: second variation (tau_x, tau_y) = parent, exactly', feq(Hx, PAR[s]['Hpr'][0]) and feq(Hy, PAR[s]['Hpr'][1]))
    check('prism: no cross term (S-forced: the form is S-even, tau_x and tau_y have opposite S-parity)', abs(Hc) < TOL)
    F = S.forcing(u, tx)
    check('prism: the S-odd forcing component vanishes (Im F along tau_x)', abs(mp.im(F)) < TOL, mp.nstr(mp.im(F), 3))
    # The field. First pass, recorded: the stated prediction, rational x sqrt(230), FAILED in both sectors (the best
    # fits had 25-digit denominators). The parity argument held (checked just below); the error was the unit tangent's
    # norm sqrt(4 + 2 x0^2) = sqrt(2) sqrt(2 + x0^2), which I had dropped. With it, F = (x0/sqrt 2) x (a rational
    # function of x0^2), and x0/sqrt(2) = sqrt(115)/10, so the field is rational x sqrt(115).
    f115, r115 = ident(F, mp.sqrt(115))
    check('prism: forcing = rational x sqrt(115), the corrected prediction', r115 < TOL_ID, f'{f115} sqrt(115)')
    check('arm: the first-pass field, rational x sqrt(230), fails', ident(F, mp.sqrt(230))[1] > TOL_ID)
    check('arm: the forcing is not a small-denominator rational x sqrt(39)', ident(F, mp.sqrt(39))[1] > TOL_ID)
    check('prism parity: F(-x0) = -F(x0)', abs(S.forcing(*prism(S, -X0)) + F) < TOL)
    vx = -mp.re(F) / Hx
    v230, rv = ident(vx, mp.sqrt(115))
    check('prism: tilt v_x = rational x sqrt(115); v_y = 0', rv < TOL_ID, f'{v230} sqrt(115)')
    f230 = f115
    lv = S.levels(u)
    xi = {n: ident(v / den(n) ** 2) for n, v in lv.items()}
    lam, rl = ident(S.lam4(u))
    check('prism: xi level norms and lambda_4 rational', all(r < TOL_ID for _, r in xi.values()) and rl < TOL_ID, str(lam))
    check('prism: lambda_4 lemma, <Phi, DN_Phi[tau_x]> = 0', abs(S.lemma(u, tx)) < TOL)
    fl = FLOAT['prism'][s]
    check('prism: the float pipeline agrees (tilt, lambda_4, forcing, level norms)', rel(vx, fl['tilt'][0]) < 1e-12 and rel(S.lam4(u), fl['lambda4']) < 1e-12
          and rel(mp.re(F), fl['forcing'][0]) < 1e-12 and all(rel(lv[n] / den(n) ** 2, fl['xi_level_norms'][str(n)]) < 1e-11 for n in lv))
    o['prism'] = dict(Q=str(PAR[s]['Qpr']), H=[str(h) for h in PAR[s]['Hpr']], forcing=f'{f230}*sqrt(115)', tilt=[f'{v230}*sqrt(115)', '0'],
                      tilt_norm2=str(v230 ** 2 * 115), lambda4=str(lam), xi={n: str(f) for n, (f, _) in xi.items()})
    print(f'  prism exact: forcing {f230} sqrt(115), tilt v_x {v230} sqrt(115), ||v||^2 {v230 ** 2 * 115}, lambda_4 {lam}')
    print('  prism xi level norms: ' + ', '.join(f'{n}: {f}' for n, (f, _) in xi.items()))

# Arm across sectors, and the prism's lambda_4 sector ratio (F2's flag).
u3 = pyramid(S0['3p'])[0]
swapped = {n: S0['4'].CJ.get(n // 2, 0) * mp.re(inner(S0['3p'].Y(n // 2, u3, u3, u3), S0['3p'].Y(n // 2, u3, u3, u3))) for n in (10, 14, 16, 18)}
check("arm: sector 4's constants on sector 3' fail #547's pyramid values", any(ident(v / den(n) ** 2)[0] != PAR['3p']['xi'][n] for n, v in swapped.items()))
l3, l4 = Fraction(OUT['3p']['prism']['lambda4']), Fraction(OUT['4']['prism']['lambda4'])
ratio = l3 / l4
print(f"\n  prism lambda_4 ratio (3'/4) = {ratio} = {float(ratio):.15f}; 1 - ratio = {1 - ratio} = {float(1 - ratio):.6e}")
OUT['prism_lambda4_ratio'] = str(ratio)
os.makedirs('out', exist_ok=True)
json.dump(OUT, open('out/exact.json', 'w'), indent=1)
print(f'\n  {count[0]} checks at {DPS} digits; ' + ('ALL PASS, ALL ARMS FIRE' if not bad else f'FAILURES: {bad}'))
sys.exit(1 if bad else 0)
