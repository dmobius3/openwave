"""Items 2-9 by finite algebra: sections as Peter-Weyl coefficient arrays, products by Clebsch-Gordan.
Usage: ./py engine.py DPS      (writes out_engine_<DPS>.json next to this file)

A function f: SU(2) -> C^d is stored as {J: array c[M+J, K+J, a]} meaning f_a(g) = sum c_{MKa} D^J_{MK}(g).
Product rule (textbook, verified numerically below):
  D^{J1}_{M1K1} D^{J2}_{M2K2} = sum_J <J1 M1;J2 M2|J M><J1 K1;J2 K2|J K> D^J_{MK}
Conjugation: conj D^J_{MK} = (-1)^{M-K} D^J_{-M,-K} (verified numerically below).
Norm (Schur orthogonality, Haar mass 1): ||f||^2 = sum_J sum |c|^2/(2J+1)."""
import sys
import json
import random
import numpy as np
import mpmath
from common import (HERE, Q5, generate_gamma, q5num, cg_exact, quat_to_su2, Dj_generic,
                    spin_matrices_generic)

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp = mpmath.mp
mp.dps = DPS
I = mp.mpc(0, 1)
ZERO = mp.mpc(0)
EPS = mp.mpf(10) ** (-(DPS - 12))
out = {"dps": DPS}


def S(x):
    return mpmath.nstr(x, DPS, strip_zeros=False) if not isinstance(x, (int,)) else str(x)


_cg = {}


def cg(j1, m1, j2, m2, J, M):
    key = (j1, m1, j2, m2, J, M)
    if key not in _cg:
        s, a = cg_exact(*key)
        _cg[key] = mp.mpf(0) if s == 0 else mp.mpf(s.numerator) / s.denominator * mp.sqrt(mp.mpf(a.numerator) / a.denominator)
    return _cg[key]


def Dmat(j, U):
    return Dj_generic(2 * j, U, lambda k: mp.factorial(k), mp.sqrt, ZERO)


def zeros(J, d):
    a = np.empty((2 * J + 1, 2 * J + 1, d), dtype=object)
    a[...] = ZERO
    return a


def cg_list(J1, J2, J):
    L = []
    for M in range(-J, J + 1):
        row = []
        for M1 in range(-J1, J1 + 1):
            M2 = M - M1
            if abs(M2) <= J2:
                c = cg(J1, M1, J2, M2, J, M)
                if c != 0:
                    row.append((M1 + J1, M2 + J2, c))
        L.append(row)
    return L


def product(Fd, Gd, Jt, mode):
    """mode 'sv': Fd scalar (d=1) times vector Gd; mode 'dot': sum_a Fd_a Gd_a (scalar result)."""
    dF = next(iter(Fd.values())).shape[2]
    dG = next(iter(Gd.values())).shape[2]
    dout = dG if mode == "sv" else 1
    res = zeros(Jt, dout)
    for J1, A in Fd.items():
        for J2, B in Gd.items():
            if not (abs(J1 - J2) <= Jt <= J1 + J2):
                continue
            L = cg_list(J1, J2, Jt)
            for iM in range(2 * Jt + 1):
                for iK in range(2 * Jt + 1):
                    acc = [ZERO] * dout
                    for (a1, a2, cl) in L[iM]:
                        for (k1, k2, cr) in L[iK]:
                            c = cl * cr
                            if mode == "sv":
                                f = A[a1, k1, 0] * c
                                if f == 0:
                                    continue
                                for a in range(dout):
                                    acc[a] += f * B[a2, k2, a]
                            else:
                                s = ZERO
                                for a in range(dF):
                                    s += A[a1, k1, a] * B[a2, k2, a]
                                acc[0] += c * s
                    for a in range(dout):
                        res[iM, iK, a] += acc[a]
    return res


def conjf(Fd):
    outd = {}
    for J, A in Fd.items():
        B = zeros(J, A.shape[2])
        n = 2 * J
        for i in range(n + 1):
            for k in range(n + 1):
                sgn = 1 if (i - k) % 2 == 0 else -1
                for a in range(A.shape[2]):
                    B[n - i, n - k, a] = sgn * mp.conj(A[i, k, a])
        outd[J] = B
    return outd


def inner(Fd, Gd):
    s = ZERO
    for J in set(Fd) & set(Gd):
        A, B = Fd[J], Gd[J]
        t = ZERO
        for x, y in zip(A.flat, B.flat):
            t += mp.conj(x) * y
        s += t / (2 * J + 1)
    return s


def norm2(Fd):
    return mp.re(inner(Fd, Fd))


def add(*terms):
    outd = {}
    for coef, Fd in terms:
        for J, A in Fd.items():
            if J in outd:
                outd[J] = outd[J] + A * coef
            else:
                outd[J] = A * coef
    return outd


# ---------------------------------------------------------------- verify the product and conjugation rules (checks that can fail)
rng = random.Random(12345)


def rand_su2():
    v = [mp.mpf(rng.gauss(0, 1)) for _ in range(4)]
    n = mp.sqrt(sum(x * x for x in v))
    v = [x / n for x in v]
    return quat_to_su2(v, mp.mpf(1), I)


g0 = rand_su2()
Dg = {j: Dmat(j, g0) for j in range(0, 10)}
dev = 0
for (J1, J2) in [(3, 3), (6, 3), (3, 6), (5, 4)]:
    for trial in range(6):
        M1, K1 = rng.randint(-J1, J1), rng.randint(-J1, J1)
        M2, K2 = rng.randint(-J2, J2), rng.randint(-J2, J2)
        lhs = Dg[J1][M1 + J1][K1 + J1] * Dg[J2][M2 + J2][K2 + J2]
        rhs = ZERO
        for J in range(abs(J1 - J2), J1 + J2 + 1):
            if J > 9 or abs(M1 + M2) > J or abs(K1 + K2) > J:
                continue
            rhs += cg(J1, M1, J2, M2, J, M1 + M2) * cg(J1, K1, J2, K2, J, K1 + K2) * Dg[J][M1 + M2 + J][K1 + K2 + J]
        if J1 + J2 <= 9:
            dev = max(dev, abs(lhs - rhs))
devc = 0
for j in range(0, 10):
    for i in range(2 * j + 1):
        for k in range(2 * j + 1):
            sgn = 1 if (i - k) % 2 == 0 else -1
            devc = max(devc, abs(mp.conj(Dg[j][i][k]) - sgn * Dg[j][2 * j - i][2 * j - k]))
print("check product rule max dev:", mpmath.nstr(dev, 5), "  conj rule max dev:", mpmath.nstr(devc, 5))
assert dev < EPS and devc < EPS
out["check_product_rule_dev"] = mpmath.nstr(dev, 5)
out["check_conj_rule_dev"] = mpmath.nstr(devc, 5)

# ---------------------------------------------------------------- group, D^3, sectors
G = generate_gamma()
D3 = []
for q in G:
    U = quat_to_su2([q5num(c, mp) for c in q], mp.mpf(1), I)
    D3.append(mp.matrix(Dmat(3, U)))
H = mp.matrix(7, 7)
for i in range(7):
    for k in range(7):
        H[i, k] = mp.mpc(mp.mpf(1 + i * 7 + 2 * k * k) / (13 + i + k), mp.mpf(i - k) / (5 + i + 2 * k))
H = (H + H.H) / 2
A = mp.matrix(7, 7)
for D in D3:
    A += D * H * D.H
A /= len(G)
E, Qm = mp.eighe(A)
ev = sorted([(E[i], i) for i in range(7)], key=lambda t: t[0])
clusters = []
for e, i in ev:
    if clusters and abs(e - clusters[-1][0][0]) < mp.mpf(10) ** (-(DPS // 2)):
        clusters[-1].append((e, i))
    else:
        clusters.append([(e, i)])
assert sorted(len(c) for c in clusters) == [3, 4]
ETA = {}
for c in clusters:
    d = len(c)
    eta = np.empty((7, d), dtype=object)
    for a, (e, i) in enumerate(c):
        for r in range(7):
            eta[r, a] = Qm[r, i]
    ETA[d] = eta

# projector cross-check against the character formula P = (d/|G|) sum chi(h) D(h), chi from item 1 (exact)
items01 = json.loads((HERE / "out_item01.json").read_text())


def parse_q5(s):
    s = s.strip("()")
    a, b = s.split(" + ")
    from fractions import Fraction as Fr
    return Q5(Fr(a), Fr(b.replace("*sqrt5", "")))


def re_q(q):
    return q[0]


cls_re = [parse_q5(c["Re_q"]) for c in items01["classes"]]
for d, key in ((3, "character_3dim"), (4, "character_4dim")):
    chi = [parse_q5(x) for x in items01[key]]
    P = mp.matrix(7, 7)
    for q, D in zip(G, D3):
        idx = cls_re.index(re_q(q))
        P += D * q5num(chi[idx], mp)
    P *= mp.mpf(d) / len(G)
    eta = ETA[d]
    Pn = mp.matrix(7, 7)
    for r in range(7):
        for s in range(7):
            Pn[r, s] = sum(eta[r, a] * mp.conj(eta[s, a]) for a in range(d))
    dv = mp.mnorm(P - Pn, 1)
    # eta^dag eta = I
    dI = max(abs(sum(mp.conj(eta[r, a]) * eta[r, b] for r in range(7)) - (1 if a == b else 0)) for a in range(d) for b in range(d))
    print(f"sector {d}: |P(eigen) - P(characters)| = {mpmath.nstr(dv, 5)}, |eta^dag eta - I| = {mpmath.nstr(dI, 5)}")
    assert dv < EPS and dI < EPS


# ---------------------------------------------------------------- vectors in V_j, coupling, Theta
def theta(u):
    """(Theta u)_m = (-1)^m conj(u_{-m}), u indexed by m + j."""
    j = (len(u) - 1) // 2
    return [(1 if (i - j) % 2 == 0 else -1) * mp.conj(u[2 * j - i]) for i in range(2 * j + 1)]


def couple(x, j1, y, j2, J):
    return [sum((cg(j1, m1, j2, M - m1, J, M) * x[m1 + j1] * y[M - m1 + j2]
                 for m1 in range(-j1, j1 + 1) if abs(M - m1) <= j2), ZERO) for M in range(-J, J + 1)]


def rho(u, K):
    return couple(u, 3, theta(u), 3, K)


def Mmap(u, K, J=3):
    return couple(rho(u, K), K, u, 3, J)


def vnorm2(u):
    return mp.re(sum(mp.conj(x) * x for x in u))


def vinner(u, w):
    return sum(mp.conj(x) * y for x, y in zip(u, w))


# ---------------------------------------------------------------- item 2: w(v) = sum_L r_L M_L(v)
out["item2"] = {}
for d, eta in ETA.items():
    cols = [[eta[r, a] for r in range(7)] for a in range(d)]
    T = {}
    rL = {}
    info = {}
    for L in range(0, 7):
        TL = [ZERO] * (2 * L + 1)
        for b in range(d):
            cb = couple(cols[b], 3, theta(cols[b]), 3, L)
            TL = [x + y for x, y in zip(TL, cb)]
        T[L] = TL
        R = [couple(TL, L, cols[a], 3, 3) for a in range(d)]      # R[a][K]
        r = sum(sum(mp.conj(cols[a][k]) * R[a][k] for k in range(7)) for a in range(d)) / d
        resid = mp.sqrt(sum(abs(R[a][k] - r * cols[a][k]) ** 2 for a in range(d) for k in range(7)))
        rL[L] = r
        info[L] = {"T_L_norm": S(mp.sqrt(vnorm2(TL))), "r_L_re": S(mp.re(r)), "r_L_im": S(mp.im(r)),
                   "R_minus_r_eta_resid": S(resid)}
        print(f"sector {d} L={L}: |T_L| = {mpmath.nstr(mp.sqrt(vnorm2(TL)), 8)}, r_L = {mpmath.nstr(r, 20)}, resid {mpmath.nstr(resid, 3)}")
    out["item2"][str(d)] = info
    # check against the full function engine for random unnormalized v (can fail)
    for trial in range(2):
        v = [mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(7)]
        Phi = {3: zeros(3, d)}
        for i in range(7):
            for k in range(7):
                for a in range(d):
                    Phi[3][i, k, a] = v[i] * eta[k, a]
        Sf = {J: product(conjf(Phi), Phi, J, "dot") for J in range(0, 7)}
        N3 = product(Sf, Phi, 3, "sv")
        w = [sum(N3[i, k, a] * mp.conj(eta[k, a]) for k in range(7) for a in range(d)) / d for i in range(7)]
        # reconstruct Pi_6 N from w and compare (checks it is of block form Phi_{sigma,w})
        blk = max(abs(N3[i, k, a] - w[i] * eta[k, a]) for i in range(7) for k in range(7) for a in range(d))
        wf = [sum(rL[L] * Mmap(v, L)[i] for L in range(7)) for i in range(7)]
        dw = max(abs(x - y) for x, y in zip(w, wf))
        w06 = [rL[0] * x + rL[6] * y for x, y in zip(Mmap(v, 0), Mmap(v, 6))]
        dw06 = max(abs(x - y) for x, y in zip(w, w06))
        m0 = Mmap(v, 0)
        dm0 = max(abs(m0[i] + vnorm2(v) * v[i] / mp.sqrt(7)) for i in range(7))
        print(f"  random v #{trial}: block-form resid {mpmath.nstr(blk, 3)}, |w - sum_L r_L M_L| {mpmath.nstr(dw, 3)}, "
              f"|w - (r0 M0 + r6 M6)| {mpmath.nstr(dw06, 3)}, |M0(v) + |v|^2 v/sqrt7| {mpmath.nstr(dm0, 3)}")
        assert blk < EPS and dw < EPS and dw06 < EPS and dm0 < EPS
    out["item2"][str(d)]["r0"] = S(mp.re(rL[0]))
    out["item2"][str(d)]["r6"] = S(mp.re(rL[6]))

# rank of the span of the maps M_0..M_6 (as cubic maps), from samples
samples = [[mp.mpc(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(7)] for _ in range(6)]
rows = []
for K in range(7):
    row = []
    for v in samples:
        mk = Mmap(v, K)
        for x in mk:
            row += [mp.re(x), mp.im(x)]
    rows.append(row)
Mt = mp.matrix(rows)
sv = mp.svd_r(Mt, compute_uv=False)
print("singular values of the 7 maps M_K (stacked on 6 random v):", [mpmath.nstr(s, 6) for s in sv])
out["item2"]["M_K_singular_values"] = [mpmath.nstr(s, 10) for s in sv]
rank = sum(1 for s in sv if s > EPS * 1e6)
out["item2"]["M_K_rank"] = rank
# null space: relations sum_K c_K M_K = 0
U_, sv2, V_ = mp.svd_r(Mt.T)
null = []
for idx in range(rank, 7):
    null.append([V_[idx, K] for K in range(7)])
out["item2"]["M_K_relations_raw"] = [[S(x) for x in vec] for vec in null]

# ---------------------------------------------------------------- rays
c13, s12 = mp.sqrt(13) / 5, 2 * mp.sqrt(3) / 5
RAYS = {}
e = lambda m: [mp.mpc(1) if i - 3 == m else ZERO for i in range(7)]
RAYS["R1"] = e(3)
RAYS["R2"] = e(0)
RAYS["R3"] = [x + y for x, y in zip(e(2), e(-2))]
RAYS["R4"] = [x + y for x, y in zip(e(3), e(-3))]
RAYS["R5"] = [c13 * x + s12 * y for x, y in zip(e(2), e(-3))]

# Casimir from spin matrices, and the Laplacian normalization check (generator of exp(t i sigma_k) = 2 i J_k)
CAS = {}
for J in range(0, 10):
    Jz, Jp, Jm = [mp.matrix(M) for M in spin_matrices_generic(J, mp.sqrt, ZERO, mp.mpc(1))]
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2 * I)
    CAS[J] = (Jx * Jx + Jy * Jy + Jz * Jz, (Jx, Jy, Jz))
t = mp.mpf(10) ** (-(DPS // 3))
pauli_i = [[[0, I], [I, 0]], [[0, 1], [-1, 0]], [[I, 0], [0, -I]]]   # i sigma_x, i sigma_y, i sigma_z
gen_dev = 0
for J in (1, 3, 6, 9):
    Jx, Jy, Jz = CAS[J][1]
    for kidx, Jk in enumerate((Jx, Jy, Jz)):
        X = mp.matrix(pauli_i[kidx])
        Up = mp.expm(X * t)
        Um = mp.expm(-X * t)
        Dp = mp.matrix(Dmat(J, [[Up[0, 0], Up[0, 1]], [Up[1, 0], Up[1, 1]]]))
        Dm_ = mp.matrix(Dmat(J, [[Um[0, 0], Um[0, 1]], [Um[1, 0], Um[1, 1]]]))
        deriv = (Dp - Dm_) / (2 * t)
        gen_dev = max(gen_dev, mp.mnorm(deriv - 2 * I * Jk, 1))
print("check d/dt D^J(exp(t i sigma_k)) = 2 i J_k : max dev", mpmath.nstr(gen_dev, 5))
out["check_generator_dev"] = mpmath.nstr(gen_dev, 5)
assert gen_dev < mp.mpf(10) ** (-(DPS // 2))

mu = {J: 4 * J * (J + 1) - 48 for J in range(0, 10)}
out["rays"] = {}
for d, eta in ETA.items():
    for name, u0 in RAYS.items():
        # normalize the section: ||Phi||^2 = d |u|^2 / 7 = 1
        sc = mp.sqrt(mp.mpf(7) / (d * vnorm2(u0)))
        u = [x * sc for x in u0]
        Phi = {3: zeros(3, d)}
        for i in range(7):
            for k in range(7):
                for a in range(d):
                    Phi[3][i, k, a] = u[i] * eta[k, a]
        nPhi = norm2(Phi)
        cPhi = conjf(Phi)
        Sf = {J: product(cPhi, Phi, J, "dot") for J in range(0, 7)}
        Snorms = {J: norm2({J: Sf[J]}) for J in Sf}
        N = {J: product(Sf, Phi, J, "sv") for J in range(0, 10)}
        Nn = {J: norm2({J: N[J]}) for J in N}
        kappa = inner(Phi, {3: N[3]})
        par = mp.sqrt(norm2(add((1, {3: N[3]}), (-kappa, Phi))))
        # xi / g  (order a^3):  mu_J Pi_J xi = -g Pi_J N,  J != 3
        xi = {J: N[J] * (mp.mpf(-1) / mu[J]) for J in N if J != 3}
        xin = {J: norm2({J: xi[J]}) for J in xi}
        # order a^3 residual with the Casimir: (-Delta - 48) xi - lambda2 Phi + g N, with g = 1, lambda2 = kappa
        def order3_resid(factor):
            r2 = ZERO
            for J in range(0, 10):
                C = CAS[J][0]
                vec = N[J].copy()
                if J == 3:
                    vec = vec - Phi[3] * kappa
                else:
                    X = xi[J]
                    for a in range(d):
                        Xa = mp.matrix([[X[i, k, a] for k in range(2 * J + 1)] for i in range(2 * J + 1)])
                        LX = C * Xa * factor - Xa * 48
                        for i in range(2 * J + 1):
                            for k in range(2 * J + 1):
                                vec[i, k, a] += LX[i, k]
                r2 += norm2({J: vec})
            return mp.sqrt(mp.re(r2))
        resid = order3_resid(4)          # -Delta = 4 (Jx^2 + Jy^2 + Jz^2) on the free index
        resid_control = order3_resid(3)  # deliberately wrong normalization: must NOT vanish
        # order a^5 projected onto the block: Pi_6 DN_Phi[xi] with xi = xi/g
        cxi = conjf(xi)
        t1 = product(Sf, xi, 3, "sv")
        s2 = {J: product(cPhi, xi, J, "dot") for J in range(0, 7)}
        s3 = {J: product(cxi, Phi, J, "dot") for J in range(0, 7)}
        t2 = product(s2, Phi, 3, "sv")
        t3 = product(s3, Phi, 3, "sv")
        X = {3: t1 + t2 + t3}
        along = inner(Phi, X)
        orth = mp.sqrt(norm2(add((1, X), (-along, Phi))))
        # fibre vector y of X = Pi_6 DN_Phi[xi] (block form X = Phi_{sigma,y}), and its part orthogonal to u
        yv = [sum(X[3][i, k, a] * mp.conj(eta[k, a]) for k in range(7) for a in range(d)) / d for i in range(7)]
        blkX = max(abs(X[3][i, k, a] - yv[i] * eta[k, a]) for i in range(7) for k in range(7) for a in range(d))
        uu = vnorm2(u)
        yperp = [y - vinner(u, yv) / uu * x for x, y in zip(u, yv)]
        lam4_formula = -3 * sum(Nn[J] / mu[J] for J in Nn if J != 3)
        # w(v) for the normalized ray and its relation to v
        rec = {
            "norm_Phi": S(nPhi),
            "S_level_norms": {str(J): S(Snorms[J]) for J in Snorms},
            "N_level_norms": {str(J): S(Nn[J]) for J in Nn},
            "kappa_re": S(mp.re(kappa)), "kappa_im": S(mp.im(kappa)),
            "parallel_resid": S(par),
            "xi_level_norm2_over_g2": {str(2 * J): S(xin[J]) for J in xin},
            "order3_casimir_resid": S(resid),
            "order3_casimir_resid_control_factor3": S(resid_control),
            "DN_along_re": S(mp.re(along)), "DN_along_im": S(mp.im(along)),
            "DN_orth_norm": S(orth),
            "DN_orth_norm2": S(orth ** 2),
            "DN_blockform_resid": S(blkX),
            "DN_fibre_y": [[S(mp.re(x)), S(mp.im(x))] for x in yv],
            "DN_fibre_yperp": [[S(mp.re(x)), S(mp.im(x))] for x in yperp],
            "u_normalized": [[S(mp.re(x)), S(mp.im(x))] for x in u],
            "lambda4_over_g2": S(mp.re(along)),
            "lambda4_formula": S(lam4_formula),
        }
        out["rays"][f"{d}:{name}"] = rec
        print(f"[d={d} {name}] |Phi|^2={mpmath.nstr(nPhi, 15)} kappa={mpmath.nstr(kappa, 25)} par_resid={mpmath.nstr(par, 3)}")
        print("    N level norms^2:", {2 * J: mpmath.nstr(Nn[J], 12) for J in Nn})
        print("    S level norms^2:", {2 * J: mpmath.nstr(Snorms[J], 12) for J in Snorms})
        print("    ||Pi_n xi||^2/g^2:", {2 * J: mpmath.nstr(xin[J], 20) for J in xin})
        print(f"    casimir resid={mpmath.nstr(resid, 3)}  DN along={mpmath.nstr(along, 25)}  orth={mpmath.nstr(orth, 3)}  "
              f"-3 sum|Pi N|^2/mu = {mpmath.nstr(lam4_formula, 25)}")

(HERE / f"out_engine_{DPS}.json").write_text(json.dumps(out, indent=1))
print("wrote", f"out_engine_{DPS}.json")
