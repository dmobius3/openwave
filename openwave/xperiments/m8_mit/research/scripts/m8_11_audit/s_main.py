"""Items 2, 4, 5, 6, 7, 8, 11 exactly (MQ arithmetic), for U1..U6 (+ U5' of item 11) and both sectors."""
import pickle, json, sys, time
from fractions import Fraction
from mq import MQ, ZERO, ONE, mq
from su2 import couple, theta
from engine import block, N, DN, compress, block_fibre, inner, norm2, fadd, fscale, vdot, vadd, vscale
from linalg_mq import nullspace, rank, solve, matvec, transpose
from checks import check, save_log

F = Fraction


def vec(**kw):
    v = [ZERO] * 7
    for k, c in kw.items():
        m = int(k[1:]) * (-1 if k[0] == "m" else 1)
        v[m + 3] = mq(c)
    return v

U = {
    "U1": vec(p3=1),
    "U2": vec(p0=1),
    "U3": vec(p2=1, m2=1),
    "U4": vec(p3=1, m3=1),
    "U5": vec(p2=MQ.sqrt(F(13, 25)), m3=MQ.sqrt(F(12, 25))),
    "U6": vec(p3=1, p0=MQ.sqrt(F(23, 10)), m3=1),
    "U5q": vec(p2=MQ.sqrt(F(3, 4)), m3=MQ.sqrt(F(1, 4))),   # item 11: sin^2 t = 1/4
}
# tangent fibre directions (unnormalized) and orientation
TANG = {
    "U5": {"e_t": vec(p2=-MQ.sqrt(F(12, 25)), m3=MQ.sqrt(F(13, 25)))},
    "U5q": {"e_t": vec(p2=-MQ.sqrt(F(1, 4)), m3=MQ.sqrt(F(3, 4)))},
}
# U6: component of v0 orthogonal to U6
u6 = U["U6"]
v0 = vec(p0=1)
tx = vadd(v0, vscale(u6, -(vdot(u6, v0) / vdot(u6, u6))))
TANG["U6"] = {"tau_x": tx}


def sqrt_mq_rational(x):
    assert list(x.t.keys()) == [1] and x.t[1][1] == 0
    return MQ.sqrt(x.t[1][0])


def real_coords(v):
    return [a.re() for a in v] + [a.im() for a in v]


def from_real(c):
    return [c[k] + MQ.I * c[k + 7] for k in range(7)]


def S(x):
    return x.pretty() if isinstance(x, MQ) else str(x)


def fl(x):
    return float(x.to_mpc().real) if x.is_real() else complex(x.to_mpc())


def mu(J):
    n = 2 * J
    return n * (n + 2) - 48


def rho(x, y, K=6):
    return couple(x, 3, theta(y, 3), 3, K)


def analyze(name, d):
    t0 = time.time()
    R = {}
    u0 = U[name]
    nu2 = vdot(u0, u0)
    uhat = vscale(u0, sqrt_mq_rational(nu2).inv())            # unit fibre vector
    s7d = MQ.sqrt(F(7, d))
    u = vscale(uhat, s7d)                                       # section-normalized fibre vector
    Phi = block(u, d)
    R["norm_Phi_sq"] = S(norm2(Phi))
    assert norm2(Phi) == ONE
    # ---------------- item 2
    Nf = compress(N(Phi))
    Fu = block_fibre(Nf, d)
    Q = vdot(u, Fu) / vdot(u, u)
    perpF = vadd(Fu, vscale(u, -Q))
    is_mult = all(c.is_zero() for c in perpF)
    Qint = inner(Phi, Nf)          # int |Phi|^4
    R["Pi6N_multiple_of_Phi"] = is_mult
    R["Q_fibre_ratio"] = S(Q)
    R["int_Phi4"] = S(Qint)
    R["Pi6N_perp_norm_sq"] = S(vdot(perpF, perpF) * F(d, 7))
    # r6hat and stationarity
    r0 = rho(uhat, uhat)
    r6 = vdot(r0, r0)
    R["r6hat"] = S(r6)
    grads = []
    for k in range(14):
        w = [ZERO] * 7
        w[k % 7] = ONE if k < 7 else MQ.I
        B = vadd(rho(uhat, w), rho(w, uhat))
        g = vdot(r0, B).re() * 2 - r6 * vdot(uhat, w).re() * 4
        grads.append(g)
    R["r6hat_stationary"] = all(g.is_zero() for g in grads)
    R["r6hat_grad_nonzero_components"] = [S(g) for g in grads if not g.is_zero()]
    # ---------------- item 4
    levels = sorted(J for J in Nf if J != 3)
    xi = {J: [(vscale(x, F(-1, mu(J))), Y) for x, Y in Nf[J]] for J in levels}
    R["levels_N"] = {2 * J: S(norm2({J: Nf[J]})) for J in sorted(Nf)}
    R["xi_level_norm_sq_over_g2"] = {2 * J: S(norm2({J: xi[J]})) for J in levels}
    # ---------------- item 5
    Gf = block_fibre(compress(DN(Phi, xi, levels={3})), d)      # fibre of Pi_6 DN_Phi[xi] / g
    along = vdot(u, Gf) * F(d, 7)                                # <Phi, DN xi>
    Gperp = vadd(Gf, vscale(u, -(vdot(u, Gf) / vdot(u, u))))
    R["DNxi_along_Phi"] = S(along)
    R["DNxi_perp_norm_sq"] = S(vdot(Gperp, Gperp) * F(d, 7))
    lam4_identity = sum((norm2({J: Nf[J]}) * F(-3, mu(J)) for J in levels), ZERO)
    R["lambda4_identity_-3sum"] = S(lam4_identity)
    # tangents (normalized sections)
    tang = {}
    if name in TANG:
        for tn, w0 in TANG[name].items():
            what = vscale(w0, sqrt_mq_rational(vdot(w0, w0)).inv())
            assert vdot(uhat, what).is_zero()
            tang[tn] = vscale(what, s7d)
        if name.startswith("U5"):
            tang["i_e_t"] = vscale(tang["e_t"], MQ.I)
        else:
            tang["tau_y"] = vscale(tang["tau_x"], MQ.I)
        R["DNxi_tangent_components"] = {tn: S(vdot(e, Gf).re() * F(d, 7)) for tn, e in tang.items()}
    # ---------------- item 6: operator L~ on R^14
    def PT(v):
        return vadd(v, vscale(u, -(vdot(u, v) / vdot(u, u))))

    def DF(v):
        return block_fibre(compress(DN(Phi, block(v, d), levels={3})), d)

    cols = []
    for k in range(14):
        e = [ZERO] * 7
        e[k % 7] = ONE if k < 7 else MQ.I
        pe = PT(e)
        cols.append(real_coords(vadd(PT(DF(pe)), vscale(pe, -Q))))
    Lm = transpose(cols)   # Lm[i][k]
    sym = all((Lm[i][k] - Lm[k][i]).is_zero() for i in range(14) for k in range(14))
    R["L_symmetric"] = sym

    def qform(e1, e2):
        # <e1, (DN - Q) e2>_R for block sections with fibres e1, e2
        return (vdot(e1, DF(e2)).re() - Q * vdot(e1, e2).re()) * F(d, 7)
    if tang:
        tn = list(tang)
        R["qform"] = {f"{a}": S(qform(tang[a], tang[a])) for a in tn}
        R["qform_cross"] = S(qform(tang[tn[0]], tang[tn[1]]))
        R["qform_cross_sym"] = S(qform(tang[tn[1]], tang[tn[0]]))
        # r6hat second derivatives along unit tangents
        def r6dd(wh):
            B = vadd(rho(uhat, wh), rho(wh, uhat))
            return (vdot(B, B) + vdot(r0, rho(wh, wh)).re() * 2 - r6 * 2) * 2
        units = {a: vscale(tang[a], s7d.inv()) for a in tn}
        R["r6hat_dd"] = {a: S(r6dd(units[a])) for a in tn}
        a, b = units[tn[0]], units[tn[1]]
        h = MQ.sqrt(F(1, 2))
        R["r6hat_dd_cross"] = S((r6dd(vscale(vadd(a, b), h)) - r6dd(vscale(vadd(a, vscale(b, -1)), h))) * F(1, 2))
    # ---------------- kernel and item 7
    ker = nullspace(Lm)
    R["dim_ker_Ltilde_R"] = len(ker)
    # rotation directions i conj(J_a) u, and phase i u
    from su2 import cg  # noqa
    Jp = [[ZERO] * 7 for _ in range(7)]
    for m in range(-3, 3):
        Jp[m + 1 + 3][m + 3] = MQ.sqrt((3 - m) * (3 + m + 1))
    Jm = transpose(Jp)
    Jx = [[(a + b) * F(1, 2) for a, b in zip(r1, r2)] for r1, r2 in zip(Jp, Jm)]
    Jy = [[(a - b) * F(1, 2) * (-MQ.I) for a, b in zip(r1, r2)] for r1, r2 in zip(Jp, Jm)]
    Jz = [[mq(m) if i == k else ZERO for k in range(7)] for i, m in enumerate(range(-3, 4))]
    conjM = lambda A: [[x.conj() for x in r] for r in A]
    rot = [real_coords(PT(vscale(matvec(conjM(Jn), u), MQ.I))) for Jn in (Jx, Jy, Jz)]
    rot_rank = rank(rot) if any(not c.is_zero() for r in rot for c in r) else 0
    R["dim_rotation_orbit_tangent"] = rot_rank
    rot_in_ker = all(all(c.is_zero() for c in matvec(Lm, r)) for r in rot)
    R["rotation_dirs_in_ker"] = rot_in_ker
    # kernel of L on the tangent space T (real dim 12) = ker L~ minus the 2 dims (u, iu)
    R["dim_ker_L_on_T"] = len(ker) - 2
    ker_plus_rot = rank(ker + rot) if rot_rank else len(ker)
    R["ker_equals_span(u,iu,rotations)"] = (len(ker) == 2 + rot_rank) and ker_plus_rot == len(ker)
    # solve L~ kappa = rhs, kappa perp ker (Euclidean in real coords == real inner product up to d/7)
    rhs = [-c for c in real_coords(Gperp)]
    nk = len(ker)
    A = [Lm[i] + [ker[j][i] for j in range(nk)] for i in range(14)] + [ker[j] + [ZERO] * nk for j in range(nk)]
    sol = solve(A, rhs + [ZERO] * nk)
    kc, muk = sol[:14], sol[14:]
    kappa = from_real(kc)
    resid_ker = [sum((muk[j] * ker[j][i] for j in range(nk)), ZERO) for i in range(14)]
    resid_vec = from_real(resid_ker)
    R["perp_eq_norm_sq_without_kappa"] = S(vdot(Gperp, Gperp) * F(d, 7))
    R["perp_eq_norm_sq_with_kappa"] = S(vdot(resid_vec, resid_vec) * F(d, 7))
    R["rhs_in_range"] = all(c.is_zero() for c in resid_ker)
    R["kappa_norm_sq_over_g2"] = S(vdot(kappa, kappa) * F(d, 7))
    assert vdot(u, kappa).is_zero()
    if tang:
        comps = {}
        rem = kappa
        for tn, e in tang.items():
            c = vdot(e, kappa).re() * F(d, 7)
            comps[tn] = S(c)
            rem = vadd(rem, vscale(e, -c))
        R["kappa_tangent_components"] = comps
        R["kappa_remainder_norm_sq"] = S(vdot(rem, rem) * F(d, 7))
    # verify the kappa equation exactly: PT(DF kappa) - Q kappa + Gperp = resid (kernel part)
    chk = vadd(vadd(PT(DF(kappa)), vscale(kappa, -Q)), Gperp)
    R["kappa_eq_residual_equals_ker_component"] = all((a - b).is_zero() for a, b in zip(chk, resid_vec))
    # ---------------- item 8
    Kf = block(kappa, d)
    DNk = block_fibre(compress(DN(Phi, Kf, levels={3})), d)
    lam4_with = (vdot(u, vadd(Gf, DNk)) - Q * vdot(u, kappa)) * F(d, 7)
    R["lambda4_over_g2_without_kappa"] = S(along)
    R["lambda4_over_g2_with_kappa"] = S(lam4_with)
    R["lambda4_sign"] = along.sign() if along.is_real() else "complex!"
    R["time_s"] = round(time.time() - t0, 1)
    store = {"u": u, "uhat": uhat, "Q": Q, "Fu": Fu, "Nf": Nf, "xi": xi, "Gf": Gf, "kappa": kappa, "ker": ker,
             "Lm": Lm, "tang": tang, "rot": rot, "grads": grads, "r6": r6, "perpF": perpF}
    return R, store


if __name__ == "__main__":
    names = sys.argv[1:] or list(U)
    allR = {}
    stores = {}
    for name in names:
        for d in (3, 4):
            R, st = analyze(name, d)
            allR[f"{name}_s{d}"] = R
            stores[(name, d)] = st
            print(f"== {name} sector {d}")
            for k, v in R.items():
                print("  ", k, ":", v)
            sys.stdout.flush()
    tag = "_".join(names) if sys.argv[1:] else "all"
    with open(f"out/main_{tag}.json", "w") as f:
        json.dump(allR, f, indent=1)
    with open(f"out/main_{tag}.pkl", "wb") as f:
        pickle.dump(stores, f)
