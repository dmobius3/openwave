"""Independent route: sections as explicit functions of a unit quaternion, NO CG coefficients.

eta from group averaging (commutant of {D^3(h)}); level projections by exact cubature on SU(2);
xi, DN_Phi[xi], lambda_4 computed from pointwise function values.  Item 9 pointwise residual with
-Delta computed two ways: (a) Casimir from spin-J matrices on the level coefficients, evaluated at
random points, (b) 9-point finite differences along the geodesics g exp(t e_a), e_a = i, j, k
(the spin-1/2 generators through the quaternion identification).
"""
import json
import sys

import numpy as np

import lib

IDENT = sys.argv[1] if len(sys.argv) > 1 else "primary"
identf = lib.su2_primary if IDENT == "primary" else lib.su2_secondary
C = lib.Checker(f"functional_{IDENT}")
rng = np.random.default_rng(2026)

G = lib.group()
Us = np.array([lib.su2_num(h, identf) for h in G])
D3G = lib.D_num(6, Us[:, 0, 0], Us[:, 0, 1], Us[:, 1, 0], Us[:, 1, 1])

# ---- sectors by commutant averaging
X = rng.normal(size=(7, 7)) + 1j * rng.normal(size=(7, 7))
X = X + X.conj().T
Cm = np.einsum("gij,jk,glk->il", D3G, X, D3G.conj()) / len(G)
ev, V = np.linalg.eigh(Cm)
gaps = np.diff(ev)
cut = int(np.argmax(gaps)) + 1
groups = {cut: V[:, :cut], 7 - cut: V[:, cut:]}
eta = {d: groups[d] for d in (3, 4)}
spread = max(np.ptp(ev[:cut]), np.ptp(ev[cut:]))
C.check(f"commutant has exactly two eigenvalue clusters (spread {spread:.1e}, gap {gaps.max():.2f})",
        lambda z: z[0] < 1e-10 and z[1] > 1e-3, (spread, gaps.max()), (np.ptp(ev), gaps.max()),
        "treat all 7 eigenvalues as one cluster")
sigma = {d: np.einsum("id,gij,je->gde", eta[d].conj(), D3G, eta[d]) for d in eta}
for d in eta:
    err = np.abs(np.einsum("gij,jd->gid", D3G, eta[d]) - np.einsum("id,gde->gie", eta[d], sigma[d])).max()
    C.check(f"D^3(h) eta = eta sigma(h) for all 120 h, d={d} (err {err:.1e})", lambda e: e < 1e-12,
            err, np.abs(np.einsum("gij,jd->gid", D3G, eta[d]) - eta[d][None]).max(),
            "replace sigma(h) by identity")

# ---- cubature on SU(2): |a|^2 = t uniform; phases trapezoid (MPH points), t Gauss-Legendre (NT)
# exact for polynomials of degree < MPH in the phases and t-degree <= 2*NT-1.
# largest integrand: N (deg 18) x conj D^10 (deg 20) = 38; DN (deg 30) x conj D^3 = 36.
MPH, NT = 40, 12
JMAX = 10           # one level beyond the largest possible level of N(Phi) (J=9), to see it vanish
xg, wg = np.polynomial.legendre.leggauss(NT)
tt = (xg + 1) / 2
wt = wg / 2
th = 2 * np.pi * np.arange(MPH) / MPH
T, T1, T2 = np.meshgrid(tt, th, th, indexing="ij")
W = (np.broadcast_to(wt[:, None, None], T.shape) / MPH ** 2).ravel()
A = (np.sqrt(T) * np.exp(1j * T1)).ravel()
B = (np.sqrt(1 - T) * np.exp(1j * T2)).ravel()
UA, UB, UC, UD = A, B, -B.conj(), A.conj()
D3 = lib.D_num(6, UA, UB, UC, UD)


def DJ_grid(J):
    return lib.D_num(2 * J, UA, UB, UC, UD)


# cubature exactness sanity within the design degree (<= 38 < MPH): Schur orthogonality for J <= 9
# (degree 4J <= 36) and the cross-level pair (9, 10) (degree 38).  J=10 with itself is degree 40,
# beyond the grid, and is deliberately not used as a check.
err = 0
for J in (0, 3, 9):
    DJJ = D3 if J == 3 else DJ_grid(J)
    Gm = np.einsum("p,pij,pkl->ijkl", W, DJJ, DJJ.conj())
    ref = np.einsum("ik,jl->ijkl", np.eye(2 * J + 1), np.eye(2 * J + 1)) / (2 * J + 1)
    err = max(err, np.abs(Gm - ref).max())
D9 = DJ_grid(9)
err = max(err, np.abs(np.einsum("p,pij,pkl->ijkl", W, D9, DJ_grid(10).conj())).max())
del D9
C.check(f"cubature reproduces Schur orthogonality up to J=10 (err {err:.1e})", lambda e: e < 1e-12,
        err, np.abs(np.einsum("p,pij,pkl->ijkl", W, D3, D3)).max() + err,
        "drop the complex conjugate in the Schur integral")


def coeffs(f, J, DJJ):
    """level-J coefficient tensor of f (points x comps): (2J+1) int f_a conj(D^J_MK)."""
    return (2 * J + 1) * np.einsum("p,pmk,pa->mka", W, DJJ.conj(), f)


def evalc(Ct, Jp):
    return np.einsum("pmk,mka->pa", Jp, Ct)


def rayv(name):
    v = np.zeros(7)
    ct, st = np.sqrt(13) / 5, 2 * np.sqrt(3) / 5
    if name == "R1":
        v[6] = 1
    elif name == "R2":
        v[3] = 1
    elif name == "R3":
        v[5] = v[1] = 1
    elif name == "R4":
        v[6] = v[0] = 1
    else:
        v[5], v[0] = ct, st
    return v.astype(complex)


mu = lambda J: 4 * J * (J + 1) - 48
out = {"identification": IDENT}
for d in (3, 4):
    et = eta[d]
    sec = {}
    for rn in ("R1", "R2", "R3", "R4", "R5"):
        v = rayv(rn)
        Phi = np.einsum("m,pmk,ka->pa", v, D3, et)
        nrm = np.sum(W * np.sum(np.abs(Phi) ** 2, 1))
        v = v / np.sqrt(nrm)      # normalisation of the SECTION by cubature (not by the 7/d formula)
        Phi = Phi / np.sqrt(nrm)
        # equivariance of the explicit section: Phi(g h) = Phi(g) sigma(h), sampled points
        idx = rng.integers(0, len(W), 5)
        e_eq = 0
        for hidx in range(len(G)):
            Uh = Us[hidx]
            Ug = np.array([[UA[idx], UB[idx]], [UC[idx], UD[idx]]]).transpose(2, 0, 1)
            Ugh = Ug @ Uh
            Dgh = lib.D_num(6, Ugh[:, 0, 0], Ugh[:, 0, 1], Ugh[:, 1, 0], Ugh[:, 1, 1])
            lhs = np.einsum("m,pmk,ka->pa", v, Dgh, et)
            e_eq = max(e_eq, np.abs(lhs - Phi[idx] @ sigma[d][hidx]).max())
            if G[hidx] == lib.Q1:
                e_mut = np.abs(lhs - Phi[idx]).max()   # mutant: drop sigma(q1)
        N = np.sum(np.abs(Phi) ** 2, 1)[:, None] * Phi
        CN = {}
        xi = np.zeros_like(Phi)
        for J in range(0, JMAX + 1):
            DJJ = D3 if J == 3 else DJ_grid(J)
            CN[J] = coeffs(N, J, DJJ)
            if J != 3:
                xi = xi - evalc(CN[J], DJJ) / mu(J)
            del DJJ
        normN2 = {J: float(np.sum(np.abs(CN[J]) ** 2).real / (2 * J + 1)) for J in CN}
        # block: C^3 = w (x) eta ?
        w = np.einsum("mka,ka->m", CN[3], et.conj()) / et.shape[1]   # sum_(K,a) |eta_Ka|^2 = d
        e_rank1 = np.abs(CN[3] - np.einsum("m,ka->mka", w, et)).max()
        mult = np.vdot(v, w) / np.vdot(v, v)
        perpw = np.linalg.norm(w - mult * v)
        dotPx = np.sum(Phi.conj() * xi, 1)
        DN = (dotPx + dotPx.conj())[:, None] * Phi + np.sum(np.abs(Phi) ** 2, 1)[:, None] * xi
        CD3 = coeffs(DN, 3, D3)
        u = np.einsum("mka,ka->m", CD3, et.conj()) / et.shape[1]
        e_rank1_dn = np.abs(CD3 - np.einsum("m,ka->mka", u, et)).max()
        along = np.sum(W * np.sum(Phi.conj() * DN, 1))
        normPi3DN2 = np.sum(np.abs(CD3) ** 2).real / 7
        perp2 = normPi3DN2 - abs(along) ** 2
        lam4_id = -3 * sum(normN2[J] / mu(J) for J in normN2 if J != 3)
        xi_norms = {2 * J: normN2[J] / mu(J) ** 2 for J in normN2 if J != 3}
        rec = {"equivariance_err": e_eq, "block_rank1_err": e_rank1, "multiple": mult.real,
               "multiple_imag": mult.imag, "w_perp_norm": perpw, "normPiN2_by_J": normN2,
               "xi_norms_by_n": xi_norms, "along": along.real, "along_imag": along.imag,
               "perp_norm2": perp2, "lam4_identity": lam4_id, "dn_rank1_err": e_rank1_dn,
               "u": [complex(x) for x in u]}
        # ---- item 9 pointwise residual at random points (not cubature nodes)
        q = rng.normal(size=(6, 4))
        q /= np.linalg.norm(q, axis=1)[:, None]

        def section_vals(qs):
            Uq = lib.su2_from_quat_float(*qs.T, ident=IDENT)
            DJq = {J: lib.D_num(2 * J, *Uq) for J in range(0, JMAX + 1)}
            Ph = np.einsum("m,pmk,ka->pa", v, DJq[3], et)
            Xi = sum(-evalc(CN[J], DJq[J]) / mu(J) for J in range(0, JMAX + 1) if J != 3)
            return Ph, Xi, DJq

        Ph, Xi, DJq = section_vals(q)
        lap_cas = 0
        for J in range(0, JMAX + 1):
            if J == 3:
                continue
            n = 2 * J
            Jz = np.diag([i - J for i in range(n + 1)]).astype(complex)
            Jp = np.zeros((n + 1, n + 1), complex)
            for i in range(n):
                m = i - J
                Jp[i + 1, i] = np.sqrt((J - m) * (J + m + 1))
            Jm = Jp.T.copy()
            Cas = Jz @ Jz + (Jp @ Jm + Jm @ Jp) / 2
            Cx = -CN[J] / mu(J)
            lap_cas = lap_cas + evalc(4 * np.einsum("mn,nka->mka", Cas, Cx), DJq[J])
        NPh = np.sum(np.abs(Ph) ** 2, 1)[:, None] * Ph
        lam2 = mult.real
        res_cas = lap_cas - 48 * Xi - lam2 * Ph + NPh
        # finite differences along g*(cos t + e sin t), e in {i, j, k}
        h = 1e-3
        cst = np.array([-1 / 560, 8 / 315, -1 / 5, 8 / 5, -205 / 72, 8 / 5, -1 / 5, 8 / 315, -1 / 560])
        lap_fd = 0
        for e in range(1, 4):
            acc = 0
            for s, cc in zip(range(-4, 5), cst):
                t = s * h
                qe = np.zeros(4)
                qe[0], qe[e] = np.cos(t), np.sin(t)
                # quaternion product q * qe (float)
                w1, x1, y1, z1 = q.T
                w2, x2, y2, z2 = qe
                qp = np.stack([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
                               w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                               w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                               w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2], 1)
                _, Xs, _ = section_vals(qp)
                acc = acc + cc * Xs
            lap_fd = lap_fd + acc / h ** 2
        res_fd = -lap_fd - 48 * Xi - lam2 * Ph + NPh
        scale = np.abs(NPh).max()
        rec["item9_residual_casimir_max"] = float(np.abs(res_cas).max())
        rec["item9_residual_fd_max"] = float(np.abs(res_fd).max())
        rec["item9_scale_maxN"] = float(scale)
        rec["item9_casimir_vs_fd_max"] = float(np.abs(lap_cas + lap_fd).max())
        tag = f"[{IDENT} d={d} {rn}]"
        C.check(f"{tag} Phi(gh) = Phi(g) sigma(h), all h (err {e_eq:.1e})", lambda e: e < 1e-12, e_eq,
                e_mut, "drop sigma(h) for h = q1 (compare Phi(g q1) with Phi(g))")
        C.check(f"{tag} level-6 part of N and of DN is (fibre) x eta (err {max(e_rank1, e_rank1_dn):.1e})",
                lambda e: e < 1e-12, max(e_rank1, e_rank1_dn),
                np.abs(CN[3] - np.einsum("m,ka->mka", w * et.shape[1], et)).max(),
                "reintroduce the factor d in the fibre extraction (the bug found in the first run)")
        sc = scale
        C.check(f"{tag} item 9 residual, Casimir (spin-J matrices), max {rec['item9_residual_casimir_max']:.1e} "
                f"(scale {sc:.2f})", lambda r: r < 1e-10 * max(sc, 1), np.abs(res_cas).max(),
                np.abs(res_cas + 1e-3 * lam2 * Ph).max(), "lambda_2 -> 1.001 lambda_2")
        C.check(f"{tag} item 9 residual, finite differences on geodesics, max {rec['item9_residual_fd_max']:.1e}",
                lambda r: r < 1e-6 * max(sc, 1), np.abs(res_fd).max(),
                np.abs(-lap_fd - 47 * Xi - lam2 * Ph + NPh).max(), "48 -> 47 in (-Delta - 48)")
        sec[rn] = rec
        print(f"[{IDENT} d={d} {rn}] equiv {e_eq:.1e} mult {mult.real:.12f} (|w_perp| {perpw:.2e}) "
              f"lam4/g^2 {along.real:.12f} (identity {lam4_id:.12f}) perp^2 {perp2:.3e} "
              f"res9 cas {rec['item9_residual_casimir_max']:.1e} fd {rec['item9_residual_fd_max']:.1e}")
        print("   ||Pi_n N||^2 by J:", {J: f"{normN2[J]:.12g}" for J in normN2})
    out[str(d)] = sec

out["checks"] = C.records
out["all_ok"] = C.all_ok()
json.dump(out, open(lib.ROOM / f"functional_{IDENT}.json", "w"), indent=1, default=lambda o: str(o))
print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
