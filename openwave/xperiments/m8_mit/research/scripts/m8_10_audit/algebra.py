"""Exact route for items 2-8 and the coefficient-space part of item 9.

Sections in a sector are handled through the isometric embedding psi -> psi eta^dagger, so the right
part of every section is a (2J+1) x 7 matrix whose columns are indexed like the columns of the exact
projector P = eta eta^dagger.  Products of D functions factor as (left CG coupling) x (right CG
coupling), so every quantity is a sum over coupling paths of (left vector) x (right matrix).
Right parts depend only on the sector; left parts only on the fibre vector.
"""
import json
import sys

import sympy as sp

import lib
import sectors

C = lib.Checker("algebra")
FLAGS = []


def canon(e):
    return sp.expand(e)


def is_zero(e):
    e = canon(e)
    if e == 0:
        return True
    val = abs(sp.N(e, 80))
    if val < sp.Float("1e-60"):
        e2 = sp.nsimplify(sp.radsimp(e))
        if e2 == 0:
            return True
        FLAGS.append(f"numerically zero but not symbolically: {e}")
        return True
    return False


def vzero(v):
    return v is None or all(is_zero(x) for x in v)


def theta_cols(Y):
    return [lib.theta(col) for col in Y]


def dot_couple(Y1, Y2, K):
    out = None
    for a in range(len(Y1)):
        z = lib.couple(Y1[a], Y2[a], K)
        if z is None:
            return None
        out = z if out is None else [p + q for p, q in zip(out, z)]
    return [canon(x) for x in out]


def sv_couple(y, Y, J):
    cols = [lib.couple(y, col, J) for col in Y]
    if cols[0] is None:
        return None
    return cols


def fro(Y1, Y2):
    return canon(sum(lib.inner(c1, c2) for c1, c2 in zip(Y1, Y2)))


def scalar_times_P(Z, Pcols, d):
    """Return c with Z == c P exactly (asserted)."""
    c = canon(fro(Pcols, Z) / d)
    for col, pc in zip(Z, Pcols):
        for z, p in zip(col, pc):
            if not is_zero(z - c * p):
                raise AssertionError("right part is not a multiple of P")
    return c


def mu_of(J):
    return 4 * J * (J + 1) - 48


# ------------------------------------------------------------------ right parts (sector only)
def right_parts(P, d):
    Pcols = [[P[i, a] for i in range(7)] for a in range(7)]
    TP = theta_cols(Pcols)
    B = {K: dot_couple(TP, Pcols, K) for K in range(7)}
    Bnz = [K for K in range(7) if not vzero(B[K])]
    R = {}
    for K in Bnz:
        for J in range(abs(K - 3), K + 4):
            R[(K, J)] = sv_couple(B[K], Pcols, J)
    cblock = {K: scalar_times_P(R[(K, 3)], Pcols, d) for K in Bnz if (K, 3) in R}
    # DN paths
    cdn = {}
    for (K, J), Y in R.items():
        if J == 3:
            continue
        TY = theta_cols(Y)
        for K2 in range(abs(J - 3), J + 4):
            if K2 > 6:
                continue
            y1 = dot_couple(TP, Y, K2)
            if y1 is not None and not vzero(y1):
                cdn[("T1", K, J, K2)] = scalar_times_P(sv_couple(y1, Pcols, 3), Pcols, d)
            y2 = dot_couple(TY, Pcols, K2)
            if y2 is not None and not vzero(y2):
                cdn[("T2", K, J, K2)] = scalar_times_P(sv_couple(y2, Pcols, 3), Pcols, d)
        for K2 in Bnz:
            if abs(K2 - J) <= 3 <= K2 + J:
                Z = sv_couple(B[K2], Y, 3)
                cdn[("T3", K, J, K2)] = scalar_times_P(Z, Pcols, d)
    return {"Pcols": Pcols, "B": B, "Bnz": Bnz, "R": R, "cblock": cblock, "cdn": cdn}


# ------------------------------------------------------------------ left parts (fibre vector only)
def left_parts(v, rp):
    tv = lib.theta(v)
    A = {K: lib.couple(tv, v, K) for K in range(7)}
    L = {(K, J): lib.couple(A[K], v, J) for (K, J) in rp["R"]}
    ldn = {}
    for key in rp["cdn"]:
        term, K, J, K2 = key
        x = L[(K, J)]
        if term == "T1":
            ldn[key] = lib.couple(lib.couple(tv, x, K2), v, 3)
        elif term == "T2":
            ldn[key] = lib.couple(lib.couple(lib.theta(x), v, K2), v, 3)
        else:
            ldn[key] = lib.couple(A[K2], x, 3)
    return {"A": A, "L": L, "ldn": ldn}


def M_K(v, K):
    return lib.couple(lib.couple(v, lib.theta(v), K), v, 3)


def nice(e):
    e = canon(e)
    try:
        e2 = sp.nsimplify(sp.radsimp(e))
        if sp.N(e2 - e, 50) == 0 or abs(sp.N(e2 - e, 50)) < 1e-45:
            return e2
    except Exception:
        pass
    return e


# ------------------------------------------------------------------ rays
ct, st = sp.sqrt(13) / 5, 2 * sp.sqrt(3) / 5


def ray(name):
    v = [sp.Integer(0)] * 7
    if name == "R1":
        v[6] = 1
    elif name == "R2":
        v[3] = 1
    elif name == "R3":
        v[5] = 1; v[1] = 1
    elif name == "R4":
        v[6] = 1; v[0] = 1
    elif name == "R5":
        v[5] = ct; v[0] = st
    return v


def main(ident_name="primary"):
    ident = lib.su2_primary if ident_name == "primary" else lib.su2_secondary
    P, info = sectors.exact_projectors(ident=ident)
    out = {"identification": ident_name}
    homdims = json.load(open(lib.ROOM / "item01.json"))["item1"]["homdims"]
    for d in (3, 4):
        rp = right_parts(P[d], d)
        sec = {"B_nonzero_K": rp["Bnz"],
               "c_block": {K: str(nice(c)) for K, c in rp["cblock"].items()}}
        print(f"\n=== sector dim {d}: K with nonzero right invariant [sum_b Theta eta_b x eta_b]_K:",
              rp["Bnz"])
        # item 2: w(v) = sum_K (-1)^K c_K M_K(v)
        coeffs = {K: nice((-1) ** K * rp["cblock"][K]) for K in rp["cblock"]}
        print("item 2: w(v) = sum_K kappa_K M_K(v), kappa =", coeffs)
        sec["item2_kappa"] = {K: str(c) for K, c in coeffs.items()}
        sec["item2_kappa_float"] = {K: float(sp.N(c)) for K, c in coeffs.items()}
        # right-part level norms r_J = |R_(6,J)|_F^2 (sector constants)
        sec["rJ"] = {J: str(nice(fro(rp["R"][(6, J)], rp["R"][(6, J)]))) for (K, J) in rp["R"] if K == 6}
        rays = {}
        for rn in ("R1", "R2", "R3", "R4", "R5"):
            v0 = ray(rn)
            n2 = lib.inner(v0, v0)
            s = sp.sqrt(sp.Rational(7, d) / n2)
            v = [canon(s * x) for x in v0]
            lp = left_parts(v, rp)
            # check L_(K,3) == (-1)^K M_K(v) (ordering identity used in item 2)
            if rn == "R5":
                # ordering identity for ALL K = 0..6 (odd K make the sign visible)
                tv_ = lib.theta(v)
                Lall = {K: lib.couple(lib.couple(tv_, v, K), v, 3) for K in range(7)}
                ok = all(vzero([a - (-1) ** K * b for a, b in zip(Lall[K], M_K(v, K))]) for K in range(7))
                C.check(f"[{d},{rn}] [[Theta v x v]_K x v]_3 == (-1)^K M_K(v), K=0..6", lambda z: z, ok,
                        all(vzero([a - b for a, b in zip(Lall[K], M_K(v, K))]) for K in range(7)),
                        "drop the (-1)^K")
            w = [0] * 7
            for K in coeffs:
                mk = M_K(v, K)
                w = [a + coeffs[K] * b for a, b in zip(w, mk)]
            w = [canon(x) for x in w]
            # direct block part from paths (independent of M_K bookkeeping)
            wdirect = [0] * 7
            for K, c in rp["cblock"].items():
                wdirect = [a + c * b for a, b in zip(wdirect, lp["L"][(K, 3)])]
            C.check(f"[{d},{rn}] w from M_K formula == w from coupling paths",
                    lambda z: vzero([a - b for a, b in zip(z, wdirect)]), w,
                    [a * 2 for a in w], "double w")
            muv = nice(lib.inner(v, w) / lib.inner(v, v))
            resid = [canon(a - muv * b) for a, b in zip(w, v)]
            is_mult = vzero(resid)
            perp_w = nice(sum(abs(sp.N(x, 40)) ** 2 for x in resid))
            # item 4: level norms
            levels = {}
            for J in range(0, 10):
                keys = [k for k in rp["R"] if k[1] == J]
                tot = 0
                for k1 in keys:
                    for k2 in keys:
                        tot += lib.inner(lp["L"][k1], lp["L"][k2]) * fro(rp["R"][k1], rp["R"][k2])
                normN2 = nice(canon(tot) / (2 * J + 1))
                levels[J] = normN2
            xi_norms = {2 * J: nice(levels[J] / mu_of(J) ** 2) for J in levels if J != 3}
            lam4_identity = nice(-3 * sum(levels[J] / mu_of(J) for J in levels if J != 3))
            lam2 = nice(sp.Rational(d, 7) * lib.inner(v, w))   # <Phi, N(Phi)> = int |Phi|^4 (g=1)
            # item 6: Pi_6 DN_Phi[xi] (g = 1 in xi)
            u = [0] * 7
            for key, c in rp["cdn"].items():
                term, K, J, K2 = key
                fac = -sp.Integer(1) / mu_of(J)
                u = [a + fac * c * b for a, b in zip(u, lp["ldn"][key])]
            u = [canon(x) for x in u]
            along = nice(sp.Rational(d, 7) * lib.inner(v, u))
            coef = canon(lib.inner(v, u) / lib.inner(v, v))
            uperp = [canon(a - coef * b) for a, b in zip(u, v)]
            perp_norm2 = nice(sp.Rational(d, 7) * lib.inner(uperp, uperp))
            C.check(f"[{d},{rn}] <Phi,DN_Phi[xi]> == 3<N,xi> == -3 sum ||Pi_n N||^2/mu_n",
                    lambda z: is_zero(z - lam4_identity), along, along + sp.Rational(1, 10**6),
                    "shift by 1e-6")
            rec = {"v_normalised": [str(nice(x)) for x in v],
                   "w": [str(nice(x)) for x in w],
                   "item3_is_multiple": bool(is_mult),
                   "item3_multiple": str(muv) if is_mult else None,
                   "item3_multiple_float": float(sp.N(muv)) if is_mult else None,
                   "item3_residual_norm2_if_not": str(perp_w),
                   "lambda2_over_g": str(lam2),
                   "normPiN2_by_J": {J: str(levels[J]) for J in levels},
                   "item4_normPi_xi2_over_g2_by_n": {n: str(x) for n, x in xi_norms.items()},
                   "item4_float": {n: float(sp.N(x)) for n, x in xi_norms.items()},
                   "item6_along": str(along), "item6_along_float": float(sp.N(along)),
                   "item6_perp_fibre": [str(nice(x)) for x in uperp],
                   "item6_perp_norm2": str(perp_norm2),
                   "item6_perp_norm2_float": float(sp.N(perp_norm2)),
                   "item7_lambda4_over_g2": str(lam4_identity),
                   "item7_float": float(sp.N(lam4_identity))}
            print(f"[{d}-dim, {rn}] multiple? {is_mult} mu={muv} ; lambda4/g^2 = {lam4_identity} "
                  f"= {float(sp.N(lam4_identity)):.12g}; perp |.|^2 = {perp_norm2}")
            print("    ||Pi_n xi||^2/g^2:", {n: str(x) for n, x in xi_norms.items()})
            # item 9 (coefficient space, exact): Casimir from spin matrices on the left index
            if True:
                maxres = 0
                for J in range(0, 10):
                    keys = [k for k in rp["R"] if k[1] == J]
                    if not keys:
                        continue
                    Jz, Jp, Jm = lib.spin_mats_sym(2 * J)
                    Cas = (Jz * Jz + (Jp * Jm + Jm * Jp) / 2).applyfunc(canon)
                    # coefficient tensors of N and xi at level J (left index M, columns (k, a))
                    for a in range(7):
                        NJ = sp.zeros(2 * J + 1, 2 * J + 1)
                        for k in keys:
                            x = sp.Matrix(lp["L"][k])
                            y = sp.Matrix(rp["R"][k][a])
                            NJ += x * y.T
                        XiJ = sp.zeros(2 * J + 1, 2 * J + 1) if J == 3 else -NJ / mu_of(J)
                        lap = 4 * Cas * XiJ                      # -Delta xi via left generators
                        lapR = 4 * XiJ * Cas.T                   # -Delta xi via right generators
                        blockPhi = sp.zeros(2 * J + 1, 2 * J + 1)
                        if J == 3:
                            blockPhi = sp.Matrix(v) * sp.Matrix(rp["Pcols"][a]).T
                        res = (lap - 48 * XiJ - (lam2 * blockPhi if is_mult else 0 * blockPhi) + NJ)
                        resR = (lapR - lap)
                        for M_ in (res, resR):
                            for e in M_:
                                if not is_zero(e):
                                    maxres = max(maxres, float(abs(sp.N(e, 30))))
                rec["item9_exact_coefficient_residual_max"] = maxres
                print(f"    item 9 exact coefficient residual max |.| = {maxres}")
            rays[rn] = rec
        sec["rays"] = rays
        out[str(d)] = sec
    return out


if __name__ == "__main__":
    ident_name = sys.argv[1] if len(sys.argv) > 1 else "primary"
    res = main(ident_name)
    res["flags"] = FLAGS
    res["checks"] = C.records
    res["all_ok"] = C.all_ok()
    json.dump(res, open(lib.ROOM / f"algebra_{ident_name}.json", "w"), indent=1, default=str)
    print("FLAGS:", FLAGS)
    print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
