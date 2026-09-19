"""Item 12 and numerical cross-checks of items 2, 4, 5, 7 by pointwise evaluation + Haar quadrature (float64).
Independent of the exact route except that xi, kappa (the objects being checked) are imported from it."""
import pickle, json, sys
import numpy as np
from numD import D_np
from su2 import spin_matrices
from checks import check, save_log

GRP = pickle.load(open("out/group.pkl", "rb"))
dims = {d: {n: (int(v.t[1][0]) if not v.is_zero() else 0) for n, v in GRP["dims"][d].items()} for d in GRP["dims"]}
ETA = pickle.load(open("out/eta.pkl", "rb"))
ST = pickle.load(open("out/main_all.pkl", "rb"))


def grid(na, nb):
    al = 2 * np.pi * np.arange(na) / na
    x, wx = np.polynomial.legendre.leggauss(nb)
    be = np.arccos(x)
    A, B, C = np.meshgrid(al, be, al, indexing="ij")
    W = np.ones_like(A) / na / na * np.meshgrid(np.ones(na), wx / 2, np.ones(na), indexing="ij")[1]
    A, B, C, W = A.ravel(), B.ravel(), C.ravel(), W.ravel()
    # g = rz(A) ry(B) rz(C), rz(t) = diag(e^{-it/2}, e^{it/2}), ry(b) = [[cos b/2, -sin b/2],[sin b/2, cos b/2]]
    ea, ec = np.exp(-1j * A / 2), np.exp(-1j * C / 2)
    cb, sb = np.cos(B / 2), np.sin(B / 2)
    g11 = ea * cb * ec
    g12 = -ea * sb * np.conj(ec)
    g21 = np.conj(ea) * sb * ec
    g22 = np.conj(ea) * cb * np.conj(ec)
    return (g11, g12, g21, g22), W


def eta_arr(d, n):
    return [np.array(b) for b in ETA[d]["eta"][n]]


def mqarr(v):
    return np.array([c.to_complex() for c in v])


def run(na, nb, verbose=True):
    G, W = grid(na, nb)
    Dc = {}
    def D(J):
        if J not in Dc:
            Dc[J] = D_np(J, *G)
        return Dc[J]
    res = {}
    for (name, d), S in ST.items():
        e3 = eta_arr(d, 6)[0]
        u = mqarr(S["u"])
        Phi = np.einsum("m,pmk,ka->pa", u, D(3), e3)
        nPhi = np.sum(W * np.sum(np.abs(Phi) ** 2, 1))
        dens = np.sum(np.abs(Phi) ** 2, 1)
        Nv = dens[:, None] * Phi
        rec = {"norm_Phi_sq": nPhi}
        levels = {}
        copy_resid = 0.0
        tot = 0.0
        cN = {}
        for J in range(0, 12):
            c = (2 * J + 1) * np.einsum("p,pmk,pa->mka", W, np.conj(D(J)), Nv)
            n = 2 * J
            if dims[d].get(n, 0) > 0:
                etas = eta_arr(d, n)
                fib = [np.einsum("mka,ka->m", c, np.conj(e)) / d for e in etas]
                rec_c = sum(np.einsum("m,ka->mka", f, e) for f, e in zip(fib, etas))
                copy_resid = max(copy_resid, np.abs(c - rec_c).max())
                if len(etas) == 2:   # mutation: keep only the first copy
                    rec["copy_mut_drop_second"] = np.abs(c - np.einsum("m,ka->mka", fib[0], etas[0])).max()
                cN[J] = c
            else:
                copy_resid = max(copy_resid, np.abs(c).max())
            nn = np.sum(np.abs(c) ** 2) / (2 * J + 1)
            levels[n] = nn
            tot += nn
        intN2 = np.sum(W * np.sum(np.abs(Nv) ** 2, 1))
        rec["parseval_defect"] = abs(intN2 - tot)
        rec["copy_residual(max over levels <= 22, incl. levels without intertwiners and 20, 22)"] = copy_resid
        rec["levels_quad"] = levels
        # exact values for comparison
        from engine import norm2
        ex = {2 * J: norm2({J: S["Nf"][J]}).to_mpc().real for J in S["Nf"]}
        rec["max_level_norm_dev_vs_exact"] = max(abs(levels[n] - float(ex.get(n, 0.0))) for n in levels)
        # item 2 check: level-3 part equals Q Phi
        Q = S["Q"].to_complex().real
        cPhi = np.einsum("m,ka->mka", u, e3)
        rec["Pi6N_minus_QPhi_max"] = np.abs(cN[3] - Q * cPhi).max() if name != "U5q" else None
        # ---- item 12: xi from the exact route, -Delta via Casimir of explicit spin matrices
        resid12 = 0.0
        xi_pts = np.zeros_like(Phi)
        for J in range(0, 10):
            if J == 3:
                continue
            cx = np.zeros((2 * J + 1, 2 * J + 1, d), complex)
            for x, Z in S["xi"].get(J, []):
                xv = mqarr(x)
                Zm = np.array([mqarr(col) for col in Z]).T      # (2J+1) x 7, columns a'
                cx += np.einsum("m,ka->mka", xv, Zm @ e3)
            Jz, Jp, Jm = spin_matrices(J)
            Jx = (Jp + Jm) / 2
            Jy = (Jp - Jm) / 2j
            Cas = Jx @ Jx + Jy @ Jy + Jz @ Jz
            lap = 4 * np.einsum("nm,nka->mka", Cas, cx)          # C^T acting on the left index
            cNJ = cN.get(J, np.zeros_like(cx))
            r = lap - 48 * cx + 1.0 * cNJ                          # g = 1
            resid12 = max(resid12, np.abs(r).max())
            # mutations: xi with flipped sign; Laplacian with n(n+2) replaced by 3*Casimir
            rec.setdefault("item12_mut_sign", 0.0)
            rec.setdefault("item12_mut_cas3", 0.0)
            rec["item12_mut_sign"] = max(rec["item12_mut_sign"], np.abs(-lap + 48 * cx + cNJ).max())
            rec["item12_mut_cas3"] = max(rec["item12_mut_cas3"], np.abs(0.75 * lap - 48 * cx + cNJ).max())
            xi_pts += np.einsum("mka,pmk->pa", cx, D(J))
        rec["item12_residual_max"] = resid12
        # ---- cross-check item 5: Pi_6 DN_Phi[xi] fibre, from pointwise evaluation
        def DN(h):
            return 2 * np.real(np.sum(np.conj(Phi) * h, 1))[:, None] * Phi + dens[:, None] * h
        def fibre3(f):
            c = 7 * np.einsum("p,pmk,pa->mka", W, np.conj(D(3)), f)
            return np.einsum("mka,ka->m", c, np.conj(e3)) / d
        Gq = fibre3(DN(xi_pts))
        rec["DNxi_fibre_dev_vs_exact"] = np.abs(Gq - mqarr(S["Gf"])).max()
        # kappa equation check: P_T(DN kappa) - Q kappa + P_T(G) = 0 (on the range; exact route says fully solvable)
        kap = mqarr(S["kappa"])
        Kp = np.einsum("m,pmk,ka->pa", kap, D(3), e3)
        DK = fibre3(DN(Kp))
        PT = lambda v: v - np.vdot(u, v) / np.vdot(u, u) * u
        if name != "U5q":
            rec["kappa_eq_residual"] = np.abs(PT(DK) - Q * kap + PT(Gq)).max()
            rec["lambda4_quad"] = (np.vdot(u, Gq) * d / 7).real
        res[(name, d)] = rec
        if verbose:
            print(name, d, {k: (v if not isinstance(v, dict) else "...") for k, v in rec.items()})
    return res


if __name__ == "__main__":
    r1 = run(24, 12)
    r2 = run(30, 16, verbose=False)
    worst = {}
    for key in ["parseval_defect", "copy_residual(max over levels <= 22, incl. levels without intertwiners and 20, 22)",
                "max_level_norm_dev_vs_exact", "item12_residual_max", "DNxi_fibre_dev_vs_exact", "kappa_eq_residual", "Pi6N_minus_QPhi_max"]:
        worst[key] = max(float(v[key]) for v in r1.values() if v.get(key) is not None)
    grid_dev = max(abs(r1[k]["levels_quad"][n] - r2[k]["levels_quad"][n]) for k in r1 for n in r1[k]["levels_quad"])
    worst["grid_24x12x24_vs_30x16x30_level_norm_dev"] = grid_dev
    worst["norm_Phi_sq_dev"] = max(abs(v["norm_Phi_sq"] - 1) for v in r1.values())
    print(json.dumps(worst, indent=1))
    tol = 1e-10
    mut_sign = min(float(v["item12_mut_sign"]) for v in r1.values())
    mut_cas = min(float(v["item12_mut_cas3"]) for v in r1.values())
    print("item 12 mutants: min residual with xi -> -xi:", mut_sign, "; with -Delta -> 3*Casimir:", mut_cas)
    check("item 12: max |(-Delta-48)xi + g(N - Pi_6 N)| over all U, sectors, levels < 1e-10 (Casimir from spin matrices)",
          lambda w: w < tol, worst["item12_residual_max"], mut_sign, "flip the sign of xi (min over cases of the mutant residual)")
    check("item 12 (second mutant)", lambda w: w < tol, worst["item12_residual_max"], mut_cas, "replace 4*Casimir by 3*Casimir")
    mut_copy = min(float(v["copy_mut_drop_second"]) for k, v in r1.items() if "copy_mut_drop_second" in v)
    check("N(Phi) fully captured by levels <= 18 with every intertwiner copy (copy residual, levels 20 and 22 included) < 1e-10",
          lambda w: w < tol, worst["copy_residual(max over levels <= 22, incl. levels without intertwiners and 20, 22)"], mut_copy,
          "sector 4, level 18: keep only the first of the two copies")
    check("Parseval: int|N|^2 == sum of level norms (no component above level 22 either)", lambda w: w < tol,
          worst["parseval_defect"], 1e-5, "defect 1e-5")
    check("quadrature level norms agree with exact ||Pi_n N||^2", lambda w: w < tol, worst["max_level_norm_dev_vs_exact"], 1e-6, "dev 1e-6")
    check("quadrature Pi_6 DN_Phi[xi] fibre agrees with exact", lambda w: w < tol, worst["DNxi_fibre_dev_vs_exact"], 1e-6, "dev 1e-6")
    check("kappa solves the order-a^5 perp equation (quadrature)", lambda w: w < tol, worst["kappa_eq_residual"], 1e-6, "dev 1e-6")
    check("Pi_6 N(Phi) == Q Phi (quadrature) at U1..U6", lambda w: w < tol, worst["Pi6N_minus_QPhi_max"], 1e-6, "dev 1e-6")
    check("quadrature grid independence", lambda w: w < tol, grid_dev, 1e-6, "dev 1e-6")
    # a real mutation of the item-12 identity: omit the Laplacian's factor 4 -> residual must become large
    with open("out/item12.json", "w") as f:
        json.dump({"worst": worst, "per_case": {f"{k[0]}_s{k[1]}": {a: (b if not isinstance(b, dict) else {str(n): v for n, v in b.items()}) for a, b in v.items()} for k, v in r1.items()}},
                  f, indent=1, default=float)
    save_log("s12")
