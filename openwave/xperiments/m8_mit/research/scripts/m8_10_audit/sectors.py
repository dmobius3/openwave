"""Exact sector projectors in V_3 from a class sum of Gamma (commutant), no tables."""
import numpy as np
import sympy as sp

import lib


def class_sum_D3(cls, ident=lib.su2_primary):
    S = sp.zeros(7, 7)
    for h in cls:
        S += lib.D_sym(6, lib.su2_sym(h, ident))
    return S.applyfunc(sp.expand)


def exact_projectors(ident=lib.su2_primary, verbose=False):
    G = lib.group()
    classes = lib.conj_classes(G)
    # choose a class whose class sum in D^3 has two distinct eigenvalues (numerically first)
    for cls in sorted(classes, key=len):
        Sn = sum(lib.D_num(6, *np.array(lib.su2_num(h, ident)).reshape(4)) for h in cls)
        ev = np.linalg.eigvalsh((Sn + Sn.conj().T) / 2)
        if np.max(ev) - np.min(ev) > 1e-6 and np.allclose(Sn, Sn.conj().T):
            break
    S = class_sum_D3(cls, ident)
    # minimal polynomial S^2 = alpha S + beta I, solved exactly from entries
    al, be = sp.symbols("al be")
    S2 = (S * S).applyfunc(sp.expand)
    eqs = []
    for i in range(7):
        for k in range(7):
            e = sp.expand(S2[i, k] - al * S[i, k] - (be if i == k else 0))
            if e != 0:
                eqs.append(e)
    sol = sp.solve(eqs[:20], [al, be], dict=True)
    assert len(sol) == 1
    alv, bev = sp.nsimplify(sol[0][al]), sp.nsimplify(sol[0][be])
    R = (S2 - alv * S - bev * sp.eye(7)).applyfunc(sp.expand)
    assert all(R[i, k] == 0 for i in range(7) for k in range(7)), "minimal polynomial not exact"
    lam = sp.solve(sp.Symbol("t") ** 2 - alv * sp.Symbol("t") - bev, sp.Symbol("t"))
    l1, l2 = [sp.nsimplify(x) for x in lam]
    inv = sp.radsimp(1 / (l1 - l2))      # rationalised: entries become canonical radical sums
    P1 = ((S - l2 * sp.eye(7)) * inv).applyfunc(lambda e: sp.expand(sp.radsimp(sp.expand(e))))
    P2 = (sp.eye(7) - P1).applyfunc(sp.expand)
    r1 = sp.nsimplify(P1.trace())
    out = {}
    for P, r in ((P1, r1), (P2, 7 - r1)):
        out[int(r)] = P
    info = {"class_size": len(cls), "class_rep": cls[0], "eigs": (l1, l2),
            "min_poly": (alv, bev)}
    if verbose:
        print("class size", len(cls), "eigenvalues", l1, l2)
    return out, info
