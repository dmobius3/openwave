"""Item 5, part b: exact arguments for the zeros of ||Pi_n xi|| at sector levels.

Structure (from the exact route): off the block, N(Phi) has only the K = 6 path, so
    ||Pi_(2J) N||^2 = l_J(v) r_J(sigma) / (2J+1),
    l_J(v) = |[[Theta v x v]_6 x v]_J|^2  (left, fibre only),  r_J = |[I_6 x eta]_J|_F^2 (sector only).
A zero at a sector level is a zero of l_J(v) (r_J != 0 at every sector level is checked from R1).

Arguments checked here, exactly:
 (A) coupling to the top spin is the product of binary forms, and coupling one below the top is
     proportional to the Jacobian (first transvectant); verified exactly on generic vectors.
     If Theta v = eps v then [Theta v x v]_6 ~ f^2 and [f^2 x f]_8 ~ (f^2, f)_1 = 2 f (f_x f_y - f_y f_x) = 0.
 (B) stabiliser semi-invariants: if D^3(h) v = chi(h) v for all h in a finite H, then
     D^J(h) l(v) = chi(h) l(v), so l(v) = 0 whenever V_J has no chi-semi-invariant vector.
     Multiplicities (1/|H|) sum conj(chi(h)) chi_J(h) from explicit matrices.
"""
import json

import numpy as np
import sympy as sp

import lib

C = lib.Checker("item5b")
res = {}
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
    else:
        v[5] = ct; v[0] = st
    return v


rays = ["R1", "R2", "R3", "R4", "R5"]
# ---- left factors l_J(v), exact, unnormalised fibre vectors
lJ = {}
for r in rays:
    v = ray(r)
    A6 = lib.couple(lib.theta(v), v, 6)
    lJ[r] = {J: sp.nsimplify(sp.expand(lib.inner(lib.couple(A6, v, J), lib.couple(A6, v, J))))
             for J in range(3, 10)}
    print(r, "l_J (J=3..9):", {J: str(x) for J, x in lJ[r].items()})
res["l_J_unnormalised"] = {r: {J: str(x) for J, x in lJ[r].items()} for r in rays}

# ---- Theta-eigen test
eps = {}
for r in rays:
    v = ray(r)
    tv = lib.theta(v)
    e = None
    for s in (1, -1):
        if all(sp.expand(a - s * b) == 0 for a, b in zip(tv, v)):
            e = s
    eps[r] = e
print("Theta v = eps v:", eps)
res["theta_eigen"] = eps

# ---- (A) forms: coordinates u_m -> f = sum u_m x^(j+m) y^(j-m) / sqrt((j+m)!(j-m)!)
x, y = sp.symbols("x y")


def form(u):
    j = (len(u) - 1) // 2
    return sp.expand(sum(u[i] * x ** i * y ** (2 * j - i) / sp.sqrt(lib.Nfac(2 * j, i))
                         for i in range(2 * j + 1)))




def const_ratio(P1, P2):
    """c with P1 == c * P2 exactly (None if not proportional): c from one coefficient, then verified."""
    t2 = sp.Poly(P2, x, y).terms()
    mon, c2 = t2[0]
    c1 = sp.Poly(P1, x, y).coeff_monomial(x ** mon[0] * y ** mon[1])
    c = sp.radsimp(c1 / c2)
    return c if sp.expand(sp.radsimp(sp.expand(P1 - c * P2))) == 0 else None

rng = np.random.default_rng(7)
def rvec(j):
    return [sp.Rational(int(k), 7) for k in rng.integers(-9, 10, 2 * j + 1)]

a3, b3 = rvec(3), rvec(3)
prod = sp.expand(form(a3) * form(b3))
top = form(lib.couple(a3, b3, 6))
ratio_top = const_ratio(top, prod)
C.check(f"[a x b]_6 form == const * f_a f_b (const {ratio_top})",
        lambda z: const_ratio(*z) is not None, (top, prod),
        (sp.expand(form(lib.couple(a3, b3, 6)) + x ** 12 / 7), prod), "add x^12/7 to the coupled form")
F6, f3 = rvec(6), rvec(3)
jac = sp.expand(sp.diff(form(F6), x) * sp.diff(form(f3), y) - sp.diff(form(F6), y) * sp.diff(form(f3), x))
c8 = form(lib.couple(F6, f3, 8))
ratio8 = const_ratio(c8, jac)
c7 = form(lib.couple(F6, f3, 7))
# degree mismatch would make the mutant trivially fail; compare instead the spin-8 coupling of a
# different pair (F6, f3') against the Jacobian of (F6, f3): not proportional unless the map is trivial
f3b = rvec(3)
C.check(f"[F x f]_8 form == const * Jacobian(F, f) (const {ratio8})",
        lambda z: const_ratio(*z) is not None, (c8, jac),
        (form(lib.couple(F6, f3b, 8)), jac), "couple F with a different f than the one in the Jacobian")
# generic-vector constants must be vector-independent: repeat on a second pair
a3b, b3b = rvec(3), rvec(3)
r2 = const_ratio(form(lib.couple(a3b, b3b, 6)), sp.expand(form(a3b) * form(b3b)))
F6b, f3c = rvec(6), rvec(3)
jb = sp.expand(sp.diff(form(F6b), x) * sp.diff(form(f3c), y) - sp.diff(form(F6b), y) * sp.diff(form(f3c), x))
r8b = const_ratio(form(lib.couple(F6b, f3c, 8)), jb)
C.check(f"constants independent of the vectors (product {r2}, Jacobian {r8b})",
        lambda z: z[0] is not None and z[1] is not None and sp.expand(z[0] - ratio_top) == 0
        and sp.expand(z[1] - ratio8) == 0, (r2, r8b), (r2, ratio_top), "compare the Jacobian constant with the product constant")
for r in rays:
    if eps[r] is not None:
        C.check(f"{r}: Theta v = eps v, hence l_8(v) = 0 (Jacobian argument); computed l_8 = {lJ[r][8]}",
                lambda z: z == 0, lJ[r][8], lJ["R1"][8], "use R1's l_8 (R1 is not a Theta eigenvector)")

# ---- (B) stabiliser multiplicities from explicit finite groups
def su2_rot(axis, ang):
    n = np.array(axis, float); n /= np.linalg.norm(n)
    sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]])
    return np.cos(ang / 2) * np.eye(2) - 1j * np.sin(ang / 2) * (n[0] * sx + n[1] * sy + n[2] * sz)


def closure_num(gens):
    els = [np.eye(2, dtype=complex)]
    frontier = list(els)
    while frontier:
        new = []
        for e in frontier:
            for g in gens:
                p = e @ g
                if not any(np.abs(p - q).max() < 1e-9 for q in els):
                    els.append(p); new.append(p)
        frontier = new
    return els


groups = {
    "R3": closure_num([su2_rot([0, 0, 1], np.pi / 2), su2_rot([1, 1, 0], np.pi / 2)]),
    "R4": closure_num([su2_rot([0, 0, 1], np.pi / 3), su2_rot([0, 1, 0], np.pi)]),
    "R5": closure_num([su2_rot([0, 0, 1], 2 * np.pi / 5)]),
}
mult = {}
for r, H in groups.items():
    v = np.array([complex(sp.N(t)) for t in ray(r)])
    chis = []
    fix_err = 0
    for h in H:
        Dh = lib.D_num(6, h[0, 0], h[0, 1], h[1, 0], h[1, 1])
        w = Dh @ v
        chi = np.vdot(v, w) / np.vdot(v, v)
        fix_err = max(fix_err, np.abs(w - chi * v).max())
        chis.append(chi)
    chis = np.array(chis)
    C.check(f"{r}: all {len(H)} elements of H fix the ray up to a phase (err {fix_err:.1e})",
            lambda e: e < 1e-10, fix_err,
            np.abs(lib.D_num(6, *su2_rot([1, 0, 0], 0.3).ravel()) @ v - v).max(),
            "test a generic rotation instead")
    mJ = {}
    for J in range(3, 10):
        trs = np.array([np.trace(lib.D_num(2 * J, h[0, 0], h[0, 1], h[1, 0], h[1, 1])) for h in H])
        m = (np.conj(chis) * trs).sum() / len(H)
        assert abs(m - round(m.real)) < 1e-9, (r, J, m)
        mJ[J] = int(round(m.real))
    mult[r] = {"order": len(H), "multiplicity_by_J": mJ}
    print(f"{r}: |H| = {len(H)}, dim of chi-semi-invariants in V_J (J=3..9):", mJ)
    for J in range(3, 10):
        if mJ[J] == 0:
            C.check(f"{r}: no chi-semi-invariant in V_{J}, and l_{J}(v) == 0 exactly",
                    lambda z: z == 0, lJ[r][J], lJ["R1"][J], "use R1's l_J")
res["stabiliser_multiplicities"] = mult
res["checks"] = C.records
res["all_ok"] = C.all_ok()
json.dump(res, open(lib.ROOM / "item5b.json", "w"), indent=1, default=str)
print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
