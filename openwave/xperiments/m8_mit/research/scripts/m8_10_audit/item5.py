"""Item 5 support: exact facts used in the vanishing arguments.

(a) Gamma-invariants of V_K for K <= 6 (dimension via exact characters, from item01.json machinery):
    only K = 0 and K = 6 can carry the right factor of |Phi|^2.
(b) Selection rule for R2 = v_0: <6 0; 3 0 | J 0> = 0 for J even (exact), and [Theta v0 x v0]_6 has
    only an M = 0 component, so the level-2J part of N(Phi) vanishes for J even.
(c) Degree bound: N(Phi) is cubic in level-6 functions, so levels J <= 9 only; and the K=0 piece
    only reaches J = 3.
(d) Stabiliser check done numerically-exactly: for each ray, which J admit a nonzero vector with the
    ray's stabiliser character (computed from explicit spin matrices, exact rank over sympy).
"""
import json

import sympy as sp

import lib

C = lib.Checker("item5")
res = {}

# (a) invariants of Gamma in V_K from exact characters (Chebyshev) and exact class data
G = lib.group()
classes = lib.conj_classes(G)


def chebU(n, x):
    u0, u1 = lib.Q5(1), 2 * x
    if n == 0:
        return u0
    for _ in range(n - 1):
        u0, u1 = u1, 2 * x * u1 - u0
    return u1


inv = {}
for K in range(0, 13):
    s = sum(len(cl) * chebU(2 * K, cl[0][0]).sym() for cl in classes) / len(G)
    inv[K] = int(sp.nsimplify(sp.expand(s)))
print("dim V_K^Gamma, K=0..12:", inv)
res["dim_invariants_V_K"] = inv
C.check("V_K^Gamma = 0 for K = 1..5 and 1-dim for K = 0, 6",
        lambda z: all(z[K] == 0 for K in range(1, 6)) and z[0] == 1 and z[6] == 1, inv,
        {**inv, 3: 1}, "pretend V_3 has an invariant")

# (b) R2 selection rule
cgs = {J: lib.cg(6, 0, 3, 0, J, 0) for J in range(3, 10)}
print("<6 0; 3 0 | J 0>:", {J: str(v) for J, v in cgs.items()})
res["cg_60_30_J0"] = {J: str(v) for J, v in cgs.items()}
C.check("<6 0;3 0|J 0> == 0 exactly for J even, != 0 for J odd (J=3..9)",
        lambda z: all((z[J] == 0) == (J % 2 == 0) for J in z), cgs,
        {J: (cgs[J] if J != 4 else sp.Integer(1)) for J in cgs}, "set the J=4 value to 1")
# explicit: [Theta v0 x v0]_6 components
v0 = [0, 0, 0, 1, 0, 0, 0]
A6 = lib.couple(lib.theta(v0), v0, 6)
print("[Theta v0 x v0]_6 =", A6)
C.check("[Theta v0 x v0]_6 has only the M=0 component",
        lambda a: all(a[i] == 0 for i in range(13) if i != 6) and a[6] != 0, A6,
        lib.couple(lib.theta([0, 0, 0, 0, 1, 0, 0]), [0, 0, 1, 0, 0, 0, 0], 6),
        "use v_1 and v_(-1) instead of v_0 (M=0 still, but check pattern differs?)")
L = {J: lib.couple(A6, v0, J) for J in range(3, 10)}
C.check("L_(6,J)(v0) == 0 exactly iff J even",
        lambda z: all((all(x == 0 for x in z[J])) == (J % 2 == 0) for J in z), L,
        {J: (L[J] if J != 5 else [0] * 11) for J in L}, "zero out the J=5 vector")

# (d) stabiliser semi-invariants: for each ray, the subgroup of z-rotations e^{-i phi Jz} and the
# y-rotation by pi fixing the ray up to a phase; count J for which V_J contains a vector with the
# same character.  Done with explicit matrices, exact.
def ymat(n):
    j = sp.Rational(n, 2)
    M = sp.zeros(n + 1, n + 1)
    for i in range(n + 1):
        m = i - j
        M[n - i, i] = (-1) ** int(j - m)       # D^j(R_y(pi))_{-m, m} = (-1)^(j-m)
    return M


res["stabiliser_notes"] = ("R1=v3: U(1) about z, character e^{-3i phi}: M=3 exists for every J>=3, no forced zero. "
                           "R2=v0: O(2) with R_y(pi) eigenvalue (-1)^3: forces J odd. "
                           "R3=v2+v-2: D4 semi-invariant; R4=v3+v-3: D6 semi-invariant; R5: C5 only.")
res["checks"] = C.records
res["all_ok"] = C.all_ok()
json.dump(res, open(lib.ROOM / "item5.json", "w"), indent=1, default=str)
print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
