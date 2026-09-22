"""Independent route: rhat6(u)*||u||^4 = <U, H U> with U = u (x) u in Sym^2(C^7).

H is a Hermitian form on the 28-dim Sym^2.  Because rhat6 is rotation
invariant, H commutes with the SO(3) action on Sym^2 = V6 + V4 + V2 + V0, so by
Schur it is a scalar h_K on each irrep.  Then

      min_K h_K  <=  rhat6(u)  <=  max_K h_K

for EVERY u, and equality holds exactly on those u whose square u(x)u lies in
the corresponding eigenspace.  This is an exact certificate, not a search.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append
rng = np.random.default_rng(11)

# ---- B_Q with  rho_Q = <u, B_Q u> ----------------------------------------
B = {}
for Q in range(-6, 7):
    M = np.zeros((7, 7), dtype=complex)
    for m1, g in C.CG_NUM[Q].items():
        m2 = Q - m1
        eps = (-1.0) ** (3 - m2)
        M[C.IDX[m1 - Q], C.IDX[m1]] = g.real * eps      # cbar index m1-Q = -m2
    B[Q] = M

u = rng.normal(size=7) + 1j * rng.normal(size=7)
chk = max(abs(np.vdot(u, B[Q] @ u) - C.rho6(u)[Q + 6]) for Q in range(-6, 7))
P(f"B_Q reproduces rho_Q            : {chk:.3e}")

# ---- H on C^7 (x) C^7:  H[(a,c),(b,d)] = sum_Q B_ab conj(B_dc) ------------
H = np.zeros((49, 49), dtype=complex)
for Q in range(-6, 7):
    Bq = B[Q]
    for a in range(7):
        for c_ in range(7):
            for b in range(7):
                for d in range(7):
                    H[a * 7 + c_, b * 7 + d] += Bq[a, b] * np.conj(Bq[d, c_])

err = 0.0
for _ in range(30):
    u = rng.normal(size=7) + 1j * rng.normal(size=7)
    U = np.kron(u, u)
    err = max(err, abs(np.vdot(U, H @ U).real - C.rhat6(u) * np.vdot(u, u).real ** 2))
P(f"<U,H U> = rhat6*||u||^4 max err : {err:.3e}")

# symmetric subspace basis
pairs = [(i, j) for i in range(7) for j in range(i, 7)]
S = np.zeros((49, len(pairs)), dtype=complex)
for k, (i, j) in enumerate(pairs):
    if i == j:
        S[i * 7 + j, k] = 1.0
    else:
        S[i * 7 + j, k] = S[j * 7 + i, k] = 1 / np.sqrt(2)
Hs = S.conj().T @ H @ S
P(f"Sym^2 dim                       : {Hs.shape[0]}")
P(f"||Hs - Hs^dagger||              : {np.abs(Hs-Hs.conj().T).max():.3e}")
Hs = (Hs + Hs.conj().T) / 2
w = np.linalg.eigvalsh(Hs)
P("")
P("eigenvalues of H on Sym^2 (times 924):")
lev, cur = [], [w[0]]
for x in w[1:]:
    if x - cur[-1] > 1e-9:
        lev.append((np.mean(cur), len(cur)))
        cur = [x]
    else:
        cur.append(x)
lev.append((np.mean(cur), len(cur)))
for v, n in lev:
    P(f"   value = {v!r}   x924 = {v*924!r}   multiplicity = {n}")
P(f"lambda_min = {w[0]!r},  lambda_max = {w[-1]!r}")
P(f"463/924    = {463/924!r}")
P(f"1/924      = {1/924!r}")
P(f"lambda_max - 463/924 = {w[-1]-463/924:.3e}    lambda_min - 1/924 = {w[0]-1/924:.3e}")

# ---- identify which irrep carries which eigenvalue, exactly ---------------
P("")
P("exact h_K by evaluating the form on the highest-weight vector of each V_K")


def irrep_vec(K, Qv):
    """|K Q> inside V3 (x) V3 as a 49-vector (exact sympy)."""
    v = sp.zeros(49, 1)
    for m1 in C.MS:
        m2 = Qv - m1
        if m2 in C.IDX:
            v[C.IDX[m1] * 7 + C.IDX[m2]] = clebsch_gordan(3, 3, K, m1, m2, Qv)
    return v


Hex = sp.zeros(49, 49)
Bex = {}
for Q in range(-6, 7):
    M = sp.zeros(7, 7)
    for m1, g in C.CG_EXACT[Q].items():
        m2 = Q - m1
        M[C.IDX[m1 - Q], C.IDX[m1]] = g * sp.Integer((-1) ** (3 - m2))
    Bex[Q] = M

for K in (0, 2, 4, 6):
    vK = irrep_vec(K, K)
    nrm = sp.simplify((vK.T * vK)[0])
    if nrm == 0:
        P(f"  K={K}: highest weight vector vanishes?!")
        continue
    # <v, H v> = sum_Q |<v, (B_Q x I?)...>| -- evaluate directly from the 49x49 form
    tot = sp.Integer(0)
    for Q in range(-6, 7):
        Bq = Bex[Q]
        s = sp.Integer(0)
        for a in range(7):
            for c_ in range(7):
                va = vK[a * 7 + c_]
                if va == 0:
                    continue
                for b in range(7):
                    for d in range(7):
                        vb = vK[b * 7 + d]
                        if vb == 0:
                            continue
                        s += sp.conjugate(va) * Bq[a, b] * sp.conjugate(Bq[d, c_]) * vb
        tot += s
    hK = sp.nsimplify(sp.simplify(tot / nrm))
    P(f"  h_{K} = {hK} = {sp.nsimplify(hK*924)}/924 = {sp.N(hK,20)}   (mult {2*K+1})")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s3.txt").write_text(txt + "\n")
np.save(HERE / "Hsym.npy", Hs)
