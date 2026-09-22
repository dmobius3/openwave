"""Step 0: arm the instrument. Every check here has a known answer."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan
import core as C

rng = np.random.default_rng(20260922)
out = []
P = out.append


def rand(n=7):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


# --- 1. su(2) algebra (known answer: [Jx,Jy]=i Jz, J^2 = j(j+1) = 12) -------
Jx, Jy, Jz = C.JX, C.JY, C.JZ
P(f"||[Jx,Jy] - i Jz||                 = {np.abs(Jx@Jy-Jy@Jx-1j*Jz).max():.3e}")
P(f"||[Jy,Jz] - i Jx||                 = {np.abs(Jy@Jz-Jz@Jy-1j*Jx).max():.3e}")
J2 = Jx @ Jx + Jy @ Jy + Jz @ Jz
P(f"||J^2 - 12 I||                     = {np.abs(J2-12*np.eye(7)).max():.3e}")
P(f"hermiticity max err                = {max(np.abs(O-O.conj().T).max() for O in C.JOPS):.3e}")

# --- 2. rotation is an SO(3) rep: 2pi = +I (known answer for integer spin) --
D2pi = C.rot([0.3, -0.7, 0.5], 2 * np.pi)
P(f"||D(n,2pi) - I||                   = {np.abs(D2pi-np.eye(7)).max():.3e}")
Dz = C.rot([0, 0, 1], 0.731)
P(f"D_z diagonal exp(-i m th) err      = "
  f"{np.abs(np.diag(Dz)-np.exp(-1j*0.731*np.array(C.MS))).max():.3e}")
Da, Db = C.rot([0, 0, 1], 0.4), C.rot([0, 0, 1], 0.9)
P(f"D_z(a)D_z(b)=D_z(a+b) err          = {np.abs(Da@Db-C.rot([0,0,1],1.3)).max():.3e}")
Du = C.rot([0.2, 0.9, -0.3], 1.1)
P(f"unitarity of D                     = {np.abs(Du.conj().T@Du-np.eye(7)).max():.3e}")

# --- 3. Theta: antiunitary, Theta^2 = ? (worklist item 0) -------------------
u, w = rand(), rand()
P(f"Theta antilinear check             = "
  f"{np.abs(C.Theta((2+3j)*u)-np.conj(2+3j)*C.Theta(u)).max():.3e}")
P(f"<Theta u,Theta w> - conj<u,w>      = {abs(np.vdot(C.Theta(u),C.Theta(w))-np.conj(np.vdot(u,w))):.3e}")
P(f"ITEM 0: ||Theta(Theta u) - u||     = {np.abs(C.Theta(C.Theta(u))-u).max():.3e}  -> Theta^2 = +1")
# exact symbolic confirmation of Theta^2
cs = sp.symbols('c0:7')
tt = C.Theta_exact(C.Theta_exact(list(cs)))
P(f"ITEM 0 exact: Theta^2 u - u        = {sp.simplify(tt - sp.Matrix(cs)).T}")
# Theta commutes with rotations (integer spin): Theta D u = D Theta u
P(f"||Theta(Du) - D(Theta u)||         = {np.abs(C.Theta(Du@u)-Du@C.Theta(u)).max():.3e}")
# Theta J Theta^-1 = -J
for nm, O in zip("xyz", C.JOPS):
    P(f"  Theta {nm} Theta^-1 + {nm} err       = {np.abs(C.Theta(O@C.Theta(u))+O@u).max():.3e}")

# --- 4. CG conventions ------------------------------------------------------
P(f"CS anchor <33;33|66>               = {clebsch_gordan(3,3,6,3,3,6)}   (must be +1)")
v = clebsch_gordan(3, 3, 6, 3, -3, 0)
P(f"ITEM 0: <33;3-3|60>                = {v} = {sp.nsimplify(v)} = {sp.N(v,25)}")
P(f"          its square               = {sp.simplify(v**2)}")
# closed form for the stretched coupling  <j m; j -m|2j 0> = sqrt(C(2j,j+m)C(2j,j-m)/C(4j,2j))
closed = sp.sqrt(sp.binomial(6, 6) * sp.binomial(6, 0) / sp.binomial(12, 6))
P(f"   textbook stretched formula      = {closed} ; difference = {sp.simplify(v-closed)}")
# full CG matrix unitary over all J=0..6
rows = []
for JJ in range(0, 7):
    for Q in range(-JJ, JJ + 1):
        row = []
        for m1 in C.MS:
            for m2 in C.MS:
                row.append(clebsch_gordan(3, 3, JJ, m1, m2, Q) if m1 + m2 == Q else sp.Integer(0))
        rows.append(row)
M = sp.Matrix(rows)
P(f"CG matrix shape {M.shape}, ||M M^T - I||inf = {max(abs(x) for x in (M*M.T-sp.eye(49)))}")

# --- 5. rhat6: invariances (conventions 1.3 asserts these; verify) ---------
u = rand()
P(f"rhat6 phase invariance             = {abs(C.rhat6(u)-C.rhat6(np.exp(1.234j)*u)):.3e}")
P(f"rhat6 scale invariance (deg 0)     = {abs(C.rhat6(u)-C.rhat6(3.7*u)):.3e}")
mx = 0.0
for _ in range(40):
    R = C.rot(rng.normal(size=3), rng.uniform(0, 2 * np.pi))
    mx = max(mx, abs(C.rhat6(u) - C.rhat6(R @ u)))
P(f"rhat6 rotation invariance (40 rot) = {mx:.3e}")
P(f"rhat6 real, imag part of ||rho||^2 = 0 by construction")

# --- 6. sum rule: sum over ALL J of ||rho_J||^2 = ||u||^2 ||Theta u||^2 ----
def rho_all(c):
    t = C.Theta(c)
    tot = 0.0
    for JJ in range(0, 7):
        for Q in range(-JJ, JJ + 1):
            s = 0j
            for m1 in C.MS:
                m2 = Q - m1
                if m2 in C.IDX:
                    s += complex(sp.N(clebsch_gordan(3, 3, JJ, m1, m2, Q), 30)) * c[C.IDX[m1]] * t[C.IDX[m2]]
            tot += abs(s) ** 2
    return tot
u = rand()
P(f"sum rule  sum_J||rho_J||^2/||u||^4 = {rho_all(u)/np.vdot(u,u).real**2:.15f}  (must be 1)")
P(f"   -> rhat6 is the J=6 weight, so 0 <= rhat6 <= 1 identically")

# --- 7. Nbar rotation covariance: Nbar(Ru) = R_SO3 Nbar(u) R_SO3^T ---------
def rodrigues(n, th):
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K
mx = 0.0
for _ in range(20):
    nv, th = rng.normal(size=3), rng.uniform(0, 2 * np.pi)
    R3, D = rodrigues(nv, th), C.rot(nv, th)
    mx = max(mx, np.abs(C.Nbar(D @ u) - R3 @ C.Nbar(u) @ R3.T).max())
P(f"Nbar covariance max err (20 rot)   = {mx:.3e}")
mx = 0.0
for _ in range(20):
    nv, th = rng.normal(size=3), rng.uniform(0, 2 * np.pi)
    R3, D = rodrigues(nv, th), C.rot(nv, th)
    mx = max(mx, np.abs(C.fvec(D @ u) - R3 @ C.fvec(u)).max())
P(f"f-vector covariance max err        = {mx:.3e}")
P(f"Tr Nbar / ||u||^2 (must be 12)     = {np.trace(C.Nbar(u))/np.vdot(u,u).real:.12f}")

# --- 8. known values ---------------------------------------------------------
v3 = C.basis(3)
P(f"Theta v3 (should be v_-3)          = {np.round(C.Theta(v3).real,12)}")
P(f"rhat6(v3) num                      = {C.rhat6(v3):.18f}")
P(f"rhat6(v3) exact                    = {C.rhat6_exact([1,0,0,0,0,0,0])}")
cat = [1, 0, 0, 0, 0, 0, 1]
P(f"rhat6((v3+v-3)/sqrt2) exact        = {C.rhat6_exact(cat)}")
P(f"   as float                        = {sp.N(C.rhat6_exact(cat),25)}")

txt = "\n".join(str(x) for x in out)
print(txt)
pathlib.Path(__file__).parent.joinpath("out_s0.txt").write_text(txt + "\n")
