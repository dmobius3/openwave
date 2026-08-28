"""Structure of the fluctuation operator for R(phi;om) = c^2 D phi - c1<phi,phi>phi + om^2 phi.

Two claims checked on a faithful finite stand-in (D = diag of Laplacian eigenvalues on a
complex level-n eigenspace; the nonlinearity is the exact one from the pinned wave_engine v_mode 1):

  C1  phi -> <phi,phi>phi is NOT complex-differentiable, so the fluctuation operator is
      real-linear on R^{2m}, not a complex m x m matrix. A complex Jacobian drops a term.
  C2  i*phi is an exact kernel vector of the true (real) Jacobian at every solution: the U(1)
      gauge mode of the standing-wave ansatz. It sits INSIDE the scored eigenspace.
"""
import numpy as np, sys
rng = np.random.default_rng(0)

c2, c1 = 1.0, 0.7
m = 6                      # complex dim of H_{rho,d_rho}: d_rho + 1, here d_rho = 5
lam = -3.0                 # Laplacian eigenvalue -n(n+2)/R^2 on this level (degenerate)

def R(phi, om2):
    return c2*lam*phi - c1*np.vdot(phi, phi).real*phi + om2*phi

# a genuine solution: on a single degenerate level, <phi,phi> = (om2 + c2*lam)/c1
om2 = 4.0
nrm2 = (om2 + c2*lam)/c1
v = rng.normal(size=m) + 1j*rng.normal(size=m)
phi = v/np.linalg.norm(v)*np.sqrt(nrm2)
print("residual at the constructed solution:", np.linalg.norm(R(phi, om2)))

def real_vec(z): return np.concatenate([z.real, z.imag])
def cplx_vec(x): return x[:m] + 1j*x[m:]

# true real-linear Jacobian, by complex-step-free finite differences on R^{2m}
J = np.zeros((2*m, 2*m))
h = 1e-7
for k in range(2*m):
    e = np.zeros(2*m); e[k] = h
    J[:, k] = (real_vec(R(cplx_vec(real_vec(phi)+e), om2))
               - real_vec(R(cplx_vec(real_vec(phi)-e), om2)))/(2*h)

# naive complex-linear Jacobian: treats <phi,phi> as a frozen scalar (drops the <dphi,phi>phi term)
Jc = (c2*lam - c1*np.vdot(phi, phi).real + om2)*np.eye(m, dtype=complex)
Jc_as_real = np.block([[Jc.real, -Jc.imag], [Jc.imag, Jc.real]])

print("\nC1  ||J_true - J_naive||_2 / ||J_true||_2 =",
      f"{np.linalg.norm(J-Jc_as_real,2)/np.linalg.norm(J,2):.3f}",
      "  (0 would mean the naive complex Jacobian is right)")
print("    eigenvalues, true real J :", np.round(np.sort(np.linalg.eigvals(J).real), 6))
print("    eigenvalues, naive       :", np.round(np.sort(np.linalg.eigvals(Jc_as_real).real), 6))

gauge = real_vec(1j*phi)
kern = np.linalg.norm(J@gauge)/np.linalg.norm(gauge)
kern_ok = kern < 1e-6                       # the label is COMPUTED, not printed unconditionally
print("\nC2  ||J @ (i*phi)|| / ||i*phi|| =", f"{kern:.3e}",
      " -> exact kernel vector" if kern_ok else " -> NOT a kernel vector, C2 FAILS")
print("    count of near-zero eigenvalues of J:", int(np.sum(np.abs(np.linalg.eigvals(J)) < 1e-6)))
print("    i*phi lies inside the scored eigenspace H_{rho,d_rho}: True by construction",
      "(complex vector space, same level, same bundle)")

# mutation arm: break the gauge symmetry, the zero mode must disappear
def R_broken(phi, om2):
    return c2*lam*phi - c1*(np.vdot(phi, phi).real + 0.3*phi[0].real)*phi + om2*phi
Jb = np.zeros((2*m, 2*m))
for k in range(2*m):
    e = np.zeros(2*m); e[k] = h
    Jb[:, k] = (real_vec(R_broken(cplx_vec(real_vec(phi)+e), om2))
                - real_vec(R_broken(cplx_vec(real_vec(phi)-e), om2)))/(2*h)
mut = np.linalg.norm(Jb@gauge)/np.linalg.norm(gauge)
mut_ok = mut > 1e3 * max(kern, 1e-16)       # red RELATIVE to the parent, per the house pattern
print("\n    mutation (gauge-breaking term): ||J_b @ (i*phi)||/||i*phi|| =",
      f"{mut:.3e}", "-> arm goes red" if mut_ok else "-> arm did NOT fire (broken arm)")
print(f"\n  VERDICT: {'CONFIRMED' if (kern_ok and mut_ok) else 'NOT CONFIRMED'}: i*phi is an exact "
      f"kernel vector of the real Jacobian and the gauge-breaking mutation destroys it")
sys.exit(0 if (kern_ok and mut_ok) else 1)
