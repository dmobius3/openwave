"""AUDIT A6: boundedness of the v4 Hamiltonian (try to break it) + A7 dimension count."""
import numpy as np
import sympy as sp

print("=== A6a: pointwise lower bound (bounded rho, compact support) ===")
phi_s, m_s, lam_s, b_s = sp.symbols("phi m lam b", real=True)
# per-point potential with worst-case linear tilt b = g*rho_max*|phi| (take phi<0 branch):
f = sp.Rational(1, 2) * m_s ** 2 * phi_s ** 2 + lam_s / 24 * phi_s ** 4 - b_s * phi_s
crit = sp.solve(sp.diff(f, phi_s), phi_s)
print("critical points exist; quartic dominates: lim f(phi->+-oo) =",
      sp.limit(f.subs({m_s: 1, lam_s: 1, b_s: 1}), phi_s, sp.oo),
      sp.limit(f.subs({m_s: 1, lam_s: 1, b_s: 1}), phi_s, -sp.oo))
# explicit finite bound ignoring the (positive) mass term:
# min over phi of (lam/24)phi^4 - b*phi = -(3/4) b (6b/lam)^(1/3)
bmin = -sp.Rational(3, 4) * b_s * (6 * b_s / lam_s) ** sp.Rational(1, 3)
chk = sp.simplify(sp.Min(*[f.subs({m_s: 0}).subs(phi_s, cp) for cp in
                           sp.solve(sp.diff(f.subs({m_s: 0}), phi_s), phi_s) if cp.is_real]) - bmin)
print("analytic pointwise min of (lam/24)phi^4 - b phi equals -(3/4) b (6b/lam)^(1/3):",
      sp.simplify(chk) == 0 or chk)
print("=> H >= -(3/4) g rho_max (6 g rho_max/lam)^(1/3) * Vol(supp rho) > -oo : BOUNDED BELOW")
print("   (pi^2, (grad phi)^2 >= 0 only add; outside supp(rho) the density is >= 0)")

print("\n=== A6b: try to break it -- point/unbounded source (rho -> A delta^3) ===")
# trial phi_h(x) = -h * exp(-(|x|/eps)^2), eps = h^-2, m=lam=g=1, A=1:
# E(h) = grad + mass + quartic - g A h ; show E -> -oo as h -> oo  (assumption necessary)
for h in [1, 10, 100, 1000]:
    eps = h ** -2.0
    rr = np.linspace(1e-9, 12 * eps, 400000)
    ph = -h * np.exp(-((rr / eps) ** 2))
    dph = np.gradient(ph, rr)
    grad = np.trapezoid(0.5 * dph ** 2 * 4 * np.pi * rr ** 2, rr)
    mass = np.trapezoid(0.5 * ph ** 2 * 4 * np.pi * rr ** 2, rr)
    quart = np.trapezoid(ph ** 4 / 24 * 4 * np.pi * rr ** 2, rr)
    E = grad + mass + quart - 1.0 * h  # interaction g*A*phi(0) = -h
    print(f"  h={h:5d}: grad={grad:.3e} mass={mass:.3e} quartic={quart:.3e} "
          f"interaction={-h:.0f}  E_total={E:.3f}")
print("=> for a DELTA-function (unbounded) rho_n the energy is UNBOUNDED below;")
print("   the 'bounded + compactly supported rho_n' assumption is load-bearing.")

print("\n=== A7: dimension counting, 3+1D, hbar=c=1 ===")
print("[L]=4; [(d phi)^2]=4 -> [phi]=1; [m^2 phi^2] -> [m]=1;")
print("[lam phi^4]=4 -> [lam]=4-4*1=0 (DIMENSIONLESS);")
print("[g phi rho_n]: [rho_n]=3 (number density, cf. v4 Fig3 rho0=0.16 fm^-3)")
print(" -> [g]=4-1-3=0 (DIMENSIONLESS).")
print("Superficial divergence of phi^4 in d=4: D = 4 - E (E=external legs);")
print("D>=0 for E=2,4 at EVERY loop order -> infinitely many divergent diagrams")
print("-> STRICTLY renormalizable (textbook phi^4), NOT superrenormalizable")
print("(superrenormalizable needs [coupling]>0, e.g. phi^4 in d=3 or phi^3 in d=4).")
