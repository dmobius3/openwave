"""AUDIT A3/A4/A5/A8/A9: independent arithmetic + integrals (own profiles, own code)."""
import numpy as np

HBARC = 197.3269788  # MeV fm (CODATA)
ALPHA = 1 / 137.035999084
me = 0.511  # MeV, as the paper uses

print("=== A3: Table 2 gap arithmetic ===")
e_nat = np.sqrt(4 * np.pi * ALPHA)
target = me / e_nat
print(f"e_nat = sqrt(4 pi alpha) = {e_nat:.6f}  (paper 0.30282)")
print(f"target H/Q = me/e_nat    = {target:.6f}  (paper 1.6875)")
for a_, b_, lbl in [(1.6969, 1.6875, "Table 2 pairing: 1.6969 vs 1.6875"),
                    (1.6918, 1.6969, "Sec 5.1 pairing: 1.6918 vs 1.6969"),
                    (1.6890, 1.6875, "Step-2 pairing: 1.6890 vs 1.6875")]:
    print(f"{lbl}: |diff|/ref = {abs(a_-b_)/b_*100:.3f}%")

print("\n=== A4: 'Q_predicted' identity ===")
Qpred = me / 1.6969
print(f"me / 1.6969 = {Qpred:.5f}  (paper prints Q_predicted = 0.3011)")
print(f"gap vs e_nat: {(e_nat-Qpred)/e_nat*100:.3f}%  (paper prints 0.56%)")
print(f"H/Q gap restated: 1.6875(target)/1.6969(code) = {1.68746/1.6969:.5f} "
      f"-> {(1-1.68746/1.6969)*100:.3f}%  (same number: Q_pred/e_nat = (me/e_nat)/(H/Q))")

print("\n=== A5: R_phys = 191 fm scale-setting ===")
H_code = me * 191 / HBARC
print(f"H_code implied by 'H_code*(hbarc/191fm)=0.511': {H_code:.5f}")
print(f"hbarc/191 fm = {HBARC/191:.5f} MeV  (claim: equals July-papers m_J = 1.033 MeV)")
for H in [1.6890, 1.6918, 1.6969, 2.1173]:
    print(f"  if H_code were {H}: R_phys = {HBARC*H/me:.1f} fm  (NOT 191)")

print("\n=== A8: Yukawa(429 fm) vs 1/r ===")
lam = HBARC / 0.460
print(f"range hbar_c/0.460 MeV = {lam:.1f} fm")
for r in [1, 5, 10, 40, 100, 400, 429]:
    dev = 1 - np.exp(-r / lam)
    print(f"  r = {r:4d} fm: relative deviation |Yukawa - 1/r|/(1/r) = {dev*100:.2f}%")

print("\n=== A9: E_int = g phi(0) A bookkeeping (own profiles) ===")
# my own choice, independent of the task's script: Woods-Saxon a=0.50 fm AND a hard
# uniform ball, field phi = exp(-r/429 fm) (any long-range profile), ratio is g,A-free
r = np.linspace(1e-4, 30, 300000)
dr = r[1] - r[0]


def ratio(rho):
    phi = np.exp(-r / lam)
    num = np.trapezoid(phi * rho * 4 * np.pi * r ** 2, r)
    den = 1.0 * np.trapezoid(rho * 4 * np.pi * r ** 2, r)  # phi(0)=1
    return num / den


for A in [1, 2, 40, 85, 184, 208]:
    R = 1.2 * A ** (1 / 3)
    ws = 1 / (1 + np.exp((r - R) / 0.50))
    ball = (r <= R).astype(float)
    print(f"  A={A:3d}: E_int/(g phi(0) A)  Woods-Saxon = {ratio(ws):.4f}   uniform ball = {ratio(ball):.4f}")
print("  (a standard per-nucleon Yukawa coupling g*sum_i phi(x_i) = g*int phi rho_n d3x"
      " is the SAME integral -> same A-linear energy)")
