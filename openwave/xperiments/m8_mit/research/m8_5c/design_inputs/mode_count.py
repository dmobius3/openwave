"""Galerkin mode count per 2I-sector on S^3 as a function of harmonic cutoff N.
Sections of E_rho at level n: (V_n (x) W_rho)^{2I} (x) V_n^*, complex dim m_rho(n)*(n+1).
Self-check: first occurrence must be n = d_rho with multiplicity 1 (M8.4 prereg s.3 table).
"""
import numpy as np
from math import pi, sin, cos

# 2I conjugacy classes: SU(2) rotation angle phi (eigenvalues e^{+-i phi}), class size
classes = [(0.0,1),(pi,1),(pi/2,30),(pi/3,20),(2*pi/3,20),(pi/5,12),(2*pi/5,12),(3*pi/5,12),(4*pi/5,12)]
assert sum(s for _,s in classes) == 120
galois = {pi/5:3*pi/5, 3*pi/5:pi/5, 2*pi/5:4*pi/5, 4*pi/5:2*pi/5}

def chiV(n, phi):
    if abs(phi) < 1e-12: return n+1
    if abs(phi-pi) < 1e-12: return (n+1)*(-1)**n
    return sin((n+1)*phi)/sin(phi)

def R5(phi):  # A5 standard 4-dim rep
    if abs(phi)<1e-12 or abs(phi-pi)<1e-12: return 4
    if abs(phi-pi/2)<1e-12: return 0
    if abs(phi-pi/3)<1e-12 or abs(phi-2*pi/3)<1e-12: return 1
    return -1

g = lambda phi: galois.get(phi, phi)
irreps = {
 'R0': lambda p: 1.0,
 'R1': lambda p: chiV(1,p),
 'R2': lambda p: chiV(1,g(p)),
 'R3': lambda p: chiV(2,p),
 'R4': lambda p: chiV(2,g(p)),
 'R5': R5,
 'R6': lambda p: chiV(3,p),
 'R7': lambda p: chiV(4,p),
 'R8': lambda p: chiV(5,p),
}
dims = {k: int(round(f(0.0))) for k,f in irreps.items()}

# self-check 1: orthonormal characters
names = list(irreps)
G = np.array([[sum(s*irreps[a](p)*irreps[b](p) for p,s in classes)/120 for b in names] for a in names])
assert np.allclose(G, np.eye(9), atol=1e-9), "character table fails orthonormality"

def mult(rho, n):
    return int(round(sum(s*chiV(n,p)*irreps[rho](p) for p,s in classes)/120))

# self-check 2: first occurrence = d_rho with multiplicity 1 (prereg s.3)
expected_d = {'R0':0,'R1':1,'R3':2,'R6':3,'R7':4,'R8':5,'R4':6,'R5':6,'R2':7}
for rho,d in expected_d.items():
    first = next(n for n in range(0,40) if mult(rho,n)>0)
    assert first == d and mult(rho,d) == 1, (rho, first, mult(rho,d))
print("self-checks: characters orthonormal; first occurrence n=d_rho, mult 1, all nine sectors")

Ns = [10,20,30,40,60,80,100]
print("\ncomplex Galerkin modes in sector rho up to cutoff N (sum_n m_rho(n)(n+1)):")
print("sector dim d_rho | " + " ".join(f"N={N:>4}" for N in Ns))
for rho in ['R0','R1','R3','R6','R7','R8','R4','R5','R2']:
    counts = [sum(mult(rho,n)*(n+1) for n in range(0,N+1)) for N in Ns]
    print(f"{rho:>5} {dims[rho]:>3} {expected_d[rho]:>5} | " + " ".join(f"{c:>6}" for c in counts))
print("\nasymptotic: ~ d_rho * N^3 / 360")
