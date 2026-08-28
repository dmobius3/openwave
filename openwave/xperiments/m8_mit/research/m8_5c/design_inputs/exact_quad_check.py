"""Verify F3/P2: for band-limited psi on S^3, every Galerkin coefficient of |psi|^2 psi
is a polynomial of degree <= 4N, so a Hopf product rule exact to that degree returns the
orthogonal projection to rounding. Mutation: a rule exact only to 2N must fail.
Harmonics are the repo's own: vec(sym_power(quat_to_su2(x), n)), degree-n polynomials in x.
"""
import numpy as np, sys, time
sys.path.insert(0, "openwave/xperiments/m8_mit/research/m8_5b/pilot")
from route_a_nonabelian import quat_to_su2, sym_power

N = 3                                   # band limit of psi
rng = np.random.default_rng(20260826)

def hopf_rule(D):
    """Product rule exact for polynomials of degree <= D on S^3 (normalized measure).
    x = (cos(eta) cos(xi1), cos(eta) sin(xi1), sin(eta) cos(xi2), sin(eta) sin(xi2)),
    measure prop. to cos(eta) sin(eta) d(eta) d(xi1) d(xi2); u = sin^2(eta) linearizes it."""
    K = D + 1                            # Fourier modes to order D in each angle
    nu = D // 2 + 1                      # u-polynomials to degree D/2, GL exact to 2*nu-1
    xs, ws = np.polynomial.legendre.leggauss(nu)
    u = (xs + 1) / 2; wu = ws / 2
    xi = 2 * np.pi * np.arange(K) / K
    ce, se = np.sqrt(1 - u), np.sqrt(u)
    X, W = [], []
    for cu, su, w in zip(ce, se, wu):
        for a in xi:
            for b in xi:
                X.append([cu*np.cos(a), cu*np.sin(a), su*np.cos(b), su*np.sin(b)])
                W.append(w)
    X = np.array(X); W = np.array(W); W /= W.sum()
    return X, W

def eval_modes(X, M):
    """Columns: all matrix-element harmonics D^n_ij, n = 0..M."""
    cols = []
    for x in X:
        U = quat_to_su2(np.asarray(x, dtype=float))
        v = [np.array([1.0+0j])]
        for n in range(1, M+1):
            v.append(sym_power(U, n).reshape(-1))
        cols.append(np.concatenate(v))
    return np.array(cols)

def project(X, W, M, psi_vals):
    """Galerkin coefficients of |psi|^2 psi against all modes up to M, via the rule."""
    Y = eval_modes(X, M)
    g = (np.abs(psi_vals)**2) * psi_vals
    return (Y.conj() * (W * g)[:, None]).sum(axis=0)

nmodes = sum((n+1)**2 for n in range(N+1))
coef = rng.standard_normal(nmodes) + 1j*rng.standard_normal(nmodes)

t0 = time.time()
X1, W1 = hopf_rule(4*N)                  # exact for the cubic Galerkin integrand
X2, W2 = hopf_rule(8*N)                  # oversampled reference
X3, W3 = hopf_rule(2*N)                  # MUTATION: exact only to 2N, must alias

psi1 = eval_modes(X1, N) @ coef
psi2 = eval_modes(X2, N) @ coef
psi3 = eval_modes(X3, N) @ coef

print(f"  N={N}, modes={nmodes};  nodes: exact-4N {len(W1)}, reference {len(W2)}, mutated {len(W3)}")

# side check: the rule reproduces Schur orthogonality (degree-2N integrands)
Y1 = eval_modes(X1, N)
G = Y1.conj().T @ (W1[:, None] * Y1)
Gex = np.zeros_like(G); k = 0
for n in range(N+1):
    d = (n+1)**2
    Gex[k:k+d, k:k+d] = np.eye(d) / (n+1)
    k += d
print(f"  Gram vs Schur orthogonality (exact-4N rule): max err {np.abs(G-Gex).max():.2e}")

F1 = project(X1, W1, N, psi1)
F2 = project(X2, W2, N, psi2)
F3 = project(X3, W3, N, psi3)
scale = np.abs(F2).max()
e_exact = np.abs(F1 - F2).max() / scale
e_mut   = np.abs(F3 - F2).max() / scale
print(f"  cubic Galerkin coefficients, exact-4N vs reference : max rel {e_exact:.2e}")
print(f"  MUTATION, rule exact only to 2N vs reference       : max rel {e_mut:.2e}")
print(f"  wall time {time.time()-t0:.1f}s")
print()
ok = e_exact < 1e-11 and e_mut > 1e3 * max(e_exact, 1e-16)
print(f"  VERDICT: {'CONFIRMED' if ok else 'NOT CONFIRMED'}: degree-4N quadrature returns the "
      f"projection to rounding, and under-resolution is DETECTED, not silent")

# The Gram mismatch above is against NAIVE Schur normalization. S1b already recorded that
# sym_power is NOT unitary in this basis (its Reynolds projector was oblique for that reason).
# The right check: the RULE's Gram vs the reference rule's Gram, and block-diagonality.
Y2 = eval_modes(X2, N)
G2 = Y2.conj().T @ (W2[:, None] * Y2)
print(f"  Gram, exact-4N rule vs oversampled reference: max diff {np.abs(G-G2).max():.2e}")
off = 0.0; k = 0
blocks = []
for n in range(N+1):
    d = (n+1)**2; blocks.append((k, k+d)); k += d
M = np.abs(G2).copy()
for (a,b) in blocks: M[a:b, a:b] = 0
print(f"  cross-level Gram entries (must vanish, levels are orthogonal): max {M.max():.2e}")
print(f"  within-level Gram is NOT the identity, consistent with S1b's finding that")
print(f"  sym_power is non-unitary; the basis is a non-orthonormal spanning set per level.")
