"""Verify Feedback-2 fix 1: the cascade monitor measures the band N < n <= 3N of |psi|^2 psi,
so its integrands reach degree 4N... no: level-a harmonics with a <= 3N against a degree-3N
integrand reach 3N + 3N = 6N. The production rule, exact only to 4N, must therefore alias the
monitor's own band. Demonstrated with a 6N-exact rule vs an oversampled reference vs the 4N rule.
"""
import numpy as np, sys, time
sys.path.insert(0, "openwave/xperiments/m8_mit/research/m8_5b/pilot")
from route_a_nonabelian import quat_to_su2, sym_power

N = 3
rng = np.random.default_rng(20260827)

def hopf_rule(D):
    K = D + 1; nu = D // 2 + 1
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
    cols = []
    for x in X:
        U = quat_to_su2(np.asarray(x, dtype=float))
        v = [np.array([1.0+0j])]
        for n in range(1, M+1):
            v.append(sym_power(U, n).reshape(-1))
        cols.append(np.concatenate(v))
    return np.array(cols)

nlow = sum((n+1)**2 for n in range(N+1))            # psi modes, levels 0..N
coef = rng.standard_normal(nlow) + 1j*rng.standard_normal(nlow)

def band_coeffs(D):
    X, W = hopf_rule(D)
    Y = eval_modes(X, 3*N)                           # all modes to 3N
    psi = Y[:, :nlow] @ coef
    g = (np.abs(psi)**2) * psi
    c = (Y.conj() * (W * g)[:, None]).sum(axis=0)
    return c[nlow:], len(W)                          # the cascade band N < n <= 3N

t0 = time.time()
c_ref, n_ref = band_coeffs(9*N)                      # oversampled reference
c_6N,  n_6N  = band_coeffs(6*N)                      # the monitor's own rule
c_4N,  n_4N  = band_coeffs(4*N)                      # the PRODUCTION rule
scale = np.abs(c_ref).max()
e6 = np.abs(c_6N - c_ref).max() / scale
e4 = np.abs(c_4N - c_ref).max() / scale
print(f"  N={N}; cascade band = levels {N+1}..{3*N}, {len(c_ref)} coefficients")
print(f"  nodes: 6N-exact {n_6N}, production 4N {n_4N}, reference {n_ref};"
      f"  node ratio 6N/4N = {n_6N/n_4N:.2f}")
print(f"  band coefficients, 6N rule vs reference        : max rel {e6:.2e}")
print(f"  band coefficients, PRODUCTION 4N rule vs ref   : max rel {e4:.2e}")
frac_ref = np.linalg.norm(c_ref); frac_4N = np.linalg.norm(c_4N)
print(f"  monitor's own reading (band norm): ref {frac_ref:.6f}, under 4N rule {frac_4N:.6f}"
      f"  ({100*abs(frac_4N-frac_ref)/frac_ref:.1f}% off)")
print(f"  wall {time.time()-t0:.1f}s")
ok = e6 < 1e-11 and e4 > 1e3*max(e6, 1e-16)
print(f"\n  VERDICT: {'CONFIRMED' if ok else 'NOT CONFIRMED'}: the monitor needs its own 6N-exact"
      f" rule; the production 4N rule aliases the very band the monitor exists to measure")
