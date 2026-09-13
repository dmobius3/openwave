"""Independent float64 cross-check of items 3, 4, 6, 7 by exact-degree quadrature on SU(2) (no Clebsch-Gordan used).
Euler angles U = Rz(al) Ry(be) Rz(ga), Rz(a) = diag(e^{-ia/2}, e^{ia/2}), Ry(b) = [[cos b/2, -sin b/2],[sin b/2, cos b/2]].
D^J_{MK}(U) = e^{-iM al} d^J_{MK}(be) e^{-iK ga}. For integer-J integrands the normalized Haar average is
  mean over al, ga in [0, 2pi) (uniform grid) and (1/2) int_0^pi sin(be) d be (Gauss-Legendre in cos be).
Grid exactness: alpha/gamma frequencies <= 18 need > 18 points (40 used); the be-integrand is a polynomial of degree <= 18 in cos be
(24 Gauss nodes are exact to degree 47)."""
import json
import math
import numpy as np
from common import HERE, generate_gamma, q5float, quat_to_su2, Dj_generic

fact = lambda k: float(math.factorial(k))


def Dnum(j, U):
    return np.array(Dj_generic(2 * j, U, fact, math.sqrt, 0j), dtype=complex)


G = generate_gamma()
D3 = [Dnum(3, quat_to_su2([q5float(c) for c in q], 1.0, 1j)) for q in G]
rng = np.random.default_rng(7)
H = rng.normal(size=(7, 7)) + 1j * rng.normal(size=(7, 7))
H = H + H.conj().T
A = sum(D @ H @ D.conj().T for D in D3) / len(G)
E, Q = np.linalg.eigh(A)
order = np.argsort(E)
E, Q = E[order], Q[:, order]
gaps = np.diff(E)
cut = int(np.argmax(gaps)) + 1
ETA = {cut: Q[:, :cut], 7 - cut: Q[:, cut:]}
assert sorted(ETA) == [3, 4], E
for d, eta in ETA.items():
    res = max(np.abs(D @ eta - eta @ (eta.conj().T @ D @ eta)).max() for D in D3)
    print(f"sector {d}: invariance residual {res:.2e}")

Na, Ng, Nb = 40, 40, 24
al = 2 * np.pi * np.arange(Na) / Na
ga = 2 * np.pi * np.arange(Ng) / Ng
xb, wb = np.polynomial.legendre.leggauss(Nb)
be = np.arccos(xb)
wb = wb / 2.0
dJ = {}
for J in range(0, 10):
    dJ[J] = np.array([Dnum(J, [[math.cos(b / 2), -math.sin(b / 2)], [math.sin(b / 2), math.cos(b / 2)]]).real for b in be])  # (Nb, 2J+1, 2J+1)
ph = lambda J, ang: np.exp(-1j * np.outer(ang, np.arange(-J, J + 1)))   # (N, 2J+1): e^{-i m ang}
W = wb[None, :, None] / (Na * Ng)     # weights on (al, be, ga)


def synth(coef, J):
    """function values f_a(al,be,ga) from coefficient array c[M,K,a] at level J."""
    return np.einsum("im,bmk,gk,mka->ibga", ph(J, al), dJ[J], ph(J, ga), coef)


def analyze(f, J):
    """c^J_{MKa} = (2J+1) int conj(D^J_{MK}) f_a."""
    return (2 * J + 1) * np.einsum("ibga,im,bmk,gk->mka", f * W[..., None], ph(J, al).conj(), dJ[J], ph(J, ga).conj())


def integ(f):
    return np.sum(f * W)


mu = {J: 4 * J * (J + 1) - 48 for J in range(10)}
e = lambda m: np.eye(7)[m + 3].astype(complex)
c13, s12 = math.sqrt(13) / 5, 2 * math.sqrt(3) / 5
RAYS = {"R1": e(3), "R2": e(0), "R3": e(2) + e(-2), "R4": e(3) + e(-3), "R5": c13 * e(2) + s12 * e(-3)}
out = {}
for d, eta in sorted(ETA.items()):
    for name, u0 in RAYS.items():
        u = u0 * math.sqrt(7 / (d * np.vdot(u0, u0).real))
        coef = np.einsum("m,ka->mka", u, eta)
        Phi = synth(coef, 3)
        s = np.sum(np.abs(Phi) ** 2, axis=-1)
        N = s[..., None] * Phi
        nPhi = integ(s).real
        kappa = integ(np.sum(Phi.conj() * N, axis=-1))
        cN = {J: analyze(N, J) for J in range(10)}
        Nn = {J: np.sum(np.abs(cN[J]) ** 2) / (2 * J + 1) for J in cN}
        par = np.sqrt(np.sum(np.abs(cN[3] - kappa * coef) ** 2) / 7)
        xi = sum(synth(-cN[J] / mu[J], J) for J in cN if J != 3)
        xin = {2 * J: Nn[J] / mu[J] ** 2 for J in cN if J != 3}
        dot = np.sum(Phi.conj() * xi, axis=-1)
        DN = s[..., None] * xi + (dot + dot.conj())[..., None] * Phi
        along = integ(np.sum(Phi.conj() * DN, axis=-1))
        c3 = analyze(DN, 3)
        orth = np.sqrt(np.sum(np.abs(c3 - along * coef) ** 2) / 7)
        lam4f = -3 * sum(Nn[J] / mu[J] for J in Nn if J != 3)
        out[f"{d}:{name}"] = {"norm_Phi": nPhi, "kappa": kappa.real, "kappa_im": kappa.imag, "parallel_resid": par,
                              "xi_level_norm2": {str(k): v for k, v in xin.items()},
                              "DN_along": along.real, "DN_along_im": along.imag, "DN_orth": orth, "lam4_formula": lam4f}
        print(f"[d={d} {name}] |Phi|^2={nPhi:.15f} kappa={kappa.real:.15f} par={par:.1e} along={along.real:.15f} "
              f"orth={orth:.1e} -3sum={lam4f:.15f}")
        print("    ||Pi_n xi||^2/g^2:", {k: f"{v:.15e}" for k, v in xin.items()})
(HERE / "out_quadcheck.json").write_text(json.dumps(out, indent=1))
