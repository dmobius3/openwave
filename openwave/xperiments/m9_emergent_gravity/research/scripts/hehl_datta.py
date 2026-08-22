#!/usr/bin/env python3
"""Einstein-Cartan + Dirac: extract the Hehl-Datta contact coefficient.

Equations (mostly-minus signature unless a Signature object says otherwise)
--------------------------------------------------------------------
eta = diag(+1, -1, -1, -1) or diag(-1, +1, +1, +1)
eps^{0123} = +1
gamma_ab = (1/2) [gamma_a, gamma_b]
gamma^5 (mostly minus) = i gamma^0 gamma^1 gamma^2 gamma^3
  (mostly plus: the unique matrix with (gamma^5)^2 = +I that anticommutes
   with every gamma^mu, oriented so eps_{0123} gamma^0 gamma^1 gamma^2 gamma^3
   matches the same Hodge dual used for s)

Hermitian Dirac:
    L_D = (i/2) e psibar gamma^mu <->D_mu psi - e m psibar psi
    D_mu = d_mu + (1/4) omega_mu^{ab} gamma_ab

omega-linear piece:
    L_D[omega] = (i/8) e omega_mu^{ab} psibar {gamma^mu, gamma_ab} psi

Canonical spin (this module's definition):
    s^{mu}_{ab} := (i/4) psibar {gamma^mu, gamma_ab} psi
    so that L_D[omega] = (1/2) e s^{mu}_{ab} omega_mu^{ab}

Palatini:
    S_g = (1/(4 kappa)) int eps_{abcd} e^a /\\ e^b /\\ R^{cd}
        = int e R / (2 kappa)
    R(omega) = R(LC) + div(K) + Q(K)
    Q = eta^{sigma nu} (K^mu_{lambda mu} K^lambda_{sigma nu}
                        - K^mu_{lambda nu} K^lambda_{sigma mu})

Totally antisymmetric K is parameterized as K_{abc} = eps_{abcr} v^r.
The on-shell v is the stationary point of L_g[Q] + L_D[K], found from
the finite-difference Hessian (the action is exactly quadratic in K).
This module returns L_on-shell / (-kappa J5·J5). It does not hard-code 3/16.

Clifford algebra is constructed explicitly. No CAS, no literature lookup.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

I2 = np.eye(2, dtype=np.complex128)
SX = np.array([[0, 1], [1, 0]], dtype=np.complex128)
SY = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
SZ = np.array([[1, 0], [0, -1]], dtype=np.complex128)
PAULI = (SX, SY, SZ)


def _block(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    return np.block([[a, b], [c, d]])


@dataclass(frozen=True)
class Signature:
    """Metric signature and the matching Dirac representation."""

    name: str
    eta: np.ndarray
    gamma: np.ndarray  # gamma[mu] is 4x4, raised index
    gamma5: np.ndarray

    @property
    def mostly_minus(self) -> bool:
        return bool(self.eta[0, 0] > 0)


def signature_mostly_minus() -> Signature:
    """Bjorken-Drell / mostly-minus Dirac representation."""
    eta = np.diag([1.0, -1.0, -1.0, -1.0])
    g0 = _block(I2, np.zeros((2, 2)), np.zeros((2, 2)), -I2)
    gi = [_block(np.zeros((2, 2)), s, -s, np.zeros((2, 2))) for s in PAULI]
    gamma = np.stack([g0, *gi], axis=0)
    # gamma5 = i gamma0 gamma1 gamma2 gamma3
    g5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
    return Signature("mostly_minus", eta, gamma, g5)


def signature_mostly_plus() -> Signature:
    """Mostly-plus Dirac matrices with a real vector current.

    gamma^0 = [[0, I], [-I, 0]] (anti-Hermitian),
    gamma^i = [[0, sigma^i], [sigma^i, 0]] (Hermitian),
    gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3.
    Adjoint matrix is A = -gamma^0 (see adjoint()).
    """
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    g0 = _block(np.zeros((2, 2)), I2, -I2, np.zeros((2, 2)))
    gi = [_block(np.zeros((2, 2)), s, s, np.zeros((2, 2))) for s in PAULI]
    gamma = np.stack([g0, *gi], axis=0)
    g5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
    return Signature("mostly_plus", eta, gamma, g5)


def levi_civita() -> np.ndarray:
    """eps[a,b,c,d] with eps[0,1,2,3] = +1 (all-lowered, Minkowski)."""
    eps = np.zeros((4, 4, 4, 4), dtype=np.float64)
    import itertools

    for perm in itertools.permutations(range(4)):
        sign = 1
        seq = list(perm)
        for i in range(4):
            for j in range(i + 1, 4):
                if seq[i] > seq[j]:
                    sign *= -1
        eps[perm] = float(sign)
    return eps


def lower_gamma(sig: Signature) -> np.ndarray:
    """gamma_mu = eta_{mu nu} gamma^nu."""
    return np.einsum("mn,nij->mij", sig.eta, sig.gamma)


def gamma_ab_lower(sig: Signature) -> np.ndarray:
    """gamma_{ab} = (1/2)[gamma_a, gamma_b], both indices lowered."""
    gl = lower_gamma(sig)
    return 0.5 * (np.einsum("aij,bjk->abik", gl, gl) - np.einsum("bij,ajk->abik", gl, gl))


def anticommutator_gamma_mu_gamma_ab(sig: Signature) -> np.ndarray:
    """{gamma^mu, gamma_{ab}} as an array [mu, a, b, i, j]."""
    gab = gamma_ab_lower(sig)
    gmu = sig.gamma
    left = np.einsum("mij,abjk->mabik", gmu, gab)
    right = np.einsum("abij,mjk->mabik", gab, gmu)
    return left + right


def check_clifford(sig: Signature, atol: float = 1e-12) -> dict[str, float]:
    """Return max residuals of the Clifford relations used downstream."""
    anticom = np.einsum("mij,njk->mnik", sig.gamma, sig.gamma) + np.einsum(
        "nij,mjk->mnik", sig.gamma, sig.gamma
    )
    target = 2.0 * np.einsum("mn,ij->mnij", sig.eta, np.eye(4))
    cliff = float(np.max(np.abs(anticom - target)))
    g5sq = float(np.max(np.abs(sig.gamma5 @ sig.gamma5 - np.eye(4))))
    anticom5 = [
        float(np.max(np.abs(sig.gamma5 @ sig.gamma[mu] + sig.gamma[mu] @ sig.gamma5)))
        for mu in range(4)
    ]
    return {
        "clifford_max": cliff,
        "gamma5_square_minus_I": g5sq,
        "gamma5_anticommutator_max": max(anticom5),
    }


def random_spinor(rng: np.random.Generator) -> np.ndarray:
    """A generic 4-component complex spinor, not normalized."""
    return rng.normal(size=4) + 1j * rng.normal(size=4)


def adjoint(sig: Signature, psi: np.ndarray) -> np.ndarray:
    """psibar = psi^dagger A.

    Mostly minus: A = gamma^0 (Hermitian). Mostly plus: A = -gamma^0
    (anti-Hermitian). Both choices make J^mu = psibar gamma^mu psi real
    and make A gamma^mu Hermitian.
    """
    a = sig.gamma[0] if sig.mostly_minus else -sig.gamma[0]
    return psi.conj() @ a


def axial_current(sig: Signature, psi: np.ndarray) -> np.ndarray:
    """J5^mu = psibar gamma^5 gamma^mu psi (raised)."""
    pbar = adjoint(sig, psi)
    return np.array(
        [np.real(pbar @ (sig.gamma5 @ sig.gamma[mu]) @ psi) for mu in range(4)],
        dtype=np.float64,
    )


def spin_tensor_from_dirac(sig: Signature, psi: np.ndarray) -> np.ndarray:
    """s^{mu}_{ab} = (i/4) psibar {gamma^mu, gamma_ab} psi.

    Returns a real array shaped (4, 4, 4) with indices (mu, a, b).
    Antisymmetric in (a, b) by construction of gamma_ab.
    """
    pbar = adjoint(sig, psi)
    anticom = anticommutator_gamma_mu_gamma_ab(sig)
    # (i/4) pbar_i anticom[mu,a,b,i,j] psi_j
    sandwiched = np.einsum("i,mabij,j->mab", pbar, anticom, psi)
    s = (1j / 4.0) * sandwiched
    if np.max(np.abs(s.imag)) > 1e-10:
        raise RuntimeError(
            f"spin tensor not real: max Im = {np.max(np.abs(s.imag)):.3e}"
        )
    return np.real(s)


def raise_last_two(sig: Signature, s_m_ab: np.ndarray) -> np.ndarray:
    """s^{mu ab} from s^{mu}_{ab}."""
    return np.einsum("mab,ac,bd->mcd", s_m_ab, sig.eta, sig.eta)


def all_raised_spin(sig: Signature, s_m_ab: np.ndarray) -> np.ndarray:
    """s^{lambda mu nu} with all indices raised, identifying (lambda, mu, nu) = (mu_form, a, b)."""
    return raise_last_two(sig, s_m_ab)


def dual_of_vector(sig: Signature, j_raised: np.ndarray) -> np.ndarray:
    """(*J)^{lambda mu nu} = eps^{lambda mu nu rho} J_rho, all raised.

    eps^{lambda mu nu rho} = det(eta) * eps_{lambda mu nu rho} in an orthonormal
    frame? In Minkowski with orthonormal indices, raised eps carries
    sign(det eta) = -1 in 4D Lorentzian signature either way? Careful:

    eps^{0123} = eta^{00} eta^{11} eta^{22} eta^{33} eps_{0123} = det(eta) eps_{0123}.
    det(eta) = -1 in both signatures, eps_{0123} = +1, so eps^{0123} = -1.

    Then (*J)^{012} = eps^{012 rho} J_rho = eps^{0123} J_3 = - J_3.
    """
    eps_low = levi_civita()
    det_eta = float(np.linalg.det(sig.eta))
    eps_high = det_eta * eps_low  # orthonormal frame
    j_low = sig.eta @ j_raised
    return np.einsum("lmnr,r->lmn", eps_high, j_low)


def claimed_spin_from_j5(sig: Signature, j5: np.ndarray) -> np.ndarray:
    """-(1/4) (*J5)^{lambda mu nu}, all raised."""
    return -0.25 * dual_of_vector(sig, j5)


def palatini_quadratic(sig: Signature, k_low: np.ndarray) -> float:
    """Q = eta^{sigma nu} (K^mu_{lambda mu} K^lambda_{sigma nu}
                         - K^mu_{lambda nu} K^lambda_{sigma mu}).

    k_low[a,b,c] = K_{a b c} with all indices lowered.
    K^mu_{lambda nu} = eta^{mu alpha} K_{alpha lambda nu}.
    Returns the scalar Q (not yet divided by 2 kappa).
    """
    eta_inv = sig.eta  # orthonormal: eta^{mu nu} = eta_{mu nu} numerically
    # K^mu_{lambda nu} = eta^{mu a} K_{a lambda nu}
    k_up_ll = np.einsum("ma,aln->mln", eta_inv, k_low)
    k_trace = np.einsum("mlm->l", k_up_ll)
    # K^lambda_{sigma nu} = eta^{lambda b} K_{b sigma nu}
    k_lam_sig_nu = np.einsum("lb,bsn->lsn", eta_inv, k_low)
    k_lam_tr = np.einsum("sn,lsn->l", eta_inv, k_lam_sig_nu)
    term1 = float(np.dot(k_trace, k_lam_tr))
    # eta^{sigma nu} K^mu_{lambda nu} K^lambda_{sigma mu}
    term2 = float(np.einsum("sn,mln,lsm->", eta_inv, k_up_ll, k_lam_sig_nu))
    return term1 - term2


def palatini_quadratic_alt(sig: Signature, k_low: np.ndarray) -> float:
    """Same Q as palatini_quadratic, explicit index sums.

    Independent of the einsum path. The previous contraction
    K^{μcν}K_{cμν} had the second factor's first two indices
    swapped and returned −Q for totally antisymmetric K.
    """
    eta = sig.eta
    kup = np.einsum("ma,aln->mln", eta, k_low)
    q = 0.0
    for si in range(4):
        for nu in range(4):
            t1 = 0.0
            t2 = 0.0
            for la in range(4):
                tr = (
                    kup[0, la, 0]
                    + kup[1, la, 1]
                    + kup[2, la, 2]
                    + kup[3, la, 3]
                )
                t1 += tr * kup[la, si, nu]
                t2 += (
                    kup[0, la, nu] * kup[la, si, 0]
                    + kup[1, la, nu] * kup[la, si, 1]
                    + kup[2, la, nu] * kup[la, si, 2]
                    + kup[3, la, nu] * kup[la, si, 3]
                )
            q += eta[si, nu] * (t1 - t2)
    return float(q)


def k_from_vector(v: np.ndarray) -> np.ndarray:
    """K_{abc} = eps_{abcr} v^r with eps_{0123} = +1 (component contraction)."""
    return np.einsum("abcr,r->abc", levi_civita(), v)


def lagrangian_of_k(
    sig: Signature, s_m_ab: np.ndarray, k_low: np.ndarray, kappa: float
) -> tuple[float, float, float, float]:
    """Return (L_int, L_g, L_D, Q) for a given contorsion."""
    q = palatini_quadratic(sig, k_low)
    l_g = q / (2.0 * kappa)
    k_mu_ab = np.einsum("ap,bq,mpq->mab", sig.eta, sig.eta, k_low)
    l_d = 0.5 * float(np.einsum("mab,mab->", s_m_ab, k_mu_ab))
    return l_g + l_d, l_g, l_d, q


def stationary_vector(sig: Signature, s_m_ab: np.ndarray, kappa: float) -> np.ndarray:
    """Solve dL/dv = 0 for the 4-vector that parameterizes totally antisymmetric K.

    Finite-difference Hessian at v=0 (the action is exactly quadratic in K).
    """
    step = 1e-6
    v0 = np.zeros(4, dtype=np.float64)

    def l_only(v: np.ndarray) -> float:
        return lagrangian_of_k(sig, s_m_ab, k_from_vector(v), kappa)[0]

    grad = np.zeros(4, dtype=np.float64)
    for i in range(4):
        vp = v0.copy()
        vp[i] += step
        vm = v0.copy()
        vm[i] -= step
        grad[i] = (l_only(vp) - l_only(vm)) / (2.0 * step)
    hess = np.zeros((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            vpp = v0.copy()
            vpp[i] += step
            vpp[j] += step
            vpm = v0.copy()
            vpm[i] += step
            vpm[j] -= step
            vmp = v0.copy()
            vmp[i] -= step
            vmp[j] += step
            vmm = v0.copy()
            vmm[i] -= step
            vmm[j] -= step
            hess[i, j] = (l_only(vpp) - l_only(vpm) - l_only(vmp) + l_only(vmm)) / (
                4.0 * step * step
            )
    return -np.linalg.solve(hess, grad)


def dual_factor(s_high: np.ndarray, dual: np.ndarray) -> float:
    """Least-squares factor a in s = a * (*J)."""
    num = float(np.vdot(dual, s_high))
    den = float(np.vdot(dual, dual))
    return num / den if abs(den) > 1e-30 else float("nan")


def contact_lagrangian(sig: Signature, psi: np.ndarray, kappa: float = 1.0) -> dict:
    """On-shell L_int from stationary K, and L_int / (-kappa J5·J5)."""
    j5 = axial_current(sig, psi)
    s_m_ab = spin_tensor_from_dirac(sig, psi)
    s_high = all_raised_spin(sig, s_m_ab)
    claimed = claimed_spin_from_j5(sig, j5)
    dual = dual_of_vector(sig, j5)

    v_star = stationary_vector(sig, s_m_ab, kappa)
    k_low = k_from_vector(v_star)
    l_int, l_g, l_d, q = lagrangian_of_k(sig, s_m_ab, k_low, kappa)
    q_alt = palatini_quadratic_alt(sig, k_low)

    jsq = float(j5 @ sig.eta @ j5)
    ratio = l_int / (-kappa * jsq) if abs(jsq) > 1e-30 else float("nan")

    spin_diff = s_high - claimed
    spin_res = float(np.max(np.abs(spin_diff))) / max(float(np.max(np.abs(j5))), 1e-30)
    factor = dual_factor(s_high, dual)

    anti = (
        s_high
        - s_high.transpose(0, 2, 1)
        - s_high.transpose(1, 0, 2)
        + s_high.transpose(1, 2, 0)
        + s_high.transpose(2, 0, 1)
        - s_high.transpose(2, 1, 0)
    ) / 6.0
    mixed_res = float(np.max(np.abs(s_high - anti))) / max(
        float(np.max(np.abs(s_high))), 1e-30
    )
    v_over_j = v_star / np.where(np.abs(j5) > 1e-12, j5, np.nan)

    return {
        "j5": j5,
        "j5_sq": jsq,
        "s_high": s_high,
        "claimed_s": claimed,
        "spin_residual": spin_res,
        "spin_dual_factor": factor,
        "v_star": v_star,
        "v_over_j5": v_over_j,
        "mixed_residual": mixed_res,
        "Q": q,
        "Q_alt": q_alt,
        "Q_alt_residual": abs(q - q_alt),
        "L_g": l_g,
        "L_D": l_d,
        "L_int": l_int,
        "ratio": ratio,
        "kappa": kappa,
    }


def mutated_ratio(
    sig: Signature,
    psi: np.ndarray,
    *,
    palatini_factor: float = 0.5,
    source: str = "dirac",
    kappa: float = 1.0,
) -> float:
    """Recompute the C1 ratio with a mutated Palatini factor or a fake source.

    palatini_factor: coefficient of Q/kappa in L_g (canonical 1/2).
    source: "dirac" (canonical) or "plus_quarter_dual" (s = +1/4 *J5).
    Stationarity is re-solved so the mutation is not a spectator.
    """
    j5 = axial_current(sig, psi)
    if source == "dirac":
        s_m_ab = spin_tensor_from_dirac(sig, psi)
    elif source == "plus_quarter_dual":
        s_high = 0.25 * dual_of_vector(sig, j5)
        s_m_ab = np.einsum("mcd,ca,db->mab", s_high, sig.eta, sig.eta)
    else:
        raise ValueError(source)

    def l_only(v: np.ndarray) -> float:
        k_low = k_from_vector(v)
        q = palatini_quadratic(sig, k_low)
        l_g = palatini_factor * q / kappa
        k_mu_ab = np.einsum("ap,bq,mpq->mab", sig.eta, sig.eta, k_low)
        l_d = 0.5 * float(np.einsum("mab,mab->", s_m_ab, k_mu_ab))
        return l_g + l_d

    step = 1e-6
    v0 = np.zeros(4, dtype=np.float64)
    grad = np.zeros(4, dtype=np.float64)
    for i in range(4):
        vp = v0.copy()
        vp[i] += step
        vm = v0.copy()
        vm[i] -= step
        grad[i] = (l_only(vp) - l_only(vm)) / (2.0 * step)
    hess = np.zeros((4, 4), dtype=np.float64)
    for i in range(4):
        for j in range(4):
            vpp = v0.copy()
            vpp[i] += step
            vpp[j] += step
            vpm = v0.copy()
            vpm[i] += step
            vpm[j] -= step
            vmp = v0.copy()
            vmp[i] -= step
            vmp[j] += step
            vmm = v0.copy()
            vmm[i] -= step
            vmm[j] -= step
            hess[i, j] = (l_only(vpp) - l_only(vpm) - l_only(vmp) + l_only(vmm)) / (
                4.0 * step * step
            )
    v_star = -np.linalg.solve(hess, grad)
    l_int = l_only(v_star)
    jsq = float(j5 @ sig.eta @ j5)
    return l_int / (-kappa * jsq)
