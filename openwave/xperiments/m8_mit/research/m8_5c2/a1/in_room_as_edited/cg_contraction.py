"""CG contraction route: P_N[|psi|^2 psi] purely in coefficient space.

Second (transform-free) route required by § 4.3. Products of representation
matrix elements are expanded via SU(2) Clebsch-Gordan coefficients; the quartic
coupling is applied per-level without going to physical space.

Key identity: the unitarized representation pi^u_n[a,b] equals the Wigner
D-matrix D^{n/2}_{n/2-a, n/2-b} after the index map m = j - a. So the CG
expansion of pi_bar^u * pi^u is the standard D-bar * D CG decomposition with
NO additional unitarization correction.

Conjugate-unconjugate product:
  pi_bar^u_{n1}[a1,b1] * pi^u_{n2}[a2,b2]
    = (-1)^{b1-a1} * sum_l C_row(a1,a2,A) * pi^u_l[A,B] * C_col(b1,b2,B)

where C_row/C_col are standard SU(2) CG coefficients with the first weight negated.
"""
import numpy as np
from functools import lru_cache
from math import lgamma, log, sqrt, exp
from build.group import multiplicity, DIMS
from build.sections import total_modes


def _log_fact(n):
    """Log-factorial via lgamma for non-negative integer n."""
    return lgamma(n + 1)


@lru_cache(maxsize=None)
def _cg(j1, m1, j2, m2, J, M):
    """Cached CG coefficient <j1 m1; j2 m2 | J M> via Racah formula."""
    if abs(m1 + m2 - M) > 1e-10:
        return 0.0
    if J < abs(j1 - j2) - 1e-10 or J > j1 + j2 + 1e-10:
        return 0.0
    if abs(M) > J + 1e-10 or abs(m1) > j1 + 1e-10 or abs(m2) > j2 + 1e-10:
        return 0.0

    tj1 = int(round(2 * j1))
    tj2 = int(round(2 * j2))
    tJ = int(round(2 * J))
    tm1 = int(round(2 * m1))
    tm2 = int(round(2 * m2))
    tM = int(round(2 * M))

    if (tj1 + tj2 + tJ) % 2 != 0:
        return 0.0
    if tm1 + tm2 != tM:
        return 0.0

    # Triangle coefficient: Δ(j1,j2,J)
    a = (tj1 + tj2 - tJ) // 2
    b = (tj1 - tj2 + tJ) // 2
    c = (-tj1 + tj2 + tJ) // 2
    d = (tj1 + tj2 + tJ) // 2 + 1
    if a < 0 or b < 0 or c < 0:
        return 0.0

    log_tri = _log_fact(a) + _log_fact(b) + _log_fact(c) - _log_fact(d)

    # Prefactor
    log_pre = (log(tJ + 1)
               + _log_fact((tj1 + tm1) // 2) + _log_fact((tj1 - tm1) // 2)
               + _log_fact((tj2 + tm2) // 2) + _log_fact((tj2 - tm2) // 2)
               + _log_fact((tJ + tM) // 2) + _log_fact((tJ - tM) // 2))

    log_sqrt_arg = log_tri + log_pre
    prefactor = sqrt(exp(log_sqrt_arg))

    # Sum over k: bounds from requiring all 6 factorial args ≥ 0
    k_min = max(0, -(tJ - tj2 + tm1) // 2, -(tJ - tj1 - tm2) // 2)
    k_max = min((tj1 + tj2 - tJ) // 2, (tj1 - tm1) // 2, (tj2 + tm2) // 2)

    if k_min > k_max:
        return 0.0

    total = 0.0
    for k in range(k_min, k_max + 1):
        f1 = k
        f2 = (tj1 + tj2 - tJ) // 2 - k
        f3 = (tj1 - tm1) // 2 - k
        f4 = (tj2 + tm2) // 2 - k
        f5 = (tJ - tj2 + tm1) // 2 + k
        f6 = (tJ - tj1 - tm2) // 2 + k

        if any(x < 0 for x in (f1, f2, f3, f4, f5, f6)):
            continue

        log_denom = (_log_fact(f1) + _log_fact(f2) + _log_fact(f3)
                     + _log_fact(f4) + _log_fact(f5) + _log_fact(f6))

        sign = (-1) ** k
        total += sign * exp(-log_denom)

    return prefactor * total


_coupling_cache = {}


def _build_coupling_tensor(n1, n2, l, conj_first):
    """Build C_row[idx1, idx2, out] and C_col[idx1, idx2, out] coupling tensors.

    When conj_first=True: conjugate-unconjugate product coupling.
    When conj_first=False: unconjugate-unconjugate product coupling.
    """
    j1 = n1 / 2.0
    j2 = n2 / 2.0
    J = l / 2.0
    d1 = n1 + 1
    d2 = n2 + 1
    dl = l + 1

    C_row = np.zeros((d1, d2, dl))
    C_col = np.zeros((d1, d2, dl))

    # Row CG
    for a1 in range(d1):
        m1 = (a1 - j1) if conj_first else (j1 - a1)
        for a2 in range(d2):
            m2 = j2 - a2
            M = m1 + m2
            A_float = J - M
            A_int = int(round(A_float))
            if A_int < 0 or A_int > l:
                continue
            C_row[a1, a2, A_int] = _cg(j1, m1, j2, m2, J, M)

    # Col CG (same structure, same sign convention)
    for b1 in range(d1):
        m1p = (b1 - j1) if conj_first else (j1 - b1)
        for b2 in range(d2):
            m2p = j2 - b2
            Mp = m1p + m2p
            B_float = J - Mp
            B_int = int(round(B_float))
            if B_int < 0 or B_int > l:
                continue
            C_col[b1, b2, B_int] = _cg(j1, m1p, j2, m2p, J, Mp)

    return C_row, C_col


def _conj_product_coupling(n1, n2, l):
    """Cached conjugate-unconjugate product coupling."""
    key = ('conj', n1, n2, l)
    if key not in _coupling_cache:
        _coupling_cache[key] = _build_coupling_tensor(n1, n2, l, conj_first=True)
    return _coupling_cache[key]


def _product_coupling(n1, n2, l):
    """Cached unconjugate-unconjugate product coupling."""
    key = ('prod', n1, n2, l)
    if key not in _coupling_cache:
        _coupling_cache[key] = _build_coupling_tensor(n1, n2, l, conj_first=False)
    return _coupling_cache[key]


def _multiply_conj_fields(per_level_1, per_level_2, basis_obj_1, basis_obj_2, N_max):
    """Compute the harmonic expansion of psi_bar_1 * psi_2 (vectorized).

    Avoids materializing the O(d^4) M tensor by factoring through C_col first.
    Returns dict: l -> rho_l[A, B] (complex array of size (l+1, l+1)).
    """
    d_rho = basis_obj_1[next(iter(basis_obj_1))].shape[1]
    levels_1 = sorted(per_level_1.keys())
    levels_2 = sorted(per_level_2.keys())

    rho = {}

    for n1 in levels_1:
        A1, c1 = per_level_1[n1]
        d1 = n1 + 1
        P1 = np.conj(np.einsum('ika,ib->kab', A1, c1))

        # phase_P1[alpha, a1, b1] = (-1)^{b1-a1} * P1[alpha, a1, b1]
        a1_idx = np.arange(d1)
        b1_idx = np.arange(d1)
        phase = (-1.0) ** (b1_idx[None, :] - a1_idx[:, None])
        phase_P1 = P1 * phase[None, :, :]

        for n2 in levels_2:
            A2, c2 = per_level_2[n2]
            d2 = n2 + 1
            P2 = np.einsum('ika,ib->kab', A2, c2)

            l_min = abs(n1 - n2)
            l_max = min(n1 + n2, 2 * N_max)

            for l in range(l_min, l_max + 1, 2):
                dl = l + 1
                if l not in rho:
                    rho[l] = np.zeros((dl, dl), dtype=complex)

                C_row, C_col = _conj_product_coupling(n1, n2, l)

                # Z[alpha, a2, b1, B] = sum_{b2} P2[alpha,a2,b2] * C_col[b1,b2,B]
                Z = np.einsum('kcb,pbB->kcpB', P2, C_col)

                # W[a1, a2, B] = sum_{alpha, b1} phase_P1[alpha,a1,b1] * Z[alpha,a2,b1,B]
                W = np.einsum('kap,kcpB->acB', phase_P1, Z)

                # rho_l[A,B] += sum_{a1,a2} C_row[a1,a2,A] * W[a1,a2,B]
                rho[l] += np.einsum('abA,abB->AB', C_row, W)

    return rho


def cg_project_cubic(coeffs, basis_obj, N, rho='R0'):
    """Compute P_N[|psi|^2 psi] via CG contraction in coefficient space.

    Vectorized two-pass algorithm:
    Pass 1: density rho_l[A,B] = expansion of |psi|^2 in harmonic basis
    Pass 2: project rho * psi onto section basis
    """
    d_rho = DIMS[rho]
    levels = sorted(basis_obj.keys())

    per_level = {}
    offset = 0
    for n in levels:
        A_all = basis_obj[n]
        mult = A_all.shape[0]
        dim_n = n + 1
        block_size = mult * dim_n
        c_block = coeffs[offset:offset + block_size].reshape(mult, dim_n)
        offset += block_size
        per_level[n] = (A_all, c_block)

    # Pass 1: density
    density = _multiply_conj_fields(per_level, per_level, basis_obj, basis_obj, N)

    # Pass 2: multiply density by psi and project onto section basis
    # Reorganize: for each (n_out, n3), accumulate T over all valid l
    result_coeffs = np.zeros_like(coeffs)

    # Precompute Q3 for each level
    Q3_all = {}
    for n3 in levels:
        A3, c3 = per_level[n3]
        Q3_all[n3] = np.einsum('ika,ib->kab', A3, c3)

    # Precompute conj(A_out) for each level
    A_out_conj_all = {}
    for n_out in levels:
        A_out_conj_all[n_out] = np.conj(per_level[n_out][0])

    # Prefilter density levels with nonzero content
    active_density = {l: rho_l for l, rho_l in density.items()
                      if np.abs(rho_l).max() > 1e-30}

    offset_out = 0
    for n_out in levels:
        m_out = per_level[n_out][0].shape[0]
        d_out = n_out + 1
        c_result = np.zeros((m_out, d_out), dtype=complex)
        A_out_conj = A_out_conj_all[n_out]

        for n3 in levels:
            d3 = n3 + 1
            Q3 = Q3_all[n3]

            U_accum = np.zeros((d_rho, d_out, d_out), dtype=complex)
            any_hit = False

            for l, rho_l in active_density.items():
                if n_out < abs(l - n3) or n_out > l + n3:
                    continue
                if (l + n3 + n_out) % 2 != 0:
                    continue

                C_row, C_col = _product_coupling(l, n3, n_out)

                # V[Q, a3, A_out] = sum_P rho_l[P,Q] * C_row[P, a3, A_out]
                V = np.einsum('PQ,PaA->QaA', rho_l, C_row)

                # W[k, b3, Q, A_out] = sum_{a3} Q3[k, a3, b3] * V[Q, a3, A_out]
                W = np.einsum('kab,QaA->kbQA', Q3, V)

                # U_l[k, A_out, B_out] = sum_{b3, Q} W[k, b3, Q, A_out] * C_col[Q, b3, B_out]
                U_accum += np.einsum('kbQA,QbB->kAB', W, C_col)
                any_hit = True

            if not any_hit:
                continue

            c_result += np.einsum('ikA,kAB->iB', A_out_conj, U_accum)

        c_result /= (n_out + 1)
        block_size = m_out * d_out
        result_coeffs[offset_out:offset_out + block_size] = c_result.ravel()
        offset_out += block_size

    return result_coeffs
