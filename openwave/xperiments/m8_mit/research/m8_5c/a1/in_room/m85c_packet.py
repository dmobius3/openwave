"""Section 4.3 preregistered manufactured field set and canonical serialization.

Generation: SINGULAR stream from numpy.random.Generator(numpy.random.PCG64(20260901)).
Traverse agreement ladder ascending {24,32,40,48}, then control-B ladder ascending
{36,44,52,60}. Within each rung: 20 scalar fields (indices 0-19), then 20 C^2
extension fields (indices 0-19). Within a field: one coefficient per basis element
in level-major, intertwiner-index, multiplet-index, component-index order (C^2
innermost). Per coefficient: real then imaginary, standard normal draws. Field
normalized to unit L2.

Serialization: little-endian IEEE-754 float64 (real, imag) pairs in C-order.
"""
import numpy as np
import hashlib
from m85c_group import multiplicity


AGREEMENT_RUNGS = [24, 32, 40, 48]
CONTROL_B_RUNGS = [36, 44, 52, 60]
ALL_RUNGS = AGREEMENT_RUNGS + CONTROL_B_RUNGS
N_SCALAR = 20
N_C2 = 20
N_PER_RUNG = N_SCALAR + N_C2
PACKET_SEED = 20260901


def r0_mode_count(N):
    """Total complex R0-invariant modes at cutoff N."""
    total = 0
    for n in range(N + 1):
        m = multiplicity("R0", n)
        if m > 0:
            total += m * (n + 1)
    return total


def r0_invariant_levels(N):
    """List of R0-invariant levels up to N."""
    return [n for n in range(N + 1) if multiplicity("R0", n) > 0]


def generate_packet():
    """Generate the full preregistered field set.
    Returns: dict {rung: {"scalar": list of arrays, "c2": list of arrays}},
    and the SHA-256 digest of the canonical serialization."""
    rng = np.random.Generator(np.random.PCG64(PACKET_SEED))
    h = hashlib.sha256()

    packet = {}
    all_rungs = AGREEMENT_RUNGS + CONTROL_B_RUNGS

    for N in all_rungs:
        M = r0_mode_count(N)
        scalar_fields = []
        c2_fields = []

        for idx in range(N_SCALAR):
            coeffs = np.empty(M, dtype=complex)
            for k in range(M):
                coeffs[k] = complex(rng.standard_normal(), rng.standard_normal())
            norm = np.sqrt(np.sum(np.abs(coeffs)**2))
            coeffs /= norm
            scalar_fields.append(coeffs)
            for c in coeffs:
                h.update(np.float64(c.real).tobytes())
                h.update(np.float64(c.imag).tobytes())

        for idx in range(N_C2):
            coeffs = np.empty(2 * M, dtype=complex)
            for k in range(2 * M):
                coeffs[k] = complex(rng.standard_normal(), rng.standard_normal())
            norm = np.sqrt(np.sum(np.abs(coeffs)**2))
            coeffs /= norm
            c2_fields.append(coeffs)
            for c in coeffs:
                h.update(np.float64(c.real).tobytes())
                h.update(np.float64(c.imag).tobytes())

        packet[N] = {"scalar": scalar_fields, "c2": c2_fields}

    return packet, h.hexdigest()


def generate_multiseed_set(N_proj=36):
    """Generate the deterministic multi-seed set for branch enumeration (§ 5).
    The N=36 rung's twenty SCALAR fields from the § 4.3 stream are PROJECTED
    onto H_{R0,12}, normalized survivors become the seed set."""
    rng = np.random.Generator(np.random.PCG64(PACKET_SEED))

    for N in AGREEMENT_RUNGS + CONTROL_B_RUNGS:
        M = r0_mode_count(N)
        for _ in range(N_SCALAR):
            for _ in range(M):
                rng.standard_normal()
                rng.standard_normal()
        for _ in range(N_C2):
            for _ in range(2 * M):
                rng.standard_normal()
                rng.standard_normal()
        if N == N_proj:
            break

    return None
