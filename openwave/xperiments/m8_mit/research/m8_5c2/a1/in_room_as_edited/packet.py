"""§ 4.3 manufactured field packet: 40 pseudorandom band-limited fields per rung,
20 scalar E_R0 + 20 manufactured E_R0 ⊗ C² extension.

Generation from a SINGULAR stream: PCG64(20260901).
Traversal order per § 4.3:
  Agreement ladder ascending {24, 32, 40, 48}, then
  Control-B ladder ascending {36, 44, 52, 60}.
  Within each rung: 20 scalar fields (indices 0-19), then 20 ⊗C² fields (0-19).
  Within a field: coefficients in level-major, intertwiner-index, multiplet-index,
    then component-index (C² component innermost, varying fastest;
    scalar fields have no component index).
  Per coefficient: real part then imaginary part, each a standard normal draw.
  Each field normalized to unit L2 norm.

Canonical serialization: little-endian IEEE-754 float64 (real, imag) pairs in C-order,
concatenated in draw order.
"""
import numpy as np
import hashlib
from build.group import multiplicity


AGREEMENT_RUNGS = [24, 32, 40, 48]
CONTROL_B_RUNGS = [36, 44, 52, 60]
RUNG_ORDER = AGREEMENT_RUNGS + CONTROL_B_RUNGS

N_SCALAR = 20
N_C2 = 20
N_FIELDS_PER_RUNG = N_SCALAR + N_C2

SEED = 20260901


def r0_modes_at_rung(N):
    """Total complex Galerkin modes for R0 at cutoff N."""
    return sum(multiplicity('R0', n) * (n + 1) for n in range(N + 1))


def generate_packet():
    """Generate the full § 4.3 field packet.
    Returns dict: rung -> {'scalar': array(20, modes), 'c2': array(20, modes, 2)}
    and the canonical serialization bytes.
    """
    rng = np.random.Generator(np.random.PCG64(SEED))
    packet = {}
    serial_parts = []

    for N in RUNG_ORDER:
        modes = r0_modes_at_rung(N)
        scalars = np.zeros((N_SCALAR, modes), dtype=complex)
        c2s = np.zeros((N_C2, modes, 2), dtype=complex)

        # 20 scalar fields
        for f in range(N_SCALAR):
            coef = np.zeros(modes, dtype=complex)
            idx = 0
            for n in range(N + 1):
                m = multiplicity('R0', n)
                if m == 0:
                    continue
                dim_n = n + 1
                for i_interp in range(m):
                    for j_mult in range(dim_n):
                        re = rng.standard_normal()
                        im = rng.standard_normal()
                        coef[idx] = re + 1j * im
                        idx += 1
            assert idx == modes
            coef /= np.linalg.norm(coef)
            scalars[f] = coef
            for z in coef:
                serial_parts.append(np.float64(z.real).tobytes())
                serial_parts.append(np.float64(z.imag).tobytes())

        # 20 C² fields
        for f in range(N_C2):
            coef = np.zeros((modes, 2), dtype=complex)
            idx = 0
            for n in range(N + 1):
                m = multiplicity('R0', n)
                if m == 0:
                    continue
                dim_n = n + 1
                for i_interp in range(m):
                    for j_mult in range(dim_n):
                        for comp in range(2):
                            re = rng.standard_normal()
                            im = rng.standard_normal()
                            coef[idx, comp] = re + 1j * im
                        idx += 1
            assert idx == modes
            coef /= np.linalg.norm(coef)
            c2s[f] = coef
            for idx_flat in range(modes):
                for comp in range(2):
                    z = coef[idx_flat, comp]
                    serial_parts.append(np.float64(z.real).tobytes())
                    serial_parts.append(np.float64(z.imag).tobytes())

        packet[N] = {'scalar': scalars, 'c2': c2s}

    serial_bytes = b''.join(serial_parts)
    return packet, serial_bytes


def packet_digest():
    """SHA-256 of the canonical serialization."""
    _, serial_bytes = generate_packet()
    return hashlib.sha256(serial_bytes).hexdigest()
