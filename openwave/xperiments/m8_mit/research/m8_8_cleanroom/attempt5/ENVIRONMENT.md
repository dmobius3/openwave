# M8.8 Environment Declaration

## Runtime

- **Language:** Python 3 (CPython)
- **Arithmetic:** `fractions.Fraction` for exact rational arithmetic throughout
- **Field extensions:** Q(φ) via custom `QG` class; Q(φ,𝕚) via custom `GC` class
- **No floating-point arithmetic** in any intermediate computation; floats appear only in final display/verification
- **No external libraries:** only Python standard library (`fractions`, `json`, `hashlib`, `math.comb`, `math.gcd`)

## Determinism

- All computations are exact over Q(φ) or Q(φ,𝕚)
- No randomness, no floating-point rounding, no numerical tolerance
- Reproducible on any platform with Python 3.8+

## Isolation

- All code runs within `/Users/blake/Desktop/OpenWave/M8.8_CLEANROOM/`
- No network access
- No imports from outside the cleanroom directory
- No git operations
- No access to M8.3 implementation, spectral data, or any external torsion tables

## Input files

- `m8_5a_packet.json` — group packet defining 2I
- `m8_8_construction_packet.json` — based chain complex over Z[2I]

## Output files

- `RAW_OUTPUT.json` — T²(ρ) for all 8 nontrivial irreps in protocol § 5.5 schema
- `METHOD_AND_GATE_MANIFEST.md` — method, conventions, gates, coverage
- `ENVIRONMENT.md` — this file
- `CONSULTED_FILES.md` — list of consulted files

## Validation artifacts

- `validate_enumeration.py` — group enumeration and SHA-256
- `validate_complex.py` — chain complex properties over Z[2I]
- `validate_saturation.py` — integral homology via saturation certificates
- `validate_representations.py` — all 9 irreps with character orthogonality
- `validate_torsion_dry.py` — twisted complex ∂∂=0 and acyclicity
- `validate_fixture.py` — convention fixture with mutation tests
- `validate_manifest.py` — registry-coverage checker

## Production code

- `compute_torsion.py` — computes T²(ρ) for all nontrivial irreps
