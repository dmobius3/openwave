# M8.8 Environment Record

## Interpreter

- Python 3.13.13

## Platform

- macOS (Darwin 24.6.0, arm64)
- Machine: Apple Silicon (M-series)

## Libraries

No external libraries used. The implementation uses only Python standard library modules:

- `json` (packet parsing, output serialization)
- `hashlib` (SHA-256 checksums)
- `sys` (exit codes)
- `fractions.Fraction` (exact rational arithmetic)
- `math.comb` (binomial coefficients for symmetric powers)
- `re` (parsing generator strings from the group packet)

No `numpy`, `scipy`, `sympy`, or any other third-party package is imported.

## Arithmetic

All arithmetic is exact. Elements of Q(phi) are represented as pairs of `Fraction` objects
(a, b) meaning a + b*phi. Elements of Q(phi, i) are represented as pairs of Q(phi)
elements (re, im). No floating-point approximation is used at any stage of the computation.

## Reproduction

To rerun: `python3 torsion.py` from the clean room directory. The script reads
`m8_5a_packet.json` and `m8_8_construction_packet.json`, and writes `RAW_OUTPUT.json`.
Exit code 0 indicates all gates passed with all mutations red.
