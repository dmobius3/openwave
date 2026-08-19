# Environment

## Interpreter

- Python 3.13.13

## Platform

- macOS (Darwin 24.6.0), arm64

## Libraries

- Standard library only: `json`, `hashlib`, `sys`, `re`, `fractions.Fraction`
- No third-party packages

## Files

| File | Purpose |
| --- | --- |
| `m88_torsion.py` | Main computation: group construction, representations, torsion, output |
| `m88_gates.py` | All 19 gates with mutation tests |

## Reproducing

```
python3 m88_gates.py    # runs all gates and mutations, writes gate_results.json
python3 m88_torsion.py  # runs torsion computation, writes RAW_OUTPUT.json
```

Both scripts must be run from this directory. No network access, no external data.
