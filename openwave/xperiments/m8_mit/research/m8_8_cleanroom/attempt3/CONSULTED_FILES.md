# Consulted Files

## Files read during this computation

| File | Purpose |
| --- | --- |
| `TASK.md` | Task specification and deliverable order |
| `PROTOCOL.md` | Governing protocol (hashes, conventions, gate requirements) |
| `m8_5a_packet.json` | Group packet (2I generators as quaternions in Q(phi)) |
| `m8_8_construction_packet.json` | Based chain complex for S^3/2I |
| `METHOD_AND_GATE_MANIFEST.md` | This reproduction's method manifest (read for SHA-256) |
| `validate_manifest.py` | Pre-implementation validation script (read for SHA-256) |
| `reproduce.py` | This production implementation (self-hash) |

## Files written during this computation

| File | Purpose |
| --- | --- |
| `reproduce.py` | Production implementation |
| `DERIVATION_ARTIFACTS.json` | Route-native intermediates (protocol section 7) |
| `RAW_OUTPUT.json` | Torsion values and gate results |
| `ENVIRONMENT.md` | Interpreter and platform information |
| `CONSULTED_FILES.md` | This file |

## External references

**None.** No network access, no external literature, no files outside this directory
were consulted. All mathematical content derives from:
1. The supplied packets (group and construction)
2. The protocol document
3. Standard algebraic facts (quaternion multiplication, symmetric power construction,
   determinant formulas, Reidemeister torsion definition) applied from first principles

This is an affirmative statement: the context firewall was maintained throughout.
