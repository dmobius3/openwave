# M8.8 Consulted Files

## Files read during this implementation

| File | Purpose |
|------|---------|
| `TASK.md` | Operational instructions for the clean room |
| `PROTOCOL.md` | Governing document: reproduction target, method class, firewall, gates, adjudication |
| `m8_5a_packet.json` | Group packet: generators and coefficient field for 2I |
| `m8_8_construction_packet.json` | Construction packet: boundary maps of the based chain complex |
| `METHOD_AND_GATE_MANIFEST.md` | The manifest produced in this session (read back for SHA-256 hashing and by validate_manifest.py) |
| `validate_pre_impl.py` | Pre-implementation validation script produced in this session (read back for reference during production implementation) |
| `validate_manifest.py` | Manifest validator produced in this session |

## Files produced during this implementation

| File | Purpose |
|------|---------|
| `validate_pre_impl.py` | Pre-implementation validation of group enumeration, boundary maps, and model gates |
| `validate_manifest.py` | Automated manifest validator |
| `METHOD_AND_GATE_MANIFEST.md` | Method-and-gate manifest (finalized before production implementation) |
| `torsion.py` | Production implementation: representation construction, torsion computation, all gates |
| `RAW_OUTPUT.json` | Raw output in the protocol's frozen schema |
| `ENVIRONMENT.md` | Environment record |
| `CONSULTED_FILES.md` | This file |

## External references

No external references were consulted. No web access, no literature lookup, no published
source of any kind was used. All computation was performed from the permitted inputs
(the two packets in this directory) and generic algebraic knowledge of group rings,
representation theory, and chain complexes. This is stated affirmatively: the implementation
was derived entirely from the supplied packets and standard mathematical technique, with
no external consultation.
