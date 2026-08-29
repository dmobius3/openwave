# M8.5-C Build Unit COMMISSIONING record

Attempt: M85C-A1
Commissioned: 2026-08-29T03:57:32Z

## 1. Provenance checks

### 1.1 Room manifest

ROOM_MANIFEST.json SHA-256: `02429a9516cb009e9d7baae9c30e1800ba91995858f4724cc51fc644c6dac815`
Handoff value: `02429a9516cb009e9d7baae9c30e1800ba91995858f4724cc51fc644c6dac815`
Result: MATCH

### 1.2 Protocol frozen region

Computed frozen region SHA-256:
`e253558b5a767084d4d7777550ac72de5b8a0591ec3d2b847108f04e17c0cc6b`

Protocol's own freeze record:
`e253558b5a767084d4d7777550ac72de5b8a0591ec3d2b847108f04e17c0cc6b`

ROOM_MANIFEST protocol_frozen_region_sha256:
`e253558b5a767084d4d7777550ac72de5b8a0591ec3d2b847108f04e17c0cc6b`

Freeze boundary marker count: 1 (exactly once)
Result: MATCH (all three agree)

### 1.3 Whitelist file hashes

All 13 files verified against ROOM_MANIFEST.json recorded hashes:

| file | manifest hash | computed hash | result |
| --- | --- | --- | --- |
| 00_COMMISSION.md | 83fb78f7…9a44c | 83fb78f7…9a44c | MATCH |
| m8_5c_protocol.md | c8d9941e…9928 | c8d9941e…9928 | MATCH |
| m8_5c_symmetry_derivation.md | 5b231af4…820d1a | 5b231af4…820d1a | MATCH |
| route_a_nonabelian.py | e5572853…931b | e5572853…931b | MATCH |
| design_inputs/README.md | 86d8b604…4654 | 86d8b604…4654 | MATCH |
| cascade_quad_check.py | 7f31bc30…6bfa | 7f31bc30…6bfa | MATCH |
| exact_quad_check.py | 51fa9b52…55ed7 | 51fa9b52…55ed7 | MATCH |
| jacobian_check.py | dce687ea…9d8 | dce687ea…9d8 | MATCH |
| mode_count.py | 794d7063…59a6 | 794d7063…59a6 | MATCH |
| right_translation_check.py | 752baa2d…8f547 | 752baa2d…8f547 | MATCH |
| controlA_references.json | dde909f6…e4e44 | dde909f6…e4e44 | MATCH |
| lattice_tables.json | efe6c2d6…b504 | efe6c2d6…b504 | MATCH |
| room_import_gate.py | 097c5366…6b9d | 097c5366…6b9d | MATCH |

### 1.4 § 15 freeze-time pins vs manifest

All files appearing in both § 15 freeze-time pins and the manifest have matching hashes:
- route_a_nonabelian.py: AGREE
- mode_count.py: AGREE
- exact_quad_check.py: AGREE
- cascade_quad_check.py: AGREE
- jacobian_check.py: AGREE
- right_translation_check.py: AGREE
- m8_5c_symmetry_derivation.md: AGREE
- lattice_tables.json: AGREE
- controlA_references.json: AGREE

### 1.5 Collision guard

ledger/: empty (no prior attempt records)
build/: empty (no prior attempt records)
Result: CLEAN

## 2. Launch gate

### 2.1 Selftest

```
RED as required: library with unresolvable import
  RED as required: script exiting 0 early, no verdict line
SELFTEST: GREEN, both arms fire
```

### 2.2 Gate

```
room import gate (§ 12 semantics: per-module subprocess + sentinel; two probe classes)
interpreter: 3.13.13 at /opt/homebrew/Caskroom/miniforge/base/bin/python3
  PASS  LIB    route_a_nonabelian           rc=0
  PASS  SCRIPT mode_count.py                
  PASS  SCRIPT jacobian_check.py            
  PASS  SCRIPT exact_quad_check.py          
  PASS  SCRIPT cascade_quad_check.py        
  PASS  SCRIPT right_translation_check.py   
  PASS  3RDPTY numpy                        version 2.5.0
  PASS  3RDPTY scipy                        version 1.18.0
GATE: GREEN, room launches
```

Exit code: 0

## 3. Environment

- Interpreter: /opt/homebrew/Caskroom/miniforge/base/bin/python3
- Python: 3.13.13
- numpy: 2.5.0
- scipy: 1.18.0
- Platform: darwin (Darwin 24.6.0, arm64)

Design-inputs README records evidence re-run under: python 3.13.13, numpy 2.5.0, scipy 1.18.0
Version difference: NONE

## 4. INPUT_MANIFEST closing hash

ledger/INPUT_MANIFEST.json SHA-256: `f4a3ddf8415c828a582fb5fcf47c3ca424fa53f29c4a635d7789bb0630f57d9a`

Contents: whitelist (13 files), arena registry (82 entries), law exception (A-R0-N36 at c1=0),
§4.3 packet digest (`2055fa05ab020deda51db11e1a34c2a8f208b8a2848f720986674c29314e9055`),
§3 basis-object hashes (72 entries: 9 sectors × 8 rungs),
controlA_references hash (`dde909f6746994468805c55c1eae192b449b2eb4862e2432f61d55a9f3be4e44`),
lattice_tables hash (`efe6c2d64ba9664edaf8b49b9026a683bf93aa78c0c994f09f41bd5a63a7b504`),
commissioning hash (`b5a0d622e00a258432cc3199f47430fb98109aff3570f74c01659fad5b6db7cd`)

Manifest closed before first GATE record.

## 5. Gate execution record

### Gates 1–3: GREEN/FIRES

All recorded in OUTPUT_LEDGER.jsonl. See ledger for full measured values.

| Gate | Parent | Mutation | Wall-clock |
| --- | --- | --- | --- |
| G1-GLIN | GREEN | FIRES | 0.1s |
| G2-GLABEL | GREEN | FIRES | 3.2s |
| G3-BASES | GREEN | FIRES | 2.0s |

### Gate 4: Projector exactness — PARENT GREEN, MUTATION DEAD

Recorded in OUTPUT_LEDGER.jsonl at 2026-08-29T06:05:17Z.

Parent check (N=24, 4N=96 rule, 461,041 nodes, 60 R0 modes):
- Gram error: 7.66e-15 (threshold: 1e-10) — PASS
- Max roundtrip error: 6.36e-15 (threshold: 1e-10) — PASS
- 5 packet fields tested, all exact to rounding

Node-drop mutation (2N=48 rule on projected cubic):
- 2N cubic disagreement: 7.44e-15 — **machine zero, does NOT fire**
- N=24 rule cubic disagreement: 1.06e-02 (supplementary)

## 6. STOP — Gate 4 node-drop mutation cannot fire as written

**Invoked per commission §6 line 104:** "If a check cannot be implemented as written,
STOP and report. Never substitute a weaker check silently."

**Also per commission §6 line 102:** "Every PASS needs an arm that can fail, run in the
same session."

### The specified check

Protocol §8 gate table (line 622): gate 4 mutation arm states
"node-drop to 2N must err O(1)."

Protocol §4.1 (line 295): "The design-input record, 1.5e-01 against 2.2e-15, was taken
at N = 3, which sits on neither production ladder."

### What was measured

| Rule | Nodes | Cubic disagreement vs 4N |
| --- | --- | --- |
| 4N (D=96) | 461,041 | — (reference) |
| 2N (D=48) | 60,025 | 7.44e-15 (machine zero) |
| N (D=24) | 8,125 | 1.06e-02 |

The 2N rule gives **identically** the same projected cubic as 4N for R0 fields.

### Root cause (mathematical, not implementation)

R0 invariant vectors under 2I have nonzero entries **only at even k indices**
(k = 0, 2, 4, …, n) in the unitarized sym_power basis. This is verified:
`build_r0_basis` returns intertwiner matrices whose columns are zero at all odd rows.

This even-k constraint propagates through the cubic |ψ|²ψ: when ψ has even-k-only
support, the products ψ·ψ* and ψ·ψ*·ψ also have even-parity Fourier content in the
Hopf angular coordinates (ξ₁, ξ₂).

The critical aliasing frequency for a degree-2N quadrature rule applied to degree-4N
polynomial content is K = 2N + 1, which is **odd**. Since R0 cubic products have
zero energy at odd frequencies, the aliasing frequency is unreachable. The 2N rule
integrates R0 cubics **exactly** — the node-drop mutation cannot fire.

Verification: at N=3 with **general** (non-R0) fields, the 2N node-drop produces
22% error, confirming the mechanism is correct and the code is sound. The R0 parity
structure was not anticipated in the design-input, which used general fields.

### What is NOT affected

1. **Projector parent check**: GREEN. The 4N production rule is exact to rounding.
   The projector works correctly.

2. **Gate 4's second arm** (protocol line 622): "route disagreement is STOP-QUAL."
   The CG contraction dual-route provides an independent falsifiability mechanism.
   It is not yet implemented (requires SU(2) Clebsch-Gordan coefficients) but the
   gate CAN go red via this route, arguably satisfying "must be able to go red."

3. **All other gates**: gates 1–3 are GREEN/FIRES. Gates 5–10 are unaffected.

### The decision I am not making

Two readings are defensible:

**(A) The dead arm kills gate 4.** The protocol says "must err O(1)"; it doesn't;
the mutation is dead; commission line 102 says every PASS needs a failing arm; gate 4
fails; STOP-QUAL per §11.

**(B) The dead arm is a finding; gate 4 survives.** The gate has a second arm
(route disagreement) that CAN fail. Commission line 102 says "an arm" (singular).
The protocol's "must be able to go red" (line 632) is satisfied by route disagreement.
The projector parent check is GREEN, and the dead mutation is an unanticipated
mathematical property of R0, not a code deficiency.

Per commission §1 (line 12): "If this commission and the protocol ever disagree, the
protocol wins; if the disagreement is material, STOP and report it instead of choosing."

The protocol's stated expectation (O(1) error at 2N) cannot materialize as written
for R0 fields. **The disposition of this finding belongs to the author, not the
Build Unit.**

### Ledger state at STOP

- Gates 1–3: GREEN/FIRES (6 records)
- Gate 4: parent GREEN, mutation DEAD (2 records)
- Gates 5–10: not reached
- Gate 11 record: not reached
- Cumulative wall-clock: ~333s of 172,800s ceiling
- CG dual-route: not implemented
- Fast Hopf transform: not implemented (needed for rungs beyond N=24)
