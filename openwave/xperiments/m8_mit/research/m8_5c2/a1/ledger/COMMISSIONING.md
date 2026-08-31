# M8.5-C2 COMMISSIONING RECORD

Attempt ID: C2-A1

## Step 1 — ROOM_MANIFEST.json integrity

```
$ shasum -a 256 ROOM_MANIFEST.json
54ffc7ea770e9d005119d6dc113549227b66a4ff9ae2efb2bce051ce3cf79720  ROOM_MANIFEST.json
```

Expected (handoff value): `54ffc7ea770e9d005119d6dc113549227b66a4ff9ae2efb2bce051ce3cf79720`

**RESULT: PASS** — hashes match; room has not changed since commissioning.

## Step 2 — Collision guard

`ledger/` — empty (no entries besides `.` and `..`)
`build/` — empty (no entries besides `.` and `..`)

No prior records found.

**RESULT: PASS** — no prior attempt records present.

## Step 4 — Protocol frozen-region integrity

Command (run from `openwave/xperiments/m8_mit/research/findings/`):
```
$ sed '/^<!-- M85C2-FREEZE-BOUNDARY -->$/,$d' m8_5c2_protocol.md | shasum -a 256
aadf4d9218bc36cf4763b42ce587725c16f673b3a9e4272f19a33b6ec0fbd5ef  -
```

Expected (freeze record below boundary): `aadf4d9218bc36cf4763b42ce587725c16f673b3a9e4272f19a33b6ec0fbd5ef`

FREEZE-BOUNDARY marker (`<!-- M85C2-FREEZE-BOUNDARY -->`) full-line occurrences: 1

**RESULT: PASS** — frozen region hash matches freeze record; COMPLETE marker occurs exactly once.

## Step 5 — File integrity: ROOM_MANIFEST hashes and § 15 freeze-time cross-check

All 13 whitelist files verified against ROOM_MANIFEST.json recorded hashes:

| file | status |
| --- | --- |
| `00_COMMISSION.md` | PASS — `4338252b…` matches manifest |
| `openwave/xperiments/m8_mit/research/findings/m8_5c2_protocol.md` | PASS — `6a58fe74…` matches manifest |
| `openwave/xperiments/m8_mit/research/findings/m8_5c_symmetry_derivation.md` | PASS — `5b231af4…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5b/pilot/route_a_nonabelian.py` | PASS — `e5572853…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/README.md` | PASS — `86d8b604…` matches manifest |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/cascade_quad_check.py` | PASS — `7f31bc30…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/exact_quad_check.py` | PASS — `51fa9b52…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/jacobian_check.py` | PASS — `dce687ea…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/mode_count.py` | PASS — `794d7063…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/design_inputs/right_translation_check.py` | PASS — `752baa2d…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/s14/raw/controlA_references.json` | PASS — `dde909f6…` matches manifest AND § 15 |
| `openwave/xperiments/m8_mit/research/m8_5c/s14/raw/lattice_tables.json` | PASS — `efe6c2d6…` matches manifest AND § 15 |
| `room_import_gate.py` | PASS — `e0a1eef8…` matches manifest |

9 of 13 files have § 15 freeze-time pins; all 9 three-way matches confirmed (file == manifest == § 15).
4 files without § 15 pins (`00_COMMISSION.md`, protocol itself, `README.md`, `room_import_gate.py`) verified against manifest only.

**RESULT: PASS** — all files verified; no mismatches.

## Step 6 — Import gate selftest

```
$ PYTHONDONTWRITEBYTECODE=1 python3 room_import_gate.py --selftest
RED as required: library with unresolvable import
  RED as required: script exiting 0 early, no verdict line
  RED as required: token printed EARLY then exit 0 above the true terminal line (the presence-test defect)
SELFTEST: GREEN, all three arms fire
```

**RESULT: PASS** — SELFTEST GREEN, all three arms fire.

## Step 7 — Import gate (full)

```
$ PYTHONDONTWRITEBYTECODE=1 python3 room_import_gate.py
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

Interpreter and library versions:
- Python 3.13.13 (conda-forge, Clang 19.1.7)
- numpy 2.5.0
- scipy 1.18.0

ROOM_MANIFEST recorded versions: python 3.13.13, numpy 2.5.0, scipy 1.18.0
Difference from design-inputs README: the README does not record interpreter/library versions; no version drift to report.

**RESULT: PASS** — GATE GREEN, room launches.

## § 3 Build — basis objects, arena registry, § 4.3 packet

All 9 sector bases built via the deterministic construction: stacked-constraint SVD
(gesvd driver) in the frozen 120-icosian element order, CG-recurrence unitarized
representation matrices for numerical stability at high n, Lowdin symmetric
orthonormalization against the per-level analytic Gram, sign fixed by the first nonzero
component. All basis objects at all levels n ≤ 180 (3 × 60) verified:
- Equivariance: `max |ρ(g) A - A π_n(g)| < 3e-15` at all levels
- Lowdin orthonormality: `max |H - I| < 1e-10` at all levels
- Multiplicity matches character route at all levels

§ 10 arena-constructor registry built per § 10's enumerated table.

§ 4.3 field packet generated from `PCG64(20260901)`: 40 fields per rung (20 scalar
E_R0 + 20 E_R0 ⊗ C²) at each of 8 rungs, draw order per § 4.3's specification.

## INPUT_MANIFEST close

`ledger/INPUT_MANIFEST.json` closed. SHA-256 of the closed manifest:

```
bfa210b5c3c78790c80e5e9249619cb22e8732d71e7a9899cad13dd2fe74175d
```

Coverage: SHA-256 of `ledger/INPUT_MANIFEST.json` bytes after writing, before any gate executes.

Contents: whitelist files (verbatim from ROOM_MANIFEST.json), arena registry (§ 10),
§ 4.3 packet digest (`1210596955a2680f62e7f2b90ce19b95d0d8f6c168cf363f8322de6f52259aaf`),
§ 3 basis-object hash per sector and rung, Control A reference and lattice-table hashes,
commission record hash. Every entry carries its § 13 coverage statement.

The INPUT_MANIFEST is now frozen and will not change.
