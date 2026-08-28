# § 14 execution package: disposition table and verdict

One pass, per the pre-recorded provenance in `00_PROVENANCE.md`. Three layers: (I) inputs
and provenance = `00_PROVENANCE.md` plus `raw/package_hashes.txt`; (II) raw execution
records = everything under `raw/`; (III) this table.

| # | obligation | class | disposition | record |
| --- | --- | --- | --- | --- |
| 1 | script 5 green at pinned hash `752baa2d…` | FALLBACK-COVERED (TRIGGER S) | GREEN: nine armed report lines, exit 0, hash matches the pin | `raw/rerun_right_translation_check.py.out` |
| 2 | the other four design-input scripts re-run | BLOCKING | GREEN: all exit 0, outputs match the pinned claims | `raw/rerun_*.out`, `raw/rerun_script_hashes.txt` |
| 3 | wall-clocks at the top rungs, scratch, feasibility-only | BLOCKING | GREEN: composed 4.35 h against the 48 h ceiling, margin 11x; licenses FEASIBILITY only, and crossing the ceiling is STOP-QUAL regardless | `raw/item3_wallclocks.json` |
| 4 | isotropy lattice tables, `l = 0..7, 12`, armed (script 6, author-side) | FALLBACK-COVERED (TRIGGER L) | GREEN: dims by two routes agree at every class; representatives fixed-by-class and moved-by-non-member; arms A (2T-for-2I), B (twist discrimination at `l = 4`, non-vacuous), C (perturbed icosahedral rep at `l = 12`, green parent first) all fire | `raw/lattice_tables.json` (sha `4e8b7f72…`) |
| 5 | Control A references at 50 dps under the frozen definitions | BLOCKING | GREEN: structured 50-dps assembly; arm P (s = 0 limit equals `diag(9 − λ)` exactly, err 0.0); arm X (generic float64 node assembly agrees, `1.2e-14`); arms M, N, S pass; separations ≈ 38 vs threshold 10; zero-count reference 0 with the smallest cluster eigenvalue O(1), fold quiescent | `raw/controlA_references.json` (sha `567fe152…`) |
| 6 | control (i) momentum-liveness on the pinned `φ` | BLOCKING | **RED. The pinned equal-weight `φ` fails its own assertion: `M₁ = 0` EXACTLY (any real `φ` gives `M₁ ≡ 0`, the generator is real antisymmetric) and `M₃ ≈ 0` by weight symmetry; only `M₂` is live.** Diagnosis and a VERIFIED amendment candidate recorded: `φ' ∝ Σ (1 + j/13) e^{ij/3} v_j` gives all three live at [+2.09e+01, +6.05e+01, −1.75e+01] | `raw/item6_momenta.json`, `raw/item6_run.out` |
| 7 | denominator preflight, three pairs, EVERY rung of both ladders | BLOCKING | GREEN: global min `D_I = 6.96e+02`, fifteen orders above machine floor; the stream machinery independently reproduced the 677 modes at `N = 60` and the multiplicity-2 level | `raw/item7_denominators.json` |
| 8 | four-site fallback cross-check | BLOCKING | GREEN: §§ 2, 5, gate 5's row, and § 14's table carry the same two degraded contracts, seven claims checked | `raw/item8_crosscheck.json` |

## Package verdict

**RED, on item 6 alone.** Seven of eight obligations are green; item 6 is BLOCKING and red,
so the package reports red per the one-pass discipline: no repair inside the package. The
red is the obligation WORKING: it exists to discover a bad control pin before the freeze
rather than at execution, and it did, with the failure mode sharpened (any real `φ` kills
`M₁` identically, so no real-weighted amendment can work) and a deterministic complex
candidate verified all-live. The required action is a PROTOCOL AMENDMENT pre-freeze: replace
§ 7's pinned `φ` with the recorded candidate (or another complex-weighted pin), then re-run
item 6 alone as the amendment's verification.

## Development record, disclosed

Item 4's first run had two defective ARMS (a vacuous twist predicate probed at `l = 2` where
all dims vanish; a perturbed-representative arm built at `l = 6` where the icosahedral
invariant does not exist); both were arm-construction bugs in the scratch tooling, fixed
before the item's final armed run, first output preserved as `raw/lattice_tables_dev1.json`.
Item 5's first implementation (generic 50-dps node quadrature, `raw/s14_item5_dev1.py`) was
killed on a multi-hour trajectory and carried a tautologous arm P; the structured
reimplementation runs in 2 s, and its arm P is the exact known limit. Item 6 had a JSON
serialization crash after the (already red) result printed; fixed, same result. One
invocation error (items 1 and 2 first run with the wrong PYTHONPATH) is recorded in
`raw/invocation_note.txt`.
