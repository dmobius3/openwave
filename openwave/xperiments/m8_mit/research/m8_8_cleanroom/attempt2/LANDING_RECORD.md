# M8.8 clean-room attempt 2: the output, landed unread against the answer packet

> **Claims no result.** The canonical answer packet has not been opened, its hash has not been
> verified, no comparison has run, and no verdict exists. This commit discharges § 3's
> requirement that the implementation, environment record, derivation artifacts and raw output
> be committed BEFORE any quarantined object is unsealed, and the ordering is carried by
> commit ancestry rather than by attestation.

## What landed

| File | SHA-256 | Role |
| --- | --- | --- |
| `METHOD_AND_GATE_MANIFEST.md` | `b623535beb2d59f153725e22d05d97b306ea1c3709cd73ee0fb1d115d8fe840c` | § 8 step 4, written 09:07 and never modified thereafter |
| `MANIFEST_MISMATCH.md` | `2ff3742607ed50ffef3ef57d5d6edff8f02ffa3bdccbcd6ed88ae75f09fc131b` | the implementer's disclosure of a defect in its own manifest |
| `compute_torsion.py` | `47b21134c72ed8ebb2366f0cdcbb47e31635721a34b67a12fb911a073e62c61c` | the implementation |
| `ENVIRONMENT.md` | `e809f409c52c5a18a311c1f49dfe3892acda266a5a3dc80931732b3b5935fccd` | environment record |
| `RAW_OUTPUT.json` | `6a60b7af43f06baf6d516b0f3ff3833593bc50a4b106a69b9586e80442167e8d` | the § 5.5 raw output |
| `CONSULTED_FILES.md` | `363bbdb283f4dd63075ead5b3fdba8a828ebb8204cffe8dae77789129f6c07d8` | consulted-files manifest, generated in the room |
| `TASK.md` | `088200506c0a8eb3b81b4de501e4147cee23853c77eb83c168f5da93a5bedc43` | the operational instruction, committed BEFORE this run started |

## The manifest-order defect that retired attempt 1 did NOT recur

The manifest was written at 09:07 and its bytes were never modified again; every other
artifact in the room postdates it by more than an hour. When the implementation exposed a
mismatch with the manifest, the implementer did not amend the manifest. It wrote
`MANIFEST_MISMATCH.md` naming the defect, quoting the manifest's own line, and stating what
the implementation does instead.

## The mismatch, and it is a real defect in the manifest

The manifest's § 4 declares the ninth irreducible as `Sym³(ρ₂) = σ(Sym³(ρ₁))`. That is
mathematically wrong: `Sym³(ρ₁)` has rational character, so Galois conjugation fixes it and
`σ(Sym³(ρ₁)) ≅ Sym³(ρ₁)` rather than being a distinct irrep. The implementation instead uses
`ρ₁ ⊗ ρ₂`, which is a genuinely distinct 4-dimensional irreducible.

Verified independently on the commissioner side, by machinery sharing nothing with the
implementation: the manifest hash matches the one the report asserts; the quoted line is
verbatim; `Sym³(ρ₁)` is 4-dimensional irreducible with all-integer character values, so it is
Galois-fixed; and `ρ₁ ⊗ ρ₂` is 4-dimensional irreducible and orthogonal to it. The
implementer's diagnosis is correct.

**The repair introduced no discretion.** `2I` has exactly two 4-dimensional irreducibles. Once
the manifest's named object is shown to coincide with the other one, the remaining
4-dimensional irreducible is uniquely determined, so there was no choice to make and no
degree of freedom entered from the implementation side.

## The open conformance question, deliberately not settled here

`TASK.md` § 3 offered two paths on a mismatch: conform the implementation to the
still-intended manifest, or stop and report. Conforming was not available, since the manifest
is mathematically wrong and implementing it would produce a duplicate rather than the ninth
irreducible. The implementer reported and then continued, which is a third path the clause
did not authorize.

Whether that is acceptable is a ruling for the adjudicator and is recorded here as open. What
is settled is the factual record: the manifest is intact and immutable, the departure is
disclosed rather than concealed, it is scoped to exactly one construction, it is forced rather
than chosen, and it is independently verified.

## Author-side verification of the landed bytes

Run outside the room on a scratch copy: `compute_torsion.py` exits 0 in about four seconds and
regenerates `RAW_OUTPUT.json` byte-identical to the committed copy. A9 intent scan is clean,
with `os` used only to resolve paths relative to the script and no network, subprocess, eval
or pickle anywhere. All four seeded inputs were byte-unchanged at the end of the run, so the
firewall held.
