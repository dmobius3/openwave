# § 14 supplement: the two REDs and one exposed bug, resolved

Scope per F1: the original one-pass package stands untouched; this supplement carries the
scoped rework only. Items 1, 2, 3, 7, 8 are unchanged from `DISPOSITIONS.md`; items 4, 5, 6
are superseded by the records below.

## Provenance reconciliation (F2's hash mismatch)

No text moved after § 14 ran. `ed712eec…` is the protocol copy F2 last received (sent at
the fresh-context close); `5b164fc1…`, the recorded pre-§ 14 hash, is that text plus
exactly the three fixes F2 itself ordered in the same review (clause-(f) indexing and the
`d_c`/`gap` generalizations), recorded in the hash chain appended to `SEEDED_CONTROL.md`
and `R2_STOPPING_RULE.md` after that send. Item 5's references were computed under
`5b164fc1…`, the current definitions. The chain since:
`5b164fc1…` → `4d6305e9…` (the § 7 `φ` amendment below) → `b88fa741…` (arm C of the
amendment caught a sign error in the amendment's own closed form; fixed). Current true
protocol: `b88fa741cc2940cb27518a32d358b19eb93dbbbf426ef0aae7e1195cd634abde`.

## Item 6: the § 7 amendment and its scoped verification, GREEN

The amendment replaces the equal-weight pin with the recorded complex candidate, every
convention explicit: weights DESCENDING `+6..−6`, `φ ∝ Σ_{j=1..13} (1 + j/13) e^{ij/3} v_j`
unit-normalized, `J₊/J₋/J₃` and `T₁/T₂/T₃` matrices pinned, `M_a = Re⟨ψ̇, T_aψ⟩`. F2's
reproducibility complaint is answered by construction: the amendment carries the index
range, the weight ordering, the normalization, and the exact expression, plus TWO analytic
identity arms of this specific pin: `|M₁/M₂| = tan(1/3)` (F2's own rotation identity,
adopted as an arm) and the `M₃` closed form. Verification (item 6 alone, per F1):
`M = [+2.093976e+01, +6.047521e+01, −1.753244e+01]`, all three live; arm R matches
`tan(1/3)` to 12 digits; arm C matches the closed form to `1e-10`, AFTER first firing on
the amendment draft itself, whose closed form carried a sign error, fixed and chained
before anything froze. `raw/item6_v2_amended.json` carries `protocol_sha`.

## Item 5: the ν-pairing bug, found by F2's field question, fixed, GREEN

F2 asked what `f64_pos_diff` measures; the honest answer was a bug. The 50-dps eigenvalue
path split blocks by SIGNED `ν`, but the conjugate (`Ju`) coupling connects `ν` to `−ν`,
so every `ν ≠ 0` block dropped its pair couplings and the shipped references were wrong at
the `1e-4` level, which is exactly what those fields were showing. v3 blocks on `|ν|`
unions; the diffs collapse from `5.6e-4` to `2.4e-15`, eigensolver rounding, and a new arm
B freezes that requirement. The references are regenerated
(`raw/controlA_references.json`, sha `0544a3b9…`) and now carry a `_field_semantics` block
distinguishing the three comparisons (50-dps record; cross-EIGENSOLVER diagnostic; arm X's
cross-ASSEMBLY comparison) and the `protocol_sha` whose definitions produced them.

## Item 4: the D-series, restored over the full classification, GREEN

F2's RED was exact: the enumeration covered A and E and silently omitted the entire
binary-dihedral D-series, in a McKay program. v3 enumerates the FULL closed-subgroup
classification of SU(2): cyclic, binary dihedral `BD_n` (n = 2..2l, all four characters
from full-commutator abelianizations), 2T, 2O, 2I, torus classes, and SU(2) itself, with
maximality decided by NUMERICAL projector containment rather than an abstract table, and
the degenerate `l = 0` case excluded outright (in the trivial rep every stabilizer is the
full group). Verified against F2's independent numbers: `Q8` at `l = 12` has dimension 4
by BOTH routes, matching the pinned analytic `(13+13+6)/8`; 2T gives 2 and 2O gives 1
there; `BD₅(χ₁)` and `BD₆(χ₁)` carry the 1-dim lines F2 predicted; 40 BD classes at
`l = 12` in all. A new ARM D tests COMPLETENESS, the axis F2 noted no arm covered:
disabling the BD family must drop classes and the pinned Q8 value must reproduce.
`raw/lattice_tables.json`, sha `24b1f741…`, supersedes both earlier versions.

## Item 3, scope note added where the margin is quoted

The 11x margin covers the three components the obligation names: the § 4.3 agreement gate,
the § 4.2 monitor, and Control B. The sector-basis build across nine sectors, gate 2 over
`n ≤ 3N`, gate 5, and the time arm sit OUTSIDE the rehearsal, and the ceiling remains
STOP-QUAL for the whole attempt regardless.

## Development record

Item 4 v3's first run failed its own rep checks at `l = 0` (the degenerate case above),
guarded and re-run; earlier dev artifacts preserved (`raw/s14_item4_dev2.py`,
`raw/lattice_tables_dev1.json`). Item 5's v2 is preserved as `raw/s14_item5_dev2.py`.
Item 6 v2's arm C fired on the amendment's sign error before anything froze; that is the
arm's job and the chain records it.

## Addendum, F2's scoped-close findings, both resolved

**Item 4, the `l = 1` `Z2` contradiction.** The maximality machinery's blind spot was
structural: numerical containment tests LISTED projectors, but the absorber of a generic
`C²` vector is a CONJUGATE torus, invisible to that test. The fix is the acts-by-scalar
rule: a class whose fixed space is ALL of `V_l` (here `π_l(−1) = −I` at odd `l`, so the
twisted `Z2` fixed space is the whole space at `l = 1, 3, 5, 7`) constrains nothing and is
absorbed. At `l = 1` the entry now reads exactly as specified: `orbit_dim 3`,
`maximal_generic False`, `absorbed_by W(|2m|=1)` up to conjugacy; at `l = 3, 5, 7` the
same rule fires with the generic-action wording. A new ARM E freezes the derivation's own
corollary as an assertion: no maximal finite-family class may exist at `l = 1`, and no
acts-by-scalar class may ship as maximal anywhere. Table superseded:
`raw/lattice_tables.json`, sha `dc99b5ae…`.

**Item 5, the two unmoved quantities.** F2's inference (a route mislabel) was reasonable
and WRONG in an instructive way; the record now proves the true explanation instead of
asserting it. All five quantities ride the 50-dps route. `min_abs_cluster` is bit-identical
across v2/v3 because it sits in the `|ν| = 0` block, the one block the signed-ν bug never
touched, and the file now RECORDS the block membership per `s` (`min_abs_block_nu = 0` at
all three). `leakage` barely moved because it is a property of the cluster SPAN, and the
dropped couplings rotated eigenvectors within that span. The file's `route_per_quantity`
block states the route of each quantity and the reason for each near-invariance.
References superseded: `raw/controlA_references.json`, sha `18414ff0…`.

The protocol was UNTOUCHED by that round; the review chain's last full-file value was
`b88fa741…`.

## Closing the chain at the shipped text

`b88fa741…` is a full-FILE hash of the pre-filing draft and resolves to no file in this
repository, because filing changed the document in exactly three mechanical ways and then
changed what is hashed. What was appended: § 15's freeze-time product hashes (this
package's own artifacts, unknown until § 14 finished), the filing banner replacing the
draft banner, the upstream-main pin, and the boundary marker with its freeze record. What
was NOT touched: no gate, no threshold, no ladder, no definition, no commitment, no
outcome sentence. The pin then moved from the whole file to the FROZEN REGION, every byte
above the boundary, which is the object § 15 and the roadmap both quote:

    cf5faa10401e02ed00bdfb92d8ee329fff5863b108f979d4c90abfe80e2b8138

superseded once more at review by the maintainer's blocking finding below. The shipped
frozen-region digest is the one in the protocol's own freeze record; it is the only hash a
reader needs, and the earlier values in this file and in `DISPOSITIONS.md` are the dated
record of how it got there.

## Maintainer review, the manifest finding, and what it cost

`raw/package_hashes.txt` shipped covering 18 of 30 tracked members, with six orphan hashes
and three mislabels: the four artifacts regenerated after their first hashing (the lattice
tables, the Control A references, and two of the item scripts) each kept a superseded hash
under their own path name, and eight members had no line at all. The cause is precise and
worth recording, because it is a variant of a failure this program has hit before: the
manifest was built by APPENDING `shasum` runs at three different times rather than being
regenerated, so every superseded value stayed under the same path and no member was ever
retired. Recomputing a hash and appending it is not the same as recomputing the manifest.
The manifest is now generated in one shot from `git ls-files` and verifies member for
member with `shasum -c`.

Three comment-tier items rode along, all of them real. `jacobian_check.py` printed two
verdict labels unconditionally, so its gauge-breaking arm could not fail; the reviewer
demonstrated it by zeroing the breaking coefficient and getting the same green. Both labels
are now computed, the mutation is scored against the parent at a `1e3` ratio in the house
pattern, and the script exits nonzero on failure; the fix was itself mutation-tested, the
neutralized arm now printing "arm did NOT fire (broken arm)" and exiting 1. Two § 14
scripts read a scratch `/tmp` file for the protocol hash they stamped into their records,
which is the mechanism that put an unresolvable hash into the shipped JSONs. The naive
repair, hashing the protocol at its repo path, does not converge, since § 15 pins those
same JSONs and each would move the other forever; the records now carry the protocol PATH,
which is stable, while the run records print the whole-file and frozen-region hashes
resolved at run time. Both halves ship; neither is circular.

