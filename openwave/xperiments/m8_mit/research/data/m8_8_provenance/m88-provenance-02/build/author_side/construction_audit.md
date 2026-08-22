# Construction-audit artifact (author side), candidate 2

**Status: outside the room, delivered to the maintainer.** This is the AUTHOR-side account
of the construction. It is delivered to the maintainer as input to the maintainer-side
construction-audit artifact that § 4.2 check 6 names, which is the maintainer's own held
record and a different document from this one. Both stay out of the clean room until
commitment and are published with the reveal; this file is not committed to the repository
and not posted publicly before then. The clean-room implementer never sees it.

| | |
| --- | --- |
| packet | `m8_8_construction_packet.json`, 2657 bytes, canonical form |
| packet SHA-256 | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| derivation source SHA-256 | `8a3a1c87f54372a446356a5c2a5ece4d9b4ba7a32367ef129b8baf18b44733f6` |
| provenance ID | `M88-CONSTR-02` |
| provenance class | `derived` (§ 4.2) |
| supersedes | `M88-CONSTR-01`, packet `2b51ce55…`, REJECTED |

## 1. What the packet is

A finite based 3-dimensional chain complex over `Z[2I]`, free ranks `[1, 2, 2, 1]` in degrees
0 through 3, representing `S³/2I`. The whole object is 20 group-ring terms: `d1` 4, `d2` 12,
`d3` 4.

## 2. Candidate 1 was rejected, and why that matters here

The first packet passed the author-side gate set 31/31 and its own mutation battery 8/8, and
was **wrong**. The maintainer-side audit (`research/scripts/m8_8_packet_audit.py`, upstream
`87cf21a9`) failed it on A7: read as a complex of free `Z`-modules, `C_*` must be the
cellular chain complex of the universal cover `S³`, so `H_*(C_* ⊗ F_p)` must be `(1,0,0,1)`
at every prime. It was `(1,0,6,7)` at 2, `(1,0,14,15)` at 3, `(1,0,16,17)` at 5,
`(1,0,11,12)` at 7, and correct only at large primes. `im ∂₃` sat at finite index inside
`ker ∂₂`, so `H₂(C_*)` carried torsion that `S³` does not have.

Root cause: the acceptance predicate accepted on rank 119 at a single prime `2⁶¹ − 1`. That
is answer-independent, which is what it was defended on, and it is **not sufficient**: it
certifies rational generation and is structurally blind to a finite index. Those are two
separate requirements and only the first was checked.

The author-side gates could not have caught it. `∂₃ → k·∂₃` leaves `∂∂ = 0`, `ε(∂₃) = 0`,
the augmented homology and per-irrep acyclicity all unchanged, because each is a rank
statement or is computed from `ε(∂₂)` alone. Confirmed directly: `∂₃ → 2·∂₃` passed all 31.

The rejected bytes are retained by the author outside this archive, with their own
rejection record. They are the evidence that the maintainer-side audit was independent and
effective, and they are not superseded quietly.

Nothing in THIS archive verifies the candidate-1 figures above (the 31/31, the 8/8, the
per-prime Betti profile, packet `2b51ce55…`). They are recorded as history; the bytes that
would verify them are retained outside the archive.

## 3. Provenance: derived, and what changed from candidate 1

**Forced, no freedom.** `∂₁` and `∂₂` are determined by Fox calculus once the presentation,
relator ordering and basing are fixed. The presentation is the geometric one,
`⟨s, t | s³ = t⁵ = (st)²⟩`, written as the symmetric relator pair `s³(st)⁻²` and `t⁵(st)⁻²`.

Candidate 1 used the asymmetric pair `s³t⁻⁵` and `t⁵(ts)⁻²`. **An earlier version of this
artifact claimed that pair yielded zero saturated generators in 4000 trials while the
geometric pair yielded one immediately. That comparison was invalid.** The two searches
tested different candidate spaces: the symmetric search enumerated the 119 saturated kernel
basis vectors first and hit at basis 6, while the asymmetric search only ever tried random
combinations and never tested a basis vector. Tested identically, the asymmetric pair yields
saturated generators at basis 7, 12 and 14, and that complex passes the universal-cover gate
at every prime with `ε(∂₂) = [[3,−5],[−2,3]]`, determinant `−1`.

How rich is each kernel in generators? The honest answer requires one care: the count "how
many BASIS vectors are saturated generators" depends on which basis, and the maintainer's
independent rebuild uses a different kernel basis than ours, so his counts and ours are counts
of different sets and are not comparable digit for digit.

Over THIS archive's own basis (the `left_kernel` elimination the frozen `search_sym.py`
performs), classified with the exact certificate and none left inconclusive: **24 of 119 on
the shipped symmetric pair, 46 of 119 on the rejected asymmetric one**, first hits at basis
6, 7, 13, 17 and 7, 12, 14, 16 respectively. The maintainer's rebuild reported 23 and 50 over
his basis. An earlier version of this artifact repeated his figures as if they were counts
over our basis, without recomputing them; that was the same verify-before-repeating failure
this project keeps writing rules about, and it is corrected here from our own recount.

The conclusion is basis-robust either way: the asymmetric kernel is roughly twice as rich in
generators (1.9x on our basis), so the 210-term generator candidate 1 carried (229 terms
for its whole complex) was an unlucky
draw from a dense space, not a property of the presentation at all.

So `∂₂` did NOT have to change and A10 could have been a prose fix. The shipped packet is
unaffected and stands: it passes A1 through A11, and the symmetric pair is the standard
geometric presentation, so it is the better object on its own merits. Only the justification
was wrong. Surfaced by the maintainer asking that a mundane explanation be excluded before
it became folklore; his specific hypothesis (an opposite-group convention slip) was not the
cause, but the challenge was.

A consequence worth recording: candidate 1 also failed A10, because its declared
`basis_order` named `t⁵(st)⁻²` while the shipped matrix was the Fox jacobian of `t⁵(ts)⁻²`.
That was going to be repaired as a prose edit, and per the paragraph above that repair WOULD
have sufficed. The rebuild resolved it a different way instead: the declared and shipped
relators now agree, verified independently by re-deriving the Fox jacobian from the packet's
own `basis_order` prose.

**Selected, one choice, and NOT a unique one.** `∂₃` is basis vector 6 of the saturated
integral basis of `ker ∂₂`, accepted as candidate 7 of 6119, inside the deterministic prefix.

The run is deterministic end to end, and the header seed is inert for the accepted result,
for two separate reasons stated precisely. First, the frozen `search_sym.py` queues the 119
kernel-basis vectors BEFORE its `random.Random(20260803)` generates any sparse candidate, and
acceptance happens inside that prefix, so the header seed feeds only candidates that are
never reached. Second, the minor-sampling RNG inside the acceptance test is seeded per
candidate by its index (`random.Random(tried)`), not by the header seed, so even the
accept/inconclusive verdicts are independent of it. Rerunning the file exactly as frozen,
seed constant included, reproduces the canonical packet byte for byte, which is the § 4.2
check 2 obligation; the seed-inertness above is a stronger structural fact established by
reading the frozen code, not by editing it. (At build time, working copies edited to seeds 1
and 999999 also reproduced basis 6; those runs required editing the file and are recorded as
build-time checks, not as something reproducible from the frozen bytes.) The published
`provenance_id` was produced under `p = 10⁹ + 7`. No uniqueness is claimed: `ker ∂₂` is cyclic but not free of rank one, since `∂₃` is
not injective and `ker ∂₃ = Z·N` for the norm element, forced by `ε(∂₃) = 0`. Structurally
acceptable generators need not be unit-related.

**Invariance under the DECLARED basis changes, argued exactly, and covering the TRIVIAL units
only.** A top-cell change by a unit `±g` multiplies the torsion by `det ρ(±g)^{±1}`. Every
element of a finite group has finite order, so `ρ(±g)` has finite order in `GL_d(C)`, so its
determinant is a root of unity and `|τ|²` is exactly unchanged. Basis-independent, and it does
not require the realization to be unitary in the chosen basis.

**Read that scope narrowly, and note what it leaves open.** It does NOT extend to a different
choice of generator. `C₃` has rank 1, so any two generators of `ker ∂₂` differ by `∂₃' = u·∂₃`
for a unit `u` of `Z[2I]/(N)`, and on every nontrivial irrep, which kills `N`, the torsion
moves by `det ρ(u)`. The units are not exhausted by `±g`: `2I` has 9 complex irreps in 7
Galois orbits with every character real, so `rank Wh(2I) = 2` and nontrivial units exist with
no reason for `|det ρ(u)|` to be 1. **Under a `derived` packet the acceptance predicate pins
the COMPLEX, not the reproduced quantity.** Raised by the maintainer, who observed that this
disclaimer previously sat two paragraphs from a passage asserting exact invariance and that an
implementer would read the stronger one. It is now item 6 of § 9's explicitly-not-verified
list, where it is unmissable.

## 4. The saturation certificate

`certify.py` establishes every premise of its own conclusion. No seed, no random choice;
re-running reproduces identical output, and `--emit` reproduces the identical JSON. It does
not write without `--emit`.

- **containment**: `∂₃∂₂ = 0`, `∂₂∂₁ = 0`, `ε∂₁ = 0`, exact over `Z[2I]`;
- **exact ranks, no Smith form at size 240**: `rank_p` is a lower bound for any `p`; the
  chain relations give matching ceilings, `rank ∂₁ ≤ 119` from `ε∂₁ = 0`,
  `rank ∂₂ ≤ 240 − rank ∂₁`, `rank ∂₃ ≤ 240 − rank ∂₂`. Measured lower bounds 119, 121, 119
  meet all three ceilings, so `rank ker ∂₂ = 119`;
- **saturation**: the gcd of all maximal minors is the product of the elementary divisors and
  divides the gcd of any subset. The deterministic column set and dropped row 0 give a
  `119 × 119` minor of determinant `−1`, so every elementary divisor is 1.

Both lattices are saturated of rank 119 and one contains the other, so they are equal. A
saturated sublattice cannot properly contain another of the same rank: the quotient would be
torsion inside a torsion-free group. **No prime list is load-bearing on the accept side.**

The search helper `sat.py` is one-sided in the safe direction. A modular rank drop rejects
definitively; a subset gcd of 1 accepts definitively; failure to reach gcd 1 returns `None`,
never a false rejection. During the search it returned `None` on kernel basis vector 4 with
gcd 19902511 rather than issuing a verdict.

Mutation-tested at build time, 5 of 5 caught, packet restored byte-identical: `∂₃ → 2·∂₃`,
`∂₃ → 6·∂₃`, a flipped `∂₃` sign, a retargeted `∂₂` element ID, and a dropped `∂₁` term. The
first of those is the exact substitution that passed all 31 gates of candidate 1. Like the
candidate-1 figures in § 2, this battery is author-recorded history: no artifact of it ships
in the archive, and nothing here verifies it. What IS verifiable in-archive is the
maintainer suite's ten mutations, transcript at `build/audit_A1_A11_output.txt`.

Two defects of the certificate generator itself are on the record, both found by the
maintainer and both fixed in the shipped `certify.py`. Its emitted JSON once recorded the
containment premises and the expected ranks as literals rather than measured outcomes, so a
red run could emit a file asserting its premises held; every flag in the emitted certificate
is now the boolean the corresponding check returned, plus a `checks_failed` list. And it once
wrote its certificate unconditionally on every run, which is how a mutation-battery run
against a deliberately broken packet left a stale `NOT CERTIFIED` certificate in the tree,
pinned to a packet that never shipped; it now writes only under `--emit`.

## 4b. Maintainer A11, and the int64 guard (2026-08-04, upstream `af178091`)

The maintainer added **A11**, an exact top-boundary certificate, and it agrees with ours
while reaching the conclusion by a DIFFERENT route. A11 reports `ranks [119,121,119]`,
`ceilings [119,121,119]`, `bounds_close: true`, `rank_ker_d2: 119`, `im_d3_saturated: true`.
It adopts the same rank-ceiling argument, then certifies saturation with 119 unit pivots
from a unimodular elimination, where ours exhibits a single `119 × 119` minor of determinant
`−1`. Two independent methods, same verdict. A11's route is also cheaper, needing no
determinant of a large minor.

He also added an **int64 overflow guard** to his `rank_mod_p`. An int64 elimination is exact
only while `(p−1)² < 2⁶³`; above that the products wrap SILENTLY and the routine returns a
WRONG rank rather than raising.

**Where that hazard lives in OUR code, corrected from an earlier version of this section.**
An earlier account here, repeated in the PR, claimed the guard had to live in `certify.py`
because frozen `sat.py` could not carry it. That rationale was false in every clause:
`sat.py` is pure-Python arbitrary-precision arithmetic, uses no numpy and no int64, and
cannot overflow at ANY prime, verified by running its elimination at `p = 2⁸⁹ − 1`. The only
int64 arithmetic that can approach overflow is `certify.py`'s own numpy `rank_p` at the
user-visible prime, so the guard sits in `certify.py` because that is where the exposure is.
(The frozen `search_sym.py` also carries an int64 screen, at frozen primes of at most 37,
where reduced entries keep products under 37² and overflow is unreachable as frozen.) `certify.py` passes the same
prime to `sat.pivot_cols` so the certificate records one `p` throughout; that is uniformity
of record, not protection.

**Our exposure, checked rather than assumed:** every rank in the published certificate ran at
`p = 10⁹ + 7`, so `(p−1)² ≈ 1.0e18` against `2⁶³ ≈ 9.22e18`, inside the bound with 9.2×
headroom. No result was affected.

What survives of the earlier section's general claim is narrower and still worth keeping: a
hash-pinned source set cannot absorb even a strictly safe fix without invalidating itself, so
any future guard or repair must land in files outside the frozen set, and those files then
have to travel with the archive for a rerun to be faithful.

## 5. Verification record

**Maintainer-side audit, A1 through A11: ALL PASS, all mutations detected.** Including A7 at
every tested prime, A10 reporting `[["s^3 (st)^-2"], ["t^5 (st)^-2"]]` with a unique match
per row, and A11's exact top-boundary certificate. The full transcript, regenerated against
the shipped auditor from inside this archive with relative paths, is
`build/audit_A1_A11_output.txt`; its header carries the packet hash it ran against. The
mutation suite includes `∂₃ → 11·∂₃`, which reddens A11 and nothing else, demonstrating that
gate is load-bearing.

**Author-side, independent of that audit:**

- `ε(∂₂) = [[1,−2],[−2,3]]`, determinant `−1`, so unimodular;
- `H_*(C_* ⊗ F_p) = (1,0,0,1)` at 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 10007;
- element IDs `s = 118`, `t = 80`, unchanged from candidate 1 and structurally incapable of
  moving, since IDs are assigned by lexicographic rank of exact coordinate tuples over the
  group packet's closure and depend on the group packet alone. Rebuilt from scratch to
  confirm rather than argued;
- `basing` fields `module_side`, `vector_convention`, `evaluation`, `augmentation` unchanged;
  only `basis_order` changed, correctly, since it is the field that names the relators.

## 6. What this does not establish

That the packet is THE intended model rather than A model with the right invariants. A7
narrows that gap considerably and the § 4.2 audit is what carries the rest. § 1 branch four
is retained, not discharged.

The sparsity of the rebuilt `∂₃`, 4 group-ring terms against candidate 1's 210 (whole
complex: 20 against 229; `∂₁` and `∂₂` are 4 and 12 here against 4 and 15 there), is
descriptive corroboration that this is a geometric 3-cell rather than a dense kernel
combination. It is not formal evidence and is not offered as any. The formal evidence is the
unimodular minor, the chain relations, and the rank equality.
