# M8.5-C2 change register: every deviation from the superseded frozen text

The draft `M8_5C2_PROTOCOL_DRAFT.md` is DERIVED from the superseded M8.5-C protocol's
frozen region, extracted byte-exact from merged main at `0b7505b9` (the extraction hashes
to the frozen digest `e253558b…`, kept beside this register as `M8_5C2_BASE_FROZEN.md`).
Every deviation is in `M8_5C2_DERIVATION.diff` (18 hunks, 208 content additions, 56 content deletions, file-header lines excluded), and every LINE of every hunk is authorized by
a class below (a diff hunk may merge adjacent deviations from more than one class, so the
authorization unit is the deviation, not the hunk; four hunks carry multi-class content), which are the #501 ruling's item 5 plus the two
mechanical classes any refiling needs. The verification a reviewer runs: diff the base
against the draft, confirm the hunk set matches this register, and confirm the base
recomputes to `e253558b…`. The base ships IN this package as `M8_5C2_BASE_FROZEN.md`, so
both halves run without any other artifact.

| class | authorization | hunks it owns |
| --- | --- | --- |
| IDENTITY | new protocol identity per the ruling | title; the banner rewritten to DRAFT status with the C2 filing path and ten-obligation state; outcome tokens `M8.5-C2-QUALIFIED` / `M8.5-C2-FAILED` (all occurrences); self-references in §§ 0, 1, 2, 3.4, 8 row 5, 11; the § 15 rulings-record row repointed to `M8_5C2_OPEN_DECISIONS.md`; the § 16 marker renamed `M85C2-FREEZE-BOUNDARY` |
| LINEAGE | the record must say what this is | § 0 lineage paragraph (supersession named as ruled, A1 pointer, no execution credit, the four-changes commitment itself); § 12 named-OUT additions (superseded protocol, A1 archive); § 12 commissioning NEW-room sentence; § 15 rows for the superseded protocol, the A1 archive, and the two at-filing placeholders; § 15 inheritance note |
| CHANGE 1: arms | both parity-dead node-drop arms replaced | § 4.1 arm block (even-`K` angular drop + cross-parity leakage with its ARM-LOCAL read construction at scalar levels `n ∈ {1,3,5}`, dimensionless scores `E_mut ≥ 1e3 × max(E_parent, 1e-16)` with the floor-coupling safety sentence, u-drop excluded on the structural per-field ground, odd-`K` REGRESSION CONTROL); § 4.2 monitor arm block (even-`K` 4N-point drop + leakage with the same read construction + regression control, injected-content arm unchanged); § 8 rows 4 and 6 cells naming draw AND read spaces |
| CHANGE 2: arm-on-arena | every arm and every design-input record demonstrated on its gate's arena before freeze | § 8 arena column (header + all 11 rows) and preamble rule; § 14 item 9 (the demonstration record, sub-items a through m) and item 10 (the audit, with the regression-control exemption); § 14 preamble and dispositions paragraph (items 1, 2, 3, 8 FRESH per the D-5 ruling, item 5's inheritance made MECHANICAL via the § 5 definitions-block digest, which recomputes unchanged at `823e9066…` on this draft); § 8 rows 3, 5, 7 naming FRESH re-run records |
| CHANGE 3: PROTOCOL-INVALID | the administrative class | § 11 paragraph, with the maintainer-independent-reproduction requirement and the licenses-nothing sentence |
| CHANGE 4: ledger schema | the three A1 lessons | § 13 schema paragraph: cumulative `cumulative_seconds`, per-gate completion states {COMPLETE, PARENT-ONLY, NOT-REACHED}, every hash names the bytes it covers, attempt identifiers namespaced by protocol identity (`C2-A1` / `C-A1`) |

Two entries forced during the D-7 round, both inside CHANGE 2's authorization: § 8 row
8's stray arm label `broken-symmetry drift` DELETED, because § 7 defines exactly one
mutation (kick-drift symplectic Euler) and the label had no operative definition and no
demonstration anywhere; recorded per the redline ruling as removal of an undefined arm
label required to make the arm-on-arena obligation executable, not as a fifth scientific
change. And § 15's definitions-block row now records the RECOMPUTATION of the sub-region
digest (string anchors, first-occurrence invariant, the exact command assembling its
anchor at run time so the row never contains it), because a digest that cannot be
recomputed by a stranger is an assertion, not a mechanism; the recorded command reproduces
`823e9066…` on this draft and the anchor-uniqueness invariants are verified.

Two further entries from the arena-conformance room (its hardened rerun; the first run's
verdict was VACATED by the seeded control), both CHANGE-2 conformance fixes confirmed
against the true draft: § 8 row 5's arena cell now names each arm's domain separately,
because its C7 arm is a linear port whose § 3.4(d) domain (the `W_ρ`-valued bases at zero
amplitude, all sectors) the cell had not named; and § 4.2 gains an explicit ALLOCATION
clause for the two new arms (production rung spaces, following the carried
injected-content allocation sentence, Control B carrying the live reading and never an
arm target), with § 8 row 6's cell stating that split instead of silently omitting
Control B.

One further entry from the scoped re-ask (Q3): row 5's Spec cell gains § 0. The arena
cell's domain restriction for the two nonlinear arms is CORRECT and comes from § 0's
containment via the Check cell, but the Spec column cited only §§ 2 and 3.4(d), so the
mechanical audit keyed to the named sections reproduced "every E_rho" and correctly
flagged the mismatch. Same authorization class and same reasoning as the gate-8
stray-label deletion: a minimal conforming edit to an inherited cell, forced by CHANGE
2's audit being executable, not a fifth change.

From the escalation reading (both redline units, rows 5 and 6): row 6 GREEN on both
readings, fix 3 holding. Row 5 took two more CHANGE-2 entries. The C7 arm's gate-5
instance is scoped to gate 5's OWN spaces (`E_R0` / `E_R0 ⊗ C²`), because § 3.4(d)'s
all-sectors zero-amplitude split is a statement about GATE 3's port arenas, gate 3's row
already carries that instance, and a gate-5 arm scoped to sectors gate 5 never evaluates
cannot make any gate-5 check red. And the arena cell is rewritten one line per arm, with
§ 14 item 9 gaining sub-item (k): the non-gradient perturbation demonstration against the
Jacobian-symmetry predicate, NEW evidence, because inspection of the pinned
`jacobian_check.py` shows its mutation is the gauge-breaking term against the kernel
predicate, so the previous cell's citation of it for gate 5 was an inferential bridge to
a false conclusion; `jacobian_check.py` leaves row 5 entirely, its records serving § 2's
design inputs and gate 7's arm, where the gate table already cites it.

Closure round of the escalation reading: item 9 now owns ALL THREE gate-5 arms on the
gate-5 arena, (k) non-gradient perturbation against the Jacobian-symmetry predicate, (l)
symmetry-breaking coupling against the equivariance identity, (m) intra-level spectrum
break against the semi-discrete action predicate, with the fresh layer-1 C3/C7 records
explicitly demoted to provenance carrying no arena credit, per § 2's own
no-inherited-credit rule; row 5's cell maps one-to-one to (k)/(l)/(m). And the § 15
recomputation row's justification is corrected to the TRUE invariant, first occurrence of
the end anchor AFTER the start anchor rather than global uniqueness, with BOTH anchors
now assembled at run time inside the recorded command so the row's own text contains
neither; the command reproduces `823e9066…` unchanged, its semantics being identical.

Amendment from § 14 item 9's execution (the C1 item-6 precedent: a pre-freeze obligation
catching a false sentence in the draft it serves), CORRECTED at redline to the
RUNG-RELATIVE form: § 4.2(b)'s read levels move from `{1, 3, 5}` to
`n ∈ {N+1, N+3, N+5}` (both units independently showed the interim literal `{25, 27, 29}`
would be dead again at every rung above 24: at fixed low read level the substituted rule
regains exactness as `N` grows, so the levels must scale with the rung), with row 6's
cell following. The drafted sentence claimed
the substitution's aliased mass lands at the low odd levels; the demonstration measured
`3e-17` there and the arithmetic agrees, because the read integrand at low levels has
degree at most `3N + 5 <= 4N` and the substituted rule integrates it EXACTLY; aliasing
requires the read to exceed the rule's exactness margin, `n > N`. Gate 4's arm at
`{1, 3, 5}` is unaffected (its mutated rule is exact only to `2N`, which the low-level
read already exceeds). CHANGE-1 authorization; the item-9 record carries the RED, the
amendment, and the re-run, per the house pattern.

§ 15's item-9 pin row is FILLED (LINEAGE/CHANGE-2 authorization, the row's own frozen
instruction): script `42ff8da8…`, record `8ce491c2…`, captured stdout `dc54f673…`, with
the stdout labeled a measurement record per the A1 idempotence lesson. The values are the
executed run-5 artifacts of the preserved five-run chain.

Freeze-mechanics fills, the last entries (IDENTITY/LINEAGE authorization, both rows'
own frozen instructions): the DRAFT banner replaced by the filing banner with the C2
review and § 14 closure state; the genesis-provenance paragraph updated to C2's own
lineage with the C1 memo inherited through the carried text; and § 15's upstream-main
row filled with the branch commit `154825a4…`. After these fills the frozen region is
final; the boundary, freeze record, and digest sit BELOW it and are not derivation
content, and the package manifest is generated LAST from the tracked files, the
explicit ordering the C1 append-versus-regenerate defect and this round's
stale-delivery finding both taught.

Nothing outside these hunks moved: §§ 2, 3 (but for the one identity token), 4.3, 4.4, 5,
6, 7, 9, 10 (but for outcome tokens), the ladders, thresholds, controls, containment,
claim ceiling, and both § 1 sentences' operative content are the superseded text verbatim.

Design inputs already on the record for § 14 item 9, all post-stop provenance, none
inherited as numbers: author even-`K` `1.8e-01`; maintainer even-`K` `1.2e-01` / `6.3e-02`
and cross-parity `7.4e-02` to `2.2e-01` (the #501 verification table); the u-drop
non-reproduction (`2.7e-15` at `N = 3` vs `5.5e-03` at `N = 4`) that excluded it.
