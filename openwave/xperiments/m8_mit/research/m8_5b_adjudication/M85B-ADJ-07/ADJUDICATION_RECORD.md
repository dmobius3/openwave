# M85B-ADJ-07: adjudication record

    Case         M85B-ADJ-07, committed by Addendum 12.4
    Outcome      rung 3a GREEN, rung 3b GREEN
    Sequence     § 4.1 steps 1 through 7, completed in order

## The ordering, checkable from this repository rather than from testimony

    26c68094   Addendum 12.4 adopted and published, carrying BOTH packet
               hashes, while both packets were still sealed and unopened
    6b5abba0   both frozen routes run on the disclosed case, both raw
               outputs committed complete and unrevised, Packet II still
               sealed and unread
    9d96f6e7   rung 3a, after Packet II was hash-verified and revealed
    1dd7451f   rung 3b, against the route outputs already committed above

Each commit descends from the one before it. The pre-reveal commitment is
therefore a fact about commit ancestry, not a claim anyone has to be
trusted about.

## What each gate established

**Commitment before reveal.** The two hashes were published in the
preregistration before either packet opened, and the published bytes were
fetched back from the public repository and both hashes verified in them
before Packet I was opened. A previous case, `M85B-ADJ-06`, was sealed
validly and then retired pre-route precisely because that ordering was not
followed; the barrier that caught it was written after that failure.

**Rung 3a, published-value adjudication.** Both frozen production routes
independently reproduced the complete sealed scalar multiplicity vector
through the packet-owned certified band at exact integer equality. All
five in-band levels participated, INCLUDING the two zero-valued reference
cells, which the contract requires to participate identically to nonzero
ones. No divergence, no cell ignored, no post-reveal reconciliation.
Packet II owned `n_max`; no route metadata could shrink or widen the band.

**Rung 3b, theorem evaluation.** The frozen evaluator ran case-specifically
at the pinned configuration `{p: 1, mapping: "corrected"}`; the
`as_printed` mapping exists only inside the evaluator's mutation battery
and was not used. The frozen adapter owned the acceptance predicate, the
exact-sector comparison at `n = k`, and the coexact aggregation
`T(M) = m_up(M-2) + m_down(M)` over `M = 2..n_max`. Both routes GREEN,
with empty exact and coexact divergence lists and no refusals.

**Reporting asymmetry, permanent and attached.** Against route (a), rung
3b is substantially method-independent. Against route (b), it is a
source-pinned theorem-consistency check that likely overlaps route (b)'s
own character machinery, and it is NOT a third independent derivation of
route (b).

## Post-adjudication publication, a reported gate

Both canonical packet byte streams are published here without
reserialization, copied byte for byte and never regenerated from parsed
objects. Recomputed from the committed files, they reproduce their
pre-reveal commitments exactly:

    packet_I.json    a2ea9172688df7c194ddf221824bf3d3fd69b462d5936f2e8efdd66b1fc4c4f2
    packet_II.json   5fed19674928c2525e0b31476529195ccc88ba5dcba8e493163c87b84e4dbfcf

GATE PASS. A later reader can therefore check the rung-3a comparison
against the very bytes the harness loaded, rather than against the
adjudicator's word.

## Process history carried here rather than in the pull request

**The sealing unit STOPPED on a genuine source ambiguity.** The cited
publication prints its analytic operator convention with one sign, at
`-k(k+2)` in its section 4 and equation (7), and its spectral convention
with the other at the commissioned loci, where section 9 and the Table 1
caption both describe `k(k+2)` as the eigenvalue. No build-time gate could
adjudicate that: V5 is a closure predicate, and a fully sign-reversed
declaration with sign-reversed coefficients satisfies it identically. The
unit ruled out a decoding artifact first, pinning the sign glyph by
mathematics rather than by font tables, then brought both readings and
resolved neither. The resolution came from Addendum 12.1.2 as already
frozen: the field is the sign of the source's PRINTED spectrum, and Packet
II transcribes that printed convention rather than an independent
reconstruction of the operator convention. A later unit re-derived the
same reading independently from the contract before proceeding.

**An earlier case was retired pre-route on a commitment-order failure.**
`M85B-ADJ-06` was sealed validly, and then Packet I was opened while its
commitment addendum was still an unadopted draft. § 4.1 requires both
packet hashes to be published before either packet opens, so the case was
retired with no scientific result, its Packet II left sealed and unread
permanently. No contamination occurred: no route ran and the opening
disclosed only case input that had already been chosen. The commissioner
side then wrote an explicit post-custody barrier, under the stop line
SEALED != CLEARED TO OPEN PACKET I, requiring the commitment to be adopted,
published, and read back from the public record before any packet opens.
That barrier is what this case's ordering satisfies.

## What this case is not

`M85B-ADJ-07` is a HOMOGENEOUS lens space, so its action pair carries
`v = 1`. It cannot satisfy the § 6.1b manufactured-pullback or
accumulation-order discrimination gates, which are pre-freeze
qualification gates for the two-sided architecture, validated before
freeze on configurations capable of discriminating them. In particular
this case does NOT exercise the route-(a) central-equivalence repair
adopted in Addendum 12.3; that repair rests on its own mutation-backed
requalification and independent maintainer attack.

The case is now spent: its identity and its answer packet have been
revealed and adjudicated, so no future blind use of it exists. The § 0
claim ceilings and the M8.5-A gate condition are unchanged.
