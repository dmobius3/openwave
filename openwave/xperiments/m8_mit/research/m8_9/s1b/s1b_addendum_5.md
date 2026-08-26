# S1b addendum 5: the root of trust for q3a's evidence

> **APPEND-ONLY.** The frozen regions of the parent rule and addenda 1 to 4 are UNCHANGED and all
> five still verify. No byte of attempt `q3a` is altered. This addendum publishes one hash and
> imposes one requirement.

## A5.1 The gap

Addendum 4 § A4.3 publishes independent hashes for `target_n12.json` and `target_n20.json`, and
those are the inputs to branch rules 2 through 7. **Rule 1 asks whether any instrument gate
failed**, and its inputs are the gate records, which addendum 4 authenticates only THROUGH
`q3a/OUTPUT_MANIFEST.json`.

That manifest had no independently published hash. So it was an unauthenticated root of trust: an
alteration to a gate record together with a matching alteration to its manifest entry would change
rule 1's input while both independently pinned target records still verified. The chain had no
anchor at its head.

## A5.2 The hash

    q3a/OUTPUT_MANIFEST.json
    f5401c5179d0cde42a8763de175542d3d04aaa316be0329ec30a42cd9d6bc3a4

Computed independently by the author and by the redline from separate copies, and confirmed equal
across `q3a/` and the adjudication room's `evidence/` copy before publication here.

## A5.3 The requirement

**Any adjudication of `q3a` verifies `OUTPUT_MANIFEST.json` against the value above BEFORE using any
record it pins, and STOPS on mismatch.** With the manifest authenticated, the chain of custody is

    frozen addendum (hash published out of band)
      -> q3a/OUTPUT_MANIFEST.json
        -> gate records and target records

**Every evidence file actually used in a branch trace must match its entry in the authenticated
manifest.** An adjudicator may not rely on a record merely because the room contains it. A record
that is present but unverifiable, whether because its manifest entry is absent from the supplied
subset or because it does not match, is reported as unverifiable and is not used; if a rule's input
is unverifiable, the rule cannot be evaluated and the adjudication STOPS.

`QUALIFICATION_NOTE.md` and `run.log` are known to differ from their manifest entries by
post-manifest appending, documented in addendum 4 § A4.3. Both are non-adjudicative: they are
forensic context, and no branch rule takes an input from either.
<!-- ADDENDUM5-BOUNDARY -->

**Freeze record, addendum 5.** SHA-256 covers every byte ABOVE the boundary comment: `7ca059f074059f925a7f231fcb1ac93932e890121ef9abd1975a4df78542b3e5`

```bash
sed '/^<!-- ADDENDUM5-BOUNDARY -->$/,$d' S1B_ADDENDUM_5.md | shasum -a 256
```

The parent rule and addenda 1 to 4 are untouched and verify independently.
