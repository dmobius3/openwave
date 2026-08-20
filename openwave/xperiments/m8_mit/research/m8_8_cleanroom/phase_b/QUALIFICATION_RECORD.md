# M8.8 Phase B Qualification Record

## What was executed

Phase B parsed 19 gate identifiers and their declared mutations from
`METHOD_AND_GATE_MANIFEST.md` section 4 (SHA-256 `8aa140e3...`), then executed
each mutation against scratch copies of the frozen Phase A machinery via
`run_phase_b.py`. No Phase A artifact was modified at any point.

All 13 Phase A artifact SHA-256 hashes were verified before and after
qualification against the table in PROTOCOL.md Addendum 1. Both public packet
hashes were verified against section 11 pins.

### Gates executed

| Gate | Domain | Mutation applied | Baseline | Mutated | Red |
|------|--------|-----------------|----------|---------|-----|
| G-M01 | dd=0 at degree 2 | d3[0][0] +1*e | zero | nonzero | Yes |
| G-M02 | dd=0 at degree 1 | d2[0][0] +1*e | zero | nonzero | Yes |
| G-M03 | Free ranks / chi | ranks[3]: 1->2 | chi=0 | chi=-1 | Yes |
| G-M04 | Augmented homology | d2[0][0] +2*e | det=-1 | det=5 | Yes |
| G-M05 | Universal-cover saturation | d2 row 0 x2 | \|det\|=1 | \|det\|>>1 | Yes |
| G-M06 | Augmentation terminal | non-augmentation | eps_d1=0 | eps_d1!=0 | Yes |
| G-M07 | Generator correspondence | swap s,t IDs | relators hold | relators fail | Yes |
| G-M08 | Per-irrep acyclicity | V1 M3 row zeroed | rank=d | rank<d | Yes |
| G-T01 | Unitarity | rho(s)[0][0] +1/10 | invariant | not invariant | Yes |
| G-T02 | Row signature | swap chi(s),chi(t) | distinct | changed | Yes |
| G-T03a | Convention: eval map | g->rho(g^-1) | dd=0 | dd!=0 | Yes |
| G-T03b | Convention: boundary dir | cochain reversal | dd=0 | dd!=0 | Yes |
| G-T03c | Convention: module side | rho(g)^T | dd=0 | dd!=0 | Yes |
| G-T03d | Convention: vector conv | GR transpose | dd=0 | dd!=0 | Yes |
| G-D01 | Twisted dd=0 | M3[0][0] +1 | zero | nonzero | Yes |
| G-D02 | Twisted ranks | M3 row zeroed | rank=d | rank<d | Yes |
| G-D03 | Det sub-matrices | minor col zeroed | det!=0 | det=0 | Yes |
| G-D04 | Galois consistency | sigma(T2_V2) for V1 | match | mismatch | Yes |
| G-D05 | Code-path dependency | identity reps | T2=8+12phi | non-acyclic | Yes |

## What this establishes

1. **All 19 preregistered mutations reddened.** Each declared mutation, applied
   to a scratch copy of the frozen Phase A machinery, caused its gate predicate
   to fail. No mutation was skipped, substituted, or narrowed.

2. **Exact registry coverage.** The set of 19 executed gate identifiers equals
   the set of 19 gate identifiers parsed from the frozen manifest. No missing
   record, no duplicate, no unregistered gate.

3. **Phase A integrity preserved.** All 13 artifact hashes match the protocol
   table both before and after qualification. Phase A was not altered.

4. **Machine-readable record.** `MUTATION_RESULTS.json` carries per-gate
   records with: gate_id, object_mutated, gate_predicate, baseline_result,
   mutated_result, and red_outcome.

## Disposition

Phase B qualification is complete. The pre-reveal gate contract from section 9
is satisfied. Outcome determination proceeds under section 8.
