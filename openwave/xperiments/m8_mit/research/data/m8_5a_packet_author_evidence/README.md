# M8.5-A clean-room packet (author-side, local staging)

The group input the M8.5-A protocol § 4 assigns to the author, plus the author-side
verification evidence. Local only; nothing here is in git.

## What goes to the maintainer, and in which role

| File | Role |
| --- | --- |
| `m8_5a_packet.json` | **the packet.** Generators and format definition only. This is what the maintainer audits, canonicalizes, and hashes, and it is the only file that enters the clean room (alongside the frozen protocol itself) |
| `verify_packet.py` | author-side verification script, exact arithmetic, with a runnable mutation suite |
| `verification_report.json` | gate outcomes, packet hash, environment. Audit evidence, NOT part of the packet |

The packet deliberately contains no irrep labels, dimensions, McKay distances, character
values, worked examples, or any target value. It is generators and the field convention,
nothing else.

## Format

Coefficients live in `Q(phi)` with `phi^2 = phi + 1`, written canonically as
`(a + b*phi)/2` with integer `a` and `b`, the denominator fixed at 2 for every component and
not reduced, matching the maintainer's stated form. Quaternion basis
is `(1, i, j, k)`. Exact throughout: no decimal approximation appears in the authoritative
packet, so there is no tolerance, no dedup ambiguity, and no closure ambiguity, and the
packet hash is meaningful. A numerical array may be generated downstream from the exact
packet, but must not be separately authored or hashed as a competing source of truth.

## Verification, as run

```
python3 verify_packet.py --packet m8_5a_packet.json --mutation-tests
```

| Gate | Result |
| --- | --- |
| P1 every supplied generator has quaternion norm exactly 1 | PASS |
| P2 the generators close under exact multiplication to exactly 120 distinct elements | PASS |
| P3 identity present, every inverse in the set, center exactly `{+1, -1}` | PASS |
| P4 the quotient by the center has order 60 and the `A_5` element-order profile (1, 15, 20, 24 at orders 1, 2, 3, 5; none at 4 or 6) | PASS |

Two generators, of orders 6 and 4. Packet SHA-256
`e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9`.

Mutation suite: 4 mutations, every gate reddened by at least one, coverage complete.
`quotient_by_trivial` reddens P4 and nothing else, which is what shows P4 carries
independent weight rather than riding on P2.

## Why P4 exists, and the limit of what the mutations demonstrate

Order and closure do not identify the group. The finite subgroups of `SU(2)` include a
binary dihedral group of order 120 and a cyclic group of order 120, so a mistranscribed
generator set could in principle close to exactly 120 elements and still not be `2I`. P4 is
the gate that separates them, since the central quotient of `2I` is `A_5` while the binary
dihedral quotient is dihedral.

**Stated honestly:** the mutation suite demonstrates that P4 reddens when the quotient
construction is corrupted. It does NOT exhibit a genuine competing order-120 group that
passes P2 and fails P4, because such a group is not constructible in this packet's field: a
cyclic or binary dihedral subgroup of that order needs `cos(pi/60)` or `cos(pi/30)`, whose
degrees over `Q` are 16 and 8, so neither lies in the degree-2 field `Q(phi)`. That is an
argument for why the competing case cannot arise here, not a mutation-demonstrated
discrimination, and it is recorded as such rather than as a stronger claim.

## Scope

This certifies the packet's group input. It computes no characters, no irreducible
dimensions, no McKay data, and no first-occurrence anything, and it says nothing about any
downstream table.
