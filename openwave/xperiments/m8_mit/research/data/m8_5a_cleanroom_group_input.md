# GROUP INPUT

Two generators of a finite subgroup of the unit quaternions, supplied by the model author.

`m8_5a_packet.json` in this directory is the authoritative form. It is what the SHA-256 below
pins, and it is what you should parse. The copy printed at the end of this file is a
convenience for reading, not a second source: if the two ever disagree, the file wins.

| Field | Value |
| --- | --- |
| file | `m8_5a_packet.json` |
| SHA-256 | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| format version | `m8_5a-packet-v1` |

## The coefficient format

Every component is an exact element of the field `Q(phi)`, where `phi` satisfies
`phi^2 = phi + 1`. Components are written `(a + b*phi)/2` with integer `a` and `b`, and the
denominator is fixed at 2 for every component and is not reduced. The quaternion basis is
`(1, i, j, k)`, in that order.

No decimal rendering of these components exists anywhere in this packet. That is deliberate:
exact arithmetic over `Q(phi)` is available to you through `fractions.Fraction`, and working
exactly avoids putting a tolerance where an equality belongs.

## What is deliberately not here

The generators, and nothing else. No group order, no element count, no labels, no dimensions,
no distances, no character values. `PROTOCOL.md` permits the group order to be supplied as a
construction input; it is withheld here on purpose, so that the gate checking it stays a check
that can fail rather than a restatement of something handed to you.

## The packet, verbatim

```json
{
  "coefficient_field": {
    "generator": "phi",
    "minimal_polynomial": "phi^2 - phi - 1"
  },
  "coefficient_form": "(a + b*phi)/2 with integer a and b; the denominator is fixed at 2 for every component and is not reduced",
  "format_version": "m8_5a-packet-v1",
  "generators": [
    [
      "(1 + 0*phi)/2",
      "(1 + 0*phi)/2",
      "(1 + 0*phi)/2",
      "(1 + 0*phi)/2"
    ],
    [
      "(0 + 0*phi)/2",
      "(1 + 0*phi)/2",
      "(-1 + 1*phi)/2",
      "(0 + 1*phi)/2"
    ]
  ],
  "quaternion_basis": [
    "1",
    "i",
    "j",
    "k"
  ]
}
```
