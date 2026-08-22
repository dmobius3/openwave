# P1 attachment bundle: the exact bytes P1 read

**HANDLING.** The three source records and both packets copied here are
public and ARE part of the P5 independent reader's declared input set. This
MANIFEST, and every other piece of author-generated packaging, is NOT:
author-side informed review only.

Public does not mean target-free. The source records STATE THE TARGET VALUES
and trigger leak-scan hits. Nothing here is a § 4.3 permitted clean-room
input, and neither are the originals it copies.

Every file below is a byte copy of the artifact P1 used, with its
SHA-256. The two packets are pinned by `adjudicates` and the
protocol sections are the governing text.

| role | file | sha256 |
| --- | --- | --- |
| group | `group__m8_5a_packet.json` | `e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9` |
| construction | `construction__m8_8_construction_packet.json` | `df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06` |
| theory | `theory__torsion-correction.test.py` | `2662e4b56571f0278c7b9223d8682886612b79873640fc4da4bc2ff5e13d033f` |
| note | `note__m8_3_method_note.md` | `3e0c1901d4089991a1de7cff0b1cde453257a29891793249ce63955f144ef06d` |
| repro | `repro__m8_3_mass_reproducer.py` | `cab9a6cd09c9943d21af636f9a1be56bf12b52f4954ecf78907fe5e1718871df` |
| protocol, full | `protocol__m8_8_reproduction_protocol.md` | `8bf24b9fe23e2c2f182c9fabb98a84680fc12a6ab7dec47ccc9a8b3d8bbe306d` |
| protocol sections | `protocol__sections_5_2_to_5_5.md` | `c63ecf2a484d3c1317c36a5799da25aea9f98bdfae6fdf9863fd24252f1be313` |

## Environment when this bundle was built

```json
{
  "mpmath": "1.3.0",
  "numpy": "2.5.0",
  "platform": "macOS-15.7.7-arm64-arm-64bit-Mach-O",
  "python": "3.13.13"
}
```

Recorded, NOT gated: P1 runs green on other runtimes. Provenance,
not a compatibility guarantee.

## External trust root

Source pins live in `../sources.pinned.json`, sha256 `5f7a6a4302ffb275451af204b53ef164b1f3a170f7d894dde735c491208aacbf`. That file is caller-selectable, so P1
cannot itself distinguish an approved repin from a synchronized
source-and-pin replacement. Pin it here and review changes to it.
