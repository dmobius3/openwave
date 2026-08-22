# M8.8 answer-packet builder, frozen bytes

> Published after § 8 step 9 ([#458](https://github.com/openwave-labs/openwave/pull/458)).
> Until this directory existed, the adjudication record cited the builder's source-domain
> pin as author-asserted. It is now checkable from source, and so is the packet itself.

## What the pin is

`p2_schema.py` line 205:

```python
INDEXING_SOURCE_DOMAIN = "packet_rows_canonical_position"
```

and `p2_emit.py` writes exactly that string into `indexing_map.source_domain`, defining it
as the packet's canonical row order. These bytes were frozen at issuance on 2026-08-10 and
are the reason the post-reveal comparator repair in
[#457](https://github.com/openwave-labs/openwave/pull/457) was an interface-literal repair
determined before reveal, not an accommodation to the revealed packet.

## Reproduce the packet from the pinned sources

The chain has two stages. P1 cross-checks the nine values, two ratios and two products
against three pinned M8.3 sources and emits a canonical artifact; P2 emits the packet from
that artifact under the frozen schema. Both are deterministic; no timestamps enter either
output.

```bash
python3 p1.py --group bundle/group__m8_5a_packet.json \
  --construction bundle/construction__m8_8_construction_packet.json \
  --theory bundle/theory__torsion-correction.test.py \
  --note bundle/note__m8_3_method_note.md \
  --repro bundle/repro__m8_3_mass_reproducer.py \
  --source-pins sources.pinned.json --emit /tmp/p1
python3 p2.py --p1-artifact /tmp/p1/p1_crosscheck.json --emit /tmp/p2
```

| Stage | Output | SHA-256 | Bytes |
| --- | --- | --- | --- |
| P1 | `p1_crosscheck.json` | `93acd8376da92687626cb6715aae7a5cd35c8adbb8c9d3eb7a0fd2ee006b3df4` | 14515 |
| P2 | `m8_8_answer_packet.CANDIDATE.json` | `744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e` | 12088 |

The P2 value is the § 11 commitment and the hash of the packet published in
`../data/m8_8_answer_packet.json`. The frozen P1 artifact this directory carries in `out/`
is the one P2 consumed at issuance; a fresh P1 run reproduces it byte for byte.

## Contents

| Path | Role |
| --- | --- |
| `p1*.py`, `qphi_exact.py` | P1: source extraction, exact `Q(φ)` arithmetic, signatures, mutation suite |
| `p2*.py` | P2: frozen schema, emitter, mutation suite |
| `p3*.py`, `p4*.py`, `p6.py` | the later stages: self-test, reference and self-check, closing checks. Not needed to regenerate the packet; part of the frozen record |
| `sources.pinned.json` | SHA-256 pins of the three M8.3 sources P1 reads |
| `bundle/` | the frozen inputs exactly as consumed, the three pinned sources among them |
| `out/p1_crosscheck.json` | the frozen P1 artifact |
| `out2/m8_8_prelock_hash_record.IMMUTABLE.json` | the issuance record, 2026-08-10 |

Machine-local runtime logs and bytecode caches are not included; each log declares itself
noncanonical. The builder stages were adversarially reviewed on
[#408](https://github.com/openwave-labs/openwave/pull/408) before issuance; the review
record is `P2_REPORT.md` on the author side and is summarized in the thread.
