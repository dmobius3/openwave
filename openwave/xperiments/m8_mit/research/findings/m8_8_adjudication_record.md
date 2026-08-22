# M8.8 § 8 steps 6 to 9: the adjudication record

> The maintainer-side record required by [`m8_8_reproduction_protocol.md` § 8](m8_8_reproduction_protocol.md)
> steps 6 to 9. Two adjudication attempts ran on the one delivered packet; both are permanent.
> The protocol is frozen ([lock manifest](m8_8_lock_manifest.md)); this file changes nothing
> in it and no § 11 pin. Task record: [`../tasks/m8_8_task_details.md`](../tasks/m8_8_task_details.md).

## Outcome

> **The M8.3 torsion closed forms were reproduced by a context-isolated independent-method
> run, from a based chain complex rather than the spectral-zeta definition. The supplied
> topological model was verified, not independently derived.**

§ 8 category **`convention difference`**, a success category disjoint from `reproduced`: the
§ 5.4 selection at `R7` resolved to the GLOBAL INVERSE, so the committed run's native
orientation is the opposite of the packet's. Under that selection the seven free forms and
`R7` are exactly equal in `Q(φ)` (`row_mismatches = []`), the four identities recomputed from
the selected rows are position-wise equal, and the two sector products cover 8 of 8
nontrivial rows. The § 2 claim ceiling applies as written; `blind` is not used.

**Provenance sentence, fixed in the [#457](https://github.com/openwave-labs/openwave/pull/457) commit message before the official run on `main` and before any category was disclosed** (the maintainer's branch-side run had already executed, since that message reports its byte-identity to the scratch control; the disposition was therefore fixed without knowledge of the category in the public record and on the author's side, not on the maintainer's): the successful
adjudication was obtained on a separately recorded rerun after the initially committed
comparator refused pre-comparison on an exact packet-domain spelling mismatch; the
post-reveal repair changed that frozen-builder literal only and altered neither the
committed reproduction output nor the comparison semantics.

## The two attempts

| | Attempt 1 | Attempt 2 |
| --- | --- | --- |
| Harness | [`scripts/m8_8_adjudication.py`](../scripts/m8_8_adjudication.py) at `ea3452d31f92c599ca0dbbdb4aafde261efc3864` (pre-reveal) | same file at `main` after [#457](https://github.com/openwave-labs/openwave/pull/457), merge `a3ae231dc6481c2127f5f1b1098f919f371185a6` |
| Step 6 | ✅ packet opened, plaintext SHA-256 `744c7f25…` = § 11 pin, 12088 bytes, canonical form unchanged | same bytes, same check, same result |
| Step 7 | ❌ REFUSED at the packet-domain gate before the indexing map was applied: `indexing_map.source_domain 'packet_rows_canonical_position' is not an implemented domain ['answer_packet.rows']` | ✅ map applied, convention map loaded and validated, orientation not yet applied |
| Step 8 | not reached: no pairing, orientation selection, row comparison, identity recomputation or category comparison occurred | ✅ `R7` selection = global inverse; 8/8 rows equal; 4/4 identities equal; sector coverage 8/8 |
| § 8 category | **`structural failure`** (a failed packet gate) | **`convention difference`** |
| Record | [`../data/m8_8_adjudication_attempt1_refusal.json`](../data/m8_8_adjudication_attempt1_refusal.json), 1210 bytes, SHA-256 `46be3d088681b37f28199081f0cfebfe3aa0d57176ba47b2a6c386cc3d121c48` | [`../data/m8_8_adjudication.json`](../data/m8_8_adjudication.json), 9688 bytes, SHA-256 `a60367443cf261f0208b6c827793a84adbb46c1948e0fbaf8774b4b36c32339f` |

**Disposition of attempt 1, as fixed in the [#457](https://github.com/openwave-labs/openwave/pull/457)
commit message before attempt 2 ran on `main`.** First adjudication attempt: STRUCTURAL
FAILURE at the step-7 packet-domain gate. No pairing, orientation selection, row comparison,
identity recomputation, or category comparison occurred. The failure was a comparator
interface defect: the pre-reveal harness used a maintainer placeholder for the packet
source-domain literal ([#453](https://github.com/openwave-labs/openwave/pull/453)), while the
exact replacement literal and its semantics were fixed by the frozen pre-reveal packet
builder. The seam belongs to the interface between the two sides: the
[#451](https://github.com/openwave-labs/openwave/pull/451) answer gave the destination literal
and only prose for the source. A separately recorded rerun was authorized after the narrow
repair; the first failure remains part of the permanent record and is not superseded or
erased.

The ORDERING RECORD of § 3 and § 8 step 8 is the commit graph plus the records above: the
harness landed at `ea3452d3` (and its two earlier readings, [#452](https://github.com/openwave-labs/openwave/pull/452)
and #453) before the packet was opened; the packet's delivery and both hash verifications
are on the [#408](https://github.com/openwave-labs/openwave/pull/408) thread; the pairing and
comparison output is the `rows`, `selection`, `orientation`, `row_mismatches`, `identities`
and `selected_table` members of the attempt-2 record, written by the harness with no manual
transcription anywhere in the path.

## Why a post-reveal comparator repair was acceptable here, and why it is the last one

The protocol does not address repair of the comparator after a post-reveal structural
refusal, in either direction. It is not amended retroactively (§ 12 admits changes only by
dated addendum, and no addendum is filed for this); the gap and the decision are recorded
here instead. The repair was accepted because three facts hold at once, each checkable:

| Fact | Evidence |
| --- | --- |
| The replacement literal was frozen independently before reveal | the author's packet builder pins `INDEXING_SOURCE_DOMAIN = "packet_rows_canonical_position"` and its emitter writes exactly that string, defining it as the packet's canonical row order (author-asserted on #457, 2026-08-22, at `p2_schema.py:205`; checkable from public objects only once the builder bytes publish, and cited as verified then, not before). The delivered packet carries the same string, which the refusal printed |
| Attempt 1 stopped before pairing or any value comparison | the refusal record: `category = "structural failure"`, `refusal` = the domain-gate message, no `rows`, `selection`, `orientation` or `identities` member |
| The patch has mutation-backed semantic nonmovement | #457 self-test 98/98; accepting the old spelling alongside the new turns C3k red (97/98, exit 1); reverting the printer turns C17 red with `KeyError: 'recomputed'` (97/98, exit 1); and the byte-identity control below |

**The byte-identity control, reproducible from public objects.**
[`../scripts/m8_8_control_byte_identity.py`](../scripts/m8_8_control_byte_identity.py) takes
the `ea3452d3` harness bytes from the commit graph, substitutes ONLY the four occurrences of
`answer_packet.rows` by `packet_rows_canonical_position`, runs that scratch copy on the
published packet, and compares the record it writes against the committed attempt-2 record:
BYTE-IDENTICAL, SHA-256 `a6036744…`, 9688 bytes. The scratch copy then exits 1 on
`KeyError: 'recomputed'`, the `ea3452d3` printer defect, after the record is written: the
printer extraction in #457 is downstream of record construction and touches no record field.
This establishes that attempt 2's result is the result of the one-string correction alone. It
does not make attempt 1 retroactively green.

**This was the last ordinary post-reveal comparator repair.** Had the rerun exposed another
packet-facing incompatibility, the run would have stopped rather than landed a third
mechanical fix; an adaptive comparator is exactly what the precommitment exists to prevent.

## Step 9: the published objects, each verified against its frozen hash

Step 9 has two layers. The maintainer layer, published here: the canonical packet and the
frozen construction-audit artifact against their pins, plus the two adjudication records.
The author layer, following before final closeout: the provenance archive plaintext and the
frozen builder behind the `p2_schema.py:205` pin, with its pinned inputs and recipe.

**Maintainer layer (this pull request)**

| Object | Path | SHA-256 | Frozen where |
| --- | --- | --- | --- |
| Canonical answer packet | [`../data/m8_8_answer_packet.json`](../data/m8_8_answer_packet.json), 12088 bytes | `744c7f25e2312d90fc356b11510da685328f05e80ae62721d0a0f418dcf9697e` | § 11 pin; [lock manifest](m8_8_lock_manifest.md); recorded by the author 2026-08-10 |
| Maintainer-side construction-audit artifact, provenance class `derived` | [`../data/m8_8_packet_audit.json`](../data/m8_8_packet_audit.json), 9463 bytes; output of [`../scripts/m8_8_packet_audit.py`](../scripts/m8_8_packet_audit.py) at `e3165904` against the construction packet `df00c022…`: 12 of 12 checks pass, 12 of 12 mutations detected | `d5bb04b9c747d3780a3e931d33d8ed9c7ab79c759a9626ccf664d88b478ef0bb` | [content-commit record](m8_8_content_commit.md) (`b1b6ce48`) + [lock manifest](m8_8_lock_manifest.md); held until now on a maintainer-side orphan commit outside the clean-room base |
| Attempt-1 refusal record | [`../data/m8_8_adjudication_attempt1_refusal.json`](../data/m8_8_adjudication_attempt1_refusal.json) | `46be3d088681b37f28199081f0cfebfe3aa0d57176ba47b2a6c386cc3d121c48` | the 8-hex prefix on the [#408 thread](https://github.com/openwave-labs/openwave/pull/408) (2026-08-22 19:15:29Z, 24 s AFTER #457 was opened); the full digest is first published here. The record itself was written before the repair existed, which the commit graph cannot show and this file does not claim |
| Attempt-2 adjudication record | [`../data/m8_8_adjudication.json`](../data/m8_8_adjudication.json) | `a60367443cf261f0208b6c827793a84adbb46c1948e0fbaf8774b4b36c32339f` | not frozen anywhere before this file: the #457 commit message asserts byte-identity to the scratch control without a digest. What binds it is reproducibility, the rerun recipe and the control below, from public objects |

**Author layer (published in [#459](https://github.com/openwave-labs/openwave/pull/459))**

| Object | Frozen where | Published |
| --- | --- | --- |
| Provenance archive plaintext (ciphertext `2ba72660…`, tag `m8.8-provenance-02`) | hash on the [#408 thread](https://github.com/openwave-labs/openwave/pull/408), verified on decryption; [content-commit record](m8_8_content_commit.md) | [`../data/m8_8_provenance/`](../data/m8_8_provenance/): tarball `4fa0228b…`, 59948 bytes, its 20 members extracted beside it; verified `--strict` at review |
| Frozen packet builder behind the `p2_schema.py:205` pin, with pinned inputs and recipe | author-asserted on #457; no pre-reveal hash of the builder bytes exists, the § 11 packet pin is the pre-reveal anchor | [`../m8_8_answer_builder/`](../m8_8_answer_builder/): P1 regenerates `93acd837…` and P2 regenerates `744c7f25…` byte for byte from the pinned sources, reproduced at review on a different runtime |

With both layers published the closeout is complete and the roadmap row reads `DONE`.

Pins consumed by the harness, all unchanged from § 11 and Addendum 1: answer packet
`744c7f25…`, raw output `1a9b56ce…`, method-and-gate manifest `8aa140e3…`, group packet
`e3b0c945…`, construction packet `df00c022…`.

## Rerun recipe

```text
python3 openwave/xperiments/m8_mit/research/scripts/m8_8_adjudication.py \
    --packet openwave/xperiments/m8_mit/research/data/m8_8_answer_packet.json \
    --json /tmp/m8_8_rerun.json          # must hash a6036744…, exit 0
python3 openwave/xperiments/m8_mit/research/scripts/m8_8_control_byte_identity.py   # exit 0
```
