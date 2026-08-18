"""M85B-ADJ-07, § 4.1 step 5: rung 3a against the revealed reference values.

Adjudication-side driver. It adds no comparison logic: staged ingestion
hash-verifies Packet II, the frozen compare_3a owns validation, the
indexing map, the band and the comparison, and this script only reports
what it returns for each route. No manual transcription, no field
selection, no reconciliation.
"""
import hashlib
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
TREE = HERE / "../../m8_5b"
for sub in ("gates", "production", "pilot", "eval3b"):
    sys.path.insert(0, str((TREE / sub).resolve()))
sys.path.insert(0, str(TREE.resolve()))

import adjudication_gates as adj                            # noqa: E402
import ingestion                                            # noqa: E402

PACKET_II = pathlib.Path(sys.argv[1])
COMMITTED = "5fed19674928c2525e0b31476529195ccc88ba5dcba8e493163c87b84e4dbfcf"

raw = PACKET_II.read_bytes()
ing = ingestion.ingest(raw, COMMITTED, require_canonical=True)
print(f"staged ingestion: hash_verified={ing.hash_verified} "
      f"canonical={ing.canonical_confirmed}")
packet = ing.data
assert packet["case_id"] == "M85B-ADJ-07"

results = {}
for label in ("a", "b"):
    art = json.loads((HERE / f"route_{label}.step3.json").read_text())
    r = adj.compare_3a(packet, art)
    results[label] = r
    print(f"\nroute ({label}): {r['outcome']}")
    if r["outcome"] == "STRUCTURAL_REFUSAL":
        for x in r["refusals"][:6]:
            print(f"   refusal: {x}")
    else:
        print(f"   compared levels: {r['compared_levels']}, "
              f"ignored above band: {r['ignored_above_band']}")
        for d in r["divergences"][:8]:
            print(f"   DIVERGENCE level {d['level']}: "
                  f"reference {d['reference']} vs observed {d['observed']} "
                  f"({d['reference_provenance']})")

outcome = ("GREEN" if all(results[k]["outcome"] == "GREEN" for k in results)
           else "RED" if any(results[k]["outcome"] == "RED" for k in results)
           else "STRUCTURAL_REFUSAL")
print(f"\nRUNG 3a: {outcome}   (both routes must pass)")
(HERE / "STEP5_RUNG3A_RESULT.json").write_text(
    json.dumps({"case_id": "M85B-ADJ-07", "rung": "3a", "outcome": outcome,
                "packet_ii_sha256": COMMITTED, "per_route": results},
               indent=2, sort_keys=True) + "\n")
