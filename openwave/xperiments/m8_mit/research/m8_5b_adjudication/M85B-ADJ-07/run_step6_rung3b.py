"""M85B-ADJ-07, § 4.1 step 6: the frozen 3b evaluator against the ALREADY
COMMITTED route outputs.

Adjudication-side driver only. It invokes the frozen evaluator at the
pinned RUN_CONFIGURATION and hands its rows to the frozen compare_3b,
which owns the acceptance predicate, the exact-sector comparison, and the
coexact aggregation T(M) = m_up(M-2) + m_down(M). This script chooses no
mapping, performs no aggregation, and reconciles nothing.
"""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
TREE = HERE / "../../m8_5b"
for sub in ("gates", "production", "pilot", "eval3b"):
    sys.path.insert(0, str((TREE / sub).resolve()))
sys.path.insert(0, str(TREE.resolve()))

import ingestion                                            # noqa: E402
from adapter_3b import compare_3b, RUN_CONFIGURATION        # noqa: E402
from lauret_evaluator import p_form_spectrum                # noqa: E402

CASE = "M85B-ADJ-07"
PKT_I  = pathlib.Path(sys.argv[1])
PKT_II = pathlib.Path(sys.argv[2])
H_I  = "a2ea9172688df7c194ddf221824bf3d3fd69b462d5936f2e8efdd66b1fc4c4f2"
H_II = "5fed19674928c2525e0b31476529195ccc88ba5dcba8e493163c87b84e4dbfcf"

pkt_i  = ingestion.ingest(PKT_I.read_bytes(),  H_I,  require_canonical=True).data
pkt_ii = ingestion.ingest(PKT_II.read_bytes(), H_II, require_canonical=True).data
assert pkt_i["case_id"] == CASE and pkt_ii["case_id"] == CASE

n_max = pkt_ii["indexing_map"]["certified_band"]["n_max"]
q, s = pkt_i["parameters"]["q"], tuple(pkt_i["parameters"]["s"])
print(f"run configuration (pinned): {RUN_CONFIGURATION}")
print(f"case input from Packet I: q={q}, s={s}; n_max from Packet II: {n_max}")

rows = p_form_spectrum(RUN_CONFIGURATION["p"], q, s, n_max,
                       mapping=RUN_CONFIGURATION["mapping"])
print(f"evaluator produced {len(rows)} rows for k=1..{n_max}")

results = {}
for label in ("a", "b"):
    art = json.loads((HERE / f"route_{label}.step3.json").read_text())
    r = compare_3b(rows, art, CASE, n_max)
    results[label] = r
    print(f"\nroute ({label}): {r['outcome']}")
    if r["outcome"] != "GREEN":
        for x in (r.get("refusals") or [])[:6]:
            print(f"   refusal: {x}")
        for d in (r["exact_divergences"] or [])[:6]:
            print(f"   exact divergence: {d}")
        for d in (r["coexact_divergences"] or [])[:6]:
            print(f"   coexact divergence: {d}")
    else:
        print(f"   exact n=1..{n_max} and coexact M=2..{n_max} agree with "
              f"the evaluator totals; rows ignored above band: "
              f"{r['evaluator_rows_ignored_above_band']}")

outcome = "GREEN" if all(results[k]["outcome"] == "GREEN" for k in results) else \
          "RED" if any(results[k]["outcome"] == "RED" for k in results) else \
          "STRUCTURAL_REFUSAL"
print(f"\nRUNG 3b: {outcome}   (both routes must pass)")
print(f"independence: {results['a']['independence_note']}")
(HERE / "STEP6_RUNG3B_RESULT.json").write_text(
    json.dumps({"case_id": CASE, "rung": "3b", "outcome": outcome,
                "run_configuration": RUN_CONFIGURATION,
                "per_route": results}, indent=2, sort_keys=True) + "\n")
