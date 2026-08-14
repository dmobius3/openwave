"""The deletion-test limb of Q4, run in a SUBPROCESS while the route-(b)
recompute surface is physically absent.

Re-ingests the synthetic Packet II by hash, loads both COMMITTED artifacts,
and runs the real 3a comparator and the real 3b adapter over the real
evaluator.  If any of that needed a route module, the import or the run dies
here and the parent records the failure.  Exit 0 only on GREEN everywhere.

argv: OUT_DIR CASE_ID N_MAX Q S1 S2 SHA_II
"""

import json
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
for _p in ("gates", "production", "eval3b"):
    sys.path.insert(0, os.path.join(ROOT, _p))
sys.path.insert(0, ROOT)

import ingestion                                                     # noqa: E402
from adjudication_gates import compare_3a                            # noqa: E402
from adapter_3b import compare_3b, RUN_CONFIGURATION                 # noqa: E402
from lauret_evaluator import p_form_spectrum                         # noqa: E402


def main():
    out, case_id, n_max, q, s1, s2, sha_ii = sys.argv[1:8]
    n_max, q, s = int(n_max), int(q), (int(s1), int(s2))

    ing = ingestion.ingest(
        open(os.path.join(out, "packet_II.json"), "rb").read(), sha_ii)
    rows = p_form_spectrum(RUN_CONFIGURATION["p"], q, s, n_max,
                           mapping=RUN_CONFIGURATION["mapping"])
    for label in ("a", "b"):
        art = json.loads(open(os.path.join(
            out, f"step3_route_{label}_{case_id}.json"), "rb").read().decode("ascii"))
        r3a = compare_3a(ing.data, art)
        r3b = compare_3b(rows, art, case_id, n_max)
        print(f"    [deletion subprocess] route {label}: "
              f"3a {r3a['outcome']}, 3b {r3b['outcome']}")
        if r3a["outcome"] != "GREEN" or r3b["outcome"] != "GREEN":
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
