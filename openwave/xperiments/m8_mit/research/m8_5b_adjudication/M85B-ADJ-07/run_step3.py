"""M85B-ADJ-07, § 4.1 step 3: run both frozen routes, commit raw outputs, stop.

Adjudication-side driver. It adds no logic: Packet I in, the two frozen
producers at their frozen defaults, the frozen step-3 runner's route-local
validation, canonical bytes out. Cross-route questions belong to step 5 and
are not asked here.
"""
import hashlib
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
TREE = HERE / "../../m8_5b"
sys.path.insert(0, str((TREE / "production").resolve()))
sys.path.insert(0, str((TREE / "gates").resolve()))
sys.path.insert(0, str((TREE / "pilot").resolve()))

import numpy as np                                          # noqa: E402
import route_a_producer as route_a                          # noqa: E402
import route_b_producer as route_b                          # noqa: E402
import step3_runner as runner                               # noqa: E402
from route_a_twosided import close_pairs                    # noqa: E402

# the rehearsal's certification resolution: top rung of the frozen ladder
CERTIFICATION_SEEDS = 480

PACKET_I = pathlib.Path(sys.argv[1])
COMMITTED = "a2ea9172688df7c194ddf221824bf3d3fd69b462d5936f2e8efdd66b1fc4c4f2"

raw = PACKET_I.read_bytes()
got = hashlib.sha256(raw).hexdigest()
assert got == COMMITTED, f"Packet I hash {got} != adopted commitment"
packet = json.loads(raw)
assert packet["case_id"] == "M85B-ADJ-07"
pairs = [(np.array(u, dtype=float), np.array(v, dtype=float))
         for u, v in packet["generators"]]

pairs_a = close_pairs(pairs)
art_a = route_a.produce(pairs_a, "adj07-step3-a", "m8_5b-v1", "M85B-ADJ-07",
                        seeds=CERTIFICATION_SEEDS, adjudication=True)
art_b = route_b.produce(pairs, "adj07-step3-b", "m8_5b-v1", "M85B-ADJ-07",
                        adjudication=True)

def writer(label, raw_bytes):
    (HERE / f"route_{label}.step3.json").write_bytes(raw_bytes)

summary = runner.run(art_a, art_b, pairs, writer=writer)
(HERE / "STEP3_RUNNER_SUMMARY.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(json.dumps(summary["committed"], indent=2, sort_keys=True))
print("step 3 complete; stopped before any cross-route question")
