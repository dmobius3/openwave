"""Q4: the integration rehearsal (addendum 12.1.6), the acceptance criterion
the gate matrix named and the burned execution never ran.

    synthetic Packet II in the REAL schema (SYN-, unsealable)
      -> real staged ingestion (hash-verify, parse, canonical check)
      -> real committed-shape Step-3 artifacts from BOTH routes
      -> real Step-3 schema + full-lattice validator BEFORE either adapter
      -> real 3a comparator (compare_3a)
      -> real 3b adapter over the real evaluator (compare_3b)
      -> route-(b) deletion test: recompute physically unavailable, the
         adjudication completes from committed artifacts alone

Case: L(3,1), a section 6.1 pilot tuning case.  Values are pilot-anchored
(TUNING_REFERENCE); nothing here is sealed, evidentiary, or a certification
claim.  Exit 0 only if every stage lands exactly as the addendum requires.
"""

import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
for _p in ("gates", "production", "pilot", "eval3b"):
    sys.path.insert(0, os.path.join(ROOT, _p))
sys.path.insert(0, ROOT)

import numpy as np                                                   # noqa: E402

import ingestion                                                     # noqa: E402
import packet_schema                                                 # noqa: E402
import step3_schema as schema                                        # noqa: E402
import step3_runner as runner                                        # noqa: E402
import route_a_producer as route_a                                   # noqa: E402
import route_b_producer as route_b                                   # noqa: E402
from route_a_twosided import close_pairs                             # noqa: E402
from adjudication_gates import compare_3a                            # noqa: E402
from adapter_3b import compare_3b, RUN_CONFIGURATION                 # noqa: E402
from synthetic_instances import packet_ii_for_case                   # noqa: E402
from lauret_evaluator import p_form_spectrum                         # noqa: E402

OUT = os.path.join(ROOT, "rehearsal")
CASE_ID = "SYN-L31-Q4"
Q, S = 3, (1, 1)
CONFIGURATION_ID = "m8_5b-v1-rehearsal"
CERTIFICATION_SEEDS = 480

# route (b) recompute surface, made physically unavailable for the deletion test
ROUTE_B_FILES = [
    "production/route_b_producer.py", "production/route_b_core.py",
    "pilot/route_b_spectral.py", "pilot/route_b_oneform.py",
    "pilot/levels_routeb.py",
]


def mkq(t, q):
    return [float(np.cos(2 * np.pi * t / q)), float(np.sin(2 * np.pi * t / q)),
            0.0, 0.0]


def build_packet_i():
    inv2 = pow(2, -1, Q)
    al = ((S[0] + S[1]) * inv2) % Q
    be = ((S[0] - S[1]) * inv2) % Q
    return {"case_id": CASE_ID, "family": "lens",
            "parameters": {"q": Q, "s": list(S)},
            "generators": [[mkq(al, Q), mkq(be, Q)]],
            "action_convention": "two_sided",
            "format_version": "m8_5b-packet-1"}


def stage(title):
    print(f"\n== {title}")


def fail(msg):
    print(f"  FAIL  {msg}")
    sys.exit(1)


def main():
    os.makedirs(OUT, exist_ok=True)
    record = {"case_id": CASE_ID, "configuration_id": CONFIGURATION_ID,
              "stages": {}}

    stage("build + validate synthetic packets (REAL schemas, SYN- unsealable)")
    pkt_i = build_packet_i()
    hits = packet_schema.judge(pkt_i, "I") + packet_schema.structural_checks(pkt_i)
    if hits:
        fail(f"Packet I not conforming: {hits}")
    pkt_ii = packet_ii_for_case(CASE_ID, "L(3,1)", 3)
    hits = packet_schema.packet_ii_gate(pkt_ii)
    if hits:
        fail(f"Packet II not conforming: {hits}")
    if not packet_schema.production_seal_refusals(pkt_ii):
        fail("SYN- packet was NOT refused by the production sealing gate")
    print("  ok    both packets conform; sealing gate refuses the SYN- packet")

    raw_i = ingestion.canonical_bytes(pkt_i)
    raw_ii = ingestion.canonical_bytes(pkt_ii)
    sha_i = hashlib.sha256(raw_i).hexdigest()
    sha_ii = hashlib.sha256(raw_ii).hexdigest()
    open(os.path.join(OUT, "packet_I.json"), "wb").write(raw_i)
    open(os.path.join(OUT, "packet_II.json"), "wb").write(raw_ii)
    record["packet_sha256"] = {"I": sha_i, "II": sha_ii}

    stage("real staged ingestion (hash -> parse -> canonical)")
    ing_i = ingestion.ingest(open(os.path.join(OUT, "packet_I.json"), "rb").read(), sha_i)
    ing_ii = ingestion.ingest(open(os.path.join(OUT, "packet_II.json"), "rb").read(), sha_ii)
    print(f"  ok    I  {sha_i[:16]}...  canonical={ing_i.canonical_confirmed}")
    print(f"  ok    II {sha_ii[:16]}...  canonical={ing_ii.canonical_confirmed}")

    stage("both routes run and commit (real producers, real frozen runner)")
    gens = [(np.asarray(u, float), np.asarray(v, float))
            for u, v in pkt_i["generators"]]
    t0 = time.time()
    pairs_a = close_pairs(gens)
    art_a = route_a.produce(pairs_a, f"q4-a-{CASE_ID}", CONFIGURATION_ID,
                            CASE_ID, seeds=CERTIFICATION_SEEDS, adjudication=True)
    print(f"  route (a): {art_a['nodes']} nodes, n_max {art_a['n_max']}, "
          f"{len(art_a['records'])} records, {time.time()-t0:.1f}s")
    t0 = time.time()
    art_b = route_b.produce(gens, f"q4-b-{CASE_ID}", CONFIGURATION_ID,
                            CASE_ID, adjudication=True)
    print(f"  route (b): group {art_b['group_order']}, n_max {art_b['n_max']}, "
          f"{len(art_b['records'])} records, {time.time()-t0:.1f}s")

    written = {}

    def writer(label, raw):
        path = os.path.join(OUT, f"step3_route_{label}_{CASE_ID}.json")
        with open(path, "wb") as fh:
            fh.write(raw)
        if open(path, "rb").read() != raw:
            raise IOError(f"{path}: written bytes differ")
        written[label] = path

    result = runner.run(art_a, art_b, gens, writer=writer)
    record["stages"]["step3"] = {k: v for k, v in result["committed"].items()}
    for label, info in sorted(result["committed"].items()):
        print(f"  committed route {label}: {info['byte_length']} bytes  "
              f"{info['sha256'][:16]}...")

    stage("Step-3 schema + FULL-LATTICE validator, before either adapter")
    arts = {}
    for label, path in sorted(written.items()):
        art = json.loads(open(path, "rb").read().decode("ascii"))
        arts[label] = art
        v = schema.validate_records(art["records"],
                                    schema.expected_cells(art["n_max"]))
        if not v["pass"]:
            fail(f"route {label} artifact failed the full-lattice validator: "
                 f"{v['problems']}")
        print(f"  ok    route {label}: {v['records']} records, lattice complete "
              f"through n_max {art['n_max']}, nulls {v['null_valued_fields']}, "
              f"zeros {v['computed_zero_fields']}")

    stage("rung-3a comparator (real path, packet-band authority)")
    for label, art in sorted(arts.items()):
        res = compare_3a(ing_ii.data, art)
        record["stages"][f"3a_route_{label}"] = res["outcome"]
        if res["outcome"] != "GREEN":
            fail(f"route {label} 3a: {res['outcome']}: "
                 f"{res.get('refusals') or res.get('divergences')}")
        ba = res["band_authority"]
        print(f"  GREEN route {label}: {res['compared_levels']} levels; band "
              f"packet={ba['packet_n_max']} claimed={ba['artifact_claimed_n_max']} "
              f"recomputed={ba['recomputed_from_artifact_cells']}")

    stage("rung-3b adapter over the real evaluator")
    n_max = pkt_ii["indexing_map"]["certified_band"]["n_max"]
    q, s = pkt_i["parameters"]["q"], tuple(pkt_i["parameters"]["s"])
    rows = p_form_spectrum(RUN_CONFIGURATION["p"], q, s, n_max,
                           mapping=RUN_CONFIGURATION["mapping"])
    for label, art in sorted(arts.items()):
        res = compare_3b(rows, art, CASE_ID, n_max)
        record["stages"][f"3b_route_{label}"] = res["outcome"]
        if res["outcome"] != "GREEN":
            fail(f"route {label} 3b: {res['outcome']}: "
                 f"{res.get('refusals') or res['exact_divergences'] or res['coexact_divergences']}")
        print(f"  GREEN route {label}: exact n=1..{n_max} and coexact "
              f"M=2..{n_max} agree with the evaluator totals")

    stage("route-(b) deletion test: recompute physically unavailable")
    held = os.path.join(OUT, "held_out")
    os.makedirs(held, exist_ok=True)
    moved = []
    try:
        for rel in ROUTE_B_FILES:
            src = os.path.join(ROOT, rel)
            if os.path.exists(src):
                dst = os.path.join(held, rel.replace("/", "__"))
                shutil.move(src, dst)
                moved.append((src, dst))
        check = subprocess.run(
            [sys.executable, os.path.join(ROOT, "rehearsal_deletion_check.py"),
             OUT, CASE_ID, str(n_max), str(q), str(s[0]), str(s[1]),
             sha_ii],
            capture_output=True, text=True, timeout=300)
        print(check.stdout.rstrip())
        if check.returncode != 0:
            print(check.stderr.rstrip())
            fail("adjudication did NOT complete with route (b) removed")
        print(f"  ok    adjudication completed from committed artifacts alone "
              f"({len(moved)} route-(b) files were absent)")
    finally:
        for src, dst in moved:
            shutil.move(dst, src)
    record["stages"]["deletion_test"] = {
        "result": "PASS",
        "files_physically_absent": sorted(os.path.relpath(src, ROOT)
                                          for src, _ in moved),
        "mechanism": ("files moved out of the tree; adjudication re-run in "
                      "subprocess rehearsal_deletion_check.py against the "
                      "committed artifacts; files restored afterwards"),
    }

    blob = (json.dumps(record, sort_keys=True, indent=2, ensure_ascii=True)
            + "\n").encode("ascii")
    open(os.path.join(OUT, "REHEARSAL_RECORD.json"), "wb").write(blob)
    print("\nQ4 REHEARSAL: PASS.  Record written to rehearsal/REHEARSAL_RECORD.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
