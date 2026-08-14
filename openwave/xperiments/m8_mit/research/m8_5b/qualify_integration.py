"""Q1/Q2/Q3/Q5: the integrated qualification battery (addendum 12.1.6).

Every line printed PASS is demonstrated to fail in its REQUIRED LAYER by a
deliberate defect: comparison RED, validator reject, or STRUCTURAL REFUSAL are
distinct outcomes and a red in the wrong layer is displaced verification.
Consumes the committed Q4 rehearsal artifacts (run rehearsal_q4.py first).
Exit is nonzero if any item misses, any layer is wrong, or any predicate goes
unexercised.
"""

import copy
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
for _p in ("gates", "production", "pilot", "eval3b"):
    sys.path.insert(0, os.path.join(ROOT, _p))
sys.path.insert(0, ROOT)

import numpy as np                                                   # noqa: E402

import ingestion                                                     # noqa: E402
import packet_schema                                                 # noqa: E402
import step3_schema as schema                                        # noqa: E402
import route_b_producer as route_b                                   # noqa: E402
from adjudication_gates import (compare_3a, G5a_map_consumption,     # noqa: E402
                                G9_transcription_mutation)
from adapter_3b import compare_3b, RUN_CONFIGURATION                 # noqa: E402
from synthetic_instances import (instance_matrix, g5a_fixture_packet,  # noqa: E402
                                 G5A_WRONG_TRANSFORM, make_scalar_artifact,
                                 packet_ii_for_case)
from lauret_evaluator import p_form_spectrum                         # noqa: E402

OUT = os.path.join(ROOT, "rehearsal")
CASE_ID = "SYN-L31-Q4"
N_MAX = 3

RESULTS = []


def check(name, layer, ok, detail=""):
    RESULTS.append({"item": name, "required_layer": layer, "pass": bool(ok),
                    "detail": detail})
    print(f"  {'PASS' if ok else 'FAIL'}  {name}  [{layer}]"
          + (f"  {detail}" if detail and not ok else ""))
    return ok


def scalar_rec(art, n):
    for r in art["records"]:
        if r["sector_id"] == "scalar" and r["harmonic_level"] == n:
            return r
    raise KeyError(n)


def rec_at(art, n, sid):
    for r in art["records"]:
        if r["sector_id"] == sid and r["harmonic_level"] == n:
            return r
    raise KeyError((n, sid))


def main():
    art_a = json.loads(open(os.path.join(
        OUT, f"step3_route_a_{CASE_ID}.json"), "rb").read().decode("ascii"))
    pkt_ii = packet_ii_for_case(CASE_ID, "L(3,1)", N_MAX)
    rows_l31 = p_form_spectrum(RUN_CONFIGURATION["p"], 3, (1, 1), N_MAX,
                               mapping=RUN_CONFIGURATION["mapping"])

    print("Q1: packet_schema validator battery (V1-V8, single-defect mutants)")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "packet_schema.py")],
                       capture_output=True, text=True)
    check("Q1 packet_schema suite exits 0 (all predicates fire, controls clean)",
          "validator", r.returncode == 0, r.stdout[-400:] if r.returncode else "")

    print("\nQ2: conforming instance matrix (pilot-anchored, SYN- unsealable)")
    for label, notes, pkt in instance_matrix():
        clean = not packet_schema.packet_ii_gate(pkt)
        refused = bool(packet_schema.production_seal_refusals(pkt))
        check(f"Q2 {label} accepted clean and seal-refused ({notes})",
              "validator", clean and refused)

    print("\nQ3: G5a rebuilt in the real schema, fault-injected arms")
    fix = g5a_fixture_packet()
    fix_raw = ingestion.canonical_bytes(fix)
    fix_ing = ingestion.ingest(fix_raw, hashlib.sha256(fix_raw).hexdigest())
    fix_obs = make_scalar_artifact("SYN-G5A-FIX", {0: 1, 1: 0, 2: 3, 3: 8})
    g5a = G5a_map_consumption(fix_ing, fix_obs, G5A_WRONG_TRANSFORM)
    check("Q3 G5a passes: clean GREEN, identity RED, wrong RED",
          "comparison", g5a["pass"], json.dumps(g5a["arms"]))
    check("Q3 no arm degraded to a refusal (injection is downstream of validation)",
          "comparison", g5a["no_arm_degraded_to_refusal"])
    check("Q3 fixture sanity: nonidentity, derangement, distinct values",
          "fixture", g5a["fixture_transform_nonidentity"]
          and g5a["fixture_is_derangement"]
          and g5a["reference_values_pairwise_distinct"])

    print("\nQ5: integrated mutation battery")

    # IM1 identity-transform injection (fixture pair; the rehearsal packet's
    # true transform IS identity, so the shift fixture carries this item)
    res = compare_3a(fix_ing.data, fix_obs,
                     transform_override={"kind": "affine", "a": 1, "b": 0})
    check("IM1 identity transform injected at application point reds",
          "comparison RED", res["outcome"] == "RED", res["outcome"])
    # IM1b the stride packet under identity injection (broader family coverage)
    stride_pkt = dict(instance_matrix()[2][2])
    stride_obs = make_scalar_artifact("SYN-Q2-C", {0: 1, 1: 0, 2: 9, 3: 0, 4: 25})
    res_clean = compare_3a(stride_pkt, stride_obs)
    res_inj = compare_3a(stride_pkt, stride_obs,
                         transform_override={"kind": "affine", "a": 1, "b": 0})
    check("IM1b stride packet: clean GREEN, identity injection RED",
          "comparison RED",
          res_clean["outcome"] == "GREEN" and res_inj["outcome"] == "RED",
          f"clean={res_clean['outcome']} injected={res_inj['outcome']}")

    # IM2 preregistered wrong transform injection
    res = compare_3a(fix_ing.data, fix_obs,
                     transform_override=G5A_WRONG_TRANSFORM)
    check("IM2 preregistered wrong transform injected reds",
          "comparison RED", res["outcome"] == "RED", res["outcome"])

    # IM3/IM4 need a branch-asymmetric case: L(7;1,2), route (b), real producer
    inv2 = pow(2, -1, 7)
    al, be = ((1 + 2) * inv2) % 7, ((1 - 2) * inv2) % 7
    mk7 = lambda t: (np.array([np.cos(2 * np.pi * t / 7),                # noqa: E731
                               np.sin(2 * np.pi * t / 7), 0.0, 0.0]))
    gens712 = [(mk7(al), mk7(be))]
    art712 = route_b.produce(gens712, "im-b-L712", "m8_5b-v1-rehearsal",
                             "SYN-L712-IM", adjudication=True)
    n712 = art712["n_max"]
    rows712 = p_form_spectrum(RUN_CONFIGURATION["p"], 7, (1, 2), n712,
                              mapping=RUN_CONFIGURATION["mapping"])
    base712 = compare_3b(rows712, art712, "SYN-L712-IM", n712)
    check("IM3/IM4 precondition: inhomogeneous L(7;1,2) clean 3b is GREEN",
          "comparison", base712["outcome"] == "GREEN",
          json.dumps(base712.get("refusals") or base712["coexact_divergences"]))
    up = {M: rec_at(art712, M - 2, "oneform_coexact_up")["quotient_multiplicity"]
          for M in range(2, n712 + 1)}
    down = {M: rec_at(art712, M, "oneform_coexact_down")["quotient_multiplicity"]
            for M in range(2, n712 + 1)}
    check("IM3/IM4 precondition: branches differ at some in-band M",
          "fixture", any(up[M] != down[M] for M in up), f"up={up} down={down}")
    # N2: assert BOTH contributor choices red, each with a nonempty divergence
    # list, so the mutation cannot be formally present yet non-discriminating on
    # a fixture where the chosen contributor happens to equal the total.
    im3 = {}
    for name, sid in (("up", "oneform_coexact_up"), ("down", "oneform_coexact_down")):
        off = 2 if name == "up" else 0
        r = compare_3b(rows712, art712, "SYN-L712-IM", n712,
                       coexact_total_fn=lambda c, M, s=sid, o=off: c[(M - o, s)])
        im3[name] = (r["outcome"], len(r["coexact_divergences"]))
    check("IM3 EACH single contributor vs evaluator total reds with a nonempty "
          "divergence list (both choices, H1)",
          "comparison RED",
          all(o == "RED" and d > 0 for o, d in im3.values()), str(im3))
    res = compare_3b(rows712, art712, "SYN-L712-IM", n712,
                     coexact_total_fn=lambda c, M: 2 * c[(M - 2, "oneform_coexact_up")])
    check("IM4 doubling shortcut reds on the branch-asymmetric case",
          "comparison RED", res["outcome"] == "RED"
          and bool(res["coexact_divergences"]))

    # IM5 missing consumed cell -> structural refusal
    mut = copy.deepcopy(art_a)
    mut["records"] = [r for r in mut["records"]
                      if not (r["sector_id"] == "scalar"
                              and r["harmonic_level"] == 1)]
    res = compare_3a(pkt_ii, mut)
    check("IM5 missing scalar cell refuses structurally (never a red)",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL"
          and not res["divergences"])

    # IM6 unknown packet field -> validator layer
    mut_pkt = copy.deepcopy(pkt_ii)
    mut_pkt["adjudicator_notes"] = "x"
    res = compare_3a(mut_pkt, art_a)
    check("IM6 unknown packet field refuses at the packet-validation layer",
          "validator reject", res["outcome"] == "STRUCTURAL_REFUSAL"
          and res.get("layer", "").startswith("packet validation"))

    # IM7 fixture-private shapes fail loudly, never coerce
    old_shape = {"case_id": CASE_ID,
                 "citation": dict(pkt_ii["citation"]),
                 "indexing_map": {"source_k": "protocol_n"},
                 "reference_values": [0, 0, 0],
                 "format_version": "m8_5b-packet-1"}
    res = compare_3a(old_shape, art_a)
    check("IM7 packet in the burned-era fixture-private shape refuses loudly",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL")
    bad_rows = [{"k": 1, "lower": {"eigenvalue": 3}}]
    res = compare_3b(bad_rows, art_a, CASE_ID, N_MAX)
    check("IM7b evaluator rows in a foreign shape refuse loudly",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL")

    # IM12: a canonical NON-OBJECT Packet II travels the REAL ingestion path
    # (hash-verified, canonical, parsed) and must be REFUSED at the
    # packet-validation layer without an exception (Redline audit blocker 1)
    im12_ok, im12_detail = True, []
    for bad_top in ([], None, 3, "x"):
        raw = ingestion.canonical_bytes(bad_top)
        ing_bad = ingestion.ingest(raw, hashlib.sha256(raw).hexdigest())
        try:
            res = compare_3a(ing_bad.data, art_a)
            ok_one = (res["outcome"] == "STRUCTURAL_REFUSAL"
                      and res.get("layer", "").startswith("packet validation"))
        except Exception as exc:
            ok_one = False
            im12_detail.append(f"{bad_top!r} CRASHED: {exc!r}")
        im12_ok &= ok_one
        if not ok_one and not im12_detail:
            im12_detail.append(f"{bad_top!r} -> {res['outcome']}")
    check("IM12 canonical non-object Packet II refuses at packet validation, "
          "never crashes, on the real ingestion path",
          "STRUCTURAL_REFUSAL", im12_ok, "; ".join(im12_detail))

    # IM8 G9 transcription mutation, schema-valid cell choice
    raw_ii = ingestion.canonical_bytes(pkt_ii)
    ing_ii = ingestion.ingest(raw_ii, hashlib.sha256(raw_ii).hexdigest())

    def perturb(d):
        d["reference_values"][2][1] = 4          # level-2 value 3 -> 4; n_max fixed

    mutant_valid = copy.deepcopy(pkt_ii)
    perturb(mutant_valid)
    check("IM8 precondition: the perturbed packet is still schema-valid",
          "fixture", not packet_schema.packet_ii_checks(mutant_valid))
    g9 = G9_transcription_mutation(
        ing_ii, lambda d, o: compare_3a(d, o)["outcome"] == "GREEN", art_a,
        perturb, description="reference level-2 multiplicity 3 -> 4")
    # G9's boolean collapses RED and STRUCTURAL_REFUSAL, so the required layer
    # is asserted DIRECTLY on the mutated comparison result (Redline audit
    # blocker 2): outcome RED, zero refusals, and exactly the injected
    # divergence.  Credit needs both this and G9's own pass.
    res_mut = compare_3a(mutant_valid, art_a)
    layer_ok = (res_mut["outcome"] == "RED"
                and not res_mut["refusals"]
                and any(d["level"] == 2 and d["reference"] == 4
                        and d["observed"] == 3 for d in res_mut["divergences"]))
    check("IM8 G9: baseline green, one perturbed cell reds the comparison, "
          "and the mutated outcome IS a comparison RED with no refusal",
          "comparison RED", g9["pass"] and layer_ok,
          f"g9={g9.get('reason', 'pass')} mutated_outcome={res_mut['outcome']} "
          f"refusals={len(res_mut['refusals'])}")

    # IM9 extra/duplicate required-surface coordinates -> structural refusal
    mut = copy.deepcopy(art_a)
    mut["records"].append(copy.deepcopy(scalar_rec(art_a, 1)))
    res = compare_3a(pkt_ii, mut)
    check("IM9 duplicate scalar record refuses structurally",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL")
    res = compare_3b(rows_l31 + [copy.deepcopy(rows_l31[0])], art_a,
                     CASE_ID, N_MAX)
    check("IM9b duplicate evaluator row refuses structurally",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL")

    # IM10 a structurally red baseline earns NO mutation credit
    mut = copy.deepcopy(art_a)
    mut["records"] = [r for r in mut["records"]
                      if not (r["sector_id"] == "scalar"
                              and r["harmonic_level"] == 1)]
    g9 = G9_transcription_mutation(
        ing_ii, lambda d, o: compare_3a(d, o)["outcome"] == "GREEN", mut,
        perturb, description="mutation over a structurally red baseline")
    check("IM10 harness refuses credit when the clean baseline is already red",
          "refusal, no credit", (not g9["pass"])
          and g9.get("reason", "").startswith("refused"), g9.get("reason", ""))

    # IM11 a required cell NEITHER adapter consumes: validator layer catches it,
    # and both adapters stay green on the same mutant (the delegation is real)
    mut = copy.deepcopy(art_a)
    mut["records"] = [r for r in mut["records"]
                      if not (r["sector_id"] == "oneform_exact"
                              and r["harmonic_level"] == 0)]
    v = schema.validate_records(mut["records"], schema.expected_cells(N_MAX))
    a3 = compare_3a(pkt_ii, mut)["outcome"]
    b3 = compare_3b(rows_l31, mut, CASE_ID, N_MAX)["outcome"]
    check("IM11 removing the level-0 exact null cell reds the artifact validator",
          "artifact validator", not v["pass"], json.dumps(v["problems"]))
    check("IM11b both adapters stay GREEN on that mutant (delegation load-bearing)",
          "adapters unaffected", a3 == "GREEN" and b3 == "GREEN",
          f"3a={a3} 3b={b3}")

    print("\nband authority (12.1.4): the three verdict rows")
    mut = copy.deepcopy(art_a)
    mut["n_max"] = 5
    res = compare_3a(pkt_ii, mut)
    check("BA1 claimed n_max contradicting own cells refuses structurally",
          "STRUCTURAL_REFUSAL", res["outcome"] == "STRUCTURAL_REFUSAL"
          and any("band authority" in x for x in res["refusals"]))
    # BA2 (coverage refusal) is the same defect IM5 injects; its verdict is
    # DERIVED from IM5's recorded result so this line can fail with it.
    im5 = next(r for r in RESULTS if r["item"].startswith("IM5"))
    check("BA2 coverage refusal exercised (derived from IM5's result)",
          "STRUCTURAL_REFUSAL", im5["pass"])
    mut = copy.deepcopy(art_a)
    scalar_rec(mut, 1)["quotient_multiplicity"] = 7.0
    scalar_rec(mut, 1)["eigenvalue_R2"] = 3.0
    mut["n_max"] = 2                      # self-consistent with its own cells now
    res = compare_3a(pkt_ii, mut)
    check("BA3 self-consistent route with divergent band: RED with the mismatch "
          "reported, never a refusal", "comparison RED",
          res["outcome"] == "RED"
          and res["band_authority"]["band_mismatch_vs_packet"]
          and any(d["level"] == 1 for d in res["divergences"]),
          json.dumps(res.get("refusals") or res.get("divergences")))
    mut = copy.deepcopy(art_a)
    scalar_rec(mut, 2)["quotient_multiplicity"] = 4.0
    res = compare_3a(pkt_ii, mut)
    check("BA4 pure in-band value divergence: RED, band consistent (the Q2 "
          "corrupted-artifact pairing)", "comparison RED",
          res["outcome"] == "RED"
          and not res["band_authority"]["band_mismatch_vs_packet"]
          and res["divergences"][0]["reference_provenance"] == "entry")

    failed = [r for r in RESULTS if not r["pass"]]
    blob = (json.dumps({"results": RESULTS, "failed": len(failed)},
                       sort_keys=True, indent=2, ensure_ascii=True) + "\n")
    open(os.path.join(OUT, "QUALIFY_RECORD.json"), "w").write(blob)
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} items PASS; "
          f"record written to rehearsal/QUALIFY_RECORD.json")
    if failed:
        print("FAILED items:")
        for f in failed:
            print(f"  {f['item']}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
