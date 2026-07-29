import json
J=json.load(open("../solverA/m8_1_1_coexact.json"))
print("TOP KEYS:", list(J.keys()))
print("CONVENTIONS:", json.dumps(J["conventions"], indent=1))
print("AUDIT:", json.dumps(J["audit"], indent=1))
print("T9b:", json.dumps(J["T9b_full_unreduced_sweep"], indent=1))
for r in J["T9"]:
    print("T9", r["group"], "m=",r["m"], "tau=",r["tau_sigma"], "dimtau=",r["dim_tau"],
          "rank=",r["svd_rank"], "csum=",r["character_sum"], "agree=",r["agree"],
          "idem_err=%.2e"%r["projector_idempotency_err"])
print()
print("GROUPS:", list(J["groups"].keys()))
for n,b in J["groups"].items():
    print("="*70)
    print(n, "order", b["order_verified"], "expected", b["order_expected"], "match", b["order_matches"],
          "exp", b["exponent"], "-I", b["minus_I_in_group"], "fallback", b["fallback_generator_list_used"],
          "explicit_agrees", b["closure_equals_explicit_spec_list"])
    print(" dims:", [c["dim"] for c in b["irreducible_characters"]])
    print(" dist:", b["distance_vector"], "diam", b["graph_diameter"])
    print(" checks:", json.dumps(b["checks"]))
    print(" A:", b["adjacency_A"])
    print(" T1 least_a:", [(t["sigma"],t["dim"],t["d"],t["least_a"]) for t in b["T1"]])
    print(" T2:", json.dumps(b["T2"]))
    print(" T4:", json.dumps(b["T4"]))
    print(" T5:", json.dumps(b["T5"]))
    print(" T6:", json.dumps(b["T6"]))
