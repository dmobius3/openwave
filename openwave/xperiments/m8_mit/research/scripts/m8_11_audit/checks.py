"""PASS/FAIL checks with a recorded mutation for every PASS line."""
import json, os

LOG = []


def check(label, pred, real, mutated, mutation):
    """pred(x) -> bool. Print PASS/FAIL on real; require pred(mutated) False and record it."""
    ok = bool(pred(real))
    mut_ok = bool(pred(mutated))
    status = "PASS" if ok else "FAIL"
    mstat = "mutation detected (check FAILS on mutant)" if not mut_ok else "MUTATION NOT DETECTED"
    print(f"{status}: {label}   [mutation: {mutation} -> {mstat}]")
    LOG.append({"label": label, "status": status, "mutation": mutation, "mutation_detected": not mut_ok})
    if mut_ok:
        raise RuntimeError("mutation not detected for check: " + label)
    return ok


def save_log(name):
    os.makedirs("out", exist_ok=True)
    with open(os.path.join("out", f"checks_{name}.json"), "w") as f:
        json.dump(LOG, f, indent=1)
