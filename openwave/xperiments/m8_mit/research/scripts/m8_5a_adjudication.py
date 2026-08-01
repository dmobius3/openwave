"""M8.5-A ADJUDICATION: the protocol section 7 comparison, run after unsealing.

WHAT THIS IS.  The maintainer-side comparison the
[M8.5-A protocol](../findings/m8_5a_reproduction_protocol.md) section 7 specifies, written
AFTER the section 9 commitment landed on main (commit dac2b6a1, merged c3dc2b5f) and
committed separately from it, so the ordering is auditable: the clean-room result was
committed before this file, or the target transcription below, existed anywhere.

THE TARGET, PINNED.  The scalar first-occurrence table of the locked pre-registration's
section 6.1, exactly as at commit ec877ee0 (content unchanged since 269456b7).  TARGET_61
below is THIS ADJUDICATOR'S OWN transcription from that pinned text, keyed by the
pre-registration's own irrep labels; it is deliberately not copied from object 2's DOC
fixture, so the two transcriptions are independent and a shared typo is not silently
self-consistent.  The section 6.1 table carries labels and McKay distances but no
dimensions, and rows are matched label-free by (dim, distance), so the label-to-dimension
map is taken from object 1's literals (LABEL_DIMS below), which is exactly the adjudication
use the protocol opens objects 1 and 2 to after commitment.  The map is checked for internal
consistency against the transcribed distances before any cell is compared.

THE COMPARISON.  Exact integer equality, cell by cell, all three columns, no tolerances.
Rows matched label-free by the (dim, mckay_distance) signature; signatures verified pairwise
distinct on each side first; row counts verified equal.  Scalar result categories, frozen in
section 7: reproduced / partial disagreement / structural failure / not completed.

THREE-WAY AGREEMENT.  Object 2's committed reconstruction JSON is read as a second reference,
and agreement of the clean-room result with both section 6.1 and object 2 is reported as
three-way agreement per section 7.  It does not raise the claim label.

THE TRANSCRIPTION MUTATION (G10).  --mutation-tests perturbs one transcribed target cell
(the doc_typo pattern from object 2) and requires the comparison to go red, so the
comparison itself is a check that can fail.  A second mutation perturbs a copy of the
result instead, so the harness is shown sensitive on both sides.

THE SECTION 6 MODULE IS NOT ADJUDICATED HERE.  Its verdict was fixed in the section 9
commitment; the standing rule says numerical agreement upgrades nothing, and the derivation
behind the implementer's `structurally derived and reproduced` verdict is examined in the
adversarial audit, not by this harness.  This file only echoes the declared verdict into the
adjudication record.

EXIT CODES.  0 = adjudication completed, every cell reproduced.  2 = adjudication completed
with disagreement (a finding, not a failure of the task).  1 = structural failure, an
incomplete run, or a mutation that fails to redden its target.

Run: python3 m8_5a_adjudication.py [--mutation-tests] [--json PATH]
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE.parent / "data" / "m8_5a_result.json"
OBJECT2_PATH = HERE.parent / "data" / "m8_2_indep_reconstruction.json"
DEFAULT_OUT = HERE.parent / "data" / "m8_5a_adjudication.json"

COMMITMENT_COMMIT = "dac2b6a1"  # the section 9 commitment this adjudication postdates
TARGET_PIN = "ec877ee0"         # the pre-registration commit section 7 pins

# ----------------------------------------------------------------------------------------
# The adjudicator's own transcription of pre-registration section 6.1 at ec877ee0.
# Column order as printed there: label, d, trivial, standard, Galois.
# ----------------------------------------------------------------------------------------
TARGET_61 = {
    #  label   d   trivial standard galois
    "R_1": (1, 1, 1, 5),
    "R_3": (2, 2, 0, 4),
    "R_6": (3, 3, 1, 3),
    "R_7": (4, 4, 2, 2),
    "R_8": (5, 5, 3, 1),
    "R_4": (6, 6, 4, 0),
    "R_5": (6, 6, 4, 2),
    "R_2": (7, 7, 5, 3),
    "R_0": (0, 0, 2, 6),
}

# Label -> dimension, from object 1's literals (labels/dims/dist arrays), used here as the
# adjudication reference the protocol permits after commitment.  Object 1 spells the labels
# without underscores; the pre-registration writes R_n.  Object 1's dist array is carried
# too, so the transcription and the reference can be cross-checked before use.
LABEL_DIMS = {
    "R_0": (1, 0),
    "R_1": (2, 1),
    "R_2": (2, 7),
    "R_3": (3, 2),
    "R_4": (3, 6),
    "R_5": (4, 6),
    "R_6": (4, 3),
    "R_7": (5, 4),
    "R_8": (6, 5),
}

COLUMNS = ("trivial", "standard", "galois")


def build_target(target_61: dict) -> dict:
    """Label-free target: (dim, distance) -> (trivial, standard, galois).

    Structural failure (raised) if a label is unknown, if object 1's distance for a label
    disagrees with the transcribed d, or if the resulting signatures collide."""
    out = {}
    for label, (d, triv, std, gal) in target_61.items():
        if label not in LABEL_DIMS:
            raise SystemExit(f"STRUCTURAL FAILURE: unknown label {label!r} in transcription")
        dim, d_ref = LABEL_DIMS[label]
        if d != d_ref:
            raise SystemExit(
                f"STRUCTURAL FAILURE: transcribed d={d} for {label} disagrees with "
                f"object 1's distance {d_ref}"
            )
        sig = (dim, d)
        if sig in out:
            raise SystemExit(f"STRUCTURAL FAILURE: duplicate target signature {sig}")
        out[sig] = (triv, std, gal)
    return out


def load_result_rows(result: dict) -> dict:
    rows = {}
    for r in result["rows"]:
        sig = (int(r["dim"]), int(r["mckay_distance"]))
        if sig in rows:
            raise SystemExit(f"STRUCTURAL FAILURE: duplicate result signature {sig}")
        vals = r["n_first"]
        rows[sig] = tuple(int(vals[c]) for c in COLUMNS)
    return rows


def load_object2_rows(path: Path) -> dict:
    """Object 2's committed table, keyed 'dim<D>_d<K>' -> [trivial, standard, galois]."""
    table = json.loads(path.read_text())["table"]
    out = {}
    for key, vals in table.items():
        dim_s, d_s = key.split("_")
        out[(int(dim_s[3:]), int(d_s[1:]))] = tuple(int(v) for v in vals)
    return out


def compare(target: dict, candidate: dict) -> tuple:
    """Returns (category, mismatches).  Exact equality only."""
    if set(target) != set(candidate):
        return "structural failure", {
            "target_only": sorted(set(target) - set(candidate)),
            "candidate_only": sorted(set(candidate) - set(target)),
        }
    mismatches = []
    for sig in sorted(target):
        for i, col in enumerate(COLUMNS):
            if target[sig][i] != candidate[sig][i]:
                mismatches.append({
                    "signature": list(sig), "column": col,
                    "target": target[sig][i], "candidate": candidate[sig][i],
                })
    return ("reproduced" if not mismatches else "partial disagreement"), mismatches


def structural_gate_check(result: dict) -> list:
    """Section 7 structural-failure conditions visible from the committed result."""
    problems = []
    for g in result.get("gates", []):
        if not g.get("pass"):
            problems.append(f"gate {g.get('gate')} recorded as failing")
    if len(result.get("rows", [])) != 9:
        problems.append(f"row count {len(result.get('rows', []))} != 9")
    for r in result.get("rows", []):
        for c in COLUMNS:
            v = r["n_first"][c]
            if not isinstance(v, int) or v < 0:
                problems.append(f"non-integer or negative cell {r['dim']},{r['mckay_distance']},{c}")
    return problems


def run_mutation_suite(target: dict, result_rows: dict) -> list:
    """Both sides of the comparison must be able to go red."""
    results = []

    # doc_typo: one transcribed target cell altered (the object 2 pattern).
    mut_target = dict(target)
    sig = sorted(mut_target)[0]
    vals = list(mut_target[sig])
    vals[0] += 1
    mut_target[sig] = tuple(vals)
    cat, mm = compare(mut_target, result_rows)
    results.append({
        "mutation": "doc_typo_target_cell",
        "expected": "comparison goes red",
        "observed_category": cat,
        "detected": cat != "reproduced" and len(mm) >= 1,
    })

    # result_typo: one candidate cell altered, same requirement from the other side.
    mut_rows = copy.deepcopy(result_rows)
    sig2 = sorted(mut_rows)[-1]
    vals2 = list(mut_rows[sig2])
    vals2[2] += 1
    mut_rows[sig2] = tuple(vals2)
    cat2, mm2 = compare(target, mut_rows)
    results.append({
        "mutation": "result_typo_candidate_cell",
        "expected": "comparison goes red",
        "observed_category": cat2,
        "detected": cat2 != "reproduced" and len(mm2) >= 1,
    })
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--mutation-tests", action="store_true")
    args = ap.parse_args()

    result = json.loads(RESULT_PATH.read_text())
    if result.get("schema_version") != "m8_5a-v1":
        print(f"STRUCTURAL FAILURE: schema {result.get('schema_version')!r} != 'm8_5a-v1'")
        return 1

    problems = structural_gate_check(result)
    target = build_target(TARGET_61)
    result_rows = load_result_rows(result)
    object2_rows = load_object2_rows(OBJECT2_PATH)

    category, mismatches = compare(target, result_rows)
    if problems:
        category = "structural failure"

    cat_o2_vs_result, mm_o2_result = compare(object2_rows, result_rows)
    cat_o2_vs_target, mm_o2_target = compare(object2_rows, target)
    three_way = (category == "reproduced" and cat_o2_vs_result == "reproduced"
                 and cat_o2_vs_target == "reproduced")

    mut = run_mutation_suite(target, result_rows) if args.mutation_tests else None
    mut_ok = all(m["detected"] for m in mut) if mut else None

    record = {
        "what": "M8.5-A section 7 adjudication",
        "target": f"pre-registration section 6.1 at {TARGET_PIN}, adjudicator's own transcription",
        "commitment_commit": COMMITMENT_COMMIT,
        "result_file": RESULT_PATH.name,
        "packet_sha256": result.get("packet_sha256"),
        "structural_problems": problems,
        "category": category,
        "mismatches": mismatches,
        "three_way": {
            "object2_vs_result": cat_o2_vs_result,
            "object2_vs_target": cat_o2_vs_target,
            "agreement": three_way,
        },
        "coexact_module": {
            "declared_verdict": result.get("coexact_module", {}).get("verdict"),
            "note": "echoed from the committed result; standing examined by the adversarial "
                    "audit, never upgraded by numerical agreement",
        },
        "mutation_suite": mut,
        "mutation_suite_all_detected": mut_ok,
    }
    args.json.write_text(json.dumps(record, indent=2) + "\n")

    print(f"target       section 6.1 at {TARGET_PIN} (own transcription, 9 rows)")
    print(f"candidate    {RESULT_PATH.name} (committed at {COMMITMENT_COMMIT})")
    print(f"category     {category.upper()}")
    if problems:
        for p in problems:
            print(f"  problem: {p}")
    if mismatches:
        for m in mismatches:
            print(f"  mismatch at {tuple(m['signature'])} {m['column']}: "
                  f"target {m['target']} vs result {m['candidate']}")
    print(f"three-way    object2 vs result: {cat_o2_vs_result} | object2 vs target: "
          f"{cat_o2_vs_target} | agreement: {three_way}")
    print(f"coexact      declared verdict: {record['coexact_module']['declared_verdict']} "
          f"(not adjudicated here)")
    if mut is not None:
        for m in mut:
            print(f"  {'DETECTED' if m['detected'] else 'MISSED  '}  {m['mutation']}: "
                  f"category under mutation = {m['observed_category']}")
    print(f"written      {args.json}")

    if category == "structural failure" or (mut_ok is False):
        return 1
    if category == "partial disagreement":
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
