"""The § 4.1 Step-3 runner: commit two independent route artifacts, and stop.

DELIBERATELY STUPID.  It orchestrates the frozen producers, checks that each
artifact is structurally sound ON ITS OWN, writes both canonically, and stops.

WHAT IT MUST NOT DO, stated so a later refactor cannot drift:

    It does NOT compare route (a) with route (b), inspect multiplicity
    agreement, reconcile n_max, apply G2-G8, compute RESOLVED or UNRESOLVED,
    determine an achieved band, invoke the 3b evaluator, ingest Packet II, or
    adjudicate anything.

Even the checks it DOES perform are strictly route-local.  It verifies that each
artifact's declared effective-group digest matches the action THAT route was
given.  It contains no rule saying A's digest must equal B's, or A's n_max must
equal B's, although qualification showed both agree.  Those are cross-route facts
and belong to step 5.  Putting them here would migrate adjudication backwards
into production.

THE CENTRAL CRITERION.  A valid Step-3 runner must be capable of faithfully
committing two MUTUALLY CONTRADICTORY route outputs.  If it cannot, adjudication
has leaked upstream.  `test_commits_contradictory_outputs` in the qualification
harness asserts exactly that.

QUALIFICATION ARTIFACTS ARE NOT STEP-3 OUTPUTS.  The tuning artifacts produced
while qualifying the producers prove the producers are Step-3-ready.  They are
not Step-3 evidence and are never copied, promoted or renamed into it.  This
runner generates fresh outputs when § 4.1 actually reaches Step 3.
"""

import sys

sys.path.insert(0, "../gates")

import route_gates as gates                                 # noqa: E402
import step3_schema as schema                               # noqa: E402

FORBIDDEN = (
    "compare routes, inspect multiplicity agreement, reconcile n_max, apply "
    "G2-G8, compute RESOLVED/UNRESOLVED, determine an achieved band, invoke the "
    "3b evaluator, ingest Packet II, adjudicate"
)


class Step3Error(Exception):
    """Raised on a structural failure.  Never swallowed, never downgraded."""


def _validate_one(artifact, action_pairs, label):
    """Route-LOCAL structural validation.  Nothing here looks at the other route."""
    problems = []

    if artifact.get("schema_version") != schema.SCHEMA_VERSION:
        problems.append(f"{label}: schema_version {artifact.get('schema_version')!r} "
                        f"!= {schema.SCHEMA_VERSION!r}")

    recs = artifact.get("records", [])
    v = schema.validate_records(recs)
    problems.extend(f"{label}: {p}" for p in v["problems"])

    # census: complete, unique, and derived from the artifact's OWN n_max
    n_max = artifact.get("n_max")
    keys = [(r.get("harmonic_level"), r.get("sector_id")) for r in recs]
    if len(keys) != len(set(keys)):
        dup = sorted({k for k in keys if keys.count(k) > 1})
        problems.append(f"{label}: duplicate coordinates {dup}")
    expected = {(n, s) for n in range(0, (n_max or 0) + 1)
                for s in ("scalar", "oneform_exact",
                          "oneform_coexact_up", "oneform_coexact_down")}
    missing = sorted(expected - set(keys))
    extra = sorted(set(keys) - expected)
    if missing:
        problems.append(f"{label}: FULL-BAND INCOMPLETE, missing {missing}")
    if extra:
        problems.append(f"{label}: coordinates outside the required band {extra}")

    # digest: this route's declaration against the action THIS route was given.
    # Not against the other route's digest.
    _, _, digest = gates.certified_effective_group(action_pairs)
    if artifact.get("effective_group_digest") != digest:
        problems.append(f"{label}: declared effective-group digest does not match "
                        f"the action supplied to this route")

    # canonical round-trip
    import json
    raw = schema.canonical_bytes(artifact)
    back = json.loads(raw.decode("ascii"))
    if schema.canonical_bytes(back) != raw:
        problems.append(f"{label}: artifact does not canonical-round-trip")

    return problems


def run(route_a_artifact, route_b_artifact, action_pairs, writer=None):
    """Validate each artifact route-locally, then commit both.  Then stop.

    `writer(label, bytes) -> None` performs the actual commit; supplying it
    keeps the runner testable without touching a repository.
    """
    problems = []
    problems += _validate_one(route_a_artifact, action_pairs, "route a")
    problems += _validate_one(route_b_artifact, action_pairs, "route b")
    if problems:
        raise Step3Error("; ".join(problems))

    committed = {}
    for label, art in (("a", route_a_artifact), ("b", route_b_artifact)):
        raw = schema.canonical_bytes(art)
        committed[label] = {"sha256": schema.artifact_sha256(art),
                            "byte_length": len(raw)}
        if writer:
            writer(label, raw)

    return {
        "step": 3,
        "committed": committed,
        "cross_route_checks_performed": "none, by design",
        "forbidden_here": FORBIDDEN,
        "note": ("both artifacts are now frozen evidence; step 4 opens Packet II "
                 "and step 5 owns every cross-route question"),
    }
