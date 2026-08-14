"""G9 plus the rung-3a comparison layer, implemented against addendum 12.1.

G9 (unchanged): the transcription mutation required of any comparison harness.
Protocol section 8: "any comparison harness against a transcribed reference
carries a transcription mutation (perturb one transcribed cell; the comparison
must go red)".  It is a meta-gate on a harness, not a check on data, and it
REFUSES to award credit when ingestion did not hash-verify or when the clean
baseline is already red (that refusal is qualification item IM10).

The comparison layer below replaces the burned execution's fixture-private
shapes.  It consumes ONLY the addendum 12.1.2 Packet II schema and the frozen
section 5 record surface, and it separates three outcomes that are never
merged:

    STRUCTURAL_REFUSAL   the inputs are malformed, incomplete, self-
                         contradictory, or fail packet validation (12.1.3
                         re-run at ingestion): no comparison verdict exists
    RED                  the comparison ran and found divergences: this is
                         adjudication content, a result
    GREEN                the comparison ran and found none

Fault injection (Q3): `transform_override` is applied at the transform-
application point, downstream of packet validation, on in-memory comparison
state that is deliberately not re-validated.  It exists for the G5a arms and
the IM1/IM2 battery items; production adjudication never passes it.

Path note: this module needs `packet_schema` (working-folder root) and is
imported next to `step3_schema` consumers; it bootstraps both paths so it
behaves identically from any caller.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, os.path.join(_ROOT, "production")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import packet_schema                                        # noqa: E402

__all__ = ["G9_transcription_mutation", "compare_3a", "G5a_map_consumption",
           "project_scalar_observed", "band_authority", "recompute_ceiling",
           "integer_near", "STEP3_SCHEMA_VERSION"]

STEP3_SCHEMA_VERSION = "m8_5b-v1"   # section 6.4, frozen
INTEGER_NEARNESS_TOL = 1e-6          # section 6.3, recomputed here, never trusted


# --- G9, unchanged from the qualified implementation --------------------------

def _first_divergence(clean, mutated):
    """A named cell where the two references differ, for the report."""
    if isinstance(clean, dict) and isinstance(mutated, dict):
        for k in sorted(set(clean) | set(mutated)):
            sub = _first_divergence(clean.get(k), mutated.get(k))
            if sub is not None:
                return (str(k),) + sub if isinstance(sub, tuple) else (str(k),)
        return None
    if isinstance(clean, list) and isinstance(mutated, list):
        for i, (a, b) in enumerate(zip(clean, mutated)):
            sub = _first_divergence(a, b)
            if sub is not None:
                return (f"[{i}]",) + sub if isinstance(sub, tuple) else (f"[{i}]",)
        return None
    return () if clean != mutated else None


def G9_transcription_mutation(ingested, compare_fn, observed,
                              perturb, description="one transcribed cell perturbed",
                              collateral_fns=None):
    """G9 on a single comparison harness.

    `compare_fn(reference_data, observed) -> bool` returns True for agreement.
    `perturb(data) -> None` mutates ONE cell of the parsed reference in place.
    `collateral_fns` is an optional mapping name -> fn(reference_data) -> bool,
    reporting honestly which other predicates the same perturbation disturbs.
    """
    if not getattr(ingested, "hash_verified", False):
        return {"gate": "G9", "pass": False,
                "reason": "refused: ingestion did not hash-verify, so a red "
                          "below would not be attributable to the comparison layer"}

    baseline_green = bool(compare_fn(ingested.data, observed))
    if not baseline_green:
        return {"gate": "G9", "pass": False,
                "reason": "refused: the clean baseline comparison is already red, "
                          "so the mutation proves nothing",
                "provenance": ingested.provenance()}

    mutant = ingested.mutated_copy(description, perturb)
    mutated_green = bool(compare_fn(mutant.data, observed))
    cell = _first_divergence(ingested.data, mutant.data)

    collateral = {}
    for name, fn in (collateral_fns or {}).items():
        try:
            collateral[name] = {"clean": bool(fn(ingested.data)),
                                "mutated": bool(fn(mutant.data))}
        except Exception as exc:                       # a crash is not a verdict
            collateral[name] = {"error": repr(exc)}

    return {
        "gate": "G9",
        "baseline_comparison_green": baseline_green,
        "mutated_comparison_green": mutated_green,
        "perturbed_cell": ".".join(cell) if cell else None,
        "mutation": description,
        "pass": bool(baseline_green and not mutated_green),
        "provenance": mutant.provenance(),
        "collateral": collateral,
        "credited": ("the comparison predicate alone; the hash gate is not "
                     "credited here and stayed green throughout"),
    }


# --- the rung-3a comparison layer (addendum 12.1.4) ---------------------------

def integer_near(x):
    """Section 6.3 integer-nearness, recomputed from the consumed value itself.

    Returns the rounded nonnegative integer, or None when the value is null,
    non-numeric, non-integral within 1e-6, or negative.  The comparator never
    trusts the artifact's own `integer_nearness_margin` field.
    """
    if isinstance(x, bool) or not isinstance(x, (int, float)):
        return None
    r = round(x)
    if abs(x - r) > INTEGER_NEARNESS_TOL or r < 0:
        return None
    return int(r)


def _scalar_metadata_refusals(rec, packet_case_id, idx):
    bad = []
    if rec.get("arena_case_id") != packet_case_id:
        bad.append(f"scalar record {idx}: arena_case_id "
                   f"{rec.get('arena_case_id')!r} != packet case_id "
                   f"{packet_case_id!r}")
    if rec.get("schema_version") != STEP3_SCHEMA_VERSION:
        bad.append(f"scalar record {idx}: schema_version "
                   f"{rec.get('schema_version')!r} != {STEP3_SCHEMA_VERSION!r}")
    if rec.get("form_degree") != 0:
        bad.append(f"scalar record {idx}: form_degree {rec.get('form_degree')!r} != 0")
    if rec.get("hodge_sector") is not None:
        bad.append(f"scalar record {idx}: hodge_sector "
                   f"{rec.get('hodge_sector')!r} != null")
    if rec.get("rung") != "3a":
        bad.append(f"scalar record {idx}: rung {rec.get('rung')!r} != '3a' "
                   "(adjudication-case scalar records carry rung 3a)")
    return bad


def project_scalar_observed(artifact, packet_case_id, n_max):
    """The 12.1.4 observed projection.  Returns a dict with `refusals`,
    `values` (level -> compared integer, only [0, n_max]), `all_levels`
    (every scalar level -> integer, for band authority), and
    `ignored_above_band` (count, reported, never a refusal by itself)."""
    refusals, values, all_levels = [], {}, {}
    ignored = 0
    for i, rec in enumerate(artifact.get("records", [])):
        if rec.get("sector_id") != "scalar":
            continue
        refusals += _scalar_metadata_refusals(rec, packet_case_id, i)
        n = rec.get("harmonic_level")
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            refusals.append(f"scalar record {i}: harmonic_level {n!r} malformed")
            continue
        v = integer_near(rec.get("quotient_multiplicity"))
        if v is None:
            refusals.append(
                f"scalar record {i} (level {n}): quotient_multiplicity "
                f"{rec.get('quotient_multiplicity')!r} is null, non-numeric, "
                "negative, or non-integral within 1e-6: the compared integer "
                "does not exist")
            continue
        if n in all_levels:
            refusals.append(f"duplicate scalar record at level {n}: the cell "
                            "value is ambiguous")
            continue
        all_levels[n] = v
        if n <= n_max:
            values[n] = v
        else:
            ignored += 1
    missing = [n for n in range(0, n_max + 1) if n not in values]
    if missing:
        refusals.append(f"scalar cells missing for levels {missing}: required "
                        "comparison support absent")
    return {"refusals": refusals, "values": values, "all_levels": all_levels,
            "ignored_above_band": ignored}


def recompute_ceiling(levels_to_int):
    """Section 6.2 applied to a route's own scalar cells: the second smallest
    n >= 1 with multiplicity > 0, or None when fewer than two exist."""
    positives = sorted(n for n, m in levels_to_int.items() if n >= 1 and m > 0)
    return positives[1] if len(positives) >= 2 else None


def band_authority(artifact, all_levels, packet_n_max):
    """The 12.1.4 band-authority verdicts.  Returns (refusals, report).

    Row 1: claimed n_max vs the ceiling recomputed from the artifact's own
    cells (or that ceiling not recomputable): STRUCTURAL REFUSAL.
    Row 3: recomputed ceiling vs the packet's n_max: NOT a refusal; it is
    reported as the consequence of an in-band divergence, which the per-level
    comparison owns."""
    refusals = []
    claimed = artifact.get("n_max")
    own = recompute_ceiling(all_levels)
    if own is None:
        refusals.append("band authority: fewer than two positive nonzero scalar "
                        "levels exist among the artifact's own cells, so its "
                        "section 6.2 ceiling cannot be recomputed: the artifact "
                        "contradicts its own band claim")
    elif claimed != own:
        refusals.append(f"band authority: the artifact claims n_max {claimed!r} "
                        f"but its own cells give {own}: the artifact contradicts "
                        "its own data")
    report = {
        "packet_n_max": packet_n_max,
        "artifact_claimed_n_max": claimed,
        "recomputed_from_artifact_cells": own,
        "band_mismatch_vs_packet": (own is not None and own != packet_n_max),
        "note": ("the packet's V7-derived n_max owns the comparison band; a "
                 "self-consistent route ceiling differing from it can arise "
                 "only through an in-band divergence, which the per-level "
                 "comparison reports as the result"),
    }
    return refusals, report


def compare_3a(packet_data, artifact, transform_override=None):
    """The rung-3a comparison, per addendum 12.1.4.

    `packet_data` is the PARSED Packet II handed over by staged ingestion (the
    hash gate lives there and is not re-credited here).  The packet is
    re-validated (12.1.3 enforcement at ingestion) BEFORE any injection;
    `transform_override` then applies at the transform-application point only.

    Returns a dict whose `outcome` is STRUCTURAL_REFUSAL, RED, or GREEN, with
    per-level divergences carrying the 12.1.2b provenance class of every
    zero-valued reference cell.
    """
    violations = packet_schema.packet_ii_gate(packet_data)
    if violations:
        return {"outcome": "STRUCTURAL_REFUSAL",
                "layer": "packet validation (12.1.3 re-run at ingestion)",
                "refusals": violations, "divergences": []}

    ref, provenance, n_max = packet_schema.filled_reference(
        packet_data, transform_override)
    proj = project_scalar_observed(artifact, packet_data["case_id"], n_max)
    band_refusals, band_report = band_authority(
        artifact, proj["all_levels"], n_max)

    refusals = proj["refusals"] + band_refusals
    if refusals:
        return {"outcome": "STRUCTURAL_REFUSAL", "layer": "observed projection",
                "refusals": refusals, "divergences": [],
                "band_authority": band_report}

    divergences = []
    for n in range(0, n_max + 1):
        if ref[n] != proj["values"][n]:
            divergences.append({"level": n, "reference": ref[n],
                                "observed": proj["values"][n],
                                "reference_provenance": provenance[n]})
    return {
        "outcome": "RED" if divergences else "GREEN",
        "refusals": [],
        "divergences": divergences,
        "band_authority": band_report,
        "ignored_above_band": proj["ignored_above_band"],
        "compared_levels": n_max + 1,
        "transform_injected": transform_override is not None,
    }


def _transform_sanity(packet_data):
    """Section 8 fixture requirements on the packet's own transform."""
    tf = packet_data["indexing_map"]["index_transform"]
    a, b = tf["a"], tf["b"]
    n_max = packet_data["indexing_map"]["certified_band"]["n_max"]
    nonidentity = not (a == 1 and b == 0)
    # derangement over K_band: no k with a*k + b = k, i.e. k*(a-1) = -b
    fixed = [k for k, _ in packet_data["reference_values"] if a * k + b == k]
    values = [m for _, m in packet_data["reference_values"]]
    distinct = len(set(values)) == len(values)
    return {"nonidentity": nonidentity, "derangement": not fixed,
            "values_pairwise_distinct": distinct, "n_max": n_max}


def G5a_map_consumption(ingested, artifact, wrong_transform):
    """G5a rebuilt per addendum 12.1 Q3: three arms, fault-injected.

    Ordering, not a validation count: the fixture packet's validation
    COMPLETES CLEAN BEFORE any fault is injected (compare_3a re-validates the
    packet on every arm, and every arm validates the same clean object, since
    the packet itself is never mutated); the identity arm and the
    preregistered-wrong arm then INJECT their transform at the application
    point via `transform_override`, downstream of that validation.  An arm
    that produces STRUCTURAL_REFUSAL fails the gate: a refusal tests the
    validator, not the comparison, and proves the wrong proposition.
    """
    if not getattr(ingested, "hash_verified", False):
        return {"gate": "G5a", "pass": False,
                "reason": "refused: ingestion did not hash-verify"}

    sanity = _transform_sanity(ingested.data)
    identity = {"kind": "affine", "a": 1, "b": 0}

    arms = {}
    for name, override in (("supplied", None), ("identity", identity),
                           ("preregistered_wrong", wrong_transform)):
        res = compare_3a(ingested.data, artifact, transform_override=override)
        div = res["divergences"]
        arms[name] = {
            "outcome": res["outcome"],
            "divergences": len(div),
            "example": (f"level {div[0]['level']}: reference {div[0]['reference']}"
                        f" vs observed {div[0]['observed']}") if div else None,
            "hash_verified_throughout": True,
            "injected": override is not None,
        }

    passed = (arms["supplied"]["outcome"] == "GREEN"
              and arms["identity"]["outcome"] == "RED"
              and arms["preregistered_wrong"]["outcome"] == "RED"
              and sanity["nonidentity"] and sanity["derangement"]
              and sanity["values_pairwise_distinct"])

    return {
        "gate": "G5a",
        "fixture_transform_nonidentity": sanity["nonidentity"],
        "fixture_is_derangement": sanity["derangement"],
        "reference_values_pairwise_distinct": sanity["values_pairwise_distinct"],
        "transform_application_proved": arms["identity"]["outcome"] == "RED",
        "transform_correctness_proved": (
            arms["supplied"]["outcome"] == "GREEN"
            and arms["preregistered_wrong"]["outcome"] == "RED"),
        "no_arm_degraded_to_refusal": all(
            a["outcome"] in ("GREEN", "RED") for a in arms.values()),
        "arms": arms,
        "pass": bool(passed),
        "credited": ("the map-processing path controls the comparison; packet "
                     "validation completed clean before every injection, the "
                     "packet itself was never mutated, and the hash gate is "
                     "not credited here"),
    }
