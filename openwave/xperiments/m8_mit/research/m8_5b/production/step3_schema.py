"""The Step-3 raw-output contract: an implementation of the FROZEN section 5 schema.

NOT A DESIGN.  Section 5 freezes the core output schema and states that both
routes emit these fields with these meanings "so the two result surfaces are
comparable by construction rather than by later reconciliation.  The pilot may
add route-specific diagnostic fields; it may not change, rename, or drop these,
and their meanings are fixed at freeze."  This module therefore transcribes the
field list and enforces it; it invents nothing.

NULL SEMANTICS (frozen, section 5), three distinct states that must never be
collapsed:

    null            NOT APPLICABLE under that rung's definition
    numeric zero    a COMPUTED MATHEMATICAL ZERO
    omission        a STRUCTURAL FAILURE

The distinction is load-bearing: outside `Gamma = 2I` the McKay signature does
not apply, so `sector_signature` carries null rather than a fabricated value,
while a genuinely empty sector carries 0.  A validator that treated missing keys
as null would destroy exactly the evidence full-band completeness exists to
preserve.

FULL-BAND COMPLETENESS (section 5, section 4).  "Every level-by-sector cell
through the common certification band exists as a record, including
zero-multiplicity cells; an absent cell is a structural failure, not an implicit
zero."

ON `gate_results`, AND WHY IT IS NOT THE EVIDENCE.  The frozen schema includes a
`gate_results` field, so the artifact carries the producer's own gate outcomes.
That field is a CROSS-CHECK, never the credited evidence.  The artifact also
carries the raw measurements every frozen gate consumes
(`integer_nearness_margin`, `eigenvalue_residual`, `subspace_residual`,
`degeneracy_splitting`, `convergence_statistic`, `measured_rank`,
`quotient_multiplicity`, ...), so the adjudicator RE-APPLIES the frozen gates to
the measurements and may then compare its verdict against the producer's
recorded one.  An adjudicator that read `gate_results` and stopped would be
letting the artifact certify itself.  `recoverable_gate_inputs` below is the
instrument that proves the measurements suffice on their own.
"""

import hashlib
import json

__all__ = ["SCHEMA_VERSION", "CORE_FIELDS", "canonical_bytes", "artifact_sha256",
           "validate_records", "expected_cells", "recoverable_gate_inputs", "GATE_INPUTS",
           "write_artifact", "verify_artifact_file", "NONDETERMINISTIC_FIELDS",
           "determinism_manifest", "project_deterministic"]

SCHEMA_VERSION = "m8_5b-v1"   # section 6.4, frozen: the output schema version string

# section 5, verbatim and in order.  Not to be changed, renamed or dropped.
CORE_FIELDS = (
    "schema_version",
    "route",
    "run_id",
    "configuration_id",
    "arena_case_id",
    "group_order",
    "rung",
    "sector_id",
    "form_degree",
    "hodge_sector",
    "connection_class",
    "harmonic_level",
    "eigenvalue_R2",
    "sector_signature",          # {dimension, mckay_distance} or null
    "restriction_multiplicity",
    "quotient_multiplicity",
    "measured_rank",
    "integer_nearness_margin",
    "eigenvalue_residual",
    "subspace_residual",
    "degeneracy_splitting",
    "convergence_statistic",
    "gate_results",
    "wall_time",
    "peak_memory",
)

# Which raw measurements each downstream gate must be able to reach WITHOUT the
# producer.  This is the deletion test expressed as data.
GATE_INPUTS = {
    "G3": ("hodge_sector", "harmonic_level", "eigenvalue_R2",
           "quotient_multiplicity", "subspace_residual"),
    "G4": ("quotient_multiplicity", "measured_rank", "integer_nearness_margin"),
    "G6": ("harmonic_level", "hodge_sector", "eigenvalue_R2",
           "quotient_multiplicity", "sector_signature", "rung"),
    "G7": ("degeneracy_splitting", "convergence_statistic", "configuration_id"),
    "G8": ("run_id", "configuration_id", "eigenvalue_R2", "measured_rank",
           "quotient_multiplicity", "convergence_statistic"),
}


def canonical_bytes(obj):
    """Canonical serialization: sorted keys, two-space indent, ASCII, LF, one
    trailing newline.  Step 3 is a commitment boundary, so two artifacts are the
    same object exactly when these bytes agree."""
    return (json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True)
            + "\n").encode("ascii")


def artifact_sha256(obj):
    return hashlib.sha256(canonical_bytes(obj)).hexdigest()


def write_artifact(path, obj):
    """Write canonical bytes and VERIFY the file back.  Returns the file's sha256.

    This exists because ad-hoc writes did not produce canonical bytes.  Every
    qualification artifact on disk was off by one trailing newline, in both
    directions: `canonical_bytes` already terminates with LF, so a script that
    appended another produced one byte too many, and a script that used
    `json.dump` without one produced one byte too few.  Either way the file hash
    was not the artifact hash, and section 5 defines artifact identity BY THOSE
    BYTES.  `step3_runner` was never affected; it hands `canonical_bytes` straight
    to its writer.

    The read-back is the point.  A write helper that only formats correctly can
    still be bypassed; one that fails closed on a mismatch cannot be bypassed
    silently.
    """
    raw = canonical_bytes(obj)
    with open(path, "wb") as fh:
        fh.write(raw)
    back = open(path, "rb").read()
    if back != raw:
        raise IOError(f"{path}: written bytes are not canonical")
    return hashlib.sha256(back).hexdigest()


def verify_artifact_file(path):
    """Do a file's bytes equal the canonical serialization of its own content?"""
    raw = open(path, "rb").read()
    canon = canonical_bytes(json.loads(raw.decode("ascii")))
    return {"path": path, "file_bytes": len(raw), "canonical_bytes": len(canon),
            "canonical": raw == canon,
            "file_sha256": hashlib.sha256(raw).hexdigest(),
            "canonical_sha256": hashlib.sha256(canon).hexdigest()}


def expected_cells(n_max):
    """The complete section 5 cell lattice through `n_max`, keyed by
    (harmonic_level, sector_id).

    Addendum 12.1 (IM11) exposed why the key must be sector_id and not
    hodge_sector: `oneform_coexact_up` and `oneform_coexact_down` share
    hodge_sector "coexact", so a (level, hodge_sector) key cannot see one of
    them missing while the other is present.  Section 5 defines the record
    surface as one record per (harmonic_level, sector_id); this helper and
    `validate_records` now use exactly that key.
    """
    sectors = ("scalar", "oneform_exact", "oneform_coexact_up",
               "oneform_coexact_down")
    return {(n, s) for n in range(0, n_max + 1) for s in sectors}


def validate_records(records, expected_cells=None):
    """Structural validation of a Step-3 record list.

    `expected_cells` is the set of (harmonic_level, sector_id) pairs the
    common certification band requires (see `expected_cells(n_max)` above).
    When given, missing cells are reported as STRUCTURAL FAILURES rather than
    treated as implicit zeros, and duplicate or out-of-lattice coordinates are
    likewise structural failures.
    """
    problems = []

    for i, rec in enumerate(records):
        missing = [f for f in CORE_FIELDS if f not in rec]
        if missing:
            problems.append(f"record {i}: OMITTED core fields {missing} "
                            f"(section 5: omission is a structural failure, "
                            f"not an implicit null)")
        extra_ok = True   # route-specific diagnostic fields are permitted
        if rec.get("schema_version") not in (None, SCHEMA_VERSION) and extra_ok:
            problems.append(f"record {i}: schema_version "
                            f"{rec.get('schema_version')!r} != {SCHEMA_VERSION!r}")

    if expected_cells is not None:
        keys = [(r.get("harmonic_level"), r.get("sector_id")) for r in records]
        present = set(keys)
        absent = sorted(set(expected_cells) - present,
                        key=lambda t: (t[0] if t[0] is not None else -1, str(t[1])))
        if absent:
            problems.append(
                f"FULL-BAND COMPLETENESS: {len(absent)} required level-by-sector "
                f"cell(s) absent, e.g. {absent[:4]}. Section 5: an absent cell is "
                f"a structural failure, not an implicit zero")
        dup = sorted({k for k in present if keys.count(k) > 1},
                     key=lambda t: (t[0] if t[0] is not None else -1, str(t[1])))
        if dup:
            problems.append(
                f"DUPLICATE CELLS: {len(dup)} level-by-sector coordinate(s) carry "
                f"more than one record, e.g. {dup[:4]}: the cell value is ambiguous")
        outside = sorted(present - set(expected_cells),
                         key=lambda t: (t[0] if t[0] is not None else -1, str(t[1])))
        if outside:
            problems.append(
                f"OUT-OF-LATTICE CELLS: {len(outside)} coordinate(s) outside the "
                f"required lattice, e.g. {outside[:4]}")

    # null vs zero must both occur somewhere for the distinction to be real;
    # a producer that emits only one of them has probably collapsed them.
    nulls = sum(1 for r in records for f in CORE_FIELDS if r.get(f, "absent") is None)
    zeros = sum(1 for r in records for f in CORE_FIELDS
                if isinstance(r.get(f), (int, float)) and r.get(f) == 0)

    return {
        "records": len(records),
        "problems": problems,
        "null_valued_fields": nulls,
        "computed_zero_fields": zeros,
        "pass": not problems,
    }


def recoverable_gate_inputs(records):
    """The deletion test as data: can each gate reach its inputs from the artifact?

    Answers, per gate, whether every field it consumes is PRESENT in every
    record.  A field carrying `null` counts as present, because null is a frozen
    meaning (not applicable) rather than an absence; a field that is missing
    does not.
    """
    out = {}
    for gate, fields in sorted(GATE_INPUTS.items()):
        missing = sorted({f for f in fields
                          for r in records if f not in r})
        out[gate] = {
            "required_fields": list(fields),
            "missing_from_some_record": missing,
            "recoverable_without_producer": not missing,
        }
    return out

# --- G8: the determinism guarantee, declared rather than assumed ---------------

# Section 8 makes bit-for-bit identity apply to "deterministic serialized
# summaries THE IMPLEMENTATION EXPLICITLY GUARANTEES".  This is that declaration.
# It travels inside the artifact so a reader of the bytes alone can see exactly
# what was promised.
#
# The guaranteed object is the artifact with these fields removed, serialized by
# `canonical_bytes`.  They are excluded because they are genuine resource
# measurements and vary between runs; emitting them as `null` would misuse
# section 5's NOT APPLICABLE, and dropping them would violate the frozen 25.
#
# G8 checks this declaration for MINIMALITY, not just sufficiency: the excluded
# set must EQUAL the set of fields that actually vary, so a producer cannot pass
# by declaring everything non-deterministic.
NONDETERMINISTIC_FIELDS = ("wall_time", "peak_memory")


def determinism_manifest():
    return {
        "guaranteed": ("bit-for-bit identity of this artifact with the listed "
                       "fields removed, under `canonical_bytes`"),
        "nondeterministic_fields": list(NONDETERMINISTIC_FIELDS),
        "reason": ("genuine resource measurements that vary between runs; "
                   "section 5 reserves null for NOT APPLICABLE, so they are "
                   "reported truthfully and excluded from the guarantee rather "
                   "than nulled or dropped"),
        "serializer": "step3_schema.canonical_bytes",
        "repeat_run_count": 3,
        "minimality_note": ("this list must EQUAL the set of fields that actually "
                            "vary across repeat runs; an over-broad declaration "
                            "is a G8 structural failure"),
    }


def project_deterministic(artifact):
    """The guaranteed summary: the artifact minus its declared varying fields."""
    import copy
    a = copy.deepcopy(artifact)
    a.pop("determinism", None)
    for r in a.get("records", []):
        for f in NONDETERMINISTIC_FIELDS:
            r.pop(f, None)
    return a
