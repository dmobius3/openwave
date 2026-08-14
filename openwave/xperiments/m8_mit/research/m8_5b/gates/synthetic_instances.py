"""Q2 synthetic Packet II instances and the G5a fixture (addendum 12.1.6).

ANCHORING (the co-design guard).  Multiplicity values for every lens-family
instance are imported at build time from the pilot's frozen tuning tables
(`pilot/route_a_quotient.py`, TUNING_REFERENCE), which reproduced published
scalar multiplicities exactly through k = 9 before the parameter freeze
(preregistration section 6.1).  The values therefore do not originate with the
unit that wrote the comparator.  This module contains no multiplicity literal
of its own.

LABELING.  Every instance carries a SYN- case_id (mechanically unsealable as a
production packet; packet_schema refuses at the sealing gate).  All instances
here are labeled CONVENTION-TRANSFORMED: their citation fields are synthetic,
non-evidentiary fixture data naming no real publication, and their transforms,
signs, and listing conventions are chosen to exercise schema degrees of
freedom.  No instance claims to represent what any published source states.
The label lives HERE and in the qualification record, never inside a packet.

The instance matrix spans: identity and non-identity transforms, a shift
b != 0, a stride a >= 2, both laplacian_sign values, zero rows listed and
unlisted, off-image classes empty and nonempty, and two distinct n_max values.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, os.path.join(_ROOT, "pilot"), os.path.join(_ROOT, "production")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from route_a_quotient import TUNING_REFERENCE               # noqa: E402
from adjudication_gates import STEP3_SCHEMA_VERSION         # noqa: E402

__all__ = ["instance_matrix", "packet_ii_for_case", "g5a_fixture_packet",
           "G5A_WRONG_TRANSFORM", "make_scalar_artifact", "SYN_CITATION"]

SYN_CITATION = {"authors": "SYN FIXTURE", "title": "CONVENTION-TRANSFORMED",
                "venue": "NONE", "year": 2026, "doi": "10.0000/syn",
                "table": "TUNING", "row": "PILOT"}


def _packet(case_id, a, b, sign, entries, n_max, off_image):
    """Assemble a conforming Packet II.  V5 coefficients are DERIVED from
    (a, b, sign) via the closure identities, which is legitimate here because
    these are convention-transformed fixtures whose 'source convention' is the
    transform itself; a production sealer transcribes A, B, C from the source
    and V5 checks them against the declared transform."""
    s = 1 if sign == "nonnegative" else -1
    A, B, C = s * a * a, s * 2 * a * (b + 1), s * b * (b + 2)
    return {
        "case_id": case_id,
        "citation": dict(SYN_CITATION),
        "indexing_map": {
            "index_transform": {"kind": "affine", "a": a, "b": b},
            "source_eigenvalue": {"form": "quadratic", "A": A, "B": B, "C": C},
            "laplacian_sign": sign,
            "radius_normalization": "unit_radius_dimensionless",
            "multiplicity_convention": {"counts": "per_protocol_level",
                                        "source_dimension_field": "real"},
            "unlisted_source_rows": "zero_multiplicity",
            "off_image_levels": off_image,
            "certified_band": {"n_max": n_max},
        },
        "reference_values": [[k, m] for k, m in entries],
        "format_version": "m8_5b-packet-II-2",
    }


def _band(case, n_max):
    """Pilot-anchored multiplicities for protocol levels 0..n_max."""
    table = TUNING_REFERENCE[case]
    if len(table) <= n_max:
        raise ValueError(f"tuning table for {case} shorter than band {n_max}")
    return table[:n_max + 1]


def instance_matrix():
    """The four Q2 instances.  Returns [(label, notes, packet), ...]."""
    l31 = _band("L(3,1)", 3)          # [1, 0, 3, 8], n_max = 3
    l21 = _band("L(2,1)", 4)          # [1, 0, 9, 0, 25], n_max = 4
    l41 = _band("L(4,1)", 4)          # [1, 0, 3, 0, 15], n_max = 4

    a_inst = _packet("SYN-Q2-A", 1, 0, "nonnegative",
                     [(n, m) for n, m in enumerate(l31)], 3, "empty")
    b_inst = _packet("SYN-Q2-B", 1, -2, "nonpositive",
                     [(n + 2, m) for n, m in enumerate(l31) if m != 0],
                     3, "empty")
    c_inst = _packet("SYN-Q2-C", 2, 0, "nonnegative",
                     [(k, l21[2 * k]) for k in range(0, 3)],
                     4, "spectrum_excludes")
    d_inst = _packet("SYN-Q2-D", 1, 0, "nonnegative",
                     [(n, m) for n, m in enumerate(l41) if m != 0],
                     4, "empty")
    return [
        ("SYN-Q2-A", "L(3,1) values; identity transform; zeros listed; "
                     "nonnegative; n_max 3", a_inst),
        ("SYN-Q2-B", "L(3,1) values; shift b=-2; nonpositive sign; zeros "
                     "unlisted; n_max 3", b_inst),
        ("SYN-Q2-C", "L(2,1) values; stride a=2; off-image class {1,3} "
                     "affirmed spectrum_excludes; n_max 4", c_inst),
        ("SYN-Q2-D", "L(4,1) values; identity transform; zeros unlisted; "
                     "n_max 4", d_inst),
    ]


def packet_ii_for_case(case_id, case="L(3,1)", n_max=3):
    """The rehearsal Packet II: SYN-Q2-A's shape bound to the rehearsal
    case_id, so the 3a comparison keys against the produced artifacts."""
    band = _band(case, n_max)
    return _packet(case_id, 1, 0, "nonnegative",
                   [(n, m) for n, m in enumerate(band)], n_max, "empty")


# --- the G5a fixture (section 8, rebuilt in the real schema per Q3) -----------

# Preregistered incorrect transform for the third arm: a wrong shift within the
# admitted family.  Injected downstream of validation (Q3 layering).
G5A_WRONG_TRANSFORM = {"kind": "affine", "a": 1, "b": -1}


def g5a_fixture_packet():
    """The fixture packet: STRUCTURED NONIDENTITY transform (shift b = -2),
    pairwise-distinct reference values, derangement over K_band.  L(3,1)
    values under the shift, zeros unlisted so every listed value is distinct."""
    l31 = _band("L(3,1)", 3)
    entries = [(n + 2, m) for n, m in enumerate(l31) if m != 0]   # [(2,1),(4,3),(5,8)]
    return _packet("SYN-G5A-FIX", 1, -2, "nonnegative", entries, 3, "empty")


def make_scalar_artifact(case_id, values_by_level, claimed_n_max=None,
                         schema_version=STEP3_SCHEMA_VERSION, rung="3a"):
    """A synthetic section 5 SCALAR artifact for harness fixtures (G5a's
    observed side).  Not a route output and never promoted to one: it carries
    a synthetic run_id and a SYN case_id.  Only the fields the 3a projection
    consumes are populated meaningfully."""
    levels = sorted(values_by_level)
    n_max = claimed_n_max if claimed_n_max is not None else max(levels)
    records = []
    for n in levels:
        m = values_by_level[n]
        records.append({
            "schema_version": schema_version,
            "route": "synthetic-fixture",
            "run_id": f"g5a-fixture-{case_id}",
            "configuration_id": "m8_5b-v1-fixture",
            "arena_case_id": case_id,
            "group_order": None,
            "rung": rung,
            "sector_id": "scalar",
            "form_degree": 0,
            "hodge_sector": None,
            "connection_class": None,
            "harmonic_level": n,
            "eigenvalue_R2": float(n * (n + 2)) if m > 0 else None,
            "sector_signature": None,
            "restriction_multiplicity": None,
            "quotient_multiplicity": float(m),
            "measured_rank": None,
            "integer_nearness_margin": 0.0,
            "eigenvalue_residual": None,
            "subspace_residual": None,
            "degeneracy_splitting": None,
            "convergence_statistic": None,
            "gate_results": None,
            "wall_time": 0.0,
            "peak_memory": 0,
        })
    return {"artifact_kind": "synthetic-g5a-observed", "n_max": n_max,
            "records": records}
