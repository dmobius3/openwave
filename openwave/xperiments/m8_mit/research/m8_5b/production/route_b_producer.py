"""Route (b) Step-3 producer: predictions only, over the same fixed coordinate band.

WHAT ROUTE (B) UNIQUELY OWES.  Its method is already productionized in
`route_b_core` and tied to the G1b-certified effective group, so this module is
mostly serialization and boundary enforcement:

    1. emit the FULL required 0..n_max coordinate surface, with computed zeros
       and inapplicable cells under the same section 5 semantics;
    2. derive its OWN n_max by the frozen character method, independently of
       route (a), because section 6.2 requires both routes to certify the same
       ceiling and that is only a check if each derives it;
    3. serialize EVERY route-(b) prediction the adjudicator will need, so no
       step-5 code ever imports or reruns `levels_routeb.py`;
    4. bind the computation to its own independently generated effective-group
       digest;
    5. populate only route-(b)-applicable fields.

WHY IT DOES NOT MANUFACTURE ROUTE-(A) FIELDS.  Section 6.1: "Route (b): closed
form, no solver."  Section 8 G7: route (b) "evaluates exact finite group sums, so
it has no convergence sequence to report and is never judged by the convergence
statistic."  So `measured_rank`, `subspace_residual`, `degeneracy_splitting` and
`convergence_statistic` carry NULL here, meaning NOT APPLICABLE, which is exactly
what section 5 reserves null for.  Filling them with plausible numbers would make
the two artifacts look comparable in places where they are not.

`eigenvalue_residual` is 0.0, a COMPUTED ZERO rather than null: route (b)
evaluates the frozen eigenvalue map in closed form, so its residual against that
map is genuinely zero, not inapplicable.

NO ADJUDICATION.  No `RESOLVED`, no achieved band, no cross-route comparison.
Section 6.3's stopping rule needs both artifacts and belongs to step 5.
"""

import sys
import time

import numpy as np

sys.path.insert(0, "../gates")

import route_b_core as rb                                   # noqa: E402
import route_gates as gates                                 # noqa: E402
import step3_schema as schema                               # noqa: E402

ROUTE = "b"
INTEGER_TOL = 1e-6          # section 6.3


def derive_n_max_character(pairs, search_to=40, tol=1e-6):
    """Section 6.2's ceiling, by route (b)'s OWN character averaging.

    Same frozen rule as route (a) uses -- the SECOND n > 0 with positive scalar
    invariant dimension -- but reached by the character method route (b) owns.
    Independence is the point: section 6.2 requires both routes to certify the
    same n_max, which is evidence only if neither took it from the other.
    """
    positive = []
    for n in range(1, search_to + 1):
        d = rb.scalar_invariant_dims(pairs, n)[n]
        if abs(d) > tol:
            positive.append(n)
            if len(positive) == 2:
                return {"n_max": positive[1], "first_two_positive_levels": positive,
                        "rule": "second n > 0 with positive scalar invariant dimension",
                        "derived_by": "route (b) two-sided character averaging"}
    raise ValueError("fewer than two positive invariant levels")


def sector_coordinates(n_max):
    """The same required coordinate set route (a) emits, derived independently
    from the frozen section 6.1 maps."""
    coords = []
    for n in range(0, n_max + 1):
        coords.append(("scalar", n, None, 0, None, n * (n + 2), True))
        coords.append(("oneform_exact", n, n, 1, "exact", n * (n + 2), n >= 1))
        coords.append(("oneform_coexact_up", n, n + 2, 1, "coexact", (n + 2) ** 2, True))
        coords.append(("oneform_coexact_down", n, n - 2, 1, "coexact", n * n, n >= 2))
    return coords


def produce(gen_pairs, run_id, configuration_id, arena_case_id, cutoff_factor=2,
            adjudication=False):
    t0 = time.time()

    verdict, elems, digest = gates.certified_effective_group(gen_pairs)
    if not verdict["pass"]:
        raise ValueError(f"G1b did not pass: {verdict}")

    # route (b) closes independently and REFUSES unless it matches G1b's group
    pairs, own_digest = rb.certified_average_group(gen_pairs, digest)

    nm = derive_n_max_character(pairs)
    n_max = nm["n_max"]

    coords = sector_coordinates(n_max)
    expected = {(c[1], c[0]) for c in coords}
    if len(expected) != len(coords):
        raise ValueError("coordinate census contains duplicates")

    records, emitted = [], set()
    for sid, n, m, deg, hodge, lam, in_range in coords:
        if not in_range:
            qm, margin, stab = None, None, None
        else:
            if sid == "scalar":
                raw = rb.scalar_invariant_dims(pairs, n)[n]
            else:
                raw = sum(rb.su2_character(u, n) * rb.su2_character(v, m)
                          for u, v in pairs) / len(pairs)
            qm = float(raw)
            margin = abs(raw - round(raw))
            # G7 for route (b): in-band cutoff stability, cutoff vs
            # cutoff_factor x cutoff, section 8 and section 4's "for every
            # observable".
            #
            # BOTH SIDES ARE RECOMPUTED THROUGH THE PIPELINE.  The previous form
            # assigned `hi = raw` for every one-form observable, so `stab` was
            # |raw - raw| and no defect in the one-form path could have moved it
            # off zero: a computed zero that had never been measured.  The scalar
            # limb was already honest and legitimately returns zero, because
            # `scalar_invariant_dims` evaluates each level independently and the
            # value at index n does not depend on the list length.  That zero is
            # earned by a real second evaluation; the one-form zero was not.
            hi_cut = cutoff_factor * (n_max + 1)
            if sid == "scalar":
                hi = rb.scalar_invariant_dims(pairs, hi_cut)[n]
            else:
                want = {"oneform_exact": "exact",
                        "oneform_coexact_up": "coexact",
                        "oneform_coexact_down": "coexact"}[sid]
                hi = next(r_["multiplicity_raw"]
                          for r_ in rb.oneform_levels(pairs, hi_cut)
                          if r_["n"] == n and r_["m"] == m and r_["sector"] == want)
            stab = float(abs(hi - raw))

        rec = {
            "schema_version": schema.SCHEMA_VERSION,
            "route": ROUTE,
            "run_id": run_id,
            "configuration_id": configuration_id,
            "arena_case_id": arena_case_id,
            "group_order": len(pairs),
            # Addendum 12.1.4: adjudication-case records carry "3a"/"3b"; the
            # hardcoded "4-6" was a disclosed frozen defect.
            "rung": ("3a" if deg == 0 else "3b") if adjudication else "4-6",
            "sector_id": sid,
            "form_degree": deg,
            "hodge_sector": hodge,
            "connection_class": None,
            "harmonic_level": n,
            "eigenvalue_R2": float(lam) if in_range else None,
            "sector_signature": None,           # McKay signature: not applicable off 2I
            "restriction_multiplicity": None,
            "quotient_multiplicity": qm,
            "measured_rank": None,              # NOT APPLICABLE: closed form, no solver
            "integer_nearness_margin": margin,
            "eigenvalue_residual": 0.0 if in_range else None,   # exact, a COMPUTED zero
            "subspace_residual": None,          # NOT APPLICABLE to route (b)
            "degeneracy_splitting": None,       # NOT APPLICABLE: no cluster
            "convergence_statistic": None,      # section 8 G7: route (b) has none
            "gate_results": None,               # adjudicator recomputes
            "wall_time": None,
            "peak_memory": None,

            # --- authorized route-(b) diagnostics ---
            "diag_route_b": {
                "m_index": m,
                "in_required_range": in_range,
                "method": rb.METHOD,
                "in_band_cutoff_stability": stab,
                "cutoff_factor": cutoff_factor,
                "effective_group_digest": own_digest,
            },
        }
        key = (n, sid)
        if key in emitted:
            raise ValueError(f"duplicate coordinate: {key}")
        emitted.add(key)
        records.append(rec)

    if emitted != expected:
        raise ValueError(f"census mismatch: missing {sorted(expected - emitted)}, "
                         f"unexpected {sorted(emitted - expected)}")

    wall = time.time() - t0
    for r in records:
        r["wall_time"] = round(wall, 4)

    return {
        "artifact_kind": "step3_raw_output",
        "route": ROUTE,
        "run_id": run_id,
        "configuration_id": configuration_id,
        "arena_case_id": arena_case_id,
        "schema_version": schema.SCHEMA_VERSION,
        "determinism": schema.determinism_manifest(),
        "effective_group_digest": own_digest,
        "group_order": len(pairs),
        "n_max": n_max,
        "n_max_derivation": nm,
        "method": rb.METHOD,
        "records": records,
        "adjudication_note": (
            "PREDICTIONS ONLY.  Every route (b) quantity the adjudicator needs is "
            "serialized here; step 5 must not import or rerun any route (b) "
            "prediction code.  RESOLVED, UNRESOLVED and the achieved band are "
            "step-5 outputs derived from BOTH committed artifacts."),
    }
