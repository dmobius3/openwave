"""Route (a) Step-3 producer: observations only, over a fixed coordinate band.

WHAT THIS IS NOT.  It is not `run_step8.py` productionized.  That pilot script
did three things forbidden at Step 3, and each is structurally impossible here:

    1. it imported `levels_routeb` and took route (b)'s predicted totals;
    2. it used the stopping rule as a computation cutoff (`stop = lam; break`);
    3. it skipped levels where route (b) predicted nothing (`if tot == 0: continue`),
       letting route (b)'s support decide which route-(a) records exist.

THE LOAD-BEARING ARCHITECTURAL RULE.  The producer computes the ENTIRE required
section 6.2 band without consulting the achieved-band stopping rule.  `RESOLVED`,
`UNRESOLVED` and the achieved-band cutoff are ADJUDICATION-layer outputs derived
only after BOTH Step-3 artifacts are committed.  Section 6.3's stopping rule has
a limb -- "its multiplicity equals the cross-route predicted multiplicity" --
that no single route can evaluate, so a producer that emitted `RESOLVED` would
either be consuming route (b) or fabricating.

Section 6.2: the required band "does not shrink in response to a numerical
failure", and "nothing in the stopping rule may redefine the required band
downward to produce a smaller passing one".

CENSUS FAILS CLOSED.  The coordinate set is computed BEFORE the solve from
n_max and the frozen per-sector eigenvalue maps.  The emitted key set must equal
it exactly; omission AND duplication both raise.  Completeness is therefore not
discovered downstream.

ROUTE-LOCALITY.  Imports: route (a)'s own geometry, route (a)'s own
representation module, and the G1/G1b gate.  Nothing from route (b), nothing
from `levels_routeb`, no packet.

`measured_rank` IS A BRANCH RANK, NOT A CLUSTER DEGENERACY COUNT.  Section 6.1
sends TWO coexact branches to one eigenvalue: `coexact_up` at level `n` and
`coexact_down` at level `n + 2` both sit at `lambda = (n+2)^2`.  Counting the
numerical cluster therefore reports the SUM of the two branches under each
branch's name.  The left-action Casimir separates them, since it acts as
`n(n+2)` on the `V_n(u)` factor that section 6.1's `n` labels, so

    measured_rank_(n, branch) = dim ker( C_L|E_lambda - n(n+2) I )

and `dim E_lambda` is recorded separately as a diagnostic.  The old behaviour is
retained as a required negative mutation; it is a real historical defect.  Scalar
and exact cells are unaffected, because no two branches ever share an eigenvalue
there, which is why the defect was invisible in those sectors.

THE DISCRETISATION IS THE ORBIT-TRANSPORT BACKEND, NOT INDEPENDENT KNN.
`route_a_twosided.operators` selects every node's stencil independently, and a
distance tie at the stencil boundary is then resolved by an implementation
choice that is not deck-equivariant.  On L(7;1,2) at the frozen 60-seed rung that
broke equivariance at 1.86e-03 for the scalar operator and 1.79e-03 for the Hodge
operator, while the cloud itself was exactly invariant.  `operators` is still
imported, ONLY so the qualification harness can run it as the structural
mutation; the producer does not call it.  See QUOTIENT_BACKEND_REPAIR.md.
"""

import sys
import time

import numpy as np

sys.path.insert(0, "../pilot")
sys.path.insert(0, "../gates")

import route_a_repn as repn                                    # noqa: E402
import cl_descent as cl                                         # noqa: E402
import route_a_oneform_repn as oneform                          # noqa: E402
import route_gates as gates                                    # noqa: E402
from route_a_twosided import cloud, relax, operators, reduce_scalar, reduce_oneform  # noqa: E402
from equivariant_oneform import build_operators_equivariant     # noqa: E402
import step3_schema as schema                                  # noqa: E402

CLUSTER_WINDOW = 0.35          # section 6.3, ASSIGNMENT ONLY
ZERO_MODE_ABS = 1e-9           # section 6.3
ROUTE = "a"


# --- the fixed coordinate band, built before any solving ---------------------

def sector_coordinates(n_max):
    """Every required (harmonic_level, sector_id) coordinate through n_max.

    `n` labels the SCALAR factor throughout (section 6.1).  Each sector derives
    its candidate eigenvalue from its OWN frozen map, because section 6.2 states
    n_max and a single eigenvalue ceiling are not interchangeable across sectors.
    """
    coords = []
    for n in range(0, n_max + 1):
        coords.append({"harmonic_level": n, "sector_id": "scalar",
                       "form_degree": 0, "hodge_sector": None,
                       "eigenvalue_R2": n * (n + 2), "in_range": True})
        coords.append({"harmonic_level": n, "sector_id": "oneform_exact",
                       "form_degree": 1, "hodge_sector": "exact",
                       "eigenvalue_R2": n * (n + 2), "in_range": n >= 1})
        coords.append({"harmonic_level": n, "sector_id": "oneform_coexact_up",
                       "form_degree": 1, "hodge_sector": "coexact",
                       "eigenvalue_R2": (n + 2) ** 2, "in_range": True})
        coords.append({"harmonic_level": n, "sector_id": "oneform_coexact_down",
                       "form_degree": 1, "hodge_sector": "coexact",
                       "eigenvalue_R2": n * n, "in_range": n >= 2})
    return coords


def _cluster(evals, lam):
    """Members within the frozen 0.35 absolute window of a candidate level."""
    return np.asarray([e for e in evals if abs(e - lam) <= CLUSTER_WINDOW])


# --- the producer ------------------------------------------------------------

def produce(pairs, run_id, configuration_id, arena_case_id,
            seeds=120, relax_iters=30, stencil=110, poly_m=7, poly_deg=4,
            rng_seed=20260811, adjudication=False):
    """Drive route (a) over the full required band and emit section 5 records."""
    t0 = time.time()

    # 1. certify the action, and bind everything downstream to that group
    verdict, elems, digest = gates.certified_effective_group(pairs)
    if not verdict["pass"]:
        raise ValueError(f"G1b did not pass on this action: {verdict}")

    # 2. n_max, route-A-owned, character-free, from the same certified action
    nm = repn.derive_n_max(pairs)
    n_max = nm["n_max"]

    # 3. the coordinate census, BEFORE the solve
    coords = sector_coordinates(n_max)
    expected_keys = {(c["harmonic_level"], c["sector_id"]) for c in coords}
    if len(expected_keys) != len(coords):
        raise ValueError("coordinate census contains duplicates")

    # 4. the numerics, over the whole band, with no early exit
    rng = np.random.default_rng(rng_seed)
    S = rng.normal(size=(seeds, 4))
    S /= np.linalg.norm(S, axis=1, keepdims=True)
    S = relax(S, pairs, relax_iters)
    X, oid, gid, M = cloud(S, pairs)
    # orbit transport, NOT independent kNN: see the header and the repair record
    Lpos, H = build_operators_equivariant(X, oid, gid, pairs,
                                          stencil, poly_m, poly_deg)
    scal = np.sort(np.linalg.eigvals(reduce_scalar(Lpos, oid, M)).real)
    H_red = reduce_oneform(H, oid, gid, pairs, M, len(X))
    onef = np.sort(np.linalg.eigvals(H_red).real)

    # The LEFT-action Casimir, reduced through the same map.  It is what separates
    # the two section 6.1 coexact branches that share one eigenvalue, and it is the
    # difference between `measured_rank` and a cluster degeneracy count.
    scal_red = reduce_scalar(Lpos, oid, M)
    CL_red = reduce_oneform(cl.build_CL(Lpos), oid, gid, pairs, M, len(X))
    sev, sV = np.linalg.eig(scal_red)
    hev, hV = np.linalg.eig(H_red)

    # 5. fill EVERY coordinate; emptiness has no control-flow power
    records, emitted = [], set()
    for c in coords:
        evals = scal if c["form_degree"] == 0 else onef
        lam = c["eigenvalue_R2"]
        members = _cluster(evals, lam) if c["in_range"] else np.asarray([])

        # `k` is the CLUSTER DEGENERACY COUNT: dim E_lambda.  It is NOT
        # `measured_rank`, and conflating the two was the defect this corrects.
        k = int(len(members))
        mean = float(members.mean()) if k else None
        lo = float(members.min()) if k else None
        hi = float(members.max()) if k else None
        if k and lam != 0:
            eig_res = abs(mean - lam) / abs(lam)
        elif k:
            eig_res = abs(mean)
        else:
            eig_res = None
        spread = ((hi - lo) / abs(mean)) if (k and mean not in (None, 0.0)
                                             and lam != 0) else None

        # invariant dimension from route (a)'s OWN character-free machinery.
        # Scalars: stacked nullity on the level carrier.  One-forms: rank INSIDE
        # the certified sector projector, never an apportioned total.
        # An in-range sector that is genuinely empty yields numeric 0, NOT null:
        # section 5 reserves null for NOT APPLICABLE.
        n_lv = c["harmonic_level"]
        if not c["in_range"]:
            qm, cond = None, {"gap": None, "state": "not applicable at this level"}
        elif c["sector_id"] == "scalar":
            dim, _, cond = repn.invariant_dim_and_basis(pairs, n_lv)
            qm = float(dim)
        else:
            m_of = {"oneform_exact": n_lv,
                    "oneform_coexact_up": n_lv + 2,
                    "oneform_coexact_down": n_lv - 2}[c["sector_id"]]
            qm = float(oneform.sector_invariant_rank(pairs, n_lv, m_of))
            cond = {"gap": None, "state": f"sector rank inside certified P_m, m={m_of}"}

        # `measured_rank`: the BRANCH rank, dim ker( C_L|E_lambda - n(n+2) I ),
        # measured on the numerical objects with no reference to `qm`.  Two
        # independent readings are computed and BOTH are recorded; they must agree,
        # and a disagreement is a convergence signal rather than something to
        # resolve by loosening a threshold.
        if not c["in_range"]:
            branch, rank = None, 0
        else:
            ev, V, C = ((sev, sV, scal_red) if c["form_degree"] == 0
                        else (hev, hV, CL_red))
            b = cl.branch_ranks(C, ev, V, float(lam))
            d = b["branches"].get(n_lv)
            rank = int(d["rank_by_eigenvalue_window"]) if d else 0
            branch = {
                "left_casimir_target": float(n_lv * (n_lv + 2)),
                "cluster_degeneracy_count": k,
                "cluster_dim_rank_revealed": b["cluster_dim"],
                "rank_by_eigenvalue_window": (d["rank_by_eigenvalue_window"] if d else 0),
                "rank_by_svd_nullity": (d["rank_by_svd_nullity"] if d else 0),
                "readings_agree": (d["agree"] if d else True),
                "band_leak": cl.band_leak(C, ev, V, float(lam)),
                "note": ("measured_rank is the BRANCH rank; cluster_degeneracy_count "
                         "is dim E_lambda.  They differ exactly where two section 6.1 "
                         "branches share one eigenvalue"),
            }

        rec = {
            "schema_version": schema.SCHEMA_VERSION,
            "route": ROUTE,
            "run_id": run_id,
            "configuration_id": configuration_id,
            "arena_case_id": arena_case_id,
            "group_order": len(elems),
            # Addendum 12.1.4: adjudication-case records carry "3a" (scalar) or
            # "3b" (one-form); the hardcoded "4-6" was a disclosed frozen defect.
            "rung": (("3a" if c["form_degree"] == 0 else "3b") if adjudication
                     else ("4-6" if c["form_degree"] in (0, 1) else None)),
            "sector_id": c["sector_id"],
            "form_degree": c["form_degree"],
            "hodge_sector": c["hodge_sector"],
            "connection_class": None,          # untwisted
            "harmonic_level": c["harmonic_level"],
            "eigenvalue_R2": mean,             # cluster MEAN, section 6.3
            "sector_signature": None,          # McKay signature: not applicable off 2I
            "restriction_multiplicity": None,
            "quotient_multiplicity": qm,
            "measured_rank": rank,             # BRANCH rank, not the cluster count
            "integer_nearness_margin": (abs(qm - round(qm)) if qm is not None else None),
            "eigenvalue_residual": eig_res,
            "subspace_residual": None,         # frozen meaning: sin(theta_max); filled by G3 stage
            "degeneracy_splitting": spread,
            "convergence_statistic": None,     # one configuration in this run
            "gate_results": None,              # CROSS-CHECK only; adjudicator recomputes
            "wall_time": None,
            "peak_memory": None,

            # --- authorized route-(a) diagnostics (section 5 permits additions) ---
            "diag_route_a": {
                "candidate_eigenvalue_R2": lam,
                "in_required_range": c["in_range"],
                "cluster_lambda_min": lo,
                "cluster_lambda_max": hi,
                "cluster_members": [float(x) for x in members],
                "cluster_window": CLUSTER_WINDOW,
                "zero_mode_abs": ([abs(float(x)) for x in members]
                                  if lam == 0 and k else None),
                "invariant_dim_conditioning": cond,
                "branch_separation": branch,
                "effective_group_digest": digest,
            },
        }
        key = (rec["harmonic_level"], rec["sector_id"])
        if key in emitted:
            raise ValueError(f"duplicate coordinate emitted: {key}")
        emitted.add(key)
        records.append(rec)

    # 6. fail closed on the census
    if emitted != expected_keys:
        raise ValueError(
            f"coordinate census mismatch: missing {sorted(expected_keys - emitted)}, "
            f"unexpected {sorted(emitted - expected_keys)}")

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
        "effective_group_digest": digest,
        "group_order": len(elems),
        "n_max": n_max,
        "n_max_derivation": nm,
        "realization": repn.REALIZATION,
        "nodes": int(len(X)),
        "records": records,
        "adjudication_note": (
            "This artifact contains OBSERVATIONS ONLY. RESOLVED, UNRESOLVED and "
            "the achieved-band cutoff are adjudication-layer outputs derived at "
            "step 5 from BOTH committed artifacts; section 6.3's stopping rule "
            "has a cross-route limb no single route can evaluate."),
    }
