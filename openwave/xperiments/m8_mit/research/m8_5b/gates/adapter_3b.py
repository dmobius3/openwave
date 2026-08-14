"""The rung-3b production adapter, implemented against addendum 12.1.5.

Bridges the frozen evaluator's output surface (rows of `k`, `lower_branch`,
`upper_branch`) to the frozen section 5 record surface, with NO adapter
discretion anywhere:

    exact    oneform_exact[n = k].quotient_multiplicity
                 <->  lower_branch.multiplicity of row k,        n = 1..n_max
    coexact  T(M) = m_up(M-2) + m_down(M)
                 <->  upper_branch.multiplicity of row k = M-1,  M = 2..n_max

The coexact comparison lives on the PHYSICAL EIGENVALUE surface lambda*R^2 =
M^2: the evaluator emits per-eigenvalue totals (its Gamma = 1 gate pins the
upper branch against the full 2(M^2 - 1) tower), while section 5 carries one
record per branch, each half of the tower.  Comparing any single branch cell
against the total is the uniform factor-2 hazard (H1); doubling one branch is
H2; both are mandatory battery mutations, injected here through
`coexact_total_fn`, which exists for the battery only and is never passed in
production.

Outcome vocabulary matches the 3a comparator and is never merged:
STRUCTURAL_REFUSAL (malformed or incomplete required support; no verdict),
RED (divergences found; adjudication content), GREEN.

The full section 5 lattice, including cells this adapter never consumes (the
level-0 exact null cell, the out-of-range coexact_down null cells, the scalar
sector), is the Step-3 artifact validator's obligation
(`step3_schema.validate_records` with `expected_cells(n_max)`), demonstrated
on the integrated path by Q4 and exercised by IM11.  This adapter neither
requires nor rejects those cells.
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (_ROOT, os.path.join(_ROOT, "production")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from adjudication_gates import integer_near, STEP3_SCHEMA_VERSION   # noqa: E402

__all__ = ["accept_evaluator_rows", "project_oneform_observed", "compare_3b",
           "RUN_CONFIGURATION"]

# 12.1.5, pinned: the step-6 runner invokes the evaluator at exactly these.
# "as_printed" exists only inside the evaluator's mutation battery.
RUN_CONFIGURATION = {"p": 1, "mapping": "corrected"}

_ROW_KEYS = {"k", "lower_branch", "upper_branch"}
_BRANCH_KEYS = {"eigenvalue", "M_Gamma_index", "multiplicity"}


def _is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def accept_evaluator_rows(rows, n_max):
    """The evaluator-side acceptance predicate (structural; violations refuse).

    Returns (refusals, consumed, ignored_above_band) where `consumed` maps
    k -> {"exact": int, "coexact_total": int} over the required support
    k = 1..n_max.
    """
    refusals, consumed = [], {}
    ignored = 0
    for i, row in enumerate(rows):
        if not isinstance(row, dict) or set(row) != _ROW_KEYS:
            refusals.append(f"evaluator row {i}: key set != {sorted(_ROW_KEYS)}")
            continue
        k = row["k"]
        if not _is_int(k):
            refusals.append(f"evaluator row {i}: k {k!r} is not a JSON integer")
            continue
        if k <= 0:
            refusals.append(f"evaluator row {i}: k = {k} <= 0 is malformed "
                            "(the evaluator's rows begin at k = 1)")
            continue
        if k > n_max:
            ignored += 1
            continue
        if k in consumed:
            refusals.append(f"duplicate evaluator row at k = {k}: the value is "
                            "ambiguous")
            continue
        branches = {}
        row_bad = False
        for bname, want_eig, want_idx in (
                ("lower_branch", k * (k + 2), 1),
                ("upper_branch", (k + 1) ** 2, 2)):
            br = row[bname]
            if not isinstance(br, dict) or set(br) != _BRANCH_KEYS:
                refusals.append(f"evaluator row k={k} {bname}: key set != "
                                f"{sorted(_BRANCH_KEYS)}")
                row_bad = True
                continue
            if br["eigenvalue"] != want_eig:
                refusals.append(
                    f"evaluator row k={k} {bname}: eigenvalue {br['eigenvalue']!r} "
                    f"!= correspondence-law value {want_eig} (this is the check "
                    "that catches adapter k-bookkeeping errors)")
                row_bad = True
            if br["M_Gamma_index"] != want_idx:
                refusals.append(
                    f"evaluator row k={k} {bname}: M_Gamma_index "
                    f"{br['M_Gamma_index']!r} != {want_idx}: the output does not "
                    "attest the corrected Theorem 3.3 mapping")
                row_bad = True
            if not (_is_int(br["multiplicity"]) and br["multiplicity"] >= 0):
                refusals.append(f"evaluator row k={k} {bname}: multiplicity "
                                f"{br['multiplicity']!r} is not a JSON integer >= 0")
                row_bad = True
            branches[bname] = br
        if row_bad:
            continue
        consumed[k] = {"exact": branches["lower_branch"]["multiplicity"],
                       "coexact_total": branches["upper_branch"]["multiplicity"]}
    missing = [k for k in range(1, n_max + 1) if k not in consumed]
    if missing:
        refusals.append(f"evaluator rows missing for k = {missing}: required "
                        "support absent")
    return refusals, consumed, ignored


_SECTOR_META = {
    "oneform_exact": "exact",
    "oneform_coexact_up": "coexact",
    "oneform_coexact_down": "coexact",
}


def project_oneform_observed(artifact, case_id, n_max):
    """The route-side acceptance and lookup (structural; violations refuse).

    Consumed cells, exactly: (n, oneform_exact) for n = 1..n_max;
    (M-2, oneform_coexact_up) and (M, oneform_coexact_down) for M = 2..n_max.
    Returns (refusals, cells) with cells keyed (harmonic_level, sector_id).
    """
    need = {(n, "oneform_exact") for n in range(1, n_max + 1)}
    need |= {(M - 2, "oneform_coexact_up") for M in range(2, n_max + 1)}
    need |= {(M, "oneform_coexact_down") for M in range(2, n_max + 1)}

    refusals, cells = [], {}
    for i, rec in enumerate(artifact.get("records", [])):
        sid = rec.get("sector_id")
        if sid not in _SECTOR_META:
            continue
        key = (rec.get("harmonic_level"), sid)
        if key not in need:
            continue                      # unconsumed lattice: the validator's job
        if rec.get("arena_case_id") != case_id:
            refusals.append(f"record {i} {key}: arena_case_id "
                            f"{rec.get('arena_case_id')!r} != {case_id!r}")
            continue
        if rec.get("schema_version") != STEP3_SCHEMA_VERSION:
            refusals.append(f"record {i} {key}: schema_version "
                            f"{rec.get('schema_version')!r} != "
                            f"{STEP3_SCHEMA_VERSION!r}")
            continue
        if rec.get("form_degree") != 1 or rec.get("rung") != "3b" \
                or rec.get("hodge_sector") != _SECTOR_META[sid]:
            refusals.append(
                f"record {i} {key}: identity metadata contradiction "
                f"(form_degree {rec.get('form_degree')!r}, rung "
                f"{rec.get('rung')!r}, hodge_sector {rec.get('hodge_sector')!r})")
            continue
        if key in cells:
            refusals.append(f"duplicate record at consumed coordinate {key}")
            continue
        v = integer_near(rec.get("quotient_multiplicity"))
        if v is None:
            refusals.append(
                f"record {i} {key}: quotient_multiplicity "
                f"{rec.get('quotient_multiplicity')!r} is null, non-numeric, "
                "negative, or non-integral within 1e-6")
            continue
        cells[key] = v
    missing = sorted(need - set(cells))
    if missing:
        refusals.append(f"consumed cells missing: {missing}")
    return refusals, cells


def compare_3b(rows, artifact, case_id, n_max, coexact_total_fn=None):
    """The rung-3b comparison per 12.1.5.  Deterministic: one result per
    conforming (evaluator output, route artifact) pair.

    `coexact_total_fn(cells, M) -> int` is the battery's fault-injection port
    for IM3/IM4 (single branch; doubling shortcut).  Production passes None and
    gets the normative T(M) = m_up(M-2) + m_down(M).
    """
    ev_refusals, ev, ev_ignored = accept_evaluator_rows(rows, n_max)
    rt_refusals, cells = project_oneform_observed(artifact, case_id, n_max)
    refusals = ev_refusals + rt_refusals
    if refusals:
        return {"outcome": "STRUCTURAL_REFUSAL", "refusals": refusals,
                "exact_divergences": [], "coexact_divergences": []}

    if coexact_total_fn is None:
        def coexact_total_fn(c, M):
            return (c[(M - 2, "oneform_coexact_up")]
                    + c[(M, "oneform_coexact_down")])
        injected = False
    else:
        injected = True

    exact_div = []
    for n in range(1, n_max + 1):
        obs = cells[(n, "oneform_exact")]
        ref = ev[n]["exact"]
        if obs != ref:
            exact_div.append({"level": n, "evaluator": ref, "observed": obs})

    coexact_div = []
    for M in range(2, n_max + 1):
        obs = coexact_total_fn(cells, M)
        ref = ev[M - 1]["coexact_total"]
        if obs != ref:
            coexact_div.append({"eigenvalue_M2": M * M, "M": M,
                                "evaluator_total": ref, "observed_total": obs})

    red = bool(exact_div or coexact_div)
    return {
        "outcome": "RED" if red else "GREEN",
        "refusals": [],
        "exact_divergences": exact_div,
        "coexact_divergences": coexact_div,
        "evaluator_rows_ignored_above_band": ev_ignored,
        "aggregation_injected": injected,
        "run_configuration_required": RUN_CONFIGURATION,
        "independence_note": ("class-2 theorem evaluation; against route (b) "
                              "this is a theorem-consistency check overlapping "
                              "its character machinery, never a third "
                              "independent derivation (section 4.2)"),
    }
