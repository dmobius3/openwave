#!/usr/bin/env python3
"""
Leanness and structural gate for the two sealed adjudication packets (M8.5-B § 4.1).

IDENTIFIERS (must match the pre-registration § 4.1 scheme exactly):
  S1..S6  structural, BUILD-TIME: free action, pair form, parameters reproduce the
          generators, closure order, freeness plus central canonicalization, and the
          element-order census.
  L1..L3  leanness, BUILD-TIME: unknown fields, prose, string length.
  G1/G1b in § 8 are RUN-TIME gates applied by each route to the action it executes.
          They overlap in what they test but are separate identifiers on purpose;
          do not renumber one to match the other.

Runs BEFORE sealing. Rejects any packet carrying a field the pre-registration does not
name, so "no working notes" is enforced mechanically rather than by inspection.

It also refuses free-text: every string value must look like data or a citation, not
prose. A note smuggled into an allowed field is exactly the failure this is for.

The mutation suite runs on SYNTHETIC fixtures with dummy values. This file never reads
the real case, so running the gate cannot expose the sealed material.
"""

import copy
import json
import math
import re
import sys

import numpy as np

SCHEMA = {
    "I": {"case_id", "family", "parameters", "generators", "action_convention",
          "format_version"},
    "II": {"case_id", "citation", "indexing_map", "reference_values",
           "format_version"},
}
CITATION_FIELDS = {"authors", "title", "venue", "year", "doi", "table", "row"}

# --- Packet II internal contract (addendum 12.1, adopted for Phase A) ---------
# Section 4.1's five-field list is frozen; everything below is internal
# structure of `indexing_map` and `reference_values` plus typing of `citation`.
PACKET_II_FORMAT_VERSION = "m8_5b-packet-II-2"
INDEXING_MAP_KEYS = {"index_transform", "source_eigenvalue", "laplacian_sign",
                     "radius_normalization", "multiplicity_convention",
                     "unlisted_source_rows", "off_image_levels", "certified_band"}
CASE_ID_RE = re.compile(r"[A-Z0-9][A-Z0-9-]{2,31}")
# Packet I `parameters`, frozen: the lens-space family carries exactly these
PARAMETER_KEYS = {"q", "s"}
SYN_PREFIX = "SYN-"      # reserved for unsealed qualification instances (Q2)

# A field may hold a short identifier or a citation string. Anything longer, or
# containing sentence punctuation, is prose and is rejected.
MAX_STR = 200
# Two rules, deliberately with different case sensitivity. The sentence rule is
# CASE-SENSITIVE: a lowercase word after a full stop is mid-sentence prose, whereas
# "A. Author" is an initial. Compiling it with re.I made [a-z] match uppercase and
# rejected every well-formed citation, which the negative control caught.
SENTENCE = re.compile(r"[.;]\s+[a-z]")
KEYWORDS = re.compile(r"\b(because|however|note that|rationale|reviewed|corrected|"
                      r"revised|round \d|see also|TODO|FIXME|deprecated)\b", re.I)


def judge(packet, which):
    """Return a list of violations. Empty list means lean enough to seal."""
    # R1/V1 totality: a packet is ONE JSON object.  A canonical array, string,
    # number, or null parses and hash-verifies fine, so the type guard must
    # live here, before any key operation, or a malformed packet crashes the
    # gate instead of being rejected (Redline audit blocker 1).
    if not isinstance(packet, dict):
        return [f"{which}: top-level value is {type(packet).__name__}, not an "
                "object: a packet is ONE JSON object (R1/V1)"]
    bad = []
    allowed = SCHEMA[which]
    keys = set(packet)
    for extra in sorted(keys - allowed):
        bad.append(f"unknown field {extra!r}: packets carry no fields beyond the manifest")
    for missing in sorted(allowed - keys):
        bad.append(f"missing required field {missing!r}")

    # Packet I's `parameters` is a nested object, and key exactness owns its key
    # set exactly as it owns the top level.  Without this, an extra key smuggled
    # inside `parameters` was accepted, and a missing `q` or `s` was caught only
    # incidentally when `structural_checks` later tried to read it.  A shape
    # defect must be refused BY the shape predicate, not by a downstream
    # IndexError: that is the same displaced-verification failure the S-battery
    # exists to refuse.
    if which == "I":
        params = packet.get("parameters")
        if not isinstance(params, dict):
            bad.append("parameters: must be an object with exactly the keys "
                       "{'q', 's'}")
        else:
            for extra in sorted(set(params) - PARAMETER_KEYS):
                bad.append(f"parameters: unknown key {extra!r}: the frozen "
                           "lens-space family carries exactly {'q', 's'}")
            for miss in sorted(PARAMETER_KEYS - set(params)):
                bad.append(f"parameters: missing required key {miss!r}")

    if which == "II" and isinstance(packet.get("citation"), dict):
        for extra in sorted(set(packet["citation"]) - CITATION_FIELDS):
            bad.append(f"unknown citation field {extra!r}")

    def scan(v, path):
        if isinstance(v, dict):
            for k, x in v.items():
                scan(x, f"{path}.{k}")
        elif isinstance(v, list):
            for i, x in enumerate(v):
                scan(x, f"{path}[{i}]")
        elif isinstance(v, str):
            if len(v) > MAX_STR:
                bad.append(f"{path}: string of {len(v)} chars exceeds {MAX_STR}, reads as prose")
            hit = SENTENCE.search(v) or KEYWORDS.search(v)
            if hit:
                bad.append(f"{path}: contains prose or a working note ({hit.group(0)!r})")
    scan(packet, which)
    return bad


# ---------------------------------------------------------------------------
# Structural checks. Separate concern from leanness: these ask whether the
# declared group IS the declared group, in the manner of the maintainers' M8.5-A
# element-order census, which separated 2I from A_5 x C_2 by counting order-4
# elements (30 against none) rather than appealing to a theorem.
#
# The B analogue is a cyclic deck group. Order alone does not pin it, and neither
# does order plus abelian-ness, so the census against phi(d) is what proves the
# group is cyclic of the declared order rather than merely the right size.
# ---------------------------------------------------------------------------

def _qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2, w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2, w1*z2 + x1*y2 - y1*x2 + z1*w2])


def _act_matrix(u, v):
    """4x4 orthogonal matrix of the isometry q -> u q v.

    Every deck group this protocol admits is represented this way, including the
    left-only (v = 1) and right-only (u = 1) cases. A single uniform representation
    is deliberate: an inhomogeneous lens space is NOT a left action, and a checker
    that assumed left multiplication would silently accept a packet its own routes
    could not consume.
    """
    E = np.eye(4)
    return np.array([_qmul(_qmul(u, E[i]), v) for i in range(4)]).T


def _close(mats, tol=1e-9, cap=4096):
    elems = [np.eye(4)]
    frontier = list(elems)
    while frontier:
        nxt = []
        for A in frontier:
            for Gm in mats:
                B = Gm @ A
                if not any(np.abs(B - E).max() < tol for E in elems):
                    elems.append(B); nxt.append(B)
                    if len(elems) > cap:
                        raise ValueError("closure cap exceeded")
        frontier = nxt
    return elems


def _order(A, tol=1e-9, cap=4096):
    B = A.copy()
    for n in range(1, cap + 1):
        if np.abs(B - np.eye(4)).max() < tol:
            return n
        B = A @ B
    return None


def canonical_action_pair(q, s):
    """The action pair [u, v] that the declared parameters (q, s) denote.

    The lens action on S^3 subset C^2 is (z1, z2) -> (zeta^s1 z1, zeta^s2 z2)
    with zeta = exp(2 pi i / q).  Write x = z1 + z2 j and take u, v as unit
    quaternions in the span of {1, i}.  Since j e^{i t} = e^{-i t} j,

        u x v = e^{i(tu + tv)} z1  +  e^{i(tu - tv)} z2 j

    so tu + tv = 2 pi s1 / q and tu - tv = 2 pi s2 / q, giving the HALF-ANGLES

        tu = pi (s1 + s2) / q            tv = pi (s1 - s2) / q

    which is the whole point.  The earlier form solved the same two equations
    with a modular inverse of 2, which exists only for ODD q, and the caller
    then SKIPPED S3 entirely for every even order: a required contract
    predicate silently uncovered over half the admissible class.  The
    half-angle form needs no inverse and is uniform over the class.

    For odd q the two agree up to the central identification [u, v] ~ [-u, -v]:
    writing s1 + s2 = 2a + m q and s1 - s2 = 2b + m' q, the angles differ by
    m pi and m' pi, and adding the two relations gives (m + m') q = 2(s1 - a - b),
    which is even, so m and m' have the same parity for odd q.  Both quaternions
    therefore flip together or not at all, and the induced 4x4 isometry is
    identical.  `s3_battery` verifies that numerically rather than asserting it,
    so nothing that passed under the old form can change verdict under this one.
    """
    tu = np.pi * (s[0] + s[1]) / q
    tv = np.pi * (s[0] - s[1]) / q
    return (np.array([np.cos(tu), np.sin(tu), 0.0, 0.0]),
            np.array([np.cos(tv), np.sin(tv), 0.0, 0.0]))


def lens_action_matrix(q, s):
    """The lens action written directly in COMPLEX coordinates, independent of
    the quaternion construction above.

    This exists so `canonical_action_pair` is checked against the contract's
    own definition of the action rather than against its own derivation. With
    `x = z1 + z2 j`, `z1 = a + b i`, `z2 = c + d i`, the quaternion basis
    `(1, i, j, k)` carries coordinates `(a, b, c, d)`, and
    `(z1, z2) -> (zeta^s1 z1, zeta^s2 z2)` is block-diagonal: a rotation by
    `2 pi s1 / q` on `(a, b)` and by `2 pi s2 / q` on `(c, d)`.
    """
    M = np.zeros((4, 4))
    for blk, sj in enumerate(s):
        t = 2 * np.pi * sj / q
        i0 = 2 * blk
        M[i0:i0 + 2, i0:i0 + 2] = np.array([[np.cos(t), -np.sin(t)],
                                            [np.sin(t), np.cos(t)]])
    return M


def _is_int_like(x):
    """An integer value, admitting the JSON float form 3.0 that json.load may
    produce, rejecting bools and genuinely fractional values."""
    if isinstance(x, bool):
        return False
    if isinstance(x, int):
        return True
    return isinstance(x, float) and float(x).is_integer()


def structural_checks(packet):
    """Return violations. Empty means the declared group is the declared group."""
    bad = []
    # S0 is a SHAPE predicate evaluated before anything indexes the parameters,
    # not an exception caught after the fact.  Key exactness is `judge`'s (see
    # PARAMETER_KEYS); S0 owns the types and arity.
    params = packet.get("parameters")
    if not isinstance(params, dict):
        return ["S0 parameters: must be an object"]
    if not (_is_int_like(params.get("q")) and int(params["q"]) >= 1):
        bad.append(f"S0 parameters.q: {params.get('q')!r} is not an integer >= 1")
    sraw = params.get("s")
    if not (isinstance(sraw, list) and len(sraw) == 2
            and all(_is_int_like(x) for x in sraw)):
        bad.append(f"S0 parameters.s: {sraw!r} is not a two-element array of "
                   "integers (the frozen family is lens spaces of S^3)")
    if bad:
        return bad
    try:
        q = int(params["q"])
        svec = [int(x) for x in sraw]
        gens = packet["generators"]
        conv = packet["action_convention"]
    except (KeyError, TypeError, ValueError) as e:
        return [f"S0 parameters/generators unreadable: {e}"]

    # S1 free action: a lens space needs gcd(s_i, q) = 1, else the action has
    # fixed points and the quotient is not a manifold.
    for i, sfac in enumerate(svec):
        if math.gcd(sfac % q, q) != 1:
            bad.append(f"S1 action not free: gcd(s[{i}]={sfac}, q={q}) != 1")

    if conv != "two_sided":
        bad.append(f"S2 action_convention {conv!r}: only 'two_sided' is admitted, "
                   "with left-only written as v = 1 (an inhomogeneous lens space is "
                   "not a left action)")
        return bad

    mats = []
    for g in gens:
        if not (isinstance(g, list) and len(g) == 2):
            bad.append("S2 each generator must be a pair [u, v] of unit quaternions")
            return bad
        u, v = np.array(g[0], float), np.array(g[1], float)
        for nm, w in (("u", u), ("v", v)):
            if abs(np.linalg.norm(w) - 1) > 1e-12:
                bad.append(f"S2 generator {nm} is not a unit quaternion")
        mats.append(_act_matrix(u, v))
    if bad:
        return bad

    # S3 the declared parameters must REPRODUCE the declared generator, so a
    # transcription slip between the two fields cannot survive.
    # TOTAL over the admissible class: no even/odd branch, no skip.  See
    # `canonical_action_pair` for why the earlier modular form covered odd
    # orders only, and `s3_battery` for the coverage evidence.
    if len(svec) != 2:
        bad.append(f"S3 parameters s has arity {len(svec)}: the frozen family "
                   "is lens spaces of S^3, where s carries exactly two entries")
    else:
        U, V = canonical_action_pair(q, svec)
        want = _act_matrix(U, V)
        if not any(np.abs(want - Gm).max() < 1e-9 for Gm in mats):
            bad.append(f"S3 generators do not reproduce parameters q={q}, s={svec}")

    # S4 closure order
    try:
        elems = _close(mats)
    except ValueError as e:
        return bad + [f"S4 {e}"]
    if len(elems) != q:
        bad.append(f"S4 closure order {len(elems)} != declared q={q}")
        return bad

    # S5a FREENESS, checked geometrically rather than inferred from parameters.
    # A 4x4 rotation fixes a point of S^3 exactly when it has eigenvalue +1, so a
    # nonidentity element with eigenvalue 1 makes the quotient an orbifold. The S1
    # gcd test only covers the declared parameters; this covers what is executed.
    for A in elems:
        if np.abs(A - np.eye(4)).max() < 1e-9:
            continue
        if np.abs(np.linalg.eigvals(A) - 1.0).min() < 1e-7:
            bad.append("S5a action not free: a nonidentity element fixes a point of S^3, "
                       "so the quotient is an orbifold rather than a manifold")
            break

    # S5b the pair representation is two-to-one onto the effective action, so a
    # packet that lists both [u,v] and [-u,-v] as distinct generators would
    # double-count the deck group while still closing correctly.
    seen = []
    for g in gens:
        uu, vv = np.array(g[0], float), np.array(g[1], float)
        for w in seen:
            if (np.abs(uu - w[0]).max() < 1e-9 and np.abs(vv - w[1]).max() < 1e-9) or \
               (np.abs(uu + w[0]).max() < 1e-9 and np.abs(vv + w[1]).max() < 1e-9):
                bad.append("S5b duplicate generator under the central identification "
                           "[u,v] ~ [-u,-v]: the deck group would be double-counted")
        seen.append((uu, vv))

    # S6 the element-order census: the check the A audit added.
    census = {}
    for A in elems:
        o = _order(A)
        census[o] = census.get(o, 0) + 1
    want = {}
    for d in range(1, q + 1):
        if q % d == 0:
            want[d] = sum(1 for t in range(1, d + 1) if math.gcd(t, d) == 1)
    if census != want:
        bad.append(f"S6 element-order census {dict(sorted(census.items()))} != "
                   f"cyclic C_{q} census {dict(sorted(want.items()))}")
    return bad


# ---------------------------------------------------------------------------
# Addendum 12.1 predicates V1-V8.  Each is a mechanical check on the parsed
# Packet II; together they are what "structurally ambiguous packets cannot be
# sealed" means.  They run at BUILD time here (sealing is blocked on any
# violation) and are re-run by the Step-5 comparator at ingestion, where a
# violation is a STRUCTURAL REFUSAL, never a rung-3a red.  One implementation
# serves both enforcement points on purpose: the contract is the addendum, the
# battery in main() is what makes this implementation falsifiable, and
# independence comes from the battery plus the Phase-B fresh sealer, not from
# duplicating the predicate code.
# ---------------------------------------------------------------------------

def _is_int(x):
    """JSON integer strictness (R4): bool is not an integer, 3.0 is not 3."""
    return isinstance(x, int) and not isinstance(x, bool)


def packet_ii_checks(packet):
    """V1-V8. Returns a list of violations; empty means structurally sound.

    Staged deliberately: shape and type predicates (V1-V3) guard the derived
    predicates (V4-V8), so one defect fires its own predicate rather than
    crashing a later one.
    """
    if not isinstance(packet, dict):
        return ["V1 II: top-level value is not a JSON object: the packet is "
                "ONE closed object shape (R1), and every non-object top-level "
                "value rejects here rather than crashing downstream"]
    bad = []

    def scan_null(v, path):
        if v is None:
            bad.append(f"V2 {path}: null is banned in packets (R5)")
        elif isinstance(v, dict):
            for k in v:
                scan_null(v[k], f"{path}.{k}")
        elif isinstance(v, list):
            for i, x in enumerate(v):
                scan_null(x, f"{path}[{i}]")
    scan_null(packet, "II")

    if packet.get("format_version") != PACKET_II_FORMAT_VERSION:
        bad.append(f"V2 format_version: must be exactly {PACKET_II_FORMAT_VERSION!r}")
    cid = packet.get("case_id")
    if not (isinstance(cid, str) and CASE_ID_RE.fullmatch(cid)):
        bad.append("V2 case_id: must match ^[A-Z0-9][A-Z0-9-]{2,31}$")

    cit = packet.get("citation")
    if not isinstance(cit, dict):
        bad.append("V1 citation: must be an object")
    else:
        for missing in sorted(CITATION_FIELDS - set(cit)):
            bad.append(f"V1 citation.{missing}: required key missing")
        for extra in sorted(set(cit) - CITATION_FIELDS):
            bad.append(f"V1 citation.{extra}: unknown key")
        if "year" in cit and not (_is_int(cit["year"]) and 1800 <= cit["year"] <= 2100):
            bad.append("V2 citation.year: JSON integer in [1800, 2100] required")
        if "doi" in cit and not (isinstance(cit["doi"], str) and cit["doi"].startswith("10.")):
            bad.append("V2 citation.doi: string beginning '10.' required")
        for f in ("table", "row"):
            if f in cit and not (isinstance(cit[f], str) and 1 <= len(cit[f]) <= 32):
                bad.append(f"V2 citation.{f}: string of 1..32 chars required")
        for f in ("authors", "title", "venue"):
            if f in cit and not isinstance(cit[f], str):
                bad.append(f"V2 citation.{f}: string required")

    im = packet.get("indexing_map")
    if not isinstance(im, dict):
        return bad + ["V1 indexing_map: must be an object"]
    for extra in sorted(set(im) - INDEXING_MAP_KEYS):
        bad.append(f"V1 indexing_map.{extra}: unknown key")
    for missing in sorted(INDEXING_MAP_KEYS - set(im)):
        bad.append(f"V1 indexing_map.{missing}: required key missing")

    rv = packet.get("reference_values")
    entries = None
    if not isinstance(rv, list) or not rv:
        bad.append("V2 reference_values: non-empty array of [k, m] pairs required")
    else:
        entries, ok_rv = [], True
        for i, e in enumerate(rv):
            if not (isinstance(e, list) and len(e) == 2
                    and _is_int(e[0]) and _is_int(e[1])):
                bad.append(f"V2 reference_values[{i}]: exactly [k, m], both JSON integers")
                ok_rv = False
                continue
            if e[1] < 0:
                bad.append(f"V2 reference_values[{i}]: m >= 0 required")
                ok_rv = False
            entries.append((e[0], e[1]))
        if ok_rv and any(entries[i][0] >= entries[i + 1][0]
                         for i in range(len(entries) - 1)):
            bad.append("V2 reference_values: k strictly increasing required")
            ok_rv = False
        if not ok_rv:
            entries = None

    a = b = None
    tf = im.get("index_transform")
    if isinstance(tf, dict):
        if set(tf) != {"kind", "a", "b"}:
            bad.append("V1 index_transform: keys exactly {kind, a, b}")
        if tf.get("kind") != "affine":
            bad.append(f"V3 index_transform.kind: {tf.get('kind')!r} not in the "
                       "closed set {'affine'}")
        elif _is_int(tf.get("a")) and _is_int(tf.get("b")):
            a, b = tf["a"], tf["b"]
            if a < 1:
                bad.append("V3 index_transform.a: integer >= 1 required")
                a = None
        else:
            bad.append("V2 index_transform: a and b must be JSON integers")
    elif "index_transform" in im:
        bad.append("V1 index_transform: must be an object")

    ABC = None
    se = im.get("source_eigenvalue")
    if isinstance(se, dict):
        if set(se) != {"form", "A", "B", "C"}:
            bad.append("V1 source_eigenvalue: keys exactly {form, A, B, C}")
        if se.get("form") != "quadratic":
            bad.append(f"V2 source_eigenvalue.form: {se.get('form')!r} not in the "
                       "closed set {'quadratic'}")
        elif all(_is_int(se.get(x)) for x in ("A", "B", "C")):
            ABC = (se["A"], se["B"], se["C"])
        else:
            bad.append("V2 source_eigenvalue: A, B, C must be JSON integers")
    elif "source_eigenvalue" in im:
        bad.append("V1 source_eigenvalue: must be an object")

    sign = im.get("laplacian_sign")
    if "laplacian_sign" in im and sign not in ("nonnegative", "nonpositive"):
        bad.append(f"V2 laplacian_sign: {sign!r} not in the closed set")
    if "radius_normalization" in im and \
            im["radius_normalization"] != "unit_radius_dimensionless":
        bad.append(f"V2 radius_normalization: {im['radius_normalization']!r} "
                   "not in the closed set")
    mc = im.get("multiplicity_convention")
    if isinstance(mc, dict):
        if set(mc) != {"counts", "source_dimension_field"}:
            bad.append("V1 multiplicity_convention: keys exactly "
                       "{counts, source_dimension_field}")
        if mc.get("counts") != "per_protocol_level":
            bad.append("V2 multiplicity_convention.counts: closed set "
                       "{'per_protocol_level'}")
        if mc.get("source_dimension_field") not in ("real", "complex"):
            bad.append("V2 multiplicity_convention.source_dimension_field: "
                       "closed set {'real', 'complex'}")
    elif "multiplicity_convention" in im:
        bad.append("V1 multiplicity_convention: must be an object")
    if "unlisted_source_rows" in im and \
            im["unlisted_source_rows"] != "zero_multiplicity":
        bad.append(f"V2 unlisted_source_rows: {im['unlisted_source_rows']!r} "
                   "not in the closed set")
    off = im.get("off_image_levels")
    if "off_image_levels" in im and off not in ("empty", "spectrum_excludes"):
        bad.append(f"V2 off_image_levels: {off!r} not in the closed set")

    n_max = None
    cb = im.get("certified_band")
    if isinstance(cb, dict):
        if set(cb) != {"n_max"}:
            bad.append("V1 certified_band: keys exactly {n_max}")
        if _is_int(cb.get("n_max")) and cb["n_max"] >= 2:
            n_max = cb["n_max"]
        else:
            bad.append("V2 certified_band.n_max: JSON integer >= 2 required")
    elif "certified_band" in im:
        bad.append("V1 certified_band: must be an object")

    if a is None or b is None or n_max is None or entries is None:
        return bad

    # V4 domain containment over K_band = {k : 0 <= a*k+b <= n_max}.
    #
    # INJECTIVITY IS NOT CHECKED HERE, because it cannot fail: V2 requires k
    # strictly increasing and V3 requires a >= 1, so k1 < k2 implies
    # a*k1 + b < a*k2 + b and two entries can never share a level.  An earlier
    # revision carried an injectivity branch and a contract sentence assigning
    # duplicates to V4; both were removed once the branch was shown unreachable.
    # A duplicate k is a strict-increase violation and is refused by V2.
    mapped, v4_bad = {}, False
    for k, m in entries:
        n = a * k + b
        if not (0 <= n <= n_max):
            bad.append(f"V4 entry k={k}: maps to n={n} outside [0, {n_max}] "
                       "(k outside K_band)")
            v4_bad = True
            continue
        mapped[n] = (k, m)

    # V5 coefficient closure: s*(A,B,C) == (a^2, 2a(b+1), b(b+2)) exactly
    if ABC is not None and sign in ("nonnegative", "nonpositive"):
        s = 1 if sign == "nonnegative" else -1
        got = (s * ABC[0], s * ABC[1], s * ABC[2])
        want = (a * a, 2 * a * (b + 1), b * (b + 2))
        if got != want:
            bad.append(f"V5 coefficient closure: s*(A,B,C) = {got} != "
                       f"(a^2, 2a(b+1), b(b+2)) = {want}")

    if v4_bad:
        return bad

    # V6 zero-level anchor
    if 0 not in mapped:
        bad.append("V6 no entry maps to level 0 (the constants)")
    elif mapped[0][1] != 1:
        bad.append(f"V6 level-0 entry carries m = {mapped[0][1]}, must be exactly 1")

    # V7 band closure over the 12.1.2b filled sequence
    filled = {n: (mapped[n][1] if n in mapped else 0) for n in range(0, n_max + 1)}
    positives = [n for n in range(1, n_max + 1) if filled[n] > 0]
    if len(positives) < 2:
        bad.append(f"V7 fewer than two positive nonzero levels in the filled band "
                   f"(found {positives})")
    elif positives[1] != n_max:
        bad.append(f"V7 second positive nonzero level is {positives[1]}, "
                   f"declared n_max is {n_max}")

    # V8 off-image declaration consistency
    off_class = [n for n in range(0, n_max + 1) if (n - b) % a != 0]
    if off == "empty" and off_class:
        bad.append(f"V8 off_image_levels = 'empty' but the off-image class is {off_class}")
    if off == "spectrum_excludes" and not off_class:
        bad.append("V8 off_image_levels = 'spectrum_excludes' but the off-image "
                   "class is empty")
    return bad


def packet_ii_gate(packet):
    """The complete build-time gate for Packet II: leanness plus V1-V8."""
    return judge(packet, "II") + packet_ii_checks(packet)


def production_seal_refusals(packet):
    """The sealing gate's SYN- refusal (12.1.2): qualification instances are
    mechanically unsealable as production packets."""
    cid = packet.get("case_id", "")
    if isinstance(cid, str) and cid.startswith(SYN_PREFIX):
        return [f"case_id {cid!r} begins with the reserved prefix {SYN_PREFIX!r}: "
                "qualification instances are never sealed as production packets"]
    return []


def filled_reference(packet, transform=None):
    """The 12.1.2b comparison surface: (filled values, provenance classes, n_max).

    With `transform=None` the packet's own validated transform applies and every
    level in [0, n_max] gets exactly one justified value with a provenance class
    from {entry, unlisted, off-image}.  A non-None `transform` is the Q3 fault
    injection: it is applied AT THIS POINT, downstream of validation, the
    corrupted state is deliberately not re-validated, entries landing outside
    the band are dropped, and every level is marked 'injected'.
    """
    im = packet["indexing_map"]
    tf = transform if transform is not None else im["index_transform"]
    a, b = tf["a"], tf["b"]
    n_max = im["certified_band"]["n_max"]
    injected = transform is not None

    landed = {}
    for k, m in packet["reference_values"]:
        n = a * k + b
        if 0 <= n <= n_max:
            landed[n] = m
    values, provenance = {}, {}
    for n in range(0, n_max + 1):
        if n in landed:
            values[n] = landed[n]
            provenance[n] = "injected" if injected else "entry"
        else:
            values[n] = 0
            if injected:
                provenance[n] = "injected"
            elif (n - b) % a != 0:
                provenance[n] = "off-image"
            else:
                provenance[n] = "unlisted"
    return values, provenance, n_max


def _packet_I(q, s, gen_s=None):
    """A Packet-I-shaped dict for exercising the structural gate.

    `gen_s` builds the GENERATOR from different parameters than the ones
    declared, which is how a parameter-to-generator mismatch is constructed
    without disturbing any other predicate: the supplied generator is still a
    free order-q action with the right census, so S1, S2, S4, S5 and S6 all
    stay clean and only S3 can fire.
    """
    u, v = canonical_action_pair(q, gen_s if gen_s is not None else s)
    return {"case_id": "SYN-S3-FIXTURE", "family": "lens",
            "parameters": {"q": q, "s": list(s)},
            "generators": [[[float(x) for x in u], [float(x) for x in v]]],
            "action_convention": "two_sided",
            "format_version": "m8_5b-packet-1"}


def _units(q):
    return [t for t in range(1, q) if math.gcd(t, q) == 1]


def s3_battery(verbose=True):
    """S3 coverage over the FULL admissible parameter class.

    S3 is the contract's parameter-to-generator consistency predicate.  It was
    covered for odd orders only, and skipped in silence for even ones, so a
    conforming-looking even-order packet could be sealed with a generator its
    own declared parameters do not denote.  This battery is the evidence that
    the gap is closed and stays closed:

        - an even-order conforming packet PASSES (homogeneous and inhomogeneous);
        - an even-order wrong generator REDS, and reds FOR S3 alone;
        - an odd-order conforming packet still PASSES;
        - an odd-order wrong generator REDS for S3 alone;
        - across a sweep of every admissible order in range, a wrong generator
          ALWAYS reds, so no parameter class silently skips the predicate;
        - on every odd order the new half-angle construction reproduces the old
          modular one as a 4x4 isometry, so no previously passing packet can
          change verdict.

    Returns (ok, rows).
    """
    rows, ok = [], True
    touched = set()          # every (q, s) this battery exercises, DERIVED

    def record(label, passed, detail=""):
        nonlocal ok
        ok &= passed
        rows.append({"item": label, "pass": bool(passed), "detail": detail})
        if verbose:
            print(f"  {'PASS' if passed else 'FAIL'}  {label}"
                  + (f"  {detail}" if detail and not passed else ""))

    def touch(q, *ss):
        for s in ss:
            touched.add((q, tuple(s)))

    def only_s3(hits):
        return bool(hits) and all(h.startswith("S3 ") for h in hits)

    for label, q, s in (("even order, homogeneous, conforming L(4;1,1)", 4, (1, 1)),
                        ("even order, inhomogeneous, conforming L(12;1,5)", 12, (1, 5))):
        touch(q, s)
        hits = structural_checks(_packet_I(q, s))
        record(label, not hits, "; ".join(hits))

    touch(7, (1, 3))
    hits = structural_checks(_packet_I(7, (1, 3)))
    record("odd order, conforming, still passes L(7;1,3)", not hits, "; ".join(hits))

    for label, q, s, gen_s in (
            ("even order, wrong generator, reds FOR S3 ALONE", 12, (1, 5), (1, 7)),
            ("odd order, wrong generator, reds FOR S3 ALONE", 7, (1, 3), (1, 5))):
        touch(q, s, gen_s)
        a = _act_matrix(*canonical_action_pair(q, s))
        b = _act_matrix(*canonical_action_pair(q, gen_s))
        distinct = np.abs(a - b).max() > 1e-9
        hits = structural_checks(_packet_I(q, s, gen_s=gen_s))
        record(label, distinct and only_s3(hits),
               f"distinct={distinct} hits={hits}")

    # no admissible class may silently skip the predicate
    swept, skipped = [], []
    for q in range(3, 17):
        us = _units(q)
        if len(us) < 2:
            continue
        s, gen_s = (1, us[0]), (1, us[1])
        touch(q, s, gen_s)
        a = _act_matrix(*canonical_action_pair(q, s))
        b = _act_matrix(*canonical_action_pair(q, gen_s))
        if np.abs(a - b).max() <= 1e-9:            # not actually a mutant here
            continue
        swept.append(q)
        if not only_s3(structural_checks(_packet_I(q, s, gen_s=gen_s))):
            skipped.append(q)
    record(f"no-skip sweep: wrong generator reds for S3 at every admissible "
           f"order in {swept[0]}..{swept[-1]} ({len(swept)} orders, "
           f"{sum(1 for q in swept if q % 2 == 0)} even)",
           bool(swept) and not skipped, f"skipped at q={skipped}")

    # the construction is judged against the CONTRACT'S action, not against its
    # own derivation: compare with the complex-coordinate form over the class
    worst, mism = 0.0, []
    for q in range(2, 17):
        for s2 in _units(q):
            touch(q, (1, s2))
            got = _act_matrix(*canonical_action_pair(q, (1, s2)))
            d = float(np.abs(got - lens_action_matrix(q, (1, s2))).max())
            worst = max(worst, d)
            if d > 1e-9:
                mism.append((q, s2))
    record("constructed pair reproduces the lens action computed independently "
           f"in complex coordinates, all orders 2..16 (worst {worst:.1e})",
           not mism, f"mismatches {mism[:4]}")

    # the (u,v) ~ (-u,-v) central lift induces the SAME isometry, so an
    # equivalent lift must be ACCEPTED: S3 judges the induced action, and the
    # half-angle form must not smuggle in a lift convention of its own
    lift_bad = []
    for q, s in ((4, (1, 1)), (12, (1, 5)), (7, (1, 3))):
        touch(q, s)
        u, v = canonical_action_pair(q, s)
        if np.abs(_act_matrix(u, v) - _act_matrix(-u, -v)).max() > 0:
            lift_bad.append((q, s, "induced actions differ"))
        pkt = _packet_I(q, s)
        pkt["generators"] = [[[float(x) for x in -u], [float(x) for x in -v]]]
        if structural_checks(pkt):
            lift_bad.append((q, s, "equivalent central lift rejected"))
    record("central lift [-u,-v] induces an identical isometry and is ACCEPTED "
           "at both parities (no new lift convention introduced)",
           not lift_bad, str(lift_bad))

    # the MIXED lift [u,-v] is a different isometry and must be refused; S3
    # must be among the violations at both parities.  At odd order S4 also
    # fires, because the flipped element generates a group of order 2q; that is
    # two predicates independently noticing a malformed packet, not a defect
    mixed_bad = []
    for q, s in ((12, (1, 5)), (7, (1, 3))):
        touch(q, s)
        u, v = canonical_action_pair(q, s)
        pkt = _packet_I(q, s)
        pkt["generators"] = [[[float(x) for x in u], [float(x) for x in -v]]]
        hits = structural_checks(pkt)
        if not any(h.startswith("S3 ") for h in hits):
            mixed_bad.append((q, s, hits))
    record("mixed lift [u,-v] is refused with S3 among the violations at both "
           "parities", not mixed_bad, str(mixed_bad))

    # backward compatibility: odd orders keep their old verdict exactly
    drift = []
    for q in range(3, 20, 2):
        i2 = pow(2, -1, q)
        for s2 in _units(q):
            s = (1, s2)
            touch(q, s)
            al, be = ((s[0] + s[1]) * i2) % q, ((s[0] - s[1]) * i2) % q
            old = _act_matrix(
                np.array([np.cos(2*np.pi*al/q), np.sin(2*np.pi*al/q), 0., 0.]),
                np.array([np.cos(2*np.pi*be/q), np.sin(2*np.pi*be/q), 0., 0.]))
            if np.abs(old - _act_matrix(*canonical_action_pair(q, s))).max() > 1e-9:
                drift.append((q, s))
    record("odd orders reproduce the previous construction as a 4x4 isometry "
           "(no previously passing packet changes verdict)", not drift,
           f"drift at {drift[:4]}")

    # COVERAGE DIAGNOSTIC, derived from what actually ran rather than declared
    # alongside it: how broad the sweep was.  The set is computed by the
    # battery itself, so it cannot drift from the code that touches it.  It
    # does NOT constrain Stage-1 case eligibility; an earlier draft made it a
    # commissioning constraint and that policy was withdrawn as unnecessary.
    # Eligibility comes from the adopted preregistration alone: non-2I,
    # distinct from the frozen pilot tuning set, informing no pilot choice.
    orders = sorted({q for q, _ in touched})
    evens = sorted({q for q, _ in touched if q % 2 == 0})
    odds = sorted({q for q, _ in touched if q % 2})
    if verbose:
        print(f"\n  COVERAGE DIAGNOSTIC: sweep touched {len(touched)} case(s) "
              f"over orders {orders[0]}..{orders[-1]}")
        print(f"    even orders touched: {evens}")
        print(f"    odd orders touched:  {odds}")
        print("    This sweep is diagnostic only and does not constrain Stage-1")
        print("    case eligibility; the full list is returned by s3_battery()")
        print("    as its third value.")
    return ok, rows, sorted(touched)


def structural_battery(verbose=True):
    """S1-S6 coverage, with each mutant scored against its TARGET predicate.

    The previous suite scored `bool(hits)`: ANY rejection counted, so a mutant
    aimed at S4 that reddened only through S3 still passed, and S5a, S5b and S6
    were never exercised by any fixture at all.  A required contract predicate
    with no falsifiability evidence is exactly the class of defect that the S3
    gap belonged to, so it is closed the same way.

    ISOLATION, and where it is provably unreachable.  Two of these predicates
    cannot fire alone, for reasons that are theorems rather than accidents:

      S5a  For S1-valid parameters the canonical generator is always FREE:
           a fixed point of x -> u x v needs u conjugate to v^-1, which for
           u = e^{i pi (s1+s2)/q}, v = e^{i pi (s1-s2)/q} forces s1 = 0 or
           s2 = 0 mod q, and gcd(s_i, q) = 1 excludes both.  So a fixed point
           implies the supplied generator is not the canonical one, and S3
           fires too.  S5a remains an independent GEOMETRIC backstop, checked
           on the action actually executed rather than inferred from the
           declared parameters (section 8, G1b).
      S6   A group of order q containing an element of order q is cyclic.  So
           whenever S3 passes (the canonical generator, of order q, is present)
           and S4 passes, the census is automatically correct, and S6 can only
           fire alongside S3.

    Where "alone" is unreachable, the item carries a DIFFERENTIAL CONTROL: a
    second packet sharing the co-firing defect but not the target property,
    which must NOT fire the target.  The difference between the two is what
    makes the target attributable rather than incidental.
    """
    rows, ok, exercised = [], True, set()

    def record(label, passed, detail=""):
        nonlocal ok
        ok &= passed
        rows.append({"item": label, "pass": bool(passed), "detail": detail})
        if verbose:
            print(f"  {'PASS' if passed else 'FAIL'}  {label}"
                  + (f"  {detail}" if detail and not passed else ""))

    def tags(hits):
        return {h.split()[0] for h in hits}

    good_I, _ = fixtures()
    u7, v7 = canonical_action_pair(7, (1, 3))

    def pkt(**over):
        p = copy.deepcopy(good_I)
        for k, val in over.items():
            if k == "generators":
                p["generators"] = [[[float(x) for x in a], [float(x) for x in b]]
                                   for a, b in val]
            elif k == "parameters":
                p["parameters"] = val
            else:
                p[k] = val
        return p

    # a wrong-but-otherwise-sound generator: the shared control for S4/S5a/S6
    u_wrong, v_wrong = canonical_action_pair(7, (1, 5))
    ctl_S3_only = pkt(generators=[(u_wrong, v_wrong)])

    I4 = np.array([0., 1, 0, 0])
    J4 = np.array([0., 0, 1, 0])
    ONE = np.array([1., 0, 0, 0])
    conj7 = np.array([u7[0], -u7[1], -u7[2], -u7[3]])
    u83, v83 = canonical_action_pair(8, (1, 3))

    cases = [
        ("S1", "declared action not free, gcd(s,q) != 1",
         pkt(parameters={"q": 7, "s": [7, 3]}), None),
        ("S2", "left-only convention, which cannot express the case",
         pkt(action_convention="left"), None),
        ("S3", "parameters do not reproduce the declared generator",
         pkt(parameters={"q": 7, "s": [1, 5]}), None),
        ("S4", "declared order disagrees with the closure",
         pkt(parameters={"q": 8, "s": [1, 3]}), ctl_S3_only),
        ("S5a", "generator [u, u^-1] fixes a point of S^3 (orbifold, not manifold)",
         pkt(generators=[(u7, conj7)]), ctl_S3_only),
        ("S5b", "both [u,v] and [-u,-v] declared: the deck group is double-counted",
         pkt(generators=[(u7, v7), (-u7, -v7)]), None),
        ("S6", "Q8 (non-cyclic, order 8, free) declared as the cyclic C8",
         pkt(parameters={"q": 8, "s": [1, 1]},
             generators=[(I4, ONE), (J4, ONE)]),
         pkt(parameters={"q": 8, "s": [1, 1]}, generators=[(u83, v83)])),
    ]

    for target, label, mutant, control in cases:
        fired = tags(structural_checks(mutant))
        hit = target in fired
        if control is None:
            passed = hit
            detail = f"fired {sorted(fired)}"
        else:
            ctl_fired = tags(structural_checks(control))
            passed = hit and target not in ctl_fired
            detail = (f"mutant fired {sorted(fired)}, differential control "
                      f"fired {sorted(ctl_fired)}")
        if hit:
            exercised.add(target)
        record(f"{target}: {label}"
               + ("" if control is None else "  [differential control]"),
               passed, detail)

    # `parameters` shape and key exactness (R5).  Two different owners, and the
    # mutants separate them: a key defect must be refused BY key exactness, not
    # incidentally by a later read, and a type defect BY S0.
    def jtags(p):
        return [h for h in judge(p, "I") if h.startswith("parameters")]

    p_extra = pkt(parameters={"q": 7, "s": [1, 3], "note": "x"})
    p_missing = pkt(parameters={"s": [1, 3]})
    p_arity = pkt(parameters={"q": 7, "s": [1, 3, 5]})
    p_type = pkt(parameters={"q": 7, "s": [1, "3"]})

    record("key exactness: an extra `parameters` key is refused by key exactness",
           bool(jtags(p_extra)) and not tags(structural_checks(p_extra)) - {"S3"},
           f"judge {jtags(p_extra)}")
    record("key exactness: a missing `parameters` key is refused by key exactness, "
           "not incidentally downstream", bool(jtags(p_missing)),
           f"judge {jtags(p_missing)}")
    record("S0: `s` of wrong arity is refused by the shape predicate",
           "S0" in tags(structural_checks(p_arity)) and not jtags(p_arity),
           f"S {tags(structural_checks(p_arity))}, judge {jtags(p_arity)}")
    record("S0: a non-integer entry in `s` is refused by the shape predicate",
           "S0" in tags(structural_checks(p_type)) and not jtags(p_type),
           f"S {tags(structural_checks(p_type))}")
    if "S0" in tags(structural_checks(p_arity)):
        exercised.add("S0")

    hits = structural_checks(good_I)
    record("conforming packet passes every structural predicate clean",
           not hits, "; ".join(hits))
    record("conforming packet passes key exactness clean", not judge(good_I, "I"),
           "; ".join(judge(good_I, "I")))

    missing = {"S0", "S1", "S2", "S3", "S4", "S5a", "S5b", "S6"} - exercised
    record(f"coverage: every structural predicate exercised by a mutant that "
           f"fires IT ({len(exercised)}/8)", not missing, f"never fired: {sorted(missing)}")

    # supporting evidence for S6: the census genuinely discriminates groups of
    # the SAME order, which is what makes it more than an order check
    def census(mats):
        els = _close(mats)
        c = {}
        for A in els:
            c[_order(A)] = c.get(_order(A), 0) + 1
        return len(els), dict(sorted(c.items()))
    nA, cA = census([_act_matrix(I4, ONE)])
    nB, cB = census([_act_matrix(I4, I4), _act_matrix(J4, J4)])
    record(f"census discriminates same-order groups: order {nA} both, "
           f"C4 {cA} vs Klein {cB}", nA == nB and cA != cB)

    return ok, rows


def fixtures():
    """Synthetic packets. Dummy values throughout; no real case data here."""
    # A NON-CASE synthetic lens example, q = 7, used only to exercise the gate.
    # This is not the adjudication case and is not derived from it.
    q, s1, s2 = 7, 1, 3
    u, v = canonical_action_pair(q, (s1, s2))       # one construction, one path
    good_I = {"case_id": "CASE-X", "family": "lens",
              "parameters": {"q": q, "s": [s1, s2]},
              "generators": [[[float(x) for x in u], [float(x) for x in v]]],
              "action_convention": "two_sided",
              "format_version": "m8_5b-packet-1"}
    # Conforming under addendum 12.1: identity transform, zeros listed,
    # off-image class empty.  SYN- case_id: unsealable as production by design.
    good_II = {"case_id": "SYN-SCHEMA-X",
               "citation": {"authors": "SYN FIXTURE", "title": "SCHEMA SUITE",
                            "venue": "NONE", "year": 2026, "doi": "10.0000/syn",
                            "table": "T0", "row": "R0"},
               "indexing_map": {
                   "index_transform": {"kind": "affine", "a": 1, "b": 0},
                   "source_eigenvalue": {"form": "quadratic", "A": 1, "B": 2, "C": 0},
                   "laplacian_sign": "nonnegative",
                   "radius_normalization": "unit_radius_dimensionless",
                   "multiplicity_convention": {"counts": "per_protocol_level",
                                               "source_dimension_field": "real"},
                   "unlisted_source_rows": "zero_multiplicity",
                   "off_image_levels": "empty",
                   "certified_band": {"n_max": 3}},
               "reference_values": [[0, 1], [1, 0], [2, 3], [3, 8]],
               "format_version": "m8_5b-packet-II-2"}
    return good_I, good_II


def family_controls():
    """Two more conforming instances proving the validator accepts the whole
    admitted family, not only the identity corner: a stride-2 source with a
    nonempty off-image class, and a shifted source stated with a nonpositive
    Laplacian.  Values are pilot-anchored (L(2,1) and L(3,1) tuning tables)."""
    stride = {"case_id": "SYN-STRIDE-X",
              "citation": {"authors": "SYN FIXTURE", "title": "SCHEMA SUITE",
                           "venue": "NONE", "year": 2026, "doi": "10.0000/syn",
                           "table": "T1", "row": "R1"},
              "indexing_map": {
                  "index_transform": {"kind": "affine", "a": 2, "b": 0},
                  "source_eigenvalue": {"form": "quadratic", "A": 4, "B": 4, "C": 0},
                  "laplacian_sign": "nonnegative",
                  "radius_normalization": "unit_radius_dimensionless",
                  "multiplicity_convention": {"counts": "per_protocol_level",
                                              "source_dimension_field": "real"},
                  "unlisted_source_rows": "zero_multiplicity",
                  "off_image_levels": "spectrum_excludes",
                  "certified_band": {"n_max": 4}},
              "reference_values": [[0, 1], [1, 9], [2, 25]],
              "format_version": "m8_5b-packet-II-2"}
    shifted = {"case_id": "SYN-SHIFT-X",
               "citation": {"authors": "SYN FIXTURE", "title": "SCHEMA SUITE",
                            "venue": "NONE", "year": 2026, "doi": "10.0000/syn",
                            "table": "T2", "row": "R2"},
               "indexing_map": {
                   "index_transform": {"kind": "affine", "a": 1, "b": -2},
                   "source_eigenvalue": {"form": "quadratic", "A": -1, "B": 2, "C": 0},
                   "laplacian_sign": "nonpositive",
                   "radius_normalization": "unit_radius_dimensionless",
                   "multiplicity_convention": {"counts": "per_protocol_level",
                                               "source_dimension_field": "real"},
                   "unlisted_source_rows": "zero_multiplicity",
                   "off_image_levels": "empty",
                   "certified_band": {"n_max": 3}},
               "reference_values": [[2, 1], [4, 3], [5, 8]],
               "format_version": "m8_5b-packet-II-2"}
    return stride, shifted


def main():
    good_I, good_II = fixtures()
    stride, shifted = family_controls()
    muts = [
        ("I", {**good_I, "notes": "kept for the adjudicator"}, "extra notes field"),
        ("I", {**good_I, "provenance": "round 3 rework"}, "extra provenance field"),
        ("II", {**good_II, "adjudication_notes": "see the log"}, "extra notes field"),
        ("II", {**good_II, "citation": {**good_II["citation"], "comment": "x"}},
         "extra citation field"),
        ("II", {**good_II, "case_id": "CASE-X. However, we corrected this later"},
         "prose smuggled into an allowed field"),
        ("II", {**good_II, "indexing_map": {**good_II["indexing_map"],
                                            "note": "n " + "x" * 250}},
         "overlong string in an allowed field"),
        ("I", {k: v for k, v in good_I.items() if k != "generators"},
         "missing required field"),
    ]
    print("structural suite S1-S6 (each mutant scored against its TARGET):")
    sok, _ = structural_battery()

    print("\nS3 full-class coverage (parameter-to-generator consistency):")
    s3_ok, _, _ = s3_battery()
    sok &= s3_ok

    if not sok:
        sys.exit("\nSTRUCTURAL GATE NOT TRUSTWORTHY. Do not seal.")
    print()

    print("leanness suite (each MUST be rejected):")
    ok = True
    for which, pkt, name in muts:
        hits = judge(pkt, which)
        print(f"  {'REJECTED' if hits else 'ACCEPTED'}  {name}")
        ok &= bool(hits)

    # ---- addendum 12.1 V-suite: every predicate fires on a single-defect
    # ---- mutant, and the whole admitted family is accepted clean
    print("\nV-suite (addendum 12.1; the TARGET predicate must fire):")

    def vmut(fn):
        m = copy.deepcopy(good_II)
        fn(m)
        return m
    vmuts = [
        ("V1", "unknown key inside indexing_map",
         vmut(lambda p: p["indexing_map"].__setitem__("extra_convention", "x"))),
        ("V1", "missing laplacian_sign",
         vmut(lambda p: p["indexing_map"].pop("laplacian_sign"))),
        ("V2", "float where integer required (a = 1.0)",
         vmut(lambda p: p["indexing_map"]["index_transform"].__setitem__("a", 1.0))),
        ("V2", "null smuggled into citation.year",
         vmut(lambda p: p["citation"].__setitem__("year", None))),
        ("V2", "negative multiplicity",
         vmut(lambda p: p["reference_values"].__setitem__(2, [2, -3]))),
        ("V2", "duplicate k: a strict-increase violation, V2's not V4's",
         vmut(lambda p: p["reference_values"].__setitem__(1, [0, 0]))),
        ("V2", "k distinct but unsorted: strict increase is V2's",
         vmut(lambda p: p.__setitem__("reference_values",
              [p["reference_values"][0], p["reference_values"][2],
               p["reference_values"][1], p["reference_values"][3]]))),
        ("V2", "wrong format_version",
         vmut(lambda p: p.__setitem__("format_version", "m8_5b-packet-1"))),
        ("V3", "transform kind outside the closed set",
         vmut(lambda p: p["indexing_map"]["index_transform"].update(kind="permutation"))),
        ("V3", "a = 0",
         vmut(lambda p: p["indexing_map"]["index_transform"].update(a=0))),
        ("V4", "entry beyond the certified band",
         vmut(lambda p: p["reference_values"].append([9, 7]))),
        ("V5", "coefficient closure broken (B = 3)",
         vmut(lambda p: p["indexing_map"]["source_eigenvalue"].update(B=3))),
        ("V6", "level-0 entry not the constants (m = 2)",
         vmut(lambda p: p["reference_values"].__setitem__(0, [0, 2]))),
        ("V7", "second positive nonzero level disagrees with declared n_max",
         vmut(lambda p: p["reference_values"].__setitem__(1, [1, 5]))),
        ("V7", "fewer than two positive nonzero levels",
         vmut(lambda p: p["reference_values"].__setitem__(2, [2, 0]))),
        ("V8", "off-image affirmation while the class is empty",
         vmut(lambda p: p["indexing_map"].update(off_image_levels="spectrum_excludes"))),
    ]
    covered = set()
    for target, name, pkt in vmuts:
        hits = packet_ii_checks(pkt)
        fired = any(h.startswith(target + " ") for h in hits)
        covered.add(target) if fired else None
        print(f"  {'RED ' + target if fired else 'MISSED  '}  {name}")
        ok &= fired
    uncovered = {"V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"} - covered
    if uncovered:
        print(f"  UNCOVERED predicates: {sorted(uncovered)}")
        ok = False
    # V8's other direction is exercised by the stride control below: it declares
    # spectrum_excludes and MUST be accepted, so 'empty'-while-nonempty is the
    # mutant and 'spectrum_excludes'-while-nonempty is the control.
    v8b = copy.deepcopy(stride)
    v8b["indexing_map"]["off_image_levels"] = "empty"
    hits = packet_ii_checks(v8b)
    fired = any(h.startswith("V8 ") for h in hits)
    print(f"  {'RED V8' if fired else 'MISSED  '}  'empty' declared while the "
          "off-image class is nonempty")
    ok &= fired

    # totality over malformed top-level values: canonical JSON that is not an
    # object must REJECT at V1/R1, never crash (Redline audit blocker 1)
    print("  totality over non-object top-level values:")
    for bad_top in ([], None, 3, "x"):
        try:
            hits = packet_ii_gate(bad_top)
            fired, crashed = any(h.startswith("V1 ") or ": top-level" in h
                                 for h in hits), False
        except Exception as exc:
            fired, crashed = False, True
            print(f"            CRASHED: {exc!r}")
        print(f"  {'RED V1' if fired and not crashed else 'MISSED  '}  "
              f"Packet II top-level {bad_top!r}")
        ok &= fired and not crashed
    try:
        hits_i = judge([], "I")
        fired, crashed = bool(hits_i), False
    except Exception as exc:
        fired, crashed = False, True
        print(f"            CRASHED: {exc!r}")
    print(f"  {'RED     ' if fired and not crashed else 'MISSED  '}  "
          "Packet I top-level [] rejected by judge without crashing")
    ok &= fired and not crashed

    print("\nnegative controls (each MUST pass clean):")
    for name, pkt in (("packet I", None), ("identity-corner packet II", good_II),
                      ("stride-2 / off-image packet II", stride),
                      ("shifted / nonpositive-sign packet II", shifted)):
        if pkt is None:
            hits = judge(good_I, "I")
        else:
            hits = packet_ii_gate(pkt)
        print(f"  {'OK      ' if not hits else 'FALSE+  '}  {name}")
        for h in hits:
            print(f"            {h}")
        ok &= not hits

    print("\nsealing gate (SYN- reservation):")
    refusals = production_seal_refusals(good_II)
    print(f"  {'REFUSED ' if refusals else 'MISSED  '}  SYN- case_id at the "
          "production sealing gate")
    ok &= bool(refusals)
    prod_like = copy.deepcopy(good_II)
    prod_like["case_id"] = "QX-CTRL-1"
    refusals = production_seal_refusals(prod_like)
    print(f"  {'OK      ' if not refusals else 'FALSE+  '}  non-SYN case_id "
          "passes the reservation check")
    ok &= not refusals

    if not ok:
        sys.exit("\nGATE NOT TRUSTWORTHY: suite misbehaved. Do not seal.")
    print("\nsuite behaves: the gate can fail.")

    # Judge packets ONLY when paths are passed explicitly (as 'I:path' or
    # 'II:path').  The automatic packets/ scan is deliberately gone: the burned
    # packets are a closed record judged by their own era's gate, and this unit
    # stays clean of Packet II's contents (HANDOFF: kept clean on purpose).
    for arg in sys.argv[1:]:
        which, _, path = arg.partition(":")
        if which not in ("I", "II") or not path:
            print(f"\nignored argument {arg!r}: expected I:path or II:path")
            continue
        pkt = json.loads(open(path).read())
        if which == "I":
            hits = judge(pkt, "I") + structural_checks(pkt)
        else:
            hits = packet_ii_gate(pkt) + production_seal_refusals(pkt)
        print(f"\n{path}: {'CLEAR TO SEAL' if not hits else 'NOT CLEAR'}")
        for h in hits:
            print(f"  {h}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
