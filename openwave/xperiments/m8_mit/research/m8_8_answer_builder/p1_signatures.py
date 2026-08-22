"""P1 step 2: the nine § 5.5 row signatures, computed fresh.

Target-free. Signatures are character values on three group elements; they carry
no torsion value.

Two things are DERIVED here rather than assumed, because both have already bitten
this build once:

  1. Irreducibles are selected by MEASUREMENT (norm 1 under the group inner
     product, then deduplicated), never by assuming which Sym^n stays irreducible.
     Assuming it once produced a duplicate ninth "irrep": Sym^3's Galois conjugate
     is a different irrep, but Sym^4's is itself.
  2. The R-label assignment is recovered as the UNIQUE isomorphism from the
     measured McKay graph onto the declared affine-E8 edge list, and uniqueness is
     checked by exhaustion rather than asserted. Labels decide which torsion value
     attaches to which signature, so a silent swap of R1 and R2 would mispair the
     packet while every count still looked right.

The generator identities come from the construction packet's `abstract_generators`
and are element IDs in the § 4.2 enumeration, so p1_enumerate must agree first.
"""

import ast
import itertools
import json
import pathlib
import re
import sys

from qphi_exact import Phi
from p1_enumerate import close_group, enumerate_group


def parse_mckay_declaration(path):
    """Pull `labels` and `edges` out of a source file without importing it."""
    txt = pathlib.Path(path).read_text()
    out = {}
    for name in ("labels", "edges"):
        m = re.search(rf"^{name}\s*=\s*(\[.*?\])\s*$", txt, re.S | re.M)
        if not m:
            raise ValueError(f"{name} not found in {path}")
        out[name] = ast.literal_eval(m.group(1))
    return out["labels"], [tuple(e) for e in out["edges"]]


def character_v1(elems):
    """chi of the defining 2-dimensional rep: twice the real part.

    Components are (A + B*phi)/2, so 2*Re is exactly A + B*phi, integral."""
    return [Phi(2 * g.c[0].p, 2 * g.c[0].q) for g in elems]


def inner(chi, psi, inv_index):
    """(1/|G|) sum_g chi(g) psi(g^-1). Exact; returns a Phi."""
    acc = Phi(0)
    for i, x in enumerate(chi):
        acc = acc + x * psi[inv_index[i]]
    from fractions import Fraction as F
    return Phi(acc.p / F(len(chi)), acc.q / F(len(chi)))


def is_one(v):
    return v.p == 1 and v.q == 0


def build_irreps(elems):
    """All irreducible characters, obtained by measurement rather than assumption.

    Seed with the Sym^n tower by McKay recursion and its Galois twists, keeping
    whatever measures irreducible. That alone is NOT enough and the shortfall is
    the instructive part: conj(Sym^3) equals Sym^3, so the Sym tower plus twists
    yields only eight of the nine irreps and misses one of the two dimension-4
    reps entirely. The count is therefore closed by tensoring known irreducibles
    and peeling off components already known; any residual of norm 1 is new.

    Completeness is asserted against sum(dim^2) = |G|, which is what makes the
    shortfall loud instead of silent."""
    n = len(elems)
    pos = {g: i for i, g in enumerate(elems)}
    inv_index = [pos[g.conj()] for g in elems]  # unit quaternions: g^-1 = conj g

    chi1 = character_v1(elems)
    sym = [[Phi(1)] * n, chi1]                       # Sym^0, Sym^1
    while len(sym) < 12:
        prev, cur = sym[-2], sym[-1]
        sym.append([chi1[i] * cur[i] - prev[i] for i in range(n)])

    irreps, seen = [], []

    def offer(name, chi):
        if not is_one(inner(chi, chi, inv_index)):
            return False
        key = tuple((x.p, x.q) for x in chi)
        if key in seen:
            return False
        seen.append(key)
        irreps.append((name, chi))
        return True

    for k, chi in enumerate(sym):
        if offer(f"Sym^{k}", chi):
            offer(f"conj(Sym^{k})", [x.conj() for x in chi])

    # Tensor closure. Bounded by |G|, and the bound is the completeness check.
    ident_i = next(i for i, g in enumerate(elems)
                   if (g.c[0].p, g.c[0].q) == (1, 0)
                   and all((c.p, c.q) == (0, 0) for c in g.c[1:]))

    def total():
        return sum(int(chi[ident_i].p) ** 2 for _, chi in irreps)

    grew = True
    while grew and total() < n:
        grew = False
        for (na, ca), (nb, cb) in itertools.product(list(irreps), repeat=2):
            prod = [ca[i] * cb[i] for i in range(n)]
            resid = prod
            for _, ck in list(irreps):
                m = inner(prod, ck, inv_index)
                if m.q != 0 or m.p.denominator != 1:
                    raise ValueError("a tensor multiplicity is not an integer")
                if m.p:
                    resid = [resid[i] - Phi(m.p) * ck[i] for i in range(n)]
            if any((x.p, x.q) != (0, 0) for x in resid):
                if offer(f"{na}(x){nb} residual", resid):
                    offer(f"conj({na}(x){nb} residual)", [x.conj() for x in resid])
                    grew = True
                    break

    if total() != n:
        raise ValueError(f"incomplete: sum(dim^2) = {total()}, expected {n}")
    return irreps, inv_index


def mckay_graph(irreps, chi1, inv_index):
    """Undirected edges {a, b} where b occurs in a (x) V1."""
    edges = set()
    for a, (_, ca) in enumerate(irreps):
        prod = [chi1[i] * ca[i] for i in range(len(ca))]
        for b, (_, cb) in enumerate(irreps):
            m = inner(prod, cb, inv_index)
            if m.q != 0 or m.p < 0 or m.p.denominator != 1:
                raise ValueError("tensor multiplicity is not a non-negative integer")
            if m.p > 0:
                edges.add(frozenset((a, b)))
    return {e for e in edges if len(e) == 2}


def assign_labels(n, edges, decl_labels, decl_edges):
    """Every graph isomorphism from the declared diagram onto the measured one.

    Exhaustive over all label assignments, so the count returned IS the
    uniqueness evidence. The declared diagram carries no dimensions, so nothing
    here is seeded with the affine marks; the graph alone must pin the labels,
    and if it does not, the caller sees more than one solution and stops."""
    sols = []
    for perm in itertools.permutations(range(n)):
        m = dict(zip(decl_labels, perm))
        if {frozenset((m[a], m[b])) for a, b in decl_edges} == edges:
            sols.append(m)
    return sols


def prepare(group_packet_path, decl_source_paths):
    """The expensive, generator-independent part: enumeration, irreps, labels.

    Split out from run() so an alternate generator can be evaluated without
    re-deriving any of it. That matters for the wrong-class diagnostic, which
    has to sweep every admissible generator rather than assume one."""
    gp = json.loads(pathlib.Path(group_packet_path).read_text())
    elems = enumerate_group(gp)
    pos = {g: i for i, g in enumerate(elems)}

    decls = [parse_mckay_declaration(p) for p in decl_source_paths]
    if len({(tuple(l), tuple(e)) for l, e in decls}) != 1:
        raise ValueError("the declared McKay diagrams disagree between sources")
    decl_labels, decl_edges = decls[0]

    irreps, inv_index = build_irreps(elems)
    chi1 = character_v1(elems)

    ident_i = next(i for i, g in enumerate(elems)
                   if (g.c[0].p, g.c[0].q) == (1, 0)
                   and all((c.p, c.q) == (0, 0) for c in g.c[1:]))
    dims = []
    for _, chi in irreps:
        d = chi[ident_i]
        if d.q != 0 or d.p.denominator != 1 or d.p <= 0:
            raise ValueError("a dimension is not a positive integer")
        dims.append(int(d.p))

    edges = mckay_graph(irreps, chi1, inv_index)
    sols = assign_labels(len(irreps), edges, decl_labels, decl_edges)

    return {"elems": elems, "pos": pos, "irreps": irreps, "dims": dims,
            "inv_index": inv_index, "ident_i": ident_i, "n_label_solutions": len(sols),
            "label_map": sols[0] if len(sols) == 1 else None}


def order_of(prep, i):
    elems, ident = prep["elems"], prep["elems"][prep["ident_i"]]
    x, k = elems[i], 1
    while k <= len(elems):
        if x == ident:
            return k
        x, k = x * elems[i], k + 1
    raise ValueError("element has no finite order")


def signatures_for(prep, s_id, t_id):
    """The nine § 5.5 signatures under a given generator pair."""
    if prep["label_map"] is None:
        return {}, None
    pos, elems, irreps, dims = prep["pos"], prep["elems"], prep["irreps"], prep["dims"]
    st_id = pos[elems[s_id] * elems[t_id]]
    rows = {}
    for lab, i in prep["label_map"].items():
        chi = irreps[i][1]
        rows[lab] = {"built_as": irreps[i][0], "dimension": dims[i],
                     "s": chi[s_id].triple(), "t": chi[t_id].triple(),
                     "st": chi[st_id].triple()}
    return rows, st_id


def sig_tuple(r):
    return (r["dimension"], r["s"], r["t"], r["st"])


def admissible_generators(prep, s_id):
    """Every t' that satisfies the presentation with this s and generates 2I.

    DERIVED, not read off a list. The relations of <s,t | s^3 = t^5 = (st)^2>
    are checked directly and generation is checked by closure, so the count of
    admissible choices is a measurement."""
    elems, pos = prep["elems"], prep["pos"]
    s = elems[s_id]
    s3 = s * s * s
    out = []
    for j, t in enumerate(elems):
        t5 = t
        for _ in range(4):
            t5 = t5 * t
        st = s * t
        if t5 == s3 and st * st == s3 and len(close_group([s, t])) == len(elems):
            out.append(j)
    return out


def conjugacy_classes(prep, ids):
    """Partition the given element IDs by actual conjugacy in the group.

    Conjugacy, not equality of character columns: the columns are a consequence,
    and using them as the definition would assume the thing being measured."""
    elems, pos = prep["elems"], prep["pos"]
    remaining, classes = list(ids), []
    while remaining:
        a = remaining[0]
        orbit = {pos[g * elems[a] * g.conj()] for g in elems}
        cls = sorted(i for i in remaining if i in orbit)
        classes.append(cls)
        remaining = [i for i in remaining if i not in orbit]
    return classes


def run(group_packet_path, construction_packet_path, decl_source_paths, prep=None):
    cp = json.loads(pathlib.Path(construction_packet_path).read_text())
    if prep is None:
        prep = prepare(group_packet_path, decl_source_paths)

    s_id, t_id = cp["abstract_generators"]["s"], cp["abstract_generators"]["t"]
    rows, st_id = signatures_for(prep, s_id, t_id)

    # § 5.5 makes separation a freeze-review obligation, so it is measured here
    # rather than inherited from the two earlier times it was checked.
    sigs = [sig_tuple(v) for v in rows.values()]
    order = sorted(rows)
    collisions = sorted(
        (a, b) for i, a in enumerate(order) for b in order[i + 1:]
        if sig_tuple(rows[a]) == sig_tuple(rows[b]))

    dims = prep["dims"]
    # When the label assignment is not unique there are no rows and no st, and
    # the run must still REJECT rather than raise: the caller's gate on
    # n_label_solutions is the thing that should speak.
    orders = {"s": order_of(prep, s_id), "t": order_of(prep, t_id),
              "st": order_of(prep, st_id) if st_id is not None else None}
    return {
        "n_irreps": len(prep["irreps"]),
        "dims": sorted(dims),
        "sum_of_squares": sum(d * d for d in dims),
        "n_label_solutions": prep["n_label_solutions"],
        "s_id": s_id, "t_id": t_id, "st_id": st_id,
        "orders": orders,
        "signatures_distinct": len(set(sigs)) == len(sigs) == 9,
        "signature_collisions": collisions,
        "rows": rows,
        "_prep": prep,
    }


def generator_class_diagnostic(prep, s_id, t_id):
    """Sweep every admissible generator and DERIVE what a wrong choice does.

    Nothing here is encoded: the candidate set, the class partition, the induced
    row correspondence and the resulting permutation are all measured. The
    caller gates the sensitivity pattern against this."""
    cands = admissible_generators(prep, s_id)
    classes = conjugacy_classes(prep, cands)
    base_rows, _ = signatures_for(prep, s_id, t_id)
    base_by_sig = {sig_tuple(v): k for k, v in base_rows.items()}
    correct = next((c for c in classes if t_id in c), None)

    # EVERY candidate is evaluated, not one representative per class. Conjugacy
    # predicts they agree within a class; that prediction is checked rather than
    # relied on, since it is the whole basis of the finding.
    per_cand = {}
    for j in cands:
        rows, _ = signatures_for(prep, s_id, j)
        same_set = {sig_tuple(v) for v in rows.values()} == set(base_by_sig)
        perm = ({lab: base_by_sig[sig_tuple(v)] for lab, v in rows.items()}
                if same_set else None)
        per_cand[j] = {"identical_rows": rows == base_rows,
                       "same_signature_set": same_set, "permutation": perm}

    out = []
    for cls in classes:
        perms = {json.dumps(per_cand[j]["permutation"], sort_keys=True) for j in cls}
        first = per_cand[cls[0]]
        out.append({
            "members": cls,
            "is_packet_class": cls is correct,
            "uniform_within_class": len(perms) == 1 and len(
                {per_cand[j]["identical_rows"] for j in cls}) == 1,
            "identical_rows": all(per_cand[j]["identical_rows"] for j in cls),
            "same_signature_set": all(per_cand[j]["same_signature_set"] for j in cls),
            "permutation": first["permutation"],
            "moved": sorted(k for k, v in (first["permutation"] or {}).items() if k != v),
        })
    return {"candidates": cands, "classes": out, "per_candidate": per_cand}
