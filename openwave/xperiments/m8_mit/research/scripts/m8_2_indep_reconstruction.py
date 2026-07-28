"""M8.2 REVIEW ARTIFACT: independent reconstruction of the first-occurrence table.

WHAT THIS IS.  The verification a maintainer ran while reviewing PR #350, landed as a
repository artifact at the author's request (PR #350 close-out, 2026-07-28) so that the
M8.5 reproduction protocol can quarantine it by name.  It is the SECOND of the three
objects that protocol distinguishes:

  1  m8_2_first_occurrence.py       the M8.2 generator (the author's method)
  2  m8_2_indep_reconstruction.py   THIS FILE, the review-time independent check
  3  the M8.5 implementation        written in a fresh context; forbidden to read (1)
                                    or (2) until its own source and raw output are
                                    committed, then may use both as adjudication
                                    references only

Agreement of all three is reported as THREE-WAY AGREEMENT.  It does not raise the claim
label: this file's author had already seen the tables, so the ceiling stays
independent-method reproduction, never blind (roadmap section CONVENTIONS).

INDEPENDENCE FROM (1).  No shared method:

  (1) builds the 2I character table from the McKay / affine-E8 recursion
      V_{n+1} = V_1 (x) V_n - V_{n-1} over a hardcoded 9-class list, and carries the
      irrep labels, the dims and the McKay distances as literals.
  (2) builds the 120 elements of 2I explicitly as unit quaternions, finds the conjugacy
      classes by brute-force conjugation, extracts the irreducible characters by Burnside
      class-sum diagonalization, and DERIVES the dims, the McKay adjacency and the graph
      distances.  Rows are matched label-free, by the (dim, distance) signature, because
      the irrep names are an input to (1) and an output here.

This file imports nothing from (1) and reads none of its fixtures.  Two values do cross
over and are declared rather than derived: the published table itself (DOC below,
transcribed from the pre-registration section 6.1), which is the object under comparison,
and the three flat connections tau = Sym^2(sigma) for sigma in {trivial, Q, Q'}, which are
a contract choice from the pre-registration section 2, not a computed result.

WHAT THIS VERIFIES.  Eight checks, each mutation-tested (`--mutation-tests`), i.e. each is
shown to go red under a deliberate defect, so no PASS line here is a line that cannot fail:

  C1  the 120 quaternions form a group (order 120, closed under multiplication)
  C2  9 conjugacy classes, sizes dividing 120 and summing to 120
  C3  the derived character table is orthonormal
  C4  the derived 2-dim rep identified as Q has chi(g) = 2 cos(theta_g)
  C5  Sym^2 Q and Sym^2 Q' are irreducible, both of dim 3, and distinct (the Galois pair)
  C6  the derived dims satisfy the affine-E8 mark condition A.dims = 2 dims
  C7  McKay distances from the trivial node are well defined, and the 9 (dim, distance)
      signatures are distinct (the precondition for label-free row matching)
  C8  all 9 rows of the section 6.1 0-form table, reproduced and compared entry by entry

WHAT THIS DOES NOT VERIFY:

  - the COEXACT 1-FORM entry rule of the M8.2 script (level 2 at d=0, 3 at d=1, else d).
    Nothing here derives it, and it has no published target, so it remains ASSERTED.
    This file reconstructs only the scalar (0-form) first-occurrence table.
  - the physical interpretation: which object each table is a fixture for (M4_int, the
    prospective M7_ad) is a contract question settled in the pre-registration, not here.
  - the author's mass-spectrum.md section 4 target table.  The M8.2 script compares against
    it; this file compares against the pre-registration section 6.1 table.  The two agreeing
    is a consequence of the M8.2 script's own cross-check, not an independent result here.
  - anything about the arena's dynamics, spectrum on the quotient, or the M8.4 Lagrangians.

Method note: ../findings/m8_2_indep_reconstruction_note.md
Raw output:  ../data/m8_2_indep_reconstruction_raw.txt  (+ .json for the structured form)

Run: python3 m8_2_indep_reconstruction.py [--json PATH] [--mutation-tests]
"""

import argparse
import functools
import itertools
import json
import pathlib
import platform
import subprocess
import sys

import numpy as np

PHI = (1.0 + 5.0**0.5) / 2.0
ORDER = 120
TOL = 1e-9

# The object under comparison: the 0-form first-occurrence table of the pre-registration
# section 6.1, keyed label-free by (irrep dim, McKay distance) -> level n for the three
# connections (trivial, standard, galois).  Transcribed by hand; nothing reads it back.
DOC = {
    (1, 0): (0, 2, 6),
    (2, 1): (1, 1, 5),
    (3, 2): (2, 0, 4),
    (4, 3): (3, 1, 3),
    (5, 4): (4, 2, 2),
    (6, 5): (5, 3, 1),
    (3, 6): (6, 4, 0),
    (4, 6): (6, 4, 2),
    (2, 7): (7, 5, 3),
}

CHECKS = {
    "C1": "120 quaternions form a group (order, closure)",
    "C2": "9 conjugacy classes, sizes divide and sum to 120",
    "C3": "derived character table orthonormal",
    "C4": "the rep identified as Q has chi = 2 cos(theta)",
    "C5": "Sym^2 Q, Sym^2 Q' irreducible dim 3 and distinct",
    "C6": "derived dims satisfy A.dims = 2 dims (E8 marks)",
    "C7": "McKay distances defined, (dim, dist) signatures distinct",
    "C8": "all 9 rows of section 6.1 reproduced",
}

# Each mutation must turn its target checks red.  Cross-effects on other checks are
# reported but not required.
MUTATIONS = {
    "phi_wrong": ("C1", "golden ratio replaced by 1.6, so the icosians leave the group"),
    "conj_without_inverse": ("C2", "conjugation g x g instead of g x g^-1"),
    "perturb_char": ("C3", "one derived character row scaled by 1.05"),
    "q_is_qprime": ("C4", "the Galois partner Q' picked as the defining rep"),
    "sym2_as_square": ("C5", "tau built as chi^2 instead of Sym^2(chi)"),
    "edge_threshold": ("C6,C7", "McKay edge accepted only at multiplicity > 1.5"),
    "bfs_forget_increment": ("C7", "BFS assigns the parent's distance, not parent + 1"),
    "chiv_offbyone": ("C8", "SU(2) character summed over n weights instead of n + 1"),
    "doc_typo": ("C8", "one transcribed target entry altered"),
}


class Halt(Exception):
    """A structural precondition failed; later checks cannot be evaluated."""


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])


def qkey(q):
    return tuple(round(c, 7) + 0.0 for c in q)


def build_2i(phi):
    """The 120 icosians: 8 units + 16 Hurwitz + 96 even permutations of the golden set."""
    elems = []
    for i in range(4):  # +-1, +-i, +-j, +-k
        for s in (1.0, -1.0):
            v = [0.0] * 4
            v[i] = s
            elems.append(tuple(v))
    for signs in itertools.product((0.5, -0.5), repeat=4):
        elems.append(tuple(signs))
    even = [
        (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2), (1, 0, 3, 2),
        (1, 2, 0, 3), (1, 3, 2, 0), (2, 0, 1, 3), (2, 1, 3, 0),
        (2, 3, 0, 1), (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0),
    ]  # fmt: skip
    base = [0.0, 0.5, 0.5 / phi, 0.5 * phi]
    for perm in even:
        for sg in itertools.product((1.0, -1.0), repeat=3):
            v = [0.0] * 4
            vals = [base[0], sg[0] * base[1], sg[1] * base[2], sg[2] * base[3]]
            for slot, val in zip(perm, vals):
                v[slot] = val
            elems.append(tuple(v))
    uniq = {}
    for q in elems:
        uniq[qkey(q)] = q
    return list(uniq.values())


@functools.cache
def provenance():
    """Environment and repository state, so a later reader can place this run exactly."""
    try:
        sha = subprocess.run(
            ["git", "-C", str(pathlib.Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        sha = "unknown"
    return {
        "commit": sha,
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "seed": 11,
    }


def chi_su2(n, thetas, mut=None):
    """SU(2) character of the spin-n/2 rep, as the explicit sum over its n+1 weights."""
    top = n if mut == "chiv_offbyone" else n + 1
    return np.array(
        [sum(np.cos((n - 2 * k) * t) for k in range(top)) for t in np.atleast_1d(thetas)]
    )


def reconstruct(mut=None):
    """Rebuild the table from scratch.  Returns {checks, values, table, comparison}."""
    checks, values = {}, {}
    table, comparison = {}, []

    def finish():
        return {
            "provenance": provenance(),
            "checks": checks,
            "values": values,
            "table": {f"dim{k[0]}_d{k[1]}": list(v) for k, v in table.items()},
            "comparison": comparison,
        }

    # ---- C1: the group ----------------------------------------------------------
    g_elems = build_2i(1.6 if mut == "phi_wrong" else PHI)
    kidx = {qkey(q): n for n, q in enumerate(g_elems)}
    closed = all(qkey(qmul(a, b)) in kidx for a in g_elems for b in g_elems)
    checks["C1"] = len(g_elems) == ORDER and closed
    values["group_order"] = len(g_elems)
    values["closed"] = bool(closed)
    if not checks["C1"]:
        raise Halt(finish())

    # ---- C2: conjugacy classes --------------------------------------------------
    unassigned, classes = set(range(len(g_elems))), []
    while unassigned:
        seed = min(unassigned)
        orbit = set()
        for g in g_elems:
            gx = qmul(g, g_elems[seed])
            right = g if mut == "conj_without_inverse" else qconj(g)
            orbit.add(kidx[qkey(qmul(gx, right))])
        classes.append(sorted(orbit))
        unassigned -= orbit
    classes.sort(key=lambda c: (len(c), -g_elems[c[0]][0]))
    nc = len(classes)
    sizes = np.array([len(c) for c in classes], float)
    checks["C2"] = nc == 9 and sizes.sum() == ORDER and all(ORDER % int(s) == 0 for s in sizes)
    values["n_classes"] = nc
    values["class_sizes"] = sizes.astype(int).tolist()
    if not checks["C2"]:
        raise Halt(finish())

    cls_of = {e: ci for ci, c in enumerate(classes) for e in c}
    # a unit quaternion is cos(t) + sin(t) n_hat, so its SU(2) trace is 2 cos(t)
    thetas = np.array([np.arccos(np.clip(g_elems[c[0]][0], -1.0, 1.0)) for c in classes])

    # ---- C3: characters by Burnside class-sum diagonalization -------------------
    cmat = np.zeros((nc, nc, nc))
    for i, ci in enumerate(classes):
        for j, cj in enumerate(classes):
            for a in ci:
                for b in cj:
                    cmat[i, j, cls_of[kidx[qkey(qmul(g_elems[a], g_elems[b]))]]] += 1.0
    cmat /= sizes[None, None, :]  # class-sum structure constants, per target element

    rng = np.random.default_rng(11)  # fixed seed: the mix only has to be generic
    mix = sum(rng.normal() * cmat[i] for i in range(nc))
    chars = []
    for vec in np.real(np.linalg.eig(mix)[1]).T:
        omega = vec / vec[0]  # central character, identity class normalized to 1
        deg = np.sqrt(ORDER / np.sum(omega**2 / sizes))
        chars.append(omega * deg / sizes)
    chars = np.array(chars)
    if mut == "perturb_char":
        chars[3] *= 1.05
    dims = np.round(chars[:, 0]).astype(int)
    order_by_dim = np.argsort(dims, kind="stable")
    chars, dims = chars[order_by_dim], dims[order_by_dim]

    gram = np.array(
        [[np.sum(sizes * chars[i] * chars[j]) / ORDER for j in range(nc)] for i in range(nc)]
    )
    ortho_err = float(np.max(np.abs(gram - np.eye(nc))))
    checks["C3"] = ortho_err < TOL
    values["orthonormality_error"] = ortho_err
    values["dims"] = sorted(dims.tolist())
    values["sum_dims_squared"] = int(np.sum(dims**2))

    # ---- C4: identify Q, the defining 2-dim rep ---------------------------------
    chi_q_target = 2.0 * np.cos(thetas)
    two_dim = [i for i in range(nc) if dims[i] == 2]
    i_q = min(two_dim, key=lambda i: np.max(np.abs(chars[i] - chi_q_target)))
    i_qp = [i for i in two_dim if i != i_q][0]
    if mut == "q_is_qprime":
        i_q, i_qp = i_qp, i_q
    q_err = float(np.max(np.abs(chars[i_q] - chi_q_target)))
    checks["C4"] = q_err < TOL
    values["defining_rep_error"] = q_err

    # ---- C5: the three connections tau = Sym^2(sigma) ---------------------------
    def chi_of_square(chi):
        return np.array(
            [chi[cls_of[kidx[qkey(qmul(g_elems[classes[c][0]], g_elems[classes[c][0]]))]]]
             for c in range(nc)]
        )  # fmt: skip

    def sym2(chi):
        if mut == "sym2_as_square":
            return chi**2
        return 0.5 * (chi**2 + chi_of_square(chi))

    tau = {"trivial": np.ones(nc), "standard": sym2(chars[i_q]), "galois": sym2(chars[i_qp])}

    def decompose(tvec):
        return np.array([np.sum(sizes * tvec * chars[i]) / ORDER for i in range(nc)])

    dec = {name: decompose(t) for name, t in tau.items()}
    irreducible = {
        name: (abs(np.sum(d**2) - 1.0) < 1e-6 and abs(tau[name][0] - 3.0) < 1e-6)
        for name, d in dec.items()
        if name != "trivial"
    }
    hit = {name: int(np.argmax(d)) for name, d in dec.items()}
    checks["C5"] = all(irreducible.values()) and hit["standard"] != hit["galois"]
    values["tau_dims"] = {name: float(t[0]) for name, t in tau.items()}
    values["tau_decomposition_dims"] = {
        name: [int(dims[i]) for i, mv in enumerate(d) if mv > 0.5] for name, d in dec.items()
    }

    # ---- C6/C7: McKay graph, derived by tensoring with Q ------------------------
    edge_min = 1.5 if mut == "edge_threshold" else 0.5
    adj = np.zeros((nc, nc))
    for i in range(nc):
        prod = chars[i] * chars[i_q]
        for j in range(nc):
            adj[i, j] = 1.0 if np.sum(sizes * prod * chars[j]) / ORDER > edge_min else 0.0
    checks["C6"] = bool(np.allclose(adj.dot(dims.astype(float)), 2.0 * dims.astype(float)))
    values["n_edges"] = int(adj.sum() / 2)

    root = int(np.where(dims == 1)[0][0])
    dist, frontier = {root: 0}, [root]
    while frontier:
        nxt = []
        for a in frontier:
            for b in range(nc):
                if adj[a, b] > 0.5 and b not in dist:
                    dist[b] = dist[a] if mut == "bfs_forget_increment" else dist[a] + 1
                    nxt.append(b)
        frontier = nxt
    signatures = [(int(dims[i]), dist.get(i)) for i in range(nc)]
    checks["C7"] = len(dist) == nc and len(set(signatures)) == nc
    values["signatures"] = [list(s) for s in signatures]
    if not checks["C7"]:
        raise Halt(finish())

    # ---- C8: the table, and the comparison against section 6.1 -----------------
    def first_level(irow, tvec, nmax=24):
        for n in range(nmax):
            v = chi_su2(n, thetas, mut)
            if np.sum(sizes * v * tvec * chars[irow]) / ORDER > 0.5:
                return n
        return None

    for i in sorted(range(nc), key=lambda i: (dist[i], dims[i])):
        table[(int(dims[i]), dist[i])] = tuple(
            first_level(i, tau[c]) for c in ("trivial", "standard", "galois")
        )

    target = dict(DOC)
    if mut == "doc_typo":
        target[(5, 4)] = (4, 2, 3)
    matched = 0
    for sig, want in sorted(target.items()):
        got = table.get(sig)
        ok = got is not None and tuple(got) == want
        matched += int(ok)
        comparison.append(
            {"dim": sig[0], "dist": sig[1], "target": list(want),
             "independent": None if got is None else list(got), "match": ok}
        )  # fmt: skip
    checks["C8"] = matched == len(target)
    values["rows_matched"] = matched
    values["rows_total"] = len(target)
    return finish()


def run(mut=None):
    try:
        return reconstruct(mut)
    except Halt as halt:
        return halt.args[0]


def label(state):
    return "PASS" if state else ("FAIL" if state is False else "n/a ")


def report(res):
    print("=" * 78)
    print("INDEPENDENT RECONSTRUCTION of the M8.2 0-form first-occurrence table")
    print("  explicit quaternions -> conjugacy classes -> Burnside class-sums")
    print("  no method, no fixture and no import shared with m8_2_first_occurrence.py")
    print("=" * 78)
    p = res["provenance"]
    print(f"  repository commit : {p['commit']}")
    print(f"  environment       : python {p['python']}, numpy {p['numpy']}, {p['platform']}")
    print(f"  seed              : {p['seed']} (fixed; the Burnside mix only has to be generic)")
    print("-" * 78)
    for cid, text in CHECKS.items():
        print(f"  {cid}  {text:<58} {label(res['checks'].get(cid))}")
    v = res["values"]
    print("\nDERIVED (nothing below is read from the M8.2 script or the document)")
    print(f"  group order / closed under multiplication : {v['group_order']} / {v['closed']}")
    print(f"  conjugacy classes                         : {v['n_classes']}")
    print(f"  class sizes                               : {v['class_sizes']}")
    print(f"  irrep dims                                : {v['dims']}")
    print(f"  sum of squared dims (= |G|)               : {v['sum_dims_squared']}")
    print(f"  orthonormality error                      : {v['orthonormality_error']:.2e}")
    print(f"  |chi_Q - 2cos(theta)|                     : {v['defining_rep_error']:.2e}")
    print(f"  McKay graph edges                         : {v['n_edges']}")
    for name, comp in v["tau_decomposition_dims"].items():
        print(f"  tau[{name:<8}] dim {v['tau_dims'][name]:.0f} -> irrep dims {comp}")

    print("\nRECONSTRUCTED TABLE, first-occurrence level n   [eigenvalue n(n+2)/R^2]")
    print(f"  {'dim':>4}{'dist':>6}{'trivial':>10}{'standard':>10}{'galois':>10}")
    for key, row in res["table"].items():
        dim, dst = key[3:].split("_d")
        print(f"  {dim:>4}{dst:>6}{row[0]:>10}{row[1]:>10}{row[2]:>10}")

    print("\nCOMPARISON against the pre-registration section 6.1")
    for row in res["comparison"]:
        verdict = "match" if row["match"] else "*** MISMATCH ***"
        print(
            f"  dim {row['dim']} d={row['dist']}: target {tuple(row['target'])}  "
            f"independent {None if row['independent'] is None else tuple(row['independent'])}"
            f"  {verdict}"
        )
    print(f"\nVERDICT: {v['rows_matched']}/{v['rows_total']} rows reproduced independently")
    print("  claim label: independent-method reproduction, NOT blind")
    print("  NOT verified here: the coexact 1-form entry rule (ASSERTED, no target)")


def mutation_tests():
    """Every PASS line above must be able to go red.  Show it, one defect at a time."""
    base = run()
    print("=" * 78)
    print("MUTATION TESTS: each check is shown to fail under a deliberate defect")
    print("=" * 78)
    print(f"  baseline: {sum(1 for s in base['checks'].values() if s)}/{len(CHECKS)} PASS")
    if not all(base["checks"].get(c) for c in CHECKS):
        print("  baseline is not all-PASS; aborting")
        return 1
    covered, failures = set(), 0
    for mut, (targets, what) in MUTATIONS.items():
        res = run(mut)
        want = targets.split(",")
        red = [c for c in want if res["checks"].get(c) is False]
        ok = len(red) == len(want)
        failures += int(not ok)
        covered.update(red)
        also = [c for c, s in res["checks"].items() if s is False and c not in want]
        print(f"\n  {mut}")
        print(f"    defect        : {what}")
        print(f"    target        : {targets}")
        print(
            f"    went red      : {','.join(red) if red else '(none)'}   {'OK' if ok else 'BROKEN'}"
        )
        if also:
            print(f"    also red      : {','.join(sorted(also))}")
    uncovered = [c for c in CHECKS if c not in covered]
    print("\n" + "=" * 78)
    print(f"  checks with a failing mutation : {len(covered)}/{len(CHECKS)}")
    print(f"  uncovered                      : {uncovered if uncovered else '(none)'}")
    print(f"  mutations that did NOT redden their target : {failures}")
    print("=" * 78)
    return 1 if (failures or uncovered) else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", metavar="PATH", help="write the structured result here")
    ap.add_argument("--mutation-tests", action="store_true", help="run the mutation suite only")
    args = ap.parse_args()
    if args.mutation_tests:
        raise SystemExit(mutation_tests())
    res = run()
    report(res)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(res, fh, indent=2)
        print(f"\nwrote {args.json}")
    raise SystemExit(0 if all(res["checks"].get(c) for c in CHECKS) else 1)


if __name__ == "__main__":
    main()
