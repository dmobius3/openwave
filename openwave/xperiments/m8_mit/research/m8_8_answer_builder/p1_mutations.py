"""P1 mutation battery: every declared P1 check must be shown able to reject.

Target-free. Mutations are described structurally and every substitute value is
read out of a source, never written here.

Two rules this battery follows, both learned the hard way:

  1. A mutation must run the SAME acceptance predicate the driver runs, which is
     `p1.assess`, and must show the TARGETED check id among the failures. An
     earlier version watched an intermediate change and called that a reddened
     gate; observing sensitivity is not observing rejection.
  2. Coverage is enforced against `p1.CHECK_IDS`. Nine named fault modes are not
     the same claim as "every gate is operative", and the earlier battery made
     the second claim while establishing only the first.

Mutations are one of two kinds, and the kind is reported:
  FILE      the source artifacts are edited on disk, in a scratch copy
  INJECTED  a derivation step is replaced in memory, to reach a check whose
            precondition cannot be produced by editing an input
"""

import contextlib
import json
import pathlib
import re
import shutil
import sys
import tempfile

import p1
import p1_enumerate
import p1_signatures
import p1_values

CASES = []


def case(kind, *targets, reason=None, require_green=()):
    """reason: substring every targeted check's detail must contain.
    require_green: check ids that must NOT reject, so a case cannot pass by
    tripping an earlier gate instead of the one it claims to test."""
    def deco(fn):
        CASES.append((fn.__name__, kind, tuple(targets), fn,
                      reason, tuple(require_green)))
        return fn
    return deco


@contextlib.contextmanager
def patched(mod, name, repl):
    orig = getattr(mod, name)
    setattr(mod, name, repl)
    try:
        yield orig
    finally:
        setattr(mod, name, orig)


def case_passed(targets, bad, raised, reason=None, details=None, require_green=()):
    """A case passes only if EVERY named target rejected FOR THE STATED REASON,
    nothing raised, and every required-green check stayed green.

    The reason clause is the important one. Making the pin an execution
    precondition silently broke three source-run cases: their mutations changed
    executable bytes without repinning, so the program was never run and the
    check failed as "source is not pinned". They still counted as red, and the
    parser paths they claimed to cover were never reached. Rejection is not
    coverage unless it is rejection by the branch under test."""
    if raised or not set(targets) <= set(bad):
        return False
    if any(g in bad for g in require_green):
        return False
    if reason is not None:
        d = details or {}
        want = reason if isinstance(reason, dict) else {t: reason for t in targets}
        return all(sub in str(d.get(t, "")) for t, sub in want.items())
    return True


def meta_self_test():
    """The meta-gate is itself a check, so it has to be able to fail.

    The first row is the case the redline used to establish the defect: two
    named targets, only one rejecting, previously reported as a success."""
    cases = [
        (("VAL_REQUIRED_CARRIERS", "VAL_FAMILIES"), {"VAL_REQUIRED_CARRIERS"}, False, False),
        (("VAL_REQUIRED_CARRIERS", "VAL_FAMILIES"),
         {"VAL_REQUIRED_CARRIERS", "VAL_FAMILIES"}, False, True),
        (("ADJ_GROUP",), {"ADJ_GROUP"}, False, True),
        (("ADJ_GROUP",), {"ADJ_CONSTRUCTION"}, False, False),
        (("ADJ_GROUP",), {"ADJ_GROUP"}, True, False),          # raised, not rejected
        (("ADJ_GROUP",), set(), False, False),
    ]
    results = [(t, b, r, want, case_passed(t, b, r)) for t, b, r, want in cases]
    agree = [x for x in results if x[3] == x[4]]
    ok = (len(agree) == len(results)
          and any(x[3] for x in results) and any(not x[3] for x in results))
    print("  meta-gate self-test:")
    for t, b, r, want, got in results:
        print(f"    [{'ok  ' if want == got else 'FAIL'}] targets={sorted(t)} "
              f"rejected={sorted(b)} raised={r} -> {got} (want {want})")
    print(f"    {len(agree)}/{len(results)} agree, both outcomes present: "
          f"{any(x[3] for x in results) and any(not x[3] for x in results)}")
    return ok


# Stops that are not per-check ids but must still be reachable as targets.
EXTRA_TARGETS = {"REGISTRY"}


def failed_ids(paths):
    """Run the driver's acceptance predicate; return everything that rejected.

    Includes non-check stops keyed by their prefix, so registry integrity is
    testable: it is fatal but it is not a check id."""
    bad, _ = failed_with_details(paths)
    return bad


def failed_with_details(paths):
    """(set of rejecting ids, {id: detail}). Details are what make a case able
    to assert WHY a check rejected."""
    try:
        checks, stops, _data = p1.assess(**paths)
    except Exception as e:                                   # noqa: BLE001
        return {"__raised__": f"{type(e).__name__}: {e}"}, {}
    bad = {cid for cid, ok, _ in checks if not ok}
    details = {cid: d for cid, ok, d in checks if not ok}
    for s in stops:
        k = s.split(":", 1)[0]
        if k in EXTRA_TARGETS:
            bad.add(k)
            details[k] = s
    return bad, details


def repin(paths, tmp, role, mutated_path):
    """A pin file that accepts the mutated source but keeps its gate inventory.

    Needed whenever a case wants to test the RUNNER: without it the pin refuses
    execution and the case measures fail-closed behaviour it did not intend."""
    pins = json.loads(pathlib.Path(p1.__file__).with_name(p1.SOURCE_PINS).read_text())
    pins["sources"][role] = p1.sha256_file(mutated_path)
    pf = tmp / f"repinned_{role}.json"
    pf.write_text(json.dumps(pins))
    return str(pf)


def _span_of_dict(txt, name):
    m = re.search(rf"\b{name}\s*=\s*\{{", txt)
    if not m:
        raise ValueError(f"{name} not found")
    i = m.end() - 1
    depth, j = 0, i
    while j < len(txt):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                return i, j + 1
        j += 1
    raise ValueError(f"{name} unterminated")


def edit_dict(txt, name, label, action):
    """action(expr) -> replacement text, or None to delete the entry."""
    i, j = _span_of_dict(txt, name)
    body = txt[i:j]
    m = re.search(rf'"{label}"\s*:\s*([^,}}]+),?', body)
    if not m:
        raise ValueError(f"{label} not inside {name}")
    repl = action(m.group(1).strip())
    new = body[:m.start()] + ("" if repl is None else f'"{label}": {repl},') + body[m.end():]
    return txt[:i] + new + txt[j:]


def copy_inputs(paths, tmp):
    out = {}
    for k, v in paths.items():
        d = tmp / k
        d.mkdir(exist_ok=True)
        p = d / pathlib.Path(v).name
        shutil.copy(v, p)
        out[k] = str(p)
    return out


# ------------------------------------------------------------ adjudicates

@case("FILE", "ADJ_GROUP")
def adjudicates_group(paths, tmp):
    m = copy_inputs(paths, tmp)
    pathlib.Path(m["group"]).write_text(pathlib.Path(m["group"]).read_text() + "\n")
    return m, {}


@case("FILE", "ADJ_CONSTRUCTION")
def adjudicates_construction(paths, tmp):
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["construction"])
    p.write_text(p.read_text() + "\n")
    return m, {}


# ------------------------------------------------------------ enumeration

@case("INJECTED", "ENUM_DIGEST", "ENUM_RANK_KEYS")
def enum_wrong_sort_key(paths, tmp):
    """The § 4.2 misread: sort by the § 4.4 value normalization instead."""
    def wrong(packet):
        elems = p1_enumerate.close_group(p1_enumerate.load_generators(packet))
        return sorted(elems, key=lambda q: tuple(t for x in q.c for t in x.triple()))
    return paths, {(p1_enumerate, "enumerate_group"): wrong}


@case("INJECTED", "ENUM_ORDER", "ENUM_BYTES")
def enum_short_closure(paths, tmp):
    """Drop an element, so the closure is not the whole group."""
    orig = p1_enumerate.enumerate_group
    return paths, {(p1_enumerate, "enumerate_group"): lambda p: orig(p)[:-1]}


# ------------------------------------------------------------- signatures

def _doctor_run(**overrides):
    orig = p1_signatures.run

    def repl(*a, **k):
        r = orig(*a, **k)
        r = dict(r)
        r.update(overrides)
        return r
    return {(p1_signatures, "run"): repl}


@case("INJECTED", "SIG_COUNT")
def sig_count(paths, tmp):
    return paths, _doctor_run(n_irreps=8)


@case("INJECTED", "SIG_SUMSQ")
def sig_sumsq(paths, tmp):
    return paths, _doctor_run(sum_of_squares=119)


@case("INJECTED", "SIG_DIMS")
def sig_dims(paths, tmp):
    return paths, _doctor_run(dims=[1, 2, 2, 3, 3, 4, 4, 5, 7])


@case("INJECTED", "SIG_LABELS_UNIQUE")
def sig_labels_ambiguous(paths, tmp):
    return paths, _doctor_run(n_label_solutions=2)


@case("INJECTED", "SIG_SEPARATION")
def sig_collision(paths, tmp):
    return paths, _doctor_run(signatures_distinct=False,
                              signature_collisions=[("Ra", "Rb")])


@case("INJECTED", "SIG_ORDERS")
def sig_wrong_order(paths, tmp):
    return paths, _doctor_run(orders={"s": 6, "t": 10, "st": 6})


@case("FILE", "SIG_LABELS_UNIQUE")
def mckay_edge_moved(paths, tmp):
    """Re-attach one McKay edge to a different node, in BOTH declaring sources.

    The edge is read out of the declaration rather than written here, so the
    diagram's shape is not restated in this file."""
    m = copy_inputs(paths, tmp)
    labels, edges = p1_signatures.parse_mckay_declaration(paths["theory"])
    a, b = edges[-1]
    other = next(x for x in labels if x not in (a, b))
    for k in ("theory", "repro"):
        p = pathlib.Path(m[k])
        p.write_text(p.read_text().replace(f'("{a}", "{b}")', f'("{other}", "{b}")', 1))
    return m, {}


# ----------------------------------------------------------------- values

@case("FILE", "VAL_PARSE")
def value_unparseable(paths, tmp):
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["note"])
    p.write_text(re.sub(r"(T\^2\(R3\)=)\S+", r"\1(4/5)psi^-2", p.read_text()))
    return m, {}


def note_labels(path):
    """Labels of the note's value table, in the order they are stated."""
    return [l for l, _ in p1_values.source_method_note(path)]


@case("FILE", "VAL_NO_INTRA_SOURCE_CONFLICT")
def value_duplicate_conflict(paths, tmp):
    """State one label twice, with different values, inside ONE source.

    The blocking case: a dict literal keeps the last duplicate key and a
    setdefault keeps the first, so before the fix neither conflict reached
    adjudication. Both the victim label and the substitute expression are
    chosen at runtime from the source, so no membership or value is written."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["repro"])
    txt = p.read_text()
    i, j = _span_of_dict(txt, "exact")
    entries = re.findall(r'"(R\d)"\s*:\s*([^,}]+)', txt[i:j])
    victim, donor = entries[0][0], next(
        e for e in entries if e[1].strip() != entries[0][1].strip())
    p.write_text(edit_dict(txt, "exact", victim,
                           lambda e: f'{e}, "{victim}": {donor[1].strip()}'))
    return m, {}


@case("FILE", "VAL_REQUIRED_CARRIERS", "VAL_FAMILIES")
def value_theory_side_missing(paths, tmp):
    """Remove a row from the THEORY source only, leaving both platform
    renderings intact. Counting renderings would still see two and pass."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    txt = p.read_text()
    i, j = _span_of_dict(txt, "tgt_int")
    victim = re.findall(r'"(R\d)"\s*:', txt[i:j])[0]
    p.write_text(edit_dict(txt, "tgt_int", victim, lambda e: None))
    return m, {}


@case("FILE", "VAL_AGREE")
def value_disagreement(paths, tmp):
    """Retarget one label in one source onto a different row's stated value.

    Both labels are picked at runtime from the table, and the substitute is
    copied from the table, so this case states neither."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["note"])
    txt = p.read_text()
    occ = p1_values.source_method_note(paths["note"])
    victim, donor = occ[0][0], next(o for o in occ if o[1] != occ[0][1])
    p.write_text(re.sub(rf"(T\^2\({victim}\)=)\S+",
                        lambda g: g.group(1) + donor[1], txt))
    return m, {}


ANCHOR = "R7"        # named in the PUBLIC protocol § 5.2 and § 5.4, not private


@case("FILE", "ANCHOR_NOT_SELF_INVERSE")
def anchor_self_inverse(paths, tmp):
    """Make the anchor its own inverse in ALL sources, so they still agree.

    The substitute is the anchor's OWN expression divided by itself, which is
    exactly one without naming any other row or any value. An earlier version
    borrowed a row known to be unit-valued, which stated a value fact in a file
    classified target-free."""
    m = copy_inputs(paths, tmp)

    note = pathlib.Path(m["note"])
    ntxt = note.read_text()
    note.write_text(re.sub(rf"(T\^2\({ANCHOR}\)=)(\S+)",
                           lambda g: f"{g.group(1)}({g.group(2)})/({g.group(2)})", ntxt))

    repro = pathlib.Path(m["repro"])
    repro.write_text(edit_dict(repro.read_text(), "exact", ANCHOR,
                               lambda e: f"({e})/({e})"))

    th = pathlib.Path(m["theory"])
    th.write_text(re.sub(rf'(T2\["{ANCHOR}"\]\s*/\s*\()(.+?)(\)\s*-\s*1)',
                         lambda g: f"{g.group(1)}({g.group(2)})/({g.group(2)}){g.group(3)}",
                         th.read_text()))
    return m, {}


# ------------------------------------------------------- source program runs

def _break_own_gate(path, dict_name):
    """Point one entry of a program's own expected-value table at another entry,
    so the program's internal gate fails. Both entries are read from the file."""
    txt = pathlib.Path(path).read_text()
    i, j = _span_of_dict(txt, dict_name)
    entries = re.findall(r'"(R\d)"\s*:\s*([^,}]+)', txt[i:j])
    victim = entries[0][0]
    donor = next(e for e in entries if e[1].strip() != entries[0][1].strip())
    return edit_dict(txt, dict_name, victim, lambda e: donor[1].strip())


@case("FILE", "SRC_RUN_THEORY", reason="exit 1",
      require_green=("SRC_PINNED",))
def src_run_theory_fails(paths, tmp):
    """The program runs and its OWN gate fails. Repinned so the runner is
    reached: without the repin this measured pin refusal instead."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    p.write_text(_break_own_gate(m["theory"], "tgt_half"))
    m["source_pins"] = repin(paths, tmp, "theory", m["theory"])
    return m, {}


@case("FILE", "SRC_RUN_PLATFORM", reason="exit 1",
      require_green=("SRC_PINNED",))
def src_run_platform_fails(paths, tmp):
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["repro"])
    p.write_text(_break_own_gate(m["repro"], "exact"))
    m["source_pins"] = repin(paths, tmp, "repro", m["repro"])
    return m, {}


# ------------------------------------------- wrong-generator-class diagnostic

def _doctor_diag(fn):
    orig = p1_signatures.generator_class_diagnostic

    def repl(*a, **k):
        return fn(orig(*a, **k))
    return {(p1_signatures, "generator_class_diagnostic"): repl}


@case("INJECTED", "GEN_CANDIDATE_SET")
def gen_candidate_count(paths, tmp):
    def f(d):
        d = dict(d)
        d["candidates"] = d["candidates"][:-1]
        return d
    return paths, _doctor_diag(f)


@case("INJECTED", "GEN_CLASS_PARTITION")
def gen_class_partition(paths, tmp):
    def f(d):
        d = dict(d)
        d["classes"] = [dict(c, members=c["members"][:2]) for c in d["classes"]]
        return d
    return paths, _doctor_diag(f)


@case("INJECTED", "GEN_UNIFORM_IN_CLASS")
def gen_not_uniform(paths, tmp):
    def f(d):
        d = dict(d)
        d["classes"] = [dict(c, uniform_within_class=False) for c in d["classes"]]
        return d
    return paths, _doctor_diag(f)


@case("INJECTED", "GEN_WRONG_CLASS_SAME_SIGSET")
def gen_sigset_differs(paths, tmp):
    def f(d):
        d = dict(d)
        d["classes"] = [dict(c, same_signature_set=c["is_packet_class"])
                        for c in d["classes"]]
        return d
    return paths, _doctor_diag(f)


@case("INJECTED", "GEN_WRONG_CLASS_IS_GALOIS_PARTNER_MAP")
def gen_permutation_not_partner(paths, tmp):
    def f(d):
        d = dict(d)
        d["classes"] = [
            c if c["is_packet_class"]
            else dict(c, permutation={k: k for k in (c["permutation"] or {})})
            for c in d["classes"]]
        return d
    return paths, _doctor_diag(f)


@case("FILE", "IDY_PARSED")
def identity_statements_removed(paths, tmp):
    """Remove the identity statements from the record entirely."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["note"])
    txt = p.read_text()
    # Replace the EXACT spans that were parsed. An earlier version used a
    # backtick-delimited pattern, which ran across the fenced value table and
    # destroyed it, so the run raised instead of rejecting on this check.
    for i in p1_values.identities_from_note(paths["note"]):
        txt = txt.replace(f"`{i['span']}`", "`removed`", 1)
    p.write_text(txt)
    return m, {}


@case("FILE", "IDY_HOLDS")
def identity_stated_value_wrong(paths, tmp):
    """Point one identity's stated result at a different identity's result.

    Both are read from the record; neither is written here."""
    m = copy_inputs(paths, tmp)
    idents = p1_values.identities_from_note(paths["note"])
    victim = idents[0]
    donor = next(i for i in idents if i["stated"] != victim["stated"])
    p = pathlib.Path(m["note"])
    p.write_text(p.read_text().replace(
        f"= {victim['stated']}`", f"= {donor['stated']}`", 1))
    return m, {}


def _first_ratio(idents):
    return next(i for i in idents if i["kind"] == "ratio")


@case("INJECTED", "IDY_PAIRWISE_SENSITIVITY")
def identity_ratio_dropped(paths, tmp):
    """Drop one ratio identity, leaving its Galois exchange bound by nothing.

    This is what the P2 recommendation rests on: NOT that a lone ratio is the
    only detection of the hazard, but that each ratio is the sole binder of its
    own pair, so removing one leaves one exchange unbound."""
    orig = p1_values.identities_from_note

    def repl(path):
        ids = orig(path)
        return [i for i in ids if i["id"] != _first_ratio(ids)["id"]]
    return paths, {(p1_values, "identities_from_note"): repl}


@case("INJECTED", "IDY_PAIRWISE_SENSITIVITY")
def identity_ratio_factors_misdeclared(paths, tmp):
    """Make both ratios declare the SAME factor pair.

    The redline asked for a case that exchanges the two ratio definitions while
    retaining both. With identities parsed and each id DERIVED from its own
    factors, an exchange is a pure relabelling and cannot change any verdict:
    a name that disagrees with what it computes is not constructible. The
    defect that IS reachable is a misdeclared factor pair, which leaves one
    exchange with no owner, and that is what this exhibits."""
    orig = p1_values.identities_from_note

    def repl(path):
        ids = orig(path)
        first = _first_ratio(ids)
        out = []
        for i in ids:
            if i["kind"] == "ratio" and i["id"] != first["id"]:
                out.append(dict(i, factors=list(first["factors"]),
                                id=i["id"] + "+misdeclared"))
            else:
                out.append(i)
        return out
    return paths, {(p1_values, "identities_from_note"): repl}


@case("INJECTED", "IDY_SHAPE")
def identity_products_removed(paths, tmp):
    """Keep both ratios, drop both products.

    The audit's finding: every remaining identity holds, every exchange still
    has its owner, and the whole run stayed green on an identity set that was
    half missing. Nothing required the declared slot count."""
    orig = p1_values.identities_from_note

    def repl(path):
        return [i for i in orig(path) if i["kind"] != "product"]
    return paths, {(p1_values, "identities_from_note"): repl}


@case("FILE", "SRC_PINNED")
def source_revision_unpinned(paths, tmp):
    """Change a source record's bytes without re-pinning it.

    Content is preserved, so every other value check still passes: the point is
    that P1 must notice it was handed a revision it was not authorized for."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["note"])
    p.write_text(p.read_text() + "\n<!-- unpinned revision -->\n")
    return m, {}


@case("FILE", "SRC_RUN_THEORY", reason="gate summaries",
      require_green=("SRC_PINNED",))
def src_run_second_summary(paths, tmp):
    """Append a SECOND, failing gate summary after the passing one.

    The audit's finding: taking the first regex match let the reassuring line
    win while the real verdict below it was never read."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    lines = p.read_text().splitlines(keepends=True)
    # Insert beside the program's OWN summary line. Appending at module level
    # does nothing here: main() ends in sys.exit, so trailing code never runs
    # and the mutation silently did not mutate.
    k = next(i for i, l in enumerate(lines) if "GATES:" in l and "print" in l)
    indent = lines[k][:len(lines[k]) - len(lines[k].lstrip())]
    lines.insert(k + 1, f'{indent}print("GATES: 0/1 pass")\n')
    p.write_text("".join(lines))
    m["source_pins"] = repin(paths, tmp, "theory", m["theory"])
    return m, {}


@case("FILE", "SRC_RUN_THEORY", reason="not a bare final line",
      require_green=("SRC_PINNED",))
def src_run_summary_not_bare(paths, tmp):
    """A summary embedded in a longer line. An unanchored search took it."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    lines = p.read_text().splitlines(keepends=True)
    k = next(i for i, l in enumerate(lines) if "GATES:" in l and "print" in l)
    lines[k] = lines[k].replace('f"GATES:', 'f"prefix GATES:').replace(
        '"GATES:', '"prefix GATES:')
    if "prefix" not in lines[k]:
        lines[k] = lines[k].rstrip("\n") + "\n"
    p.write_text("".join(lines))
    m["source_pins"] = repin(paths, tmp, "theory", m["theory"])
    return m, {}


@case("FILE", "SRC_RUN_THEORY", reason="summary counts",
      require_green=("SRC_PINNED",))
def src_run_count_mismatch(paths, tmp):
    """Two gate ids reported, summary says 1/1. Comparing numerator to
    denominator, and separately ids to pins, left this green."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    p.write_text('print("PASS  [a] ok")\nprint("PASS  [b] ok")\n'
                 'print("GATES: 1/1 pass")\n')
    pins = json.loads(pathlib.Path(p1.__file__).with_name(p1.SOURCE_PINS).read_text())
    pins["sources"]["theory"] = p1.sha256_file(m["theory"])
    pins["gate_ids"]["theory"] = ["a", "b"]
    pf = tmp / "count.json"
    pf.write_text(json.dumps(pins))
    m["source_pins"] = str(pf)
    return m, {}


@case("INJECTED", "REGISTRY")
def registry_check_dropped(paths, tmp):
    """Drop a check from what production registry validation actually sees.

    The earlier version of this case removed a check AFTER validation had run
    and then appended the stop it expected, so it tested itself rather than the
    production logic. Now `validate_registry` is a function and the mutation
    drives its INPUT, so the stop is produced by the real code path."""
    orig = p1.validate_registry

    def repl(emitted):
        return orig([c for c in emitted][:-1])
    return paths, {(p1, "validate_registry"): repl}


@case("FILE", "SRC_RUN_THEORY", reason="summary counts",
      require_green=("SRC_PINNED",))
def src_gate_inventory_gutted(paths, tmp):
    """A repinned revision that prints a valid summary and runs almost nothing.

    The audit's finding: `GATES: 1/1 pass` satisfies numerator == denominator,
    so a source could shed its entire gate inventory and stay green. The pin is
    updated too, so this is exactly the legitimate-repin scenario rather than a
    hash mismatch in disguise."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    p.write_text('print("GATES: 1/1 pass")\n')
    pins = json.loads(pathlib.Path(p1.__file__).with_name(p1.SOURCE_PINS).read_text())
    pins["sources"]["theory"] = p1.sha256_file(str(p))
    pf = tmp / "repinned.json"
    pf.write_text(json.dumps(pins))
    m["source_pins"] = str(pf)
    return m, {}


@case("FILE", "SRC_RUN_THEORY", "SRC_PINNED",
      reason={"SRC_RUN_THEORY": "refusing to execute", "SRC_PINNED": "unpinned"})
def unpinned_source_refused(paths, tmp):
    """An unpinned source must be REFUSED, not run and then reported.

    Rejection is what this case establishes. That the code never ran is a
    separate, stronger claim and needs a side effect to witness it, so it has
    its own experiment: fail_closed_self_test."""
    m = copy_inputs(paths, tmp)
    p = pathlib.Path(m["theory"])
    p.write_text(p.read_text() + "\n# unpinned revision\n")
    return m, {}


def fail_closed_self_test(paths):
    """Prove the unpinned source is never EXECUTED, not merely marked failed.

    The audit's finding had a security edge: P1 recorded the pin failure and ran
    the code anyway, so untrusted code had already executed by the time the
    verdict stopped. A marker file witnesses execution directly; asserting
    refusal without one would be taking the fix on trust."""
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="p1fc-"))
    try:
        m = copy_inputs(paths, tmp)
        marker = tmp / "EXECUTED"
        p = pathlib.Path(m["theory"])
        p.write_text(f'open({str(marker)!r}, "w").write("ran")\n' + p.read_text())
        bad = failed_ids(m)
        ran = marker.exists()
        rejected = isinstance(bad, set) and {"SRC_PINNED", "SRC_RUN_THEORY"} <= bad
        print("  fail-closed self-test (unpinned source):")
        print(f"    SRC_PINNED and SRC_RUN_THEORY both rejected: {rejected}")
        print(f"    marker file created, i.e. the code RAN:      {ran}")
        return rejected and not ran
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


@case("INJECTED", "IDY_PARSED", "IDY_PAIRWISE_SENSITIVITY")
def identity_factor_unknown_row(paths, tmp):
    """An identity naming a row that does not exist.

    Made IDY_PARSED false and then raised KeyError in the sensitivity loop, so
    the run lost the verdict the failing check was about to give."""
    orig = p1_values.identities_from_note

    def repl(path):
        ids = [dict(i) for i in orig(path)]
        ids[0]["factors"] = [("R9", 1)] + list(ids[0]["factors"][1:])
        return ids
    return paths, {(p1_values, "identities_from_note"): repl}


@case("INJECTED", "GEN_WRONG_CLASS_IS_GALOIS_PARTNER_MAP")
def partner_map_not_total(paths, tmp):
    """A row whose Galois conjugate signature matches nothing.

    Left a None in the map; pair construction then tried to order None against
    a string and raised, losing the complete verdict over a diagnostic."""
    orig = p1.galois_partner_map

    def repl(rows):
        m = dict(orig(rows))
        m[sorted(m)[0]] = None
        return m
    return paths, {(p1, "galois_partner_map"): repl}


def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--meta-self-test":
        ok = meta_self_test()
        print(f"  meta-gate {'PROVEN' if ok else 'BROKEN'}")
        return 0 if ok else 1
    if sys.argv[1:2] == ["--counts"]:
        pass
    elif len(sys.argv) != 6 and not (len(sys.argv) == 7 and sys.argv[1] == "--fail-closed"):
        print("usage: p1_mutations.py <group> <construction> <theory> <note> <repro>",
              file=sys.stderr)
        print("       p1_mutations.py --meta-self-test", file=sys.stderr)
        print("       p1_mutations.py --fail-closed <group> <construction> <theory> "
              "<note> <repro>", file=sys.stderr)
        return 2
    if sys.argv[1] == "--counts":
        # Totals GENERATED, so the report never carries a stale hard-coded number.
        kinds = {}
        for _n, k, _ts, _f, _r, _g in CASES:
            kinds[k] = kinds.get(k, 0) + 1
        named = {t for _n, _k, ts, _f, _r, _g in CASES for t in ts}
        with_reason = sum(1 for c in CASES if c[4] is not None)
        print(f"declared_checks={len(p1.CHECK_IDS)}")
        print(f"mutation_cases={len(CASES)}")
        for k in sorted(kinds):
            print(f"cases_{k.lower()}={kinds[k]}")
        print(f"cases_with_causal_reason={with_reason}")
        print(f"checks_named_by_a_case={len(named & set(p1.CHECK_IDS))}")
        return 0
    if sys.argv[1] == "--fail-closed":
        paths = dict(zip(("group", "construction", "theory", "note", "repro"), sys.argv[2:]))
        ok = fail_closed_self_test(paths)
        print(f"  fail-closed {'PROVEN' if ok else 'BROKEN'}")
        return 0 if ok else 1
    base = dict(zip(("group", "construction", "theory", "note", "repro"), sys.argv[1:]))

    covered = {t for _, _, ts, _, _, _ in CASES for t in ts}
    uncovered = [c for c in p1.CHECK_IDS if c not in covered]
    unknown = sorted(covered - set(p1.CHECK_IDS) - EXTRA_TARGETS)

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="p1mut-"))
    print("P1 mutation battery: each case must make its TARGET check reject\n")
    results = []
    try:
        for name, kind, targets, fn, reason, req_green in CASES:
            d = tmp / name
            d.mkdir()
            mpaths, injections = fn(base, d)
            with contextlib.ExitStack() as st:
                for (mod, attr), repl in injections.items():
                    st.enter_context(patched(mod, attr, repl))
                res = failed_with_details(mpaths)
                bad, details = res
            raised = isinstance(bad, dict) and "__raised__" in bad
            required = set(targets)
            # EVERY named target must reject, not any one of them. Intersecting
            # let a two-target case certify both when only one rejected, and
            # coverage was then counted from declarations rather than from what
            # was observed to fail. The predicate is factored out so it can be
            # self-tested; see meta_self_test.
            ok = case_passed(required, bad if not raised else set(), raised,
                             reason=reason, details=details, require_green=req_green)
            missing = sorted(required - bad) if not raised else sorted(required)
            results.append((name, kind, targets, ok, bad, missing))
            if reason and not raised:
                want = reason if isinstance(reason, dict) else {x: reason for x in required}
                bad_reason = [x for x, sub in want.items() if sub not in str(details.get(x, ""))]
            else:
                bad_reason = []
            green_violation = [g for g in req_green if not raised and g in bad]
            mark = "RED " if ok else "MISS"
            print(f"  {mark} [{kind:8s}] {name}")
            print(f"         targets: {', '.join(targets)}")
            if raised:
                print(f"         RAISED instead of rejecting: {bad['__raised__']}")
            elif ok:
                extra = sorted(bad - required)
                why = f' via "{reason}"' if reason else ""
                print(f"         rejected: {sorted(required)}{why}"
                      + (f"  (also: {extra})" if extra else ""))
            elif bad_reason:
                print(f"         WRONG REASON: {bad_reason} rejected, but not via "
                      f'"{reason}". Actual: '
                      + "; ".join(f"{x}={details.get(x)!r}" for x in bad_reason))
            elif green_violation:
                print(f"         TRIPPED AN EARLIER GATE instead: {green_violation}")
            else:
                print(f"         STILL GREEN: {missing}; failures were {sorted(bad)}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # Coverage is computed from checks OBSERVED failing in a successful case.
    observed = set()
    for _n, _k, targets, ok, _b, _m in results:
        if ok:
            observed |= set(targets)
    never = [c for c in p1.CHECK_IDS if c not in observed]

    red = sum(1 for *_x, ok, _b, _m in results if ok)
    # FILE and INJECTED support DIFFERENT claims and are reported apart. A FILE
    # case shows a realistic defect propagating through the whole derivation.
    # An INJECTED case replaces a derived field, so it shows the comparator is
    # wired and the id reddens, but says nothing about upstream propagation.
    # Summing them into one number overstates end-to-end fault coverage.
    for kind in ("FILE", "INJECTED"):
        sel = [r for r in results if r[1] == kind]
        print(f"  {sum(1 for r in sel if r[3])}/{len(sel)} {kind} cases "
              f"({'end-to-end through the derivation' if kind == 'FILE' else 'comparator wiring only'})")
    print(f"\n  {red}/{len(results)} cases made EVERY named target reject")
    print(f"  coverage: {len(observed & set(p1.CHECK_IDS))}/{len(p1.CHECK_IDS)} "
          f"declared checks OBSERVED rejecting")
    if unknown:
        print(f"  CASES TARGET UNDECLARED IDS: {unknown}")
    if uncovered:
        print(f"  CHECKS WITH NO MUTATION DECLARED: {uncovered}")
    if never:
        print(f"  CHECKS NEVER OBSERVED REJECTING: {never}")
    ok = red == len(results) and not never and not unknown
    print(f"  {'EVERY DECLARED P1 CHECK IS OPERATIVE' if ok else 'COVERAGE INCOMPLETE'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
