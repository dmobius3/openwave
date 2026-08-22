"""P1 driver: assemble the answer-packet data and report every check.

This SCRIPT is target-free. Its OUTPUT is answer-bearing and is written to
answer/build/out/ with a quarantine banner. Between P0a and P0b there is no
populated leak gate, so P1 output is quarantined by default rather than scanned.

Every check has a STABLE ID. `assess()` is the single acceptance predicate: it
returns the checks and the stop list, and the mutation battery calls it rather
than observing intermediates of its own. A mutation that only shows some
internal value moved does not establish that a gate rejects anything.
"""

import argparse
import hashlib
import json
import re
import pathlib
import sys

import p1_enumerate
import p1_signatures
import p1_values
from qphi_exact import Phi

PINNED = {
    "group_packet": "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9",
    "construction_packet": "df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06",
}

BANNER = (
    "QUARANTINED, ANSWER-BEARING.\n"
    "ONE authorized destination: an author-side informed redline session, which\n"
    "is permanently author-side once it receives this. Nothing else.\n"
    "BARRED from the PR, commit messages, the protocol, the sidecar, the public\n"
    "ciphertext README, the P5 independent reader, and any § 4.3 permitted input.\n"
    "Do not paste, quote or commit it anywhere.\n"
)

# The artifact is a P1 CROSS-CHECK TABLE. It is NOT a § 4.4 answer packet and
# must never be promoted, renamed or loaded as one: the packet is closed at
# seven keys and is built and validated at P2, from this data. The marker below
# exists so an answer-packet loader can reject this file positively rather than
# by noticing absent fields.
ARTIFACT_KIND = "m8_8-p1-crosscheck"

# The identity set is PARSED from the record by p1_values.identities_from_note.
# It is deliberately not written here: hard-coded memberships are quarantined
# STRUCTURE sitting in a file classified target-free, and the numeric leak gate
# cannot see structure. Anything that hardcodes a membership is answer-bearing
# no matter what the scanner says about it.

IDENTITY_SLOTS = 4          # § 4.4: the packet carries four `identities` slots

CHECK_IDS = [
    "ADJ_GROUP", "ADJ_CONSTRUCTION", "SRC_PINNED",
    "ENUM_ORDER", "ENUM_DIGEST", "ENUM_BYTES", "ENUM_RANK_KEYS",
    "SIG_COUNT", "SIG_SUMSQ", "SIG_DIMS", "SIG_LABELS_UNIQUE",
    "SIG_SEPARATION", "SIG_ORDERS",
    "SRC_RUN_THEORY", "SRC_RUN_PLATFORM",
    "VAL_PARSE", "VAL_NO_INTRA_SOURCE_CONFLICT", "VAL_REQUIRED_CARRIERS",
    "VAL_FAMILIES", "VAL_AGREE",
    "ANCHOR_NOT_SELF_INVERSE",
    "GEN_CANDIDATE_SET", "GEN_CLASS_PARTITION", "GEN_UNIFORM_IN_CLASS",
    "GEN_WRONG_CLASS_SAME_SIGSET", "GEN_WRONG_CLASS_IS_GALOIS_PARTNER_MAP",
    "IDY_PARSED", "IDY_SHAPE", "IDY_HOLDS", "IDY_PAIRWISE_SENSITIVITY",
]

SOURCE_PINS = "sources.pinned.json"     # written by p1_snapshot.py, gated here


def environment():
    """Recorded with every run: identical inputs on a different numerical stack
    are not the same evidence."""
    import platform
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for mod in ("numpy", "mpmath"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:                                     # noqa: BLE001
            out[mod] = "absent"
    return out

_RUN_CACHE = {}


def sha256_file(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def gate_ids_of(out):
    """The stable gate ids a source program reported, in order."""
    return re.findall(r"^\s*(?:PASS|FAIL)\s+\[([^\]]+)\]", out, re.M)


def run_source_program(path, expected_ids):
    """Execute a source artifact IN A SANDBOX and check its gate inventory.

    Three things this does not do, each learned from an audit:

      * It does not accept `n/n` as evidence. A program that exits early after
        printing `GATES: 1/1 pass` satisfies numerator == denominator while
        having run almost nothing, so a repinned revision could silently shed
        its gates. The EXACT expected gate-id set must match.
      * It does not run outside a scratch tree. The platform program writes an
        artifact next to its source, which made a verification check depend on
        write access and left side effects in a review workspace.
      * It is never reached for an unpinned source; see assess().

    Cached on (content, expected ids), so the battery pays once per distinct
    source rather than once per case."""
    import shutil
    import subprocess
    import tempfile
    key = (sha256_file(path), tuple(expected_ids or ()))
    if key in _RUN_CACHE:
        return _RUN_CACHE[key]
    p = pathlib.Path(path)
    sand = None
    try:
        sand = pathlib.Path(tempfile.mkdtemp(prefix="p1src-"))
        (sand / "scripts").mkdir()
        shutil.copy(p, sand / "scripts" / p.name)
        r = subprocess.run([sys.executable, p.name], cwd=str(sand / "scripts"),
                           capture_output=True, text=True, timeout=900,
                           env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin"})
        out = r.stdout + r.stderr
        # EXACTLY one summary, and it must be the last non-blank line. Taking the
        # first match accepted a program that printed a passing summary followed
        # by a failing one: the reassuring line wins and the real verdict is
        # never read.
        ms = re.findall(r"GATES:\s*(\d+)\s*/\s*(\d+)\s*pass", out)
        lines = [l for l in out.splitlines() if l.strip()]
        last = lines[-1].strip() if lines else ""
        # FULL-LINE match. An unanchored search accepted `prefix GATES: 1/1 pass
        # suffix`, so a summary embedded in other text counted as the verdict.
        last_is_summary = re.fullmatch(r"GATES:\s*(\d+)\s*/\s*(\d+)\s*pass", last)
        got = gate_ids_of(out)
        want = list(expected_ids or ())

        if r.returncode != 0:
            ok, detail = False, f"exit {r.returncode}"
        elif len(ms) != 1:
            ok, detail = False, f"{len(ms)} gate summaries, expected exactly 1"
        elif not last_is_summary:
            ok, detail = False, "gate summary is not a bare final line"
        elif ms[0][0] != ms[0][1] or "FAIL" in out:
            ok, detail = False, f"GATES: {ms[0][0]}/{ms[0][1]}, or a FAIL was printed"
        elif not want:
            ok, detail = False, "no expected gate inventory pinned for this source"
        elif len(set(got)) != len(got):
            ok, detail = False, f"duplicate gate ids: {sorted({g for g in got if got.count(g) > 1})}"
        # The summary COUNT must equal the inventory it claims to summarise.
        # Comparing numerator to denominator, and separately ids to pins, left
        # `GATES: 1/1 pass` green alongside two reported gate ids.
        elif int(ms[0][1]) != len(got) or len(got) != len(want):
            ok, detail = False, (f"summary counts {ms[0][1]} but {len(got)} gate ids ran "
                                 f"and {len(want)} are pinned")
        elif got != want:
            miss, extra = sorted(set(want) - set(got)), sorted(set(got) - set(want))
            ok, detail = False, (f"gate inventory differs: {len(got)} ran, {len(want)} pinned"
                                 + (f", missing {miss}" if miss else "")
                                 + (f", unexpected {extra}" if extra else "")
                                 + ("" if miss or extra else ", order changed"))
        else:
            ok, detail = True, f"GATES: {ms[0][0]}/{ms[0][1]} pass, {len(got)} pinned ids matched"
    except Exception as e:                                    # noqa: BLE001
        ok, detail = False, f"{type(e).__name__}: {e}"
    finally:
        if sand is not None:
            shutil.rmtree(sand, ignore_errors=True)
    _RUN_CACHE[key] = (ok, detail)
    return _RUN_CACHE[key]


def eval_identity(ident, T):
    acc = Phi(1)
    for lab, e in ident["factors"]:
        acc = acc * (T[lab] ** e)
    return acc


def as_phi(triple):
    a, b, c = triple
    return Phi(a, b) / Phi(c)


def conj_triple(t):
    """Galois conjugate of a normalized triple: phi -> 1 - phi."""
    return as_phi(t).conj().triple()


def galois_partner_map(rows):
    """label -> the label whose signature is the Galois conjugate of its own.

    Derived from the signatures. This is what makes the wrong-class permutation
    checkable without writing the permutation down."""
    by_sig = {p1_signatures.sig_tuple(v): k for k, v in rows.items()}
    out = {}
    for lab, v in rows.items():
        want = (v["dimension"], conj_triple(v["s"]),
                conj_triple(v["t"]), conj_triple(v["st"]))
        out[lab] = by_sig.get(want)
    return out


def validate_registry(emitted):
    """Registry integrity, as a function so the battery can drive its INPUT.

    Factored out because the mutation that claimed to test this used to drop a
    check AFTER production validation had run and then append the stop it
    expected, which tested the mutation rather than the production logic."""
    missing = [c for c in CHECK_IDS if c not in emitted]
    dupes = sorted({c for c in emitted if emitted.count(c) > 1})
    undeclared = sorted(set(emitted) - set(CHECK_IDS))
    stops = []
    if missing:
        stops.append(f"REGISTRY: declared checks not executed: {missing}")
    if dupes:
        stops.append(f"REGISTRY: check ids emitted more than once: {dupes}")
    if undeclared:
        stops.append(f"REGISTRY: undeclared check ids emitted: {undeclared}")
    return stops, {"declared": len(CHECK_IDS), "emitted": len(emitted),
                   "missing": missing, "duplicated": dupes, "undeclared": undeclared}


ENUM_CHECKS = ["ENUM_ORDER", "ENUM_DIGEST", "ENUM_BYTES", "ENUM_RANK_KEYS"]
SIG_CHECKS = ["SIG_COUNT", "SIG_SUMSQ", "SIG_DIMS", "SIG_LABELS_UNIQUE",
              "SIG_SEPARATION", "SIG_ORDERS"]
GEN_CHECKS = ["GEN_CANDIDATE_SET", "GEN_CLASS_PARTITION", "GEN_UNIFORM_IN_CLASS",
              "GEN_WRONG_CLASS_SAME_SIGSET", "GEN_WRONG_CLASS_IS_GALOIS_PARTNER_MAP"]


def assess(group, construction, theory, note, repro, source_pins=None):
    """The acceptance predicate. Returns (checks, stops, data).

    Every stage sits behind an exception boundary that converts a failure into
    BLOCKED checks. A stage that raises would otherwise take the whole run down
    with a traceback, and a traceback is not a verdict: the caller cannot tell a
    rejected input from a broken harness, and the registry ends up incomplete."""
    checks, stops, data = [], [], {}

    def add(cid, ok, detail):
        checks.append((cid, ok, detail))
        if not ok:
            stops.append(f"{cid}: {detail}")
        return ok

    def blocked(cids, e):
        for c in cids:
            add(c, False, f"BLOCKED: {type(e).__name__}: {e}")

    # adjudicates -----------------------------------------------------------
    for cid, name, path in (("ADJ_GROUP", "group_packet", group),
                            ("ADJ_CONSTRUCTION", "construction_packet", construction)):
        try:
            got = sha256_file(path)
        except Exception as e:                                # noqa: BLE001
            add(cid, False, f"BLOCKED: {type(e).__name__}: {e}")
            continue
        add(cid, got == PINNED[name],
            got[:16] if got == PINNED[name] else f"{got[:16]} != § 11 pin")

    # the three source records, pinned rather than merely recorded ----------
    # Without this P1 establishes that a SUPPLIED source set agrees internally,
    # not that it is the canonical revision. The pin file is deliberate state:
    # changing a source means re-pinning it on purpose.
    pins_path = source_pins or str(pathlib.Path(__file__).with_name(SOURCE_PINS))
    pinned_ok = {"theory": False, "repro": False}
    gate_pins = {}
    try:
        doc = json.loads(pathlib.Path(pins_path).read_text())
        pins, gate_pins = doc["sources"], doc.get("gate_ids", {})
        actual = {"theory": sha256_file(theory), "note": sha256_file(note),
                  "repro": sha256_file(repro)}
        bad = {k: v[:12] for k, v in actual.items() if pins.get(k) != v}
        for k in pinned_ok:
            pinned_ok[k] = pins.get(k) == actual[k]
        extra_keys = sorted(set(pins) - set(actual))
        missing_keys = sorted(set(actual) - set(pins))
        # The detail has to describe the reason the check failed. Reporting
        # "3/3 match" while failing on an unexpected pin key is a green-sounding
        # message on a red check.
        if bad:
            det = f"unpinned: {bad}"
        elif extra_keys or missing_keys:
            det = f"pin key set differs: unexpected {extra_keys}, missing {missing_keys}"
        else:
            det = "3/3 source records match their pins"
        add("SRC_PINNED", not bad and not extra_keys and not missing_keys, det)
        data["source_pins"] = pins
    except Exception as e:                                    # noqa: BLE001
        add("SRC_PINNED", False, f"BLOCKED: {type(e).__name__}: {e}")

    # enumeration -----------------------------------------------------------
    try:
        _, enum_results, _ = p1_enumerate.check(group)
        by_name = {n: (e, a) for n, e, a in enum_results}
        add("ENUM_ORDER", by_name["order 120"][0] == by_name["order 120"][1],
            by_name["order 120"][1])
        add("ENUM_DIGEST", by_name["digest"][0] == by_name["digest"][1],
            by_name["digest"][1][:16])
        add("ENUM_BYTES", by_name["byte length"][0] == by_name["byte length"][1],
            by_name["byte length"][1])
        rk = [(n, e, a) for n, e, a in enum_results if n.startswith("rank ")]
        add("ENUM_RANK_KEYS", all(e == a for _, e, a in rk),
            f"{sum(e == a for _, e, a in rk)}/{len(rk)} pinned rank keys")
    except Exception as e:                                    # noqa: BLE001
        blocked(ENUM_CHECKS, e)

    # signatures ------------------------------------------------------------
    prep = sig = None
    try:
        prep = p1_signatures.prepare(group, [theory, repro])
        sig = p1_signatures.run(group, construction, [theory, repro], prep=prep)
    except Exception as e:                                    # noqa: BLE001
        blocked(SIG_CHECKS, e)
    data["sig"] = sig
    if sig is not None:
        add("SIG_COUNT", sig["n_irreps"] == 9, sig["n_irreps"])
        add("SIG_SUMSQ", sig["sum_of_squares"] == 120, sig["sum_of_squares"])
        add("SIG_DIMS", sig["dims"] == [1, 2, 2, 3, 3, 4, 4, 5, 6], sig["dims"])
        add("SIG_LABELS_UNIQUE", sig["n_label_solutions"] == 1, sig["n_label_solutions"])
        add("SIG_SEPARATION", sig["signatures_distinct"],
            "9 distinct" if sig["signatures_distinct"] else sig["signature_collisions"])
        add("SIG_ORDERS",
            (sig["orders"]["s"], sig["orders"]["t"], sig["orders"]["st"]) == (6, 10, 4),
            sig["orders"])

    # the two source programs, executed here rather than trusted ------------
    # FAIL CLOSED. The pin is an EXECUTION PRECONDITION, not a parallel opinion.
    # Previously P1 recorded SRC_PINNED as failed and then executed the source
    # anyway: the final verdict stopped, but unpinned code had already run and
    # could have done anything before it did.
    for cid, role, path in (("SRC_RUN_THEORY", "theory", theory),
                            ("SRC_RUN_PLATFORM", "repro", repro)):
        if not pinned_ok[role]:
            add(cid, False, "BLOCKED: source is not pinned; refusing to execute it")
            continue
        ok, detail = run_source_program(path, gate_pins.get(role))
        add(cid, ok, detail)

    # values ----------------------------------------------------------------
    # A malformed source must REJECT, not raise. The extractors are strict by
    # design (they refuse an ambiguous or incomplete value table), and that
    # strictness has to surface as a failing check rather than as a traceback.
    try:
        parsed, errors, conflicts = p1_values.collect(theory, note, repro)
    except Exception as e:                                    # noqa: BLE001
        parsed, errors, conflicts = {}, [("extract", "-", "-", str(e))], []
    vrows = p1_values.adjudicate(parsed)
    data["vrows"] = vrows
    add("VAL_PARSE", not errors, f"{len(errors)} unevaluable")
    add("VAL_NO_INTRA_SOURCE_CONFLICT", not conflicts,
        conflicts[0] if conflicts else "none")
    add("VAL_REQUIRED_CARRIERS", all(not r["missing_required"] for r in vrows),
        "A and B2 carry all 9" if all(not r["missing_required"] for r in vrows)
        else [r["label"] for r in vrows if r["missing_required"]])
    add("VAL_FAMILIES",
        all(set(r["families"]) >= {"theory", "platform"} for r in vrows),
        "both families on all 9")
    add("VAL_AGREE", all(r["agree"] for r in vrows),
        f"{sum(r['agree'] for r in vrows)}/9")

    T = {r["label"]: as_phi(r["triple"]) for r in vrows if r["triple"]}

    # anchor ----------------------------------------------------------------
    if "R7" in T:
        si = (T["R7"] * T["R7"]).triple() == (1, 0, 1)
        add("ANCHOR_NOT_SELF_INVERSE", not si,
            "not self-inverse" if not si else "SELF-INVERSE, § 8 invalid-anchor branch")

    # wrong-generator-class diagnostic, all of it derived --------------------
    partner, pairs = {}, []
    try:
        if sig is None or prep is None:
            raise RuntimeError("signature stage did not complete")
        diag = p1_signatures.generator_class_diagnostic(prep, sig["s_id"], sig["t_id"])
        data["diag"] = diag
        wrong = [c for c in diag["classes"] if not c["is_packet_class"]]
        sizes = sorted(len(c["members"]) for c in diag["classes"])
        add("GEN_CANDIDATE_SET", len(diag["candidates"]) == 6,
            f"{len(diag['candidates'])} admissible t derived from the presentation")
        add("GEN_CLASS_PARTITION", sizes == [3, 3], f"class sizes {sizes}")
        add("GEN_UNIFORM_IN_CLASS", all(c["uniform_within_class"] for c in diag["classes"]),
            "every member of a class behaves identically")
        add("GEN_WRONG_CLASS_SAME_SIGSET",
            bool(wrong) and all(c["same_signature_set"] for c in wrong),
            "identical signature set, so a join by signature still succeeds")

        partner = galois_partner_map(sig["rows"])
        data["partner"] = partner
        # The map must be TOTAL and involutive before anything sorts it. A row
        # whose conjugate signature matched nothing left a None in the map, and
        # pair construction then tried to order None against a string and raised,
        # losing the whole verdict over a diagnostic.
        total = all(v is not None for v in partner.values())
        involutive = total and all(partner[partner[k]] == k for k in partner)
        add("GEN_WRONG_CLASS_IS_GALOIS_PARTNER_MAP",
            bool(wrong) and total and involutive
            and all(c["permutation"] == partner for c in wrong),
            f"moved {sorted(k for k, v in partner.items() if k != v)}" if total
            else f"partner map is not total: {sorted(k for k, v in partner.items() if v is None)}")
        if not (total and involutive):
            partner = {}
    except Exception as e:                                    # noqa: BLE001
        blocked(GEN_CHECKS, e)
        partner = {}

    # identities, parsed from the record -------------------------------------
    try:
        idents = p1_values.identities_from_note(note)
        parse_err = None
    except Exception as e:                                    # noqa: BLE001
        idents, parse_err = [], f"{type(e).__name__}: {e}"
    data["identities"] = idents
    add("IDY_PARSED", bool(idents) and not parse_err
        and all(l in T for i in idents for l, _ in i["factors"]),
        parse_err or f"{len(idents)} identities, all factors are known rows")

    # COMPLETENESS. Without this, removing both sector products leaves the two
    # ratios, every remaining identity holds, every exchange still has its owner,
    # and all checks stay green on an identity set that is half missing. The
    # protocol's § 4.4 fixes the slot count; kinds and uniqueness are checked
    # alongside it so the shape is pinned, not just the total.
    kinds = sorted(i["kind"] for i in idents)
    unique_ids = {i["id"] for i in idents}
    shape_ok = (len(idents) == IDENTITY_SLOTS
                and len(unique_ids) == IDENTITY_SLOTS
                and kinds == ["product", "product", "ratio", "ratio"])
    add("IDY_SHAPE", shape_ok,
        f"{len(idents)} identities ({len(unique_ids)} unique), kinds {kinds}; "
        f"§ 4.4 declares {IDENTITY_SLOTS} slots")

    holds = []
    for i in idents:
        if all(l in T for l, _ in i["factors"]):
            try:
                got = eval_identity(i, T)
                want = p1_values.closed_form(i["stated"])
                holds.append((i["id"], got.triple() == want.triple()))
            except Exception as e:                            # noqa: BLE001
                # An unevaluable stated form is a FAILED identity, not a crash.
                holds.append((i["id"] + f" [{type(e).__name__}]", False))
    # Every declared identity must be EVALUATED, not just every evaluated one
    # hold. Skipping an identity whose factors name an unknown row and then
    # reporting "3/3" reads green while a quarter of the set was never checked.
    add("IDY_HOLDS",
        len(holds) == len(idents) == IDENTITY_SLOTS and all(ok for _, ok in holds),
        f"{sum(ok for _, ok in holds)} hold, {len(holds)} evaluated, "
        f"{len(idents)} parsed, {IDENTITY_SLOTS} declared")

    # which identity binds which exchange, one exchange at a time ------------
    # Only identities whose factors are ALL present are evaluated. A parsed
    # factor naming an unknown row made IDY_PARSED false and then raised KeyError
    # in this loop, which lost the verdict the failing check was about to give.
    evaluable = [i for i in idents if all(l in T for l, _ in i["factors"])]
    pairs, sens = [], []
    sens_err = None
    try:
        pairs = sorted({tuple(sorted((k, v))) for k, v in partner.items()
                        if v is not None and k != v})
        if len(T) == 9:
            for pair in pairs:
                swapped = dict(T)
                swapped[pair[0]], swapped[pair[1]] = T[pair[1]], T[pair[0]]
                for i in evaluable:
                    sens.append({"pair": list(pair), "identity": i["id"],
                                 "detects": eval_identity(i, T).triple()
                                 != eval_identity(i, swapped).triple()})
    except Exception as e:                                    # noqa: BLE001
        sens_err = f"{type(e).__name__}: {e}"
        pairs, sens = [], []
    data["sensitivity"] = sens
    data["pairs"] = [list(p) for p in pairs]

    # The FULL correspondence, not "some galois-looking thing detects it".
    # Expected binder for a pair is DERIVED: the ratio identity whose factor set
    # is exactly that pair. Checking only that one binder exists and its name
    # starts with "galois" would pass a set in which neither ratio binds its own
    # pair, because names are not evidence about what a function computes.
    def det(pair, iid):
        return any(s["detects"] for s in sens
                   if s["pair"] == list(pair) and s["identity"] == iid)

    expected = {}
    for pair in pairs:
        owners = [i["id"] for i in evaluable
                  if i["kind"] == "ratio" and {l for l, _ in i["factors"]} == set(pair)]
        expected[pair] = owners
    want_ok = (sens_err is None and len(evaluable) == len(idents)
               and bool(pairs) and all(len(v) == 1 for v in expected.values()))
    detail = [] if sens_err is None else [f"BLOCKED: {sens_err}"]
    if len(evaluable) != len(idents):
        detail.append(f"{len(idents) - len(evaluable)} identity(ies) name unknown rows")
    for pair in pairs:
        binders = sorted(i["id"] for i in evaluable if det(pair, i["id"]))
        detail.append(f"{pair[0]}<->{pair[1]} bound by {binders or ['NOTHING']}")
        if binders != sorted(expected[pair]):
            want_ok = False
    add("IDY_PAIRWISE_SENSITIVITY", want_ok, "; ".join(detail))

    # REGISTRY INTEGRITY, and it is fatal. Previously a missing check produced a
    # printed note and a green verdict with exit 0, so a harness that silently
    # dropped a check still read as a full pass. Missing, duplicated and
    # undeclared ids are all stops.
    reg_stops, data["registry"] = validate_registry([c for c, _, _ in checks])
    stops += reg_stops
    data["environment"] = environment()

    return checks, stops, data


def main():
    ap = argparse.ArgumentParser()
    for k in ("group", "construction", "theory", "note", "repro"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--emit", metavar="DIR")
    ap.add_argument("--source-pins", metavar="FILE")
    a = ap.parse_args()

    checks, stops, data = assess(a.group, a.construction, a.theory, a.note, a.repro,
                                 source_pins=a.source_pins)

    # Invalidate any prior artifact BEFORE the checks can stop the run, so a
    # failure can never leave the previous green output standing.
    if a.emit:
        prior = pathlib.Path(a.emit)
        for name in ("p1_crosscheck.json", "p1_runtime.log.json"):
            f = prior / name
            if f.exists():
                f.rename(f.with_suffix(f.suffix + ".superseded"))

    print("P1 data assembly\n")
    w = max(len(c) for c, _, _ in checks)
    for cid, ok, detail in checks:
        print(f"  {'ok  ' if ok else 'FAIL'} {cid:{w}s}  {detail}")

    env = data["environment"]
    print(f"\n  environment: python {env['python']}, numpy {env['numpy']}, "
          f"mpmath {env['mpmath']}")

    print("\n  identity binding, one Galois exchange at a time:")
    for s in data.get("sensitivity", []):
        print(f"    {'binds ' if s['detects'] else 'blind '} "
              f"{s['pair'][0]}<->{s['pair'][1]}  {s['identity']}")

    print()
    if stops:
        print("  STOP. P1 does not pass to P2.")
        for s in stops:
            print(f"    - {s}")
        if a.emit:
            e = pathlib.Path(a.emit)
            for f in e.glob("*.superseded"):
                f.unlink()
            (e / "P1_FAILED.txt").write_text(
                "P1 STOPPED. No canonical artifact was produced by this run, and any\n"
                "artifact from an earlier run has been removed rather than left to be\n"
                "mistaken for current output.\n\n" + "\n".join(stops) + "\n")
            print(f"    (prior artifacts removed; {e / 'P1_FAILED.txt'} written)")
        return 1
    print("  P1 GREEN on all declared checks.")

    if a.emit:
        out = pathlib.Path(a.emit)
        out.mkdir(parents=True, exist_ok=True)
        # Publication is all-or-nothing. A failed run used to return before
        # touching --emit, leaving the PREVIOUS green artifact in place with its
        # original hash, so a reviewer or downstream step could consume a stale
        # green answer-bearing file after the current run had stopped.
        sig, vrows = data["sig"], data["vrows"]
        table = {
            "_banner": BANNER,
            "_artifact_kind": ARTIFACT_KIND,
            "_not_an_answer_packet": (
                "P1 cross-check table. NOT the § 4.4 answer packet, which is "
                "closed at seven keys and is constructed and validated at P2."),
            "adjudicates": {"group_packet": sha256_file(a.group),
                            "construction_packet": sha256_file(a.construction)},
            # ROLE AND CONTENT HASH ONLY. Absolute paths went first; basenames
            # had to follow, because the same bytes under a different filename
            # still produced a different artifact hash. Environment goes too:
            # embedding platform and library versions in a file whose selling
            # point is cross-machine byte identity was self-contradictory.
            # All of it lives in the noncanonical runtime log instead.
            "sources": {k: {"sha256": sha256_file(v)}
                        for k, v in (("A_theory", a.theory), ("B1_method_note", a.note),
                                     ("B2_reproducer", a.repro))},
            "registry": data["registry"],
            "signatures": sig["rows"],
            "generator_ids": {"s": sig["s_id"], "t": sig["t_id"], "st": sig["st_id"]},
            "orders": sig["orders"],
            "values": [{"label": r["label"], "triple": r["triple"],
                        "sources": r["sources"], "families": r["families"],
                        "as_written": r["exprs"]} for r in vrows],
            "generator_class_diagnostic": {
                "admissible_t": data["diag"]["candidates"],
                "classes": [{k: v for k, v in c.items() if k != "permutation"}
                            for c in data["diag"]["classes"]],
                "galois_partner_map": data["partner"],
                "exchanged_pairs": data["pairs"],
            },
            # The STRUCTURED identity set, not just its derived ids. P2 emits the
            # packet's `identities` key from this; without it P2 would have to
            # parse factors back out of an id string, which is the transcription
            # pattern this build removed everywhere else. Found by building P2,
            # which is the only thing that could have found it.
            "identities": [{"id": i["id"], "kind": i["kind"],
                            "factors": [[l, e] for l, e in i["factors"]],
                            "stated": i["stated"]}
                           for i in data.get("identities", [])],
            "identity_pairwise_sensitivity": data["sensitivity"],
            "checks": [{"id": c, "ok": o, "detail": str(d)} for c, o, d in checks],
        }
        p = out / "p1_crosscheck.json"
        tmp = out / "p1_crosscheck.json.partial"
        tmp.write_text(json.dumps(table, indent=2, sort_keys=True) + "\n")
        tmp.replace(p)
        for f in out.glob("*.superseded"):
            f.unlink()
        (out / "P1_FAILED.txt").unlink(missing_ok=True)
        # Machine-local detail lives here, outside the canonical artifact, so
        # the artifact stays byte-identical across machines given equal inputs.
        (out / "p1_runtime.log.json").write_text(json.dumps({
            "_note": "NONCANONICAL. Machine-local paths; not part of the artifact.",
            "paths": {k: str(pathlib.Path(v).resolve()) for k, v in
                      (("group", a.group), ("construction", a.construction),
                       ("theory", a.theory), ("note", a.note), ("repro", a.repro))},
            "basenames": {k: pathlib.Path(v).name for k, v in
                          (("group", a.group), ("construction", a.construction),
                           ("theory", a.theory), ("note", a.note), ("repro", a.repro))},
            "environment": data["environment"],
        }, indent=2, sort_keys=True) + "\n")
        print(f"\n  wrote {p}  (ANSWER-BEARING, quarantined)")
        print(f"  wrote {out / 'p1_runtime.log.json'}  (noncanonical)")
        print(f"  canonical artifact sha256: {sha256_file(p)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
