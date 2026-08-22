"""P4 self-check: the adjudicated regression set, rerun against this tree.

Runs from `build/`, mutates only throwaway COPIES of the tree, and writes
`P4_SELFCHECK.txt` beside the report. It never modifies the tree it reads.

WHAT THIS IS. Engineering evidence, in the same standing as `p3_selftest.py`
and NOT in the standing of the battery itself. `p4.py` is the formal evidence
that every P3 check can be made to reject; this file is the record that the
specific defects earlier redlines found stay fixed, and that the two bounded
claims say the same thing in the battery, the report and the review brief.

Every line here must be able to fail. Each check either mutates something and
requires the battery to notice, or requires an exact string in the battery's
real output. No line asserts a property by restating it. Two lines in an
earlier version did not meet that standard and were replaced: one injected an
unused assignment instead of the substring-matching regression it claimed to
test, and one hollowed a P3 check with an anchor that matched only the
single-line spelling of `add(...)`, so it silently disabled an exception branch
and reported a false survivor.

ANCHORS LOCATE CODE; THEY MUST NOT LEARN FROM IT. A mutation has to find its
target somehow, so every anchor here is a literal string matched against the
file. What makes that safe is the failure mode: `patch` requires the anchor to
appear exactly once and raises otherwise, so a refactor that moves the target
ERRORS rather than quietly testing nothing. The two inventories that are read
from the code under test, the P3 check list and `p4.BOUNDS`, are pinned to
counts on this side. Verified by refactoring the tree four ways: renaming a
targeted function and rewrapping the loader both raise here; renaming a P3 check
and renaming a BOUNDS key both fail checks. None of the four passed silently.
"""
import sys

# Before any local import. This file reads `p4.BOUNDS`, and importing p4
# left `build/__pycache__/p4.cpython-*.pyc` in the tree it is auditing.
sys.dont_write_bytecode = True

import ast
import hashlib
import pathlib
import re
import shutil
import subprocess
import tempfile

TREE = pathlib.Path(__file__).resolve().parent.parent
ARGS = ["out/p1_crosscheck.json",
        "out2/m8_8_answer_packet.CANDIDATE.json",
        "out2/m8_8_prelock_hash_record.CANDIDATE.json",
        "bundle/group__m8_5a_packet.json",
        "bundle/construction__m8_8_construction_packet.json"]
RESULTS = []


def record(ok, title, detail):
    RESULTS.append((ok, title, detail))
    print(f"  {'ok  ' if ok else 'FAIL'} {title}\n         {detail}", flush=True)


def scratch_copy():
    d = pathlib.Path(tempfile.mkdtemp(prefix="p4regress-"))
    shutil.copytree(TREE / "build", d / "build",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    return d


TREE_ROOT = TREE          # the real tree, for the writes-nothing check only


def battery(root, timeout=300):
    r = subprocess.run([sys.executable, "p4.py", *ARGS], cwd=root / "build",
                       capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout + r.stderr


def line_for(out, token):
    for ln in out.splitlines():
        if token in ln:
            return " ".join(ln.split())[:112]
    return "(no line found)"


def caught_by(out, *names):
    """True when the battery reddened AND one of `names` is what reported it.

    `rc != 0` alone is not evidence: a mutation that raises inside the gate also
    exits non-zero, and this file has already thrown out one line for scoring a
    `NameError` as a catch. A regression that cannot tell "the check I named
    rejected it" from "something exploded" is testing the interpreter.

    TWO GUARDS, COVERING DIFFERENT SHAPES, and both are load-bearing:

      the traceback guard   covers a CRASH, which exits non-zero and produces no
                            verdict line at all. It works only because
                            `battery()` returns stdout AND stderr: a raise inside
                            `p3.gate` puts the traceback on stderr with nothing
                            on stdout, so a guard reading stdout alone would miss
                            it entirely.
      the named-reporter    covers a red FOR THE WRONG REASON, which does not
      guard                 crash and would satisfy the first guard completely.
                            This is the stronger of the two.
    """
    if "Traceback" in out or "NameError" in out:
        return False
    return any(f"FAIL {n}" in out or f"MISS  {n}" in out for n in names)


def patch(root, rel, old, new, all_sites=False):
    p = root / rel
    t = p.read_text()
    n = t.count(old)
    assert n == 1 or (all_sites and n >= 1), \
        f"{rel}: anchor appears {n}x: {old[:70]!r}"
    p.write_text(t.replace(old, new))
    ast.parse((root / rel).read_text())          # never ship a syntax-error pass
    return n


def mutate(rel, old, new, all_sites=False):
    """Apply one mutation in a copy and return (exit code, output, sites)."""
    d = scratch_copy()
    try:
        n = patch(d, rel, old, new, all_sites)
        rc, out = battery(d)
        return rc, out, n
    finally:
        shutil.rmtree(d, ignore_errors=True)


REF = "build/p4_reference.py"

# 1. baseline ---------------------------------------------------------------
base = scratch_copy()
rc, BASE_OUT = battery(base)
# The counts are FROZEN here, not read from the battery. Adding four content
# cases for the fifth review made this line fail, which is the point: a baseline
# that tracked the battery's own totals would have said nothing.
record(rc == 0 and "P4 GREEN. Content: 28/28, coverage 12/12. Ingestion: 8/8." in BASE_OUT,
       "baseline 28/28 content, 12/12 P3 checks, 8/8 ingestion",
       line_for(BASE_OUT, "P4 GREEN"))

# 2. hollow every P3 check --------------------------------------------------
P3SRC = (TREE / "build/p3.py").read_text()
CHECKS = re.findall(r'add\(\s*"([A-Z][A-Z0-9_]*)"\s*,', P3SRC)
CHECKS = list(dict.fromkeys(CHECKS))
survivors, sites = [], {}
for name in CHECKS:
    # `add("X", cond, ...)` is written both on one line and wrapped after the
    # comma. An anchor that only matched one shape silently hollowed the
    # exception branch instead of the check, and reported a false survivor.
    pat = re.compile(r'add\(\s*"%s"\s*,\s*' % re.escape(name))
    d = scratch_copy()
    try:
        src = (d / "build/p3.py").read_text()
        src, n = pat.subn(f'add("{name}", True or ', src)
        assert n >= 1, f"no recording site for {name}"
        (d / "build/p3.py").write_text(src)
        ast.parse(src)
        rc2, out2 = battery(d)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    sites[name] = n
    if rc2 == 0 and "P4 GREEN" in out2:
        survivors.append(name)
# 19 is FROZEN here, like the check count and the bound count. A wrong site
# total is the signal that the sweep itself is broken, so reporting it without
# pinning it was the one inventory this file only printed.
record(len(CHECKS) == 12 and sum(sites.values()) == 19 and not survivors,
       f"all {len(CHECKS)} P3 checks hollow-tested, none survives",
       f"{len(CHECKS) - len(survivors)}/{len(CHECKS)} caught at "
       f"{sum(sites.values())} recording sites; "
       f"survivors: {', '.join(survivors) if survivors else 'none'}")

# 3. F1 address the committed rows by the join, discarding the map ----------
rc3, out3, _ = mutate(
    REF, 'if selected(committed_rows[corr[i]]["value"]) != tuple(rows[i]["value"])]',
    'if selected(committed_rows[(corr[i], join[signature_key(rows[i])][1])[1]]'
    '["value"]) != tuple(rows[i]["value"])]')
record(rc3 != 0 and "FAIL SYNTHETIC_NONIDENTITY_MAP" in out3,
       "F1 index-and-discard then align through the join: caught",
       line_for(out3, "SYNTHETIC_NONIDENTITY_MAP"))

# 4. F2 parse before the hash is checked ------------------------------------
rc4, out4, _ = mutate(
    REF, "def load_hash_verified(path, expected_sha256, parse=json.loads):",
    "def load_hash_verified(path, expected_sha256, parse=json.loads):\n"
    "    json.JSONDecoder().decode(open(path, 'rb').read().decode())")
record(rc4 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in out4,
       "F2 parse before the hash check, entering by scan_once: caught",
       line_for(out4, "WRONG_HASH_BEFORE_PARSE"))

# 5. F3 revert declaration validation to substring matching ----------------
rc5, out5, _ = mutate(REF, "        if block[leaf] != expected:",
                   "        if expected not in block[leaf]:")
record(rc5 != 0 and "FAIL CONVENTION_MAP_VALIDATED" in out5,
       "F3 substring validation, so a leaf may keep the wording and negate it: caught",
       line_for(out5, "CONVENTION_MAP_VALIDATED"))

# 6. F4 constant selector, and the inverse branch disabled ------------------
rc6a, out6a, _ = mutate(
    REF, "    x, r = tuple(anchor_committed), tuple(anchor_reference)",
    "    return \"committed\"\n"
    "    x, r = tuple(anchor_committed), tuple(anchor_reference)")
rc6b, out6b, _ = mutate(REF, '    if inverted:\n        return "inverse"',
                     '    if inverted:\n        return "committed"')
record(caught_by(out6a, "ORIENTATION_SELECTION_EXERCISED")
       and caught_by(out6b, "ORIENTATION_SELECTION_EXERCISED"),
       "F4 constant selector and disabled inversion: both caught",
       f"constant selector {'caught' if rc6a else 'SURVIVED'}, "
       f"{line_for(out6a, 'ORIENTATION_SELECTION')[:40]}...; "
       f"inversion disabled {'caught' if rc6b else 'SURVIVED'}")

# 7. F5 stage-two identity recomputation removed ---------------------------
rc7, out7, _ = mutate(REF, '    for slot, ident in enumerate(packet["identities"]):',
                   '    for slot, ident in []:')
record(caught_by(out7, "DOWNSTREAM_CELL_MUTATION"),
       "F5 § 5.3 stage-two identity recomputation removed: caught",
       line_for(out7, "DOWNSTREAM_CELL_MUTATION"))

# --- the four shortcutes the external redline found, as standing regressions ---
rcE1, outE1, _ = mutate(
    REF, 'if selected(committed_rows[corr[i]]["value"]) != tuple(rows[i]["value"])]',
    'if selected(next(r for j, r in enumerate(committed_rows)\n'
    '                      if j == join[signature_key(rows[i])][1])["value"])\n'
    '                 != tuple(rows[i]["value"])]')
record(rcE1 != 0 and "FAIL SYNTHETIC_NONIDENTITY_MAP" in outE1,
       "E1 committed rows reached by ITERATION, never by index: caught",
       line_for(outE1, "SYNTHETIC_NONIDENTITY_MAP"))

LOADER = "    return parse(p4_verify.verified_bytes(path, expected_sha256))"
VER = "build/p4_verify.py"

rcE2, outE2, _ = mutate(
    REF, LOADER,
    '    json._default_decoder.scan_once(pathlib.Path(path).read_bytes()'
    '.decode("ascii"), 0)\n' + LOADER)
record(rcE2 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outE2,
       "E2 parse through the module's pre-bound default decoder: caught",
       line_for(outE2, "WRONG_HASH_BEFORE_PARSE"))

rcE3, outE3, _ = mutate(
    REF, LOADER,
    '    return parse(p4_verify.verified_bytes(path, expected_sha256)'
    '.decode("ascii") + " ")')
record(rcE3 != 0 and "FAIL SAME_PARSER_PATH" in outE3,
       "E3 verify the canonical bytes, parse a re-rendering: caught",
       line_for(outE3, "SAME_PARSER_PATH"))

# The sixth failure of the parse-order control, and the three limbs that replaced
# measurement with structure. Each limb has to be able to fail on its own.
rcE5, outE5, _ = mutate(
    REF, LOADER,
    '    _E = json.JSONDecoder()\n'
    '    _E.scan_once(pathlib.Path(path).read_bytes().decode("ascii"), 0)\n' + LOADER)
record(rcE5 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outE5,
       "E5 a SEPARATELY pre-bound decoder parsing before the hash: caught",
       line_for(outE5, "WRONG_HASH_BEFORE_PARSE"))

rcE6, outE6, _ = mutate(VER, "import hashlib\nimport pathlib",
                        "import hashlib\nimport os\nimport pathlib")
record(rcE6 != 0 and "AND SHOULD NOT" in outE6,
       "E6 ISOLATION limb, imports: an extra import is caught",
       line_for(outE6, "WRONG_HASH_BEFORE_PARSE"))

# Imports were never the whole property. Built-ins need no import, and `eval`
# parses the packet with nothing imported at all.
rcE6b, outE6b, _ = mutate(
    VER, "    got = hashlib.sha256(raw).hexdigest()",
    '    eval(raw.decode("ascii").replace("true", "True")\n'
    '         .replace("false", "False").replace("null", "None"),\n'
    '         {"__builtins__": {}})\n'
    "    got = hashlib.sha256(raw).hexdigest()")
record(rcE6b != 0 and "NO IT DOES NOT" in outE6b,
       "E6b ISOLATION limb, body: a built-in parse before the hash is caught",
       line_for(outE6b, "WRONG_HASH_BEFORE_PARSE"))

rcE7, outE7, _ = mutate(VER, "    return raw", "    return raw[:-1]")
record(caught_by(outE7, "WRONG_HASH_BEFORE_PARSE"),
       "E7 EXECUTION limb: verification returns other than the bytes: caught",
       line_for(outE7, "WRONG_HASH_BEFORE_PARSE"))

rcE8, outE8, _ = mutate(REF, LOADER, "    _ = 1\n" + LOADER)
record(caught_by(outE8, "WRONG_HASH_BEFORE_PARSE"),
       "E8 COMPOSITION limb: a statement before the verification call: caught",
       line_for(outE8, "WRONG_HASH_BEFORE_PARSE"))

rcE9, outE9, _ = mutate("build/p4.py", "sys.meta_path.insert(0, _RefuseJson())",
                        "pass  # import guard disabled")
record(rcE9 != 0 and "proves nothing" in outE9,
       "E9 POSITIVE CONTROL: the no-JSON import guard disabled, probe discards itself",
       line_for(outE9, "WRONG_HASH_BEFORE_PARSE"))

rcE4, outE4, _ = mutate(
    REF, "    owning = [r for r in packet[\"rows\"] if r.get(\"evidentiary_class\") == ANCHOR_CLASS]",
    "    owning = [r for r in packet[\"rows\"]\n"
    "              if signature_key(r) == declared]  # believe the claimant")
record(caught_by(outE4, "ORIENTATION_SELECTION_EXERCISED",
                 "SYNTHETIC_NONIDENTITY_MAP"),
       "E4 anchor validation reverted to believing the packet: caught",
       line_for(outE4, "ORIENTATION_SELECTION") if "ORIENTATION_SELECTION"
       in outE4 else line_for(outE4, "FAIL"))

# The third external redline: an argument is an expression, a map check that
# never had to fire, and an identity claim about an object nobody read.
rcG1, outG1, _ = mutate(
    REF, LOADER,
    "    return parse(\n        p4_verify.verified_bytes(\n"
    "            (json.JSONDecoder().scan_once(pathlib.Path(path).read_text(), 0),\n"
    "             path)[1],\n            expected_sha256,\n        )\n    )")
record(rcG1 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outG1,
       "G1 a parse placed in an ARGUMENT of the verification call: caught",
       line_for(outG1, "WRONG_HASH_BEFORE_PARSE"))

rcG2, outG2, _ = mutate(REF, "        check_indexing_map(packet, committed_rows, join)",
                        "        pass  # check_indexing_map removed")
record(rcG2 != 0 and "FAIL SYNTHETIC_NONIDENTITY_MAP" in outG2,
       "G2 check_indexing_map removed outright: caught",
       line_for(outG2, "SYNTHETIC_NONIDENTITY_MAP"))

rcG3, outG3, _ = mutate(REF, '    rows = packet["rows"]\n    try:',
                        '    packet = dict(packet)\n    rows = packet["rows"]\n    try:')
record(rcG3 != 0 and "FAIL SAME_PARSER_PATH" in outG3,
       "G3 compare() rebinds the packet to a shallow copy: caught",
       line_for(outG3, "SAME_PARSER_PATH"))

BRANCHES = [
    ("join-agreement", "        if joined != dst:", "        if False and joined != dst:"),
    ("duplicate-destination", "        if dst in seen:", "        if False and dst in seen:"),
    ("coverage", "    if len(seen) != len(committed_rows):",
     "    if False and len(seen) != len(committed_rows):"),
    ("range", "        if not (0 <= src < len(rows) and 0 <= dst < len(committed_rows)):",
     "        if False and not (0 <= src < len(rows) and 0 <= dst < len(committed_rows)):"),
]
survived = []
for label, old_b, new_b in BRANCHES:
    rcB, outB, _ = mutate(REF, old_b, new_b)
    if rcB == 0 or "FAIL SYNTHETIC_NONIDENTITY_MAP" not in outB:
        survived.append(label)
record(not survived,
       "G4-G7 four check_indexing_map branches hollowed; the header, type and "
       "destination arms are covered at J1, J2 and N4, not here",
       f"{len(BRANCHES) - len(survived)}/{len(BRANCHES)} branches caught; "
       f"survivors: {', '.join(survived) if survived else 'none'}")

# The fourth external redline: qualify the receiver, and prove causal use rather
# than a warm-up touch.
# `evil` must be a WORKING delegate, or the mutation crashes at runtime and
# tests nothing about the structural check it is aimed at. A first version left
# the name undefined and reddened P4 with a NameError, which is a pass for the
# wrong reason: the whole point is that the AST test must reject the receiver
# even when the code runs fine.
_DELEGATE = ("import json\nimport types\n\nimport p4_verify\n\n"
         "def _delegating_verified_bytes(path, expected_sha256):\n"
         "    json.JSONDecoder().scan_once(pathlib.Path(path).read_text(), 0)\n"
         "    return p4_verify.verified_bytes(path, expected_sha256)\n\n"
         "other = types.SimpleNamespace(verified_bytes=_delegating_verified_bytes)\n")
d_ = scratch_copy()
try:
    patch(d_, REF, "import json\n\nimport p4_verify\n", _DELEGATE)
    patch(d_, REF, LOADER, "    return parse(other.verified_bytes(path, expected_sha256))")
    rcH1, outH1 = battery(d_)
finally:
    shutil.rmtree(d_, ignore_errors=True)
record(rcH1 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outH1
       and "NameError" not in outH1,
       "H1 a DIFFERENT receiver, delegating to the real verifier: caught",
       line_for(outH1, "WRONG_HASH_BEFORE_PARSE"))

BODY = '    rows = packet["rows"]\n    try:'
COPIES = [
    ("touch-then-copy",
     '    for _k in ("rows", "identities", "indexing_map", "convention_map",\n'
     '               "adjudicates"):\n        packet[_k]\n    packet = dict(packet)\n'),
    ("plain shallow copy", "    packet = dict(packet)\n"),
    ("**unpacking", "    packet = {**packet}\n"),
]
survived_copy = []
for label, prefix in COPIES:
    rcC, outC, _ = mutate(REF, BODY, prefix + BODY)
    if rcC == 0 or "FAIL SAME_PARSER_PATH" not in outC:
        survived_copy.append(label)
record(not survived_copy,
       "H2 three ways of comparing a COPY instead of the parsed object: all caught",
       f"{len(COPIES) - len(survived_copy)}/{len(COPIES)} caught; "
       f"survivors: {', '.join(survived_copy) if survived_copy else 'none'}")

# The negative control this whole mechanism needs: a consumer that legitimately
# binds its data THROUGH the parsed object must still pass, or the check is
# rejecting correct behaviour rather than incorrect behaviour.
rcH3, outH3, _ = mutate(
    REF, BODY, '    _r = packet["rows"]\n    _i = packet["identities"]\n' + BODY)
record(rcH3 == 0 and "ok   SAME_PARSER_PATH" in outH3,
       "H3 NEGATIVE CONTROL: binding through the parsed object stays green",
       line_for(outH3, "SAME_PARSER_PATH"))

# The fifth external redline: header declarations the consumer assumed, and a
# check that named a source it never read.
rcJ1, outJ1, _ = mutate(
    REF, "    for key, expected in INDEXING_MAP_HEADER.items():",
    "    for key, expected in {}.items():")
record(caught_by(outJ1, "SYNTHETIC_NONIDENTITY_MAP"),
       "J1 indexing-map header validation removed: caught",
       line_for(outJ1, "SYNTHETIC_NONIDENTITY_MAP"))

rcJ2, outJ2, _ = mutate(
    REF, "            if isinstance(v, bool) or not isinstance(v, (int, _IntLike)):",
    "            if False:")
record(caught_by(outJ2, "SYNTHETIC_NONIDENTITY_MAP"),
       "J2 position type-check removed, so JSON false indexes as 0: caught",
       line_for(outJ2, "SYNTHETIC_NONIDENTITY_MAP"))

# The fifth external review: a bound asserted over a stage the fixture never
# reached, a structural argument about text that nothing tied to the running
# objects, and two P3 branches no content case touched.
P3F = "build/p3.py"

rcK1, outK1, _ = mutate(
    REF, 'committed_rows[corr[fidx]]["value"]',
    'committed_rows[join[signature_key(rows[fidx])][1]]["value"]')
record(rcK1 != 0 and "FAIL SYNTHETIC_NONIDENTITY_MAP" in outK1,
       "K1 § 5.3 identity stage aligns through the join, not the map: caught",
       line_for(outK1, "SYNTHETIC_NONIDENTITY_MAP"))

DELEGATE = """import pathlib
import json as _j
import types as _t

import p4_verify as _real_verify

_EARLY = _j.JSONDecoder().scan_once


def _delegated(path, expected_sha256):
    _EARLY(pathlib.Path(path).read_text(), 0)
    return _real_verify.verified_bytes(path, expected_sha256)


p4_verify = _t.SimpleNamespace(verified_bytes=_delegated,
                               HashMismatch=_real_verify.HashMismatch,
                               __file__=_real_verify.__file__)"""
rcK2, outK2, _ = mutate(REF, "import p4_verify\nimport pathlib", DELEGATE)
record(rcK2 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outK2,
       "K2 a DELEGATE under the name `p4_verify`, real `__file__`: caught",
       line_for(outK2, "WRONG_HASH_BEFORE_PARSE"))

DECOR = """def _pre_parse(fn):
    import json as _j
    _early = _j.JSONDecoder().scan_once

    def wrapper(path, expected_sha256, parse=json.loads):
        _early(pathlib.Path(path).read_text(), 0)
        return fn(path, expected_sha256, parse)
    return wrapper


@_pre_parse
def load_hash_verified("""
rcK3, outK3, _ = mutate(REF, "def load_hash_verified(", DECOR)
record(rcK3 != 0 and "FAIL WRONG_HASH_BEFORE_PARSE" in outK3,
       "K3 a DECORATOR, loader body byte-identical, runs first: caught",
       line_for(outK3, "WRONG_HASH_BEFORE_PARSE"))

BRANCHES_P3 = [
    ("adjudicates construction, declared",
     '        if adj.get("construction_packet_sha256") != p1_construction:',
     '        if False and adj.get("construction_packet_sha256") != p1_construction:',
     "wrong_adjudicates_construction"),
    ("adjudicates construction, supplied",
     "        if sha256_bytes(construction_bytes) != p1_construction:",
     "        if False and sha256_bytes(construction_bytes) != p1_construction:",
     "substituted_construction_packet"),
    ("both-families",
     '            if not {"theory", "platform"} <= fams:',
     '            if False and not {"theory", "platform"} <= fams:',
     "one_family_missing"),
    ("families cross-check",
     '            if fams != set(rec.get("families") or []):',
     '            if False and fams != set(rec.get("families") or []):',
     "families_field_desynchronised"),
]
survived_p3 = []
for label, old_b, new_b, want_case in BRANCHES_P3:
    rcB, outB, _ = mutate(P3F, old_b, new_b)
    # The NAMED case must be the one that goes MISS. Requiring only a red let a
    # `NameError` inside the gate score as a catch, which is this file's own
    # recorded mistake at a new site.
    if not caught_by(outB, want_case):
        survived_p3.append(f"{label} ({'no red' if rcB == 0 else 'red, but not by ' + want_case})")
record(not survived_p3,
       "K4-K7 four ADJUDICATES and FAMILIES branches, each caught BY ITS OWN CASE",
       f"{len(BRANCHES_P3) - len(survived_p3)}/{len(BRANCHES_P3)} caught by the "
       f"case that names them; survivors: "
       f"{'; '.join(survived_p3) if survived_p3 else 'none'}")

# The verifier's permitted-call set is a PROHIBITION now, not a permission:
# deleting the dead `bytes` entry shape let `isinstance` come out of it, so an
# `isinstance` reappearing in that file has to redden the limb.
rcL1, outL1, _ = mutate(
    VER, "    raw = pathlib.Path(path).read_bytes()",
    "    raw = pathlib.Path(path).read_bytes() if not isinstance(path, bytes) else path")
record(caught_by(outL1, "WRONG_HASH_BEFORE_PARSE"),
       "L1 an `isinstance` reintroduced into the verifier: caught",
       line_for(outL1, "WRONG_HASH_BEFORE_PARSE"))

# The seventh review: a file is not a function, and three controls whose success
# condition was satisfied by doing nothing at all.
SAME_FILE_REBIND = """    return raw


def _verified_bytes_impl(path, expected_sha256):
    raw = pathlib.Path(path).read_bytes()
    _verified_bytes_impl.parsed_before_hash_check = eval(
        raw.decode("ascii").replace("true", "True")
        .replace("false", "False").replace("null", "None"),
        {"__builtins__": {}})
    got = hashlib.sha256(raw).hexdigest()
    if got != expected_sha256:
        raise HashMismatch("refused before parsing")
    return raw


verified_bytes = _verified_bytes_impl
"""
rcM1, outM1, _ = mutate(VER, "    return raw\n", SAME_FILE_REBIND)
record(caught_by(outM1, "WRONG_HASH_BEFORE_PARSE"),
       "M1 a SECOND def in the same file, rebound to the name: caught",
       line_for(outM1, "WRONG_HASH_BEFORE_PARSE"))

# Three vacuity mutations. Each empties an inventory whose control asserted only
# that NOTHING FAILED, which trying nothing satisfies perfectly.
VACUOUS = [
    ("ORIENTATION_SELECTION_EXERCISED, alternate anchors",
     '    for row in pkt["rows"]:', "    for row in []:",
     "ORIENTATION_SELECTION_EXERCISED"),
    ("SYNTHETIC_NONIDENTITY_MAP, malformed-map branches",
     "    for label, fn, want in BRANCHES:", "    for label, fn, want in []:",
     "SYNTHETIC_NONIDENTITY_MAP"),
    ("SAME_PARSER_PATH, required-key inventory",
     '    NEEDED = {"rows", "identities", "indexing_map", "convention_map", "adjudicates"}',
     "    NEEDED = set()", "SAME_PARSER_PATH"),
]
vacuous_survived = []
for label, old_v, new_v, want in VACUOUS:
    rcV, outV, _ = mutate("build/p4.py", old_v, new_v)
    if not caught_by(outV, want):
        vacuous_survived.append(label)
record(not vacuous_survived,
       "M2-M4 three inventories emptied; each control must count, not just not-fail",
       f"{len(VACUOUS) - len(vacuous_survived)}/{len(VACUOUS)} caught; "
       f"survivors: {'; '.join(vacuous_survived) if vacuous_survived else 'none'}")

# The eighth review: the gate's rejection channel, the globals a pinned body
# calls into, a negative that accepted any reason, an unreached arm, a scanner
# with no positive control, and a schema version nothing read.
rcN1, outN1, _ = mutate(
    P3F, "        if not ok:\n            stops.append",
    "        if False:\n            stops.append")
record(rcN1 != 0, "N1 p3's rejection channel disabled; per-check booleans intact: caught",
       line_for(outN1, "Content:") or "battery reddened")

SHIM = """HashMismatch = p4_verify.HashMismatch

_EARLY_DECODER = json.JSONDecoder()


class _ParsingPathlib:
    class Path:
        def __init__(self, p):
            self._p = pathlib.Path(p)

        def read_bytes(self):
            raw = self._p.read_bytes()
            _EARLY_DECODER.scan_once(raw.decode("ascii"), 0)
            return raw


p4_verify.pathlib = _ParsingPathlib"""
rcN2, outN2, _ = mutate(REF, "HashMismatch = p4_verify.HashMismatch", SHIM)
record(caught_by(outN2, "WRONG_HASH_BEFORE_PARSE"),
       "N2 `p4_verify.pathlib` swapped for a parsing shim: caught",
       line_for(outN2, "WRONG_HASH_BEFORE_PARSE"))

rcN3, outN3, _ = mutate(REF, '        elif corr_mode == "join":', "        elif False:")
record(caught_by(outN3, "SYNTHETIC_NONIDENTITY_MAP"),
       "N3 join mode raises instead of disagreeing; the negative must name it: caught",
       line_for(outN3, "SYNTHETIC_NONIDENTITY_MAP"))

rcN4, outN4, _ = mutate(
    REF, '        for field in ("source_position", "destination_position"):',
    '        for field in ("source_position",):')
record(caught_by(outN4, "SYNTHETIC_NONIDENTITY_MAP"),
       "N4 the integer guard narrowed to the source arm: caught",
       line_for(outN4, "SYNTHETIC_NONIDENTITY_MAP"))

rcN5, outN5, _ = mutate(
    "build/p4.py", "def _triple_shaped_literals(src):",
    "def _triple_shaped_literals(src):\n    return []\n\n\ndef _unused_tsl(src):")
record(caught_by(outN5, "NO_MANUAL_TRANSCRIPTION"),
       "N5 the transcription scanner made blind; positive control fires: caught",
       line_for(outN5, "NO_MANUAL_TRANSCRIPTION"))

rcN6, outN6, _ = mutate(
    REF, '    if packet.get("format_version") != FORMAT_VERSION:', "    if False:")
record(caught_by(outN6, "CONVENTION_MAP_VALIDATED"),
       "N6 the fixture's format_version check removed: caught",
       line_for(outN6, "CONVENTION_MAP_VALIDATED"))

# 8. both halves of the adjudicates binding --------------------------------
rc8, out8, _ = mutate(
    REF, "def check_adjudicates(packet, group_bytes, construction_bytes):",
    "def check_adjudicates(packet, group_bytes, construction_bytes):\n    return")
adj = line_for(BASE_OUT, "ADJUDICATES_BINDING")
record(rc8 != 0 and "ok   ADJUDICATES_BINDING" in BASE_OUT and "construction" in adj,
       "both adjudicates halves; removing the binding is caught", adj)

# 9-10. controls that report their own internals ---------------------------
cm = line_for(BASE_OUT, "CONVENTION_MAP_VALIDATED")
# 11, down from 19. An independent read adjudicated the Q(phi) encoding out of
# the convention map, so its four leaves left the relied-on inventory and took
# eight mutations with them. Not coverage lost: coverage with nothing left to
# point at. The frozen number caught the drop on the run that caused it, which is
# what freezing it here is for.
record("ok   CONVENTION_MAP_VALIDATED" in BASE_OUT and "11 mutations" in cm,
       "11 convention mutations, inventory frozen control-side", cm)
record(all(w in BASE_OUT for w in ("INVERSE orientation", "INVALID_ANCHOR", "DISAGREEMENT")),
       "all three § 5.4 outcomes exercised",
       line_for(BASE_OUT, "ORIENTATION_SELECTION_EXERCISED"))


# 11. the bounded claims, word for word, in the battery's own output -------
def norm(t):
    return " ".join(re.sub(r'(?m)^\s*>\s?', ' ', t).replace('`', '').split())


# The bounds are read FROM THE CODE, not restated here. A self-check holding its
# own copy is a fifth hand-maintained copy of the thing it is checking for drift.
sys.path.insert(0, str(TREE / "build"))
import p4 as _P4                                                  # noqa: E402
BOUNDS = dict(_P4.BOUNDS)
assert len(BOUNDS) == 2, f"expected two bounded controls, found {sorted(BOUNDS)}"
nb = norm(BASE_OUT)
cnt = {cid: nb.count(norm(b)) for cid, b in BOUNDS.items()}
record(all(v >= 2 for v in cnt.values()),
       "each bound is printed VERBATIM, with its own result and in the scope block",
       "exact-string counts in stdout: " + ", ".join(f"{k} {v}x" for k, v in cnt.items()))

# 12. the same sentences in report, docstring and brief; no wider claim ----
def all_docstrings(path):
    """EVERY docstring in a module, not just the module's own.

    The self-check used to read only the module docstring. `control_2`'s own
    docstring carried an exclusion the printed bound did not state, and this
    check reported the bounds in agreement anyway: it was looking at the one
    place the discrepancy was not. An external redline found that; the fix is to
    read all of them."""
    out = []
    for node in ast.walk(ast.parse(pathlib.Path(path).read_text())):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            d = ast.get_docstring(node)
            if d:
                out.append(d)
    return out


DOC_SOURCES = ("build/p4.py", "build/p4_reference.py", "build/p4_verify.py")
DOCSTRINGS = {f"{f.split('/')[-1]}#{i}": norm(d)
              for f in DOC_SOURCES for i, d in enumerate(all_docstrings(TREE / f))}
# Docstrings are deliberately NOT here: after the BOUNDS refactor they refer to
# the bounds rather than restating them, and the paraphrase check below is what
# holds them to that. These are the two places that legitimately QUOTE a bound.
# BASE_OUT belongs here. The recorded detail said this compared "battery output
# + report" while `WHERE` held only the report, so a `print` of an overclaim in
# the battery itself went unseen. A check that names a source it never read is
# the plainest form of the defect this whole file exists to catch.
WHERE = {"battery output": nb, "report": norm((TREE / "P4_REPORT.md").read_text())}
# The review brief is written OUTSIDE the bundle on purpose, so it is present
# author-side and absent from a reviewer's extraction. Say which sources were
# actually compared rather than quietly comparing fewer.
brief = TREE / "P4_REDLINE_PROMPT.md"
if brief.exists():
    WHERE["brief"] = norm(brief.read_text())
short = [f"{k}: {cid}" for k, v in WHERE.items() for cid, b in BOUNDS.items()
         if norm(b) not in v]
WIDER = ["§ 4.4's eight", "eight § 4.4", "all eight § 4.4", "§ 4.4 requires eight",
         "§ 4.4 assigns eight controls", "all § 4.4 controls"]
# A phrase may appear as a QUOTATION of wording that is then rejected. Require
# every occurrence to be followed closely by the rejection, and count the ones
# that stand on their own as claims.
REJECTS = r'was wrong|would not be|narrowed|an earlier wording|no longer'
hits, quoted = [], []
for w in WIDER:
    for where, v in WHERE.items():
        for m in re.finditer(re.escape(w), v):
            # A rejection reads either way round: "X was wrong" and "narrowed
            # from X" both disown the phrase, so look on both sides of it.
            window = v[max(0, m.start() - 70):m.end() + 70]
            (quoted if re.search(REJECTS, window) else hits).append(f"{where}: {w}")
record(not short and not hits,
       "bounds agree wherever they are stated; no text claims § 4.4 assigns eight",
       f"compared: {', '.join(WHERE)}"
       + ("" if "brief" in WHERE else " (the review brief ships outside the bundle "
          "and is not in this tree)")
       + f"; missing the sentences: {short or 'none'}; "
       f"wider-claim phrases found among {len(WIDER)} searched: {hits or 'none'}"
       + (f"; cleared as quotations that the surrounding text disowns: "
          f"{quoted}" if quoted else ""))

# 13. no docstring anywhere states a BOUND the battery does not print --------
# A docstring that begins to restate a bound must restate it EXACTLY. The
# defect this replaces was a docstring carrying an exclusion the printed bound
# did not, so a near-copy is the shape to catch, not an unrelated sentence.
drift = []
for where, d in DOCSTRINGS.items():
    for cid, b in BOUNDS.items():
        head = norm(b)[:44]
        if head in d and norm(b) not in d:
            drift.append(f"{where} starts {cid}'s bound and diverges")
record(not drift,
       "no docstring paraphrases a bound; each is held once in p4.BOUNDS",
       f"{len(DOCSTRINGS)} docstrings across "
       f"{', '.join(f.split('/')[-1] for f in DOC_SOURCES)} scanned against "
       f"{len(BOUNDS)} bounds; divergent restatements: {drift or 'none'}")

# 14. the brief's quoted scope block is what the battery actually prints ---
BRIEF = TREE / "P4_REDLINE_PROMPT.md"
if BRIEF.exists():
    q = norm(BRIEF.read_text())
    a, b = q.find("Scope: SIX of the eight"), q.find("here claims that harness is correct.")
    segs = [s.strip() for s in q[a:b + 36].split("[...]") if s.strip()]
    absent = [s[:44] + "..." for s in segs if s not in nb]
    record(len(segs) == 2 and not absent,
           "the brief quotes the battery's scope block, not a paraphrase of it",
           f"{len(segs)} quoted segment(s) either side of the elision; "
           f"absent from the battery's output: {absent or 'none'}")
else:
    record(True, "the brief quotes the battery's scope block, not a paraphrase of it",
           "SKIPPED: the review brief ships outside the bundle and is not in this tree")

# 15. the battery writes nothing -------------------------------------------
# __pycache__ is INCLUDED. Excluding it made this check unable to fail on the
# one thing the claim is actually about: the battery gained a subprocess, the
# subprocess imported a module, and the interpreter wrote a .pyc that this
# snapshot could not see.
def snap(root):
    """Every path AND the content hash of every file.

    Paths alone could not see an overwrite: rewriting an existing file left the
    listing identical, so a battery that appended to a source in place would have
    passed this check. "Writes nothing" is a claim about bytes, not about the
    directory listing."""
    out = {}
    for q in sorted(root.rglob("*")):
        out[str(q.relative_to(root))] = (
            "d" if q.is_dir() else hashlib.sha256(q.read_bytes()).hexdigest())
    return out
# Run in THE REAL TREE, not in a scratch copy. Snapshotting TREE while running
# the battery somewhere else made this check unable to observe a write at all:
# it reported "added none" while a run in TREE was leaving .pyc files behind.
# The claim is that a normal run writes nothing, so a normal run is what has to
# happen here. Nothing else in this file touches TREE.
before = snap(TREE)
rc13, _ = battery(TREE_ROOT)
after = snap(TREE)
changed = sorted(k for k in set(before) & set(after) if before[k] != after[k])
record(before == after and rc13 == 0, "filesystem unchanged after a normal run",
       f"{len(before)} paths before, {len(after)} after in {TREE.name}/, "
       f"__pycache__ included and every file content-hashed; added "
       f"{sorted(set(after) - set(before)) or 'none'}, rewritten in place "
       f"{changed or 'none'}. The absolute count is whatever tree this ran in, "
       f"which is smaller from a bundle extraction than author-side and larger "
       f"again if a stale __pycache__ is present; the before-and-after "
       f"comparison is the claim, not the number")
shutil.rmtree(base, ignore_errors=True)

# --------------------------------------------------------------------------
n_ok = sum(1 for ok, _, _ in RESULTS if ok)
print(f"\n  {n_ok}/{len(RESULTS)} self-checks pass")
lines = ["P4 FINAL SELF-CHECK, against the regenerated tree", ""]
for ok, title, detail in RESULTS:
    lines += [f"  {'ok  ' if ok else 'FAIL'} {title}", f"         {detail}"]
lines += ["", f"  {n_ok}/{len(RESULTS)} self-checks pass", ""]
(TREE / "P4_SELFCHECK.txt").write_text("\n".join(lines))
sys.exit(0 if n_ok == len(RESULTS) else 1)
