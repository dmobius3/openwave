"""P1 step 3: the nine torsion values, extracted from every source that states
them, converted to exact Q(phi) triples, and cross-checked.

Target-free BY CONSTRUCTION, and that is the point: no value is written here.
Every closed form is PARSED out of a source file and evaluated exactly. A value
typed into this file would be a transcription that the cross-check could not
catch, since it would be comparing a source against itself.

Sources, per the plan's data inventory:
  A   theory side, torsion-correction.test.py, the corrected record
  B1  platform side, m8_3_method_note.md § 1.4 (the plan says § 4; § 4 is
      "Not computed" and states no values, so the plan's citation is wrong)
  B2  platform side, m8_3_mass_reproducer.py, the Part E gate targets

A and B are NOT independent mathematical evidence and this file does not claim
they are: both recompute the same spectral-zeta construction and gate against
closed forms with shared provenance. What agreement establishes is that the
packet faithfully captures the settled M8.3 targets. M8.8 is the independent
test of those targets, which is why it exists.
"""

import ast
import pathlib
import re
import sys
from fractions import Fraction as F

from qphi_exact import Phi

LABELS = [f"R{i}" for i in range(9)]


# ---------------------------------------------------------------- evaluation

def _eval(node):
    """Evaluate a parsed closed form exactly in Q(phi)."""
    if isinstance(node, ast.Expression):
        return _eval(node.body)
    if isinstance(node, ast.Constant):
        if isinstance(node.value, int):
            return Phi(node.value)
        raise ValueError(f"non-integer literal in a closed form: {node.value!r}")
    if isinstance(node, ast.Name):
        if node.id == "phi":
            return Phi(0, 1)
        raise ValueError(f"unknown name in a closed form: {node.id}")
    if isinstance(node, ast.Call):
        # mp.mpf(x) and mpf(x) are exact-integer wrappers in these sources
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", None)
        if name in ("mpf", "mpmathify") and len(node.args) == 1:
            return _eval(node.args[0])
        raise ValueError(f"unsupported call in a closed form: {name}")
    if isinstance(node, ast.UnaryOp):
        v = _eval(node.operand)
        if isinstance(node.op, ast.USub):
            return -v
        if isinstance(node.op, ast.UAdd):
            return v
        raise ValueError("unsupported unary operator")
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Pow):
            e = node.right
            neg = isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub)
            base = e.operand if neg else e
            if not (isinstance(base, ast.Constant) and isinstance(base.value, int)):
                raise ValueError("non-integer exponent")
            k = base.value * (-1 if neg else 1)
            return _eval(node.left) ** k
        a, b = _eval(node.left), _eval(node.right)
        if isinstance(node.op, ast.Add):
            return a + b
        if isinstance(node.op, ast.Sub):
            return a - b
        if isinstance(node.op, ast.Mult):
            return a * b
        if isinstance(node.op, ast.Div):
            return a / b
        raise ValueError("unsupported binary operator")
    raise ValueError(f"unsupported syntax: {type(node).__name__}")


def closed_form(expr):
    """Evaluate a closed-form string exactly. Accepts '^' for '**' and the
    implicit multiplication used in the method note's display block."""
    s = expr.strip().rstrip(".,;")
    s = s.replace("^", "**")
    s = re.sub(r"(\d|\))\s*(?=phi\b)", r"\1*", s)   # (4/5)phi -> (4/5)*phi
    s = re.sub(r"(?<=phi)\s*(?=\()", "*", s)         # phi(...) -> phi*(...)
    return _eval(ast.parse(s, mode="eval"))


# ------------------------------------------------------------------ sources

# Extractors return a LIST of (label, expression) occurrences, never a dict.
# Folding to a dict inside an extractor destroys conflicts before adjudication:
# a dict literal silently keeps the last duplicate key and a setdefault silently
# keeps the first. Both were doing exactly that here. "Stop on any disagreement"
# has to include disagreement WITHIN a source, so occurrences are carried
# through and collapsed only after they have been compared.

def _dict_literal(txt, name):
    """Extract `name = { ... }` as a list of (key, expression) occurrences."""
    m = re.search(rf"\b{name}\s*=\s*\{{", txt)
    if not m:
        return []
    i = m.end() - 1
    depth, j = 0, i
    while j < len(txt):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    body = txt[i:j + 1]
    tree = ast.parse(body, mode="eval").body
    out = []
    for k, v in zip(tree.keys, tree.values):
        if isinstance(k, ast.Constant) and k.value in LABELS:
            out.append((k.value, ast.get_source_segment(body, v) or ast.unparse(v)))
    return out


def source_theory(path):
    txt = pathlib.Path(path).read_text()
    got = []
    for name in ("tgt_int", "tgt_half"):
        got += _dict_literal(txt, name)
    got += [("R7", m.group(1)) for m in
            re.finditer(r'T2\["R7"\]\s*/\s*\((.+?)\)\s*-\s*1', txt)]
    got += [("R0", m.group(1).strip()) for m in
            re.finditer(r'T2\["R0"\]\s*=\s*(.+)', txt)]
    return got


def source_reproducer(path):
    txt = pathlib.Path(path).read_text()
    got = _dict_literal(txt, "exact")
    got += [("R0", m.group(1).strip()) for m in
            re.finditer(r'T2\["R0"\]\s*=\s*(.+)', txt)]
    return got


NOTE_SECTION = r"###\s*1\.4\b"


def source_method_note(path):
    """The § 1.4 VALUE TABLE, and only that.

    Scoped twice over. Document-wide would take whichever occurrence came first
    and ignore a conflicting restatement elsewhere. Section-wide is still wrong:
    the same section states the identities, and a ratio of two rows matches a
    naive value pattern, yielding a spurious value for the row in the
    denominator. So the extraction is confined to the fenced block that carries
    the table, and the block must be unique and complete."""
    txt = pathlib.Path(path).read_text()
    m = re.search(NOTE_SECTION, txt)
    if not m:
        raise ValueError(f"method note § 1.4 not found in {path}")
    nxt = re.search(r"^###?\s", txt[m.end():], re.M)
    body = txt[m.end():m.end() + (nxt.start() if nxt else len(txt))]

    blocks = re.findall(r"```[a-z]*\n(.*?)```", body, re.S)
    tables = [b for b in blocks if re.search(r"T\^2\(R0\)", b)]
    if len(tables) != 1:
        raise ValueError(f"expected exactly one § 1.4 value table, found {len(tables)}")
    occ = [(g.group(1), g.group(2)) for g in
           re.finditer(r"T\^2\((R\d)\)\s*=\s*(\S+)", tables[0])]
    if sorted(l for l, _ in occ) != LABELS:
        raise ValueError(
            f"§ 1.4 value table does not state each label exactly once: "
            f"{sorted(l for l, _ in occ)}")
    return occ


# ----------------------------------------------------------------- checking

def identities_from_note(path):
    """Parse the § 1.4 identity statements into structured factor form.

    Returns [{"id", "kind", "factors": [(label, exponent)...], "stated": expr}].

    Memberships are PARSED, never written here. Hard-coding which rows appear in
    which identity puts quarantined STRUCTURE into a file classified target-free,
    and the numeric leak gate cannot see structure: it matches renderings of
    values, so a membership list passes it silently. The gate's own SCOPE notice
    says exactly this. Parsing removes the question instead of arguing it.

    The id is derived from the factors, so an identity cannot carry a name that
    disagrees with what it computes. That divergence is not guarded against
    elsewhere because it is not constructible here."""
    txt = pathlib.Path(path).read_text()
    m = re.search(NOTE_SECTION, txt)
    if not m:
        raise ValueError(f"method note § 1.4 not found in {path}")
    nxt = re.search(r"^###?\s", txt[m.end():], re.M)
    body = txt[m.end():m.end() + (nxt.start() if nxt else len(txt))]
    body = re.sub(r"```[a-z]*\n.*?```", "", body, flags=re.S)   # drop the value table

    out = []
    for span in re.findall(r"`([^`]+)`", body):
        if "T^2(" not in span or "=" not in span:
            continue
        # Only spans whose every T^2 argument is a row label. The section also
        # states the general propagation formula, whose arguments are symbols.
        args = re.findall(r"T\^2\(([^)]*)\)", span)
        if not args or any(not re.fullmatch(r"R\d", x.strip()) for x in args):
            continue
        lhs, rhs = span.split("=", 1)
        parts = lhs.split("/")
        factors = []
        for k, part in enumerate(parts):
            labs = re.findall(r"T\^2\((R\d)\)", part)
            if not labs:
                raise ValueError(f"identity side with no rows: {span!r}")
            factors += [(l, -1 if k else 1) for l in labs]
        if len(factors) < 2:
            continue
        kind = ("ratio" if len(parts) > 1 else "product")
        out.append({
            "id": f"{kind}:" + ",".join(f"{l}^{e}" for l, e in factors),
            "kind": kind,
            "factors": factors,
            "stated": rhs.strip().rstrip(".,;`"),
            "span": span,
        })
    if not out:
        raise ValueError(f"no identities parsed from § 1.4 of {path}")
    return out


# Source FAMILIES, which is what the plan's evidentiary claim is about. B1 and B2
# are two renderings of the same platform-side record, so counting renderings
# instead of families lets the entire theory side vanish while every row stays
# green. The count was doing exactly that.
SOURCES = {
    "A  theory/torsion-correction": ("theory", source_theory),
    "B1 platform/method-note-1.4": ("platform", source_method_note),
    "B2 platform/reproducer-partE": ("platform", source_reproducer),
}
REQUIRED = ("A  theory/torsion-correction", "B2 platform/reproducer-partE")


def collect(theory, note, reproducer):
    paths = {"A  theory/torsion-correction": theory,
             "B1 platform/method-note-1.4": note,
             "B2 platform/reproducer-partE": reproducer}
    parsed, errors, conflicts = {}, [], []
    for sname, (_family, fn) in SOURCES.items():
        occurrences = {}
        for lab, expr in fn(paths[sname]):
            # Evaluate BEFORE touching the accumulator. setdefault-then-append
            # creates the list first, so a raising evaluation left an empty list
            # behind and the later occ[0] raised IndexError instead of the run
            # rejecting cleanly on VAL_PARSE.
            try:
                triple = closed_form(expr).triple()
            except Exception as e:                       # noqa: BLE001
                errors.append((sname, lab, expr, str(e)))
                continue
            occurrences.setdefault(lab, []).append((expr, triple))
        parsed[sname] = {}
        for lab, occ in occurrences.items():
            distinct = {t for _, t in occ}
            if len(distinct) > 1:
                conflicts.append(
                    f"{sname} states {lab} more than once with different values: "
                    f"{sorted(distinct)}")
                continue
            if len(occ) > 1:
                conflicts.append(
                    f"{sname} states {lab} {len(occ)} times (agreeing, but this "
                    f"location is declared one row per label)")
            parsed[sname][lab] = occ[0]
    return parsed, errors, conflicts


def adjudicate(parsed):
    """Per label: the distinct triples across sources, coverage, and families."""
    rows = []
    for lab in LABELS:
        have = {s: v[lab] for s, v in parsed.items() if lab in v}
        triples = {v[1] for v in have.values()}
        rows.append({
            "label": lab,
            "sources": sorted(have),
            "n_sources": len(have),
            "families": sorted({SOURCES[s][0] for s in have}),
            "missing_required": [s for s in REQUIRED if s not in have],
            "agree": len(triples) == 1,
            "triple": next(iter(triples)) if len(triples) == 1 else None,
            "distinct": sorted(triples),
            "exprs": {s: have[s][0] for s in sorted(have)},
        })
    return rows


def verdict(rows, errors, conflicts=()):
    """STOP conditions: any disagreement between sources, any disagreement
    WITHIN a source, any missing required carrier, any unevaluable form.

    Coverage is by FAMILY and by required carrier, not by a count of renderings.
    Both the theory-side record and the independently executed platform
    reproducer must carry every row; a count of two would be satisfied by the
    two platform renderings alone."""
    problems = list(conflicts)
    if errors:
        problems.append(f"{len(errors)} closed form(s) failed to evaluate")
    for r in rows:
        if r["missing_required"]:
            problems.append(
                f"{r['label']} is not carried by {', '.join(r['missing_required'])}")
        elif {"theory", "platform"} - set(r["families"]):
            problems.append(f"{r['label']} lacks a source family: {r['families']}")
        elif not r["agree"]:
            problems.append(f"{r['label']} DISAGREES across sources: {r['distinct']}")
    return problems
