"""Lint MODELS.md: the 55-word cell budget + matrix-icon sync.

Two checks, both defined in MODELS.md "Cell format (the 55-word rule)":
1. Budget: every summary cell in the per-model RESULTS BY MODEL tables carries
   at most LIMIT words of prose. Prose = cell text minus the status tag
   ("<icon> [status]"), minus every markdown link, minus "<br>-> ..." pointer
   tail segments (mid-prose arrows still count).
2. Sync: the at-a-glance matrix carries icons only, and each icon equals the
   status tag icon of the same criterion in the model's own table.

Usage: python3 dev_docs/check_models_md.py [limit]   (default 55)
Exit 0 = clean, 1 = violations (listed on stdout).
"""

import re
import sys
from pathlib import Path

PATH = Path(__file__).resolve().parent.parent / "MODELS.md"
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 55

LINK = re.compile(r"\[[^\]]*\]\([^)]*\)")
TAG = re.compile(r"^\s*[^\[]*\[[^\]]*\]")  # leading icon + [status]
ICON = re.compile(r"(✅|❌|⚠️|🔶|🚧)")
MODEL_HEAD = re.compile(r"^### .*\((M\d)\)")
MATRIX_COL = re.compile(r"\[(M\d)\]")


def split_row(line):
    parts = re.split(r"(?<!\\)\|", line)
    return [p.strip() for p in parts][1:-1]


def prose_words(cell):
    segs = [s for s in cell.split("<br>") if not s.strip().startswith("→")]
    body = TAG.sub("", " ".join(segs), count=1)
    body = LINK.sub(" ", body)
    body = re.sub(r"[`*]", "", body)
    return [w for w in body.split() if re.search(r"[0-9A-Za-zα-ωΑ-Ω]", w)]


def icon_of(cell):
    m = ICON.search(cell)
    return m.group(1) if m else None


def is_divider(cells):
    return set("".join(cells).strip()) <= set("- ") or cells[0].startswith("**")


matrix = {}  # (criterion, model) -> icon
sections = {}  # (criterion, model) -> (icon, wordcount, lineno)
errors = []

lines = PATH.read_text().splitlines()
in_matrix = False
model = None
header_models = []

for i, line in enumerate(lines, 1):
    if line.startswith("### Status at a glance"):
        in_matrix, model = True, None
        continue
    m = MODEL_HEAD.match(line)
    if m:
        in_matrix, model = False, m.group(1)
        continue
    if line.startswith("## ") and not line.startswith("## COVERAGE"):
        in_matrix, model = False, None
        continue
    if not line.startswith("| "):
        continue
    cells = split_row(line)
    if in_matrix:
        if cells[0] == "Criteria":
            header_models = [MATRIX_COL.search(c).group(1) for c in cells[1:]]
            continue
        if is_divider(cells):
            continue
        for mod, cell in zip(header_models, cells[1:]):
            matrix[(cells[0], mod)] = icon_of(cell)
            extra = prose_words(cell)
            if extra:
                errors.append(f"L{i} matrix cell ({cells[0]}, {mod}) carries prose: {extra}")
    elif model:
        if cells[0] == "Criteria" or is_divider(cells):
            continue
        sections[(cells[0], model)] = (icon_of(cells[1]), len(prose_words(cells[1])), i)

if not matrix:
    errors.append("no at-a-glance matrix found")
if not sections:
    errors.append("no RESULTS BY MODEL tables found")

for key, (icon, n, ln) in sorted(sections.items(), key=lambda kv: kv[1][2]):
    crit, mod = key
    if n > LIMIT:
        errors.append(f"L{ln} ({crit}, {mod}) over budget: {n} > {LIMIT} words")
    if key not in matrix:
        errors.append(f"L{ln} ({crit}, {mod}) missing from the at-a-glance matrix")
    elif matrix[key] != icon:
        errors.append(f"L{ln} ({crit}, {mod}) icon mismatch: matrix {matrix[key]} vs section {icon}")

for key in sorted(matrix):
    if key not in sections:
        errors.append(f"matrix entry ({key[0]}, {key[1]}) has no row in the model's table")

counts = sorted(n for _, n, _ in sections.values())
print(f"summary cells: {len(sections)} | matrix icons: {len(matrix)} | limit {LIMIT}")
if counts:
    print(f"max {counts[-1]}  median {counts[len(counts) // 2]}")
if errors:
    print(f"\n{len(errors)} violation(s):")
    for e in errors:
        print(f"  {e}")
    sys.exit(1)
print("clean")
