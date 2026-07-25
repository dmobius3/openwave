"""Lint MODELS.md: the 55-word cell budget + status/score-board consistency.

Three checks. The first two are defined in MODELS.md "Cell format (the 55-word
rule)"; the third keeps the derived score-board honest:
1. Budget: every summary cell in the per-model RESULTS BY MODEL tables carries
   at most LIMIT words of prose. Prose = cell text minus the status tag
   ("<icon> [status]"), minus every markdown link, minus "<br>-> ..." pointer
   tail segments (mid-prose arrows still count).
2. Sync: the summary-status table carries icons only, and each icon equals the
   status tag icon of the same criterion in the model's own table. Missing rows
   on either side are reported too.
3. Score-board: every count in the SCORE-BOARD table equals the tally of that
   icon over that model's rows, each total equals the criteria count, and every
   icon actually used has a score-board row (so adding a status back, e.g. the
   retired 🔶, cannot go uncounted).

Tables are found by SHAPE, not by heading text, so renaming a section does not
silently disable a check:
  - score-board table    = header "| **SCORE-BOARD** | [.. (M5)](..) | ... |"
  - summary-status table = header "| Criteria | [M5](..) | [M7](..) | ... |"
  - per-model table      = header "| Criteria | Status + result summary |"
    under the nearest "### <name> (M<n>)" heading.

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
MODEL_HEAD = re.compile(r"^#{2,4} .*\((M\d)\)\s*$")
MODEL_COL = re.compile(r"^\[(M\d)\]\(#")
BOARD_COL = re.compile(r"\((M\d)\)\]\(")  # "[Liquid Crystal<br>(M5)](path)"
NUM = re.compile(r"\d+")


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


def is_group_header(cells):
    """A domain divider row: "| **PARTICLES** | | | ... |"."""
    return cells[0].startswith("**") and not "".join(cells[1:]).strip()


def is_rule(cells):
    return set("".join(cells).strip()) <= set("- :")


status = {}  # (criterion, model) -> icon, from the summary-status table
summary = {}  # (criterion, model) -> (icon, wordcount, lineno), per-model tables
board = {}  # (icon or "total", model) -> (count, lineno), from the score-board
errors = []

mode = None  # "status" | "model" | "board" | None
model = None  # current model id when mode == "model"
model_cols = []  # column -> model id when mode in ("status", "board")
seen_model = None  # nearest "### ... (Mn)" heading above

for i, line in enumerate(PATH.read_text().splitlines(), 1):
    head = MODEL_HEAD.match(line)
    if head:
        seen_model = head.group(1)
    if not line.startswith("|"):
        if line.startswith("#"):
            mode = None  # a heading always closes the previous table
        continue

    cells = split_row(line)
    if is_rule(cells):
        continue

    if cells[0] == "Criteria":  # a header row: decide what table this is
        cols = [MODEL_COL.match(c) for c in cells[1:]]
        if len(cols) > 1 and all(cols):
            mode, model_cols = "status", [c.group(1) for c in cols]
        elif len(cells) == 2:
            mode, model = "model", seen_model
            if model is None:
                errors.append(f"L{i} per-model table with no '### ... (Mn)' heading above")
        else:
            mode = None
        continue

    cols = [BOARD_COL.search(c) for c in cells[1:]]
    if "SCORE-BOARD" in cells[0] and len(cols) > 1 and all(cols):
        mode, model_cols = "board", [c.group(1) for c in cols]
        continue

    if mode is None or is_group_header(cells):
        continue

    if mode == "board":
        icon = icon_of(cells[0]) or ("total" if "total" in cells[0].lower() else None)
        if icon is None:
            errors.append(f"L{i} score-board row '{cells[0]}' is neither an icon nor a total")
            continue
        for mod, cell in zip(model_cols, cells[1:]):
            n = NUM.search(cell)
            if not n:
                errors.append(f"L{i} score-board cell ({cells[0]}, {mod}) is not a number: {cell!r}")
            else:
                board[(icon, mod)] = (int(n.group()), i)
    elif mode == "status":
        for mod, cell in zip(model_cols, cells[1:]):
            status[(cells[0], mod)] = icon_of(cell)
            extra = prose_words(cell)
            if extra:
                errors.append(f"L{i} status cell ({cells[0]}, {mod}) carries prose: {extra}")
    elif mode == "model" and model:
        summary[(cells[0], model)] = (icon_of(cells[1]), len(prose_words(cells[1])), i)

if not status:
    errors.append("no summary-status table found (header '| Criteria | [M5](#..) | ... |')")
if not summary:
    errors.append("no per-model tables found (header '| Criteria | Status + result summary |')")

for key, (icon, n, ln) in sorted(summary.items(), key=lambda kv: kv[1][2]):
    crit, mod = key
    if n > LIMIT:
        errors.append(f"L{ln} ({crit}, {mod}) over budget: {n} > {LIMIT} words")
    if icon is None:
        errors.append(f"L{ln} ({crit}, {mod}) has no status icon")
    if key not in status:
        errors.append(f"L{ln} ({crit}, {mod}) missing from the summary-status table")
    elif status[key] != icon:
        errors.append(f"L{ln} ({crit}, {mod}) icon mismatch: status {status[key]} vs cell {icon}")

for key in sorted(status):
    if key not in summary:
        errors.append(f"status entry ({key[0]}, {key[1]}) has no row in that model's table")

if board:
    for (icon, mod), (claimed, ln) in sorted(board.items(), key=lambda kv: kv[1][1]):
        rows = [i for (c, m), i in status.items() if m == mod]
        actual = len(rows) if icon == "total" else rows.count(icon)
        if actual != claimed:
            label = "total criteria" if icon == "total" else icon
            errors.append(f"L{ln} score-board ({label}, {mod}): claims {claimed}, rows give {actual}")
    for mod in {m for _, m in status}:
        if not any(m == mod for _, m in board):
            errors.append(f"score-board has no column for {mod}")
        for icon in {i for (c, m), i in status.items() if m == mod}:
            if (icon, mod) not in board:
                errors.append(f"{mod} uses {icon} in its rows but the score-board has no {icon} row")
else:
    errors.append("no SCORE-BOARD table found (header '| **SCORE-BOARD** | [.. (M5)](..) | ... |')")

counts = sorted(n for _, n, _ in summary.values())
models = sorted({m for _, m in summary})
print(f"summary cells: {len(summary)} | status icons: {len(status)} | models: {','.join(models)}")
if counts:
    print(f"limit {LIMIT} | max {counts[-1]}  median {counts[len(counts) // 2]}")
if errors:
    print(f"\n{len(errors)} violation(s):")
    for e in errors:
        print(f"  {e}")
    sys.exit(1)
print("clean")
