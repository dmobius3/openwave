"""Lint MODELS.md: the 55-word cell budget + status/score-board consistency.

Six checks. The first two are defined in MODELS.md "Cell format (the 55-word
rule)"; the rest keep the derived tables honest:
1. Budget: every summary cell in the per-model RESULTS BY MODEL tables carries
   at most LIMIT words of prose. Prose = cell text minus the status tag
   ("<icon> [status]"), minus every markdown link, minus "<br>-> ..." pointer
   tail segments (mid-prose arrows still count). Those exclusions are
   deliberate: the cell exists to point at the record, so taxing its pointers
   would discourage linking. They also mean LIMIT is NOT comparable to the
   65-word Description cap in ROADMAP_STANDARDS.md, whose rule counts link
   labels; 55 here and 65 there already render the same size of cell. Measure
   before changing either (dev_docs/utils/models_cell_stats.py, T3).
2. Sync: the summary-status table carries icons only, and each icon equals the
   status tag icon of the same criterion in the model's own table. Missing rows
   on either side are reported too.
3. Score-board: every count in the SCORE-BOARD table equals the tally of that
   icon over that model's rows, each total equals the criteria count, and every
   icon actually used has a score-board row (so adding a status back, e.g. the
   retired 🔶, cannot go uncounted).
4. Regime: if the summary-status table carries a `regime` column, every
   criterion row declares one of REGIMES. That column is a property of the
   criterion rather than of any model (MODELS.md "Summary Status"), so it is
   matched by NAME and never counted as a model column.
5. Simplest test: the simplest-test companion table (header "| Criteria |
   simplest test |") must fill a non-empty test for every criterion, and its
   criteria set must match the summary-status table's, both directions: a row
   with no named test (or a criterion missing from either side) violates the
   "every row names its simplest passing test" rule that table exists to carry.
   A `simplest test` column carried by the summary-status table itself (the
   pre-companion layout) is checked the same way.
6. Shape: every data row carries exactly as many cells as its table's header.
   This is what an unescaped "|" inside a cell looks like from here, and it is
   worth its own message because such a row otherwise parses by position and
   silently reports the wrong column's contents.

Tables are found by SHAPE, not by heading text, so renaming a section does not
silently disable a check:
  - score-board table    = header "| **SCORE-BOARD** | [.. (M5)](..) | ... |"
  - summary-status table = header "| Criteria | [M5](..) | ... |", two or more
    model-link columns; any further column is matched by name (see check 4)
  - simplest-test table  = header "| Criteria | simplest test |" (two columns,
    the second named exactly that)
  - per-model table      = any other two-column header starting "| Criteria |"
    under the nearest "### <name> (M<n>)" heading.

Model columns are located by INDEX in the header, not by assuming they run to
the end of the row, so adding a criterion-level column cannot disable checks
2 through 4 (it did, when `regime` was introduced).

Usage: python3 dev_docs/utils/check_models_md.py [limit]   (default 55)
Exit 0 = clean, 1 = violations (listed on stdout).
"""

import re
import sys
from pathlib import Path

PATH = Path(__file__).resolve().parents[2] / "MODELS.md"
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 55

LINK = re.compile(r"\[[^\]]*\]\([^)]*\)")
TAG = re.compile(r"^\s*[^\[]*\[[^\]]*\]")  # leading icon + [status]
ICON = re.compile(r"(✅|❌|⚠️|🔶|🚧)")
MODEL_HEAD = re.compile(r"^#{2,4} .*\((M\d)\)\s*$")
MODEL_COL = re.compile(r"^\[(M\d)\]\(#")
BOARD_COL = re.compile(r"\((M\d)\)\]\(")  # "[Liquid Crystal<br>(M5)](path)"
NUM = re.compile(r"\d+")
REGIMES = {"static", "dynamic", "both"}  # MODELS.md "Summary Status" legend


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
tests = {}  # criterion -> lineno, from the simplest-test table (non-empty rows)
errors = []

mode = None  # "status" | "model" | "board" | "tests" | None
model = None  # current model id when mode == "model"
model_cols = []  # (column index, model id) pairs when mode in ("status", "board")
regime_col = None  # column index of the status table's `regime`, if it has one
test_col = None  # column index of the status table's `simplest test`, if it has one
width = None  # cell count of the current table's header row
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
        cols = [(j, m.group(1)) for j, c in enumerate(cells) if j and (m := MODEL_COL.match(c))]
        if len(cols) > 1:
            mode, model_cols, width = "status", cols, len(cells)
            named = {c.strip("* `").lower(): j for j, c in enumerate(cells) if j}
            regime_col = named.get("regime")
            test_col = named.get("simplest test")
            spare = [j for j in range(1, len(cells)) if j not in {k for k, _ in cols}]
            unknown = [cells[j] for j in spare if j not in (regime_col, test_col)]
            if unknown:
                errors.append(
                    f"L{i} summary-status header has unrecognized column(s): {unknown}."
                    " A criterion-level column must be registered here before use"
                    " (see `regime`), so that adding one cannot silently skip a check"
                )
        elif len(cells) == 2 and cells[1].strip("* `").lower() == "simplest test":
            mode, width = "tests", len(cells)
        elif len(cells) == 2:
            mode, model, width = "model", seen_model, len(cells)
            if model is None:
                errors.append(f"L{i} per-model table with no '### ... (Mn)' heading above")
        else:
            mode = None
        continue

    cols = [(j, m.group(1)) for j, c in enumerate(cells) if j and (m := BOARD_COL.search(c))]
    if "SCORE-BOARD" in cells[0] and len(cols) > 1:
        mode, model_cols, width = "board", cols, len(cells)
        continue

    if mode is None:
        continue

    if width is not None and len(cells) != width:
        errors.append(
            f"L{i} row '{cells[0][:40]}' has {len(cells)} cells, header has {width}"
            r" (an unescaped '|' inside a cell does this; write it as '\|')"
        )
        continue  # parsing by column index would report the wrong cells

    if is_group_header(cells):
        continue

    if mode == "board":
        icon = icon_of(cells[0]) or ("total" if "total" in cells[0].lower() else None)
        if icon is None:
            errors.append(f"L{i} score-board row '{cells[0]}' is neither an icon nor a total")
            continue
        for ci, mod in model_cols:
            cell = cells[ci]
            n = NUM.search(cell)
            if not n:
                errors.append(
                    f"L{i} score-board cell ({cells[0]}, {mod}) is not a number: {cell!r}"
                )
            else:
                board[(icon, mod)] = (int(n.group()), i)
    elif mode == "status":
        for ci, mod in model_cols:
            cell = cells[ci]
            status[(cells[0], mod)] = icon_of(cell)
            extra = prose_words(cell)
            if extra:
                errors.append(f"L{i} status cell ({cells[0]}, {mod}) carries prose: {extra}")
        if regime_col is not None:
            reg = cells[regime_col].strip("` ").lower()
            if reg not in REGIMES:
                errors.append(
                    f"L{i} regime cell ({cells[0]}) is {cells[regime_col]!r},"
                    f" expected one of {sorted(REGIMES)}"
                )
        if test_col is not None and not cells[test_col].strip("* `"):
            errors.append(
                f"L{i} simplest-test cell ({cells[0]}) is empty:"
                " every criterion row must name its simplest passing test"
            )
    elif mode == "tests":
        if not cells[1].strip("* `"):
            errors.append(
                f"L{i} simplest-test cell ({cells[0]}) is empty:"
                " every criterion row must name its simplest passing test"
            )
        else:
            tests[cells[0]] = i
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
            errors.append(
                f"L{ln} score-board ({label}, {mod}): claims {claimed}, rows give {actual}"
            )
    for mod in {m for _, m in status}:
        if not any(m == mod for _, m in board):
            errors.append(f"score-board has no column for {mod}")
        for icon in {i for (c, m), i in status.items() if m == mod}:
            if (icon, mod) not in board:
                errors.append(
                    f"{mod} uses {icon} in its rows but the score-board has no {icon} row"
                )
else:
    errors.append(
        "no SCORE-BOARD table found (header '| **SCORE-BOARD** | [.. (M5)](..) | ... |')"
    )

crit_set = {c for c, _ in status}
if tests:
    for c in sorted(crit_set - set(tests)):
        errors.append(f"criterion ({c}) missing from the simplest-test table")
    for c in sorted(set(tests) - crit_set):
        errors.append(
            f"L{tests[c]} simplest-test row ({c}) has no criterion in the summary-status table"
        )
elif status and test_col is None:
    errors.append(
        "no simplest-test table found (header '| Criteria | simplest test |')"
        " and the summary-status table carries no `simplest test` column"
    )

counts = sorted(n for _, n, _ in summary.values())
models = sorted({m for _, m in summary})
print(
    f"summary cells: {len(summary)} | status icons: {len(status)}"
    f" | tests: {len(tests)} | models: {','.join(models)}"
)
if counts:
    print(f"limit {LIMIT} | max {counts[-1]}  median {counts[len(counts) // 2]}")
if errors:
    print(f"\n{len(errors)} violation(s):")
    for e in errors:
        print(f"  {e}")
    sys.exit(1)
print("clean")
