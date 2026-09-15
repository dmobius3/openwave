"""
One-shot encoding fixer. Run once from inside pipeline_engine/.

Recurses into subdirectories. Reads each .py file (UTF-8 strict, fall back
to cp1252), replaces non-ASCII characters with ASCII equivalents, writes
back as UTF-8 without BOM.

Safe to run multiple times.
"""

from pathlib import Path

REPLACEMENTS = [
    ("\u2014", "--"),  # em-dash
    ("\u2013", "-"),  # en-dash
    ("\u2192", "->"),  # right arrow
    ("\u2190", "<-"),  # left arrow
    ("\u00b5", "u"),  # micro sign
    ("\u00d7", "x"),  # multiplication sign
    ("\u00b1", "+/-"),  # plus-minus
    ("\u2026", "..."),  # ellipsis
    ("\u201c", '"'),
    ("\u201d", '"'),
    ("\u2018", "'"),
    ("\u2019", "'"),
]


def is_ascii(s):
    try:
        s.encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


def fix(path):
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
        source = "utf-8"
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
        source = "cp1252"

    original = text
    for old, new in REPLACEMENTS:
        text = text.replace(old, new)

    remaining = [c for c in text if not is_ascii(c)]
    if remaining:
        chars = "".join(sorted(set(remaining)))
        print(f"  ! {path}: still non-ASCII: {chars!r}")

    path.write_text(text, encoding="utf-8", newline="\n")
    changed = "yes" if text != original else "no"
    rel = path.relative_to(path.parents[0])
    print(f"  {rel}: read={source}, replaced={changed}, wrote=utf-8")


def main():
    here = Path(__file__).resolve().parent
    print(f"Fixing files under: {here}")
    files = [p for p in here.rglob("*.py") if "__pycache__" not in p.parts]
    for p in sorted(files):
        fix(p)
    print(f"Done. {len(files)} files processed.")


if __name__ == "__main__":
    main()
