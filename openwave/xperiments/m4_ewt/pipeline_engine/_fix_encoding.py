"""
One-shot encoding fixer. Run once from inside pipeline_engine/.

Reads each .py file (UTF-8 strict, fall back to cp1252), replaces non-ASCII
characters with ASCII equivalents, writes back as UTF-8 without BOM.

Safe to run multiple times.
"""
from pathlib import Path

FILES = [
    "__init__.py",
    "context.py",
    "pipeline.py",
    "sinks.py",
    "loggers.py",
    "runner.py",
    "_smoke_test.py",
]

# Order matters: longer patterns first.
REPLACEMENTS = [
    ("—", "--"),    # em-dash
    ("–", "-"),     # en-dash
    ("→", "->"),    # arrow
    ("←", "<-"),
    ("µ", "u"),     # micro sign (us)
    ("×", "x"),
    ("±", "+/-"),
    ("…", "..."),
    ("“", '"'),
    ("”", '"'),
    ("‘", "'"),
    ("’", "'"),
]


def is_ascii(s: str) -> bool:
    try:
        s.encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


def fix(path: Path) -> None:
    raw = path.read_bytes()

    # Try UTF-8 strict first.
    try:
        text = raw.decode("utf-8")
        source = "utf-8"
    except UnicodeDecodeError:
        # Fall back to cp1252 (Windows default).
        text = raw.decode("cp1252")
        source = "cp1252"

    original = text
    for old, new in REPLACEMENTS:
        text = text.replace(old, new)

    # Report any remaining non-ASCII so nothing slips through.
    remaining = [c for c in text if not is_ascii(c)]
    if remaining:
        chars = "".join(sorted(set(remaining)))
        print(f"  ! {path.name}: still has non-ASCII after replacement: {chars!r}")

    path.write_text(text, encoding="utf-8", newline="\n")

    changed = "yes" if text != original else "no"
    print(f"  {path.name}: read={source}, replaced={changed}, wrote=utf-8")


def main() -> None:
    here = Path(__file__).resolve().parent
    print(f"Fixing files in: {here}")
    for name in FILES:
        p = here / name
        if not p.exists():
            print(f"  - {name}: not found, skipped")
            continue
        fix(p)
    print("Done.")


if __name__ == "__main__":
    main()