#!/usr/bin/env python3
"""Regenerate the _DATASETS.md manifest for a research data/ folder.

Policy (2026-07-20, supersedes the "delete raw data > 1 MB" rule): heavy binary
arrays under any data/ folder are LOCAL-ONLY (gitignored), never deleted. The
tracked evidence stays in git: summary .json / .csv / .txt, plots, scripts. This
manifest is the index of what lives locally and how to rebuild it.

Grouping: files are bucketed by task prefix (m5_21_6, t12, ...). For each group
the producing script and the task record are DERIVED by convention, not guessed:
a script matching the prefix under ../scripts/, a task doc matching the prefix
under ../tasks/. Anything not resolvable is printed as "not resolved" so the gap
is visible instead of invented.

Usage:
    python3 dev_docs/utils/gen_datasets_manifest.py <path/to/data> [--write]

Without --write it prints the manifest to stdout (diff-friendly dry run).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Binary array formats that are local-only under the policy.
BLOB_SUFFIXES = {".npz", ".npy", ".h5", ".hdf5", ".pkl", ".pt", ".mat"}

# Prefix = the task id at the head of the filename: m5_21_6, m5_8_2cb, t12.
# Leading underscores are scratch markers and are stripped first.
PREFIX_RE = re.compile(r"^(m\d+(?:_\d+)*[a-z]*|t\d+)")

# Invocation printed into the generated manifest header. Repo-relative on
# purpose: a clone can run it, and no absolute local path reaches the repo.
SELF = "dev_docs/utils/gen_datasets_manifest.py"


def prefix_of(stem: str) -> str:
    m = PREFIX_RE.match(stem.lstrip("_"))
    return m.group(1) if m else "(unprefixed)"


def human_mb(n: int) -> str:
    return f"{n / 1048576:.2f} MB"


def resolve_producers(prefix: str, scripts_dir: Path) -> list[str]:
    if not scripts_dir.is_dir():
        return []
    hits = sorted(p.name for p in scripts_dir.glob(f"{prefix}*.py"))
    if hits:
        return hits
    # Fall back to progressively shorter prefixes (m5_21_6 -> m5_21 -> m5).
    parts = prefix.split("_")
    while len(parts) > 1:
        parts.pop()
        stub = "_".join(parts)
        hits = sorted(p.name for p in scripts_dir.glob(f"{stub}*.py"))
        if hits:
            return hits
    return []


def resolve_record(prefix: str, data_dir: Path) -> tuple[str, str] | None:
    """Return (relative_path_to_dir, filename) for the doc recording this group's runs.

    Two layouts are probed: tasks/ and findings/ beside data/ (the
    research/{data,tasks,findings} layout), and the same one level up.
    """
    candidates_dirs = [
        ("..", data_dir.parent / "tasks"),
        ("..", data_dir.parent / "findings"),
        ("../..", data_dir.parent.parent / "tasks"),
        ("../..", data_dir.parent.parent / "findings"),
    ]
    stubs = [prefix]
    parts = prefix.split("_")
    while len(parts) > 1:
        parts.pop()
        stubs.append("_".join(parts))
    for stub in stubs:
        for rel, d in candidates_dirs:
            if not d.is_dir():
                continue
            for cand in (
                f"{stub}_task_details.md",
                f"{stub}_task.md",
                f"{stub}_findings.md",
                f"{stub}_note.md",
            ):
                if (d / cand).exists():
                    return f"{rel}/{d.name}", cand
    return None


def build(data_dir: Path, shown_path: str) -> str:
    scripts_dir = data_dir.parent / "scripts"

    groups: dict[str, list[Path]] = {}
    for f in sorted(data_dir.iterdir()):
        if f.is_file() and f.suffix.lower() in BLOB_SUFFIXES:
            groups.setdefault(prefix_of(f.stem), []).append(f)

    total_files = sum(len(v) for v in groups.values())
    total_bytes = sum(f.stat().st_size for v in groups.values() for f in v)

    lines: list[str] = []
    lines.append("# Local-only datasets manifest")
    lines.append("")
    lines.append(
        f"> AUTO-GENERATED, do not hand-edit the table: "
        f"`python3 {SELF} {shown_path} --write`"
    )
    lines.append("")
    lines.append(
        "Heavy binary arrays in this folder are **local-only**: gitignored, never deleted "
        "(policy 2026-07-20, which supersedes the earlier \"delete raw data > 1 MB\" rule). "
        "They stay on the working machine so later tasks can consume them directly, and they "
        "stay OUT of the repo so clones stay light. What IS tracked in git and readable on "
        "GitHub: the summary `.json` / `.csv` / `.txt` in this same folder, the plots, and the "
        "scripts that rebuild everything here."
    )
    lines.append("")
    lines.append(
        f"**Inventory**: {total_files} local-only files, {human_mb(total_bytes)}, "
        f"in {len(groups)} task groups."
    )
    lines.append("")
    lines.append("| Task group | Files | Size | Producing script(s) | Record (regen commands + context) |")
    lines.append("| --- | --- | --- | --- | --- |")

    for prefix in sorted(groups):
        files = groups[prefix]
        size = human_mb(sum(f.stat().st_size for f in files))
        producers = resolve_producers(prefix, scripts_dir)
        if producers:
            prod = " · ".join(f"[`{p}`](../scripts/{p})" for p in producers[:3])
            if len(producers) > 3:
                prod += f" (+{len(producers) - 3} more)"
        else:
            prod = "⚠️ not resolved"
        found = resolve_record(prefix, data_dir)
        rec = f"[`{found[1]}`]({found[0]}/{found[1]})" if found else "⚠️ not resolved"
        lines.append(f"| `{prefix}` | {len(files)} | {size} | {prod} | {rec} |")

    lines.append("")
    lines.append(
        "**Regeneration**: the exact command + runtime per dataset lives in the task record "
        "linked on its row (the task_details / findings doc), which is where the run "
        "configuration is already written down. Runs are deterministic from their fixed "
        "seeds and configs, so a regenerated array reproduces the original bit-for-bit at "
        "the stored precision."
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    if not data_dir.is_dir():
        print(f"not a directory: {data_dir}", file=sys.stderr)
        return 1

    text = build(data_dir, str(args.data_dir))
    if args.write:
        out = data_dir / "_DATASETS.md"
        out.write_text(text)
        print(f"wrote {out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
