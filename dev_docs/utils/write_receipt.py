"""Write a validation receipt for a reproduction or verification run (platform task T10).

A receipt records WHAT RAN, in a fixed shape two receipts can be diffed on; it is never evidence of
correctness (the verdict words stay in the task record with their gates). Shape: METHOD_NOTE.md
§ The validation receipt.

Usage:
  write_receipt.py --task M5.35 --base <commit> --package <commit|doi|url> --logs <dir>
                   [--compare <comparator.txt>] [--machine workstation|laptop|container]
                   [--verdict "<words>"] --out research/data/m5_35_receipt.json
Logs: one file per job with a `START <hh:mm:ss> <command>` line and an `END <hh:mm:ss> exit=<n>
wall_s=<n>` line (the pool launchers of M5.33 to M5.35 write this shape). Comparator lines:
`<file>.json: <n> leaves, max rel diff <x> at <leaf> (<a> vs <b>)`.
"""

import argparse, datetime, json, os, re, sys


def jobs_from_logs(d):
    out = {}
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".log") or fn.startswith("_"):
            continue
        txt = open(os.path.join(d, fn), errors="replace").read()
        st = re.findall(r"^START \S+ (.*)$", txt, re.M)
        # the interpreter is recorded by name, never by its path (a receipt carries no machine paths)
        st = [re.sub(r"^\S*/", "", c, count=1) for c in st]
        en = re.findall(r"^END \S+ exit=(\d+) wall_s=(\d+)", txt, re.M)
        out[fn[:-4]] = {
            "command": st[-1] if st else None,
            "exit": int(en[-1][0]) if en else None,
            "wall_s": int(en[-1][1]) if en else None,
            "log": fn,
        }
    return out


def comparator_from(p):
    out = {}
    for line in open(p, errors="replace"):
        m = re.match(
            r"(?:\S+ )?(\S+\.json): (\d+) leaves, max rel diff (\S+) at (.*?) \(([^()]*)\)$",
            line.strip(),
        )
        if m:
            out[m.group(1)] = {
                "leaves": int(m.group(2)),
                "max_rel_diff": float(m.group(3)),
                "leaf": m.group(4),
            }
    return out


def versions():
    v = {"python": sys.version.split()[0]}
    for name in ("numpy", "scipy", "sympy", "mpmath"):
        try:
            v[name] = __import__(name).__version__
        except Exception:
            v[name] = None
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--package", required=True)
    ap.add_argument("--logs", required=True)
    ap.add_argument("--compare")
    ap.add_argument(
        "--machine", default="workstation", choices=["workstation", "laptop", "container"]
    )
    ap.add_argument("--verdict", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    jobs = jobs_from_logs(a.logs)
    rec = {
        "task": a.task,
        "written": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "base": a.base,
        "package": a.package,
        "machine": a.machine,
        "versions": versions(),
        "jobs": jobs,
        "comparator": comparator_from(a.compare) if a.compare else {},
        "verdict": a.verdict,
    }
    bad = [k for k, j in jobs.items() if j["exit"] not in (0, None)]
    missing = [k for k, j in jobs.items() if j["exit"] is None]
    rec["summary"] = {
        "jobs": len(jobs),
        "exit_nonzero": bad,
        "no_end_line": missing,
        "wall_total_s": sum(j["wall_s"] or 0 for j in jobs.values()),
        "files_compared": len(rec["comparator"]),
    }
    json.dump(rec, open(a.out, "w"), indent=1)
    print(
        f"wrote {a.out}: {len(jobs)} jobs, {len(bad)} non-zero, {len(missing)} without END, "
        f"{len(rec['comparator'])} files compared"
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
