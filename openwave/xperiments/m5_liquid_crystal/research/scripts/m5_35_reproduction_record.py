"""M5.35: compile the reproduction record of the scorecard v1.2 delta (review_checks/, V40 to V51 from the
v1.2 copy, V52 to V57) from the pool logs and the comparator output into data/m5_35_reproduction.json.

Inputs (paths as arguments): the logs folder (one <lane>_<name>.log per job with START/END lines and
the script's own prints), the comparator text (one line per regenerated JSON: leaves, max relative
difference, its leaf). Output: per job exit code, wall seconds, the README's one-core time where
stated, the comparator's max relative difference and its leaf; plus the headline values read from
the regenerated JSONs (E15 of V44, V48, V52, V54; the DMC energies from the logs).
Usage: python3 m5_35_reproduction_record.py <logs_dir> <compare_txt> <run_dir> [--write]
"""
import json, os, re, sys

README_S = {  # the README's one-core times (s) for V40 to V57 and review_checks (2-core container, Python 3.13)
    "r1c2_V40": 39, "r1c2_V40b": 109, "r1c2_V41": 2, "r1c2_V42": 0.5, "r1c2_V43": 22, "r1c2_V44_0": 137,
    "r1c2_V44": 421, "r1c2_V44b": 28, "r1c2_V45": 3, "r1c2_V47": 1, "r1c2_V46": 4, "r1c2_V48_ell": 549,
    "r1c2_V48_fixed": 484, "r1c2_V48b": 49, "r1c2_V49": 1, "r1c2_V51": 42, "r1c2_V50": 11,
    "r1r_frozen_bag": 15, "r1r_svm": 540, "r1r_cib": 120, "r1r_v11_v5check": 44, "r1r_v11_v5free": 44,
    "r1r_v3_hess": 1, "r1r_v3_inertia": 1, "r1r_v3_stress3d": 2, "r1r_v11_nums": 0.1,
}


def parse_logs(d):
    out = {}
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".log") or fn.startswith("_"): continue
        txt = open(os.path.join(d, fn)).read()
        ends = re.findall(r"^END (\S+) exit=(\d+) wall_s=(\d+)", txt, re.M)
        starts = re.findall(r"^START (\S+) (.*)$", txt, re.M)
        name = fn[:-4]
        rec = {"command": starts[-1][1] if starts else None,
               "exit": int(ends[-1][1]) if ends else None, "wall_s": int(ends[-1][2]) if ends else None,
               "readme_s": README_S.get(name)}
        m = re.search(r"^(refit|mtv) N=(\d) .*E = (\S+) \+- (\S+)", txt, re.M)
        if m: rec["dmc"] = {"case": m.group(1), "N": int(m.group(2)), "E": float(m.group(3)), "err": float(m.group(4))}
        out[name] = rec
    return out


def parse_compare(p):
    out = {}
    for line in open(p):
        m = re.match(r"(\S+\.json): (\d+) leaves, max rel diff (\S+) at (.*?) \(([^()]*)\)$", line.strip())
        if m: out[m.group(1)] = {"leaves": int(m.group(2)), "max_rel_diff": float(m.group(3)), "leaf": m.group(4), "values": m.group(5)}
    return out


def headlines(run):
    def J(f): return json.load(open(os.path.join(run, f)))
    h = {}
    h["V44_E15_by_drag"] = {k: v["E15"] for k, v in J("V44_comoving_clock_refit.json").items() if isinstance(v, dict) and "E15" in v}
    for tag in ("ell", "fixed"):
        h[f"V48_{tag}_E15_by_W"] = {k: v["E15"] for k, v in J(f"V48_covariant_omega_refit_{tag}.json").items()}
    rb = J("V48b_readout.json")
    h["V48b"] = {k: {kk: v[kk] for kk in ("E15", "M", "Delta", "gA", "rMV2") if kk in v} for k, v in rb.items() if isinstance(v, dict)}
    for t in ("n4", "lock"): h[f"V52_{t}_E15"] = J(f"V52_refit_{t}.json")["5.8"]["E15"]
    for t in ("nopipi", "consistent", "pole", "poleguard", "vbp"): h[f"V54_{t}_E15"] = J(f"V54_refit_{t}.json")["5.8"]["E15"]
    h["V54b_E15"] = {k: {kk: v[kk] for kk in ("E15", "E_own_scoring") if kk in v} for k, v in J("V54b_replay_fine.json").items() if isinstance(v, dict)}
    v57 = J("V57_delta_level_vs_pole.json")
    h["V57"] = {k: v for k, v in v57.items() if not isinstance(v, (dict, list))}
    return h


def main():
    logs, cmp, run = sys.argv[1:4]
    rec = {"task": "M5.35", "versions": {"python": "3.12.14", "numpy": "2.5.3", "scipy": "1.18.1", "sympy": "1.14.0", "mpmath": "1.3.0"},
           "jobs": parse_logs(logs), "comparator": parse_compare(cmp), "headlines": headlines(run)}
    rec["summary"] = {"jobs": len(rec["jobs"]), "exit_nonzero": [k for k, v in rec["jobs"].items() if v["exit"] not in (0, None)],
                      "wall_total_s": sum(v["wall_s"] or 0 for v in rec["jobs"].values()),
                      "regenerated_json": len(rec["comparator"]),
                      "max_rel_diff_over_1e-3": {k: v["max_rel_diff"] for k, v in rec["comparator"].items() if v["max_rel_diff"] > 1e-3}}
    if "--write" in sys.argv:
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "m5_35_reproduction.json")
        json.dump(rec, open(out, "w"), indent=1, default=float); print("wrote", os.path.normpath(out))
    print(json.dumps(rec["summary"], indent=1, default=float))


if __name__ == "__main__":
    main()
