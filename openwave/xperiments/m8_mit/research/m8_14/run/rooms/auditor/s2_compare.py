"""Stage 2: script / JSON / report consistency for both rooms.

Run after s2_rerun.py, as `./py s2_compare.py`.  Writes s2_compare.json.

  C1  every numeric leaf of the rerun's results.json equals the room's results.json:
      |x - y| / (|x| + 1e-3) < 1e-9, i.e. relative 1e-9 with an absolute floor of 1e-12
      for derived small differences (ARPACK may differ in the last digits between runs)
  C2  every T4 number printed in RETURN.md (levels, p, extrapolation, error) equals the room's
      results.json at the printed precision
  C3  every T3 lowest-six value and T5 lowest-Jacobi value printed in RETURN.md equals the JSON
      at the printed precision
Planted twins: one JSON leaf altered by 1e-6 relative (C1); one printed eigenvalue altered in
the last printed digit (C2); one printed T3 value altered (C3).  They must print FAIL.
"""
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
lines, out = [], {}
LABS = ["1/4", "1/2", "1", "7/5", "3/2", "pi/2", "17/10"]


def report(name, ok, planted=False):
    s = f"{'[planted] ' if planted else ''}{'PASS' if ok else 'FAIL'}  {name}"
    print(s, flush=True)
    lines.append(s)
    return ok


def leaves(o, path=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, f"{path}/{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, f"{path}[{i}]")
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        yield path, float(o)


def c1(room, planted=False):
    a = dict(leaves(json.load(open(os.path.join(HERE, f"room_{room}", "results.json")))))
    b = dict(leaves(json.load(open(os.path.join(HERE, f"rerun_{room}", "results.json")))))
    if planted:
        k = next(p for p in a if "levels" in p)
        a[k] *= 1 + 1e-6
    keys_ok = set(a) == set(b)
    worst, wk = 0.0, None
    for k in set(a) & set(b):
        x, y = a[k], b[k]
        if math.isnan(x) and math.isnan(y):
            continue
        # relative to the leaf, with an absolute floor of 1e-12 for derived differences
        # (e.g. extrap - exact ~ 5e-8, where 1e-14 of run-to-run ARPACK noise is 1e-7 relative)
        d = abs(x - y) / (max(abs(x), abs(y)) + 1e-3)
        if d > worst:
            worst, wk = d, k
    report(f"C1 room {room}: rerun JSON vs return JSON, {len(a)} leaves, same keys {keys_ok}, "
           f"max rel {worst:.1e} at {wk}", keys_ok and worst < 1e-9, planted)
    return {"n_leaves": len(a), "same_keys": keys_ok, "max_rel": worst, "where": wk}


out["C1_a"], out["C1_b"] = c1("a"), c1("b")
out["C1_planted_fired"] = c1("a", True)["max_rel"] >= 1e-9


def num(s):
    return float(s)


def ndig(s):
    """half a unit in the last printed place of s (fixed or exponent form)."""
    m = re.match(r"^[+\-−]?(\d*)\.?(\d*)(?:e([+\-]?\d+))?$", s)
    dec = len(m.group(2))
    ex = int(m.group(3)) if m.group(3) else 0
    return 0.5 * 10 ** (ex - dec) * 1.0000001


def sections(text, pat):
    parts = re.split(pat, text)
    return {parts[i].strip(): parts[i + 1] for i in range(1, len(parts) - 1, 2)}


def c2(room, planted=False):
    text = open(os.path.join(HERE, f"room_{room}", "RETURN.md")).read()
    js = json.load(open(os.path.join(HERE, f"room_{room}", "results.json")))["t4"]
    t4 = text[text.index("T4"):]
    if room == "a":
        secs = sections(t4, r"#### W/R = (\S+)\n")
        rowpat = r"^\| (\d) \| \d+ \| \d+ \| \d+ \| (.*) \|$"
        estpat = r"p = ([\d.]+), extrapolated bottom = ([\d.]+), error estimate = ([\d.e+\-]+)"
    else:
        secs = sections(t4[:t4.index("**Extra (labelled):**")], r"### W = (\S+)\n")
        rowpat = r"^\| (\d) \| (.*) \|$"
        estpat = r"p = ([\d.]+), extrap = ([\d.]+), err = ([\d.e+\-]+)"
    bad, n = [], 0
    for lab in LABS:
        sec = secs[lab if room == "a" else lab.replace("pi", "π")]
        rows = re.findall(rowpat, sec, flags=re.M)
        if len(rows) != 5:
            bad.append((lab, "rows", len(rows)))
            continue
        for lev, body in rows:
            vals = [v.strip() for v in body.split("|")]
            for q, v in enumerate(vals):
                if planted and lab == "7/5" and lev == "4" and q == 0:
                    v = v[:-1] + str((int(v[-1]) + 1) % 10)
                n += 1
                if abs(num(v) - js[lab]["levels"][int(lev)][q]) > ndig(v):
                    bad.append((lab, lev, q, v, js[lab]["levels"][int(lev)][q]))
        m = re.search(estpat, sec)
        for s, key in zip(m.groups(), ("p", "extrap", "err")):
            n += 1
            if abs(num(s) - js[lab][key]) > ndig(s):
                bad.append((lab, key, s, js[lab][key]))
    report(f"C2 room {room}: {n} printed T4 numbers vs JSON, mismatches {bad[:3]}", not bad, planted)
    return {"n": n, "mismatches": bad}


out["C2_a"], out["C2_b"] = c2("a"), c2("b")
out["C2_planted_fired"] = bool(c2("a", True)["mismatches"])


def c3(room, planted=False):
    text = open(os.path.join(HERE, f"room_{room}", "RETURN.md")).read()
    js = json.load(open(os.path.join(HERE, f"room_{room}", "results.json")))
    bad, n = [], 0
    for lab0 in LABS:
        lab = lab0 if room == "a" else lab0.replace("pi", "π")
        if room == "a":
            m = re.search(r"^\| " + re.escape(lab) + r" \| [\d.]+ \| ([^|]+) \| \(1", text, flags=re.M)
            vals = [v.strip() for v in m.group(1).split(",")]
            want = js["t3_closed_form"][lab0]["lowest6"]
            m5 = re.search(r"^\| " + re.escape(lab) + r" \| ([+\-]\d+\.\d+) \| [+0−] \|", text, flags=re.M)
            jv, jw = m5.group(1), js["t5"][lab0]["jacobi_lowest"]
        else:
            m = re.search(r"^\| " + re.escape(lab) + r" \| [^|]+ \| (.*) \|$", text, flags=re.M)
            vals = re.findall(r"([\d.]+) \(\d+,\d+\)", m.group(1))
            want = [e["value"] for e in js["t3_t5_exact"][lab0]["eigs"]]
            m5 = re.search(r"^\| " + re.escape(lab) + r" \| [^|]+ \| ([−\-]?[\d.]+) \| ", text[text.index("## T5"):], flags=re.M)
            jv, jw = m5.group(1).replace("−", "-"), js["t3_t5_exact"][lab0]["J_lowest"]
        if planted and lab == "1":
            vals[0] = vals[0][:-1] + str((int(vals[0][-1]) + 1) % 10)
        for v, w in zip(vals, want):
            n += 1
            tol = ndig(v) if "." in v else 1e-9
            if abs(num(v) - w) > tol:
                bad.append((lab, v, w))
        n += 1
        if abs(num(jv) - jw) > (ndig(jv) if "." in jv else 1e-9):
            bad.append((lab, "J", jv, jw))
    report(f"C3 room {room}: {n} printed T3/T5 numbers vs JSON, mismatches {bad[:3]}", not bad, planted)
    return {"n": n, "mismatches": bad}


out["C3_a"], out["C3_b"] = c3("a"), c3("b")
out["C3_planted_fired"] = bool(c3("a", True)["mismatches"])
out["lines"] = lines
with open(os.path.join(HERE, "s2_compare.json"), "w") as fh:
    json.dump(out, fh, indent=1)
