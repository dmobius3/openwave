"""M8.14 compare adapter: score a room's T4 block by the frozen rules, every field derived from the return.

Usage: python compare.py <results.json> [--auditor stage1.json] [--selftest]

Per width, from the room's own "t4" block:
  - the room's reported p, extrap, err, and the same three recomputed from its own levels by the
    protocol (p from the last three bottoms, Richardson from the last two, err = |extrap - finest}),
    with both readings of the Richardson exponent (estimated p, and p = 2);
  - the frozen per-width outcome on the room's REPORTED numbers;
  - the seam diagnostic: the finest level's six against F1b's six, and against the auditor's six.
Then F2 over the five widths and the two control gates.
"""

import copy
import json
import math
import sys

from frozen_check import F2_WIDTHS, WIDTHS, f1b_lowest, score_controls, score_f2, score_width

ALL = list(WIDTHS)


def recompute(levels):
    b = [lv[0] for lv in levels]
    if len(b) != 5:
        return None
    d1, d2 = b[2] - b[3], b[3] - b[4]
    p = math.log2(d1 / d2) if d1 * d2 > 0 else float("nan")
    ext_p = b[4] + (b[4] - b[3]) / (2 ** p - 1) if math.isfinite(p) and p > 0 else float("nan")
    ext_2 = b[4] + (b[4] - b[3]) / 3.0
    return {"p": p, "ext_p": ext_p, "err_p": abs(ext_p - b[4]), "ext_2": ext_2, "err_2": abs(ext_2 - b[4]),
            "finest": b[4]}


def compare(ret, aud=None, quiet=False):
    t4 = ret.get("t4", {})
    out = {"missing": [n for n in ALL if n not in t4], "widths": {}}
    for n in ALL:
        if n not in t4:
            continue
        r = t4[n]
        lv = r.get("levels")
        rc = recompute(lv) if lv else None
        rec = {"reported": {k: r.get(k) for k in ("p", "extrap", "err")}, "recomputed": rc}
        if rc:
            rec["reading"] = ("estimated p" if abs(rc["ext_p"] - r["extrap"]) <= 1e-9 * abs(r["extrap"])
                              else "p = 2" if abs(rc["ext_2"] - r["extrap"]) <= 1e-9 * abs(r["extrap"])
                              else "NEITHER")
            rec["p_consistent"] = abs(rc["p"] - r["p"]) <= 1e-6
            rec["err_consistent"] = (abs(rc["err_p"] - r["err"]) <= 1e-9 * abs(r["extrap"]) or
                                     abs(rc["err_2"] - r["err"]) <= 1e-9 * abs(r["extrap"]))
            fin = lv[-1][:6]
            ref = [float(v) for v in f1b_lowest(WIDTHS[n])]
            rec["seam_vs_f1b"] = [(a - b) / b for a, b in zip(fin, ref)]
            if aud and n in aud.get("v2", {}):
                al = aud["v2"][n]["lowest"][:6]
                rec["seam_vs_auditor"] = [(a - b) / b for a, b in zip(fin, al)]
        if n in F2_WIDTHS:
            rec["outcome"] = score_width(n, r["extrap"], r["err"], r["p"])
        out["widths"][n] = rec
    if not out["missing"]:
        out["F2"] = score_f2(t4)[0]
        n1, n2, jb = score_controls(t4)
        out["N1"], out["N2"], out["N2_J_bottom"] = n1, n2, jb
    else:
        out["F2"] = "INCOMPLETE RETURN"
    if not quiet:
        print(json.dumps(out, indent=1, default=float))
    return out


def selftest(ret):
    ok = True

    def arm(label, cond):
        nonlocal ok
        print(("FIRES " if cond else "DEAD  ") + label)
        ok &= cond

    base = compare(ret, quiet=True)
    arm("base return scores Reproduced, N1 and N2 pass",
        base["F2"] == "Reproduced" and base["N1"] and base["N2"])
    m = copy.deepcopy(ret); del m["t4"]["7/5"]
    arm("missing width -> INCOMPLETE RETURN, never a pass", compare(m, quiet=True)["F2"] == "INCOMPLETE RETURN")
    m = copy.deepcopy(ret); m["t4"]["1"]["extrap"] *= 1.01
    r = compare(m, quiet=True)
    arm("1% shift of the reported extrap -> Contradicted and reading NEITHER",
        r["F2"] == "Contradicted" and r["widths"]["1"]["reading"] == "NEITHER")
    m = copy.deepcopy(ret); m["t4"]["1"]["p"] = 1.5
    arm("reported p not from its levels -> p_consistent False", not compare(m, quiet=True)["widths"]["1"]["p_consistent"])
    m = copy.deepcopy(ret); m["t4"]["1"]["levels"][-1][1] = m["t4"]["1"]["levels"][-1][0]
    r = compare(m, quiet=True)
    arm("seamless-like second level (repeat of the bottom) shows in seam_vs_f1b[1]",
        abs(r["widths"]["1"]["seam_vs_f1b"][1]) > 0.3)
    m = copy.deepcopy(ret); m["t4"]["1"]["err"] *= 3
    arm("reported err not from its levels -> err_consistent False", not compare(m, quiet=True)["widths"]["1"]["err_consistent"])
    print("selftest:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    ret = json.load(open(sys.argv[1]))
    aud = None
    if "--auditor" in sys.argv:
        aud = json.load(open(sys.argv[sys.argv.index("--auditor") + 1]))
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest(ret) else 1)
    compare(ret, aud)
