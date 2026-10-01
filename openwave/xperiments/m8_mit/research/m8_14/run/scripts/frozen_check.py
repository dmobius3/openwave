"""M8.14 maintainer route, run before the rooms open.

1. Parse the FROZEN VALUES table from the pinned task doc and recompute every entry at 30 digits.
2. The F1b reference for the seam diagnostic: the lowest six of (m a0 + l)(m a0 + l + 1), m >= 1, l >= 0.
3. score_width / score_f2 / score_controls: the frozen outcome rules, the control gates.
4. Mutation arms: every scoring line must be able to fail.

Usage: python frozen_check.py <task_doc.md> <fe_out.json>
"""

import json
import re
import sys

import mpmath as mp

mp.mp.dps = 30
TOL = mp.mpf("1e-3")
WIDTHS = {"1/4": mp.mpf(1) / 4, "1/2": mp.mpf(1) / 2, "1": mp.mpf(1), "7/5": mp.mpf(7) / 5,
          "3/2": mp.mpf(3) / 2, "pi/2": mp.pi / 2, "17/10": mp.mpf(17) / 10}
F2_WIDTHS = ["1/4", "1/2", "1", "7/5", "3/2"]
fails = []


def check(label, ok):
    print(("PASS " if ok else "FAIL ") + label)
    if not ok:
        fails.append(label)
    return ok


def alpha0(W):
    return mp.pi / (2 * W)


def target(W):
    a = alpha0(W)
    return a * (a + 1)


def f1b_lowest(W, k=6):
    a = alpha0(W)
    vals = sorted((m * a + l) * (m * a + l + 1) for m in range(1, 40) for l in range(0, 40))
    return vals[:k]


def parse_frozen(doc):
    sec = doc[doc.index("## FROZEN VALUES"):doc.index("## FEASIBILITY")]
    rows = {}
    for line in sec.splitlines():
        m = re.match(r"\| (\S+)(?: \((N1|N2)\))? \| ([\d.]+) \| ([\d.]+) \| ([+-][\d.]+) \|", line)
        if m:
            rows[m.group(1).replace("π", "pi")] = (m.group(3), m.group(4), m.group(5))
    return rows


def frozen_ok(rows):
    ok = len(rows) == 7 and set(rows) == set(WIDTHS)
    for name, (a_s, b_s, j_s) in rows.items():
        W = WIDTHS[name]
        a = alpha0(W)
        ok &= abs(mp.mpf(a_s) - a) <= mp.mpf("5e-7")
        ok &= abs(mp.mpf(b_s) - a * (a + 1)) <= mp.mpf("5e-7")
        ok &= abs(mp.mpf(j_s) - (a - 1) * (a + 2)) <= mp.mpf("5e-7")
        ok &= abs(a * (a + 1) - 2 - (a - 1) * (a + 2)) < mp.mpf("1e-25")
    return ok


def score_width(name, extrap, err, p):
    """Frozen outcome at one width: Reproduced / Contradicted / Unresolved."""
    t = target(WIDTHS[name])
    extrap, err = mp.mpf(extrap), mp.mpf(err)
    miss = abs(extrap - t)
    if p is None or not mp.isfinite(mp.mpf(p)) or mp.mpf(p) <= 0:
        converged = False
    else:
        converged = True
    if converged and miss > TOL * t and miss > 10 * err:
        return "Contradicted"
    if converged and miss <= TOL * t and err <= TOL * t:
        return "Reproduced"
    return "Unresolved"


def score_f2(ret):
    per = {n: score_width(n, ret[n]["extrap"], ret[n]["err"], ret[n]["p"]) for n in F2_WIDTHS}
    if any(v == "Contradicted" for v in per.values()):
        return "Contradicted", per
    if all(v == "Reproduced" for v in per.values()):
        return "Reproduced", per
    return "Unresolved", per


def score_controls(ret):
    n1 = ret["pi/2"]
    n1_ok = abs(mp.mpf(n1["extrap"]) - 2) <= TOL * 2
    n2 = ret["17/10"]
    jb = mp.mpf(n2["extrap"]) - 2
    n2_ok = jb < 0 and abs(jb) > 10 * mp.mpf(n2["err"])
    return n1_ok, n2_ok, float(jb)


def fe_to_ret(fe):
    return {n: {"extrap": fe[n]["p_est"]["extrap"], "err": fe[n]["p_est"]["err"], "p": fe[n]["p"],
                "levels": fe[n]["levels"]} for n in fe}


def main():
    doc = open(sys.argv[1]).read()
    rows = parse_frozen(doc)
    check(f"frozen table parsed, {len(rows)} rows", len(rows) == 7)
    check("frozen table recomputes at 30 digits (alpha0, bottom, J, J = bottom - 2)", frozen_ok(rows))
    # arm: a one-unit change in the sixth decimal must fail
    bad = dict(rows)
    a_s, b_s, j_s = bad["7/5"]
    bad["7/5"] = (a_s, b_s[:-1] + str((int(b_s[-1]) + 1) % 10), j_s)
    check("ARM frozen: sixth-decimal change at 7/5 fires", not frozen_ok(bad))

    print("\nF1b lowest six (seam-diagnostic reference):")
    ref = {}
    for n, W in WIDTHS.items():
        ref[n] = [float(v) for v in f1b_lowest(W)]
        print(f"  {n:6s}", " ".join(f"{v:.6f}" for v in ref[n]))
    check("F1b bottom equals alpha0(alpha0+1) at every width",
          all(abs(ref[n][0] - float(target(WIDTHS[n]))) < 1e-12 for n in WIDTHS))
    check("N1 reference is the hemisphere's 2, 6, 6, 12, 12, 12",
          max(abs(a - b) for a, b in zip(ref["pi/2"], [2, 6, 6, 12, 12, 12])) < 1e-12)

    fe = json.load(open(sys.argv[2]))
    ret = fe_to_ret(fe)
    f2, per = score_f2(ret)
    print("\nreviewer instrument, per width:", per)
    check("reviewer instrument scores F2 Reproduced", f2 == "Reproduced")
    n1, n2, jb = score_controls(ret)
    check(f"reviewer instrument N1 within tolerance", n1)
    check(f"reviewer instrument N2 J bottom {jb:.6f} negative beyond 10x err", n2)

    print("\nscoring arms (each must fire):")
    import copy
    m = copy.deepcopy(ret); m["1"]["extrap"] = float(target(WIDTHS["1"])) * 1.01
    check("ARM 1% miss at W=1 with small err -> Contradicted", score_f2(m)[0] == "Contradicted")
    m = copy.deepcopy(ret); m["1"]["extrap"] = float(target(WIDTHS["1"])) * 1.01; m["1"]["err"] = 0.01
    check("ARM 1% miss within 10x err -> Unresolved", score_f2(m)[0] == "Unresolved")
    m = copy.deepcopy(ret); m["3/2"]["err"] = 0.01
    check("ARM err estimate above tolerance -> Unresolved", score_f2(m)[0] == "Unresolved")
    m = copy.deepcopy(ret); m["1/4"]["p"] = -0.5
    check("ARM nonpositive p -> Unresolved", score_f2(m)[0] == "Unresolved")
    m = copy.deepcopy(ret); m["1/2"]["extrap"] = float(target(WIDTHS["1/2"])) * (1 + 9e-4)
    check("ARM miss 9e-4 relative stays Reproduced (boundary inside)", score_f2(m)[0] == "Reproduced")
    m = copy.deepcopy(ret); m["pi/2"]["extrap"] = 2.01
    check("ARM N1 off by 0.01 -> N1 fails", not score_controls(m)[0])
    m = copy.deepcopy(ret); m["17/10"]["extrap"] = 2.05
    check("ARM N2 positive J bottom -> N2 fails", not score_controls(m)[1])
    m = copy.deepcopy(ret); m["17/10"]["err"] = 0.05
    check("ARM N2 err 0.05 (10x > |J|) -> N2 fails", not score_controls(m)[1])

    json.dump({"f1b_lowest6": ref, "targets": {n: float(target(W)) for n, W in WIDTHS.items()}},
              open("maintainer_reference.json", "w"), indent=1)
    print(f"\n{len(fails)} failed")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
