"""Assemble stage1.json (every number) and table_stage1.md from the result
files; then check that STAGE1.md contains the table verbatim."""
import os

from common import WIDTHS, NEIG, load_json, path, save_json


def main():
    RA = load_json("results_A.json")
    RB = load_json("results_B.json")
    RC = load_json("results_C.json")
    CK = load_json("results_checks.json")

    v2 = {}
    for label, W in WIDTHS:
        b = RB["widths"][label]
        v2[label] = {"W_over_R": W, "lowest": b["best"], "bottom": b["best"][0], "err": b["err"][0],
                     "err_all": b["err"], "labels_k_m": b["levels"][-1]["labels_k_m"],
                     "source": "method B (best), cross-checked by method A"}
    out = {"R": 1, "v2": v2, "method_A": RA, "method_B": RB,
           "closed_form_crosscheck": RC, "checks": CK}
    save_json("stage1.json", out)

    lines = ["| W/R | method | λ1 (bottom) | λ2 | λ3 | λ4 | λ5 | λ6 | err(bottom) |",
             "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for label, W in WIDTHS:
        b = RB["widths"][label]
        a = RA["widths"][label]["richardson"]
        c = RC[label]
        lines.append("| %s | **B (best)** | %s | %.1e |" % (
            label, " | ".join("%.12f" % x for x in b["best"]), b["err"][0]))
        lines.append("| | A extrap. | %s | %.1e |" % (
            " | ".join("%.10f" % r["extrap"] for r in a), a[0]["err_combined"]))
        lines.append("| | A p (last 3) | %s | |" % " | ".join("%.3f" % r["p"] for r in a))
        lines.append("| | closed form (check only) | %s | |" % " | ".join("%.12f" % x for x in c["lowest"]))
        lines.append("| | labels (k,m) | %s | |" % " | ".join("(%d,%d)" % tuple(x) for x in c["labels_k_m"]))
    table = "\n".join(lines) + "\n"
    with open(path("table_stage1.md"), "w") as f:
        f.write(table)
    print(table)

    s1 = path("STAGE1.md")
    if os.path.exists(s1):
        with open(s1) as f:
            txt = f.read()
        print(("PASS" if table in txt else "FAIL") + "  STAGE1.md contains table_stage1.md verbatim")
    else:
        print("FAIL  STAGE1.md not found")


if __name__ == "__main__":
    main()
