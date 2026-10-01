"""Format results.json into the markdown tables of RETURN.md (writes t4_tables.md). No new numbers."""
import os, json

here = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(here, "results.json")))
t3, t4, t5 = r["t3_closed_form"], r["t4"], r["t5"]
L = []
for lab, v in t4.items():
    L.append(f"#### W/R = {lab}\n")
    L.append("| level | N per side | cells across | ndof | λ1 | λ2 | λ3 | λ4 | λ5 | λ6 |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for i, row in enumerate(v["levels"]):
        L.append(f"| {i} | {16 * 2**i} | {8 * 2**i} | {v['ndof'][i]} | " + " | ".join(f"{x:.10f}" for x in row) + " |")
    L.append(f"\np = {v['p']:.6f}, extrapolated bottom = {v['extrap']:.10f}, error estimate = {v['err']:.3e}\n")
L.append("#### Summary\n")
L.append("| W/R | p | extrapolated bottom | error estimate | closed form ν1(ν1+1) (T3) | extra: extrap − closed form |")
L.append("| --- | --- | --- | --- | --- | --- |")
for lab, v in t4.items():
    L.append(f"| {lab} | {v['p']:.6f} | {v['extrap']:.10f} | {v['err']:.3e} | {t3[lab]['lowest6'][0]:.10f} | "
             f"{r['extra_t4_vs_closed_form'][lab]['extrap_minus_exact']:+.3e} |")
L.append("\n#### Closed-form lowest six (T3), labels (k, n)\n")
L.append("| W/R | ν1 | lowest six λ_{k,n} | labels |")
L.append("| --- | --- | --- | --- |")
for lab, v in t3.items():
    L.append(f"| {lab} | {v['nu1']:.12f} | " + ", ".join(f"{x:.10f}" for x in v["lowest6"]) + " | "
             + " ".join(f"({k},{n})" for k, n in v["lowest6_labels"]) + " |")
L.append("\n#### T5 table\n")
L.append("| W/R | lowest Jacobi eigenvalue ν1(ν1+1) − 2 | sign | from T4 extrap − 2 | negative eigenvalues | nullity |")
L.append("| --- | --- | --- | --- | --- | --- |")
for lab, v in t5.items():
    s = "0" if v["nullity"] and abs(v["jacobi_lowest"]) < 1e-15 else ("+" if v["jacobi_lowest"] > 0 else "−")
    L.append(f"| {lab} | {v['jacobi_lowest']:+.12f} | {s} | {v['jacobi_lowest_from_t4_extrap']:+.9f} | "
             f"{v['negative_count']} | {v['nullity']} |")
open(os.path.join(here, "t4_tables.md"), "w").write("\n".join(L) + "\n")
print("wrote t4_tables.md")
