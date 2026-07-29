"""P8 general rule for e(sigma) and q(tau), plus the truncation boundary probe."""
import json
M=json.load(open("_sweep_raw.json"))
bad=[]
for name,g in M.items():
    n=len(g["dims"]); br=g["branching"]
    for s in range(n):
        occ={aa for aa in range(33) if br[aa][s]}
        pred=next((m_ for m_ in range(2,31) if (m_ in occ) or (m_-2 in occ)),None)
        if pred!=g["T6"][s]["e"]: bad.append(("e",name,s,pred,g["T6"][s]["e"]))
    # trivial
    occ={aa for aa in range(33) if br[aa][0]}
    pred=next((m_ for m_ in range(2,31) if (m_ in occ) or (m_-2 in occ)),None)
    if pred!=g["T5"]["q"]: bad.append(("q_triv",name,pred,g["T5"]["q"]))
print("P8 general rule  e(tau) = min{m>=2 : tau* in V_m or tau* in V_(m-2)}  violations:",bad)
# where does d(sigma) exceed agent A's JMAX=14 / MMAX=12 ?
print()
print("max distance per family (agent A: a<=14, m<=12):")
for name,g in sorted(M.items()):
    print("  %-6s diam=%-3d max_d=%-3d max_least_a=%-3d max_e=%-3d"%(
      name,g["diameter"],max(g["distance_vector"]),
      max(t["least_a"] for t in g["T1"]),max(t["e"] for t in g["T6"])))
