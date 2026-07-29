from mutate import run
import json
M=[]
def add(tag,desc,patches,groups=("2T",)):
    r=run(tag,patches,groups)
    M.append({"tag":tag,"description":desc,"result":r[:340]})
    print("\n[%s] %s\n   -> %s"%(tag,desc,r[:340]))

add("M1","corrupt one branching multiplicity Mult[5][0]+=1 (dim_check / gram_matches must go red)",
    [("    G.Mult = Mult\n","    G.Mult = Mult\n    Mult[5][0] += 1\n")])

add("M2","drop the complex conjugate in the isotypic projector (assert rk==ds must fire)",
    [("        cs = np.array([complex(mp.conj(G.chartab[s][cls_of[i]])) for i in range(N)])",
      "        cs = np.array([complex(G.chartab[s][cls_of[i]]) for i in range(N)])")])

add("M3","corrupt one adjacency entry A[1][2]+=1 (A_symmetric / A_rowsum_ok must go red)",
    [("    G.A = A\n","    G.A = A\n    A[1][2] += 1\n")])

add("M4","truncate the closed group by one element (does closed_under_multiplication ever report False?)",
    [("        if len(elems) > cap:\n            raise RuntimeError(\"closure blew past cap\")\n    return elems",
      "        if len(elems) > cap:\n            raise RuntimeError(\"closure blew past cap\")\n    return elems[:-1]")])

add("M5","replace chi_(E_m) by a NON-self-dual class function (does the convA==convB check notice?)",
    [("            chiE = [(m_ - 1) * chi[m_, i] + (m_ + 1) * chi[m_ - 2, i] for i in range(nc)]",
      "            chiE = [((m_ - 1) * chi[m_, i] + (m_ + 1) * chi[m_ - 2, i]) * G.chartab[1][i] for i in range(nc)]")])

add("M6","swap the tau/tau* convention in mu (convention A now uses conj(chi_tau)): does ANY check go red?",
    [("                sa += G.class_sizes[i] * chiE[i] * ch[i]\n                sb += G.class_sizes[i] * chiE[i] * mp.conj(ch[i])",
      "                sa += G.class_sizes[i] * chiE[i] * mp.conj(ch[i])\n                sb += G.class_sizes[i] * chiE[i] * ch[i]")])

add("M7","T4: use Lambda^2 instead of Sym^2 for S(rho) (assert sum(dim*mult)==3 must fire)",
    [("        Schar = [(G.chartab[s][i] ** 2 + G.chartab[s][G.class_sq[i]]) / 2",
      "        Schar = [(G.chartab[s][i] ** 2 - G.chartab[s][G.class_sq[i]]) / 2")])

add("M8","loosen SVD_TOL from 1e-8 to 1e-20 (does any check catch the inflated ranks?)",
    [("SVD_TOL = 1e-8           # stated tolerance for the numerical-rank decision",
      "SVD_TOL = 1e-20          # MUTATED")])

add("M9","tighten SVD_TOL to 0.5 (rank decision should be unaffected: sv are 1 vs ~1e-15)",
    [("SVD_TOL = 1e-8           # stated tolerance for the numerical-rank decision",
      "SVD_TOL = 0.5            # MUTATED")])

add("M10","off-by-one in the symmetric-power level: chi_a built as Sym^(a+1) (many checks must go red)",
    [("            for k in range(a + 1):\n                s += lam ** (a - 2 * k)",
      "            for k in range(a + 2):\n                s += lam ** (a + 1 - 2 * k)")])
json.dump(M,open("_mutations.json","w"),indent=1)
