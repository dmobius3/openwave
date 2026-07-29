import subprocess, sys, textwrap, json
BASE=open("mut/base.py").read()
R = textwrap.dedent('''
    import json, traceback
    import %(mod)s as S
    log=[]
    try:
        gs={n:(g,o) for n,g,o in S.build_groups()}
        out={}
        for nm in %(groups)r:
            G=S.analyse_group(nm,gs[nm][0],gs[nm][1],log)
            out[nm]={"irrep_char_err":G.irrep_char_err,"irrep_hom_err":G.irrep_hom_err,
                     "irrep_unitary_err":G.irrep_unitary_err,
                     "mu_agree_all":all(r["agree"] for t in G.mu.values() for r in t["rows"]),
                     "mu_rows_disagreeing":[(l,r["m"],r["char_tau"],r["char_tau_conj"],r["svd_rank"])
                        for l,t in G.mu.items() for r in t["rows"] if not r["agree"]][:6],
                     "T4":[(r["rho_sigma"],r["q"],r["q_squared"]) for r in G.T4],
                     "T5q":G.T5["q"]}
        print("RESULT "+json.dumps(out))
    except BaseException as e:
        print("EXCEPTION "+type(e).__name__+": "+str(e)[:250])
''')
def run(tag,patches,groups):
    src=BASE
    for o,n in patches:
        assert o in src, (tag,o[:50]); src=src.replace(o,n,1)
    m="mut2_%s"%tag
    open("mut/%s.py"%m,"w").write(src); open("mut/r2_%s.py"%tag,"w").write(R%{"mod":m,"groups":list(groups)})
    r=subprocess.run([sys.executable,"mut/r2_%s.py"%tag],capture_output=True,text=True)
    ls=[l for l in (r.stdout+r.stderr).splitlines() if l.startswith(("RESULT","EXCEPTION"))]
    return ls[-1] if ls else "NO OUTPUT: "+ (r.stdout+r.stderr)[-250:]

O=[]
def add(t,d,p,g=("2T",)):
    res=run(t,p,g); O.append({"tag":t,"description":d,"result":res[:500]})
    print("\n[%s] %s\n   -> %s"%(t,d,res[:500]))

add("M11","loosen SVD_TOL ONLY inside the T3 mu block (irrep extraction untouched)",
    [("                sv = np.linalg.svd(P, compute_uv=False)\n                r_ = int(np.sum(sv > SVD_TOL))",
      "                sv = np.linalg.svd(P, compute_uv=False)\n                r_ = int(np.sum(sv > 1e-20))")],("2T",))
add("M12","corrupt ONE svd rank in the T3 mu block (does row['agree'] go red?)",
    [("                rk += rep * r_","                rk += rep * r_ + (1 if (m_==5 and a==m_) else 0)")],("2T",))
add("M13","is irrep_matrix_character_error ever asserted? mutate the extracted matrices by a phase",
    [("        mats = np.array([B.conj().T @ U[found][i] @ B for i in range(N)])",
      "        mats = np.array([(1j) * (B.conj().T @ U[found][i] @ B) for i in range(N)])")],("2T",))
json.dump(O,open("_mutations2.json","w"),indent=1)
