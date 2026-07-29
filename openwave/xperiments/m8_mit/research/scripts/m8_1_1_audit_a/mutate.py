"""Mutation harness: patch a copy of agent A's script, run one group, report whether the check goes red."""
import subprocess, sys, os, textwrap
BASE=open("mut/base.py").read()
RUNNER = textwrap.dedent('''
    import traceback, json
    import %(mod)s as S
    log=[]
    try:
        gs={n:(g,o) for n,g,o in S.build_groups()}
        out={}
        for nm in %(groups)r:
            G=S.analyse_group(nm,gs[nm][0],gs[nm][1],log)
            out[nm]={"dim_check":G.dim_check,"gram_matches":G.gram_matches,
                     "A_symmetric":G.A_symmetric,"A_rowsum_ok":G.A_rowsum_ok,
                     "sum_dim_sq_ok":G.sum_dim_sq_ok,"closed":G.closed,
                     "order_ok":G.order_ok,
                     "mu_agree_all":all(r["agree"] for t in G.mu.values() for r in t["rows"]),
                     "mu_convA_vs_convB_all_equal":all(r["char_tau"]==r["char_tau_conj"]
                                                       for t in G.mu.values() for r in t["rows"]),
                     "T5":G.T5,"T4":[(r["rho_sigma"],r["q"],r["q_squared"]) for r in G.T4]}
        print("RESULT "+json.dumps(out))
    except BaseException as e:
        print("EXCEPTION "+type(e).__name__+": "+str(e)[:300])
''')

def run(tag, patches, groups=("2T",), expect_red=None):
    src=BASE
    for old,new in patches:
        assert old in src, (tag,"patch anchor not found: "+old[:60])
        src=src.replace(old,new,1)
    mod="mut_%s"%tag
    open("mut/%s.py"%mod,"w").write(src)
    open("mut/run_%s.py"%tag,"w").write(RUNNER%{"mod":mod,"groups":list(groups)})
    r=subprocess.run([sys.executable,"mut/run_%s.py"%tag],capture_output=True,text=True,cwd=".")
    out=(r.stdout+r.stderr).strip().splitlines()
    line=[l for l in out if l.startswith(("RESULT","EXCEPTION"))]
    return line[-1] if line else ("NO OUTPUT / rc=%d :: "%r.returncode)+ " | ".join(out[-3:])

if __name__=="__main__":
    print("M0 CONTROL (unmutated, 2T):")
    print("   ",run("m0",[])[:400])
