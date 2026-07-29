import subprocess, sys, textwrap, json
BASE=open("mut/base.py").read()
R=textwrap.dedent('''
    import json
    import %(mod)s as S
    log=[]
    try:
        gs={n:(g,o) for n,g,o in S.build_groups()}
        out={}
        for nm in %(groups)r:
            G=S.analyse_group(nm,gs[nm][0],gs[nm][1],log)
            out[nm]={"gram_matches":G.gram_matches,"dim_check":G.dim_check,
                     "sum_dim_sq_ok":G.sum_dim_sq_ok,"A_symmetric":G.A_symmetric,
                     "A_rowsum_ok":G.A_rowsum_ok,
                     "dist":G.dist,"T1_least_a":[t["least_a"] for t in G.T1],
                     "T6_e":[t["e"] for t in G.T6]}
        print("RESULT "+json.dumps(out))
    except BaseException as e:
        print("EXCEPTION "+type(e).__name__+": "+str(e)[:250])
''')
def run(tag,patches,groups):
    src=BASE
    for o,n in patches:
        assert o in src,(tag,o[:60]); src=src.replace(o,n,1)
    m="mut3_%s"%tag
    open("mut/%s.py"%m,"w").write(src); open("mut/r3_%s.py"%tag,"w").write(R%{"mod":m,"groups":list(groups)})
    r=subprocess.run([sys.executable,"mut/r3_%s.py"%tag],capture_output=True,text=True)
    ls=[l for l in (r.stdout+r.stderr).splitlines() if l.startswith(("RESULT","EXCEPTION"))]
    return ls[-1] if ls else "NO OUTPUT: "+(r.stdout+r.stderr)[-250:]
O=[]
def add(t,d,p,g=("2T",)):
    res=run(t,p,g); O.append({"tag":t,"description":d,"result":res[:600]}); print("\n[%s] %s\n   -> %s"%(t,d,res[:600]))

# D7: does the exact cyclotomic Gram check notice a SWAP of two irreducible labels?
add("M14","swap two columns of the multiplicity matrix, i.e. relabel two irreducibles "
          "(does the exact cyclotomic Gram check gram_matches go red?)",
    [("    G.Mult = Mult\n",
      "    G.Mult = Mult\n    for _r in Mult: _r[1], _r[2] = _r[2], _r[1]\n")])

# D5: does agent A's script survive a group whose distances exceed JMAX=14?
add("M15","run agent A's analyse_group on C_30 (max d = 15 > JMAX = 14)",
    [("    for n in range(1, 11):","    for n in (30,):")],("C_30",))
add("M16","run agent A's analyse_group on C_26 (max d = 13 > MMAX = 12, but <= JMAX)",
    [("    for n in range(1, 11):","    for n in (26,):")],("C_26",))
add("M17","run agent A's analyse_group on BD_13 (max d = 13 > MMAX = 12)",
    [("    for n in range(2, 7):","    for n in (13,):")],("BD_13",))
json.dump(O,open("_mutations3.json","w"),indent=1)

add("M14b","2O: swap Mult columns 2,3 -- two 2-dim irreps at DIFFERENT distances (5 and 1). "
          "Does gram_matches stay green while T1 least_a stops equalling d?",
    [("    G.Mult = Mult\n",
      "    G.Mult = Mult\n    for _r in Mult: _r[2], _r[3] = _r[3], _r[2]\n")],("2O",))
add("M18","agent A's script on BD_16 (2-dim irreps out to d = 15 > JMAX = 14)",
    [("    for n in range(2, 7):","    for n in (16,):")],("BD_16",))
add("M19","agent A's script on BD_15 (2-dim irreps out to d = 15 > JMAX = 14)",
    [("    for n in range(2, 7):","    for n in (15,):")],("BD_15",))
json.dump(O,open("_mutations3.json","w"),indent=1)
