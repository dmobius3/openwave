import sys, json, functools, random, time; sys.path.insert(0,'.')
import numpy as np
from complex import *
from kernel import left_kernel
import sat
pres=json.load(open('pres.json')); s,t = pres['s_id'], pres['t_id']
R1 = [(0,1)]*3 + [(1,-1),(0,-1)]*2          # s^3 (st)^-2
R2 = [(1,1)]*5 + [(1,-1),(0,-1)]*2          # t^5 (st)^-2
M2=[[fox(r,k,[s,t]) for k in (0,1)] for r in (R1,R2)]
def rowM(x,M):
    return [functools.reduce(radd,(rmul(x[i],M[i][j]) for i in range(len(M))),ZR)
            for j in range(len(M[0]))]
def flat(x):
    v=[0]*(2*N)
    for j in range(2):
        for g,c in x[j].items(): v[j*N+g]+=c
    return v
def unflat(v):
    return [{g:c for g,c in enumerate(v[:N]) if c},{g:c for g,c in enumerate(v[N:]) if c}]
A=[]
for i in range(2):
    for g in range(N):
        x=[ZR,ZR]; x[i]={g:1}
        A.append(flat(rowM(x,M2)))
KER,rk = left_kernel(A)
print(f"  symmetric presentation: rank {rk}, saturated kernel basis {len(KER)}, "
      f"max entry {max(abs(v) for r in KER for v in r)}", flush=True)

def rank_p_np(rows,p):
    A=np.array(rows,dtype=np.int64)%p; R,C=A.shape; rk=0
    for c in range(C):
        nz=np.nonzero(A[rk:,c])[0]
        if not nz.size: continue
        r=rk+nz[0]; A[[rk,r]]=A[[r,rk]]
        A[rk]=(A[rk]*pow(int(A[rk,c]),p-2,p))%p
        col=A[:,c].copy(); col[rk]=0
        A=(A-np.outer(col,A[rk]))%p; rk+=1
        if rk==R: break
    return rk
SCREEN=(2,3,5,7,11,13,17,19,23,29,31,37)
def screen(x):
    rows=sat.zmat(x,MUL)
    for p in SCREEN:
        if rank_p_np(rows,p)<119: return None
    return rows

t0=time.time(); tried=0; found=None
cands=[]
for i,b in enumerate(KER): cands.append(('basis %d'%i, b))
rng=random.Random(20260803)
for k in range(6000):
    supp=rng.sample(range(len(KER)), rng.randint(2,5))
    v=[0]*(2*N)
    for j in supp:
        c=rng.choice((-1,1))
        for i in range(2*N): v[i]+=c*KER[j][i]
    cands.append((f'sparse {k}', v))
for name,v in cands:
    tried+=1
    rows=screen(unflat(v))
    if rows is None: continue
    print(f"    {name}: passed the {len(SCREEN)}-prime screen", flush=True)
    verdict,why = sat.saturated(rows, random.Random(tried))
    print(f"      -> saturated = {verdict}: {why}", flush=True)
    if verdict is True:
        x=unflat(v)
        assert all(not p for p in rowM(x,M2)), "not in ker d2"
        json.dump({"M3":[{str(a):b for a,b in x[0].items()},{str(a):b for a,b in x[1].items()}],
                   "s_id":s,"t_id":t,"relators":["s^3 (st)^-2","t^5 (st)^-2"],
                   "found_as":name,"seed":20260803,
                   "acceptance":"im(d3) = ker(d2) as integral lattices; gcd of maximal minors = 1"},
                  open("m3_v2.json","w"))
        found=name; break
print(f"\n  tried {tried} candidates in {time.time()-t0:.1f}s -> "
      f"{'FOUND: '+found if found else 'none saturated'}")
