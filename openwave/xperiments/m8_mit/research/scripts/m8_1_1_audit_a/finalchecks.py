import json, math, itertools
import numpy as np
# --- independent float check of agent A's min pairwise element separation + 2I closure claim
PHI=(1+math.sqrt(5))/2
def q2m(a,b,c,d): return np.array([[a+b*1j,c+d*1j],[-c+d*1j,a-b*1j]],dtype=complex)
def close(gens):
    els=[np.eye(2,dtype=complex)]; fr=list(els)
    while fr:
        nw=[]
        for M in fr:
            for g in gens:
                P=M@g
                if not any(np.max(np.abs(P-E))<1e-9 for E in els+nw): nw.append(P)
        els+=nw; fr=nw
    return els
g2I=[q2m(.5,.5,.5,.5), q2m(PHI/2,(1/PHI)/2,.5,0)]
E=close(g2I)
print("2I closure from the SPEC generators (independent float run):",len(E),"elements")
# explicit spec list
def parity(pm):
    return sum(1 for i in range(4) for j in range(i+1,4) if pm[i]>pm[j])%2
qs=[]
for s in (1,-1):
    for pos in range(4):
        v=[0.]*4; v[pos]=float(s); qs.append(tuple(v))
for sg in itertools.product((1,-1),repeat=4): qs.append(tuple(x*0.5 for x in sg))
base=[0.,1.,1/PHI,PHI]
for pm in [pp for pp in itertools.permutations(range(4)) if parity(pp)==0]:
    for s1,s2,s3 in itertools.product((1,-1),repeat=3):
        vals=[0.,s1*base[1],s2*base[2],s3*base[3]]; q=[0.]*4
        for src in range(4): q[pm[src]]=vals[src]/2
        qs.append(tuple(q))
ded=[]
for q in qs:
    if not any(max(abs(q[i]-r[i]) for i in range(4))<1e-9 for r in ded): ded.append(q)
X=[q2m(*q) for q in ded]
matched=sum(1 for m in X if any(np.max(np.abs(m-e))<1e-9 for e in E))
print("explicit spec 120-list: %d elements, %d of them found in the closure -> identical set: %s"%(
      len(X),matched,matched==len(X)==len(E)))
mins={}
def grp_min(els):
    return min(float(np.max(np.abs(els[i]-els[j]))) for i in range(len(els)) for j in range(i+1,len(els)))
print("min pairwise element separation for 2I:", grp_min(E))
allmin=grp_min(E)
for n in range(2,11):
    z=np.exp(2j*np.pi/n); els=[np.diag([z**k,np.conj(z)**k]) for k in range(n)]
    allmin=min(allmin,grp_min(els))
for n in range(2,7):
    z=np.exp(1j*np.pi/n)
    els=close([np.diag([z,np.conj(z)]),np.array([[0,1],[-1,0]],dtype=complex)])
    allmin=min(allmin,grp_min(els))
els=close([q2m(.5,.5,.5,.5),q2m(0,1,0,0)]); allmin=min(allmin,grp_min(els))
els=close([q2m(.5,.5,.5,.5),q2m(0,1,0,0),q2m(1/math.sqrt(2),1/math.sqrt(2),0,0)])
allmin=min(allmin,grp_min(els))
print("min pairwise separation over ALL 18 spec groups (mine):",allmin)
print("agent A reported                                       : 0.43701602444882115")

# --- T7 relation over the widened set
M=json.load(open("_sweep_raw.json"))
bad1=[];bad2=[]
for n,g in M.items():
    for r in g["T4"]:
        if r["T7_k_k_plus_2"]!=r["q"]**2+2*r["q"]: bad1.append((n,r["q"],r["T7_k_k_plus_2"]))
        if r["T7_k_k_plus_2"]==r["q_squared"]: bad2.append((n,"q2 == k(k+2)"))
    t=g["T5"]
    if t["T7_k_k_plus_2"]==t["q_squared"]: bad2.append((n,"trivial q2 == k(k+2)"))
print()
print("T4 twists: k(k+2) == q^2 + 2q  violations:",bad1)
print("any case where q^2 == k(k+2):",bad2)
print("T5 trivial: q^2=%s vs k(k+2)=%s for every group -> never equal"%(
    {g["T5"]["q_squared"] for g in M.values()},{g["T5"]["T7_k_k_plus_2"] for g in M.values()}))
