"""S14 item 3: measured wall-clocks at the TOP RUNG of each ladder, author-side SCRATCH,
feasibility-only. What is timed: the load-bearing primitives of the separable-Hopf transform
path at production sizes, composed per the protocol's own operation counts. This licenses
FEASIBILITY of the 48-hour ceiling, never a prediction of the Build Unit's in-room cost,
and crossing the ceiling is STOP-QUAL regardless of these numbers.
"""
import numpy as np, time, json, sys

def t_of(f, reps=3):
    ts=[]
    for _ in range(reps):
        t0=time.perf_counter(); f(); ts.append(time.perf_counter()-t0)
    return min(ts)

out={}
# grids: production 4N at N=48 and N=60; monitor 6N at N=60
GRIDS={"agree_N48_4N":(97,193,193), "ctrlB_N60_4N":(121,241,241), "monitor_N60_6N":(181,361,361)}
for name,(nu,k1,k2) in GRIDS.items():
    A=(np.random.standard_normal((nu,k1,k2))+1j*np.random.standard_normal((nu,k1,k2)))
    fft2=t_of(lambda: np.fft.fft2(A,axes=(1,2)))
    cubic=t_of(lambda: (np.abs(A)**2)*A)
    U=np.random.standard_normal((nu,nu))
    umat=t_of(lambda: np.einsum('uv,vjk->ujk',U,A))
    per_eval = 2*fft2 + cubic + 2*umat        # synthesis + analysis + pointwise cubic
    out[name]=dict(nodes=nu*k1*k2, fft2_s=fft2, cubic_s=cubic, umatmul_s=umat,
                   per_projected_cubic_s=per_eval)
    print(f"  {name}: nodes {nu*k1*k2:,}  fft2 {fft2:.2f}s  cubic {cubic:.2f}s  "
          f"u-matmul {umat:.2f}s  per-evaluation {per_eval:.2f}s")

# Control B eigensolve + LDL at top rung: 1354-dim real symmetric
n=1354
S=np.random.standard_normal((n,n)); S=(S+S.T)/2
eigh=t_of(lambda: np.linalg.eigh(S), reps=2)
import scipy.linalg as sla
ldl=t_of(lambda: sla.ldl(S), reps=2)
out["ctrlB_eigh_1354"]=eigh; out["ctrlB_ldl_1354"]=ldl
print(f"  Control B: eigh(1354) {eigh:.2f}s, ldl(1354) {ldl:.2f}s")

# composition, stated as the protocol counts them:
# agreement gate: 40 fields x 2 routes at each of 8 rungs, bounded by top-rung cost
agree = 40*2*8*out["ctrlB_N60_4N"]["per_projected_cubic_s"]
# monitor: one 6N reading per field per rung, bounded by the N=60 monitor cost
monitor = 40*8*out["monitor_N60_6N"]["per_projected_cubic_s"]
# Control B: 4 eta x ~12 lattice classes x <=30 Newton iters, each ~ one projected cubic
# + one dense solve; eigen at 4x12 accepted points; inertia 2 LDL per point
ctrlB = 4*12*30*(out["ctrlB_N60_4N"]["per_projected_cubic_s"] + 0.3) + 4*12*(eigh+2*ldl)
total = agree+monitor+ctrlB
out["composed_hours"]=dict(agreement=agree/3600, monitor=monitor/3600,
                           controlB=ctrlB/3600, total=total/3600, ceiling=48.0)
print(f"  composed (bounding everything by top-rung cost): agreement {agree/3600:.2f} h, "
      f"monitor {monitor/3600:.2f} h, Control B {ctrlB/3600:.2f} h")
print(f"  TOTAL {total/3600:.2f} h against the 48 h ceiling "
      f"(margin {48/(total/3600):.0f}x); FEASIBILITY-ONLY")
json.dump(out,open("raw/item3_wallclocks.json","w"),indent=1,default=float)
ok = total/3600 < 24
print(f"  ITEM 3 VERDICT: {'GREEN' if ok else 'RED'}")
sys.exit(0 if ok else 1)
