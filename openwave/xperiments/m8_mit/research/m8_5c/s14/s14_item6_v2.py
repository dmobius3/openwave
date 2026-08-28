"""S14 item 6 v2: the AMENDED phi pin, verified alone (scoped re-run per F1). Conventions
exactly as the amended S 7 pins them; two analytic identity arms per the amendment."""
import numpy as np, json, sys
s=6
m=np.arange(s,-s-1,-1,dtype=float)                 # weights DESCENDING +6..-6
J3=np.diag(m)
ap=np.sqrt((s-m[1:])*(s+m[1:]+1))
Jp=np.zeros((13,13)); Jp[np.arange(12),np.arange(1,13)]=ap
Jm=Jp.T.copy()
T1=0.5*(Jp-Jm)+0j; T2=0.5j*(Jp+Jm); T3=1j*J3
omega=np.sqrt(168.0)
j=np.arange(1,14)
w=(1.0+j/13.0)*np.exp(1j*j/3.0)
phi=w/np.linalg.norm(w)
psid=1j*omega*phi
S_cons=np.linalg.norm(phi)*np.linalg.norm(psid)+np.linalg.norm(phi)**2
thr=1e-6*S_cons
M=[float(np.real(np.vdot(psid,T@phi))) for T in (T1,T2,T3)]
live=[bool(abs(x)>=thr) for x in M]
ratio=abs(M[0]/M[1]); tan13=np.tan(1/3)
armR=abs(ratio-tan13)<1e-12
M3_closed=+omega*float(np.sum(m*np.abs(1.0+j/13.0)**2))/float(np.sum(np.abs(1.0+j/13.0)**2))
armC=abs(M[2]-M3_closed)<1e-10
print(f"  amended phi: M = [{M[0]:+.6e}, {M[1]:+.6e}, {M[2]:+.6e}]")
print(f"  liveness (thr {thr:.3e}): {live}")
print(f"  arm R: |M1/M2| = {ratio:.12f} vs tan(1/3) = {tan13:.12f}  "
      f"{'PASS' if armR else 'FAIL'}")
print(f"  arm C: M3 = {M[2]:+.10f} vs closed form {M3_closed:+.10f}  "
      f"{'PASS' if armC else 'FAIL'}")
res=dict(M=M,live=live,S_cons=float(S_cons),threshold=float(thr),
         ratio=float(ratio),tan_one_third=float(tan13),M3_closed=float(M3_closed),
         protocol_sha=open("/tmp/h5.txt").read().strip())
json.dump(res,open("raw/item6_v2_amended.json","w"),indent=1)
ok=all(live) and armR and armC
print(f"  ITEM 6 v2 VERDICT: {'GREEN' if ok else 'RED'}")
sys.exit(0 if ok else 1)
