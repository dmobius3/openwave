from driver import analyse
from part1 import zeta, quat, inv2, p
import part1
# BD_n for n=13..16 and C_n for n=26..32 : where do least_a / e exceed agent A's 14 / 12 ?
rows=[]
for n in (13,14,15,16):
    z=zeta(2*n); gens=[(z,0,0,pow(z,2*n-1,p)),(0,1,p-1,0)]
    _,r,*_=analyse("BD_%d"%n,gens,4*n,verbose=False)
    rows.append(("BD_%d"%n,max(t["least_a"] for t in r["T1"]),max(t["e"] for t in r["T6"]),r["diameter"]))
for n in (26,28,30,32):
    z=zeta(n); gens=[(z,0,0,pow(z,n-1,p))]
    _,r,*_=analyse("C_%d"%n,gens,n,verbose=False)
    rows.append(("C_%d"%n,max(t["least_a"] for t in r["T1"]),max(t["e"] for t in r["T6"]),r["diameter"]))
print("%-8s %-12s %-8s %s"%("group","max least_a","max e","diameter"))
for g,la,e,d in rows:
    flag=[]
    if la>14: flag.append("least_a EXCEEDS agent A JMAX=14")
    if e>12: flag.append("e EXCEEDS agent A MMAX=12")
    print("%-8s %-12d %-8d %-8d %s"%(g,la,e,d,"; ".join(flag)))
