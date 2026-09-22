import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np, sympy as sp, core as C
from grad import f_and_grad_real
HERE = pathlib.Path(__file__).parent
vals = np.load(HERE/"vals_MAX.npy"); pts = np.load(HERE/"pool_MAX.npy")
print("level            count   max|tangential grad|   exact guess")
for lev in sorted(set(np.round(vals,9))):
    sel = np.where(np.abs(vals-lev)<1e-8)[0]
    if len(sel)<5: continue
    gn = 0.0
    for i in sel[:200]:
        y = np.concatenate([pts[i].real, pts[i].imag]); y/=np.linalg.norm(y)
        v,g = f_and_grad_real(y)
        g = g - np.dot(g,y)*y          # project out the radial direction
        gn = max(gn, np.linalg.norm(g))
    r = sp.nsimplify(sp.Rational(round(lev*924*35), 924*35))
    print(f"{lev!r:22} {len(sel):5d}   {gn:.3e}   {r} = {sp.N(r,18)}  match={abs(float(r)-lev)<1e-12}")
