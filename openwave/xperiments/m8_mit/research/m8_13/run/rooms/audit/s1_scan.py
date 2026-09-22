"""Global scan of rhat6 over the unit sphere of C^7 (14 real dims).

Attacks C1 (max value), C3 (min value) and produces the candidate pool for the
orbit tests of C2/C3.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
from scipy.optimize import minimize
import core as C
from grad import f_and_grad_real

HERE = pathlib.Path(__file__).parent
rng = np.random.default_rng(20260922)
NSTART = 4000


def optimize(sign, nstart=NSTART):
    def fun(y):
        v, g = f_and_grad_real(y)
        return sign * v, sign * g

    vals, pts = [], []
    for k in range(nstart):
        r = k % 5
        if r == 0:
            x = rng.normal(size=14)
        elif r == 1:
            x = rng.normal(size=14) * (rng.random(14) < 0.35)
        elif r == 2:
            x = np.concatenate([rng.normal(size=7), np.zeros(7)])
        elif r == 3:
            x = rng.normal(size=14) * 10 ** rng.uniform(-3.0, 1.0, size=14)
        else:
            x = np.zeros(14)
            idx = rng.choice(14, size=rng.integers(1, 4), replace=False)
            x[idx] = rng.normal(size=len(idx))
        if np.linalg.norm(x) < 1e-9:
            x = rng.normal(size=14)
        x /= np.linalg.norm(x)
        res = minimize(fun, x, jac=True, method="L-BFGS-B",
                       options=dict(maxiter=3000, ftol=1e-18, gtol=1e-14))
        y = res.x / np.linalg.norm(res.x)
        # polish with a few extra restarts from the result
        for _ in range(2):
            res = minimize(fun, y, jac=True, method="L-BFGS-B",
                           options=dict(maxiter=3000, ftol=1e-18, gtol=1e-14))
            y = res.x / np.linalg.norm(res.x)
        vals.append(f_and_grad_real(y)[0])
        pts.append(y[:7] + 1j * y[7:])
    return np.array(vals), np.array(pts)


out = []
P = out.append

for sign, name in ((-1, "MAX"), (1, "MIN")):
    vals, pts = optimize(sign)
    best = vals.max() if sign < 0 else vals.min()
    P(f"--- {name}: {NSTART} starts ---")
    P(f"best found            = {best!r}")
    P(f"463/924               = {463/924!r}")
    P(f"1/924                 = {1/924!r}")
    P(f"best - target         = {best - (463/924 if sign<0 else 1/924):.3e}")
    # distinct converged levels
    sv = np.sort(vals)
    levels, cur = [], [sv[0]]
    for v in sv[1:]:
        if v - cur[-1] > 1e-7:
            levels.append((np.mean(cur), len(cur)))
            cur = [v]
        else:
            cur.append(v)
    levels.append((np.mean(cur), len(cur)))
    P(f"distinct converged levels ({len(levels)}):")
    for lv, n in (levels if sign < 0 else levels):
        P(f"    {lv!r}   count={n}")
    np.save(HERE / f"pool_{name}.npy", pts)
    np.save(HERE / f"vals_{name}.npy", vals)

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s1.txt").write_text(txt + "\n")
