"""Item 1 (explicit intertwiners): eta at every level of each sector, mpmath 50 digits; residual checks.
The one matrix representation sigma(h) := eta3^dag D^3(h) eta3 (eta3 = orthonormal basis of range P_sigma)."""
import pickle, json, random
import mpmath as mp
from math import comb, factorial
from fractions import Fraction
from checks import check, save_log

mp.mp.dps = 50
with open("out/group.pkl", "rb") as f:
    GRP = pickle.load(f)
G, P, dims = GRP["G"], GRP["P"], GRP["dims"]
q1, q2 = GRP["q1"], GRP["q2"]
from su2 import qkey
keys = [qkey(g) for g in G]
i1, i2 = keys.index(qkey(q1)), keys.index(qkey(q2))


def su2_mp(q):
    w, x, y, z = [c.to_mpc().real for c in q]
    I = mp.mpc(0, 1)
    return [[w - I * z, -y - I * x], [y - I * x, w + I * z]]


def D_mp(j, g):
    (g11, g12), (g21, g22) = g
    n = 2 * j
    D = mp.matrix(n + 1, n + 1)
    wm = lambda m: mp.mpf(factorial(j + m) * factorial(j - m)) / factorial(2 * j)
    for m in range(-j, j + 1):
        for k1 in range(j + m + 1):
            c1 = comb(j + m, k1) * g11 ** k1 * g21 ** (j + m - k1)
            for k2 in range(j - m + 1):
                c2 = comb(j - m, k2) * g12 ** k2 * g22 ** (j - m - k2)
                mpr = k1 + k2 - j
                D[mpr + j, m + j] += c1 * c2
    for a in range(-j, j + 1):
        for b in range(-j, j + 1):
            D[a + j, b + j] *= mp.sqrt(wm(a) / wm(b))
    return D


def dagm(A):
    return A.transpose_conj()


def maxabs(A):
    return max(abs(A[i, j]) for i in range(A.rows) for j in range(A.cols))


gs = [su2_mp(g) for g in G]
Dcache = {}
def Dall(j):
    if j not in Dcache:
        Dcache[j] = [D_mp(j, g) for g in gs]
    return Dcache[j]

out = {}
resid = {}
random.seed(12345)
for d in (3, 4):
    Pm = mp.matrix([[P[d][i][k].to_mpc() for k in range(7)] for i in range(7)])
    # orthonormal basis of range P: Gram-Schmidt on columns of P
    cols = []
    for k in range(7):
        v = Pm[:, k]
        for c in cols:
            v = v - c * (dagm(c) * v)[0]
        nv = mp.sqrt(sum(abs(v[i]) ** 2 for i in range(7)))
        if nv > mp.mpf(10) ** -30:
            cols.append(v / nv)
    assert len(cols) == d
    eta3 = mp.matrix(7, d)
    for a, c in enumerate(cols):
        for k in range(7):
            eta3[k, a] = c[k]
    D3 = Dall(3)
    sig = [dagm(eta3) * D3[i] * eta3 for i in range(len(G))]
    out[d] = {"sigma_q1": sig[i1], "sigma_q2": sig[i2], "eta": {}}
    for twoj in range(0, 19, 2):
        j = twoj // 2
        dim = int(dims[d][twoj].t[1][0]) if not dims[d][twoj].is_zero() else 0
        if dim == 0:
            continue
        Dj = Dall(j)
        basis = []
        tries = 0
        while len(basis) < dim:
            tries += 1
            X = mp.matrix(2 * j + 1, 7)
            for r in range(2 * j + 1):
                for c in range(7):
                    X[r, c] = mp.mpf(random.randint(-9, 9)) / 7 + mp.mpc(0, 1) * random.randint(-9, 9) / 11
            E = mp.matrix(2 * j + 1, 7)
            for i in range(len(G)):
                E += Dj[i] * X * dagm(D3[i])
            E = E / len(G)
            eta = E * eta3
            for b in basis:
                ip = sum((dagm(b) * eta)[a, a] for a in range(d)) / d
                eta = eta - b * ip
            nn = mp.sqrt(sum(abs(eta[r, c]) ** 2 for r in range(eta.rows) for c in range(eta.cols)) / d)
            if nn > mp.mpf(10) ** -20:
                basis.append(eta / nn)
            assert tries < 20
        # schur: eta^dag eta must be scalar; normalized so eta^dag eta = I
        r_int, r_norm, r_orth = mp.mpf(0), mp.mpf(0), mp.mpf(0)
        for b in basis:
            for i in (i1, i2):
                r_int = max(r_int, maxabs(Dj[i] * b - b * sig[i]))
            r_norm = max(r_norm, maxabs(dagm(b) * b - mp.eye(d)))
        if len(basis) == 2:
            r_orth = maxabs(dagm(basis[0]) * basis[1])
        resid[(d, twoj)] = (r_int, r_norm, r_orth)
        print(f"sector {d} level n={twoj}: dim {dim}  max|D(h)eta - eta sigma(h)| (h=q1,q2) = {mp.nstr(r_int, 3)}, "
              f"max|eta^dag eta - I| = {mp.nstr(r_norm, 3)}" + (f", max|eta^dag eta'| = {mp.nstr(r_orth, 3)}" if len(basis) == 2 else ""))
        out[d]["eta"][twoj] = basis

tol = mp.mpf(10) ** -40
allres = [max(v) for v in resid.values()]
check("all intertwiner residuals (intertwining at q1,q2; eta^dag eta = I; eta^dag eta' = 0) < 1e-40 at 50 digits",
      lambda rs: max(rs) < tol, allres, allres + [mp.mpf(10) ** -30], "append a residual 1e-30")
# the intertwining check itself must detect a wrong sigma: use sigma(q2) in place of sigma(q1)
d = 3; b = out[3]["eta"][6][0]
bad = maxabs(Dall(3)[i1] * b - b * out[3]["sigma_q2"])
check("intertwining test detects a wrong group element (sector 3, level 6)", lambda r: r < tol,
      maxabs(Dall(3)[i1] * b - b * out[3]["sigma_q1"]), bad, "compare D(q1) eta with eta sigma(q2)")

# save as complex numbers (float) and mp strings
ser = {}
for d in out:
    ser[d] = {"sigma_q1": [[str(out[d]["sigma_q1"][i, k]) for k in range(d)] for i in range(d)],
              "sigma_q2": [[str(out[d]["sigma_q2"][i, k]) for k in range(d)] for i in range(d)],
              "eta": {n: [[[complex(b[r, c]) for c in range(d)] for r in range(b.rows)] for b in basis]
                      for n, basis in out[d]["eta"].items()}}
with open("out/eta.pkl", "wb") as f:
    pickle.dump(ser, f)
with open("out/item1_residuals.json", "w") as f:
    json.dump({f"sector{d}_n{n}": [mp.nstr(x, 3) for x in v] for (d, n), v in resid.items()}, f, indent=1)
save_log("s02")
