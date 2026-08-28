"""S14 item 5 (v2): Control A reference values at s in {0.1, 0.3, 0.5}, 50 dps.
Structured assembly: in the clean Hopf convention (a = cos(eta) e^{i xi1}, b = sin(eta)
e^{i xi2}), every monomial sym-power entry factors as e^{i xi1 mu} e^{i xi2 nu} f(u) with
mu = l - r - k, nu = r - k, so four-product integrals reduce to ANALYTIC phase selection
times a 13-node Gauss-Legendre radial sum at 50 dps. The convention differs from the repo's
by an axis relabeling, an SU(2) rotation, so by the protocol's own equivariance theorem all
reference eigenquantities are identical; the record carries that argument. J splits into
(+nu, -nu)-paired real blocks; 50-dps eigenquantities come from mp.eigsy per block.
ARMS: (P) at s = 0 the assembled J equals diag(omega^2 - lambda) exactly (known limit);
(X) generic float64 node-quadrature assembly agrees with the structured route to ~1e-12;
(M) an s-perturbation moves every reference; (N) splitting and leakage nonzero (G-NULL-c);
(S) the separation check >= 0.25 * gap_A = 10 at every s.
"""
import numpy as np, mpmath as mp, json, hashlib, time, sys
from math import comb
mp.mp.dps = 50
t0 = time.time()
L = [2, 6]; OM2 = 9; C1 = 1

# ---- entries: for level l, entry (r,k): phase (mu, nu) = (l-r-k, r-k); radial poly in
# (c, s) with c^2 = 1-u, s^2 = u: sum over i, j of coeff * c^{A} s^{B}, A+B = l
def radial_poly(l, r, k):
    """returns list of (apow_c, bpow_s, coeff) for entry (r,k) of the monomial sym power,
    with a = c e^{i xi1}, b = s e^{i xi2}: term a^{l-k-i} b^i (-b bar)^{k-j} (a bar)^j,
    i + j = r."""
    terms = []
    for i in range(min(r, l - k) + 1):
        j = r - i
        if j < 0 or j > k: continue
        coeff = comb(l - k, i) * comb(k, j) * (-1) ** (k - j)
        cpow = (l - k - i) + j          # powers of c from a and a-bar
        spow = i + (k - j)              # powers of s from b and b-bar
        terms.append((cpow, spow, coeff))
    return terms

# 13-node GL on u in [0,1] at 50 dps
NU = 13
pts = mp.polyroots(mp.taylor(lambda t: mp.legendre(NU, t), 0, NU)[::-1], maxsteps=300, extraprec=120)
GLX, GLW = [], []
for rt in sorted([mp.re(r) for r in pts]):
    dP = mp.diff(lambda t: mp.legendre(NU, t), rt)
    GLX.append((rt + 1) / 2); GLW.append(1 / ((1 - rt ** 2) * dP ** 2))
CS = [(mp.sqrt(1 - u), mp.sqrt(u)) for u in GLX]

BASIS = []      # (l, r, k, mu, nu, radial values at nodes [mp], norm)
for l in L:
    for r in range(l + 1):
        for k in range(l + 1):
            vals = [sum(cf * c ** cp * s ** sp for cp, sp, cf in radial_poly(l, r, k))
                    for (c, s) in CS]
            BASIS.append([l, r, k, l - r - k, r - k, vals, None])
def rint(vals):     # integral over S3 of a pure-radial function: sum w * v (angular = 1)
    return sum(w * v for w, v in zip(GLW, vals))
for e in BASIS:     # norms: <|entry|^2> has phase 0 automatically
    e[6] = mp.sqrt(rint([v * v for v in e[5]]))
IDX = {(e[0], e[1], e[2]): i for i, e in enumerate(BASIS)}
n = len(BASIS); lam = [e[0] * (e[0] + 2) for e in BASIS]
iv2, iv6 = IDX[(2, 0, 0)], IDX[(6, 0, 0)]

def four_int(i, cci, j, ccj, p, ccp, q, ccq):
    """integral of e_i^(c) e_j^(c) e_p^(c) e_q^(c), cc = +1 plain, -1 conjugated; zero
    unless both phase sums vanish; else radial 13-node sum of the normalized entries."""
    ei, ej, ep, eq = BASIS[i], BASIS[j], BASIS[p], BASIS[q]
    mu = cci * ei[3] + ccj * ej[3] + ccp * ep[3] + ccq * eq[3]
    nu = cci * ei[4] + ccj * ej[4] + ccp * ep[4] + ccq * eq[4]
    if mu != 0 or nu != 0: return mp.mpf(0)
    vals = [a * b * c * d for a, b, c, d in zip(ei[5], ej[5], ep[5], eq[5])]
    return rint(vals) / (ei[6] * ej[6] * ep[6] * eq[6])

def assemble(s_val):
    s = mp.mpf(s_val)
    cA = {iv2: s, iv6: s * s / 2}
    Jc = mp.zeros(n, n); Ju = mp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            acc_c = mp.mpf(0); acc_u = mp.mpf(0)
            for p, cp in cA.items():
                for q, cq in cA.items():
                    acc_c += cp * cq * four_int(i, -1, j, +1, p, +1, q, -1)
                    acc_u += cp * cq * four_int(i, -1, j, -1, p, +1, q, +1)
            Jc[i, j] = 2 * C1 * acc_c
            Ju[i, j] = C1 * acc_u
    A = mp.zeros(2 * n, 2 * n)
    for i in range(n):
        for j in range(n):
            c = Jc[i, j]; u = Ju[i, j]        # both REAL here (real coefficients, real radials)
            A[i, j] += -(c + u)
            A[n + i, n + j] += -(c - u)
        A[i, i] += OM2 - lam[i]; A[n + i, n + i] += OM2 - lam[i]
    return A

def references(s_val):
    A = assemble(s_val)
    sym = max(abs(A[i, j] - A[j, i]) for i in range(0, 2 * n, 17) for j in range(0, 2 * n, 13))
    Af = np.array([[float(A[i, j]) for j in range(2 * n)] for i in range(2 * n)])
    evf, Vf = np.linalg.eigh((Af + Af.T) / 2)
    order = np.argsort(np.abs(evf)); d_c = 18
    sep = float(abs(evf[order[d_c]]) - abs(evf[order[d_c - 1]]))
    # 50-dps eigenvalues by nu-PAIRED blocks. Jc couples nu_i = nu_j, but Ju (the
    # conjugate term) couples nu_j = -nu_i, so the invariant blocks are the UNIONS
    # {i: |nu_i| = v} in each of the re and im sectors. The first implementation split by
    # SIGNED nu and dropped the Ju couplings for nu != 0; the f64 diff fields exposed it
    # (F2), and arm B below now requires those diffs at eigensolver rounding.
    mus_mp = []
    vecs = np.zeros((2 * n, d_c))
    for sector in (0, 1):
        for v in sorted({abs(e[4]) for e in BASIS}):
            ids = [i for i, e in enumerate(BASIS) if abs(e[4]) == v]
            off = sector * n
            sub = mp.matrix([[A[off + a, off + b] for b in ids] for a in ids])
            E, Q = mp.eigsy(sub)
            for c in range(len(ids)):
                mus_mp.append((E[c], [(off + a) for a in ids], [Q[r, c] for r in range(len(ids))]))
    for idx in range(len(mus_mp)):
        val, ids, q = mus_mp[idx]
        nu_of_block = abs(BASIS[ids[0] % n][4])
        mus_mp[idx] = (val, ids, q, nu_of_block)
    mus_mp.sort(key=lambda t: abs(t[0]))
    cluster = mus_mp[:d_c]
    pos = sum(m[0] for m in cluster) / d_c
    spl = max(m[0] for m in cluster) - min(m[0] for m in cluster)
    sep50 = float(abs(mus_mp[d_c][0]) - abs(mus_mp[d_c - 1][0]))
    for cidx, (val, ids, q, _nb) in enumerate(cluster):
        for a, qa in zip(ids, q): vecs[a, cidx] = float(qa)
    free = [IDX[(2, r, k)] for r in range(3) for k in range(3)]
    free_idx = free + [n + f for f in free]
    Qc, _ = np.linalg.qr(vecs)
    P = np.zeros((2 * n, 2 * n)); P[free_idx, free_idx] = 1
    svals = np.linalg.svd(P @ Qc, compute_uv=False)[:d_c]
    leak = float(np.sqrt(max(0.0, 1 - min(svals) ** 2)))
    f64_pos = float(np.mean(evf[order[:d_c]])); f64_spl = float(evf[order[:d_c]].max() - evf[order[:d_c]].min())
    return dict(s=s_val, sym=float(sym), sep=sep50, position=mp.nstr(pos, 30),
                splitting=mp.nstr(spl, 30), leakage=leak, zero_count_ref=0,
                min_abs_cluster=float(abs(mus_mp[0][0])),
                min_abs_block_nu=int(mus_mp[0][3]),
                route_per_quantity=dict(
                    position="50-dps |nu|-paired block eigenVALUES",
                    splitting="50-dps |nu|-paired block eigenVALUES",
                    min_abs_cluster="50-dps |nu|-paired block eigenVALUES; sits in the "
                        "|nu| block recorded in min_abs_block_nu, and where that is a "
                        "block the signed-nu bug never touched, its v2 bit-identity is "
                        "explained, not suspicious",
                    leakage="the SPAN of the 50-dps cluster eigenVECTORS against the free "
                        "block; a span is insensitive to rotations within the cluster, "
                        "which is why the signed-nu bug moved it only at 1e-11",
                    sep="50-dps |nu|-paired block eigenVALUES",
                    f64_pos_diff_and_spl_diff="cross-EIGENSOLVER diagnostic only"),
                f64_pos_diff=abs(float(pos) - f64_pos), f64_spl_diff=abs(float(spl) - f64_spl))

# ---- ARM P: s = 0 limit
A0 = assemble(0)
p_err = max(abs(A0[i, j] - ((OM2 - lam[i % n]) if i == j else 0))
            for i in range(0, 2 * n, 7) for j in range(0, 2 * n, 11))
print(f"  arm P (s=0 gives exactly diag(9 - lambda)): max err {mp.nstr(p_err,3)} "
      f"-> {'PASS' if p_err < mp.mpf('1e-45') else 'FAIL'}", flush=True)

# ---- ARM X: generic float64 node assembly vs structured route at s = 0.3
def generic_f64(s):
    K = 25; nuq = 13
    xs, ws = np.polynomial.legendre.leggauss(nuq)
    uu = (xs + 1) / 2; wu = ws / 2
    xi = 2 * np.pi * np.arange(K) / K
    nodes = []; W = []
    for u0, w0 in zip(uu, wu):
        c0, s0 = np.sqrt(1 - u0), np.sqrt(u0)
        for a1 in xi:
            for b1 in xi:
                nodes.append((c0 * np.exp(1j * a1), s0 * np.exp(1j * b1))); W.append(w0 / K ** 2)
    W = np.array(W)
    B = np.zeros((len(nodes), n), dtype=complex)
    for col, e in enumerate(BASIS):
        l, r, k = e[0], e[1], e[2]
        for row, (a, b) in enumerate(nodes):
            val = 0
            for cp, sp, cf in radial_poly(l, r, k):
                val += cf * abs(a) ** cp * abs(b) ** sp
            B[row, col] = val * np.exp(1j * (np.angle(a) * e[3] + np.angle(b) * e[4]))
        B[:, col] /= float(e[6])
    coef = np.zeros(n, dtype=complex); coef[iv2] = s; coef[iv6] = s * s / 2
    phi = B @ coef
    w2 = 2 * np.abs(phi) ** 2; u2 = phi ** 2
    Jc = (B.conj().T * (W * w2)) @ B
    Ju = (B.conj().T * (W * u2)) @ B.conj()
    Af = np.zeros((2 * n, 2 * n))
    for i in range(n):
        for j in range(n):
            Af[i, j] += -np.real(Jc[i, j] + Ju[i, j])
            Af[n + i, n + j] += -np.real(Jc[i, j] - Ju[i, j])
            Af[i, n + j] += -(-np.imag(Jc[i, j]) + np.imag(Ju[i, j])) * (-1)
            Af[n + i, j] += -(np.imag(Jc[i, j]) + np.imag(Ju[i, j]))
        Af[i, i] += OM2 - lam[i]; Af[n + i, n + i] += OM2 - lam[i]
    return Af
Ax = assemble(mp.mpf("0.3"))
Axf = np.array([[float(Ax[i, j]) for j in range(2 * n)] for i in range(2 * n)])
Ag = generic_f64(0.3)
x_err = float(np.abs(Axf - Ag).max())
print(f"  arm X (structured 50-dps vs generic float64 assembly at s=0.3): "
      f"max abs diff {x_err:.2e} -> {'PASS' if x_err < 1e-10 else 'FAIL'}", flush=True)

PROTO_SHA = open("/tmp/h5.txt").read().strip()
results = {"_field_semantics": {
  "position/splitting": "50-dps nu-paired block eigenVALUES under the frozen S5 definitions; the object of record",
  "f64_pos_diff/f64_spl_diff": "cross-EIGENSOLVER diagnostic: 50-dps block route vs one float64 LAPACK eigh of the full matrix; must sit at eigensolver rounding; its earlier 1e-4 values were the bug detector",
  "arm_X": "ASSEMBLY comparison: structured 50-dps matrix vs generic float64 node-quadrature matrix, entrywise; a different object from the eigen diagnostic",
  "protocol_sha": PROTO_SHA}}
for s in (0.1, 0.3, 0.5):
    results[str(s)] = references(s)
    r = results[str(s)]
    print(f"  s={s}: sep {r['sep']:.3f} (>=10), position {r['position'][:22]}, "
          f"splitting {r['splitting'][:22]}, leakage {r['leakage']:.3e}, "
          f"J-sym {r['sym']:.1e}, f64 diffs {r['f64_pos_diff']:.1e}/{r['f64_spl_diff']:.1e}", flush=True)
rp = references(0.3 + 1e-3)
dmove = abs(float(mp.mpf(rp["splitting"])) - float(mp.mpf(results["0.3"]["splitting"])))
armM = dmove > 1e-7
armN = all(float(mp.mpf(results[k]["splitting"])) > 1e-8 and results[k]["leakage"] > 1e-12
           for k in results if not k.startswith("_"))
armS = all(results[k]["sep"] >= 10 for k in results if not k.startswith("_"))
armX = x_err < 1e-10
armP = p_err < mp.mpf("1e-45")
armB = all(results[k]["f64_pos_diff"] < 1e-9 and results[k]["f64_spl_diff"] < 1e-9
           for k in results if not k.startswith("_"))
print(f"  arm B (nu-paired blocks vs float64 full-matrix eigenvalues): "
      f"{'PASS' if armB else 'FAIL'}")
print(f"  arm M (s-perturbation moves splitting by {dmove:.2e}): {'PASS' if armM else 'FAIL'}")
print(f"  arm N (splitting, leakage nonzero at every s): {'PASS' if armN else 'FAIL'}")
json.dump(results, open("raw/controlA_references.json", "w"), indent=1)
h = hashlib.sha256(open("raw/controlA_references.json", "rb").read()).hexdigest()
print(f"  controlA_references.json sha256 {h}")
print(f"  wall {time.time()-t0:.0f}s")
ok = armP and armX and armM and armN and armS and armB
print(f"  ITEM 5 VERDICT: {'GREEN' if ok else 'RED'}")
sys.exit(0 if ok else 1)
