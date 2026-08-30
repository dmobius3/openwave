"""S14-C2 item 9: the arm-on-arena demonstration record, sub-items (a) through (m).

One armed script. Field-space demonstrations run on fields drawn from the arena at rung
N = 24 (the smallest agreement rung). SELF-ARMED: section 0's arms must all fire or the
script exits red before any demonstration runs; every demonstration is a parent-plus-
mutation pair with computed verdicts and the exit code carries the result. Scores use the
protocol's ratio form E_mut >= 1e3 * max(E_parent, 1e-16); regression CONTROLS are not
arms and their green condition is machine zero. Record: raw/c2_item9_record.json,
rewritten after every section so a dying run leaves its progress. Large node sets are
STREAMED: no reference-rule basis is ever materialized whole.
"""
import numpy as np, json, time, sys, os
from math import comb

sys.path.insert(0, "/Users/blake/Desktop/MIT/OpenWave/openwave/openwave/xperiments/m8_mit/research/m8_5b/pilot")
from route_a_nonabelian import quat_to_su2, sym_power

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
REC_PATH = os.path.join(HERE, "raw", "c2_item9_record.json")
RECORD = {"item": "S14-C2 item 9", "rung": 24, "seed": 20260830, "sections": {}}
GREEN = True
SMOKE = os.environ.get("C2_SMOKE") == "1"

def save(name, d):
    global GREEN
    RECORD["sections"][name] = d
    GREEN = GREEN and bool(d.get("green", False))
    json.dump(RECORD, open(REC_PATH, "w"), indent=1, default=float)
    print(f"[{time.time()-T0:8.1f}s] {name}: {'GREEN' if d.get('green') else 'RED'}  {d.get('note','')}", flush=True)

def ratio_fire(e_mut, e_parent):
    return e_mut >= 1e3 * max(e_parent, 1e-16)

# ------------------------------------------------------------ section 0: primitives + arms
def hopf_nodes(K, nu):
    xs, ws = np.polynomial.legendre.leggauss(nu)
    u = (xs + 1) / 2; wu = ws / 2
    xi = 2 * np.pi * np.arange(K) / K
    ce, se = np.sqrt(1 - u), np.sqrt(u)
    ca, sa = np.cos(xi), np.sin(xi)
    X = np.empty((nu * K * K, 4)); W = np.empty(nu * K * K)
    idx = 0
    for i in range(nu):
        x0 = np.repeat(ce[i] * ca, K); x1 = np.repeat(ce[i] * sa, K)
        x2 = np.tile(se[i] * ca, K);  x3 = np.tile(se[i] * sa, K)
        X[idx:idx+K*K] = np.stack([x0, x1, x2, x3], axis=1)
        W[idx:idx+K*K] = wu[i]; idx += K*K
    W /= W.sum()
    return X, W

def exact_rule(D):
    return hopf_nodes(D + 1, D // 2 + 1)

def su2_arrays(X):
    # quat_to_su2's own convention: [[w+ix, y+iz], [-y+iz, w-ix]]
    a = X[:, 0] + 1j * X[:, 1]; b = X[:, 2] + 1j * X[:, 3]
    c = -X[:, 2] + 1j * X[:, 3]; d = X[:, 0] - 1j * X[:, 1]
    return a, b, c, d

def mono_table(n, C):
    """Monomial table for f_p(x) = sum_jk C[p,j,k] S^n[j,k](x), S per the pinned formula."""
    monos = {}; cols = []
    for k in range(n + 1):
        for j in range(n + 1):
            for m in range(max(0, j - (n - k)), min(k, j) + 1):
                key = ((n - k) - (j - m), j - m, k - m, m)
                w = comb(n - k, j - m) * comb(k, m)
                if key not in monos:
                    monos[key] = len(monos); cols.append(np.zeros(C.shape[0], dtype=complex))
                cols[monos[key]] += w * C[:, j, k]
    keys = np.array(sorted(monos, key=lambda kk: monos[kk]), dtype=int)
    return keys, np.stack(cols, axis=1)

def eval_tab(tab, X):
    keys, M = tab
    a, b, c, d = su2_arrays(X)
    nmax = int(keys.max()) if keys.size else 0
    pa = [np.ones(X.shape[0], dtype=complex)]; pb = [pa[0].copy()]
    pc = [pa[0].copy()]; pd = [pa[0].copy()]
    for _ in range(nmax):
        pa.append(pa[-1]*a); pb.append(pb[-1]*b); pc.append(pc[-1]*c); pd.append(pd[-1]*d)
    V = np.empty((keys.shape[0], X.shape[0]), dtype=complex)
    for i, (ea, ec, eb, ed) in enumerate(keys):
        V[i] = pa[ea] * pc[ec] * pb[eb] * pd[ed]
    return M @ V

def build_icosians():
    phi = (1 + np.sqrt(5)) / 2
    G = []
    for i in range(4):
        for s in (1, -1):
            q = np.zeros(4); q[i] = s; G.append(q)
    for signs in range(16):
        G.append(np.array([(1 if signs >> k & 1 else -1) * 0.5 for k in range(4)]))
    base = np.array([1.0, phi, 1/phi, 0.0]) / 2
    evens = [(0,1,2,3),(0,2,3,1),(0,3,1,2),(1,0,3,2),(1,2,0,3),(1,3,2,0),
             (2,0,1,3),(2,1,3,0),(2,3,0,1),(3,0,2,1),(3,1,0,2),(3,2,1,0)]
    for p in evens:
        for s0 in (1,-1):
            for s1 in (1,-1):
                for s2 in (1,-1):
                    q = np.zeros(4)
                    q[p[0]] = s0*base[0]; q[p[1]] = s1*base[1]; q[p[2]] = s2*base[2]
                    G.append(q)
    return np.unique(np.round(np.array(G), 12), axis=0)

def qmul(p, q):
    w1,x1,y1,z1 = p; w2,x2,y2,z2 = q
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])

def group_arms(G):
    if G.shape[0] != 120: return False, {}
    S = {tuple(np.round(q, 9)) for q in G}
    closed = all(tuple(np.round(qmul(G[i], G[j]), 9)) in S
                 for i in range(0, 120, 7) for j in range(0, 120, 11))
    orders = {}
    for q in G:
        p = q.copy(); k = 1
        while not np.allclose(p, [1,0,0,0], atol=1e-9):
            p = qmul(p, q); k += 1
            if k > 60: break
        orders[k] = orders.get(k, 0) + 1
    return closed and orders == {1:1,2:1,3:20,4:30,5:24,6:20,10:24}, orders

def reynolds_basis(G, n):
    dim = n + 1
    P = np.zeros((dim*dim, dim*dim), dtype=complex)
    for q in G:
        qi = np.array([q[0], -q[1], -q[2], -q[3]])
        A = sym_power(quat_to_su2(qi), n).T
        P += np.kron(A, np.eye(dim))   # row-major vec: C' = A C is kron(A, I); kron(I, A) built the RIGHT quotient (the l/m failure)
    P /= len(G)
    w, V = np.linalg.eig(P)
    B = V[:, np.abs(w - 1) < 1e-8]
    B, _ = np.linalg.qr(B)
    return B.T.reshape(-1, dim, dim)

print("== section 0: primitives and self-arms ==", flush=True)
rng0 = np.random.default_rng(11)
werr = 0.0
for n in (3, 12, 24):
    for _ in range(15):
        q = rng0.standard_normal(4); q /= np.linalg.norm(q)
        C1 = np.zeros((1, n+1, n+1), dtype=complex)
        jj, kk = rng0.integers(0, n+1, 2); C1[0, jj, kk] = 1.0
        val = eval_tab(mono_table(n, C1), q[None, :])
        werr = max(werr, abs(val[0,0] - sym_power(quat_to_su2(q), n)[jj, kk]))
G120 = build_icosians()
armC, census = group_arms(G120)
armD = reynolds_basis(G120, 6).shape[0] == 0
save("0-selfarms", {"green": bool(werr < 1e-11 and armC and armD),
    "batched_vs_sym_power": werr, "group_census": {str(k): v for k, v in census.items()},
    "reynolds_negative_control_n6_dim": 0 if armD else -1,
    "note": "pinned-formula arm; 120-closure and order-census arm; n=6 empty-invariant arm"})
if not GREEN: sys.exit(1)

# ------------------------------------------------------------ bases at N = 24
print("== bases at N=24 ==", flush=True)
N = 24
LEVELS = [0, 12, 20, 24]
CINV = {0: np.ones((1,1,1), dtype=complex)}
for n in LEVELS[1:]:
    CINV[n] = reynolds_basis(G120, n)
dims = {n: CINV[n].shape[0] for n in LEVELS}
dims_ok = dims == {0:1, 12:13, 20:21, 24:25}
TAB = {n: mono_table(n, CINV[n]) for n in LEVELS}

def arena_raw(X):
    return np.vstack([eval_tab(TAB[n], X) for n in LEVELS])

X4, W4 = exact_rule(4*N)
Braw = arena_raw(X4)
Gm = (Braw * W4) @ Braw.conj().T
ev, Uv = np.linalg.eigh(Gm)
LwT = (Uv @ np.diag(ev**-0.5) @ Uv.conj().T).conj().T     # orthonormalizer (left-multiply)
Bo = LwT @ Braw
# SYMMETRIC refinement of the SAME Lowdin map (per the redline ruling: S 3.3 freezes
# symmetric orthonormalization, canonical and ordering-independent, so no triangular
# factor may become the effective basis map). The raw basis is ill-conditioned (the S1b
# non-unitarity fact) and a single eigh-based inverse square root leaves ~1e-10 residual;
# a second symmetric inverse-square-root factor of the RESIDUAL Gram converges to the
# same canonical G^{-1/2} (exactly it in exact arithmetic, since the residual Gram is I),
# evaluated more accurately. No ordering enters at any step.
G2 = (Bo * W4) @ Bo.conj().T
ev2, Uv2 = np.linalg.eigh(G2)
S2 = Uv2 @ np.diag(ev2**-0.5) @ Uv2.conj().T
Bo = S2 @ Bo
LwT = S2 @ LwT
gram_err = np.abs((Bo * W4) @ Bo.conj().T - np.eye(60)).max()
# direct left-invariance arm: for a group element h, the level-12 invariant coefficient
# matrices must satisfy D(h^{-1})^T C = C exactly; the kron-convention bug produced
# right-invariant matrices and this check is the one that fires on that mistake directly.
qh = G120[37]
qhi = np.array([qh[0], -qh[1], -qh[2], -qh[3]])
Dh = sym_power(quat_to_su2(qhi), 12).T
linv_err = max(np.abs(Dh @ CINV[12][i] - CINV[12][i]).max() for i in range(13))
rng = np.random.default_rng(20260830)
coef = rng.standard_normal(60) + 1j*rng.standard_normal(60); coef /= np.linalg.norm(coef)
psi = coef @ Bo
rt_err = np.abs(((Bo * W4) @ psi.conj()).conj() - coef).max()
save("B-bases", {"green": bool(dims_ok and gram_err < 1e-11 and rt_err < 1e-11 and linv_err < 1e-11),
    "dims": {str(k): v for k, v in dims.items()}, "gram_err": gram_err, "roundtrip_err": rt_err,
    "left_invariance_err": linv_err,
    "note": "Reynolds dims {1,13,21,25}; symmetric Lowdin (refined, still the canonical map); Gram, roundtrip, and DIRECT left-invariance at rounding"})
if not GREEN or SMOKE:
    sys.exit(0 if GREEN else 1)

def coeffs_of(Bmat, Wt, vals):
    return ((Bmat * Wt) @ vals.conj()).conj()

def cubic_vals(vals):
    return (np.abs(vals)**2) * vals

ref4 = coeffs_of(Bo, W4, cubic_vals(psi))
Eref = np.abs(ref4).max()

def stream_reads(D_or_rule, read_tabs, orth_read=None, chunk=200000):
    """Stream a rule's nodes: synthesize psi from the arena basis, form the cubic, and
    accumulate projections onto each read table (raw) and optionally the orthonormal
    arena basis. Nothing large is materialized. The chunk is CAPPED by the largest read
    table's monomial count so the per-chunk value array stays near 1.5 GB: at the
    high-odd read levels the default chunk allocated ~16 GB per level and run 3 thrashed
    on a 24 GB machine. Pure batching; no computed value depends on the chunk."""
    X, W = D_or_rule
    maxmono = max([tab[0].shape[0] for tab in read_tabs] + [4000])
    chunk = min(chunk, max(20000, int(1.2e8 / maxmono)))
    accs = [None]*len(read_tabs); acc_o = None
    for s in range(0, X.shape[0], chunk):
        Xc = X[s:s+chunk]; Wc = W[s:s+chunk]
        Bc = LwT @ arena_raw(Xc)
        gc = cubic_vals(coef @ Bc)
        for i, tab in enumerate(read_tabs):
            Rc = eval_tab(tab, Xc)
            p = ((Rc * Wc) @ gc.conj()).conj()
            accs[i] = p if accs[i] is None else accs[i] + p
        if orth_read:
            p = ((Bc * Wc) @ gc.conj()).conj()
            acc_o = p if acc_o is None else acc_o + p
    return accs, acc_o

# ------------------------------------------------------------ (a)(b)(c): gate 4
print("== (a)(b)(c) gate-4 arms ==", flush=True)
_, ref8 = stream_reads(exact_rule(8*N), [], orth_read=True)
E_parent = np.abs(ref4 - ref8).max() / np.abs(ref8).max()

Xe, We = hopf_nodes(2*N, 2*N + 1)
Be = LwT @ arena_raw(Xe)
r_evenK = coeffs_of(Be, We, cubic_vals(coef @ Be))
E_mut_a = np.abs(r_evenK - ref8).max() / np.abs(ref8).max()
save("a-evenK-drop", {"green": bool(E_parent < 1e-11 and ratio_fire(E_mut_a, E_parent)),
    "E_parent": E_parent, "E_mut": E_mut_a, "drop_nodes": int(Xe.shape[0]),
    "note": "even-K angular drop (K=2N, u exact) on the arena-drawn field; parent = 4N-vs-8N agreement"})

LEAK = [1, 3, 5]
def full_level_tab(n):
    dim = n + 1
    C = np.zeros((dim*dim, dim, dim), dtype=complex)
    for i in range(dim):
        for j in range(dim):
            C[i*dim + j, i, j] = 1.0
    return mono_table(n, C)
LTAB = {n: full_level_tab(n) for n in LEAK}
Lk4 = np.vstack([eval_tab(LTAB[n], X4) for n in LEAK])
leak_parent = np.abs(coeffs_of(Lk4, W4, cubic_vals(psi))).max() / Eref
Xo, Wo = exact_rule(2*N)
BoK = LwT @ arena_raw(Xo)
psiO = coef @ BoK
LkO = np.vstack([eval_tab(LTAB[n], Xo) for n in LEAK])
offs = np.cumsum([0] + [(n+1)**2 for n in LEAK])
lm_all = coeffs_of(LkO, Wo, cubic_vals(psiO))
leak_per_level = {str(n): float(np.abs(lm_all[offs[i]:offs[i+1]]).max() / Eref) for i, n in enumerate(LEAK)}
leak_mut = max(leak_per_level.values())
save("b-crossparity-leakage", {"green": bool(leak_parent < 1e-11 and ratio_fire(leak_mut, leak_parent)),
    "E_parent": leak_parent, "E_mut": leak_mut, "read_levels": LEAK, "per_level": leak_per_level,
    "note": "arm-local scalar read at n in {1,3,5}, truth zero by parity; odd-K rule fires, exact rule at rounding"})

r_odd = coeffs_of(BoK, Wo, cubic_vals(psiO))
E_reg = np.abs(r_odd - ref8).max() / np.abs(ref8).max()
save("c-oddK-regression-CONTROL", {"green": bool(E_reg < 1e-11), "E": E_reg,
    "note": "superseded odd-K drop DEAD on arena fields; green condition IS machine zero (control, not an arm)"})

# ------------------------------------------------------------ (d): gate 1
lam = np.concatenate([[0.0], np.full(13, 12*14.0), np.full(21, 20*22.0), np.full(25, 24*26.0)])
A0 = np.diag(lam.astype(complex))
K0 = np.linalg.norm(A0 - A0.conj().T, 2); J0 = float(np.abs(np.linalg.eigvals(A0).imag).max())
gam = 1e-3 * np.linalg.norm(A0, 2)
Ai = A0.copy(); Ai[3, 7] += gam; Ai[7, 3] -= gam   # SAME level (both lambda=168): distinct-level coupling keeps eigenvalues real and J cannot fire
Ki = np.linalg.norm(Ai - Ai.conj().T, 2); Ji = float(np.abs(np.linalg.eigvals(Ai).imag).max())
thr = 100 * np.finfo(float).eps * np.linalg.norm(A0, 2)
save("d-gate1-injection", {"green": bool(K0 <= thr and J0 <= thr and Ki > 1e3*thr and Ji > 1e3*thr),
    "K_parent": K0, "J_parent": J0, "K_inj": Ki, "J_inj": Ji, "threshold": thr,
    "note": "built N=24 R0 diagonal; anti-Hermitian injection between two retained modes"})

# ------------------------------------------------------------ (e): gate 2 mislabel
def fd_casimir(idx, nsample=140, eps=1e-4):
    rs = np.random.default_rng(5)
    pick = rs.choice(X4.shape[0], nsample, replace=False)
    num = den = 0.0
    for p in pick:
        x = X4[p]
        Tfr = np.linalg.svd(x[None, :])[2][1:]
        lap = 0.0
        for t in Tfr:
            xp = (np.cos(eps)*x + np.sin(eps)*t)[None, :]
            xm = (np.cos(eps)*x - np.sin(eps)*t)[None, :]
            fp = (LwT @ arena_raw(xp))[idx, 0]
            fm = (LwT @ arena_raw(xm))[idx, 0]
            lap += (fp + fm - 2*Bo[idx, p]) / eps**2
        num += (np.conj(Bo[idx, p]) * lap).real
        den += abs(Bo[idx, p])**2
    return -num / den

meas = fd_casimir(5)                       # a built level-12 element
true_ev, fake_ev = 12*14.0, 20*22.0
save("e-gate2-mislabel", {"green": bool(abs(meas-true_ev)/true_ev < 1e-3 and abs(meas-fake_ev)/fake_ev > 0.1),
    "measured": meas, "true": true_ev, "mislabel_claim": fake_ev,
    "note": "FD Casimir through the sampling map on a BUILT element; relabel 12->20 fires"})

# ------------------------------------------------------------ (k)(l)(m): gate 5
print("== (k)(l)(m) gate-5 arms ==", flush=True)
absq = np.abs(psi)**2; sq = psi**2
J = np.zeros((120, 120))
for p in range(60):
    for reim in (0, 1):
        dvec = Bo[p] * (1.0 if reim == 0 else 1j)
        dc = coeffs_of(Bo, W4, 2*absq*dvec + sq*np.conj(dvec))
        J[:60, p+60*reim] = dc.real; J[60:, p+60*reim] = dc.imag
symrel = np.linalg.norm(J - J.T, 2) / np.linalg.norm(J, 2)
L = np.zeros((120, 120)); L[3, 47] = 0.05 * np.linalg.norm(J, 2)
Jp = J + L
symrel_mut = np.linalg.norm(Jp - Jp.T, 2) / np.linalg.norm(Jp, 2)
save("k-nongradient", {"green": bool(symrel <= 1e-12 and ratio_fire(symrel_mut, symrel)),
    "parent_rel_asym": symrel, "mut_rel_asym": symrel_mut,
    "note": "cubic Jacobian on the arena-drawn E_R0 state; parent within the frozen 1e-12 bound; a non-gradient term reds exactly that predicate"})

qg = np.array([0.3, 0.5, -0.2, 0.78]); qg /= np.linalg.norm(qg)
Xg = np.stack([qmul(x, qg) for x in X4])
Bg = LwT @ arena_raw(Xg)                    # values of v_q at X4·g == values of R_g v_q at X4
Agm = (Bo * W4) @ Bg.conj().T               # conj(<v_p, R_g v_q>)
unit_err = np.abs(Agm.conj().T @ Agm - np.eye(60)).max()
coef2 = rng.standard_normal((2, 60)) + 1j*rng.standard_normal((2, 60))
coef2 /= np.linalg.norm(coef2)
psi2 = coef2 @ Bo
dens = (np.abs(psi2)**2).sum(axis=0)
Ncub = dens * psi2
psi2g = coef2 @ Bg                          # values of R_g psi at X4
densg = (np.abs(psi2g)**2).sum(axis=0)
Ng = densg * psi2g
NgA = np.stack([coeffs_of(Bo, W4, Ng[i]) for i in range(2)])
NA  = np.stack([coeffs_of(Bo, W4, Ncub[i]) for i in range(2)])
NA_g = NA @ Agm.conj().T   # coeff of v_p in R_g v_q is conj(Agm[p,q]); the missing .T was masked in run 1 by the near-diagonal wrong-quotient Agm
equiv_rel = np.abs(NgA - NA_g).max() / np.abs(NA).max()
w1 = 1 + 0.3 * np.real(eval_tab(LTAB[1], X4)[0])
w1g = 1 + 0.3 * np.real(eval_tab(LTAB[1], Xg)[0])
Nb  = np.stack([coeffs_of(Bo, W4, w1 * Ncub[i]) for i in range(2)])
Nbg = np.stack([coeffs_of(Bo, W4, w1g * Ng[i]) for i in range(2)])
equiv_mut = np.abs(Nbg - Nb @ Agm.conj().T).max() / np.abs(Nb).max()
save("l-symbreak-coupling", {"green": bool(equiv_rel < 1e-10 and ratio_fire(equiv_mut, equiv_rel) and unit_err < 1e-10),
    "parent": equiv_rel, "mut": equiv_mut, "represented_action_unitarity_err": unit_err,
    "note": "E_R0 (x) C^2 arena-drawn state; cubic equivariance under a right translation; a non-invariant weight reds it"})

Dl = np.diag(lam)
comm = np.linalg.norm(Dl @ Agm - Agm @ Dl, 2) / np.linalg.norm(Dl, 2)
lam_b = lam.copy(); lam_b[1:14] += np.linspace(0, 30, 13)
Db = np.diag(lam_b)
comm_mut = np.linalg.norm(Db @ Agm - Agm @ Db, 2) / np.linalg.norm(Db, 2)
save("m-intralevel-break", {"green": bool(comm < 1e-10 and ratio_fire(comm_mut, comm)),
    "parent": comm, "mut": comm_mut,
    "note": "level-diagonal commutes with the measured right action; intra-level spread reds it"})

del Bg, psi2, psi2g, Ncub, Ng

# ------------------------------------------------------------ (f): gate 6 monitor arms
print("== (f) gate-6 monitor arms (streamed) ==", flush=True)
C30 = reynolds_basis(G120, 30)
T30 = mono_table(30, C30)
(bandacc,), _ = stream_reads(exact_rule(6*N), [T30])
(band10, leak10), _ = stream_reads(exact_rule(10*N), [T30, LTAB[1]])
E_par_f = np.abs(bandacc - band10).max() / np.abs(band10).max()
(band_drop,), _ = stream_reads(hopf_nodes(4*N, 3*N + 1), [T30])
E_mut_f = np.abs(band_drop - band10).max() / np.abs(band10).max()
Eband = np.abs(band10).max()
# gate-6 leakage reads HIGH odd levels {25,27,29}: the read must push the integrand past
# the substitution rule's exactness (72 + n > 96 needs n >= 25); at {1,3,5} the 4N rule is
# EXACT on the read integrand and the arm cannot fire, the defect this run caught in the
# draft's own S 4.2(b) sentence.
LEAK6 = [N + 1, N + 3, N + 5]          # RUNG-RELATIVE per the amended S 4.2(b)
LT6 = {n: full_level_tab(n) for n in LEAK6}
l6 , _ = stream_reads(exact_rule(6*N), [LT6[n] for n in LEAK6])
leak_par_f = max(np.abs(v).max() for v in l6) / Eband
s6, _ = stream_reads((X4, W4), [LT6[n] for n in LEAK6])
leak6_per_level = {str(n): float(np.abs(s6[i]).max() / Eband) for i, n in enumerate(LEAK6)}
leak_mut_f = max(leak6_per_level.values())
(band_sub,), _ = stream_reads((X4, W4), [T30])
E_reg_f = np.abs(band_sub - band10).max() / np.abs(band10).max()
v30 = eval_tab(T30, X4)[0]
n30 = np.sqrt(np.real((np.abs(v30)**2 * W4).sum()))
cinj = 2 * 0.1 * np.linalg.norm(band10) / n30       # inject band content at 2x the 0.1 threshold scale
inj_read = coeffs_of(eval_tab(T30, X4), W4, cubic_vals(psi) + cinj * v30)
inj_gain = np.abs(inj_read - band_sub).max()
fire_inj = inj_gain > 0.5 * cinj * n30**2
save("f-monitor-arms", {"green": bool(E_par_f < 1e-10 and ratio_fire(E_mut_f, E_par_f)
        and leak_par_f < 1e-10 and ratio_fire(leak_mut_f, leak_par_f)
        and E_reg_f < 1e-10 and fire_inj),
    "band_level": 30, "E_parent_drop": E_par_f, "E_mut_drop": E_mut_f,
    "leak_parent": leak_par_f, "leak_mut": leak_mut_f, "leak_per_level": leak6_per_level, "read_levels": LEAK6,
    "regression_CONTROL_substitution": E_reg_f, "injected_gain": inj_gain,
    "note": "band content at n=30: even-K monitor drop fires; leakage under the 4N-for-6N substitution fires; the substitution itself is DEAD on the band read (control); injected high-band content at 2x threshold is read"})

# ------------------------------------------------------------ (g): Control A
print("== (g) Control A ==", flush=True)
XA, WA = exact_rule(24)
BA = np.vstack([eval_tab(full_level_tab(n), XA) for n in (2, 6)])
GA = (BA * WA) @ BA.conj().T
evA, UA = np.linalg.eigh(GA)
BAo = (UA @ np.diag(evA**-0.5) @ UA.conj().T).conj().T @ BA
lamA = np.concatenate([np.full(9, 8.0), np.full(49, 48.0)])
om2 = 9.0
def RA(phi):
    vals = phi @ BAo
    return (lamA - om2) * phi + coeffs_of(BAo, WA, cubic_vals(vals))
v2 = np.zeros(58, dtype=complex); v2[0] = 1
v6 = np.zeros(58, dtype=complex); v6[9] = 1
v22 = np.zeros(58, dtype=complex); v22[1] = 1
s = 0.3
phistar = s*v2 + (s**2/2)*v6
Gforce = -RA(phistar)
def F(phi, breakop=False):
    r = RA(phi) + Gforce
    return r + 0.05*np.roll(phi, 3) if breakop else r
def newton(phi0, breakop=False):
    phi = phi0.copy()
    for _ in range(30):
        r = F(phi, breakop)
        if np.linalg.norm(r) < 8e-13:
            break
        Jr = np.zeros((116, 116))
        for p in range(58):
            for reim in (0, 1):
                d = np.zeros(58, dtype=complex); d[p] = 1.0 if reim == 0 else 1j
                h = 1e-7
                dr = (F(phi + h*d, breakop) - F(phi - h*d, breakop)) / (2*h)
                Jr[:58, p+58*reim] = dr.real; Jr[58:, p+58*reim] = dr.imag
        st = np.linalg.solve(Jr, np.concatenate([r.real, r.imag]))
        phi = phi - (st[:58] + 1j*st[58:])
    return phi
ret = np.linalg.norm(newton(phistar + 1e-3*np.linalg.norm(phistar)*v22) - phistar) / np.linalg.norm(phistar)
retb = np.linalg.norm(newton(phistar + 1e-3*np.linalg.norm(phistar)*v22, True) - phistar) / np.linalg.norm(phistar)
save("g-controlA-return", {"green": bool(ret < 1e-10 and ratio_fire(retb, max(ret, 1e-16))),
    "return_parent": ret, "return_mut": retb, "s": s,
    "note": "perturbed seed returns to the exact path; a mutated operator does not"})

# ------------------------------------------------------------ (h): kick-drift
def integrate(dt, T, kind):
    q, p = 1.0, 0.0
    n = int(round(T/dt)); out = np.empty(n)
    for i in range(n):
        if kind == "kick":
            p -= dt * q**3; q += dt * p
        else:
            p -= 0.5*dt * q**3; q += dt * p; p -= 0.5*dt * q**3
        out[i] = q
    return out
Th = 6.0; base = 2e-5
refq = integrate(base, Th, "leap")
def err(dt, kind):
    tr = integrate(dt, Th, kind)
    step = int(round(dt/base))
    return np.abs(tr - refq[step-1::step][:len(tr)]).max()
dts = [0.02, 0.01, 0.005]
e_leap = [err(d, "leap") for d in dts]
e_kick = [err(d, "kick") for d in dts]
c_leap = min(e_leap[i]/e_leap[i+1] for i in range(2))
c_kick = min(e_kick[i]/e_kick[i+1] for i in range(2))
save("h-kickdrift", {"green": bool(c_leap >= 3.0 and c_kick < 3.0),
    "contraction_parent_leapfrog": c_leap, "contraction_kickdrift": c_kick,
    "errors_leapfrog": e_leap, "errors_kickdrift": e_kick,
    "note": "control (ii): parent passes the 3x per-rung gate, first-order kick-drift fails it"})

# ------------------------------------------------------------ (i): section-9 tail
floor = 2.22e-14
def conv(seq):
    okc = True
    for r in range(len(seq)-1):
        if seq[r] > floor:
            okc &= (seq[r+1] <= seq[r]/3)
    return okc
save("i-noncontracting-tail", {"green": bool(conv([1e-3,3e-4,9e-5,2.6e-5]) and not conv([1e-3,3e-4,2.5e-4,2.4e-4])),
    "note": "the frozen e/3 rule: contracting parent GREEN, injected non-contracting tail RED"})

# ------------------------------------------------------------ (j): stub registry
REG = {"A-R0-N24": True, "A-R0C2-N24": True, "A-CTRLA": True, "A-SECTOR-R3-N24": False}
log = []
def nl_eval(arena):
    if not REG[arena]:
        log.append(("REFUSED", arena)); return "REFUSED"
    log.append(("RAN", arena)); return "RAN"
r1 = nl_eval("A-SECTOR-R3-N24")
for a in ("A-R0-N24", "A-R0C2-N24"):      # A-CTRLA deliberately skipped
    nl_eval(a)
ran = {a for v, a in log if v == "RAN"}
cov_trips = not {a for a, okf in REG.items() if okf}.issubset(ran)
save("j-stub-registry", {"green": bool(r1 == "REFUSED" and cov_trips),
    "log": log, "note": "out-of-registry nonlinear call REFUSED and logged; skipped arena trips the two-sided coverage check"})

print(f"\nITEM 9 VERDICT: {'GREEN, all thirteen demonstrations' if GREEN else 'RED'}  wall {time.time()-T0:.1f}s")
RECORD["verdict"] = "GREEN" if GREEN else "RED"
RECORD["wall_seconds"] = time.time() - T0
json.dump(RECORD, open(REC_PATH, "w"), indent=1, default=float)
sys.exit(0 if GREEN else 1)
