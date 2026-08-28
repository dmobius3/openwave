"""Arm the IMPLEMENTATION bridges of the M8.5-C symmetry derivation (see
M8_5C_SYMMETRY_DERIVATION.md). The mathematics is DERIVED there; each check arms one
theorem-to-code bridge, separately, with its own green parent and its own mutation. No omnibus
PASS: a red names the broken bridge.

SCOPE (two layers, do not conflate): a green here means the REPO SCALAR PRIMITIVES
(sym_power/quat_to_su2 harmonics, full-S^3 toy projector) realize the derived right action.
The future W_rho-valued sector bases and the real Galerkin system are qualified LATER, by
gates 3 and 5 of the protocol; this script gives them no pre-emptive credit.

  C1  right-action realization: pi_l is a homomorphism, l = 1..7 and 12
  C2  projector commutation: assembled R_num matches the analytic coefficient rep
  C3  cubic equivariance in projected coefficients (mutation under the 6N rule)
  C4  complete-level necessity: parent = complete level, mutation = dropped column
  C5  multiplicity: Molien half two routes; first occurrence DISCOVERED, not transcribed
  C6  2I-commutant census: dim = <chi_l, chi_l>, commutative iff multiplicity-free
  C7  level-diagonal spectrum vs the assembled right action (the assembled object is the
      action; the spectrum is analytic n(n+2), and C2 already implies this parent)

Preflights P0 (group forensics) and P1 (character table) are ARMED checks like the rest,
routed through report() so the output ledger records structure, never a traceback.

Dependency: route_a_nonabelian (repo, m8_5b/pilot) via PYTHONPATH.
"""
import numpy as np, sys, time, hashlib, itertools
import route_a_nonabelian as RAN
from route_a_nonabelian import quat_to_su2, sym_power
import scipy
from scipy.linalg import block_diag

# ---------- provenance: a green must name the implementation that produced it ----------
_src = open(RAN.__file__, "rb").read()
print(f"  provenance: route_a_nonabelian = {RAN.__file__}")
print(f"              sha256 {hashlib.sha256(_src).hexdigest()}")
print(f"              numpy {np.__version__}, scipy {scipy.__version__}")

rng = np.random.default_rng(20260827)
t0 = time.time()
results = []
def report(tag, green, red_mut, note=""):
    results.append(green and red_mut)
    print(f"  {tag}: parent {'PASS' if green else 'FAIL'}, "
          f"mutation {'RED (good)' if red_mut else 'NOT RED (broken arm)'}  {note}")

# ---------- the group: all 120 unit icosians, explicit, fully checked ----------
TAU = (1 + np.sqrt(5)) / 2
def icosians():
    Q = []
    for i in range(4):
        for s in (1.0, -1.0):
            q = np.zeros(4); q[i] = s; Q.append(q)
    for signs in range(16):
        Q.append(np.array([0.5 if (signs >> k) & 1 == 0 else -0.5 for k in range(4)]))
    base = np.array([0.0, 1.0, 1/TAU, TAU]) / 2
    evens = [p for p in itertools.permutations(range(4))
             if sum(1 for a in range(4) for b in range(a+1,4) if p[a] > p[b]) % 2 == 0]
    for p in evens:
        for signs in range(8):
            sb = np.array([1.0, (-1)**(signs & 1), (-1)**((signs >> 1) & 1),
                           (-1)**((signs >> 2) & 1)])
            v = base * sb
            Q.append(np.array([v[p.index(k)] for k in range(4)]))
    return np.array(Q)

def qmul(a, b):
    w1,x1,y1,z1 = a; w2,x2,y2,z2 = b
    return np.array([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2])

G120 = icosians()
def order_of(q):
    p = q.copy()
    for k in range(1, 12):
        if np.abs(p - np.array([1,0,0,0])).sum() < 1e-9: return k
        p = qmul(p, q)
    return -1
def closure_residual(G):
    prods = np.einsum("iab,jb->ija",
                      np.array([[[ q[0],-q[1],-q[2],-q[3]],
                                 [ q[1], q[0],-q[3], q[2]],
                                 [ q[2], q[3], q[0],-q[1]],
                                 [ q[3],-q[2], q[1], q[0]]] for q in G]), G)
    return np.abs(prods[:, :, None, :] - G[None, None, :, :]).sum(axis=3).min(axis=2).max()
def P0():
    census = {}
    for q in G120: census[order_of(q)] = census.get(order_of(q), 0) + 1
    d = closure_residual(G120)
    green = (G120.shape == (120, 4)
             and len({tuple(np.round(q, 9)) for q in G120}) == 120
             and np.allclose(np.einsum("ij,ij->i", G120, G120), 1.0)
             and d < 1e-9
             and census == {1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24})
    G_bad = G120.copy()
    v = rng.standard_normal(4); G_bad[7] = v / np.linalg.norm(v)   # a non-icosian unit quaternion
    d_bad = closure_residual(G_bad)
    report("P0 group forensics", green, d_bad > 1e-3,
           f"(closure {d:.1e}, census {census == {1:1,2:1,3:20,4:30,5:24,6:20,10:24}}; "
           f"non-icosian substitution mutation {d_bad:.1e})")
    return green
P0()

def pi(n, q):
    return np.array([[1.0+0j]]) if n == 0 else sym_power(quat_to_su2(np.asarray(q, float)), n)
def rand_su2_quat():
    q = rng.standard_normal(4); return q / np.linalg.norm(q)

# ---------- C1: right-action realization, l = 1..7 and 12 ----------
def C1():
    Ls = [1,2,3,4,5,6,7,12]
    xs = [rand_su2_quat() for _ in range(3)]
    gs = [rand_su2_quat(), G120[17], G120[93]]
    err = max(np.abs(pi(l, qmul(x, g)) - pi(l, x) @ pi(l, g)).max()
              for l in Ls for x in xs for g in gs)
    errL = max(np.abs(pi(l, qmul(g, x)) - pi(l, g) @ pi(l, x)).max()
               for l in Ls for x in xs for g in gs)
    green = err < 1e-9 and errL < 1e-9
    err_m = max(np.abs(pi(12, qmul(x, g)) - pi(12, x) @ pi(12, g).T).max()
                for x in xs for g in gs[:2])
    report("C1 right-action realization", green, err_m > 1e-3,
           f"(l in {{1..7,12}}: homom {err:.1e}/{errL:.1e}; transpose mutation {err_m:.1e})")
C1()

# ---------- quadrature machinery ----------
def hopf_rule(D):
    K = D + 1; nu = D // 2 + 1
    xs, ws = np.polynomial.legendre.leggauss(nu)
    u = (xs + 1) / 2; wu = ws / 2
    xi = 2*np.pi*np.arange(K)/K
    ce, se = np.sqrt(1-u), np.sqrt(u)
    X, W = [], []
    for cu, su, w in zip(ce, se, wu):
        for a in xi:
            for b in xi:
                X.append([cu*np.cos(a), cu*np.sin(a), su*np.cos(b), su*np.sin(b)]); W.append(w)
    X = np.array(X); W = np.array(W); W /= W.sum()
    return X, W

N = 3
def modes(X, M, drop_cols=None):
    cols = []
    for x in X:
        v = [np.array([1.0+0j])]
        for n in range(1, M+1):
            v.append(pi(n, x).reshape(-1))
        cols.append(np.concatenate(v))
    Y = np.array(cols)
    return np.delete(Y, drop_cols, axis=1) if drop_cols is not None else Y

X4, W4 = hopf_rule(4*N); X6, W6 = hopf_rule(6*N)
Y4 = modes(X4, N);       Y6 = modes(X6, N)
G4 = np.linalg.pinv(Y4.conj().T @ (W4[:, None] * Y4))
G6 = np.linalg.pinv(Y6.conj().T @ (W6[:, None] * Y6))
nmodes = Y4.shape[1]
lam = np.concatenate([[n*(n+2)]*(n+1)**2 for n in range(N+1)]).astype(float)

def right_translate_nodes(X, g):
    return np.array([qmul(x, g) for x in X])
def coeffs4(vals): return G4 @ (Y4.conj().T @ (W4 * vals))
def coeffs6(vals): return G6 @ (Y6.conj().T @ (W6 * vals))

def right_rep_on_coeffs(g, swap_kron=False):
    """Field coefficients under R_g: basis entries B_ij(xg) = sum_m B_im(x) P_mj regroup the
    field as C' = C P^T, whose row-major vec is (I kron P) vec(C). swap_kron is C2's mutation."""
    blocks = [np.array([[1.0+0j]])]
    for n in range(1, N+1):
        P = pi(n, g)
        blocks.append(np.kron(P, np.eye(n+1)) if swap_kron else np.kron(np.eye(n+1), P))
    return block_diag(*blocks)

g0 = rand_su2_quat()
Y4g = modes(right_translate_nodes(X4, g0), N)
Y6g = modes(right_translate_nodes(X6, g0), N)
R_num = G4 @ (Y4.conj().T @ (W4[:, None] * Y4g))     # assembled action, columns = basis images

# ---------- C2: assembled vs analytic coefficient rep ----------
def C2():
    R_an = right_rep_on_coeffs(g0)
    err = np.abs(R_num - R_an).max() / np.abs(R_an).max()
    c = rng.standard_normal(nmodes) + 1j*rng.standard_normal(nmodes)
    lhs = coeffs4(Y4g @ c); rhs = R_an @ c
    err2 = np.abs(lhs - rhs).max() / np.abs(rhs).max()
    green = err < 1e-10 and err2 < 1e-10
    R_mut = right_rep_on_coeffs(g0, swap_kron=True)
    err_m = np.abs(R_num - R_mut).max() / np.abs(R_num).max()
    report("C2 projector commutation", green, err_m > 1e-3,
           f"(assembled-vs-analytic {err:.1e}, field test {err2:.1e}; kron-swap mutation {err_m:.1e})")
C2()

# ---------- C3: cubic equivariance; mutation with CORRECT sides under the 6N rule ----------
def C3():
    cube = lambda v: (np.abs(v)**2) * v
    c = rng.standard_normal(nmodes) + 1j*rng.standard_normal(nmodes)
    vals, vals_g = Y4 @ c, Y4g @ c
    lhs = coeffs4(cube(vals_g))                       # N(R_g psi), degree 3N: exact at 4N
    rhs = right_rep_on_coeffs(g0) @ coeffs4(cube(vals))
    err = np.abs(lhs - rhs).max() / np.abs(rhs).max()
    green = err < 1e-9
    # mutation: N_w(psi) = w(x)|psi|^2 psi with w(x) = 1 + 3 x_0 breaks right-equivariance.
    #   N_w(R_g psi)(x) = w(x) |psi(xg)|^2 psi(xg)      -> w AT THE ORIGINAL POINT
    #   R_g N_w(psi)(x) = w(xg)|psi(xg)|^2 psi(xg)
    # Integrand degree <= 4N + 1 <= 6N: run BOTH sides under the 6N rule so any red is the
    # symmetry break, never aliasing.
    w6  = 1.0 + 3.0 * X6[:, 0]
    v6, v6g = Y6 @ c, Y6g @ c
    lhs_m = coeffs6(w6 * ((np.abs(v6g)**2) * v6g))                    # N_w(R_g psi)
    rhs_m = right_rep_on_coeffs(g0) @ coeffs6(w6 * ((np.abs(v6)**2) * v6))  # R_g N_w(psi)
    err_m = np.abs(lhs_m - rhs_m).max() / np.abs(rhs_m).max()
    report("C3 cubic equivariance", green, err_m > 1e-3,
           f"(parent {err:.1e}; corrected-sides mutation under 6N rule {err_m:.1e})")
C3()

# ---------- C4: complete-level necessity, real parent and real mutation ----------
def C4():
    g = rand_su2_quat()
    def resid(drop):
        Yd  = modes(X4, N, drop_cols=drop)
        Ydg = modes(right_translate_nodes(X4, g), N, drop_cols=drop)
        Gd  = np.linalg.pinv(Yd.conj().T @ (W4[:, None] * Yd))
        cd  = rng.standard_normal(Yd.shape[1]) + 1j*rng.standard_normal(Yd.shape[1])
        proj = Gd @ (Yd.conj().T @ (W4 * (Ydg @ cd)))
        return np.abs(Yd @ proj - Ydg @ cd).max() / np.abs(Yd @ cd).max()
    err_green = resid(None)                # complete levels: translation stays in the span
    err_mut   = resid([nmodes - 1])        # one dropped column: it leaks
    report("C4 complete-level necessity", err_green < 1e-10, err_mut > 1e-3,
           f"(complete {err_green:.1e}, dropped-column {err_mut:.1e})")
C4()

# ---------- per-element 2I characters (mode_count's construction, inlined) ----------
CLASS_ANGLES = np.array([0, np.pi, np.pi/2, np.pi/3, 2*np.pi/3,
                         np.pi/5, 2*np.pi/5, 3*np.pi/5, 4*np.pi/5])
GALOIS = {5:7, 7:5, 6:8, 8:6}              # indices into CLASS_ANGLES: pi/5<->3pi/5, 2pi/5<->4pi/5
def snap(q):
    th = np.arccos(np.clip(q[0], -1, 1))
    k = int(np.argmin(np.abs(CLASS_ANGLES - th)))
    assert abs(CLASS_ANGLES[k] - th) < 1e-9, "angle off every 2I class"
    return k
def chiV(n, k):
    th = CLASS_ANGLES[k]
    if k == 0: return float(n + 1)
    if k == 1: return float((n+1) * (-1)**n)
    return np.sin((n+1)*th) / np.sin(th)
IRREPS8 = {"R0": lambda k: 1.0,          "R1": lambda k: chiV(1, k),
           "R2": lambda k: chiV(1, GALOIS.get(k, k)), "R3": lambda k: chiV(2, k),
           "R4": lambda k: chiV(2, GALOIS.get(k, k)),
           "R6": lambda k: chiV(3, k),   "R7": lambda k: chiV(4, k),
           "R8": lambda k: chiV(5, k)}
IRREP_NAMES = ["R0","R1","R2","R3","R4","R5","R6","R7","R8"]
D_RHO = {"R1":1, "R3":2, "R6":3, "R7":4, "R8":5, "R4":6, "R5":6, "R2":7}
def P1():
    """Eight characters from the SU(2) formulas; R5 DERIVED in-room, never typed: column
    orthogonality gives |chi(g)|^2 = |C(g)| - sum of the other eight, and a sign search over
    the ambiguous classes has a UNIQUE row orthonormal to all eight formula-built characters."""
    global KLASS, CHARTAB
    try:
        KLASS = [snap(q) for q in G120]
    except AssertionError as e:
        report("P1 character table", False, False, f"(class snapping failed: {e})")
        return False
    import collections
    csize = collections.Counter(KLASS)
    cent = {k: 120 // csize[k] for k in csize}
    absq = {k: cent[k] - sum(IRREPS8[r](k)**2 for r in IRREPS8) for k in csize}
    amb = [k for k in sorted(csize) if absq[k] > 0.5 and k != 0]
    others_vec = {r: np.array([IRREPS8[r](k) for k in KLASS]) for r in IRREPS8}
    hits = []
    for bits in range(2 ** len(amb)):
        row = {k: float(np.sqrt(max(absq[k], 0.0))) for k in csize}
        for j, k in enumerate(amb):
            if (bits >> j) & 1: row[k] = -row[k]
        vec = np.array([row[k] for k in KLASS])
        if (all(abs(np.mean(vec * others_vec[r])) < 1e-9 for r in IRREPS8)
                and abs(np.mean(vec * vec) - 1.0) < 1e-9):
            hits.append(vec)
    unique = (len(hits) == 1)
    CHARTAB = {r: others_vec[r] for r in IRREPS8}
    CHARTAB["R5"] = hits[0] if unique else np.zeros(120)
    dev = np.abs(np.array([[np.mean(CHARTAB[a]*CHARTAB[b]) for b in IRREP_NAMES]
                           for a in IRREP_NAMES]) - np.eye(9)).max()
    bad = {r: v.copy() for r, v in CHARTAB.items()}
    bad["R5"][3] += 1.0                                  # perturb one element value of R5
    dev_m = np.abs(np.array([[np.mean(bad[a]*bad[b]) for b in IRREP_NAMES]
                             for a in IRREP_NAMES]) - np.eye(9)).max()
    report("P1 character table", unique and dev < 1e-9, dev_m > 1e-3,
           f"(R5 derived, sign solution unique over {2**len(amb)} candidates; "
           f"orthonormality dev {dev:.1e}; perturbed-R5 mutation {dev_m:.1e})")
    return unique and dev < 1e-9
P1()

# ---------- C5: multiplicity one at every scored level, BOTH halves ----------
def C5():
    chi_level = lambda n: np.array([chiV(n, k) for k in KLASS])
    molien_levels = [0, 12, 20, 24, 30, 32, 36]
    char_mult = {n: round(float(np.mean(chi_level(n)))) for n in molien_levels}
    reyn_mult = {}
    for n in molien_levels:
        s = np.linalg.svd(np.mean([pi(n, q) for q in G120], axis=0), compute_uv=False)
        reyn_mult[n] = int((s > 0.5).sum())
        assert s[s <= 0.5].max(initial=0.0) < 1e-6, "Reynolds spectrum not cleanly split"
    half1 = all(char_mult[n] == reyn_mult[n] == 1 for n in molien_levels)
    SCAN = 10                                            # fixed independent range, all d_rho < 10
    disc = {}
    for r in D_RHO:
        ms = [round(float(np.mean(chi_level(n) * CHARTAB[r]))) for n in range(SCAN + 1)]
        found = next((n for n, m in enumerate(ms) if m > 0), None)
        disc[r] = (found, ms[found] if found is not None else 0)
    half2 = all(disc[r] == (D_RHO[r], 1) for r in D_RHO)
    green = half1 and half2
    G24 = G120[:24]
    ok24 = all(np.abs(G24 - qmul(G24[i], G24[j])).sum(axis=1).min() < 1e-9
               for i in range(24) for j in range(24))
    mult12_2T = round(float(np.mean([chiV(12, snap(q)) for q in G24])))   # 2T subset of 2I
    report("C5 multiplicity both halves", green, ok24 and mult12_2T != 1,
           f"(Molien half {char_mult}; first occurrence DISCOVERED over n = 0..10, equals the "
           f"pinned d_rho with mult 1, all eight; 2T mutation at n=12 -> {mult12_2T})")
C5()

# ---------- generators of 2I, for the commutant census ----------
def generate(seed_idx):
    S = {tuple(np.round(G120[i], 9)) for i in seed_idx}
    frontier = [G120[i] for i in seed_idx]
    while frontier:
        nxt = []
        for a in frontier:
            for i in seed_idx:
                p = qmul(a, G120[i]); tp = tuple(np.round(p, 9))
                if tp not in S: S.add(tp); nxt.append(p)
        frontier = nxt
    return S
ord10s = [i for i, q in enumerate(G120) if order_of(q) == 10]
ord4s  = [i for i, q in enumerate(G120) if order_of(q) == 4]
GENS = next(([a, b] for a in ord10s[:3] for b in ord4s
             if len(generate([a, b])) == 120), None)

# ---------- C6: 2I-commutant census (replaces the Schur-trivial SU(2) version) ----------
def commutant(mats):
    d = mats[0].shape[0]
    rows = [np.kron(np.eye(d), A) - np.kron(A.T, np.eye(d)) for A in mats]
    M = np.vstack(rows)
    s = np.linalg.svd(M, compute_uv=False)
    nullity = int((s < 1e-8 * s.max()).sum())
    _, _, Vh = np.linalg.svd(M)
    basis = [Vh[-(k+1)].reshape(d, d) for k in range(nullity)]
    comm = max((np.abs(A @ B - B @ A).max() for A in basis for B in basis), default=0.0)
    return nullity, comm
def C6():
    if GENS is None:
        report("C6 2I-commutant census", False, False, "(no generating pair found)")
        return
    expect = {1:1, 2:1, 3:1, 4:1, 5:1, 6:2, 7:2, 12:4}
    dims, comms, chis = {}, {}, {}
    for l in expect:
        mats = [pi(l, G120[i]) for i in GENS]
        dims[l], comms[l] = commutant(mats)
        chis[l] = round(float(np.mean(np.array([chiV(l, k) for k in KLASS])**2)))
    green = all(dims[l] == expect[l] == chis[l] for l in expect) \
            and max(comms.values()) < 1e-8
    # mutation: pi_1 + pi_1 of 2I has commutant M_2(C): dim 4 and NONcommutative, so the
    # "commutative iff multiplicity-free" reading must fire on it.
    mats_m = [block_diag(pi(1, G120[i]), pi(1, G120[i])) for i in GENS]
    dim_m, comm_m = commutant(mats_m)
    report("C6 2I-commutant census", green, dim_m == 4 and comm_m > 1e-3,
           f"(dims {dims} = <chi,chi> {chis}, max commutator {max(comms.values()):.1e}; "
           f"pi1+pi1 mutation dim {dim_m}, commutator {comm_m:.1e})")
C6()

# ---------- C7: assembled Laplacian commutes with the assembled right action ----------
def C7():
    Lam = np.diag(lam)
    err = np.abs(Lam @ R_num - R_num @ Lam).max() / max(lam.max(), 1.0)
    green = err < 1e-10
    lam_m = lam.copy()
    top = slice(nmodes - (N+1)**2, nmodes)            # perturb WITHIN the top level
    lam_m[top] = lam_m[top] + np.linspace(0, 0.5*lam.max(), (N+1)**2)
    Lm = np.diag(lam_m)
    err_m = np.abs(Lm @ R_num - R_num @ Lm).max() / max(lam.max(), 1.0)
    report("C7 level-diagonal spectrum vs assembled right action", green, err_m > 1e-3,
           f"(commutator {err:.1e}; intra-level-broken mutation {err_m:.1e})")
C7()

print(f"\n  wall {time.time()-t0:.1f}s")
if all(results):
    print("  VERDICT: ALL BRIDGES ARMED. The repo scalar primitives realize the derived right")
    print("  action; every mutation went red on its own check. The W_rho sector bases and the")
    print("  real Galerkin system remain to be qualified by protocol gates 3 and 5.")
    sys.exit(0)
print("  VERDICT: AT LEAST ONE BRIDGE FAILED OR AN ARM DID NOT FIRE. See lines above.")
sys.exit(1)
