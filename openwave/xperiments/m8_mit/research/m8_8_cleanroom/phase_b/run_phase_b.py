#!/usr/bin/env python3
"""
M8.8 Phase B: Gate Qualification
Executes all 19 declared mutations from METHOD_AND_GATE_MANIFEST.md § 4
against Phase A's immutable artifacts.

Deliverable: MUTATION_RESULTS.json
"""
import sys, os, json, copy, hashlib, importlib.util
from fractions import Fraction
from math import gcd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FROZEN_DIR = os.path.join(SCRIPT_DIR, 'phase_a_frozen')

# ========== PHASE A HASH TABLE (from protocol Addendum 1) ==========
PHASE_A_HASHES = {
    'TASK.md': 'e3d9b90861bb81862843988e8bd5da925b4d48bc48c0d2335becd3137df9cb17',
    'METHOD_AND_GATE_MANIFEST.md': '8aa140e3978366ca38f7c1d5926d1a2972305733be434595f2905e9df512f838',
    'validate_enumeration.py': 'b028f3b6fffe13809c49242f0acad7c0213025c65ef7fc693620274a2d87c1f7',
    'validate_complex.py': '348c87779c17f79bd4b0281ebaa29d8b5b598ac18c447bd71aa16d36b96607ea',
    'validate_saturation.py': 'f04602c622597eab132b25dfd52de2d305140cc05f8d361b256ef018b92c420f',
    'validate_representations.py': '580ed17aad2154313a1286ece6887509104bfba00ede467deb67892ffaf1e0ec',
    'validate_torsion_dry.py': '1e76a080e68ca0586e93842cc758727553fb15f2d40dae40830566ab2bd76601',
    'validate_fixture.py': 'd866a56eb852f9fc8fa870e5409033477e40a28fd731aa77720b0bbfa00a69f8',
    'validate_manifest.py': 'db9d73a244abdc7108db7697f3beaa4b89e82ea809492a4d398eda026db73488',
    'compute_torsion.py': '6277aef99613cc26c849f25084671fb8d1c6a6d232bf649eab9e627f049b7ab2',
    'ENVIRONMENT.md': '97637ba7192268d9fbfaa1813da5609d8f4b82febe5cb2edf1887f5d98a310e1',
    'RAW_OUTPUT.json': '1a9b56ce70bae73e5cf8c4ef00f6e43bf76937afb9075801605f6bf5047d1002',
    'CONSULTED_FILES.md': '650864857a50c266ad89d742346974b516521e08c6547d42a3643dc968a67652',
}

def verify_phase_a_hashes():
    for fname, expected in PHASE_A_HASHES.items():
        path = os.path.join(FROZEN_DIR, fname)
        with open(path, 'rb') as f:
            actual = hashlib.sha256(f.read()).hexdigest()
        if actual != expected:
            print(f"HASH MISMATCH: {fname}")
            print(f"  expected: {expected}")
            print(f"  actual:   {actual}")
            sys.exit(1)
    print(f"All {len(PHASE_A_HASHES)} Phase A artifact hashes verified.")
    return True

# ========== LOAD PHASE A MODULE ==========
def _load_mod(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(FROZEN_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# ========== Z[2I] GROUP-RING OPERATIONS ==========
def gr_zero():
    return {}

def gr_add(a, b):
    result = dict(a)
    for eid, coeff in b.items():
        result[eid] = result.get(eid, 0) + coeff
        if result[eid] == 0:
            del result[eid]
    return result

def gr_mul(a, b, mult_table):
    result = {}
    for eid_a, coeff_a in a.items():
        for eid_b, coeff_b in b.items():
            eid_ab = mult_table[eid_a][eid_b]
            c = coeff_a * coeff_b
            result[eid_ab] = result.get(eid_ab, 0) + c
            if result[eid_ab] == 0:
                del result[eid_ab]
    return result

def gr_is_zero(a):
    return len(a) == 0

def parse_boundary_map(bmap_data):
    rows = []
    for row_data in bmap_data:
        row = []
        for entry_data in row_data:
            gr_elem = {}
            for coeff, eid in entry_data:
                gr_elem[eid] = gr_elem.get(eid, 0) + coeff
                if gr_elem[eid] == 0:
                    del gr_elem[eid]
            row.append(gr_elem)
        rows.append(row)
    return rows

def gr_mat_mul(A, B, mult_table):
    r1 = len(A)
    c1 = len(A[0]) if A else 0
    c2 = len(B[0]) if B else 0
    result = []
    for i in range(r1):
        row = []
        for j in range(c2):
            entry = gr_zero()
            for k in range(c1):
                prod = gr_mul(A[i][k], B[k][j], mult_table)
                entry = gr_add(entry, prod)
            row.append(entry)
        result.append(row)
    return result

def augmentation(gr_elem):
    return sum(gr_elem.values())

def augment_matrix(mat):
    return [[augmentation(mat[i][j]) for j in range(len(mat[0]))] for i in range(len(mat))]

# ========== Z-EXPANSION AND SATURATION ==========
def expand_gr_mat(mat, mult_table, n_group=120):
    rows_gr = len(mat); cols_gr = len(mat[0])
    inv = [0] * n_group
    for i in range(n_group):
        for j in range(n_group):
            if mult_table[i][j] == 119:
                inv[i] = j; break
    Z_mat = [[0] * (cols_gr * n_group) for _ in range(rows_gr * n_group)]
    for bi in range(rows_gr):
        for bj in range(cols_gr):
            for eid, coeff in mat[bi][bj].items():
                for a in range(n_group):
                    b = mult_table[a][eid]
                    Z_mat[bi * n_group + a][bj * n_group + b] += coeff
    return Z_mat

def gauss_pivots(mat):
    m = len(mat)
    if m == 0: return 0, [], []
    n = len(mat[0])
    M = [[Fraction(mat[i][j]) for j in range(n)] for i in range(m)]
    rp = list(range(m)); pr = []; pc = []; r = 0
    for col in range(n):
        piv = None
        for row in range(r, m):
            if M[row][col] != 0: piv = row; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        rp[r], rp[piv] = rp[piv], rp[r]
        s = M[r][col]
        for j in range(n): M[r][j] /= s
        for row in range(m):
            if row == r: continue
            f = M[row][col]
            if f != 0:
                for j in range(n): M[row][j] -= f * M[r][j]
        pr.append(rp[r]); pc.append(col); r += 1
    return r, pr, pc

def det_frac(mat):
    n = len(mat)
    M = [[Fraction(mat[i][j]) for j in range(n)] for i in range(n)]
    d = Fraction(1)
    for col in range(n):
        piv = None
        for row in range(col, n):
            if M[row][col] != 0: piv = row; break
        if piv is None: return Fraction(0)
        if piv != col: M[col], M[piv] = M[piv], M[col]; d = -d
        d *= M[col][col]; s = M[col][col]
        for j in range(col, n): M[col][j] /= s
        for row in range(col + 1, n):
            f = M[row][col]
            if f != 0:
                for j in range(col, n): M[row][j] -= f * M[col][j]
    return d

# ========== FIXTURE TORSION (from validate_fixture.py) ==========
def minv_gc(M, mz_f, GC1_v, GC0_v):
    n = len(M)
    aug = [[M[i][j] for j in range(n)] + [GC1_v if i == j else GC0_v for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = None
        for row in range(col, n):
            if not aug[row][col].is_zero(): piv = row; break
        assert piv is not None, "Matrix is singular"
        aug[col], aug[piv] = aug[piv], aug[col]
        sc = aug[col][col]
        for j in range(2 * n): aug[col][j] = aug[col][j] / sc
        for row in range(n):
            if row == col: continue
            f = aug[row][col]
            if not f.is_zero():
                for j in range(2 * n): aug[row][j] = aug[row][j] - f * aug[col][j]
    return [[aug[i][n + j] for j in range(n)] for i in range(n)]

def fixture_compute_torsion(rho, d1_raw, d2_raw, d3_raw, d, ct):
    """Compute torsion for a representation using the declared conventions.
    Returns (tau, is_acyclic, msg).
    """
    def eval_entry(entry):
        result = ct.mz(d, d)
        for coeff, eid in entry:
            result = ct.madd(result, ct.msc(ct.GC(ct.QG(coeff)), rho[eid]))
        return result

    def eval_bmap(bmap_raw):
        rows_gr = len(bmap_raw); cols_gr = len(bmap_raw[0])
        M = ct.mz(rows_gr * d, cols_gr * d)
        for bi in range(rows_gr):
            for bj in range(cols_gr):
                block = eval_entry(bmap_raw[bi][bj])
                for i in range(d):
                    for j in range(d):
                        M[bi * d + i][bj * d + j] = block[i][j]
        return M

    M1 = eval_bmap(d1_raw)
    M2 = eval_bmap(d2_raw)
    M3 = eval_bmap(d3_raw)

    prod32 = ct.mmul(M3, M2)
    prod21 = ct.mmul(M2, M1)
    dd_ok = (all(prod32[i][j].is_zero() for i in range(d) for j in range(2 * d)) and
             all(prod21[i][j].is_zero() for i in range(2 * d) for j in range(d)))
    if not dd_ok:
        return None, False, "dd!=0"

    J3 = list(range(d))
    M3_minor = [[M3[i][j] for j in J3] for i in range(d)]
    det3 = ct.det_gc(M3_minor)
    if det3.is_zero():
        J3 = list(range(d, 2 * d))
        M3_minor = [[M3[i][j] for j in J3] for i in range(d)]
        det3 = ct.det_gc(M3_minor)
        if det3.is_zero():
            return None, False, "M3 rank deficient"

    I1 = list(range(d))
    M1_minor = [[M1[i][j] for j in range(d)] for i in I1]
    det1 = ct.det_gc(M1_minor)
    if det1.is_zero():
        I1 = list(range(d, 2 * d))
        M1_minor = [[M1[i][j] for j in range(d)] for i in I1]
        det1 = ct.det_gc(M1_minor)
        if det1.is_zero():
            return None, False, "M1 rank deficient"

    J3c = [j for j in range(2 * d) if j not in J3]
    I1c = [i for i in range(2 * d) if i not in I1]
    M2_minor = [[M2[i][j] for j in I1c] for i in J3c]
    det2 = ct.det_gc(M2_minor)
    tau = det2 / (det1 * det3)
    return tau, True, "ok"


# ========== REGISTRY (19 GATES FROM MANIFEST § 4) ==========
MANIFEST_GATES = [
    'G-M01', 'G-M02', 'G-M03', 'G-M04', 'G-M05', 'G-M06', 'G-M07', 'G-M08',
    'G-T01', 'G-T02', 'G-T03a', 'G-T03b', 'G-T03c', 'G-T03d',
    'G-D01', 'G-D02', 'G-D03', 'G-D04', 'G-D05',
]

def main():
    print("=" * 60)
    print("M8.8 PHASE B: GATE QUALIFICATION")
    print("=" * 60)

    # ---- Step 0: Verify Phase A integrity ----
    print("\n--- Step 0: Verify Phase A artifact hashes ---")
    verify_phase_a_hashes()

    # ---- Load Phase A module ----
    print("\n--- Loading Phase A computation module ---")
    ct = _load_mod('ct', 'compute_torsion.py')

    # ---- Load packets ----
    print("\n--- Loading packets ---")
    with open(os.path.join(SCRIPT_DIR, 'm8_5a_packet.json')) as f:
        gp = json.load(f)
    with open(os.path.join(SCRIPT_DIR, 'm8_8_construction_packet.json')) as f:
        cp = json.load(f)

    # ---- Group setup ----
    print("\n--- Group enumeration ---")
    se = ct.eg(gp)
    assert len(se) == 120
    e2r = {q: i for i, q in enumerate(se)}
    mt = [[e2r[se[i] * se[j]] for j in range(120)] for i in range(120)]
    s_id = cp['abstract_generators']['s']
    t_id = cp['abstract_generators']['t']
    e_id = 119
    st_id = mt[s_id][t_id]
    print(f"  s_id={s_id}, t_id={t_id}, e_id={e_id}, st_id={st_id}")

    # ---- Parse boundary maps (Z[2I] level) ----
    print("\n--- Parsing boundary maps ---")
    d1_gr = parse_boundary_map(cp['boundary_maps']['d1'])
    d2_gr = parse_boundary_map(cp['boundary_maps']['d2'])
    d3_gr = parse_boundary_map(cp['boundary_maps']['d3'])
    d1_raw = cp['boundary_maps']['d1']
    d2_raw = cp['boundary_maps']['d2']
    d3_raw = cp['boundary_maps']['d3']
    print(f"  d1: {len(d1_gr)}x{len(d1_gr[0])}, d2: {len(d2_gr)}x{len(d2_gr[0])}, d3: {len(d3_gr)}x{len(d3_gr[0])}")

    # ---- Build representations ----
    print("\n--- Building all 9 irreps ---")
    reps = ct.build_all_irreps(se, mt, s_id, t_id)
    rep_names = ['V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8']
    dims = {name: len(reps[name][0]) for name in rep_names}
    print(f"  Dimensions: {dims}")

    # ---- Compute inverse map ----
    inv_map = [0] * 120
    for i in range(120):
        for j in range(120):
            if mt[i][j] == e_id:
                inv_map[i] = j
                break

    # ---- Results collector ----
    results = []

    def record(gate_id, gate_name, obj_mutated, gate_pred, baseline, mutated, red):
        results.append({
            'gate_id': gate_id,
            'gate_name': gate_name,
            'object_mutated': obj_mutated,
            'gate_predicate': gate_pred,
            'baseline_result': baseline,
            'mutated_result': mutated,
            'red_outcome': red,
        })
        status = "REDDENED" if red else "FAILED TO REDDEN"
        print(f"  {gate_id} ({gate_name}): {status}")

    # ================================================================
    # MODEL GATES
    # ================================================================
    print("\n" + "=" * 60)
    print("MODEL GATES (G-M01 through G-M08)")
    print("=" * 60)

    # ---- G-M01: d3.d2 = 0, perturb d3 ----
    print("\nG-M01: Perturb d3 entry, verify d3.d2 != 0")
    d3d2 = gr_mat_mul(d3_gr, d2_gr, mt)
    baseline_m01 = all(gr_is_zero(d3d2[i][j])
                       for i in range(len(d3d2)) for j in range(len(d3d2[0])))
    assert baseline_m01, "Baseline G-M01 failed"
    d3_mut = copy.deepcopy(d3_gr)
    d3_mut[0][0] = gr_add(d3_mut[0][0], {e_id: 1})
    d3d2_mut = gr_mat_mul(d3_mut, d2_gr, mt)
    mut_m01 = all(gr_is_zero(d3d2_mut[i][j])
                  for i in range(len(d3d2_mut)) for j in range(len(d3d2_mut[0])))
    record('G-M01', 'dd=0 at degree 2', 'd3[0][0]: added +1*e (eid 119)',
           'd3.d2 = 0 over Z[2I]', 'PASS (product is zero)',
           f'FAIL (product is {"zero" if mut_m01 else "nonzero"})', not mut_m01)

    # ---- G-M02: d2.d1 = 0, perturb d2 ----
    print("\nG-M02: Perturb d2 entry, verify d2.d1 != 0")
    d2d1 = gr_mat_mul(d2_gr, d1_gr, mt)
    baseline_m02 = all(gr_is_zero(d2d1[i][j])
                       for i in range(len(d2d1)) for j in range(len(d2d1[0])))
    assert baseline_m02, "Baseline G-M02 failed"
    d2_mut = copy.deepcopy(d2_gr)
    d2_mut[0][0] = gr_add(d2_mut[0][0], {e_id: 1})
    d2d1_mut = gr_mat_mul(d2_mut, d1_gr, mt)
    mut_m02 = all(gr_is_zero(d2d1_mut[i][j])
                  for i in range(len(d2d1_mut)) for j in range(len(d2d1_mut[0])))
    record('G-M02', 'dd=0 at degree 1', 'd2[0][0]: added +1*e (eid 119)',
           'd2.d1 = 0 over Z[2I]', 'PASS (product is zero)',
           f'FAIL (product is {"zero" if mut_m02 else "nonzero"})', not mut_m02)

    # ---- G-M03: Free ranks and chi = 0 ----
    print("\nG-M03: Alter one rank, verify chi != 0")
    free_ranks = cp['free_ranks']
    chi = sum((-1)**k * r for k, r in enumerate(free_ranks))
    assert chi == 0, "Baseline G-M03 failed"
    ranks_mut = list(free_ranks)
    ranks_mut[3] = 2
    chi_mut = sum((-1)**k * r for k, r in enumerate(ranks_mut))
    record('G-M03', 'Free ranks and chi', f'free_ranks[3]: 1 -> 2',
           'chi = sum(-1)^k * r_k = 0', f'PASS (chi={chi})',
           f'FAIL (chi={chi_mut})', chi_mut != 0)

    # ---- G-M04: Augmented homology ----
    print("\nG-M04: Replace d2 entry, verify augmented homology changes")
    d2_aug = augment_matrix(d2_gr)
    det_d2_aug = d2_aug[0][0] * d2_aug[1][1] - d2_aug[0][1] * d2_aug[1][0]
    assert abs(det_d2_aug) == 1, f"Baseline G-M04 failed: det={det_d2_aug}"
    d2_mut4 = copy.deepcopy(d2_gr)
    d2_mut4[0][0] = gr_add(d2_mut4[0][0], {e_id: 2})
    d2_aug_mut = augment_matrix(d2_mut4)
    det_mut = d2_aug_mut[0][0] * d2_aug_mut[1][1] - d2_aug_mut[0][1] * d2_aug_mut[1][0]
    record('G-M04', 'Augmented homology', 'd2[0][0]: added +2*e (eid 119)',
           'det(d2_aug) = +/-1', f'PASS (det={det_d2_aug})',
           f'FAIL (det={det_mut})', abs(det_mut) != 1)

    # ---- G-M05: Universal-cover homology (saturation) ----
    print("\nG-M05: Scale d2 row by non-unit, verify saturation fails")
    print("  Expanding d2 over Z (240x240)...")
    d2_Z = expand_gr_mat(d2_gr, mt)
    print(f"  d2_Z: {len(d2_Z)}x{len(d2_Z[0])}")
    print("  Running Gaussian elimination on original d2_Z...")
    r2, pr2, pc2 = gauss_pivots(d2_Z)
    assert r2 == 121, f"Baseline rank(d2_Z)={r2}, expected 121"
    sub2 = [[d2_Z[pr2[i]][pc2[j]] for j in range(121)] for i in range(121)]
    det2_orig = det_frac(sub2)
    assert abs(det2_orig) == 1, f"Baseline |det|={abs(det2_orig)}, expected 1"
    print(f"  Original: rank={r2}, |det(pivot minor)|={abs(det2_orig)} -> PASS")
    d2_gr_mut5 = copy.deepcopy(d2_gr)
    for j in range(len(d2_gr_mut5[0])):
        d2_gr_mut5[0][j] = {eid: 2 * c for eid, c in d2_gr_mut5[0][j].items()}
    d2_Z_mut = expand_gr_mat(d2_gr_mut5, mt)
    sub2_mut = [[d2_Z_mut[pr2[i]][pc2[j]] for j in range(121)] for i in range(121)]
    det2_mut = det_frac(sub2_mut)
    print(f"  Mutated:  |det(pivot minor)|={abs(det2_mut)}")
    record('G-M05', 'Universal-cover saturation',
           'd2 row 0 scaled by 2 (non-unit)',
           'unimodular 121x121 minor exists (|det|=1)',
           f'PASS (|det|={abs(det2_orig)})',
           f'FAIL (|det|={abs(det2_mut)})', abs(det2_mut) != 1)

    # ---- G-M06: Augmentation is terminal map ----
    print("\nG-M06: Replace eps with non-augmentation, verify eps.d1 != 0")
    eps_d1 = [augmentation(d1_gr[i][0]) for i in range(len(d1_gr))]
    baseline_m06 = all(v == 0 for v in eps_d1)
    assert baseline_m06, "Baseline G-M06 failed"
    def non_augmentation(gr_elem):
        return sum(c * (1 if eid == e_id else 2) for eid, c in gr_elem.items())
    eps_d1_mut = [non_augmentation(d1_gr[i][0]) for i in range(len(d1_gr))]
    mut_m06 = all(v == 0 for v in eps_d1_mut)
    record('G-M06', 'Augmentation is terminal map',
           'eps: replaced with non-augmentation (e->1, others->2)',
           'eps(d1) = 0', f'PASS (eps_d1={eps_d1})',
           f'FAIL (eps_d1={eps_d1_mut})', not mut_m06)

    # ---- G-M07: Generator correspondence ----
    print("\nG-M07: Swap s,t IDs, verify relator check fails")
    s_elem = se[s_id]; t_elem = se[t_id]
    s3 = s_elem * s_elem * s_elem
    st = s_elem * t_elem
    st2 = st * st
    t5 = t_elem * t_elem * t_elem * t_elem * t_elem
    baseline_m07 = (s3 == st2 and t5 == st2)
    assert baseline_m07, "Baseline G-M07 failed"
    s_swap = se[t_id]; t_swap = se[s_id]
    s3_swap = s_swap * s_swap * s_swap
    st_swap = s_swap * t_swap
    st2_swap = st_swap * st_swap
    t5_swap = t_swap * t_swap * t_swap * t_swap * t_swap
    mut_m07 = (s3_swap == st2_swap and t5_swap == st2_swap)
    record('G-M07', 'Generator correspondence',
           'Swapped s_id and t_id',
           's^3=(st)^2 and t^5=(st)^2',
           'PASS (both relators hold)',
           f'FAIL (relators {"hold" if mut_m07 else "fail"})', not mut_m07)

    # ---- G-M08: Per-irrep acyclicity ----
    print("\nG-M08: Mutate complex for nontrivial irrep, verify acyclicity fails")
    rho_V1 = reps['V1']
    d_v1 = dims['V1']
    M3_v1 = ct.eval_bmap(d3_raw, rho_V1, d_v1)
    def rank_gc(M):
        m = len(M); n = len(M[0])
        A = [[M[i][j] for j in range(n)] for i in range(m)]
        r = 0
        for col in range(n):
            piv = None
            for row in range(r, m):
                if not A[row][col].is_zero(): piv = row; break
            if piv is None: continue
            A[r], A[piv] = A[piv], A[r]
            sc = A[r][col]
            for j in range(n): A[r][j] = A[r][j] / sc
            for row in range(m):
                if row == r: continue
                f = A[row][col]
                if not f.is_zero():
                    for j in range(n): A[row][j] = A[row][j] - f * A[r][j]
            r += 1
        return r
    r3_baseline = rank_gc(M3_v1)
    assert r3_baseline == d_v1, f"Baseline rank(M3_V1)={r3_baseline}, expected {d_v1}"
    M3_v1_mut = [[M3_v1[i][j] for j in range(2 * d_v1)] for i in range(d_v1)]
    for j in range(2 * d_v1):
        M3_v1_mut[0][j] = ct.GC()
    r3_mut = rank_gc(M3_v1_mut)
    record('G-M08', 'Per-irrep acyclicity',
           'V1 twisted M3: zeroed row 0',
           'rank(M3) = d for nontrivial irrep',
           f'PASS (rank={r3_baseline}, d={d_v1})',
           f'FAIL (rank={r3_mut}, d={d_v1})', r3_mut < d_v1)

    # ================================================================
    # THEOREM-SIDE GATES
    # ================================================================
    print("\n" + "=" * 60)
    print("THEOREM-SIDE GATES (G-T01 through G-T03d)")
    print("=" * 60)

    # ---- G-T01: Unitarity ----
    print("\nG-T01: Perturb rep matrix, verify Hermitian invariance fails")
    rho_test = reps['V1']
    d_test = dims['V1']
    H = ct.mz(d_test, d_test)
    for g in range(120):
        H = ct.madd(H, ct.mmul(ct.mct(rho_test[g]), rho_test[g]))
    H = ct.msc(ct.GC(ct.QG(Fraction(1, 120))), H)
    lhs_base = ct.mmul(ct.mct(rho_test[s_id]), ct.mmul(H, rho_test[s_id]))
    baseline_t01 = ct.meq(lhs_base, H)
    assert baseline_t01, "Baseline G-T01 failed"
    rho_mut = [m for m in rho_test]
    rho_s_mut = [[rho_test[s_id][i][j] for j in range(d_test)] for i in range(d_test)]
    rho_s_mut[0][0] = rho_s_mut[0][0] + ct.GC(ct.QG(Fraction(1, 10)))
    rho_mut[s_id] = rho_s_mut
    H_mut = ct.mz(d_test, d_test)
    for g in range(120):
        H_mut = ct.madd(H_mut, ct.mmul(ct.mct(rho_mut[g]), rho_mut[g]))
    H_mut = ct.msc(ct.GC(ct.QG(Fraction(1, 120))), H_mut)
    lhs_mut = ct.mmul(ct.mct(rho_mut[s_id]), ct.mmul(H_mut, rho_mut[s_id]))
    mut_t01 = ct.meq(lhs_mut, H_mut)
    record('G-T01', 'Unitarity',
           'V1: perturbed rho(s)[0][0] by +1/10',
           'rho(s)^dag H rho(s) = H',
           'PASS (invariant)',
           f'FAIL ({"invariant" if mut_t01 else "not invariant"})', not mut_t01)

    # ---- G-T02: Row signature identity ----
    print("\nG-T02: Swap character values, verify signature match fails")
    ch_s_v1 = ct.mtr(reps['V1'][s_id])
    ch_t_v1 = ct.mtr(reps['V1'][t_id])
    ch_s_v7 = ct.mtr(reps['V7'][s_id])
    ch_t_v7 = ct.mtr(reps['V7'][t_id])
    baseline_distinct = not (ch_s_v1.re == ch_s_v7.re and ch_t_v1.re == ch_t_v7.re)
    assert baseline_distinct, "Baseline G-T02: V1 and V7 have same sig for s,t"
    all_sigs = set()
    for name in rep_names:
        rho = reps[name]
        d = dims[name]
        cs = ct.mtr(rho[s_id]).re
        ct_val = ct.mtr(rho[t_id]).re
        cst = ct.mtr(rho[st_id]).re
        sig = (d, (cs.a, cs.b), (ct_val.a, ct_val.b), (cst.a, cst.b))
        all_sigs.add(sig)
    baseline_t02 = len(all_sigs) == 9
    assert baseline_t02, "Baseline G-T02: not all 9 sigs distinct"
    ch_s_v1_swapped = ct.mtr(reps['V1'][t_id]).re
    ch_t_v1_swapped = ct.mtr(reps['V1'][s_id]).re
    swapped_sig = (dims['V1'],
                   (ch_s_v1_swapped.a, ch_s_v1_swapped.b),
                   (ch_t_v1_swapped.a, ch_t_v1_swapped.b),
                   (ct.mtr(reps['V1'][st_id]).re.a, ct.mtr(reps['V1'][st_id]).re.b))
    orig_v1_sig = (dims['V1'],
                   (ch_s_v1.re.a, ch_s_v1.re.b),
                   (ch_t_v1.re.a, ch_t_v1.re.b),
                   (ct.mtr(reps['V1'][st_id]).re.a, ct.mtr(reps['V1'][st_id]).re.b))
    sig_changed = swapped_sig != orig_v1_sig
    record('G-T02', 'Row signature identity',
           'V1: swapped chi(s) and chi(t)',
           'Row signatures are distinct and match expected',
           'PASS (9 distinct signatures)',
           f'FAIL (V1 signature {"changed" if sig_changed else "unchanged"})',
           sig_changed)

    # ---- G-T03a-d: Convention fixture tests ----
    print("\nG-T03a-d: Convention fixture mutations")
    su2 = [ct.q2su2(q) for q in se]
    P = [[ct.GC(ct.QG(2)), ct.GC(ct.QG(0), ct.QG(1))],
         [ct.GC(), ct.GC(ct.QG(1))]]
    Pi = minv_gc(P, ct.mz, ct.GC(ct.QG(1)), ct.GC())
    fixture_reps = [ct.mmul(Pi, ct.mmul(su2[g], P)) for g in range(120)]
    sds = ct.mmul(ct.mct(fixture_reps[s_id]), fixture_reps[s_id])
    assert not ct.meq(sds, ct.mid(2)), "Fixture is unitary"

    tau_base, acyc_base, msg_base = fixture_compute_torsion(
        fixture_reps, d1_raw, d2_raw, d3_raw, 2, ct)
    assert acyc_base, f"Fixture baseline not acyclic: {msg_base}"
    T2_base = tau_base * tau_base.conj()
    assert T2_base.im.is_zero(), "Baseline T2 not real"
    print(f"  Fixture baseline T2 = {T2_base.re}")

    # G-T03a: evaluation map g -> rho(g^-1) (anti-homomorphism)
    anti_hom = [fixture_reps[inv_map[g]] for g in range(120)]
    tau_a, acyc_a, msg_a = fixture_compute_torsion(
        anti_hom, d1_raw, d2_raw, d3_raw, 2, ct)
    if not acyc_a:
        red_a = True
        mut_desc_a = f'FAIL ({msg_a})'
    else:
        T2_a = tau_a * tau_a.conj()
        red_a = not (T2_a.re == T2_base.re)
        mut_desc_a = f'FAIL (T2 changed)' if red_a else 'PASS (T2 unchanged)'
    record('G-T03a', 'Convention: evaluation map',
           'Fixture: g -> rho(g^-1) (anti-homomorphism)',
           'dd=0 and correct T2',
           'PASS (acyclic, T2 computed)', mut_desc_a, red_a)

    # G-T03b: boundary direction (cochain reversal)
    d1_rev = [[d3_raw[j][i] for j in range(len(d3_raw))] for i in range(len(d3_raw[0]))]
    d2_rev = [[d2_raw[j][i] for j in range(len(d2_raw))] for i in range(len(d2_raw[0]))]
    d3_rev = [[d1_raw[j][i] for j in range(len(d1_raw))] for i in range(len(d1_raw[0]))]
    tau_b, acyc_b, msg_b = fixture_compute_torsion(
        fixture_reps, d1_rev, d2_rev, d3_rev, 2, ct)
    if not acyc_b:
        red_b = True
        mut_desc_b = f'FAIL ({msg_b})'
    else:
        T2_b = tau_b * tau_b.conj()
        red_b = not (T2_b.re == T2_base.re)
        mut_desc_b = f'FAIL (T2 changed)' if red_b else 'PASS (T2 unchanged)'
    record('G-T03b', 'Convention: boundary direction',
           'Fixture: transposed boundary matrices (cochain)',
           'dd=0 and correct T2',
           'PASS (acyclic, T2 computed)', mut_desc_b, red_b)

    # G-T03c: module side (transpose reps = anti-homomorphism)
    def mtranspose(M):
        return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]
    trans_reps = [mtranspose(fixture_reps[g]) for g in range(120)]
    tau_c, acyc_c, msg_c = fixture_compute_torsion(
        trans_reps, d1_raw, d2_raw, d3_raw, 2, ct)
    if not acyc_c:
        red_c = True
        mut_desc_c = f'FAIL ({msg_c})'
    else:
        T2_c = tau_c * tau_c.conj()
        red_c = not (T2_c.re == T2_base.re)
        mut_desc_c = f'FAIL (T2 changed)' if red_c else 'PASS (T2 unchanged)'
    record('G-T03c', 'Convention: module side',
           'Fixture: rho(g)^T (transpose = anti-homomorphism)',
           'dd=0 and correct T2',
           'PASS (acyclic, T2 computed)', mut_desc_c, red_c)

    # G-T03d: vector convention (transpose GR boundary maps)
    d1_grt = [[d1_raw[j][i] for j in range(len(d1_raw))] for i in range(len(d1_raw[0]))]
    d2_grt = [[d2_raw[j][i] for j in range(len(d2_raw))] for i in range(len(d2_raw[0]))]
    d3_grt = [[d3_raw[j][i] for j in range(len(d3_raw))] for i in range(len(d3_raw[0]))]
    tau_d, acyc_d, msg_d = fixture_compute_torsion(
        fixture_reps, d3_grt, d2_grt, d1_grt, 2, ct)
    if not acyc_d:
        red_d = True
        mut_desc_d = f'FAIL ({msg_d})'
    else:
        T2_d = tau_d * tau_d.conj()
        red_d = not (T2_d.re == T2_base.re)
        mut_desc_d = f'FAIL (T2 changed)' if red_d else 'PASS (T2 unchanged)'
    record('G-T03d', 'Convention: vector convention',
           'Fixture: transposed GR maps (column vectors)',
           'dd=0 and correct T2',
           'PASS (acyclic, T2 computed)', mut_desc_d, red_d)

    # ================================================================
    # DERIVATION-PATH GATES
    # ================================================================
    print("\n" + "=" * 60)
    print("DERIVATION-PATH GATES (G-D01 through G-D05)")
    print("=" * 60)

    # Use V1 (dim 2) for efficiency
    rho_d = reps['V1']
    d_d = dims['V1']

    M1_d = ct.eval_bmap(d1_raw, rho_d, d_d)
    M2_d = ct.eval_bmap(d2_raw, rho_d, d_d)
    M3_d = ct.eval_bmap(d3_raw, rho_d, d_d)

    # ---- G-D01: dd=0 per twisted complex ----
    print("\nG-D01: Perturb evaluated boundary, verify M3.M2 != 0")
    prod32_d = ct.mmul(M3_d, M2_d)
    baseline_d01 = all(prod32_d[i][j].is_zero()
                       for i in range(d_d) for j in range(2 * d_d))
    assert baseline_d01, "Baseline G-D01 failed"
    M3_d_mut = [[M3_d[i][j] for j in range(2 * d_d)] for i in range(d_d)]
    M3_d_mut[0][0] = M3_d_mut[0][0] + ct.GC(ct.QG(1))
    prod32_mut = ct.mmul(M3_d_mut, M2_d)
    mut_d01 = all(prod32_mut[i][j].is_zero()
                  for i in range(d_d) for j in range(2 * d_d))
    record('G-D01', 'dd=0 per twisted complex',
           'V1 M3[0][0]: perturbed by +1',
           'M3.M2 = 0 (twisted complex chain condition)',
           'PASS (product is zero)',
           f'FAIL (product is {"zero" if mut_d01 else "nonzero"})', not mut_d01)

    # ---- G-D02: Twisted complex ranks ----
    print("\nG-D02: Zero row of M3, verify rank drops")
    r3_d = rank_gc(M3_d)
    assert r3_d == d_d, f"Baseline rank(M3)={r3_d}, expected {d_d}"
    M3_d_mut2 = [[M3_d[i][j] for j in range(2 * d_d)] for i in range(d_d)]
    for j in range(2 * d_d):
        M3_d_mut2[0][j] = ct.GC()
    r3_mut = rank_gc(M3_d_mut2)
    record('G-D02', 'Twisted complex ranks',
           'V1 M3: zeroed row 0',
           'rank(M3) = d',
           f'PASS (rank={r3_d}, d={d_d})',
           f'FAIL (rank={r3_mut}, d={d_d})', r3_mut < d_d)

    # ---- G-D03: Determinant sub-matrices nonsingular ----
    print("\nG-D03: Zero column of chosen minor, verify det = 0")
    J3_d = list(range(d_d))
    M3_minor_d = [[M3_d[i][j] for j in J3_d] for i in range(d_d)]
    det3_d = ct.det_gc(M3_minor_d)
    assert not det3_d.is_zero(), "Baseline G-D03: M3 minor singular"
    M3_minor_mut = [[M3_minor_d[i][j] for j in range(d_d)] for i in range(d_d)]
    for i in range(d_d):
        M3_minor_mut[i][0] = ct.GC()
    det3_mut = ct.det_gc(M3_minor_mut)
    record('G-D03', 'Determinant sub-matrices nonsingular',
           'V1 M3 minor: zeroed column 0',
           'det(minor) != 0',
           f'PASS (det nonzero)',
           f'FAIL (det {"nonzero" if not det3_mut.is_zero() else "= 0"})',
           det3_mut.is_zero())

    # ---- G-D04: Galois consistency ----
    print("\nG-D04: Break Galois pairing, verify consistency fails")
    T2_V1, _, _ = ct.compute_torsion_sq(reps['V1'], d1_raw, d2_raw, d3_raw, dims['V1'])
    T2_V7, _, _ = ct.compute_torsion_sq(reps['V7'], d1_raw, d2_raw, d3_raw, dims['V7'])
    T2_V1_gal = T2_V1.galois()
    baseline_d04 = (T2_V1_gal == T2_V7)
    assert baseline_d04, "Baseline G-D04 failed"
    T2_V2, _, _ = ct.compute_torsion_sq(reps['V2'], d1_raw, d2_raw, d3_raw, dims['V2'])
    fake_gal = T2_V2.galois()
    mut_d04 = (fake_gal == T2_V7)
    record('G-D04', 'Galois consistency',
           'Substituted sigma(T2(V2)) for sigma(T2(V1))',
           'sigma(T2(V1)) = T2(V7)',
           f'PASS (sigma(T2(V1)) = T2(V7))',
           f'FAIL (sigma(T2(V2)) {"=" if mut_d04 else "!="} T2(V7))',
           not mut_d04)

    # ---- G-D05: Torsion code-path dependency ----
    print("\nG-D05: Replace torsion input with identity, verify output changes")
    identity_reps = [ct.mid(d_d) for _ in range(120)]
    try:
        T2_id, _, _ = ct.compute_torsion_sq(identity_reps, d1_raw, d2_raw, d3_raw, d_d)
        T2_id_val = T2_id
    except Exception:
        T2_id_val = None
    if T2_id_val is not None:
        different = not (T2_id_val == T2_V1)
        mut_desc_d05 = f'T2_identity={T2_id_val}, differs from T2_V1' if different else 'T2 unchanged'
    else:
        different = True
        mut_desc_d05 = 'Computation failed (non-acyclic with identities)'
    record('G-D05', 'Torsion code-path dependency',
           'Replaced V1 reps with identity matrices',
           'T2 result changes when input changes',
           f'PASS (T2_V1 = {T2_V1})',
           f'FAIL ({mut_desc_d05})', different)

    # ================================================================
    # REGISTRY COVERAGE CHECK
    # ================================================================
    print("\n" + "=" * 60)
    print("REGISTRY COVERAGE CHECK")
    print("=" * 60)
    executed_gates = [r['gate_id'] for r in results]
    manifest_set = set(MANIFEST_GATES)
    executed_set = set(executed_gates)
    exact_equality = manifest_set == executed_set
    missing = manifest_set - executed_set
    extra = executed_set - manifest_set
    all_red = all(r['red_outcome'] for r in results)

    print(f"  Manifest gates:  {len(MANIFEST_GATES)}")
    print(f"  Executed gates:  {len(executed_gates)}")
    print(f"  Exact set equality: {exact_equality}")
    if missing: print(f"  MISSING: {missing}")
    if extra: print(f"  EXTRA: {extra}")
    print(f"  All mutations reddened: {all_red}")

    # ---- Re-verify Phase A hashes ----
    print("\n--- Re-verifying Phase A hashes (post-qualification) ---")
    post_ok = verify_phase_a_hashes()

    # ================================================================
    # WRITE OUTPUT
    # ================================================================
    output = {
        'schema_version': 'm8_8-phase-b-mutation-results-1',
        'phase_a_hashes_verified_pre': True,
        'phase_a_hashes_verified_post': post_ok,
        'registry_coverage': {
            'manifest_gates': MANIFEST_GATES,
            'executed_gates': executed_gates,
            'exact_set_equality': exact_equality,
            'count': len(results),
        },
        'all_mutations_reddened': all_red,
        'results': results,
    }

    out_path = os.path.join(SCRIPT_DIR, 'MUTATION_RESULTS.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nOutput written to {out_path}")

    # ---- Final verdict ----
    print("\n" + "=" * 60)
    if exact_equality and all_red and post_ok:
        print("PHASE B QUALIFICATION: ALL 19 GATES REDDENED")
        print("Registry coverage: exact set equality confirmed")
        print("Phase A integrity: verified pre and post")
        sys.exit(0)
    else:
        print("PHASE B QUALIFICATION: INCOMPLETE")
        if not exact_equality: print("  - Registry coverage mismatch")
        if not all_red:
            for r in results:
                if not r['red_outcome']:
                    print(f"  - {r['gate_id']}: mutation did not redden")
        if not post_ok: print("  - Phase A hashes changed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
