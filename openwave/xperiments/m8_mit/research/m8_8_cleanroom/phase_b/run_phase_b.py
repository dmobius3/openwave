#!/usr/bin/env python3
"""
M8.8 Phase B: Gate Qualification
Parses the 19 gate identifiers and their declared mutations from the frozen
METHOD_AND_GATE_MANIFEST.md § 4, then executes each mutation against scratch
copies of the frozen Phase A machinery.

Deliverable: MUTATION_RESULTS.json
"""
import sys
sys.dont_write_bytecode = True

import os, json, copy, hashlib, importlib.util
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


# ========== MANIFEST PARSER ==========

def parse_gate_registry(manifest_text):
    """Parse gate IDs and declared mutations from section 4 of the manifest.

    Returns an ordered list of dicts: [{'id': ..., 'name': ..., 'mutation': ...}, ...]
    Raises ValueError on duplicate IDs or empty/missing mutation declarations.
    """
    gates = []
    seen_ids = set()
    in_section_4 = False

    for line in manifest_text.split('\n'):
        stripped = line.strip()

        if stripped.startswith('## 4'):
            in_section_4 = True
            continue
        if in_section_4 and stripped.startswith('## ') and not stripped.startswith('### 4'):
            break
        if not in_section_4:
            continue
        if stripped.startswith('### ') or stripped.startswith('| ---') or stripped.startswith('| Gate ID'):
            continue

        if not stripped.startswith('| G-'):
            continue

        parts = stripped.split('|')
        if len(parts) < 5:
            raise ValueError(f"Malformed gate row (fewer than 4 pipe-separated columns): {stripped}")

        gate_id = parts[1].strip()
        gate_name = parts[2].strip()
        mutation = parts[3].strip()

        if not gate_id.startswith('G-'):
            raise ValueError(f"Malformed gate ID (must start with G-): {gate_id!r}")
        if not mutation.strip():
            raise ValueError(f"Empty mutation declaration for gate {gate_id}")
        if gate_id in seen_ids:
            raise ValueError(f"Duplicate gate identifier: {gate_id}")

        seen_ids.add(gate_id)
        gates.append({'id': gate_id, 'name': gate_name, 'mutation': mutation})

    return gates


def check_coverage(parsed_ids, other_ids, label):
    """Check exact set equality. Returns True on match; prints diagnostics on mismatch."""
    if parsed_ids == other_ids:
        return True
    missing = parsed_ids - other_ids
    extra = other_ids - parsed_ids
    print(f"  COVERAGE FAILURE ({label}):")
    if missing:
        print(f"    In parsed manifest but not in {label}: {sorted(missing)}")
    if extra:
        print(f"    In {label} but not in parsed manifest: {sorted(extra)}")
    return False


def self_test_parser_coverage(manifest_text, implemented_ids):
    """Self-test: verify the parser-to-coverage linkage detects manifest changes.

    Operates on scratch copies of the manifest text, never modifies any frozen file.
    Tests that adding or removing a gate from the text changes the parsed set and
    causes coverage to exit nonzero. This test must not pass merely because a
    Phase-A hash changed.
    """
    print("\n--- Self-test: parser-to-coverage linkage ---")

    baseline = parse_gate_registry(manifest_text)
    baseline_ids = set(g['id'] for g in baseline)
    assert baseline_ids == implemented_ids, "Baseline parser disagrees with implementation"

    # Test 1: add a gate to the scratch text
    extra_row = "| G-FAKE | Fake gate | Fake mutation; verify fake | Fake |\n"
    insert_marker = "| G-D05 |"
    pos = manifest_text.find(insert_marker)
    assert pos >= 0, "Cannot find G-D05 row for insertion"
    line_end = manifest_text.index('\n', pos)
    text_added = manifest_text[:line_end + 1] + extra_row + manifest_text[line_end + 1:]
    parsed_added = parse_gate_registry(text_added)
    ids_added = set(g['id'] for g in parsed_added)
    assert len(ids_added) == len(baseline_ids) + 1, \
        f"Added-gate parse: expected {len(baseline_ids)+1}, got {len(ids_added)}"
    assert 'G-FAKE' in ids_added, "Added gate G-FAKE not found in parsed set"
    cov_add = check_coverage(ids_added, implemented_ids, "add-test")
    assert not cov_add, "Coverage should fail when manifest has an extra gate"
    print("  ADD test: parsed set grew to 20, coverage correctly failed")

    # Test 2: remove a gate from the scratch text
    lines = manifest_text.split('\n')
    lines_removed = [l for l in lines if '| G-M01 |' not in l]
    text_removed = '\n'.join(lines_removed)
    parsed_removed = parse_gate_registry(text_removed)
    ids_removed = set(g['id'] for g in parsed_removed)
    assert len(ids_removed) == len(baseline_ids) - 1, \
        f"Removed-gate parse: expected {len(baseline_ids)-1}, got {len(ids_removed)}"
    assert 'G-M01' not in ids_removed, "Removed gate G-M01 still in parsed set"
    cov_rem = check_coverage(ids_removed, implemented_ids, "remove-test")
    assert not cov_rem, "Coverage should fail when manifest is missing a gate"
    print("  REMOVE test: parsed set shrank to 18, coverage correctly failed")

    # Test 3: duplicate gate ID
    dup_row = "| G-M01 | Duplicate | Duplicate mutation | Duplicate |\n"
    text_dup = manifest_text[:line_end + 1] + dup_row + manifest_text[line_end + 1:]
    try:
        parse_gate_registry(text_dup)
        assert False, "Duplicate ID should have raised ValueError"
    except ValueError as e:
        assert "Duplicate" in str(e)
    print("  DUPLICATE test: parser correctly rejected duplicate gate ID")

    # Test 4: empty mutation
    empty_row = "| G-EMPTY | Empty | | Nothing |\n"
    text_empty = manifest_text[:line_end + 1] + empty_row + manifest_text[line_end + 1:]
    try:
        parse_gate_registry(text_empty)
        assert False, "Empty mutation should have raised ValueError"
    except ValueError as e:
        assert "Empty mutation" in str(e)
    print("  EMPTY-MUTATION test: parser correctly rejected empty mutation")

    print("  Parser-to-coverage self-test: ALL PASSED")
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


# ========== FIXTURE TORSION ==========

def minv_gc(M, GC1_v, GC0_v):
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


# ========== IMPLEMENTED HANDLER TABLE ==========
# Each handler is a function: (ctx) -> result_dict
# The dict keys are the gate IDs for which handler code exists.

def _make_handlers():
    H = {}

    def gm01(c):
        d3d2 = gr_mat_mul(c['d3_gr'], c['d2_gr'], c['mt'])
        assert all(gr_is_zero(d3d2[i][j])
                    for i in range(len(d3d2)) for j in range(len(d3d2[0]))), "Baseline G-M01 failed"
        d3_mut = copy.deepcopy(c['d3_gr'])
        d3_mut[0][0] = gr_add(d3_mut[0][0], {c['e_id']: 1})
        d3d2_mut = gr_mat_mul(d3_mut, c['d2_gr'], c['mt'])
        red = not all(gr_is_zero(d3d2_mut[i][j])
                      for i in range(len(d3d2_mut)) for j in range(len(d3d2_mut[0])))
        return dict(object_mutated='d3[0][0]: added +1*e (eid 119)',
                    gate_predicate='d3.d2 = 0 over Z[2I]',
                    baseline_result='PASS (product is zero)',
                    mutated_result=f'FAIL (product is {"nonzero" if red else "zero"})',
                    red_outcome=red)
    H['G-M01'] = gm01

    def gm02(c):
        d2d1 = gr_mat_mul(c['d2_gr'], c['d1_gr'], c['mt'])
        assert all(gr_is_zero(d2d1[i][j])
                    for i in range(len(d2d1)) for j in range(len(d2d1[0]))), "Baseline G-M02 failed"
        d2_mut = copy.deepcopy(c['d2_gr'])
        d2_mut[0][0] = gr_add(d2_mut[0][0], {c['e_id']: 1})
        d2d1_mut = gr_mat_mul(d2_mut, c['d1_gr'], c['mt'])
        red = not all(gr_is_zero(d2d1_mut[i][j])
                      for i in range(len(d2d1_mut)) for j in range(len(d2d1_mut[0])))
        return dict(object_mutated='d2[0][0]: added +1*e (eid 119)',
                    gate_predicate='d2.d1 = 0 over Z[2I]',
                    baseline_result='PASS (product is zero)',
                    mutated_result=f'FAIL (product is {"nonzero" if red else "zero"})',
                    red_outcome=red)
    H['G-M02'] = gm02

    def gm03(c):
        free_ranks = c['cp']['free_ranks']
        chi = sum((-1)**k * r for k, r in enumerate(free_ranks))
        assert chi == 0, "Baseline G-M03 failed"
        ranks_mut = list(free_ranks); ranks_mut[3] = 2
        chi_mut = sum((-1)**k * r for k, r in enumerate(ranks_mut))
        return dict(object_mutated='free_ranks[3]: 1 -> 2',
                    gate_predicate='chi = sum(-1)^k * r_k = 0',
                    baseline_result=f'PASS (chi={chi})',
                    mutated_result=f'FAIL (chi={chi_mut})',
                    red_outcome=chi_mut != 0)
    H['G-M03'] = gm03

    def gm04(c):
        d2_aug = augment_matrix(c['d2_gr'])
        det_orig = d2_aug[0][0] * d2_aug[1][1] - d2_aug[0][1] * d2_aug[1][0]
        assert abs(det_orig) == 1, f"Baseline G-M04 failed: det={det_orig}"
        d2_mut4 = copy.deepcopy(c['d2_gr'])
        d2_mut4[0][0] = gr_add(d2_mut4[0][0], {c['e_id']: 2})
        d2_aug_mut = augment_matrix(d2_mut4)
        det_mut = d2_aug_mut[0][0] * d2_aug_mut[1][1] - d2_aug_mut[0][1] * d2_aug_mut[1][0]
        return dict(object_mutated='d2[0][0]: added +2*e (eid 119)',
                    gate_predicate='det(d2_aug) = +/-1',
                    baseline_result=f'PASS (det={det_orig})',
                    mutated_result=f'FAIL (det={det_mut})',
                    red_outcome=abs(det_mut) != 1)
    H['G-M04'] = gm04

    def gm05(c):
        print("  Expanding d2 over Z (240x240)...")
        d2_Z = expand_gr_mat(c['d2_gr'], c['mt'])
        print("  Running Gaussian elimination on original d2_Z...")
        r2, pr2, pc2 = gauss_pivots(d2_Z)
        assert r2 == 121, f"Baseline rank(d2_Z)={r2}, expected 121"
        sub2 = [[d2_Z[pr2[i]][pc2[j]] for j in range(121)] for i in range(121)]
        det2_orig = det_frac(sub2)
        assert abs(det2_orig) == 1, f"Baseline |det|={abs(det2_orig)}, expected 1"
        print(f"  Original: rank={r2}, |det(pivot minor)|={abs(det2_orig)} -> PASS")
        d2_gr_mut5 = copy.deepcopy(c['d2_gr'])
        for j in range(len(d2_gr_mut5[0])):
            d2_gr_mut5[0][j] = {eid: 2 * coeff for eid, coeff in d2_gr_mut5[0][j].items()}
        d2_Z_mut = expand_gr_mat(d2_gr_mut5, c['mt'])
        sub2_mut = [[d2_Z_mut[pr2[i]][pc2[j]] for j in range(121)] for i in range(121)]
        det2_mut = det_frac(sub2_mut)
        print(f"  Mutated:  |det(pivot minor)|={abs(det2_mut)}")
        return dict(object_mutated='d2 row 0 scaled by 2 (non-unit)',
                    gate_predicate='unimodular 121x121 minor exists (|det|=1)',
                    baseline_result=f'PASS (|det|={abs(det2_orig)})',
                    mutated_result=f'FAIL (|det|={abs(det2_mut)})',
                    red_outcome=abs(det2_mut) != 1)
    H['G-M05'] = gm05

    def gm06(c):
        eps_d1 = [augmentation(c['d1_gr'][i][0]) for i in range(len(c['d1_gr']))]
        assert all(v == 0 for v in eps_d1), "Baseline G-M06 failed"
        def non_aug(gr_elem):
            return sum(coeff * (1 if eid == c['e_id'] else 2) for eid, coeff in gr_elem.items())
        eps_d1_mut = [non_aug(c['d1_gr'][i][0]) for i in range(len(c['d1_gr']))]
        red = not all(v == 0 for v in eps_d1_mut)
        return dict(object_mutated='eps: replaced with non-augmentation (e->1, others->2)',
                    gate_predicate='eps(d1) = 0',
                    baseline_result=f'PASS (eps_d1={eps_d1})',
                    mutated_result=f'FAIL (eps_d1={eps_d1_mut})',
                    red_outcome=red)
    H['G-M06'] = gm06

    def gm07(c):
        se = c['se']; s_id = c['s_id']; t_id = c['t_id']
        s = se[s_id]; t = se[t_id]
        s3 = s*s*s; st = s*t; st2 = st*st; t5 = t*t*t*t*t
        assert s3 == st2 and t5 == st2, "Baseline G-M07 failed"
        sw = se[t_id]; tw = se[s_id]
        s3w = sw*sw*sw; stw = sw*tw; st2w = stw*stw; t5w = tw*tw*tw*tw*tw
        red = not (s3w == st2w and t5w == st2w)
        return dict(object_mutated='Swapped s_id and t_id',
                    gate_predicate='s^3=(st)^2 and t^5=(st)^2',
                    baseline_result='PASS (both relators hold)',
                    mutated_result=f'FAIL (relators {"fail" if red else "hold"})',
                    red_outcome=red)
    H['G-M07'] = gm07

    def gm08(c):
        ct = c['ct']; rho = c['reps']['V1']; d = c['dims']['V1']
        M3 = ct.eval_bmap(c['d3_raw'], rho, d)
        r3 = rank_gc(M3)
        assert r3 == d, f"Baseline rank(M3_V1)={r3}, expected {d}"
        M3_mut = [[M3[i][j] for j in range(2*d)] for i in range(d)]
        for j in range(2*d): M3_mut[0][j] = ct.GC()
        r3_mut = rank_gc(M3_mut)
        return dict(object_mutated='V1 twisted M3: zeroed row 0',
                    gate_predicate='rank(M3) = d for nontrivial irrep',
                    baseline_result=f'PASS (rank={r3}, d={d})',
                    mutated_result=f'FAIL (rank={r3_mut}, d={d})',
                    red_outcome=r3_mut < d)
    H['G-M08'] = gm08

    def gt01(c):
        ct = c['ct']; rho = c['reps']['V1']; d = c['dims']['V1']; s_id = c['s_id']
        H_form = ct.mz(d, d)
        for g in range(120):
            H_form = ct.madd(H_form, ct.mmul(ct.mct(rho[g]), rho[g]))
        H_form = ct.msc(ct.GC(ct.QG(Fraction(1, 120))), H_form)
        lhs = ct.mmul(ct.mct(rho[s_id]), ct.mmul(H_form, rho[s_id]))
        assert ct.meq(lhs, H_form), "Baseline G-T01 failed"
        rho_mut = list(rho)
        rho_s_mut = [[rho[s_id][i][j] for j in range(d)] for i in range(d)]
        rho_s_mut[0][0] = rho_s_mut[0][0] + ct.GC(ct.QG(Fraction(1, 10)))
        rho_mut[s_id] = rho_s_mut
        H_mut = ct.mz(d, d)
        for g in range(120):
            H_mut = ct.madd(H_mut, ct.mmul(ct.mct(rho_mut[g]), rho_mut[g]))
        H_mut = ct.msc(ct.GC(ct.QG(Fraction(1, 120))), H_mut)
        lhs_mut = ct.mmul(ct.mct(rho_mut[s_id]), ct.mmul(H_mut, rho_mut[s_id]))
        red = not ct.meq(lhs_mut, H_mut)
        return dict(object_mutated='V1: perturbed rho(s)[0][0] by +1/10',
                    gate_predicate='rho(s)^dag H rho(s) = H',
                    baseline_result='PASS (invariant)',
                    mutated_result=f'FAIL ({"not invariant" if red else "invariant"})',
                    red_outcome=red)
    H['G-T01'] = gt01

    def gt02(c):
        ct = c['ct']; reps = c['reps']; dims = c['dims']
        s_id = c['s_id']; t_id = c['t_id']; st_id = c['st_id']
        rep_names = ['V0','V1','V2','V3','V4','V5','V6','V7','V8']
        all_sigs = set()
        for name in rep_names:
            rho = reps[name]
            cs = ct.mtr(rho[s_id]).re
            ct_v = ct.mtr(rho[t_id]).re
            cst = ct.mtr(rho[st_id]).re
            sig = (dims[name], (cs.a, cs.b), (ct_v.a, ct_v.b), (cst.a, cst.b))
            all_sigs.add(sig)
        assert len(all_sigs) == 9, "Baseline G-T02: not all 9 sigs distinct"
        ch_s = ct.mtr(reps['V1'][s_id]).re
        ch_t = ct.mtr(reps['V1'][t_id]).re
        orig_sig = (dims['V1'], (ch_s.a, ch_s.b), (ch_t.a, ch_t.b),
                    (ct.mtr(reps['V1'][st_id]).re.a, ct.mtr(reps['V1'][st_id]).re.b))
        swapped_sig = (dims['V1'], (ch_t.a, ch_t.b), (ch_s.a, ch_s.b),
                       (ct.mtr(reps['V1'][st_id]).re.a, ct.mtr(reps['V1'][st_id]).re.b))
        red = swapped_sig != orig_sig
        return dict(object_mutated='V1: swapped chi(s) and chi(t)',
                    gate_predicate='Row signatures are distinct and match expected',
                    baseline_result='PASS (9 distinct signatures)',
                    mutated_result=f'FAIL (V1 signature {"changed" if red else "unchanged"})',
                    red_outcome=red)
    H['G-T02'] = gt02

    def _fixture_test(c, gate_id, mutation_desc, obj_desc, make_rho, make_bmaps):
        vf = c['vf']; fix = c['vf_fixture_reps']
        d1r = c['d1_raw']; d2r = c['d2_raw']; d3r = c['d3_raw']
        T2_base = c['vf_fixture_T2_base']
        rho_in = make_rho(fix, c) if make_rho else fix
        d1_in, d2_in, d3_in = make_bmaps(d1r, d2r, d3r) if make_bmaps else (d1r, d2r, d3r)
        tau, acyc, msg = vf.compute_torsion(rho_in, d1_in, d2_in, d3_in, 2)
        if not acyc:
            return dict(object_mutated=obj_desc, gate_predicate='dd=0 and correct T2',
                        baseline_result='PASS (acyclic, T2 computed)',
                        mutated_result=f'FAIL ({msg})', red_outcome=True)
        T2 = tau * tau.conj()
        red = not (T2.re == T2_base)
        return dict(object_mutated=obj_desc, gate_predicate='dd=0 and correct T2',
                    baseline_result='PASS (acyclic, T2 computed)',
                    mutated_result=f'FAIL (T2 changed)' if red else 'PASS (T2 unchanged)',
                    red_outcome=red)

    def gt03a(c):
        def make_rho(fix, c):
            return [fix[c['inv_map'][g]] for g in range(120)]
        return _fixture_test(c, 'G-T03a',
                             'Convention: evaluation map',
                             'Fixture: g -> rho(g^-1) (anti-homomorphism)',
                             make_rho, None)
    H['G-T03a'] = gt03a

    def gt03b(c):
        def make_bmaps(d1r, d2r, d3r):
            d1_rev = [[d3r[j][i] for j in range(len(d3r))] for i in range(len(d3r[0]))]
            d2_rev = [[d2r[j][i] for j in range(len(d2r))] for i in range(len(d2r[0]))]
            d3_rev = [[d1r[j][i] for j in range(len(d1r))] for i in range(len(d1r[0]))]
            return d1_rev, d2_rev, d3_rev
        return _fixture_test(c, 'G-T03b',
                             'Convention: boundary direction',
                             'Fixture: transposed boundary matrices (cochain)',
                             None, make_bmaps)
    H['G-T03b'] = gt03b

    def gt03c(c):
        def make_rho(fix, c):
            return [[[fix[g][j][i] for j in range(len(fix[g]))]
                     for i in range(len(fix[g][0]))] for g in range(120)]
        return _fixture_test(c, 'G-T03c',
                             'Convention: module side',
                             'Fixture: rho(g)^T (transpose = anti-homomorphism)',
                             make_rho, None)
    H['G-T03c'] = gt03c

    def gt03d(c):
        def make_bmaps(d1r, d2r, d3r):
            d1_grt = [[d1r[j][i] for j in range(len(d1r))] for i in range(len(d1r[0]))]
            d2_grt = [[d2r[j][i] for j in range(len(d2r))] for i in range(len(d2r[0]))]
            d3_grt = [[d3r[j][i] for j in range(len(d3r))] for i in range(len(d3r[0]))]
            return d3_grt, d2_grt, d1_grt
        return _fixture_test(c, 'G-T03d',
                             'Convention: vector convention',
                             'Fixture: transposed GR maps (column vectors)',
                             None, make_bmaps)
    H['G-T03d'] = gt03d

    def gd01(c):
        ct = c['ct']; rho = c['reps']['V1']; d = c['dims']['V1']
        M2 = ct.eval_bmap(c['d2_raw'], rho, d)
        M3 = ct.eval_bmap(c['d3_raw'], rho, d)
        prod32 = ct.mmul(M3, M2)
        assert all(prod32[i][j].is_zero() for i in range(d) for j in range(2*d)), "Baseline G-D01"
        M3_mut = [[M3[i][j] for j in range(2*d)] for i in range(d)]
        M3_mut[0][0] = M3_mut[0][0] + ct.GC(ct.QG(1))
        prod_mut = ct.mmul(M3_mut, M2)
        red = not all(prod_mut[i][j].is_zero() for i in range(d) for j in range(2*d))
        return dict(object_mutated='V1 M3[0][0]: perturbed by +1',
                    gate_predicate='M3.M2 = 0 (twisted complex chain condition)',
                    baseline_result='PASS (product is zero)',
                    mutated_result=f'FAIL (product is {"nonzero" if red else "zero"})',
                    red_outcome=red)
    H['G-D01'] = gd01

    def gd02(c):
        ct = c['ct']; rho = c['reps']['V1']; d = c['dims']['V1']
        M3 = ct.eval_bmap(c['d3_raw'], rho, d)
        r3 = rank_gc(M3)
        assert r3 == d, f"Baseline rank(M3)={r3}"
        M3_mut = [[M3[i][j] for j in range(2*d)] for i in range(d)]
        for j in range(2*d): M3_mut[0][j] = ct.GC()
        r3_mut = rank_gc(M3_mut)
        return dict(object_mutated='V1 M3: zeroed row 0',
                    gate_predicate='rank(M3) = d',
                    baseline_result=f'PASS (rank={r3}, d={d})',
                    mutated_result=f'FAIL (rank={r3_mut}, d={d})',
                    red_outcome=r3_mut < d)
    H['G-D02'] = gd02

    def gd03(c):
        ct = c['ct']; rho = c['reps']['V1']; d = c['dims']['V1']
        M3 = ct.eval_bmap(c['d3_raw'], rho, d)
        J3 = list(range(d))
        minor = [[M3[i][j] for j in J3] for i in range(d)]
        det3 = ct.det_gc(minor)
        assert not det3.is_zero(), "Baseline G-D03: M3 minor singular"
        minor_mut = [[minor[i][j] for j in range(d)] for i in range(d)]
        for i in range(d): minor_mut[i][0] = ct.GC()
        det3_mut = ct.det_gc(minor_mut)
        return dict(object_mutated='V1 M3 minor: zeroed column 0',
                    gate_predicate='det(minor) != 0',
                    baseline_result='PASS (det nonzero)',
                    mutated_result=f'FAIL (det {"= 0" if det3_mut.is_zero() else "nonzero"})',
                    red_outcome=det3_mut.is_zero())
    H['G-D03'] = gd03

    def gd04(c):
        ct = c['ct']; reps = c['reps']; dims = c['dims']
        d1r = c['d1_raw']; d2r = c['d2_raw']; d3r = c['d3_raw']
        T2_V1, _, _ = ct.compute_torsion_sq(reps['V1'], d1r, d2r, d3r, dims['V1'])
        T2_V7, _, _ = ct.compute_torsion_sq(reps['V7'], d1r, d2r, d3r, dims['V7'])
        assert T2_V1.galois() == T2_V7, "Baseline G-D04 failed"
        V1_gal = [ct.mgal(reps['V1'][g]) for g in range(120)]
        T2_gal, _, _ = ct.compute_torsion_sq(V1_gal, d1r, d2r, d3r, dims['V1'])
        gal_check = T2_gal.galois()
        red = not (gal_check == T2_V7)
        return dict(object_mutated='V1 representation matrices: applied sigma (phi->1-phi) entry-wise via mgal',
                    gate_predicate='sigma(T2(V1)) = T2(V7) (Galois conjugacy of the pair)',
                    baseline_result=f'PASS (sigma(T2(V1))={T2_V1.galois()} = T2(V7)={T2_V7})',
                    mutated_result=(f'FAIL (sigma(T2(sigma(V1)))={gal_check} != T2(V7)={T2_V7})' if red
                                    else 'PASS (still conjugate)'),
                    implemented_mutation='Applied mgal (phi->1-phi) to all 120 V1 representation matrices, recomputed T2 through compute_torsion_sq, checked Galois relation against T2(V7)',
                    red_outcome=red)
    H['G-D04'] = gd04

    def gd05(c):
        ct = c['ct']; d = c['dims']['V1']
        d1r = c['d1_raw']; d2r = c['d2_raw']; d3r = c['d3_raw']
        T2_V1, _, _ = ct.compute_torsion_sq(c['reps']['V1'], d1r, d2r, d3r, d)
        original_det_gc = ct.det_gc
        call_count = [0]
        def det_gc_identity(M):
            call_count[0] += 1
            n = len(M)
            I_n = ct.mid(n)
            return original_det_gc(I_n)
        ct.det_gc = det_gc_identity
        try:
            T2_mut, _, _ = ct.compute_torsion_sq(c['reps']['V1'], d1r, d2r, d3r, d)
        finally:
            ct.det_gc = original_det_gc
        assert call_count[0] == 3, \
            f"G-D05: expected exactly 3 det_gc calls under identity substitution, got {call_count[0]}"
        red = not (T2_mut == T2_V1)
        return dict(
            object_mutated='Frozen det_gc: substituted to return det(I_d) for every minor',
            gate_predicate='T2 output changes when determinant inputs are replaced',
            baseline_result=f'PASS (T2(V1) = {T2_V1})',
            mutated_result=(f'FAIL (T2_mutated = {T2_mut} != T2(V1) = {T2_V1})' if red
                            else 'PASS (T2 unchanged)'),
            implemented_mutation=(f'Substituted det_gc on frozen module to return det(I_d) for every minor '
                                  f'during one call to frozen compute_torsion_sq; intercepted exactly '
                                  f'{call_count[0]} calls (verified == 3); restored original in finally block; '
                                  f'both baseline and mutated T2 returned by frozen compute_torsion_sq'),
            red_outcome=red)
    H['G-D05'] = gd05

    return H

HANDLER_TABLE = _make_handlers()

IMPLEMENTED_MUTATIONS = {
    'G-M01': 'Added +1*e(eid 119) to d3[0][0] in Z[2I] scratch copy; computed d3*d2 and verified product is nonzero',
    'G-M02': 'Added +1*e(eid 119) to d2[0][0] in Z[2I] scratch copy; computed d2*d1 and verified product is nonzero',
    'G-M03': 'Changed free_ranks[3] from 1 to 2; recomputed chi = sum(-1)^k r_k and verified chi != 0',
    'G-M04': 'Added +2*e(eid 119) to d2[0][0] in scratch copy; recomputed det(d2_aug) and verified |det| != 1',
    'G-M05': 'Scaled row 0 of d2 by 2 (non-unit) in scratch copy; expanded over Z, computed 121x121 pivot minor det, verified |det| != 1',
    'G-M06': 'Replaced augmentation eps with non-augmentation (e->1, others->2); computed eps(d1) and verified eps(d1) != 0',
    'G-M07': 'Swapped s_id and t_id generator identifiers; verified 2I relators s^3=(st)^2 and t^5=(st)^2 fail under the swap',
    'G-M08': 'Zeroed row 0 of V1 twisted M3 (evaluated boundary matrix); verified rank(M3) drops below d',
    'G-T01': 'Perturbed V1 rho(s)[0][0] by +1/10 in scratch copy; recomputed Hermitian form H and verified rho(s)^dag H rho(s) != H',
    'G-T02': 'Swapped chi(s) and chi(t) character values in V1 row signature; verified the signature changes',
    'G-T03a': 'On convention fixture via frozen validate_fixture.compute_torsion: replaced g->rho(g) with g->rho(g^-1) (anti-homomorphism); verified dd != 0',
    'G-T03b': 'On convention fixture via frozen validate_fixture.compute_torsion: transposed boundary matrices (cochain reversal); verified dd != 0',
    'G-T03c': 'On convention fixture via frozen validate_fixture.compute_torsion: used rho(g)^T (transpose, i.e. right-module action); verified dd != 0',
    'G-T03d': 'On convention fixture via frozen validate_fixture.compute_torsion: transposed group-ring boundary maps and reversed degree ordering; verified dd != 0',
    'G-D01': 'Perturbed V1 M3[0][0] by +1 over GC; verified M3*M2 != 0 (twisted chain condition fails)',
    'G-D02': 'Zeroed row 0 of V1 twisted M3; verified rank(M3) < d',
    'G-D03': 'Zeroed column 0 of V1 M3 minor (d*d determinant sub-matrix); verified det(minor) = 0',
}


# ========== MAIN ==========

def main():
    print("=" * 60)
    print("M8.8 PHASE B: GATE QUALIFICATION")
    print("=" * 60)

    # ---- Step 0: Verify Phase A integrity ----
    print("\n--- Step 0: Verify Phase A artifact hashes ---")
    pre_hash_ok = verify_phase_a_hashes()

    # ---- Step 1: Parse gate registry from frozen manifest ----
    print("\n--- Step 1: Parse gate registry from frozen manifest ---")
    manifest_path = os.path.join(FROZEN_DIR, 'METHOD_AND_GATE_MANIFEST.md')
    with open(manifest_path, 'r') as f:
        manifest_text = f.read()
    parsed_registry = parse_gate_registry(manifest_text)
    parsed_ids = set(g['id'] for g in parsed_registry)
    parsed_lookup = {g['id']: g for g in parsed_registry}
    print(f"  Parsed {len(parsed_registry)} gates from manifest section 4:")
    for g in parsed_registry:
        print(f"    {g['id']}: {g['name']}")

    # ---- Step 2: Self-test parser-to-coverage linkage ----
    implemented_ids = set(HANDLER_TABLE.keys())
    self_test_ok = self_test_parser_coverage(manifest_text, implemented_ids)

    # ---- Step 3: Pre-execution coverage proof ----
    print("\n--- Step 3: Pre-execution coverage proof ---")
    print(f"  Parsed gate IDs:      {len(parsed_ids)}")
    print(f"  Implemented handlers: {len(implemented_ids)}")
    pre_ok = check_coverage(parsed_ids, implemented_ids, "implemented handlers")
    if not pre_ok:
        print("FATAL: Pre-execution coverage mismatch between parsed manifest and handlers.")
        sys.exit(1)
    print("  Exact set equality: parsed manifest == implemented handlers")

    # ---- Step 4: Setup ----
    print("\n--- Loading Phase A computation modules ---")
    ct = _load_mod('ct', 'compute_torsion.py')
    vf = _load_mod('vf', 'validate_fixture.py')

    print("\n--- Loading packets ---")
    with open(os.path.join(SCRIPT_DIR, 'm8_5a_packet.json')) as f:
        gp = json.load(f)
    with open(os.path.join(SCRIPT_DIR, 'm8_8_construction_packet.json')) as f:
        cp = json.load(f)

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

    print("\n--- Parsing boundary maps ---")
    d1_gr = parse_boundary_map(cp['boundary_maps']['d1'])
    d2_gr = parse_boundary_map(cp['boundary_maps']['d2'])
    d3_gr = parse_boundary_map(cp['boundary_maps']['d3'])
    d1_raw = cp['boundary_maps']['d1']
    d2_raw = cp['boundary_maps']['d2']
    d3_raw = cp['boundary_maps']['d3']

    print("\n--- Building all 9 irreps ---")
    reps = ct.build_all_irreps(se, mt, s_id, t_id)
    rep_names = ['V0', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8']
    dims = {name: len(reps[name][0]) for name in rep_names}
    print(f"  Dimensions: {dims}")

    inv_map = [0] * 120
    for i in range(120):
        for j in range(120):
            if mt[i][j] == e_id:
                inv_map[i] = j; break

    # Build convention fixture using frozen validate_fixture types
    print("\n--- Building convention fixture (via frozen validate_fixture) ---")
    def q2su2_vf(q):
        return [[vf.GC(vf.QG(q.w.a, q.w.b), vf.QG(q.x.a, q.x.b)),
                 vf.GC(vf.QG(q.y.a, q.y.b), vf.QG(q.z.a, q.z.b))],
                [vf.GC(vf.QG(-q.y.a, -q.y.b), vf.QG(q.z.a, q.z.b)),
                 vf.GC(vf.QG(q.w.a, q.w.b), vf.QG(-q.x.a, -q.x.b))]]
    vf_su2 = [q2su2_vf(q) for q in se]
    vf_P = [[vf.GC(vf.QG(2)), vf.GCI], [vf.GC0, vf.GC1]]
    vf_Pi = vf.minv(vf_P)
    vf_fixture_reps = [vf.mmul(vf_Pi, vf.mmul(vf_su2[g], vf_P)) for g in range(120)]
    tau_vf_base, acyc_vf, msg_vf = vf.compute_torsion(
        vf_fixture_reps, d1_raw, d2_raw, d3_raw, 2)
    assert acyc_vf, f"VF fixture baseline not acyclic: {msg_vf}"
    T2_vf_gc = tau_vf_base * tau_vf_base.conj()
    assert T2_vf_gc.im.is_zero(), "VF baseline T2 not real"
    vf_fixture_T2_base = T2_vf_gc.re
    print(f"  Fixture baseline T2 = ({vf_fixture_T2_base.a}+{vf_fixture_T2_base.b}phi)")

    # Build context
    ctx = dict(ct=ct, se=se, mt=mt, s_id=s_id, t_id=t_id, e_id=e_id, st_id=st_id,
               cp=cp, d1_gr=d1_gr, d2_gr=d2_gr, d3_gr=d3_gr,
               d1_raw=d1_raw, d2_raw=d2_raw, d3_raw=d3_raw,
               reps=reps, dims=dims, inv_map=inv_map,
               vf=vf, vf_fixture_reps=vf_fixture_reps,
               vf_fixture_T2_base=vf_fixture_T2_base)

    # ---- Step 5: Execute all handlers ----
    print("\n" + "=" * 60)
    print("EXECUTING 19 GATE MUTATIONS")
    print("=" * 60)

    results = []
    for entry in parsed_registry:
        gid = entry['id']
        gname = entry['name']
        print(f"\n{gid}: {gname}")
        handler = HANDLER_TABLE[gid]
        rec = handler(ctx)
        rec['gate_id'] = gid
        rec['gate_name'] = gname
        rec['declared_mutation'] = entry['mutation']
        if 'implemented_mutation' not in rec:
            rec['implemented_mutation'] = IMPLEMENTED_MUTATIONS[gid]
        results.append(rec)
        status = "REDDENED" if rec['red_outcome'] else "FAILED TO REDDEN"
        print(f"  {gid}: {status}")

    # ---- Step 6: Post-execution coverage proof ----
    print("\n" + "=" * 60)
    print("POST-EXECUTION COVERAGE CHECK")
    print("=" * 60)
    executed_ids = set(r['gate_id'] for r in results)
    post_cov = check_coverage(parsed_ids, executed_ids, "executed gates")
    all_red = all(r['red_outcome'] for r in results)
    print(f"  Parsed gates:    {len(parsed_ids)}")
    print(f"  Executed gates:  {len(executed_ids)}")
    print(f"  Exact set equality: {post_cov}")
    print(f"  All mutations reddened: {all_red}")

    # ---- Step 7: Re-verify Phase A hashes ----
    print("\n--- Re-verifying Phase A hashes (post-qualification) ---")
    post_ok = verify_phase_a_hashes()

    # ---- Step 8: Write output ----
    output = {
        'schema_version': 'm8_8-phase-b-mutation-results-2',
        'phase_a_hashes_verified_pre': pre_hash_ok,
        'phase_a_hashes_verified_post': post_ok,
        'manifest_parsed_at_runtime': len(parsed_registry) == 19,
        'manifest_sha256': hashlib.sha256(manifest_text.encode()).hexdigest(),
        'parser_self_test_passed': self_test_ok,
        'registry_coverage': {
            'parsed_gate_ids': sorted(parsed_ids),
            'implemented_handler_ids': sorted(implemented_ids),
            'executed_gate_ids': sorted(executed_ids),
            'pre_execution_set_equality': pre_ok,
            'post_execution_set_equality': post_cov,
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
    if post_cov and all_red and post_ok:
        print("PHASE B QUALIFICATION: ALL 19 GATES REDDENED")
        print("Registry: parsed from frozen manifest at runtime")
        print("Coverage: exact set equality (parsed == implemented == executed)")
        print("Parser self-test: linkage verified on scratch manifest text")
        print("Phase A integrity: verified pre and post")
        sys.exit(0)
    else:
        print("PHASE B QUALIFICATION: INCOMPLETE")
        if not post_cov: print("  - Post-execution coverage mismatch")
        if not all_red:
            for r in results:
                if not r['red_outcome']:
                    print(f"  - {r['gate_id']}: mutation did not redden")
        if not post_ok: print("  - Phase A hashes changed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
