#!/usr/bin/env python3
"""
Pre-implementation validation script for M8.8 clean-room reproduction.

This script validates:
1. Canonical enumeration of the 120 elements of 2I (SHA-256 check)
2. Abstract generator correspondence: s (ID 118), t (ID 80) satisfy relators
3. d_n d_{n-1} = 0 over Z[2I]
4. Free ranks and Euler characteristic
5. Augmentation: eps(d1) = 0
6. Group packet SHA-256 matches construction packet declaration

This is a VALIDATION ARTIFACT, not production implementation.
It does NOT compute torsion values or populate RAW_OUTPUT.json.
"""

import json
import hashlib
from fractions import Fraction
import sys

# ============================================================
# Q(phi) arithmetic: elements are (a, b) meaning a + b*phi
# where phi = (1+sqrt(5))/2, phi^2 = phi + 1
# a, b are Fraction (exact rationals)
# ============================================================

def qphi(a, b=0):
    return (Fraction(a), Fraction(b))

def qphi_add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def qphi_sub(x, y):
    return (x[0] - y[0], x[1] - y[1])

def qphi_mul(x, y):
    a, b = x
    c, d = y
    return (a*c + b*d, a*d + b*c + b*d)

def qphi_neg(x):
    return (-x[0], -x[1])

def qphi_eq(x, y):
    return x[0] == y[0] and x[1] == y[1]

ZERO = qphi(0)
ONE = qphi(1)

# ============================================================
# Quaternion arithmetic over Q(phi)
# q = (q0, q1, q2, q3) where each qi is a Q(phi) element
# representing q0 + q1*i + q2*j + q3*k
# ============================================================

def quat_mul(p, q):
    a0, a1, a2, a3 = p
    b0, b1, b2, b3 = q
    c0 = qphi_sub(qphi_sub(qphi_sub(qphi_mul(a0, b0), qphi_mul(a1, b1)), qphi_mul(a2, b2)), qphi_mul(a3, b3))
    c1 = qphi_add(qphi_add(qphi_sub(qphi_mul(a0, b1), qphi_mul(a3, b2)), qphi_mul(a1, b0)), qphi_mul(a2, b3))
    c2 = qphi_add(qphi_add(qphi_sub(qphi_mul(a0, b2), qphi_mul(a1, b3)), qphi_mul(a2, b0)), qphi_mul(a3, b1))
    c3 = qphi_add(qphi_add(qphi_sub(qphi_mul(a0, b3), qphi_mul(a2, b1)), qphi_mul(a1, b2)), qphi_mul(a3, b0))
    return (c0, c1, c2, c3)

def quat_inv(q):
    """Inverse of a unit quaternion: conjugate (since |q|=1)."""
    q0, q1, q2, q3 = q
    return (q0, qphi_neg(q1), qphi_neg(q2), qphi_neg(q3))

def quat_eq(p, q):
    return all(qphi_eq(p[i], q[i]) for i in range(4))

def quat_identity():
    return (ONE, ZERO, ZERO, ZERO)

def quat_neg(q):
    return tuple(qphi_neg(c) for c in q)

# ============================================================
# Parse group packet
# ============================================================

def parse_component(s):
    """Parse '(a + b*phi)/2' into Q(phi) value = (a/2) + (b/2)*phi."""
    s = s.strip()
    # Format: "(A + B*phi)/2"
    import re
    m = re.match(r'\((-?\d+)\s*\+\s*(-?\d+)\*phi\)/2', s)
    if not m:
        raise ValueError(f"Cannot parse: {s}")
    A, B = int(m.group(1)), int(m.group(2))
    return qphi(Fraction(A, 2), Fraction(B, 2))

def parse_generator(gen_list):
    """Parse a list of 4 component strings into a quaternion."""
    return tuple(parse_component(s) for s in gen_list)

# ============================================================
# Canonical ID computation
# ============================================================

def to_8int_key(q):
    """Convert quaternion to the 8-integer sort key per protocol §4.2.
    Components in quaternion_basis order [1, i, j, k].
    Each component is (A + B*phi)/2, key is (A, B) per component."""
    key = []
    for comp in q:
        a, b = comp
        A = int(a * 2)
        B = int(b * 2)
        if Fraction(A, 2) != a or Fraction(B, 2) != b:
            raise ValueError(f"Component {comp} does not have integer (A,B)/2 form")
        key.extend([A, B])
    return tuple(key)

# ============================================================
# Main validation
# ============================================================

def main():
    results = {}

    # Load packets
    with open("m8_5a_packet.json") as f:
        group_pkt_bytes = f.read()
    group_pkt = json.loads(group_pkt_bytes)

    with open("m8_8_construction_packet.json") as f:
        constr_pkt_bytes = f.read()
    constr_pkt = json.loads(constr_pkt_bytes)

    # ---- Check 1: Group packet SHA-256 ----
    gp_canonical = json.dumps(group_pkt, sort_keys=True, indent=2,
                               separators=(',', ': '), ensure_ascii=True) + '\n'
    gp_hash = hashlib.sha256(gp_canonical.encode('ascii')).hexdigest()
    print(f"Group packet SHA-256: {gp_hash}")
    expected_gp_hash = "e3b0c945bbbb15b4549fa641234c9461062c2337b3d1e372af621b614d4883a9"
    gp_hash_ok = (gp_hash == expected_gp_hash)
    print(f"  Matches protocol pin: {gp_hash_ok}")
    results['group_packet_hash'] = gp_hash_ok

    # Also check construction packet declares the right group packet hash
    cp_gp_hash = constr_pkt["group_packet_sha256"]
    print(f"  Construction packet declares: {cp_gp_hash}")
    print(f"  Match: {cp_gp_hash == expected_gp_hash}")
    results['cp_gp_hash_match'] = (cp_gp_hash == expected_gp_hash)

    # ---- Check 2: Construction packet SHA-256 ----
    cp_canonical = json.dumps(constr_pkt, sort_keys=True, indent=2,
                               separators=(',', ': '), ensure_ascii=True) + '\n'
    cp_hash = hashlib.sha256(cp_canonical.encode('ascii')).hexdigest()
    print(f"\nConstruction packet SHA-256: {cp_hash}")
    expected_cp_hash = "df00c0222f98c481eb56b882cd867a6c3a4f8604b8633e81dec0cce1f8460a06"
    cp_hash_ok = (cp_hash == expected_cp_hash)
    print(f"  Matches protocol pin: {cp_hash_ok}")
    results['construction_packet_hash'] = cp_hash_ok

    # ---- Generate the group ----
    print("\n=== Generating 2I from group packet generators ===")
    g0 = parse_generator(group_pkt["generators"][0])
    g1 = parse_generator(group_pkt["generators"][1])
    print(f"Generator 0: scalar part = {g0[0]}")
    print(f"Generator 1: scalar part = {g1[0]}")

    # Generate by BFS
    identity = quat_identity()
    elements = {to_8int_key(identity): identity}
    frontier = [identity, g0, g1, quat_inv(g0), quat_inv(g1)]
    # Add initial elements
    for e in frontier:
        k = to_8int_key(e)
        elements[k] = e

    changed = True
    while changed:
        changed = False
        keys = list(elements.keys())
        for k1 in keys:
            e1 = elements[k1]
            for gen in [g0, g1, quat_inv(g0), quat_inv(g1)]:
                prod = quat_mul(e1, gen)
                pk = to_8int_key(prod)
                if pk not in elements:
                    elements[pk] = prod
                    changed = True
                prod2 = quat_mul(gen, e1)
                pk2 = to_8int_key(prod2)
                if pk2 not in elements:
                    elements[pk2] = prod2
                    changed = True

    print(f"Generated {len(elements)} elements")
    results['group_order'] = len(elements)
    if len(elements) != 120:
        print("ERROR: Expected 120 elements!")
        results['group_order_ok'] = False
    else:
        print("  Order 120: OK")
        results['group_order_ok'] = True

    # Sort by canonical key and assign IDs
    sorted_keys = sorted(elements.keys())
    id_to_quat = {}
    key_to_id = {}
    for rank, key in enumerate(sorted_keys):
        id_to_quat[rank] = elements[key]
        key_to_id[key] = rank

    # ---- Check 3: Canonical enumeration SHA-256 ----
    print("\n=== Verifying canonical enumeration hash ===")
    enum_array = [list(k) for k in sorted_keys]
    enum_json = json.dumps(enum_array, separators=(',', ':'))
    enum_hash = hashlib.sha256(enum_json.encode('ascii')).hexdigest()
    enum_bytes = len(enum_json.encode('ascii'))
    print(f"Enumeration SHA-256: {enum_hash}")
    print(f"Enumeration bytes: {enum_bytes}")
    expected_enum_hash = "27ff780d28d5d854d464ead87e8fc20244fac8334bda9f0600c6ee1b3c89561e"
    enum_hash_ok = (enum_hash == expected_enum_hash)
    print(f"  Matches protocol: {enum_hash_ok}")
    results['enum_hash'] = enum_hash_ok

    # Check specific ranks from protocol
    print(f"  Rank 0: {list(sorted_keys[0])}")
    print(f"  Rank 118: {list(sorted_keys[118])}")
    print(f"  Rank 119: {list(sorted_keys[119])}")
    results['rank0_check'] = (list(sorted_keys[0]) == [-2,0,0,0,0,0,0,0])
    results['rank119_check'] = (list(sorted_keys[119]) == [2,0,0,0,0,0,0,0])
    results['rank118_check'] = (list(sorted_keys[118]) == [1,0,1,0,1,0,1,0])
    print(f"  Rank 0 = [-2,0,0,0,0,0,0,0]: {results['rank0_check']}")
    print(f"  Rank 118 = [1,0,1,0,1,0,1,0]: {results['rank118_check']}")
    print(f"  Rank 119 = [2,0,0,0,0,0,0,0]: {results['rank119_check']}")

    # Helper: get element ID from quaternion
    def get_id(q):
        return key_to_id[to_8int_key(q)]

    # Helper: multiply elements by ID
    def mul_id(a_id, b_id):
        prod = quat_mul(id_to_quat[a_id], id_to_quat[b_id])
        return get_id(prod)

    # Helper: inverse by ID
    def inv_id(a_id):
        return get_id(quat_inv(id_to_quat[a_id]))

    # ---- Check 4: Identity and negation ----
    identity_id = get_id(identity)
    neg_one_id = get_id(quat_neg(identity))
    print(f"\n=== Group structure checks ===")
    print(f"Identity element ID: {identity_id}")
    print(f"  Expected 119: {identity_id == 119}")
    print(f"-1 element ID: {neg_one_id}")
    print(f"  Expected 0: {neg_one_id == 0}")
    results['identity_at_119'] = (identity_id == 119)
    results['neg_one_at_0'] = (neg_one_id == 0)

    # ---- Check 5: Abstract generators satisfy relators ----
    print("\n=== Checking relators: s^3 = t^5 = (st)^2 ===")
    s_id = constr_pkt["abstract_generators"]["s"]
    t_id = constr_pkt["abstract_generators"]["t"]
    print(f"s = element {s_id}, t = element {t_id}")

    s = id_to_quat[s_id]
    t = id_to_quat[t_id]

    # s^3
    s2 = quat_mul(s, s)
    s3 = quat_mul(s2, s)
    s3_id = get_id(s3)

    # t^5
    t2 = quat_mul(t, t)
    t3 = quat_mul(t2, t)
    t4 = quat_mul(t3, t)
    t5 = quat_mul(t4, t)
    t5_id = get_id(t5)

    # (st)^2
    st = quat_mul(s, t)
    st2 = quat_mul(st, st)
    st2_id = get_id(st2)

    # For the balanced presentation <s,t | s^3 = t^5 = (st)^2>,
    # these should all be equal (= the central element -1, since
    # in 2I the standard balanced presentation gives s^3 = t^5 = (st)^2 = -1)
    print(f"s^3 = element {s3_id}")
    print(f"t^5 = element {t5_id}")
    print(f"(st)^2 = element {st2_id}")
    relators_ok = (s3_id == t5_id == st2_id)
    print(f"All equal: {relators_ok}")
    results['relators_equal'] = relators_ok

    # Check that s^3 is the central element -1 (ID 0)
    print(f"s^3 = -1 (ID 0): {s3_id == 0}")
    results['s3_is_neg1'] = (s3_id == 0)

    # Check orders
    s6 = quat_mul(s3, s3)
    s6_id = get_id(s6)
    print(f"s^6 = element {s6_id} (should be 119=identity): {s6_id == 119}")
    results['s_order_divides_6'] = (s6_id == 119)

    t10 = quat_mul(t5, t5)
    t10_id = get_id(t10)
    print(f"t^10 = element {t10_id} (should be 119=identity): {t10_id == 119}")
    results['t_order_divides_10'] = (t10_id == 119)

    st_id = get_id(st)
    st4 = quat_mul(st2, st2)
    st4_id = get_id(st4)
    print(f"(st)^4 = element {st4_id} (should be 119=identity): {st4_id == 119}")
    results['st_order_divides_4'] = (st4_id == 119)

    # Verify s and t generate the full group
    generated = {identity_id}
    frontier_ids = [s_id, t_id, inv_id(s_id), inv_id(t_id)]
    for fid in frontier_ids:
        generated.add(fid)
    changed = True
    while changed:
        changed = False
        current = list(generated)
        for a in current:
            for b in [s_id, t_id, inv_id(s_id), inv_id(t_id)]:
                p = mul_id(a, b)
                if p not in generated:
                    generated.add(p)
                    changed = True
    print(f"\n<s,t> generates {len(generated)} elements")
    results['st_generates'] = (len(generated) == 120)
    print(f"  Generates full group: {results['st_generates']}")

    # ---- Check 6: Boundary maps d*d = 0 ----
    print("\n=== Checking d_n d_{n-1} = 0 over Z[2I] ===")

    # Parse boundary maps
    # Group ring element: list of [coefficient, element_id] pairs
    def parse_gr_element(terms):
        """Parse group ring element from construction packet format."""
        return [(c, eid) for c, eid in terms]

    def gr_mul(gr_a, gr_b):
        """Multiply two group ring elements."""
        result = {}
        for ca, ga in gr_a:
            for cb, gb in gr_b:
                prod_id = mul_id(ga, gb)
                coeff = ca * cb
                result[prod_id] = result.get(prod_id, 0) + coeff
        return [(c, g) for g, c in result.items() if c != 0]

    def gr_add(gr_a, gr_b):
        """Add two group ring elements."""
        result = {}
        for c, g in gr_a:
            result[g] = result.get(g, 0) + c
        for c, g in gr_b:
            result[g] = result.get(g, 0) + c
        return [(c, g) for g, c in result.items() if c != 0]

    def gr_is_zero(gr):
        return all(c == 0 for c, _ in gr)

    # Parse boundary maps from construction packet
    d1_raw = constr_pkt["boundary_maps"]["d1"]  # 2x1 matrix
    d2_raw = constr_pkt["boundary_maps"]["d2"]  # 2x2 matrix
    d3_raw = constr_pkt["boundary_maps"]["d3"]  # 1x2 matrix

    def parse_matrix(raw):
        """Parse a matrix of group ring elements."""
        rows = len(raw)
        cols = len(raw[0])
        mat = [[parse_gr_element(raw[i][j]) for j in range(cols)] for i in range(rows)]
        return mat, rows, cols

    d1, d1_rows, d1_cols = parse_matrix(d1_raw)
    d2, d2_rows, d2_cols = parse_matrix(d2_raw)
    d3, d3_rows, d3_cols = parse_matrix(d3_raw)

    print(f"d1: {d1_rows}x{d1_cols} matrix over Z[2I]")
    print(f"d2: {d2_rows}x{d2_cols} matrix over Z[2I]")
    print(f"d3: {d3_rows}x{d3_cols} matrix over Z[2I]")

    # Check dimensions
    print(f"\nDimension checks (free ranks [1,2,2,1]):")
    print(f"  d1 is {d1_rows}x{d1_cols}, expected 2x1: {d1_rows==2 and d1_cols==1}")
    print(f"  d2 is {d2_rows}x{d2_cols}, expected 2x2: {d2_rows==2 and d2_cols==2}")
    print(f"  d3 is {d3_rows}x{d3_cols}, expected 1x2: {d3_rows==1 and d3_cols==2}")
    results['d1_dims'] = (d1_rows==2 and d1_cols==1)
    results['d2_dims'] = (d2_rows==2 and d2_cols==2)
    results['d3_dims'] = (d3_rows==1 and d3_cols==2)

    # Matrix multiply over Z[2I]: C = A * B where A is rxm, B is mxn
    def mat_mul_gr(A, A_rows, A_cols, B, B_rows, B_cols):
        assert A_cols == B_rows
        C = [[[] for _ in range(B_cols)] for _ in range(A_rows)]
        for i in range(A_rows):
            for j in range(B_cols):
                entry = []
                for k in range(A_cols):
                    prod = gr_mul(A[i][k], B[k][j])
                    entry = gr_add(entry, prod)
                C[i][j] = entry
        return C, A_rows, B_cols

    # Check d2 * d1 = 0 (2x2 * 2x1 = 2x1, should be zero matrix)
    d2d1, _, _ = mat_mul_gr(d2, d2_rows, d2_cols, d1, d1_rows, d1_cols)
    d2d1_zero = all(gr_is_zero(d2d1[i][j]) for i in range(2) for j in range(1))
    print(f"\nd2 * d1 = 0: {d2d1_zero}")
    results['d2d1_zero'] = d2d1_zero

    # Check d3 * d2 = 0 (1x2 * 2x2 = 1x2, should be zero matrix)
    d3d2, _, _ = mat_mul_gr(d3, d3_rows, d3_cols, d2, d2_rows, d2_cols)
    d3d2_zero = all(gr_is_zero(d3d2[i][j]) for i in range(1) for j in range(2))
    print(f"d3 * d2 = 0: {d3d2_zero}")
    results['d3d2_zero'] = d3d2_zero

    # ---- Check 7: Augmentation ----
    print("\n=== Checking augmentation ===")
    # eps sends every group element to 1
    # eps(d1) should be 0, i.e., for each row of d1, sum of all coefficients = 0
    def augment(gr_elem):
        """Apply augmentation: sum of all coefficients."""
        return sum(c for c, _ in gr_elem)

    eps_d1 = [[augment(d1[i][j]) for j in range(d1_cols)] for i in range(d1_rows)]
    eps_d1_zero = all(eps_d1[i][j] == 0 for i in range(d1_rows) for j in range(d1_cols))
    print(f"eps(d1) = 0: {eps_d1_zero}")
    for i in range(d1_rows):
        for j in range(d1_cols):
            print(f"  eps(d1[{i},{j}]) = {eps_d1[i][j]}")
    results['eps_d1_zero'] = eps_d1_zero

    # ---- Check 8: Free ranks and Euler characteristic ----
    print("\n=== Free ranks and Euler characteristic ===")
    ranks = constr_pkt["free_ranks"]
    print(f"Free ranks: {ranks}")
    chi = sum((-1)**i * ranks[i] for i in range(len(ranks)))
    print(f"Euler characteristic: {chi}")
    results['euler_char_zero'] = (chi == 0)
    print(f"  chi = 0: {results['euler_char_zero']}")

    # ---- Check 9: d1 matches generator correspondence ----
    print("\n=== Checking d1 generator correspondence ===")
    # d1 is 2x1. Row 0 corresponds to e_s, row 1 to e_t.
    # d1[0,0] should be s - 1 (i.e., element s minus identity)
    # d1[1,0] should be t - 1

    d1_00 = d1[0][0]
    d1_10 = d1[1][0]
    print(f"d1[0,0] terms: {d1_00}")
    print(f"d1[1,0] terms: {d1_10}")

    # Expected: d1[0,0] = s - e = [(1, 118), (-1, 119)]
    expected_d1_00 = {(1, s_id), (-1, identity_id)}
    actual_d1_00 = {(c, g) for c, g in d1_00}
    d1_00_ok = (actual_d1_00 == expected_d1_00)
    print(f"d1[0,0] = s - 1: {d1_00_ok}")

    # Expected: d1[1,0] = t - e = [(1, 80), (-1, 119)]
    expected_d1_10 = {(1, t_id), (-1, identity_id)}
    actual_d1_10 = {(c, g) for c, g in d1_10}
    d1_10_ok = (actual_d1_10 == expected_d1_10)
    print(f"d1[1,0] = t - 1: {d1_10_ok}")
    results['d1_gen_correspondence'] = (d1_00_ok and d1_10_ok)

    # ---- Summary ----
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    all_pass = True
    for k, v in results.items():
        status = "PASS" if v else "FAIL"
        if not v:
            all_pass = False
        print(f"  {k}: {status}")
    print(f"\nOverall: {'ALL PASS' if all_pass else 'SOME FAILURES'}")
    return 0 if all_pass else 1

if __name__ == "__main__":
    sys.exit(main())
