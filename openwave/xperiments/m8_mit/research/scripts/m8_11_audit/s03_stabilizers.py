"""Item 3 (and data for items 4, 9): stabilizers of the lines [u], characters, character-subspace dimensions.
Action on sections: (r.Phi)(g) = Phi(r^-1 g)  <=>  fibre u -> conj(D(r)) u (all levels)."""
import json, pickle
import numpy as np
from math import factorial
from numD import D_of, su2_axis_angle
from checks import check, save_log

np.set_printoptions(precision=4, suppress=True)
s = np.sqrt
U = {
    "U1": {3: 1.0}, "U2": {0: 1.0}, "U3": {2: 1.0, -2: 1.0}, "U4": {3: 1.0, -3: 1.0},
    "U5": {2: s(13 / 25), -3: s(12 / 25)}, "U6": {3: 1.0, 0: s(23 / 10), -3: 1.0},
}


def uvec(dct):
    v = np.zeros(7, complex)
    for m, c in dct.items():
        v[m + 3] = c
    return v / np.linalg.norm(v)


def act(J, g, x):
    return np.conj(D_of(J, g)) @ x


def roots_points(u):
    # p(x,y) = sum u_m x^{3+m} y^{3-m} / sqrt(w_m), t = x/y
    w = [factorial(3 + m) * factorial(3 - m) / factorial(6) for m in range(-3, 4)]
    c = [u[k] / np.sqrt(w[k]) for k in range(7)]   # coefficient of t^{k}
    coeffs = np.array(c[::-1])                     # highest power first
    nz = np.nonzero(np.abs(coeffs) > 1e-12)[0]
    lead = nz[0]
    rts = np.roots(coeffs[lead:]) if len(coeffs[lead:]) > 1 else np.array([])
    pts = []
    for t in rts:
        a = abs(t) ** 2
        pts.append(np.array([2 * t.real, 2 * t.imag, a - 1]) / (a + 1))
    pts += [np.array([0, 0, 1.0])] * lead           # roots at infinity
    return pts


def multiset(pts, tol=1e-6):
    groups = []
    for p in pts:
        for gp in groups:
            if np.linalg.norm(gp[0] - p) < tol:
                gp[1] += 1
                break
        else:
            groups.append([p, 1])
    return groups


def frame(a, b):
    e1 = a
    e2 = b - (a @ b) * a
    e2 /= np.linalg.norm(e2)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def rot_to_su2(R):
    """robust rotation-matrix -> unit quaternion (w,x,y,z) -> g = w - i (x sx + y sy + z sz)."""
    tr = np.trace(R)
    cands = [1 + tr, 1 + R[0, 0] - R[1, 1] - R[2, 2], 1 - R[0, 0] + R[1, 1] - R[2, 2], 1 - R[0, 0] - R[1, 1] + R[2, 2]]
    k = int(np.argmax(cands))
    if k == 0:
        w = np.sqrt(cands[0]) / 2
        x, y, z = (R[2, 1] - R[1, 2]) / (4 * w), (R[0, 2] - R[2, 0]) / (4 * w), (R[1, 0] - R[0, 1]) / (4 * w)
    elif k == 1:
        x = np.sqrt(cands[1]) / 2
        w, y, z = (R[2, 1] - R[1, 2]) / (4 * x), (R[0, 1] + R[1, 0]) / (4 * x), (R[0, 2] + R[2, 0]) / (4 * x)
    elif k == 2:
        y = np.sqrt(cands[2]) / 2
        w, x, z = (R[0, 2] - R[2, 0]) / (4 * y), (R[0, 1] + R[1, 0]) / (4 * y), (R[1, 2] + R[2, 1]) / (4 * y)
    else:
        z = np.sqrt(cands[3]) / 2
        w, x, y = (R[1, 0] - R[0, 1]) / (4 * z), (R[0, 2] + R[2, 0]) / (4 * z), (R[1, 2] + R[2, 1]) / (4 * z)
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]], complex); sz = np.array([[1, 0], [0, -1]], complex)
    return w * np.eye(2) - 1j * (x * sx + y * sy + z * sz)


def stabilizer_finite(u):
    grp = multiset(roots_points(u))
    pts = [g[0] for g in grp]
    mult = [g[1] for g in grp]
    Rs = []
    # choose a reference pair (a,b) distinct, not antipodal
    pair = None
    for i in range(len(pts)):
        for j in range(len(pts)):
            if i != j and abs(pts[i] @ pts[j]) < 1 - 1e-9:
                pair = (i, j); break
        if pair: break
    assert pair is not None
    i, j = pair
    F = frame(pts[i], pts[j])
    for k in range(len(pts)):
        for l in range(len(pts)):
            if k == l or mult[k] != mult[i] or mult[l] != mult[j]:
                continue
            if abs(pts[k] @ pts[l] - pts[i] @ pts[j]) > 1e-9:
                continue
            R = frame(pts[k], pts[l]) @ F.T
            ok = True
            for p, mlt in zip(pts, mult):
                q = R @ p
                idx = [n for n in range(len(pts)) if np.linalg.norm(pts[n] - q) < 1e-7]
                if len(idx) != 1 or mult[idx[0]] != mlt:
                    ok = False; break
            if ok and not any(np.linalg.norm(R - R2) < 1e-7 for R2 in Rs):
                Rs.append(R)
    return Rs, grp


out = {}
data = {}
for name, dct in U.items():
    u = uvec(dct)
    rec = {}
    grp = multiset(roots_points(u))
    rec["majorana_points"] = [[list(np.round(p, 6)), m] for p, m in grp]
    # Lie algebra stabilizer: X in su(2) with (i conj J_a) u in C u
    Jz = np.diag(np.arange(-3, 4)).astype(complex)
    Jp = np.zeros((7, 7), complex)
    for m in range(-3, 3):
        Jp[m + 4, m + 3] = np.sqrt((3 - m) * (3 + m + 1))
    Jx = (Jp + Jp.T) / 2; Jy = (Jp - Jp.T) / 2j
    Pp = np.eye(7) - np.outer(u, u.conj())
    M = np.column_stack([Pp @ (1j * np.conj(Jn) @ u) for Jn in (Jx, Jy, Jz)])
    Mr = np.vstack([M.real, M.imag])
    sv = np.linalg.svd(Mr, compute_uv=False)
    lie_dim = int(np.sum(sv < 1e-10))
    rec["lie_stabilizer_dim"] = lie_dim
    if lie_dim == 0:
        Rs, _ = stabilizer_finite(u)
        els = []
        for R in Rs:
            g = rot_to_su2(R)
            v = act(3, g, u)
            lam = np.vdot(u, v)
            res = np.linalg.norm(v - lam * u)
            els.append((g, lam, res))
        rec["order_SO3"] = len(Rs)
        rec["order_SU2_lift"] = 2 * len(Rs)
        rec["max_residual_line_fixed"] = float(max(e[2] for e in els))
        # rotation orders and characters
        orders = []
        for R in Rs:
            k = 1; M_ = R.copy()
            while np.linalg.norm(M_ - np.eye(3)) > 1e-7:
                M_ = M_ @ R; k += 1
            orders.append(k)
        rec["SO3_element_orders"] = sorted(orders)
        chars = [complex(np.round(e[1], 10)) for e in els]
        rec["character_values_on_Phi"] = sorted(set((round(c.real, 8), round(c.imag, 8)) for c in chars))
        # chi-subspace dims in V_J (J = 0..9) under the SU(2) lifts (integer J: -1 acts trivially)
        dims = {}
        for J in range(0, 10):
            Pj = sum(np.conj(e[1]) * np.conj(D_of(J, e[0])) for e in els) / len(els)
            dims[2 * J] = int(round(np.trace(Pj).real))
        rec["chi_subspace_dim_by_level"] = dims
        rec["block_chi_subspace_dim"] = dims[6]
        data[name] = {"elements": [(e[0], e[1]) for e in els]}
    else:
        rec["note"] = "continuous stabilizer; see analytic treatment"
    out[name] = rec
    print(name, json.dumps(rec, default=str))

# ---- continuous cases, analytic + numeric verification
th = 0.7318
gT = np.diag([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
w = su2_axis_angle([0, 1, 0], np.pi)       # pi rotation about y (Weyl element)
u1 = uvec(U["U1"]); u2 = uvec(U["U2"])
lam1 = np.vdot(u1, act(3, gT, u1)); lam2T = np.vdot(u2, act(3, gT, u2)); lam2w = np.vdot(u2, act(3, w, u2))
print("U1 torus: conj(D(exp(-i th Jz))) v3 =", np.round(lam1, 12), " expected e^{3 i th} =", np.round(np.exp(3j * th), 12))
print("U2 torus:", np.round(lam2T, 12), " Weyl (pi about y):", np.round(lam2w, 12))
check("U1: torus acts on Phi by e^{3 i theta}", lambda x: abs(x - np.exp(3j * th)) < 1e-12, lam1, np.exp(-3j * th), "compare with e^{-3i theta}")
check("U2: torus acts trivially on Phi", lambda x: abs(x - 1) < 1e-12, lam2T, lam2T * np.exp(1j * th), "multiply by e^{i theta}")
check("U2: Weyl element acts on Phi by -1", lambda x: abs(x + 1) < 1e-12, lam2w, -lam2w, "negate")
for name in ("U1", "U2"):
    dims = {}
    for J in range(10):
        if name == "U1":
            dims[2 * J] = 1 if J >= 3 else 0            # weight +3 (on conj-action: character e^{3i th}) exists iff J >= 3
        else:
            dims[2 * J] = 1 if J % 2 == 1 else 0        # weight 0 vector v_0, Weyl eigenvalue (-1)^J must equal -1
    # numeric verification of these by torus-averaging and Weyl
    for J in range(10):
        thetas = np.linspace(0, 2 * np.pi, 41)[:-1]
        if name == "U1":
            Pj = sum(np.exp(-3j * t) * np.conj(D_of(J, np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)]))) for t in thetas) / len(thetas)
        else:
            PT = sum(np.conj(D_of(J, np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)]))) for t in thetas) / len(thetas)
            Pj = PT @ (np.eye(2 * J + 1) - np.conj(D_of(J, w))) / 2
        assert int(round(np.trace(Pj).real)) == dims[2 * J], (name, J)
    out[name]["chi_subspace_dim_by_level"] = dims
    out[name]["block_chi_subspace_dim"] = dims[6]
    print(name, "chi-subspace dims by level:", dims)

# which stabilizer groups: element-order census identifies the SO(3) image
check("U3 stabilizer image has 24 elements (octahedral), orders census", lambda o: len(o) == 24 and o.count(4) == 6 and o.count(3) == 8,
      out["U3"]["SO3_element_orders"], out["U4"]["SO3_element_orders"], "use U4's census")
check("U4 stabilizer image has 12 elements with a 6-fold rotation (dihedral D6)", lambda o: len(o) == 12 and o.count(6) == 2 and o.count(2) == 7,
      out["U4"]["SO3_element_orders"], out["U3"]["SO3_element_orders"], "use U3's census")
check("U5 stabilizer image cyclic of order 5", lambda o: o == [1, 5, 5, 5, 5], out["U5"]["SO3_element_orders"], out["U6"]["SO3_element_orders"], "use U6's census")
check("U6 stabilizer image dihedral D3 (order 6: 2 threefold, 3 twofold)", lambda o: o == [1, 2, 2, 2, 3, 3], out["U6"]["SO3_element_orders"], out["U5"]["SO3_element_orders"], "use U5's census")
for name in ("U3", "U4", "U5", "U6"):
    check(f"{name}: every stabilizer element fixes the line of u (residual < 1e-12)", lambda r: r < 1e-12, out[name]["max_residual_line_fixed"], 1e-3, "residual 1e-3")

with open("out/stab.pkl", "wb") as f:
    pickle.dump(data, f)


# ---- explicit generators and their characters on Phi; block chi-subspace bases
def rz(t):
    return np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])


gens = {
    "U1": {"rz(theta)": rz(th)},
    "U2": {"rz(theta)": rz(th), "ry(pi)": su2_axis_angle([0, 1, 0], np.pi)},
    "U3": {"rz(pi/2)": rz(np.pi / 2), "r_(sqrt2,0,1)(2pi/3)": su2_axis_angle([np.sqrt(2), 0, 1], 2 * np.pi / 3),
           "rx(pi)": su2_axis_angle([1, 0, 0], np.pi), "r_(1,1,0)(pi)": su2_axis_angle([1, 1, 0], np.pi)},
    "U4": {"rz(pi/3)": rz(np.pi / 3), "ry(pi)": su2_axis_angle([0, 1, 0], np.pi), "rx(pi)": su2_axis_angle([1, 0, 0], np.pi)},
    "U5": {"rz(2pi/5)": rz(2 * np.pi / 5)},
    "U6": {"rz(2pi/3)": rz(2 * np.pi / 3), "rx(pi)": su2_axis_angle([1, 0, 0], np.pi)},
}
genres = {}
for name, gd in gens.items():
    u = uvec(U[name])
    genres[name] = {}
    for gn, g in gd.items():
        v = act(3, g, u)
        lam = np.vdot(u, v)
        res = np.linalg.norm(v - lam * u)
        genres[name][gn] = (complex(np.round(lam, 12)), float(res))
    print(name, "generator characters (value on Phi, residual):", genres[name])
    out[name]["generators"] = {k: [str(v[0]), v[1]] for k, v in genres[name].items()}
check("U3 generators rz(pi/2), r_(sqrt2,0,1)(2pi/3), rx(pi), r_(1,1,0)(pi) all fix [u] (residual < 1e-12)",
      lambda d: all(v[1] < 1e-12 for v in d.values()), genres["U3"],
      {**genres["U3"], "rz(pi/2)": (0, 1.0)}, "set a residual to 1")
check("U6 generators rz(2pi/3), rx(pi) fix [u]", lambda d: all(v[1] < 1e-12 for v in d.values()), genres["U6"],
      {**genres["U6"], "rx(pi)": (0, 0.5)}, "set a residual to 0.5")
# block chi-subspace for U5, U6: orthonormal basis
for name in ("U5", "U6"):
    els = data[name]["elements"]
    P3 = sum(np.conj(l) * np.conj(D_of(3, g)) for g, l in els) / len(els)
    wv, V = np.linalg.eigh((P3 + P3.conj().T) / 2)
    B = V[:, wv > 0.5]
    print(name, "block chi-subspace basis (columns, index m=-3..3):\n", np.round(B, 6))
    out[name]["block_chi_basis_support_m"] = sorted({m - 3 for m in range(7) if np.linalg.norm(B[m]) > 1e-9})

# ---- item 4: xi transforms with the same character (numerical check on the exact xi)
st = pickle.load(open("out/main_all.pkl", "rb"))
worst = 0.0
for (name, d), S_ in st.items():
    if name not in gens:
        continue
    for gn, g in gens[name].items():
        lam = genres[name][gn][0]
        for J, terms in S_["xi"].items():
            Cx = sum(np.outer(np.array([c.to_complex() for c in x]),
                              np.array([[c.to_complex() for c in col] for col in Y]).reshape(-1)) for x, Y in terms)
            lhs = np.conj(D_of(J, g)) @ Cx
            worst = max(worst, np.abs(lhs - lam * Cx).max())
print("max |r.xi - chi(r) xi| over all U, sectors, levels, generators:", worst)
check("xi transforms by the same character as Phi under all listed stabilizer generators (max dev < 1e-10)",
      lambda w_: w_ < 1e-10, worst, 0.1, "deviation 0.1")
out["xi_character_max_dev"] = worst
with open("out/item3.json", "w") as f:
    json.dump(out, f, indent=1, default=str)
save_log("s03")
