"""Stage 2, task 1 (item 3): test the solver's closed-form stabilizer generators and characters
with MY D-matrix code (numD.py) and my fibre vectors; close the generated groups and count their
orders; compare the characters with my stage-1 generator characters (out/item3.json)."""
import json
import numpy as np
from numD import D_of, su2_axis_angle

s = np.sqrt
U = {"U1": {3: 1.0}, "U2": {0: 1.0}, "U3": {2: 1.0, -2: 1.0}, "U4": {3: 1.0, -3: 1.0},
     "U5": {2: s(13 / 25), -3: s(12 / 25)}, "U6": {3: 1.0, 0: s(23 / 10), -3: 1.0}}


def uvec(d):
    v = np.zeros(7, complex)
    for m, c in d.items():
        v[m + 3] = c
    return v / np.linalg.norm(v)


def act(g, u):
    return np.conj(D_of(3, g)) @ u


def rz(a):
    return np.diag([np.exp(-1j * a / 2), np.exp(1j * a / 2)])


pi = np.pi
# the solver's closed-form generators and claimed characters (common.stabilizer_generators)
gens = {
    "U1": [(rz(1.0), np.exp(3j))],
    "U2": [(rz(1.0), 1.0), (su2_axis_angle([1, 0, 0], pi), -1.0)],
    "U3": [(rz(pi / 2), -1.0), (su2_axis_angle([1, 1, 0], pi / 2), -1.0)],
    "U4": [(rz(pi / 3), -1.0), (su2_axis_angle([1, 0, 0], pi), -1.0)],
    "U5": [(rz(2 * pi / 5), np.exp(4j * pi / 5))],
    "U6": [(rz(2 * pi / 3), 1.0), (su2_axis_angle([1, 0, 0], pi), -1.0)],
}


def close(gs, cap=400):
    els = [(np.eye(2, dtype=complex), 1.0 + 0j)]
    fr = list(els)
    while fr:
        new = []
        for e, c in fr:
            for g, ch in gs:
                h, ch2 = e @ g, c * ch
                if not any(np.abs(h - x).max() < 1e-9 for x, _ in els):
                    els.append((h, ch2)); new.append((h, ch2))
                    if len(els) > cap:
                        return els
        fr = new
    return els


out = {}
for k, gl in gens.items():
    u = uvec(U[k])
    rec = {}
    worst = 0.0
    for g, chi in gl:
        v = act(g, u)
        lam = np.vdot(u, v)
        worst = max(worst, np.linalg.norm(v - chi * u))
    rec["max |r.u - chi u| over solver generators (my D)"] = worst
    if k not in ("U1", "U2"):
        els = close(gl)
        rec["order of group generated (my closure)"] = len(els)
        # character is a homomorphism: every product acts by the product of generator characters
        rec["max |h.u - chi(h) u| over the whole group"] = max(np.linalg.norm(act(h, u) - ch * u) for h, ch in els)
    out[k] = rec
    print(k, rec)

# my stage-1 generator characters vs the solver's description on the same elements
mine = json.load(open("out/item3.json"))
extra = {
    "U2 Rx(pi) (solver) vs ry(pi) (mine)": (su2_axis_angle([1, 0, 0], pi), su2_axis_angle([0, 1, 0], pi), "U2"),
    "U3 rz(pi/2)": (rz(pi / 2), None, "U3"),
    "U4 ry(pi) (solver says +1)": (su2_axis_angle([0, 1, 0], pi), None, "U4"),
    "U6 rx(pi)": (su2_axis_angle([1, 0, 0], pi), None, "U6"),
}
for lab, (g1, g2, k) in extra.items():
    u = uvec(U[k])
    l1 = np.vdot(u, act(g1, u))
    l2 = np.vdot(u, act(g2, u)) if g2 is not None else None
    print(lab, np.round(l1, 12), "" if l2 is None else np.round(l2, 12))
    out[lab] = [str(np.round(l1, 12)), None if l2 is None else str(np.round(l2, 12))]
print("my stage-1 generator characters:", {k: mine[k].get("generators") for k in ("U2", "U3", "U4", "U5", "U6")})
json.dump(out, open("out2/item3_crosscheck.json", "w"), indent=1, default=str)
