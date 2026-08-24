"""E_ρ bundle Laplacian on S³/2I via orbit transport.

P0.1 DERIVATION (filed here as the implementation):

The associated vector bundle E_ρ = (S³ × W_ρ) / 2I, where the equivalence
is (x, w) ~ (γx, ρ(γ)w) for γ ∈ 2I.

A section ψ of E_ρ lifts to an equivariant function ψ̃: S³ → W_ρ satisfying

    ψ̃(γ · x) = ρ(γ) ψ̃(x)     for all γ ∈ 2I, x ∈ S³.

DERIVATION: An element of E_ρ at [x] ∈ S³/2I is an equivalence class
[(x, w)]. A section σ([x]) = [(x, ψ̃(x))] defines a lift ψ̃. Well-definedness
requires [(γx, ψ̃(γx))] = [(x, ψ̃(x))], i.e. (γx, ψ̃(γx)) ~ (x, ψ̃(x)).
Using the equivalence relation with group element γ:
    (x, ψ̃(x)) ~ (γx, ρ(γ) ψ̃(x))
so ψ̃(γx) = ρ(γ) ψ̃(x).  This is the LEFT equivariance law.

Not ρ(γ)⁻¹, not ρ(γ)†, not ρ(γ)ᵀ: the equivalence relation dictates ρ(γ).
The convention is set by the bundle construction, not by choice.

The fibre transport law for the quotient operator is then:
    (L_ρ)_{a,b} = Σ_{j ∈ stencil(a), orbit(j)=b} w_j · ρ(gid(j))

where w_j are the SCALAR RBF-FD weights (geometry, shared), and ρ(gid(j))
is the fibre map for the group element labelling node j within orbit b.
"""

import numpy as np
from scipy.spatial import cKDTree
from .algebra import qmul
from .rbffd import rbf_row


def orbit_stencils(X, oid, gid, elems, k=110):
    """Orbit-transported stencil selection (shared geometry).

    Returns (node_of, plan) where plan[orbit] = (rep_idx, stencil_idx, moved).
    """
    G = len(elems)
    N = len(X)
    node_of = {}
    for i in range(N):
        node_of[(int(oid[i]), int(gid[i]))] = i

    mult = _multiplication_table_1sided(elems)

    tree = cKDTree(X)
    plan = {}
    for orbit in sorted(set(int(o) for o in oid)):
        rep = node_of[(orbit, 0)]
        _, idx = tree.query(X[rep], k=k)
        idx = np.asarray(idx)
        moved = []
        for g in range(G):
            transported = np.array([
                node_of[(int(oid[s]), mult[g][int(gid[s])])]
                for s in idx])
            moved.append(transported)
        plan[orbit] = (rep, idx, moved)
    return node_of, mult, plan


def _multiplication_table_1sided(elems, tol=1e-8):
    """Multiplication table for 2I elements (exact matching)."""
    coords = np.array(elems)
    G = len(elems)

    def find(p):
        diffs = np.max(np.abs(coords - p), axis=1)
        best = int(np.argmin(diffs))
        if diffs[best] < tol:
            return best
        raise ValueError(f"product not found, best dist {diffs[best]:.2e}")

    return [[find(qmul(elems[i], elems[j]))
             for j in range(G)]
            for i in range(G)]


def build_L_scalar(X, oid, gid, elems, k=110, m=7, p=4):
    """Scalar Laplacian by orbit transport (trivial fibre, for P0.2)."""
    node_of, mult, plan = orbit_stencils(X, oid, gid, elems, k)
    G = len(elems)
    N = len(X)
    L = np.zeros((N, N))

    for orbit, (rep, idx, moved) in plan.items():
        w = rbf_row(X[rep], X[idx], m, p)
        for g in range(G):
            tgt = node_of[(orbit, g)]
            L[tgt, moved[g]] = w
    return -L


def build_L_bundle(X, oid, gid, elems, rho, k=110, m=7, p=4,
                   mode="correct"):
    """Bundle Laplacian L_ρ by orbit transport with fibre map ρ(γ).

    The quotient operator on the orbit representatives (seeds) is:
        (L_ρ)_{a,b} = Σ_{j: orbit(j)=b} w_j · F(gid(j))

    where F is the fibre map selected by `mode`:
        "correct":   ρ(γ)           — the derived law
        "inverse":   ρ(γ)⁻¹ = ρ(γ)† (unitary)
        "transpose": ρ(γ)ᵀ
        "conjugate": conj(ρ(γ))
        "omitted":   I              — no fibre map
        "wrong_g":   ρ(γ') for γ' ≠ γ (cyclic shift)

    Only "correct" is admissible; the rest are mutation arms.
    """
    node_of, mult, plan = orbit_stencils(X, oid, gid, elems, k)
    N_seeds = len(plan)
    d_rho = rho[0].shape[0]
    G = len(elems)

    def fibre_map(gi):
        R = rho[gi]
        if mode == "correct":
            return R
        elif mode == "inverse":
            return R.conj().T
        elif mode == "transpose":
            return R.T
        elif mode == "conjugate":
            return R.conj()
        elif mode == "omitted":
            return np.eye(d_rho, dtype=complex)
        elif mode == "wrong_g":
            return rho[(gi + 1) % G]
        else:
            raise ValueError(f"unknown mode {mode!r}")

    seed_orbits = sorted(plan.keys())
    seed_map = {o: i for i, o in enumerate(seed_orbits)}

    L_q = np.zeros((N_seeds * d_rho, N_seeds * d_rho), dtype=complex)

    for orbit, (rep, idx, moved) in plan.items():
        w = rbf_row(X[rep], X[idx], m, p)
        a = seed_map[orbit]
        for s in range(len(idx)):
            j = idx[s]
            b = seed_map[int(oid[j])]
            gj = int(gid[j])
            F = fibre_map(gj)
            L_q[a*d_rho:(a+1)*d_rho, b*d_rho:(b+1)*d_rho] += w[s] * F

    return -L_q, seed_orbits


def oracle_row(X, oid, gid, elems, rho, node_idx, k=110, m=7, p=4):
    """P0.4 independent oracle: recompute RBF-FD weights at node_idx directly.

    Does NOT use orbit transport. Computes the stencil at node_idx by
    independent nearest-neighbour selection and RBF-FD, then returns
    the d_ρ-vector result of applying the operator to a test section.

    Returns (stencil_indices, weights, fibre_blocks) where fibre_blocks[s]
    is the d_ρ × d_ρ block ρ(gid(j_s)) for each stencil node j_s.
    """
    tree = cKDTree(X)
    _, idx = tree.query(X[node_idx], k=k)
    idx = np.asarray(idx)
    w = rbf_row(X[node_idx], X[idx], m, p)

    d_rho = rho[0].shape[0]
    blocks = []
    for s in range(len(idx)):
        gj = int(gid[idx[s]])
        blocks.append(rho[gj])

    return idx, w, blocks


TRANSPORT_MODES = {
    "correct": "ρ(γ), the derived law from E_ρ = (S³ × W_ρ)/2I",
    "inverse": "ρ(γ)⁻¹ = ρ(γ)†, the inverse/adjoint",
    "transpose": "ρ(γ)ᵀ",
    "conjugate": "conj(ρ(γ))",
    "omitted": "I, fibre map dropped entirely",
    "wrong_g": "ρ(γ'), wrong group element (cyclic shift)",
}
