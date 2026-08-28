"""M5.32 rung R0: INDEPENDENT adversarial audit of the six quadratic F-invariants.

Built from the term DEFINITIONS only (the producer script was not read).

EQUATIONS
=========

Field and jet
-------------
    M(x)      4x4 real symmetric, transformation  M -> Lambda M Lambda^T,  x -> Lambda x
              => the raw entries are CONTRAVARIANT:  M^{ab}
    A_mu      = d_mu M                (derivative slot covariant, internal slots contravariant)
    A'_mu     = (Lambda^{-1})^nu_mu  Lambda A_nu Lambda^T
    eta       = diag(-1, 1, 1, 1);  Lambda^T eta Lambda = eta   (so Lambda^{-1} eta Lambda^{-T} = eta)

Curvature
---------
    F_{mu nu}^{ab} = (A_mu eta A_nu)^{ab} - (A_nu eta A_mu)^{ab}
                   = A_mu^{ac} eta_{cd} A_nu^{db} - (mu <-> nu)
    antisymmetric in (mu, nu) and in (a, b) (the second follows from A symmetric).
    F'_{mu nu} = (Lambda^{-1})^rho_mu (Lambda^{-1})^sigma_nu  Lambda F_{rho sigma} Lambda^T

Contraction rule that follows from the index positions
------------------------------------------------------
    derivative-derivative : eta^{mu nu}   (two lower indices)
    internal-internal     : eta_{ab}      (two upper indices)
    derivative-internal   : delta^mu_a    (one lower, one upper)
Storage: F[mu, nu, a, b]; a delta contraction is a plain index identification.

The six candidate scalars (audit definitions)
---------------------------------------------
    I1 = sum_{mu<nu} eta^{mu mu} eta^{nu nu} tr(eta F_{mu nu} eta F_{mu nu}^T)
       = (1/2) F_{mu nu}^{ab} F^{mu nu}_{ab}
    I2 = F_{mu nu}^{ab} F_{ab}^{mu nu}                 (four mixed pairs, all delta)
    I3 = F_{mu nu}^{ab} F^{mu}_{a}{}^{nu}_{b}           (one eta^DD, one eta_II, two deltas)
    Rt^b_nu = F_{mu nu}^{mu b}   (the ONLY mixed trace, up to sign)
    I4 = Rt_{nu b} Rt^{nu b}   = eta^{nu nu} eta_{bb} Rt[nu,b]^2
    I5 = Rt^b_nu Rt^nu_b       = Rt[nu,b] Rt[b,nu]
    I6 = (Rt^mu_mu)^2
The producer's R_ac = sum_mu F[mu, c, a, mu] = -Rt[c, a]; I4..I6 are sign-blind.

Static 3x3 sector
-----------------
    M = diag(-s g) (+) B(x),  B 3x3 symmetric, no time dependence:
    A_0 = 0, A_i = 0 (+) b_i  =>  F_{ij} lives in the 3x3 block, eta restricted = 1,
    F_{ij}^{kl} = eps_{ijm} eps_{kln} G_{mn}  =>  every quadratic SO(3) scalar is one of
    tr(G G^T), tr(G^2), (tr G)^2  =>  rank of {I1..I6} on that sector is at most 3.

omega scan
----------
    A_0 = omega a0  =>  F_{0i} linear in omega, F_{ij} independent  =>  every I_k is a
    polynomial of degree <= 2 in omega:  I = a + b omega + c omega^2  (b = the odd piece).
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETA_D = np.array([-1.0, 1.0, 1.0, 1.0])


# --------------------------------------------------------------------------- jets
def rand_sym(rng, n=4):
    x = rng.standard_normal((n, n))
    return 0.5 * (x + x.T)


def rand_jet(rng):
    """(M, A) with A[mu] symmetric 4x4 (M is not used by F but carried for the law)."""
    return rand_sym(rng), np.stack([rand_sym(rng) for _ in range(4)])


def curvature(A):
    """F[mu, nu, a, b] = (A_mu eta A_nu - A_nu eta A_mu)^{ab}."""
    T = np.einsum("mac,cd,ndb->mnab", A, ETA, A)
    return T - np.transpose(T, (1, 0, 2, 3))


# --------------------------------------------------------------------------- scalars
def inv_I1(F):
    s = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[mu, nu]
            s += ETA_D[mu] * ETA_D[nu] * np.trace(ETA @ Fm @ ETA @ Fm.T)
    return float(s)


def inv_I2(F):
    return float(np.einsum("mnab,abmn->", F, F))


def inv_I3(F):
    return float(np.einsum("m,b,mnab,manb->", ETA_D, ETA_D, F, F))


def mixed_trace(F):
    """Rt[nu, b] = sum_mu F[mu, nu, mu, b]."""
    return np.einsum("mnmb->nb", F)


def inv_I4(F):
    R = mixed_trace(F)
    return float(np.einsum("n,b,nb,nb->", ETA_D, ETA_D, R, R))


def inv_I5(F):
    R = mixed_trace(F)
    return float(np.einsum("nb,bn->", R, R))


def inv_I6(F):
    R = mixed_trace(F)
    return float(np.trace(R) ** 2)


INVS = {"I1": inv_I1, "I2": inv_I2, "I3": inv_I3, "I4": inv_I4, "I5": inv_I5, "I6": inv_I6}


def all_invs(F):
    return np.array([INVS[k](F) for k in INVS])


# the WRONG rule: eta on every raw index of the mixed pairings
def inv_I2_eta_all(F):
    return float(np.einsum("m,n,a,b,mnab,abmn->", ETA_D, ETA_D, ETA_D, ETA_D, F, F))


def inv_I3_eta_all(F):
    return float(np.einsum("m,n,a,b,mnab,manb->", ETA_D, ETA_D, ETA_D, ETA_D, F, F))


def mixed_trace_eta(F):
    """Rt_eta[nu, b] = sum_mu eta_mumu F[mu, nu, mu, b]  (the eta-everywhere mixed trace)."""
    return np.einsum("m,mnmb->nb", ETA_D, F)


def inv_I5_eta_all(F):
    R = mixed_trace_eta(F)
    return float(np.einsum("n,b,nb,bn->", ETA_D, ETA_D, R, R))


def inv_I6_eta_all(F):
    R = mixed_trace_eta(F)
    return float(np.einsum("n,nn->", ETA_D, R) ** 2)


# --------------------------------------------------------------------------- Lorentz
def boost(rng):
    v = rng.standard_normal(3)
    v *= rng.uniform(0.1, 0.8) / np.linalg.norm(v)
    b2 = v @ v
    gam = 1.0 / np.sqrt(1.0 - b2)
    L = np.eye(4)
    L[0, 0] = gam
    L[0, 1:] = -gam * v
    L[1:, 0] = -gam * v
    L[1:, 1:] = np.eye(3) + (gam - 1.0) * np.outer(v, v) / b2
    return L


def rotation(rng):
    q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    L = np.eye(4)
    L[1:, 1:] = q
    return L


def rand_lorentz(rng):
    L = rotation(rng) @ boost(rng) @ rotation(rng)
    assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-10)
    return L


def transform_jet_contra(A, L):
    """A'_mu = (L^{-1})^nu_mu  L A_nu L^T   (M -> L M L^T law)."""
    Linv = np.linalg.inv(L)
    return np.einsum("nm,ac,ncd,bd->mab", Linv, L, A, L)


def transform_jet_cov(A, L):
    """A'_mu = (L^{-1})^nu_mu  L^{-T} A_nu L^{-1}   (M -> L^{-T} M L^{-1} law, for contrast)."""
    Linv = np.linalg.inv(L)
    return np.einsum("nm,ca,ncd,db->mab", Linv, Linv, A, Linv)


# --------------------------------------------------------------------------- helpers
def numeric_rank(V, rtol=1e-9):
    """Rank of a k x N value matrix with the singular spectrum and the null vector."""
    Vn = V / np.linalg.norm(V, axis=1, keepdims=True)
    s = np.linalg.svd(Vn, compute_uv=False)
    rank = int(np.sum(s > rtol * s[0]))
    _, _, vt = np.linalg.svd(Vn.T, full_matrices=False)
    return rank, s.tolist(), vt[-1].tolist()


def static_jet(rng):
    """A_0 = 0, A_i = 0 (+) b_i with b_i random symmetric 3x3 (the gradient of B(x))."""
    A = np.zeros((4, 4, 4))
    for i in range(1, 4):
        A[i, 1:, 1:] = rand_sym(rng, 3)
    return A


def static_G(F):
    """G_{mn} = (1/4) eps_{ijm} eps_{kln} F_{ij}^{kl}."""
    eps = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        eps[p] = np.linalg.det(np.eye(3)[list(p)])
    Fs = F[1:, 1:, 1:, 1:]
    return 0.25 * np.einsum("ijm,kln,ijkl->mn", eps, eps, Fs)


def enumerate_pairings():
    """Every perfect matching of the 8 slots of F (x) F, kept if no pair is a same-tensor
    derivative-derivative or internal-internal pair (those vanish by antisymmetry).
    Slot labels: 0,1 = F1 derivative; 2,3 = F1 internal; 4,5 = F2 derivative; 6,7 = F2 internal.
    """
    kind = {0: "D", 1: "D", 2: "I", 3: "I", 4: "D", 5: "D", 6: "I", 7: "I"}
    tensor = {s: (0 if s < 4 else 1) for s in range(8)}

    def matchings(slots):
        if not slots:
            yield []
            return
        a = slots[0]
        for b in slots[1:]:
            rest = [s for s in slots if s not in (a, b)]
            for m in matchings(rest):
                yield [(a, b)] + m

    out = []
    for m in matchings(list(range(8))):
        ok = True
        for a, b in m:
            if tensor[a] == tensor[b] and kind[a] == kind[b]:
                ok = False
        if ok:
            out.append(m)
    return out


def eval_pairing(F, m):
    """Evaluate a pairing with eta on DD and II pairs and delta on DI pairs."""
    letters = "abcdefgh"
    lab = list(letters)
    weights = []
    kind = "DDIIDDII"
    for a, b in m:
        lab[b] = lab[a]
        if kind[a] == kind[b]:
            weights.append(lab[a])
    sub = "".join(lab[:4]) + "," + "".join(lab[4:])
    ops = [F, F]
    for w in weights:
        sub = w + "," + sub
        ops.insert(0, ETA_D)
    return float(np.einsum(sub + "->", *ops))


# --------------------------------------------------------------------------- audit
def main():
    t0 = time.time()
    rng = np.random.default_rng(20260827)
    out = {"seed": 20260827, "claims": {}}
    names = list(INVS)

    # ---------------- C1: covariance drift
    n_l = 40
    drift = {k: 0.0 for k in names}
    drift_eta_all = {"I2": 0.0, "I3": 0.0, "I5": 0.0, "I6": 0.0}
    drift_cov_law = {k: 0.0 for k in names}
    drift_cov_law_eta_all = {"I2": 0.0, "I3": 0.0, "I5": 0.0, "I6": 0.0}
    bad = {"I2": inv_I2_eta_all, "I3": inv_I3_eta_all, "I5": inv_I5_eta_all, "I6": inv_I6_eta_all}
    for _ in range(n_l):
        _, A = rand_jet(rng)
        L = rand_lorentz(rng)
        F0, F1 = curvature(A), curvature(transform_jet_contra(A, L))
        F2 = curvature(transform_jet_cov(A, L))
        for k in names:
            v0, v1, v2 = INVS[k](F0), INVS[k](F1), INVS[k](F2)
            drift[k] = max(drift[k], abs(v1 - v0) / max(1.0, abs(v0)))
            drift_cov_law[k] = max(drift_cov_law[k], abs(v2 - v0) / max(1.0, abs(v0)))
        for k in bad:
            v0, v1, v2 = bad[k](F0), bad[k](F1), bad[k](F2)
            drift_eta_all[k] = max(drift_eta_all[k], abs(v1 - v0) / max(1.0, abs(v0)))
            drift_cov_law_eta_all[k] = max(drift_cov_law_eta_all[k], abs(v2 - v0) / max(1.0, abs(v0)))
    # equivalence of the two bookkeepings: M_cov = eta M eta
    _, A = rand_jet(rng)
    Fa = curvature(A)
    Fb = curvature(np.einsum("ac,ncd,db->nab", ETA, A, ETA))
    equiv = {k: abs(bad[k](Fb) - INVS[k](Fa)) for k in bad}
    out["claims"]["C1"] = {
        "n_lorentz": n_l,
        "drift_delta_rule_under_contra_law": drift,
        "drift_eta_all_rule_under_contra_law": drift_eta_all,
        "drift_delta_rule_under_cov_law": drift_cov_law,
        "drift_eta_all_rule_under_cov_law": drift_cov_law_eta_all,
        "eta_all_on_(eta M eta)_minus_delta_on_M": equiv,
    }

    # ---------------- C2: traces
    n = 30
    maxDD = maxII = 0.0
    asym = []
    mixed_forms = []
    for _ in range(n):
        _, A = rand_jet(rng)
        F = curvature(A)
        maxDD = max(maxDD, np.abs(np.einsum("m,mmab->ab", ETA_D, F)).max())
        maxII = max(maxII, np.abs(np.einsum("a,mnaa->mn", ETA_D, F)).max())
        R = mixed_trace(F)
        asym.append(np.linalg.norm(R - R.T) / np.linalg.norm(R))
        # the four mixed traces (slot0-slot2, slot0-slot3, slot1-slot2, slot1-slot3)
        t = [np.einsum("mnmb->nb", F), np.einsum("mnam->na", F),
             np.einsum("mnnb->mb", F), np.einsum("mnan->ma", F)]
        mixed_forms.append([np.linalg.norm(t[0] + t[1]), np.linalg.norm(t[0] + t[2]),
                            np.linalg.norm(t[0] - t[3])])
    out["claims"]["C2"] = {
        "max_abs_derivative_pair_trace": float(maxDD),
        "max_abs_internal_pair_trace": float(maxII),
        "R_asymmetry_ratio_min": float(min(asym)),
        "R_asymmetry_ratio_max": float(max(asym)),
        "mixed_traces_reduce_to_one_up_to_sign_max_residual": float(np.max(mixed_forms)),
    }

    # ---------------- C3: static sector
    N = 300
    Vs = np.zeros((6, N))
    Gq = np.zeros((3, N))
    frob_vs_eta = 0.0
    ratios = {k: [] for k in names[1:]}
    for j in range(N):
        A = static_jet(rng)
        F = curvature(A)
        Vs[:, j] = all_invs(F)
        G = static_G(F)
        Gq[:, j] = [np.sum(G * G), np.trace(G @ G), np.trace(G) ** 2]
        # I1 with Frobenius internal metric
        s = 0.0
        for mu in range(4):
            for nu in range(mu + 1, 4):
                s += ETA_D[mu] * ETA_D[nu] * np.sum(F[mu, nu] ** 2)
        frob_vs_eta = max(frob_vs_eta, abs(s - Vs[0, j]) / max(1.0, abs(s)))
        for i, k in enumerate(names[1:], start=1):
            ratios[k].append(Vs[i, j] / Vs[0, j])
    rank_s, sv_s, null_s = numeric_rank(Vs)
    rank_s5, sv_s5, null_s5 = numeric_rank(Vs[1:])
    # express each I_k as a fixed combination of the three G-invariants
    coeff = {}
    resid = {}
    for i, k in enumerate(names):
        c, *_ = np.linalg.lstsq(Gq.T, Vs[i], rcond=None)
        coeff[k] = c.tolist()
        resid[k] = float(np.linalg.norm(Gq.T @ c - Vs[i]) / np.linalg.norm(Vs[i]))
    out["claims"]["C3"] = {
        "n_fields": N,
        "I1_frobenius_vs_eta_max_rel_diff": float(frob_vs_eta),
        "ratio_to_I1_min_max": {k: [float(min(v)), float(max(v))] for k, v in ratios.items()},
        "ratio_to_I1_std": {k: float(np.std(v)) for k, v in ratios.items()},
        "identically_zero": {k: bool(np.abs(Vs[i]).max() < 1e-12) for i, k in enumerate(names)},
        "rank_I1_to_I6_static": rank_s,
        "singular_values_I1_to_I6_static": sv_s,
        "null_vector_I1_to_I6_static_normalized_rows": null_s,
        "rank_I2_to_I6_static": rank_s5,
        "singular_values_I2_to_I6_static": sv_s5,
        "null_vector_I2_to_I6_static_normalized_rows": null_s5,
        "G_basis": ["tr(G G^T)", "tr(G^2)", "(tr G)^2"],
        "coeff_on_G_basis": coeff,
        "fit_residual_on_G_basis": resid,
    }

    # ---------------- C4 + C6: generic jets
    N = 400
    Vg = np.zeros((6, N))
    pairings = enumerate_pairings()
    Vp = np.zeros((len(pairings), N))
    for j in range(N):
        _, A = rand_jet(rng)
        F = curvature(A)
        Vg[:, j] = all_invs(F)
        for p, m in enumerate(pairings):
            Vp[p, j] = eval_pairing(F, m)
    rank_g, sv_g, null_g = numeric_rank(Vg)
    rank_g5, sv_g5, null_g5 = numeric_rank(Vg[1:])
    rank_p, sv_p, _ = numeric_rank(Vp)
    rank_all, sv_all, _ = numeric_rank(np.vstack([Vg, Vp]))
    # express each pairing on I1..I6
    pair_fit = []
    for p, m in enumerate(pairings):
        c, *_ = np.linalg.lstsq(Vg.T, Vp[p], rcond=None)
        r = float(np.linalg.norm(Vg.T @ c - Vp[p]) / max(1e-300, np.linalg.norm(Vp[p])))
        pair_fit.append({"pairing": m, "coeff_on_I1..I6": np.round(c, 6).tolist(),
                         "residual": r, "max_abs_value": float(np.abs(Vp[p]).max())})
    out["claims"]["C6"] = {
        "n_jets": N,
        "rank_I1_to_I6_generic": rank_g,
        "singular_values_I1_to_I6_generic": sv_g,
        "null_vector_generic": null_g,
        "rank_I2_to_I6_generic": rank_g5,
        "singular_values_I2_to_I6_generic": sv_g5,
        "n_pairings_enumerated_nonvanishing_type": len(pairings),
        "rank_of_all_pairings": rank_p,
        "rank_of_pairings_plus_I1_to_I6": rank_all,
        "pairing_fits": pair_fit,
    }

    # C4: omega polynomial + odd piece
    n = 40
    Bmax = {k: 0.0 for k in names}
    Brel = {k: 0.0 for k in names}
    cubic_resid = {k: 0.0 for k in names}
    omegas = np.array([-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0])
    for _ in range(n):
        _, A = rand_jet(rng)
        a0 = rand_sym(rng)
        vals = np.zeros((6, len(omegas)))
        for i, w in enumerate(omegas):
            Aw = A.copy()
            Aw[0] = w * a0
            vals[:, i] = all_invs(curvature(Aw))
        V = np.vander(omegas, 4, increasing=True)
        for i, k in enumerate(names):
            c, *_ = np.linalg.lstsq(V, vals[i], rcond=None)
            cubic_resid[k] = max(cubic_resid[k], abs(c[3]) / max(1.0, np.abs(c).max()))
            Bmax[k] = max(Bmax[k], abs(c[1]))
            Brel[k] = max(Brel[k], abs(c[1]) / max(1e-300, np.abs(c[:3]).max()))
    # the odd piece with a STATIC-compatible base (time row of A_i zero): does it survive?
    Bstatic = {k: 0.0 for k in names}
    for _ in range(n):
        A = static_jet(rng)
        a0 = rand_sym(rng)
        vals = np.zeros((6, len(omegas)))
        for i, w in enumerate(omegas):
            Aw = A.copy()
            Aw[0] = w * a0
            vals[:, i] = all_invs(curvature(Aw))
        V = np.vander(omegas, 4, increasing=True)
        for i, k in enumerate(names):
            c, *_ = np.linalg.lstsq(V, vals[i], rcond=None)
            Bstatic[k] = max(Bstatic[k], abs(c[1]))
    out["claims"]["C4"] = {
        "n_jets": n,
        "omega_grid": omegas.tolist(),
        "cubic_coefficient_rel_max": cubic_resid,
        "B_omega_odd_abs_max": Bmax,
        "B_omega_odd_rel_max": Brel,
        "B_omega_odd_abs_max_static_base": Bstatic,
    }

    out["runtime_s"] = time.time() - t0
    return out


if __name__ == "__main__":
    res = main()
    if "--json" in sys.argv:
        Path("data/m5_32_r0_audit.json").write_text(json.dumps(res, indent=1, default=float))
        print("wrote data/m5_32_r0_audit.json")
    else:
        for c, v in res["claims"].items():
            print("==", c)
            for kk, vv in v.items():
                if kk != "pairing_fits":
                    print("  ", kk, ":", vv)
            if "pairing_fits" in v:
                for pf in v["pairing_fits"]:
                    print("   ", pf)
        print("runtime", res["runtime_s"])
