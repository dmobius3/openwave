"""S14 item 6: control (i)'s nonzero-momenta setup assertion, RUN pre-freeze.
The pinned phi = normalized sum of the 13 canonical basis vectors of H_{R0,12}.
Canonical basis construction (the deterministic S 3 endpoint on this space): the R0
intertwiner at level 12 is the 1-dim fixed space of {pi_12(gamma)^T}; sections
f_j(x) = v_inv^T pi_12(x) e_j; in the unitary D-basis Schur makes the f_j orthogonal with
equal norms, so Lowdin = normalization and the canonical multiplet index j IS the weight
index of V_12^* in row-major (canonical) order. Right generators on the multiplet: the
spin-6 angular momentum matrices. Momenta for the linear standing wave psi = e^{i omega t}
phi: M_a = Re<psi_dot, T_a psi> with T_a the anti-Hermitian generators; S_cons per S 7.
Assertion: |M_a| >= 1e-6 * S_cons for a = 1,2,3. This item REPORTS; a red is a protocol
misconfiguration caught pre-freeze, which is the reason the obligation exists.
"""
import numpy as np, json, hashlib, sys

s = 6                                     # spin of V_12^*
m = np.arange(s, -s-1, -1, dtype=float)   # weights +6..-6, canonical (row-major) order
J3 = np.diag(m)
ap = np.sqrt((s - m[1:]) * (s + m[1:] + 1))      # J+ |m> = a |m+1>: subdiag placement
Jp = np.zeros((13,13)); Jp[np.arange(12), np.arange(1,13)] = ap
Jm = Jp.T.copy()
# su(2) check
comm = Jp@Jm - Jm@Jp - 2*J3
assert np.abs(comm).max() < 1e-12
T1 = 0.5*(Jp - Jm)                        # anti-Hermitian (real antisymmetric here)
T2 = 0.5j*(Jp + Jm)                       # anti-Hermitian
T3 = 1j*J3                                # anti-Hermitian
omega = np.sqrt(168.0)

phi = np.ones(13, dtype=complex); phi /= np.linalg.norm(phi)
psid = 1j*omega*phi
S_cons = np.linalg.norm(phi)*np.linalg.norm(psid) + np.linalg.norm(phi)**2
thr = 1e-6*S_cons
M = [float(np.real(np.vdot(psid, T@phi))) for T in (T1,T2,T3)]
live = [abs(x) >= thr for x in M]
print(f"  pinned phi = equal sum of the 13 canonical (weight) basis vectors")
print(f"  S_cons = {S_cons:.6f}, threshold = {thr:.3e}")
for a,(v,ok) in enumerate(zip(M,live),1):
    print(f"  M_{a} = {v:+.6e}   {'LIVE' if ok else 'ZERO / BELOW THRESHOLD'}")
# diagnosis, sharpened: any REAL phi has M_1 = 0 identically (T1 is real antisymmetric,
# so <i omega phi, T1 phi> = Re(i * real) = 0), and any weight-symmetric phi has M_3 = 0.
# A viable pin needs COMPLEX weights. Deterministic candidate for the protocol amendment:
w = (1.0 + np.arange(1,14)/13.0) * np.exp(1j*np.arange(1,14)/3.0)
phi2 = w/np.linalg.norm(w)
psid2 = 1j*omega*phi2
M2 = [float(np.real(np.vdot(psid2, T@phi2))) for T in (T1,T2,T3)]
live2 = [bool(abs(x) >= thr) for x in M2]
print(f"  candidate amendment phi' ~ sum (1 + j/13) e^(i j/3) v_j: M = "
      f"[{M2[0]:+.3e}, {M2[1]:+.3e}, {M2[2]:+.3e}] -> {['LIVE' if k else 'zero' for k in live2]}")
res = dict(M=M, live=[bool(x) for x in live], S_cons=float(S_cons), threshold=float(thr),
           diagnosis="pinned equal-weight phi: M1 = 0 exactly (real phi), M3 = 0 (weight symmetry)",
           amendment_candidate="normalized sum of (1 + j/13) exp(i j/3) v_j",
           amendment_candidate_M=M2, amendment_live=live2)
json.dump(res, open("raw/item6_momenta.json","w"), indent=1)
ok = all(live)
print(f"  ITEM 6 VERDICT: {'GREEN' if ok else 'RED, protocol phi pin fails its own liveness assertion'}")
sys.exit(0 if ok else 1)
