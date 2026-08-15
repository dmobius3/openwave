# M9.28: first law at fixed \(H\) --- 1d enclosed, 3d balls CHM

> Paper 36: a Hamiltonian point source is the wrong variation
> (\(\delta S\neq\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\)).
> This note never changes \(H\). The instrument now passes.
> The kernel does not: 1d still prefers enclosed energy.
> 3d balls prefer CHM. That is not Einstein.

## Equations

Hop Hamiltonian \(H\) is fixed. Occupied packet \(L\) is the
occupied projection of a Gaussian. Unoccupied packet \(R\) is
the unoccupied projection of the bipartite stagger
(\(\Gamma\psi(x)=(-1)^x\psi(x)\) in 1d,
\((-1)^{x+y+z}\) in 3d). State

\[
C(\alpha)=C_0+\alpha\bigl(\lvert R\rangle\langle R\rvert
-\lvert L\rangle\langle L\rvert\bigr).
\]

Site energy of the same \(H\)

\[
e_i=\sum_j H_{ij}C_{ij},\qquad
P_w=\sum_{i\in A}w_i\,\delta e_i.
\]

1d: \(w=(L/2)^2-x^2\) or \(w=1\).
3d ball: \(w=R^2-r^2\), \(R-r\), or \(1\).

Peschel \(S\) and \(K=\log((1-C)/C)\). Instrument identity

\[
\delta S\;\stackrel{?}{=}\;\mathrm{Tr}(K_{\mathrm{vac}}\Delta C).
\]

Discarded, not scored: coherent particle-hole rotation
(first-order \(\langle H\rangle\) vanishes by bipartite
selection; \(C_{\mathrm{vac}}\) fails independently of
\(\theta\)); first-order \(dC\) from a potential
(\(\mathrm{Tr}(H_0\,dC)=0\) by Hellmann).

## 1d gates (\(N=200\), \(L=16\), \(\sigma=2\), \(\alpha=0.02\))

| Gate | Lock | Result |
| --- | --- | --- |
| C_vac | \(\lvert\rho(\delta S,\mathrm{Tr}(K_{\mathrm{vac}}\Delta C))\rvert>0.95\) | PASS \(0.957\) |
| C0 | \(\max\|\delta S\|>10^{-6}\) | PASS \(0.196\) |
| C1 | Pearson \(\delta S(\alpha),\delta S(2.5\alpha)>0.95\) | PASS \(0.9998\) |
| C2 PRIMARY | \(R_{\mathrm{CHM}}<R_{\mathrm{flat}}\) | **FAIL** \(0.266>0.060\) |
| C4 | \(\lvert\rho_{\mathrm{CHM}}\rvert>0.60\) | PASS \(0.964\) |

\(\rho_{\mathrm{flat}}=0.998\). Verdict `1D_FIXEDH_FLAT`.

Auditor \(N=160\), \(L=12\), \(\sigma=2.5\), \(\alpha=0.03\):
C_vac **CONFIRMED** (\(0.975\)). C2 **REFUTED**
(\(R=0.199>0.038\)).

## 3d gates (\(N=12\), \(R=2\), \(512\) balls, \(\sigma=1.5\))

| Gate | Result |
| --- | --- |
| C_vac | PASS \(0.999\) |
| C0 / C1 | PASS \(0.118\) / \(0.9995\) |
| C2 CHM \(<\) flat | **PASS** \(0.018<0.102\) |
| C3 CHM \(<\) linear | PASS \(0.018<0.029\) |
| C4 \(\lvert\rho_{\mathrm{CHM}}\rvert\) | PASS \(0.9998\) |

Verdict `3D_FIXEDH_CHM_WINS`. Auditor \(N=10\), \(216\) balls:
C_vac **CONFIRMED** (\(0.997\)). C2 **CONFIRMED**
(\(0.034<0.123\)).

## What this is

The vacuum first law now holds for this variation. In 1d,
where Casini--Huerta--Myers is a theorem for the *operator*
\(K\), the pairing of this excitation against \(K\) is
whether the packet sits in the interval, not the CHM
weight. On 3d balls of radius \(2\), the same pairing
selects \(R^2-r^2\) over flat, with a small residual gap
over \(R-r\).

That is a linear functional of local energy, with a
kernel that depends on the region. It is not
\(G_{\mu\nu}=8\pi G T_{\mu\nu}\). It is not FGHMV in AdS.
It is not de Sitter. Paper 25 (shape of \(K\)) is
untouched. Paper 35 remains an instrument failure.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(C(\alpha)\), \(e_i\), \(S\), \(K\), 1d gates | `scripts/m9_28_1d_state.py` |
| 1d adversary | `scripts/m9_28_audit_1d.py` |
| 3d balls | `scripts/m9_28_3d_state.py` |
| 3d adversary | `scripts/m9_28_audit_3d.py` |

Paper: [`../latex/37_FixedH_First_Law.tex`](../latex/37_FixedH_First_Law.tex).
