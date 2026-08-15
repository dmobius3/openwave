# M9.23: Bloch dimer vs CHM --- complete covering

> The guess: modular hop tracks a Bloch 3-sphere, not CHM.
> Test: every NN dimer in the ball. One eigh. No subsample.

## Bloch coordinates

For the real \(2\times 2\) correlator of a dimer \((i,j)\):

\[
r_x=2C_{ij},\qquad r_z=C_{ii}-C_{jj}.
\]

Pauli coordinates of a correlator, not a qubit state. Not
identified with curvature.

## Covering

| Grid | Dimers | \(\rho(K,r_x)\) | \(R(K,r_x)\) | \(\rho(K,w_{\mathrm{CHM}})\) | \(R(K,w)\) |
| --- | --- | --- | --- | --- | --- |
| 1d \(L=32\) | 31 | \(0.091\) | \(0.996\) | \(-0.942\) | \(0.334\) |
| \(N=16,R=5\) | **1302** | \(0.034\) | \(0.999\) | \(-0.987\) | \(0.162\) |
| \(N=12,R=4\) | **624** | \(-0.033\) | \(0.999\) | \(-0.988\) | \(0.156\) |
| auditor 1d | 23 | \(0.106\) | \(0.994\) | --- | \(0.177\) |
| auditor \(N=14,R=5\) | **1302** | \(0.031\) | \(1.000\) | --- | \(0.163\) |

\(r_x\) is almost constant (\(\mathrm{std}/\mathrm{mean}\sim 1\%\)).
\(\langle r_z\rangle=0\) to \(10^{-17}\).

## Verdict

C2 PRIMARY **FAIL**. Auditor **REFUTED**. CHM predicts \(K\)
(\(R\sim 0.16\)). Bloch \(r_x\) does not (\(R\sim 1\),
\(\lvert\rho\rvert\sim 0.03\)).

The Bloch 3-sphere guess, as a description of the modular hop
on the diamond, is **closed**. More dimers will not change a
flat \(r_x\). Curvature-as-axis is still a guess and was not
tested as physics.

Scripts: `m9_23_bloch.py`, `m9_23_audit_bloch.py`.
Paper: [`../latex/32_Bloch_Not_CHM.tex`](../latex/32_Bloch_Not_CHM.tex).
