# M9.12: required A2 tests (1d scalar, 2d Dirac, 3d Dirac)

> Same locked gate as M9.10: \(R(m)/R(0)<2\) for \(0<mL\le 8\).
> Local ansatz = on-site + lattice nearest neighbours.

## Verdicts

| Field | C2 (primary) | Notes |
| --- | --- | --- |
| 1d staggered fermion (M9.10) | PASS (ratios \(\le 1.32\)) | already recorded |
| 1d massive scalar | **no verdict** | instrument rejected: \(\xi\) pinned at \(1/2+\varepsilon\), \(R_\pi(0)=0.95\). Not a physics fail |
| 2d staggered Dirac | PASS | ratios 1.19, 1.29, 1.38, 1.38. Mutation fires. Auditor (N=28, L=10): 1.35, 1.29 CONFIRMED |
| 3d staggered Dirac | PASS | ratios 0.45, 0.47, 0.44 (mass *shortens* range). Grid \(10^3\), region \(6^3\) |

A2 is **not refuted** on 1d/2d/3d free lattice Dirac. Remainder in 2d still *grows* with \(m\) (as in 1d) and stays under the locked factor of 2.

## What this is not

Not a 4d continuum diamond. Not the Standard Model. Not a proof of Jacobson 2016. Not FGHMV in de Sitter.

Scripts: `m9_12_A2_scalar_1d.py` (invalidated), `m9_12_A2_dirac_2d.py`,
`m9_12_A2_dirac_3d.py`, `m9_12_audit_A2.py`.
