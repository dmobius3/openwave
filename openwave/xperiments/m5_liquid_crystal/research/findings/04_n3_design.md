# N3 design , the parameter+potential search (closed-loop field theory -> PMNS)

Go: 2026-06-21 21:27 EDT (01:27 UTC 06-22). Task: OpenWave #236 N3 (HELD, local). Resume ping NOT
armed (the maintainer: 11pm reset is the cap RELEASE not a deadline, won't hit it, one completion ping).

## The job (from 10a § "Sub-tasks" + the TBM gate)

Make the PMNS matrix `U` EMERGE from the closed-loop LC field theory, not be hardcoded (N2 hardcoded
`U_TBM` + a phenomenological two-level theta13). The TBM gate's literal bar: the 3 TBM angles come
"from the closed-loop FIELD DYNAMICS, not the group argument". Then derive theta13 (the SO(3)-breaking)
and resolve the central tension (theta13=8.5deg needs delta~0.15 OR a ~1.5e9 resonance, not delta~1e-10).

## The bridge (the central idea N3 rests on)

The three neutrino FLAVOUR states = three closed-loop configurations (Z3/S3-related). Define the
effective 3x3 mass matrix in FLAVOUR space from the loop field energetics:

```text
  M_ab = E_overlap[ loop_a , loop_b ]        (a,b in {e,mu,tau})
  M_aa = self-energy (the loop "mass")        (diagonal)
  M_ab = loop-loop coupling / field overlap   (off-diagonal, a != b)
  U    = eigenvectors of M_mass  (orthogonal; the columns = mass eigenstates)
  angles <- pmns_angles(U)   (reuse N2's PDG extraction)
```

This is the standard neutrino model-building structure made FIELD-THEORETIC: TBM <=> M_mass has the
"magic" (all row-sums equal -> (1,1,1)/sqrt3 is an eigenvector -> trimaximal solar col, sin^2 th12=1/3)
AND mu-tau (2<->3) symmetry. delta (the LC SO(3)-breaking eigenvalue) breaks them -> theta13 != 0.

The NOVELTY (publishable): show the closed-loop LC overlaps NATURALLY produce the magic+mu-tau matrix
(via a geometric Z3 arrangement of the three loops + mu-tau reflection), and that the LC delta sources
theta13 of the right size/sign.

## Plan (de-risk -> full search -> theta13/tension -> scorecard)

| Stage | File | What | Gate |
| --- | --- | --- | --- |
| D (de-risk) | `m5_11_n3_derisk.py` | (D1) verify U=eigvecs(M) bridge on a KNOWN TBM mass matrix round-trips the angles; (D2) democratic Z3 matrix -> trimaximal (sin^2 th12=1/3); (D3) magic+mu-tau -> exact TBM; (D4) which breakings turn on theta13, sign/size | scaffold correct before touching loops |
| S1 (machine) | `m5_11_n3_mass_matrix.py` | loop-overlap field theory: seed 3 flavour loops (Z3 geometric), compute M_mass from LC energy overlaps (N1 precision-safe when dressed), diagonalize -> U -> angles. Convention index-0. | M_mass real-symmetric, sensible |
| S2 (gate) | `m5_11_n3_search.py` | SEARCH (g, delta, LdG a/b/c, R, loop arrangement, core, dressing b*) for U -> TBM. Grid + refine, 16-core parallel. ★ TBM GATE | 3 TBM angles within tol from dynamics |
| S3 (crux) | `m5_11_n3_theta13.py` | turn on delta at the TBM point; measure theta13(delta); resolve tension (is eff delta ~1e-10 -> need resonance, or larger?); test the near-degenerate-gap resonance amplification | theta13 derived or honest gap reported |
| S4 (deliverable) | findings + plots | full scorecard vs NuFIT 6.0 (th12/th23/th13/delta_CP), document every param, reproducibility | peer-review-grade record |

## Honest FAIL path (gate clause)
If no (g,delta,potential,arrangement) reproduces TBM from the dynamics, that is a reportable finding
(closed-loop SO(3) dynamics don't yield TBM -> different structure needed) and redirects effort.

## Discipline
Headless numpy f64 (the f32 engine can't carry 1e10/1e-10; N0 proved the port is faithful). Checkpoint
each sub-result here on arrival. No GitHub, no git. Total invisibility (public OpenWave physics only).
Cross-link findings <-> 10a <-> #199.
