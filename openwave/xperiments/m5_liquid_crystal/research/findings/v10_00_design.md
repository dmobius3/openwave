# N1+N2 foundation build , design + progress log (LOCAL, do not publish)

Go-time: 2026-06-21 19:00 EDT. Reset 11:00pm EDT (resume ping armed 03:05 UTC Jun 22).
HOLD: no GitHub #236 posting until the whole N-program (N1-N5) finishes (the maintainer's call,
`research/10a_neutrino_oscillations.md § N1 + N2 foundation scope` WORKFLOW NOTE). Everything LOCAL. Progress logged here.

Scope source: `research/10a_neutrino_oscillations.md § N1 + N2 foundation scope`. Master plan: `research/10a_neutrino_oscillations.md`.

## The two deliverables

| Phase | Deliverable | Script |
| --- | --- | --- |
| N1 | precision-safe numerical method (non-dim + perturbative-delta), Duda index-0 convention, cancellation test passes | `m5_11_n1_precision_method.py` |
| N2 | closed-vortex-loop seeder + 3 flavour states + SO(3) rotation + mixing observables (in the N1 formulation) | `m5_11_n2_closed_loop.py` |

## The precision problem (from Duda round-2, made concrete by reading the engine)

`M = O * diag(g, 1, delta, 0) * O^T`, `g ~ 1e10`, `delta ~ 1e-10`. The energy functional
(read from `engine3_observables.compute_energyH_density_M` + `engine2_pde.V_M`):

```text
H = 1/2 ||Mdot||^2  +  c^2 * 4 * sum_{mu<nu} ||[d_mu M, d_nu M]||^2  +  (V_M(M) - v0)
V_M = a*Tr(M^2) - b*Tr(M^3) + c*(Tr(M^2))^2      (on the NON-g block only)
```

Load-bearing facts verified from the engine:
- `V_M` acts on the NON-g 3x3 block ONLY (the g-axis is excluded: g^2~1e20 would inflate V).
  => the potential's delta-sensitivity sits at the `1 + delta^2` level (delta^2 lost vs 1 in f64).
- The CURVATURE includes the g-axis (full 4x4 commutators with the time axis). When the frame
  mixes the g-axis into the spatial structure (boost dressing), the commutators carry `g^2~1e20`
  AND `delta*g`, `delta^2` pieces in the SAME f64 sum. f64 ULP at 1e20 is ~1e4, so EVERY term
  below ~1e4 (i.e. the entire delta-sector: the SO(3)-breaking, the `theta_13` channel) is
  ANNIHILATED in a naive sum. This is the core error Duda flagged.

## N1 method (precision-safe)

1. Non-dimensionalization: measure energies in units of g (the fixed background scale) so the
   dynamical quantities are O(1) and g enters as a fixed multiplicative constant.
2. Perturbative-delta: `D(delta) = D_0 + delta*D_1` (linear), so `M(delta) = M_0 + delta*M_1`
   EXACTLY (`M_1 = O*D_1*O^T`). The polynomial energy expands exactly: `E(delta) = sum_k delta^k E_k`.
   `E_0` = the SO(3)-symmetric (TBM) structure; `E_1` = the leading SO(3)-BREAKING (-> theta_13).
   Each order computed in its own O(1) arithmetic, never summed against g^2.
3. Cancellation test: extract the O(delta) coefficient (the breaking) two ways:
   - NAIVE: [E(delta_h) - E(0)] / delta_h , both ~1e20 -> catastrophic, ~1e4 noise -> garbage.
   - PERTURBATIVE: the analytic directional derivative dE/ddelta, single pass, no big subtraction.
   - REFERENCE: complex-step derivative Im[E(i*delta_h)]/delta_h (no cancellation) + an mpmath
     arbitrary-precision check. PASS if perturbative matches the reference to ~machine precision
     while naive is wrong by many orders.
4. Convention: Duda index-0, `D = diag(g,1,delta,0)`, `eta = diag(-1,1,1,1)` (g + minus both index 0).
   (Engine uses the equivalent index-3 ordering; sandbox_v10 uses index-0 to match Duda.)

## N2 (closed loop)

- Closed disclination-loop seeder: circular loop radius R; the eigen-frame O(x) winds around the
  loop core. New geometry (engine has point hedgehogs + straight disclinations only).
- 3 flavour states = 3 loop/frame configs related by an SO(3) rotation in flavour space (per #199:
  oscillation = a single SO(3) rotation; TBM = a specific element).
- SO(3) rotation dynamics + mixing observables: the overlap (PMNS-like) matrix U between flavour
  and mass triads -> the standard angle extraction (sin^2 th13 = |U_e3|^2, etc.), reproducing #199's
  TBM angles as the delta=0 consistency check, and exposing the O(delta) -> theta_13 channel.
- Stability of a closed disclination loop = the honest engineering unknown (energy-relax + report).
- Does NOT claim any NuFIT value (that is N3 search + N4 theta_13).

## Sub-task status

| # | Sub-task | Status |
| --- | --- | --- |
| N1.1 | index-0 convention + energy functional in numpy f64 | 🔶 building |
| N1.2 | cancellation test (naive vs perturbative vs complex-step vs mpmath) | 🔶 building |
| N1.3 | non-dim + delta-expansion (M_0, M_1; E_0, E_1, E_2) | 🔶 building |
| N2.1 | closed-loop seeder (order-parameter field on a grid) | 🚧 next |
| N2.2 | 3 flavour states + SO(3) rotation + mixing observables (angle extraction) | 🚧 next |
| N2.3 | loop energy + relaxation (stability note) | 🚧 next |
| , | findings doc + plots + summary JSON | 🚧 next |

## Env

python 3.12.13 (openwave312), numpy 2.4.3, matplotlib 3.10.8 (Agg/headless), mpmath 1.3.0.
