# Conventions, and a worklist

Everything below happens in one seven-dimensional complex space. No quotient, no sections and no bundles are needed, and none are supplied.

## 1. Setup

### 1.1 The space

`V₃ = ℂ⁷` with basis `v₃, v₂, v₁, v₀, v₋₁, v₋₂, v₋₃`, written `u = Σ_m c_m v_m`. The Hermitian product is `⟨u, w⟩ = Σ_m conj(c_m) d_m`, and the real inner product is `Re⟨u, w⟩`.

### 1.2 The rotation action

`J_z v_m = m·v_m`, and `J_± v_m = √(12 − m(m ± 1))·v_{m±1}`, with `J_x = (J_+ + J_−)/2` and `J_y = (J_+ − J_−)/2i`. A rotation by angle `θ` about a unit axis `n` acts by `D³(n, θ) = exp(−iθ·(n_x J_x + n_y J_y + n_z J_z))`. This is a representation of `SO(3)`: the rotation by `2π` acts as the identity.

### 1.3 Time reversal and the quartic

`(Θu)_m = (−1)^{3−m}·conj(c_{−m})`.

With the Clebsch-Gordan coefficients `⟨3 m₁; 3 m₂ | 6 Q⟩` in the Condon-Shortley convention, fixed by `⟨3 3; 3 3 | 6 6⟩ = +1`, define, for `Q = −6, …, 6`,

`ρ₆(u)_Q = Σ_{m₁} ⟨3 m₁; 3 (Q − m₁) | 6 Q⟩ · c_{m₁} · (Θu)_{Q − m₁}`,

and then

`r̂₆(u) = Σ_Q |ρ₆(u)_Q|² / ‖u‖⁴`.

`r̂₆` is real, invariant under `u ↦ e^{iφ}u` and under every rotation, and homogeneous of degree 0.

## 2. Worklist

Answer in order. Where an item asks for an argument, give it as prose with every hypothesis listed, and for each step say what would fail if that step were omitted. A numerical search or scan is welcome as corroboration, and you must say what it can and cannot establish; on its own it establishes neither a maximum, nor a minimum, nor that a set of extremizers is complete.

**0.** Report `⟨3 3; 3 −3 | 6 0⟩`, and `Θ(Θu)` in terms of `u`, and say what precision your route carries.

**1.** Find the maximum of `r̂₆` over the unit sphere of `V₃`, exactly, with an argument.

**2.** Find the complete set of unit vectors at which that maximum is attained, exactly, and give an argument that the set you report is complete. Describe the structure of that argument: the cases or branches it splits into and how they cover everything. For each case, say what would fail if it were omitted. If any step rests on a solver, say which, and say what would fail if you relied on it alone.

**3.** Answer items 1 and 2 for the minimum of `r̂₆` over the unit sphere.

**4.** For every quantity in items 1 to 3 that vanishes or is found equal to a bound, say whether it was computed and found so. A value your computation did not reach is reported as missing, never as zero. For each computed equality, give an exact reason.

## 3. What to return

One markdown file with:

- the maximum and the minimum, exactly, each with its argument and hypotheses;
- the complete set of maximizers and of minimizers, exactly, each with its completeness argument and the "what would fail if omitted" answer for every case;
- for every numerical route, the precision it carries and what it can and cannot establish;
- the tools you used, and anything you looked up rather than derived.

Report disagreement with your own earlier steps rather than smoothing it. If an item cannot be completed, say which and why; an incomplete item reported as incomplete is worth more than a filled one that is not verified.
