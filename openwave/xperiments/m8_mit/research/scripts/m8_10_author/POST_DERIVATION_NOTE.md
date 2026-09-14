# M8.10: post-derivation observation (not scored, not used in the pre-registered computation)

Recorded 2026-09-12, after all eight values were known. By ruling it stays out of M8.10's scored claims, and it is to be derived before it is stated as more than an observation.

**Observation.** For each sector σ and each off-block level n = 2J, the level norm factorizes as

‖Πₙ N(Φ_{σ,v})‖² = C_{σ,n} · ‖[ρ₆(v) ⊗ v]_J‖²,

with C_{σ,n} independent of the ray. At the four rays the ratio is flat to about 10⁻¹⁵ (the diagnostic line of `m810_zeros.py`), including level 18 in sector 4, where the norm is taken over the whole two-copy multiplicity space. Values of C_{σ,n}: sector 3′, levels 10, 14, 16, 18: 1.318681×10⁻¹, 8.144796×10⁻², 5.042017×10⁻², 5.715647×10⁻²; sector 4, levels 8, 12, 14, 16, 18: 1.318681×10⁻¹, 1.318681×10⁻¹, 5.042017×10⁻², 8.144796×10⁻², 7.471167×10⁻². None has been identified exactly.

**Why it should hold (a sketch, not a proof).** By the paper's Lemma 4.1, R_K(P) = 0 for K = 1 to 5, so the density is a constant plus a rank-6 term, |Φ|² = c₀ + d₁₂. The constant multiplies Φ and stays at level 6, so every off-block level comes from d₁₂Φ. Its fibre factor couples ρ₆(v) with v, and spin 6 ⊗ spin 3 contains each J once, so that factor is [ρ₆(v) ⊗ v]_J up to a constant; the other factor depends on σ and n but not on v. A derivation also needs the Peter-Weyl structure of that other factor and, where the multiplicity is two, the multiplicity-space norm.

**What it would explain, and where it could lead.** The ray dependence of λ₄ would enter only through the fibre invariants ‖[ρ₆(v) ⊗ v]_J‖², and the sector dependence only through C_{σ,n}, which is why the sector ratios vary by ray instead of sharing one multiplier. If the fibre invariants and the constants are shown to be rational, the eight fractions become exact results rather than candidates. Blind agents who find the same structure without being told would make it a post-run result.

**Also recorded: a defect found and fixed before freeze.** `m810_exact.py` sorted 2I elements into conjugacy classes with a fixed tolerance of 10⁻⁴⁰, which silently finds no class below about 42 digits. The tolerance is now 10^−(dps/2); the classes are at least 0.19 apart. Runs at 60 and 100 digits were unaffected. The shape is worth naming: a fixed tolerance that fails only at a precision nobody had reason to run.
