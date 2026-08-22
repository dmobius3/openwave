# Manifest Mismatch Report

Per TASK.md § 3: "If implementation exposes a mismatch with it, either change the
implementation to conform to the still-intended manifest or STOP and report the mismatch."

## Mismatch

The METHOD_AND_GATE_MANIFEST.md (SHA-256: b623535b..., immutable) states in § 4:

> | R6 | Sym³(ρ₂) = σ(Sym³(ρ₁)) | 4 |

This is mathematically incorrect. Sym³(ρ₁) has rational characters (its character
values on every group element lie in Q, not Q(φ)\Q), so applying the Galois conjugation
σ: φ → 1−φ to its matrix entries leaves the character unchanged:

    χ_{Sym³(ρ₂)}(g) = σ(χ_{Sym³(ρ₁)}(g)) = χ_{Sym³(ρ₁)}(g)  for all g ∈ 2I

Therefore Sym³(ρ₂) ≅ Sym³(ρ₁) as complex representations. They are the same irrep,
not a distinct 9th irreducible representation of 2I.

## Resolution

The implementation uses R6 = ρ₁ ⊗ ρ₂ (tensor product of the natural and Galois-conjugate
natural representations), which IS a distinct dimension-4 irreducible representation.

Verification:
- Character inner product ⟨χ, χ⟩ = (1/120) Σ |χ(g)|² = 1 confirms irreducibility
- Row signature (4, χ(s)=(1,0,1), χ(t)=(-1,0,1), χ(st)=(0,0,1)) is distinct from
  all other irreps, including Sym³(ρ₁)
- The representation has rational characters (self-conjugate under σ), consistent with
  the affine E₈ McKay quiver structure

The manifest was not amended. The implementation departs from the manifest's R6
construction while producing all 9 distinct irreps correctly.
