# M5.27 convo record

> Companion to [`m5_27_task_details.md`](m5_27_task_details.md) (the background-scalar time sector, closed 2026-07-24 with a structural null). This is an **author-channel record filed on the M5.27 ID**, not a record of the M5.27 run: the exchange is public, it carries the author's own statement of the three approximation levels and the `λ` / `Λ` / `O` notation convention, and its level-3 clause bears directly on the M5.27 verdict (§ Bearing on the time sector below).

## 2026-07-25: the three-level ladder + the λ / Λ / O notation convention (public, GitHub Discussions [#334](https://github.com/openwave-labs/openwave/discussions/334))

**Context**: the M8 (Mode Identity Theory) column author opened a Q&A discussion asking whether M5's spatial frame should be identified (soldered) with the MIT adjoint holonomy `Ad(σ): 2I → SO(3)`, since `M`'s transformation law `M → Λ⁻ᵀ M Λ⁻¹` is a frame law rather than a separate internal gauge index. Both branches were pre-registered in the opening post (solder it under a separately-named track, or leave it native), and the no-reply default was declared: native and untwisted, three-connection compatibility recorded "not applicable" as a neutral status, per [`dev_docs/CROSS_MODEL_TESTING.md`](../../../../../dev_docs/CROSS_MODEL_TESTING.md) § 1 and § 4. The maintainer side posted nothing substantive: the question is author-gated (§ 6 routing), and there is no implementation layer we hold facts about. The author replied the same night with notation rather than a yes/no.

Verbatim:

> Hi, seeing "M -> Lambda^-T M Lambda^-1" we might have notation differences, so let me briefly elaborate - liquid crystal approach allows to work in 2 approximated modes (1,2) and complete (3):
>
> 1. EM only: S^2 vacuum, Faber, uniaxial/Oseen-Frank: vector field n with Higgs-like potential preferring unitary n, getting EM with topological charge quantization and running coupling effect (e.g. <https://journals.aps.org/prd/abstract/10.1103/8zn4-rwth> ),
> 2. EM+QM: SO(3) vacuum, 3D biaxial/Landau-de Gennes: 3x3 real symmetric matrix field M with potential preferring Lambda=(1,delta,0) eigenspectrum of M for delta~10^-10, allowing to recognize low energy twists of Faber's unitary vectors n, interpreted as quantum phase (getting Klein-Gordon-like equation),
> 3. EM+QM+GEM: SO(1,3) vacuum, 4D liquid crystal/Einstein's teleparallelism: 4x4 real symmetric tensor field, with potential preferring Lambda=(g,1,delta,0) eigenspectrum of M for g~10^10, adding boosts getting 2nd set of Maxwell equations for gravitoelectromagnetism (GEM), also making that energy minimization prefers nonzero time derivatives: for e.g. electron's angular momentum, neutrino oscillations.
>
> So personally I denote lambda as current M eigenspectrum, Lambda as (vacuum) eigenspectrum preferred by potential, and O as SO(3) rotation in level 1,2, and SO(1,3) with boosts in level 3.
>
> You can use M5 as it suits you best (not modifying it), but please be careful regarding notation differences.

### The ladder, level by level

The hierarchy itself is not new to the program (the EM ≫ QM ≫ GEM ordering and the teleparallel 4D extension were recorded from the 2026-05-12 call, [`m5_4a_convo_2026.05.12.md`](m5_4a_convo_2026.05.12.md), and the ladder is taught in [`__M5_course.md`](../../__M5_course.md) L3). What is new is that all three levels are now stated together, in the author's voice, on a public record, each with its vacuum manifold and its named liquid-crystal framework.

| Level | Vacuum + framework | Where our work sits |
| --- | --- | --- |
| 1, EM only | `S²` vacuum, uniaxial / Frank-Oseen, vector field `n` with a Higgs-like potential; gives EM with topological charge quantization + the running-coupling effect | the M5.6 Faber port (`Λ = q₀⁶/r₀⁴`, mass pinned `E ∝ 1/r₀`) and the M5.17 two-charge running-coupling read live here |
| 2, EM+QM | `SO(3)` vacuum, 3D biaxial Landau-de Gennes, 3×3 real symmetric `M`, `Λ = (1, δ, 0)` at `δ ~ 1e-10`; low-energy twists of Faber's `n` read as quantum phase, giving a Klein-Gordon-like equation | the certified 3D statics instrument ([`m5_21_2b_note.md`](../findings/m5_21_2b_note.md)) and the whole lepton census. Note `O ∈ SO(3)` here is the **spatial** rotation group, so `M`'s indices are frame indices, not a separate internal gauge index: the reason the soldering question partly dissolves rather than resolving yes or no |
| 3, EM+QM+GEM | `SO(1,3)` vacuum, 4D liquid crystal / Einstein teleparallelism, 4×4 real symmetric tensor, `Λ = (g, 1, δ, 0)` at `g ~ 1e10`; boosts give the second set of Maxwell equations (GEM), and energy minimization "prefers nonzero time derivatives" | the 4×4 lift ([`m5_21_3_note.md`](../findings/m5_21_3_note.md), [`m5_27_note.md`](../findings/m5_27_note.md)). The closing clause is the time-sector claim the program has been measuring against: see below |

### The notation convention

⚠️ **The collision that prompted the reply.** The opening post wrote `Λ` for the Lorentz transform acting on `M` (`M → Λ⁻ᵀ M Λ⁻¹`); the author writes `Λ` for the vacuum eigenspectrum preferred by the potential. Same letter, different object, which is what "we might have notation differences" points at. In the author's convention the transform is `O`.

| Symbol | The author's convention | Our docs |
| --- | --- | --- |
| `λ` | the **current** eigenspectrum of `M` | the eigenvalue reads (the `λ₂ − λ₃` uniaxial-escape split, the biaxial census) |
| `Λ` | the **vacuum** eigenspectrum preferred by the potential | written `D` in [`__M5_model_briefing.md`](../../__M5_model_briefing.md) and as the `M_vac` spectrum in [`m5_theory_canonical.md`](../m5_theory_canonical.md): same object, different letter |
| `O` | the rotation: `SO(3)` at levels 1-2, `SO(1,3)` with boosts at level 3 | `O` in `M = O·D·Oᵀ`: same symbol, same meaning |

⚠️ **The η caution, which is where a real notation difference bites.** The author writes `Λ` as "eigenspectrum of M" with a **positive** leading `g`. Our implemented vacuum **tensor** is `M_vac = diag(−g, 1, δ, 0)` with the time-time entry NEGATIVE, and `diag(+g, 1, δ, 0)` is **not** a vacuum (`V ≈ 1.05e6` at g = 8, machine-exact, [`m5_18_verification_note.md`](../findings/m5_18_verification_note.md) § 5), because spectrum statements in this program always mean the spectrum of `η·M` with `η = diag(−1, 1, 1, 1)`. Anyone reading `D` or `Λ` as the literal diagonal of `M` and building a background from it writes the wrong sign at `[0,0]`. This entry does **not** resolve which sign is physical: that is [Q27](../m5_question_tracker.md), answered-as-open on 2026-07-17 ("I am not certain if g should be positive or negative"), so the theory does not select the branch and the branch-qualification policy stands. The exchange is the reason the briefing's Substrate and Vacuum rows now carry the η convention explicitly.

The physical scales quoted here (`g ~ 1e10`, `δ ~ 1e-10`) match the locked-era values already on record; the dynamics-era working values are `g = 8`, `δ = 0.3`, bridged by the measured scaling laws (transverse split ∝ `δ^1.034`, vacuum stiff mode ∝ `g^2.992`), per [`m5_theory_canonical.md`](../m5_theory_canonical.md).

### Bearing on the time sector (why this is filed on the M5.27 ID)

| The level-3 closing clause | Status against measurement |
| --- | --- |
| "energy minimization prefers nonzero time derivatives: for e.g. electron's angular momentum, neutrino oscillations" | The claim is re-asserted, the mechanism is not supplied, and the measured gap is unchanged. [M5.21.8](m5_21_8_task_details.md) found **no finite nonzero ω in the author's own formulas** (`ω* = 0` or `−∞`; the `Hm` `ω²` coefficient is negative for all `g > 1`, `δ ≤ 1.9`, including `(1e10, 1e-10)` exactly), and the lattice agreed. [M5.27](m5_27_task_details.md) then measured a **structural null** for the background-scalar route: the drive owns eigenVALUES and exerts zero torque on the eigenFRAME that carries the clock. The live successor lane is a NON-COMMUTING coupling (the mixed `(0,i)` block or derivative/frame couplings), which is Lagrangian-level work, carried by the [M5.22](m5_22_task_details.md) hand-off |

Read together with the standing time-sector record: the fixed-J isorotation state remains the clock of record, and nothing in this entry changes that.

### Doc consequences (applied 2026-07-25)

| Change | Where |
| --- | --- |
| The η convention stated on the two rows a cross-model reader would build a background from | [`__M5_model_briefing.md`](../../__M5_model_briefing.md) Substrate + Vacuum rows, pointing at the canonical row and the M5.18 measurement |
| The level-1 citation the author linked is a **new publication**: Faber & Golubich, *Phys. Rev. D* **114**, 014510, published 2026-07-16, [DOI 10.1103/8zn4-rwth](https://doi.org/10.1103/8zn4-rwth), the peer-reviewed version of arXiv:2604.12021 (the arXiv record carries no journal-ref yet). Citation sweep applied across both models | [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md), [`__M5_model_briefing.md`](../../__M5_model_briefing.md), M7 briefing + citations + `m7_0_bootstrap.md`, repo [`README.md`](../../../../../README.md) |
| Two citation errors found in the same sweep and fixed: `Universe` **11**(4):113 is Wabnig, Resch, Theuerkauf, Anmasser & Faber (previously miscredited here as "Faber & Golubich" and as "Manfried Faber et al."), and the stale "Duda 2026?, running coupling, MDPI 2076-3417" entry was a mis-transcription resolved to Faber's two papers | [`m5_5a_lagrangian_evolution.md`](m5_5a_lagrangian_evolution.md), [`m5_1a_lagrangian_framework.md`](m5_1a_lagrangian_framework.md) |
| `faber_universe_2025.pdf` is a filename misnomer: its content is the arXiv:2604.12021v1 preprint of the dipole paper, not a Universe 2025 article (title page read). The phantom "Faber universe note (author draft)" bibliography row was removed and the local-corpus mapping made explicit | [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md) |

### Cross-model outcome

No M5-side change was requested or made. The reply grants the § 1 native default plus permission for a separately-named track ("You can use M5 as it suits you best (not modifying it)"), which is exactly the branch the opening post pre-registered as the no-reply default. The companion M4 thread ([#333](https://github.com/openwave-labs/openwave/discussions/333)) closed the same way: native M4 untwisted, three-connection rows "not applicable", the model retained as a dynamical comparison column.
