# M5.23 convo record

> Companion to [`m5_23_task_details.md`](m5_23_task_details.md) (the VIZ.5 "ellipsoid" visualization). Thread lineage: the viz spec first landed in [`m5_21_convo.md § 2026-07-18 19:44`](m5_21_convo.md) ("one value per 3D angle / one value per 2D angle").

## 2026-07-19: the author's ellipsoid-viz feedback (single ellipsoid per angle + the density warning + the split-vortex ask)

**Context**: the author's reply to a visualization share, forwarded mid-M5.23 run; it sharpens the render spec and names the follow-up animation arc.

Verbatim:

> Looks very nice, but to be more readable I would suggest to use single elipsoid per angle, like at the top of <https://community.wolfram.com/groups/-/m/t/3398814> - where I have emphasized single vortex (because split was too difficult for me - I have seen it some years ago, but somehow rejected mentally as too scary to include) - would be great to make similar animations for split vortices, also for muon and taon.
>
> I suspect its four 1/2-vortices are usually connected forming two loops - hence usually emitting two neutrinos in muon/taon decay ... would be great if making such animation, especially from simulation. Again for visibility it cannot be too dense: ellipsoid field only emphasizing central behavior and vortices.

Decode:

| Point | The author's words | Design consequence |
| --- | --- | --- |
| Single ellipsoid per angle + the author's reference render | "use single elipsoid per angle, like at the top of" the Wolfram post (the author's "Time crystal ϕ⁴ kink" Featured-Contributor post) | CONFIRMS the M5.23 shell design (one glyph per direction û); the reference look is individually distinguishable ellipsoids at moderate density, never a continuous field |
| The density warning | "for visibility it cannot be too dense: ellipsoid field only emphasizing central behavior and vortices" | Default direction count LOWERED to 162 (icosphere subdiv 2) with a density option (42 / 162 / 642); shells stay pinned to defect centers + vortex lines, never space-filling |
| Split vortices, muon / taon | "would be great to make similar animations for split vortices, also for muon and taon" | The Stage C vortex arm is reinforced; the μ/τ SPLIT-vortex animation from simulation is the follow-up arc riding the μ/τ decay task ([M5.21.6](m5_21_6_task_details.md)), not this run |
| The two-loop conjecture | "I suspect its four 1/2-vortices are usually connected forming two loops - hence usually emitting two neutrinos in muon/taon decay" | A pre-registerable read for the μ/τ program (loop-count on the simulated split-vortex state); logged here for [M5.21.6](m5_21_6_task_details.md) planning |

## 2026-07-19 evening: the VIZ.5 static-hedgehog video SENT

The rendered static biaxial hedgehog (the VIZ.5 ellipsoid shell + disclination rods delivered at [M5.23](m5_23_task_details.md), running on the [M5.24](m5_24_task_details.md) canonical stack) shared as a short video, <https://youtu.be/7LKbE_QYEd8>. Framing carried in the message: the dynamics stack is already ported and gated, the spinning view arrives with the fixed-J run, and no hand-drawn or keyframed geometry is involved (every shape, size and orientation comes from the simulated field tensor, the standing fidelity directive in [`m5_visualization.md`](../m5_visualization.md)).

## 2026-07-20 01:16 + 01:33: his reply to the video (the split-vortex render ask + the isosurface spec)

Verbatim (two messages, 17 minutes apart). One editorial substitution: the author names the specific AI assistant he consulted, written here as `[his AI]` (this record is public and the attribution is his to make, not ours):

> Thank you, looks very vice - I understand it is for this basic ansatz w with two 1-vortces.
>
> But [his AI] said, as in liquid crystal e.g. Alexander article, that both should split into two 1/2-vortices - would be great to visualize them, also for muon and taon.
>
> Maybe simplified as point with two loops - for muon/taon decay the central field should quickly rotate, releasing these two loops with large velocities - as neutrinos.

> ps. isosurfaces from these simulations would be great - visualized like e.g. <https://reference.wolfram.com/language/ref/ContourPlot3D.html> , or uniformly covered with ellipsoids.
>
> Especially for energy density - optimizing surface both around vortices, and central charge.

Decode:

| Point | The author's words | Routing |
| --- | --- | --- |
| The seed reading is correct | "it is for this basic ansatz with two 1-vortices" | ✅ Accurate: the video shows the SEEDED ansatz. The ½-split is a property of the RELAXED field, and the launcher has carried the FIRE relaxer since [M5.24](m5_24_task_details.md), so relax-then-trace shows the split in the same live tool |
| Visualize the ½-vortices | "both should split into two 1/2-vortices - would be great to visualize them" | Already MEASURED and receipted by the author ([Q31](../m5_question_tracker.md#q31-detail): TWO +½ through-lines, net +2 half-units on every plane at ρ ≥ 10, pair azimuth BRAIDING 72° → 169° → 135°); never RENDERED. This is exactly [M5.23.2](m5_23_2_task_details.md) arm (1), the tracer, ungated |
| Also for muon and taon | "also for muon and taon" | The converged μ/τ fields are ON DISK (the T2 pinned 32³ B and C endpoints from [M5.21.6](m5_21_6_task_details.md)); what is missing is a loader that seeds the launcher from a research endpoint: [M5.23.2](m5_23_2_task_details.md) arm (3) |
| Isosurfaces of energy density | "isosurfaces ... like ContourPlot3D, or uniformly covered with ellipsoids ... especially for energy density ... around vortices, and central charge" | Energy density has been live in production since [M5.24](m5_24_task_details.md) (true-zero vacuum floor) and the VIZ.5 mesh machinery is already in `medium.py`: [M5.23.2](m5_23_2_task_details.md) arm (4). His "covered with ellipsoids" variant = the existing ellipsoid samples moved from the S² shell onto the isosurface |
| Point with two loops, decay releasing them as neutrinos | "Maybe simplified as point with two loops - for muon/taon decay the central field should quickly rotate, releasing these two loops ... as neutrinos" | ⚠️ NOT drawn. The rotation-dominant transition IS measured at dynamics grade ([M5.21.6](m5_21_6_task_details.md) finding 7: ledger closes to 3 decimals, 23% of start energy radiated out), but the released count at descent grade was **ONE equatorial ring, not two loops**, and loop CLOSURE is unmeasured (in-box lines run boundary to boundary). A "simplified point with two loops" would be his conjecture drawn as a picture, which the standing directive forbids: we render what the tracer finds |
| The Alexander reference | "as in liquid crystal e.g. Alexander article" | AI-relayed by his side, so evidence-not-resolution per [`AI_HYGIENE.md`](../../../../../AI_HYGIENE.md); the underlying claim is independently ours (Q31 measured it before this message). The citation itself is unchased: chase it when the split arm writes up |
