# T1: MODELS.md criteria housekeeping: split the coarse rows, name the simplest passing test per row

**Status**: 🚧 STAGED 2026-07-27 (user decision, same day as the author proposal). Origin: the author's 2026-07-27 reply on the go-pack thread ([`m5_22_convo.md § 2026-07-27 afternoon`](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_convo.md)), two structural asks about the [`MODELS.md`](../../MODELS.md) score-board. **Expanded 2026-07-28**: the author posted a priority-ranked properties-to-test slide to the group ([`t1_convo.md`](t1_convo.md)), with the MODELS.md matrix embedded in the public talk deck; the per-row test list this task was staged to build is now author-proposed almost in full (§ The author's test map below). Platform-wide housekeeping: the criteria rows are shared by ALL model columns, so every edit here touches the whole matrix, not just M5. Docs-only, no runs.

## The two asks

| Ask | Content | Source |
| --- | --- | --- |
| Split the coarse criteria | Some rows bundle results of very different difficulty, so a single icon is dishonest in both directions. Named examples: **Weak force** contains muon decay (the simplest example, measured for M5) and beta decay (much harder, needs a proper neutron); **Baryons (p, n)** contains why the neutron is heavier, the positive core / negative shell, exact masses, and beta decay | author, 2026-07-27 (summarized; 1:1 channel) |
| Every row names its simplest passing experiment | "There should be specified some simplest e.g. experiment used to test given property." For "chiral SU(2)" the author's own question is what the test would even be. For **Magnetic force**: Coulomb + Lorentz covariance jointly imply magnetism, but **Lorentz covariance is itself missing as a criterion**; the closest direct test is **Larmor precession** (electron magnetic dipole in an external field); ortho- vs para-positronium (annihilation lifetimes, 2 vs 3 photons) is the harder follow-on | same |

## The author's test map (2026-07-28 group post)

The "Properties to test for SM + gravity?" slide ([`t1_convo.md`](t1_convo.md); local snapshot [`theory/duda_2026-07-28_properties_to_test_slide.png`](../../openwave/xperiments/m5_liquid_crystal/theory/duda_2026-07-28_properties_to_test_slide.png), deck [public](https://th.if.uj.edu.pl/~dudaj/AIphysics.pdf) + [local](../../openwave/xperiments/m5_liquid_crystal/theory/duda_2026-07-28_AIphysics_talk_slides.pdf)). Shading is the author's explicit priority code: **strong green = the most crucial property to test in a given type**, light green = secondary, no shading = long-horizon. Thin rules split the list into five blocks (charge/EM, spin/leptons, baryons/nuclei, SM parameters, gravity), mirroring the matrix sections. Transcribed from the emailed PNG; the deck snapshot regenerated 44 min later carries two line edits (noted in the convo record).

| Slide line (shading of each phrase) | MODELS.md home today | The named test, decoded |
| --- | --- | --- |
| **Charge quantization** (strong): why not e → e/2 + e/2 (the split highlighted as the thing to rule out) | Charge quantization | Show the topological charge cannot halve; the Gauss-law degree instrument ([M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) Q5 default) is the measurement |
| **Coulomb, Lorentz covariance** (strong) → magnetism | Electric force (Coulomb 1/r) + Magnetic force; Lorentz covariance row MISSING | Confirms the 2026-07-27 ask: covariance is its own criterion, and magnetism is graded as the derived consequence of Coulomb + covariance |
| running coupling (light), two linked papers (Universe, Phys. Rev. D) | none today | Candidate secondary row: running of the coupling with scale |
| ½ spin: same field by π rotation (light) | Spin-½ statistics (720° double cover) | The row's simplest test named: the field configuration returns to itself under a π spatial rotation |
| "time crystals": **(3) neutrino oscillations** (strong) | Neutrinos + de Broglie clock (Zitterbewegung) | The crucial neutrino test: 3 flavors + oscillations from the time-crystal mechanism |
| **electron: angular momentum** (strong) (→ pilot wave; deck edit: Klein-Gordon) | Angular momentum J (spin ℏ/2) | Field-carried angular momentum of the electron soliton, crucial tier |
| **magnetic dipole, Larmor precession** (strong) | Magnetic moment μ (g ≈ 2) + Magnetic force | Larmor precession confirmed as THE magnetic-force test (2026-07-27 ask, now priority-shaded) |
| **3 leptons** (strong), muon/taon decay → neutrinos (light) | Lepton mass spectrum (μ, τ) + Weak force | The 3-family spectrum is crucial; decays with neutrino release = the muon-decay sub-row of the Weak force split |
| para/ortho-positronium → 2/3 photons (light) | Antimatter + annihilation | The harder follow-on, confirmed secondary |
| **m_prot < m_neu** (strong), positive core/negative shell (light) | Baryons (p, n) | The mass-ordering + core/shell sub-row of the Baryons split = [M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) rungs 1-2 pre-registered reads |
| **beta decay** (strong), Cornell potential +~1 GeV/fm | Weak force + Strong force / confinement | The beta-decay sub-row; the ~1 GeV/fm Cornell scale anchors the [M5.21.4](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_21_4_task_details.md) string-tension instrument ([Q38](../../openwave/xperiments/m5_liquid_crystal/research/m5_question_tracker.md#q38-detail)) |
| **m_d < m_p + m_n** (strong), quadrupole electric mom. (light) | none today | Candidate NEW row: deuteron binding + electric quadrupole moment; already M5.22 rung 3 + pre-registered read (iv) |
| larger nuclei, levels, decays e.g. beta, halos (light) | none today | Candidate NEW row: nuclear structure; [M5.22.1](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_1_task_details.md) scope (the deck's halo-nuclei slide names Boron-8, Ne-17, Li-11, Borromean/Efimov targets) |
| pion, kaon, decays, Λ⁰ → p⁺ + π⁻ (light) | Mesons (π, K) | Strange-decay channels as the mesons test |
| masses, CKM, PMNS, SM parameters (no shading) | spread across rows | The long-horizon tier: derive the SM parameter set; no priority shading |
| **Newton** (strong) (→ GEM), light bend., time dilation (light) | Gravity | The gravity row's tests named: Newton limit via GEM crucial, light bending + time dilation secondary |

What the slide changes for the scope:

| Change | Content |
| --- | --- |
| The test column is author-proposed | Most rows now have their simplest test named by the author with a priority tier attached; the design work shifts from inventing tests to fitting them into the matrix (55-word budget, width) and marking the priority tier honestly |
| Confirmations of the 2026-07-27 asks | Lorentz covariance row, Larmor as the magnetic-force test, positronium as the harder follow-on, the Baryons mass-ordering/core-shell + beta-decay sub-rows: all reappear priority-shaded on the slide |
| Three candidate NEW rows beyond the 2026-07-27 asks | Running coupling (secondary); deuteron binding + quadrupole moment; larger nuclei / halos. The nuclear rows enter as 🚧 not yet tested until [M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) / [M5.22.1](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_1_task_details.md) produce evidence; adding a row with an honest 🚧 is in scope, upgrading an icon is not |
| The chiral-SU(2) question stands | The slide does not resolve which test the chiral-SU(2) row would get; the 2026-07-27 position (name a test or rephrase the row as testable) is unchanged |
| The stakes rose | The matrix is embedded in a public talk deck and the deck names OpenWave on its opening slide; score-board integrity now faces an external physics audience, and the author's open group question ("do you agree/disagree?") makes the T1 output the natural reply vehicle (reply user-gated) |

## Scope

| Piece | Content |
| --- | --- |
| Row-split proposal | A proposed new criteria set: at minimum Weak force → (muon decay \| neutron beta decay) and Baryons → (mass ordering + core/shell profile \| exact masses \| beta decay), plus a sweep of the other 20 rows for the same bundling smell (candidates on sight: Neutrinos, Gravity, Strong force). Each split row keeps a defensible icon PER MODEL from the evidence already linked, no new runs |
| The test column | A per-row "simplest passing test" entry (design decision: a new column vs a first line of the row label; weigh against the 55-word cell budget and the matrix width). Larmor precession lands as the magnetic-force test; the chiral-SU(2) row gets its test named or the row gets rephrased as testable |
| The Lorentz covariance row | New criterion row, missing today and load-bearing (it is what upgrades Coulomb into magnetism) |
| Score-board integrity | The summary counts (`Total criteria` = 22 today) and every per-model detail table re-derived to match the new row set; `python3 dev_docs/utils/check_models_md.py` green (55-word cells, icon/status sync) |
| Editing rules | The standing MODELS.md rules apply: cells are state of the art, no names/dates/endorsements in cells, verifiable refs only. The author-proposal provenance lives HERE and in the convo record, not in the matrix |

## What this task is not

No physics runs, no new evidence, no icon upgrades. Splitting a row may DOWNGRADE the apparent coverage of a model (a bundled ⚠️ can become one ✅ sub-row and one 🚧 sub-row); that is the point, and the score-board counts move accordingly. Any run the new rows motivate (Larmor precession on the 3×3/4×4 stack) is staged separately.

**Gated by**: user "go" (docs-only; can interleave between physics tasks). **Run-order promoted ahead of [M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) (user decision 2026-07-28)**: no dependency on census results (missing-evidence rows enter 🚧), and freezing the criteria set first means the census files into stable rows and [T2](t2_task_details.md) reads a stable N. Awaiting the go.

## Cross-links

| Doc | Why |
| --- | --- |
| [`m5_22_convo.md`](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_convo.md) | The 2026-07-27 proposal record + routing |
| [`t1_convo.md`](t1_convo.md) | The 2026-07-28 group post: the properties-to-test slide, verbatim + decode |
| [`MODELS.md`](../../MODELS.md) | The target doc |
| [M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) | The thread the 2026-07-27 proposal arrived on; its results will feed the split Baryons rows and the candidate nuclear rows |
| [`theory/_CITATIONS.md`](../../openwave/xperiments/m5_liquid_crystal/theory/_CITATIONS.md) | The deck citation (public URL + local snapshot) |
| [`platform_roadmap.md`](../platform_roadmap.md) | The live roadmap row (T1) + the migration change-log |
| [`m5_roadmap.md`](../../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) | Provenance change-log entries (staged as M5.29; input #2; the migration note) |
