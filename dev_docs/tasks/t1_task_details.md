# T1: MODELS.md criteria housekeeping: split the coarse rows, name the simplest passing test per row

**Status**: ✅ DONE 2026-07-28 (go 10:26 EDT, review approved same morning; staged 🚧 2026-07-27, user decision, same day as the author proposal). Origin: the author's 2026-07-27 reply on the go-pack thread ([`m5_22_convo.md § 2026-07-27 afternoon`](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_convo.md)), two structural asks about the [`MODELS.md`](../../MODELS.md) score-board. **Expanded 2026-07-28**: the author posted a priority-ranked properties-to-test slide to the group ([`t1_convo.md`](t1_convo.md)), with the MODELS.md matrix embedded in the public talk deck; the per-row test list this task was staged to build is now author-proposed almost in full (§ The author's test map below). Platform-wide housekeeping: the criteria rows are shared by ALL model columns, so every edit here touches the whole matrix, not just M5. Docs-only, no runs.

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

## DEVIATIONS LOG (2026-07-28 run)

| Deviation | What and why |
| --- | --- |
| Beta decay gets ONE home | The staging scope listed beta decay in BOTH splits (Weak and Baryons). One experiment = one row, or the score-board double-counts it; it lands under Weak force (the author's 2026-07-27 framing: the hard sub-row that needs a proper neutron), and the Baryons rows do not repeat it |
| Baryons split 3-way, not 2-way | The author's decomposition (ordering, core/shell, exact masses) leaves no honest home for M4's existing evidence (the K = 10 tetrahedron is a bound-state existence result, neither ordering nor masses). A `bound state` sub-row preserves it; the author's two items become the other two rows |
| Priority marking kept OUT of the matrix | The slide carries three tiers (strong green / light green / unshaded). A bold crucial-tier marking was tried at execute and removed at review (user call): the test column stays plain, and the full tier decode lives in § The author's test map above |
| Running coupling = one row, not two | The slide names running coupling in the charge/EM block AND the M5 strong-force cell carried a running-onset result. One criterion row (`Running coupling`, FORCES) serves both; the strong row narrows to confinement |
| Neutrinos + Gravity split by the sweep | The staged sweep candidates confirmed: Neutrinos → neutral states \| oscillations (PMNS), Gravity → Newton limit (GEM) \| metric phenomena. Strong/confinement handled via the running-coupling extraction. M8's Λ evidence lands on the metric row, M5's boost-tilt evidence on the Newton row |
| Sawada target relocated | M6's Sawada long-range nuclear anomaly (v(r) ~ −C/r⁶) is an internucleon-potential falsifier: it moves from the old Strong-force cell to the new Nuclear-structure row |

## FINDINGS (2026-07-28)

The criteria set went **22 → 31 rows**, every count re-derived, [`check_models_md.py`](../utils/check_models_md.py) green (155 cells, first try).

| Deliverable | Landed as |
| --- | --- |
| Row splits | Weak force → muon decay ⚠️(M5) \| beta decay 🚧; Baryons → bound state \| mass ordering + charge profile \| exact masses; Neutrinos → neutral states \| oscillations (PMNS); Gravity → Newton limit (GEM) \| metric phenomena; Strong → confinement \| Running coupling |
| The test table | Every one of the 31 criteria names its cheapest passing experiment in a dedicated two-column companion table right below the matrix (`§ Simplest Test per Criterion`), keeping the matrix itself to icons; the linter enforces non-empty tests and a criteria set synced with the matrix, both directions. Cells stay unmarked; the priority-tier decode lives in § The author's test map |
| Lorentz covariance row | Added to FORCES: M5 ⚠️ (c tilt modes + relativistic KG dispersion = linearized sector measured; boosted defects untested), M7 ⚠️ (exact KG branches, harmonic-frame + tachyon caveat), M4/M6/M8 🚧 |
| Chiral SU(2) | Resolved by test-naming rather than a row: the beta-decay row's test is `n → p + e + ν̄, parity-violating`, parity violation being the experimental face of the chiral question |
| Candidate NEW rows | All three entered 🚧 across every model: Running coupling (M5 ⚠️: onset at r₀ measured), Deuteron (binding + quadrupole; M5.22 rung 3), Nuclear structure (levels, halos; M5.22.1 scope, M6's Sawada target) |
| Score-board | M5 9✅ 10⚠️ 1❌ 11🚧 · M7 0/10/0/21 · M4 0/8/3/20 · M6 3/3/3/22 · M8 0/1/0/30; column order unchanged (✅+⚠️: 19/10/8/6/1) |
| Honest downgrades (the point of the split) | M5's bundled Weak ⚠️ now shows beta decay 🚧; M5's bundled Gravity ⚠️ now shows metric phenomena 🚧; M5's bundled Baryons 🚧 stays 🚧 on all three sub-rows |
| Linter | Docstring now six checks; `simplest test` registered as a named criterion-level column (the unknown-column guard still fires for anything unregistered) and an empty test cell is a violation |

Icon provenance: no icon was upgraded anywhere; every ⚠️/✅/❌ on a split or new row is carried by evidence already linked in the pre-split cell (verifiable via git diff of this commit).

## TASK REVIEW (2026-07-28)

Task Duration: 00:15 (from 10:26 to 10:41 EDT)
Usage Cap Triggered: NO

Approved by the user same morning, with two adjustments at review: the bold crucial-tier marking was removed from the test cells (plain text; the tier decode lives in this doc only), and the test column moved OUT of the summary-status table into a dedicated companion table below it (`§ Simplest Test per Criterion`), so the matrix stays lean; the linter now parses that table and syncs its criteria set with the matrix both directions (mutation-tested: empty cell, missing criterion, orphan row all fail). Results: all deliverables ✅ (31-row criteria set, the simplest-test table, Lorentz covariance row, three candidate rows entered 🚧, score-board re-derived, no icon upgraded; `check_models_md.py` green first try at 155 cells, doc checker + roadmap linter clean). Issues: none blocking; one judgement call flagged (M7's Lorentz ⚠️ rests on the measured exact KG branches while the harmonic frame fixes ω by construction; 🚧 also defensible under a stricter read). Deviations: six, logged above. Follow-ons user-gated: a short reply to the author announcing the update was drafted at review (user-voice, user sends; drafted text stays out of repo docs) and the deck's "21 criterions" vs the live 31 is left to the author.

**Findings**: splitting the coarse rows surfaced three hidden 🚧 inside M5's former ⚠️ bundles (beta decay, metric gravity) while every model kept the evidence it had actually earned; all 31 rows now name the cheapest experiment that would settle them, with parity violation answering the chiral-SU(2) test question.

**Research docs created/updated**: this file (deviations log, findings, review) · [`MODELS.md`](../../MODELS.md) (the restructured matrix) · [`check_models_md.py`](../utils/check_models_md.py) (6th check) · [`platform_roadmap.md`](../platform_roadmap.md) (T1 → Done).

## Post-close follow-up (2026-07-29)

The author replied 1:1 with one question on the new Lorentz covariance row: why test covariance in solutions rather than only in equations. The exchange record (the question, the sent reply, and what the answer decides) is in [`t1_convo.md § 2026-07-29`](t1_convo.md#2026-07-29-0354-edt-11-the-lorentz-covariance-question). No task action until the answer arrives; the possible edit (reword the simplest-test row, or split the criterion into equation-level and boosted-solution rows) is user-gated.

## Cross-links

| Doc | Why |
| --- | --- |
| [`m5_22_convo.md`](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_convo.md) | The 2026-07-27 proposal record + routing |
| [`t1_convo.md`](t1_convo.md) | The thread record: the 2026-07-28 slide post (verbatim + decode) + the 2026-07-29 Lorentz covariance follow-up (question, reply summary, pending decision) |
| [`MODELS.md`](../../MODELS.md) | The target doc |
| [M5.22](../../openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_task_details.md) | The thread the 2026-07-27 proposal arrived on; its results will feed the split Baryons rows and the candidate nuclear rows |
| [`theory/_CITATIONS.md`](../../openwave/xperiments/m5_liquid_crystal/theory/_CITATIONS.md) | The deck citation (public URL + local snapshot) |
| [`platform_roadmap.md`](../platform_roadmap.md) | The live roadmap row (T1) + the migration change-log |
| [`m5_roadmap.md`](../../openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) | Provenance change-log entries (staged as M5.29; input #2; the migration note) |
