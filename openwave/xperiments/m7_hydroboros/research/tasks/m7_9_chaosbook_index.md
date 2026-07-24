# M7.9 ChaosBook index (working set + where every benchmark value lives)

> E1 deliverable of [M7.9](m7_9_chaosbook.md). Source: Cvitanović, Artuso, Mainieri, Tanner, Vattay, *Chaos: Classical and Quantum* ("ChaosBook"), webbook at [chaosbook.org](https://chaosbook.org/), downloaded per chapter from the version 17 stable page ([version17/paper.shtml](https://chaosbook.org/version17/paper.shtml)); the PDFs self-identify as edition 16.0 text. Local copies: `theory/chaosbook/` (gitignored; re-fetch per chapter if absent). Citation entry: [`theory/_CITATIONS.md`](../../theory/_CITATIONS.md). The context rule: NEVER load the whole book (~1000 pages); this index maps the small working set, load one section or exercise at a time.

## 1. Working set (9 files, 3.4 MB, 174 pages total)

| File | Chapter | Pages | Why it is in the set |
| --- | --- | --- | --- |
| `contents.pdf` | front matter | 15 | the map of the tome; use to decide what else to fetch |
| `intro.pdf` | ch. 1 Overture | 35 | the author's suggested reading dose for the maintainer (E5 digest source) |
| `flows.pdf` | ch. 2 Go with the flow | 28 | Rössler flow defined (2.28); equilibria formula + values (2.29); exercises 2.7, 2.8 |
| `maps.pdf` | ch. 3 Discrete time dynamics | 18 | Poincaré sections + return maps; Hénon map (3.17)/(3.18); Lozi map (3.19) |
| `stability.pdf` | ch. 4 Local stability | 21 | stability matrix, exponents; Rössler equilibrium eigenvalues (4.36), example 4.5 |
| `invariants.pdf` | ch. 5 Cycle stability | 14 | Floquet multipliers; map-cycle Jacobian product (example 5.2); exercise 5.1 (analytic Floquet) |
| `cycles.pdf` | ch. 16 Fixed points, and how to get them | 13 | Newton for flows (16.3); multipoint shooting; exercise 16.9 = the Rössler cycle table to length 7 |
| `count.pdf` | ch. 18 Counting | 26 | topological trace formula; tables 18.1-18.3 (periodic points + prime-cycle counts) |
| `appendCycle.pdf` | appendix A16 Finding cycles | 4 | the multi-shooting Newton linear system (A16.1) implemented in the toolkit |

## 2. Where every benchmark value comes from (the hygiene ledger)

Every published number asserted in [`scripts/m7_9_benchmark.py`](../scripts/m7_9_benchmark.py) traces to a page in the working set; nothing is supplied from model memory (task doc § 4 E1 hygiene rule). Note: the book **truncates** (not rounds) its printed decimals, so comparison gates are 1 ULP of the printed last digit.

| Benchmark | Published value(s) | Source (file, page, label) |
| --- | --- | --- |
| B1 Rössler equilibria | formula `x± = (1/2 ± 1/2 sqrt(1−4ab/c²)) (c, −c/a, c/a)`; `(0.0070, −0.0351, 0.0351)`, `(5.6929, −28.464, 28.464)` | `flows.pdf` p. 24, eq. (2.29), example 2.3 |
| B2 equilibrium exponents | `(−5.686, 0.0970 ± i 0.9951)` inner; `(0.1929, −4.596e−6 ± i 5.428)` outer | `stability.pdf` p. 17, eq. (4.36), example 4.5 |
| B3 analytic limit-cycle Floquet | system `q̇ = p + q(1−q²−p²), ṗ = −q + p(1−q²−p²)`; the answer `{1, e^(−4π)}`, `T = 2π` is OUR hand derivation (the exercise asks for it; not printed in the working set) | `invariants.pdf` p. 14, exercise 5.1 (G. Bard Ermentrout) |
| B4 Rössler cycle table | 14 rows to topological length 7: itinerary, `(0, y_p, z_p)`, `Λ_e` (e.g. `1`: 6.091768, 1.299732, −2.403953) | `cycles.pdf` p. 12, exercise 16.9 table (J. Mathiesen, G. Simon, A. Basu) |
| B5 cycle counts | `N_n = 2, 4, 8, ..., 1024`; `M_n(2) = 2, 1, 2, 3, 6, 9, 18, 30, 56, 99`; `M_n(3)`, `M_n(4)` columns; formulas (18.26) | `count.pdf` p. 5 table 18.2, p. 14 table 18.3 |
| Rössler parameters | `a = b = 0.2, c = 5.7` | `flows.pdf` p. 24, eq. (2.28) |
| Hénon map | `x' = 1 − ax² + by, y' = x`, `(a, b) = (1.4, 0.3)` | `maps.pdf` p. 15, example 3.6, eq. (3.17) |
| Multi-shooting Newton system | block linear system (A16.1) | `appendCycle.pdf` p. 1, A16.1.2 |
| Map-cycle Floquet matrix | ordered product `M_p = Π_k J(x_k)` | `invariants.pdf` p. 13, example 5.2 |

Values NOT found in the working set (and therefore not benchmarked from the book): a Hénon periodic-orbit table (Hénon runs as an analytic pre-book gate instead, [`scripts/m7_9_gates.py`](../scripts/m7_9_gates.py) G2); the 3-disk pinball cycle tables (ch. 9 `billiards.pdf` not downloaded; the counting benchmark B5 uses the complete binary + N-alphabet tables 18.2/18.3, which the working set does print).

## 3. Candidate extensions (not needed for M7.9 DoD)

| Chapter | Would add |
| --- | --- |
| ch. 6 `Lyapunov.pdf` | Lyapunov exponents (exercise 6.4 evaluates them for Rössler) |
| ch. 9 `billiards.pdf` | the 3-disk pinball (exercise 16.8's analytic 10-cycle stability formula (16.14)) |
| ch. 14-15 `knead.pdf`, `smale.pdf` | symbolic dynamics, kneading theory (itinerary labels from first principles) |
| ch. 30 `PDEs.pdf` | Kuramoto-Sivashinsky: the bridge from small ODEs to field flows (closest to the M7.10-M7.14 use case) |

## 4. Cross-refs

[M7.9 task doc](m7_9_chaosbook.md) · [benchmark report](m7_9_benchmark_report.md) · [toolkit module](../scripts/m7_9_orbits.py) · [pre-book gates](../scripts/m7_9_gates.py) · [benchmark driver](../scripts/m7_9_benchmark.py) · [roadmap](../m7_roadmap.md) · consumer: [M7.10](m7_10_pure_maxwell.md)
