# M9.10--M9.11: A2 then A1 (1d free fermion)

> Pre-registered. A2 is a necessary-condition test, not a 4d proof.

## A2 verdict

`A2_PASS` on the locked gate: \(R(m)/R(0)<2\) for \(0<mL\le 8\).
Solver ratios 1.07, 1.16, 1.32. Auditor (N=180, L=24) 1.26, 1.69.
Mutation C3 fires (diagonal mass visible). Remainder *grows* with
\(m\); it does not cross the locked factor of 2.

Scope: 1d staggered-mass fermion, interval. Not a 4d diamond. Not
the SM. The local-\(X\) ansatz is **not refuted** here.

## A1 verdict

`A1_PASS`. \(\alpha(0)=0.323\approx 1/3\). UV relative drift
\(\le 0.072\). IR mutation \(\alpha\to 0.010\). Auditor
\(\alpha(0)=0.325\), UV rel \(0.032\).

Scope: 1d leading log. Not a 4d area law. Not \(\eta=1/4G\).

## Tags

| Claim | Tag |
| --- | --- |
| A2 C2 on this lattice fermion | *computed* (threshold 2, locked) |
| A2 in 4d / for the SM | *unresolved* |
| A1 UV \(\alpha\) IR-independent, this field | *computed* |
| Jacobson 2016 nonconformal in 4d | still not `[P]` |

Scripts: `m9_10_A2_modular.py`, `m9_10_audit_A2.py`,
`m9_11_A1_uv.py`, `m9_11_audit_A1.py`.
