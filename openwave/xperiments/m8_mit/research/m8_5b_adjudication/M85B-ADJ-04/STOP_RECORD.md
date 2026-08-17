# M85B-ADJ-04: STOP / PRE-REVEAL ROUTE-A STRUCTURAL FAILURE

    Case        M85B-ADJ-04, committed under Addendum 12.2
    Outcome     retired as a pre-reveal bug-discovery case under
                Addendum 12.3. No adjudication occurred
    Packet II   sealed and unread throughout; its commitment stands

## What happened

Packet I was opened at § 4.1 step 2 and its bytes verified against the
adopted commitment. During step 3, route (a) refused before producing an
output: its multiplication-table lookup required sign-exact equality of
action-pair representatives, while the effective-group closure had chosen
representatives modulo the diagonal central kernel `(u, v) ~ (-u, -v)`, so
products landing on a central partner were reported as absent. Route (b)
was not run, no route output was committed, and no external reference value
was consulted or compared. The external test is therefore not scored: no
prediction ever met the reference.

The defect is latent, not case-induced: the § 6.1 tuning case `L(7;1,2)`
fails by the identical mechanism on the same path. It survived
qualification because every case that previously reached route (a)'s
stencil path was sign-exact closed by accident of selection. Addendum 12.3
records the repair authorization and the requalification obligations;
`route_a_closure_battery.py` in the qualified tree now covers the missing
regime, with a mutation arm demonstrating the pre-repair lookup fails.

## Files here

`run_step3.py` is the step-3 driver exactly as executed at the STOP. It
produced no artifacts. This directory deliberately contains no route
output and no packet bytes.
