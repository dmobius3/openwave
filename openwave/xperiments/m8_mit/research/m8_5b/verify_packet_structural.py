#!/usr/bin/env python3
"""Packet-schema structural qualification: coverage and non-vacuity.

Establishes three things about the Packet-I structural surface and the
rung-3a/3b index-transform surface, each printed with its own evidence:

    A  the structural battery, scored by TARGET predicate, 8 of 8 covered
    B  non-vacuity: suppress each predicate in turn, its own item must red
    C  the S3 full-class battery over both parities of group order

Exit 0 only if every section passes.  Invoked by `run_qualification.py`; also
runnable on its own from the tree root:

    python3 verify_packet_structural.py
"""

import sys

import packet_schema as ps

PREDICATES = ("S0", "S1", "S2", "S3", "S4", "S5a", "S5b", "S6")
ok = True


def head(t):
    print(f"\n{'=' * 72}\n{t}\n{'=' * 72}")


# --- C ----------------------------------------------------------------------
head("C. the new structural battery, scored by TARGET predicate")
c_ok, rows = ps.structural_battery(verbose=True)
print(f"\n   {sum(r['pass'] for r in rows)}/{len(rows)} PASS")
ok &= c_ok

# --- D ----------------------------------------------------------------------
head("D. non-vacuity: suppress each predicate, its own item must red")
orig = ps.structural_checks
for tgt in PREDICATES:
    ps.structural_checks = (lambda p, t=tgt:
                            [h for h in orig(p) if not h.startswith(t + " ")])
    d_ok, d_rows = ps.structural_battery(verbose=False)
    red = sorted({r["item"].split(":")[0] for r in d_rows
                  if not r["pass"] and ":" in r["item"]})
    own = tgt in red
    print(f"   suppress {tgt:4} -> battery {'RED ' if not d_ok else 'PASS'}   "
          f"items red: {red}   own item red: {own}")
    ok &= (not d_ok) and own
ps.structural_checks = orig
r_ok, r_rows = ps.structural_battery(verbose=False)
print(f"   restored        -> {sum(r['pass'] for r in r_rows)}/{len(r_rows)} PASS")
ok &= r_ok

# --- C ----------------------------------------------------------------------
head("C. S3 full-class battery, both parities of group order")
s_ok, s_rows, touched = ps.s3_battery(verbose=True)
print(f"\n   {sum(r['pass'] for r in s_rows)}/{len(s_rows)} PASS; "
      f"{len(touched)} fixture cases swept (diagnostic only, constrains nothing)")
ok &= s_ok

head("VERDICT")
print(f"   {'ALL SECTIONS PASS' if ok else 'FAILURE, see above'}")
sys.exit(0 if ok else 1)
