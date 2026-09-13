#!/bin/sh
# Reproduce every number of the stage-1 audit from a clean directory (outputs: *.json, *.log).
# Run from inside the room. Sequential on purpose (one single-threaded process at a time).
set -e
cd "$(dirname "$0")"
./py check_cg_D.py      > check_cg_D.log 2>&1         # CG (Racah vs lowering), D^j construction, product formula
./py item01.py          > item01.log 2>&1             # items 0 and 1
./py item5.py           > item5.log 2>&1              # item 5: invariants of V_K, R2 selection rule
./py algebra.py primary   > algebra_primary.log 2>&1  # items 2-8 exact, item 9 coefficient residual
./py algebra.py secondary > algebra_secondary.log 2>&1
./py functional.py primary   > functional_primary.log 2>&1   # CG-free explicit-function route + item 9 pointwise
./py functional.py secondary > functional_secondary.log 2>&1
./py item5b.py          > item5b.log 2>&1             # item 5: Jacobian and stabiliser arguments
# order: item5b before report (report reads item5b.json)
./py report.py          > report.log 2>&1             # cross-route comparison, item 8, audit_results.json
echo done
