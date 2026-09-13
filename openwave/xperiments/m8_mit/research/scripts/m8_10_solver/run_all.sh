#!/bin/sh
# Reproduce every number from a clean room directory. Run from inside the room: sh run_all.sh
set -e
cd "$(dirname "$0")"
./py item0_1.py        > log_item01.txt       # items 0, 1 (exact group theory) -> out_item01.json
./py engine.py 60      > log_engine_60.txt    # items 2-9, CG engine, 60 digits -> out_engine_60.json
./py engine.py 110     > log_engine_110.txt   # same, 110 digits
./py engine.py 160     > log_engine_160.txt   # same, 160 digits (item 6 orthogonal norm identification)
./py quadcheck.py      > log_quadcheck.txt    # independent float64 quadrature check -> out_quadcheck.json
./py relations.py 60   > log_relations_60.txt # item 2 relations among the M_K
./py relations.py 110  > log_relations_110.txt
./py exact_parts.py    > log_exact.txt        # exact symbolic route (items 2-5, 7, 8) -> out_exact.json
./py item5_zeros.py    > log_item5.txt        # item 5 arguments, stabilizer multiplicities
./py identify.py       > log_identify.txt     # two-precision identification + cross-checks -> results.json
echo done
