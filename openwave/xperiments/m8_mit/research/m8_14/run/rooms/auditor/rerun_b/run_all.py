"""Run every script in order; planted-defect runs are labelled.

Invoke as `./py run_all.py`. The children are started with sys.executable,
i.e. the interpreter that ./py already launched (calling ./py again from
inside ./py is refused by the room's sandbox).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
steps = [
    ["t1_t3_symbolic.py"],
    ["t1_t3_symbolic.py", "--plant"],
    ["t3_t5_exact.py"],
    ["t3_t5_exact.py", "--plant"],
    ["t4_fem.py", "--checks"],
    ["t4_fem.py", "--checks", "--plant"],
    ["t4_fem.py"],
    ["collect.py"],
]
for s in steps:
    print(f"\n===== {' '.join(s)} =====", flush=True)
    r = subprocess.run([PY] + s, cwd=HERE)
    if r.returncode != 0:
        sys.exit(r.returncode)
