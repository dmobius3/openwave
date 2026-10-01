"""Stage 2: rerun each room's run_all.py on a copy inside this room, then compare.

Copies room_a/ -> rerun_a/ and room_b/ -> rerun_b/ (scripts only; the outputs are
regenerated), runs `run_all.py` there with this interpreter (./py cannot be nested),
and writes the logs to rerun_a.log / rerun_b.log.  The comparison is done by
s2_compare.py.  Relative paths only.
"""
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = {
    "a": ["run_all.py", "t1_second_variation.py", "t3_spectrum.py", "t4_fem.py", "make_tables.py"],
    "b": ["run_all.py", "t1_t3_symbolic.py", "t3_t5_exact.py", "t4_fem.py", "collect.py"],
}
which = sys.argv[1:] or ["a", "b"]
for r in which:
    src, dst = os.path.join(HERE, f"room_{r}"), os.path.join(HERE, f"rerun_{r}")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    for f in SCRIPTS[r]:
        shutil.copy(os.path.join(src, f), os.path.join(dst, f))
    t0 = time.time()
    p = subprocess.run([sys.executable, "run_all.py"], cwd=dst, capture_output=True, text=True)
    with open(os.path.join(HERE, f"rerun_{r}.log"), "w") as fh:
        fh.write((p.stdout + p.stderr).replace(HERE + os.sep, ""))
    print(f"room {r}: run_all.py exit {p.returncode}, {time.time() - t0:.0f} s", flush=True)
