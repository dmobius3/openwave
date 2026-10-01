"""Run every stage-1 script in order.  Invoke as `./py run_all.py`; the scripts
run in-process (runpy) under that one ./py interpreter, because ./py cannot be
nested (its sandbox refuses a second sandbox-exec)."""
import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["fem2d_q2.py", "fem1d_sectors.py", "closed_form.py", "checks.py", "make_report.py"]

os.chdir(HERE)
sys.path.insert(0, HERE)
for s in SCRIPTS:
    print("=== %s" % s, flush=True)
    runpy.run_path(os.path.join(HERE, s), run_name="__main__")
