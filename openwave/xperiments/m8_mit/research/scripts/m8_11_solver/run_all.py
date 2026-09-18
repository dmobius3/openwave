"""Run every script in order.  Launch as `./py run_all.py`: the room interpreter cannot be nested
inside its own sandbox (sandbox-exec refuses), so each script is started with the same interpreter
and flags as ./py (sys.executable -S); the child inherits ./py's sandbox, environment (one thread,
PYTHONPATH) and priority.  Stops at the first script that exits with an error; reports the number
of PASS and FAIL lines."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ['s00_group.py', 's01_levels.py', 's02_su2_exact.py', 's03_stabilizers.py', 's04_pipeline.py',
           's09_zeros.py', 's10_existence.py', 's12_quadrature.py', 's13_identification_B.py',
           's14_negative_controls.py', 's99_assemble.py']
npass = nfail = 0
for s in SCRIPTS:
    t0 = time.time()
    print('=' * 20, s, flush=True)
    p = subprocess.run([sys.executable, '-S', s], cwd=HERE, capture_output=True, text=True)
    sys.stdout.write(p.stdout)
    sys.stderr.write(p.stderr)
    npass += sum(1 for l in p.stdout.splitlines() if l.startswith('PASS'))
    nfail += sum(1 for l in p.stdout.splitlines() if l.startswith('FAIL'))
    print('-- %s finished in %.0fs (exit %d)' % (s, time.time() - t0, p.returncode), flush=True)
    if p.returncode != 0:
        sys.exit('script %s failed' % s)
print('TOTAL: %d PASS, %d FAIL' % (npass, nfail))
