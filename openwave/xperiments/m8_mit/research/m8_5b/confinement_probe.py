"""Subprocess confinement probe.

`run_qualification.py` launches this THROUGH its `run()` helper, the same entry
point every battery uses, so what it reports is the environment the batteries
actually receive rather than one reconstructed alongside them. It reports its
own view and a grandchild's, because the confinement claim covers nested
children and no parent-side inspection can reach them.

The caller deliberately poisons `PYTHONPATH` and `PYTHONNOUSERSITE` before
launching this. Without that, the probe could not fail: with nothing hostile to
inherit, dropping the confinement would look exactly like enforcing it.
"""
import json
import os
import subprocess
import sys

MARKER = "CONFINEMENT_PROBE_JSON "

_GRANDCHILD = (
    "import json, os, sys;"
    "print(json.dumps({"
    "'PYTHONPATH': os.environ.get('PYTHONPATH'),"
    "'PYTHONNOUSERSITE': os.environ.get('PYTHONNOUSERSITE'),"
    "'sys_path': sys.path}))"
)


def view():
    return {"PYTHONPATH": os.environ.get("PYTHONPATH"),
            "PYTHONNOUSERSITE": os.environ.get("PYTHONNOUSERSITE"),
            "sys_path": sys.path}


# No `env=` on this call, deliberately: the grandchild must inherit whatever
# this process actually holds, which is what makes nested inheritance
# observable rather than asserted.
_g = subprocess.run([sys.executable, "-c", _GRANDCHILD],
                    capture_output=True, text=True)

print(MARKER + json.dumps({"child": view(),
                           "grandchild": json.loads(_g.stdout)}))
