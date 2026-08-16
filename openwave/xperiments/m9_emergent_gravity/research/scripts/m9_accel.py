"""Prefer a Nuitka-built helper tree when present.

Source ``.py`` files remain the record. After
``./compile_nuitka.sh``, ``_nuitka/*.so`` is imported
first. Same equations, compiled mpmath loops.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACCEL = os.path.join(HERE, "_nuitka")
if os.path.isdir(ACCEL) and ACCEL not in sys.path:
    sys.path.insert(0, ACCEL)
